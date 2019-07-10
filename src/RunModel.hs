{-# LANGUAGE BangPatterns              #-}
{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE MultiParamTypeClasses     #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TupleSections             #-}
{-# LANGUAGE TypeFamilies              #-}

module RunModel
  ( runModel, X.TDigest, latextable, posteriorMatrices
  ) where

import           Control.Foldl                 (FoldM (..))
import qualified Control.Foldl                 as F
import           Control.Lens
import           Control.Monad                 (when)
import           Data.Foldable                 (fold)
import qualified Data.HashMap.Strict           as M
import           Data.List                     (intercalate, intersperse, nub,
                                                sort)
import           Data.Maybe                    (fromMaybe)
import           Data.Monoid                   ((<>))
import           Data.Reflection               (Reifies)
import           Data.TDigest                  as X
import qualified Data.Text                     as T
import           Data.Traversable              (mapAccumL)
import           Data.Vector                   (Vector)
import qualified Data.Vector                   as V
import           InMatrix
import           Linear.Matrix
import           MarkovChain
import           Matrix
import           Model
import           Numeric.AD
import           Numeric.AD.Internal.Reverse   (Reverse, Tape)
import           Numeric.AD.Mode.Reverse       as Reverse (gradWith')
import           Pipes
import qualified Pipes.Prelude                 as P
import qualified Probability                   as P
import           System.IO                     (BufferMode (..), IOMode (..),
                                                hPutStrLn, hSetBuffering,
                                                stdout, withFile)
import           System.Random.MWC.Probability
import           Text.Printf                   (PrintfArg, printf)



toError :: Either String c -> c
toError = either error id

type HMT = M.HashMap T.Text

runModel
  :: (Maybe (Int, Double))
  -> Int
  -> String
  -> V.Vector Int
  -> Model Double
  -> (forall a. Floating a => V.Vector a -> a)
  -> HMT (ModelParam Double)
  -> IO (HMT (Maybe Double, TDigest 3), M.HashMap (T.Text, T.Text) Double)
runModel hamParams nsamps outfile dataH model' logReg modelparams = do
  let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
      start = _mpInitialValue <$> mps
      priors = _mpPrior <$> mps

      variations :: (Mode c, Scalar c ~ Double) => V.Vector (ModelVar c)
      variations = fmap auto . _mpVariation <$> mps

      model :: forall c. (Mode c, Scalar c ~ Double) => Model c
      model = auto <$> model'

      logLH
        :: forall c. (Ord c, Floating c, Mode c, Scalar c ~ Double)
        => V.Vector c -> c
      logLH =
        toError
        . modelLogPosterior
            dataH
            model
            variations
            (ppToFunc . fmap auto <$> priors)
            logReg

      gLogLH = grad logLH

  hSetBuffering stdout LineBuffering
  putStrLn "finding best-fit starting parameters"
  putStrLn "starting from here:"
  print $ V.zip mpnames start
  putStrLn "\nwith a reco prediction of:"
  print . toError $ prediction =<< appVars variations start model
  putStrLn "\ncompared to the data:"
  print dataH
  putStrLn "\nlog-likelihood of starting params:"
  print $ logLH start
  putStrLn ""


  -- if we want to print the full likelihood
  -- putStrLn "likelihood:"
  -- print . logLH $ fmap (var . T.unpack) mpnames

  let xs' = take (100 * length start) $ gradientAscent' logLH start
      xs = last $ start : xs'

  start' <-
    if any isNaN xs || isNaN (logLH xs)
      then do
        putStrLn "warning: could not find a likelihood maximum"
        putStrLn "based on your model and initial parameter values."
        putStrLn "using initial values to determine hessian:"
        putStrLn "this could be extremely inefficient."
        putStrLn ""
        -- putStrLn "the 'best' signal histogram appears to be"
        -- print $ signalVars dataH model
        return start
      else
        return xs

  putStrLn "best-fit params:"
  print $ V.zip mpnames start'
  putStrLn "\ngradient of log-likelihood given best-fit params:"
  print . V.zip mpnames $ gLogLH start'
  putStrLn "\nlog-likelihood of best-fit params:"
  print $ logLH start'
  putStrLn ""

  -- invert hessian -> covariance matrix
  -- then find the transform from canonical variables
  -- (with mu ~0 and sigma ~1) to the original variables.
  let potential
        :: forall c. (Ord c, Floating c, Mode c, Scalar c ~ Double)
        => V.Vector c -> c
      potential = negate . logLH

      hess' = hessian potential start'
      cov = toError $ invM hess'
      t = cholM cov

      transform'
        :: forall c. (Ord c, Floating c, Mode c, Scalar c ~ Double)
        => V.Vector c -> Vector c
      transform' v = ((fmap auto <$> t) !* v) ^+^ (auto <$> start')

      potential'
        :: forall c. (Ord c, Floating c, Mode c, Scalar c ~ Double)
        => V.Vector c -> c
      potential' = potential . transform'

      predstart = toError $ prediction =<< appVars variations start' model
      rref' = toError $ inM rref hess'

      covlist = V.toList $ V.toList <$> cov
      varlist = V.toList . getDiag . fromLists $ covlist
      absuncerts = zipWith (\x v -> abs . (/sqrt x) <$> v) varlist covlist
      reluncerts = zipWith (\x v -> (/x) <$> v) (V.toList start') absuncerts

      mpnameslist = V.toList mpnames
      matToTable = latextable' mpnameslist

  putStrLn "prediction given starting params:"
  print predstart
  putStrLn ""
  putStrLn "actual data:"
  print dataH
  putStrLn ""

  putStrLn "hessian matrix:"
  print . fst $ toMatrix hess'
  putStrLn ""
  putStrLn "hessian matrix (rref):"
  print . fst $ toMatrix rref'
  putStrLn ""
  putStrLn "covariance matrix:"
  print . fst $ toMatrix cov
  putStrLn ""
  putStr . T.unpack $ matToTable covlist
  putStrLn ""
  putStrLn "absolute uncertainties:"
  putStr . T.unpack $ matToTable absuncerts
  putStrLn ""
  putStrLn "relative uncertainties:"
  putStr . T.unpack $ matToTable reluncerts
  putStrLn ""
  putStrLn "transformation matrix:"
  print . fst $ toMatrix t

  when
    (anyOf (traverse.traverse) isNaN t)
    . error
      $ "we have NaNs in the transformation matrix."
        ++ "\nthis might mean the hessian matrix is singular"
        ++ "\n(e.g. you may need more constraints in your model)."
        ++ "\nexiting."

  -- need an RNG...
  g <- createSystemRandom

  let zeros = const 0 <$> start'

      chain :: PrimMonad m => Producer (V.Vector Double, Double) (Prob m) ()
      chain =
        case hamParams of
          Just (steps, eps) ->
            hamiltonian
              steps
              eps
              potential'
              (grad potential')
              (zeros, potential' zeros)
            >-> P.map (\(v, p) -> (v, negate p))

          Nothing ->
            metropolis
              (traverse (`normal` sig))
              (logLH . transform')
              (zeros, logLH start')

      -- optimal metropolis size taken from
      -- Gelman et al. (1996)
      sig = 2.38 / (sqrt . fromIntegral . length $ start')

  case hamParams of
    Just (steps, eps) -> do
      putStrLn ""
      putStrLn "hamiltonian sampling parameters:"
      putStrLn $ "leapfrog steps: " ++ show steps
      putStrLn $ "epsilon: " ++ show eps

    Nothing -> do
      putStrLn ""
      putStrLn "metropolis step size:"
      print sig

  -- write the walk locations to file.
  withFile outfile WriteMode $ \f -> do
    hSetBuffering f LineBuffering

    let binnames = iover traversed (\i _ -> "recobin" <> T.pack (show i)) predstart
        normnames = ("norm" <>) <$> V.filter (T.isInfixOf "truthbin") mpnames
        names = V.cons "llh" $ mpnames V.++ binnames V.++ normnames

        allParams (ps, llh) =
          let theseparams = transform' ps
              normalize' ts =
                let (s', ts') = mapAccumL (\s u -> (s+u, u/s')) 0 ts
                in ts'
              thispred =
                toError $ prediction =<< appVars variations theseparams model
              normparams =
                normalize' . fmap snd . V.filter (T.isInfixOf "truthbin" . fst)
                $ V.zip mpnames theseparams
          in V.cons llh $ theseparams V.++ thispred V.++ normparams

        showMC ps = mconcat . intersperse ", " . V.toList $ show <$> ps

        toHandle h =
          FoldM
            (\() u -> liftIO . hPutStrLn h $ showMC u)
            (return ())
            return

        vectorize :: Int -> F.Fold a b -> F.Fold (Vector a) (Vector b)
        vectorize n (F.Fold h u r) =
          F.Fold (\vxs -> forceV . V.zipWith h vxs) (V.replicate n u) (fmap r)
          where
            forceV v = foldr seq () v `seq` v

        tdigestf :: F.Fold Double (TDigest 3)
        tdigestf = F.Fold (flip insert) mempty id

        pairwise :: [a] -> [(a, a)]
        pairwise ys@(y:ys') = fmap (y,) ys ++ pairwise ys'
        pairwise _          = []

        vcovf n =
          F.premap (V.fromList . pairwise . V.toList)
          $ vectorize ((n+1)*n) P.covf

        folder =
          F.impurely P.foldM printAndStore
          $ chain >-> P.take nsamps >-> P.map allParams


        printAndStore
          :: MonadIO m
          => FoldM m (Vector Double) (Vector (Maybe Double, TDigest 3), Vector Double)
        printAndStore =
          const
          <$> F.generalize
            ((,)
              <$> vectorize (V.length names) ((,) <$> F.head <*> tdigestf)
              <*> vcovf (V.length names)
            )
          <*> toHandle f

    hPutStrLn f . mconcat . intersperse ", " $ T.unpack <$> V.toList names

    (vval, vcov) <- sample folder g

    let hmval = M.fromList . V.toList $ V.zip names vval
        covnames = V.fromList . pairwise $ V.toList names
        hmcov = symmetrize . M.fromList . V.toList $ V.zip covnames vcov

        symmetrize = M.fromList . go . M.toList
          where
            go []                   = []
            go (((a1, a2), b):rest) = ((a1, a2), b) : ((a2, a1), b) : go rest

    return (hmval, hmcov)




posteriorMatrices
  :: Ord a
  => M.HashMap T.Text (Maybe a, TDigest comp)
  -> M.HashMap (T.Text, T.Text) Double
  -> M.HashMap (T.Text, T.Text) (Double, Double, Double)
posteriorMatrices params covariances =
  let params' =
        fromMaybe (error "error getting mode or quantiles") . sequence
        $ quant <$> params

      quant (mx, y) = do
        x <- mx
        q16 <- quantile 0.16 y
        q50 <- quantile 0.50 y
        q84 <- quantile 0.84 y
        return (x, (q16, q50, q84))

      symmetrize hm =
        let insert1 h (x, y) = M.insert (y, x) (h M.! (x, y)) h
        in foldl insert1 hm $ M.keys hm


      uncerts =
        flip M.mapMaybeWithKey (symmetrize covariances) $ \(name, name') cov ->
          let var =
                fromMaybe (error "missing variance")
                $ M.lookup (name, name) covariances

              var' =
                fromMaybe (error "missing variance")
                $ M.lookup (name', name') covariances

              mean1 =
                fromMaybe (error "missing best fit value") $ do
                  (_, (_, q50, _)) <- M.lookup name params'
                  return q50

              corr = cov / sqrt var / sqrt var'
              absuncert = abs $ cov / sqrt var'
              reluncert = absuncert / mean1


          in
            if T.isPrefixOf "normtruthbin" name
                && not (T.isPrefixOf "truthbin" name')
                && not (T.isPrefixOf "recobin" name')
                && name' /= "llh"
              then Just (absuncert, reluncert, corr)
              else Nothing

  in uncerts



eitherA :: (a -> Bool) -> a -> Either a a
eitherA f x = if f x then Left x else Right x

collapseEither :: Either a a -> a
collapseEither (Left x)  = x
collapseEither (Right x) = x

latextable :: PrintfArg a => M.HashMap (T.Text, T.Text) a -> String
latextable m =
  let ks = M.keys m

      poinames = fmap collapseEither . sort . fmap (eitherA (T.isInfixOf "truthbin")) . nub $ fst <$> ks
      npnames = fmap collapseEither . sort . fmap (eitherA (T.isInfixOf "truthbin")) . nub $ snd <$> ks
      fmtLine npname =
        T.unpack (paramToName npname)
        ++ " & "
        ++ intercalate " & "
            ( printf "%.3f" . (M.!) m . (,npname)
              <$> poinames
            )
        ++ " \\\\"

  in unlines $
    [ "\\begin{tabular}{ l " ++ fold (replicate (length poinames) "| r ") ++ "}"
    , " & " ++ intercalate " & " (T.unpack . paramToName <$> poinames) ++ " \\\\"
    , "\\hline"
    ]
    ++ (fmtLine <$> npnames)
    ++ ["\\end{tabular}"]




-- taken directly from the `ad` package, but now testing for NaNs.

-- | The 'gradientDescent' function performs a multivariate
-- optimization, based on the naive-gradient-descent in the file
-- @stalingrad\/examples\/flow-tests\/pre-saddle-1a.vlad@ from the
-- VLAD compiler Stalingrad sources.  Its output is a stream of
-- increasingly accurate results.  (Modulo the usual caveats.)
--
-- It uses reverse mode automatic differentiation to compute the gradient.
gradientDescent'
  :: (Traversable f, RealFloat a, Ord a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> Reverse s a)
  -> f a
  -> [f a]
gradientDescent' f x0 = go x0 fx0 xgx0 0.1 (0 :: Int)
  where
    (fx0, xgx0) = Reverse.gradWith' (,) f x0
    go x fx xgx !eta !i
      | eta == 0              = [] -- step size is 0
      | fx1 > fx || isNaN fx1 = go x fx xgx (eta/2) 0 -- we stepped too far
      | zeroGrad xgx          = [] -- gradient is 0
      | otherwise             =
        x1 : if i == 10
          then go x1 fx1 xgx1 (eta*2) 0
          else go x1 fx1 xgx1 eta (i+1)

      where
        zeroGrad = all (\(_,g) -> g == 0)
        x1 = fmap (\(xi,gxi) -> xi - eta * gxi) xgx
        (fx1, xgx1) = Reverse.gradWith' (,) f x1

gradientAscent'
  :: (Traversable f, RealFloat a, Ord a)
  => (forall s. Reifies s Tape => f (Reverse s a) -> Reverse s a)
  -> f a
  -> [f a]
gradientAscent' f = gradientDescent' (negate . f)


latextable' :: [T.Text] -> [[Double]] -> T.Text
latextable' names mat =
  let fmtLine npname vals =
        paramToName npname
        <> " & "
        <> T.intercalate " & " (T.pack . printf "%.3f" <$> vals)
        <> " \\\\"

  in T.unlines $
      [ "\\begin{tabular}{ l " <> fold (replicate (length names) "| c ") <> "}"
      , " & " <> T.intercalate " & " (paramToName <$> names) <> " \\\\"
      , "\\hline"
      ]
      ++ zipWith fmtLine names mat
      ++ ["\\end{tabular}"]


paramToName :: T.Text -> T.Text
paramToName s =
  if "normtruthbin" `T.isPrefixOf` s
    then
      T.replace "normtruthbin" "\\sigma_"
      . T.cons '$'
      $ T.snoc s '$'
    else
      T.replace "__1" ""
      . T.replace "v2trk" "trk"
      $ T.replace "_" " " s

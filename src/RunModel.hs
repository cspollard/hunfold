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
  ( runModel, X.TDigest
  ) where

import           Control.DeepSeq
import           Control.Foldl                 (FoldM (..))
import qualified Control.Foldl                 as F
import           Control.Lens
import           Control.Monad                 (when)
import qualified Data.HashMap.Strict           as M
import           Data.List                     (intersperse)
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
import           System.IO                     (BufferMode (..), IOMode (..),
                                                hFlush, hPutStrLn,
                                                hSetBuffering, stdout, withFile)
import           System.Random.MWC.Probability

toError :: Either String c -> c
toError = either error id

type HMT = M.HashMap T.Text

runModel
  :: Int
  -> String
  -> V.Vector Int
  -> Model Double
  -> HMT (ModelParam Double)
  -> IO (HMT (Maybe Double, TDigest 3), M.HashMap (T.Text, T.Text) Double)
runModel nsamps outfile dataH model' modelparams = do
  let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
      start = _mpInitialValue <$> mps
      priors = _mpPrior <$> mps
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

      gLogLH = grad logLH

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
  let hess' = hessian (negate . logLH) start'
      cov = toError $ invM hess'
      t = cholM cov
      it = toError $ invM t
      transform' v = (t !* v) ^+^ start'
      invtransform' v' = it !* (v' ^-^ start')
      predstart = toError $ prediction =<< appVars variations start' model
      rref' = toError $ inM rref hess'

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

  let chain :: PrimMonad m => Producer (V.Vector Double, Double) (Prob m) ()
      chain =
        metropolis
          (traverse (`normal` sig))
          (logLH . transform')
          (invtransform' start', logLH start')

      -- optimal metropolis size taken from
      -- Gelman et al. (1996)
      sig = 2.38 / (sqrt . fromIntegral . length $ start')

  putStrLn ""
  putStrLn "metropolis step size:"
  print sig
  hFlush stdout

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
        vectorize n fol =
          case fol of
            (F.Fold h u r) -> F.Fold (\vxs -> forceV . V.zipWith h vxs) (V.replicate n u) (fmap r)
          where
            forceV v = foldr seq () v `seq` v

        tdigestf :: F.Fold Double (TDigest 3)
        tdigestf = F.Fold (flip insert) mempty id

        -- fold for calculating the covariance between two parameters
        covf :: F.Fold (Double, Double) Double
        covf = F.Fold step (0, 0, 0, 0 :: Int) done
          where
            step (c, meanx, meany, n) (x, y) =
              let dx = x - meanx
                  dy = y - meany
                  n' = n+1
                  nd = fromIntegral n'
                  meanx' = meanx + dx/nd
                  meany' = meany + dy/nd
                  c' = c + dx*dy
              in force (c', meanx', meany', n')

            done (c, _, _, n)
              | n >= 2 = c / fromIntegral (n-1)
              | otherwise = 1.0/0.0

        pairwise :: [a] -> [(a, a)]
        pairwise (x:ys) = fmap (x,) ys ++ pairwise ys
        pairwise _      = []

        vcovf n =
          F.premap (V.fromList . pairwise . V.toList) $ vectorize (n*(n-1)) covf

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
        hmcov = M.fromList . V.toList $ V.zip covnames vcov

    return (hmval, hmcov)


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

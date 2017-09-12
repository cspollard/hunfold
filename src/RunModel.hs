{-# LANGUAGE BangPatterns              #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE MultiParamTypeClasses     #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TypeFamilies              #-}

module RunModel where

import           Control.Lens
import           Control.Monad                 (when)
import qualified Data.HashMap.Strict           as M
import           Data.List                     (intersperse)
import           Data.Monoid                   ((<>))
import           Data.Reflection               (Reifies)
import qualified Data.Text                     as T
import qualified Data.Vector                   as V
import           InMatrix                      hiding (trace, transpose, zero)
import           Linear.Matrix                 hiding (trace)
import           MarkovChain
import           Matrix                        hiding (trace)
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

runModel
  :: Int
  -> String
  -> V.Vector Int
  -> Model Double
  -> M.HashMap T.Text (ModelParam Double)
  -> IO ()
runModel nsamps outfile dataH model' modelparams = do
  let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
      start = _mpInitialValue <$> mps
      priors = _mpPrior <$> mps
      variations = fmap auto . _mpVariation <$> mps
      model :: forall c. (Mode c, Scalar c ~ Double) => Model c
      model = auto <$> model'

      logLH
        :: forall c. (Floating c, Mode c, Scalar c ~ Double)
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
  putStrLn "starting from here..."
  print mpnames
  print start


  -- if we want to print the full likelihood
  -- putStrLn "likelihood:"
  -- print . logLH $ fmap (var . T.unpack) mpnames

  let xs = take (10 * length start) $ gradientAscent' logLH start
      x = last xs

  start' <-
    if any isNaN x || isNaN (logLH x)
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
        return x

  putStrLn "starting params:"
  print start'
  putStrLn ""
  putStrLn "gradient of log-likelihood given starting params:"
  print $ gLogLH start'
  putStrLn ""
  putStrLn "log-likelihood of starting params:"
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
        showMC (ps, llh) =
          let theseparams = transform' ps
              thispred =
                toError $ prediction =<< appVars variations theseparams model
          in
            mconcat . intersperse ", " . V.toList
            $ show <$> V.cons llh theseparams V.++ thispred

    hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
      $ "llh" : V.toList mpnames ++ V.toList binnames

    let effect :: (PrimMonad m, MonadIO m) => Effect (Prob m) ()
        effect = chain >-> P.take nsamps >-> P.map showMC >-> P.toHandle f

    sample (runEffect effect) g


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

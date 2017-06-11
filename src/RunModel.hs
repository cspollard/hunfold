{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}

module RunModel where

import           Control.Lens
import           Control.Monad                 (when)
import           Data.List                     (intersperse)
import qualified Data.Map.Strict               as M
import           Data.Monoid                   ((<>))
import qualified Data.Text                     as T
import qualified Data.Vector                   as V
import           InMatrix                      hiding (transpose, zero)
import           Linear.Matrix
import           MarkovChain
import           Matrix
import           Model
import           NDSL
import           Numeric.AD
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
  -> M.Map T.Text (ModelParam Double)
  -> IO ()
runModel nsamps outfile dataH model modelparams = do
  let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
      start = fmap _mpInitialValue mps
      priors = fmap _mpPrior mps
      variations = fmap _mpVariation mps

      logLH
        :: forall c. (Floating c, Mode c, Scalar c ~ Double)
        => V.Vector c -> c
      logLH =
        toError
        . modelLogPosterior
            dataH
            (fmap auto model)
            (fmap (fmap auto) variations)
            (fmap (ppToFunc . fmap auto) priors)

      gLogLH = grad logLH

  putStrLn ""

  putStrLn "likelihood:"
  print . logLH $ fmap (var . T.unpack) mpnames

  -- TODO
  -- what is best value?!
  -- find the maximum likelihood starting location
  let xs = take (length start) $ conjugateGradientAscent logLH start
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
  print mpnames
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

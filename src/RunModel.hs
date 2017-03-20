{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}

module RunModel where

import           Control.Lens
import           Control.Monad    (when)
import           Data.List        (intersperse)
import qualified Data.Map.Strict  as M
import           Data.Monoid      ((<>))
import qualified Data.Text        as T
import qualified Data.Vector      as V
import           InMatrix         hiding (transpose, zero)
import           Linear.Matrix
import qualified List.Transformer as LT
import           MarkovChain
import           Matrix
import           Model
import           Numeric.AD
import           Numeric.MCMC
import           System.IO        (BufferMode (..), IOMode (..), hFlush,
                                   hPutStr, hPutStrLn, hSetBuffering, stdout,
                                   withFile)

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

      -- I'm not sure why we need an explicit type here.
      -- probably because of the RankNType going on here
      logLH
        :: forall c. (Floating c, Ord c, Mode c, Scalar c ~ Double)
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

  -- finally, build the chain, metropolis transition, and the MCMC walk
  let c =
        Chain
          (Target (logLH . transform') $ Just (gLogLH . transform'))
          (logLH start')
          (invtransform' start')
          Nothing

      -- optimal metropolis size taken from
      -- Gelman et al. (1996)
      sig = 2.38 / (sqrt . fromIntegral . length $ start')
      trans = metropolis sig

      walk = runMC trans c g

  putStrLn ""
  putStrLn "metropolis step size:"
  print sig
  hFlush stdout

  -- write the walk locations to file.
  withFile outfile WriteMode $ \f -> do
    hSetBuffering f LineBuffering
    let binnames = iover traversed (\i _ -> "recobin" <> T.pack (show i)) predstart

    hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
      $ "llh" : V.toList mpnames ++ V.toList binnames

    LT.runListT . LT.take nsamps $ do
      Chain{..} <- walk
      LT.liftIO $ do
        hPutStr f $ show chainScore ++ ", "
        -- TODO
        -- this is very wasteful: we are already getting prediction
        -- to calculate the log likelihood...
        let theseparams = transform' chainPosition
            thispred =
              toError $ prediction =<< appVars variations theseparams model
        hPutStrLn f
          . mconcat . intersperse ", " . V.toList
          $ show <$> (theseparams V.++ thispred)

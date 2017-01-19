{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}
{-# LANGUAGE FlexibleContexts             #-}

module Dagostini where

import           Control.Monad.Trans.Class     (lift)
import qualified Data.Vector.Fixed.Cont        as V
import           Likelihood
import qualified List.Transformer              as LT
import           Matrix
import           System.IO                     (IOMode (..), hPutStrLn,
                                                openFile)
import           System.Random.MWC.Probability


logPostLH :: Num a => (b -> a) -> (b -> a) -> b -> a
logPostLH logLikelhood logPrior sigmas = logLikelhood sigmas + logPrior sigmas

logPostLHNIter :: Num a => Int -> (b -> a) -> (b -> a) -> b -> a
-- TODO
-- which is better here?
-- logPostLHNIter n logCondProb = head . drop n . iterate (logPostLH logCondProb)
logPostLHNIter n logLikelhood logPrior sigmas = fromIntegral n * logLikelhood sigmas + logPrior sigmas


-- test case:
-- 3 effect bins
-- 2 cause bins

type NE = ToPeano 3
type NC = ToPeano 2

test :: Int -> Int -> Int -> IO ()
test niters nextsamps nintsamps = do
  f <- openFile "test.dat" WriteMode
  hPutStrLn f "bin1, bin2"
  _ <- withSystemRandom . samples nextsamps
      $ LT.runListT . LT.take nintsamps
      $ do
        [x1, x2] <- V.toList <$> sampledPosterior
        LT.liftIO . hPutStrLn f
          $ show x1 ++ ", " ++ show x2

  return ()

  where
    logLikelihood dat bkgs smear lumi sigmas =
      -- xs: the prediction based on the given cross sections
      -- NB: the first row of the migration matrix should be the inefficiency of
      --     each cause bin.
      let xs = fmap (*lumi) $ (+) <$> tailV (multMV smear sigmas) <*> bkgs
      in sum $ logPoissonP <$> dat <*> xs

    -- the probability distribution of log posterior functions given the data
    -- and distribution of "worlds"
    logPosteriors :: Prob IO (Vec NC Double -> Double)
    logPosteriors = do
      lumi <- cutOffNormal 1 0.3
      bkgs <- myBkgs
      smears <- mySmears
      let llhf = logLikelihood myData bkgs smears lumi
      -- NB: flat spectrum prior
      return $ logPostLHNIter niters llhf (const 0)

    sampledPosterior :: ListT (Prob IO) (Vec NC Double)
    sampledPosterior = fmap fst $ do
      logPost <- lift logPosteriors
      start <- lift myStart
      g <- LT.liftIO createSystemRandom
      runLLH (metropolis 0.1) logPost start g


cutOffNormal :: PrimMonad m => Double -> Double -> Prob m Double
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x


myData :: Vec NE Int
myData = [0, 1, 0]

myBkgs :: PrimMonad m => Prob m (Vec NE Double)
myBkgs = do
  norm <- normal 1 0.2
  return . fmap (*norm) $ [0.5, 0.0, 0.5]

mySmears :: PrimMonad m => Prob m (Mat NC (S NE) Double)
mySmears = do
  let effs1 = [0.0, 0.5, 0.0, 0.5]
  let effs2 = [0.0, 0.0, 1.0, 0.0]
  return . transpose
    $ [effs1, effs2]

myStart :: PrimMonad m => Prob m (Vec NC Double)
myStart = sequence [uniform, uniform]


{-

-- Ci: cause bin i
-- si: events in cause bin i due to process in question (what we are measuring)
-- sij: events in cause bin i and effect bin j
-- Ej: effect bin j
-- xj: events in effect bin j (what we observe)

-- P(Ej|Ci):
--   probability that an event lands in effect bin j given an event in cause
--   bin i, ie. the smearing matrix
-- P(Ci|Ej):
--   probability that an event lands in cause bin i given an event in effect
--   bin j

-- Bayes' theorem
-- P(x) = P(x|y) P(y) / P(y|x)
-- ie. P(Ci) = P(Ci|Ej) P(Ej) / P(Ej|Ci)

-- so P(Ci|Ej) = P(Ej|Ci) P(Ci) / P(Ej)
--   where
--     P(Ej) = sum_i P(Ej|Ci) P(Ci)

-- so sij = P(Ci|Ej) xj
-- si = (sum_j P(Ci|Ej) xj) / (sum_j P(Ej|Ci))
--   where sum_j P(Ej|Ci) is the observation efficiency of cause bin i

-- we observe a spectrum xj (integers) in n effect bins.
-- we hope to measure differential cross sections si in m cause bins.
-- we have m*(n+1) P(Ej|Ci) terms corresponding to the smearing matrix.
--   (n+1 takes into account the inefficiency)
-- we need to marginalize over uncertainties in our backgrounds.
-- we need to marginalize over uncertainties in our signal model, P(Ej|Ci).
-- we need to deal with the prior distribution P(Ci) (iterative?)
-- given a particular P(Ej|Ci), a background prediction bj and the obvervations xj,
--   we can obtain P(xj|si)

-- we want P(si|xj) = P(xj|si) P(si) / P(xj)
--   - where P(xj) = bj + sum_i P(Ej|Ci) P(si)

-- P(sigma)

-}

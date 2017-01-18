{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}

module Dagostini where

import           Control.Monad.Trans.Class     (lift)
import qualified Data.Vector.Fixed.Cont        as V
import           Likelihood
import qualified List.Transformer              as LT
import           Matrix
import           System.Random.MWC.Probability


logPostLH :: Num a => (b -> a) -> (b -> a) -> b -> a
logPostLH logCondProb logPrior sigmas = logCondProb sigmas + logPrior sigmas

logPostLHNIter :: Num a => Int -> (b -> a) -> (b -> a) -> b -> a
logPostLHNIter n logCondProb = head . drop n . iterate (logPostLH logCondProb)


-- test case:
-- 3 data bins
-- 2 cause bins

type NE = ToPeano 3
type NC = ToPeano 2

test :: Int -> Int -> Int -> IO ()
test niters nintsamps nextsamps = do
    writeFile "test.dat" "bin1, bin2\n"
    _ <- withSystemRandom . samples nextsamps
        $ LT.runListT . LT.take nintsamps $
        do
          [x1, x2] <- V.toList <$> chain
          LT.liftIO . appendFile "test.dat"
            $ show x1 ++ ", " ++ show x2 ++ "\n"

    return ()



  where
    -- NB: the first row of the migration matrix should be the inefficiency of
    --     each cause bin.
    logProb dat bkgs smear lumi sigmas =
      -- xs: the prediction based on the given cross sections
      let xs = fmap (*lumi) $ (+) <$> tailV (multMV smear sigmas) <*> bkgs
      in sum $ logPoissonP <$> dat <*> xs

    -- the "bank" of (logs of) posterior distributions iterated n times
    logPosts :: Prob IO (Vec NC Double -> Double)
    logPosts = do
      logCondProb <- logCondProbs
      return $ iterate logCondProb myLogPrior !! niters

    -- the "bank" of (logs of) conditional probabilities functions to marginalize over "the world"
    logCondProbs :: Prob IO ((Vec NC Double -> Double) -> Vec NC Double -> Double)
    logCondProbs = do
      lumi <- cutOffNormal 1 0.2
      return $ logPostLH (logProb myData myBkgs mySmears lumi)

    cutOffNormal mu s = do
      x <- normal mu s
      if x < 0 then cutOffNormal mu s else return x

    chain :: ListT (Prob IO) (Vec NC Double)
    chain = fmap fst $ do
      logPost <- lift logPosts
      g <- LT.liftIO createSystemRandom
      runLLH (metropolis 0.1) logPost myStart g

myData :: Vec NE Int
myData = V.fromList' [1, 0, 1]

myBkgs :: Fractional a => Vec NE a
myBkgs = V.fromList' [0.5, 0.0, 0.5]

mySmears :: Fractional a => Mat NC (S NE) a
mySmears = V.fromList'
  [ V.fromList' [0.0, 0.0]
  , V.fromList' [0.5, 0.0]
  , V.fromList' [0.0, 1.0]
  , V.fromList' [0.5, 0.0]
  ]

myLogPrior :: Num a => Vec NC a -> a
myLogPrior = const 0

myStart :: Num a => ContVec NC a
myStart = V.fromList' [1, 1]

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

{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}

module Dagostini where

import qualified Data.Vector.Fixed.Cont as V
import           Likelihood
import qualified List.Transformer       as LT
import           Matrix



postLLH :: (Integral b, Floating a, Arity neffect, Arity ncause)
  => Vec neffect b
  -> Vec neffect a
  -> Mat ncause (S neffect) a
  -> a
  -> (Vec ncause a -> a)
  -> Vec ncause a -> a
postLLH dat bkgs smear lumi logPrior sigmas =
  conditionalLogProb + logPrior sigmas
  where
    -- NB: the first row of the migration matrix should be the inefficiency of
    --     each cause bin.
    -- xs: the prediction based on the given cross sections
    xs = fmap (*lumi) $ (+) <$> tailV (multMV smear sigmas) <*> bkgs

    -- the log conditional probability of the data bins given the sigmas
    conditionalLogProb = sum $ logPoissonP <$> dat <*> xs


-- test case:
-- 3 data bins
-- 2 cause bins

type NE = ToPeano 3
type NC = ToPeano 2

-- TODO
-- odd numbers are broken for some reason
test :: Int -> IO ()
test n = LT.runListT $ LT.liftIO . print =<< c

  where
    f :: Floating a => Vec NC a -> a
    f = iterate (postLLH myData myBkgs mySmears 1.0) myLogPrior !! n
    c =
      LT.take 100
        $ runLLH (metropolis 1.0) f (V.fromList' [1, 1])
          =<< LT.liftIO createSystemRandom

myData :: Vec NE Int
myData = V.fromList' [1, 2, 1]

myBkgs :: Fractional a => Vec NE a
myBkgs = V.fromList' [0.0, 0.0, 0.0]

mySmears :: Fractional a => Mat NC (S NE) a
mySmears = V.fromList'
  [ V.fromList' [0.0, 0.0]
  , V.fromList' [0.5, 0.0]
  , V.fromList' [0.0, 1.0]
  , V.fromList' [0.5, 0.0]
  ]

myLogPrior :: Num a => Vec NC a -> a
myLogPrior = const 0

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

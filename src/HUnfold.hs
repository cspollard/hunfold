{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE ScopedTypeVariables       #-}

module HUnfold where

import           Data.List              (intersperse)
import qualified Data.Vector.Fixed.Cont as V
import           Likelihood
import qualified List.Transformer       as LT
import           Matrix
import           Probability
import           System.IO              (IOMode (..), hClose, hPutStrLn,
                                         openFile)



logPostLH :: Num a => (b -> a) -> (b -> a) -> b -> a
logPostLH logLikelhood logPrior sigmas = logLikelhood sigmas + logPrior sigmas

-- test case

type NE = ToPeano 14
type NC = ToPeano 10

test :: Int -> Int -> IO ()
test nburn nsamps = do
  f <- openFile "test.dat" WriteMode

  let n = arity (undefined :: NC) - 1
  hPutStrLn f
    $ mconcat . intersperse ", bin" $ "bin0" : fmap show [1..n]

  _ <- withSystemRandom . asGenIO
      $ \g -> LT.runListT . LT.take nsamps
        $ do
          (xs :: T (ContVec NC Double) Double) <- gen g
          LT.liftIO . hPutStrLn f
            $ mconcat . intersperse ", " . V.toList $ show <$> fstT xs

  hClose f

  where
    logLikelihood dat bkgs smear lumi sigmas =
      -- xs: the prediction based on the given cross sections
      -- NB: the first row of the migration matrix should be the inefficiency of
      --     each cause bin.
      let xs = fmap (*lumi) $ (+) <$> tailV (multMV smear sigmas) <*> bkgs
      in sum $ logPoissonP <$> dat <*> xs

    logPost :: Floating a => ContVec NC a -> a
    logPost = logLikelihood myData myBkgs mySmears 1.0

    nonNegLogPrior :: (Foldable t1, Ord a, Num a, Fractional t) => t1 a -> t
    nonNegLogPrior xs = if any (< 0) xs then neginf else 0
      where neginf = negate $ 1/0

    gen :: (InvErf b, Variate b, PrimMonad m, RealFloat b)
        => Gen (PrimState m) -> ListT m (T (ContVec NC b) b)
    gen g =
      do
        let start = pure $ (fromIntegral (sum myData) - sum myBkgs) / fromIntegral nC
        -- let start = findMaximum logPost initial
        LT.drop nburn
          $ runMC (prop (logPostLH nonNegLogPrior logPost)) (T start $ logPost start) g

    prop = weightedProposal (metropolis 100)


cutOffNormal :: (InvErf b, Variate b, PrimMonad m, Ord b) => b -> b -> Prob m b
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x

myBkgs :: Num a => Vec NE a
myBkgs = pure 0

mySmears :: Fractional a => Mat NC (S NE) a
mySmears = dagSmear2

{-
mySmears = transpose
  [ [0.25, 0.5, 0.25]
  , [0.1, 0.2, 0.7]
  ]
-}

nE :: Int
nE = arity (undefined :: NE)

nC :: Int
nC = arity (undefined :: NC)

myData :: Vec NE Int
-- myData = [1, 3]
-- myData = fmap floor . tailV . multMV dagSmear2 . fmap (*150). fmap dagFun1 $ [1..10]
myData = fmap floor . tailV . multMV dagSmear2 . fmap dagFun2 $ [(1 :: Int)..10]

dagSmear2 :: Fractional a => Mat NC (S NE) a
dagSmear2 = transpose
  [ [0.1, 0.045, 0.09, 0.945, 0, 0, 0, 0, 0, 0, 0, 0.18, 0.27, 0.27, 0.09]
  , [0.1, 0, 0.054, 0.072, 0.054, 0, 0, 0, 0.045, 0.135, 0.36, 0.18, 0, 0, 0]
  , [0.1, 0, 0, 0.045, 0.045, 0, 0, 0.135, 0.27, 0.27, 0.45, 0, 0, 0.045, 0.045]
  , [0.1, 0, 0, 0, 0.45,0.18, 0.315, 0.27, 0.045, 0, 0, 0, 0.045, 0, 0]
  , [0.1, 0, 0, 0, 0.225, 0.36, 0.18, 0.045, 0, 0.045, 0.045, 0, 0, 0, 0]
  , [0.1, 0, 0, 0, 0.18, 0.45, 0.18, 0, 0, 0.045, 0.045, 0, 0, 0, 0]
  , [0.1, 0, 0, 0, 0.36, 0.315, 0.225, 0, 0.045, 0, 0, 0, 0, 0, 0]
  , [0.1, 0, 0, 0.225, 0.405, 0.225, 0, 0.045, 0, 0, 0, 0, 0, 0, 0]
  , [0.1, 0, 0.225, 0.36, 0.225, 0, 0.045, 0, 0.045, 0, 0, 0, 0, 0, 0]
  , [0.1, 0.225, 0.405, 0.225, 0, 0, 0, 0, 0, 0.045, 0, 0, 0, 0, 0]
  ]

dagFun1 :: Integral a => a -> Double
dagFun1 x = (*) 150 $ 11 - fromIntegral x

dagFun2 :: Integral a => a -> Double
dagFun2 x = (*) 25 $ fromIntegral $ x*x

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

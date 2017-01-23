{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE ScopedTypeVariables       #-}

module HUnfold where

import           Data.List              (intersperse)
import qualified Data.Vector.Fixed.Cont as V
import qualified List.Transformer       as LT
import           System.IO              (IOMode (..), hClose, hPutStrLn,
                                         openFile)

import           MarkovChain
import           Matrix
import           Metropolis
import           Probability



logPostLH :: Num a => (b -> a) -> (b -> a) -> b -> a
logPostLH logLikelhood logPrior sigmas = logLikelhood sigmas + logPrior sigmas

-- test case

type NE = ToPeano 4
type NC = ToPeano 5

test :: Int -> Int -> IO ()
test nburn nsamps = do
  f <- openFile "test.dat" WriteMode

  hPutStrLn f
    $ mconcat . intersperse ", bin" $ "llh" : fmap show [0..(nC-1)]

  _ <- withSystemRandom . asGenIO
      $ \g -> LT.runListT . LT.take nsamps
        $ do
          (T xs llh :: T (ContVec NC Double) Double) <- gen g
          let ys = V.cons llh xs
          LT.liftIO . hPutStrLn f
            $ mconcat . intersperse ", " . V.toList $ show <$> ys
  hClose f

  where
    logLikelihood dat bkgs smear lumi sigmas =
      -- xs: the prediction based on the given cross sections
      -- NB: the first row of the migration matrix should be the inefficiency of
      --     each cause bin.
      let xs = fmap (*lumi) $ (+) <$> multMV smear sigmas <*> bkgs
      in sum $ logPoissonP <$> dat <*> xs

    logPost :: Floating a => ContVec NC a -> a
    logPost = logLikelihood zeeData zeeBkgs zeeSmear zeeLumi

    nonNegLogPrior :: (Foldable t1, Ord a, Num a, Fractional t) => t1 a -> t
    nonNegLogPrior xs = if any (< 0) xs then neginf else 0
      where neginf = negate $ 1/0

    gen :: (InvErf b, Variate b, PrimMonad m, RealFloat b)
        => Gen (PrimState m) -> ListT m (T (ContVec NC b) b)
    gen g =
      do
        let start = pure $ (fromIntegral (sum zeeData) - sum zeeBkgs * zeeLumi) / fromIntegral (nC * zeeLumi)
        LT.drop nburn
          $ runMC (prop (logPostLH nonNegLogPrior logPost)) (T start $ logPost start) g

    prop = weightedProposal . dynamicMetropolis $ exp <$> uniformR (-7, -2)


cutOffNormal :: (InvErf b, Variate b, PrimMonad m, Ord b) => b -> b -> Prob m b
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x

zeeBkgs :: Fractional a => Vec NE a
zeeBkgs = [0.32, 0.12, 0.015, 0.0026]

mySmears :: Fractional a => Mat NC NE a
mySmears = zeeSmear

zeeLumi :: Num t => t
zeeLumi = 34000

zeeSmear :: Fractional a => Mat NC NE a
zeeSmear = transpose
  [ [0.05, 0.00, 0.0, 0.0]
  , [0.84, 0.16, 0.0, 0.0]
  , [0.07, 0.87, 0.06, 0.0]
  , [0.0, 0.08, 0.84, 0.08]
  , [0.0, 0.0, 0.05, 0.95]
  ]

zeeData :: Vec NE Int
zeeData = floor . (* 34) . (/ 3.2) <$> [6766, 2309, 288, 54 :: Double]

nE :: Int
nE = arity (undefined :: NE)

nC :: Int
nC = arity (undefined :: NC)

{-
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
-}

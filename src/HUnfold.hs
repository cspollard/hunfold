{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TemplateHaskell           #-}
{-# LANGUAGE TypeFamilies              #-}


module HUnfold where

import           Control.Applicative (liftA2)
import           Control.Lens
import           Data.List           (intersperse)
import           Data.Map.Strict     (Map)
import qualified Data.Map.Strict     as M
import           Data.Text           (Text)
import qualified Data.Text           as T
import qualified List.Transformer    as LT
import           Numeric.AD
import           System.IO           (IOMode (..), hPutStr, hPutStrLn, withFile)

import           Hamiltonian
import           MarkovChain
import           Matrix
import           Probability



-- test case

type NE = ToPeano 4
type NC = ToPeano 5


cutOffNormal :: (InvErf b, Variate b, PrimMonad m, Ord b) => b -> b -> Prob m b
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x

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

-- what do I need from model?
--

type Hist = Vec
type Param n m a = (a -> (Model n m a -> Model n m a, a))


data Model n m a =
  Model
    { _mBkgs   :: Map Text (Hist m a)
    , _mSigs   :: Hist n a
    , _mSmears :: Mat n m a
    , _mLumi   :: a
    }

makeLenses ''Model


myModel :: Fractional a => Model NC NE a
myModel =
  Model
    (M.singleton "ttbar" [0.32, 0.12, 0.015, 0.0026])
    (pure 1)
    zeeSmear
    34000


myModelParams :: Floating a => Map Text (Param NC NE a)
myModelParams = M.fromList
  [ ("ttbarnorm", \x -> (over (mBkgs.ix "ttbar") (fmap (*x)), logLogNormalP 1 0.2 x))
  , ("sigma0", \x -> (set (mSigs.element 0) x, 0))
  , ("sigma1", \x -> (set (mSigs.element 1) x, 0))
  , ("sigma2", \x -> (set (mSigs.element 2) x, 0))
  , ("sigma3", \x -> (set (mSigs.element 3) x, 0))
  , ("sigma4", \x -> (set (mSigs.element 4) x, 0))
  , ("lumi", \x -> (over mLumi (*x), logLogNormalP 1 0.2 x))
  ]


myInitialParams :: Fractional a => Map Text a
myInitialParams = M.fromList
  [ ("ttbarnorm", 1)
  , ("sigma0", 1)
  , ("sigma1", 1)
  , ("sigma2", 0.1)
  , ("sigma3", 0.01)
  , ("sigma4", 0.001)
  , ("lumi", 1)
  ]


modelPred :: (Arity n, Arity m, Num a) => Model n m a -> Hist m a
modelPred (Model bkgs sigs smears lumi) =
  let bkgTot = foldl (liftA2 (+)) (pure 0) bkgs
  in fmap (*lumi) $ (+) <$> multMV smears sigs <*> bkgTot


appParams :: Num a
          => Map Text (Param n m a)
          -> Model n m a
          -> Map Text a
          -> (Model n m a, a)
appParams fs m ps =
  let fs' = M.intersectionWith ($) fs ps
  in foldr (\(f, p) (model, prior) -> (f model, prior+p)) (m, 0) fs'


modelLogPosterior :: (Arity m, Arity n, Integral a, Floating b)
                  => Map Text (Param n m b)
                  -> Model n m b
                  -> Hist m a
                  -> Map Text b
                  -> b
modelLogPosterior fs model mdata ps =
  let (model', prior) = appParams fs model ps
      preds = modelPred model'
      logLH = sum $ logPoissonP <$> mdata <*> preds
  in logLH + prior


test :: Int -> Int -> IO ()
test nburn nsamps = do

  let llh = modelLogPosterior myModelParams myModel zeeData
  -- TODO
  -- remove test
  -- let llh m = let x = m M.! "sigma0" in logLogNormalP 0 1 x
  let prop = hamiltonian 3 0.0001 llh (grad llh)
              -- weightedProposal $ metropolis 0.01
              -- $ dynamicMetropolis $ exp <$> uniformR (-7, -3)

  g <- createSystemRandom

  let l = LT.drop nburn $ runMC prop (T myInitialParams $ llh myInitialParams) g

  withFile "test.dat" WriteMode
    $ \f -> do
      hPutStrLn f
        $ mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys myInitialParams

      LT.runListT . LT.take nsamps
        $ do
          (T xs llh :: T (Map Text Double) Double) <- l
          LT.liftIO . hPutStr f $ show llh ++ ", "
          LT.liftIO . hPutStrLn f
            $ mconcat . intersperse ", " . M.elems $ show <$> xs

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

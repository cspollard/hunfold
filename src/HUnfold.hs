{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TemplateHaskell           #-}


module HUnfold where

import           Control.Applicative (liftA2)
import           Control.Lens
import           Data.List           (intersperse)
import           Data.Map.Strict     (Map)
import qualified Data.Map.Strict     as M
import           Data.Text           (Text)
import qualified Data.Text           as T
import qualified List.Transformer    as LT
import           System.IO           (IOMode (..), hClose, hPutStr, hPutStrLn,
                                      openFile)

import           MarkovChain
import           Matrix
import           Metropolis
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
type Param n m = (Double -> (Model n m -> Model n m, Double))

data Model n m =
  Model
    { _mBkgs   :: Map Text (Hist m Double)
    , _mSigs   :: Hist n Double
    , _mSmears :: Mat n m Double
    , _mLumi   :: Double
    }

makeLenses ''Model

myModel :: Model NC NE
myModel =
  Model
    (M.singleton "ttbar" [0.32, 0.12, 0.015, 0.0026])
    (pure 1)
    zeeSmear
    34000

myModelParams :: Map Text (Param NC NE)
myModelParams = M.fromList
  [ ("ttbarnorm", \x -> (over (mBkgs.ix "ttbar") (fmap (*(1+x))), logNormalP 0 0.2 x))
  , ("sigma0", \x -> (set (mSigs.element 0) x, 0))
  , ("sigma1", \x -> (set (mSigs.element 1) x, 0))
  , ("sigma2", \x -> (set (mSigs.element 2) x, 0))
  , ("sigma3", \x -> (set (mSigs.element 3) x, 0))
  , ("sigma4", \x -> (set (mSigs.element 4) x, 0))
  , ("lumi", \x -> (over mLumi (*(1+x)), logNormalP 0 0.2 x))
  ]

myInitialParams :: Map Text Double
myInitialParams = M.fromList
  [ ("ttbarnorm", 0.0)
  , ("sigma0", 1.0)
  , ("sigma1", 1.0)
  , ("sigma2", 1.0)
  , ("sigma3", 1.0)
  , ("sigma4", 1.0)
  , ("lumi", 0.0)
  ]



modelPred :: (Arity n, Arity m) => Model n m -> Hist m Double
modelPred (Model bkgs sigs smears lumi) =
  let bkgTot = foldl (liftA2 (+)) (pure 0) bkgs
  in fmap (*lumi) $ (+) <$> multMV smears sigs <*> bkgTot


appParams :: Map Text (Param n m)
          -> Model n m
          -> Map Text Double
          -> (Model n m, Double)
appParams fs m ps =
  let fs' = M.intersectionWith ($) fs ps
  in foldr (\(f, p) (model, prior) -> (f model, prior+p)) (m, 0) fs'


modelLogPosterior :: (Arity m, Arity n, Integral a)
                  => Map Text (Param n m)
                  -> Model n m
                  -> Hist m a
                  -> Map Text Double
                  -> Double
modelLogPosterior fs model mdata ps =
  let (model', prior) = appParams fs model ps
      preds = modelPred model'
      logLH = sum $ logPoissonP <$> mdata <*> preds
  in logLH + prior


-- hamiltonian MCMC?

test :: Int -> Int -> IO ()
test nburn nsamps = do
  f <- openFile "test.dat" WriteMode

  hPutStrLn f
    $ mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys myInitialParams

  _ <- withSystemRandom . asGenIO
      $ \g -> LT.runListT . LT.take nsamps
        $ do
          (T xs llh) <- gen g
          LT.liftIO . hPutStr f $ show llh ++ ", "
          LT.liftIO . hPutStrLn f
            $ mconcat . intersperse ", " . M.elems $ show <$> xs
  hClose f

  where
    gen g =
      do
        let f = modelLogPosterior myModelParams myModel zeeData
        let prop = weightedProposal . dynamicMetropolis
                    $ exp <$> uniformR (-7, -3)
        LT.drop nburn
          $ runMC (prop f) (T myInitialParams $ f myInitialParams) g


{-

type Model f n a = f a -> (Hist n a, a)

modelLLH :: (Foldable t, Applicative t, Integral c, Floating b) => (a -> (t b, b)) -> t c -> a -> b
modelLLH model mdata ps =
  let (mpred, mprior) = model ps
  in mprior + sum (logPoissonP <$> mdata <*> mpred)


addBkg :: (Arity n, Num a) => Model f n a -> Model f n a -> Model f n a
addBkg bkgModel model fx =
  let (bkgHist, bkgPrior) = bkgModel fx
      (totHist, totPrior) = model fx
  in ((+) <$> bkgHist <*> totHist, bkgPrior + totPrior)

ttbarBkg :: (Arity n, Floating a) => Model (Map Text) n a
ttbarBkg ps =
  let norm = ps M.! "ttbarNorm"
      prior = logNormalP 1 0.2 norm
  in ((*norm) <$> , prior)

signal :: (Arity n, Num a) => Model (Map Text) n a
signal ps =
  let h = (ps M.!) <$> ["sigma1", "sigma2", "sigma3", "sigma4"]
  in (h, 1)
-}



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

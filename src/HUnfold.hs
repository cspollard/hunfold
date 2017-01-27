{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TemplateHaskell           #-}
{-# LANGUAGE TypeFamilies              #-}


module HUnfold where

import           Control.Applicative (liftA2)
import           Control.Arrow       ((&&&))
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
import           Metropolis
import           Probability



-- test case

type NE = ToPeano 4
type NC = ToPeano 4


cutOffNormal :: (InvErf b, Variate b, PrimMonad m, Ord b) => b -> b -> Prob m b
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x

zeeSmear :: Fractional a => Mat NC NE a
zeeSmear = transpose
  [ [0.84, 0.16, 0.0, 0.0]
  , [0.07, 0.87, 0.06, 0.0]
  , [0.0, 0.08, 0.84, 0.08]
  , [0.0, 0.0, 0.05, 0.95]
  ]

zeeData :: Vec NE Int
zeeData = [6766, 2309, 288, 54]

nE :: Int
nE = arity (undefined :: NE)

nC :: Int
nC = arity (undefined :: NC)

-- TODO
-- move this to hmatrix form

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
    3200


myModelParams :: (Floating a, Ord a) => Map Text (Param NC NE a)
myModelParams = M.fromList
  [ ("ttbarnorm", \x -> (over (mBkgs.ix "ttbar") (fmap (*x)), logLogNormalP 0 0.2 x))
  , ("sigma0", set (mSigs.element 0) &&& nonNegPrior)
  , ("sigma1", set (mSigs.element 1) &&& nonNegPrior)
  , ("sigma2", set (mSigs.element 2) &&& nonNegPrior)
  , ("sigma3", set (mSigs.element 3) &&& nonNegPrior)
  , ("lumi", \x -> (over mLumi (*x), logLogNormalP 0 0.2 x))
  ]

  where
    nonNegPrior x
      | x < 0 = negate $ 1/0
      | otherwise = 0


myInitialParams :: Fractional a => Map Text a
myInitialParams = M.fromList
  [ ("ttbarnorm", 1)
  , ("sigma0", 2)
  , ("sigma1", 0.5)
  , ("sigma2", 0.1)
  , ("sigma3", 0.01)
  , ("lumi", 1)
  ]

myParamRadii :: Fractional a => Map Text a
myParamRadii = M.fromList
  [ ("ttbarnorm", 0.1)
  , ("sigma0", 0.1)
  , ("sigma1", 0.01)
  , ("sigma2", 0.01)
  , ("sigma3", 0.001)
  , ("lumi", 0.1)
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


modelLogPosterior :: (Arity m, Arity n, Integral a, Floating b, Ord b)
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
test nskip nsamps = do

  let llh = modelLogPosterior myModelParams myModel zeeData
      start = last . take 100 $ conjugateGradientAscent llh myInitialParams
      prop = weightedProposal (multibandMetropolis myParamRadii) llh
              -- weightedProposal (dynamicMetropolis $ exp <$> uniformR (-7, -3)) llh
              -- hamiltonian 3 0.005 llh (grad llh)
              -- weightedProposal $ metropolis 0.01


  g <- createSystemRandom

  let takeEvery n l = ListT $ do
        c <- next $ LT.drop n l
        case c of
          Cons x l' -> return . Cons x $ takeEvery n l'
          Nil       -> return Nil

  let chain = takeEvery nskip $ runMC prop (T start $ llh start) g

  withFile "test.dat" WriteMode
    $ \f -> do
      hPutStrLn f
        $ mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys start

      LT.runListT . LT.take nsamps
        $ do
          (T xs llhxs :: T (Map Text Double) Double) <- chain
          LT.liftIO . hPutStr f $ show llhxs ++ ", "
          LT.liftIO . hPutStrLn f
            $ mconcat . intersperse ", " . M.elems $ show <$> xs

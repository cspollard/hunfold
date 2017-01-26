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

import           Hamiltonian         (gzipWith)
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
zeeData = [6766, 2309, 288, 54]

nE :: Int
nE = arity (undefined :: NE)

nC :: Int
nC = arity (undefined :: NC)


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
  , ("sigma0", \x -> (set (mSigs.element 0) x, nonNegPrior x))
  , ("sigma1", \x -> (set (mSigs.element 1) x, nonNegPrior x))
  , ("sigma2", \x -> (set (mSigs.element 2) x, nonNegPrior x))
  , ("sigma3", \x -> (set (mSigs.element 3) x, nonNegPrior x))
  , ("sigma4", \x -> (set (mSigs.element 4) x, nonNegPrior x))
  , ("lumi", \x -> (over mLumi (*x), logLogNormalP 0 0.2 x))
  ]

  where
    nonNegPrior x = if x < 0 then (-1)/0 else 0


myInitialParams :: Fractional a => Map Text a
myInitialParams = M.fromList
  [ ("ttbarnorm", 1)
  , ("sigma0", 2)
  , ("sigma1", 2)
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
test nskip nsamps = do

  let llh = modelLogPosterior myModelParams myModel zeeData

  let start = conjugateGradientAscent llh myInitialParams !! 500

  -- TODO
  -- look at adaptive algorithms.
  -- covariance is not covering space very well.
  
  print "orig params"
  print myInitialParams

  print "orig llh:"
  print $ llh myInitialParams

  print "best params"
  print start

  print "best llh:"
  print . llh $ start

  let h = hessian (negate . llh) start
      cov = newtonInverse 50 h (goodSeed h)
      radii = imap (\i x -> x ^?! ix i * 10) cov

  print "hessian:"
  print h

  print "covariance:"
  print cov

  print "hess*cov:"
  print $ h `multM` cov

  print "radii:"
  print radii

  print "grad:"
  print $ grad llh start

  -- TODO
  -- remove test
  -- let llh m = let x = m M.! "sigma0" in logLogNormalP 0 1 x
  let prop = weightedProposal (multibandMetropolis radii) llh
              -- hamiltonian 3 0.0001 llh (grad llh)
              -- weightedProposal $ metropolis 0.01
              -- $ dynamicMetropolis $ exp <$> uniformR (-7, -3)

  g <- createSystemRandom

  let takeEvery n l = ListT $ do
        c <- next $ LT.drop n l
        case c of
          Cons x l' -> return . Cons x $ LT.drop n l'
          Nil       -> return Nil


  let l = takeEvery nskip $ runMC prop (T start $ llh start) g

  withFile "test.dat" WriteMode
    $ \f -> do
      hPutStrLn f
        $ mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys start

      LT.runListT . LT.take nsamps
        $ do
          (T xs llhxs :: T (Map Text Double) Double) <- l
          LT.liftIO . hPutStr f $ show llhxs ++ ", "
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


type M = Map Text (Map Text Double)

transposeM :: M -> M
transposeM m = imap (\i x -> imap (\j _ -> (m ^?! ix j) ^?! ix i) x) m

dotV :: (FunctorWithIndex (Index s) t, Ixed s, Foldable t, Num (IxValue s)) => t (IxValue s) -> s -> IxValue s
dotV x y = sum $ gzipWith (*) x y

multM :: M -> M -> M
multM mm nn = M.map (\x -> dotV x <$> transposeM nn) mm

addM :: M -> M -> M
addM = gzipWith (gzipWith (+))

goodSeed :: M -> M
goodSeed m =
  let mt = transposeM m
      a1 = maximum . fmap (sum . fmap abs) $ m
      a2 = maximum . fmap (sum . fmap abs) $ mt

  in (fmap.fmap) (/ (a1*a2)) mt

newtonInverse :: Int -> M -> M -> M
newtonInverse n _ seed | n < 1 = seed
newtonInverse n m seed =
  let sub = gzipWith (gzipWith (-))
      m' = seed `multM` m `multM` seed
      m'' = (fmap.fmap) (*2) seed `sub` m'

  in newtonInverse (n-1) m m''

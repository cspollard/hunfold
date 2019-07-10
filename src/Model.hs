{-# LANGUAGE DataKinds                  #-}
{-# LANGUAGE DeriveFunctor              #-}
{-# LANGUAGE DeriveGeneric              #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE NoMonomorphismRestriction  #-}
{-# LANGUAGE OverloadedStrings          #-}
{-# LANGUAGE RankNTypes                 #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE ScopedTypeVariables        #-}
{-# LANGUAGE TemplateHaskell            #-}
{-# LANGUAGE TypeFamilies               #-}


module Model
  ( Model(..)
  , mBkgs, mSig, mMig, mLumi
  , ModelVar(..)
  , mvBkgs, mvSig, mvMig, mvLumi
  , ModelParam(..)
  , mpInitialValue, mpPrior, mpVariation
  , ParamPrior(..)
  , ppToFunc
  , modelLogLikelihood
  , modelLogPosterior
  , varDiff
  , appVars
  , prediction
  ) where

import           Control.Applicative ((<|>))
import           Control.Lens
import           Control.Monad       (foldM, join)
import           Data.Aeson
import           Data.Aeson.Types
import           Data.HashMap.Strict
import           Data.Text           (Text)
import           Data.Vector         (Vector)
import           GHC.Generics
import           GHC.TypeLits

import           Matrix
import           Probability

-- TODO
-- TODO!
-- update this to work with Vars from atlas.git!!
data Model a =
  Model
    { _mBkgs :: HashMap Text (Vector a)
    , _mSig  :: Vector a
    , _mMig  :: Vector (Vector a)
    , _mLumi :: a
    } deriving (Show, Generic, Functor)

makeLenses ''Model

instance FromJSON a => FromJSON (Model a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 2}


addM :: Num a => Model a -> Model a -> Model a
addM (Model b s m l) (Model b' s' m' l') =
  Model
    (liftU2 (^+^) b b')
    (s ^+^ s')
    (liftU2 (^+^) m m')
    (l + l')


data ModelVar a =
  ModelVar
    { _mvBkgs :: Maybe (HashMap Text (Vector a))
    , _mvSig  :: Maybe (Vector a)
    , _mvMig  :: Maybe (Vector (Vector a))
    , _mvLumi :: Maybe a
    } deriving (Show, Generic, Functor)

makeLenses ''ModelVar


instance FromJSON a => FromJSON (ModelVar a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 3}


data ParamPrior a =
  Flat
  | NonNegative
  | Normal a a
  | LogNormal a a
  deriving (Generic, Functor, Show)


instance FromJSON a => FromJSON (ParamPrior a) where
  parseJSON (String "Flat") = return Flat
  parseJSON (String "NonNegative") = return NonNegative

  parseJSON (Object o) =
      (o .: "Normal" >>= p Normal)
      <|> (o .: "LogNormal" >>= p LogNormal)
    where
      p f (Object m) = f <$> m .: "Mu" <*> m .: "Sigma"
      p _ invalid    = typeMismatch "ParamPrior" invalid

  parseJSON invalid = typeMismatch "ParamPrior" invalid


ppToFunc :: (Ord a, Floating a) => ParamPrior a -> (a -> a)
ppToFunc Flat            = const 0
ppToFunc NonNegative = \x ->
  if x < 0 then negate (1/0) else 0
ppToFunc (Normal m s)    = logNormalP m s
ppToFunc (LogNormal m s) = logLogNormalP m s


data ModelParam a =
  ModelParam
    { _mpInitialValue :: a
    , _mpPrior        :: ParamPrior a
    , _mpVariation    :: ModelVar a
    } deriving (Generic, Functor)

makeLenses ''ModelParam

instance FromJSON a => FromJSON (ModelParam a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 3}


modelLogPosterior
  :: (Ord a, Floating a, Integral b)
  => Vector b
  -> Model a
  -> Vector (ModelVar a)
  -> Vector (a -> a)
  -> (Vector a -> a)
  -> Vector a
  -> Either String a
modelLogPosterior dats model mps logPriors logReg ps =
  reifyVector ps $ \ps' -> do
    mps' <- toEither "incorrect length of mps" $ fromVector mps
    logPriors' <- toEither "incorrect length of logPriors" $ fromVector logPriors
    model' <- appVars mps' ps' model
    logLike <- toEither "incorrect length of data" $ modelLogLikelihood dats model'

    let logPrior = sum $ logPriors' <*> ps'

    return $ logLike + logPrior + logReg (_mSig model')



prediction'
  :: (KnownNat n, KnownNat m, Num a)
  => HashMap Text (V m a) -> V n a -> M n m a -> a -> V m a
prediction' bkgs sig mig lumi =
  let bkgtot = foldl (^+^) zero bkgs
      sigtot = sig *! mig
  in lumi *^ (sigtot ^+^ bkgtot)


prediction :: Num a => Model a -> Either String (Vector a)
prediction Model{..} =
  toEither "backgrounds have incorrect length." . join
    $ reifyVectorNat _mSig
      $ \sig -> reifyMatrix1 _mMig
        $ \mig -> do
          bkgs <- traverse fromVector _mBkgs
          return . toVector $ prediction' bkgs sig mig _mLumi


modelLogLikelihood
  :: (Integral b, Floating a)
  => Vector b -> Model a -> Maybe a
modelLogLikelihood dats Model{..} =
  join $ reifyVectorNat _mSig
    $ \sig -> reifyMatrix1 _mMig
      $ \mig -> do
        bkgs <- traverse fromVector _mBkgs
        dat <- fromVector dats
        return . sum $ logPoissonP <$> dat <*> prediction' bkgs sig mig _mLumi


appVars
  :: (Additive t, Foldable t, Num a)
  => t (ModelVar a) -> t a -> Model a -> Either String (Model a)
appVars mvs ps m = foldM (fmap . addM) m ms
  where ms = ($ m) <$> liftI2 varDiff mvs ps


bkgDiff
  :: (KnownNat n, Num a)
  => a -> HashMap Text (V n a) -> HashMap Text (V n a) -> HashMap Text (V n a)
bkgDiff x = liftU2 (sigDiff x)


sigDiff :: (KnownNat n, Num a) => a -> V n a -> V n a -> V n a
sigDiff x s = (x *^) . (^-^ s)


migDiff
  :: (KnownNat n, KnownNat m, Num a)
  => a -> M n m a -> M n m a -> M n m a
migDiff x m = (x *!!) . (!-! m)


lumiDiff :: Num c => c -> c -> c -> c
lumiDiff x l = ((x-1) *) . flip (-) l


-- the outgoing "Model" is really a difference between the nominal model and
-- the varied model.
varDiff :: Num a => ModelVar a -> a -> Model a -> Either String (Model a)
varDiff ModelVar{..} x Model{..} = toEither "failed to apply model variation." $
  join $ reifyVectorNat _mSig
    $ \(sig :: V n a) -> reifyMatrix1 _mMig
      $ \(mig :: M n m a) -> do

        (bkgs :: HashMap Text (V m a)) <- traverse fromVector _mBkgs

        -- this is a bit tricky
        -- it is Nothing if the vector conversion failed.
        -- it is Just Nothing if there was Nothing in there to begin with.
        vBkgs <- sequenceA $ traverse fromVector <$> _mvBkgs
        vSig <- sequenceA $ fromVector <$> _mvSig
        vMig <- sequenceA $ fromVectorM <$> _mvMig

        return
          $ Model
            (toVector <$> maybe zero (bkgDiff x bkgs) vBkgs)
            (toVector $ maybe zero (sigDiff x sig) vSig)
            (toVectorM $ maybe zero (migDiff x mig) vMig)
            (maybe 0 (lumiDiff x _mLumi) _mvLumi)


toEither :: a -> Maybe b -> Either a b
toEither x = maybe (Left x) Right

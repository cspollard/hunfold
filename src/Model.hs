{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}


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
  , appVar
  , appVars
  , prediction
  ) where

import Control.Lens
import GHC.TypeLits
import Control.Monad (join, foldM)
import GHC.Generics
import Data.Map (Map)
import Data.Vector (Vector)
import Data.Text (Text)
import Data.Aeson
import Data.Aeson.Types
import Control.Applicative ((<|>))

import Probability
import Matrix

data Model a =
  Model
    { _mBkgs :: Map Text (Vector a)
    , _mSig:: Vector a
    , _mMig:: Vector (Vector a)
    , _mLumi  :: a
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
    { _mvBkgs :: Maybe (Map Text (Vector a))
    , _mvSig:: Maybe (Vector a)
    , _mvMig:: Maybe (Vector (Vector a))
    , _mvLumi  :: Maybe a
    } deriving (Show, Generic, Functor)

makeLenses ''ModelVar


instance FromJSON a => FromJSON (ModelVar a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 3}


data ParamPrior a =
  Flat
  | Normal a a
  | LogNormal a a
  deriving (Generic, Functor, Show)

instance FromJSON a => FromJSON (ParamPrior a) where
  parseJSON (String "Flat") = return Flat

  parseJSON (Object o) =
    (o .: "Normal" >>= p Normal)
      <|> (o .: "LogNormal" >>= p LogNormal)
    where
      p f (Object m) = f <$> m .: "Mu" <*> m .: "Sigma"
      p _ invalid    = typeMismatch "ParamPrior" invalid

  parseJSON invalid = typeMismatch "ParamPrior" invalid


ppToFunc :: Floating a => ParamPrior a -> (a -> a)
ppToFunc Flat = const 0
ppToFunc (Normal m s) = logNormalP m s
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
  :: (Floating a, Ord a, Integral b)
  => Vector b
  -> Model a
  -> Vector (ModelVar a)
  -> Vector (a -> a)
  -> Vector a
  -> Either String a
modelLogPosterior dats model mps logPriors ps =
  reifyVector ps $ \ps' -> do
    mps' <- toEither "incorrect length of mps" $ fromVector mps
    logPriors' <- toEither "incorrect length of logPriors" $ fromVector logPriors
    model' <- appVars mps' ps' model
    logLike <- toEither "incorrect length of dats" $ modelLogLikelihood dats model'

    let logPrior = sum $ logPriors' <*> ps'

    return $ logLike + logPrior



prediction'
  :: (KnownNat n, KnownNat m, Num a)
  => Map Text (V m a) -> V n a -> M n m a -> a -> V m a
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
  :: (Integral b, Floating a, Ord a)
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
  where ms = ($ m) <$> liftI2 appVar mvs ps


mergeBkgs :: (KnownNat n, Num a) => a -> Map Text (V n a) -> Map Text (V n a) -> Map Text (V n a)
mergeBkgs x b = fmap (x *^) . liftU2 (^-^) b


mergeSig :: (KnownNat n, Num a) => a -> V n a -> V n a -> V n a
mergeSig x s = (x *^) . (^-^ s)


mergeMig :: (KnownNat n, KnownNat m, Num a) => a -> M n m a -> M n m a -> M n m a
mergeMig x m = (x *!!) . (!-! m)


mergeLumi :: Num c => c -> c -> c -> c
mergeLumi x l = (x *) . flip (-) l


appVar :: Num a => ModelVar a -> a -> Model a -> Either String (Model a)
appVar ModelVar{..} x Model{..} = toEither "appVar failed" $
  join $ reifyVectorNat _mSig
    $ \(sig :: V n a) -> reifyMatrix1 _mMig
      $ \(mig :: M n m a) -> do

        (bkgs :: Map Text (V m a)) <- traverse fromVector _mBkgs

        -- this is a bit tricky
        -- it is Nothing if the vector conversion failed.
        -- it is Just Nothing if there was Nothing in there to begin with.
        vBkgs <- sequenceA $ traverse fromVector <$> _mvBkgs
        vSig <- sequenceA $ fromVector <$> _mvSig
        vMig <- sequenceA $ fromVectorM <$> _mvMig

        return
          $ Model
            (toVector <$> maybe zero (mergeBkgs x bkgs) vBkgs)
            (toVector $ maybe zero (mergeSig x sig) vSig)
            (toVectorM $ maybe zero (mergeMig x mig) vMig)
            (maybe 0 (mergeLumi x _mLumi) _mvLumi)


toEither :: a -> Maybe b -> Either a b
toEither x = maybe (Left x) Right

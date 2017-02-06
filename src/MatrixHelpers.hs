{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise #-}
{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}

module MatrixHelpers
  ( module X
  , M
  , toVectorM
  , fromVectorM
  , reifyMatrix
  , reifyMatrix1
  , reifyMatrix2
  , reifySqMatrix
  ) where

import           Data.Proxy
import           Data.Vector  (Vector)
import qualified Data.Vector  as V
import           GHC.TypeLits
import           Linear       as X
import           Linear.V     as X

type M (n :: Nat) (m :: Nat) a = V n (V m a)

instance KnownNat n => Trace (V n)

fromVectorM :: (Dim n, Dim m) => Vector (Vector a) -> Maybe (M n m a)
fromVectorM m = do
  vs <- traverse fromVector m
  fromVector vs

toVectorM :: M n m a -> Vector (Vector a)
toVectorM = toVector . fmap toVector


reifyMatrix :: Vector (Vector a) -> (forall n m. (KnownNat n, KnownNat m) => M n m a -> r) -> Maybe r
reifyMatrix m f = do
  let l1 = V.length m
      l2 = if l1 == 0 then 0 else V.length (m V.! 0)
  reifyDimNat l1 $
    \(Proxy :: Proxy k1) -> reifyDimNat l2 $
      \(Proxy :: Proxy k2) -> do
        (mV :: M k2 k1 a) <- fromVectorM m
        return $ f mV

reifySqMatrix
  :: forall a r. Vector (Vector a) -> (forall n. KnownNat n => M n n a -> r) -> Maybe r
reifySqMatrix m f = do
  let l1 = V.length m
  reifyDimNat l1 $
    \(Proxy :: Proxy n) -> do
      (mV :: M n n a) <- fromVectorM m
      return $ f mV


reifyMatrix1
  :: forall n a r. KnownNat n
  => Vector (Vector a) -> (forall m. KnownNat m => M n m a -> r) -> Maybe r
reifyMatrix1 m f = do
  let l1 = reflectDim (Proxy :: Proxy n)
      l2 = if l1 == 0 then 0 else V.length (m V.! 0)
  reifyDimNat l2 $
    \(Proxy :: Proxy k2) -> do
      (mV :: M n k2 a) <- fromVectorM m
      return $ f mV


reifyMatrix2
  :: forall m a r. KnownNat m
  => Vector (Vector a) -> (forall n. KnownNat n => M n m a -> r) -> Maybe r
reifyMatrix2 m f = do
  let l1 = V.length m
  reifyDimNat l1 $
    \(Proxy :: Proxy k1) -> do
      (mV :: M k1 m a) <- fromVectorM m
      return $ f mV

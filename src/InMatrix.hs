{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TupleSections             #-}
{-# LANGUAGE TypeFamilies              #-}

module InMatrix
  ( module X
  , toMatrix, inM, inM2
  , invM, cholM
  ) where

import           Control.Lens hiding (Indexable)
import           Data.Matrix  as X

-- NB
-- Matrix is indexed starting at 1
toMatrix :: (Traversable t, Traversable t1)
         => t (t1 b) -> (Matrix b -> t (t1 b), Matrix b)
toMatrix xs =
  let back = iover traversed (const . (+1))
      bb = iover traversed (\i x -> (i+1,) <$> back x) xs
      f m = over (traverse.traverse) (m !) bb
      ls = toListOf (traverse . to (toListOf traverse)) xs
  in (f, fromLists ls)

inM
  :: (Functor f, Traversable t, Traversable t1)
  => (Matrix b -> f (Matrix b)) -> t (t1 b) -> f (t (t1 b))
inM f xs =
  let (b, m) = toMatrix xs
  in b <$> f m

inM2
  :: (Functor f, Traversable t, Traversable t1)
  => (Matrix a -> Matrix a -> f (Matrix a)) -> t (t1 a) -> t (t1 a) -> f (t (t1 a))
inM2 f xs ys =
  let (bx, mx) = toMatrix xs
      (_, my) = toMatrix ys
  in bx <$> f mx my

invM :: (Traversable t1, Traversable t, Fractional b, Eq b) => t (t1 b) -> Either String (t (t1 b))
invM = inM inverse

cholM :: (Traversable t1, Traversable t, Floating b) => t (t1 b) -> t (t1 b)
cholM = runIdentity . inM (Identity . cholDecomp)

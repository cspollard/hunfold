{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE TupleSections             #-}
{-# LANGUAGE TypeFamilies             #-}
{-# LANGUAGE RankNTypes             #-}

module InMatrix
  ( module X
  , toMatrix, inM, inM2
  , toVector, inV
  , invM, eigM, cholM
  ) where

import           Control.Lens          hiding (Indexable)
import           Numeric.LinearAlgebra as X

-- TODO
-- these could probably be rewritten in terms of the Vector equivalents below.
toMatrix :: (Container Vector b, Num b, Traversable t, Traversable t1)
         => t (t1 b) -> (Matrix b -> t (t1 b), Matrix b)
toMatrix xs =
  let back = iover traversed const
      bb = iover traversed (\i x -> (i,) <$> back x) xs
      f m = over (traverse.traverse) (atIndex m) bb
      ls = toListOf (traverse . to (toListOf traverse)) xs
  in (f, fromLists ls)

inM :: (Container Vector b, Num b, Traversable t, Traversable t1) => (Matrix b -> Matrix b) -> t (t1 b) -> t (t1 b)
inM f xs =
  let (b, m) = toMatrix xs
  in b $ f m

inM2 :: (Container Vector a, Traversable t, Traversable t1, Num a) => (Matrix a -> Matrix a -> Matrix a) -> t (t1 a) -> t (t1 a) -> t (t1 a)
inM2 f xs ys =
  let (bx, mx) = toMatrix xs
      (_, my) = toMatrix ys
  in bx $ f mx my

invM :: (Field b, Traversable t1, Traversable t) => t (t1 b) -> t (t1 b)
invM = inM (pinvTol 1e-3)

cholM :: (Field b, Traversable t1, Traversable t) => t (t1 b) -> t (t1 b)
cholM = inM (chol . sym)

eigM :: (Traversable t1, Traversable t) => t (t1 (Complex Double)) -> t (t1 (Complex Double))
eigM = inM (snd . eig)


toVector :: (Container Vector b, Traversable t)
         => t b -> (Vector b -> t b, Vector b)
toVector xs =
  let back = iover traversed const xs
      f v = over traverse (atIndex v) back
      vs = toListOf traverse xs
  in (f, fromList vs)

inV :: (Container Vector b, Traversable t) => (Vector b -> Vector b) -> t b -> t b
inV f xs =
  let (b, v) = toVector xs
  in b $ f v

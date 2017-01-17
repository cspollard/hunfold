{-# LANGUAGE DataKinds #-}

module Dagostini where

import           Control.Monad.Primitive
import           Data.Foldable
import           Data.Monoid
import           Data.Traversable
import           Data.Vector.Fixed.Cont        (Arity (..), ContVec (..),
                                                ToPeano (..))
import qualified Data.Vector.Fixed.Cont        as V
import           System.Random.MWC.Probability hiding (dirichlet)

-- | The dirichlet distribution.
dirichlet :: (Traversable f, PrimMonad m) => f Double -> Prob m (f Double)
dirichlet as = do
  zs <- mapM (`gamma` 1) as
  return $ let (s, zs') = mapAccumL (\x z -> (x+z, z/s)) 0 zs in zs'
{-# INLINABLE dirichlet #-}


type Vec = ContVec

instance (Show a, Arity n) => Show (ContVec n a) where
  show = show . toList

-- rows are innermost
type Mat n m a = Vec m (Vec n a)

transpose :: (Arity n, Arity m) => Mat n m a -> Mat m n a
transpose = sequenceA

inner :: (Monoid c, Foldable t, Applicative t)
      => (a -> b -> c) -> t a -> t b -> c
inner f v v' = fold $ f <$> v <*> v'

dot :: (Num a, Foldable t, Applicative t)
    => t a -> t a -> a
dot v v' = getSum $ inner (\x y -> Sum $ x * y) v v'

multMM :: (Traversable t1, Traversable t, Applicative t1, Applicative f, Num b)
       => t (t1 b) -> t1 (f b) -> f (t b)
multMM m m' = outer dot m $ sequenceA m'

multMV :: (Traversable t, Foldable f, Applicative f, Num b)
       => t (f b) -> f b -> t b
multMV = traverse dot

multVM :: (Traversable t1, Traversable t, Applicative t1, Applicative t, Num b)
       => t1 b -> t1 (t b) -> t b
multVM v m = traverse dot (sequenceA m) v

-- arbitrary cartesian product
infix 5 ><
(><) :: (Traversable t, Functor f) => t a -> f b -> f (t (a, b))
(><) = outer (,)

outer :: (Traversable t, Functor f)
      => (a1 -> a -> b) -> t a1 -> f a -> f (t b)
outer f v v' = traverse f v <$> v'


-- test case:
-- 3 data bins
-- 2 cause bins
-- 1 trash bin

myData :: Vec (ToPeano 3) Int
myData = V.fromList' [1, 2, 3]

myLambdas :: Mat (ToPeano 2) (ToPeano 4) Double
myLambdas = V.fromList'
  [ V.fromList' [ 0.25, 0.25 ]
  , V.fromList' [ 0.25, 0.50 ]
  , V.fromList' [ 0.00, 0.25 ]
  , V.fromList' [ 0.50, 0.00 ]
  ]

myPrior :: Vec (ToPeano 2) Double
myPrior = V.fromList' [1, 1000]

myThetas :: Mat (ToPeano 4) (ToPeano 2) Double
-- TODO
-- this line could be cleaned I'm sure
myThetas = (\x -> (/) <$> x <*> denom) <$> numer
  where
    ljis = transpose myLambdas
    pis = myPrior
    -- TODO
    -- this line could be cleaned I'm sure
    numer = (\x y -> fmap (*y) x) <$> ljis <*> pis
    denom = multVM pis ljis

myEfficiency :: Vec (ToPeano 2) Double
myEfficiency = (-) <$> pure 1 <*> V.head myLambdas

myMeas :: Vec (ToPeano 2) Double
myMeas = (/) <$> multMV smear dat <*> myEfficiency
  where
    smear = (transpose . V.tail . transpose) myThetas
    dat = fromIntegral <$> myData



{-

-- Ci: cause bin i
-- xi: events in cause bin i (what we are measuring)
-- xij: events in cause bin i and effect bin j
-- Ej: effect bin j
-- yj: events in effect bin j (what we observe)

-- P(Ej|Ci):
--   probability that an event lands in effect bin j given an event in cause
--   bin i, ie. the smearing matrix
-- P(Ci|Ej):
--   probability that an event lands in cause bin i given an event in effect
--   bin j

-- Bayes' theorem
-- P(x) = P(x|y) P(y) / P(y|x)
-- ie. P(Ci) = P(Ci|Ej) P(Ej) / P(Ej|Ci)

-- so P(Ci|Ej) = P(Ej|Ci) P(Ci) / P(Ej)
--   where
--     P(Ej) = sum_i P(Ej|Ci) P(Ci)

-- so xij = P(Ci|Ej) yj
-- xi = (sum_j P(Ci|Ej) yj) / (sum_j P(Ej|Ci))
--   where sum_j P(Ej|Ci) is the observation efficiency of cause bin i

-- we observe a spectrum yj (integers) in n effect bins.
-- we hope to measure differential cross sections sigma_i in m cause bins.
-- we have m*(n+1) P(Ej|Ci) terms corresponding to the smearing matrix.
--   (n+1 takes into account the inefficiency)
-- we need to marginalize over uncertainties in our backgrounds.
-- we need to marginalize over uncertainties in our signal model, P(Ej|Ci).
-- we need to deal with the prior distribution P(Ci) (iterative?)
-- given a particular P(Ej|Ci), a background prediction bj and the obvervations yj,
--   we can obtain P(mui|yj) = P(yj|mui) P(mui) / P(yj)

-}

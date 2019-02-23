{-# LANGUAGE DeriveFunctor              #-}
{-# LANGUAGE GADTs                      #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE NoMonomorphismRestriction  #-}
{-# LANGUAGE TemplateHaskell            #-}
{-# LANGUAGE TypeFamilies               #-}

module NDSL where

import           Control.Monad.Free
import           Numeric.AD
import           Text.Show.Deriving

data NDSLFun a where
  Var :: String -> NDSLFun a
  Times :: a -> a -> NDSLFun a
  Plus :: a -> a -> NDSLFun a
  Abs :: a -> NDSLFun a
  Signum :: a -> NDSLFun a
  Neg :: a -> NDSLFun a
  Recip :: a -> NDSLFun a
  Exp :: a -> NDSLFun a
  Log :: a -> NDSLFun a
  Sqrt :: a -> NDSLFun a
  Sin :: a -> NDSLFun a
  Cos :: a -> NDSLFun a
  ASin :: a -> NDSLFun a
  ACos :: a -> NDSLFun a
  ATan :: a -> NDSLFun a
  Sinh :: a -> NDSLFun a
  Cosh :: a -> NDSLFun a
  ASinh :: a -> NDSLFun a
  ACosh :: a -> NDSLFun a
  ATanh :: a -> NDSLFun a
  deriving (Eq, Ord, Show, Read, Functor)


deriveShow1 ''NDSLFun


var :: String -> NDSL a
var = wrap . Var

newtype NDSL a = NDSL { runNDSL :: Free NDSLFun a }
  deriving (Functor, Applicative, Monad, MonadFree NDSLFun, Show)

instance Mode a => Mode (NDSL a) where
  type Scalar (NDSL a) = a
  auto = pure

instance Num a => Num (NDSL a) where
  fromInteger = pure . fromInteger
  x * y =  wrap $ Times x y
  x + y =  wrap $ Plus x y
  negate = wrap . Neg
  abs = wrap . Abs
  signum = wrap . Signum

instance Fractional a => Fractional (NDSL a) where
  fromRational = pure . fromRational
  recip = wrap . Recip

instance Floating a => Floating (NDSL a) where
  pi = pure pi
  exp = wrap . Exp
  log = wrap . Log
  sqrt = wrap . Sqrt
  sin = wrap . Sin
  cos = wrap . Cos
  asin = wrap . ASin
  acos = wrap . ACos
  atan = wrap . ATan
  sinh = wrap . Sinh
  cosh = wrap . Cosh
  asinh = wrap . ASinh
  acosh = wrap . ACosh
  atanh = wrap . ATanh

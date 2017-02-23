{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

{-
 - note that parts of this module are taken
 - directly from the source code of the mwc-probability
 - package with thanks to the authors.
 -}

module Probability
  ( module X
  , logNormalP, logLogNormalP, logPoissonP
  , standard, normal, exponential, truncatedExp
  , gamma, chiSquare, bernoulli, dirichlet
  , T(..)
  ) where

import           Control.Monad.Primitive       as X (PrimMonad, PrimState)
import           Data.Monoid                   (Sum (..), (<>))
import           Data.Number.Erf               as X
import           Data.Traversable              (mapAccumL)
import           System.Random.MWC             as X (Gen, asGenIO,
                                                     createSystemRandom,
                                                     withSystemRandom)
import           System.Random.MWC.Probability as X (Prob, Variate, sample,
                                                     samples, uniform, uniformR)

lnsqrt2pi
    :: Fractional a
    => a
lnsqrt2pi =
    0.9189385332046727417803297364056176398613974736377834128171
{-# INLINABLE lnsqrt2pi #-}
{-# SPECIALIZE lnsqrt2pi :: Double #-}


logNormalP
    :: Floating a
    => a -> a -> a -> a
logNormalP m s x =
    let e = negate . (/ 2) . raising 2 $ (x - m) / s
    in e - log s - lnsqrt2pi
{-# INLINABLE logNormalP #-}
{-# SPECIALIZE logNormalP :: Double -> Double -> Double -> Double #-}


logLogNormalP
    :: Floating a
    => a -> a -> a -> a
logLogNormalP m s x =
    let x' = log x
        e = negate . (/ 2) . raising 2 $ (x' - m) / s
    in e - log s - x' - lnsqrt2pi
{-# INLINABLE logLogNormalP #-}
{-# SPECIALIZE logLogNormalP :: Double -> Double -> Double -> Double #-}

logPoissonP
    :: (Integral a, Floating b, Ord b)
    => a -> b -> b
logPoissonP k l
  | k < 0 = error "logPoissonP: negative counts"
  | otherwise = fromIntegral k * log l - l - logFactorial k
{-# INLINABLE logPoissonP #-}
{-# SPECIALIZE logPoissonP :: Integral a => a -> Double -> Double #-}
{-# SPECIALIZE logPoissonP :: Int -> Double -> Double #-}

logFactorial
    :: (Integral a, Floating b)
    => a -> b
logFactorial k
  | k < 0 = error "logFactorial: negative integer input!"
  | k <= 14 = log . fromIntegral $ factorial k
  | otherwise =
      (x - 0.5) * log x - x + lnsqrt2pi + 1 / (12 * x) -
      1 / (360 * raising 3 x) +
      1 / (1260 * raising 5 x) -
      1 / (1680 * raising 7 x)
  where
    x = fromIntegral (k + 1)
    factorial 0 = 1
    factorial 1 = 1
    factorial n = n * factorial (n - 1)
{-# INLINABLE logFactorial #-}
{-# SPECIALIZE logFactorial :: Integral a => a -> Double #-}
{-# SPECIALIZE logFactorial :: Int -> Double #-}

raising :: Num a => Int -> a -> a
raising y x = x ^ y
{-# INLINABLE raising #-}
{-# SPECIALIZE raising :: Int -> Double -> Double #-}


-- strict 2-tuple
data T a b =
  T
    { fstT :: !a
    , sndT :: !b
    } deriving Show

standard :: (Variate b, PrimMonad m, InvErf b)
         => Prob m b
standard = do
  x <- uniform
  return $! inverf (2*x-1) * (1 - sqrt 2)

normal :: (Variate b, PrimMonad m, InvErf b)
       => b -> b -> Prob m b
normal m s = do
  x <- standard
  return $! x*s + m

exponential :: (Variate b, PrimMonad m, Floating b)
            => b -> Prob m b
exponential b = do
  x <- uniform
  return $! - log x / b

truncatedExp :: (Variate b, PrimMonad m, Floating b) => b -> (b, b) -> Prob m b
truncatedExp scale (a,b) = do
  -- We shift a to 0 and then generate distribution truncated to [0,b-a]
  -- It's easier
  let delta = b - a
  p <- uniform
  return $! a - log ( (1 - p) + p*exp(-scale*delta)) / scale

pkgError :: String -> String -> a
pkgError func msg = error $ "System.Random.MWC.Distributions." ++ func ++
                            ": " ++ msg

sqr :: Num a => a -> a
sqr x = x*x
{-# INLINABLE sqr #-}
{-# SPECIALIZE sqr :: Double -> Double #-}

gamma :: (Variate a, PrimMonad m, InvErf a, Ord a) => a -> a -> Prob m a
gamma a b
  | a <= 0    = pkgError "gamma" "negative alpha parameter"
  | otherwise = mainloop
    where
      mainloop = do
        T x v <- innerloop
        u     <- uniform
        let cont =  u > 1 - 0.331 * sqr (sqr x)
                 && log u > 0.5 * sqr x + a1 * (1 - v + log v) -- Rarely evaluated
        case () of
          _| cont      -> mainloop
           | a >= 1    -> return $! a1 * v * b
           | otherwise -> do y <- uniform
                             return $! y ** (1 / a) * a1 * v * b
      -- inner loop
      innerloop = do
        x <- standard
        case 1 + a2*x of
          v | v <= 0    -> innerloop
            | otherwise -> return $! T x (v*v*v)
      -- constants
      a' = if a < 1 then a + 1 else a
      a1 = a' - 1/3
      a2 = 1 / sqrt(9 * a1)

chiSquare :: (Variate b, PrimMonad m, InvErf b, Ord b, Integral a) => a -> Prob m b
chiSquare n
  | n <= 0    = pkgError "chiSquare" "number of degrees of freedom must be positive"
  | otherwise = do
    x <- gamma (0.5 * fromIntegral n) 1
    return $! 2 * x

bernoulli :: (Variate a, PrimMonad m, Ord a) => a -> Prob m Bool
bernoulli p = (< p) <$> uniform

dirichlet :: (Traversable f, Variate a, Fractional a, InvErf a, Ord a, PrimMonad m) => f a -> Prob m (f a)
dirichlet as = do
  zs <- traverse (`gamma` 1) as
  return $ zs `normTo` 1

-- use Sum here to avoid mixups
normTo :: (Fractional a, Traversable f) => f a -> a -> f a
normTo xs y = let (Sum s, xs') = mapAccumL (\sofar x -> (Sum x <> sofar, x*y/s)) 0 xs in xs'

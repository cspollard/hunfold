{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module Likelihood
  ( module X
  , logNormalP
  , logLogNormalP
  , logPoissonP
  , runMC, findMaximum
  , weightedProposal, metropolis
  , standard, normal
  , exponential, truncatedExp
  , gamma
  , chiSquare
  , testMetro
  , T(..)
  ) where

import           Control.Monad.Primitive       as X (PrimMonad, PrimState)
import           Data.Functor.Identity
import           Data.Number.Erf               as X
import           List.Transformer              as X (ListT (..), Step (..),
                                                     liftIO, runListT)
import           Numeric.AD
import           System.Random.MWC             as X (Gen, asGenIO,
                                                     withSystemRandom)
import           System.Random.MWC.Probability as X (Prob, Variate, sample,
                                                     samples, uniform)

import           Numeric.AD.Internal.Forward   (Forward)
import           Numeric.AD.Internal.Kahn
import           Numeric.AD.Internal.On
import           Numeric.AD.Internal.Or        hiding (T)

lnsqrt2pi
    :: Fractional a
    => a
lnsqrt2pi =
    0.9189385332046727417803297364056176398613974736377834128171


logNormalP
    :: Floating a
    => a -> a -> a -> a
logNormalP m s x =
    let e = negate . (/ 2) . raising 2 $ (x - m) / s
    in e - log s - lnsqrt2pi

logLogNormalP
    :: Floating a
    => a -> a -> a -> a
logLogNormalP m s x =
    let x' = log x
        e = negate . (/ 2) . raising 2 $ (x' - m) / s
    in e - log s - x' - lnsqrt2pi

logPoissonP
    :: (Integral a, Floating b)
    => a -> b -> b
logPoissonP k l = fromIntegral k * log l - l - logFactorial k

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

raising :: Num a => Int -> a -> a
raising y x = x ^ y

-- ~from https://hackage.haskell.org/package/declarative-0.2.1/docs/src/Numeric-MCMC.html#mcmc
-- a Markov chain driven by an arbitrary transition operator.
runMC :: PrimMonad m
      => (a -> Prob m a) -> a -> Gen (PrimState m) -> ListT m a
runMC t s g = ListT $ do
  s' <- sample (t s) g
  return . Cons s' $ runMC t s' g

testMetro :: IO ()
testMetro = do
  withSystemRandom . asGenIO $
    \g -> runListT $
      do
        x <- runMC prop (T (Identity 0) 1) g
        liftIO $ print (x :: T (Identity Double) Double)

  return ()

  where prop = weightedProposal (metropolis 1) (const 1)

weightedProposal :: (RealFloat b, Variate b, PrimMonad m)
                 => (a -> Prob m a) -> (a -> b) -> T a b -> Prob m (T a b)
weightedProposal t logLH (T x y) =
  do
    prop <- t x
    let score = logLH prop
    -- TODO
    -- what happens with NaNs?
    let prob = exp . min 0 $ score - y
    accept <- bernoulli prob
    return $ if accept then T prop score else T x y



-- TODO
-- if we knew f was a "zipping" type then we could have many different radii
-- should we pick a particular type here?
-- TODO
-- does this actually sample normal n times or only once?
metropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t) => b -> t b -> Prob m (t b)
metropolis r = traverse (\x -> (+x) <$> normal 0 r)


findMaximum :: (Traversable f, Ord a, Fractional a)
            => (forall s. Chosen s => f (Or s (On (Forward (Forward a))) (Kahn a)) -> Or s (On (Forward (Forward a))) (Kahn a))
            -> f a -> f a
findMaximum llhf start = last . take n . conjugateGradientAscent llhf $ start
  where n = length start + 1

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

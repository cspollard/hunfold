{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , runMC, findMaximum
  , weightedProposal
  ) where

import           List.Transformer              as X (ListT (..), Step (..),
                                                     liftIO, runListT)
import           Numeric.AD
import           Numeric.AD.Internal.Forward   (Forward)
import           Numeric.AD.Internal.Kahn
import           Numeric.AD.Internal.On
import           Numeric.AD.Internal.Or        hiding (T)
import           Probability
import           System.Random.MWC             as X (Gen, asGenIO,
                                                     withSystemRandom)
import           System.Random.MWC.Probability as X (Prob, Variate, sample,
                                                     samples, uniform)

runMC :: PrimMonad m
      => (a -> Prob m a) -> a -> Gen (PrimState m) -> ListT m a
runMC t s g = ListT $ do
  s' <- sample (t s) g
  return . Cons s' $ runMC t s' g


-- a proposal weighted by a log likelihood function
weightedProposal :: (Floating b, Ord b, Variate b, PrimMonad m)
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


findMaximum :: (Traversable f, Ord a, Fractional a)
            => (forall s. Chosen s => f (Or s (On (Forward (Forward a))) (Kahn a)) -> Or s (On (Forward (Forward a))) (Kahn a))
            -> f a -> f a
findMaximum llhf start = last . take n . conjugateGradientAscent llhf $ start
  where n = length start + 1

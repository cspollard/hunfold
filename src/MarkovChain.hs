{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , runMC
  , weightedProposal
  ) where

import           List.Transformer              as X (ListT (..), Step (..),
                                                     liftIO, runListT)
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
weightedProposal :: (Floating b, Variate b, PrimMonad m, RealFloat b)
                 => (a -> Prob m a) -> (a -> b) -> T a b -> Prob m (T a b)
weightedProposal t logLH (T x y) =
  do
    prop <- t x
    let score = logLH prop
    -- TODO
    -- what happens with NaNs?
    let prob = exp . min 0 $ score - y
    accept <- if isNaN prob then return False else bernoulli prob
    return $ if accept then T prop score else T x y

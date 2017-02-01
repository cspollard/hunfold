{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , runMC
  ) where

import           List.Transformer                 as X (ListT (..), Step (..),
                                                        liftIO, runListT)
import           Numeric.MCMC

import           Control.Monad.Trans.State.Strict
import           System.Random.MWC                as X (Gen, asGenIO,
                                                        withSystemRandom)
import           System.Random.MWC.Probability    as X (Prob, Variate, sample,
                                                        samples, uniform)

runMC :: PrimMonad m
      => Transition m a -> a -> Gen (PrimState m) -> ListT m a
runMC t s g = ListT $ do
  s' <- sample (execStateT t s) g
  return . Cons s' $ runMC t s' g

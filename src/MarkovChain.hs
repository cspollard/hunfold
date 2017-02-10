{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , runMC
  ) where

import           Control.Monad.Primitive          as X (PrimMonad (..))
import           Control.Monad.Trans.State.Strict
import           Data.Sampling.Types              as X (Transition)
import           List.Transformer                 as X (ListT (..), Step (..),
                                                        liftIO, runListT)
import           System.Random.MWC                as X (Gen, asGenIO,
                                                        withSystemRandom)
import           System.Random.MWC.Probability    as X (sample)

runMC :: PrimMonad m
      => Transition m a -> a -> Gen (PrimState m) -> ListT m a
runMC t s g = ListT $ do
  s' <- sample (execStateT t s) g
  return . Cons s $ runMC t s' g

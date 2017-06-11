{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , markovChain, metropolis
  ) where

import           Control.Monad.Primitive       as X (PrimMonad (..))
import           Pipes
import           Pipes.Prelude                 as P
import           System.Random.MWC             as X (Gen, asGenIO,
                                                     withSystemRandom)
import           System.Random.MWC.Probability

markovChain :: Monad m => (t -> m t) -> t -> Producer t m r
markovChain step start = do
  -- if we don't yield the starting location it's not part of the chain!
  yield start
  P.unfoldr (\x -> Right <$> (fmap dup . step) x) start

  where dup x = (x, x)


metropolis
  :: PrimMonad m
  => (t -> Prob m t)
  -> (t -> Double)
  -> (t, Double)
  -> Producer (t, Double) (Prob m) r
metropolis prop logLH = markovChain $ metropolisStep prop logLH


metropolisStep
  :: PrimMonad m
  => (t -> Prob m t)
  -> (t -> Double)
  -> (t, Double)
  -> Prob m (t, Double)
metropolisStep proposal logLH (currloc, currllh) = do
  nextloc <- proposal currloc
  let nextllh = whenNaN negInf $ logLH nextloc
  accept <- bernoulli . exp . min 0 $ nextllh - currllh
  if accept
    then return (nextloc, nextllh)
    else return (currloc, currllh)

  where
    negInf = negate $ 1 / 0
    whenNaN x y = if isNaN y then x else y

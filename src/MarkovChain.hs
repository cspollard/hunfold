{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}

module MarkovChain
  ( module X
  , markovChain, metropolis, hamiltonian
  ) where

import           Control.Monad.Primitive       as X (PrimMonad (..))
import           Linear.Vector
import           Pipes
import qualified Pipes.Prelude                 as P
import           System.Random.MWC             as X (Gen, asGenIO,
                                                     withSystemRandom)
import           System.Random.MWC.Probability


type Transition m t = t -> m t
type MarkovChain m t = forall r. t -> Producer t m r


markovChain :: Monad m => Transition m t -> MarkovChain m t -- (t -> m t) -> t -> Producer t m r
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


hamiltonian
  :: (PrimMonad m, Traversable f, Additive f, Show (f Double))
  => Int                        -- ^ the number of leapfrog steps per sample
  -> Double                     -- ^ epsilon in the leapfrog steps
  -> (f Double -> Double)       -- ^ the potential energy function
  -> (f Double -> f Double)     -- ^ the derivative of the potential energy function
  -> MarkovChain (Prob m) (f Double, Double)
hamiltonian steps eps u du =
  markovChain $ hamiltonianStep steps eps u du


hamiltonianStep
  :: (PrimMonad m, Traversable f, Additive f, Show (f Double))
  => Int
  -> Double
  -> (f Double -> Double)
  -> (f Double -> f Double)
  -> Transition (Prob m) (f Double, Double)
hamiltonianStep steps eps u du start@(q0s, pe0) = do
  p0s <- traverse (const $ normal 0 1) q0s

  let (qfs, pfs) = leapfrog steps eps du (q0s, p0s)
      pef = u qfs
      ke xs = sum $ (\p -> p*p / 2.0) <$> xs
      kef = ke pfs


  if isNaN pef || isNaN kef
    then hamiltonianStep steps eps u du start
    else do
      step <- bernoulli . min 1 . exp $ (ke p0s + pe0) - (kef + pef)

      if step
        then return (qfs, pef)
        else return start


leapfrog
  :: (Floating a, RealFloat a, Additive f)
  => Int
  -> a
  -> (f a -> f a)
  -> (f a, f a)
  -> (f a, f a)
leapfrog steps eps du (q0s, p0s) =
  let p1s = dienan "p1s" <$> ((-eps/2) *^ du q0s ^+^ p0s)

      its =
        flip iterate (q0s, p1s) $ \(qs, ps) ->
          let qs' = dienan "qs'" <$> (eps *^ ps ^+^ qs)
              ps' = dienan "ps'" <$> ((-eps) *^ du qs' ^+^ ps)
          in (qs', ps')

      (q2s, p2s) = its !! (steps-1)
      qfs = eps *^ p2s ^+^ q2s
      pfs = (-eps/2) *^ du qfs ^+^ p2s

  in (qfs, pfs)

  where
    dienan s x = if isNaN x then error s else x

{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards           #-}

module Metropolis
  ( module X
  , metropolis
  , adaptiveMetropolis
  , AMInfo(..)
  , updateAMInfo
  ) where

import           Control.Monad.Trans.State.Strict
import           Linear
import           Numeric.MCMC                     (Transition)

import           Probability                      as X


-- dynamicMetropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t)
--                   => Prob m b -> t b -> Prob m (t b)
-- dynamicMetropolis pr tb =
--   forM tb $ \x -> do
--     r <- normal 0 =<< pr
--     return $ x + r
--
--
-- multibandMetropolis
--   :: (Additive f, InvErf a, Variate a, PrimMonad m, Traversable f)
--   => f a -> f a -> Prob m (f a)
-- multibandMetropolis radii pos = do
--   move <- traverse (normal 0) radii
--   return $ move ^+^ pos

data AMInfo f a =
  AMInfo
    { amt      :: a
    , ampost   :: f a
    , amavgt   :: f a
    , amavgtm1 :: f a
    , amcovt   :: f (f a)
    }

metropolis
  :: ( Additive t1, InvErf b, Variate b, Variate t
     , PrimMonad m, Traversable t1, RealFloat t )
  => b -> (t1 b -> t) -> Transition m (T t (t1 b))
metropolis r llhf = StateT $
  \(T llh xs) -> do
    move <- traverse (const standard) xs
    let prop = (r *^ move) ^+^ xs
        llh' = llhf prop
        prob = exp . min 0 $ llh' - llh

    accept <- if isNaN prob then return False else bernoulli prob

    let next = if accept then prop else xs
        newllh = if accept then llh' else llh

    return ((), T newllh next)


adaptiveMetropolis
  :: ( Additive t1, InvErf b, Ord b, Variate b, Variate t
     , PrimMonad m, Traversable t1, Applicative t1, RealFloat t )
  => b -> t1 (t1 b) -> b -> b -> (t1 b -> t) -> Transition m (T t (AMInfo t1 b))
adaptiveMetropolis t0 cov0 sd eps llhf = StateT $
  \tt@(T llh ami@AMInfo{..}) -> do
    move <- traverse (const standard) ampost
    let cov = if amt > t0 then amcovt else cov0
        prop = (cov !* move) ^+^ ampost
        llh' = llhf prop
        prob = exp . min 0 $ llh' - llh

    accept <- if isNaN prob then return False else bernoulli prob

    -- do not update if we haven't moved and aren't past t0
    if accept || amt > t0
      then return ((), T llh' $ updateAMInfo sd eps ami prop)
      else return ((), tt)



updateAMInfo
  :: (Additive f, Applicative f, Traversable f, Fractional a)
  => a -> a -> AMInfo f a -> f a -> AMInfo f a
updateAMInfo sd eps AMInfo{..} next =
  AMInfo
    tp1
    next
    (((amt *^ amavgt) ^+^ next) ^/ tp1)
    amavgt
    (covtp1 sd eps amt amavgtm1 amavgt amcovt ampost)

  where tp1 = amt + 1


covtp1
  :: (Additive n, Traversable n, Applicative n, Fractional a)
  => a -> a -> a -> n a -> n a -> n (n a) -> n a -> n (n a)
covtp1 sd eps t avgtm1 avgt covt post = x !+! y !-! z !+! w !+! u
  where
    x = (t-1) / t *!! covt
    y = sd *!! outer avgtm1 avgtm1
    z = (sd * (t+1) / t) *!! outer avgt avgt
    w = (sd / t) *!! outer post post
    u = (sd*eps) *!! identity

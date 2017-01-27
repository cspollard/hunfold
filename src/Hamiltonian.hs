{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TypeFamilies              #-}

module Hamiltonian where

import Control.Monad.IO.Class (MonadIO(..))
import           Control.Lens
import           Control.Monad (forM)

import           MarkovChain
import           Probability


-- TODO
-- far too many "negate"s in this function.
hamiltonian :: ( MonadIO m, PrimMonad m, Ord a, Variate a, InvErf a, Ixed (f a),
                TraversableWithIndex (Index (f a)) f, IxValue (f a) ~ a )
            => Int
            -> a
            -> (f a -> a)
            -> (f a -> f a)
            -> T (f a) a
            -> Prob m (T (f a) a)
hamiltonian n eps logLLHF dlogLLHF t@(T qs logLLHqs) = do

  ps <- forM qs (const standard)

  let u = negate . logLLHF
      du = fmap negate . dlogLLHF
      (ps', qs') = iterate (leapfrog eps du) (ps, qs) !! n
      uqs' = u qs'
      prob = min 1 . exp $ negate logLLHqs - uqs' + kinetic ps - kinetic ps'
  move <- bernoulli prob

  return $ if move then T qs' (negate uqs') else t

  where
    kinetic p = (*0.5) . sum $ gzipWith (*) p p


leapfrog :: (Fractional a, Ixed (f a), TraversableWithIndex (Index (f a)) f, IxValue (f a) ~ a)
         => a
         -> (f a -> f a)
         -> (f a, f a)
         -> (f a, f a)
leapfrog eps du (ps, qs) = (ps'', qs')
  where
    update f e x dx = f x $ e * dx
    ps' = gzipWith (update (-) (eps/2)) ps $ du qs
    qs' = gzipWith (update (+) eps) qs ps'
    ps'' = gzipWith (update (-) (eps/2)) ps' $ du qs'


gzipWith
  :: (FunctorWithIndex (Index s) f, Ixed s)
  => (a -> IxValue s -> b) -> f a -> s -> f b
gzipWith f xs ys = imap (\i x -> f x $ ys ^?! ix i) xs

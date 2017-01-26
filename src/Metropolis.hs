{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module Metropolis where

import           Control.Lens
import           Control.Monad (forM)

import           Probability

metropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t)
           => b -> t b -> Prob m (t b)
metropolis r = traverse (\x -> (+x) <$> normal 0 r)


dynamicMetropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t)
                  => Prob m b -> t b -> Prob m (t b)
dynamicMetropolis pr tb =
  forM tb $ \x -> do
    r <- normal 0 =<< pr
    return $ x + r


multibandMetropolis :: ( Ixed s, TraversableWithIndex (Index s) t,
                         InvErf (IxValue s), Variate (IxValue s), PrimMonad m )
                    => s -> t (IxValue s) -> Prob m (t (IxValue s))
multibandMetropolis tr tb = do
  let f i x = (+x) <$> normal 0 (tr ^?! ix i)
  imapM f tb

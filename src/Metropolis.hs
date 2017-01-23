{-# LANGUAGE NoMonomorphismRestriction #-}

module Metropolis where

import Probability

metropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t)
           => b -> t b -> Prob m (t b)
metropolis r = traverse (\x -> (+x) <$> normal 0 r)

{-# LANGUAGE NoMonomorphismRestriction #-}

module Metropolis where

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


-- ASSUMES ZIPABLE TYPES
multibandMetropolis :: (Variate b, InvErf b, PrimMonad m, Traversable t, Applicative t)
                  => t b -> t b -> Prob m (t b)
multibandMetropolis tr tb = do
  rs <- sequence $ normal 0 <$> tr
  return $ (+) <$> tb <*> rs

{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE RecordWildCards           #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}


module Main where

import           Control.Lens
import           Control.Monad        (when)
import           Data.Aeson
import           Data.Aeson.Types     (Parser, parseEither)
import qualified Data.ByteString.Lazy as BS
import           Data.List            (intersperse)
import           Data.Map.Strict      (Map)
import qualified Data.Map.Strict      as M
import           Data.Text            (Text)
import qualified Data.Text            as T
import           Data.Vector          (Vector)
import qualified Data.Vector          as V
import           Linear.Matrix
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
import           Numeric.AD
import           Numeric.MCMC
import           Options.Applicative  hiding (Parser, auto)
import qualified Options.Applicative  as OA
import           System.IO            (BufferMode (..), IOMode (..), hPutStr,
                                       hPutStrLn, hSetBuffering, stdout,
                                       withFile)

import           InMatrix             hiding (transpose, zero)
import           MarkovChain
import           Matrix
import           Model

data InArgs =
    InArgs
        { nburn   :: Int
        , nskip   :: Int
        , nsamps  :: Int
        , outfile :: String
        , infile  :: String
        }

inArgs :: OA.Parser InArgs
inArgs = InArgs
    <$> option OA.auto
      ( long "burn"
      <> help "number of MCMC steps to burn before recording"
      )
    <*> option OA.auto
      ( long "skip"
      <> help "number of MCMC steps to skip between each recording"
      )
    <*> option OA.auto
      ( long "samples"
      <> help "number of samples to record"
      )
    <*> strOption
      ( long "outfile"
      <> help "text file to record to"
      )
    <*> strOption
      ( long "infile"
      <> help "json file to read model from"
      )

opts :: ParserInfo InArgs
opts = info (helper <*> inArgs) fullDesc

main :: IO ()
main = do
  InArgs {..} <- execParser opts

  -- write a line as soon as it comes...
  hSetBuffering stdout LineBuffering

  -- parse the JSON file to Aeson.Values first
  values <- eitherDecode' <$> BS.readFile infile

  -- then try to parsse to our data, Model, and ModelParams
  -- NB: need to give explicit types here so the parser knows what to look for.
  case parseEither parseModel =<< values of
    Left err -> error err
    Right
      ( dataH :: Vector Int
      , model :: Model Double
      , modelparams :: Map Text (ModelParam Double)
      ) -> do

      let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
          start = fmap _mpInitialValue mps
          priors = fmap _mpPrior mps
          variations = fmap _mpVariation mps

          -- I'm not sure why we need an explicit type here.
          -- probably because of the RankNType going on here
          logLH
            :: forall a. (Floating a, Ord a, Mode a, Scalar a ~ Double)
            => Vector a -> a
          logLH =
            toError
            . modelLogPosterior
                dataH
                (fmap auto model)
                (fmap (fmap auto) variations)
                (fmap (ppToFunc . fmap auto) priors)

          gLogLH = grad logLH

      putStrLn ""

      -- find the maximum likelihood starting location
      let xs = take 100 $ conjugateGradientAscent logLH start
          x = last xs

      start' <-
        if any isNaN x || isNaN (logLH x)
          then do
            putStrLn "warning: could not find a likelihood maximum"
            putStrLn "based on your model and initial parameter values."
            putStrLn "using initial values to determine hessian:"
            putStrLn "this could be extremely inefficient."
            putStrLn ""
            -- putStrLn "the 'best' signal histogram appears to be"
            -- print $ signalVars dataH model
            return start
          else
            return x

      putStrLn "starting params:"
      print $ V.zip mpnames start'
      putStrLn "log-likelihood of starting params:"
      print $ logLH start'
      putStrLn ""
      putStrLn "data prediction given starting params:"
      print . toError $ prediction =<< appVars variations start' model
      putStrLn ""
      putStrLn "actual data:"
      print dataH
      putStrLn ""

      -- invert hessian -> covariance matrix
      -- then find the transform from canonical variables
      -- (with mu ~0 and sigma ~1) to the original variables.
      let hess' = hessian (negate . logLH) start'
          cov = toError $ invM hess'
          t = cholM cov
          it = toError $ invM t
          transform' v = (t !* v) ^+^ start'
          invtransform' v' = it !* (v' ^-^ start')

      putStrLn "hessian matrix:"
      print . fst $ toMatrix hess'
      putStrLn ""
      putStrLn "covariance matrix:"
      print . fst $ toMatrix cov
      putStrLn ""
      putStrLn "transformation matrix:"
      print . fst $ toMatrix t

      when
        (anyOf (traverse.traverse) isNaN t)
        $ error "we have NaNs in the transformation matrix; exiting."

      -- need an RNG...
      g <- createSystemRandom

      -- finally, build the chain, metropolis transition, and the MCMC walk
      let c =
            Chain
              (Target (logLH . transform') $ Just (gLogLH . transform'))
              (logLH start')
              (invtransform' start')
              Nothing

          trans = metropolis (1 / fromIntegral nskip)
          walk = takeEvery nskip . LT.drop nburn $ runMC trans c g


      -- write the walk locations to file.
      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
          $ "llh" : V.toList mpnames

        LT.runListT . LT.take nsamps $ do
          Chain{..} <- walk
          LT.liftIO $ do
            hPutStr f $ show chainScore ++ ", "
            hPutStrLn f
              . mconcat . intersperse ", " . V.toList
              $ show <$> transform' chainPosition


everyLT :: Monad m => Int -> (a -> m ()) -> ListT m a -> ListT m a
everyLT n f = go 0
  where
    go m ll = ListT $ do
      c <- next ll
      case c of
        Cons x lll ->
          if m >= n
            then f x >> (return . Cons x $ go 0 lll)
            else return . Cons x $ go (m+1) lll
        Nil -> return Nil


takeEvery :: Monad m => Int -> ListT m a -> ListT m a
takeEvery n l = ListT $ do
  c <- next $ LT.drop (n-1) l
  case c of
    Cons x l' -> return . Cons x $ takeEvery n l'
    Nil       -> return Nil


dropWhileL :: Monad m => (a -> Bool) -> ListT m a -> ListT m a
dropWhileL f l = ListT $ do
  c <- next l
  case c of
    Cons x l' ->
      if f x
        then next $ dropWhileL f l'
        else return c
    Nil -> return Nil


parseModel
  :: (FromJSON a, FromJSON b)
  => Value -> Parser (Vector b, Model a, Map Text (ModelParam a))
parseModel = withObject "error: parseModel was not given a json object" $
  \o -> do
    d <- o .: "Data"
    m <- o .: "Nominal"
    mps <- o .: "ModelVars"

    return (d, m, mps)


signalVars
  :: forall a b. (Ord a, Floating a, Integral b)
  => Vector b
  -> Model a
  -> (Vector Text, Vector (ParamPrior a), Vector (ModelVar a), Vector a)
signalVars dataH model@Model{..} =
  let names = iover traversed (\i _ -> "sigma" <> T.pack (show i)) _mSig
      priors = fmap (const Flat) _mSig
      variations =
        _mSig
          & iover traversed
            (\i _ -> ModelVar Nothing (Just $ iover traversed (\j _ -> if i == j then 1 else 0) _mSig) Nothing Nothing)
      ones = over traverse (const 1) _mSig
      logLH
        :: forall c. (Floating c, Ord c, Mode c, Scalar c ~ a)
        => Vector c -> c
      logLH =
        toError
          . modelLogPosterior
            dataH
            (fmap auto model)
            (fmap (fmap auto) variations)
            (fmap (ppToFunc . fmap auto) priors)
      xs = take 100 $ gradientAscent logLH ones
      best = last xs
  in (names, priors, variations, best)


toError :: Either String c -> c
toError = either error id

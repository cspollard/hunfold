{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}


module Main where

import           Data.Aeson
import           Data.Aeson.Types     (Parser, parseEither)
import qualified Data.ByteString.Lazy as BS
import           Data.Map.Strict      (Map)
import           Data.Text            (Text)
import           Data.Vector          (Vector)
import           Model
import           Options.Applicative  hiding (Parser, auto)
import qualified Options.Applicative  as OA
import           RunModel
import           System.IO            (BufferMode (..), hSetBuffering, stdout)

data InArgs =
  InArgs
    { nsamps  :: Int
    , outfile :: String
    , infile  :: String
    }

inArgs :: OA.Parser InArgs
inArgs =
  InArgs
  <$> option OA.auto
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
    Right (dataH, model, modelparams)
      -> runModel nsamps outfile dataH model modelparams


parseModel
  :: (FromJSON a, FromJSON b)
  => Value -> Parser (Vector b, Model a, Map Text (ModelParam a))
parseModel = withObject "error: parseModel was not given a json object" $
  \o -> do
    d <- o .: "Data"
    m <- o .: "Nominal"
    mps <- o .: "ModelVars"

    return (d, m, mps)


-- everyLT :: Monad m => Int -> (a -> m ()) -> ListT m a -> ListT m a
-- everyLT n f = go 0
--   where
--     go m ll = ListT $ do
--       c <- next ll
--       case c of
--         Cons x lll ->
--           if m >= n
--             then f x >> (return . Cons x $ go 0 lll)
--             else return . Cons x $ go (m+1) lll
--         Nil -> return Nil
--
--
-- takeEvery :: Monad m => Int -> ListT m a -> ListT m a
-- takeEvery n l = ListT $ do
--   c <- next $ LT.drop (n-1) l
--   case c of
--     Cons x l' -> return . Cons x $ takeEvery n l'
--     Nil       -> return Nil
--
--
-- dropWhileL :: Monad m => (a -> Bool) -> ListT m a -> ListT m a
-- dropWhileL f l = ListT $ do
--   c <- next l
--   case c of
--     Cons x l' ->
--       if f x
--         then next $ dropWhileL f l'
--         else return c
--     Nil -> return Nil

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}


module Main where

import           Control.Lens
import           Data.Aeson
import           Data.Aeson.Types     (Parser, parseEither)
import qualified Data.ByteString.Lazy as BS
import           Data.HashMap.Strict  (HashMap)
import           Data.Monoid          ((<>))
import           Data.Text            (Text)
import           Data.Vector          (Vector)
import           Model
import           Options.Applicative  hiding (Parser, auto)
import qualified Options.Applicative  as OA
import           RunModel
import           System.IO            (BufferMode (..), IOMode (..), hPutStrLn,
                                       hSetBuffering, stdout, withFile)

data InArgs =
  InArgs
    { nsamps      :: Int
    , outfile     :: String
    , tablesfile  :: String
    , infile      :: String
    , hamiltonian :: Maybe (Int, Double)
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
    <> help "text file to record samples to"
    )
  <*> strOption
    ( long "tablesfile"
    <> help "text file to write tables to"
    )
  <*> strOption
    ( long "infile"
    <> help "json file to read model from"
    )
  <*> optional
    ( option OA.auto
      ( long "hamiltonian"
      <> help "use hamiltonian sampling rather than metropolis-hastings"
      )
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
  (params, covariances) <-
    case parseEither parseModel =<< values of
      Left err -> error err
      Right (dataH, model, modelparams)
        -> runModel hamiltonian nsamps outfile dataH model (const 0) modelparams

  let uncerts = posteriorMatrices params covariances

  withFile tablesfile WriteMode $ \h -> do
      hPutStrLn h "absolute uncertainties:"
      hPutStrLn h . latextable $ view _1 <$> uncerts

      hPutStrLn h ""
      hPutStrLn h "relative uncertainties:"
      hPutStrLn h . latextable $ view _2 <$> uncerts

      hPutStrLn h ""
      hPutStrLn h "correlations:"
      hPutStrLn h . latextable $ view _3 <$> uncerts



parseModel
  :: (FromJSON a, FromJSON b)
  => Value -> Parser (Vector b, Model a, HashMap Text (ModelParam a))
parseModel = withObject "error: parseModel was not given a json object" $
  \o -> do
    d <- o .: "Data"
    m <- o .: "Nominal"
    mps <- o .: "ModelVars"

    return (d, m, mps)

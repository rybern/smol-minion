{-# LANGUAGE OverloadedLists, DeriveGeneric, TypeApplications #-}
module Plot where

import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import qualified Data.Map as M
import Control.Monad
import Data.Maybe (fromMaybe)
import Data.Csv hiding ((.=))
import qualified Data.ByteString.Lazy.Char8 as BS
import GHC.Generics
import qualified Data.Vector as V

signal :: [Double] -> [(Double,Double)]
signal xs = [ (x,(sin (x*3.14159/45) + 1) / 2 * (sin (x*3.14159/5))) | x <- xs ]

plotLinesFromFile :: FilePath
                  -> (Double -> Double)
                  -> FilePath
                  -> String
                  -> String
                  -> String
                  -> (Double, Double)
                  -> IO ()
plotLinesFromFile inFile fn outFile title xaxis yaxis (lowerY, upperY) = do
  bsContent <- BS.readFile inFile
  let (Right lines) = decode NoHeader bsContent
      ls = V.map (\(title:pts) -> (map (fn . read @Double) pts, title)) lines
      xs = [1..length (fst (V.head ls))]

  toFile def outFile $ do
    layout_title .= title
    forM_ ls $ \(ys, name) -> do
      plot (line name [(zip xs ys)])

plotPoints :: FilePath
           -> [(Double, Double)]
           -> String
           -> IO ()
plotPoints fp points title = toFile def fp $ do
    layout_title .= title
    plot (line "" [points])

plotLines :: (Enum n, Num n, Num m, PlotValue n, PlotValue m)
          => FilePath
          -> String
          -> String
          -> String
          -> (n, n)
          -> (m, m)
          -> [n]
          -> [([m], String)]
          -> IO ()
plotLines fp title xaxis yaxis (lowerX, upperX) (lowerY, upperY) xs lines = toFile def fp $ do
    layout_title .= title
    forM_ lines $ \(ys, name) -> do
      plot (line name [(zip xs ys)])

plotProbLines :: (Enum n, Num n, Num m, PlotValue n, PlotValue m)
              => FilePath
              -> String
              -> String
              -> String
              -> (n, n)
              -> (m, m)
              -> [(M.Map n m, String)]
              -> IO ()
plotProbLines fp title xaxis yaxis (lowerX, upperX) (lowerY, upperY) maps = toFile def fp $ do
    layout_title .= title
    forM_ maps $ \(m, name) -> do
      plot (line name [(interpolate (lowerX, upperX) 0 m)])

interpolate :: (Enum n, Ord n) => (n, n) -> m -> M.Map n m -> [(n, m)]
interpolate (lower, upper) def m = map (\k -> (k, fromMaybe def (M.lookup k m))) [lower..upper]


plotHist :: FilePath
         -> (Double, Double)
         -> Int
         -> [Double]
         -> String
         -> IO ()
plotHist outputFile bounds nBuckets vals title = toFile def outputFile $ do
  let points = histPoints bounds nBuckets vals
  layout_title .= title
  plot (line "" [points])

histPoints :: (Double, Double) -> Int -> [Double] -> [(Double, Double)]
histPoints (lower, upper) nBuckets vals =
  zipWith (\(lower, upper) n -> ((upper + lower) / 2, fromIntegral n / total )) buckets bucketed
  where width = (upper - lower) / fromIntegral nBuckets
        boundaries = [lower, lower+width .. upper]
        buckets = zipWith (,) boundaries (tail boundaries)
        inBuckets = map (\(lower, upper) -> \x -> x >= lower && x <= upper) buckets
        bucketed = map (\pred -> length $ filter pred vals) inBuckets
        total = fromIntegral $ length vals

histFile :: FilePath -> FilePath -> String -> IO ()
histFile inputFile outputFile title = do
  vals <- (map read . lines) <$> readFile inputFile
  plotHist outputFile (0, 1) 300 vals title

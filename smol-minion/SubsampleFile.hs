module SubsampleFile
  ( subsampleFile
  , countFileLines
  ) where

import qualified Data.ByteString.Lazy.Char8 as BS

countFileLines :: FilePath -> IO Int
countFileLines inFile = (return . length . BS.lines) =<< BS.readFile inFile

subsampleFile :: (Int, Int) -> FilePath -> FilePath -> IO ()
subsampleFile bounds inFile outFile = BS.writeFile outFile =<< (return . takeLines bounds) =<< BS.readFile inFile

takeLines :: (Int, Int) -> BS.ByteString -> BS.ByteString
takeLines (start, length)
  | start < 0 = takeLines (0, length + start)
  | otherwise = BS.unlines . take length . drop start . BS.lines

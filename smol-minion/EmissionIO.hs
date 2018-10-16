module EmissionIO where

import Data.Vector (Vector)
import qualified Data.Vector as V
import Control.Monad (liftM)
import Data.Csv
import Sequence.Types
import Sequence (VecMat)
import qualified Data.ByteString.Lazy as BS

encodeMat :: VecMat -> BS.ByteString
encodeMat = encode . V.toList

writeMat :: VecMat -> FilePath -> IO ()
writeMat emissions filepath = BS.writeFile filepath (encodeMat emissions)

decodeMat :: BS.ByteString -> Either String VecMat
decodeMat = decode NoHeader

readMat :: FilePath -> IO (Either String VecMat)
readMat = liftM decodeMat . BS.readFile

decodeViterbi :: BS.ByteString -> Either String (Vector Int)
decodeViterbi = liftM V.head . decode NoHeader

readViterbi :: FilePath -> IO (Either String (Vector Int))
readViterbi = liftM decodeViterbi . BS.readFile

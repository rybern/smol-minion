module EmissionIO where

import Data.Vector (Vector)
import qualified Data.Vector as V
import Control.Monad (liftM)
import Data.Csv
import Sequence.Types
import qualified Data.ByteString.Lazy as BS

type Emissions = Vector (Vector Prob)
type DEmissions = Vector (Vector Double)

toDecimal :: Emissions -> DEmissions
toDecimal = id -- V.map (V.map fromRational)

fromDecimal :: DEmissions -> Emissions
fromDecimal = id --V.map (V.map toRational)

encodeEmissions :: Emissions -> BS.ByteString
encodeEmissions = encode . V.toList . toDecimal

writeEmissions :: Emissions -> FilePath -> IO ()
writeEmissions emissions filepath = BS.writeFile filepath (encodeEmissions emissions)

decodeEmissions :: BS.ByteString -> Either String Emissions
decodeEmissions = liftM fromDecimal . decode NoHeader

readEmissions :: FilePath -> IO (Either String Emissions)
readEmissions = liftM decodeEmissions . BS.readFile

decodeViterbi :: BS.ByteString -> Either String (Vector Int)
decodeViterbi = liftM V.head . decode NoHeader

readViterbi :: FilePath -> IO (Either String (Vector Int))
readViterbi = liftM decodeViterbi . BS.readFile

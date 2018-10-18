{-# LANGUAGE OverloadedLists, ViewPatterns #-}
module EmissionPermutation where

import Data.Map (Map)
import qualified Data.Map as Map
import Data.Vector (Vector)
import qualified Data.Vector as V
import SMoL
import Data.List
import MinION

permuteVector :: Int -> V.Vector String -> V.Vector (V.Vector String)
permuteVector 0 _ = V.singleton []
permuteVector n alpha = V.concatMap (\a -> V.map (a `V.cons`) smaller) alpha
  where smaller = permuteVector (pred n) alpha

buildIndexMap :: (Ord s) => V.Vector s -> Map s Int
buildIndexMap v = V.ifoldl' (\m ix val -> Map.insert val ix m) Map.empty v

kmerString :: V.Vector String -> String
--kmerString = concat . V.toList
kmerString = (\s -> "(" ++ s ++ ")") . intercalate "," . V.toList

minionIndexMap :: Map String Int
minionIndexMap = buildIndexMap . V.map kmerString . permuteVector 5 . V.fromList . map return $ keyOrder

writeEmissionPerm :: Map String Int -> FilePath -> MatSeq String -> IO ()
writeEmissionPerm m file = writeFile file . intercalate "," . map show . V.toList . buildEmissionPerm m

{-# LANGUAGE OverloadedLists, ViewPatterns #-}
module EmissionPermutation where

import Data.Map (Map)
import qualified Data.Map as Map
import Data.Vector (Vector)
import qualified Data.Vector as V
import Sequence
import Data.List
import MinION

permuteVector :: Int -> V.Vector String -> V.Vector (V.Vector String)
permuteVector 0 _ = V.singleton []
permuteVector n alpha = V.concatMap (\a -> V.map (a `V.cons`) smaller) alpha
  where smaller = permuteVector (pred n) alpha

indexMap :: (Ord s) => V.Vector s -> Map s Int
indexMap v = V.ifoldl' (\m ix val -> Map.insert val ix m) Map.empty v

kmerString :: V.Vector String -> String
--kmerString = concat . V.toList
kmerString = (\s -> "(" ++ s ++ ")") . intercalate "," . V.toList

minionIndexMap = indexMap . V.map kmerString . permuteVector 5 . V.fromList . map return $ keyOrder

minionKMerIndex :: (Ord a) => Map a Int -> a -> Int
minionKMerIndex m k = let (Just ix) = Map.lookup k m in ix

buildEmissionPerm :: (Ord a) => Map a Int -> MatSeq a -> V.Vector Int
buildEmissionPerm m = V.map (minionKMerIndex m . fst) . stateLabels

writeEmissionPerm :: Map String Int -> FilePath -> MatSeq String -> IO ()
writeEmissionPerm m file = writeFile file . intercalate "," . map show . V.toList . buildEmissionPerm m

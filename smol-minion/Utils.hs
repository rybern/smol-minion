{-# LANGUAGE OverloadedLists #-}
module Utils where

import Data.Vector (Vector)
import qualified Data.Vector as V
import qualified Data.Map as M
import SMoL
import Data.Function
import Data.List
import EmissionIO
import System.IO.Temp
import System.IO

--normalize :: (Fractional a, Foldable t, Functor t) => t a -> t a
--normalize v = let s = sum v in (/ s) <$> v

groupOn :: (Foldable t, Ord b) => (a -> Maybe b) -> t a -> M.Map b [a]
groupOn f t = foldl
              (\m a -> maybe
                m
                (\b -> M.insertWith (++) b [a] m)
                (f a))
              []
              t

withTempFiles :: Int -> ([(FilePath, Handle)] -> IO r) -> IO r
withTempFiles n = go [] (map show [1..n])
  where go paths (basename:rest) f = withSystemTempFile basename $ \fp h ->
          go (paths ++ [(fp, h)]) rest f
        go paths [] f = f paths

type StateDist s = V.Vector (Prob, (s, StateTag))

pathProbs :: V.Vector (s, StateTag) -> VecMat -> V.Vector (StateDist s)
pathProbs stateLabels = V.map (flip V.zip stateLabels)

sumDist :: (Foldable t, Functor t) => t (Prob, (s, StateTag)) -> Prob
sumDist = sum . fmap fst

distOf :: (Ord a) => [a] -> M.Map a Prob
distOf = normalize . foldl (\m k -> M.insertWith (+) k 1 m) M.empty

avgDist :: (Ord a) => [M.Map a Prob] -> M.Map a Prob
avgDist = normalize . foldl1 (M.unionWith (+))

distOver :: (Ord b) => ((Prob, (s, StateTag)) -> Maybe b) -> StateDist s -> M.Map b Prob
distOver f dist = M.map sumDist (groupOn f dist)

distOverLabel :: (Ord s) => StateDist s -> M.Map s Prob
distOverLabel = distOver (Just . fst . snd)

removeImpossible :: StateDist s -> StateDist s
removeImpossible = V.filter ((/=0) . fst)

maxPosterior :: StateDist s -> (Prob, (s, StateTag))
maxPosterior = maximumBy (compare `on` fst) . V.toList

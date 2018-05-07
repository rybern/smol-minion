{-# LANGUAGE FlexibleInstances, DeriveFoldable, DeriveFunctor #-}
module MinION where

import Sequence
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.List
import Data.Foldable

class Joinable a where
  join :: (Foldable f) => f a -> a

data StateTree a = Collapsed [StateTree a]
                 | Alone a
                 deriving (Eq, Show, Functor, Foldable)

instance Joinable (StateTree a) where
  join = Collapsed . toList

instance Joinable String where
  join = (\s -> "(" ++ s ++ ")") . intercalate "," . toList

keyOrder :: [Char]
keyOrder = ['A', 'C', 'G', 'T']

joinState :: Vector (StateTree a) -> StateTree a
joinState = Collapsed . V.toList

ntTreeToString :: StateTree Char -> String
ntTreeToString (Collapsed sts) = (\s -> "(" ++ s ++ ")") . intercalate "," . map ntTreeToString $ sts
ntTreeToString (Alone a) = [a]

-- 1958 aligned read length, 3698 events
-- [0.529, _]
-- [0.529, _]
{- move counts from read data processes with metrichor:
0: 42603113
1: 13178240
2: 792009
3: 273
4: 877
5: 2658

let counts = [42603113, 13178240, 792009, 273, 877, 2658]
proportional:
[0.7530089080100684,0.23292504732916122,1.3998738360366912e-2,4.825267859809884e-6,1.5500952062466184e-5,4.6980080481225905e-5]
rounded:
[0.753009,0.232925,0.013999,0.000005,0.000016,0.000047]

this gives an expected nt per event of 0.26. this is corroborated by the average number of events per sequence by poretools, but there seems to be a wide variance, maybe event bimodal. many files i've examined seem to be around 0.4.
-}

minionSkipDist :: [Double]
minionSkipDist = [0.753009,0.232925,0.013999,0.000005,0.000016,0.000047]
--minionSkipDist = [0.4589887640449438,0.5334269662921348,7.584269662921348e-3]
--minionSkipDist = [0.44892576515193144,0.5217319632946865,9.780757184460684e-3,9.780757184460684e-3,9.780757184460684e-3]


minion :: (Joinable a) => ProbSeq a -> ProbSeq a
minion = skipDist minionSkipDist
       . collapse undefined join 5

calcAvgEventsPerNT :: Double
calcAvgEventsPerNT = recip . sum $ zipWith (*) [0..] minionSkipDist

avgEventsPerNT :: Double
avgEventsPerNT = 3.827941677480602  --1.639992790844372 --1.822836661546339 -- 0.427 



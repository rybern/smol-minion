{-# LANGUAGE OverloadedLists, ViewPatterns, RecordWildCards #-}
module MinIONRepeat where

import qualified Data.Map as M
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.List
import Control.Monad
import Data.Monoid
import Data.Maybe
import Sequence
import Data.Function

import EmissionIO
import Utils
import Pomegranate

import Plot

maxPriorReps = 20
trueReps = 15

main = do
  print . V.length . stateLabels $ priorSeq
  print "done building"
  agg <- runFeatures 50
  --print agg
  print "done running"
  print (numRun agg)
  plotProbLines "minion_gaussian1.png" "MinIONPlot 1" "Repeats" "Probability" (0, maxPriorReps) (0, 1) [
       --(reformatDists $ naiveDist agg, "Naive estimate")
       (reformatDists $ viterbiNDist agg, "Viterbi")
     , (reformatDists $ postNDist agg, "Gaussian - Max posterior")
     , (reformatDists $ avgPostDist agg, "Gaussian - Average posterior density")
     ]

reformatDists :: M.Map Int Prob -> M.Map Int Double
reformatDists = M.map fromRational . M.mapKeys succ

data SmallRepeatFeatures = SmallRepeatFeatures {
    trueN :: Int
  , viterbiN :: Int
  , postN :: Int
  , naiveN :: Int
  , postDist :: M.Map Int Prob
  } deriving Show

data SmallRepeatAggregate = SmallRepeatAggregate {
    viterbiNDist :: M.Map Int Prob
  , postNDist :: M.Map Int Prob
  , trueNDist :: M.Map Int Prob
  , avgPostDist :: M.Map Int Prob
  , naiveDist :: M.Map Int Prob
  , numRun :: Int
  } deriving Show

aggregateFeatures :: [SmallRepeatFeatures] -> SmallRepeatAggregate
aggregateFeatures fs = SmallRepeatAggregate {
    viterbiNDist = distOf (map viterbiN fs)
  , postNDist = distOf (map postN fs)
  , trueNDist = distOf (map trueN fs)
  , avgPostDist = avgDist (map postDist fs)
  , naiveDist = distOf (map naiveN fs)
  , numRun = length fs
  }

runFeatures :: Int -> IO SmallRepeatAggregate
runFeatures n = (aggregateFeatures . catMaybes . map (smallRepeatFeatures getRepeat naiveCountRepeats)) <$> replicateM n (runHMMSample genSeq priorSeq obsProb)

smallRepeatFeatures :: ((Prob, (s, StateTag)) -> Maybe Int)
                    -> (Vector s -> Int)
                    -> (PomResults s, Vector (s, Int), MatSeq s, MatSeq s)
                    -> Maybe SmallRepeatFeatures
smallRepeatFeatures f naiveFn (PomResults {..}, V.unzip -> (sample, sampleIxs), priorSeq, genSeq) = case keep of
  True -> Just $ SmallRepeatFeatures {
    trueN = trueN
  , viterbiN = viterbiN
  , postN = postN
  , postDist = postDist
  , naiveN = naiveN
  }
  False -> Nothing
  where keep = trueN == trueReps - 1

        tags = pathProbs (stateLabels priorSeq) posterior
        postDist = normalize . M.filterWithKey (\k _ -> k > 5) . distOver f $ V.last tags

        trueN = (fromMaybe 0 $ f (undefined, stateLabels genSeq V.! V.last sampleIxs))
        viterbiN = (fromMaybe 0 $ f (undefined, stateLabels priorSeq V.! V.last viterbi))
        postN = fst . maximumBy (compare `on` snd) . M.toList $ postDist
        naiveN = naiveFn sample

naiveCountRepeats :: Vector String -> Int
naiveCountRepeats = const 1

obsProb :: (Eq s) => MatSeq s -> s -> IO (Vector Prob)
obsProb seq i = return
              . normalize
              . V.map (\(i', _) -> if i == i'
                                   then advantage / denom
                                   else 1 / denom)
              . stateLabels $ seq
  where n = 3 ^ 3
        advantage = 5
        denom = n + advantage - 1

getRepeat :: ((Prob, (s, StateTag)) -> Maybe Int)
getRepeat (_, (_, StateTag 0 [StateTag 0 (last -> StateTag 1 [StateTag n _])])) = Just n
getRepeat _ = Nothing

genSeq = buildMatSeq $ minionSeq (repeatSequence trueReps)
--priorSeq = buildMatSeq $ minionSeq (uniformDistRepeat maxPriorReps)
priorSeq = buildMatSeq $ minionSeq (finiteDistRepeat gaussian15)
gaussian15 = (normalize $ replicate 11 0 ++ gaussian9)
gaussian9 = [0.000088,0.002289,0.023205,0.092566,0.146634,0.092566,0.023205,0.002289,0.000088]

minionSeq :: (ProbSeq String -> ProbSeq String) -> ProbSeq String
minionSeq reps = skipDist skipD
               . collapse undefined concat 3
               $ series' [
                   geomUniform
                 , reps . series' . map (\c -> state [c]) $ "AGTC"
                 ]
  where alphabet = map state ["G", "A", "T", "C"]
        geomUniform = geometricRepeat 0.95
          . uniformDistOver
          $ alphabet

skipD = [0.3, 0.3, 0.2, 0.2]

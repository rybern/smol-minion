{-# LANGUAGE OverloadedLists, ViewPatterns, RecordWildCards #-}
module SmallRepeat where

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

maxPriorReps = 50
minPriorReps = 20
trueReps = 35

main = do
  let uniformPriorSeq = buildMatSeq $ andThen
        (repeatSequence' minPriorReps periodSeq)
        (uniformDistRepeat' (maxPriorReps - minPriorReps) periodSeq)

  --(Just s) <- (smallRepeatFeatures getPriorRepeat getGenRepeat naiveCountRepeats) <$> (runHMMSample genSeq uniformPriorSeq obsProb)
  --print $ (postDist s M.! (trueReps -1))
  --print $ trueReps -1
  --print $ postDist s
  --print $ ((postDist s M.! (trueReps - 1)) /= 0)
  agg <- runFeatures genSeq uniformPriorSeq 200
  print "done running"
  --print (avgPostDist agg)
  print (numRun agg)
  let gaussianPriorSeq = buildMatSeq $ andThen
        (repeatSequence' (minPriorReps) periodSeq)
        (finiteDistRepeat gaussian29Tighter periodSeq)
  agg2 <- runFeatures genSeq gaussianPriorSeq 200
  print "done running"
  print (numRun agg2)
  plotProbLines "test3.png" "Simulation: Simple [1,2] repeats" "Repeats" "Probability" (0, maxPriorReps) (0, 1) [
       (reformatDists $ naiveDist agg, "Naive estimate")
     , (reformatDists $ viterbiNDist agg, "Uniform - Viterbi")
     , (reformatDists $ postNDist agg, "Uniform - Max posterior")
     , (reformatDists $ avgPostDist agg, "Uniform - Average posterior density")
     , (reformatDists $ postNDist agg2, "Gaussian - Max posterior")
     , (reformatDists $ avgPostDist agg2, "Gaussian - Average posterior density")
     ]

gaussian9 :: [Prob]
gaussian9 = [0.000088,0.002289,0.023205,0.092566,0.146634,0.092566,0.023205,0.002289,0.000088]
gaussian29 :: [Prob]
gaussian29 = [0.000128,0.00022,0.000362,0.000573,0.000871,0.001272,0.001785,0.002407,0.003119,0.003884,0.004647,0.005343,0.005903,0.006266,0.006393,0.006266,0.005903,0.005343,0.004647,0.003884,0.003119,0.002407,0.001785,0.001272,0.000871,0.000573,0.000362,0.00022,0.000128]
gaussian29Tighter :: [Prob]
gaussian29Tighter = [0.000005,0.000014,0.000038,0.000096,0.000224,0.000484,0.000964,0.00177,0.002999,0.004684,0.006746,0.00896,0.010973,0.012393,0.012905,0.012393,0.010973,0.00896,0.006746,0.004684,0.002999,0.00177,0.000964,0.000484,0.000224,0.000096,0.000038,0.000014,0.000005]
gaussian29Tighter' :: [Prob]
gaussian29Tighter' = [0,0,0,0,0,0.000002,0.000015,0.000097,0.000474,0.001825,0.005495,0.012949,0.023883,0.034481,0.038972,0.034481,0.023883,0.012949,0.005495,0.001825,0.000474,0.000097,0.000015,0.000002,0,0,0,0,0]
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
  , avgPostDist :: M.Map Int Prob
  , naiveDist :: M.Map Int Prob
  , numRun :: Int
  } deriving Show

aggregateFeatures :: [SmallRepeatFeatures] -> SmallRepeatAggregate
aggregateFeatures fs = SmallRepeatAggregate {
    viterbiNDist = distOf (map viterbiN fs)
  , postNDist = distOf (map postN fs)
  , avgPostDist = avgDist (map postDist fs)
  , naiveDist = distOf (map naiveN fs)
  , numRun = length fs
  }

runFeatures :: MatSeq Int -> MatSeq Int -> Int -> IO SmallRepeatAggregate
runFeatures genSeq priorSeq n = (aggregateFeatures . catMaybes . map (smallRepeatFeatures getPriorRepeat getGenRepeat naiveCountRepeats)) <$> replicateM n (runHMMSample genSeq priorSeq obsProb)

getPriorRepeat :: ((Prob, (s, StateTag)) -> Maybe Int)
getPriorRepeat (_, (_, StateTag 1 [StateTag n _])) = Just (n + minPriorReps)
getPriorRepeat _ = Nothing

getGenRepeat :: ((Prob, (s, StateTag)) -> Maybe Int)
getGenRepeat (_, (_, StateTag n _)) = Just n

smallRepeatFeatures :: ((Prob, (s, StateTag)) -> Maybe Int)
                    -> ((Prob, (s, StateTag)) -> Maybe Int)
                    -> (Vector s -> Int)
                    -> (PomResults s, Vector (s, Int), MatSeq s, MatSeq s)
                    -> Maybe SmallRepeatFeatures
smallRepeatFeatures priorFn genFn naiveFn (PomResults {..}, V.unzip -> (sample, sampleIxs), priorSeq, genSeq) =
  case keep of
  True -> Just $ SmallRepeatFeatures {
    trueN = trueN
  , viterbiN = viterbiN
  , postN = postN
  , postDist = postDist
  , naiveN = naiveN
  }
  False -> Nothing
  where keep = (trueN == trueReps - 1) && (postDist M.! (trueReps - 1) /= 0) -- (postDist M.! trueReps-1 > postDist M.! trueReps-2 || postDist M.! trueReps-1 > postDist M.! trueReps)

        tags = pathProbs (stateLabels priorSeq) posterior
        postDist = normalize . M.filterWithKey (\k _ -> k > 5) . distOver priorFn $ V.last tags

        trueN = (fromMaybe 0 $ genFn (undefined, stateLabels genSeq V.! V.last sampleIxs))
        viterbiN = (fromMaybe 0 $ priorFn (undefined, stateLabels priorSeq V.! V.last viterbi))
        postN = fst . maximumBy (compare `on` snd) . M.toList $ postDist
        naiveN = naiveFn sample

naiveCountRepeats :: Vector Int -> Int
naiveCountRepeats = go . V.toList
  where go (1:2:rest) = 1 + go rest
        go (1:rest) = 1 + go rest
        go (2:rest) = 1 + go rest
        go [] = -1

obsProb :: (Eq s) => MatSeq s -> s -> IO (Vector Prob)
obsProb seq i = return . normalize . V.map (\(i', _) -> if i == i' then 1.0 else 0.0) . stateLabels $ seq

genSeq = buildMatSeq $ repeatSequence trueReps periodSeq

periodSeq :: ProbSeq Int
periodSeq = series' . map (\v -> andThen (state v) skipDSeq) $
  [ 1, 2 ]

skipD :: [Prob]
skipD = [0.5, 1 - head skipD]

skipDSeq :: ProbSeq a
skipDSeq = finiteDistRepeat skipD $ skip 1

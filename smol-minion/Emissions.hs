{-# LANGUAGE OverloadedLists #-}
module Sequence.Matrix.Emissions
  (
    stateSequenceProbability
  , endAfterStepsProbability
  , sequencePrefixProbability
  , sequenceSuffixProbability
  ) where

import Sequence.Matrix.Types
import Sequence.Matrix.ProbSeqMatrixUtils
import Sequence.Matrix.Operations.AndThen
import Sequence.Matrix.Operations.Deterministic
import Control.Monad
import qualified Data.Vector as V
import qualified Math.LinearAlgebra.Sparse as M

endAfterStepsProbability :: MatSeq s -> Int -> Int -> Prob
endAfterStepsProbability seq n end = endAfterStepsProbability' (squareMain, ends) n end 1
  where (main, ends) = splitEnds (trans seq)
        squareMain = addStartColumn main

endAfterStepsProbability' :: (Trans, Trans) -> Int -> Int -> M.Index -> Prob
endAfterStepsProbability' _ 0 _ _ = 0.0
endAfterStepsProbability' (main, ends) n end current = thisRowEnd + sum nextRowEnds
  where thisRowEnd = M.row ends current M.! (end+1)
        nexts = M.row main current
        nextRowEnds = map (\(ix, p) -> p * endAfterStepsProbability' (main, ends) (n-1) end ix)
                    . tail
                    . M.vecToAssocList
                    $ nexts

sequencePrefixProbability :: (Eq s) => V.Vector s -> MatSeq s -> Prob
sequencePrefixProbability path seq = sum . pathProbs (getNormalTransWithEnds seq) . V.toList . V.map (stateIxs seq) $ path

sequenceSuffixProbability :: (Eq s) => Int -> (V.Vector s, Int) -> MatSeq s -> Prob
sequenceSuffixProbability skipped (seq, nSkip) m =
  stateSequenceProbability (seq, nSkip) $ skip skipped `andThen` m

stateSequenceProbability :: (Eq s) => (V.Vector s, Int) -> MatSeq s -> Prob
stateSequenceProbability (path, skip) seq = sum . pathProbs (getNormalTransWithEnds seq) . stateSequenceIxs seq $ (path, skip)

stateSequenceIxs :: (Eq s) => MatSeq s -> (V.Vector s, Int) -> [V.Vector M.Index]
stateSequenceIxs seq (path, skip) = V.toList . addEnd . V.map (stateIxs seq) $ path
  where endIx = V.length (stateLabels seq) + 2
        addEnd = (`V.snoc` [endIx + skip])

stateIxs :: (Eq s) => MatSeq s -> s -> V.Vector M.Index
stateIxs (MatSeq {stateLabels = stateLabels}) label =
  flip V.imapMaybe stateLabels $ \ix (label', _) ->
                                   if label' == label
                                   then Just (ix + 2)
                                   else Nothing

pathProbs :: Trans -> [V.Vector M.Index] -> V.Vector Prob
pathProbs = pathProbs' 1

pathProbs' :: M.Index -> Trans -> [V.Vector M.Index] -> V.Vector Prob
pathProbs' _ _ [] = [1.0]
pathProbs' current trans (nexts:rest) = do
  next <- nexts

  let transProb = trans M.# (current, next)
  guard $ transProb > 0

  nextPathProb <- pathProbs' next trans rest

  return $ transProb * nextPathProb

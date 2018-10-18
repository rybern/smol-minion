{-# LANGUAGE DeriveFunctor, OverloadedLists, TupleSections, BangPatterns, RecordWildCards, ViewPatterns #-}
module GeneralizedSNP where

import SNP
import Data.List
import Data.Function
import Data.Set (Set)
import qualified Data.Set as Set
import Data.Vector (Vector)
import qualified Data.Vector as V
import SMoL
import MinION

import SMoL.Matrix.ProbSeqMatrixUtils
import SMoL.Matrix.IO.StateLabels
import Data.Csv
import qualified Data.ByteString.Lazy.Char8 as BS
--import Math.LinearAlgebra.Sparse.Matrix hiding (trans)
import SparseMatrix hiding (trans)

data GenSNPState = NoiseKey Int
                 | LeftFlank Int
                 | Ref
                 | Alt
                 | RightFlank Int
                 | Untracked
                 deriving (Eq, Show)

flankSize = 50
noiseSize = 4
softFlankSize = 20

genSite :: Int -> Site GenSNPState
genSite flankSize = Site {
    pos = 0
  , alleles = [(Ref, 0.5), (Alt, 0.5)]
  , leftFlank = map LeftFlank [0..flankSize-1]
  , rightFlank = map RightFlank [0..flankSize-1]
  }

specifySNPState :: [a] -> Site a -> GenSNPState -> a
specifySNPState keyOrder _ (NoiseKey i) = keyOrder !! i
specifySNPState _ site (LeftFlank i) = leftFlank site !! i
specifySNPState _ site Ref = fst $ alleles site !! 0
specifySNPState _ site Alt = fst $ alleles site !! 1
specifySNPState _ site (RightFlank i) = rightFlank site !! i

  {- can make this fully usable outside of flank overlap by adding exactly one frame of noise -}

  -- roughly 30 seconds
genMatSeq :: MatSeq (StateTree GenSNPState)
(genMatSeq, genMatIxs@[[refIxs, altIxs]]) = snpsMatSeq Alone (map NoiseKey [0..noiseSize-1]) (-softFlankSize) softFlankSize [genSite flankSize]

genMatSeq' flankSize = snpsMatSeq'
  Alone
  True
  (Alone Untracked)
  (- (length (leftFlank (genSite flankSize))) - 1)
  (length (rightFlank (genSite flankSize)) + 1)
  [(genSite flankSize)]

isRef Ref = True
isRef _ = False
isAlt Alt = True
isAlt _ = False

genMatSeqIxs' :: Int -> (GenSNPState -> Bool) -> Vector Int
genMatSeqIxs' flankSize pred =
    V.map fst
  . V.filter (any pred . snd)
  . V.imap (\i v -> (i, stateLabel v))
  $ stateLabels (genMatSeq' flankSize)

genMatSeqAltIxs' flankSize = genMatSeqIxs' flankSize isAlt
genMatSeqRefIxs' flankSize = genMatSeqIxs' flankSize isRef

endNoiseMatSeq :: String -> Int -> MatSeq [NT]
endNoiseMatSeq token ((* avgEventsPerNT) . fromIntegral -> expected) = buildMatSeq
  $ geometricRepeat (1 / (1+expected)) (symbol token)

snpsNTMatSeq'' :: Int -> String -> Int -> Int -> [Site NT] -> (MatSeq [NT], [[Vector Int]])
snpsNTMatSeq'' flankSize token readStart readEnd sites =
  (matSeq, ixs)
  where liftSite :: Site a -> Site [a]
        liftSite = fmap (\c -> [c])
        --sites = map liftSite sites'
        gms = genMatSeq' flankSize
        specified = map (specifyGenMatSeqNT' flankSize gms) sites
        regions :: [SNPCallerRegion [NT]]
        regions = snpRegions readStart readEnd (map (fmap (\c -> [c])) sites)
        nStates matSeq = V.length (stateLabels matSeq)
        -- fold with (compiled sites, compiled sites+regions, site offsets, current offset)
        (_, reverse -> seqs, reverse -> siteIxs, _) = foldl' takeRegion (specified, [], [], 0) regions
        takeRegion (spec, matSeqs, sofar, n) (Noise exp) =
          let ms = endNoiseMatSeq token exp
          in ( spec
             , ms : matSeqs
             , sofar
             , nStates ms + n )
        takeRegion (spec, matSeqs, sofar, n ) noise@(Flank _) = ( spec, matSeqs, sofar, n )
        takeRegion (next:rest, matSeqs, sofar, n) noise@(SNP _) =
          ( rest, next : matSeqs, n : sofar, nStates next + n )

        matSeq = buildMatSeq . series $ map matrixForm seqs

        ixs = map (\siteOffset -> [ (siteOffset+) <$> genMatSeqRefIxs' flankSize
                                  , (siteOffset+) <$> genMatSeqAltIxs' flankSize]) siteIxs

specifyGenMatSeqNT' :: Int -> MatSeq (StateTree GenSNPState)
                   -> Site Char
                   -> MatSeq [NT]
specifyGenMatSeqNT' flankSize genMatSeq site =
  mapStates ntTreeToString $ specifyGenMatSeq (genMatSeq' flankSize) genMatSeqIxSets' keyOrder site

setStateContains matSeq nt =
  Set.fromList . map fst . filter (any (== nt) . stateLabel . snd) . zip [1..] . V.toList . stateLabels $ matSeq

genMatSeqIxSets' = let refSet = setStateContains (genMatSeq' flankSize) Ref
                       altSet = setStateContains (genMatSeq' flankSize) Alt
                   in (refSet, altSet, refSet `Set.union` altSet)

genMatSeqIxSets = let refSet = setStateContains genMatSeq Ref
                      altSet = setStateContains genMatSeq Alt
                  in (refSet, altSet, refSet `Set.union` altSet)

genMatSeqAlleleIxs :: [(Int, Int)]
genMatSeqAlleleIxs = allelePairs (genMatSeq, [[refIxs, altIxs]])

{-
You should only need to find the ix pairs once per genSite

test this by inspection of tags for genMatSeq
-}
specifyGenMatSeq :: MatSeq (StateTree GenSNPState)
                 -> (Set Int, Set Int, Set Int)
                 -> [a]
                 -> Site a
                 -> MatSeq (StateTree a)
specifyGenMatSeq genMs sets keyOrder site = updateProbs . updateStateLabels $ genMs
  where updateStateLabels = mapStates (specifySNPState keyOrder site <$>)
        updateProbs = updateSinkRatio sets (let ((_, p):_) = alleles site in p)

updateSinkRatio :: (Set Int, Set Int, Set Int)
                -> Prob -- ref prob
                -> MatSeq s
                -> MatSeq s
updateSinkRatio (refSet, altSet, eitherSet) p ms = ms { trans = mapWithIxs update (trans ms) }
  where update (from, to) val =
          if (from - 1) `Set.member` eitherSet
          then val
          else case (to `Set.member` refSet, to `Set.member` altSet) of
               (True, False) -> p * val * 2
               (False, True) -> (1-p) * val * 2
               (False, False) -> val
               (True, True) -> error $ "Index " ++ show to ++ " is in both altSet and refSet!"

-- mapWithIxs :: (Num a, Eq a) => ((M.Index, M.Index) -> a -> a) -> M.SparseMatrix a -> M.SparseMatrix a
-- use symmetry, don't need pairs
-- let f (p, q) (a, b) = (q * 2 * a, p * 2 * b)

allelePairs :: (MatSeq (StateTree GenSNPState), [[V.Vector Int]]) -> [(Int, Int)]
allelePairs (genMs, [[refIxs, altIxs]]) = ixPairs
  where ixStateTag ix = (stateTag . (stateLabels genMs V.!) $ ix, ix)
        ixStateTags = V.toList . V.map ixStateTag
        tagPairs = pairOff (tagPair `on` fst) (ixStateTags refIxs) (ixStateTags altIxs)
        ixPairs = map (\((_, i1), (_, i2)) -> (i1, i2)) tagPairs

tagPair :: StateTag -> StateTag -> Bool
tagPair (StateTag n1 [StateTag 0 []]) (StateTag n2 [StateTag 0 []]) = True
tagPair (StateTag _ []) _ = False
tagPair _ (StateTag _ []) = False
tagPair (StateTag n1 l1) (StateTag n2 l2) = n1 == n2 && and (zipWith tagPair l1 l2)

pairOff :: (Eq b) => (a -> b -> Bool) -> [a] -> [b] -> [(a, b)]
pairOff pred (a:as) bs = case find (pred a) bs of
  Nothing -> pairOff pred as bs
  Just b' -> (a, b') : pairOff pred as (delete b' bs)
pairOff pred [] bs = []

specifyGenMatSeqNT :: MatSeq (StateTree GenSNPState)
                   -> Site Char
                   -> (MatSeq [NT], [[Vector Int]])
specifyGenMatSeqNT genMatSeq site =
  ( mapStates ntTreeToString $ specifyGenMatSeq genMatSeq genMatSeqIxSets keyOrder site
  , genMatIxs)

exSite p = Site {
    pos = 0
  , leftFlank =  "abcdefghij"
  , alleles = [('0', p), ('1', 1-p)]
  , rightFlank = "ABCDEFGHIJ"
  }
site = exSite 0.35

check ms = putStrLn $ show (trans ms # (100, 100))

(ms, ixs) = specifyGenMatSeqNT genMatSeq site


{-

testing code


testMatSeq :: (Joinable b, Eq b) => (a -> b) -> [a] -> [Site a] -> MatSeq b
testMatSeq toJoinable (map toJoinable -> keyOrder) (map (toJoinable <$>) -> sites) =
  buildMatSeq . minion . regionsProbSeq keyOrder $ snpRegions sites

testA :: MatSeq (StateTree GenSNPState)
testA = testMatSeq Alone (map NoiseKey [0..noiseSize-1]) [genSite]
testA' = mapStates (specifySNPState keyOrder site <$>) testA
testA'' = updateSinkRatio genMatSeqIxSets (let ((_, p):_) = alleles site in p) testA'
testA''' = mapStates ntTreeToString $ testA''

testC :: MatSeq [NT]
testC = fst $ specifyGenMatSeqNT site

testB :: MatSeq [NT]
testB = testMatSeq return keyOrder [site]

--specifyGenMatSeq

--(testA, aIxs) = specifyGenMatSeqNT (exSite 0.25)
--(testB, bIxs) = snpsNTMatSeq [exSite 0.35]

eqMat :: Trans
eqMat = (trans testA''' - trans testB)
eq = isZeroMx eqMat

eqMatList :: [((Int, Int), Prob)]
eqMatList = toAssocList eqMat

{-
setStateContains matSeq nt =
  Set.fromList . map fst . filter (any (== nt) . fst . snd) . zip [0..] . V.toList . stateLabels $ matSeq
set0 = setStateContains testA '0'
set1 = setStateContains testA '1'
-}

inSet nt =
  filter (\e -> Set.member (fromIntegral $ e - 1) (setStateContains testA nt)) . map (fst . fst) $ toAssocList eqMat

-}

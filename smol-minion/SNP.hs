{-# LANGUAGE DeriveFunctor, TupleSections, BangPatterns, RecordWildCards, ViewPatterns #-}
module SNP where

import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.Set (Set)
import qualified Data.Set as S
import Data.Maybe
import Sequence
import MinION

type NT = Char

data Site a = Site {
    pos :: Int
  , alleles :: [(a, Prob)]
  , leftFlank :: [a]
  , rightFlank :: [a]
  } deriving (Show, Functor)

data SNPCallerRegion a = SNP [(a, Prob)]
                       | Flank [a]
                       | Noise Int
                       deriving Show

regionProbSeq :: [a] -> SNPCallerRegion a -> ProbSeq a
regionProbSeq keyOrder (SNP nts) = finiteDistOver . map (\(nt, p) -> (state nt, p)) $ nts
regionProbSeq keyOrder (Flank nts) = series . map state $ nts
regionProbSeq keyOrder (Noise expected) = geometricRepeat (1 / (1 + fromIntegral expected)) noise
  where noise = uniformDistOver (map state keyOrder)

regionsProbSeq :: [a] -> [SNPCallerRegion a] -> ProbSeq a
regionsProbSeq keyOrder = series . map (regionProbSeq keyOrder)

-- TODO: estimate starting and ending positions
-- need to start and end with noise of at least length 2.
--  geometric with length 1 does not accept strings of longer than 1
snpRegions :: Int -> Int -> [Site a] -> [SNPCallerRegion a]
snpRegions readStart readEnd sites =
  let end = last sites
      first = head sites
  in flankBetween True readStart Nothing (Just $ leftFlank first) (pos first)
     ++ go sites
     ++ flankBetween True (pos end) (Just $ rightFlank end) Nothing readEnd
  where go (a:b:rest) = siteSNP a : siteFlankBetween a b ++ go (b:rest)
        go [a] = [siteSNP a]
        go [] = []

siteSNP :: Site a -> SNPCallerRegion a
siteSNP a = SNP (alleles a)

siteFlankBetween :: Site a -> Site a -> [SNPCallerRegion a]
siteFlankBetween a b = flankBetween False (pos a) (Just $ rightFlank a) (Just $ leftFlank b) (pos b)

  -- assume a is earlier than b!
flankBetween :: Bool -> Int -> Maybe [a] -> Maybe [a] -> Int -> [SNPCallerRegion a]
flankBetween forceNoise aLoc maybeAFlank maybeBFlank bLoc = regions
  where regions = if noiseDist' > 0
                  then catMaybes [Flank <$> maybeAFlank, Just $ Noise noiseDist', Flank <$> maybeBFlank]
                  else [Flank (take dist aFlank ++ drop (length bFlank - (dist - length aFlank)) bFlank)]
        aFlank = fromMaybe [] maybeAFlank
        bFlank = fromMaybe [] maybeBFlank
        dist = abs $ bLoc - aLoc - 1
        sumFlankLen = length aFlank + length bFlank
        noiseDist' = if forceNoise
                     then max 2 $ noiseDist aLoc maybeAFlank maybeBFlank bLoc
                     else noiseDist aLoc maybeAFlank maybeBFlank bLoc

noiseDist :: Int -> Maybe [a] -> Maybe [a] -> Int -> Int
noiseDist aLoc maybeAFlank maybeBFlank bLoc = if noiseDist' > 0
                                              then max noiseDist' 2
                                              else noiseDist'
  where aFlank = fromMaybe [] maybeAFlank
        bFlank = fromMaybe [] maybeBFlank
        dist = abs $ bLoc - aLoc - 1
        sumFlankLen = length aFlank + length bFlank
        noiseDist' = dist - sumFlankLen


siteDistance :: Site a -> Site a -> Int
siteDistance a b = pos b - pos a - 1

ntComplement :: NT -> NT
ntComplement 'A' = 'T'
ntComplement 'T' = 'A'
ntComplement 'G' = 'C'
ntComplement 'C' = 'G'

siteComplement :: Site NT -> Site NT
siteComplement s@(Site {..}) = s {
    alleles = map (\(n, p) -> (ntComplement n, p)) alleles
  , leftFlank = map ntComplement leftFlank
  , rightFlank = map ntComplement rightFlank
  }

{-
snpsMatSeqRC :: [NT] -> [Site NT] -> (MatSeq [NT], [[Vector Int]])
snpsMatSeqRC keyOrder (map siteComplement -> sites) = (matSeq, ixPairs)
  where regions = snpRegions sites
        matSeq = buildMatSeq . minion . reverseSequence . regionsProbSeq keyOrder $ regions

        isSNP (SNP vars) = Just (length vars)
        isSNP _ = Nothing
        snpIxs = catMaybes . zipWith (\i m -> (,i) <$> m) [0..] . map isSNP $ regions

        tags = V.map snd . stateLabels $ matSeq
        tagSNP snp allele (StateTag 0 [StateTag _ [StateTag 0 seqStates]]) = flip any seqStates $ \seqState ->
          case seqState of
            (StateTag n [StateTag x _])  -> x == allele && n == snp
            --(StateTag 1 [StateTag 0 [StateTag n [StateTag x _]]]) -> x == allele && n == snp
            _ -> False
        tagSNP _ _ _ = False

        allIxs = V.enumFromTo 0 (V.length tags - 1)
        ixPairs = flip map snpIxs $ \(nAlleles, snp) ->
          flip map [0..nAlleles-1] $ \allele ->
          V.filter (\i -> tagSNP snp allele $ tags V.! i) allIxs
-}

snpsNTMatSeq :: Int -> Int -> [Site NT] -> (MatSeq [NT], [[Vector Int]])
snpsNTMatSeq = snpsMatSeq (\a -> [a]) keyOrder

-- make generic module, type for each location in snp site region, make site-swapping half of the specifying function

regionProbSeq' :: a -> SNPCallerRegion a -> ProbSeq a
regionProbSeq' token (SNP nts) = finiteDistOver . map (\(nt, p) -> (state nt, p)) $ nts
regionProbSeq' token (Flank nts) = series . map state $ nts
regionProbSeq' token (Noise expected) = geometricRepeat (1 / (1 + avgEventsPerNT * fromIntegral expected)) (state token)

isNoise (Noise _) = True
isNoise _ = False

regionsProbSeq' :: (Joinable a) => a -> [SNPCallerRegion a] -> ProbSeq a
regionsProbSeq' token = series . regionsProbSeq
  where regionsProbSeq regions = case break isNoise regions of
          ([], []) -> []
          ([], (noise@(Noise _)):rest) -> regionProbSeq' token noise : regionsProbSeq rest
          (nonNoise, rest) -> (minion . series . map (regionProbSeq' token) $ nonNoise) : regionsProbSeq rest


snpsMatSeq' :: (Joinable b, Eq b) => (a -> b) -> Bool -> b -> Int -> Int -> [Site a] -> MatSeq b
snpsMatSeq' toJoinable noNoise token readStart readEnd (map (toJoinable <$>) -> sites) =
  matSeq
  where regions = let regions' = snpRegions readStart readEnd sites
                  in if noNoise
                     then tail . init $ regions'
                     else regions'
        matSeq = buildMatSeq $ regionsProbSeq' token regions

snpsMatSeq :: (Joinable b, Eq b) => (a -> b) -> [a] -> Int -> Int -> [Site a] -> (MatSeq b, [[Vector Int]])
snpsMatSeq toJoinable (map toJoinable -> keyOrder) readStart readEnd (map (toJoinable <$>) -> sites) =
  (matSeq, ixPairs)
  where regions = snpRegions readStart readEnd sites
        matSeq = buildMatSeq . minion . regionsProbSeq keyOrder $ regions

        isSNP (SNP vars) = Just (length vars)
        isSNP _ = Nothing
        snpIxs = catMaybes . zipWith (\i m -> (,i) <$> m) [0..] . map isSNP $ regions

        tags = V.map stateTag . stateLabels $ matSeq
        tagSNP snp allele (StateTag 0 [StateTag _ seqStates]) = flip any seqStates $ \seqState ->
          case seqState of
            (StateTag n [StateTag x _])  -> x == allele && n == snp
            --(StateTag 1 [StateTag 0 [StateTag n [StateTag x _]]]) -> x == allele && n == snp
            _ -> False
        tagSNP _ _ _ = False

        allIxs = V.enumFromTo 0 (V.length tags - 1)
        ixPairs = flip map snpIxs $ \(nAlleles, snp) ->
          flip map [0..nAlleles-1] $ \allele ->
          V.filter (\i -> tagSNP snp allele $ tags V.! i) allIxs


{-
snpMatSeq :: Site a -> (MatSeq [a], Vector Int, Vector Int)
snpMatSeq (Site {..}) = (matSeq, refIxs, altIxs)
  where maf = snd $ alleles !! 1
        alt = fst $ alleles !! 1
        ref = fst $ alleles !! 0
        matSeq = buildMatSeq . minion $
          series [noise, series [left, eitherOr maf altState refState, right] , noise]
        altState = state [alt]
        refState = state [ref]
        left = series . map (state . return) $ leftFlank
        right = series . map (state . return) $ rightFlank
        noise = geometricRepeat 0.99 (uniformDistOver (map state keyOrder))

        tags = V.map snd . stateLabels $ matSeq
        tagSNP snp (StateTag 0 [StateTag 0 seqStates]) = flip any seqStates $ \seqState ->
          case seqState of
            (StateTag 1 [StateTag 1 [StateTag x _]]) -> x == snp
            _ -> False
        tagSNP _ _ = False

        allIxs = V.enumFromTo 0 (V.length tags - 1)
        altIxs = V.filter (\i -> tagSNP 0 $ tags V.! i) allIxs
        refIxs = V.filter (\i -> tagSNP 1 $ tags V.! i) allIxs
-}

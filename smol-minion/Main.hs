{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedLists, ViewPatterns #-}
module Main where

import Data.Map (Map)
import qualified Data.Map as Map
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IntMap
import Data.List
import qualified Data.ByteString.Lazy.Char8 as BS
import Data.Csv
import Control.Monad
import Control.Monad.Loops
import Data.Monoid
import Data.Maybe
import Data.Char
import SMoL
import Data.Function
import System.Environment
import System.IO
import System.IO.Temp

import Plot
import SNP
import GeneralizedSNP
import EmissionPermutation
import EmissionIO
import Utils
import MinION
import SubsampleFile
import SMoL.Inference.SHMM
import qualified SMoL.Inference as IF
import SMoL.Tags
import Data.Hashable

import SMoL.Matrix.ProbSeqMatrixUtils

--import Math.LinearAlgebra.Sparse.Matrix hiding (trans)
import SparseMatrix hiding (trans)

setupSiteAnalysis :: String -> String -> [String] -> FilePath
                  -> IO (FilePath, MatSeq (StateTree GenSNPState), Int, Int, Int, Emissions [NT], [Site NT])
setupSiteAnalysis regionFile emissionsFile siteLines outputFile = do
  (_, emissionsStart, emissionsEnd) <- readRegionFile regionFile
  numEvents <- countFileLines emissionsFile
  let sites = map parseSite siteLines
  let genMS = genMatSeq
  (Right emissions) <- readMinIONEmissions emissionsFile
  return (outputFile, genMS, numEvents, emissionsStart, emissionsEnd, emissions, sites)

setupSiteAnalysisArgs :: [String]
                      -> IO (FilePath, MatSeq (StateTree GenSNPState), Int, Int, Int, (Emissions [NT]), [Site NT])
setupSiteAnalysisArgs args = do
  (regionFile, emissionsFile, siteLines, outputFile) <- case args of
    [regionFile, emissionsFile, sitesFile, outputFile] -> do
      siteLines <- lines <$> readFile sitesFile
      hPutStrLn stderr $ "Found " ++ show (length siteLines) ++ " sites"
      return (regionFile, emissionsFile, siteLines, outputFile)
    otherwise -> do
      hPutStrLn stderr "Arguments invalid, using test arguments"
      --return ("test_data/region.csv", "test_data/minion_post.csv", testLines, "test_data/snp-calling-test-output.csv")
      return ( "test_data/region.csv"
             , "test_data/minion_post.csv"
             , take 3 testLines
             , "test_data/snp-calling-test-output.csv")
  setupSiteAnalysis regionFile emissionsFile siteLines outputFile

readReference :: Bool -> FilePath -> IO [NT]
readReference rc = (((if rc then reverse . map complement else id) . filter (\x -> any (x==) keyOrder) . map toUpper) <$>) . readFile

complement :: NT -> NT
complement 'G' = 'C'
complement 'A' = 'T'
complement 'T' = 'A'
complement 'C' = 'G'

writePosts :: (FilePath, MatSeq (StateTree GenSNPState), Int, Int, Int, (Emissions [NT]), [Site NT])
           -> IO ()
writePosts (outputFile, genMS, numEvents, emissionsStart, emissionsEnd, emissions, sites) = do
  forM_ sites $ \site -> do
    print $ siteEmissionsBounds softFlank numEvents emissionsStart emissionsEnd site
    let locusStr = show (pos site)
        baseFile = "post_deamer_" ++ locusStr
        csvFile = baseFile ++ ".csv"
        pngFile = baseFile ++ ".png"
        f = id -- (/ log 1000000) . log . (+ 1) . (1000000 *)
        emissions' = siteEmissions softFlank numEvents emissionsStart emissionsEnd emissions site
    writeSNPPost genMS numEvents emissionsStart emissionsEnd emissions (head sites) csvFile
    plotLinesFromFile csvFile f pngFile ("Lg Probability of SNP at locus " ++ locusStr) "locus" "probability" (0, 1)
    hPutStrLn stderr $ "deamer " ++ show (pos site) ++ ".[csv,png]"
  hPutStrLn stderr $ "done writing post"

seqModel :: Int -> [NT] -> ProbSeq [NT]
seqModel flankSize ref = series [ flank
                                , minion . series . map (symbol . return) $ trimmedRef
                                , flank
                                ]
  where refLength = length ref
        trimmedRef = drop flankSize . take (refLength - flankSize) $ ref
        flank =  geometricRepeat (recip $ fromIntegral flankSize * avgEventsPerNT) (symbol "_")

noiseModel :: Int -> ProbSeq [NT]
noiseModel size = geometricRepeat (recip . (avgEventsPerNT *) . fromIntegral $ size) (symbol "_")

main = main''' >> return ()

-- distribution -> p * uniform + (1-p) * distribution, where p = factor * (1-p)
regularize :: (Functor t, Foldable t, Fractional a) => a -> t a -> t a
regularize factor v = fmap (\val -> (1 - pUniform) * val + pUniform * uniform) v
  where pUniform = factor / (1 + factor)
        uniform = 1 / fromIntegral (length v)

--working config: token (1e-6), x 100, flank 50
--working config: token (1e-10), x 10, flank 50
--working config: token (1e-5), x 1, flank 50

main' :: IO ()
main' = do
  let dir = "test_data3/"
      regionFile = dir ++ "region.csv"
      emissionsFile = dir ++ "minion_post.csv"
      referenceFile = dir ++ "reference.txt"
      siteLines = testLines
      outputFile = dir ++ "output.csv"

  env@(outputFile, genMS, numEvents, emissionsStart, emissionsEnd, ems, sites) <-
    setupSiteAnalysis regionFile emissionsFile siteLines outputFile
  let x = 1
      emissions' = ems { emissions = V.map (normalize . regularize x) (emissions ems) }
      emissions'' = addToken (1e-6) "_" emissions'

  ref <- readReference False referenceFile

  let flankSize = 50
      models = [
          (seqModel flankSize ref)
        , (seqModel flankSize (reverse ref))
        , (seqModel flankSize (map complement ref))
        , (seqModel flankSize (reverse (map complement ref)))
        ]
      matSeq = buildMatSeq $ uniformDistOver models
      nLabels = V.length . stateLabels $ matSeq
      nModel = nLabels `div` length models
      ixs = [[[i*nModel .. (i+1)*nModel-1] | i <- [0..length models-1]]]

  [densities] <- runHMMSum' ixs emissions'' matSeq
  let probs = map (/ sum densities) densities
  print densities
  print probs

alleleRefModelIxs :: Int -> Int -> Int -> (Vector Int, Vector Int)
alleleRefModelIxs edgeFlankSize loc refLength = (alleleLabels 'B', alleleLabels 'D')
  where ixLabels = V.imap (,) . V.map stateLabel . stateLabels . buildMatSeq $ alleleRefModel edgeFlankSize loc 'D' (replicate loc 'A' ++ ['B'] ++ replicate (refLength - (loc + 1)) 'C')
        alleleLabels a = V.map fst $ V.filter (\(ix, nts) -> any (\nt -> nt == a) nts) ixLabels

alleleRefModel :: Int -> Int -> NT -> [NT] -> ProbSeq [NT]
alleleRefModel edgeFlankSize loc alt ref = series [
    (noiseModel edgeFlankSize)
  , minion . series $ [
        series . map (symbol . return) $ first
      , eitherOr 0.5 (symbol [ref !! loc]) (symbol [alt])
      , series . map (symbol . return) $ last
      ]
  , (noiseModel edgeFlankSize)
  ]
  where first = drop edgeFlankSize . take loc $ ref
        last = drop (loc + 1) . take (length ref - edgeFlankSize) $ ref

main'' = do
  let dir = "test_data_fwd1/"
      regionFile = dir ++ "region.csv"
      emissionsFile = dir ++ "minion_post.csv"
      referenceFile = dir ++ "reference.txt"
      siteLines = testLines
      outputFile = dir ++ "output.csv"

  env@(outputFile, genMS, numEvents, emissionsStart, emissionsEnd, emissions, sites) <-
    setupSiteAnalysis regionFile emissionsFile siteLines outputFile

  let emissions' = addToken (1e-6) "_" emissions

  ref <- readReference False referenceFile

  let locs = [100, 200 .. 1800] :: [Int]
  probs <- mapM (probRef emissions' ref) locs
  print probs

main''' = do
  args <- getArgs

  let (flipped, regionFile, emissionsFile, referenceFile, outputFile) = case args of
        [flipped, regF, emiF, refF, outF] -> (flipped == "--flip", regF, emiF, refF, outF)
        _ -> let dir = "/mnt/ubuntu/home/ryan1/documents/smol/snp-calling/10-16/"
                 regionFile = dir ++ "scrappie_events.region"
                 emissionsFile = dir ++ "scrappie_events.post"
                 referenceFile = dir ++ "scrappie_events.ref"
                 outputFile = dir ++ "scrappie_events.probs"
             in (True, regionFile, emissionsFile, referenceFile, outputFile)
 --flip albacore.region nanonet.post albacore.ref nanonet.probs
  (_, emissionsStart, emissionsEnd) <- readRegionFile regionFile
  (Right ems) <- readMinIONEmissions emissionsFile

  let emissions' = addToken (1e-6) "_" (expEmissions ems)

  ref <- readReference False referenceFile

  print (length ref)
  print (V.length . emissions $ emissions')
  print (emissionsStart, emissionsEnd)

  let edgeFlankSize = 50
      sites = testSNPs edgeFlankSize 100 emissionsStart ref
      model = fullSNPModel flipped (emissionsStart, emissionsEnd) edgeFlankSize ref sites
      matSeq = buildMatSeq model

      --ixs = fullSNPModelIxs flipped (emissionsStart, emissionsEnd) edgeFlankSize ref sites
      ixs = fullSNPModelIxs' flipped (V.map stateTag . stateLabels $ matSeq)

      x = 1
      emissions'' = emissions' { emissions = V.map (normalize . regularize x) (emissions emissions') }
      --emissions'' = V.map (normalize . regularize x) emissions'

  print "bf"
  densities <- runHMMSum' ixs emissions'' matSeq
  print (length densities)
  print "af"
  let probs = map head $ map (\v -> map (/ sum v) v) densities
      outputStr = unlines . map show $ probs

  putStrLn outputStr
  writeFile outputFile outputStr

expEmissions :: Emissions a -> Emissions a
expEmissions ems = ems { emissions = expVecMat (emissions ems) }
  where expVecMat :: VecMat -> VecMat
        expVecMat vm = V.map (V.map exp) vm


fullSNPModel :: Bool -> (Int, Int) -> Int -> [NT] -> [(Int, NT, Prob)] -> ProbSeq [NT]
fullSNPModel False (startIx, endIx) flankSize ref sites =
  fullSNPModelOriented flankSize ref [(locus - startIx, nt, p) | (locus, nt, p) <- sites]
fullSNPModel True (startIx, endIx) flankSize ref sites = fullSNPModelOriented
  flankSize
  (reverse (map complement ref))
  [(length ref - 1 - (locus - startIx), complement nt, p) | (locus, nt, p) <- sites]

fullSNPModelOriented :: Int -> [NT] -> [(Int, NT, Prob)] -> ProbSeq [NT]
fullSNPModelOriented flankSize ref sites = series [
    noiseModel flankSize
  , minion . series $ regions
  , noiseModel flankSize
  ]
  where lastFirstSites = sortBy (compare `on` (\(a,_,_) -> -a)) sites
        ref' = drop flankSize . take (length ref - flankSize) $ ref
        (before, regions') = foldl' nextRegion (map (symbol . return) ref', []) lastFirstSites
        nextRegion (reference, sofar) (locus, nt, p) =
          let (before, ref:after) = splitAt (locus - flankSize) reference
          in (before, eitherOr (1-p) ref (symbol (return nt)) : series after : sofar)
        regions = series before : regions'

fullSNPModelIxs' :: Bool -> Vector StateTag -> [[Vector Int]]
fullSNPModelIxs' flipped =
    (if flipped then reverse else id)
  . toIxs
  . catMaybes
  . foldl' (\sofar (ix, tag) -> takeLabel ix tag : sofar) []
  . V.toList
  . V.imap (,)
  where takeLabel :: Int -> StateTag -> Maybe (Int, Int, Int)
        takeLabel stateIx (StateTag 1 [ -- series[1]
                              StateTag 0 [StateTag 0 -- minion
                                          collapsedTags]]) =
          getFirstL . map (takeRegionLabel stateIx) $ collapsedTags
        takeLabel _ _ = Nothing
        takeRegionLabel :: Int -> StateTag -> Maybe (Int, Int, Int)
        takeRegionLabel stateIx (StateTag regionTag rest) =
          if even regionTag
          then Nothing
          else case rest of
            [StateTag n _] -> Just ((regionTag - 1) `div` 2, n, stateIx)
            _ -> Nothing
        toIxs :: [(Int, Int, Int)] -> [[Vector Int]]
        toIxs = map (map (V.fromList . map thd3))
          . map (groupBy ((==) `on` snd3))
          . groupBy ((==) `on` fst3)
          . sort

fullSNPModelIxs :: Bool -> (Int, Int) -> Int -> [NT] -> [(Int, NT, Prob)] -> [[Vector Int]]
fullSNPModelIxs False (startIx, endIx) flankSize ref sites =
  fullSNPModelOrientedIxs flankSize ref [(locus - startIx, nt, p) | (locus, nt, p) <- sites]
fullSNPModelIxs True (startIx, endIx) flankSize ref sites = fullSNPModelOrientedIxs
  flankSize
  (reverse (map complement ref))
  [(length ref - 1 - (locus - startIx), nt, p) | (locus, nt, p) <- sites]

fullSNPModelOrientedIxs :: Int -> [NT] -> [(Int, NT, Prob)] -> [[Vector Int]]
fullSNPModelOrientedIxs flankSize ref sites =
  [ [presentIxs (charPair nt), presentIxs nt]
  | (locus, nt, _) <- sites']
  where sites' = zipWith (\c (loc, _, p) -> (loc, c, p)) (take (length sites) [minBound..]) sites
        charPair c = toEnum $ (fromEnum (maxBound :: Char) - (fromEnum c))
        lastFirstSites = sortBy (compare `on` (\(a,_,_) -> -a)) sites'
        ref' = drop flankSize . take (length ref - flankSize) $ ref
        (before, regions') = foldl' nextRegion (map (symbol . return) ref', []) lastFirstSites
        nextRegion (reference, sofar) (locus, nt, p) =
          let (before, ref:after) = splitAt (locus - flankSize) reference
          in (before, eitherOr (1-p) (symbol (return (charPair nt))) (symbol (return nt)) : series after : sofar)
        regions = series before : regions'
        matSeq = buildMatSeq $ series [
            noiseModel flankSize
          , minion . series $ regions
          , noiseModel flankSize
          ]
        labels = V.map stateLabel . stateLabels $ matSeq
        presentIxs nt = V.map fst . V.filter (\(ix, label) -> any (==nt) label) . V.imap (,) $ labels

        
testSNPs :: Int -> Int -> Int -> [NT] -> [(Int, NT, Prob)]
testSNPs flank buffer emissionsStart ref = [ (emissionsStart + locus, complement (ref !! locus), 0.5)
                                           | locus <- loci]
  where loci = [flank + buffer, flank + 2 * buffer .. length ref - flank - buffer]

probRef :: Emissions [NT] -> [NT] -> Int -> IO Double
probRef emissions' ref loc = do
  let edgeFlankSize = 50
      alt = complement (ref !! loc)
      model = alleleRefModel edgeFlankSize loc alt ref
      matSeq = buildMatSeq model
      nLabels = V.length . stateLabels $ matSeq
      (refIxs, altIxs) = alleleRefModelIxs edgeFlankSize loc (length ref)
      ixs = [[refIxs, altIxs]]

      nEmissions = V.length . V.head . emissions $ emissions'
      --emissions'' = V.map (V.map (const (1 / fromIntegral nEmissions))) emissions'
      x = 1
      emissions'' = emissions' { emissions = V.map (normalize . regularize x) (emissions emissions') }

  [densities] <- runHMMSum' ixs emissions'' matSeq
  let probs = map (/ sum densities) densities
  return $ head probs


--main'' = do
  --post <- runHMM minionIndexMap' emissions' matSeq
  --let sums = V.map (\sums -> map (map (sum . V.map (sums V.!))) ixs) post
  --mapM print sums


  --let sites' = reverseSites sites

  --mapM_ print $ snpRegions emissionsStart emissionsEnd sites'

  --res <- callSNPs' genMS numEvents emissionsStart emissionsEnd emissions sites'

  --mapM_ putStrLn
    -- . map (\(Site {..}, pAlt) -> show pos ++ "\t " ++ show (snd . (!! 1) $ alleles) ++ ": \t" ++ show pAlt )
    -- $ zip sites' res

  --writePosts env

  {-

  probAlts <- (concat <$>) . forM sites $ \site -> do
    hPutStrLn stderr $ "running site " ++ show (pos site)
    probAlt <- callSNP genMS numEvents emissionsStart emissionsEnd emissions site
    hPutStrLn stderr $ "found p(alt) " ++ show probAlt
    return probAlt

  h <- openFile outputFile WriteMode
  forM (zip probAlts sites) $ \(p, site) ->
    hPutStrLn h $ show (pos site) ++ "," ++ show (snd $ alleles site !! 1) ++ "," ++ show p
  hClose h
-}
  -- running out of memory!
  -- calls <- forM sites $ \site -> do
    --hPutStrLn stderr "attempting site"
    --p <- callSNP emissionsFile site
    --let res = show (pos site) ++ "," ++ show (maf site) ++ "," ++ show p
    --hPutStrLn h res

  --writeFile outputFile $ unlines calls

readRegionFile :: FilePath -> IO (String, Int, Int)
readRegionFile regionFile = do
  (Right [triple]) <- decode NoHeader <$> BS.readFile regionFile
  return triple

probOfAlt :: [Prob] -> Prob
probOfAlt [refP, altP] = (altP / (refP + altP))

{- CHECK THE INDICES! COULD BE OFF BY 1! -}

softFlank = 50

siteEmissions :: Int -> Int -> Int -> Int -> Emissions [NT] -> Site a -> Emissions [NT]
siteEmissions softFlank numEvents emissionsStart emissionsEnd ems site =
  ems { emissions = V.slice start length (emissions ems)}
  where (start, length) = siteEmissionsBounds softFlank numEvents emissionsStart emissionsEnd site

siteEmissionsBounds :: Int -> Int -> Int -> Int -> Site a -> (Int, Int)
siteEmissionsBounds softFlank numEvents emissionsStart emissionsEnd site =
  (center - softFlank, 2 * softFlank + 1)
  where proportion = (fromIntegral $ pos site - emissionsStart) / (fromIntegral $ emissionsEnd - emissionsStart)
        center = round (fromIntegral numEvents * proportion) - softFlank


writeSNPPost :: MatSeq (StateTree GenSNPState) -> Int -> Int -> Int -> (Emissions [NT]) -> Site NT -> FilePath -> IO ()
writeSNPPost genMS numEvents emissionsStart emissionsEnd ems site outFile = do
  -- let (matSeq, [[refIxs, altIxs]]) = specifyGenMatSeqNT genMS site
  let (matSeq, [[refIxs, altIxs]]) = snpsNTMatSeq emissionsStart emissionsEnd [site]  -- snpsNTMatSeq sites

  post <- runHMM ems matSeq

  let refVec = V.map (\row -> sum $ V.map (row V.!) refIxs) post
      altVec = V.map (\row -> sum $ V.map (row V.!) altIxs) post
      refLine = ("ref," ++ ) . intercalate "," . map show . V.toList $ refVec
      altLine = ("alt," ++ ) . intercalate "," . map show . V.toList $ altVec
      content = unlines [refLine, altLine]

  putStrLn $ "ref loci " ++ progressString emissionsStart (pos site) emissionsEnd
  putStrLn $ "ref ixs " ++ progressString 0 (V.maxIndex refVec) (V.length refVec)
  putStrLn $ "alt ixs " ++ progressString 0 (V.maxIndex altVec) (V.length refVec)
  putStrLn $ "p(alt) " ++ show (sum altVec / (sum altVec + sum refVec))

  writeFile outFile content

progressString :: Int -> Int -> Int -> String
progressString start pnt end =
    show start ++
    "/(" ++ show pnt ++
    ", " ++ show (100 * fromIntegral (pnt - start) / fromIntegral (end - start)) ++
    "%)/" ++ show end



testCallSNP :: MatSeq (StateTree GenSNPState) -> Int -> Int -> Int -> (Emissions [NT]) -> Site NT -> IO [Prob]
testCallSNP genMS numEvents emissionsStart emissionsEnd emissions site = do
  let (matSeq, ixs) = specifyGenMatSeqNT genMS site  -- snpsNTMatSeq sites

  hPutStr stderr $ "evaluating genMS: "
  hPutStrLn stderr $ show (trans genMS # (100, 100))
  hPutStr stderr $ "evaluating siteMS: "
  hPutStrLn stderr $ show (trans matSeq # (100, 100))

  hPutStrLn stderr $ "subsampling emissions file with region size " ++ show (2 * softFlank + 1)
  --subsampledEmissionsFile <- emptySystemTempFile emissionsFile
  --putStrLn $ "using temporary subsampled emissions file: " ++ subsampledEmissionsFile

  let subEmissions = siteEmissions softFlank numEvents emissionsStart emissionsEnd emissions site

  ps <- iterateUntilDiffM $ do
    ps <- runHMMSum' ixs subEmissions matSeq
    hPutStrLn stderr $ "results: " ++ intercalate "," (map show ps)
    return ps

  return $ map probOfAlt ps

iterateUntilDiffM :: (Monad m, Eq a) => m a -> m a
iterateUntilDiffM action = do
  comp <- action
  iterateWhile (/= comp) action

callSNPs' :: MatSeq (StateTree GenSNPState) -> Int -> Int -> Int -> (Emissions [NT]) -> [Site NT] -> IO [Prob]
callSNPs' genMS numEvents emissionsStart emissionsEnd emissions sites = do
  --let (matSeq, ixs) = specifyGenMatSeqNT genMS site  -- snpsNTMatSeq sites
  let (matSeq, ixs) = snpsNTMatSeq'' 10 "_" emissionsStart emissionsEnd sites  -- snpsNTMatSeq sites

  hPutStr stderr $ "evaluating siteMS: "
  hPutStrLn stderr $ show (trans matSeq # (100, 100))

  --hPutStrLn stderr $ "subsampling emissions file with region size " ++ show (2 * softFlank + 1)
  --subsampledEmissionsFile <- emptySystemTempFile emissionsFile
  --putStrLn $ "using temporary subsampled emissions file: " ++ subsampledEmissionsFile

  --let subEmissions = siteEmissions softFlank numEvents emissionsStart emissionsEnd emissions site

  let emissions' = addToken 0.1 "_" emissions

  ps <- runHMMSum' ixs emissions' matSeq

  return $ map probOfAlt ps

addToken :: (Ord a) => Prob -> a -> Emissions a -> Emissions a
addToken p token (Emissions {..})= Emissions emissions' (Map.insert token tokenIx indexMap)
  where (tokenIx, emissions') = addEmissionsToken p emissions

addTokenToIndexMap :: (Ord a) => a -> Map a Int -> Map a Int
addTokenToIndexMap tokenKey indexMap = Map.insert tokenKey (Map.size indexMap) indexMap

addEmissionsToken :: Prob -> VecMat -> (Int, VecMat)
addEmissionsToken p emissions = ( V.length (V.head emissions)
                                , V.map ((`V.snoc` p) . (((1-p)*) <$>)) emissions)



callSNPs :: MatSeq (StateTree GenSNPState) -> Int -> Int -> Int -> Emissions [NT] -> [Site NT] -> IO [Prob]
callSNPs genMS numEvents emissionsStart emissionsEnd emissions sites = do
  --let (matSeq, ixs) = specifyGenMatSeqNT genMS site  -- snpsNTMatSeq sites
  let (matSeq, ixs) = snpsNTMatSeq emissionsStart emissionsEnd sites  -- snpsNTMatSeq sites

  hPutStr stderr $ "evaluating siteMS: "
  hPutStrLn stderr $ show (trans matSeq # (100, 100))

  hPutStrLn stderr $ "subsampling emissions file with region size " ++ show (2 * softFlank + 1)
  --subsampledEmissionsFile <- emptySystemTempFile emissionsFile
  --putStrLn $ "using temporary subsampled emissions file: " ++ subsampledEmissionsFile

  --let subEmissions = siteEmissions softFlank numEvents emissionsStart emissionsEnd emissions site

  ps <- runHMMSum' ixs emissions matSeq

  hPutStrLn stderr $ "results: " ++ intercalate "," (map show ps)

  return $ map probOfAlt ps


callSNP :: MatSeq (StateTree GenSNPState) -> Int -> Int -> Int -> (Emissions [NT])-> Site NT -> IO [Prob]
callSNP genMS numEvents emissionsStart emissionsEnd emissions site = do
  --let (matSeq, ixs) = specifyGenMatSeqNT genMS site  -- snpsNTMatSeq sites
  let (matSeq, ixs) = snpsNTMatSeq emissionsStart emissionsEnd [site]  -- snpsNTMatSeq sites

  hPutStr stderr $ "evaluating genMS: "
  hPutStrLn stderr $ show (trans genMS # (100, 100))
  hPutStr stderr $ "evaluating siteMS: "
  hPutStrLn stderr $ show (trans matSeq # (100, 100))

  hPutStrLn stderr $ "subsampling emissions file with region size " ++ show (2 * softFlank + 1)
  --subsampledEmissionsFile <- emptySystemTempFile emissionsFile
  --putStrLn $ "using temporary subsampled emissions file: " ++ subsampledEmissionsFile

  let subEmissions = siteEmissions softFlank numEvents emissionsStart emissionsEnd emissions site

  ps <- runHMMSum' ixs emissions matSeq

  hPutStrLn stderr $ "results: " ++ intercalate "," (map show ps)

  return $ map probOfAlt ps

instance Hashable a => Hashable (Vector a) where
  hashWithSalt = hashUsing V.toList

runHMMSum' :: [[Vector Int]]
           -> Emissions [NT]
           -> MatSeq String
           -> IO [[Prob]]
runHMMSum' ixs emissions priorSeq = do
  let (Posterior sumsMap) = posteriorSHMM
        emissions
        priorSeq
      sums = V.generate (IntMap.size sumsMap) (fromJust . flip IntMap.lookup sumsMap)
  return $ map (map (sum . V.map (sums V.!))) ixs

-- posteriorSHMM :: (Ord s, Show s) => Emissions s -> MatSeq s -> Posterior

runHMM :: Emissions String
       -> MatSeq String
       -> IO VecMat
runHMM emissions priorSeq = undefined
  --SHMM.shmmFull (nStates (trans priorSeq)) triples emissions permutation
  --where permutation = buildEmissionPerm indexMap priorSeq
        --triples = matSeqTriples priorSeq

{-
callSNPs :: FilePath -> [Site NT] -> IO [Prob]
callSNPs emissionsFile sites = do
  let (matSeq, ixs) = specifyGenMatSeqNT (head sites)  -- snpsNTMatSeq sites

  putStrLn $ "done building " ++ show (V.length (stateLabels matSeq))
  post <- SHMM.runHMM minionIndexMap emissionsFile matSeq
  let ps = map (map (\arr -> stateProbs arr post)) ixs
  putStrLn $ "results: " ++ intercalate "," (map show ps)

  return $ map probOfAlt ps
-}

parseSite :: String -> Site NT
parseSite row = Site {
    alleles = [(ref, 1-maf), (alt, maf)]
  , pos = pos
  , leftFlank = leftFlank
  , rightFlank = rightFlank
  }
  where
    [read -> pos, [ref], [alt], (read::String->Double) -> maf, leftFlank, [ref'], rightFlank] = words row

{-
callSNP :: FilePath
        -> Site NT
        -> IO Prob
callSNP emissionsFile site = do
  let (ms, refIxs, altIxs) = snpMatSeq site
  --write this as a configuration file
  --print =<< sampleSeq vecDist ms
  --writeEmissionPerm minionIndexMap "emission_perm.csv" ms
  --writeSTFile ms "test_ex.st"
  post <- SHMM.runHMM minionIndexMap emissionsFile ms
  let refP = stateProbs refIxs post
      altP = stateProbs altIxs post
  return (altP / (refP + altP))
-}

stateProbs :: Vector Int -> Emissions String -> Prob
stateProbs ixs = sum . V.map (\row -> sum . V.map (row V.!) $ ixs) . emissions

--main = compareSmall

  {-
compareMicrosatellites :: IO ()
compareMicrosatellites = do
  let satellite :: ProbSeq String
      satellite = series . map (\c -> state [c]) $ "ATTTA"
      genSeq = buildMatSeq . minion $ repeatSequence 15 satellite
      priorSeq = buildMatSeq . minion $ andThen
                    (repeatSequence 10 satellite)
                    (uniformDistRepeat 20 satellite)
      getRepeat :: ((Prob, (s, StateTag)) -> Maybe Int)
      getRepeat (_, (_, StateTag 1 [StateTag n _])) = Just n
      getRepeat _ = Nothing
  print =<< sampleSeq vecDist priorSeq
  print $ V.head $ stateLabels priorSeq

  compareHMM genSeq priorSeq getRepeat

compareSmall :: IO ()
compareSmall = do
  let genSeq = buildMatSeq $ repeatSequence 15 periodSeq
      priorSeq = buildMatSeq $ andThen
                    (repeatSequence 0 periodSeq)
        (uniformDistRepeat 20 periodSeq)
      f = getRepeat
  compareHMM genSeq priorSeq f

getRepeat :: ((Prob, (s, StateTag)) -> Maybe Int)
getRepeat (_, (_, StateTag 1 [StateTag n _])) = Just n
getRepeat _ = Nothing

compareHMM :: (Show s, Eq s)
           => MatSeq s
           -> MatSeq s
           -> ((Prob, (s, StateTag)) -> Maybe Int)
           -> IO ()
compareHMM genSeq priorSeq f = do
  (V.unzip -> (sample, sampleIxs), posterior, viterbi, forward, backward) <- runHMM genSeq priorSeq obsProb

  putStrLn $ "state labels:" ++ show sample
  putStrLn $ "state generating index:" ++ show sampleIxs

  putStrLn $ "truth: " ++ show (fromMaybe 0 $ f (undefined, stateLabels priorSeq V.! V.last sampleIxs))
  putStrLn $ "viterbi: " ++ show (fromMaybe 0 $ f (undefined, stateLabels priorSeq V.! V.last viterbi))
  let tags = pathProbs (stateLabels priorSeq) posterior
      post_dist = (distOver f $ V.last tags)
      max_post = fst $ maximumBy (compare `on` snd) post_dist
  putStrLn $ "max posterior: " ++ show max_post
  putStrLn "posterior dist: "
  mapM_ (\(ix, p) -> putStrLn $ show ix ++ ": " ++ show (fromRational p)) (distOver f $ V.last tags)
  let tags = pathProbs (stateLabels priorSeq) forward
      post_dist = (distOver f $ V.last tags)
      max_post = fst $ maximumBy (compare `on` snd) post_dist
  --putStrLn $ "max forward: " ++ show max_post
  --putStrLn "forward dist: "
  --mapM_ (\(ix, p) -> putStrLn $ show ix ++ ": " ++ show (fromRational p)) (distOver f $ V.last tags)
  let tags = pathProbs (stateLabels priorSeq) backward
      post_dist = (distOver f $ V.last tags)
      max_post = fst $ maximumBy (compare `on` snd) post_dist
  --putStrLn $ "max backward: " ++ show max_post
  --putStrLn "backward dist: "
  --mapM_ (\(ix, p) -> putStrLn $ show ix ++ ": " ++ show (fromRational p)) (distOver f $ V.last tags)
  --print priorSeq

  return ()
-}

obsProb :: (Eq s) => MatSeq s -> s -> IO (Vector Prob)
obsProb seq i = return . normalize . V.map (\i' -> if i == i' then 1.0 else 0.0) . V.map stateLabel . stateLabels $ seq

sampleCycleMatIxs :: IO (Vector Int)
sampleCycleMatIxs = fst <$> sampleSeqIxs vecDist cycleMatSeq

sampleCycleMatIxs' :: IO (Vector Int)
sampleCycleMatIxs' = fst <$> sampleSeqIxs vecDist cycleMatSeq'

sampleCycleMat :: IO (Vector Int)
sampleCycleMat = fst <$> sampleSeq vecDist cycleMatSeq

cycleMatSeq :: MatSeq Int
cycleMatSeq = buildMatSeq cycleSeq

cycleMatSeq' :: MatSeq Int
cycleMatSeq' = buildMatSeq cycleSeq'

cycleSeq' :: ProbSeq Int
cycleSeq' = repeatSequence n periodSeq

n = 15

cycleSeq :: ProbSeq Int
cycleSeq = andThen
  (repeatSequence nStart periodSeq)
  (uniformDistRepeat nEnd periodSeq)
  where nStart = 0
        nEnd = 20

periodSeq :: ProbSeq Int
periodSeq = series' . map (\v -> andThen (symbol v) skipDSeq) $
  [ 2, 1 ]

skipD :: [Prob]
skipD = [0.5, 1 - head skipD]

skipDSeq :: ProbSeq a
skipDSeq = finiteDistRepeat skipD $ skip 1

  -- runs out of memory between 8 and 16
  -- "62525746 G C 0.002396 CAGGAGCACC G GCCGCAGAGG"
testLines :: [String]
testLines = [
    "62525667 T C 0.3954 CCTTGGATGC T ACTGGGTTTG" -- potential
  , "62525746 G C 0.002396 CAGGAGCACC G GCCGCAGAGG"
  , "62525767 T A 0.0009984 TCTGGGAGCT T CTAGGATGGG"
  , "62525777 G C 0.09844 TCTAGGATGG G AAGTGGCCCA"
  , "62525782 G A 0.001198 GATGGGAAGT G GCCCAGGCAG"
  , "62525813 A G 0.0003994 GCAGGCCGTC A GTGAGTGGCG"
  , "62525990 G C 0.007788 TGGACAGGAA G GAAGGAAGGA"
  , "62526024 T G 0.4101 CGCTCCAGGG T TGAGCAAATG" -- potential
  , "62526065 G A 0.03035 CCGGGTGGGG G CGGGGGCGAC"
  , "62526134 T G 0.003994 ACGATTACGT T TTCTCAGTCT"
  , "62526155 G T 0.003994 TACTTAAAGC G CTGAGTAAAC"
  , "62526212 C A 0.005192 CTGCGCGGTT C CCCGCAGCAC"
  , "62526224 T C 0.03195 CCGCAGCACA T GGCGTGTCCA"
  , "62526334 C T 0.01418 CAGGGCGCCT C GGCCCCGGGC"
  , "62526352 C G 0.0003994 GGCTGTCACT C GGGACTCCGC"
  , "62526369 A G 0.000599 CCGCCCCTTC A TGGACGGAGC"
  --, "62526373 A G 0.01078 CCCTTCATGG A CGGAGCCTCC" -- can't deal with this. flank overlaps with another snp! wat!
  , "62526500 T T 0.00599 CGTGCTCGTC T CCGCTGCCGC"
  , "62526547 T T 0.08926 CCGCGCCCTC T GCCGCCGCCG"
  -- , "62526696 G A 0.04034 CTGCGGGTCG G GCGGGCGGAT" -- this one gives NaN
  , "62526719 G A 0.08187 GCCCACGTCA G GCCCGGGCAG"
  , "62526801 G C 0.09784 CCCCCGGGCC G GGGCTGCGCG"
  , "62526815 G T 0.001398 CTGCGCGGGC G CTCGGGGCCG"
  , "62526818 C T 0.0007987 CGCGGGCGCT C GGGGCCGGAG"
  , "62526898 G A 0.002796 TGCGGGAGCC G GGCCGGGCCG"
  , "62527012 C T 0.01018 GCGGCCGCCC C CAACCCCCCG"
  , "62527231 C T 0.02157 TTTATAAAAA C ATTTGAAGCC"
  , "62527305 G A 0.09006 CTAACTTGTT G GTGTTAAGTG"
  , "62527322 T A 0.002596 AGTGTCTGGA T TAAAGACTCT"
  , "62527414 A C 0.000599 TTTCAGCTTT A TTTTTGTTTT"
  , "62527427 G A 0.001997 TTTGTTTTCC G GCTTAGGCTT"
  , "62527467 C T 0.002396 TGGTTAGACA C GTCTGCCCTT"
  ]

reverseSite :: Site a -> Site a
reverseSite site = site {
    leftFlank = reverse (leftFlank site)
  , rightFlank = reverse (rightFlank site)
  }

reverseSites :: [Site a] -> [Site a]
reverseSites = reverse . map reverseSite

getFirstL :: [Maybe a] -> Maybe a
getFirstL = getFirst . mconcat . map First

fst3 :: (a, b, c) -> a
fst3 (a, _, _) = a

snd3 :: (a, b, c) -> b
snd3 (_, b, _) = b

thd3 :: (a, b, c) -> c
thd3 (_, _, c) = c

readMinIONEmissions :: FilePath -> IO (Either String (Emissions String))
readMinIONEmissions = (((\ems -> Emissions ems minionIndexMap) <$>) <$>) . readMat

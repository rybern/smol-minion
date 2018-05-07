{-# LANGUAGE OverloadedLists, RecordWildCards #-}
module Tester where

import Sequence
import EmissionIO
import Data.Scientific
import Data.List
import Data.Maybe
import Data.Vector (Vector)
import qualified Data.Vector as V
import Control.Monad.Random
import qualified SHMM as SHMM
import Sequence.Matrix.ProbSeqMatrixUtils
import Math.LinearAlgebra.Sparse.Matrix hiding (trans)
import Data.Set (Set)
import qualified Data.Set as Set
import Data.Map (Map)
import qualified Data.Map as Map

import Plot
import MinION
import GeneralizedSNP
import SNP
import EmissionPermutation

type StateIx = Int
type RMSE = Double

data SimulationResult s = SimRes {
    states :: Vector StateIx
  , rmse :: RMSE
  , misclassify :: Double
  , post :: Emissions
  , emissions :: Emissions
  , maxLabels :: Vector s
  , trueLabels :: Vector s
  , indexMap :: Map s Int
  } deriving Show

prettyPrintSimulationResult :: (Show s) => Bool -> SimulationResult s -> IO ()
prettyPrintSimulationResult showPost (SimRes {..}) = do
  putStrLn $ "True state labels:"
  putStrLn $ "\t<" ++ intercalate "," (V.toList . V.map show $ trueLabels) ++ ">"
  putStrLn $ "Max state labels:"
  putStrLn $ "\t<" ++ intercalate "," (V.toList . V.map show $ maxLabels) ++ ">"
  putStrLn $ "RMSE: " ++ show rmse
  putStrLn $ "% misclassified: " ++ show (100*misclassify)
  when showPost $ do
    putStrLn $ "posterior:"
    prettyPrintEmissions post

prettyPrintEmissions :: Emissions -> IO ()
prettyPrintEmissions = mapM_ putStrLn . prettyPrintEmissionsLines

prettyPrintEmissions2 :: [Emissions] -> IO ()
prettyPrintEmissions2 = mapM_ putStrLn . map (intercalate " || ") . transpose . map prettyPrintEmissionsLines

prettyPrintEmissionsLines :: Emissions -> [String]
prettyPrintEmissionsLines ems = V.toList . flip V.map ems $ ('\t':) . intercalate "\t" . V.toList . V.map prettyPrintFloat

prettyPrintFloat :: (RealFloat a) => a -> String
prettyPrintFloat = formatScientific Generic (Just 3) . fromFloatDigits

simulationSimple :: (Ord b) => Map b Int -> MatSeq b -> IO (SimulationResult b)
simulationSimple = simulation simpleNoise

simulationUniform :: (Ord b) => Prob -> Map b Int -> MatSeq b -> IO (SimulationResult b)
simulationUniform p = simulation (uniformNoise p)

simulation :: (Ord b) => (Int -> StateIx -> Vector Prob) -> Map b Int -> MatSeq b -> IO (SimulationResult b)
simulation noise indexMap ms = do
  (observations, truth) <- simulate indexMap noise ms
  estimated <- posterior indexMap ms observations
  return $ SimRes {
      states = truth
    , rmse = scoreRMSE truth estimated
    , misclassify = scoreMisclassify truth estimated
    , post = estimated
    , emissions = observations
    , maxLabels = maxLabelPath ms estimated
    , trueLabels = V.map (fst . (stateLabels ms V.!)) truth
    , indexMap = indexMap
    }

simulation2 :: (Ord b, Show b)
            => (Int -> StateIx -> Vector Prob) -> Map b Int -> MatSeq b -> MatSeq b -> IO (SimulationResult b)
simulation2 noise indexMap msSim msModel = do
  (observations, truth) <- simulate indexMap noise msSim
  estimated <- posterior indexMap msModel observations
  return $ SimRes {
      states = truth
    , rmse = scoreRMSE truth estimated
    , misclassify = scoreMisclassify truth estimated
    , post = estimated
    , emissions = observations
    , maxLabels = maxLabelPath msModel estimated
    , trueLabels = V.map (fst . (stateLabels msSim V.!)) truth
    , indexMap = indexMap
    }

posterior :: (Ord a) => Map a Int -> MatSeq a -> Emissions -> IO Emissions
posterior indexMap ms observations = do
  let ns = nStates (trans ms)
      triples = matSeqTriples ms
      permutation = buildEmissionPerm indexMap ms
  print $ "ns: " ++ show ns
  print $ "n triples: " ++ show (length triples)
  print $ "max row triples: " ++ show (maximum (map (\(x, _, _) -> x) triples))
  print $ "max col triples: " ++ show (maximum (map (\(_, x, _) -> x) triples))
  print $ "permutation length: " ++ show (V.length permutation)
  print $ "permutation min: " ++ show (V.minimum permutation)
  print $ "permutation max: " ++ show (V.maximum permutation)
  print $ "observations rows: " ++ show (V.length observations)
  print $ "observations min cols: " ++ show (V.minimum (V.map V.length observations))
  print $ "observations max cols: " ++ show (V.maximum (V.map V.length observations))
  post <- SHMM.shmmFull ns triples observations permutation
  print $ "called shmm got " ++ show (V.length post)
  let post' = V.map (V.tail . V.init) $ post
  return post

defaultIndexMap :: (Ord a) => MatSeq a -> Map a Int
defaultIndexMap = Map.fromList . flip zip [0..] . nub . map fst . V.toList . stateLabels

matSeqTriples :: MatSeq a
              -> [(Int, Int, Double)]
matSeqTriples = map (\((r, c), p) -> (r - 1, c - 1, p)) . tail . toAssocList . cleanTrans . trans

simulate :: (Ord a, MonadRandom m)
         => Map a Int -> (Int -> StateIx -> Vector Prob) -> MatSeq a -> m (Emissions, Vector StateIx)
simulate indexMap observe ms = do
  --let ixs = [1367,1370,1374,1376,1380,1383,1384,1389,1406,1480,2992,2861,1854,1918,3200,3202,3205,3211,3213,3218,3239,3303,4568,4523,4177,4220,3884,3854,3714,4259,4310,4514,4305,4495,4268,4347]
  --let ixs = [7,15,17,21,23,24,25,29,43,326,760,536,663]
  --let ixs = [1367,1370,1372,1374,1378,1381,1389,1708,2680,2203,2067,2072,2155,2421,2461,2619,2227,2087]
  (ixs, _) <- sampleSeqIxs vecDist ms
  let emIxs = V.map (fromJust . flip Map.lookup indexMap . fst . (stateLabels ms V.!)) ixs
      len = Map.size indexMap
      observations = V.map (observe len) emIxs
  return (observations, ixs)

uniformNoise :: Double -> Int -> StateIx -> Vector Prob
uniformNoise coef len ix = V.replicate len base V.// [(ix, coef * base)]
  where base = 1 / (fromIntegral len - 1 + coef)

simpleNoise :: Int -> StateIx -> Vector Prob
simpleNoise = onehot

onehot :: Int -> Int -> Vector Prob
onehot len ix = V.replicate len 0 V.// [(ix, 1)]

meanScoreByRow :: (StateIx -> Vector Prob -> Double) -> Vector StateIx -> Emissions -> Double
meanScoreByRow byRow ixs obs = (\v -> sum v / fromIntegral (V.length v)) . V.map (uncurry byRow) $ V.zip ixs obs

scoreMisclassify :: Vector StateIx -> Emissions -> Double
scoreMisclassify = meanScoreByRow score
  where score ix ob = if V.maxIndex ob == ix then 0.0 else 1.0

scoreRMSE :: Vector StateIx -> Emissions -> RMSE
scoreRMSE truth = sqrt . meanScoreByRow score truth
  where score ix = sum
                 . V.imap (\ix' v -> let t = if ix == ix' then 1.0 else 0.0
                                     in (v - t) ** 2)

maxIxPath :: Emissions -> Vector StateIx
maxIxPath = V.map V.maxIndex

maxLabelPath :: MatSeq a -> Emissions -> Vector a
maxLabelPath ms = V.map (fst . (stateLabels ms V.!)) . maxIxPath

basic :: MatSeq Int
basic = buildMatSeq $ uniformDistOver [ series . map state $ [1..5]
                                      , series . map state $ [6..10]
                                      , series . map state $ [11..15]]

window :: MatSeq String
window = buildMatSeq . minion . series $ [
    series' "hello "
  , eitherOr 0.5 (series' "minion") (series' "world")
  , series' "!!"
  ]
  where series' = series . map (state . return)


testCalling :: IO ()
testCalling = do
  let site = Site {pos = 50, alleles = [('X', 0.5), ('Y', 0.5)] , leftFlank = "ABCDEF" , rightFlank = "fedcba"}
  let (ms, [[refIxs, altIxs]]) = snpsNTMatSeq 0 100 [site]
      regions = snpRegions 0 100 [site]
      indexMap = defaultIndexMap ms

  print regions

  res <- simulation2 (uniformNoise 2) indexMap ms ms

  let post' = post res
      altVec = V.map (\row -> sum $ V.map (row V.!) altIxs) post'
      refVec = V.map (\row -> sum $ V.map (row V.!) refIxs) post'

  putStrLn $ "ref ixs " ++ progressString 0 (V.maxIndex refVec) (V.length refVec)
  putStrLn $ "alt ixs " ++ progressString 0 (V.maxIndex altVec) (V.length refVec)
  putStrLn $ "ref: " ++ show (sum refVec)
  putStrLn $ "alt: " ++ show (sum altVec)
  putStrLn $ "p(alt): " ++ show (sum altVec / (sum refVec + sum altVec))
  putStrLn $ show (trueLabels res)
  putStrLn $ show refVec
  putStrLn $ show altVec
  prettyPrintSimulationResult False res

progressString :: Int -> Int -> Int -> String
progressString start pnt end =
    show start ++
    "/(" ++ show pnt ++
    ", " ++ show (100 * fromIntegral (pnt - start) / fromIntegral (end - start)) ++
    "%)/" ++ show end

ms1 = buildMatSeq $ series [
    geometricRepeat (1 - 1/10) (uniformDistOver [state "A", state "B"])
  , eitherOr 0.5 (series [state "C", state "D"]) (series [state "E", state "D"])
  , geometricRepeat (1 - 1/10) (uniformDistOver [state "A", state "B"])
  ]

ms2 = buildMatSeq $ series [
    geometricRepeat (1 - 1/10) (state "_")
  , eitherOr 0.5 (series [state "C", state "D"]) (series [state "E", state "D"])
  , geometricRepeat (1 - 1/10) (state "_")
  ]

noiseWithToken :: Prob -> (Int -> StateIx -> Vector Prob) -> Int -> StateIx -> Vector Prob
noiseWithToken pToken noise len ix = (((1 - pToken) *) <$> noise (pred len) ix) `V.snoc` pToken

addTokenToIndexMap :: (Ord a) => a -> Map a Int -> Map a Int
addTokenToIndexMap tokenKey indexMap = Map.insert tokenKey (Map.size indexMap) indexMap

testSplit :: IO ()
testSplit = do
  res <- simulation2 (noiseWithToken 0.2 (uniformNoise 2)) (addTokenToIndexMap "_" $ defaultIndexMap ms1) ms1 ms2
  prettyPrintSimulationResult False res

{-
object is the size of the indexMap

*** Exception: ./Data/Vector/Generic/Mutable.hs:845 (update): index out of bounds (1570,1462)
CallStack (from HasCallStack):
  error, called at ./Data/Vector/Internal/Check.hs:87:5 in vector-0.12.0.0-DGP48m6Q6tLL7qSICSsTxc:Data.Vector.Internal.Check

probably something to do with indexMap

-}
n = 2
repeatN :: Int -> [a] -> [a]
repeatN n = concat . replicate n
flankSize = n * 10

site1 = Site {
    pos = n*50
  , alleles = [('W', 0.5), ('X', 0.5)]
  , leftFlank = repeatN n "ABCDEFGHIJ"
  , rightFlank = repeatN n "jihgfedcba"
  }
site2 = Site {
    pos = n*70
  , alleles = [('Y', 0.5), ('Z', 0.5)]
  , leftFlank = repeatN n "ABCDEFGHIJ"
  , rightFlank = repeatN n "jihgfedcba"
  }

reverseSite :: Site a -> Site a
reverseSite site = site {
    leftFlank = reverse (leftFlank site)
  , rightFlank = reverse (rightFlank site)
  }

reverseSites :: [Site a] -> [Site a]
reverseSites = reverse . map reverseSite



sites = [site1]
(cms1, ixs1) = snpsNTMatSeq (n*0) (n*100) sites

token = "_"
(cms2, ixs2) = snpsNTMatSeq'' Tester.flankSize token (n*0) (n*100) sites

cIndexMap = addTokenToIndexMap token $ defaultIndexMap cms1

avgSimulationSize :: MatSeq s -> IO Double
avgSimulationSize ms = (\xs -> fromIntegral (sum xs) / fromIntegral (length xs)) <$> replicateM 10000 (length . fst <$> sampleSeq vecDist ms)

showIxs :: MatSeq [NT] -> [[Vector Int]] -> IO ()
showIxs ms ixs = forM_ ixs $ \[refIxs, altIxs] -> do
  putStrLn "ref:"
  forM_ refIxs $ \ix -> do
    putStrLn $ "\t\t" ++ show (stateLabels ms V.! ix)
  putStrLn ""
  putStrLn "alt:"
  forM_ altIxs $ \ix -> do
    putStrLn $ "\t" ++ show (stateLabels ms V.! ix)
  putStrLn ""

extractSets :: Emissions -> [[Vector Int]] -> [[Emissions]]
extractSets ems = map (map (extractSet ems))
-- correctly more likely once aligned
-- not aligning very well
extractSet :: Emissions -> Vector Int -> Emissions
extractSet ems ixs = V.map (\row -> V.map (row V.!) ixs) ems

testCalling2 :: Int -> Double -> IO Double --(SimulationResult [NT])
testCalling2 flankSize p = do
  let (cms2, ixs2) = snpsNTMatSeq'' flankSize token (n*0) (n*100) sites
  res <- simulation2 (noiseWithToken 0.05 (uniformNoise p)) cIndexMap cms1 cms2

  ps <- forM (zip3 ixs1 ixs2 sites) $ \([ref1, alt1], [ref2, alt2], site) -> do
    let isRef = (V.any (flip V.elem ref1) (states res))
        trueAllele = let [ref, alt] = map return . map fst . alleles $ site
                     in if isRef then ref else alt
    putStrLn $ "truly allele is " ++ (if isRef then "ref" else "alt") ++ " (" ++ trueAllele ++ ")"

    let altVec = extractSet (post res) alt2
        refVec = extractSet (post res) ref2
        altP = sum (V.map sum altVec)
        refP = sum (V.map sum refVec)
        pPref = altP / (altP + refP)

        --isAlt = V.any (any (== 'W')) . trueLabels $ res
        pTrue = if isRef then 0.0 else 1.0
        err = abs (pPref - pTrue)

    putStrLn $ "post p(alt) is " ++ show pPref
    putStrLn $ "error is " ++ show err
    return err
  return (mean ps)

mean :: Floating a => [a] -> a
mean xs = sum xs / fromIntegral (length xs)

xticks :: [Double]
xticks = [2, 5, 10, 15, 20, 25, 30, 40, 50]
repeats = 10
errPoint :: Double -> IO Double
errPoint = (mean <$>) . replicateM repeats . testCalling2 1

errPoints :: IO [(Double, Double)]
errPoints = zip xticks <$> mapM errPoint xticks

flankTicks :: [Int]
flankTicks = [10,20..50]
errPoint' :: Int -> IO Double
errPoint' flankSize = (mean <$>) . replicateM repeats $ testCalling2 flankSize 30

errPoints' :: IO [(Double, Double)]
errPoints' = zip (map fromIntegral flankTicks) <$> mapM errPoint' flankTicks

pts :: [(Double, Double)]
pts = [(2.0,0.4981651132396747),(5.0,0.4896712778189328),(10.0,0.46877311841114044),(15.0,0.4472552843617038),(20.0,0.43254715918499764),(25.0,0.3678726335210324),(30.0,0.26440086960393167),(40.0,0.2206910165519158),(50.0,0.10865170115692621)]

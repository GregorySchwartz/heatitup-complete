{- Merge
Gregory W. Schwartz

Collects the functions pertaining to the merging of headers with abundances.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Merge
    ( mergeAbundance
    , mergeMates
    , bamToFasta
    ) where

-- Standard
import Data.Function (on)
import Data.List
import Data.Maybe
import Data.Monoid
import qualified Data.IntMap.Strict as IMap
import qualified Data.Map.Strict as Map

-- Cabal
import Control.Lens
import qualified Data.Fasta.Text as F
import qualified Data.Text as T
import qualified Data.Text.Read as T
import Safe
import TextShow

-- Local
import Types
import Utility

mergeAbundance :: FrequencyMap -> DuplicationRow -> DuplicationRow
mergeAbundance (FrequencyMap m) (DuplicationRow row) =
    DuplicationRow
        . Map.insert "frequency" (showt . findWithError rowAccession $ m)
        $ row
  where
    rowAccession =
        flip (Safe.at) 1 . T.splitOn "|" . findWithError "fHeader" $ row

-- | Get the position map of a sequence.
getPositionMap :: Position -> Sequence -> PositionMap
getPositionMap (Position p) =
    PositionMap . IMap.fromList . zip [p..] . T.unpack . unSequence

-- | Fill the gaps in a position sequence list with a character.
fillGapsWith :: Fill
             -> [(Position, Char)]
             -> [(Position, Char)]
             -> [(Position, Char)]
fillGapsWith _ !acc [] = reverse acc
fillGapsWith fill [] (x:xs) = fillGapsWith fill [x] xs
fillGapsWith
    fill
    acc@((Position !previousPos, _):_)
    allNext@((Position !currentPos, !currentChar):xs) =
    if abs (currentPos - previousPos) > 1
        then fillGapsWith
                fill
                ((Position (previousPos + 1), unFill fill) : acc)
                allNext
        else fillGapsWith
                fill
                ((Position currentPos, currentChar) : acc)
                xs

-- | Merge two sequences by their positions.
mergeSequences :: Fill
               -> (Position, Sequence)
               -> (Position, Sequence)
               -> Sequence
mergeSequences fill left =
    Sequence
        . T.pack
        . fmap snd
        . fillGapsWith fill []
        . fmap (over _1 Position)
        . IMap.toAscList
        . IMap.union (unPositionMap . uncurry getPositionMap $ left)
        . unPositionMap
        . uncurry getPositionMap

-- | Merge a mate pair. Assumes that the lines are sorted by mate pair already
-- (samtools sort -n). Also assumes that the left mate appears first. Removes
-- supplementary alignments.
mergePair :: Fill -> [BamRow] -> Maybe BamRow
mergePair fill xs =
    case filter valid xs of
        []  -> Nothing
        [x] -> Just x
        [BamRow left, BamRow right] ->
            Just
                . BamRow
                . set (ix 9) ( unSequence
                            $ mergeSequences fill (posSeq left) (posSeq right)
                            )
                $ left
        xs -> error ("Too many mate pairs (even after supplementary removal) for: " <> show xs)
  where
    posSeq xs = ( either error (Position . fst) . T.decimal $ xs !! 3
                , Sequence . revComplCheck (BamRow xs) $ xs !! 9
                )
    revComplCheck row fSeq =
      if checkSamFlag 16 row
        then (F.fastaSeq . F.revCompl $ F.FastaSequence "" fSeq)
        else fSeq
    valid x = not (checkSamFlag 2048 x) && not (checkSamFlag 256 x)  -- Make sure not supplementary nor not primary

-- | Merge all mate pairs. Assumes that the lines are sorted by mate pair already
-- (samtools sort -n).
mergeMates :: Fill -> [BamRow] -> [BamRow]
mergeMates fill =
    catMaybes . fmap (mergePair fill) . groupBy ((==) `on` (head . unBamRow))

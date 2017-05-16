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
import Data.Monoid
import qualified Data.IntMap.Strict as IMap
import qualified Data.Map.Strict as Map

-- Cabal
import Control.Lens
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
-- (samtools sort -n). Also assumes that the left mate appears first.
mergePair :: Fill -> [BamRow] -> BamRow
mergePair _ []  = error "Empty grouping. Impossible if used with mergeMates\
                      \ unless the input is empty."
mergePair _ [x] = x
mergePair fill [BamRow left, BamRow right] =
    BamRow
        . set (ix 9) ( unSequence
                     $ mergeSequences fill (posSeq left) (posSeq right)
                     )
        $ left
  where
    posSeq xs = ( either error (Position . fst) . T.decimal $ xs !! 3
                , Sequence $ xs !! 9
                )
mergeMate xs = error ("Too many mate pairs for: " <> show xs)

-- | Merge all mate pairs. Assumes that the lines are sorted by mate pair already
-- (samtools sort -n).
mergeMates :: Fill -> [BamRow] -> [BamRow]
mergeMates fill = fmap (mergePair fill) . groupBy ((==) `on` (head . unBamRow))

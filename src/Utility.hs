{- Utility
Gregory W. Schwartz

Collects the functions pertaining to utility helper functions.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Utility
    ( findWithError
    , bamToFasta
    , appendHeader
    , fastaToMap
    ) where

-- Standard
import Data.Bool
import Data.Char
import Data.Maybe
import Data.Monoid
import qualified Data.Map.Strict as Map
import qualified Data.ByteString.Char8 as B

-- Cabal
import qualified Data.Text as T
import Data.Fasta.Text
import Safe
import Turtle.Line

-- Local
import Types

-- | Return a specific error when the text key is not found.
findWithError :: (Eq a, Ord a, Show a) => a -> Map.Map a b -> b
findWithError x = Map.findWithDefault (error (show x <> " not found.")) x

-- | Whether to ignore the read based on the cigar (0 is ignore 1 is don't
-- ignore).
isIgnore :: T.Text -> Bool
isIgnore = (== "M") . T.filter isAlpha

-- | Convert a line in a bam file to a fasta format for duplication input.
bamToFasta :: Maybe AccessionMap
           -> T.Text
           -> Maybe Headers
           -> BamRow
           -> [T.Text]
bamToFasta aMap fileHeader headers (BamRow row) =
    [ ">" <> header
    , maybe (row `at` 9) (lookupErr (row `at` 2) . unAccessionMap) aMap
    ]
  where
    header = T.intercalate "|"
           . catMaybes
           $ [ Just fileHeader
             , Just (row `at` 0)
             , Just (row `at` 3)
             , Just (bool "1" "0" . isIgnore $ row `at` 5)
             , fmap unHeaders headers
             ]

-- | Append a header to the end of the header line in a fasta row.
appendHeader :: Maybe Headers -> Line -> Line
appendHeader Nothing row                              = row
appendHeader (Just (Headers h)) row@((T.take 1 . lineToText) -> ">") =
    unsafeTextToLine $ lineToText row <> "|" <> h
appendHeader _ row                                    = row

-- | Create a mapping from accession number to sequence. Trinity outputs ' ' as
-- a separator.
fastaToMap :: Char -> [FastaSequence] -> AccessionMap
fastaToMap sep =
    AccessionMap . Map.fromList . fmap (\ !x -> (getField 1 sep x, fastaSeq x))

-- | Get a mapping from accession to if the cigar has all matching.
getAllMatchMap :: [BamRow] -> MatchMap
getAllMatchMap = MatchMap
               . Map.fromListWith (&&)
               . fmap (\ (BamRow !x) -> (x `at` 3, isIgnore $ x `at` 5))

-- | Better lookup error message.
lookupErr :: (Ord k, Show k) => k -> Map.Map k a -> a
lookupErr k = fromMaybe (error $ "Cannot find: " <> (show k)) . Map.lookup k

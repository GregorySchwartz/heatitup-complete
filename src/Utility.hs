{- Utility
Gregory W. Schwartz

Collects the functions pertaining to utility helper functions.
-}

{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Utility
    ( findWithError
    , bamToFasta
    , appendHeader
    ) where

-- Standard
import Data.Bool
import Data.Char
import Data.Maybe
import Data.Monoid
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Text as T
import Turtle.Line

-- Local
import Types

-- | Return a specific error when the text key is not found.
findWithError :: (Eq a, Ord a, Show a) => a -> Map.Map a b -> b
findWithError x = Map.findWithDefault (error (show x <> " not found.")) x
              
-- | Whether to ignore the read based on the cigar (0 is ignore 1 is don't
-- ignore).
cigarToIgnore :: T.Text -> T.Text
cigarToIgnore = bool "0" "1" . (== "M") . T.filter isAlpha

-- | Convert a line in a bam file to a fasta format for duplication input.
bamToFasta :: T.Text -> Maybe Headers -> BamRow -> [T.Text]
bamToFasta fileHeader headers (BamRow row) = [">" <> header, row !! 9]
  where
    header = T.intercalate "|"
           . catMaybes
           $ [ Just fileHeader
             , Just (row !! 0)
             , Just (row !! 3)
             , Just (cigarToIgnore $ row !! 5)
             , fmap unHeaders headers
             ]

-- | Append a header to the end of the header line in a fasta row.
appendHeader :: Maybe Headers -> Line -> Line
appendHeader Nothing row                              = row
appendHeader (Just (Headers h)) row@((T.take 1 . lineToText) -> ">") =
    unsafeTextToLine $ lineToText row <> "|" <> h
appendHeader _ row                                    = row


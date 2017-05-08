{- Parse
Gregory W. Schwartz

Collects the functions pertaining to the parsing of the abundance output.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Parse
    ( getAbundanceMap
    , getFrequencyMap
    , getAccSet
    , parseDuplications
    , parseBAM
    ) where

-- Standard
import Data.Char
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set

-- Cabal
import Control.Lens
import Data.Csv
import Data.Text.Read
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Text as T
import qualified Data.Vector as V
import Safe

-- Local
import Types
import Utility

-- | Get the abundance map from the contents of the abundance output.
getAbundanceMap :: B.ByteString -> AbundanceMap
getAbundanceMap =
    AbundanceMap
        . Map.unions
        . fmap (\ !m -> Map.singleton
                            (findWithError ("target_id" :: T.Text) m)
                            (getDouble . findWithError ("est_counts" :: T.Text) $ m)
               )
        . V.toList
        . either error snd
        . decodeByNameWith decodeOpts
  where
    decodeOpts = defaultDecodeOptions { decDelimiter = fromIntegral (ord '\t') }
    getDouble = either error fst . double

-- | Convert an abundance map to a frequency map.
getFrequencyMap :: AbundanceMap -> FrequencyMap
getFrequencyMap (AbundanceMap m) = FrequencyMap . Map.map (/ totalCount) $ m
  where
    totalCount = Map.foldl (+) 0 m

-- | Get the list of accessions from a duplication row post-Trinity.
getAccSet :: [DuplicationRow] -> AccSet
getAccSet = AccSet
          . Set.fromList
          . fmap ( flip (Safe.at) 1
                 . T.splitOn "|"
                 . findWithError "fHeader"
                 . unDuplicationRow
                 )
    
-- | Get the duplications in an easy to read format.
parseDuplications :: B.ByteString -> (Header, [DuplicationRow])
parseDuplications =
    over _2 (fmap DuplicationRow . V.toList) . either error id . decodeByName

-- | Parse bam rows.
parseBAM :: T.Text -> [BamRow]
parseBAM = fmap (BamRow . T.splitOn "\t") . T.lines

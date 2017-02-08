{- Parse
Gregory W. Schwartz

Collects the functions pertaining to the parsing of the abundance output.
-}

{-# LANGUAGE BangPatterns #-}

module Parse
    ( getAbundanceMap
    , getFrequencyMap
    , parseDuplications
    , parseBAM
    ) where

-- Standard
import Data.Char
import qualified Data.Map.Strict as Map

-- Cabal
import Control.Lens
import qualified Data.ByteString.Lazy.Char8 as B
import Data.Csv
import Data.Text.Read
import qualified Data.Vector as V

-- Local
import Types
import Utility

-- | Get the abundance map from the contents of the abundance output.
getAbundanceMap :: B.ByteString -> AbundanceMap
getAbundanceMap =
    AbundanceMap
        . Map.unions
        . fmap (\ !m -> Map.singleton
                            (findWithError "target_id" m)
                            (getDouble . findWithError "est_counts" $ m)
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

-- | Get the duplications in an easy to read format.
parseDuplications :: B.ByteString -> (Header, [DuplicationRow])
parseDuplications = 
    over _2 (fmap DuplicationRow . V.toList) . either error id . decodeByName

-- | Parse bam rows.
parseBAM :: B.ByteString -> [BamRow]
parseBAM =
    fmap BamRow . V.toList . either error id . decodeWith decodeOpts NoHeader
  where
    decodeOpts = defaultDecodeOptions { decDelimiter = fromIntegral (ord '\t') }

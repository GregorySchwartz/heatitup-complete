{- Blast
Gregory W. Schwartz

Collects the functions pertaining to Blasting sequences.
-}

{-# LANGUAGE OverloadedStrings #-}

module Blast
    ( filterBlastBamRows
    , blastBamToFasta
    ) where

-- Standard
import Data.Maybe
import qualified Data.Set as Set

-- Cabal
import qualified Control.Foldl as Fold
import qualified Data.Text as T
import Safe
import Turtle
import Turtle.Line

-- Local
import Types

-- | Get the accessions that successfully align to the reference.
getBlastSet :: Command -> Maybe Args -> Subject -> Query -> Shell BlastSet
getBlastSet (Command cmd) args (Subject sub) query = do
    let queryIn (QueryFile x) = input x
        queryIn (QueryShell x) = x
        argsIn = maybe [] (fmap T.pack . words . unArgs) args
        blastOutput = fmap fromJust
                    . mfilter isJust
                    . fmap (headMay . cut "\t" . lineToText)
                    . inproc cmd ( [ "-subject", format fp sub
                                   , "-outfmt", "6"
                                   , "-max_target_seqs", "1"
                                   ]
                                <> argsIn
                                 )
                    . queryIn
                    $ query
    blastSet <- fmap (BlastSet . Set.fromList) . fold blastOutput $ Fold.list

    return blastSet

-- | Filter out BamRows that do not successfully align to the reference.
filterBlastBamRows :: Int -> Command -> Maybe Args -> Subject -> Query -> [BamRow] -> Shell [BamRow]
filterBlastBamRows accField cmd args sub query rows = do
    blastSet <- getBlastSet cmd args sub $ query

    let newRows =
            filter
                (flip Set.member (unBlastSet blastSet) . (`at` accField) . unBamRow)
                rows

    return newRows

-- | Convert a bam row to fasta for a quick blast check.
blastBamToFasta :: BamRow -> [T.Text]
blastBamToFasta (BamRow row) =
    [ ">" <> (fromMaybe (error "No accession in bam.") . headMay $ row)
    , fromMaybe (error "No sequence in bam.") $ row `atMay` 9
    ]

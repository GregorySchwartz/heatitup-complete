{- find-duplication-complete
Gregory W. Schwartz

Find duplications a sequence with assembly. Requires all to be in path:
find-duplication, Trinity, samtools (if input is a bam), rsem-prepare-reference,
rsem-calculate-expression, and bowtie2. Cleans up at the end by deleting all
files with INPUT.* pattern, the new period being the main difference.
-}

{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Main where

-- Standard
import Data.Bool
import Data.Maybe

-- Cabal
import qualified Control.Foldl as Fold
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Csv as CSV
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Options.Applicative as O
import Turtle
import Turtle.Line

-- Local
import Types
import Parse
import Merge

-- | Command line arguments
data Options = Options { input            :: String
                       , refInput         :: Maybe String
                       , blacklistInput   :: Maybe String
                       , output           :: Maybe String
                       , outputPlot       :: Maybe String
                       , outputLabel      :: String
                       , refField         :: Int
                       , posField         :: Maybe Int
                       , minSize          :: Int
                       , gaussWindow      :: Int
                       , gaussTime        :: Double
                       , gaussThreshold   :: Double
                       , minMut           :: Maybe Int
                       , distance         :: Int
                       , refBlacklistFlag :: Bool
                       , minRichness      :: Int
                       , preprocessType   :: Preprocess
                       , trinityCommand   :: String
                       , trinityArgs      :: String
                       , trinityAbundance :: String
                       , dirtyFlag        :: Bool
                       , headers          :: Maybe String
                       }

-- | Command line options
options :: O.Parser Options
options = Options
      <$> O.strOption
          ( O.long "input"
         <> O.short 'i'
         <> O.metavar "FILE"
         <> O.help "The input file"
          )
      <*> O.optional ( O.strOption
          ( O.long "reference-input"
         <> O.short 'I'
         <> O.metavar "[Nothing] | FILE"
         <> O.help "The input file containing the reference sequences to compare\
                 \ to. The first entry in the field must be the accession and\
                 \ match the requested field from reference-field for the\
                 \ input. With no supplied file, no spacer will be annotated."
          )
        )
      <*> O.optional ( O.strOption
          ( O.long "blacklist-input"
         <> O.short 'b'
         <> O.metavar "[Nothing] | FILE"
         <> O.help "The input fasta file containing possible false positives --\
                 \ sequences which may be duplicate nucleotides in the\
                 \ reference sequence."
          )
        )
      <*> O.optional ( O.strOption
          ( O.long "output"
         <> O.short 'o'
         <> O.metavar "FILE"
         <> O.help "The output file"
          )
        )
      <*> O.optional ( O.strOption
          ( O.long "output-plot"
         <> O.short 'O'
         <> O.metavar "FILE"
         <> O.help "The output file for the plot. Each new plot gets a new number\
                 \ on it: output_1.svg, output_2.svg, etc. Each plot uses the\
                 \ first entry in the fasta header as the label."
          )
        )
      <*> O.strOption
          ( O.long "label"
         <> O.short 'l'
         <> O.metavar "FILE"
         <> O.value ""
         <> O.help "The label to use in the label column for the output"
          )
      <*> O.option O.auto
          ( O.long "reference-field"
         <> O.short 'f'
         <> O.metavar "[1] | INT"
         <> O.value 1
         <> O.help "The field in each input header that contains the reference\
                 \ accession number to compare to."
          )
      <*> O.optional ( O.option O.auto
          ( O.long "position-field"
         <> O.short 'p'
         <> O.metavar "[Nothing] | INT"
         <> O.help "The field in each input header that contains the starting\
                 \ position of the read. Added to the annotations."
          )
        )
      <*> O.option O.auto
          ( O.long "min-size"
         <> O.short 's'
         <> O.metavar "[15] | INT"
         <> O.value 15
         <> O.help "The minimum size of a duplication"
          )
      <*> O.option O.auto
          ( O.long "gaussian-window"
         <> O.short 'w'
         <> O.metavar "[4] | Double"
         <> O.value 3
         <> O.help "The window for the discrete gaussian kernel atypical spacer\
                 \ determination"
          )
      <*> O.option O.auto
          ( O.long "gaussian-time"
         <> O.short 't'
         <> O.metavar "[2] | Double"
         <> O.value 2
         <> O.help "The time for the discrete gaussian kernel atypical spacer\
                 \ determination"
          )
      <*> O.option O.auto
          ( O.long "gaussian-threshold"
         <> O.short 'T'
         <> O.metavar "[0.4] | Double"
         <> O.value 0.4
         <> O.help "The cutoff to be considered a mutation for the discrete\
                 \ gaussian kernel atypical spacer determination"
          )
      <*> O.optional ( O.option O.auto
          ( O.long "min-mutations"
         <> O.short 'm'
         <> O.metavar "INT"
         <> O.help "The minimum number of nucleotides between mutations"
          )
        )
      <*> O.option O.auto
          ( O.long "levenshtein-distance"
         <> O.short 'L'
         <> O.metavar "[2] | INT"
         <> O.value 2
         <> O.help "The minimum Levenshtein distance to the false positive\
                 \ checker. If the distance to the false positive string\
                 \ is less than or equal to this number,\
                 \ the duplication is considered\
                 \ a false positive. Compares candidates against each sequence\
                 \ in --blacklist-input"
          )
      <*> O.switch
          ( O.long "reference-blacklist"
         <> O.short 'r'
         <> O.help "Whether to use the reference as a blacklist in addition to\
                 \ the supplied blacklist. That is, the reference sequences\
                 \ are inputed with the same parameters (except distance, which\
                 \ here is 0)\
                 \ to the duplication finder, and those duplications found are\
                 \ added to the blacklist. This process is recursive, executed\
                 \ until no more duplications are found in the reference.\
                 \ Beware, too many blacklist entries can slow down the finder\
                 \ significantly, as each blacklist entry is compared with each\
                 \ candidate."
          )
      <*> O.option O.auto
          ( O.long "min-richness"
         <> O.short 'R'
         <> O.metavar "[1] | INT"
         <> O.value 0
         <> O.help "The minimum nucleotide richness (number of different types of\
                 \ nucleotides) allowed in the duplication to be considered\
                 \ real. Useful if the user knows that a sequence like\
                 \ \"TTTTTTTTCTTTTTTTTC\" is not likely to be real."
          )
      <*> O.option O.auto
          ( O.long "type"
         <> O.short 'P'
         <> O.metavar "Assembly | NonAssembly CHAR"
         <> O.help "The type of preprocessing before duplication finding.\
                   \ Assembly would be used on exome sequencing, RNA-seq, etc.\
                   \ while NonAssembly would be for certain paired end\
                   \ sequencing like amplicon based. Basically, are the reads\
                   \ fragmented across a location (Assembly) or are they\
                   \ all piled up (NonAssembly)? Paired end is required for\
                   \ NonAssembly, otherwise just use the duplication finding\
                   \ program directly on the reads for best results.\
                   \ NonAssembly additionally requires a character to use as\
                   \ filler between non-overlapping mate pairs, as the\
                   \ program will remove duplications containing that\
                   \ character. Input would look like \"NonAssembly 'X'\"."
          )
      <*> O.strOption
          ( O.long "trinity-command"
         <> O.short 'c'
         <> O.metavar "[Trinity] | PATH"
         <> O.value "Trinity"
         <> O.help "The command used for Trinity. Useful if not in path."
          )
      <*> O.strOption
          ( O.long "trinity-args"
         <> O.short 'C'
         <> O.metavar "[--genome_guided_max_intron 10000\
                      \ --max_memory 10G\
                      \ --full_cleanup\
                      \ --no_version_check\
                      \ --genome_guided_bam]\
                      \ | PATH"
         <> O.value "--genome_guided_max_intron 10000\
                    \ --max_memory 10G\
                    \ --full_cleanup\
                    \ --no_version_check\
                    \ --genome_guided_bam"
         <> O.help "The arguments used for Trinity.\
                 \ Defaults to genome guided bams. Make sure the input argument\
                 \ is last and points to nothing (like in the default)."
          )
      <*> O.strOption
          ( O.long "trinity-abundance"
         <> O.short 'a'
         <> O.metavar "[align_and_estimate_abundance.pl] | PATH"
         <> O.value "align_and_estimate_abundance.pl"
         <> O.help "The command used for align_and_estimate_abundance.pl\
                 \ in Trinity's util folder. Useful if not in path. Make sure\
                 \ kallisto is in the path."
          )
      <*> O.switch
          ( O.long "dirty"
         <> O.short 'd'
         <> O.help "Leave behind the INPUT.* files at the end (but the trinity\
                   \ output is still deleted)."
          )
      <*> O.optional ( O.strOption
          ( O.long "header"
         <> O.short 'H'
         <> O.metavar "[Nothing] | STRING"
         <> O.help "The headers used when converting bam to fasta.\
                   \ For Assembly, the order in the resulting header is, for example\
                   \ \"TRINITY_DN1000|c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]|HEADER\".\
                   \ For NonAssembly, the order in the resulting header is\
                   \ \">FILE|ACCESSION|POSITION|HEADER\". So take that order\
                   \ into account for field options with positions and the\
                   \ rest. Also, make sure this string has fields separated by\
                   \ a pipe \"|\" character. So if you have HEADER as the\
                   \ \"ENSE00000SOMETHING\" reference accession that\
                   \ agrees with --input-reference, in Assembly that would be\
                   \ field 3 while in NonAssembly that would be field 4."
          )
        )

-- | Run the trinity command.
runTrinity :: Options -> Turtle.FilePath -> IO ()
runTrinity opts tempDir = do
    let modifiedArgs = (T.words .  T.pack $ trinityArgs opts)
                    <> [ T.pack . Main.input $ opts
                       , "--output"
                       , format fp tempDir
                       ]

    stderr . inproc (T.pack $ trinityCommand opts) modifiedArgs $ mempty

    return ()

-- | Run samtools if the input is a bam file.
runSamtools :: Options -> Turtle.FilePath -> IO Turtle.FilePath
runSamtools opts tempFasta = do
    let inputFile = T.pack . Main.input $ opts

    fastaFile <- if hasExtension (fromText inputFile) "bam"
                    then strict $ do
                        Turtle.output tempFasta
                            . inproc "samtools" ["fasta", inputFile]
                            $ mempty
                        return . head . textToLines . format fp $ tempFasta
                    else
                        return inputFile

    return . fromText $ fastaFile

-- | Run the RSEM abundance command.
runAbundance :: Options -> Turtle.FilePath -> Turtle.FilePath -> Shell ()
runAbundance opts tempDir trinityFastaFile = do
    tempFasta <- mktempfile tempDir "trinity.fasta"

    inputFasta <- liftIO . runSamtools opts $ tempFasta

    stderr
        . inproc
            (T.pack $ trinityAbundance opts)
            [ "--transcripts", format fp trinityFastaFile
            , "--seqType", "fa"
            , "--single", format fp inputFasta
            , "--est_method", "kallisto"
            , "--trinity_mode"
            , "--prep_reference"
            , "--output_dir", format fp tempDir
            ]
        $ mempty

    return ()

-- | Find duplications in the Trinity output.
runFindDuplication :: Options -> Shell Line -> Shell Line
runFindDuplication opts streamIn =
    inproc "find-duplication" commandList streamIn
  where
    commandList =
        fmap T.pack
            . concat
            . catMaybes
            $ [ fmap (("--reference-input" :) . (:[])) . refInput $ opts
              , fmap (("--blacklist-input" :) . (:[])) . blacklistInput $ opts
              , fmap (("--output" :) . (:[])) . Main.output $ opts
              , fmap (("--output-plot" :) . (:[])) . outputPlot $ opts
              , Just ["--label", outputLabel opts]
              , Just ["--reference-field", show . refField $ opts]
              , fmap (("--position-field" :) . (:[]) . show) . posField $ opts
              , Just ["--min-size", show . minSize $ opts]
              , Just ["--gaussian-window", show . gaussWindow $ opts]
              , Just ["--gaussian-time", show . gaussTime $ opts]
              , Just ["--gaussian-threshold", show . gaussThreshold $ opts]
              , fmap (("--min-mutations" :) . (:[]) . show) . minMut $ opts
              , Just ["--levenshtein-distance", show . distance $ opts]
              , bool Nothing (Just ["--reference-blacklist"])
              . refBlacklistFlag
              $ opts
              , Just ["--min-richness", show . minRichness $ opts]
              ]


-- | Cleanup the intermediate files.
cleanup :: Options -> Shell ()
cleanup opts = sh $ do
    removable <-
        fold ( grep (has . text $ T.pack (Main.input opts) <> ".")
             . join
             . fmap (select . textToLines . format fp)
             . ls
             $ "."
             )
             Fold.list

    mapM_ (rm . fromText . lineToText) removable

    return ()

assembly :: Options -> IO ()
assembly opts = sh $ do
    tempDir <- mktempdir "." "trinity"

    liftIO . runTrinity opts $ tempDir
    trinityFastaFile <-
        fold (find (has "/Trinity" <> suffix ".fasta") tempDir) Fold.head

    unless
        (isNothing trinityFastaFile)
        . runAbundance opts tempDir . fromJust $ trinityFastaFile

    abundances        <-
        strict
            . maybe mempty (Turtle.input . const (tempDir </> "abundance.tsv"))
            $ trinityFastaFile
    duplicationOutput <- strict
                       . runFindDuplication opts
                       . maybe mempty Turtle.input
                       $ trinityFastaFile

    unless (dirtyFlag opts) . cleanup $ opts

    let abundanceMap         = getAbundanceMap . B.pack . T.unpack $ abundances
        frequencyMap         = getFrequencyMap abundanceMap
        (dupHeader, dupRows) =
            parseDuplications . B.pack . T.unpack $ duplicationOutput
        newRows = fmap (mergeAbundance frequencyMap) dupRows
        finalOutput = CSV.encodeByName (V.snoc dupHeader "frequency")
                    . fmap unDuplicationRow
                    $ newRows

    liftIO . B.putStrLn $ finalOutput

-- | Remove rows where the duplication contains the filler.
removeFillDups :: Fill -> WithHeader Line -> Maybe Line
removeFillDups _ (Header x)  = Just x
removeFillDups (Fill fill) (Row h row) =
    if T.isInfixOf (T.singleton fill)
        . fromMaybe (error "dSubstring column not found.")
        . lookup "dSubstring"
        . zip (getCols h)
        . getCols
        $ row
        then Nothing
        else Just row
  where
    getCols = T.splitOn "," . lineToText

nonAssembly :: Options -> Fill -> IO ()
nonAssembly opts fill = sh $ do
    bamOutput <- strict
               . inproc "samtools" ["view", "-"] -- No headers
               . inproc "samtools" [ "sort"
                                   , "-n"
                                   , "-O", "sam"
                                   , T.pack $ Main.input opts
                                   ]
               $ mempty

    let bamRows     = parseBAM . B.pack . T.unpack $ bamOutput
        mergedMates = mergeMates fill bamRows
        fastaOutput = select
                    . concatMap textToLines
                    . concatMap ( bamToFasta
                                  (T.pack $ Main.input opts)
                                  (fmap (Headers . T.pack) . headers $ opts)
                                )
                    $ mergedMates

    stdout
        . fmap fromJust
        . mfilter isJust
        . fmap (removeFillDups fill)
        . header
        . runFindDuplication opts
        $ fastaOutput

mainFunc :: Options -> IO ()
mainFunc opts@(preprocessType -> Assembly)         = assembly opts
mainFunc opts@(preprocessType -> NonAssembly fill) =
    nonAssembly opts (Fill fill)

main :: IO ()
main = O.execParser opts >>= mainFunc
  where
    opts = O.info (O.helper <*> Main.options)
      ( O.fullDesc
     <> O.progDesc "Finds duplications in a sequence with assembly"
     <> O.header "find-duplication-complete, Gregory W. Schwartz" )

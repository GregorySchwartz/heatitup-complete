{- heatitup-complete
Gregory W. Schwartz

Find duplications a sequence with assembly. Requires all to be in path:
heatitup, Trinity, samtools (if input is a bam), rsem-prepare-reference,
rsem-calculate-expression, and bowtie2. Cleans up at the end by deleting all
files with INPUT.* pattern, the new period being the main difference.
-}

{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Main where

-- Standard
import Data.Bool
import Data.Maybe
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set

-- Cabal
import Data.Fasta.Text
import Pipes
import qualified Pipes.Prelude as P
import qualified Pipes.Text as PT
import qualified Pipes.Text.IO as PT
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
import Blast
import Utility

-- | Command line arguments
data Options = Options { input                 :: String
                       , spacerFlag            :: Bool
                       , refInput              :: Maybe String
                       , blacklistInput        :: Maybe String
                       , output                :: Maybe String
                       , outputPlot            :: Maybe String
                       , outputLabel           :: String
                       , refField              :: Int
                       , posField              :: Maybe Int
                       , ignoreField           :: Maybe Int
                       , inputMinSize          :: Int
                       , gaussWindow           :: Int
                       , gaussTime             :: Double
                       , gaussThreshold        :: Double
                       , inputMinMut           :: Maybe Int
                       , inputDistance         :: Int
                       , refCheckBlacklistFlag :: Bool
                       , refRecBlacklistFlag   :: Bool
                       , minRichness           :: Int
                       , inputColorDupL        :: String
                       , inputColorDupR        :: String
                       , inputColorMut         :: String
                       , inputColorSpacer      :: String
                       , inputColorBackground  :: String
                       , inputColorForeground  :: String
                       , preprocessType        :: Preprocess
                       , blastCommand          :: Maybe String
                       , blastArgs             :: Maybe String
                       , trinityCommand        :: String
                       , trinityArgs           :: TrinityArgs
                       , trinityAbundance      :: String
                       , dirtyFlag             :: Bool
                       , cigarFlag             :: Bool
                       , headers               :: Maybe String
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
      <*> O.switch
          ( O.long "spacer"
         <> O.short 's'
         <> O.help "Whether to characterize the spacer. Requires\
                   \ reference-input and matching labels for the input and\
                   \ reference-input to compare."
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
                   \ accession number to compare to. Results in an out of bounds\
                   \ if this field does not exist."
          )
      <*> O.optional ( O.option O.auto
          ( O.long "position-field"
         <> O.short 'p'
         <> O.metavar "[Nothing] | INT"
         <> O.help "The field in each input header that contains the starting\
                   \ position of the read. Added to the annotations. Results\
                   \ in out of bounds if this field does not exist."
          )
        )
      <*> O.optional ( O.option O.auto
          ( O.long "ignore-field"
         <> O.short 'g'
         <> O.metavar "[Nothing] | INT"
         <> O.help "The field in each input header that contains a 0 or a 1:\
                   \ 0 means to ignore this read (assign as Normal) and 1\
                   \ means to find a duplication in this read.\
                   \ Used for reads where there is known to be no duplication\
                   \ and thus helps remove false positives."
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
          ( O.long "reference-check-blacklist"
         <> O.short 'c'
         <> O.help "Whether to use the reference as a blacklist in addition to\
                   \ the supplied blacklist. That is, we check if the duplication\
                   \ can be found twice or more in the reference input."
          )
      <*> O.switch
          ( O.long "reference-recursive-blacklist"
         <> O.short 'r'
         <> O.help "Whether to use the reference as a recursive\
                   \ blacklist in addition to\
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
         <> O.value 1
         <> O.help "The minimum nucleotide richness (number of different types of\
                 \ nucleotides) allowed in the duplication to be considered\
                 \ real. Useful if the user knows that a sequence like\
                 \ \"TTTTTTTTCTTTTTTTTC\" is not likely to be real."
          )
      <*> O.strOption
          ( O.long "color-left-duplication"
         <> O.metavar "[#458588] | COLOR"
         <> O.help "The color of the left side of the repeated sequence."
         <> O.value "#458588"
          )
      <*> O.strOption
          ( O.long "color-right-duplication"
         <> O.metavar "[#b16286] | COLOR"
         <> O.help "The color of the right side of the repeated sequence."
         <> O.value "#b16286"
          )
      <*> O.strOption
          ( O.long "color-difference"
         <> O.metavar "[#cc241d] | COLOR"
         <> O.help "The color of discrepancies between the left and right side of\
                   \ the duplication."
         <> O.value "#cc241d"
          )
      <*> O.strOption
          ( O.long "color-spacer"
         <> O.metavar "[#689d6a] | COLOR"
         <> O.help "The color of the spacer."
         <> O.value "#689d6a"
          )
      <*> O.strOption
          ( O.long "color-background"
         <> O.metavar "[#ebdbb2] | COLOR"
         <> O.help "The color of the background."
         <> O.value "#ebdbb2"
          )
      <*> O.strOption
          ( O.long "color-foreground"
         <> O.metavar "[#282828] | COLOR"
         <> O.help "The color of the foreground."
         <> O.value "#282828"
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
      <*> O.optional ( O.strOption
          ( O.long "blast-command"
         <> O.short 'B'
         <> O.metavar "[Nothing] | PATH"
         <> O.help "The command used for blastn. Useful if not in path.\
                   \ Used for filtering out reads from irrelevant locations.\
                   \ If using small sequences, be sure to set blast-args\
                   \ to \"-task blastn-short\"."
          )
        )
      <*> O.optional ( O.strOption
          ( O.long "blast-args"
         <> O.short 'z'
         <> O.metavar "[Nothing] | STRING"
         <> O.help "The additional arguments used for blastn. Separated by space.\
                   \ Use \"-task blastn-short\" for small sequences."
          )
        )
      <*> O.strOption
          ( O.long "trinity-command"
         <> O.short 'C'
         <> O.metavar "[Trinity] | PATH"
         <> O.value "Trinity"
         <> O.help "The command used for Trinity. Useful if not in path."
          )
      <*> O.option O.auto
          ( O.long "trinity-args"
         <> O.short 'A'
         <> O.metavar "[TrinityBase] | TrinityGenome | TrinityCustom STRING"
         <> O.value TrinityBase
         <> O.help "The arguments used for Trinity.\
                 \ TrinityBase is --seqType fq --run_as_paired --max_memory 10G --no_version_check --single,\
                 \ TrinityGenome is --genome_guided_max_intron 10000 --max_memory 10G --no_version_check --genome_guided_bam,\
                 \ and TrinityCustom STRING is STRING.\
                 \ Make sure the input argument\
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
      <*> O.switch
          ( O.long "cigar-filter"
         <> O.short 'G'
         <> O.help "Skip the CIGAR based filtering, that is, look at all reads\
                   \ and not just reads without all M's in the CIGAR."
          )
      <*> O.optional ( O.strOption
          ( O.long "header"
         <> O.short 'H'
         <> O.metavar "[Nothing] | STRING"
         <> O.help "The headers used when converting bam to fasta.\
                   \ The order in the resulting header is\
                   \ \">FILE|ACCESSION|POSITION|IGNORE|HEADER\".\
                   \ So take that order\
                   \ into account for field options with positions and the\
                   \ rest. Also, make sure this string has fields separated by\
                   \ a pipe \"|\" character. So if you have HEADER as the\
                   \ \"ENSE00000SOMETHING\" reference accession that\
                   \ agrees with --input-reference, that would be\
                   \ field 5."
          )
        )

-- | Get the Trinity arguments.
defTrinityArgs :: TrinityArgs -> T.Text
defTrinityArgs TrinityBase       = T.pack "--seqType fq\
                                          \ --run_as_paired\
                                          \ --max_memory 10G\
                                          \ --no_version_check\
                                          \ --single"
defTrinityArgs TrinityGenome     = T.pack "--genome_guided_max_intron 10000\
                                          \ --max_memory 10G\
                                          \ --no_version_check\
                                          \ --genome_guided_bam"
defTrinityArgs (TrinityCustom x) = T.pack x

-- | Run the trinity command.
runTrinity :: Options -> Turtle.FilePath -> Turtle.FilePath -> IO ()
runTrinity opts tempFile tempDir = do
    let modifiedArgs = (T.words .  defTrinityArgs . trinityArgs $ opts)
                    <> [ format fp tempFile
                       , "--output"
                       , format fp tempDir
                       ]

    stderr . inproc (T.pack $ trinityCommand opts) modifiedArgs $ mempty

    return ()

-- | Run samtools if the input is a bam file.
runSamtools :: Options -> Turtle.FilePath -> IO Turtle.FilePath
runSamtools opts tempFasta = do
    let inputFile = T.pack . Main.input $ opts
        outputType = case preprocessType opts of
                        Assembly        -> "fastq"
                        (NonAssembly _) -> "fasta"

    fastFile <- if hasExtension (fromText inputFile) "bam"
                    then strict $ do
                        Turtle.output tempFasta
                            . inproc "samtools" [outputType, inputFile]
                            $ mempty
                        return . head . textToLines . format fp $ tempFasta
                    else
                        return inputFile

    return . fromText . T.strip $ fastFile

-- | Run the abundance command.
runAlignContig :: Options -> Turtle.FilePath -> Turtle.FilePath -> Shell Line
runAlignContig opts tempDir trinityFastaFile = do
    tempFasta <- mktempfile tempDir "trinity.fasta"

    inputFasta <- liftIO . runSamtools opts $ tempFasta
    fastaMap <-
        fmap (fastaToMap ' ') -- Trinity outputs ' ' as a separator.
            . liftIO
            . PT.runSafeT
            . runEffect
            . P.toListM
            $ pipesFasta (PT.readFile . T.unpack . format fp $ trinityFastaFile)

    bamOutput <- strict
               . inproc "samtools" ["view", "-F", "4", "-"] -- No headers and unmapped reads
               . inproc
                   (T.pack $ trinityAbundance opts)
                   [ "--transcripts", format fp trinityFastaFile
                   , "--seqType", "fa"
                   , "--single", format fp inputFasta
                   , "--est_method", "kallisto"
                   , "--kallisto_add_opts", "--pseudobam"
                   , "--trinity_mode"
                   , "--prep_reference"
                   , "--output_dir", format fp tempDir
                   ]
               $ mempty

    let bamRows' = fmap (BamRow . T.splitOn "\t") . T.lines $ bamOutput

    bamRows <- case blastCommand opts of
                (Just cmd) -> filterBlastBamRows
                                2
                                (Command . T.pack $ cmd)
                                (fmap Args . blastArgs $ opts)
                                ( Subject
                                . fromText
                                . T.pack
                                . fromMaybe (error "Need --reference-input for blast.")
                                . refInput
                                $ opts
                                )
                                (QueryFile trinityFastaFile)
                                bamRows'
                Nothing    -> return bamRows'

    let matchMap = getMatchMap bamRows

    select
        . concatMap textToLines
        . concat
        . nub'
        . fmap
            ( bamToFasta
                (Just fastaMap)
                (bool (Just matchMap) Nothing . cigarFlag $ opts)
                (T.pack $ Main.input opts)
                (fmap (Headers . T.pack) . headers $ opts)
            )
        $ bamRows

-- | Find duplications in the Trinity output.
runFindDuplication :: Options -> Shell Line -> Shell Line
runFindDuplication opts streamIn =
    inproc "heatitup" commandList streamIn
  where
    commandList =
        fmap T.pack
            . concat
            . catMaybes
            $ [ fmap (("--reference-input" :) . (:[])) . refInput $ opts
              , bool Nothing (Just ["--spacer"]) . spacerFlag $ opts
              , fmap (("--blacklist-input" :) . (:[])) . blacklistInput $ opts
              , fmap (("--output" :) . (:[])) . Main.output $ opts
              , fmap (("--output-plot" :) . (:[])) . outputPlot $ opts
              , Just ["--label", outputLabel opts]
              , Just ["--reference-field", show . refField $ opts]
              , fmap (("--position-field" :) . (:[]) . show) . posField $ opts
              , fmap (("--ignore-field" :) . (:[]) . show) . ignoreField $ opts
              , Just ["--min-size", show . inputMinSize $ opts]
              , Just ["--gaussian-window", show . gaussWindow $ opts]
              , Just ["--gaussian-time", show . gaussTime $ opts]
              , Just ["--gaussian-threshold", show . gaussThreshold $ opts]
              , fmap (("--min-mutations" :) . (:[]) . show) . inputMinMut $ opts
              , Just ["--levenshtein-distance", show . inputDistance $ opts]
              , bool Nothing (Just ["--reference-check-blacklist"])
              . refCheckBlacklistFlag
              $ opts
              , bool Nothing (Just ["--reference-recursive-blacklist"])
              . refRecBlacklistFlag
              $ opts
              , Just ["--min-richness", show . minRichness $ opts]
              , Just ["--color-left-duplication", inputColorDupL $ opts]
              , Just ["--color-right-duplication", inputColorDupR $ opts]
              , Just ["--color-difference", inputColorMut $ opts]
              , Just ["--color-spacer", inputColorSpacer $ opts]
              , Just ["--color-background", inputColorBackground $ opts]
              , Just ["--color-foreground", inputColorForeground $ opts]
              ]


-- | Cleanup the intermediate files.
cleanup :: Options -> Shell ()
cleanup opts = sh $ do
    removableHere <-
        fold ( grep (has . text $ T.pack (Main.input opts) <> ".")
             . join
             . fmap (select . textToLines . format fp)
             . ls
             $ "."
             )
             Fold.list
    removableThere <-
        fold ( grep (has . text $ T.pack (Main.input opts) <> ".")
             . join
             . fmap (select . textToLines . format fp)
             . ls
             . directory
             . fromText
             . T.pack
             . Main.input
             $ opts
             )
             Fold.list

    mapM_ (rm . fromText . lineToText) $ removableHere <> removableThere

    return ()

assembly :: Options -> IO ()
assembly opts = sh $ do
    tempDir  <- mktempdir "." "trinity"
    tempFile' <- mktempfile "." "tmp.fq"
    tempFile <-
        case trinityArgs opts of
            TrinityBase       -> liftIO . runSamtools opts $ tempFile'
            TrinityGenome     ->
                return . fromText . T.pack . Main.input $ opts
            (TrinityCustom _) ->
                return . fromText . T.pack . Main.input $ opts

    liftIO . runTrinity opts tempFile $ tempDir
    trinityFastaFile <-
        fold (find (has "/Trinity" <> suffix ".fasta") tempDir) Fold.head

    _ <- if isNothing trinityFastaFile
            then error "Fasta not generated."
            else return ()

    duplicationOutput <- strict
                       . runFindDuplication opts
                       . runAlignContig opts tempDir
                       . fromJust
                       $ trinityFastaFile

    unless (dirtyFlag opts) . cleanup $ opts

    abundances        <-
        strict
            . maybe mempty (Turtle.input . const (tempDir </> "abundance.tsv"))
            $ trinityFastaFile

    unless (dirtyFlag opts) . cleanup $ opts

    let abundanceMap'        = getAbundanceMap . B.pack . T.unpack $ abundances
        (dupHeader, dupRows) =
            parseDuplications . B.pack . T.unpack $ duplicationOutput
        accSet               = getAccSet dupRows
        abundanceMap         = --We only want the frequency of valid contigs
            AbundanceMap
                . Map.filterWithKey (\k _ -> Set.member k . unAccSet $ accSet)
                . unAbundanceMap
                $ abundanceMap'
        frequencyMap         = getFrequencyMap abundanceMap
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

    let bamRows' = parseBAM bamOutput
        query    = QueryShell
                 . select
                 . concatMap textToLines
                 . concatMap blastBamToFasta
                 $ bamRows'

    -- Get rid of sequences that don't blast to the reference.
    bamRows <- case blastCommand opts of
                (Just cmd) -> filterBlastBamRows
                                0
                                (Command . T.pack $ cmd)
                                (fmap Args . blastArgs $ opts)
                                ( Subject
                                . fromText
                                . T.pack
                                . fromMaybe (error "Need --reference-input for blast.")
                                . refInput
                                $ opts
                                )
                                query
                                bamRows'
                Nothing    -> return bamRows'

    let mergedMates = mergeMates fill bamRows
        fastaOutput = select
                    . concatMap textToLines
                    . concatMap ( bamToFasta
                                    Nothing
                                    Nothing
                                    (T.pack $ Main.input opts)
                                    (fmap (Headers . T.pack) . headers $ opts)
                                )
                    $ mergedMates

    stdout
        . fmap fromJust
        . mfilter isJust
        . fmap (removeFillDups fill)
        . Turtle.header
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
     <> O.header "heatitup-complete, Gregory W. Schwartz" )

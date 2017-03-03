{- Types
Gregory W. Schwartz

Collects the types used in the program
-}

module Types where

-- Standard
import qualified Data.IntMap.Strict as IMap
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Text as T

-- Local


-- Algebraic
data Preprocess = Assembly | NonAssembly Char deriving (Read, Show)
data TrinityArgs = TrinityBase | TrinityGenome | TrinityCustom String deriving (Read, Show)

-- Basic
newtype Position = Position
    { unPosition :: Int
    } deriving (Show)
newtype Sequence = Sequence
    { unSequence :: T.Text
    } deriving (Show)
newtype Fill = Fill
    { unFill :: Char
    } deriving (Show)

-- Advanced
newtype BamRow = BamRow
    { unBamRow :: [T.Text]
    } deriving (Read, Show)
newtype Headers = Headers
    { unHeaders :: T.Text
    }
newtype AbundanceMap = AbundanceMap
    { unAbundanceMap :: Map.Map T.Text Double
    }
newtype FrequencyMap = FrequencyMap
    { unFrequencyMap :: Map.Map T.Text Double
    }
newtype DuplicationRow = DuplicationRow
    { unDuplicationRow :: Map.Map T.Text T.Text
    }
newtype PositionMap = PositionMap
    { unPositionMap :: IMap.IntMap Char
    }
newtype AccessionMap = AccessionMap
    { unAccessionMap :: Map.Map T.Text T.Text
    } deriving (Read, Show)
newtype MatchMap = MatchMap
    { unMatchMap :: Map.Map T.Text Bool
    }

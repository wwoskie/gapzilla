"""
Data models (dataclasses) used in module
"""

from dataclasses import dataclass, field, InitVar
from Bio.SeqFeature import SeqFeature


@dataclass
class IntervaledFeature:
    """
    A class to represent a genomic interval with features.

    Attributes
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    feature_list : list of SeqFeature
        List of features within the interval.
    feature_num : int
        Number of features within the interval.
    feature_lengths : list of int
        List of lengths of the features within the interval.

    Parameters
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    feature : SeqFeature
        Initial feature, removed after class instance creation.
    """

    start: int
    """Start position on the genome"""
    end: int
    """End position on the genome"""
    feature: InitVar[SeqFeature]
    """Initial feature, removed after class instance creation"""
    feature_list: list[SeqFeature] = field(init=False)
    """List of features"""
    feature_num: int = field(init=False)
    """Number of features"""
    feature_lengths: list[int] = field(init=False)
    """list of features lengths"""

    def __post_init__(self, feature):
        self.feature_list = [feature]
        self.feature_num = len(self.feature_list)
        self.feature_lengths = [feature.length]


@dataclass
class IntervaledGap:
    """
    A class to represent a gap interval on the genome with associated features on either side.

    Attributes
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    features_left : list
        List of features on the left side of the gap.
    features_right : list
        List of features on the right side of the gap.

    Parameters
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    features_left : list
        List of features on the left side of the gap.
    features_right : list
        List of features on the right side of the gap.
    """

    start: int
    """Start position on the genome"""
    end: int
    """End position on the genome"""
    features_left: list
    """List"""
    features_right: list


@dataclass
class InsertionSite:
    """
    A class to represent an insertion site on the genome.

    Attributes
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    score : int
        Score associated with the insertion site.

    Parameters
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    score : int
        Score associated with the insertion site.
    """

    start: int
    """Start position on the genome"""
    end: int
    """End position on the genome"""
    score: int


@dataclass
class Hairpin:
    """
    A class to represent a hairpin structure on the genome.

    Attributes
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    sequence : str
        RNA sequence for the hairpin.
    structure : str
        Hairpin structure in bracket notation.
    mfe : float, optional
        Minimal free energy (MFE) of the hairpin structure (default is None).

    Parameters
    ----------
    start : int
        Start position on the genome.
    end : int
        End position on the genome.
    sequence : str
        RNA sequence for the hairpin.
    structure : str
        Hairpin structure in bracket notation.
    mfe : float, optional
        Minimal free energy (MFE) of the hairpin structure (default is None).
    """

    start: int
    """Start position on the genome"""
    end: int
    """End position on the genome"""
    sequence: str
    """RNA sequence for hairpin"""
    structure: str
    """Hairpin structure in bracket representation"""
    mfe: float = None
    """Minimal free energy (MFE) of hairpin structure"""

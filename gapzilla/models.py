from dataclasses import dataclass, field
from Bio.SeqFeature import SeqFeature


# class IntervaledFeature:
#     def __init__(self, interval, feature):
#         self.interval = interval
#         self.feature_list = [feature]
#         self.feature_num = len(self.feature_list)
#         self.feature_lengths = [feature.length]


#     def __repr__(self):
#         return f"IntervaledFeature(interval={self.interval}, feature_list={self.feature_list}, feature_num={self.feature_num}, feature_lengths={self.feature_lengths})"
@dataclass
class IntervaledFeature:
    start: int
    end: int
    feature: SeqFeature
    feature_list: list[SeqFeature] = field(init=False)
    feature_num: int = field(init=False)
    feature_lengths: list[int] = field(init=False)

    def __post_init__(self):
        self.feature_list = [self.feature]
        self.feature_num = len(self.feature_list)
        self.feature_lengths = [self.feature.length]


@dataclass
class IntervaledGap:
    start: int
    end: int
    features_left: list
    features_right: list


@dataclass
class InsertionSite:
    start: int
    end: int
    score: int


@dataclass
class Hairpin:
    start: int
    """Start position on the genome"""
    end: int
    """End position on the genome"""
    sequence: str
    """RNA sequence for hairpin"""
    structure: str
    """Hairpin structure in bracket representation"""
    mfe: float
    """Minimal free energy (MFE) of hairpin structure"""

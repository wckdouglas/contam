from enum import Enum

from pydantic.dataclasses import dataclass


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    stop: int


class Genotype(Enum):
    """
    enum for genotype of the variant call
    """

    HET = "HET"  # heterozygous
    HOM = "HOM"  # homozygous


class VariantType(Enum):
    """
    enum for variant type of the variant
    """

    SNV = "SNV"
    INDEL = "INDEL"

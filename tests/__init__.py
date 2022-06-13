from unittest.mock import MagicMock

from diploid_contam.contamination_estimation import VariantPosition
from diploid_contam.models import Genotype, VariantType


class PysamFakeVcf:
    """
    mock pysam.VariantFile object
    """

    def __init__(self, variants):
        """
        a mock object that mimics the pysam.VariantFile object
        :param List variants: reads of the mock vcf file
        """
        self.variants = variants

    def __iter__(self):
        return iter(self.variants)

    def fetch(self, contig, start, end):
        for variant in self.variants:
            if variant.contig == contig and start < variant.start < end:
                yield variant

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


MOCK_VCF_RECORDS = [
    MagicMock(
        contig="chr1",
        start=100,
        filter=["PASS"],
        alleles=("CT", "C"),
        samples=[{"DP": 1152, "AD": (460, 552), "GT": (0, 1)}],
    ),
    MagicMock(
        contig="chr1",
        start=200,
        filter=["PASS"],
        alleles=("CTCCTCTCCTTCCTCT", "C"),
        samples=[{"DP": 1122, "AD": (647, 469), "GT": (0, 1)}],
    ),
    MagicMock(
        contig="chr1",
        start=300,
        filter=["PASS"],
        alleles=("CTCTCCT", "C"),
        samples=[{"DP": 1044, "AD": (604, 437), "GT": (0, 1)}],
    ),
    MagicMock(
        contig="chr1",
        start=1300,
        filter=["RefCall"],
        alleles=("T", "C"),
        samples=[{"DP": 1044, "AD": (554, 466), "GT": (None, None)}],
    ),
    MagicMock(
        contig="chr1",
        start=1600,
        filter=["PASS"],
        alleles=("C", "T"),
        samples=[{"DP": 1131, "AD": (583, 517), "GT": (1, 1)}],
    ),
]
MOCK_VARIANT_POSITIONS = [
    VariantPosition(1152, 552, Genotype.HET, VariantType.INDEL),
    VariantPosition(1122, 469, Genotype.HET, VariantType.INDEL),
    VariantPosition(1044, 437, Genotype.HET, VariantType.INDEL),
    VariantPosition(1044, 466, Genotype.HET, VariantType.SNV),
    VariantPosition(1131, 517, Genotype.HOM, VariantType.SNV),
]

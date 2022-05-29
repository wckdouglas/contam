from unittest.mock import MagicMock, call, patch

import numpy as np
import pytest

from contamination_estimation import (
    CONTAMINATION_RANGE,
    VariantPosition,
    estimate_contamination,
    estimate_vcf_contamination_level,
    maximum_likelihood_contamination,
)


@pytest.mark.parametrize(
    "test_case,total_depth,alt_depth,variant_type,contam_level,expect_log_prob",
    [
        ("HET no contam (0.1)", 100, 50, "HET", 0.1, -3.03339),
        ("HET no contam (0.0)", 100, 50, "HET", 0.0, -2.5308),
        ("HET 10% contam HET-called (0.1)", 100, 10, "HET", 0.1, -2.0259),
        ("HET 10% contam HET-called (0.0)", 100, 10, "HET", 0.0, -38.8323),
        ("HET 10% contam (0.0)", 100, 40, "HET", 0.0, -4.52415),
        ("HET 10% contam (0.1)", 100, 40, "HET", 0.1, -3.01996),
        ("HET 10% contam looks like ALT (0.0)", 100, 60, "HET", 0.0, -4.52415),
        ("HET 10% contam looks like ALT (0.1)", 100, 60, "HET", 0.1, -2.51060),
        ("HET 10% contam a HOM ALT (0.1)", 100, 90, "HET", 0.1, -2.025973),
        ("HET 10% contam a HOM ALT (0.0)", 100, 90, "HET", 0.0, -38.83239),
        ("HOM no contam (0.1)", 100, 100, "HOM", 0.1, -10.536051),
        ("HOM no contam (0.0)", 100, 100, "HOM", 0.0, 0.0),
        ("HOM 10% contam (0.1)", 100, 90, "HOM", 0.1, -2.02597),
        ("HOM 10% contam (0.0)", 100, 90, "HOM", 0.0, -np.inf),
    ],
)
def test_variant_position(test_case, total_depth, alt_depth, variant_type, contam_level, expect_log_prob):
    """
    test function
    """
    vp = VariantPosition(total_depth, alt_depth, variant_type)
    assert np.isclose(vp.log_contam_probability(contam_level), expect_log_prob, rtol=1e-4), f"Failed {test_case}"


def test_estimate_contamination():
    """
    test function
    """
    variant_position1 = MagicMock()
    variant_position2 = MagicMock()
    variant_position1.log_contam_probability.return_value = 1
    variant_position2.log_contam_probability.return_value = 2
    variant_positions = [variant_position1, variant_position2]
    possible_contamination_level = np.arange(CONTAMINATION_RANGE[0], CONTAMINATION_RANGE[1], 0.001)

    out = estimate_contamination(variant_positions)

    for contam in possible_contamination_level:
        assert out[contam] == 1 + 2

    for variant_position in variant_positions:
        for contam in possible_contamination_level:
            variant_position.log_contam_probability.assert_has_calls([call(contam)])


def test_maximum_likelihood_contamination():
    """
    test_function
    """
    with patch("contamination_estimation.estimate_contamination") as mock_estimate_contamination:
        mock_estimate_contamination.return_value = {0.1: 0.1, 0.3: 0.3, 0.2: 0.2}
        result = maximum_likelihood_contamination(["mock", "mock"])
        assert result == 0.3


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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


def test_estimate_vcf_contamination_level(tmp_path):
    """
    test function
    """
    MOCK_VCF_RECORDS = [
        MagicMock(filter=["PASS"], alleles=("CT", "C"), samples=[{"DP": 1152, "AD": (460, 552), "GT": (0, 1)}]),
        MagicMock(
            filter=["PASS"], alleles=("CTCCTCTCCTTCCTCT", "C"), samples=[{"DP": 1122, "AD": (647, 469), "GT": (0, 1)}]
        ),
        MagicMock(filter=["PASS"], alleles=("CTCTCCT", "C"), samples=[{"DP": 1044, "AD": (604, 437), "GT": (0, 1)}]),
        MagicMock(filter=["RefCall"], alleles=("C", "T"), samples=[{"DP": 1043, "AD": (600, 442), "GT": (None, None)}]),
        MagicMock(filter=["PASS"], alleles=("T", "C"), samples=[{"DP": 978, "AD": (534, 441), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("T", "C"), samples=[{"DP": 1038, "AD": (536, 493), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("C", "A"), samples=[{"DP": 1051, "AD": (550, 487), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("T", "C"), samples=[{"DP": 1053, "AD": (547, 476), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("T", "A"), samples=[{"DP": 1054, "AD": (563, 483), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("CTT", "C"), samples=[{"DP": 1065, "AD": (577, 445), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("CT", "C"), samples=[{"DP": 1058, "AD": (597, 459), "GT": (0, 1)}]),
        MagicMock(filter=["PASS"], alleles=("C", "T"), samples=[{"DP": 1053, "AD": (569, 466), "GT": (0, 1)}]),
        MagicMock(filter=["RefCall"], alleles=("T", "C"), samples=[{"DP": 1044, "AD": (554, 466), "GT": (None, None)}]),
        MagicMock(filter=["PASS"], alleles=("TTCC", "T"), samples=[{"DP": 1068, "AD": (600, 453), "GT": (0, 1)}]),
        MagicMock(
            filter=["PASS"],
            alleles=("C", "CCCTCCCCTTCTCCTTCCTCCCCTTCTT"),
            samples=[{"DP": 1100, "AD": (679, 374), "GT": (0, 1)}],
        ),
        MagicMock(filter=["PASS"], alleles=("C", "T"), samples=[{"DP": 1131, "AD": (583, 517), "GT": (0, 1)}]),
    ]

    temp_vcf = tmp_path / "test.vcf"
    temp_vcf.write_text("content")
    with patch("pysam.VariantFile") as mock_variants:
        mock_vcf = PysamFakeVcf(variants=MOCK_VCF_RECORDS)
        mock_variants.return_value.__enter__.return_value = mock_vcf

        assert np.isclose(0.086, estimate_vcf_contamination_level(temp_vcf, snv_only=True))

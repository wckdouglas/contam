import sys
from operator import itemgetter
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import numpy as np
import pytest

from tests import MOCK_VARIANT_POSITIONS, MOCK_VCF_RECORDS, PysamFakeVcf

LIB_PATH = (Path(__file__).absolute().parents[1] / "src").as_posix()
sys.path.append(LIB_PATH)

from contamination_estimation import (
    CONTAMINATION_RANGE,
    VariantPosition,
    collect_variants_from_vcf,
    estimate_contamination,
    estimate_vcf_contamination_level,
    maximum_likelihood_contamination,
)
from models import Interval


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


def test_variant_position_exception():
    """
    test function
    """
    with pytest.raises(ValueError, match="total_depth must be > alt_depth"):
        vp = VariantPosition(10, 11, "HET")


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


@pytest.mark.parametrize(
    "test_case,interval, expected_variants_idx",
    [
        ("without interval", None, [0, 1, 2, 4]),
        ("with 1 interval", [Interval("chr1", 1250, 1701)], [4]),
        ("with 2 intervals", [Interval("chr1", 90, 101), Interval("chr1", 1250, 1701)], [0, 4]),
        ("Interval (RefCall)", [Interval("chr1", 1299, 1301)], None),
    ],
)
def test_collect_variants_from_vcf(tmp_path, test_case, interval, expected_variants_idx):
    """
    test function
    """

    temp_vcf = tmp_path / "test.vcf"
    temp_vcf.write_text("content")
    expected_variants = itemgetter(*expected_variants_idx)(MOCK_VARIANT_POSITIONS) if expected_variants_idx else ()
    if not isinstance(expected_variants, tuple):
        expected_variants = (expected_variants,)
    with patch("contamination_estimation.pysam.VariantFile") as mock_variants, patch(
        "contamination_estimation.pysam.tabix_index"
    ):
        mock_vcf = PysamFakeVcf(variants=MOCK_VCF_RECORDS)
        mock_variants.return_value.__enter__.return_value = mock_vcf
        out_variants = collect_variants_from_vcf(temp_vcf, intervals=interval)
        assert len(expected_variants) == len(out_variants), f"Failed for {test_case}"
        for expected_variant in expected_variants:
            assert expected_variant in out_variants, f"Failed for {test_case}"


@pytest.mark.parametrize(
    "test_case,snv,intervals, expected_variants_idx",
    [
        ("with interval, snv only", True, [Interval("chr1", 1250, 1701)], [4]),
        ("with no interval, snv only", True, None, [4]),
        ("with interval, snv+indel", False, [Interval("chr1", 1250, 1701)], [4]),
        ("with no interval, snv+indel", False, None, [0, 1, 2, 4]),
        ("with no interval, snv+indel", False, [Interval("chr2", 1, 2)], []),
    ],
)
def test_estimate_vcf_contamination_level(tmp_path, test_case, snv, intervals, expected_variants_idx):
    temp_vcf = tmp_path / "test.vcf"
    temp_vcf.write_text("content")
    expected_variants = itemgetter(*expected_variants_idx)(MOCK_VARIANT_POSITIONS) if expected_variants_idx else ()
    if not isinstance(expected_variants, tuple):
        expected_variants = (expected_variants,)
    with patch("contamination_estimation.pysam.VariantFile") as mock_variants, patch(
        "contamination_estimation.pysam.tabix_index"
    ), patch("contamination_estimation.maximum_likelihood_contamination") as mock_func, patch(
        "contamination_estimation.logging"
    ) as mock_log:
        mock_vcf = PysamFakeVcf(variants=MOCK_VCF_RECORDS)
        mock_variants.return_value.__enter__.return_value = mock_vcf
        contam_level = estimate_vcf_contamination_level(temp_vcf, intervals=intervals, snv_only=snv)
        if len(expected_variants_idx) > 0:
            mock_func.assert_called_once_with(list(expected_variants))
        else:
            mock_log.warning.assert_called_once_with("No variants found in vcf file %s", temp_vcf)

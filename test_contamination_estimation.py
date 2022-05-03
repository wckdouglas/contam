from unittest.mock import MagicMock, patch, call

import numpy as np
import pytest

from contamination_estimation import (
    CONTAMINATION_RANGE,
    VariantPosition,
    estimate_contamination,
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
        ("HET 10% contam looks like ALT (0.1)", 100, 60, "HET", 0.1, -3.01996),
        ("HET 10% contam a HOM ALT (0.1)", 100, 90, "HET", 0.1, -2.025973),
        ("HET 10% contam a HOM ALT (0.1)", 100, 90, "HET", 0.0, -38.83239),
        ("HOM no contam (0.1)", 100, 100, "HOM", 0.1, -10.536051),
        ("HOM no contam (0.0)", 100, 100, "HOM", 0.0, 0.0),
        ("HOM 10% contam (0.1)", 100, 90, "HOM", 0.1, -2.02597),
        ("HOM 10% contam (0.0)", 100, 90, "HOM", 0.0, -np.inf),
    ],
)
def test_VariantPosition(test_case, total_depth, alt_depth, variant_type, contam_level, expect_log_prob):
    vp = VariantPosition(total_depth, alt_depth, variant_type)
    assert np.isclose(vp.log_contam_probability(contam_level), expect_log_prob, rtol=1e-4), f"Failed {test_case}"


def test_estimate_contamination():
    variant_position1 = MagicMock()
    variant_position2 = MagicMock()
    variant_positions = [variant_position1, variant_position2]
    possible_contamination_level = np.arange(CONTAMINATION_RANGE[0], CONTAMINATION_RANGE[1], 0.001)

    estimate_contamination(variant_positions)

    for variant_position in variant_positions:
        variant_position.log_contam_probability.assert_has_calls(
            [
                call(contam) for contam in possible_contamination_level
            ]
        )



def test_maximum_likelihood_contamination():
    with patch("contamination_estimation.estimate_contamination") as mock_estimate_contamination:
        mock_estimate_contamination.return_value = {0.1: 0.1, 0.3: 0.3, 0.2: 0.2}
        result = maximum_likelihood_contamination(["mock", "mock"])
        assert result == 0.3

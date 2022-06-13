from __future__ import annotations

import logging
from enum import Enum
from typing import TYPE_CHECKING, Optional

import numpy as np
import pysam
from pydantic import FilePath, conint, validate_arguments, validator

if TYPE_CHECKING:
    from dataclasses import dataclass
else:
    from pydantic.dataclasses import dataclass

from scipy.stats import binom

# https://github.com/liguowang/dcon/blob/master/lib/DconModule/utils.py

CONTAMINATION_RANGE = (0, 0.4)


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


@dataclass(frozen=True)
class VariantPosition:
    """
    data class for a variant

    :param int total_depth: total depth of the variant
    :param int alt_depth: depth of the variant in the ALT allele
    :param Genotype genotype: genotype of the variant (HET or HOM)
    :param Genotype variant_type: type of the variant (SNV or INDEL)
    """

    total_depth: conint(ge=0)  # type: ignore
    alt_depth: conint(ge=0)  # type: ignore
    genotype: Genotype
    variant_type: Optional[VariantType] = None

    @validator("alt_depth")
    def depth_validation(cls, v, values, **kwargs):  # type: ignore
        """
        an extra validation to ensure that the alt_depth is less than the total_depth

        see pydantic documentation for more details
        """
        if values["total_depth"] < v:
            raise ValueError("total_depth must be > alt_depth")
        return v

    def homozygous_log_prob(self, contam_level: float) -> float:
        """
        return log probability of a homozygous variant for a given contamination level

        :param float contam_level: hypothesized contamination level
        :return: log probability of a homozygous variant
        :rtype: float
        """
        expected_alt_fraction = 1 - contam_level  # low AF in ALT because of contam
        max_log_prob = binom.logpmf(
            k=self.alt_depth,
            n=self.total_depth,
            p=expected_alt_fraction,
        )
        return max_log_prob

    def heterozygous_log_prob(self, contam_level: float) -> float:
        """
        return log probability of a heterozygous variant for a given contamination level

        for heterozygous variant, it is a little more complex, because it could be due to:
        1. contamination that doesn't look like the HET ALT allele: we expect lower HET alt allele frequency
        2. contamination that doesn't look like the HOM ALT allele: we expect High HET alt allele frequency
        3. contamination that looks like the ALT allele: we expect higher alt allele frequency
        4. contamination being called as ALT

        :param float contam_level: hypothesized contamination level
        :return: log probability of a heterozygous variant
        :rtype: float
        """
        possibe_expected_alt_fraction = [
            (1 - contam_level) / 2,  # low AF in HET ALT because of contam
            (1 - contam_level),  # this is when a HOM being called as HET because of contam
            (0.5 + contam_level),  # this is when contam looks like ALT
            contam_level,  # this is when the contam is being called as het
        ]
        max_log_prob = max(
            binom.logpmf(
                k=self.alt_depth,
                n=self.total_depth,
                p=expected_alt_fraction,
            )
            for expected_alt_fraction in possibe_expected_alt_fraction
        )
        return max_log_prob

    def log_contam_probability(self, contam_level: float) -> float:
        """
        Estimate the log probabibily of a given contamination level

        :param float contam_level: A
        :return: log likelihood of the contamination level give what we observed at this
            variant position
        :rtype: float
        """
        if not (0 <= contam_level < 1):
            raise ValueError(f"Contamination level must be between 0 and 1: {contam_level}")

        estimator = self.homozygous_log_prob if self.genotype == Genotype.HOM else self.heterozygous_log_prob

        return estimator(contam_level)


def estimate_contamination(
    variant_positions: list[VariantPosition],
) -> dict[float, float]:
    """
    Given a list of SNV position, we will calculate the likelihood at different contamination
    level

    :param List[VariantPosition] variant_positions: a list of VariantPosition object
    :return: a dictionary of contamination level and it's log likeliehood value
    :rtype: dict[float, float]
    """
    possible_contamination_level = np.arange(CONTAMINATION_RANGE[0], CONTAMINATION_RANGE[1], 0.001)
    result = {}

    for contam_level in possible_contamination_level:
        result[contam_level] = sum(
            variant_position.log_contam_probability(contam_level) for variant_position in variant_positions
        )

    return result


def maximum_likelihood_contamination(
    variant_positions: list[VariantPosition],
) -> float:
    """
    Given a list of SNV position, we will estimate the most probable contamination level

    :param List[VariantPosition] variant_positions: a list of VariantPosition object
    :return: a dictionary of contamination level and it's log likeliehood value
    :rtype: float
    """

    likelihoods: dict[float, float] = estimate_contamination(variant_positions)
    sorted_likelihoods: list[tuple[float, float]] = sorted(likelihoods.items(), key=lambda k: k[1])  # ascending sort
    return sorted_likelihoods[-1][0]  # the key of the last item is the max likelihood


@validate_arguments
def collect_variants_from_vcf(vcf_file: FilePath) -> list[VariantPosition]:
    """
    extract variant info from vcf file

    :param FilePath vcf_file: vcf file path
    :return: a list of VariantPosition object
    :rtype: list[VariantPosition]
    """
    variants: list[VariantPosition] = []
    with pysam.VariantFile(vcf_file.as_posix()) as vcf:  # type: ignore
        for variant in vcf:
            if "PASS" in variant.filter or len(variant.filter) == 0:
                total_depth: int = variant.samples[0]["DP"]
                alt: int = variant.samples[0]["GT"][1]
                variant_depth: int = variant.samples[0]["AD"][alt]
                variant_type: VariantType = (
                    VariantType.SNV if len(variant.alleles[1]) == len(variant.alleles[0]) else VariantType.INDEL  # type: ignore
                )
                gt_field: tuple[int, int] = variant.samples[0]["GT"]  # e.g. (0, 1)
                genotype: Genotype = Genotype.HET if len(set(gt_field)) > 1 else Genotype.HOM

                variants.append(
                    VariantPosition(
                        total_depth=total_depth, alt_depth=variant_depth, genotype=genotype, variant_type=variant_type
                    )
                )
    return variants


@validate_arguments
def estimate_vcf_contamination_level(vcf_file: FilePath, snv_only: bool = True) -> float:
    """
    estimate contamination level of a vcf file

    :param Path vcf_file: a vcf file for a sample
    :return: contamination level
    :rtype: float
    """
    variants = collect_variants_from_vcf(vcf_file)
    if snv_only:
        variants = [variant for variant in variants if variant.variant_type == VariantType.SNV]
    logging.debug(f"Processing {len(variants)} variants")
    return maximum_likelihood_contamination(variants)

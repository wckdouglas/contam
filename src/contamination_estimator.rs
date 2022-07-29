use crate::statrs::distribution::{Binomial, Discrete};

use std::vec::Vec;

use crate::model::{VariantPosition, VariantType, Zygosity};

/// Calculate log probability of seeing a number of alt calls
/// at some read depth for a given contamination level
/// # Arguments
/// * `variant_position` - a VariantPosition instance
/// * hypothetical_contamination_level - a conamination level to test for
///
/// # Returns
/// * log probability of seeing the given number of alt calls
///
/// # Examples
///
/// ```
/// let variant = VariantPosition::new(
///     "chr1", 1, 50, 25, VariantType::SNV, Zygosity::HETEROZYGOUS
/// );
/// let log_prob = calc_loglik_for_hypothetical_contam_level(&variant, 0.2);
/// assert_approx_eq!(-3.2073523851, log_prob);
/// ```
fn calc_loglik_for_hypothetical_contam_level(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> f64 {
    let binom = Binomial::new(
        hypothetical_contamination_level,
        variant_position.total_read_depth as u64,
    )
    .unwrap();
    let log_prob = binom.ln_pmf(variant_position.alt_depth as u64);
    log_prob
}

/// Return log probability of a heterozygous variant for a given contamination level
/// for heterozygous variant, it is a little more complex, because it could be due to:
/// 1. contamination that doesn't look like the HET ALT allele: we expect lower HET alt allele frequency
/// 2. contamination that doesn't look like the HOM ALT allele: we expect High HET alt allele frequency
/// 3. contamination that looks like the ALT allele: we expect higher alt allele frequency
/// 4. contamination that looks like the REF allele: we expect lower alt allele frequency
/// 5. contamination being called as ALT
///
/// # Arguments
/// * variant_position - the positional data of the heterozygous variant
/// * hypothetical_contamination_level: hypothetical contamination level
///
/// # Returns
/// * maximum log probability of seeing the given alt depth (across all the tested hypotheses)
fn calc_loglik_for_hypothetical_contam_level_heterozygous(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> f64 {
    let possibe_expected_alt_fraction: Vec<f64> = vec![
        (1.0 - hypothetical_contamination_level) / 2.0, // low AF in HET ALT because of contam doesn't look like ref or alt
        (1.0 - hypothetical_contamination_level), // this is when a HOM being called as HET because of contam
        (0.5 + hypothetical_contamination_level), // this is when contam looks like ALT
        (0.5 - hypothetical_contamination_level), // this is when contam looks like REF
        hypothetical_contamination_level,         // this is when the contam is being called as het
    ];
    let max_log_prob = possibe_expected_alt_fraction
        .into_iter()
        .map(|contam_level| {
            calc_loglik_for_hypothetical_contam_level(variant_position, contam_level)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    max_log_prob
}

/// Helper function to calculate the log probability of a given
/// variant
///
/// # Arguments
/// --------
/// * variant_position: A VariantPosition object to be evaluated
/// * hypothetical_contamination_level: hypothetical contamination level to test
///
/// # Returns
///
/// log probability of seeing the given variant alt coutn
fn calaulate_loglik_for_variant_position(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> f64 {
    let mut log_prob: f64 = 0.0;
    if variant_position.variant_type == VariantType::SNV {
        log_prob = match variant_position.zygosity {
            Zygosity::HOMOZYGOUS => calc_loglik_for_hypothetical_contam_level(
                variant_position,
                1.0 - hypothetical_contamination_level,
            ),
            Zygosity::HETEROZYGOUS => calc_loglik_for_hypothetical_contam_level_heterozygous(
                variant_position,
                hypothetical_contamination_level,
            ),
        }
    }
    log_prob
}

/// Given a list of variant position and a hypothetical contamination level
/// we calculate the log probabilty of seeing the given numbers of alt base across all positions
///
/// # Arguments
///
/// * `variant_vector` - a list of VariantPosition
/// * `hypothetical_contamination_level` - the hypthetical contamination level
///
/// # Returns
///
/// * the sum of log probabilty of seeing the given list of variants at the given contam level
///
/// # Examples
///
/// ```
/// ```
pub fn calculate_contam_hypothesis(
    variant_vector: &Vec<VariantPosition>,
    hypothetical_contamination_level: f64,
) -> f64 {
    if !(0.0..1.0).contains(&hypothetical_contamination_level) {
        panic!("Contamination level must be > 0 and <= 1");
    }

    let log_prob_sum: f64 = variant_vector
        .iter()
        .map(|variant_position| {
            calaulate_loglik_for_variant_position(
                variant_position,
                hypothetical_contamination_level,
            )
        })
        .sum::<f64>();
    log_prob_sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use rstest::*;

    #[rstest]
    #[case(50, 25, 0.2, -3.20735238519, Zygosity::HETEROZYGOUS)]
    #[case(50, 25, 0.3, -4.54456950896, Zygosity::HETEROZYGOUS)]
    #[case(50, 10, 0.3, -1.96740651296, Zygosity::HETEROZYGOUS)]
    #[case(50, 40, 0.1, -4.18755689231, Zygosity::HOMOZYGOUS)] // homozygous
    #[case(50, 40, 0.1, -4.18755689231, Zygosity::HETEROZYGOUS)] // case 2 in HET
    #[case(50, 30, 0.1, -2.16666920827, Zygosity::HETEROZYGOUS)] // case 3 in HET
    #[case(50, 20, 0.1, -2.16666920827, Zygosity::HETEROZYGOUS)] // case 3 in HET
    #[case(50, 5, 0.1,  -1.68780709970, Zygosity::HETEROZYGOUS)] // case 4 in HET
    /// SUT:  calaulate_loglik_for_variant_position
    /// Collaborators:
    ///     - calc_loglik_for_hypothetical_contam_level
    ///     - calc_loglik_for_hypothetical_contam_level_heterozygous
    ///
    fn test_calaulate_loglik_for_variant_position(
        #[case] total_read_depth: usize,
        #[case] alt_depth: usize,
        #[case] hypothetical_contamination_level: f64,
        #[case] expected_out: f64,
        #[case] zygosity: Zygosity,
    ) {
        let variant = VariantPosition::new(
            "X",
            1,
            total_read_depth,
            alt_depth,
            VariantType::SNV,
            zygosity,
        );
        let p = calaulate_loglik_for_variant_position(&variant, hypothetical_contamination_level);
        assert_approx_eq!(p, expected_out);
    }

    #[rstest]
    #[case(0.0,  -2.5308764039)]
    #[case(0.1,  -13.569444762)]
    #[case(0.3,  -42.913904771)]
    fn test_calculate_contam_hypothesis(#[case] contam_level: f64, #[case] expected_log_prob: f64) {
        let variant_list = vec![
            VariantPosition::new("X", 1, 100, 50, VariantType::SNV, Zygosity::HETEROZYGOUS),
            VariantPosition::new("X", 1, 100, 100, VariantType::SNV, Zygosity::HOMOZYGOUS),
        ];
        let log_prob = calculate_contam_hypothesis(&variant_list, contam_level);
        assert_approx_eq!(log_prob, expected_log_prob)
    }
}

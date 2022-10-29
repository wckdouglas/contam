use crate::model::{Hypothesis, VariantPosition, Zygosity};
use lazy_static::lazy_static;
use rayon::prelude::*;
use statrs::distribution::{Binomial, Discrete};
use std::vec::Vec;

lazy_static! {
    static ref HYPOTHESES: [&'static str; 6] = [
        "homozygous",
        "contam is not ref nor alt",
        "contam is called as alt",
        "contam looks like het-alt at hom-alt position",
        "contam comes from a ref-allele",
        "contam looks like het-alt at hom-ref position",
    ];
}

/// Calculate log probability of seeing a number of alt calls
/// at some read depth for a given contamination level
/// # Arguments
/// * `variant_position`: a VariantPosition instance
/// * `hypothetical_contamination_level`: a conamination level to test for
///
/// # Returns
/// * log probability of seeing the given number of alt calls
///
/// # Examples
///
/// ```
/// use assert_approx_eq::assert_approx_eq;
/// use diploid_contam_estimator::contamination_estimator::calc_loglik_for_hypothetical_contam_level;
/// use diploid_contam_estimator::model::{VariantPosition, VariantType, Zygosity};
/// let variant = VariantPosition::new(
///     "chr1", 1, 50, 25, VariantType::SNV, Zygosity::HETEROZYGOUS
/// ).unwrap();
/// let log_prob = calc_loglik_for_hypothetical_contam_level(&variant, 0.2).unwrap();
/// assert_approx_eq!(-13.343980087895773, log_prob);
/// ```
pub fn calc_loglik_for_hypothetical_contam_level(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> Result<f64, String> {
    let binom = Binomial::new(
        hypothetical_contamination_level,
        variant_position.total_read_depth as u64,
    )
    .map_err(|e| e.to_string())?;
    Ok(binom.ln_pmf(variant_position.alt_depth as u64))
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
/// * `variant_position`: the positional data of the heterozygous variant
/// * `hypothetical_contamination_level`: hypothetical contamination level
///
/// # Returns
/// * the best hypothesis with the maximum log probability of seeing the given alt depth (across all the tested hypotheses)
fn calc_loglik_for_hypothetical_contam_level_heterozygous(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> Result<Hypothesis, String> {
    let mut contamination_hypotheses: Vec<Hypothesis> = vec![
        Hypothesis::new(
            HYPOTHESES[1].to_string(),
            (1.0 - hypothetical_contamination_level) / 2.0,
        )?, // low AF in HET ALT because of contam doesn't look like ref or alt
        Hypothesis::new(
            HYPOTHESES[2].to_string(),
            1.0 - hypothetical_contamination_level,
        )?, // this is when a HOM being called as HET because of contam
        Hypothesis::new(
            HYPOTHESES[3].to_string(),
            0.5 + hypothetical_contamination_level,
        )?, // this is when contam looks like ALT
        Hypothesis::new(
            HYPOTHESES[4].to_string(),
            0.5 - hypothetical_contamination_level,
        )?, // this is when contam looks like REF
        Hypothesis::new(HYPOTHESES[5].to_string(), hypothetical_contamination_level)?, // this is when the contam is being called as low vaf het
    ];
    let best_hypothesis = contamination_hypotheses
        .iter_mut()
        .map(|contam_hypothesis| {
            let loglik = calc_loglik_for_hypothetical_contam_level(
                variant_position,
                contam_hypothesis.variant_fraction,
            )
            .unwrap();
            contam_hypothesis.set_loglik(loglik);
            contam_hypothesis
        })
        .max_by(|a, b| a.loglik.partial_cmp(&b.loglik).unwrap())
        .ok_or("MAX is not found in the loglik calculation")?;

    Ok(best_hypothesis.clone())
}

/// Helper function to calculate the log probability of a given
/// variant
///
/// # Arguments
///
/// * `variant_position`: A VariantPosition object to be evaluated
/// * `hypothetical_contamination_level`: hypothetical contamination level to test
///
/// # Returns
///
/// the best hypothesis with the highest log probability of seeing the given variant alt coutn
fn calaulate_loglik_for_variant_position(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> Result<Hypothesis, String> {
    match variant_position.zygosity {
        Zygosity::HOMOZYGOUS => {
            let variant_fraction = 1.0 - hypothetical_contamination_level;
            let loglik =
                calc_loglik_for_hypothetical_contam_level(variant_position, variant_fraction)?;

            let mut best_hypothesis = Hypothesis::new(HYPOTHESES[0].to_string(), variant_fraction)?;
            best_hypothesis.set_loglik(loglik);
            Ok(best_hypothesis.clone())
        }
        Zygosity::HETEROZYGOUS => {
            let best_hypothesis = calc_loglik_for_hypothetical_contam_level_heterozygous(
                variant_position,
                hypothetical_contamination_level,
            )?;
            Ok(best_hypothesis)
        }
    }
}

/// Given a list of variant position and a hypothetical contamination level
/// we calculate the log probabilty of seeing the given numbers of alt base across all positions
///
/// # Arguments
///
/// * `variant_list` - a list of VariantPosition
/// * `hypothetical_contamination_level` - the hypthetical contamination level
///
/// # Returns
///
/// * the sum of log probabilty of seeing the given list of variants at the given contam level
///
/// # Examples
///
/// ```        
/// use assert_approx_eq::assert_approx_eq;
/// use diploid_contam_estimator::contamination_estimator::calculate_contam_hypothesis;
/// use diploid_contam_estimator::model::{VariantPosition, VariantType, Zygosity};
/// let contam_level: f64 = 0.0;
/// let expected_log_prob: f64 = -2.5308764039;
/// let mut variant_list: Vec<VariantPosition> = vec![
///     VariantPosition::new("X", 1, 100, 50, VariantType::SNV, Zygosity::HETEROZYGOUS).unwrap(),
///     VariantPosition::new("X", 1, 100, 100, VariantType::SNV, Zygosity::HOMOZYGOUS).unwrap(),
/// ];
/// let log_prob: f64 = calculate_contam_hypothesis(&mut variant_list, contam_level).unwrap();
/// assert_approx_eq!(log_prob, expected_log_prob)
/// ```
pub fn calculate_contam_hypothesis(
    variant_list: &mut Vec<VariantPosition>,
    hypothetical_contamination_level: f64,
) -> Result<f64, String> {
    if !(0.0..1.0).contains(&hypothetical_contamination_level) {
        return Err("Contamination level must be > 0 and <= 1".to_string());
    }

    // parallel processing of the variant list
    let log_prob_sum = variant_list
        .par_iter_mut()
        .map(|variant_position| {
            // first calculate the best hypothesis and it's respective log likelihood
            // for the given contamination level
            let hyp = calaulate_loglik_for_variant_position(
                variant_position,
                hypothetical_contamination_level,
            )?;
            // transferring the contamination label to the VariantPosition object
            variant_position.set_contamination_label(hyp.label);
            hyp.loglik
                .ok_or_else(|| "loglik not calculated".to_string())
        })
        .sum::<Result<f64, String>>()?;
    Ok(log_prob_sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::VariantType;
    use assert_approx_eq::assert_approx_eq;
    use rstest::*;

    #[rstest]
    #[case(50, 25, 0.2, -3.20735238519, Zygosity::HETEROZYGOUS, HYPOTHESES[1].to_string())]
    #[case(50, 25, 0.3, -4.54456950896, Zygosity::HETEROZYGOUS, HYPOTHESES[1].to_string())]
    #[case(50, 10, 0.3, -1.96740651296, Zygosity::HETEROZYGOUS, HYPOTHESES[4].to_string())]
    #[case(50, 40, 0.1, -4.18755689231, Zygosity::HOMOZYGOUS, HYPOTHESES[0].to_string())] // homozygous
    #[case(50, 40, 0.1, -4.18755689231, Zygosity::HETEROZYGOUS, HYPOTHESES[2].to_string())] // case 2 in HET
    #[case(50, 30, 0.1, -2.16666920827, Zygosity::HETEROZYGOUS, HYPOTHESES[3].to_string())] // case 3 in HET
    #[case(50, 20, 0.1, -2.16666920827, Zygosity::HETEROZYGOUS, HYPOTHESES[4].to_string())] // case 3 in HET
    #[case(50, 5, 0.1,  -1.68780709970, Zygosity::HETEROZYGOUS, HYPOTHESES[5].to_string())] // case 4 in HET
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
        #[case] label: String,
    ) {
        let variant = VariantPosition::new(
            "X",
            1,
            total_read_depth,
            alt_depth,
            VariantType::SNV,
            zygosity,
        )
        .unwrap();
        let p = calaulate_loglik_for_variant_position(&variant, hypothetical_contamination_level)
            .unwrap();
        assert_approx_eq!(p.loglik.unwrap(), expected_out);
        assert_eq!(p.label, label);
    }

    #[rstest]
    #[case(0.0,  -2.5308764039)]
    #[case(0.1,  -13.569444762)]
    #[case(0.3,  -42.913904771)]
    fn test_calculate_contam_hypothesis(#[case] contam_level: f64, #[case] expected_log_prob: f64) {
        let mut variant_list = vec![
            VariantPosition::new("X", 1, 100, 50, VariantType::SNV, Zygosity::HETEROZYGOUS)
                .unwrap(),
            VariantPosition::new("X", 1, 100, 100, VariantType::SNV, Zygosity::HOMOZYGOUS).unwrap(),
        ];
        let log_prob = calculate_contam_hypothesis(&mut variant_list, contam_level);
        assert_approx_eq!(log_prob.unwrap(), expected_log_prob)
    }
}

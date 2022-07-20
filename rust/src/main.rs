extern crate probability;
use crate::probability::distribution::Discrete;
use probability::distribution::Binomial;
use std::vec::Vec;

#[derive(Eq, PartialEq)]
pub enum VariantType {
    SNV,
    INDEL,
}

#[derive(Eq, PartialEq)]
pub enum Zygosity {
    HOMOZYGOUS,
    HETEROZYGOUS,
}

pub struct VariantPosition {
    total_read_depth: usize,
    alt_depth: usize,
    variant_type: VariantType,
    zygosity: Zygosity,
}

/// calculate log probability of seeing a number of alt calls
/// at some read depth for a given contamination level
/// # Arguments
/// * `variant_position` - a VariantPosition instance
/// * hypothetical_contamination_level - a conamination level to test for
///
/// # Returns
/// * log probability of seeing the given number of alt calls
fn calc_loglik_for_hypothetical_contam_level(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> f64 {
    let binom = Binomial::new(
        variant_position.total_read_depth,
        hypothetical_contamination_level.into(),
    );
    let p = Discrete::mass(&binom, variant_position.alt_depth);
    return p.ln();
}


/// return log probability of a heterozygous variant for a given contamination level
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
        .map(|contam| calc_loglik_for_hypothetical_contam_level(variant_position, contam))
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    return max_log_prob;
}

fn calaulate_loglik_for_variant_position(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f64,
) -> f64 {
    let mut p: f64 = 0.0;
    if variant_position.variant_type == VariantType::SNV {
        p = match variant_position.zygosity {
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
    return p;

}


/// Given a list of variant position and a hypothetical contamination level
/// we calculate the log probabilty of seeing the given numbers of alt base across all positions
///
/// # Arguments
///
/// * `variant_vector` - a list of VariantPosition
/// * `hypothetical_contamination_level` - the hypthetical contamination level
fn calculate_contam_hypothesis(
    variant_vector: &Vec<VariantPosition>,
    hypothetical_contamination_level: f64,
) -> f64 {
    if hypothetical_contamination_level < 0.0 || hypothetical_contamination_level >= 1.0 {
        panic!("Contamination level must be > 0 and <= 1");
    }

    let log_prob_sum: f64 = variant_vector
        .into_iter()
        .map(|variant_position| {calaulate_loglik_for_variant_position(variant_position, hypothetical_contamination_level)} )
        .sum::<f64>();
    return log_prob_sum;
}

fn main() {
    let mut variant_vector: Vec<VariantPosition> = Vec::new();
    variant_vector.push(VariantPosition {
        total_read_depth: 100,
        alt_depth: 10,
        variant_type: VariantType::INDEL,
        zygosity: Zygosity::HETEROZYGOUS,
    });
    variant_vector.push(VariantPosition {
        total_read_depth: 100,
        alt_depth: 20,
        variant_type: VariantType::SNV,
        zygosity: Zygosity::HOMOZYGOUS,
    });

    for hypothetical_contamination_level in (1..400).map(|x| x as f64 * 0.001) {
        let p = calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);
        println!("contam {}: {}", hypothetical_contamination_level, p)
    }
}

extern crate probability;
use crate::probability::distribution::Discrete;
use probability::distribution::Binomial;
use std::vec::Vec;

struct VariantPosition {
    total_read_depth: usize,
    alt_depth: usize,
}

/// calculate log probability of seeing a number of alt calls at some read depth for a given contamination level
fn calc_loglik_for_hypothetical_contam_level(
    variant_position: &VariantPosition,
    hypothetical_contamination_level: f32,
) -> f64 {
    let binom = Binomial::new(
        variant_position.total_read_depth,
        hypothetical_contamination_level.into(),
    );
    let p = Discrete::mass(&binom, variant_position.alt_depth);
    return p.ln();
}

fn calculate_contam_hypothesis(
    variant_vector: &Vec<VariantPosition>,
    hypothetical_contamination_level: f32,
) -> f64 {
    let mut log_prob_sum: f64 = 0.0;
    for variant_position in variant_vector.iter() {
        let p = calc_loglik_for_hypothetical_contam_level(
            variant_position,
            hypothetical_contamination_level,
        );
        log_prob_sum += p;
    }
    return log_prob_sum;
}

fn main() {
    let mut variant_vector: Vec<VariantPosition> = Vec::new();
    variant_vector.push(VariantPosition {
        total_read_depth: 100,
        alt_depth: 10,
    });
    variant_vector.push(VariantPosition {
        total_read_depth: 100,
        alt_depth: 20,
    });

    for hypothetical_contamination_level in (1..400).map(|x| x as f32 * 0.001) {
        let p = calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);
        println!("contam {}: {}", hypothetical_contamination_level, p)
    }
}

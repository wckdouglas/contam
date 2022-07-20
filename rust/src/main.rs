use std::vec::Vec;
use std::fs::File;
use std::io::Write;

extern crate rust_htslib;
extern crate serde;
extern crate serde_json;
extern crate statrs;

pub mod model;
use model::VariantPosition;
use serde::{Deserialize, Serialize};

mod contamination_estimator;
use contamination_estimator::calculate_contam_hypothesis;

mod vcfreader;
use vcfreader::build_variant_list;

#[derive(Serialize, Deserialize, Debug)]
struct Output {
    pub contamination_level: f64,
    pub probability: f64,
    pub log_likelihood: f64,
}


fn main() {
    let vcf_file: &str = &"/Users/douglas.wu/code/contam/data/test.vcf";
    let output_json: &str = &"a.json";
    let variant_vector: Vec<VariantPosition> = build_variant_list(vcf_file);
    let mut result_vector: Vec<Output> = Vec::new();
    let mut best_guess_contam_level: f64 = 0.0;
    let mut max_log_likelihood: f64 = 1.0;
    for hypothetical_contamination_level in (1..300).map(|x| x as f64 * 0.001) {
        let p: f64 = calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);
        let output: Output = Output {
            contamination_level: hypothetical_contamination_level,
            probability: p.exp(),
            log_likelihood: p,
        };
        result_vector.push(output);

        if max_log_likelihood > 0.0 || max_log_likelihood < p {
            best_guess_contam_level = hypothetical_contamination_level;
            max_log_likelihood = p;
        }
    }

    let json_string = serde_json::to_string_pretty(&result_vector).unwrap();
    let mut output_file = File::create(output_json).unwrap();
    write!(output_file, "{}", json_string).unwrap();
    println!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );
}

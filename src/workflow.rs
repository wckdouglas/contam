use std::fs::File;
use std::io::Write;
use std::vec::Vec;

extern crate rust_htslib;
extern crate serde;
extern crate serde_json;
extern crate statrs;

use crate::contamination_estimator::calculate_contam_hypothesis;
use crate::model::{ContamProbResult, VariantPosition};
use crate::vcfreader::build_variant_list;

use log::info;

const MAX_CONTAM: usize = 300; // should be 0.399 because we divide 1000

/// write string to file
///
/// Arguments:
/// - filename: the file name of the new file to be written to
/// - json_string: String to be written to the file
pub fn write_json(filename: &str, json_string: String) {
    let mut output_file = File::create(filename).unwrap();
    write!(output_file, "{}", json_string).unwrap();
    info!("Written debug file at: {}", filename)
}

/// the actual workflow to takes in a variant vcf file and calcualte the
/// contamination level
///
/// Arguments:
///
/// - vcf_file: the file path to the input vcf file for the analysis
/// - snv_only_flag: boolean flag indicating whether we should only look at SNV instead of both SNV and indel
/// - depth_threshold: removing all variants with read depth below this threshold
/// - prob_json: for debug, a json file name for writing the contam level and the
///              respecitive log likelihoos into ("_no_file" will turn off writing a file)
/// - prob_json: for debug, a json file name for writing the list of variants that are being
///              used for the contam level compuatation ("_no_file" will turn off writing a file)
pub fn workflow(
    vcf_file: &str,
    snv_only_flag: bool,
    depth_threshold: usize,
    prob_json: &str,
    variant_json: &str,
) -> f64 {
    // collect varaints
    let variant_vector: Vec<VariantPosition> =
        build_variant_list(&*vcf_file, snv_only_flag, depth_threshold);

    // using variants as input to estimate contamination
    let mut result_vector: Vec<ContamProbResult> = Vec::with_capacity(MAX_CONTAM); // initialize a result array to store all result
    let mut best_guess_contam_level: f64 = 0.0; // initialize the final contamination level variable
    let mut max_log_likelihood: f64 = 1.0; // initialize something so that we can track the running max log prob
    for hypothetical_contamination_level in (1..MAX_CONTAM).map(|x| x as f64 * 0.001) {
        // loop over the hypothetical contamination level
        // and calculate the log likelihood
        let log_prob: f64 =
            calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);

        // store them into a result object
        let output: ContamProbResult = ContamProbResult {
            contamination_level: hypothetical_contamination_level,
            log_likelihood: log_prob,
        };
        // and put them in to a result array
        result_vector.push(output);

        if max_log_likelihood > 0.0 || max_log_likelihood < log_prob {
            // if there's a high likelihood contamination level,
            // keep it!
            best_guess_contam_level = hypothetical_contamination_level;
            max_log_likelihood = log_prob;
        }
    }

    if prob_json.ne("no_file") {
        // write result json file
        let json_string = serde_json::to_string_pretty(&result_vector).unwrap();
        write_json(prob_json, json_string)
    }

    if variant_json.ne("no_file") {
        // write variant json file
        let json_string = serde_json::to_string_pretty(&variant_vector).unwrap();
        write_json(variant_json, json_string)
    }

    return best_guess_contam_level;
}

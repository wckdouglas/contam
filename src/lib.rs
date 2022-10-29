pub mod bedreader;
pub mod cli;
pub mod contamination_estimator;
pub mod model;
pub mod vcfreader;

use bedreader::read_bed;
use contamination_estimator::calculate_contam_hypothesis;
use log::info;
use model::{ContamProbResult, VariantPosition};
use std::fs::File;
use std::io::Write;
use std::option::Option;
use std::string::String;
use std::vec::Vec;
use vcfreader::build_variant_list;

const MAX_CONTAM: usize = 400; // should be 0.399 because we divide 1000
const DECIMAL_PLACE: f64 = 0.001; // how precise we want for the contamination level

/// write string to file
///
/// # Arguments:
/// * `filename`: the file name of the new file to be written to
/// * `json_string`: String to be written to the file
pub fn write_json(filename: &str, json_string: String) -> Result<(), String> {
    let mut output_file = File::create(filename).map_err(|e| e.to_string())?;
    write!(output_file, "{}", json_string).unwrap();
    info!("Written debug file at: {}", filename);
    Ok(())
}

/// the actual workflow to takes in a variant vcf file and calcualte the
/// contamination level
///
/// # Arguments:
///
/// * `vcf_file`: the file path to the input vcf file for the analysis
/// * `snv_only_flag`: boolean flag indicating whether we should only look at SNV instead of both SNV and indel
/// * `depth_threshold`: removing all variants with read depth below this threshold
/// * `prob_json`: for debug, a json file name for writing the contam level and the
///              respecitive log likelihoos into ("_no_file" will turn off writing a file)
/// * `prob_json`: for debug, a json file name for writing the list of variants that are being
///              used for the contam level compuatation ("_no_file" will turn off writing a file)
///
/// # Return:
/// * the contamination level with the highest log likelihood
///
/// # Examples:
///
/// ```
/// use diploid_contam_estimator::run;
/// let best_guess_contam = run("data/test.vcf", None, true, 100, Some("prob.json"), Some("variant.json")).unwrap();
/// assert_eq!(best_guess_contam, 0.046);
/// ```
pub fn run(
    vcf_file: &str,
    loci_bed: Option<&str>,
    snv_only_flag: bool,
    depth_threshold: usize,
    prob_json: Option<&str>,
    variant_json: Option<&str>,
) -> Result<f64, String> {
    // collect varaints
    let regions: Vec<String> = match loci_bed {
        Some(bed) => read_bed(bed)?,
        _ => vec![],
    };
    let mut variant_vector: Vec<VariantPosition> =
        build_variant_list(vcf_file, snv_only_flag, depth_threshold, regions)?;

    // using variants as input to estimate contamination
    let mut result_vector: Vec<ContamProbResult> = Vec::with_capacity(MAX_CONTAM); // initialize a result array to store all result
    let mut best_guess: Option<ContamProbResult> = None;
    let contamination_range_to_evaluate = (1..MAX_CONTAM).map(|x| x as f64 * DECIMAL_PLACE);
    for hypothetical_contamination_level in contamination_range_to_evaluate {
        // loop over the hypothetical contamination level
        // and calculate the log likelihood
        let log_prob: f64 =
            calculate_contam_hypothesis(&mut variant_vector, hypothetical_contamination_level)?;

        // store them into a result object
        let output: ContamProbResult = ContamProbResult {
            contamination_level: hypothetical_contamination_level,
            log_likelihood: log_prob,
        };
        // and put them in to a result array
        result_vector.push(output);

        // evaluate whether the newly computed result
        // is better than the previous best one?
        // We will always keep the better guess
        match best_guess {
            None => {
                best_guess = Some(output);
            }
            Some(bg) => {
                if output.log_likelihood > bg.log_likelihood {
                    best_guess = Some(output);
                }
            }
        }
    }
    let best_guess_contam_level = best_guess
        .ok_or("No best guess contam object")?
        .contamination_level;

    // just writing out the result/intermediate files
    if prob_json.is_some() {
        // write result json file
        let json_string =
            serde_json::to_string_pretty(&result_vector).map_err(|e| e.to_string())?;
        write_json(prob_json.ok_or("No prob json name found")?, json_string)?
    }

    if variant_json.is_some() {
        // recalculate loglik
        info!("Adding labels to the variant data json");
        calculate_contam_hypothesis(&mut variant_vector, best_guess_contam_level)?;
        // write variant json file
        let json_string =
            serde_json::to_string_pretty(&variant_vector).map_err(|e| e.to_string())?;
        write_json(
            variant_json.ok_or("No variant json name found")?,
            json_string,
        )?
    }

    Ok(best_guess_contam_level)
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use rstest::*;
    use serde_json::Value;
    use std::io::Read;

    #[rstest]
    #[case(
        false,
        true,
        1000,
        Some("prob.json"),
        Some("variants.json"),
        0.046,
        None
    )]
    #[case(false, true, 1000, None, None, 0.046, None)]
    #[case(false, true, 10, None, None, 0.046, None)]
    #[case(false, true, 10, None, None, 0.046, None)]
    #[case(false, false, 1100, None, None, 0.399, None)]
    #[case(false, true, 1100, None, None, 0.043, None)]
    #[case(true, true, 200, None, None, 0.001, Some("data/test.bed"))] // fetch region from bed
    #[case(true, true, 200, None, None, 0.046, None)] // fetch region from bed
    fn test_run(
        #[case] gz_input: bool,
        #[case] snv_only_flag: bool,
        #[case] depth_threshold: usize,
        #[case] prob_json: Option<&str>,
        #[case] variant_json: Option<&str>,
        #[case] expected_out: f64,
        #[case] bed_file: Option<&str>,
    ) {
        // this is an end to end testing to test everything in
        // the workflow
        let vcf_file = match gz_input {
            false => "data/test.vcf",
            true => "data/test.vcf.gz",
        };
        let best_guess_contam_level: f64 = run(
            vcf_file,
            bed_file,
            snv_only_flag,
            depth_threshold,
            prob_json,
            variant_json,
        )
        .unwrap();
        assert_approx_eq!(best_guess_contam_level, expected_out);
    }

    #[test]
    #[should_panic(expected = "Fetching bed loci from non bgzipped")]
    fn test_workflow_exception() {
        run(
            "data/test.vcf",
            Some("data/test.bed"),
            true,
            100,
            None,
            None,
        )
        .unwrap();
    }

    #[test]
    fn test_write_json() {
        let json_string = "{\"data/test.vcf\":0.046 }";
        write_json("test.json", json_string.to_string()).unwrap();

        let mut file = File::open("test.json").unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();

        let json_data: Value = serde_json::from_str(&data).expect("Bad json data?");
        assert_eq!(json_data["data/test.vcf"], 0.046);
    }
}

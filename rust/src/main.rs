use std::fs::File;
use std::io::Write;
use std::vec::Vec;

use clap::{App, Arg, ArgMatches};

extern crate rust_htslib;
extern crate serde;
extern crate serde_json;
extern crate statrs;

pub mod model;
use model::{ContamProbResult, VariantPosition};

mod contamination_estimator;
use contamination_estimator::calculate_contam_hypothesis;

mod vcfreader;
use vcfreader::build_variant_list;

const PROGRAM_DESC: &'static str = 
    "Estimating contamination level from a diploid VCF file\n\n
    The program assume we are dealing with a diploid genome, and using the 
    deviation of allelic balance from the expected allelic frequence for homozygous
    or heterozygous variant calls to compute a contamination value.

    For homozygous variants, we deviation from allelic frequency of 1 is all introduced by contaminaion.

    For heterozygous variants, it is a little more complex, because it could be due to: 
        1. contamination that doesn't look like the HET ALT allele: we expect lower HET alt allele frequency 
        2. contamination that doesn't look like the HOM ALT allele: we expect High HET alt allele frequency 
        3. contamination that looks like the ALT allele: we expect higher alt allele frequency 
        4. contamination that looks like the REF allele: we expect lower alt allele frequency
        5. contamination being called as ALT
";
const PROGRAM_NAME: &'static str = "diploid-contam-estimator";
const MAX_CONTAM: usize = 300; // should be 0.399 because we divide 1000

/// arg parser to get input from command line
fn parse_args() -> ArgMatches {
    let matches: ArgMatches = App::new(PROGRAM_NAME)
        .version("0.1.0")
        .author("Douglas Wu <wckdouglas@gmail.com>")
        .about(PROGRAM_DESC)
        .arg(
            Arg::with_name("in_vcf")
                .short('i')
                .long("in-vcf")
                .takes_value(true)
                .required(true)
                .help("A diploid vcf file for estimating contamination"),
        )
        .arg(
            Arg::with_name("debug_json")
                .short('d')
                .long("debug-json")
                .takes_value(true)
                .required(false)
                .help("A json output file for storing all intermediate log prob"),
        )
        .arg(
            Arg::with_name("debug_variant_json")
                .short('v')
                .long("debug-variant-json")
                .takes_value(true)
                .required(false)
                .help("A json output file for storing all input variants used for calculation"),
        )
        .arg(
            Arg::with_name("snv_only")
                .long("snv-only")
                .takes_value(false)
                .help("Only use SNV (ignore indel) for contamination estimations"),
        )
        .arg(
            Arg::with_name("depth_threshold")
                .short('m')
                .long("min-depth")
                .takes_value(true)
                .default_value("0")
                .help("Minimum depth for a variant to be considered (i.e. DP tag; default: 0"),
        )
        .get_matches();
    return matches;
}

fn main() {
    // parse cli argumnets
    let args = parse_args();
    let vcf_file: &str = args.value_of::<&str>("in_vcf").unwrap();
    let output_json: &str = args.value_of::<&str>("debug_json").unwrap_or("no_file");
    let output_variant_json: &str = args.value_of::<&str>("debug_variant_json").unwrap_or("no_file");
    let depth_threshold: usize = args.value_of::<&str>("depth_threshold").unwrap_or("0").to_string().parse::<usize>().unwrap();
    let snv_only_flag: bool = args.is_present("snv_only");

    // collect varaints
    let variant_vector: Vec<VariantPosition> = build_variant_list(&*vcf_file, snv_only_flag, depth_threshold);

    // using variants as input to estimate contamination
    let mut result_vector: Vec<ContamProbResult> = Vec::with_capacity(MAX_CONTAM); // initialize a result array to store all result
    let mut best_guess_contam_level: f64 = 0.0; // initialize the final contamination level variable
    let mut max_log_likelihood: f64 = 1.0; // initialize something so that we can track the running max log prob
    for hypothetical_contamination_level in (1..MAX_CONTAM).map(|x| x as f64 * 0.001) {
        // loop over the hypothetical contamination level
        // and calculate the log likelihood
        let log_prob: f64 = calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);

        // store them into a result object
        let output: ContamProbResult = ContamProbResult {
            contamination_level: hypothetical_contamination_level,
            log_likelihood: log_prob,
        };
        result_vector.push(output); // and put them in to a result array

        if max_log_likelihood > 0.0 || max_log_likelihood < log_prob {
            // if there's a high likelihood contamination level,
            // keep it!
            best_guess_contam_level = hypothetical_contamination_level;
            max_log_likelihood = log_prob;
        }
    }

    if output_json.ne("no_file") {
        // write result json file
        let json_string = serde_json::to_string_pretty(&result_vector).unwrap();
        let mut output_file = File::create(output_json).unwrap();
        write!(output_file, "{}", json_string).unwrap();
        println!("Written debug file at: {}", output_json)
    }

    if output_variant_json.ne("no_file") {
        // write variant json file
        let json_string = serde_json::to_string_pretty(&variant_vector).unwrap();
        let mut output_file = File::create(output_variant_json).unwrap();
        write!(output_file, "{}", json_string).unwrap();
        println!("Written debug file at: {}", output_variant_json)
    }

    // this is the resultant number that we want!
    println!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );
}

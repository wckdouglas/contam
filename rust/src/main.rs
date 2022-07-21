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
        .get_matches();
    return matches;
}

fn main() {
    //let vcf_file: &str = &"/Users/wckdouglas/code/contam/data/test.vcf";
    let args = parse_args();
    let vcf_file: &str = args.value_of::<&str>("in_vcf").unwrap();
    let output_json: &str = args.value_of::<&str>("debug_json").unwrap_or("no_file");
    let output_variant_json: &str = args.value_of::<&str>("debug_variant_json").unwrap_or("no_file");
    let snv_only_flag: bool = args.is_present("snv_only");
    let variant_vector: Vec<VariantPosition> = build_variant_list(&*vcf_file, snv_only_flag);
    let mut result_vector: Vec<ContamProbResult> = Vec::new();
    let mut best_guess_contam_level: f64 = 0.0;
    let mut max_log_likelihood: f64 = 1.0;
    for hypothetical_contamination_level in (1..300).map(|x| x as f64 * 0.001) {
        let p: f64 = calculate_contam_hypothesis(&variant_vector, hypothetical_contamination_level);
        let output: ContamProbResult = ContamProbResult {
            contamination_level: hypothetical_contamination_level,
            log_likelihood: p,
        };
        result_vector.push(output);

        if max_log_likelihood > 0.0 || max_log_likelihood < p {
            best_guess_contam_level = hypothetical_contamination_level;
            max_log_likelihood = p;
        }
    }

    if output_json.ne("no_file") {
        let json_string = serde_json::to_string_pretty(&result_vector).unwrap();
        let mut output_file = File::create(output_json).unwrap();
        write!(output_file, "{}", json_string).unwrap();
        println!("Written debug file at: {}", output_json)
    }

    if output_variant_json.ne("no_file") {
        let json_string = serde_json::to_string_pretty(&variant_vector).unwrap();
        let mut output_file = File::create(output_variant_json).unwrap();
        write!(output_file, "{}", json_string).unwrap();
        println!("Written debug file at: {}", output_variant_json)
    }

    println!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );
}

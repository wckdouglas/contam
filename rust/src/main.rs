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

const PROGRAM_DESC: &'static str = "Estimating contamination level from a diploid VCF file";
const PROGRAM_NAME: &'static str = "diploid-contam";

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
        .get_matches();
    return matches;
}

fn main() {
    //let vcf_file: &str = &"/Users/wckdouglas/code/contam/data/test.vcf";
    let args = parse_args();
    let vcf_file: &str = args.value_of::<&str>("in_vcf").unwrap();
    let output_json: &str = args.value_of::<&str>("debug_json").unwrap_or("no_file");
    let variant_vector: Vec<VariantPosition> = build_variant_list(&*vcf_file);
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
    println!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );
}

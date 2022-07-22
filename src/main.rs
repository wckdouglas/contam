use clap::{App, Arg, ArgMatches};

extern crate rust_htslib;
extern crate serde;
extern crate serde_json;
extern crate statrs;

mod contamination_estimator;
mod model;
mod vcfreader;
mod workflow;

use log::info;

use workflow::{workflow, write_json};

const PROGRAM_DESC: &str = "Estimating contamination level from a diploid VCF file\n\n
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
const PROGRAM_NAME: &str = "diploid-contam-estimator";

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
            Arg::with_name("out_json")
                .short('o')
                .long("out-json")
                .takes_value(true)
                .required(false)
                .help("A json output file for storing the maximum likelihood contam level for the vcf file"),
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
                .help("Minimum depth for a variant to be considered (i.e. DP tag)"),
        )
        .get_matches();
    matches
}

fn main() {
    // parse cli argumnets
    env_logger::init();
    let args = parse_args();
    let vcf_file: &str = args.value_of::<&str>("in_vcf").unwrap();
    let prob_json: &str = args.value_of::<&str>("debug_json").unwrap_or("no_file");
    let out_json: &str = args.value_of::<&str>("out_json").unwrap_or("no_file");
    let variant_json: &str = args
        .value_of::<&str>("debug_variant_json")
        .unwrap_or("no_file");
    let depth_threshold: usize = args
        .value_of::<&str>("depth_threshold")
        .unwrap_or("0")
        .to_string()
        .parse::<usize>()
        .unwrap();
    let snv_only_flag: bool = args.is_present("snv_only");

    let best_guess_contam_level: f64 = workflow(
        vcf_file,
        snv_only_flag,
        depth_threshold,
        prob_json,
        variant_json,
    );

    // this is the resultant number that we want!
    info!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );

    if out_json.ne("no_file") {
        let json_string = format!("{{\n  \"{}\": {}\n}}\n", vcf_file, best_guess_contam_level);
        write_json(out_json, json_string)
    }
}

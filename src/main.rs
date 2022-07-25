extern crate noodles_bgzf;
extern crate noodles_tabix;
extern crate noodles_vcf;
extern crate serde;
extern crate serde_json;
extern crate statrs;

mod bed;
mod cli;
mod contamination_estimator;
mod model;
mod vcfreader;
mod workflow;

use cli::parse_args;
use log::info;
use std::option::Option;
use workflow::{workflow, write_json};

fn main() {
    // parse cli argumnets
    env_logger::init();
    let args = parse_args();
    let vcf_file: &str = args.value_of::<&str>("in_vcf").unwrap();
    let prob_json: Option<&str> = args.value_of::<&str>("debug_json");
    let out_json: Option<&str> = args.value_of::<&str>("out_json");
    let variant_json: Option<&str> = args.value_of::<&str>("debug_variant_json");
    let loci_bed: Option<&str> = args.value_of::<&str>("loci_bed");
    let depth_threshold: usize = args
        .value_of::<&str>("depth_threshold")
        .unwrap_or("0")
        .to_string()
        .parse::<usize>()
        .unwrap();
    let snv_only_flag: bool = args.is_present("snv_only");

    let best_guess_contam_level: f64 = workflow(
        vcf_file,
        loci_bed,
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

    if out_json.is_some() {
        let json_string = format!(
            "{{\n  \"vcf_file\":\"{}\",\n\"contamination_level\": {}\n}}\n",
            vcf_file,
            best_guess_contam_level * 100.0
        );
        write_json(out_json.unwrap(), json_string);
        info!("Written result json at: {}", out_json.unwrap());
    }
}

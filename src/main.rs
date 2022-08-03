use diploid_contam_estimator::cli::parse_args;
use diploid_contam_estimator::{run, write_json};
use log::info;
use serde_json::json;
use std::option::Option;

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

    let best_guess_contam_level: f64 = run(
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
        let json_data = json!(
            {
                "vcf_file": vcf_file,
                "contamination_percentage": best_guess_contam_level * 100.0,
            }
        );
        write_json(
            out_json.unwrap(),
            serde_json::to_string_pretty(&json_data).unwrap(),
        );
        info!("Written result json at: {}", out_json.unwrap());
    }
}

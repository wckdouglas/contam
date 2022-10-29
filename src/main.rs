use diploid_contam_estimator::cli::parse_args;
use diploid_contam_estimator::{run, write_json};
use log::info;
use serde_json::json;

pub fn wrapper() -> Result<i8, String> {
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
    )?;
    info!(
        "Maximum likelihood contamination level: {}",
        best_guess_contam_level
    );

    if out_json.is_some() {
        let out_json_file = out_json.ok_or("JSON filename is not given")?;
        let json_data = json!(
            {
                "vcf_file": vcf_file,
                "contamination_percentage": best_guess_contam_level * 100.0,
            }
        );
        write_json(
            out_json_file,
            serde_json::to_string_pretty(&json_data).map_err(|e| e.to_string())?,
        )?;
        info!("Written result json at: {}", out_json_file);
    }
    Ok(0)
}

fn main() {
    // parse cli argumnets
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let result = wrapper();
    if let Err(e) = result {
        println!("Error: {}", e)
    }
}

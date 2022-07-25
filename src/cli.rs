use clap::{App, Arg, ArgMatches};

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
pub fn parse_args() -> ArgMatches {
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
        .arg(
            Arg::with_name("loci_bed")
                .short('b')
                .long("bed")
                .takes_value(true)
                .required(false)
                .help("bed file containing loci for extracting variants"),
        )
        .get_matches();
    matches
}

use crate::model::{VariantPosition, VariantType, Zygosity};

use log::info;
use noodles_vcf as vcf;
use noodles_vcf::header::format::Key;
use noodles_vcf::record::filters::Filters;
use noodles_vcf::record::genotypes::genotype::field::Value::{Integer, IntegerArray};
use noodles_vcf::record::Record;
use std::fs::File;
use std::io::BufReader;
use std::vec::Vec;

/// Colelcting variants from a vcf file
///
/// Arguments:
/// - vcf_file: file path to the vcf file we want to parse
/// - snv_only_flag: boolean flag indicating whether we shopuld only look at SNV instead of both SNV and indel
/// - depth_threshold: we will skip any variants with DP tag lower than this threshold
///
/// Returns:
/// - a list of varaints that passed the given filters
pub fn build_variant_list(
    vcf_file: &str,
    snv_only_flag: bool,
    depth_threshold: usize,
) -> Vec<VariantPosition> {
    let mut reader = File::open(vcf_file)
        .map(BufReader::new)
        .map(vcf::Reader::new)
        .unwrap();
    let raw_header = reader.read_header().expect("Error reading header");
    let header = raw_header.parse().unwrap();

    let mut variants: Vec<VariantPosition> = Vec::new();
    for result in reader.records(&header) {
        let record: Record = result.expect("Cannot read vcf record");

        if record.filters().unwrap().eq(&Filters::Pass) {
            // only look at pass filter variants

            let sample_genotype = record.genotypes().get(0).expect("Error out Alelle 1");
            let read_depth = match sample_genotype[&Key::ReadDepth].value().expect("DP tag") {
                Integer(n) => *n,
                _ => -1,
            };

            if read_depth >= depth_threshold as i32 {
                let bad_vec = &vec![None];
                let allele_depths = match sample_genotype[&Key::ReadDepths].value().expect("AD tag")
                {
                    IntegerArray(n) => n,
                    _ => bad_vec,
                };
                // Genotyping sample
                let gt = sample_genotype.genotype().unwrap().unwrap();
                let ref_genotype = gt[0].position().unwrap();
                let alt_genotype = gt[1].position().unwrap();

                let mut zygosity = Zygosity::HOMOZYGOUS;
                if ref_genotype != alt_genotype {
                    zygosity = Zygosity::HETEROZYGOUS
                }
                // assume theres only one sample in the vcf file hence:  get(0)
                // and diploid call (2nd genotype is non-ref), hence: [1].index
                let ref_base = record.reference_bases();
                let alt_base = &record.alternate_bases()[alt_genotype - 1];
                let alt_depth = allele_depths[alt_genotype]
                    .expect("Alt allele depth is unavaliable (AD tag)")
                    as usize;

                let mut variant_type: VariantType = VariantType::INDEL;
                if ref_base.to_string().len() == alt_base.to_string().len() {
                    // this should be testing the len of Vec<u8> where
                    // each item represents a base
                    // only if they are the same length, it's a SNV
                    variant_type = VariantType::SNV;
                }

                if !snv_only_flag || (snv_only_flag && variant_type == VariantType::SNV) {
                    // whether we want snv-only or not
                    // make a new VariantPosition here and put into the list
                    variants.push(VariantPosition::new(
                        &record.chromosome().to_string(),
                        usize::try_from(record.position()).unwrap(),
                        read_depth as usize, // only sample in the vcf
                        alt_depth,
                        variant_type,
                        zygosity,
                    ));
                }
            }
        }
    }
    info!("Collected {} variants from {}", variants.len(), vcf_file);
    return variants;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(false, 0, 14)] // all variants
    #[case(true, 0, 7)] // all SNV
    #[case(true, 1000, 6)] // all high depth SNV
    #[case(true, 1100, 1)] // all high depth SNV
    #[case(true, 1200, 0)] // all high depth SNV
    fn test_build_variant_list(
        #[case] snv_only_flag: bool,
        #[case] depth_threshold: usize,
        #[case] expected_number_variants: usize,
    ) {
        let vcf_file = "data/test.vcf";
        let variant_list = build_variant_list(&vcf_file, snv_only_flag, depth_threshold);
        assert_eq!(variant_list.len(), expected_number_variants);
    }

    #[rstest]
    #[case(
        0,
        Zygosity::HETEROZYGOUS,
        552,
        1152,
        38144667,
        "X",
        VariantType::INDEL
    )]
    #[case(
        1,
        Zygosity::HETEROZYGOUS,
        469,
        1122,
        38145129,
        "X",
        VariantType::INDEL
    )]
    #[case(3, Zygosity::HETEROZYGOUS, 441, 978, 38145492, "X", VariantType::SNV)]
    #[case(13, Zygosity::HOMOZYGOUS, 517, 1131, 38145911, "X", VariantType::SNV)]
    fn test_build_variant_list_constructed_variant_position(
        #[case] record_idx: usize,
        #[case] zygosity: Zygosity,
        #[case] alt_depth: usize,
        #[case] total_read_depth: usize,
        #[case] position: usize,
        #[case] contig: String,
        #[case] variant_type: VariantType,
    ) {
        let vcf_file = "data/test.vcf";
        let variant_list = build_variant_list(&vcf_file, false, 0);
        let record = &variant_list[record_idx];
        assert_eq!(record.zygosity, zygosity);
        assert_eq!(record.alt_depth, alt_depth);
        assert_eq!(record.total_read_depth, total_read_depth);
        assert_eq!(record.position, position);
        assert_eq!(record.contig, contig);
        assert_eq!(record.variant_type, variant_type);
    }
}

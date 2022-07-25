use crate::model::{VariantPosition, VariantType, Zygosity};

use log::{info, warn};
use noodles_bgzf as bgzf;
use noodles_tabix as tabix;
use noodles_vcf as vcf;
use noodles_vcf::header::format::Key;
use noodles_vcf::record::filters::Filters;
use noodles_vcf::record::genotypes::genotype::field::Value::{Integer, IntegerArray};
use noodles_vcf::record::Record;
use std::fs::{metadata, File};
use std::io::BufReader;
use std::vec::Vec;

/// Evaluate a vcf record and determine whether it should be collected
/// for estimating contamination
///
/// Arguments
/// --------
/// - record: a vcf record from noodles_vcf
/// - variants: a mutable list containing the accepted variants
/// - depth_threshold: if the variant has DP tag lower than this, it will be rejected
/// - snv_only_flag: boolean flag indicating whether we should skip all InDel variants
fn filter_variants(
    record: &Record,
    depth_threshold: usize,
    snv_only_flag: bool,
) -> Option<VariantPosition> {
    if record.filters().is_none() || record.filters().unwrap().eq(&Filters::Pass) {
        // only look at pass filter variants

        let sample_genotype = record.genotypes().get(0).expect("Error out Alelle 1");
        let read_depth = match sample_genotype[&Key::ReadDepth].value().expect("DP tag") {
            Integer(n) => *n,
            _ => 0,
        };

        if read_depth >= depth_threshold as i32 {
            let bad_vec = &vec![None];
            let allele_depths = match sample_genotype[&Key::ReadDepths].value().expect("AD tag") {
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
                filtered_out = false;
                return VariantPosition::new(
                    &record.chromosome().to_string(),
                    usize::try_from(record.position()).unwrap(),
                    read_depth as usize, // only sample in the vcf
                    alt_depth,
                    variant_type,
                    zygosity,
                );
            }
        }
    }
}

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
    regions: Vec<String>,
) -> Vec<VariantPosition> {
    let mut variants: Vec<VariantPosition> = Vec::new();
    let is_gz_input = vcf_file.ends_with(".gz");
    let is_fetch: bool = regions.len() > 0;
    match (is_gz_input, is_fetch) {
        (true, true) => {
            let vcf_file_idx_fn = format!("{}.tbi", vcf_file);
            if metadata(&vcf_file_idx_fn).is_ok() {
                let index = tabix::read(vcf_file_idx_fn).expect(
                    "Error reading tabix file, tabix file is needed if a bed file is supplied",
                );

                let mut reader = File::open(vcf_file)
                    .map(bgzf::Reader::new)
                    .map(vcf::Reader::new)
                    .unwrap();

                let raw_header = reader.read_header().expect("Error reading header");
                let header = raw_header.parse().unwrap();
                for region in regions.iter() {
                    info!("Fetching {}", region);
                    let query = reader.query(
                        &header,
                        &index,
                        &region.parse().expect("Error parsing locus string"),
                    );

                    if query.is_ok() {
                        // it will not spit error only when there are record in the vcf file within the given locus
                        for record in query
                            .unwrap()
                            .map(|result| result.expect("Cannot read vcf record"))
                        {
                            filter_variants(&record, &mut variants, depth_threshold, snv_only_flag);
                        }
                    }
                }
            } else {
                panic!("Missing tabix file {}", vcf_file_idx_fn);
            }
        }
        (true, false) => {
            // in the case of gzip-ed vcf input, we will read all variants
            let mut reader = File::open(vcf_file)
                .map(bgzf::Reader::new)
                .map(BufReader::new)
                .map(vcf::Reader::new)
                .unwrap();

            let raw_header = reader.read_header().expect("Error reading header");
            let header = raw_header.parse().unwrap();
            for record in reader
                .records(&header)
                .map(|result| result.expect("Cannot read vcf record"))
            {
                filter_variants(&record, &mut variants, depth_threshold, snv_only_flag);
            }
        }
        _ => {
            if is_fetch {
                warn!("Because the input vcf is not indexed, will use all the variants!!");
            }
            // in the case of non gz vcf input, we will read all variants
            let mut reader = File::open(vcf_file)
                .map(BufReader::new)
                .map(vcf::Reader::new)
                .unwrap();

            let raw_header = reader.read_header().expect("Error reading header");
            let header = raw_header.parse().unwrap();
            for record in reader
                .records(&header)
                .map(|result| result.expect("Cannot read vcf record"))
            {
                filter_variants(&record, &mut variants, depth_threshold, snv_only_flag);
            }
        }
    };

    info!("Collected {} variants from {}", variants.len(), vcf_file);
    return variants;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(false, 0, 14, vec![])] // all variants
    #[case(true, 0, 7, vec![])] // all SNV
    #[case(true, 1000, 6, vec![])] // all high depth SNV
    #[case(true, 1100, 1, vec![])] // all high depth SNV
    #[case(true, 1200, 0, vec![])] // all high depth SNV
    fn test_build_variant_list(
        #[case] snv_only_flag: bool,
        #[case] depth_threshold: usize,
        #[case] expected_number_variants: usize,
        #[case] regions: Vec<String>,
    ) {
        let vcf_file = "data/test.vcf";
        let variant_list = build_variant_list(&vcf_file, snv_only_flag, depth_threshold, regions);
        assert_eq!(variant_list.len(), expected_number_variants);
    }

    #[rstest]
    #[case(false, 0, 14, vec![])] // all variants
    #[case(true, 0, 7, vec![])] // all SNV
    #[case(true, 1000, 6, vec![])] // all high depth SNV
    #[case(true, 1100, 1, vec![])] // all high depth SNV
    #[case(true, 1200, 0, vec![])] // all high depth SNV
    #[case(false, 10, 1, vec![String::from("X:38144665-38144669")])] // test fetch
    #[case(false, 10, 7, vec![String::from("X:38145491-38145540")])] // test fetch
    #[case(false, 10, 8, vec![String::from("X:38144665-38144669"), String::from("X:38145491-38145540")])] // test fetch
    #[case(false, 10, 0, vec![String::from("1:38145491-38145540")])] // test fetch somewhere else, no error out
    fn test_build_variant_list_from_vcf_gz(
        #[case] snv_only_flag: bool,
        #[case] depth_threshold: usize,
        #[case] expected_number_variants: usize,
        #[case] regions: Vec<String>,
    ) {
        let vcf_file = "data/test.vcf.gz";
        let variant_list = build_variant_list(&vcf_file, snv_only_flag, depth_threshold, regions);
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
    #[case(12, Zygosity::HOMOZYGOUS, 374, 1100, 38145619, "X", VariantType::INDEL)]
    #[case(
        11,
        Zygosity::HETEROZYGOUS,
        453,
        1068,
        38145582,
        "X",
        VariantType::INDEL
    )]
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
        let variant_list = build_variant_list(&vcf_file, false, 0, vec![]);
        let record = &variant_list[record_idx];
        assert_eq!(record.zygosity, zygosity);
        assert_eq!(record.alt_depth, alt_depth);
        assert_eq!(record.total_read_depth, total_read_depth);
        assert_eq!(record.position, position);
        assert_eq!(record.contig, contig);
        assert_eq!(record.variant_type, variant_type);
    }
}

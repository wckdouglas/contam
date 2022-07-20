use crate::model::{VariantPosition, VariantType, Zygosity};
use crate::rust_htslib::bcf::{Reader, Read, Record};
use crate::rust_htslib::bcf::record::FilterId;
use std::vec::Vec;

pub fn build_variant_list(vcf_file: &str) -> Vec<VariantPosition> {
    let mut bcf: Reader = Reader::from_path(vcf_file).expect("Error opening file.");
    let mut variants: Vec<VariantPosition> = Vec::new();
    for (_i, record_result) in bcf.records().enumerate() {
        let record: Record = record_result.expect("Fail to read record");

        if record.filters().map(|filter| filter.is_pass() ).all(|x| !!x) {
            // only look at pass filter variants
            
            let read_depth = record.format(b"DP").integer().ok().expect("Error reading DP field.");
            let allele_depth = record.format(b"AD").integer().ok().expect("Error reading AD field.");
            let gts = record.genotypes().expect("Error reading genotypes");
            let alt_call = gts.get(0)[1].index().unwrap() as usize; // from the only sample, and diploid call (2nd genotype is non-ref)

            variants.push(
                VariantPosition {
                    total_read_depth: read_depth[0][0] as usize, // only sample in the vcf
                    alt_depth: allele_depth[0][alt_call] as usize,
                    variant_type: VariantType::SNV,
                    zygosity: Zygosity::HOMOZYGOUS,
                }
            );
        }
    }
    eprintln!("Collected {} variants from {}", variants.len(), vcf_file);
    return variants;
}

use crate::model::{VariantPosition, VariantType, Zygosity};
use crate::rust_htslib::bcf::record::FilterId;
use crate::rust_htslib::bcf::{Read, Reader, Record};
use std::str::from_utf8;
use std::vec::Vec;

pub fn build_variant_list(vcf_file: &str, snv_only_flag: bool) -> Vec<VariantPosition> {
    let mut vcf: Reader = Reader::from_path(vcf_file).expect("Error opening file.");

    let mut variants: Vec<VariantPosition> = Vec::new();
    for (_i, record_result) in vcf.records().enumerate() {
        let record: Record = record_result.expect("Fail to read record");

        if record.filters().map(|filter| filter.is_pass()).all(|x| !!x) {
            // only look at pass filter variants

            let read_depth = record
                .format(b"DP")
                .integer()
                .ok()
                .expect("Error reading DP field.");
            let allele_depth = record
                .format(b"AD")
                .integer()
                .ok()
                .expect("Error reading AD field.");
            let gts = record.genotypes().expect("Error reading genotypes");
            let alt_call = gts.get(0)[1].index().unwrap() as usize; // from the only sample, and diploid call (2nd genotype is non-ref)

            let mut variant_type: VariantType = VariantType::INDEL;
            if record.alleles()[0].len() == record.alleles()[1].len() {
                variant_type = VariantType::SNV;
            }

            let mut zygosity = Zygosity::HOMOZYGOUS;
            if gts.get(0)[0] != gts.get(0)[1] {
                zygosity = Zygosity::HETEROZYGOUS;
            }

            if !snv_only_flag || (snv_only_flag && variant_type == VariantType::SNV) {
                variants.push(VariantPosition::new(
                    from_utf8(record.header().rid2name(record.rid().unwrap()).unwrap()).unwrap(),
                    record.pos(),
                    read_depth[0][0] as usize, // only sample in the vcf
                    allele_depth[0][alt_call] as usize,
                    variant_type, // TODO: fix this
                    zygosity,
                ));
            }
        }
    }
    eprintln!("Collected {} variants from {}", variants.len(), vcf_file);
    return variants;
}

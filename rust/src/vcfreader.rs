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
            // assume theres only one sample in the vcf file hence:  get(0)
            // and diploid call (2nd genotype is non-ref), hence: [1].index
            let alt_call = gts.get(0)[1].index().unwrap() as usize;

            let mut variant_type: VariantType = VariantType::INDEL;
            if record.alleles()[0].len() == record.alleles()[1].len() {
                // this should be testing the len of Vec<u8> where
                // each item represents a base
                // only if they are the same length, it's a SNV
                variant_type = VariantType::SNV;
            }

            let mut zygosity = Zygosity::HOMOZYGOUS;
            if gts.get(0)[0] != gts.get(0)[1] {
                // if the first genotype != the second genotype
                // e.g. 0/1, 0/2
                // then it's a heterozygous
                zygosity = Zygosity::HETEROZYGOUS;
            }

            if !snv_only_flag || (snv_only_flag && variant_type == VariantType::SNV) {
                // whether we want snv-only or not
                // make a new VariantPosition here and put into the list
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

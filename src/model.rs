/// Defining models for the code
///
use serde::{Deserialize, Serialize};
use std::string::String;

#[derive(Eq, PartialEq, Debug, Serialize, Deserialize)]
/// variant types in the vcf file
pub enum VariantType {
    /// a single nucleotide variant
    SNV,
    /// insertion or deletion
    INDEL,
}

#[derive(Eq, PartialEq, Debug, Serialize, Deserialize)]
/// zygostiy of a variant
pub enum Zygosity {
    /// a homozygous variant
    HOMOZYGOUS,
    /// a heterozygous variant
    HETEROZYGOUS,
}

#[derive(Serialize, Deserialize, Debug)]
/// A struct to hold the contamination estimation result
pub struct ContamProbResult {
    /// for the given contamination level
    pub contamination_level: f64,
    /// what is the log_likelihood given the observed variants?
    pub log_likelihood: f64,
}

#[derive(Serialize, Deserialize, Debug)]
/// data structure for  a variant position
pub struct VariantPosition {
    /// contig name for where the variant is located at
    pub contig: String,
    /// genomic position of the variant
    pub position: usize,
    /// total read depth at the variant position
    pub total_read_depth: usize,
    /// total read that showed alt alleles in the variant position
    pub alt_depth: usize,
    /// is it a indel or snv?
    pub variant_type: VariantType,
    /// the zygosity of the variant
    pub zygosity: Zygosity,
}

impl VariantPosition {
    /// Create a VariantPosition object
    ///
    /// Example::
    ///
    /// ```
    /// let variant = VariantPosition::new(
    ///     "chr1", 1, 100, 50, VariantType::SNV, Zygosity::HETEROZYGOUS
    /// );
    /// ```
    pub fn new(
        contig: &str,
        position: usize,
        total_read_depth: usize,
        alt_depth: usize,
        variant_type: VariantType,
        zygosity: Zygosity,
    ) -> Self {
        if total_read_depth < alt_depth || total_read_depth < 1 {
            // validation of the input
            panic!("Total read depth should be >= alt depth and positive")
        }
        Self {
            contig: contig.to_string(),
            position: position,
            total_read_depth: total_read_depth,
            alt_depth: alt_depth,
            variant_type: variant_type,
            zygosity: zygosity,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Total read depth should be >= alt dept")]
    fn test_variant_position_exception() {
        let _vp = VariantPosition::new(
            "chrX",
            1,
            100,
            101,
            VariantType::SNV,
            Zygosity::HETEROZYGOUS,
        );
    }

    #[test]
    #[should_panic(expected = "and positive")]
    fn test_variant_position_exception_0_depth() {
        let depth: i32 = 0;
        let _vp = VariantPosition::new(
            "chrX",
            1,
            depth as usize,
            101,
            VariantType::SNV,
            Zygosity::HETEROZYGOUS,
        );
    }
}

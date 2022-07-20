/// Defining models for the code
/// 
#[derive(Eq, PartialEq)]
pub enum VariantType {
    /// variant types in the vcf file
    SNV,
    INDEL,
}

#[derive(Eq, PartialEq)]
pub enum Zygosity {
    /// zygostiy of a variant
    HOMOZYGOUS,
    HETEROZYGOUS,
}

pub struct VariantPosition {
    /// data structure for  a variant position
    pub total_read_depth: usize, // total read depth at the variant position
    pub alt_depth: usize, // total read that showed alt alleles in the variant position
    pub variant_type: VariantType, // is it a indel or snv?
    pub zygosity: Zygosity, // the zygosity of the variant
}

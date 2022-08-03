extern crate noodles_bed;

use log::info;
use std::fs::File;
use std::io::BufReader;
use std::string::String;
use std::vec::Vec;

use noodles_bed as bed;

/// Reading bed file and output a list of regions string (contig:start-end)
///
/// # Arguments:
/// * `bed_file`: bed file path pointing to the data to be parsed
///
///
/// # Return:
/// * List of region string
///
/// # Examples
///
/// ```
/// use diploid_contam_estimator::bed::read_bed;
/// let region_list = read_bed("data/test.bed");
/// assert_eq!(region_list.len(), 4);
/// assert_eq!(region_list[0], "X:2-5");
/// ```
pub fn read_bed(bed_file: &str) -> Vec<String> {
    let mut region_list: Vec<String> = vec![];
    let mut reader = File::open(bed_file)
        .map(BufReader::new)
        .map(bed::Reader::new)
        .expect("Error reading bed file");

    for record in reader.records::<3>() {
        let bed_record = record.unwrap();
        let start: usize = usize::try_from(bed_record.start_position()).unwrap() - 1;
        let stop: usize = usize::try_from(bed_record.end_position()).unwrap();
        let contig = bed_record.reference_sequence_name();
        let region_string: String = format!("{}:{}-{}", contig, start, stop);
        region_list.push(region_string);
    }
    info!("Collected {} loci from {}", region_list.len(), bed_file);
    return region_list;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_variant_list() {
        let bed_file = "data/test.bed";
        let region_list = read_bed(bed_file);
        assert_eq!(region_list.len(), 4);
        assert_eq!(region_list[0], "X:2-5");
        assert_eq!(region_list[1], "X:3-7");
        assert_eq!(region_list[3], "1:3-7");
    }
}

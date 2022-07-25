extern crate noodles_bed;

use log::info;
use std::fs::File;
use std::io::BufReader;
use std::option::Option;
use std::string::String;
use std::vec::Vec;

use noodles_bed as bed;

pub fn read_bed(bed_file: Option<&str>) -> Vec<String> {
    let mut region_list: Vec<String> = vec![];
    if bed_file.is_some() {
        let mut reader = File::open(bed_file.unwrap())
            .map(BufReader::new)
            .map(bed::Reader::new)
            .unwrap();

        for record in reader.records::<3>() {
            let bed_record = record.unwrap();
            let start: usize = usize::try_from(bed_record.start_position()).unwrap() - 1;
            let stop: usize = usize::try_from(bed_record.end_position()).unwrap();
            let contig = bed_record.reference_sequence_name();
            let region_string: String = format!("{}:{}-{}", contig, start, stop);
            region_list.push(region_string);
        }
        info!(
            "Collected {} loci from {}",
            region_list.len(),
            bed_file.unwrap()
        );
    }
    return region_list;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_variant_list() {
        let bed_file = "data/test.bed";
        let region_list = read_bed(Some(bed_file));
        assert_eq!(region_list.len(), 3);
        assert_eq!(region_list[0], "X:2-5");
        assert_eq!(region_list[1], "X:3-7");
    }
}

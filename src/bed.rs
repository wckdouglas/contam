extern crate noodles_bed;

use std::{env, fs::File, io::BufReader};

use noodles_bed as bed;

fn main() {
    let mut reader = File::open("./test.bed")
        .map(BufReader::new)
        .map(bed::Reader::new)
        .unwrap();

    for record in reader.records::<3>() {
        let bed_record = record.unwrap();
        let start: usize = usize::try_from(bed_record.start_position()).unwrap() - 1;
        let stop: usize = usize::try_from(bed_record.end_position()).unwrap();
        let contig = bed_record.reference_sequence_name();
        println!("{}:{}-{}", contig, start, stop);
    }
}

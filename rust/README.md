# Rust #

Here's a rust version of the code.

To run the code:

```{bash}
$ cargo run -i input_vcf -d debug_json
```

or:

```
$ diploid-contam-estimator  -h
diploid-contam-estimator 0.1.0
Douglas Wu <wckdouglas@gmail.com>
Estimating contamination level from a diploid VCF file


    The program assume we are dealing with a diploid genome, and using the
    deviation of allelic balance from the expected allelic frequence for homozygous
    or heterozygous variant calls to compute a contamination value.

    For homozygous variants, we deviation from allelic frequency of 1 is all introduced by
contaminaion.

    For heterozygous variants, it is a little more complex, because it could be due to:
        1. contamination that doesn't look like the HET ALT allele: we expect lower HET alt allele
frequency
        2. contamination that doesn't look like the HOM ALT allele: we expect High HET alt allele
frequency
        3. contamination that looks like the ALT allele: we expect higher alt allele frequency
        4. contamination that looks like the REF allele: we expect lower alt allele frequency
        5. contamination being called as ALT

USAGE:
    diploid-contam-estimator [OPTIONS] --in-vcf <in_vcf>

OPTIONS:
    -d, --debug-json <debug_json>    A json output file for storing all intermediate log prob
    -h, --help                       Print help information
    -i, --in-vcf <in_vcf>            A diploid vcf file for estimating contamination
    -V, --version                    Print version information
```
# contam-estimator #

[![poetry CI](https://github.com/wckdouglas/contam/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wckdouglas/contam/actions/workflows/ci.yml)


A module to estimate contamination level from diploid variant calls. This is heavily inspired by [Dcon](https://github.com/liguowang/dcon/blob/master/lib/DconModule/utils.py).

Basically, we are hypothesizing the contamination level would be between $0-0.4$, and at each hypothesized contamination level, we'll calculate the likelihood of the observed variant-allele-frequency assuming the given contamination level is true with the for each variant in the VCF file using a binomial model.

For a homozygous variant, the expected varaint-allele-frequency is:

$$ P(\text{expected vaf}) = 1 - cl  $$ 

where $ cl \in [0,0.4] $ is the contamination level.


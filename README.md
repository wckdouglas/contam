# contam-estimator #

[![poetry CI](https://github.com/wckdouglas/contam/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wckdouglas/contam/actions/workflows/ci.yml)


A module to estimate contamination level from diploid variant calls. This is heavily inspired by [Dcon](https://github.com/liguowang/dcon/blob/master/lib/DconModule/utils.py).

Basically, we are hypothesizing the contamination level would be between $0-0.4$, and at each hypothesized contamination level, we'll calculate the likelihood observing the observed variant-allele-frequency assuming the given contamination level is true for each variant in the VCF file using a binomial model. We then sum the log likelihood for all variants and pick the maximum likelihood contamination level as the final call.

For a homozygous variant, the expected variant-allele-frequency is:

$$ P(\text{expected vaf}) = 1 - cl  $$ 

where $ cl \in [0,0.4] $ is the contamination level.

For a heterozygous variant, the expected variant-allele-frequency can either be:

1. low ALT allele frequency because of the contamination:  $ (1 - cl) / 2 $
2. a HOM being called as HET because of the contamination: $ (1 - cl) $
3. the contamination looks like ALT allele, such that the ALT allele frequency is higher than expected: $ (0.5 + cl) $
4. the contam is being called as HET variant: $ cl $

We will evaluate all these cases and pick the highest probability event when summing the likelihoods for the given contamination level.
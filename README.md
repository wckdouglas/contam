# contam-estimator #

[![poetry CI](https://github.com/wckdouglas/contam/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wckdouglas/contam/actions/workflows/ci.yml)


A module to estimate contamination level from **diploid** variant calls. This is heavily inspired by [Dcon](https://github.com/liguowang/dcon/blob/master/lib/DconModule/utils.py).

# Background #

Basically, we are hypothesizing the contamination level would be between $0-0.4$, and at each hypothesized contamination level, we'll calculate the likelihood observing the observed variant-allele-frequency assuming the given contamination level is true for each variant in the VCF file using a binomial model. We then sum the log likelihood for all variants and pick the maximum likelihood contamination level as the final call.

Pseudo code:
```
n = total_count
x = alt_allele_count
contam = 0
max_log_likelihood = -inf
for contam_level in all_contamination_level:
    p = expected_alt_fraction_for_the_given_contamination_level
    log_likelihood = sum(binom_loglik(n, x, p) for all_variants)
    if log_likelihood > max_log_likelihood:
        contam = contam_level
```

A simulated study at [here](https://github.com/wckdouglas/contam/blob/main/notebooks/contam_simulator.ipynb).


### Homozygous variants

For a homozygous variant, the probability of observing the expected variant-allele-count ($x$) with a read depth $n$ at a given contamination level $ c \in [0,0.4] $ is :

$$ P(X=x) = \binom{n}{x}p^x(1 - p)^{n-x}  $$ 

where $ p = (1-c) $ in all homozygous variants

![](https://github.com/wckdouglas/contam/blob/main/img/hom.png?raw=true)

### Heterozygous variants

For a heterozygous variant, the probablity of observing the expected variant-allele-count ($x$) with a read depth $n$ at a given contamination level $ c \in [0,0.4] $ will be following the above binomial distribution but $p$ can either be:


1. $ (1 - c) / 2 $, when low alternate allele frequency because of the contamination
2. $ (1 - c) $, when a homozygous variant being called as heterozygous because of the contamination
3. $ (0.5 + c) $, when the contamination looks like alternate allele, such that the alternate allele frequency is higher than expected
4. $ c $, when the contam is being called as HET variant

We will evaluate all these cases and pick the highest probability event when summing the log likelihoods for the given contamination level.

![](https://github.com/wckdouglas/contam/blob/main/img/het.png?raw=true)

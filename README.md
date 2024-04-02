
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProbBreed

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/ProbBreed)](https://CRAN.R-project.org/package=ProbBreed)
[![ProbBreed status
badge](https://saulo-chaves.r-universe.dev/badges/ProbBreed)](https://saulo-chaves.r-universe.dev/ProbBreed)
<!-- badges: end -->

ProbBreed employs Bayesian statistics to analyse multi-environment
trials’ data, and uses its outputs to estimate the marginal and pairwise
probabilities of superior performance and superior stability of the
genotypes, as well as their conditional probability of superior
performance. The method is thoroughly described at
<https://doi.org/10.1007/s00122-022-04041-y>.

## Installation

You can install the CRAN version of `ProbBreed` using the following
command:

``` r
install.packages("ProbBreed")
```

Alternatively, you can install the development version of `ProbBreed`
from [GitHub](https://github.com/saulo-chaves/ProbBreed) with:

``` r
# install.packages("devtools")
devtools::install_github("saulo-chaves/ProbBreed")
```

## Usage

A basic workflow using the available data is:

``` r
library(ProbBreed)

mod = bayes_met(data = maize,
                gen = "Hybrid",
                loc = "Location",
                repl = c("Rep", "Block"),
                year = NULL,
                reg = "Region",
                res.het = F,
                trait = "GY",
                iter = 6000, cores = 4, chains = 4)

outs = extr_outs(data = maize, trait = "GY", model = mod,
                 probs = c(0.05, 0.95),
                 check.stan.diag = FALSE, 
                 verbose = TRUE)

results = prob_sup(data = maize, trait = "GY", gen = "Hybrid", loc = "Location",
                   mod.output = outs, reg = 'Region', year = NULL, int = .2,
                   increase = TRUE, save.df = FALSE, interactive = FALSE, 
                   verbose = TRUE)
```

## Citation

For citing the package, use:

``` r
citation('ProbBreed')
```

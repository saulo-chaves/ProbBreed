
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

If the levels of the factors (Genotype, Environment and Region) in your
dataset are not alphanumeric, you can use the `recod` function to recode
it:

``` r
toy = data.frame(
  gen = rep(1:5, each =4),
  env = rep(1:4, times = 5),
  y = rnorm(20, mean = 20, sd = 5)
)
new.toy = recod(toy, name.fact = 'gen', cod = 'G')
```

A basic workflow using the available data is:

``` r
library(ProbBreed)

mod = bayes_met(data = soy,
                gen = "Gen",
                env = "Env",
                repl = NULL,
                reg = "Reg",
                res.het = F,
                trait = "Y",
                iter = 2000, cores = 1, chains = 4)

outs = extr_outs(data = soy, trait = "Y", gen = "Gen", model = mod,
                 effects = c('l','g','gl','m','gm'),
                 nenv = length(unique(soy$Env)),
                 probs = c(0.05, 0.95), check.stan.diag = FALSE, 
                 verbose = TRUE)

results = prob_sup(data = soy, trait = "Y", gen = "Gen", env = "Env",
                   mod.output = outs, reg = 'Reg', int = .2,
                   increase = TRUE, save.df = FALSE, interactive = FALSE, 
                   verbose = TRUE)
```

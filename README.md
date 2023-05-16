
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProbBreed

<!-- badges: start -->
<!-- badges: end -->

ProbBreed employs Bayesian statistics to analyse multi-environment
trials’ data, and uses its outputs to estimate the marginal, pairwise,
and conditional probability of superior performance of the genotypes.
The method is thoroughly described at
<https://doi.org/10.1007/s00122-022-04041-y>.

## Installation

You can install the development version of ProbBreed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("saulo-chaves/ProbBreed")
```

## Usage

A basic workflow using the available data is:

``` r
library(ProbBreed)

mod = bayes_met(data = soy,
                gen = "Gen",
                env = "Env", 
                rept = NULL,
                reg = "Reg",
                res.het = F,
                trait = "Y", 
                iter = 2000, cores = 1, chain = 1)

outs = extr_outs(data = soy, 
                 trait = "Y", 
                 gen = "Gen", 
                 model = mod, 
                 effects = c('l','g','gl','m','gm'),
                 nenv = length(unique(soy$Env)), 
                 res.het = FALSE,
                 probs = c(0.05, 0.95)
                 check.stan.diag = TRUE)

margs = marg_prob(data = soy, 
                  trait = "Y",
                  gen = "Gen", 
                  env = "Env",
                  extr_outs = outs, 
                  int = 0.2,
                  increase = TRUE, 
                  save.df = FALSE, 
                  interactive = FALSE)

conds = cond_prob(data = soy, 
                  trait = "Y",
                  gen = "Gen", 
                  env = "Env",
                  reg = "Reg",
                  extr_outs = outs, 
                  int = 0.2,
                  increase = TRUE, 
                  save.df = FALSE, 
                  interactive = FALSE)
```

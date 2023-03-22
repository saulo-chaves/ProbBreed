
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProbBreed

<!-- badges: start -->
<!-- badges: end -->

ProbBreed employs Bayesian statistics to analyse multi-environment
trials’ data, and uses its outputs to estimate the marginal, pairwise,
and conditional probability of superior performance of the genotypes.
The method is thoroughly described at
<https://doi.org/10.1007/s00122-022-04041-y>. This is a work in
progress.

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
                gen = c("Gen", "normal", "cauchy"),
                env = c("Env", "normal", "cauchy"),
                rept = NULL,
                reg = list(c("Reg", "normal", "cauchy"), c("normal", "cauchy")),
                res.het = F,
                sigma.dist = c("cauchy", "cauchy"),
                mu.dist = c("normal", "cauchy"),
                gei.dist = c("normal", "normal"),
                trait = "eBLUE", 
                hyperparam = "default",
                iter = 100, cores = 4, chain = 4)
# You may want to increase the number of iterations, cores and chains

outs = extr_outs(data = soy, 
                 trait = "eBLUE", 
                 gen = "Gen", 
                 model = mod, 
                 effects = c('l','g','gl','m','gm'),
                 nenv = length(unique(soy$Env)), 
                 res.het = FALSE,
                 probs = c(0.05, 0.95)
                 check.stan.diag = TRUE)

margs = marg_prob(data = soy, 
                  trait = "eBLUE",
                  gen = "Gen", 
                  env = "Env",
                  extr_outs = outs, 
                  int = 0.2,
                  increase = TRUE, 
                  save.df = FALSE, 
                  interactive = FALSE)

conds = cond_prob(data = soy, 
                  trait = "eBLUE",
                  gen = "Gen", 
                  env = "Env",
                  reg = "Reg",
                  extr_outs = outs, 
                  int = 0.2,
                  increase = TRUE, 
                  save.df = FALSE, 
                  interactive = FALSE)
```

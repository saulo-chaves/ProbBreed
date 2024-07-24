
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProbBreed

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/ProbBreed)](https://CRAN.R-project.org/package=ProbBreed)
[![ProbBreed status
badge](https://saulo-chaves.r-universe.dev/badges/ProbBreed)](https://saulo-chaves.r-universe.dev/ProbBreed)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/ProbBreed)](https://cran.r-project.org/package=ProbBreed)
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

``` r
library(ProbBreed)
```

Currently, `ProbBreed` has nine available models implemented in the
`bayes_met` function. See `?bayes_met` for more details. An examples
using the `soy` example dataset is described below:

``` r
model = bayes_met(data = soy,
                gen = "Gen",
                loc = "Loc",
                repl = NULL,
                year = NULL,
                reg = NULL,
                res.het = TRUE,
                trait = 'Y',
                iter = 6000, cores = 4, chains = 4)
```

`gen`, `loc`, `repl`, `year` and `reg` are all column names that contain
information on genotypes, environments (locations), replicates, years
(or seasons) and regions (or mega-environments). The `soy` dataset has
only adjusted means, so only `gen` and `loc` have values. `res.het`
indicates if a per-environmental residual variance should be estimated.
`trait` is the column in `data` that contain the phenotypic
observations. The other arguments are specifications for model fitting:
the number of iterations, cores and chains. Feel free to customize these
and other options according to your necessity.

The output of this function will be an object of class `stanfit`, which
should be used in the `extr_outs` function for further processing before
computing the probabilities per se. This function also provides some
useful diagnostics. Here is how to use it:

``` r
outs = extr_outs(model = model,
                 probs = c(0.05, 0.95),
                 verbose = TRUE)
```

The object of class `extr` provided by this function contains the
effects’ posterior and maximum posterior, the models’ variance
components and some posterior predictive checks. Here are them:

``` r
outs$variances
#>         effect     var      sd naive.se HPD_0.05 HPD_0.95
#> 1          Gen   3.381   1.295    0.012    1.622    5.794
#> 2          Loc 247.136 131.148    1.197  116.979  474.143
#> 3   error_env1  10.060   2.716    0.025    6.350   15.090
#> 4   error_env2  28.760   7.162    0.065   19.053   41.710
#> 5   error_env3  11.773   3.436    0.031    7.307   18.012
#> 6   error_env4  18.556   4.928    0.045   12.003   27.531
#> 7   error_env5  51.396  12.677    0.116   34.404   75.112
#> 8   error_env6  15.398   3.986    0.036   10.018   22.784
#> 9   error_env7  19.359   4.845    0.044   12.683   28.259
#> 10  error_env8  21.205   5.438    0.050   13.893   31.346
#> 11  error_env9  14.971   3.853    0.035    9.719   22.103
#> 12 error_env10  13.426   5.648    0.052    6.945   23.900
#> 13 error_env11  22.025   9.007    0.082   11.667   38.359
#> 14 error_env12   7.681   3.405    0.031    3.692   14.143
#> 15 error_env13  14.728   6.150    0.056    7.566   26.174
#> 16 error_env14  10.867   4.490    0.041    5.600   19.145
outs$ppcheck
#>                   Diagnostics
#> p.val_max              0.9190
#> p.val_min              0.3497
#> p.val_median           0.7192
#> p.val_mean             0.5058
#> p.val_sd               0.5553
#> Eff_No_parameters     27.5098
#> WAIC2               2715.7241
#> mean_Rhat              1.0004
#> Eff_sample_size        0.6785
```

You can also the `plot` S3 method for some useful visualizations. For
e.g., the comparison between the empirical and sampled phenotype
illustrates the model’s convergence:

``` r
plot(outs)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

See `?plot.extr` for more details and further options.

After these two steps, everything is set to compute the probabilities.
This can be done using the function `prob_sup`:

A basic workflow using the available data is:

``` r
results = prob_sup(extr = outs, 
                   int = .2,
                   increase = TRUE, 
                   save.df = FALSE, 
                   verbose = TRUE)
```

This function generates an object of class `probsup`, which contains two
lists: `across` and `within`. As their names suggest, the `across` list
has the across-environments probabilities, and is suitable for a broader
recommendation. Conversely, the `within` results are more appropriate to
specific recommendations. For example, here are some probability of
superior performances across and within environments:

``` r
head(results$across$perfo)
#>     ID      prob
#> 36 G36 0.9855000
#> 9  G09 0.8345833
#> 20 G20 0.8215833
#> 38 G38 0.7213333
#> 31 G31 0.6745000
#> 1  G01 0.5281667
head(results$within$perfo$gl)
#>   gen         E01         E02        E03         E04         E05         E06
#> 1 G01 0.528166667 0.528166667 0.66775000 0.528166667 0.528166667 0.528166667
#> 2 G02 0.006250000 0.006250000 0.01300000 0.006250000 0.006250000 0.006250000
#> 3 G03 0.162083333 0.162083333 0.24758333 0.162083333 0.162083333 0.162083333
#> 4 G04 0.014166667 0.014166667 0.02766667 0.014166667 0.014166667 0.014166667
#> 5 G05 0.001333333 0.001333333 0.00275000 0.001333333 0.001333333 0.001333333
#> 6 G06 0.199083333 0.199083333 0.28308333 0.199083333 0.199083333 0.199083333
#>           E07         E08         E09       E10       E11       E12       E13
#> 1 0.528166667 0.528166667 0.528166667 0.9563333 0.9563333 0.9563333 0.9563333
#> 2 0.006250000 0.006250000 0.006250000        NA        NA        NA        NA
#> 3 0.162083333 0.162083333 0.162083333        NA        NA        NA        NA
#> 4 0.014166667 0.014166667 0.014166667        NA        NA        NA        NA
#> 5 0.001333333 0.001333333 0.001333333        NA        NA        NA        NA
#> 6 0.199083333 0.199083333 0.199083333        NA        NA        NA        NA
#>         E14
#> 1 0.9563333
#> 2        NA
#> 3        NA
#> 4        NA
#> 5        NA
#> 6        NA
```

The S3 method `plot` is also available for `probsup` objects. Here are
some of them:

- Probability of superior performance across environments

``` r
plot(results)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

- Pairwise probability of superior performance across environments

``` r
plot(results, category = "pair_perfo", level = "across")
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

- Probability of superior stability

``` r
plot(results, category = "stabi")
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

- Probability of superior performance within environments

``` r
plot(results, category = "perfo", level = "within")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

The grey cells are environments where a given genotype was not assessed.
See more options at `?plot.probsup`.

The estimation of these probabilities are strictly related to some key
questions that constantly arises in plant breeding, like:

- **What is the risk of recommending a selection candidate for a target
  population of environments?**

- **What is the probability of a given selection candidate having good
  performance if recommended to a target population of environments? And
  for a specific environment?**

- **What is the probability of a given selection candidate having better
  performance than a cultivar check in the target population of
  environments? And in specific environments?**

- **How probable is it that a given selection candidate performs
  similarly across environments?**

- **What are the chances that a given selection candidate is more stable
  than a cultivar check in the target population of environments?**

- **What is the probability that a given selection candidate having a
  superior and invariable performance across environments?**

For a more detailed tutorial, see
<https://saulo-chaves.github.io/ProbBreed_site/>.

## Citation

Thank you for using `ProbBreed`! Please, do not forget to cite:

``` r
citation('ProbBreed')
```

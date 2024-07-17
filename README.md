
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
`bayes_met` function. See `?bayes_met` for more details. Using the `soy`
example dataset, a basic usage is as follows:

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
#> 1          Gen   3.383   1.283    0.012    1.655    5.758
#> 2          Loc 246.354 124.591    1.137  115.921  470.127
#> 3   error_env1  10.077   2.679    0.024    6.418   15.058
#> 4   error_env2  28.924   7.296    0.067   19.029   42.483
#> 5   error_env3  11.620   3.372    0.031    7.210   17.970
#> 6   error_env4  18.458   4.732    0.043   12.092   27.309
#> 7   error_env5  51.193  12.312    0.112   34.576   73.887
#> 8   error_env6  15.524   4.085    0.037   10.051   23.073
#> 9   error_env7  19.385   5.003    0.046   12.723   28.555
#> 10  error_env8  21.217   5.390    0.049   13.971   31.041
#> 11  error_env9  14.818   3.862    0.035    9.658   21.890
#> 12 error_env10  13.633   5.612    0.051    7.024   24.165
#> 13 error_env11  22.330   9.116    0.083   11.656   39.590
#> 14 error_env12   7.657   3.377    0.031    3.693   14.020
#> 15 error_env13  14.604   5.965    0.054    7.650   25.520
#> 16 error_env14  10.983   4.557    0.042    5.678   19.356
outs$ppcheck
#>                   Diagnostics
#> p.val_max              0.9218
#> p.val_min              0.3498
#> p.val_median           0.7276
#> p.val_mean             0.5012
#> p.val_sd               0.5630
#> Eff_No_parameters     27.4466
#> WAIC2               2715.4466
#> mean_Rhat              1.0002
#> Eff_sample_size        0.8112
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
#> 36 G36 0.9845833
#> 9  G09 0.8344167
#> 20 G20 0.8135833
#> 38 G38 0.7305833
#> 31 G31 0.6705000
#> 1  G01 0.5160000
head(results$within$perfo$gl)
#>   gen          E01          E02         E03          E04          E05
#> 1 G01 0.5160000000 0.5160000000 0.656916667 0.5160000000 0.5160000000
#> 2 G02 0.0062500000 0.0062500000 0.013916667 0.0062500000 0.0062500000
#> 3 G03 0.1638333333 0.1638333333 0.247416667 0.1638333333 0.1638333333
#> 4 G04 0.0138333333 0.0138333333 0.029250000 0.0138333333 0.0138333333
#> 5 G05 0.0008333333 0.0008333333 0.002583333 0.0008333333 0.0008333333
#> 6 G06 0.1968333333 0.1968333333 0.290083333 0.1968333333 0.1968333333
#>            E06          E07          E08          E09       E10       E11
#> 1 0.5160000000 0.5160000000 0.5160000000 0.5160000000 0.9541667 0.9541667
#> 2 0.0062500000 0.0062500000 0.0062500000 0.0062500000        NA        NA
#> 3 0.1638333333 0.1638333333 0.1638333333 0.1638333333        NA        NA
#> 4 0.0138333333 0.0138333333 0.0138333333 0.0138333333        NA        NA
#> 5 0.0008333333 0.0008333333 0.0008333333 0.0008333333        NA        NA
#> 6 0.1968333333 0.1968333333 0.1968333333 0.1968333333        NA        NA
#>         E12       E13       E14
#> 1 0.9541667 0.9541667 0.9541667
#> 2        NA        NA        NA
#> 3        NA        NA        NA
#> 4        NA        NA        NA
#> 5        NA        NA        NA
#> 6        NA        NA        NA
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

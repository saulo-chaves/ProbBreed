# Bayesian Probabilistic Selection Index (BPSI)

This function estimates the genotype's merit for multiple traits using
the probabilities of superior performance across environments.

## Usage

``` r
bpsi(problist, increase = NULL, lambda = NULL, int, save.df = FALSE)
```

## Arguments

- problist:

  A list of object of class `probsup`, obtained from the
  [prob_sup](https://saulo-chaves.github.io/ProbBreed/reference/prob_sup.md)
  function

- increase:

  Optional logical vector with size corresponding to the number of
  traits of `problist`, in the same order.`TRUE` if the selection is for
  increasing the trait value, `FALSE` otherwise. If not declared, `bpsi`
  will consider the information provided in
  [prob_sup](https://saulo-chaves.github.io/ProbBreed/reference/prob_sup.md)

- lambda:

  A numeric representing the weight of each trait. Defaults to 1 (equal
  weights). The trait with more economic interest should be greater.

- int:

  A numeric representing the selection intensity (between 0 and 1),
  considering the selection index.

- save.df:

  Logical. Should the data frames be saved in the work directory? `TRUE`
  for saving, `FALSE` (default) otherwise.

## Value

The function returns an object of class `bpsi`, which contains two
lists, one with the BPSI- Bayesian Probabilistic Selection Index, and
another with the original `data`- with across-environments probabilities
of superior performance for each trait.

## Details

- Bayesian Probabilistic Selection Index

\$\$BPSI_i = \sum\_{m=1}^{t} \frac{\gamma\_{pt} -\gamma\_{it}
}{(1/\lambda_t)}\$\$

where \\\gamma_p\\ is the probability of superior performance of the
worst genotype for the trait \\t\\, \\\gamma\\ is the probability of
superior performance of genotype \\i\\ for trait \\t\\, \\t\\ is the
total number of traits evaluated, \\\left(m = 1, 2, ..., t \right)\\,
and \\\lambda\\ is the weight for each trait \\t\\.

More details about the usage of `bpsi` can be found at
<https://tiagobchagas.github.io/BPSI/>.

## References

Chagas, J. T. B., Dias, K. O. G., Carneiro, V. Q., Oliveira, L. M. C.,
Nunes, N. X., Pereira Júnior, J. D., Carneiro, P. C. S., & Carneiro, J.
E. S. (2025). Bayesian probabilistic selection index in the selection of
common bean families. *Crop Science*, 65(3).
[doi:10.1002/CSC2.70072](https://doi.org/10.1002/CSC2.70072)

## See also

[`plot.bpsi()`](https://saulo-chaves.github.io/ProbBreed/reference/plot.bpsi.md)

## Author

José Tiago Barroso Chagas

## Examples

``` r
# \donttest{


met_df <-
read.table("https://raw.githubusercontent.com/tiagobchagas/BPSI/refs/heads/main/Data/blues_long.txt",header = TRUE)

mod = bayes_met(data = met_df,
                gen = "gen",
                loc = "env",
                repl = NULL,
                trait = "PH",
                reg = NULL,
                year = NULL,
                res.het = TRUE,
                iter = 2000, cores = 2, chain = 4)
#> in 2:                21.52 seconds (Total)
#> Chain 2: 
#> Warning: There were 24 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess


mod2 = bayes_met(data = met_df,
                 gen = "gen",
                 loc = "env",
                 repl = NULL,
                 trait = "GY",
                 reg = NULL,
                 year = NULL,
                 res.het = TRUE,
                 iter = 2000, cores = 2, chain = 4)
#> :                21.52 seconds (Total)
#> Chain 2: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess

mod3 = bayes_met(data = met_df,
                 gen = "gen",
                 loc = "env",
                 repl =  NULL,
                 trait = "NDM",
                 reg = NULL,
                 year = NULL,
                 res.het = TRUE,
                 iter = 2000, cores = 2, chain = 4)
#> hain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Warning: There were 41 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess



models=list(mod,mod2,mod3)
names(models) <- c("PH","GY","NDM")
increase = c(FALSE,TRUE,FALSE)
names(increase) <- names(models)

probs = list()
for (i in names(models)) {
  outs <- extr_outs(model = models[[i]],
                    probs = c(0.05, 0.95),
                    verbose = TRUE)
  probs[[i]] <- prob_sup(
    extr = outs,
    int = .2,
    increase = increase[[i]],
    save.df = FALSE,
    verbose = TRUE
  )

}
#> -> Posterior effects extracted
#> -> Variances extracted
#> -> Maximum posterior values extracted
#> -> Posterior predictive checks computed
#> 24 of 4000 iterations ended with a divergence (0.6%).
#> Try increasing 'adapt_delta' to remove the divergences.
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
#> E-BFMI indicated no pathological behavior.
#> -> Probability of superior performance estimated
#> -> Pairwise probability of superior performance estimated
#> Process completed!
#> -> Posterior effects extracted
#> -> Variances extracted
#> -> Maximum posterior values extracted
#> -> Posterior predictive checks computed
#> 0 of 4000 iterations ended with a divergence.
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
#> E-BFMI indicated no pathological behavior.
#> -> Probability of superior performance estimated
#> -> Pairwise probability of superior performance estimated
#> Process completed!
#> -> Posterior effects extracted
#> -> Variances extracted
#> -> Maximum posterior values extracted
#> -> Posterior predictive checks computed
#> 41 of 4000 iterations ended with a divergence (1.025%).
#> Try increasing 'adapt_delta' to remove the divergences.
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
#> E-BFMI indicated no pathological behavior.
#> -> Probability of superior performance estimated
#> -> Pairwise probability of superior performance estimated
#> Process completed!

index = bpsi(
  problist = probs,
  increase = increase,
  int = 0.1,
  lambda = c(1, 2, 1),
  save.df = FALSE
)
# }
```

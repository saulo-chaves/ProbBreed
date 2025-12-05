# Probabilities of superior performance and stability

This function estimates the probabilities of superior performance and
stability across environments, and probabilities of superior performance
within environments.

## Usage

``` r
prob_sup(extr, int, increase = TRUE, save.df = FALSE, verbose = FALSE)
```

## Arguments

- extr:

  An object of class `extr`, obtained from the
  [extr_outs](https://saulo-chaves.github.io/ProbBreed/reference/extr_outs.md)
  function

- int:

  A numeric representing the selection intensity (between 0 and 1)

- increase:

  Logical.`TRUE` (default) if the selection is for increasing the trait
  value, `FALSE` otherwise.

- save.df:

  Logical. Should the data frames be saved in the work directory? `TRUE`
  for saving, `FALSE` (default) otherwise.

- verbose:

  A logical value. If `TRUE`, the function will indicate the completed
  steps. Defaults to `FALSE`.

## Value

The function returns an object of class `probsup`, which contains two
lists, one with the `across`-environments probabilities, and another
with the `within`-environments probabilities. If an entry-mean model was
used in
[`ProbBreed::bayes_met`](https://saulo-chaves.github.io/ProbBreed/reference/bayes_met.md),
only the `across` list will be available.

The `across` list has the following elements:

- `g_hpd`: Highest posterior density (HPD) of the posterior genotypic
  main effects.

- `perfo`: the probabilities of superior performance.

- `pair_perfo`: the pairwise probabilities of superior performance.

- `stabi`: a list with the probabilities of superior stability. It
  contains the data frames `gl`, `gm` (when `reg` is not `NULL`) and
  `gt` (when `year` is not `NULL`). Unavailable if an entry-mean model
  was used in `bayes_met`.

- `pair_stabi`: a list with the pairwise probabilities of superior
  stability. It contains the data frames `gl`, `gm` (when `reg` is not
  `NULL`) and `gt` (when `year` is not `NULL`). Unavailable if an
  entry-mean model was used in `bayes_met`.

- `joint_prob`: the joint probabilities of superior performance and
  stability. Unavailable if an entry-mean model was used in `bayes_met`.

The `within` list has the following elements:

- `perfo`: a list of data frames containing the probabilities of
  superior performance within locations (`gl`), regions (`gm`) and years
  (`gt`).

- `pair_perfo`: lists with the pairwise probabilities of superior
  performance within locations (`gl`), regions (`gm`) and years (`gt`).

## Details

Probabilities provide the risk of recommending a selection candidate for
a target population of environments or for a specific environment.
`prob_sup` computes the probabilities of superior performance and the
probabilities of superior stability:

- Probability of superior performance

Let \\\Omega\\ represent the subset of selected genotypes based on their
performance across environments. A given genotype \\j\\ will belong to
\\\Omega\\ if its genotypic marginal value (\\\hat{g}\_j\\) is high or
low enough compared to its peers. `prob_sup` leverages the Monte Carlo
discretized sampling from the posterior distribution to emulate the
occurrence of \\S\\ trials. Then, the probability of the \\j^{th}\\
genotype belonging to \\\Omega\\ is the ratio of success (\\\hat{g}\_j
\in \Omega\\) events and the total number of sampled events, as follows:

\$\$Pr\left(\hat{g}\_j \in \Omega \vert y \right) =
\frac{1}{S}\sum\_{s=1}^S{I\left(\hat{g}\_j^{(s)} \in \Omega \vert
y\right)}\$\$

where \\S\\ is the total number of samples \\\left(s = 1, 2, ..., S
\right)\\, and \\I\left(g_j^{(s)} \in \Omega \vert y\right)\\ is an
indicator variable that can assume two values: (1) if \\\hat{g}\_j^{(s)}
\in \Omega\\ in the \\s^{th}\\ sample, and (0) otherwise. \\S\\ is
conditioned to the number of iterations and chains previously set at
[bayes_met](https://saulo-chaves.github.io/ProbBreed/reference/bayes_met.md).

Similarly, the within-environment probability of superior performance
can be applied to individual environments. Let \\\Omega_k\\ represent
the subset of superior genotypes in the \\k^{th}\\ environment, so that
the probability of the \\j^{th} \in \Omega_k\\ can calculated as
follows:

\$\$Pr\left(\hat{g}\_{jk} \in \Omega_k \vert y\right) = \frac{1}{S}
\sum\_{s=1}^S I\left(\hat{g}\_{jk}^{(s)} \in \Omega_k \vert y\right)\$\$

where \\I\left(\hat{g}\_{jk}^{(s)} \in \Omega_k \vert y\right)\\ is an
indicator variable mapping success (1) if \\\hat{g}\_{jk}^{(s)}\\ exists
in \\\Omega_k\\, and failure (0) otherwise, and \\\hat{g}\_{jk}^{(s)} =
\hat{g}\_j^{(s)} + \widehat{ge}\_{jk}^{(s)}\\. Note that when computing
within-environment probabilities, we are accounting for the interaction
of the \\j^{th}\\ genotype with the \\k^{th}\\ environment.

The pairwise probabilities of superior performance can also be
calculated across or within environments. This metric assesses the
probability of the \\j^{th}\\ genotype being superior to another
experimental genotype or a commercial check. The calculations are as
follows, across and within environments, respectively:

\$\$Pr\left(\hat{g}\_{j} \> \hat{g}\_{j^\prime} \vert y\right) =
\frac{1}{S} \sum\_{s=1}^S I\left(\hat{g}\_{j}^{(s)} \>
\hat{g}\_{j^\prime}^{(s)} \vert y\right)\$\$

or

\$\$Pr\left(\hat{g}\_{jk} \> \hat{g}\_{j^\prime k} \vert y\right) =
\frac{1}{S} \sum\_{s=1}^S I\left(\hat{g}\_{jk}^{(s)} \>
\hat{g}\_{j^\prime k}^{(s)} \vert y\right)\$\$

These equations are set for when the selection direction is positive. If
`increase = FALSE`, \\\>\\ is simply switched by \\\<\\.

- Probability of superior stability

This probability makes a direct analogy with the method of Shukla
(1972): a stable genotype is the one that has a low
genotype-by-environment interaction variance \\\[var(\widehat{ge})\]\\.
Using the same probability principles previously described, the
probability of superior stability is given as follows:

\$\$Pr \left\[var \left(\widehat{ge}\_{jk}\right) \in \Omega \vert y
\right\] = \frac{1}{S} \sum\_{s=1}^S I\left\[var
\left(\widehat{ge}\_{jk}^{(s)} \right) \in \Omega \vert y \right\]\$\$

where \\I\left\[var \left(\widehat{ge}\_{jk}^{(s)} \right) \in \Omega
\vert y \right\]\\ indicates if
\\var\left(\widehat{ge}\_{jk}^{(s)}\right)\\ exists in \\\Omega\\ (1) or
not (0). Pairwise probabilities of superior stability are also possible
in this context:

\$\$Pr \left\[var \left(\widehat{ge}\_{jk} \right) \<
var\left(\widehat{ge}\_{j^\prime k} \right) \vert y \right\] =
\frac{1}{S} \sum\_{s=1}^S I \left\[var \left(\widehat{ge}\_{jk}
\right)^{(s)} \< var \left(\widehat{ge}\_{j^\prime k} \right)^{(s)}
\vert y \right\]\$\$

Note that \\j\\ will be superior to \\j^\prime\\ if it has a **lower**
variance of the genotype-by-environment interaction effect. This is true
regardless if `increase` is set to `TRUE` or `FALSE`.

The joint probability independent events is the product of the
individual probabilities. The estimated genotypic main effects and the
variances of GEI effects are independent by design, thus the joint
probability of superior performance and stability as follows:

\$\$Pr \left\[\hat{g}\_j \in \Omega \cap var \left(\widehat{ge}\_{jk}
\right) \in \Omega \right\] = Pr \left(\hat{g}\_j \in \Omega \right)
\times Pr \left\[var \left(\widehat{ge}\_{jk} \right) \in \Omega
\right\]\$\$

The estimation of these probabilities are strictly related to some key
questions that constantly arises in plant breeding:

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

More details about the usage of `prob_sup`, as well as the other
function of the `ProbBreed` package can be found at
<https://saulo-chaves.github.io/ProbBreed_site/>.

## References

Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães,
L. J. M., Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging
probability concepts for cultivar recommendation in multi-environment
trials. *Theoretical and Applied Genetics*, 133(2):443-455.
[doi:10.1007/s00122-022-04041-y](https://doi.org/10.1007/s00122-022-04041-y)

Shukla, G. K. (1972) Some statistical aspects of partioning genotype
environmental componentes of variability. *Heredity*, 29:237-245.
[doi:10.1038/hdy.1972.87](https://doi.org/10.1038/hdy.1972.87)

## See also

[plot.probsup](https://saulo-chaves.github.io/ProbBreed/reference/plot.probsup.md)

## Examples

``` r
# \donttest{
mod = bayes_met(data = maize,
                gen = "Hybrid",
                loc = "Location",
                repl = c("Rep","Block"),
                trait = "GY",
                reg = "Region",
                year = NULL,
                res.het = TRUE,
                iter = 2000, cores = 2, chain = 4)
#>              557.127 seconds (Total)
#> Chain 1: 
#> Warning: There were 4 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.09, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

outs = extr_outs(model = mod,
                 probs = c(0.05, 0.95),
                 verbose = TRUE)
#> -> Posterior effects extracted
#> -> Variances extracted
#> -> Maximum posterior values extracted
#> -> Posterior predictive checks computed
#> 4 of 4000 iterations ended with a divergence (0.1%).
#> Try increasing 'adapt_delta' to remove the divergences.
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
#> E-BFMI indicated possible pathological behavior:
#>   Chain 1: E-BFMI = 0.091
#> E-BFMI below 0.2 indicates you may need to reparameterize your model.

results = prob_sup(extr = outs,
                   int = .2,
                   increase = TRUE,
                   save.df = FALSE,
                   verbose = FALSE)
# }
```

# Bayesian model for multi-environment trials

Fits a Bayesian multi-environment model using `rstan`, the `R` interface
to `Stan`.

## Usage

``` r
bayes_met(
  data,
  gen,
  loc,
  repl,
  trait,
  reg = NULL,
  year = NULL,
  res.het = FALSE,
  iter = 2000,
  cores = 1,
  chains = 4,
  pars = NA,
  warmup = floor(iter/2),
  thin = 1,
  seed = sample.int(.Machine$integer.max, 1),
  init = "random",
  verbose = FALSE,
  algorithm = c("NUTS", "HMC", "Fixed_param"),
  control = NULL,
  include = TRUE,
  show_messages = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame in which to interpret the variables declared in the other
  arguments.

- gen, loc:

  A string. The name of the columns that contain the evaluated
  candidates and locations (or environments, if you are working with
  factor combinations), respectively.

- repl:

  A string, a vector, or `NULL`. If the trial is randomized in complete
  blocks design, `repl` will be a string representing the name of the
  column that corresponds to the blocks. If the trial is randomized in
  incomplete blocks design, `repl` will be a string vector containing
  the name of the columns that correspond to the replicate and block
  effects on the first and second positions, respectively
  (`c(replicate, block)`). If the data does not have replicates, `repl`
  must be `NULL`.

- trait:

  A string. The analysed variable. Currently, only single-trait models
  are fitted.

- reg:

  A string or NULL. The name of the column that contain information on
  regions or mega-environments. `NULL` (default) if not applicable.

- year:

  A string or NULL. The name of the column that contain information on
  years (or seasons). `NULL` (default) if not applicable.

- res.het:

  Should the model consider heterogeneous residual variances? Defaults
  for `FALSE`. If `TRUE`, the model will estimate one residual variance
  per location (or environmnet).

- iter:

  A positive integer specifying the number of iterations for each chain
  (including warmup). The default is 2000.

- cores:

  Number of cores to use when executing the chains in parallel, which
  defaults to 1 but we recommend setting the `mc.cores` option to be as
  many processors as the hardware and RAM allow (up to the number of
  chains).

- chains:

  A positive integer specifying the number of Markov chains. The default
  is 4.

- pars:

  A vector of character strings specifying parameters of interest. The
  default is `NA` indicating all parameters in the model. If
  `include = TRUE`, only samples for parameters named in `pars` are
  stored in the fitted results. Conversely, if `include = FALSE`,
  samples for all parameters *except* those named in `pars` are stored
  in the fitted results.

- warmup:

  A positive integer specifying the number of warmup (aka burnin)
  iterations per chain. If step-size adaptation is on (which it is by
  default), this also controls the number of iterations for which
  adaptation is run (and hence these warmup samples should not be used
  for inference). The number of warmup iterations should be smaller than
  `iter` and the default is `iter/2`.

- thin:

  A positive integer specifying the period for saving samples. The
  default is 1, which is usually the recommended value.

- seed:

  The seed for random number generation. The default is generated from 1
  to the maximum integer supported by R on the machine. Even if multiple
  chains are used, only one seed is needed, with other chains having
  seeds derived from that of the first chain to avoid dependent samples.
  When a seed is specified by a number, `as.integer` will be applied to
  it. If `as.integer` produces `NA`, the seed is generated randomly. The
  seed can also be specified as a character string of digits, such as
  `"12345"`, which is converted to integer.

- init:

  Initial values specification. See the detailed documentation for the
  init argument in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- verbose:

  `TRUE` or `FALSE`: flag indicating whether to print intermediate
  output from Stan on the console, which might be helpful for model
  debugging.

- algorithm:

  One of sampling algorithms that are implemented in Stan. Current
  options are `"NUTS"` (No-U-Turn sampler, Hoffman and Gelman 2011,
  Betancourt 2017), `"HMC"` (static HMC), or `"Fixed_param"`. The
  default and preferred algorithm is `"NUTS"`.

- control:

  A named `list` of parameters to control the sampler's behavior. See
  the details in the documentation for the `control` argument in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- include:

  Logical scalar defaulting to `TRUE` indicating whether to include or
  exclude the parameters given by the `pars` argument. If `FALSE`, only
  entire multidimensional parameters can be excluded, rather than
  particular elements of them.

- show_messages:

  Either a logical scalar (defaulting to `TRUE`) indicating whether to
  print the summary of Informational Messages to the screen after a
  chain is finished or a character string naming a path where the
  summary is stored. Setting to `FALSE` is not recommended unless you
  are very sure that the model is correct up to numerical error.

- ...:

  Additional arguments can be `chain_id`, `init_r`, `test_grad`,
  `append_samples`, `refresh`, `enable_random_init`. See the
  documentation in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

## Value

An object of S4 class `stanfit` representing the fitted results. Slot
`mode` for this object indicates if the sampling is done or not.

## Details

The function has nine available models, which will be fitted according
to the options set in the arguments:

1.  Entry-mean model : fitted when `repl = NULL`, `reg = NULL` and
    `year = NULL`: \$\$y = \mu + g + l + \varepsilon\$\$ Where \\y\\ is
    the phenotype, \\\mu\\ is the intercept, \\g\\ is the genotypic
    effect, \\l\\ is the location (or environment) effect, and
    \\\varepsilon\\ is the error (which contains the
    genotype-by-location interaction, in this case).

2.  Randomized complete blocks design : fitted when `repl` is a single
    string. It will fit different models depending if `reg` and `year`
    are `NULL`:

    - `reg = NULL` and `year = NULL` : \$\$y = \mu + g + l + gl + r +
      \varepsilon\$\$ where \\gl\\ is the genotype-by-location effect,
      and \\r\\ is the replicate effect.

    - `reg = "reg"` and `year = NULL` : \$\$y = \mu + g + m + l + gl +
      gm + r + \varepsilon\$\$ where \\m\\ is the region effect, and
      \\gm\\ is the genotype-by-region effect.

    - `reg = NULL` and `year = "year"` : \$\$y = \mu + g + t + l + gl +
      gt + r + \varepsilon\$\$ where \\t\\ is the year effect, and
      \\gt\\ is the genotype-by-year effect.

    - `reg = "reg"` and `year = "year"` : \$\$y = \mu + g + m + t + l +
      gl + gm + gt + r + \varepsilon\$\$

3.  Incomplete blocks design : fitted when `repl` is a string vector of
    size 2. It will fit different models depending if `reg` and `year`
    are `NULL`:

    - `reg = NULL` and `year = NULL` : \$\$y = \mu + g + l + gl + r +
      b + \varepsilon\$\$ where \\b\\ is the block within replicates
      effect.

    - `reg = "reg"` and `year = NULL` : \$\$y = \mu + g + m + l + gl +
      gm + r + b + \varepsilon\$\$

    - `reg = NULL` and `year = "year"` : \$\$y = \mu + g + t + l + gl +
      gt + r + b + \varepsilon\$\$

    - `reg = "reg"` and `year = "year"` : \$\$y = \mu + g + m + t + l +
      gl + gm + gt + r + b + \varepsilon\$\$

The models described above have predefined priors: \$\$x \sim
\mathcal{N} \left( 0, S^{\[x\]} \right)\$\$ \$\$\sigma \sim
\mathcal{HalfCauchy}\left( 0, S^{\[\sigma\]} \right)\$\$ where \\x\\ can
be any effect but the error, and \\\sigma\\ is the standard deviation of
the likelihood. If `res.het = TRUE`, then \\\sigma_k \sim
\mathcal{HalfCauchy}\left( 0, S^{\left\[ \sigma_k \right\]} \right)\\.
The hyperpriors are set as follows: \$\$S^{\[x\]} \sim
\mathcal{HalfCauchy}\left( 0, \phi \right)\$\$ where \\\phi\\ is the
known global hyperparameter defined such as \\\phi = max(y) \times 10\\.

More details about the usage of `bayes_met` and other functions of the
`ProbBreed` package can be found at
<https://saulo-chaves.github.io/ProbBreed_site/>. Solutions to
convergence or mixing issues can be found at
<https://mc-stan.org/misc/warnings.html>.

## Methods

- `sampling`:

  `signature(object = "stanmodel")` Call a sampler (NUTS, HMC, or
  Fixed_param depending on parameters) to draw samples from the model
  defined by S4 class `stanmodel` given the data, initial values, etc.

## See also

rstan::sampling,
[rstan::stan](https://mc-stan.org/rstan/reference/stan.html),
[rstan::stanfit](https://mc-stan.org/rstan/reference/stanfit-class.html)

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
#>      174.745 seconds (Sampling)
#> Chain 2:                393.713 seconds (Total)
#> Chain 2: 
#> Warning: There were 2 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.1, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
# }
```

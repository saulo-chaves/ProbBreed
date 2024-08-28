##' @title Bayesian model for multi-environment trials
##'
##' @description
##' Fits a Bayesian multi-environment model using `rstan`, the `R` interface to `Stan`.
##'
##' @param data A data frame in which to interpret the variables declared in the other arguments.
##' @param gen,loc  A string. The name of the columns that contain the evaluated
##' candidates and locations (or environments, if you are working with factor combinations), respectively.
##' @param repl  A string, a vector, or `NULL`. If the trial is randomized in complete blocks design,
##' `repl` will be a string representing the name of the column that corresponds to the blocks.
##' If the trial is randomized in incomplete blocks design, `repl` will be a string vector
##' containing the name of the columns that correspond to the replicate and block effects on
##' the first and second positions, respectively (`c(replicate, block)`).
##' If the data does not have replicates, `repl` will be `NULL`.
##' @param trait A string. The analysed variable. Currently, only single-trait models are fitted.
##' @param reg A string or NULL. The name of the column that contain information on
##' regions or mega-environments. `NULL` (default) if not applicable.
##' @param year A string or NULL. The name of the column that contain information on
##' years (or seasons). `NULL` (default) if not applicable.
##' @param res.het Should the model consider heterogeneous residual variances?
##' Defaults for `FALSE`. If `TRUE`, the model will estimate one
##' residual variance per location (or environmnet). If `repl = NULL`, `res.het` must be `TRUE`.
##' @inheritParams rstan::sampling
##'
##' @inheritSection rstan::sampling Methods
##' @inherit rstan::sampling return
##'
##' @details
##' The function has nine available models, which will be fitted according to the
##' options set in the arguments:
##' \enumerate{
##' \item{Entry-mean model} : fitted when `repl = NULL`, `reg = NULL` and `year = NULL`:
##' \deqn{y = \mu + g + l + \varepsilon}
##' Where \eqn{y} is the phenotype, \eqn{\mu} is the intercept, \eqn{g} is the genotypic
##' effect, \eqn{l} is the location (or environment) effect, and \eqn{\varepsilon} is
##' the residue (which contains the genotype-by-location interaction, in this case).
##'
##' \item{Randomized complete blocks design} : fitted when `repl` is a single string.
##' It will fit different models depending if `reg` and `year` are `NULL`:
##'   \itemize{
##'     \item{`reg = NULL` and `year = NULL`} :
##'       \deqn{y = \mu + g + l + gl + r + \varepsilon}
##'       where \eqn{gl} is the genotype-by-location effect, and \eqn{r} is the replicate effect.
##'     \item{`reg = "reg"` and `year = NULL`} :
##'       \deqn{y = \mu + g + m + l + gl + gm + r + \varepsilon}
##'       where \eqn{m} is the region effect, and \eqn{gm} is the genotype-by-region effect.
##'     \item{`reg = NULL` and `year = "year"`} :
##'       \deqn{y = \mu + g + t + l + gl + gt + r + \varepsilon}
##'       where \eqn{t} is the year effect, and \eqn{gt} is the genotype-by-year effect.
##'     \item{`reg = "reg"` and `year = "year"`} :
##'       \deqn{y = \mu + g + m + t + l + gl + gm + gt + r + \varepsilon}
##'   }
##'
##' \item{Incomplete blocks design} : fitted when `repl` is a string vector of size 2.
##' It will fit different models depending if `reg` and `year` are `NULL`:
##'   \itemize{
##'     \item{`reg = NULL` and `year = NULL`} :
##'       \deqn{y = \mu + g + l + gl + r + b + \varepsilon}
##'       where \eqn{b} is the block within replicates effect.
##'     \item{`reg = "reg"` and `year = NULL`} :
##'       \deqn{y = \mu + g + m + l + gl + gm + r + b + \varepsilon}
##'     \item{`reg = NULL` and `year = "year"`} :
##'       \deqn{y = \mu + g + t + l + gl + gt + r + b + \varepsilon}
##'     \item{`reg = "reg"` and `year = "year"`} :
##'       \deqn{y = \mu + g + m + t + l + gl + gm + gt + r + b + \varepsilon}
##'   }
##'
##' }
##' The models described above have predefined priors:
##' \deqn{x \sim \mathcal{N} \left( 0, S^{[x]} \right)}
##' \deqn{\sigma \sim \mathcal{HalfCauchy}\left( 0, S^{[\sigma]} \right)}
##' where \eqn{x} can be any effect but the error, and \eqn{\sigma} is the standard
##' deviation of the likelihood. If `res.het = TRUE`, then \eqn{\sigma_k \sim \mathcal{HalfCauchy}\left( 0, S^{\left[ \sigma_k \right]} \right)}.
##' The hyperpriors are set as follows:
##' \deqn{S^{[x]} \sim \mathcal{HalfCauchy}\left( 0, \phi \right)}
##' where \eqn{\phi} is the known global hyperparameter defined such as \eqn{\phi = max(y) \times 10}.
##'
##' More details about the usage of `bayes_met` and other functions of
##' the `ProbBreed` package can be found at \url{https://saulo-chaves.github.io/ProbBreed_site/}.
##' Solutions to convergence or mixing issues can be found at
##' \url{https://mc-stan.org/misc/warnings.html}.
##'
##' @seealso [rstan::sampling], [rstan::stan], [rstan::stanfit]
##'
##' @import rstan
##' @importFrom stats density median model.matrix na.exclude quantile reorder sd var
##' @importFrom utils globalVariables
##'
##' @examples
##' \donttest{
##' mod = bayes_met(data = maize,
##'                 gen = "Hybrid",
##'                 loc = "Location",
##'                 repl = c("Rep","Block"),
##'                 trait = "GY",
##'                 reg = "Region",
##'                 year = NULL,
##'                 res.het = TRUE,
##'                 iter = 2000, cores = 2, chain = 4)
##' }
##' @export

bayes_met = function(data, gen, loc, repl, trait, reg = NULL, year = NULL,
                     res.het = FALSE, iter = 2000, cores = 1, chains = 4,
                     pars = NA, warmup = floor(iter/2), thin = 1,
                     seed = sample.int(.Machine$integer.max, 1),
                     init = 'random', verbose = FALSE,
                     algorithm = c("NUTS", "HMC", "Fixed_param"),
                     control = NULL, include = TRUE, show_messages = TRUE, ...)
  {

  requireNamespace('rstan')

  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data

  stopifnot("gen is not in the data" = gen %in% colnames(data))
  stopifnot("loc is not in the data" = loc %in% colnames(data))
  stopifnot("Please, specify the trait" = trait %in% colnames(data))
  if(is.null(repl)){
    stopifnot("By default, a model with no replicates can only be fitted with heterogeneous residual variances (res.het = TRUE)" = res.het)
    stopifnot("'reg' must be NULL to fit an entry-mean model (G + L + error)" = is.null(reg))
    stopifnot("'year' must be NULL to fit an entry-mean model (G + L + error)" = is.null(year))
  }

  if(!all(grepl('[A-Za-z]', data[, gen]))){data[,gen] = paste("G", data[,gen], sep = "_")}
  if(!all(grepl('[A-Za-z]', data[, loc]))){data[,loc] = paste("E", data[,loc], sep = "_")}

  if(res.het){
  # Heterogeneous residual variances -----------------
    if(is.null(year)){
      # No year effect --------------------------
      if(is.null(reg)){
        # No region effect ---------------------------
        if(is.null(repl)){
          # Only means --------------------------------
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          # Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          # p5 <- ncol(Z5)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n,
                         p3 = p3,
                         p4 = p4,
                         # p5 = p5,
                         Z3 = Z3,
                         Z4 = Z4,
                         #Z5 = Z5,
                         index = index,
                         y = y,
                         phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET: entry-mean")

          Model = rstan::sampling(object = stanmodels$entrymean, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = TRUE

        }else if(length(repl) == 1){
          # RCB ------------------------
          stopifnot("repl is not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4,
                         p5 = p5, Z1 = Z1, Z3 = Z3, Z4 = Z4,Z5 = Z5,
                         index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD1, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2){
          # incomplete blocks ------------------------
          stopifnot("repl effects are not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                         Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5,
                         index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD1, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }
      }else{
        # With region effect ------------------------
        if(length(repl) == 1){
          # RCB -------------------------
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, Z1 = Z1, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                         Z7 = Z7, index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD2, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2){
          # incomplete blocks ------------------------
          stopifnot("repl is are not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                         Z7 = Z7, index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD2, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE
        }
      }

    }else # With year effect --------------------------
      {
      stopifnot("year is not in the data" = year %in% colnames(data))

        if(!all(grepl('[A-Za-z]', data[, year]))){data[,year] = paste("Y", data[, year], sep = "_")}

      if(is.null(reg)) # No region effect ---------------------------
        {
        if(length(repl) == 1) # RCBD ------------------------
          {
          stopifnot("repl is not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p8 = p8,
                         p9 = p9,  Z1 = Z1, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z8 = Z8,
                         Z9 = Z9, index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD3, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2) # Incomplete blocks ------------------------
          {
          stopifnot("repl are not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                         p8 = p8, p9 = p9, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4,
                         Z5 = Z5, Z8 = Z8, Z9 = Z9, index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD3, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }
      }else # With region effect ------------------------
        {
          stopifnot("reg is not in the data" = reg %in% colnames(data))

          if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

        if(length(repl) == 1) # RCB -------------------------
          {
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, p8 = p8, p9 = p9, Z1 = Z1, Z3 = Z3, Z4 = Z4,
                         Z5 = Z5, Z6 = Z6, Z7 = Z7, Z8 = Z8, Z9 = Z9,
                         index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD4, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2) # incomplete blocks ------------------------
          {
          stopifnot("repl is are not in the data" = repl %in% colnames(data))
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, p8 = p8, p9 = p9, Z1 = Z1, Z2 = Z2, Z3 = Z3,
                         Z4 = Z4, Z5 = Z5, Z6 = Z6, Z7 = Z7, Z8 = Z8, Z9 = Z9,
                         index = index, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD4, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE
        }
      }
    }

}else # Homogeneous variances ------------------------------
  {
  if(is.null(year)) # No year effect ------------------------
    {
    if(is.null(reg)) # No region information -----------------------------
      {
      if(length(repl) == 1) # RCB ---------------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,loc] = as.factor(data[,loc])
        data[,repl] = as.factor(data[,repl])
        n = nrow(data)
        Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        p1 <- ncol(Z1)
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4,
                       p5 = p5, Z1 = Z1, Z3 = Z3, Z4 = Z4,Z5 = Z5, y = y, phi = phi)

        # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(object = stanmodels$RCBD5, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages,...)
        attr(Model, 'modmean') = FALSE

      }else if(length(repl) == 2) # incomplete blocks --------------------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,loc] = as.factor(data[,loc])
        data[,repl[1]] = as.factor(data[,repl[1]])
        data[,repl[2]] = as.factor(data[,repl[2]])
        n = nrow(data)
        Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
        Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        p1 <- ncol(Z1)
        p2 <- ncol(Z2)
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                       Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5, y = y, phi = phi)

        # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(object = stanmodels$IBD5, data = df_stan,iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)
        attr(Model, 'modmean') = FALSE

      }
    }else # With region information -------------------------
      {
        if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}
      if(length(repl) == 1) # RCDB --------------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,reg] = as.factor(data[,reg])
        data[,loc] = as.factor(data[,loc])
        data[,repl] = as.factor(data[,repl])
        n = nrow(data)
        Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        Z6 = model.matrix(~-1 + data[,reg])
        Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
        p1 <- ncol(Z1)
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        p6 <- ncol(Z6)
        p7 <- ncol(Z7)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                       p7 = p7, Z1 = Z1, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                       Z7 = Z7, y = y, phi = phi)

        # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(object = stanmodels$RCBD6, data = df_stan,iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)
        attr(Model, 'modmean') = FALSE

      }else if(length(repl) == 2) # Incomplete blocks --------------------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,reg] = as.factor(data[,reg])
        data[,loc] = as.factor(data[,loc])
        data[,repl[1]] = as.factor(data[,repl[1]])
        data[,repl[2]] = as.factor(data[,repl[2]])
        n = nrow(data)
        Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
        Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        Z6 = model.matrix(~-1 + data[,reg])
        Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
        p1 <- ncol(Z1)
        p2 <- ncol(Z2)
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        p6 <- ncol(Z6)
        p7 <- ncol(Z7)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                       p7 = p7, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                       Z7 = Z7, y = y, phi = phi)

        # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(object = stanmodels$IBD6, data = df_stan,iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)
        attr(Model, 'modmean') = FALSE
      }
    }
  }else # With year effect ------------------
    {
      if(!all(grepl('[A-Za-z]', data[, year]))){data[,year] = paste("Y", data[, year], sep = "_")}

      if(is.null(reg)) # No region information -----------------------------
      {
        if(length(repl) == 1) # RCB ---------------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p8 = p8,
                         p9 = p9, Z1 = Z1, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z8 = Z8,
                         Z9 = Z9, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD7, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2) # incomplete blocks --------------------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                         p8 = p8, p9 = p9, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4,
                         Z5 = Z5, Z8 = Z8, Z9 = Z9, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD7, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }
      }else # With region information -------------------------
      {
        if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

        if(length(repl) == 1) # RCDB --------------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl] = as.factor(data[,repl])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, p8 = p8, p9 = p9, Z1 = Z1, Z3 = Z3, Z4 = Z4,
                         Z5 = Z5, Z6 = Z6, Z7 = Z7, Z8 = Z8, Z9 = Z9, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$RCBD8, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE

        }else if(length(repl) == 2) # Incomplete blocks --------------------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          data[,repl[1]] = as.factor(data[,repl[1]])
          data[,repl[2]] = as.factor(data[,repl[2]])
          n = nrow(data)
          Z1 = model.matrix(~-1 + data[,repl[1]]:data[,loc])
          Z2 = model.matrix(~-1 + data[,repl[2]]:data[,loc])
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p1 <- ncol(Z1)
          p2 <- ncol(Z2)
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                         p6 = p6, p7 = p7, p8 = p8, p9 = p9, Z1 = Z1, Z2 = Z2,
                         Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6, Z7 = Z7, Z8 = Z8,
                         Z9 = Z9, y = y, phi = phi)

          # stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(object = stanmodels$IBD8, data = df_stan,iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
          attr(Model, 'modmean') = FALSE
        }
      }
  }
  }

  attr(Model, "declared_eff") = data.frame(trait = trait, gen = gen,
                                           repl = ifelse(is.null(repl), 0, repl[1]),
                                           blk = ifelse(length(repl) == 2, repl[2],0),
                                           loc = loc, reg = ifelse(is.null(reg), 0, reg),
                                           year = ifelse(is.null(year), 0, year),
                                           res.het = res.het)
  attr(Model, "data") = data
  return(Model)
}


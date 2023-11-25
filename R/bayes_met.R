##' @title Bayesian model for multi-environment trials
##'
##' @description
##' This function runs a Bayesian model for analyzing data from
##' Multi-environment trials using `rstan`, the `R` interface to `Stan`.
##'
##'
##' @param data  A data frame containing the observations.
##' @param gen,loc  A string. The name of the
##' column that corresponds to the evaluated genotype and location, respectively. If
##' the environment is a combination of other factors (for instance, location-year),
##' the name of the column that contains this information must be attributed to `loc`.
##' @param repl  A string, a vector, or `NULL`. If the trial is randomized in complete blocks,
##' `repl` will be a string representing the name of the column
##' that corresponds to the blocks. If the trial is randomized in incomplete blocks design,
##' `repl` will be a string vector containing the name of the column that corresponds to
##' the replicate and block effects on the first and second positions, respectively.
##' If the data do not have replicates, `repl` will be `NULL`.
##' @param trait A string. The name of the column that corresponds to the analysed variable.
##' @param reg A string or NULL. If the data has information of regions,
##' `reg` will be a string with the name of the column that corresponds to the
##' region information. Otherwise, `reg = NULL` (default).
##' @param year A string or NULL. If the data set has information of time-related
##' environmental factors (years, seasons...), `year` will be a string with the
##' name of the column that corresponds to the time information. Otherwise, `year = NULL` (default).
##' @param res.het Logical, indicating if the model should consider heterogeneous
##' residual variances. Default is `FALSE`. If `TRUE`, the model will estimate one
##' residual variance per location.
##' @inheritParams rstan::sampling
##'
##' @inheritSection rstan::sampling Methods
##' @inherit rstan::sampling return
##'
##' @details
##' More details about the usage of `bayes_met` and other function of
##' the `ProbBreed` package can be found at \url{https://saulo-chaves.github.io/ProbBreed_site/}.
##' Information on solutions to solve convergence or mixing issue can be found at
##' \url{https://mc-stan.org/misc/warnings.html}.
##'
##' @seealso [rstan::sampling()], [rstan::stan()], [rstan::stanfit()]
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
##'                 repl = c("Rep", "Block"),
##'                 year = NULL,
##'                 reg = 'Region',
##'                 res.het = FALSE,
##'                 trait = 'GY',
##'                 iter = 6000, cores = 4, chains = 4)
##'                 }
##' @export

bayes_met = function(data, gen, loc, repl, trait, reg = NULL, year = NULL,
                     res.het = FALSE, iter = 2000, cores = 2, chains = 4,
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
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p3 = p3, p4 = p4,
                         p5 = p5, Z3 = Z3, Z4 = Z4,Z5 = Z5,
                         index = index, y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;
    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }else if(length(repl) == 2){
          # incomplete blocks ------------------------
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;
    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }
      }else{
        # With region effect ------------------------
        if(is.null(repl)){

          if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

          # Only means ------------------------
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                         Z7 = Z7, index = index, y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        } else if(length(repl) == 1){
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
        }
      }

    }else # With year effect --------------------------
      {
      stopifnot("year is not in the data" = year %in% colnames(data))

        if(!all(grepl('[A-Za-z]', data[, year]))){data[,year] = paste("Y", data[, year], sep = "_")}

      if(is.null(reg)) # No region effect ---------------------------
        {
        if(is.null(repl)) # Only means --------------------------------
          {
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p3 = ncol(Z3)
          p4 = ncol(Z4)
          p5 = ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          index = rep(1:nlevels(data[,loc]), times = as.numeric(table(data[,loc])))
          phi = max(y) * 10
          df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p8 = p8, p9 = p9,
                         Z3 = Z3, Z4 = Z4, Z5 = Z5, Z8 = Z8, Z9 = Z9,
                         index = index, y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Year parameter/hyperparameters
    real<lower=0> s_t
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }else if(length(repl) == 1) # RCBD ------------------------
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] l;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }
      }else # With region effect ------------------------
        {
          stopifnot("reg is not in the data" = reg %in% colnames(data))

          if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

        if(is.null(repl)) # Only means ------------------------
          {

          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
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
          df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                         p7 = p7, p8 = p8, p9 = p9, Z3 = Z3, Z4 = Z4, Z5 = Z5,
                         Z6 = Z6, Z7 = Z7, Z8 = Z8, Z9 = Z9, index = index,
                         y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genotype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        } else if(length(repl) == 1) # RCB -------------------------
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genotype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

    // Indexation vector for sigma
   int index[n];
  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genptype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genptype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Declaring the vector with the residual standard deviations
    vector<lower=0>[n] sigma_vec;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma_vec);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
        }
      }
    }

}else # Homogeneous variances ------------------------------
  {
  if(is.null(year)) # No year effect ------------------------
    {
    if(is.null(reg)) # No region information -----------------------------
      {
      if(is.null(repl)) # Only-means ---------------------------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,loc] = as.factor(data[,loc])
        n = nrow(data)
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p3 = p3, p4 = p4,
                       p5 = p5, Z3 = Z3, Z4 = Z4,Z5 = Z5, y = y, phi = phi)

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);

    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)

      }else if(length(repl) == 1) # RCB ---------------------
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

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;
    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)

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

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)

      }
    }else # With region information -------------------------
      {
        if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}
      if(is.null(repl)) # Only means --------------
        {
        data[,gen] = as.factor(data[,gen])
        data[,reg] = as.factor(data[,reg])
        data[,loc] = as.factor(data[,loc])
        n = nrow(data)
        Z3 = model.matrix(~-1 + data[,gen])
        Z4 = model.matrix(~-1 + data[,loc])
        Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
        Z6 = model.matrix(~-1 + data[,reg])
        Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
        p3 <- ncol(Z3)
        p4 <- ncol(Z4)
        p5 <- ncol(Z5)
        p6 <- ncol(Z6)
        p7 <- ncol(Z7)
        y = data[,trait]
        phi = max(y) * 10
        df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p6 = p6,
                       p7 = p7, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                       Z7 = Z7, y = y, phi = phi)

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)

      } else if(length(repl) == 1) # RCDB --------------------
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

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)

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

        stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);

    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "
        stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

        Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                cores = cores, chains = chains, pars = pars,
                                warmup = warmup, thin = thin, seed = seed,
                                init = init, verbose = verbose,
                                algorithm = algorithm, control = control,
                                include = include, show_messages = show_messages, ...)
      }
    }
  }else # With year effect ------------------
    {
      if(!all(grepl('[A-Za-z]', data[, year]))){data[,year] = paste("Y", data[, year], sep = "_")}

      if(is.null(reg)) # No region information -----------------------------
      {
        if(is.null(repl)) # Only-means ---------------------------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p8 = p8, p9 = p9,
                         Z3 = Z3, Z4 = Z4,Z5 = Z5, Z8 = Z8, Z9 = Z9, y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Hybrid by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);

    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }else if(length(repl) == 1) # RCB ---------------------
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        }
      }else # With region information -------------------------
      {
        if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

        if(is.null(repl)) # Only means --------------
        {
          data[,gen] = as.factor(data[,gen])
          data[,reg] = as.factor(data[,reg])
          data[,loc] = as.factor(data[,loc])
          n = nrow(data)
          Z3 = model.matrix(~-1 + data[,gen])
          Z4 = model.matrix(~-1 + data[,loc])
          Z5 = model.matrix(~-1 + data[,gen]:data[,loc])
          Z6 = model.matrix(~-1 + data[,reg])
          Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
          Z8 = model.matrix(~-1 + data[,year])
          Z9 = model.matrix(~-1 + data[,gen]:data[,year])
          p3 <- ncol(Z3)
          p4 <- ncol(Z4)
          p5 <- ncol(Z5)
          p6 <- ncol(Z6)
          p7 <- ncol(Z7)
          p8 = ncol(Z8)
          p9 = ncol(Z9)
          y = data[,trait]
          phi = max(y) * 10
          df_stan = list(n = n, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7,
                         p8 = p8, p9 = p9, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
                         Z7 = Z7, Z8 = Z8, Z9 = Z9, y = y, phi = phi)

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genotype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

        } else if(length(repl) == 1) # RCDB --------------------
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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genotype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);
    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "

          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)

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

          stan_df = "
  data{
    // Number of observations
    int<lower=1> n;

    // Number of parameters
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;
    int<lower=1> p8;
    int<lower=1> p9;

    // Designs matrices
    matrix[n, p1] Z1;
    matrix[n, p2] Z2;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;
    matrix[n, p8] Z8;
    matrix[n, p9] Z9;

    // Phenotype vector
    real y[n];

    // Global hyperparameter
    real phi;

  }
    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p1] r;

    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p2] b;

    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p3] g;

    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p4] l;

    // Genotype by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p5] gl;

    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p6] m;

    // Genotype by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p7] gm;

    // year parameter/hyperparameters
    real<lower=0> s_t;
    vector[p8] t;

    // Genotype by year parameter/hyperparameters
    real<lower=0> s_gt;
    vector[p9] gt;

    // Defining variable to generate data from the model
    real y_gen[n];
  }

  transformed parameters{

    // Declaring variables to receive input
    vector[n] expectation;

    // Computing the expectation of the likelihood function
    expectation = mu + Z1*r + Z2*b + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm + Z8*t + Z9*gt;

  }
  model{

    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);

    // Conditional prior probabilities distributions for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    // Conditional prior probabilities distributions  for year
    s_t ~ cauchy(0, phi);
    t ~ normal(0, s_t);

    // Conditional prior probabilities distributions  for genotype by year
    s_gt ~ cauchy(0, phi);
    gt ~ normal(0, s_gt);

    // Specifying the likelihood
    y ~ normal(expectation, sigma);

    // Generating data from the model
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

} "
          stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

          Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter,
                                  cores = cores, chains = chains, pars = pars,
                                  warmup = warmup, thin = thin, seed = seed,
                                  init = init, verbose = verbose,
                                  algorithm = algorithm, control = control,
                                  include = include, show_messages = show_messages, ...)
        }
      }
  }
}
return(Model)
}


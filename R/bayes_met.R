## Function bayes_met
##
##' @title
##' Bayesian model for multi-environment trials
##'
##' @description
##' This function runs a Bayesian model for analyzing data from
##' Multi-environment trials using `rstan`, the `R` interface to `Stan`.
##'
##' @details
##' More details about the usage of `bayes_met`, as well as the other function of
##' the `ProbBreed` package can be found at \url{https://saulo-chaves.github.io/ProbBreed_site/}.
##'
##' @param data  A data frame containing the observations
##' @param gen,env  A string. The name of the
##' column that corresponds to the evaluated genotype or environment.
##' @param repl  A string, a vector, or `NULL`. If the trial is randomized in complete blocks,
##' `repl` will be a string representing the name of the column
##' that corresponds to the blocks. If the trial is randomized in incomplete blocks design,
##' `repl` will be a string vector containing the name of the column that corresponds to
##' the replicate and block effects on the first and second positions, respectively.
##' If the data do not have replicates, `repl` will be `NULL`.
##' @param trait A string. The name of the column that corresponds to the analysed variable.
##' @param reg A string or NULL. If the data set has information about regions,
##' `reg` will be a string with the name of the column that corresponds to the
##' region information. Otherwise, `reg = NULL` (default).
##' @param res.het Logical, indicating if the model should consider heterogeneous
##' residual variances. Default is `FALSE`.
##' @param chains Inherited from [rstan::sampling()].
##' A positive integer specifying the number of Markov chains. The default is 4.
##' @param iter Inherited from [rstan::sampling()].
##' A positive integer specifying the number of iterations for each chains
##' (including warmup). The default is 2000.
##' @param cores Inherited from [rstan::sampling()].
##' A positive integer specifying the number of cores to use when executing the
##' chains in parallel (defaults to 1).
##' @param ... Additional arguments passed to the [rstan::sampling()] function
##' (for instance, to change the thin number, or to set a specific seed).
##' For more information, see [rstan::sampling()] manual
##' @return The function returns an object of S4 class stanfit containing the
##' fitted results
##'
##' @seealso [rstan::sampling()]
##'
##' @import rstan
##' @importFrom stats density median model.matrix na.exclude quantile reorder sd var
##' @importFrom utils globalVariables
##'
##'
##' @examples
##' \donttest{
##' mod = bayes_met(data = soy,
##'                 gen = "Gen",
##'                 env = "Env",
##'                 repl = NULL,
##'                 reg = "Reg",
##'                 res.het = FALSE,
##'                 trait = "Y",
##'                 iter = 2000, cores = 1, chains = 4)
##'                 }
##' @export

bayes_met = function(data, gen, env, repl, trait, reg = NULL, res.het = FALSE,
                    iter = 2000, cores = 1, chains = 4,...){

  requireNamespace('rstan')

  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data

  stopifnot("gen is not in the data" = gen %in% colnames(data))
  stopifnot("env is not in the data" = env %in% colnames(data))

  if(res.het){
  # Heterogeneous residual variances
  if(is.null(reg)){
    # no region effect
    if(is.null(repl)){
      # Only means
      data[,gen] = as.factor(data[,gen])
      data[,env] = as.factor(data[,env])
      n = nrow(data)
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

    }else if(length(repl) == 1){
      # RCB
      stopifnot("repl is not in the data" = repl %in% colnames(data))
      data[,gen] = as.factor(data[,gen])
      data[,env] = as.factor(data[,env])
      data[,repl] = as.factor(data[,repl])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,repl]:data[,env])
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
      p1 <- ncol(Z1)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

    }else if(length(repl) == 2){
      # incomplete blocks
      stopifnot("repl are not in the data" = repl %in% colnames(data))
      data[,gen] = as.factor(data[,gen])
      data[,env] = as.factor(data[,env])
      data[,repl[1]] = as.factor(data[,repl[1]])
      data[,repl[2]] = as.factor(data[,repl[2]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,repl[1]]:data[,env])
      Z2 = model.matrix(~-1 + data[,repl[2]]:data[,env])
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
      p1 <- ncol(Z1)
      p2 <- ncol(Z2)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

    }
  }else{
    # With region effect
    if(is.null(repl)){
      # Only means
      data[,gen] = as.factor(data[,gen])
      data[,reg] = as.factor(data[,reg])
      data[,env] = as.factor(data[,env])
      n = nrow(data)
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
      Z6 = model.matrix(~-1 + data[,reg])
      Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      p6 <- ncol(Z6)
      p7 <- ncol(Z7)
      y = data[,trait]
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

    } else if(length(repl) == 1){
      # RCB
      data[,gen] = as.factor(data[,gen])
      data[,reg] = as.factor(data[,reg])
      data[,env] = as.factor(data[,env])
      data[,repl] = as.factor(data[,repl])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,repl]:data[,env])
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
      Z6 = model.matrix(~-1 + data[,reg])
      Z7 = model.matrix(~-1 + data[,gen]:data[,reg])
      p1 <- ncol(Z1)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      p6 <- ncol(Z6)
      p7 <- ncol(Z7)
      y = data[,trait]
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

    }else if(length(repl) == 2){
      # incomplete blocks
      stopifnot("repl is are not in the data" = repl %in% colnames(data))
      data[,gen] = as.factor(data[,gen])
      data[,reg] = as.factor(data[,reg])
      data[,env] = as.factor(data[,env])
      data[,repl[1]] = as.factor(data[,repl[1]])
      data[,repl[2]] = as.factor(data[,repl[2]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,repl[1]]:data[,env])
      Z2 = model.matrix(~-1 + data[,repl[2]]:data[,env])
      Z3 = model.matrix(~-1 + data[,gen])
      Z4 = model.matrix(~-1 + data[,env])
      Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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
      index = rep(1:nlevels(data[,env]), times = as.numeric(table(data[,env])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)
    }
  }
}else{
# Homogeneous variances
if(is.null(reg)){
  # No region information
  if(is.null(repl)){
    # Only-means
    data[,gen] = as.factor(data[,gen])
    data[,env] = as.factor(data[,env])
    n = nrow(data)
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

  }else if(length(repl) == 1){
  # RCB
    data[,gen] = as.factor(data[,gen])
    data[,env] = as.factor(data[,env])
    data[,repl] = as.factor(data[,repl])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,repl]:data[,env])
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

  }else if(length(repl) == 2){
  # incomplete blocks
    data[,gen] = as.factor(data[,gen])
    data[,env] = as.factor(data[,env])
    data[,repl[1]] = as.factor(data[,repl[1]])
    data[,repl[2]] = as.factor(data[,repl[2]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,repl[1]]:data[,env])
    Z2 = model.matrix(~-1 + data[,repl[2]]:data[,env])
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

  }
}else{
  # With region information
  if(is.null(repl)){
    data[,gen] = as.factor(data[,gen])
    data[,reg] = as.factor(data[,reg])
    data[,env] = as.factor(data[,env])
    n = nrow(data)
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

  } else if(length(repl) == 1){
  # RCDB
    data[,gen] = as.factor(data[,gen])
    data[,reg] = as.factor(data[,reg])
    data[,env] = as.factor(data[,env])
    data[,repl] = as.factor(data[,repl])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,repl]:data[,env])
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)

  }else if(length(repl) == 2){
  # incomplete blocks
    data[,gen] = as.factor(data[,gen])
    data[,reg] = as.factor(data[,reg])
    data[,env] = as.factor(data[,env])
    data[,repl[1]] = as.factor(data[,repl[1]])
    data[,repl[2]] = as.factor(data[,repl[2]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,repl[1]]:data[,env])
    Z2 = model.matrix(~-1 + data[,repl[2]]:data[,env])
    Z3 = model.matrix(~-1 + data[,gen])
    Z4 = model.matrix(~-1 + data[,env])
    Z5 = model.matrix(~-1 + data[,gen]:data[,env])
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chains = chains, ...)
    }
  }
}
return(Model)

}


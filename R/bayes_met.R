## Function bayes_met
##
##' @title
##' Bayesian model for multi-environment trials
##'
##' @description
##' This function runs a Bayesian model for analysing data from
##' Multi-environment trials using `rstan`, the `R` interface to `Stan`.
##'
##'
##' @param data  A data frame containing the observations
##' @param gen,env  A string vector. The first element is the name of the
##' column that corresponds to the evaluated genotype or environment.
##' @param rept  A vector, a list or `NULL`. If the trial is randomized in complete blocks,
##' `rept` will be a vector with the first element representing the name of the column
##' that corresponds to the blocks. If the trial is randomized in a lattice,
##' `rept` will be a list containing two vectors with the same structure previously
##' mentioned: a vector for replicate effects, and a vector for block effects.
##' If the data do not have replicates, `rept` will be `NULL`.
##' @param trait A character representing the name of the column that
##' corresponds to the analysed variable.
##' @param reg A list containing a string vector with the name
##' of the column that contain the information of region in the first position.
##' If there are no information about regions in
##' the data set, `reg = NULL` (default).
##' @param res.het Logical, indicating if the model should consider heterogeneous
##' residual variances. Default is `FALSE`.
##' @param chain Inherited from [rstan::sampling()].
##' A positive integer specifying the number of Markov chains. The default is 4.
##' @param iter Inherited from [rstan::sampling()].
##' A positive integer specifying the number of iterations for each chain
##' (including warmup). The default is 2000.
##' @param cores Inherited from [rstan::sampling()].
##' A positive integer specifying the number of cores to use when executing the
##' chains in parallel (defaults to 1).
##' @param ... Additional arguments passed to the [rstan::sampling()] function.
##' For more information, see [rstan::sampling()] manual
##' @return The function returns an object of S4 class stanfit containing the
##' fitted results
##'
##' @seealso [rstan::sampling()]
##'
##'
##' @import rstan
##' @importFrom stats density median model.matrix na.exclude quantile reorder sd var
##' @importFrom utils globalVariables
##'
##'
##' @examples
##' \dontrun{
##' # Model
##' mod = bayes_met(data = soy,
##'                 gen = "Gen",
##'                 env = "Env",
##'                 rept = NULL,
##'                 reg = "Reg",
##'                 res.het = F,
##'                 trait = "eBLUE", 
##'                 iter = 100, cores = 4, chain = 4)
##'                 #You may want to increase the number of iterations, cores and chains
##'                 }
##' @export

bayes_met = function(data, gen, env, rept, trait, reg = NULL, res.het = F,
                    iter = 2000, cores, chain,...){

  requireNamespace('rstan')

  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data

if(res.het){
  if(is.null(reg)){
    # no region effect
    if(is.null(rept)){
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,env[1]] = as.factor(data[,env[1]])
      n = nrow(data)
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

    }else if(!is.list(rept)){
      # DBC
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,env[1]] = as.factor(data[,env[1]])
      data[,rept[1]] = as.factor(data[,rept[1]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,rept[1]]:data[,env[1]])
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      p1 <- ncol(Z1)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

} ")

      stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

    }else if(is.list(rept)){
      # Látice
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,env[1]] = as.factor(data[,env[1]])
      data[,rept[[1]][1]] = as.factor(data[,rept[[1]][1]])
      data[,rept[[2]][1]] = as.factor(data[,rept[[2]][1]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,rept[[1]][1]]:data[,env[1]])
      Z2 = model.matrix(~-1 + data[,rept[[2]][1]]:data[,env[1]])
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      p1 <- ncol(Z1)
      p2 <- ncol(Z2)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

} ")

      stan_df_comp = rstan::stan_model(model_code = stan_df, model_name = "BayesMET")

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

    }
  }else{
    # Com região
    if(is.null(rept)){
      # média
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
      data[,env[1]] = as.factor(data[,env[1]])
      n = nrow(data)
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
      Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      p6 <- ncol(Z6)
      p7 <- ncol(Z7)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

    } else if(!is.list(rept)){
      # DBC
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
      data[,env[1]] = as.factor(data[,env[1]])
      data[,rept[1]] = as.factor(data[,rept[1]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,rept[1]]:data[,env[1]])
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
      Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
      p1 <- ncol(Z1)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      p6 <- ncol(Z6)
      p7 <- ncol(Z7)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

    }else if(is.list(rept)){
      # Látice
      data[,gen[1]] = as.factor(data[,gen[1]])
      data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
      data[,env[1]] = as.factor(data[,env[1]])
      data[,rept[[1]][1]] = as.factor(data[,rept[[1]][1]])
      data[,rept[[2]][1]] = as.factor(data[,rept[[2]][1]])
      n = nrow(data)
      Z1 = model.matrix(~-1 + data[,rept[[1]][1]]:data[,env[1]])
      Z2 = model.matrix(~-1 + data[,rept[[2]][1]]:data[,env[1]])
      Z3 = model.matrix(~-1 + data[,gen[1]])
      Z4 = model.matrix(~-1 + data[,env[1]])
      Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
      Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
      Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
      p1 <- ncol(Z1)
      p2 <- ncol(Z2)
      p3 <- ncol(Z3)
      p4 <- ncol(Z4)
      p5 <- ncol(Z5)
      p6 <- ncol(Z6)
      p7 <- ncol(Z7)
      y = data[,trait]
      index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

      Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)
    }
  }
}else{

if(is.null(reg)){
  if(is.null(rept)){
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,env[1]] = as.factor(data[,env[1]])
    n = nrow(data)
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

  }else if(!is.list(rept)){
  # RCBD
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,env[1]] = as.factor(data[,env[1]])
    data[,rept[1]] = as.factor(data[,rept[1]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,rept[1]]:data[,env[1]])
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    p1 <- ncol(Z1)
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

  }else if(is.list(rept)){
  # IBD
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,env[1]] = as.factor(data[,env[1]])
    data[,rept[[1]][1]] = as.factor(data[,rept[[1]][1]])
    data[,rept[[2]][1]] = as.factor(data[,rept[[2]][1]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,rept[[1]][1]]:data[,env[1]])
    Z2 = model.matrix(~-1 + data[,rept[[2]][1]]:data[,env[1]])
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    p1 <- ncol(Z1)
    p2 <- ncol(Z2)
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

  }
}else{
  if(is.null(rept)){
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
    data[,env[1]] = as.factor(data[,env[1]])
    n = nrow(data)
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
    Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    p6 <- ncol(Z6)
    p7 <- ncol(Z7)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

  } else if(!is.list(rept)){
  # RCDB
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
    data[,env[1]] = as.factor(data[,env[1]])
    data[,rept[1]] = as.factor(data[,rept[1]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,rept[1]]:data[,env[1]])
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
    Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
    p1 <- ncol(Z1)
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    p6 <- ncol(Z6)
    p7 <- ncol(Z7)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)

  }else if(is.list(rept)){
  # IBD
    data[,gen[1]] = as.factor(data[,gen[1]])
    data[,reg[[1]][1]] = as.factor(data[,reg[[1]][1]])
    data[,env[1]] = as.factor(data[,env[1]])
    data[,rept[[1]][1]] = as.factor(data[,rept[[1]][1]])
    data[,rept[[2]][1]] = as.factor(data[,rept[[2]][1]])
    n = nrow(data)
    Z1 = model.matrix(~-1 + data[,rept[[1]][1]]:data[,env[1]])
    Z2 = model.matrix(~-1 + data[,rept[[2]][1]]:data[,env[1]])
    Z3 = model.matrix(~-1 + data[,gen[1]])
    Z4 = model.matrix(~-1 + data[,env[1]])
    Z5 = model.matrix(~-1 + data[,gen[1]]:data[,env[1]])
    Z6 = model.matrix(~-1 + data[,reg[[1]][1]])
    Z7 = model.matrix(~-1 + data[,gen[1]]:data[,reg[[1]][1]])
    p1 <- ncol(Z1)
    p2 <- ncol(Z2)
    p3 <- ncol(Z3)
    p4 <- ncol(Z4)
    p5 <- ncol(Z5)
    p6 <- ncol(Z6)
    p7 <- ncol(Z7)
    y = data[,trait]
    index = rep(1:nlevels(data[,env[1]]), times = as.numeric(table(data[,env[1]])))
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

    Model = rstan::sampling(stan_df_comp, data = df_stan, iter = iter, cores = cores, chain = chain, ...)
    }
  }
}
return(Model)

}


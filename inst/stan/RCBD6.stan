// RCBD6: no year, with region - Homogeneous variances

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
    array[n] real y;

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
    array[n] real y_gen;
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
    array[n] real y_log_like;
      for (j in 1:n) {
      // Computing log-likelihood of the observed data:
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

}


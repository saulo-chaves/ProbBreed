// IBD7 - with year, no region - Homogeneous variances

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

}

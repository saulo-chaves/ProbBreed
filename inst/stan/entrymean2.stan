// Entry-mean - Homogeneous variances

data{
    int<lower=1> n;

    int<lower=1> p3;
    int<lower=1> p4;

    matrix[n, p3] Z3;
    matrix[n, p4] Z4;

    array[n] real y;

    real phi;

  }
    parameters{
    real<lower=0> s_sigma;
    real<lower=0> sigma;

    real<lower=0> s_mu;
    real mu;

    real<lower=0> s_g;
    vector[p3] g;

    real<lower=0> s_l;
    vector[p4] l;

    array[n] real y_gen;
  }

  transformed parameters{

    vector[n] expectation;

    expectation = mu + Z3*g + Z4*l;

  }
  model{

    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    y ~ normal(expectation, sigma);
    y_gen ~ normal(expectation, sigma);

  }

  generated quantities {
    array[n] real y_log_like;
      for (j in 1:n) {
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }

}

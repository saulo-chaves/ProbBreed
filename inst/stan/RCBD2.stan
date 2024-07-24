// RCBD2: No year, with region - Heterogeneous variances

data{
    int<lower=1> n;

    int<lower=1> p1;
    int<lower=1> p3;
    int<lower=1> p4;
    int<lower=1> p5;
    int<lower=1> p6;
    int<lower=1> p7;

    matrix[n, p1] Z1;
    matrix[n, p3] Z3;
    matrix[n, p4] Z4;
    matrix[n, p5] Z5;
    matrix[n, p6] Z6;
    matrix[n, p7] Z7;

    real y[n];

    real phi;

   int index[n];
  }
    parameters{
    real<lower=0> s_sigma;
    real<lower=0> sigma[p4];

    real<lower=0> s_mu;
    real mu;

    real<lower=0> s_r;
    vector[p1] r;

    real<lower=0> s_g;
    vector[p3] g;

    real<lower=0> s_l;
    vector[p4] l;

    real<lower=0> s_gl;
    vector[p5] gl;

    real<lower=0> s_m;
    vector[p6] m;

    real<lower=0> s_gm;
    vector[p7] gm;

    real y_gen[n];
  }

  transformed parameters{

    vector[n] expectation;

    vector<lower=0>[n] sigma_vec;

    expectation = mu + Z1*r + Z3*g + Z4*l + Z5*gl + Z6*m + Z7*gm;

    sigma_vec = to_vector(sigma[index]);
  }
  model{

    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma);

    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);

    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);

    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);

    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);

    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);

    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);

    y ~ normal(expectation, sigma_vec);
    y_gen ~ normal(expectation, sigma_vec);

  }

  generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma_vec[j]);
      }

}


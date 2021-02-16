// generated with brms 2.14.0
functions {
  /* turn a vector into a matrix of defined dimension 
   * stores elements in row major order
   * Args: 
   *   X: a vector 
   *   N: first dimension of the desired matrix
   *   K: second dimension of the desired matrix 
   * Returns: 
   *   a matrix of dimension N x K 
   */ 
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
    }
    return Y; 
  } 
 /* compute correlated group-level effects
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  /* zero-inflated negative binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_neg_binomial_lpmf(int y, real mu, real phi, 
                                       real zi) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                         neg_binomial_2_lpmf(0 | mu, phi)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             neg_binomial_2_lpmf(y | mu, phi); 
    } 
  } 
  /* zero-inflated negative binomial log-PDF of a single response 
   * logit parameterization of the zero-inflation part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_neg_binomial_logit_lpmf(int y, real mu, 
                                             real phi, real zi) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                         bernoulli_logit_lpmf(0 | zi) + 
                         neg_binomial_2_lpmf(0 | mu, phi)); 
    } else { 
      return bernoulli_logit_lpmf(0 | zi) +  
             neg_binomial_2_lpmf(y | mu, phi); 
    } 
  }
  /* zero-inflated negative binomial log-PDF of a single response 
   * log parameterization for the negative binomial part
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial distribution 
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_neg_binomial_log_lpmf(int y, real eta, 
                                           real phi, real zi) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                         neg_binomial_2_log_lpmf(0 | eta, phi)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             neg_binomial_2_log_lpmf(y | eta, phi); 
    } 
  } 
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * logit parameterization of the zero-inflation part
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial distribution 
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_neg_binomial_log_logit_lpmf(int y, real eta, 
                                                 real phi, real zi) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                         bernoulli_logit_lpmf(0 | zi) + 
                         neg_binomial_2_log_lpmf(0 | eta, phi)); 
    } else { 
      return bernoulli_logit_lpmf(0 | zi) +  
             neg_binomial_2_log_lpmf(y | eta, phi); 
    } 
  }
  // zero_inflated negative binomial log-CCDF and log-CDF functions
  real zero_inflated_neg_binomial_lccdf(int y, real mu, real phi, real hu) { 
    return bernoulli_lpmf(0 | hu) + neg_binomial_2_lccdf(y | mu, phi);
  }
  real zero_inflated_neg_binomial_lcdf(int y, real mu, real phi, real hu) { 
    return log1m_exp(zero_inflated_neg_binomial_lccdf(y | mu, phi, hu));
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  vector[N] offsets;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
  real<lower=0,upper=1> zi;  // zero-inflation probability
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + offsets;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    for (n in 1:N) {
      target += zero_inflated_neg_binomial_log_lpmf(Y[n] | mu[n], shape, zi);
    }
  }
  // priors including all constants
  target += normal_lpdf(b | 0, 1);
  target += normal_lpdf(Intercept | 0, 10);
  target += gamma_lpdf(shape | 0.01, 0.01);
  target += beta_lpdf(zi | 1, 1);
  target += student_t_lpdf(sd_1 | 3, 0, 3.6)
    - 2 * student_t_lccdf(0 | 3, 0, 3.6);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
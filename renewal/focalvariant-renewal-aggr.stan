/* implementation of a 2-variant renewal eqn in Stan.
 * Fit both case counts and variant sequence counts to infer
 * the selective advantage of one of the variants over the other.
 * We have N observations per region. For this, we need N-1
 * reproduction numbers and N-2 Gaussian RVs
 *
 * We allow for aggregation of the data to e.g. per-week
 * to reduce week-day patterns
 *
 * The generation-time distribution can be given by a negative binomial
 * or a (discretized) gamma distribution.
 */

functions {
  /* Generate random samples from the discretized gamma distribution.
   * This RNG function is mostly here for diagnostic purposes.
   */
  int gamma_discrete_rng(real alpha, real beta) {
    real x = gamma_rng(alpha, beta);
    int n = 0;
    while ( x > n ) {
      n += 1;
    }
    return n-1;
  }
  
  /* Convert the Gaussian random walk (parametertized by iid increments Z)
   * into a time-dependent reproduction number R_t. 
   * WARNING: the vector Times contains more time points than required 
   * for a single region. Infer the length of the R_t timeseries from Z.
   * N is the number of modeled days (without aggregation).
   */
  vector compute_Rt(data vector Times, vector Z, real Rt0, real tau) {
    int N = num_elements(Z) + 2;
    vector[N-1] Rt = Rt0 * exp(tau * cumulative_sum(append_row(0, Z)) - 0.5 * tau^2 * Times[1:N-1]);
    return Rt;
  }
  
  /* The renewal model depends on the entire history of incidence. However,
   * we truncate to a window size W. Still, to compute the first few incidences,
   * we need initial incidence "data" for time points 1-W,...,0.
   * For this, we assume that Rt was constant (and equal to Rt0) before the first 
   * time point. We then use the relation 1/R = E[exp(-r T_G)] to compute the 
   * exponential growth rate of incidence under the assumption of a constant
   * reproduction number. This leads to the equations below (using the MGF).
   */ 
  vector compute_dI_hist(real Rt0, real dI0, real mu, real phi, int W, int gentime_distr) {
    vector[W] dI_hist;
    if ( gentime_distr == 0 ) { // neg-binom
      real logx = log1p(phi/mu * (1-Rt0^(1/phi)));
      dI_hist = dI0 * exp(logx * linspaced_vector(W, 0, W-1));
    } else if ( gentime_distr == 1 ) { // gamma
      dI_hist = dI0 * exp(-phi/mu * (Rt0^(1/phi)-1) * linspaced_vector(W, 0, W-1));
    } else {
      reject("invalid generation time distribution");
    }
    return reverse(dI_hist);
  }
  
  /* Compute the incidence (dI) given the time-dependent reproduction number Rt,
   * and the generation time "Kernel". We also have to compute the history of 
   * the incidence before time 0. For this, we need the parameters of the 
   * generation-time distribution.
   */
  vector compute_dI(vector Rt, real dI0, data row_vector Kernel, real mu, real phi, int gentime_distr) {
    int N = num_elements(Rt) + 1;
    int W = num_elements(Kernel);
    vector[N-1+W] dI;
    // wirst W elements of dI will contain the history
    dI[1:W] = compute_dI_hist(Rt[1], dI0, mu, phi, W, gentime_distr);
    // then we use the renewal equation and Rt to compute the rest
    for ( n in 1:N-1 ) {
      dI[W+n] = Rt[n] * Kernel * reverse(dI[n:W+n-1]); // renewal equation
    }
    return dI;
  }
  
  /* Using the incidence and a "Kernel", we compute the 
   * expected number of cases at times 1,..,N. 
   * The kernel gives the probability that and individual infecte i days ago
   * is counted as positive today.
   */ 
  vector compute_I(vector dI, data row_vector Kernel) {
    int W = num_elements(Kernel);
    int N = num_elements(dI) - W + 1;
    vector[N] I;
    for ( n in 1:N ) {
      I[n] = Kernel * reverse(dI[n:n+W-1]);
    }
    return I;
  }
}


data {
  int R; // regions
  array[R] int N; // modeled days per region
  array[R] int NumObs; // num obs per region (NumObs <= N allowing for aggregation)
  
  array[max(NumObs),R] real<lower=0> Cases; // padded array
  array[max(NumObs),R] int<lower=0> TotalSeq; // total sequences collected
  array[max(NumObs),R] int<lower=0> VariantSeq; // variant sequences
  array[max(NumObs),R] int<lower=0> ObsTimes; // day numbers of aggregated observations
  
  int W; // window for kernel
  real mu_kernel_wt; // mean generation time E[T_G]
  real mu_kernel_mt; // idem for the mutant's generation time distribution
  real phi_kernel; // shape parameter of T_G
  int<lower=0, upper=1> gentime_distr; // 0 = neg-binom, 1 = gamma
  
  // covariates used for mixed-effects model
  int<lower=0> NumCovariates;
  matrix[R, NumCovariates] Covariates; // NB: should be standardized (mean = 0, sd = 1)
}




transformed data {
  row_vector<lower=0>[W] TransKernelWT; // transmission kernel
  row_vector<lower=0>[W] TransKernelMut; // idem for the mutant
  // a shared vector of time points used to compute Rt
  vector[max(N)] Time = linspaced_vector(max(N), 0, max(N)-1);

  for ( n in 1:W ) {
    // pre-compute kernel for transmission
    if ( gentime_distr == 0 ) { // neg-binom
      TransKernelWT[n] = exp(neg_binomial_2_lpmf(n-1 | mu_kernel_wt, phi_kernel));
      TransKernelMut[n] = exp(neg_binomial_2_lpmf(n-1 | mu_kernel_mt, phi_kernel));
    } else if ( gentime_distr == 1 ) { // gamma
      // gamma is parameterized with scale (phi) and rate (beta) parameters, and mu = phi / beta
      real beta_wt = phi_kernel / mu_kernel_wt;
      real beta_mt = phi_kernel / mu_kernel_mt;
      // discretize the distribution by integrating area between integers
      TransKernelWT[n] = exp(gamma_lcdf(n | phi_kernel, beta_wt)) 
          - exp(gamma_lcdf(n-1 | phi_kernel, beta_wt));
      TransKernelMut[n] = exp(gamma_lcdf(n | phi_kernel, beta_mt)) 
          - exp(gamma_lcdf(n-1 | phi_kernel, beta_mt));
    } else {
      reject("invalid generation time distribution");
    }
  }
  // rescale kernels to 1
  TransKernelWT /= sum(TransKernelWT);
  TransKernelMut /= sum(TransKernelMut);
}




parameters {
  
  vector[sum(N)-2*R] Z_concat; // random vector for constructing Rt
  //array[R] real<lower=0> tau; // volatility
  //real mu_tau;
  //real<lower=0> sd_tau;
  real<lower=0> tau_pop;
  
  /* we have a log-normal model for the case counts and estimate the 
   * standard deviation. WARNING: this overparameterized model has an unbounded 
   * likelihood due to "Neal's Funnel". We can fix this by setting a non-zero
   * lower bound for sigma
   */
   
  //array[R] real<lower=0> sigma_cases;
  //real<lower=0> sigma_cases_mu;
  //real<lower=0> sigma_cases_sd;
  real<lower=0.05> sigma_cases_pop; // non-zero lower bound

  array[R] real<lower=0> phi_inv_seq; // overdispersion sequencing
  real<lower=0> phi_inv_seq_mu;
  
  vector<lower=0>[R] dI0; // initial condition for incidence
  real mu_dI0;
  real<lower=0> sd_dI0;
  
  vector<lower=0, upper=1>[R] b0; // initial mutant fraction
  vector<lower=0>[R] Rt0; // initial reproduction number
  vector<lower=0>[R] fitness; // relative fitness of the variant (fitness = s+1)
  
  // hierarchical model parameters for fitness
  real mu_fitness;
  real<lower=0> sd_fitness;
  vector[NumCovariates] covariate_weights; // covariate weights
}




transformed parameters {
  // reproduction number (of wt)
  matrix[max(N)-1, R] Rt = rep_matrix(0.0, max(N)-1, R);
  
  // incidence 
  matrix[max(N)-1+W, R] dI_w = rep_matrix(0.0, max(N)-1+W, R); // wild-type
  matrix[max(N)-1+W, R] dI_m = rep_matrix(0.0, max(N)-1+W, R); // mutant

  // sampling-corrected prevalence
  matrix[max(N), R] I_w = rep_matrix(0.0, max(N), R); // wild-type
  matrix[max(N), R] I_m = rep_matrix(0.0, max(N), R); // mutant

  // compute Rt, dI, I
  for ( r in 1:R ) {
    // get right Z from Z_concat for current region
    int i1 = sum(N[:r-1])-2*r+3; // 1, n_1-1, n_2-3, ...
    int i2 = sum(N[:r])-2*r; // n_1-2, n_2-4, n_3-6, ..., n_R-2*R
    vector[N[r]-2] Z = Z_concat[i1:i2];
    
    // compute incidence
    Rt[1:N[r]-1, r] = compute_Rt(Time, Z, Rt0[r], tau_pop);
    dI_w[1:N[r]-1+W, r] = compute_dI(Rt[1:N[r]-1, r], (1-b0[r])*dI0[r], 
        TransKernelWT, mu_kernel_wt, phi_kernel, gentime_distr);
    dI_m[1:N[r]-1+W, r] = compute_dI(fitness[r]*Rt[1:N[r]-1, r], b0[r]*dI0[r], 
        TransKernelMut, mu_kernel_mt, phi_kernel, gentime_distr);
    
    // compute sampling-prevalence (useing the same kernel as for generation time)
    I_w[1:N[r], r] = compute_I(dI_w[1:N[r]-1+W, r], TransKernelWT);
    I_m[1:N[r], r] = compute_I(dI_m[1:N[r]-1+W, r], TransKernelMut);
  }
}


model {
  for ( r in 1:R ) {
    // predicted daily cases and mutant freq.
    vector[N[r]] I = I_w[1:N[r], r] + I_m[1:N[r], r];
    vector[N[r]] p = I_m[1:N[r], r] ./ I;
    
    for ( n in 1:NumObs[r] ) {
      // find indices of predicted daily observations we have to aggregate
      int a = ( n == 1 ? 1 : ObsTimes[n-1, r] + 1 );
      int b = ObsTimes[n, r];
      
      // for cases, we can just sum the expected observations for the week
      real Iaggr = sum(I[a:b]);
      // for variant freqs, we take a weighted average of daily freqs
      real paggr = dot_product(p[a:b], I[a:b]) / Iaggr;

      // likelihood cases      
      Cases[n, r] ~ lognormal(log(Iaggr), sigma_cases_pop);
      // likelihood sequences
      VariantSeq[n, r] ~ beta_binomial(
        TotalSeq[n, r], 
        paggr/phi_inv_seq[r], 
        (1-paggr)/phi_inv_seq[r]
      );
    }
  }
  
  // priors
  Rt0 ~ lognormal(0.0, 1.0);
  b0 ~ exponential(10);
  Z_concat ~ normal(0, 1);
  
  // prior for tau: do we want to make this region specific???
  
  //tau ~ lognormal(mu_tau, sd_tau);
  //mu_tau ~ normal(0, 10);
  //sd_tau ~ exponential(0.1);
  tau_pop ~ exponential(1);
  
  //sigma_cases ~ lognormal(sigma_cases_mu, sigma_cases_sd);
  //sigma_cases_mu ~ normal(0, 1);
  //sigma_cases_sd ~ exponential(1);
  sigma_cases_pop ~ exponential(1);
  
  phi_inv_seq ~ exponential(1/phi_inv_seq_mu);
  phi_inv_seq_mu ~ exponential(1e2);
  
  dI0 ~ lognormal(mu_dI0, sd_dI0);
  mu_dI0 ~ normal(0, 10);
  sd_dI0 ~ exponential(0.1);
  
  // hierarchical prior for fitness
  for ( r in 1:R ) {
    fitness[r] ~ lognormal(mu_fitness + Covariates[r,:] * covariate_weights, sd_fitness);
  }
  mu_fitness ~ normal(0, 10);
  sd_fitness ~ exponential(0.1);
  covariate_weights ~ normal(0, 10);
}

generated quantities {
  // sample generation intervals
  int gentime_wt;
  int gentime_mt;
  // quantities for posterior predictive checks (fits)
  matrix[max(N), R] I = I_w + I_m; // daily predicted cases
  matrix[max(NumObs), R] HatCases = rep_matrix(0.0, max(NumObs), R); // aggregated cases
  matrix[max(N), R] pt = rep_matrix(0.0, max(N), R); // daily mutant freq
  matrix[max(NumObs), R] HatVariantFreq = rep_matrix(0.0, max(NumObs), R); // aggregated
  // simulations for HPD
  matrix[max(NumObs), R] SimCases = rep_matrix(0, max(NumObs), R);
  matrix[max(NumObs), R] SimVariantFreq = rep_matrix(0.0, max(NumObs), R);
  // likelihood of freq
  vector[sum(NumObs)] loglik; // concatenated log-likelihoods
  
  // convenient transformations of parameter-of-interest
  vector[R] s = fitness - 1; // selection coefficient
  real s_mean = expm1(mu_fitness + 0.5 * sd_fitness^2);
  
  for ( r in 1:R ) {
    pt[1:N[r], r] = I_m[1:N[r],r] ./ I[1:N[r],r];
    
    for ( i in 1:NumObs[r] ) {
      int a = ( i == 1 ? 1 : ObsTimes[i-1, r]+1 );
      int b = ObsTimes[i, r];
      real Iaggr = sum(I[a:b, r]);
      real paggr = dot_product(pt[a:b, r], I[a:b, r]) / Iaggr;
      HatCases[i, r] = Iaggr;
      HatVariantFreq[i, r] = paggr;
      /* simulation with RNG functions */
      // observed cases
      SimCases[i, r] = lognormal_rng(log(Iaggr), sigma_cases_pop);
      // observed mutant freqs
      if ( TotalSeq[i, r] > 0 ) { // don't divide by zero
        SimVariantFreq[i, r] = beta_binomial_rng(TotalSeq[i, r], 
            paggr/phi_inv_seq[r], (1-paggr)/phi_inv_seq[r]) * inv(TotalSeq[i,r]);
      }
      // log-likelihood of freqs
      loglik[sum(NumObs[:r-1]) + i] = beta_binomial_lpmf(VariantSeq[i, r] | TotalSeq[i, r], 
          paggr/phi_inv_seq[r], (1-paggr)/phi_inv_seq[r]);
    }
  }
  
  // sample generation intervals (for diagnosics)
  if ( gentime_distr == 0 ) { // neg-binom
    gentime_wt = neg_binomial_2_rng(mu_kernel_wt, phi_kernel);
    gentime_mt = neg_binomial_2_rng(mu_kernel_mt, phi_kernel);
  } else if ( gentime_distr == 1) { // gamma
    gentime_wt = gamma_discrete_rng(phi_kernel, phi_kernel / mu_kernel_wt);
    gentime_mt = gamma_discrete_rng(phi_kernel, phi_kernel / mu_kernel_mt);
  } else {
    reject("invalid generation time distribution");
  }
}


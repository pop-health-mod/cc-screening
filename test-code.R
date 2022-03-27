cc_trends_gr <- '
data {
  int<lower = 1> N; //number of observations
  
  int<lower = 1> S; //number of surveys
  int<lower = 1> survey[N];
  
  int<lower = 1> C; //number of countries
  int<lower = 1> country[N];
  
  int<lower = 1> R; //number of regions
  int<lower = 1> region[N];
  
  int<lower = 1> A; //number of age groups - 1
  matrix[A, N] age;
  
  int<lower = 1> G; //number of gni groups - 1
  matrix[G, N] gni;
  
  //vector[N] pgm; //binary variable if national screening program exists
  
  int<lower = 1> n_hiv; // number of hiv groups
  matrix[n_hiv, N] prv_wgt;
  
  //int<lower = 0, upper = 1> whs[N];
  
  vector[N] time;
  matrix[A + 1, N] time_length;
  
  //grouped regression specific inputs
  int<lower = 0> num[N]; // numerator
  int<lower = 1> den[N]; // denominator
}

parameters {
  real alpha; //intercept
  
  vector[C] re_country; 
  real<lower = 0, upper = 10> sd_re_country;
  
  vector[S] re_survey; 
  real<lower = 0, upper = 10> sd_re_survey;
  
  //real beta_time;
  vector[R] rs_time;
  real<lower = 0, upper = 10> sd_rs_region;
  real slope_time; 
  
  //vector[C] rs_pgm; 
  //real<lower = 0, upper = 10> sd_rs_country_pgm;
  //real slope_pgm;
  
  row_vector[A] beta_age;
  row_vector[G] beta_gni;
  
  real<lower = -5, upper = 5> beta_hiv;
  //real beta_recall_hiv;
  vector[C] rs_hiv;
  real<lower = 0, upper = 10> sd_rs_country_hiv;
  
  //real beta_whs;
  //real beta_pgm;
  
  //real beta_time_length;
  row_vector[A + 1] beta_time_length;
}

transformed parameters{
  vector[N] logit_prd;
  vector[N] logit_prob;  
  vector[n_hiv] hiv_effect;
  matrix[N, n_hiv] prd;
  
  hiv_effect[1] = 0;
  hiv_effect[2] = 1;
  //hiv_effect[2] = beta_hiv;
  for(n in 1:N) {
    logit_prd[n] = alpha + re_country[country[n]] 
    + re_survey[survey[n]] 
    + rs_time[region[n]] * time[n] 
    //+ rs_pgm[country[n]] * pgm[n]
    //+ beta_time * time[n] 
    + beta_age * age[, n]
    + beta_gni * gni[, n]
    //+ beta_pgm * pgm[n]
    //+ beta_whs * whs[n]
    + beta_time_length * time_length[, n];
    for (h in 1:n_hiv) {
      prd[n, h] = inv_logit(logit_prd[n] + rs_hiv[country[n]]*hiv_effect[h]);
      //+ beta_recall_hiv * time_length[n] * hiv_effect[h]; 
    }
    // if we have no information on HIV, we do a weighted average
    logit_prob[n] = logit(prd[n, ] * prv_wgt[, n]);
  }
}

model {
  // priors
  
  // Overall Intercept
  alpha ~ normal(0, 10);
  
  // Random Intercepts by survey
  re_survey ~ normal(0, sd_re_survey); 
  sd_re_survey ~ cauchy(0, 2) T[0, 10]; // half-Chauchy prior
  
  // Random Intercepts by country 
  re_country ~ normal(0, sd_re_country); 
  sd_re_country ~ cauchy(0, 2) T[0, 10]; // half-Chauchy prior
  
  // Random Slopes
  rs_time ~ normal(slope_time, sd_rs_region); 
  slope_time ~ normal(0, 10);
  sd_rs_region ~ cauchy(0, 2) T[0, 10];
  
  rs_hiv ~ normal(beta_hiv, sd_rs_country_hiv); 
  beta_hiv ~ normal(0, 10);
  sd_rs_country_hiv ~ cauchy(0, 2) T[0, 10];
  
  //rs_pgm ~ normal(slope_time, sd_rs_country_pgm); 
  //slope_pgm ~ normal(0, 10);
  //sd_rs_country_pgm ~ cauchy(0, 1) T[0, 10];
  
  // Fixed Effects
  beta_age ~ normal(0, 5);
  beta_time_length ~ normal(0, 2);
  //beta_time ~ normal(0, 2);
  //beta_whs ~ normal(0, 5);
  //beta_recall_hiv ~ normal(0, 0.1);
  beta_gni ~ normal(0, 5);
  //beta_hiv ~ normal(0, 5);
  //beta_pgm ~ normal(0, 5);
  
  //likelihood
  num ~ binomial_logit(den, logit_prob);  
}

generated quantities {
  vector[N] log_lik;
  vector[N] pred_num;
  vector[N] pred_prob;
  
  for (n in 1:N) { //can try to vectorize all - didnt work
    //generates log-likelihood for model validation
    log_lik[n] = binomial_logit_lpmf(num[n] | den[n], logit_prob[n]); 
    
    //for posterior predictive checks
    pred_num[n] = binomial_rng(den[n], inv_logit(logit_prob[n]));
    pred_prob[n] = pred_num[n]/den[n];
  }
}
'
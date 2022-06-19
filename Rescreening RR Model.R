#import
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)
library(tidyr)
library(rstan)
library(readxl)
library(matrixStats)
library(haven)
library(survey)
library(stringr)
surveys <- read_excel("surveys.xlsx", sheet = 1)
det_avail <- read_excel("surveys.xlsx", sheet = 3)
eco_det <- read_excel("surveys.xlsx", sheet = 4)
screening_program <- read_excel("surveys.xlsx", sheet = 5)

# stan model
library(rstan)
cc_mod_count_combine <- '
functions {
  // function to run the renewal model for a given entity (e.g., RSS)
  matrix cc_model(int n_obs,
                  vector svy_evr,
                  vector lambda,
                  real beta,
                  //real age_eff,
                  //vector age,
                  real recall_period) {
    
    matrix[n_obs, 2] out;
    
    vector[n_obs] nvr = rep_vector(0.0, n_obs);
    vector[n_obs] pyr = rep_vector(0.0, n_obs);
    vector[n_obs] evr = rep_vector(0.0, n_obs);
    
    vector[n_obs] svy_nvr = rep_vector(0.0, n_obs);
    svy_nvr = 1 - svy_evr; // never tested
    
    for (i in 1:n_obs) {
      evr[i] = svy_evr[i] + lambda[i] * svy_nvr[i];
      pyr[i] = (lambda[i] * svy_nvr[i] + beta * lambda[i] * svy_evr[i]) * recall_period;
    }            
    out[, 1] = evr;
    out[, 2] = pyr;
    return(out);
  }
}
data {
      int<lower = 1> n_count_ecw;   // number of countries
      int<lower = 1> n_obs_ecw;   // number of observations for ecw
      int svy_evr_num_ecw[n_obs_ecw, n_count_ecw];  // survey estimates of % ever screend
      int svy_pyr_num_ecw[n_obs_ecw, n_count_ecw];  // survey estimates of % screend pyr
      int svy_evr_den_ecw[n_obs_ecw, n_count_ecw];  // survey estimates of % ever screend
      int svy_pyr_den_ecw[n_obs_ecw, n_count_ecw];  // survey estimates of % screend pyr
      int ind_obs_ecw[n_obs_ecw];      // indices for observations to compare with data
      matrix[n_obs_ecw, n_count_ecw] svy_ever_ecw;  // distribution of 15 years old
      
      int<lower = 1> n_count_s;   // number of countries
      int<lower = 1> n_obs_s;   // number of observations for s
      int svy_evr_num_s[n_obs_s, n_count_s];  // survey estimates of % ever screend
      int svy_pyr_num_s[n_obs_s, n_count_s];  // survey estimates of % screend pyr
      int svy_evr_den_s[n_obs_s, n_count_s];  // survey estimates of % ever screend
      int svy_pyr_den_s[n_obs_s, n_count_s];  // survey estimates of % screend pyr
      int ind_obs_s[n_obs_s];      // indices for observations to compare with data
      matrix[n_obs_s, n_count_s] svy_ever_s;  // distribution of 15 years old
      
      real recall_period;      // if there is telescoping bias
      //vector[n_obs] age; // vector of age values if the respondants age is younger than 30 or older than 30
    }
    
parameters {
      real<upper = log(0.5)> ln_lambda_ecw[n_obs_ecw, n_count_ecw];
      real<upper = log(0.5)> ln_lambda_s[n_obs_s, n_count_s];
      
      //modelling one beta per country
      real ln_beta_ecw[n_count_ecw];
      real ln_beta_s[n_count_s];
      real ln_beta_reg[2];
      real ln_beta_overall;
      real<lower = 0> sd_beta_reg[2];
      real<lower = 0, upper = 2> sd_beta_overall;
      
      //effect of age
      //real ln_age_eff[n_count];
      //real<lower = 0, upper = 1> sd_age_count;
    }
    
transformed parameters {
      matrix[n_obs_ecw, n_count_ecw] lambda_ecw;
      matrix[n_obs_s, n_count_s] lambda_s;
      vector[n_count_ecw] beta_ecw;
      vector[n_count_s] beta_s;
      vector[2] beta_reg;
      
      lambda_ecw = exp(to_matrix(ln_lambda_ecw));
      lambda_s = exp(to_matrix(ln_lambda_s));
      beta_ecw = exp(to_vector(ln_beta_ecw));
      beta_s = exp(to_vector(ln_beta_s));
      beta_reg = exp(to_vector(ln_beta_reg));
      real beta_overall = exp(ln_beta_overall);
    }

model {
      for(j in 1:n_count_ecw){ 
         matrix[n_obs_ecw, 2] model_prd_ecw = cc_model(n_obs_ecw, svy_ever_ecw[,j], lambda_ecw[,j], beta_ecw[j], recall_period);
         
         for (i in 1:n_obs_ecw) {
          svy_evr_num_ecw[i, j] ~ binomial(svy_evr_den_ecw[i, j], model_prd_ecw[i, 1]); 
          svy_pyr_num_ecw[i, j] ~ binomial(svy_pyr_den_ecw[i, j], model_prd_ecw[i, 2]);
         }
      }
      for(j in 1:n_count_s){
         matrix[n_obs_s, 2] model_prd_s = cc_model(n_obs_s, svy_ever_s[,j], lambda_s[,j], beta_s[j], recall_period);
         
         for (i in 1:n_obs_s) {
          svy_evr_num_s[i, j] ~ binomial(svy_evr_den_s[i, j], model_prd_s[i, 1]); 
          svy_pyr_num_s[i, j] ~ binomial(svy_pyr_den_s[i, j], model_prd_s[i, 2]);
        }
      }
      // prior for re-testing rate (beta)
      ln_beta_ecw ~ normal(ln_beta_reg[1], sd_beta_reg[1]); 
      ln_beta_s ~ normal(ln_beta_reg[2], sd_beta_reg[2]);
      ln_beta_reg ~ normal(ln_beta_overall, sd_beta_overall);
      ln_beta_overall ~ normal(0, 5);
        
      sd_beta_reg ~ normal(0, 15);
      sd_beta_overall ~ normal(0, 5) T[0, 100];
  
      for(j in 1:n_count_ecw){
        ln_lambda_ecw[,j] ~ normal(log(0.005), 5);
      }
      for(j in 1:n_count_s){
        ln_lambda_s[,j] ~ normal(log(0.005), 5);
      }
  }
'
cc_stan_count_combine <- stan_model(model_code = cc_mod_count_combine)

#fit data to model
dat <- readRDS("dat_retest.rds"); nrow(dat)
dat$den_previous <- NA
dat$ever_previous <- NA
dat$surveyid <- paste0(dat$surveyid, dat$Year)
uid <- unique(dat$surveyid)
df <- NULL
for (i in 1:length(uid)) {
  df_i <- subset(dat, surveyid == uid[i])
  df_i$den_previous <- c(NA, df_i$den[-nrow(df_i)])
  df_i$ever_previous <- c(NA, df_i$ever[-nrow(df_i)])
  df <- rbind(df, df_i)
}
dat$den_previous <- df$den_previous
dat$ever_previous <- df$ever_previous

region_country <- unique(surveys[,c("Region", "Country")])
dat$Country <- as.character(dat$Country)
dat <- left_join(dat, region_country)
n_survey <- length(unique(dat$surveyid))
list_surveys <- vector(mode = "list", length = n_survey)

dat$age_eff[dat$age < 10] <- 0
dat$age_eff[dat$age >= 10] <- 1

dat_ecw <- filter(dat, Region != "Southern")
dat_s <- filter(dat, Region == "Southern")
n_survey_ecw <- length(unique(dat_ecw$surveyid))
n_survey_s <- length(unique(dat_s$surveyid))
list_surveys_ecw <- vector(mode = "list", length = n_survey_ecw)
list_surveys_s <- vector(mode = "list", length = n_survey_s)

for(i in 1:n_survey_ecw){
  df <- filter(dat_ecw, surveyid == unique(dat_ecw$surveyid)[i] & age > 18 & age <= 49)
  list_surveys_ecw[[i]] <- df
}

for(i in 1:n_survey_s){
  df <- filter(dat_s, surveyid == unique(dat_s$surveyid)[i] & age > 18 & age <= 29)
  list_surveys_s[[i]] <- df
}

n_obs_ecw <- nrow(list_surveys_ecw[[1]])
ind_obs_ecw <- 1:n_obs_ecw
n_count_ecw <- n_survey_ecw

n_obs_s <- nrow(list_surveys_s[[1]])
ind_obs_s <- 1:n_obs_s
n_count_s <- n_survey_s
recall_period <- 1

svy_evr_num_ecw <- svy_evr_den_ecw <- svy_pyr_num_ecw <- svy_pyr_den_ecw <- svy_ever_ecw <- matrix(0, n_obs_ecw, n_count_ecw)
svy_evr_num_s <- svy_evr_den_s <- svy_pyr_num_s <- svy_pyr_den_s <- svy_ever_s <- matrix(0, n_obs_s, n_count_s)

for(i in 1:n_count_ecw){
  svy_evr_num_ecw[,i] <- list_surveys_ecw[[i]]$num_ever[ind_obs_ecw]
  svy_evr_den_ecw[,i] <- list_surveys_ecw[[i]]$den[ind_obs_ecw]
  svy_pyr_num_ecw[,i] <- list_surveys_ecw[[i]]$num_1y[ind_obs_ecw]
  svy_pyr_den_ecw[,i] <- list_surveys_ecw[[i]]$den[ind_obs_ecw]
  svy_ever_ecw[,i] <- list_surveys_ecw[[i]]$ever_previous
  
  init_val <- c(1 - list_surveys_ecw[[i]]$ever[1], list_surveys_ecw[[i]]$ever[1] - list_surveys_ecw[[i]]$within_1y[1], list_surveys_ecw[[i]]$within_1y[1])
  if (!identical(sum(init_val), 1)) { "stop - initial values inconsistent" }
}

for(i in 1:n_count_s){
  svy_evr_num_s[,i] <- list_surveys_s[[i]]$num_ever[ind_obs_s]
  svy_evr_den_s[,i] <- list_surveys_s[[i]]$den[ind_obs_s]
  svy_pyr_num_s[,i] <- list_surveys_s[[i]]$num_1y[ind_obs_s]
  svy_pyr_den_s[,i] <- list_surveys_s[[i]]$den[ind_obs_s]
  svy_ever_s[,i] <- list_surveys_s[[i]]$ever_previous
  
  init_val <- c(1 - list_surveys_s[[i]]$ever[1], list_surveys_s[[i]]$ever[1] - list_surveys_s[[i]]$within_1y[1], list_surveys_s[[i]]$within_1y[1])
  if (!identical(sum(init_val), 1)) { "stop - initial values inconsistent" }
}

data_stan <- list(ind_obs_ecw = ind_obs_ecw,
                  n_obs_ecw = n_obs_ecw,
                  n_count_ecw = n_count_ecw,
                  svy_evr_num_ecw = matrix(svy_evr_num_ecw, ncol = n_count_ecw),
                  svy_evr_den_ecw = matrix(svy_evr_den_ecw, ncol = n_count_ecw),
                  svy_pyr_num_ecw = matrix(svy_pyr_num_ecw, ncol = n_count_ecw),
                  svy_pyr_den_ecw = matrix(svy_pyr_den_ecw, ncol = n_count_ecw),
                  svy_ever_ecw = svy_ever_ecw,
                  ind_obs_s = ind_obs_s,
                  n_obs_s = n_obs_s,
                  n_count_s = n_count_s,
                  svy_evr_num_s = matrix(svy_evr_num_s, ncol = n_count_s),
                  svy_evr_den_s = matrix(svy_evr_den_s, ncol = n_count_s),
                  svy_pyr_num_s = matrix(svy_pyr_num_s, ncol = n_count_s),
                  svy_pyr_den_s = matrix(svy_pyr_den_s, ncol = n_count_s),
                  svy_ever_s = svy_ever_s,
                  recall_period = recall_period)
rstan_options(auto_write = TRUE)
#' -----------------------------------
# ---- Step 2: fitting the model ----
#' -----------------------------------
# we fit the model 
options(mc.cores = parallel::detectCores())
fit <- sampling(cc_stan_count_combine, data = data_stan, iter = 6000, chains = 4, refresh = 100,
                warmup = 3000, thin = 1, control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20))
#save(list_surveys, data_stan, fit, cc_stan_count_combine, file = "rescreen_1.5")
View(data.frame(rstan::summary(fit, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary))
rstan::summary(fit, pars = c("beta_ecw", "beta_s", "beta_overall", "beta_reg"), probs = c(0.025, 0.5, 0.975))$summary
rstan::stan_trace(fit, pars = c("sd_beta_reg", "sd_beta_overall"))
#fit_beta <- data.frame(rstan::summary(fit, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary)

rstan::summary(fit, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary

rstan::extract(fit_bndhs)$beta
fit_beta <- rbind(fit_beta_bndhs, fit_beta_cpvsteps, fit_beta_ethphia, fit_beta_ghsage, fit_beta_ls14dhs, fit_beta_ls9dhs, fit_beta_mlphia, 
                  fit_beta_rwphia, fit_beta_sabssm, fit_beta_sasage, fit_beta_tzphia, fit_beta_zmphia, fit_beta_zwdhs, fit_beta_zwphia, 
                  fit_beta_s, fit_beta_ecw)
fit_beta$df <- c(s_id, "S Overall", "ECW Overall")
#save(fit_beta, file = "fit_beta")

rstan::stan_trace(fit, pars = c("lambda"))
rstan::stan_trace(fit, pars = c("beta"))
rstan::stan_trace(fit, pars = c("sd_rw"))
rstan::summary(fit, pars = c("lambda"), probs = c(0.025, 0.5, 0.975))$summary
View(data.frame(rstan::summary(fit, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary))
rstan::summary(fit, pars = c("age_eff"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("ln_beta_count", "sd_beta_count"), probs = c(0.025, 0.5, 0.975))$summary

#' -----------------------------------
# ---- Step 3: View data ----
#' -----------------------------------
#beta values for new combined multi-level model
load("rescreen_1")
beta_ecw <- rstan::summary(fit, pars = c("beta_ecw"), probs = c(0.025, 0.5, 0.975))$summary
beta_s <- rstan::summary(fit, pars = c("beta_s"), probs = c(0.025, 0.5, 0.975))$summary
beta_reg <- rstan::summary(fit, pars = c("beta_reg"), probs = c(0.025, 0.5, 0.975))$summary
beta_overall <- rstan::summary(fit, pars = c("beta_overall"), probs = c(0.025, 0.5, 0.975))$summary

id_s <- unique(dat_s$surveyid)
id_ecw <- unique(dat_ecw$surveyid)

overall <- data.frame(id = "Overall", Region = "Overall", screening = paste0(round(beta_overall[,"50%"], 1), " (", round(beta_overall[,"2.5%"], 1), ", ", round(beta_overall[,"97.5%"], 1), ")"))
s <- data.frame(id = id_s, Region = "Southern", screening = paste0(round(beta_s[,"50%"], 1), " (", round(beta_s[,"2.5%"], 1), ", ", round(beta_s[,"97.5%"], 1), ")"))
ecw <- data.frame(id = id_ecw, Region = "Western/Central/Eastern", screening = paste0(round(beta_ecw[,"50%"], 1), " (", round(beta_ecw[,"2.5%"], 1), ", ", round(beta_ecw[,"97.5%"], 1), ")"))
region <- data.frame(id = "Overall", Region = c("Western/Central/Eastern", "Southern"), screening = paste0(round(beta_reg[,"50%"], 1), " (", round(beta_reg[,"2.5%"], 1), ", ", round(beta_reg[,"97.5%"], 1), ")"))

df <- rbind(overall, 
            region[region$Region == "Western/Central/Eastern",], ecw,
            region[region$Region == "Southern",], s)
View(df)


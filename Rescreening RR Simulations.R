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
library(ggplot2)
library(gridExtra)
library(ggsci)
surveys <- read_excel("surveys.xlsx", sheet = 1)
det_avail <- read_excel("surveys.xlsx", sheet = 3)
eco_det <- read_excel("surveys.xlsx", sheet = 4)
screening_program <- read_excel("surveys.xlsx", sheet = 5)

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
      real<lower = 0> sd_beta_overall;
      
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
        
      sd_beta_reg ~ normal(0, 1);
      sd_beta_overall ~ normal(0, 3);
  
      for(j in 1:n_count_ecw){
        ln_lambda_ecw[,j] ~ normal(log(0.005), 5);
      }
      for(j in 1:n_count_s){
        ln_lambda_s[,j] ~ normal(log(0.005), 5);
      }
  }
generated quantities {
    matrix[n_obs_ecw, n_count_ecw] prd_pyr_num_ecw;
    matrix[n_obs_ecw, n_count_ecw] prd_evr_num_ecw;   
    matrix[n_obs_ecw, n_count_ecw] prd_evr_ecw;  
    matrix[n_obs_ecw, n_count_ecw] prd_pyr_ecw;
    
    matrix[n_obs_s, n_count_s] prd_pyr_num_s;
    matrix[n_obs_s, n_count_s] prd_evr_num_s;   
    matrix[n_obs_s, n_count_s] prd_evr_s;  
    matrix[n_obs_s, n_count_s] prd_pyr_s;
    
    for(j in 1:n_count_ecw){
      matrix[n_obs_ecw, 2] local_ecw = cc_model(n_obs_ecw, svy_ever_ecw[,j], lambda_ecw[,j], beta_ecw[j], recall_period);
    
      for (i in 1:n_obs_ecw) {
        // we predict the observations
        prd_evr_num_ecw[i, j] = binomial_rng(svy_evr_den_ecw[i, j], local_ecw[ind_obs_ecw[i], 1]);   
        prd_pyr_num_ecw[i, j] = binomial_rng(svy_pyr_den_ecw[i, j], local_ecw[ind_obs_ecw[i], 2]);
        prd_evr_ecw[i, j] = prd_evr_num_ecw[i, j] / svy_evr_den_ecw[i, j];
        prd_pyr_ecw[i, j] = prd_pyr_num_ecw[i, j] / svy_pyr_den_ecw[i, j];
      }
    }
    
    for(j in 1:n_count_s){
      matrix[n_obs_s, 2] local_s = cc_model(n_obs_s, svy_ever_s[,j], lambda_s[,j], beta_s[j], recall_period);
    
      for (i in 1:n_obs_s) {
        // we predict the observations
        prd_evr_num_s[i, j] = binomial_rng(svy_evr_den_s[i, j], local_s[ind_obs_s[i], 1]);   
        prd_pyr_num_s[i, j] = binomial_rng(svy_pyr_den_s[i, j], local_s[ind_obs_s[i], 2]);
        prd_evr_s[i, j] = prd_evr_num_s[i, j] / svy_evr_den_s[i, j];
        prd_pyr_s[i, j] = prd_pyr_num_s[i, j] / svy_pyr_den_s[i, j];
      }
    }
  }  
'
cc_stan_count_combine <- stan_model(model_code = cc_mod_count_combine)

n_age_ecw <- 34
n_age_s <- 34
lambda_version <- 1
if (lambda_version == 1) {
  lambda_ecw0 <- exp(log(0.001) + c(1:n_age_ecw) * 0.22 - (c(1:n_age_ecw)^2) * 0.006)
  lambda_ecw1 <- exp(log(0.001) + c(1:n_age_ecw) * 0.22 - (c(1:n_age_ecw)^2) * 0.006) 
  # lambda_e0 <- exp(log(0.001) + c(1:n_age_e) * 0.22 - (c(1:n_age_e)^2) * 0.006)
  # lambda_e1 <- exp(log(0.001) + c(1:n_age_e) * 0.22 - (c(1:n_age_e)^2) * 0.006) 
  lambda_s0 <- exp(log(0.006) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006)
  lambda_s1 <- exp(log(0.006) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006) }
if (lambda_version == 2) {
  lambda_ecw0 <- exp(log(0.001) + c(1:n_age_ecw) * 0.22 - (c(1:n_age_ecw)^2) * 0.006)
  lambda_ecw1 <- exp(log(0.00115) + c(1:n_age_ecw) * 0.22 - (c(1:n_age_ecw)^2) * 0.006) 
  # lambda_e0 <- exp(log(0.001) + c(1:n_age_e) * 0.22 - (c(1:n_age_e)^2) * 0.006)
  # lambda_e1 <- exp(log(0.00115) + c(1:n_age_e) * 0.22 - (c(1:n_age_e)^2) * 0.006) 
  lambda_s0 <- exp(log(0.0055) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006)
  lambda_s1 <- exp(log(0.0065) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006)}
if (lambda_version == 3) {
  lambda_ecw0 <- exp(log(0.0005) + c(1:n_age_ecw) * 0.2 - (c(1:n_age_ecw)^2) * 0.006)
  lambda_ecw1 <- exp(log(0.0015) + c(1:n_age_ecw) * 0.2 - (c(1:n_age_ecw)^2) * 0.006)
  # lambda_e0 <- exp(log(0.001) + c(1:n_age_e) * 0.19 - (c(1:n_age_e)^2) * 0.006)
  # lambda_e1 <- exp(log(0.003) + c(1:n_age_e) * 0.12 - (c(1:n_age_e)^2) * 0.003)
  lambda_s0 <- exp(log(0.003) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006)
  lambda_s1 <- exp(log(0.009) + c(1:n_age_s) * 0.2 - (c(1:n_age_s)^2) * 0.006) }
  
beta_ecw <- c(30, 40, 50)
# beta_e <- c(15, 20, 25)
beta_s <- c(1, 5, 10)
recall_period <- 1
dt <- 1
n_yr <- 30
time <- seq(0, n_yr, by = dt)
niter <- (n_yr) / dt + 1

lambda_ecw <- array(data = NA, dim = c(niter, n_age_ecw))
# lambda_e <- array(data = NA, dim = c(niter, n_age_e))
lambda_s <- array(data = NA, dim = c(niter, n_age_s))
for (a in 1:n_age_ecw) { 
  lambda_ecw[, a] <- seq(lambda_ecw0[a], lambda_ecw1[a], length.out = niter)
}
for (a in 1:n_age_s) {
  lambda_s[, a] <- seq(lambda_s0[a], lambda_s1[a], length.out = niter)
}

a <- n_yr - 5
lambda_ecw[a:(a+2),] <- lambda_ecw[a:(a+2),]*2
lambda_s[a:(a+2),] <- lambda_s[a:(a+2),]*2

library(ggplot2)
time_yr <- time[time %in% seq(1, n_yr, 1)]
lambda_yr <- lambda_ecw[time %in% seq(1, n_yr, 1), ]
df <- expand.grid(year = time_yr, age_n = 1:n_age_ecw)
df$age <- df$age_n + 15
df$lambda <- NA
for (i in 1:nrow(df)) { df$lambda[i] <- lambda_ecw[which(time == df$year[i]), df$age_n[i]] }
pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
heat <- ggplot(df, aes(year, age, fill = lambda)) + 
  geom_tile() +
  scale_fill_gradientn(name = "Initial Screening Rate", colours = pal, breaks = c(0.002, 0.004, 0.006, 0.008)) +
  scale_y_continuous(name = "Age", breaks = seq(15, 50, 5)) +
  scale_x_continuous(name = "Years") +
  theme_minimal() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5),
        strip.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20));heat

X_ecw <- array(data = 0, dim = c(niter, n_age_ecw, 3, length(beta_ecw)), 
           dimnames = list(c(1:niter), c(1:n_age_ecw + 14), 
                           c("nvr", "evr", "pyr")))  
X_ecw[1, , "nvr",] <- 1

X_s <- array(data = 0, dim = c(niter, n_age_s, 3, length(beta_s)), 
               dimnames = list(c(1:niter), c(1:n_age_s + 14), 
                               c("nvr", "evr", "pyr"))) 
X_s[1, , "nvr",] <- 1

X <- list(X_ecw, X_s)
lambda <- list(lambda_ecw, lambda_s)
beta <- list(beta_ecw, beta_s)
n_age <- list(n_age_ecw, n_age_s)

for(r in 1:length(X)){
  for (i in 2:niter) {
    for (a in 1:n_age[[r]]) {
      X[[r]][i, a, "nvr",] <- X[[r]][i - 1, a, "nvr",] + dt * (- lambda[[r]][i, a] * X[[r]][i - 1, a, "nvr",])   
      X[[r]][i, a, "evr",] <- X[[r]][i - 1, a, "evr",] + dt * (+ lambda[[r]][i, a] * X[[r]][i - 1, a, "nvr",])  
      X[[r]][i, a, "pyr",] <- X[[r]][i - 1, a, "pyr",] + dt * (+ lambda[[r]][i, a] * X[[r]][i - 1, a, "nvr",] 
                                                     + beta[[r]] * lambda[[r]][i, a] * X[[r]][i - 1, a, "evr",] 
                                                     - 1 / recall_period * X[[r]][i - 1, a, "pyr",])  
    }
    # we age the population using a RAS
    if ((i * dt) %in% c(1:niter)) {
      PX <- X[[r]][i, , ,]
      X[[r]][i, 1, "nvr",] <- 1
      X[[r]][i, 1, "evr",] <- 0
      X[[r]][i, 1, "pyr",] <- 0
      for (a in 2:n_age[[r]]) {
        X[[r]][i, a, "nvr",] <- PX[a - 1, "nvr",]
        X[[r]][i, a, "evr",] <- PX[a - 1, "evr",]
        X[[r]][i, a, "pyr",] <- PX[a - 1, "pyr",]
      }
    }
  }
}

X_ras <- X
list_sim_country_ecw <- NULL
# list_sim_country_e <- NULL
list_sim_country_s <- NULL
for(r in 1:length(X)){
  for(j in 1:length(beta[[r]])){
    sim <- as.data.frame(X[[r]][niter, , ,j])
    sim$den <- round(runif(nrow(sim), min = 10000, max = 10000), 0)
    sim$num_evr <- rbinom(nrow(sim), size = sim$den, prob = sim$evr)
    #  the num for pyr should be lower than the num for ever, we sampled sequentially.
    prb_pyr <- ifelse(sim$pyr == 0, 0, sim$pyr / sim$evr)
    sim$num_pyr <- rbinom(nrow(sim), size = sim$num_evr, prob = prb_pyr)
    sim$est_evr <- sim$num_evr / sim$den
    sim$est_pyr <- sim$num_pyr / sim$den
    sim$ever_previous <- c(NA, sim$est_evr[1:(nrow(sim) - 1)])
    sim$age <- seq(15, n_age[[r]] + 14, by = 1)
    sim$lci_evr <- sim$uci_evr <- NA
    sim$lci_pyr <- sim$uci_pyr <- NA
    for (i in 1:nrow(sim)) {
      res <- binom.test(sim$num_evr[i], sim$den[i])$conf.int
      sim$lci_evr[i] <- res[1]
      sim$uci_evr[i] <- res[2]
      res <- binom.test(sim$num_pyr[i], sim$den[i])$conf.int
      sim$lci_pyr[i] <- res[1]
      sim$uci_pyr[i] <- res[2]
    }
    if(r == 1){
      list_sim_country_ecw[[j]] <- sim
    # } else if(r == 2){
    #   list_sim_country_e[[j]] <- sim
    } else{
      list_sim_country_s[[j]] <- sim
    }
  }
}

#figures for screening in past year and lifetime
par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
age <- c(15:48) + 0.1
plot(list_sim_country_ecw[[1]]$est_evr ~ age, type = "p", pch = 19, ylim = c(0, 1), col = "royalblue3",
     yaxt = "n", xlab = "Age", ylab = "Proportion Screened in Lifetime", lwd = 1.5, axes = T, cex.axis = 1.5, cex.lab = 1.5)
points(list_sim_country_s[[1]]$est_evr ~ age, pch = 19, col = "violetred3")
lines(list_sim_country_ecw[[1]]$evr ~ age, col = "royalblue3")
lines(list_sim_country_s[[1]]$evr ~ age, col = "violetred3")
axis(2, at = seq(0, 1, by = 0.1), cex.axis = 1.5)
legend("topleft", legend = c("True", "Region 1", "Region 2"),
       pch = c(NA, rep(21, 7)), lwd = c(1.5, rep(NA, 7)), 
       col = c("cornsilk4", "royalblue3", "violetred3"), cex = 1.5,
       pt.bg = c(NA, "royalblue3", "violetred3"), bty = "n")

plot(list_sim_country_ecw[[1]]$pyr ~ age, type = "l", pch = 19, ylim = c(0, 0.3), col = "royalblue4", 
     yaxt = "n", xlab = "Age", ylab = "Proportion Screened in the Past Year", lwd = 1.5, axes = T, cex.axis = 1.5, cex.lab = 1.5)
lines(list_sim_country_ecw[[2]]$pyr ~ age, col = "steelblue2", lwd = 1.5)
lines(list_sim_country_ecw[[3]]$pyr ~ age, col = "royalblue1", lwd = 1.5)
lines(list_sim_country_s[[1]]$pyr ~ age, col = "violetred4", lwd = 1.5)
lines(list_sim_country_s[[2]]$pyr ~ age, col = "tomato3", lwd = 1.5)
lines(list_sim_country_s[[3]]$pyr ~ age, col = "violetred1", lwd = 1.5)
points(list_sim_country_ecw[[1]]$est_pyr ~ list_sim_country_ecw[[1]]$age, pch = 19 , col = "royalblue4")
points(list_sim_country_ecw[[2]]$est_pyr ~ list_sim_country_ecw[[2]]$age, pch = 19 , col = "steelblue2")
points(list_sim_country_ecw[[3]]$est_pyr ~ list_sim_country_ecw[[3]]$age, pch = 19 , col = "royalblue1")
points(list_sim_country_s[[1]]$est_pyr ~ list_sim_country_s[[1]]$age, pch = 19 , col = "violetred4")
points(list_sim_country_s[[2]]$est_pyr ~ list_sim_country_s[[2]]$age, pch = 19 , col = "tomato3")
points(list_sim_country_s[[3]]$est_pyr ~ list_sim_country_s[[3]]$age, pch = 19 , col = "violetred1")
axis(2, at = seq(0, 0.3, by = 0.05), cex.axis = 1.5)
legend("topleft", legend = c("True", paste("Region 1: Phi", beta_ecw[1]), paste("Region 1: Phi", beta_ecw[2]), paste("Region 1: Phi", beta_ecw[3]), paste("Region 2: Phi", beta_s[1]), paste("Region 2: Phi", beta_s[2]), paste("Region 2: Phi", beta_s[3])),
       pch = c(NA, rep(21, 7)), lwd = c(1.5, rep(NA, 7)), 
       col = c("cornsilk4", "royalblue4", "steelblue2", "royalblue1", "violetred4", "tomato3", "violetred1"), cex = 1.5,
       pt.bg = c(NA, "royalblue4", "steelblue2", "royalblue1", "violetred4", "tomato3", "violetred1"), bty = "n")

ind_obs_ecw <- which(!is.na(list_sim_country_ecw[[1]]$ever_previous)) - 1
n_obs_ecw <- length(ind_obs_ecw)
n_count_ecw <- length(beta_ecw)
svy_evr_num_ecw <- svy_evr_den_ecw <- svy_pyr_num_ecw <- svy_pyr_den_ecw <- svy_ever_ecw <- matrix(0, n_obs_ecw, n_count_ecw)

ind_obs_s <- which(!is.na(list_sim_country_s[[1]]$ever_previous)) - 1
n_obs_s <- length(ind_obs_s)
n_count_s <- length(beta_s)
svy_evr_num_s <- svy_evr_den_s <- svy_pyr_num_s <- svy_pyr_den_s <- svy_ever_s <- matrix(0, n_obs_s, n_count_s)

for(i in 1:n_count_ecw){
  list_sim_country_ecw[[i]] <- na.omit(list_sim_country_ecw[[i]])
  svy_evr_num_ecw[,i] <- list_sim_country_ecw[[i]]$num_evr[ind_obs_ecw]
  svy_evr_den_ecw[,i] <- list_sim_country_ecw[[i]]$den[ind_obs_ecw]
  svy_pyr_num_ecw[,i] <- list_sim_country_ecw[[i]]$num_pyr[ind_obs_ecw]
  svy_pyr_den_ecw[,i] <- list_sim_country_ecw[[i]]$den[ind_obs_ecw]
  svy_ever_ecw[,i] <- list_sim_country_ecw[[i]]$ever_previous
  
  init_val <- c(1 - list_sim_country_ecw[[i]]$num_evr[1], list_sim_country_ecw[[i]]$num_evr[1] - list_sim_country_ecw[[i]]$est_pyr[1], list_sim_country_ecw[[i]]$est_pyr[1])
  if (!identical(sum(init_val), 1)) { "stop - initial values inconsistent" }
}

for(i in 1:n_count_s){
  list_sim_country_s[[i]] <- na.omit(list_sim_country_s[[i]])
  svy_evr_num_s[,i] <- list_sim_country_s[[i]]$num_evr[ind_obs_s]
  svy_evr_den_s[,i] <- list_sim_country_s[[i]]$den[ind_obs_s]
  svy_pyr_num_s[,i] <- list_sim_country_s[[i]]$num_pyr[ind_obs_s]
  svy_pyr_den_s[,i] <- list_sim_country_s[[i]]$den[ind_obs_s]
  svy_ever_s[,i] <- list_sim_country_s[[i]]$ever_previous
  
  init_val <- c(1 - list_sim_country_s[[i]]$num_evr[1], list_sim_country_s[[i]]$num_evr[1] - list_sim_country_s[[i]]$est_pyr[1], list_sim_country_s[[i]]$est_pyr[1])
  if (!identical(sum(init_val), 1)) { "stop - initial values inconsistent" }
}

data_stan <- list(ind_obs_ecw = ind_obs_ecw,
                  n_obs_ecw = n_obs_ecw,
                  n_count_ecw = n_count_ecw,
                  svy_evr_num_ecw = svy_evr_num_ecw,
                  svy_evr_den_ecw = svy_evr_den_ecw,
                  svy_pyr_num_ecw = svy_pyr_num_ecw,
                  svy_pyr_den_ecw = svy_pyr_den_ecw,
                  svy_ever_ecw = svy_ever_ecw,
                  ind_obs_s = ind_obs_s,
                  n_obs_s = n_obs_s,
                  n_count_s = n_count_s,
                  svy_evr_num_s = svy_evr_num_s,
                  svy_evr_den_s = svy_evr_den_s,
                  svy_pyr_num_s = svy_pyr_num_s,
                  svy_pyr_den_s = svy_pyr_den_s,
                  svy_ever_s = svy_ever_s,
                  #age = age,
                  recall_period = recall_period)

rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())
fit <- sampling(cc_stan_count_combine, data = data_stan, iter = 7000, chains = 4, refresh = 100,
                warmup = 3500, thin = 1, control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15))
#save(X, lambda, beta, n_age, data_stan, fit, cc_stan_count_combine, file = "retest_sim_three_reg")
rstan::summary(fit, pars = c("beta_ecw", "beta_s", "beta_overall", "beta_reg"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("sd_beta_overall", "sd_beta_reg"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, probs = c(0.025, 0.5, 0.975))$summary
rstan::stan_trace(fit, pars = c("beta_ecw", "beta_s", "beta_overall", "beta_reg"))
rstan::stan_trace(fit, pars = c("sd_beta_overall", "sd_beta_reg"))
View(rstan::summary(fit, probs = c(0.025, 0.5, 0.975))$summary)
rstan::summary(fit, pars = c("beta", "lambda", "sd_rw"), probs = c(0.025, 0.5, 0.975))$summary

draws <- rstan::extract(fit)
prd_evr_ecw <- prd_pyr_ecw <- lambda_est_ecw <- rate_ecw <- NULL
prd_evr_s <- prd_pyr_s <- lambda_est_s <- rate_s <- NULL
rate_reg <- NULL
re_test_ecw <- as.data.frame(matrix(data = NA, nrow = n_count_ecw, ncol = 4))
re_test_s <- as.data.frame(matrix(data = NA, nrow = n_count_s, ncol = 4))
re_test_reg <- as.data.frame(matrix(data = NA, nrow = 2, ncol = 4))

for(i in 1:n_count_ecw){
  prd_evr_ecw[[i]] <- as.data.frame(colQuantiles(draws$prd_evr_ecw[,,i], probs = c(0.025, 0.5, 0.975)))
  prd_evr_ecw[[i]]$age <- c(min(list_sim_country_ecw[[i]]$age ):max(list_sim_country_ecw[[i]]$age)) + 0.25
  prd_evr_ecw[[i]]$Country <- i
  prd_pyr_ecw[[i]] <- as.data.frame(colQuantiles(draws$prd_pyr_ecw[,,i], probs = c(0.025, 0.5, 0.975)))
  prd_pyr_ecw[[i]]$age <- c(min(list_sim_country_ecw[[i]]$age ):max(list_sim_country_ecw[[i]]$age)) + 0.25
  prd_pyr_ecw[[i]]$Country <- i
  
  lambda_est_ecw[[i]] <- as.data.frame(colQuantiles(draws$lambda_ecw[,,i], probs = c(0.025, 0.5, 0.975)))
  lambda_est_ecw[[i]]$Country <- i
  
  re_test_ecw[i,] <- round(data.frame(t(quantile(draws$beta_ecw[,i], probs = c(0.025, 0.5, 0.975))), Country = i), 2)
  rate_ecw[i] <- paste0(re_test_ecw[i,2], " (", re_test_ecw[i,1], "-", re_test_ecw[i,3], ")")
}

for(i in 1:n_count_s){
  prd_evr_s[[i]] <- as.data.frame(colQuantiles(draws$prd_evr_s[,,i], probs = c(0.025, 0.5, 0.975)))
  prd_evr_s[[i]]$age <- c(min(list_sim_country_s[[i]]$age ):max(list_sim_country_s[[i]]$age)) + 0.25
  prd_evr_s[[i]]$Country <- i
  prd_pyr_s[[i]] <- as.data.frame(colQuantiles(draws$prd_pyr_s[,,i], probs = c(0.025, 0.5, 0.975)))
  prd_pyr_s[[i]]$age <- c(min(list_sim_country_s[[i]]$age ):max(list_sim_country_s[[i]]$age)) + 0.25
  prd_pyr_s[[i]]$Country <- i
  
  lambda_est_s[[i]] <- as.data.frame(colQuantiles(draws$lambda_s[,,i], probs = c(0.025, 0.5, 0.975)))
  lambda_est_s[[i]]$Country <- i
  
  re_test_s[i,] <- round(data.frame(t(quantile(draws$beta_s[,i], probs = c(0.025, 0.5, 0.975))), Country = i), 2)
  rate_s[i] <- paste0(re_test_s[i,2], " (", re_test_s[i,1], "-", re_test_s[i,3], ")")
}

for(i in 1:2){
  re_test_reg[i,] <- round(data.frame(t(quantile(draws$beta_reg[,i], probs = c(0.025, 0.5, 0.975))), Country = i), 2)
  rate_reg[i] <- paste0(re_test_reg[i,2], " (", re_test_reg[i,1], "-", re_test_reg[i,3], ")")
}

df_prd_evr_ecw <- prd_evr_ecw[[1]]
df_prd_pyr_ecw <- prd_pyr_ecw[[1]]
df_lambda_est_ecw <- lambda_est_ecw[[1]]
df_prd_evr_s <- prd_evr_s[[1]]
df_prd_pyr_s <- prd_pyr_s[[1]]
df_lambda_est_s <- lambda_est_s[[1]]
for(i in 2:n_count_ecw){
  df_prd_evr_ecw <- rbind(df_prd_evr_ecw, prd_evr_ecw[[i]])
  df_prd_pyr_ecw <- rbind(df_prd_pyr_ecw, prd_pyr_ecw[[i]])
  df_lambda_est_ecw <- rbind(df_lambda_est_ecw, lambda_est_ecw[[i]])
}
for(i in 2:n_count_s){
  df_prd_evr_s <- rbind(df_prd_evr_s, prd_evr_s[[i]])
  df_prd_pyr_s <- rbind(df_prd_pyr_s, prd_pyr_s[[i]])
  df_lambda_est_s <- rbind(df_lambda_est_s, lambda_est_s[[i]])
}

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
x_label <- seq(15, 30, 5)
age1 <- list_sim_country_ecw[[1]]$age 
age2 <- list_sim_country_ecw[[1]]$age + 0.1
age3 <- list_sim_country_ecw[[1]]$age + 0.2
age4 <- list_sim_country_ecw[[1]]$age + 0.3
age5 <- list_sim_country_ecw[[1]]$age + 0.4
age6 <- list_sim_country_ecw[[1]]$age + 0.5

plot(lambda_ecw[niter - 1 / dt, ] ~ I(c(16:49)), type = "l", ylim = c(0, 0.06), yaxt = "n", 
     xlab = "Age", ylab = "Initial Screening Rate", col = "lightblue", lwd = 3, axes = T, 
     cex.axis = 1.5, cex.lab = 1.5)
lines(lambda_s[niter - 1 / dt, ] ~ I(c(16:49)), col = "lightpink", lwd = 3)
points(lambda_est_ecw[[1]]$`50%` ~ age1, pch = 19, col = "royalblue4", cex = 0.5)
points(lambda_est_ecw[[2]]$`50%` ~ age2, pch = 19, col = "steelblue3", cex = 0.5)
points(lambda_est_ecw[[3]]$`50%` ~ age3, pch = 19, col = "royalblue1", cex = 0.5)
points(lambda_est_s[[1]]$`50%` ~ age4, pch = 19, col = "violetred4", cex = 0.5)
points(lambda_est_s[[2]]$`50%` ~ age5, pch = 19, col = "tomato3", cex = 0.5)
points(lambda_est_s[[3]]$`50%` ~ age6, pch = 19, col = "violetred1", cex = 0.5)
axis(2, at = seq(0, 0.06, by = 0.01), cex.axis = 1.5)
segments(x0 = age1, x1 = age1, 
         y0 = lambda_est_ecw[[1]]$`2.5%`, y1 = lambda_est_ecw[[1]]$`97.5%`, col = "royalblue4")
segments(x0 = age2, x1 = age2, 
         y0 = lambda_est_ecw[[2]]$`2.5%`, y1 = lambda_est_ecw[[2]]$`97.5%`, col = "steelblue3")
segments(x0 = age3, x1 = age3, 
         y0 = lambda_est_ecw[[3]]$`2.5%`, y1 = lambda_est_ecw[[3]]$`97.5%`, col = "royalblue4")
segments(x0 = age4, x1 = age4, 
         y0 = lambda_est_s[[1]]$`2.5%`, y1 = lambda_est_s[[1]]$`97.5%`, col = "violetred4")
segments(x0 = age5, x1 = age5, 
         y0 = lambda_est_s[[2]]$`2.5%`, y1 = lambda_est_s[[2]]$`97.5%`, col = "tomato3")
segments(x0 = age6, x1 = age6, 
         y0 = lambda_est_s[[3]]$`2.5%`, y1 = lambda_est_s[[3]]$`97.5%`, col = "violetred1")
legend("topleft", legend = c("True", paste("Phi", beta_ecw[1]), paste("Phi", beta_ecw[2]), paste("Phi", beta_ecw[3]), paste("Phi", beta_s[1]), paste("Phi", beta_s[2]), paste("Phi", beta_s[3])),
       pch = c(NA, rep(21, 7)), lwd = c(1.5, rep(NA, 7)), 
       col = c("cornsilk4", "royalblue4", "steelblue2", "royalblue1", "violetred4", "tomato3", "violetred1"), cex = 1.5,
       pt.bg = c(NA, "royalblue4", "steelblue2", "royalblue1", "violetred4", "tomato3", "violetred1"), bty = "n")

sim_beta_ecw <- re_test_ecw[,2]
sim_beta_ecw_l <- re_test_ecw[,1]
sim_beta_ecw_u <- re_test_ecw[,3]
df_beta_ecw <- data.frame(Country = rep(c(paste("Phi", beta_ecw[1]), paste("Phi", beta_ecw[2]), paste("Phi", beta_ecw[3])), 2),
                      name = c(rep("True", 3), rep("Modelled", 3)), 
                      val = c(beta_ecw, sim_beta_ecw),
                      lower = c(rep(NA, 3), sim_beta_ecw_l),
                      upper = c(rep(NA, 3), sim_beta_ecw_u))
df_beta_ecw$Country = factor(df_beta_ecw$Country, levels = c(paste("Phi", beta_ecw[1]), paste("Phi", beta_ecw[2]), paste("Phi", beta_ecw[3])))
df_beta_ecw$Region <- "Region 1"

sim_beta_s <- re_test_s[,2]
sim_beta_s_l <- re_test_s[,1]
sim_beta_s_u <- re_test_s[,3]
df_beta_s <- data.frame(Country = rep(c(paste("Phi", beta_s[1]), paste("Phi", beta_s[2]), paste("Phi", beta_s[3])), 2),
                          name = c(rep("True", 3), rep("Modelled", 3)), 
                          val = c(beta_s, sim_beta_s),
                          lower = c(rep(NA, 3), sim_beta_s_l),
                          upper = c(rep(NA, 3), sim_beta_s_u))
df_beta_s$Country = factor(df_beta_s$Country, levels = c(paste("Phi", beta_s[1]), paste("Phi", beta_s[2]), paste("Phi", beta_s[3])))
df_beta_s$Region <- "Region 2"

df_beta <- rbind(df_beta_s, df_beta_ecw)

ggplot(df_beta, aes(x = Country, y = val, fill = name)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.5), width = 0.4) +
  ylab("Rate Ratio")  +
  ylim(0, 85) +
  #scale_y_continuous(breaks = seq(0, 90, by = 20)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        strip.text = element_text(size = 15)) +
  scale_fill_nejm(name = "Type") +
  facet_grid(~Region, space = "free", scales = "free")

beta_ecw; beta_s
rate_ecw; rate_s



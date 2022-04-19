#Life table function 
# ---- Import ----
library(tidyr)
library(dplyr)
library(knitr)
library(kableExtra)
library(rstan)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(GGally)
library(bayesplot)
library(rstanarm)
library(splines)
library(lme4)
library(jtools)
library(bayesplot)
library(readxl)
library(matrixStats)
library(haven)
library(survey)
library(stringr)
library(stats)
library(diagis)
library(loo); options(mc.cores = 6)
library(ggh4x)
surveys <- read_excel("surveys.xlsx", sheet = 1)
det_avail <- read_excel("surveys.xlsx", sheet = 3)
eco_det <- read_excel("surveys.xlsx", sheet = 4)
screening_program <- read_excel("surveys.xlsx", sheet = 5)
region_country <- read_excel("surveys.xlsx", sheet = 7)
load("Cleaned Pooled Surveys")

source("predict_poststrat_func.R")  #will need to update these functions based on the covariates I add
source("life_tbl_func.R")

# ----  I) Predictions - to use for life table ----
##find weights by country
df_weights_hiv <- weights(agegr = c("25-29", "30-34", "35-39", "40-44", "45-49"), hiv = T, years = c(2000:2020))

##multiply predictions by country weights
load("gr_combined70")
draws <- rstan::extract(fit)
#df_stan_org <- df_stan
df_stan$predicted <- "Y"
region_country$Region[region_country$Region == "Central" | region_country$Region == "Western"] <- "Central/Western"
df_stan_all_count <- left_join(region_country, df_stan)
df_stan_all_count$predicted[is.na(df_stan_all_count$predicted)] <- "N"
source("predict_poststrat_func.R")

pred_prob_weighted_3049 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(5:20), agegr = 2:5, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)

pred_prob_weighted_2529 <- pred_weighted(draws, df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(5), agegr = 1, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)

pred_prob_weighted_4549 <- pred_weighted(draws, df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(20), agegr = 5, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)

pred_prob_weighted_4044 <- pred_weighted(draws, df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(19), agegr = 4, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)

prob_by_count_age_year_3y <- prob_by_var(pred_prob_weighted_3049[[1]], pred_prob_weighted_3049[[3]], 3, "Country", "agegr", "Year")
prob_by_count_age_year_ever <- prob_by_var(pred_prob_weighted_3049[[2]], pred_prob_weighted_3049[[3]], 3, "Country", "agegr", "Year")
prob_by_count_year_ever_2529 <- prob_by_var(pred_prob_weighted_2529[[2]], pred_prob_weighted_2529[[3]], 2, "Country", "Year")
prob_by_count_year_ever_4044 <- prob_by_var(pred_prob_weighted_4044[[2]], pred_prob_weighted_4044[[3]], 2, "Country", "Year")
prob_by_count_year_ever_4549 <- prob_by_var(pred_prob_weighted_4549[[2]], pred_prob_weighted_4549[[3]], 2, "Country", "Year")

prob_by_hiv_count_age_year_3y <- prob_by_var(pred_prob_weighted_3049[[1]], pred_prob_weighted_3049[[3]], 4, "Country", "agegr", "Year", "hivstat")
prob_by_hiv_count_age_year_ever <- prob_by_var(pred_prob_weighted_3049[[2]], pred_prob_weighted_3049[[3]], 4, "Country", "agegr", "Year", "hivstat")
prob_by_hiv_count_year_ever_2529 <- prob_by_var(pred_prob_weighted_2529[[2]], pred_prob_weighted_2529[[3]], 3, "Country", "Year", "hivstat")
prob_by_hiv_count_year_ever_4549 <- prob_by_var(pred_prob_weighted_4549[[2]], pred_prob_weighted_4549[[3]], 3, "Country", "Year", "hivstat")

prob_hiv_ever_4549 <- prob_final(pred_prob_weighted_4549[[2]], dat_pred = pred_prob_weighted_4549[[3]], 2, "Year", "hivstat")
prob_ever_4549 <- prob_final(pred_prob_weighted_4549[[2]], dat_pred = pred_prob_weighted_4549[[3]], 1, "Year")
View(prob_ever_4549[, c("2.5%", "50%", "97.5%")])

View(prob_by_count_year_ever_4044[, c("Country", "2.5%", "50%", "97.5%")])

save(pred_prob_weighted_3049, 
     pred_prob_weighted_2529, pred_prob_weighted_4044, #pred_prob_weighted_4549,
     prob_by_count_age_year_3y,
     prob_by_count_age_year_ever,
     prob_by_count_year_ever_2529,
     prob_by_count_year_ever_4044,
     #prob_by_count_year_ever_4549,
     # prob_by_hiv_count_age_year_3y,
     # prob_by_hiv_count_age_year_ever,
     # prob_by_hiv_count_year_ever_2529,
     # prob_by_hiv_count_year_ever_4549,
     file = "pred_prob_for_lt_70")

#load("pred_prob_for_lt_70")
# ----  J) Life table values ----
# load("retest_combined_1949_17_ecw")
# fit_ecw <- fit
# load("retest_combined_1929_16_s")
# fit_s <- fit
# load("retest_combined_1949_one_beta_19_ecw")
# fit_one_beta_ecw <- fit
# load("retest_combined_1929_one_beta_18_s")
# fit_one_beta_s <- fit

# summary(fit_ecw, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary
# summary(fit_s, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary
# summary(fit_one_beta_ecw, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary
# summary(fit_one_beta_s, pars = c("beta"), probs = c(0.025, 0.5, 0.975))$summary
# rescreen_draws_ecw <- rstan::extract(fit_ecw) 
# rescreen_draws_s <- rstan::extract(fit_s) 
# rescreen_draws_one_beta_ecw <- rstan::extract(fit_one_beta_ecw) 
# rescreen_draws_one_beta_s <- rstan::extract(fit_one_beta_s) 

# beta_ecw <- data.frame(Median = summary(fit_ecw, pars = c("beta"), probs = c(0.5))$summary[,c("50%")])
# beta_s <- data.frame(Median = summary(fit_s, pars = c("beta"), probs = c(0.5))$summary[,c("50%")])
# beta_one_beta_ecw <- data.frame(Median = summary(fit_one_beta_ecw, pars = c("beta"), probs = c(0.5))$summary[,c("50%")])
# beta_one_beta_s <- data.frame(Median = summary(fit_one_beta_s, pars = c("beta"), probs = c(0.5))$summary[,c("50%")])

# beta_ecw$df <- c("bndhs", "cpvsteps", "ethphia", "ghsage", "mlphia", "rwphia", "tzphia", "zmphia")
# beta_s$df <- c("ls14dhs", "ls9dhs", "sabssm", "sasage", "zwdhs", "zwphia")
# region_count_df <- unique(surveys[,c("Region", "Country", "df")])
# beta_ecw <- left_join(beta_ecw, region_count_df)
# beta_s <- left_join(beta_s, region_count_df)

# load("fit_beta")
# fit_beta
# names(fit_beta)[names(fit_beta) == "X50."] <- "Median"
# ecw_count <- c("bndhs", "cpvsteps", "ethphia", "ghsage", "mlphia", "rwphia", "tzphia", "zmphia")
# s_count <- c("ls14dhs", "ls9dhs", "sabssm", "sasage", "zwdhs", "zwphia")
# 
# beta_ecw <- fit_beta[fit_beta$df %in% ecw_count, c("Median", "df")]
# beta_s <- fit_beta[fit_beta$df %in% s_count, c("Median", "df")]
# beta_one_beta_ecw <- fit_beta[fit_beta$df == "ECW Overall", c("Median", "df")]
# beta_one_beta_s <- fit_beta[fit_beta$df == "S Overall", c("Median", "df")]
# 
# region_count_df <- unique(surveys[,c("Region", "Country", "df")])
# beta_ecw <- left_join(beta_ecw, region_count_df)
# beta_s <- left_join(beta_s, region_count_df)
# 
# #beta_ecw$Median <- beta_one_beta_ecw$Median
# #beta_s$Median <- beta_one_beta_s$Median
# 
# beta <- rbind(beta_ecw, beta_s)
# table(beta$Country)
# #only want one for each country --> will choose the latest value
# beta <- filter(beta, df != "ls9dhs" & df != "sasage" & df != "zwdhs")
# 
# val <- 50
# beta$Median <- val
# beta_one_beta_ecw$Median <- val
# beta_one_beta_s$Median <- val

load("rescreen_1")
rescreen_draws <- rstan::extract(fit)
rstan::summary(fit, pars = c("beta_ecw", "beta_s", "beta_overall", "beta_reg"), probs = c(0.025, 0.5, 0.975))$summary

size = 1000
set.seed(321)
beta_ecw <- data.frame(matrix(NA, nrow = ncol(rescreen_draws$beta_ecw), ncol = size))
beta_s <- data.frame(matrix(NA, nrow = ncol(rescreen_draws$beta_s), ncol = size))
beta_reg <- data.frame(matrix(NA, nrow = ncol(rescreen_draws$beta_reg), ncol = size))
samp = T
if(samp == T){
  for(i in 1:ncol(rescreen_draws$beta_ecw)){
    beta_ecw[i,] <- sample(rescreen_draws$beta_ecw[,i], size = size)
  }
  for(i in 1:ncol(rescreen_draws$beta_s)){
    beta_s[i,] <- sample(rescreen_draws$beta_s[,i], size = size)
  }
  for(i in 1:ncol(rescreen_draws$beta_reg)){
    beta_reg[i,] <- sample(rescreen_draws$beta_reg[,i], size = size)
  }
} else{
  for(i in 1:ncol(rescreen_draws$beta_ecw)){
    beta_ecw[i,] <- median(rescreen_draws$beta_ecw[,i])
  }
  for(i in 1:ncol(rescreen_draws$beta_s)){
    beta_s[i,] <- median(rescreen_draws$beta_s[,i])
  }
  for(i in 1:ncol(rescreen_draws$beta_reg)){
    beta_reg[i,] <- median(rescreen_draws$beta_reg[,i])
  }
  colnames(beta_ecw) <- "beta"
  colnames(beta_s) <- "beta"
  colnames(beta_reg) <- "beta"
}

#median(as.numeric(beta_reg[2,1:1000]))
#median(beta_ecw[,i])
beta_ecw$df <- c("bndhs", "cpvsteps", "ethphia", "ghsage", "mlphia", "rwphia", "tzphia", "zmphia")
beta_s$df <- c("ls14dhs", "ls9dhs", "sabssm", "sasage", "zwdhs", "zwphia")
region_count_df <- unique(surveys[,c("Region", "Country", "df")])
beta_ecw <- left_join(beta_ecw, region_count_df)
beta_s <- left_join(beta_s, region_count_df)

beta_reg$Region <- c("ECW", "Southern")

beta <- rbind(beta_ecw, beta_s)
table(beta$Country)
#only want one for each country --> will choose the latest value
beta <- filter(beta, df != "ls9dhs" & df != "sasage" & df != "zwdhs")

col <- ncol(prob_by_count_age_year_3y)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1

rates_count <- risk_to_rate(prob_by_count_age_year_3y[,index:col], 3)
rates_count <- cbind(prob_by_count_age_year_3y[,1:3], rates_count)
region_country <- read_excel("surveys.xlsx", sheet = 7)
rates_count <- left_join(rates_count, region_country)

#indicate countries with their own beta
beta_countries <- unique(beta$Country)
rates_count$has_beta[rates_count$Region == "Southern"] <- "N;S"
rates_count$has_beta[rates_count$Region != "Southern"] <- "N;ECW"
rates_count$has_beta[rates_count$Country %in% beta_countries] <- "Y"

#merge rates_count dataset with beta dataset
indv_beta = 2
if(indv_beta == 1){
  #first split the dataset into two (one that will use individual betas and one for overall betas)
  rates_count_1 <- rates_count[rates_count$has_beta == "Y", ]
  rates_count_2 <- rates_count[rates_count$has_beta != "Y", ]
  rates_count_2$Region[rates_count_2$Region != "Southern"] <- "ECW"
  
  rates_count_1 <- left_join(rates_count_1, select(beta, -df))
  rates_count_2 <- left_join(rates_count_2, beta_reg)
  
  rates_count <- rbind(rates_count_1, rates_count_2)
} else if(indv_beta == 2){
  rates_count$Region[rates_count$Region != "Southern"] <- "ECW"
  rates_count <- left_join(rates_count, beta_reg, by = "Region")
} else if(indv_beta == 3){
  beta_manual <- data.frame(Region = c("Central", "Western", "Eastern", "Southern"),
                            beta = c(100, 100, 18, 8))
  rates_count <- left_join(rates_count, beta_manual, by = "Region")
}

#finding initial screening rates
col <- ncol(prob_by_count_age_year_3y)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1
beta_index <- ncol(rates_count) - size + 1
beta_col <- ncol(rates_count)

init_rates_count <- rates_count[,index:col]/(1 - prob_by_count_age_year_ever[,index:col] + rates_count[,beta_index:beta_col] * prob_by_count_age_year_ever[,index:col])
init_rates_count <- cbind(rates_count[,1:3], init_rates_count)

# init_rates_count <- rates_count[,5:3604]/(1 - prob_by_hiv_count_age_year_ever[,5:3604] + rates_count$Median * prob_by_hiv_count_age_year_ever[,5:3604])
# init_rates_count <- cbind(rates_count[,1:4], beta = rates_count$Median, init_rates_count)

source("life_tbl_func.R")

col <- ncol(init_rates_count)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1
prob_once_by_count <- data.frame(matrix(data = 0, nrow = length(unique(init_rates_count$Country)), ncol = iter))#*samp_size)
prob_twice_by_count <- data.frame(matrix(data = 0, nrow = length(unique(init_rates_count$Country)), ncol = iter))#*samp_size)

m = 1
hiv = F
for(c in 1:length(unique(init_rates_count$Country))){
  count <- unique(init_rates_count$Country)[c]
  beta_count <- unique(rates_count$Country)[c]
  lower_age = 30
  upper_age = 45
  if(hiv == T){
    for(h in 1:2){
      dat_rates <- filter(init_rates_count, Country == count & hivstat == (h-1))
      n = 0
      rates = data.frame(rep(0, iter))
      for(a in 2:as.numeric(max(unique(init_rates_count$agegr)))){
        r <- filter(dat_rates, agegr == a)
        r1 <- t(filter(r, Year == n)[,index:col])
        r2 <- t(filter(r, Year == (n+1))[,index:col])
        r3 <- t(filter(r, Year == (n+2))[,index:col])
        r4 <- t(filter(r, Year == (n+3))[,index:col])
        r5 <- t(filter(r, Year == (n+4))[,index:col])
        n = n+5
        rates = cbind(rates, r1, r2, r3, r4, r5)
      }
      rates = rates[,-1]
      df_beta <- data.frame(t(rep((unique(dat_rates$beta)), upper_age - lower_age)))
      
      wom_once <- filter(prob_by_hiv_count_year_ever_2529, Country == count & hivstat == (h-1))[,4:1003] #will need to change if niter changes for these values
      wom_never <- 1 - wom_once
      
      val <- life_table_ras2(rates = rates, beta = df_beta,
                             lower_age = lower_age,
                             upper_age = upper_age,
                             dt = 0.1,
                             wom_once = wom_once, 
                             wom_never = wom_never)
      
      prob_once_by_count[m, 1:iter] <- as.numeric(val[1,], 1)
      prob_twice_by_count[m, 1:iter] <- as.numeric(val[2,], 1)
      
      prob_once_by_count[m,iter+1] <- count
      prob_once_by_count[m,iter+2] <- h-1
      
      prob_twice_by_count[m,iter+1] <- count
      prob_twice_by_count[m,iter+2] <- h-1
      
      print(count)
      m = m + 1
    }
  } else{
      dat_rates <- filter(init_rates_count, Country == count)
      n = 0
      rates = data.frame(rep(0, iter))
      for(a in 2:as.numeric(max(unique(init_rates_count$agegr)))){
        r <- filter(dat_rates, agegr == a)
        r1 <- t(filter(r, Year == n)[,index:col])
        r2 <- t(filter(r, Year == (n+1))[,index:col])
        r3 <- t(filter(r, Year == (n+2))[,index:col])
        r4 <- t(filter(r, Year == (n+3))[,index:col])
        r5 <- t(filter(r, Year == (n+4))[,index:col])
        n = n+5
        
        rates = cbind(rates, r1, r2, r3, r4, r5)
      }
      rates = rates[,-1]
      df_beta <- NULL
      if(samp == T){
        dat_beta <- data.frame(t(unique(filter(rates_count, Country == count)[,beta_index:beta_col])))
        df_beta <- cbind(dat_beta, rep(dat_beta[1], upper_age - lower_age))
      }else{
        dat_beta <- unique(filter(rates_count, Country == count)[,beta_index:beta_col])
        df_beta <- data.frame(t(rep(dat_beta, upper_age - lower_age)))
      }
      wom_once <- filter(prob_by_count_year_ever_2529, Country == count)[,(index - 1):(col - 1)] #will need to change if niter changes for these values
      wom_never <- 1 - wom_once
      
      val <- life_table_ras2(rates = rates, beta = df_beta,
                             lower_age = lower_age,
                             upper_age = upper_age,
                             dt = 0.1,
                             wom_once = wom_once, 
                             wom_never = wom_never)
      
      prob_once_by_count[m, 1:iter] <- as.numeric(val[1,], 1)
      prob_twice_by_count[m, 1:iter] <- as.numeric(val[2,], 1)
      
      prob_once_by_count[m,iter+1] <- count
      
      prob_twice_by_count[m,iter+1] <- count
      
      print(count)
      m = m + 1
  }
}
# names(prob_once_by_count)[(iter+1):(iter+2)] <- c("Country", "hivstat")
# names(prob_twice_by_count)[(iter+1):(iter+2)] <- c("Country", "hivstat")
# store1 <- prob_once_by_count #store in another variable just in case
# store2 <- prob_twice_by_count #store in another variable just in case
# prob_once_by_count <- store1
# prob_twice_by_count <- store2

# prob_once_by_count_b50 <- prob_once_by_count
# prob_twice_by_count_b50 <- prob_twice_by_count
# prob_once_by_count <- prob_once_by_count_b50 
# prob_twice_by_count <- prob_twice_by_count_b50 

names(prob_once_by_count)[(iter+1)] <- c("Country")
names(prob_twice_by_count)[(iter+1)] <- c("Country")
region_country <- read_excel("surveys.xlsx", sheet = 7)

prob_once_by_count <- left_join(prob_once_by_count, region_country)
prob_twice_by_count <- left_join(prob_twice_by_count, region_country)

#can also add method predictions for screening ever
#prob once by count should be using those 50 only

prob_once_by_count <- prob_once_by_count[,c((iter+1), 1:iter)]
prob_twice_by_count <- prob_twice_by_count[,c((iter+1), 1:iter)]
# prob_once_by_count <- prob_once_by_count[,c((iter+1):(iter+2), 1:iter)]
# prob_twice_by_count <- prob_twice_by_count[,c((iter+1):(iter+2), 1:iter)]
prob_twice_among_once_by_count <- prob_twice_by_count[,2:(iter+1)]/prob_once_by_count[,2:(iter+1)]

prob_once_by_count <- cbind(Country = prob_once_by_count[,1], row_sumstats(prob_once_by_count[,2:(iter+1)], prob = c(0.5, 0.025, 0.975)))
names(prob_once_by_count)[names(prob_once_by_count) == "50%"] <- "median"
names(prob_once_by_count)[names(prob_once_by_count) == "2.5%"] <- "lower"
names(prob_once_by_count)[names(prob_once_by_count) == "97.5%"] <- "upper"

prob_twice_by_count<- cbind(Country = prob_twice_by_count[,1], row_sumstats(prob_twice_by_count[,2:(iter+1)], prob = c(0.5, 0.025, 0.975)))
names(prob_twice_by_count)[names(prob_twice_by_count) == "50%"] <- "median"
names(prob_twice_by_count)[names(prob_twice_by_count) == "2.5%"] <- "lower"
names(prob_twice_by_count)[names(prob_twice_by_count) == "97.5%"] <- "upper"

prob_twice_among_once_by_count<- cbind(Country = prob_twice_by_count[,1], row_sumstats(prob_twice_among_once_by_count, prob = c(0.5, 0.025, 0.975)))
names(prob_twice_among_once_by_count)[names(prob_twice_among_once_by_count) == "50%"] <- "median"
names(prob_twice_among_once_by_count)[names(prob_twice_among_once_by_count) == "2.5%"] <- "lower"
names(prob_twice_among_once_by_count)[names(prob_twice_among_once_by_count) == "97.5%"] <- "upper"

#View(prob_twice_among_once_by_count[prob_twice_among_once_by_count$Country %in% countries_2p, c("Country", "median", "lower", "upper")])
a <- paste0(round(prob_twice_among_once_by_count[,"median"]*100, 0), "% (", round(prob_twice_among_once_by_count[,"lower"]*100, 0), "-", round(prob_twice_among_once_by_count[,"upper"]*100, 0), "%)")
b <- data.frame(Country = prob_twice_among_once_by_count[,"Country"], a)
#View(b[b$Country %in% countries_2p,])

col <- ncol(prob_by_count_year_ever_4044)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1
prob_by_ever <- cbind(Country = prob_by_count_year_ever_4044[,"Country"], row_sumstats(data.frame(prob_by_count_year_ever_4044[,index:col]), prob = c(0.5, 0.025, 0.975)))
#prob_by_ever <- cbind(Country = prob_by_hiv_count_year_ever_4549[,"Country"], hivstat = prob_by_hiv_count_year_ever_4549[,"hivstat"], row_sumstats(data.frame(prob_by_hiv_count_year_ever_4549[,4:3603]), prob = c(0.5, 0.025, 0.975)))
names(prob_by_ever)[names(prob_by_ever) == "50%"] <- "median"
names(prob_by_ever)[names(prob_by_ever) == "2.5%"] <- "lower"
names(prob_by_ever)[names(prob_by_ever) == "97.5%"] <- "upper"

#View(prob_twice_by_count)
df_weights <- weights(agegr = c("40-44"), hiv = F, years = 2020)
a <- left_join(prob_twice_by_count, df_weights)
b <- left_join(prob_once_by_count, df_weights)

#by region
prob_once_by_count_ecw <- filter(prob_once_by_count, Region != "Southern")
prob_once_by_count_s <- filter(prob_once_by_count, Region == "Southern")
prob_twice_by_count_ecw <- filter(prob_twice_by_count, Region != "Southern")
prob_twice_by_count_s <- filter(prob_twice_by_count, Region == "Southern")

a_ecw <- left_join(prob_twice_by_count_ecw, df_weights)
a_s <- left_join(prob_twice_by_count_s, df_weights)
b_ecw <- left_join(prob_once_by_count_ecw, df_weights)
b_s <- left_join(prob_once_by_count_s, df_weights)

#for general
col <- ncol(prob_by_count_age_year_3y)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1
prob_once_by_count_weighted <- sweep(b[,index:col], 1, b$poststrat_weight, "*")
prob_twice_by_count_weighted <- sweep(a[,index:col], 1, a$poststrat_weight, "*")
prop_twice_among_once_weighted <- sweep(a[,index:col]/b[,index:col], 1, b$poststrat_weight, "*")

prob_twice_by_count_weighted <- sweep(a[,index:col], 1, a$poststrat_weight, "*")
pooled_estimate_once <- quantile(colSums(prob_once_by_count_weighted)/sum(b$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice <- quantile(colSums(prob_twice_by_count_weighted)/sum(a$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice_among_once <- quantile(colSums(prop_twice_among_once_weighted)/sum(b$poststrat_weight), probs = c(0.025, 0.5, 0.975))

pooled_estimate_once_clean <- data.frame(Country = "Pooled",
                                         median = as.numeric(pooled_estimate_once[2]),
                                         lower = as.numeric(pooled_estimate_once[1]),
                                         upper = as.numeric(pooled_estimate_once[3]))

pooled_estimate_twice_clean <- data.frame(Country = "Pooled",
                                         median = as.numeric(pooled_estimate_twice[2]),
                                         lower = as.numeric(pooled_estimate_twice[1]),
                                         upper = as.numeric(pooled_estimate_twice[3]))

#from model
prob_ever_4044 <- prob_final(pred_prob_weighted_4044[[2]], dat_pred = pred_prob_weighted_4044[[3]], 1, "Year")

pooled_estimate_once_mod_clean <- data.frame(Country = "Pooled",
                                           median = as.numeric(prob_ever_4044["50%"]),
                                           lower = as.numeric(prob_ever_4044["2.5%"]),
                                           upper = as.numeric(prob_ever_4044["97.5%"]))

#by region
col <- ncol(prob_by_count_age_year_3y)
#iter <- nrow(draws$alpha)
iter <- 1000
index <- col - iter + 1
prob_once_by_count_ecw_weighted <- sweep(b_ecw[,index:col], 1, b_ecw$poststrat_weight, "*")
prob_twice_by_count_ecw_weighted <- sweep(a_ecw[,index:col], 1, a_ecw$poststrat_weight, "*")
prop_twice_among_once_ecw_weighted <- sweep(a_ecw[,index:col]/b_ecw[,index:col], 1, b_ecw$poststrat_weight, "*")

prob_twice_by_count_ecw_weighted <- sweep(a_ecw[,index:col], 1, a_ecw$poststrat_weight, "*")
pooled_estimate_once_ecw <- quantile(colSums(prob_once_by_count_ecw_weighted)/sum(b_ecw$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice_ecw <- quantile(colSums(prob_twice_by_count_ecw_weighted)/sum(a_ecw$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice_ecw_among_once <- quantile(colSums(prop_twice_among_once_ecw_weighted)/sum(b_ecw$poststrat_weight), probs = c(0.025, 0.5, 0.975))

pooled_estimate_once_ecw_clean <- data.frame(Country = "Pooled",
                                         median = as.numeric(pooled_estimate_once_ecw[2]),
                                         lower = as.numeric(pooled_estimate_once_ecw[1]),
                                         upper = as.numeric(pooled_estimate_once_ecw[3]))

pooled_estimate_twice_ecw_clean <- data.frame(Country = "Pooled",
                                          median = as.numeric(pooled_estimate_twice_ecw[2]),
                                          lower = as.numeric(pooled_estimate_twice_ecw[1]),
                                          upper = as.numeric(pooled_estimate_twice_ecw[3]))

prob_once_by_count_s_weighted <- sweep(b_s[,index:col], 1, b_s$poststrat_weight, "*")
prob_twice_by_count_s_weighted <- sweep(a_s[,index:col], 1, a_s$poststrat_weight, "*")
prop_twice_among_once_s_weighted <- sweep(a_s[,index:col]/b_s[,index:col], 1, b_s$poststrat_weight, "*")

prob_twice_by_count_s_weighted <- sweep(a_s[,index:col], 1, a_s$poststrat_weight, "*")
pooled_estimate_once_s <- quantile(colSums(prob_once_by_count_s_weighted)/sum(b_s$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice_s <- quantile(colSums(prob_twice_by_count_s_weighted)/sum(a_s$poststrat_weight), probs = c(0.025, 0.5, 0.975))
pooled_estimate_twice_s_among_once <- quantile(colSums(prop_twice_among_once_s_weighted)/sum(b_s$poststrat_weight), probs = c(0.025, 0.5, 0.975))

pooled_estimate_once_s_clean <- data.frame(Country = "Pooled",
                                             median = as.numeric(pooled_estimate_once_s[2]),
                                             lower = as.numeric(pooled_estimate_once_s[1]),
                                             upper = as.numeric(pooled_estimate_once_s[3]))

pooled_estimate_twice_s_clean <- data.frame(Country = "Pooled",
                                              median = as.numeric(pooled_estimate_twice_s[2]),
                                              lower = as.numeric(pooled_estimate_twice_s[1]),
                                              upper = as.numeric(pooled_estimate_twice_s[3]))

#from model
prob_ever_4044_by_reg <- prob_final(pred_prob_weighted_4044[[2]], dat_pred = pred_prob_weighted_4044[[3]], 2, "Region", "Year")

pooled_estimate_once_mod_clean_by_reg <- data.frame(Country = c("Western/Central/Eastern Pooled", "Southern Pooled"),
                                             median = as.numeric(prob_ever_4044_by_reg["50%"]),
                                             lower = as.numeric(prob_ever_4044_by_reg["2.5%"]),
                                             upper = as.numeric(prob_ever_4044_by_reg["97.5%"]))
#for HIV
# a_hivn <- filter(a, hivstat == 0)
# a_hivp <- filter(a, hivstat == 1)
# prob_twice_by_count_weighted_hivn <- sweep(a_hivn[,3:3601], 1, a_hivn$poststrat_weight, "*")
# prob_twice_by_count_weighted_hivp <- sweep(a_hivp[,3:3601], 1, a_hivp$poststrat_weight, "*")
# pooled_estimate_twice_hivn <- quantile(colSums(prob_twice_by_count_weighted_hivn)/sum(a_hivn$poststrat_weight), probs = c(0.025, 0.5, 0.975))
# pooled_estimate_twice_hivp <- quantile(colSums(prob_twice_by_count_weighted_hivp)/sum(a_hivp$poststrat_weight), probs = c(0.025, 0.5, 0.975))
# 
# b_hivn <- filter(b, hivstat == 0)
# b_hivp <- filter(b, hivstat == 1)
# prob_once_by_count_weighted_hivn <- sweep(b_hivn[,3:3601], 1, b_hivn$poststrat_weight, "*")
# prob_once_by_count_weighted_hivp <- sweep(b_hivp[,3:3601], 1, b_hivp$poststrat_weight, "*")
# pooled_estimate_once_hivn <- quantile(colSums(prob_once_by_count_weighted_hivn)/sum(b_hivn$poststrat_weight), probs = c(0.025, 0.5, 0.975))
# pooled_estimate_once_hivp <- quantile(colSums(prob_once_by_count_weighted_hivp)/sum(b_hivp$poststrat_weight), probs = c(0.025, 0.5, 0.975))
# 
# #add for pooled estimates
# pooled_estimate_once_hiv <- data.frame(Country = c("Pooled", "Pooled"),
#                                        hivstat = c(0, 1),
#                                        median = c(pooled_estimate_once_hivn[2], pooled_estimate_once_hivp[2]),
#                                        lower = c(pooled_estimate_once_hivn[1], pooled_estimate_once_hivp[1]),
#                                        upper = c(pooled_estimate_once_hivn[3], pooled_estimate_once_hivp[3]))
# 
# pooled_estimate_twice_hiv <- data.frame(Country = c("Pooled", "Pooled"),
#                                        hivstat = c(0, 1),
#                                        median = c(pooled_estimate_twice_hivn[2], pooled_estimate_twice_hivp[2]),
#                                        lower = c(pooled_estimate_twice_hivn[1], pooled_estimate_twice_hivp[1]),
#                                        upper = c(pooled_estimate_twice_hivn[3], pooled_estimate_twice_hivp[3]))
# #from model
# prob_hiv_ever_4549 <- prob_final(pred_prob_weighted_4549[[2]], dat_pred = pred_prob_weighted_4549[[3]], 2, "Year", "hivstat")
# prob_ever_4549 <- prob_final(pred_prob_weighted_4549[[2]], dat_pred = pred_prob_weighted_4549[[3]], 1, "Year")
# #View(prob_hiv_ever_4549[, c("2.5%", "50%", "97.5%")])
# 
# pooled_estimate_once_mod_hiv <- data.frame(Country = c("Pooled", "Pooled"),
#                                            hivstat = c(0, 1),
#                                            median = c(prob_hiv_ever_4549[1, "50%"], prob_hiv_ever_4549[2, "50%"]),
#                                            lower = c(prob_hiv_ever_4549[1, "2.5%"], prob_hiv_ever_4549[2, "2.5%"]),
#                                            upper = c(prob_hiv_ever_4549[1, "97.5%"], prob_hiv_ever_4549[2, "97.5%"]))

# ----  K) Plots - Life table values ----
#create data frame with both prop screening twice and prop screening once
prop_twice <- rbind(prob_twice_by_count[,c("Country", "median", "lower", "upper")], pooled_estimate_twice_clean)
prop_once <- rbind(prob_once_by_count[,c("Country", "median", "lower", "upper")], pooled_estimate_once_clean)
prop_once_mod <- rbind(prob_by_ever[,c("Country", "median", "lower", "upper")], pooled_estimate_once_mod_clean)

# prop_twice <- rbind(prob_twice_by_count[,c("Country", "hivstat", "median", "lower", "upper")], pooled_estimate_twice_hiv)
# prop_once <- rbind(prob_once_by_count[,c("Country", "hivstat", "median", "lower", "upper")], pooled_estimate_once_hiv)
# prop_once_mod <- rbind(prob_by_ever[,c("Country", "hivstat", "median", "lower", "upper")], pooled_estimate_once_mod_hiv)
  
prop_twice$num_screening <- "From life table method: Screening At Least Twice Ever"
prop_once$num_screening <- "From life table method: Screening At Least Once Ever"
prop_once_mod$num_screening <- "From model: Screening Ever"
prop <- rbind(prop_twice, prop_once, prop_once_mod)
prop$Country <- factor(prop$Country, levels = unique(prop$Country), labels = unique(prop$Country))
View(prop[prop$Country == "Pooled", ])

region_country <- surveys[,c("Country", "Region")]
prop_twice <- left_join(prop_twice, region_country)
prop_twice$Country[prop_twice$Country == "Tanzania"] <- "United Republic of Tanzania"
prop_twice$Country[prop_twice$Country == "Congo (Rep)"] <- "Republic of the Congo"
prop_twice$Region[is.na(prop_twice$Region)] = "All"
prop_twice$Region[prop_twice$Region == "Central" | prop_twice$Region == "Western"] = "Central/Western"
#countries_2p <- c("Benin", "Burkina Faso", "Cote d'Ivoire", "Eswatini", "Ethiopia", "Ghana", "Kenya", "Lesotho", "Malawi", "Namibia", "Senegal", "South Africa", "Zambia", "Zimbabwe")
countries_2p <- c("Benin", "Cape Verde", "Ethiopia", "Ghana", "Malawi", "Rwanda", "Tanzania", "Zambia", "Lesotho", "South Africa", "Zimbabwe")
p_twice <- ggplot(prop_twice[prop_twice$Country %in% countries_2p | prop_twice$Country == "Pooled", ], aes(x = Country, y = median)) +
  geom_bar(stat = "identity", position = "dodge", fill = "salmon2") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9))  +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  ylab("Proportion") +
  xlab("Country") +
  # labs(title = "Proportion of women aged 49 who have been screened twice lifetime in 2020",
  #      caption = "Note: denominator is among those who were never screened at the age of 30") +
  theme_bw()  +
  #scale_fill_nejm() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90, vjust=0.5, hjust = 1),
        plot.caption = element_text(hjust = 0),
        legend.position = "bottom") +
  facet_grid(.~Region, scales = "free", space = "free") +
  ylim(0,1);p_twice

region_country <- surveys[,c("Country", "Region")]
prop <- left_join(prop, region_country)
prop$Country[prop$Country == "Tanzania"] <- "United Republic of Tanzania"
prop$Country[prop$Country == "Congo (Rep)"] <- "Republic of the Congo"
prop$Region[is.na(prop$Region)] = "All"
prop$Region[prop$Region == "Central" | prop$Region == "Western"] = "Central/Western"

p_twice_check <- ggplot(prop[prop$Country %in% countries_2p | prop$Country == "Pooled", ], aes(x = Country, y = median, fill = num_screening)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9))  +
  # labs(fill = "Number of Screenings",
  #      title = "Robustness check for proportion of women who are 49 who have been screened twice lifetime in 2020\nRecall period = 1",
  #      caption = "Note there is a difference in denominators: \n'From life table methods' is among those who were were never screened at the age of 30; \n'From model' is among everyone") +
  ylab("Proportion") +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  theme_bw()  +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90, vjust=0.5, hjust = 1),
        plot.caption = element_text(hjust = 0),
        legend.position = "bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  scale_fill_brewer(palette = "Set2", name = "Recall Period", labels = c("From life table method: Lifetime screening at least once", "From life table method: Lifetime screening at least twice", "From time trends model: Lifetime screening at least once")) +
  facet_grid(.~Region, scales = "free", space = "free") +
  ylim(0,1);p_twice_check

prop$hivstat <- factor(prop$hivstat, levels = c("0", "1"), labels = c("HIV-", "HIV+"))
prop$Country <- factor(prop$Country, levels = unique(prop$Country), labels = unique(prop$Country))
p_hiv_twice_check <- ggplot(prop, aes(x = Country, y = median, fill = num_screening)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9))  +
  labs(fill = "Number of Screenings",
       title = "Robustness check for proportion of women who are 49 who have been screened twice lifetime in 2020\nRecall period = 1",
       caption = "Note there is a difference in denominators: \n'From life table methods' is among those who were were never screened at the age of 30; \n'From model' is among everyone") +
  ylab("Proportion") +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  theme_minimal()  +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text.x = element_text(size = 12, angle = 60, hjust=0.95, vjust=1.1,),
        plot.caption = element_text(hjust = 0),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set2") +
  ylim(0,1) +
  facet_wrap(~hivstat);p_hiv_twice_check

save(rates_count, init_rates_count,
     prob_twice_by_count, prob_once_by_count, prob_by_ever,
     prob_by_count_year_ever_4044,
     p_twice, p_twice_check,
     pooled_estimate_once, pooled_estimate_twice,
     # p_hiv_twice_check,
     # prob_by_hiv_count_year_ever_4549,
     # pooled_estimate_once_hivn, pooled_estimate_once_hivp,
     # pooled_estimate_twice_hivn, pooled_estimate_twice_hivp,
     file = "prob_twice_final_only_reg")
#load("prob_twice_new_14")
# prop_twice$Country_fr <- factor(prop_twice$Country, labels = count_fr)
# prop$num_screening_fr <- factor(prop$num_screening, labels = c("À partir de la méthode des tables de mortalité : dépistage au moins une fois",
#                                                                "À partir de la méthode des tables de mortalité : dépistage au moins deux fois",
#                                                                "Du modèle : dépistage jamais"))

# p_twice_check_fr <- ggplot(prop, aes(x = Country, y = median, fill = num_screening_fr)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9))  + 
#   labs(fill = "Nombre de dépistages",
#        title = "Contrôle de robustesse pour la proportion de femmes de 49 ans qui ont été dépistées deux fois au cours de leur vie en 2020",
#        caption = "Notez qu'il y a une différence dans les dénominateurs :\nLes « méthodes des tables de mortalité » font partie de ceux qui n'ont jamais été dépistés à l'âge de 30 ans ; \n'Du modèle' est parmi tout le monde") +
#   ylab("Proportion") +
#   xlab("Pays") +
#   geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
#   theme_minimal()  +
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=17),
#         axis.text.x = element_text(size = 12, angle = 60, hjust=0.95, vjust=1.25),
#         plot.caption = element_text(hjust = 0),
#         legend.position = "bottom") +
#   scale_fill_brewer(palette = "Set2") +
#   ylim(0,1);p_twice_check_fr

# set.seed(123)
# samp_size = 15
# beta_ecw <- data.frame(matrix(0, nrow = ncol(rescreen_draws_ecw$beta), ncol = samp_size))
# for(i in 1:ncol(rescreen_draws_ecw$beta)){
#   beta_ecw[i,] <- sample(rescreen_draws_ecw$beta[,i], size = samp_size, replace = F)
# }
# rowMeans(beta_ecw)
# 
# beta_s <- data.frame(matrix(0, nrow = ncol(rescreen_draws_s$beta), ncol = samp_size))
# for(i in 1:ncol(rescreen_draws_s$beta)){
#   beta_s[i,] <- sample(rescreen_draws_s$beta[,i], size = samp_size, replace = F)
# }
# rowMeans(beta_s)
# 
# beta_one_beta_ecw <- sample(rescreen_draws_one_beta_ecw$beta, size = samp_size, replace = F)
# beta_one_beta_s <- sample(rescreen_draws_one_beta_s$beta, size = samp_size, replace = F)
# 
# beta_ecw$df <- c("bndhs", "cpvsteps", "ethphia", "ghsage", "mlphia", "rwphia", "tzphia", "zmphia", "zwdhs", "zwphia")
# beta_s$df <- c("ls14dhs", "ls9dhs", "sabssm", "sasage")
# region_count_df <- unique(surveys[,c("Region", "Country", "df")])
# beta_ecw <- left_join(beta_ecw, region_count_df)
# beta_s <- left_join(beta_s, region_count_df)

# beta <- rbind(beta_ecw, beta_s)
# table(beta$Country)
# #only want one for each country --> will choose the latest value
# beta <- filter(beta, df != "ls9dhs" & df != "sasage" & df != "zwdhs")
# table(beta$Country)
# 
# #naming the beta iterations
# colnames <- paste0("beta", c(1:samp_size))
# names(beta)[1:samp_size] <- colnames
# names(beta_one_beta_ecw)[1:samp_size] <- colnames
# names(beta_one_beta_s)[1:samp_size] <- colnames
# 
# rates_count <- risk_to_rate(prob_by_count_age_year_3y[,4:3603], 3)
# rates_count <- cbind(prob_by_count_age_year_3y[,1:3], rates_count)
# region_count <- unique(region_count_df[,c("Country", "Region")])
# rates_count <- left_join(rates_count, region_count)
# 
# rates_count_list <- vector(mode = "list", length = samp_size)
# init_rates_count_list <- vector(mode = "list", length = samp_size)
# for(i in 1:samp_size){
#   beta_i <- paste0("beta", i)
#   rates_count_list[[i]] <- left_join(rates_count, beta[,c(beta_i, "Country", "Region", "df") ])
#   names(rates_count_list[[i]])[names(rates_count_list[[i]]) == beta_i] <- "beta"
#   
#   for(j in 1:nrow(rates_count_list[[i]])){
#     if(is.na(rates_count_list[[i]]$df[j])){
#       if(rates_count_list[[i]]$Region[j] == "Southern" ){
#         rates_count_list[[i]]$beta[j] <- beta_one_beta_s[beta_i]
#       }else{
#         rates_count_list[[i]]$beta[j] <- beta_one_beta_ecw[beta_i]
#       }
#     }
#   }
#   init_rates_count_list[[i]] <- rates_count_list[[i]][,4:3603]/(1 - prob_by_count_age_year_ever[,4:3603] + rates_count_list[[i]]$beta * prob_by_count_age_year_ever[,4:3603])
#   init_rates_count_list[[i]] <- cbind(rates_count_list[[i]][,1:3], beta = rates_count_list[[i]]$beta, init_rates_count_list[[i]])
#   print(i)
# }
# save(rates_count_list, init_rates_count_list,
#      file = "rates_list_2")
# load("rates_list_2")
# rates <- array(as.numeric(unlist(rates_count_list)), dim = c(nrow(rates_count_list[[1]]), ncol(rates_count_list[[1]]), samp_size))
# b <- df_beta
# beta <- list(b,b)
# View(rates_count_list[[i]][,c("Country", "beta", "df")])

#  dat_rates <- vector(mode = "list", length = samp_size)
#  rates_list <- vector(mode = "list", length = samp_size)
#  df_beta_list <- vector(mode = "list", length = samp_size)
# # a <- array(data = 0, dim = c(2240, 3604, samp_size))
#  for(i in 1:samp_size){
#    dat_rates[[i]] <- filter(init_rates_count_list[[i]], Country == unique(init_rates_count_list[[i]]$Country)[c])
#    n = 0
#    rates <- data.frame(rep(0, iter))
#    for(a in 2:as.numeric(max(unique(init_rates_count_list[[i]]$agegr)))){
#      r <- filter(dat_rates[[i]], agegr == a)
#      r1 <- t(filter(r, Year == n)[,index:col])
#      r2 <- t(filter(r, Year == (n+1))[,index:col])
#      r3 <- t(filter(r, Year == (n+2))[,index:col])
#      r4 <- t(filter(r, Year == (n+3))[,index:col])
#      r5 <- t(filter(r, Year == (n+4))[,index:col])
#      n = n+5
#      rates = cbind(rates, r1, r2, r3, r4, r5)
#    }
#    rates = rates[,-1]
#    rates_list[[i]] <- rates
#    df_beta_list[[i]] <- data.frame(t(rep(unique(dat_rates[[i]]$beta), upper_age - lower_age)))
#    print(i)
#  }

# r <- array(as.numeric(unlist(rates_list)), dim = c(nrow(rates_list[[1]]), ncol(rates_list[[1]]), samp_size))
# b <- array(as.numeric(unlist(df_beta_list)), dim = c(nrow(df_beta_list[[1]]), ncol(df_beta_list[[1]]), samp_size))
#   
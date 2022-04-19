#Bayesian mod to investigate linkage to care
library(tidyr)
library(dplyr)
library(rstan)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(GGally)
library(readxl)
library(matrixStats)
surveys <- read_excel("surveys.xlsx", sheet = 1)

#' -----------------------------------------------------------------
# ---- Prepare/Clean data ----
#' -----------------------------------------------------------------

##I first load the pooled surveys that I've previously created
load("Cleaned Pooled Surveys")

##Prepare datasets to investigate those who received their last cc test results
dat_resultreceived <- filter(pooled_surveys, !is.na(last_cc_test_result_received))

##Prepare datasets to investigate those who received a follow up appointment among those with abnormal/positive test result
dat_fu <- filter(pooled_surveys, !is.na(last_cc_test_fu))

##Prepare datasets to investigate those who've received treatment among those with abnormal/positive test result
dat_treated <- filter(pooled_surveys, !is.na(cc_treated))
sum(table(dat_treated$psu))

##Prepare datasets to investigate those who've received treatment on the same day among those with abnormal/positive test result
dat_treatedsameday <- filter(pooled_surveys, !is.na(treated_sameday))
dat_treatedsameday$treated_sameday[dat_treatedsameday$cc_treated == 0] <- 0

#created grouped data set
dat_treated <- filter(dat_treated, weight != 0)
dat_group <- expand.grid(df = unique(dat_treated$df)); nrow(dat_group)

#obtain the weighted proportions of ever testing for each survey 
options(survey.lonely.psu = "adjust")
for(i in 1:length(dat_group$df)){
  a <- filter(dat_treated, df == dat_group$df[i])
  if(nrow(a) <= 1){
    dat_treated$prop[dat_treated$df == dat_group$df[i]] <- NA
  } else{
    des <- svydesign(data = a, id = ~psu, weights = ~weight)
    prop <- as.numeric(svyciprop(~cc_treated, design = des)[1])
    dat_treated$prop[dat_treated$df == dat_group$df[i]] <- prop
  }
  print(i)
}
design <- 2.0
dat_treated_gr <- NULL
pb <- txtProgressBar(1, nrow(dat_group), style = 3)
for (i in 1:nrow(dat_group)) {
  setTxtProgressBar(pb, i)
  dat_group_i <- data.frame(df = dat_group[i, ])
  dat_i <- dat_treated[dat_treated$df == dat_group_i$df, ]
  dat_group_i$Country <- dat_i$Country[1]
  dat_group_i$den <- round(nrow(dat_i) / design)
  dat_group_i$num <- round(dat_i$prop[1]*dat_group_i$den) 
  dat_group_i$denominator_original <- nrow(dat_i)
  dat_group_i$numerator_original <- sum(dat_i$cc_treated)
  dat_treated_gr <- rbind(dat_treated_gr, dat_group_i)    
}

cc_ltc_mod_re_gr <- '
data{
  int<lower = 1> N; //number of observations
  int<lower = 1> C; //number of countries
  int<lower = 1> country[N];
  
  //grouped regression specific inputs
  int<lower = 0> num[N]; // numerator
  int<lower = 1> den[N]; // denominator
  
  //predictions
  int<lower = 1> P; //number of predictions
  int<lower = 1> pred_country[P];
}

parameters{
  real alpha; //intercept
  vector[C] re_country; 
  real<lower = 0, upper = 10> sd_re_country;
}

transformed parameters{
  vector[N] logit_prob;
  for(n in 1:N) {
   logit_prob[n] = alpha + re_country[country[n]]; 
  }
}

model{
  // priors
  // Overall Intercept
  alpha ~ normal(0, 10);
  
  // Random Intercepts by country
  re_country ~ normal(0, sd_re_country); // we center the random effect around the intercept (mean = alpha, and not 0)
  sd_re_country ~ cauchy(0, 5) T[0, 10]; // half-Chauchy prior
  
  //likelihood
  num ~ binomial_logit(den, logit_prob); 
}

generated quantities{
  vector[P] pred_logit_prob;
  for(p in 1:P){
    pred_logit_prob[p] = alpha + re_country[pred_country[p]];
  }
}
'

# Compile the model
cc_ltc_stan_re_gr <- rstan::stan_model(model_code = cc_ltc_mod_re_gr)

#Fit the model
df_stan <- dat_treated_gr

#prep variables for stan model
N <- nrow(df_stan)
C <- length(unique(df_stan$Country))

df_stan$count_num <- as.numeric(factor(df_stan$Country))
country <- df_stan$count_num

num <- df_stan$num
den <- df_stan$den

#prediction data frame
df_pred <- expand.grid(Country = 1:C)
P = nrow(df_pred)
pred_country = df_pred$Country

# the dataset that contains observations
data_stan <- list(N = N,
                  C = C,
                  num = num,
                  den = den,
                  country = country,
                  P = P,
                  pred_country = pred_country)

#' -----------------------------------------------------------------
# ---- Run model ----
#' -----------------------------------------------------------------
options(mc.cores = parallel::detectCores())
fit <- rstan::sampling(cc_ltc_stan_re_gr, data = data_stan, iter = 5000, chains = 6, refresh = 50,
                       warmup = 2500, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))
save(data_stan, fit, df_stan, cc_ltc_stan_re_gr, file = "ltc_m5")
rstan::rstan_options(auto_write = TRUE)

#load("ltc_m3")
#model parameters outputs/diagnostics
rstan::stan_trace(fit, pars = c("alpha", "re_country", "sd_re_country"))
rstan::summary(fit, pars = c("alpha", "re_country", "sd_re_country"), probs = c(0.025, 0.5, 0.975))$summary

#extract predicted probabilities
inv_logit <- function(x){
  a <- exp(x)/(1+exp(x))
  return(a)
}

draws <- rstan::extract(fit)
pred_prob <- t(inv_logit(draws$pred_logit_prob))
pred_prob_sum <- data.frame(rowQuantiles(pred_prob, probs = c(0.025, 0.5, 0.975)))
pred_prob_sum$Country <- sort(df_stan$Country)

pred_prob_total <- data.frame(t(quantile(pred_prob, probs = c(0.025, 0.5, 0.975))))
pred_prob_total$Country <- "Pooled"

pred_prob_df <- rbind(pred_prob_sum, pred_prob_total)
colnames(pred_prob_df) <- c("Lower", "Prop", "Upper", "Country")

#values obtained from estimates
load("CC Treatment")
raw_prob <- cctreated_df[cctreated_df$Age.Groups == "All" &
                           !is.na(cctreated_df$Country), ]
names(raw_prob)[names(raw_prob) == "Mean"] <- "Prop"
# names(raw_prob)[names(raw_prob) == "Lower"] <- "raw_lower"
# names(raw_prob)[names(raw_prob) == "Upper"] <- "raw_upper"
names(raw_prob)[names(raw_prob) == "Age.Groups"] <- "age.group"
raw_prob$Country <- factor(raw_prob$Country, levels = c("Cape Verde", "Malawi", "Tanzania", "Zambia"),
                           labels = c("Cape Verde \n(2019)","Malawi \n(2017-18)", "Tanzania \n(2018-19)", "Zambia \n(2016-17)"))

## Plots
pred_prob_df$Country <- factor(pred_prob_df$Country, levels = c("Cape Verde", "Malawi", "Tanzania", "Zambia", "Pooled"),
                                  labels = c("Cape Verde \n(2019)", "Malawi \n(2017-18)", "Tanzania \n(2018-19)", "Zambia \n(2016-17)", "Pooled"))

raw_prob$Source <- "Survey Data"
pred_prob_df$Source <- "Modelled Estimates"

df_combine <- NULL
paper = T
if(paper == T){
   df_combine <- rbind(pred_prob_df[pred_prob_df$Country == "Pooled",], raw_prob[raw_prob$Country != "Pooled", c("Country", "Prop", "Lower", "Upper", "Source")])
} else{
  df_combine <- rbind(raw_prob[,c("Country", "Prop", "Lower", "Upper", "Source")], pred_prob_df[,c("Country", "Prop", "Lower", "Upper", "Source")])
}

p <- ggplot(df_combine, aes(x = Country, y = Prop, fill = "red")) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, position = position_dodge(0.9)) +
  #geom_point(data = raw_prob, aes(y = raw_mean), colour = "darkgrey", size = 3) +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  #geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper), colour = "darkgrey", size = 0.5) +
  ylim(0, 1) +
  ylab("Proportion") +
  theme_bw() +
  scale_fill_nejm() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "none");p
  #scale_fill_brewer(palette = "Set3", labels = c("Cape Verde", "Malawi", "Tanzania", "Zambia", "Pooled")) +
  # ggtitle("Proportion treated for CC among those who received a positive CC screening") +
  # labs(subtitle = "Pooled value is weighted by the female population of the countries");p

save(p, pred_prob_df, raw_prob, file = "pred_prob_ltc_m5")
#load("pred_prob_ltc_m3")

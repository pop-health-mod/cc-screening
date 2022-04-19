#Grouped logistic regression
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
det_avail <- read_excel("surveys.xlsx", sheet = 4)
eco_det <- read_excel("surveys.xlsx", sheet = 5)
screening_program <- read_excel("surveys.xlsx", sheet = 6)
region_country <- read_excel("surveys.xlsx", sheet = 7)
load("Cleaned Pooled Surveys")

source("predict_poststrat_func.R")  #will need to update these functions based on the covariates I add

# ---- A) Prep ----
##I then create datasets that remove NAs for each outcome variable (i.e. either pap smear ever, or cc test ever)
dat_pap <- filter(pooled_surveys, !is.na(ever_had_pap_smear) & !is.na(pap_smear_past_3y))
dat_cctest <- filter(pooled_surveys, !is.na(ever_cc_tested) & !is.na(cc_test_past_3y))

##I also create new dataset which has a variable for either pap smear or a general cc test (the surveys asked questions regarding pap smears only to many of the southern african countries - this approach will assume that pap smears were the only/primary form of cc testing at time of survey)
dat_anytest <- pooled_surveys

#create variable for testing ever and testing in the past 3 years
dat_anytest$any_test_ever <- NA
dat_anytest$any_test_ever[dat_anytest$ever_cc_tested == 0 | dat_anytest$ever_had_pap_smear == 0] <- 0 #if pap smear or general cc test is 0 then will code as 0
dat_anytest$any_test_ever[dat_anytest$ever_cc_tested == 1 | dat_anytest$ever_had_pap_smear == 1] <- 1 #this line allows for variable to be coded as 1 if at least one of pap smear or gen cc test was 1. If previously coded as 0, will replace

dat_anytest$any_test_past_3y <- NA
dat_anytest$any_test_past_3y[dat_anytest$cc_test_past_3y == 0 | dat_anytest$pap_smear_past_3y == 0] <- 0 #if pap smear or general cc test is 0 then will code as 0
dat_anytest$any_test_past_3y[dat_anytest$cc_test_past_3y == 1 | dat_anytest$pap_smear_past_3y == 1] <- 1 #this line allows for variable to be coded as 1 if at least one of pap smear or gen cc test was 1. If previously coded as 0, will replace

##df to have all surveys double counted - creating dataframes with ONLY ever or 3y
dat_anytest_ever <- filter(subset(dat_anytest, select = -c(any_test_past_3y)), !is.na(any_test_ever))
dat_anytest_3y <- filter(subset(dat_anytest, select = -c(any_test_ever)), !is.na(any_test_past_3y))

#add indicator variables for type of test - ever or 3y
dat_anytest_ever$time_length <- "ever"
dat_anytest_3y$time_length <- "3y"

dat_anytest_ever$y <- dat_anytest_ever$any_test_ever
dat_anytest_3y$y <- dat_anytest_3y$any_test_past_3y

#Create combined dataframe
dat_anytest_all <- merge(dat_anytest_3y, dat_anytest_ever, all = T)

dat_anytest_all$q_type <- NA
dat_anytest_all$q_type[dat_anytest_all$ever_had_pap_smear == 0 | dat_anytest_all$ever_had_pap_smear == 1 | dat_anytest_all$pap_smear_past_3y == 0 | dat_anytest_all$pap_smear_past_3y == 1] <- 0
dat_anytest_all$q_type[dat_anytest_all$ever_cc_tested == 0 | dat_anytest_all$ever_cc_tested == 1 | dat_anytest_all$cc_test_past_3y == 0 | dat_anytest_all$cc_test_past_3y == 1] <- 1
#View(unique(dat_anytest_all[,c("df", "Country", "Region", "Survey")]))

#create age groups
age.groups = data.frame(l_age = seq(15, 65, 5), u_age = c(seq(19, 64, 5), 1000))
A <- nrow(age.groups)

for(n in 1:nrow(dat_anytest_all)){
  for(a in 1:A){
    if(!is.na(dat_anytest_all$age[n]) & age.groups$l_age[a] <= dat_anytest_all$age[n] & age.groups$u_age[a] >= dat_anytest_all$age[n]){
      agegr <- paste0(age.groups$l_age[a], "-", age.groups$u_age[a])
      dat_anytest_all$agegr[n] <- agegr
    }
  }
}

#calculating design effect
run_def <- TRUE
if (run_def) {
  # alternative way
  def <- NULL
  dat_anytest_all$psu_num <- as.numeric(as.factor(dat_anytest_all$psu))
  for (i in unique(dat_anytest_all$df)) {
    dat_i <- dat_anytest_all[dat_anytest_all$df == i, ]
    avg_psu_i <- mean(aggregate(dat_i$psu, by = list(dat_i$psu), 
                                FUN = function(x) sum(!is.na(x)))$x)
    rho_i <- samplesize4surveys::ICC(y = dat_i$y, cl = dat_i$psu_num)$ICC
    def_i <- 1 + (avg_psu_i - 1) * rho_i
    def <- rbind(def, data.frame(surveyid = i, def = def_i))
    print(i)
  }
  def
  summary(def$def)
}

design <- 2.0

##we create a dataset where women with the same covariate patterns are grouped together.- note that if you want to add more variables, you need to include them here (e.g. hiv)
# remove values with NA and 0 weights
###check nm13dhs isn't being double counted
dat_anytest_all <- filter(dat_anytest_all, !is.na(agegr) & weight != 0)
dat_group <- expand.grid(df = unique(dat_anytest_all$df),
                         time_length = unique(dat_anytest_all$time_length),
                         hivstat = unique(dat_anytest_all$hivstat),
                         agegr = unique(dat_anytest_all$agegr)); nrow(dat_group)
dat_group$id <- with(dat_group, paste(df, 
                                      time_length, 
                                      hivstat, 
                                      agegr))
dat_anytest_all$id <- with(dat_anytest_all, paste(df, 
                                                  time_length, 
                                                  hivstat, 
                                                  agegr))

# to speed computations, we only select covariate patterns that are found in the data
id_to_select <- which(dat_group$id %in% dat_anytest_all$id)
dat_group <- dat_group[id_to_select , ]; nrow(dat_group)

#obtain the weighted proportions of ever testing for each survey 
for(i in 1:length(dat_group$id)){
  a <- filter(dat_anytest_all, id == dat_group$id[i])
  if(nrow(a) <= 1){
    dat_anytest_all$prop[dat_anytest_all$id == dat_group$id[i]] <- NA
  } else{
    des <- svydesign(data = a, id = ~psu, weights = ~weight)
    prop <- as.numeric(svyciprop(~y, design = des)[1])
    dat_anytest_all$prop[dat_anytest_all$id == dat_group$id[i]] <- prop
  }
  print(i)
}
# unique(dat_anytest_all[,c("id", "prop")])
# table(is.na(dat_anytest_all$prop))
View(unique(dat_anytest_all[!is.na(dat_anytest_all$hivstat), c("df", "hivstat", "prop", "time_length", "agegr", "id")]))

# this takes a long time to perform, set to TRUE if you want to re-run (takes 5-10 min)
to_update <- TRUE
if (to_update) { 
  dat_timetrends_gr <- NULL
  pb <- txtProgressBar(1, nrow(dat_group), style = 3)
  for (i in 1:nrow(dat_group)) {
    setTxtProgressBar(pb, i)
    dat_group_i <- dat_group[i, ]
    dat_i <- dat_anytest_all[dat_anytest_all$id == dat_group_i$id, ]
    dat_group_i$Country <- dat_i$Country[1]
    dat_group_i$agegr <- dat_i$agegr[1]
    dat_group_i$den <- round(nrow(dat_i) / design)
    dat_group_i$num <- round(dat_i$prop[1]*dat_group_i$den) 
    dat_group_i$denominator_original <- nrow(dat_i)
    dat_group_i$numerator_original <- sum(dat_i$y)
    dat_timetrends_gr <- rbind(dat_timetrends_gr, dat_group_i)    
  }
  
  #add values from burkina faso STEPS tabulations
  bfsteps <- c(0.079, 0.056, 0.102) #estimate for ages 30-49 will assume its the same for all age groups between those ages
  
  p_bf <- bfsteps[1]
  lci_bf <- bfsteps[2]
  uci_bf <- bfsteps[2]
  
  bfsteps_den <- round(p_bf*(1-p_bf)/((p_bf - lci_bf)/1.96)^2)
  bfsteps_num <- round(bfsteps_den*p_bf)
  
  df_bfsteps <- data.frame(agegr = c("30-34", "35-39", "40-44", "45-49", "30-49"),
                           prob = rep(bfsteps[1], 5),
                           denominator_original = c(rep(bfsteps_den/4, 4), bfsteps_den),
                           numerator_original = c(rep(bfsteps_num/4, 4), bfsteps_num))
  #here i will apply the design effect to the denominators - should double check if this is the appropriate method
  df_bfsteps$den <- round(df_bfsteps$denominator_original/design)
  df_bfsteps$num <- round(df_bfsteps$den * df_bfsteps$prob)
  
  #add values from Mozambique STEPS tabulations
  mzsteps_3055 <- 0.035
  mzsteps_den <- rep(871/5, 5)
  
  df_mzsteps <- data.frame(agegr = c("30-34", "35-39", "40-44", "45-49", "30-49"),
                           prob = rep(mzsteps_3055, 5),
                           denominator_original = c(rep(871/5, 4), 871),
                           numerator_original = c(rep(871/5*mzsteps_3055, 4), 871*mzsteps_3055))
  df_mzsteps$den <- round(df_mzsteps$denominator_original/design)
  df_mzsteps$num <- round(df_mzsteps$den * df_mzsteps$prob)
  
  # sasage <- c(0.358, 0.292, 0.32, 0.104)
  # sasage_den <- c(963, 623, 311, 110)/2
  
  #add values from zw20phia
  zw20phia <- c(0.017, 0.193, 0.27, 0.364, 0.325, 0.382, 0.429, 0.473, 0.389, 0.354, 0.186)
  zw20phia_den <- c(66, 117, 168, 257, 318, 307, 389, 176, 162, 115, 76)
  zw20phia_num <- zw20phia*zw20phia_den
  
  df_zw20phia <- data.frame(agegr = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-1000"),
                           prob = zw20phia,
                           denominator_original = zw20phia_den,
                           numerator_original = zw20phia_num)
  df_zw20phia$den <- round(df_zw20phia$denominator_original/design)
  df_zw20phia$num <- round(df_zw20phia$den * df_zw20phia$prob)
  
  #add extra columns to data frames to match dat_timetrends_gr
  df_bfsteps$df <- "bfsteps"
  df_bfsteps$time_length <- "ever"
  df_bfsteps$Country <- "Burkina Faso"
  df_bfsteps$hivstat <- NA
  
  df_mzsteps$df <- "mzsteps"
  df_mzsteps$time_length <- "ever"
  df_mzsteps$Country <- "Mozambique"
  df_mzsteps$hivstat <- NA
  
  df_zw20phia$df <- "zw20phia"
  df_zw20phia$time_length <- "ever"
  df_zw20phia$Country <- "Zimbabwe"
  df_zw20phia$hivstat <- 1

  df_to_add <- rbind(df_bfsteps, df_mzsteps, df_zw20phia)
  
  df_to_add$id <- with(df_to_add, paste(df, time_length, agegr))
  df_to_add <- select(df_to_add, -prob)
  
  dat_timetrends_gr <- rbind(dat_timetrends_gr, df_to_add)
  
  close(pb)      
  nrow(dat_timetrends_gr)
  
  #add year, regions and survey type to the df
  x <- unique(surveys[,c("Country", "Region", "Survey", "df", "Year", "ISO3")])
  dat_timetrends_gr <- left_join(dat_timetrends_gr, x)
  
  #add hiv prevalence data 
  #right now I'll just add the 15+
  f_hivprev_15 <- read.csv("f_hivprev_15+.csv")
  toMatch <- c("Footnote", "lower", "upper")
  f_hivprev_15 <- f_hivprev_15[,-grep(paste(toMatch, collapse = "|"), colnames(f_hivprev_15))]
  colnames(f_hivprev_15) <- str_replace(colnames(f_hivprev_15), "X", "")
  f_hivprev_15 <- pivot_longer(f_hivprev_15, cols = -1, names_to = "Year", values_to = "hivprev15" )
  f_hivprev_15$Country[f_hivprev_15$Country == "United Republic of Tanzania"] <- "Tanzania"
  f_hivprev_15$Country[f_hivprev_15$Country == "Congo"] <- "Congo (Rep)"
  f_hivprev_15$Country[f_hivprev_15$Country == "Cabo Verde"] <- "Cape Verde"
  f_hivprev_15$hivprev15[f_hivprev_15$hivprev15 == "<0.1 "] <- 0.1
  f_hivprev_15$Year <- as.numeric(f_hivprev_15$Year)
  f_hivprev_15$hivprev15 <- as.numeric(f_hivprev_15$hivprev15)
  f_hivprev_15$hivprev15 <- f_hivprev_15$hivprev15/100
  
  #15-24 and 15-49 values
  f_hivprev_1549 <- read.csv("f_hivprev_1549.csv")
  toMatch <- c("Footnote", "lower", "upper")
  f_hivprev_1549 <- f_hivprev_1549[,-grep(paste(toMatch, collapse = "|"), colnames(f_hivprev_1549))]
  colnames(f_hivprev_1549) <- str_replace(colnames(f_hivprev_1549), "X", "")
  f_hivprev_1549 <- pivot_longer(f_hivprev_1549, cols = -1, names_to = "Year", values_to = "hivprev1549" )
  f_hivprev_1549$Country[f_hivprev_1549$Country == "United Republic of Tanzania"] <- "Tanzania"
  f_hivprev_1549$Country[f_hivprev_1549$Country == "Congo"] <- "Congo (Rep)"
  f_hivprev_1549$Country[f_hivprev_1549$Country == "Cabo Verde"] <- "Cape Verde"
  f_hivprev_1549$hivprev1549[f_hivprev_1549$hivprev1549 == "<0.1 "] <- 0.1
  f_hivprev_1549$Year <- as.numeric(f_hivprev_1549$Year)
  f_hivprev_1549$hivprev1549 <- as.numeric(f_hivprev_1549$hivprev1549)
  f_hivprev_1549$hivprev1549 <- f_hivprev_1549$hivprev1549/100
  
  f_hivprev_1524 <- read.csv("f_hivprev_1524.csv")
  toMatch <- c("Footnote", "lower", "upper")
  f_hivprev_1524 <- f_hivprev_1524[,-grep(paste(toMatch, collapse = "|"), colnames(f_hivprev_1524))]
  colnames(f_hivprev_1524) <- str_replace(colnames(f_hivprev_1524), "X", "")
  f_hivprev_1524 <- pivot_longer(f_hivprev_1524, cols = -1, names_to = "Year", values_to = "hivprev1524" )
  f_hivprev_1524$Country[f_hivprev_1524$Country == "United Republic of Tanzania"] <- "Tanzania"
  f_hivprev_1524$Country[f_hivprev_1524$Country == "Congo"] <- "Congo (Rep)"
  f_hivprev_1524$Country[f_hivprev_1524$Country == "Cabo Verde"] <- "Cape Verde"
  f_hivprev_1524$hivprev1524[f_hivprev_1524$hivprev1524 == "<0.1 "] <- 0.1
  f_hivprev_1524$Year <- as.numeric(f_hivprev_1524$Year)
  f_hivprev_1524$hivprev1524 <- as.numeric(f_hivprev_1524$hivprev1524)
  f_hivprev_1524$hivprev1524 <- f_hivprev_1524$hivprev1524/100
  
  #trying to obtain values for 25-49 using using 15-24 and 25-49
  pop_age <- read_excel("pop_age0020.xlsx", sheet = 2)
  colnames(pop_age) <- pop_age[1,]
  pop_age <- pop_age[-c(1:4), -c(1, 3:5)]
  pop_age[,-1] <- lapply(pop_age[,-1], as.numeric)
  
  f_hivprev <- left_join(f_hivprev_15, f_hivprev_1549)

  #add hiv prevalence to surveys that didn't have hiv biomarker
  dat_timetrends_gr <- left_join(dat_timetrends_gr, det_avail[,c("df", "HIV Biomarker")])

  df_imp <- dat_timetrends_gr[dat_timetrends_gr$'HIV Biomarker' == "N" | is.na(dat_timetrends_gr$hivstat),] #includes indv from surveys with biomarkers but don't have data on their biomarker
  df_imp <- (left_join(df_imp, f_hivprev))

  df_noimp <- dat_timetrends_gr[dat_timetrends_gr$'HIV Biomarker' != "N" & !is.na(dat_timetrends_gr$hivstat),]
  df_noimp$hivprev15 <- df_noimp$hivstat
  df_noimp$hivprev1549 <- df_noimp$hivstat

  df_cmb <- rbind(df_imp, df_noimp)
  
  saveRDS(df_cmb, file = "hiv_dat_timetrends_gr.rds")
}

# read in the data
hiv_df_gr <- readRDS("hiv_dat_timetrends_gr.rds")
#df_gr <- readRDS("dat_timetrends_gr.rds")
sum(hiv_df_gr$denominator_original)
# ---- B) Model ----
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
  
  //int<lower = 1> G; //number of gni groups - 1
  //matrix[G, N] gni;
  
  //vector[N] pgm; //binary variable if national screening program exists
  
  int<lower = 1> n_hiv; // number of hiv groups
  matrix[n_hiv, N] prv_wgt;
  
  //int<lower = 0, upper = 1> whs[N];
  
  vector[N] time;
  matrix[A + 1, N] time_length;
  
  //grouped regression specific inputs
  int<lower = 0> num[N]; // numerator
  int<lower = 1> den[N]; // denominator
  
  //for correlated intercept and slopes
  //int K;
  //matrix[N, K] time;
}

parameters {
  real alpha; //intercept
  
  vector[C] re_country; 
  real<lower = 0, upper = 10> sd_re_country;
  vector[S] re_survey; 
  real<lower = 0, upper = 10> sd_re_survey;
  vector[R] re_region; 
  real<lower = 0, upper = 10> sd_re_region;
  
  //real beta_time;
  vector[R] rs_time;
  real<lower = 0, upper = 10> sd_rs_region;
  real slope_time; 
  
  //vector[C] rs_pgm; 
  //real<lower = 0, upper = 10> sd_rs_country_pgm;
  //real slope_pgm;
  
  row_vector[A] beta_age;
  //row_vector[G] beta_gni;
  
  //real beta_recall_hiv;
  
  real<lower = -5, upper = 5> beta_hiv;
  vector[C] rs_hiv_country;
  vector[R] rs_hiv_region;
  real<lower = 0, upper = 10> sd_rs_country_hiv;
  real<lower = 0, upper = 10> sd_rs_region_hiv;
  
  //real beta_whs;
  //real beta_pgm;
  
  //real beta_time_length;
  row_vector[A + 1] beta_time_length;
  
  //for varying slopes and correlation
  //real mu_a;
  //real mu_b;
  //corr_matrix[K] Omega;
  //vector<lower=0>[K] tau;
}

transformed parameters{
  vector[N] logit_prd;
  vector[N] logit_prob;  
  vector[n_hiv] hiv_effect;
  matrix[N, n_hiv] prd;
  
  //vector[2] B[C];
  //vector[2] Mu;
  //Mu[1] = mu_a;
  //Mu[2] = mu_b;
  
  hiv_effect[1] = 0;
  hiv_effect[2] = 1;
  //hiv_effect[2] = beta_hiv;
  for(n in 1:N) {
    logit_prd[n] = alpha + re_region[region[n]] + re_country[country[n]] + re_survey[survey[n]] 
    + rs_time[region[n]] * time[n] 
    //+ rs_pgm[country[n]] * pgm[n]
    //+ beta_time * time[n] 
    + beta_age * age[, n]
    //+ beta_gni * gni[, n]
    //+ beta_pgm * pgm[n]
    //+ beta_whs * whs[n]
    + beta_time_length * time_length[, n];
    for (h in 1:n_hiv) {
      prd[n, h] = inv_logit(logit_prd[n] + (beta_hiv + rs_hiv_region[region[n]] + rs_hiv_country[country[n]]) * hiv_effect[h]);
      //+ beta_recall_hiv * time_length[n] * hiv_effect[h]; 
    }
    // if we have no information on HIV, we do a weighted average
    logit_prob[n] = logit(prd[n, ] * prv_wgt[, n]);
  }
}

model {
  // priors
  
  // Varying-slope and intercept by country
  //tau ~ cauchy(0,2.5);
  //Omega ~ lkj_corr(2);
  //mu_a ~ normal(0,10);
  //mu_b ~ normal(0,10);
  //B ~ multi_normal(Mu, quad_form_diag(Omega, tau));
  
  // Overall Intercept
  alpha ~ normal(0, 10);
  
  // Random Intercepts by survey
  re_survey ~ normal(0, sd_re_survey); 
  sd_re_survey ~ cauchy(0, 2) T[0, 10]; // half-Chauchy prior
  
  // Random Intercepts by country 
  re_country ~ normal(0, sd_re_country); 
  sd_re_country ~ cauchy(0, 2) T[0, 10]; // half-Chauchy prior
  
    // Random Intercepts by region
  re_region ~ normal(0, sd_re_region); 
  sd_re_region ~ cauchy(0, 2) T[0, 10]; // half-Chauchy prior
  
  // Random Slopes
  rs_time ~ normal(slope_time, sd_rs_region); 
  slope_time ~ normal(0, 10);
  sd_rs_region ~ cauchy(0, 5) T[0, 10];
  
  beta_hiv ~ normal(0, 10);
  rs_hiv_region ~ normal(0, sd_rs_region_hiv); 
  sd_rs_region_hiv ~ cauchy(0, 5) T[0, 10];
  
  rs_hiv_country ~ normal(0, sd_rs_country_hiv); 
  sd_rs_country_hiv ~ cauchy(0, 5) T[0, 10];
  
  //rs_pgm ~ normal(slope_time, sd_rs_country_pgm); 
  //slope_pgm ~ normal(0, 10);
  //sd_rs_country_pgm ~ cauchy(0, 1) T[0, 10];
  
  // Fixed Effects
  beta_age ~ normal(0, 5);
  beta_time_length ~ normal(0, 2);
  //beta_time ~ normal(0, 2);
  //beta_whs ~ normal(0, 5);
  //beta_recall_hiv ~ normal(0, 0.1);
  //beta_gni ~ normal(0, 5);
  //beta_hiv ~ normal(0, 5);
  //beta_pgm ~ normal(0, 5);
  
  //likelihood
  num ~ binomial_logit(den, logit_prob);  
}

generated quantities {
  vector[N] log_lik;
  vector[N] pred_num;
  vector[N] pred_prob;
  //matrix[K, K] Sigma;
  //Sigma = quad_form_diag(Omega, tau);
  
  for (n in 1:N) { //can try to vectorize all - didnt work
    //generates log-likelihood for model validation
    log_lik[n] = binomial_logit_lpmf(num[n] | den[n], logit_prob[n]); 
    
    //for posterior predictive checks
    pred_num[n] = binomial_rng(den[n], inv_logit(logit_prob[n]));
    pred_prob[n] = pred_num[n]/den[n];
  }
}
'
#try a model without the hiv standardization 
# We compile the model
cc_trends_gr_stan <- rstan::stan_model(model_code = cc_trends_gr)
# ----  C) Fit ----
df_stan <- select(hiv_df_gr, -Region)

#change regions from un classification to gbd
GBD = T
if(GBD == T){
  region_country <- unique(surveys[,c("Country", "Region")]) #either un or gbd classification was chosen straight in the excel document
  df_stan <- left_join(df_stan, region_country)
  #df_stan <- filter(df_stan, Country != "Mauritius")
}

#remove values with >20% missing data
# missing <- read_excel("CC Analysis.xlsx", sheet = "Missing")
# missing$agegr[grepl("whs", missing$df) & missing$agegr == "65-69"] <- "65-1000" #this is done bc question only asked to those 15-69
# df_stan <- left_join(df_stan, missing)
# df_stan$Exclude20[is.na(df_stan$Exclude20)] <- "N" #for the added surveys with tabulations
# df_stan$Exclude15[is.na(df_stan$Exclude15)] <- "N" #for the added surveys with tabulations
# #seeing what model looks like if I exclude WHS surveys for those 40+
# df_stan$Exclude20[df_stan$Survey == "WHS" & df_stan$agegr != "15-19" & df_stan$agegr != "20-24" & df_stan$agegr != "25-29" & df_stan$agegr != "30-34" & df_stan$agegr != "35-39"] <- "Y"
# df_stan <- filter(df_stan, Exclude20 == "N")

#restrict to only those with hiv biomarkers
#df_stan <- filter(df_stan, !is.na(hivstat))

#restrict analysis by age
#15-29
#df_stan <- filter(df_stan, agegr == "15-19" | agegr == "20-24" | agegr == "25-29")
#25-49
df_stan <- filter(df_stan, agegr == "25-29" | agegr == "30-34" | agegr == "35-39" | agegr == "40-44" | agegr == "45-49")
#50-65+
#df_stan <- filter(df_stan, agegr == "50-54" | agegr == "55-59" | agegr == "60-64" | agegr == "65-1000")
#15-49
#df_stan <- filter(df_stan, agegr != "50-54" & agegr != "55-59" & agegr != "60-64" & agegr != "65-1000")

#restrict analysis by region
#df_stan <- filter(df_stan, Region == "Central" | Region == "Western")
#remove values with a denominator of 0
df_stan <- filter(df_stan, den != 0)

#restrict to countries with more than one survey
# unique(df_stan[order(df_stan$Country),c("Country", "Year")])
# df_stan <- filter(df_stan,
#                   Country == "Benin" |
#                     Country == "Burkina Faso" |
#                     Country == "Cote d'Ivoire" |
#                     Country == "Ethiopia" |
#                     Country == "Kenya" |
#                     Country == "Lesotho" |
#                     Country == "Malawi" |
#                     Country == "Namibia" |
#                     Country == "Senegal" |
#                     Country == "South Africa" |
#                     Country == "Swaziland" |
#                     Country == "Zambia" |
#                     Country == "Zimbabwe")

# add gdp per capita + national screening pgm to df
df_stan <- left_join(df_stan, eco_det)

#prep variables for stan model
N <- nrow(df_stan)

S <- length(unique(df_stan$df))
survey <- as.numeric(factor(df_stan$df))

C <- length(unique(df_stan$Country))
df_stan$count_num <- as.numeric(factor(df_stan$Country))
country <- df_stan$count_num #levels will be in alphabetical order

df_stan$Region[df_stan$Region == "Central" | df_stan$Region == "Western"] <- "Central/Western"
#df_stan$Region[df_stan$Country == "Zambia" | df_stan$Country == "Zimbabwe"] <- "Southern"
R <- length(unique(df_stan$Region))
df_stan$reg_num <- as.numeric(factor(df_stan$Region))
region <- df_stan$reg_num 

# df_stan$whs[df_stan$Survey == "WHS"] <- 1
# df_stan$whs[df_stan$Survey != "WHS"] <- 0
# whs <- df_stan$whs

#create age matrix
n_age <- length(unique(df_stan$agegr))
df_stan$agegr_num <- as.numeric(factor(df_stan$agegr))
for(a in 2:n_age){
  df_stan$x[df_stan$agegr_num == a] <- 1
  df_stan$x[df_stan$agegr_num != a] <- 0
  name <- paste0("agegr", a)
  names(df_stan)[names(df_stan) == "x"] <- name
}

age = t(as.matrix(df_stan[,c("agegr2", "agegr3", "agegr4" , "agegr5")]))
                          #,"agegr6" , "agegr7")]))
                          # , "agegr8", "agegr9", "agegr10", "agegr11")]))
A <- n_age - 1 #subtract one for the model as there will be n-1 parameters

#center year around 2000 (lowest survey year)
df_stan$year_cent <- df_stan$Year - 2000
time <- df_stan$year_cent

#create binary variable for time length, 1 = 3y, 0 = ever - will keep it at 1 and 2 for one run bc it was ran at 1 and 2 to obtain the odds ratio
df_stan$time_length2[df_stan$time_length == "3y"] <- 1
df_stan$time_length2[df_stan$time_length == "ever"] <- 0
#time_length <- df_stan$time_length2

#creating time_length matrix if trying to have a different odds ratio for screening interval by age
n_age <- length(unique(df_stan$agegr))
for(a in 1:n_age){
  df_stan$x[df_stan$time_length2 == 1 & df_stan$agegr_num == a] <- 1
  df_stan$x[df_stan$time_length2 == 1 & df_stan$agegr_num != a] <- 0
  df_stan$x[df_stan$time_length2 == 0] <- 0
  name <- paste0("time_length_agegr", a)
  names(df_stan)[names(df_stan) == "x"] <- name
}
#View(df_stan[,c("agegr_num", "time_length", "time_length_agegr1", "time_length_agegr2", "time_length_agegr3", "time_length_agegr4", "time_length_agegr5")])

time_length = t(as.matrix(df_stan[,c("time_length_agegr1", "time_length_agegr2", "time_length_agegr3", "time_length_agegr4", "time_length_agegr5")]))

# #add gdp to df_stan
# gni_raw <- read.csv("gni.csv")
# gni <- pivot_longer(gni_raw, cols = -c(1:4), names_to = "Year", values_to = "gni")
# gni$Year <- as.numeric(sub(".", "", gni$Year))
# names(gni)[names(gni) == "Country.Name"] <- "Country"
# gni <- gni[,-c(2:4)]
# gni$Country[gni$Country == "Cabo Verde"] <- "Cape Verde"
# gni$Country[gni$Country == "Congo, Rep."] <- "Congo (Rep)"
# df_stan <- left_join(df_stan, gni)
# #create income group variable
# gni_thresholds <- read_excel("gni_thresholds.xlsx")
# gni_thresholds <- pivot_wider(gni_thresholds, names_from = Class, values_from = Max_Threshold)
# df_stan <- left_join(df_stan, gni_thresholds)
# df_stan$income_group[df_stan$gni < df_stan$Low_LowerMiddle] <- "Low"
# df_stan$income_group[df_stan$gni >= df_stan$Low_LowerMiddle & df_stan$gni <= df_stan$LowerMiddle_UpperMiddle] <- "LowerMiddle"
# df_stan$income_group[df_stan$gni > df_stan$LowerMiddle_UpperMiddle & df_stan$gni <= df_stan$UpperMiddle_High] <- "UpperMiddle"
# df_stan$income_group[df_stan$gni > df_stan$UpperMiddle_High] <- "High"
# df_stan$income_group_num <- as.numeric(factor(df_stan$income_group, levels = c("Low", "LowerMiddle", "UpperMiddle", "High")))
# 
# # #create income group matrix
# G <- length(unique(df_stan$income_group_num))
# for(g in 2:G){
#   df_stan$x[df_stan$income_group_num == g] <- 1
#   df_stan$x[df_stan$income_group_num != g] <- 0
#   name <- paste0("income_group", g)
#   names(df_stan)[names(df_stan) == "x"] <- name
# }
# 
# G <- G-1 #subtract one for the model as there will be n-1 parameters
# 
# gni <- t(as.matrix(df_stan[,c("income_group2", "income_group3")]))

#national screening program
# df_stan$pgm[df_stan$`National Screening Program` == "N"] <- 0
# df_stan$pgm[df_stan$`National Screening Program` == "Y"] <- 1
# any_scr_pgm1 <- filter(screening_program, pgm_type == "Any")
# for(i in 1:nrow(df_stan)){
#   c <- df_stan$Country[i]
#   pgm1 <- filter(any_scr_pgm1, Country == c)
#   if(nrow(pgm1) == 1){
#     if(df_stan$Year[i] >= pgm1$screen_start & df_stan$Year[i] <= pgm1$screen_end){
#       df_stan$pgm1[i] <- 1
#     }else{
#       df_stan$pgm1[i] <- 0
#     }
#   }else if(nrow(pgm1) == 2){
#     if(df_stan$Year[i] >= pgm1$screen_start[1] & df_stan$Year[i] <= pgm1$screen_end[1]){
#       df_stan$pgm1[i] <- 1
#     }else if(df_stan$Year[i] >= pgm1$screen_start[2] & df_stan$Year[i] <= pgm1$screen_end[2]){
#       df_stan$pgm1[i] <- 1
#     }else{
#       df_stan$pgm1[i] <- 0
#     }
#   }else{
#     df_stan$pgm1[i] <- 0
#   }
# }
#screening program time lag
# time_lag_start <- 1
# #time_lag_end <- 3
# for(i in 1:nrow(df_stan)){ #to intro a time lag to when a screening program was implemented
#   c <- df_stan$Country[i]
#   pgm <- filter(any_scr_pgm, Country == c)
#   if(nrow(pgm) > 0){
#     for(j in 1:nrow(pgm)){
#       if(df_stan$Year[i] >= pgm$screen_start[j] & df_stan$Year[i] < (pgm$screen_start[j] + time_lag_start)){
#         df_stan$pgm[i] <- 0
#       }
#     }
#   }
# }
# pgm <- as.numeric(df_stan$pgm)
# View(df_stan[,c("Country","Year", "National Screening Program","pgm")])

#hiv values
#hiv <- df_stan$hivstat

#group logistic reg specific variables
num <- df_stan$num
den <- df_stan$den

#variables for hiv prevalence
prv_wgt <- matrix(data = NA, nrow = N, ncol = 2)
prv_wgt[, 1] <- 1 - df_stan$hivprev1549 # hiv neg
prv_wgt[, 2] <- df_stan$hivprev1549     # hiv pos
n_hiv <- ncol(prv_wgt)

data_stan <- list(N = N,
                  S = S,
                  survey = survey,
                  C = C,
                  country = country,
                  R = R,
                  region = region,
                  A = A,
                  # gni = gni,
                  # G = G,
                  # pgm = pgm,
                  # whs = whs,
                  age = age,
                  time = time,
                  time_length = time_length,
                  n_hiv = n_hiv,
                  prv_wgt = t(prv_wgt),
                  # hiv = hiv,
                  num = num,
                  den = den)

# ----  D) Run ----
options(mc.cores = parallel::detectCores())
fit <- rstan::sampling(cc_trends_gr_stan, 
                       data = data_stan, 
                       iter = 3000, chains = 6, refresh = 50,
                       warmup = 1500, thin = 1,
                       control = list(adapt_delta = 0.99, max_treedepth = 15))
save(data_stan, fit, df_stan, cc_trends_gr_stan, file = "gr_combined73")
rstan::rstan_options(auto_write = TRUE) 

rstan::stan_trace(fit, pars = c("alpha", "slope_time", "rs_time", "sd_rs_region", "re_country", "sd_re_country"))
rstan::stan_trace(fit, pars = c("beta_age"))
rstan::stan_trace(fit, pars = c("beta_time_length"))
rstan::stan_trace(fit, pars = c("rs_hiv_country"))
rstan::summary(fit, pars = c("alpha", "slope_time", "rs_time", "sd_rs_region", "re_country", "sd_re_country"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("beta_age"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("rs_hiv"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("beta_time_length"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("re_survey"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("re_country"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("rs_hiv"), probs = c(0.025, 0.5, 0.975))$summary
rstan::summary(fit, pars = c("rs_pgm"), probs = c(0.025, 0.5, 0.975))$summary
#rstan::summary(fit, pars = c("beta_gni"), probs = c(0.025, 0.5, 0.975))$summary
#load("gr_combined73")
#odds ratios for the effects of HIV
draws <- rstan::extract(fit)
rs_hiv_overall <- draws$beta_hiv
rs_hiv_region <- draws$rs_hiv_region
rs_hiv_country <- draws$rs_hiv_country

exp(quantile(rs_hiv_overall, probs = c(0.025, 0.5, 0.975)))

hiv_region_or <- matrix(data = NA, ncol = 3, nrow = 3)
for(i in 1:ncol(rs_hiv_region)){
  a <- rs_hiv_region[i] + rs_hiv_overall
  hiv_region_or[i,] <- exp(quantile(a, probs = c(0.025, 0.5, 0.975)))
}

df_rs_hiv <- data.frame(Country = unique(df_stan$Country), 
                        count_num = unique(df_stan$count_num),
                        t(rs_hiv_country[unique(df_stan$count_num),]))

df_rs_hiv <- left_join(df_rs_hiv, region_country)
df_rs_hiv$Region[df_rs_hiv$Region == "Central" | df_rs_hiv$Region == "Western"] <- "Central/Western"
df_rs_hiv$reg_num <- as.numeric(factor(df_rs_hiv$Region))

hiv_country_or <- matrix(data = NA, ncol = 3, nrow = 28)
for(i in 1:nrow(df_rs_hiv)){
  a <- rs_hiv_overall + rs_hiv_region[,df_rs_hiv$reg_num[i]] + rs_hiv_country[,df_rs_hiv$count_num[i]]
  hiv_country_or[i,] <- exp(quantile(a, probs = c(0.025, 0.5, 0.975)))
}

hiv_country_or <- data.frame(hiv_country_or)
hiv_country_or$Country <- df_rs_hiv$Country
hiv_country_or$Region <- df_rs_hiv$Region
hiv_countries <- unique(df_stan[!is.na(df_stan$hivstat),]$Country)

View(hiv_country_or[hiv_country_or$Country %in% hiv_countries,])

hiv_odds <- paste0(round(df_rs_hiv_only[,2], 1), " (", round(df_rs_hiv_only[,1], 1), "-", round(df_rs_hiv_only[,3], 1), ")")
hiv_output <- data.frame(df_rs_hiv_only[,4:5], hiv_odds)
hiv_output <- hiv_output[order(hiv_output$Country),]
names(hiv_output)[2] <- "Odds Ratio (95% CrI)"
rownames(hiv_output) <- NULL

#odds ratios for screening intervals
exp(rstan::summary(fit, pars = c("beta_time_length"), probs = c(0.025, 0.5, 0.975))$summary)

hiv_output %>%
  kbl(caption = "Table 2. Lifetime screening odds ratios among women living and not living with HIV", align = "c") %>%
  kable_classic() %>%
  row_spec(0, bold = T)

# ----  E) Predictions ----
#df_stan_org <- df_stan
df_stan$predicted <- "Y"
region_country$Region[region_country$Region == "Central" | region_country$Region == "Western"] <- "Central/Western"
df_stan_all_count <- left_join(region_country, df_stan)
df_stan_all_count$predicted[is.na(df_stan_all_count$predicted)] <- "N"
source("predict_poststrat_func.R")

#finding Countries that have x or more surveys
# a <- unique(df_stan[,c("Country", "df")])
# count_2p <- rownames(filter(data.frame(num = table(a$Country) >= 2), num == T)) #2 or more surveys
# count_3p <- rownames(filter(data.frame(num = table(a$Country) >= 3), num == T)) #3 or more surveys

##find weights by country
df_weights_hiv <- weights(agegr = c("25-29", "30-34", "35-39", "40-44", "45-49"), hiv = T, years = c(2000:2020))

#Obtain unweighted predicted probabilities for each iteration and category. This can be done in generated quantities however I'm going to do it outside generated quantities 
#create data set to predict
#load("gr_combined63")
draws <- rstan::extract(fit)

countries_2p <- c("Benin", "Burkina Faso", "Cote d'Ivoire", "Eswatini", "Ethiopia", "Ghana", "Kenya", "Lesotho", "Malawi", "Namibia", "Senegal", "South Africa", "Zambia", "Zimbabwe")

pred_prob_weighted_2549 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(0, 5, 10, 15, 20), agegr = 1:count(!is.na(unique(df_stan$agegr))), 
                                    hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)
pred_prob_weighted_3049 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(0, 5, 10, 15, 20), agegr = 2:5, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)
# pred_prob_weighted_2549_2 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(20), agegr = 1:count(!is.na(unique(df_stan$agegr))), 
#                                          hivstat_log = T, gni_log = F, pgm_log = T, multi_time_length = T, df_weights_hiv)
pred_prob_weighted_3034 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(20), agegr = 2, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)
#store<-pred_prob_weighted_2549
#now that we a matrix of weighted values we can then group the data by their characteristics and then obtain weighted means/quantiles
prob_by_year_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 1, "Year")
prob_by_year_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 1, "Year")
prob_by_reg_year_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Region", "Year")
prob_by_reg_year_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Region", "Year")
prob_by_count_year_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Country", "Year")
prob_by_count_year_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Country", "Year")

#by hiv status
prob_by_hiv_year_3y <- prob_final(pred_prob_weighted_2549[[1]], dat_pred = pred_prob_weighted_2549[[3]], 2, "hivstat", "Year", count_restrict = F, countries = hiv_countries)
prob_by_hiv_year_ever <- prob_final(pred_prob_weighted_2549[[2]], dat_pred = pred_prob_weighted_2549[[3]], 2, "hivstat", "Year", count_restrict = F, countries = hiv_countries)
prob_by_hiv_reg_year_3y <- prob_final(pred_prob_weighted_2549[[1]], dat_pred = pred_prob_weighted_2549[[3]], 3, "hivstat", "Region", "Year", count_restrict = F, countries = hiv_countries)
prob_by_hiv_reg_year_ever <- prob_final(pred_prob_weighted_2549[[2]], dat_pred = pred_prob_weighted_2549[[3]], 3, "hivstat", "Region", "Year", count_restrict = F, countries = hiv_countries)
prob_by_hiv_count_year_3y <- prob_final(pred_prob_weighted_2549[[1]], dat_pred = pred_prob_weighted_2549[[3]], 3, "hivstat", "Country", "Year")
prob_by_hiv_count_year_ever <- prob_final(pred_prob_weighted_2549[[2]], dat_pred = pred_prob_weighted_2549[[3]], 3, "hivstat", "Country", "Year")

#by age group
prob_by_year_age_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))
prob_by_year_age_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 2, "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))
prob_by_reg_year_age_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 3, "Region", "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))
prob_by_reg_year_age_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 3, "Region", "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))
prob_by_count_year_age_3y <- prob_final(pred_prob_weighted_3049[[1]], dat_pred = pred_prob_weighted_3049[[3]], 3, "Country", "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))
prob_by_count_year_age_ever <- prob_final(pred_prob_weighted_3049[[2]], dat_pred = pred_prob_weighted_3049[[3]], 3, "Country", "Year", "agegr", n_agegr = unique(pred_prob_weighted_3049[[3]]$agegr), agegr_labels = c("30-34", "35-39", "40-44", "45-49"))

prob_by_hiv_reg_year_age_3y <- prob_final(pred_prob_weighted_3y, dat_pred, 4, "hivstat", "Region", "Year", "agegr")
prob_by_hiv_reg_year_age_ever <- prob_final(pred_prob_weighted_ever, dat_pred, 4, "hivstat", "Region", "Year", "agegr")
prob_by_hiv_count_year_age_3y <- prob_final(pred_prob_weighted_3y, dat_pred, 4, "hivstat", "Country", "Year", "agegr")
prob_by_hiv_count_year_age_ever <- prob_final(pred_prob_weighted_ever, dat_pred, 4, "hivstat", "Country", "Year", "agegr")

#combined ever and past 3y dataframes
prob_by_year <- dat_combine(prob_by_year_ever, prob_by_year_3y)
prob_by_reg_year <- dat_combine(prob_by_reg_year_ever, prob_by_reg_year_3y)
prob_by_count_year <- dat_combine(prob_by_count_year_ever, prob_by_count_year_3y)

prob_by_hiv_year <- dat_combine(prob_by_hiv_year_ever, prob_by_hiv_year_3y)
prob_by_hiv_reg_year <- dat_combine(prob_by_hiv_reg_year_ever, prob_by_hiv_reg_year_3y)
prob_by_hiv_count_year <- dat_combine(prob_by_hiv_count_year_ever, prob_by_hiv_count_year_3y)

prob_by_reg_year_age <- dat_combine(prob_by_reg_year_age_ever, prob_by_reg_year_age_3y)
prob_by_count_year_age <- dat_combine(prob_by_count_year_age_ever, prob_by_count_year_age_3y)
prob_by_hiv_reg_year_age <- dat_combine(prob_by_hiv_reg_year_age_ever, prob_by_hiv_reg_year_age_3y)
prob_by_hiv_count_year_age <- dat_combine(prob_by_hiv_count_year_age_ever, prob_by_hiv_count_year_age_3y)

#create summary dataframes
year = 2020
agegr = "45-49"
hivstat = "HIV+"

a <- prob_by_year_age_3y[prob_by_year_age_3y$Year == year & prob_by_year_age_3y$agegr == agegr,]
b <- prob_by_count_year_age_3y[prob_by_count_year_age_3y$Year == year & prob_by_count_year_age_3y$Region == "Southern" & prob_by_count_year_age_3y$agegr == agegr,]
c <- prob_by_count_year_age_3y[prob_by_count_year_age_3y$Year == year & prob_by_count_year_age_3y$Region == "Central/Western" & prob_by_count_year_age_3y$agegr == agegr ,]
d <- prob_by_count_year_age_3y[prob_by_count_year_age_3y$Year == year & prob_by_count_year_age_3y$Region == "Eastern" & prob_by_count_year_age_3y$agegr == agegr,]
e <- prob_by_reg_year_age_3y[prob_by_reg_year_age_3y$Year == year & prob_by_reg_year_age_3y$agegr == agegr,]

overall <- data.frame(Country = "Overall", Region = "Overall", screening = paste0(round(a[,"50%"]*100, 0), "% (", round(a[,"2.5%"]*100, 0), "-", round(a[,"97.5%"]*100, 0), "%)"))
sa <- data.frame(Country = b[,"Country"], Region = "Southern", screening = paste0(round(b[,"50%"]*100, 0), "% (", round(b[,"2.5%"]*100, 0), "-", round(b[,"97.5%"]*100, 0), "%)"))
cwa <- data.frame(Country = c[ ,"Country"], Region = "Central/Western", screening = paste0(round(c[ ,"50%"]*100, 0), "% (", round(c[ ,"2.5%"]*100, 0), "-", round(c[ ,"97.5%"]*100, 0), "%)"))
ea <- data.frame(Country = d[ ,"Country"], Region = "Eastern", screening = paste0(round(d[ ,"50%"]*100, 0), "% (", round(d[ ,"2.5%"]*100, 0), "-", round(d[ ,"97.5%"]*100, 0), "%)"))
region <- data.frame(Country = "Overall", Region = e[ ,"Region"], screening = paste0(round(e[ ,"50%"]*100, 0), "% (", round(e[ ,"2.5%"]*100, 0), "-", round(e[ ,"97.5%"]*100, 0), "%)"))

df <- rbind(overall, 
            region[region$Region == "Central/Western",], cwa,
            region[region$Region == "Eastern",], ea,
            region[region$Region == "Southern",], sa)

countries_2p <- c("Benin", "Burkina Faso", "Cote d'Ivoire", "Eswatini", "Ethiopia", "Ghana", "Kenya", "Lesotho", "Malawi", "Namibia", "Senegal", "South Africa", "Zambia", "Zimbabwe")
View(df[df$Country %in% countries_2p | df$Country == "Overall", ])

##raw probabilities from survey
raw_prob_ever <- survey_prob(timing = "ever", gen = T, hiv = T, df_to_add)
raw_prob_3y <- survey_prob(timing = "3y", gen = T, hiv = T, df_to_add)

raw_prob <- dat_combine(raw_prob_ever[[1]], raw_prob_3y[[1]])
hiv_raw_prob <- dat_combine(raw_prob_ever[[2]], raw_prob_3y[[2]])

agegr <- c("30-34", "35-39", "40-44", "45-49")
raw_prob_age <- raw_prob[raw_prob$agegr %in% agegr, ]
  
raw_prob <- raw_prob[raw_prob$agegr == "30-49",]
hiv_raw_prob <- hiv_raw_prob[hiv_raw_prob$agegr == "25-49",]
# ----  F) Plots ----
#p_overall <- 
  ggplot(prob_by_year, aes(x = Year, y = mean)) +
  #geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  ylim(0, 1) +
  ylab("Proportion") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  scale_color_nejm(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  scale_fill_nejm(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime"))

#p_hiv_overall <- 
ggplot(prob_by_hiv_year[prob_by_hiv_year$time_length == "ever", ], aes(x = Year, y = mean)) +
  geom_pointrange(data = hiv_raw_prob[hiv_raw_prob$time_length == "ever", ], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = hivstat), size = 0.2) +
  geom_line(aes(colour = hivstat)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = hivstat), alpha = 0.25) +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  ylim(0, 1) +
  ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
  scale_fill_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF"))

p_count_3y <- 
  ggplot(prob_by_count_year_3y, aes(x = Year, y = mean)) +
  geom_pointrange(data = raw_prob[raw_prob$time_length == "3y", ], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = Region), size = 0.2) +
  geom_line(aes(colour = Region)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = Region), alpha = 0.25) +
  ylim(0, 1) +
  ylab("Proportion") +
  ggtitle("Time trends for CC screening in the past 3 years for those aged 30-49 by country") +
  labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
       caption = "*Points represent weighted estimates from survey data") +
  theme_minimal() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
        strip.text = element_text(size = 12)) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Country)

p_count_ever <- 
  ggplot(prob_by_count_year_ever, aes(x = Year, y = mean)) +
  geom_pointrange(data = raw_prob[raw_prob$time_length == "ever", ], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = Region), size = 0.2) +
  geom_line(aes(colour = Region)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = Region), alpha = 0.25) +
  ylim(0, 1) +
  ylab("Proportion") +
  ggtitle("Time trends for CC screening ever for for those aged 30-49 by country") +
  labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
       caption = "*Points represent weighted estimates from survey data") +
  theme_minimal() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
        strip.text = element_text(size = 12)) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Country)

p_count <- 
  ggplot(prob_by_count_year, aes(x = Year, y = mean)) +
  geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  ylim(0, 1) +
  ylab("Proportion") +
  ggtitle("Time trends for CC screening for for those aged 30-49 by country") +
  # labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
  #     caption = "*Points represent weighted estimates from survey data") +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=17),
        axis.text.x = element_text(size = 12, angle = 90, hjust=0.95, vjust=1.25),
        strip.text = element_text(size = 12)) +
  scale_colour_brewer(name = "Recall Period", palette = "Set2", labels = c("3Y", "Ever")) +
  scale_fill_brewer(name = "Recall Period", palette = "Set2", labels = c("3Y", "Ever")) +
  facet_wrap(~Country)

#countries with 2 or more surveys
countries_2p <- c("Benin", "Burkina Faso", "Cote d'Ivoire", "Eswatini", "Ethiopia", "Ghana", "Kenya", "Lesotho", "Malawi", "Namibia", "Senegal", "South Africa", "Zambia", "Zimbabwe")
#p_count_2p <- 
  ggplot(prob_by_count_year[prob_by_count_year$Country %in% countries_2p & prob_by_count_year$Region == "Southern",], aes(x = Year, y = mean)) +
  geom_pointrange(data = raw_prob[raw_prob$Country %in% countries_2p & raw_prob$Region == "Southern",], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  ylim(0, 1) +
  ylab("Proportion") +
  #ggtitle("Time trends for CC screening for for those aged 30-49 by country") +
  #labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
       #caption = "*Points represent weighted estimates from survey data") +
  theme_bw() +
  theme(axis.text=element_text(size = 15),
        axis.title=element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        strip.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  scale_color_nejm(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  scale_fill_nejm(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  facet_grid(Region~Country)

#p_reg <- 
  ggplot(prob_by_reg_year, aes(x = Year, y = mean)) +
  #geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
  ylim(0, 1) +
  ylab("Proportion") +
  # ggtitle("Time trends for CC screening for for those aged 30-49 by region") +
  # labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
  #      caption = "*Points represent weighted estimates from survey data") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        strip.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  scale_colour_nejm(name = "Recall Period", labels = c("Within the past 3Y", "Lifetime")) +
  scale_fill_nejm(name = "Recall Period", labels = c("Within the past 3Y", "Lifetime")) +
  facet_wrap(~Region)
  
  
hiv_countries_ever <- c("Cameroon","Cote d'Ivoire", "Ethiopia", "Kenya", "Lesotho", "Malawi", "Namibia", "Rwanda", "South Africa", "Tanzania", "Zambia", "Zimbabwe")
hiv_countries_3y <- c("Cameroon", "Ethiopia", "Lesotho", "Malawi", "Rwanda", "South Africa", "Tanzania", "Zambia", "Zimbabwe")

hiv_count_year_ever <- prob_by_hiv_count_year_ever[prob_by_hiv_count_year_ever$Country %in% hiv_countries_ever,]
hiv_count_year_3y <- prob_by_hiv_count_year_3y[prob_by_hiv_count_year_3y$Country %in% hiv_countries_3y,]
# hiv_count_year_age_ever <- prob_by_hiv_count_year_age_ever[prob_by_hiv_count_year_age_ever$Country %in% hiv_countries_ever,]
# hiv_count_year_age_3y <- prob_by_hiv_count_year_age_3y[prob_by_hiv_count_year_age_3y$Country %in% hiv_countries_3y,]

#p_hiv_count_ever <-
  ggplot(hiv_count_year_ever[hiv_count_year_ever$Region == "Southern",], aes(x = Year, y = mean)) +
    geom_line(aes(colour = hivstat)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = hivstat), alpha = 0.25) +
    #geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
    geom_pointrange(data = hiv_raw_prob[hiv_raw_prob$time_length == "ever" & hiv_raw_prob$Region == "Southern",], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = hivstat), size = 0.2) +
    ylim(0, 1) +
    ylab("Proportion") +
  #   ggtitle("Time trends for CC screening ever by Country and HIV Status for those aged 25-49") +
  # labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
  #      caption = "*Points represent weighted estimates from survey data") +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
          strip.text = element_text(size = 15),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.text=element_text(size=15),
          legend.title=element_text(size=15)) +
    scale_colour_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
    scale_fill_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
    facet_grid(Region~Country)
  
p_hiv_count_3y <-
  ggplot(hiv_count_year_3y, aes(x = Year, y = mean)) +
    geom_line(aes(colour = hivstat)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = hivstat), alpha = 0.25) +
    geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
    geom_pointrange(data = hiv_raw_prob[hiv_raw_prob$time_length == "3y",], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = hivstat), size = 0.2) +
    ylim(0, 1) +
    ylab("Proportion") +
    ggtitle("Time trends for CC screening in the last 3y by Country and HIV Status for those aged 25-49") +
    labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
         caption = "*Points represent weighted estimates from survey data") +
    theme_bw() +
    theme(#axis.text=element_text(size=15),
          #axis.title=element_text(size=17),
          axis.text.x = element_text(angle = 90, hjust=0.95, vjust=1.25),
          #strip.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    scale_colour_manual(name = "HIV Status", values = c("dodgerblue3", "lightsalmon1")) +
    scale_fill_manual(name = "HIV Status", values = c("dodgerblue3", "lightsalmon1")) +
    facet_wrap(~Country)

  p_hiv_reg_ever <-
  ggplot(prob_by_hiv_reg_year[prob_by_hiv_reg_year$time_length == "ever", ], aes(x = Year, y = mean)) +
    geom_line(aes(colour = hivstat)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = hivstat), alpha = 0.25) +
    geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
    geom_pointrange(data = hiv_raw_prob[hiv_raw_prob$time_length == "ever",], aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = hivstat), size = 0.2) +
    ylim(0, 1) +
    ylab("Proportion") +
    #   ggtitle("Time trends for CC screening ever by Country and HIV Status for those aged 25-49") +
    # labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
    #      caption = "*Points represent weighted estimates from survey data") +
    theme_bw() +
    theme(#axis.text=element_text(size=15),
      #axis.title=element_text(size=17),
      axis.text.x = element_text(angle = 90, hjust=0.95, vjust=1.25),
      #strip.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      legend.position = "bottom") +
    scale_colour_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
    scale_fill_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
    facet_wrap(~Region)
#p_hiv_count_age_ever <-
  # ggplot(hiv_count_year_age_ever, aes(x = Year, y = `50%`)) +
  # geom_line(aes(colour = hivstat)) +
  # geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = hivstat), alpha = 0.25) +
  # geom_pointrange(data = hiv_age_raw_prob_ever, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = hivstat), size = 0.2) +
  # ylim(0, 1) +
  # ylab("Proportion") +
  # ggtitle("Time trends for CC screening in by Country, age group and HIV status for those aged 15-49") +
  # theme_minimal() +
  # theme(axis.text=element_text(size=15),
  #       axis.title=element_text(size=17),
  #       axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
  #       strip.text = element_text(size = 12)) +
  # scale_colour_brewer(name = "HIV Status", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
  # scale_fill_brewer(name = "HIV Status", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
  # facet_grid(Country~agegr)

p_count_age <-
  ggplot(prob_by_count_year_age, aes(x = Year, y = `50%`)) +
  geom_line(aes(colour = time_length)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  geom_pointrange(data = raw_prob_age, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  ylim(0, 1) +
  ylab("Proportion") +
    ggtitle("Time trends for CC screening in by Country and age group for those aged 30-49") +
    theme_minimal() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=17),
          axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
          strip.text = element_text(size = 12)) +
    scale_colour_brewer(name = "Recall Period", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
    scale_fill_brewer(name = "Recall Period", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
  facet_grid(Country~agegr)

p_reg_age <-
  ggplot(prob_by_reg_year_age, aes(x = Year, y = `50%`)) +
  geom_line(aes(colour = time_length)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
  geom_pointrange(data = raw_prob_age, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  ylim(0, 1) +
  ylab("Proportion") +
  ggtitle("Time trends for CC screening in by Region and age group for those aged 30-49") +
  theme_minimal() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
        strip.text = element_text(size = 12)) +
  scale_colour_brewer(name = "Recall Period", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
  scale_fill_brewer(name = "Recall Period", palette = "Set2", labels = c("Within the past 3Y", "Ever")) +
  facet_grid(Region~agegr)

p_count_3y; p_count_ever
p_count; p_count_2p; p_reg
p_hiv_count_3y; p_hiv_count_ever
p_count_age; p_reg_age
# ----  G) Save ----
save(prob_by_count_year_3y, prob_by_count_year_ever,
     prob_by_reg_year_3y, prob_by_reg_year_ever,
     prob_by_hiv_count_year_3y, prob_by_hiv_count_year_ever,
     p_overall, p_hiv_overall,
     p_count_3y, p_count_ever, p_count, p_count_2p, p_reg,
     #p_count_age, p_reg_age,
     #p_count3049, p_reg3049,
     #p_hiv_count,
     p_hiv_count_3y, p_hiv_count_ever,
     #p_hiv_count_age_ever,
     #pred_prob_3y, pred_prob_ever,
     pred_prob_weighted_2549,
     pred_prob_weighted_3049,
     #pred_prob_weighted_2549_2, 
     #pred_prob_weighted_3034,
     file = "pred_prob_gr_combined73_all_ssa")
# load("gr_combined60")
# load("pred_prob_gr_combined73_all_ssa")

# ----  H) Model validation ----
## Model validation 
#load("gr_combined73")
#WAIC
ll <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = T)
waic_gr_combined73 <- waic(ll)
#loo
ll_2 <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = F)
r_eff <- relative_eff(exp(ll_2), cores = 6) 
loo_gr_combined73 <- loo(ll_2, r_eff = r_eff, cores = 6)
print(waic_gr_combined73); print(loo_gr_combined73)

# model_val <- data.frame(Random_intercepts = rep(" Country, Survey", 6),
#                         Random_slopes = c(rep("TimeXRegion", 3), rep("TimeXRegion, National Screening ProgramXCountry", 3)),
#                         Fixed_effects = c(rep("Age, National Screening Program", 3), rep("Age", 3)),
#                         Age = c(rep("25-49", 6)),
#                         Time_lag = rep(c("0y", "1y", "3y"), 2),
#                         WAIC = c(waic_gr_combined45$estimates[3,1], waic_gr_combined46$estimates[3,1], waic_gr_combined47$estimates[3,1], waic_gr_combined48$estimates[3,1],
#                                  waic_gr_combined49$estimates[3,1], waic_gr_combined50$estimates[3,1]),
#                         LOOIC = c(loo_gr_combined45$estimates[3,1], loo_gr_combined46$estimates[3,1], loo_gr_combined47$estimates[3,1], loo_gr_combined48$estimates[3,1],
#                                   loo_gr_combined49$estimates[3,1], loo_gr_combined50$estimates[3,1]))
# 
# kbl(model_val) %>%
#   kable_styling()

##In sample comparisons
#obtain values from model
draws <- rstan::extract(fit)
pred_prob <- data.frame(t(draws$pred_prob))
pred_prob <- row_sumstats(pred_prob, prob = c(0.5, 0.025, 0.975))
pred_prob <- cbind(df_stan, pred_prob[,c("50%", "2.5%", "97.5%")])
pred_prob$hivstat <- factor(pred_prob$hivstat, levels = c(0, 1), labels = c("HIV-", "HIV+"))

#values from survey
source("predict_poststrat_func.R")
list_raw_prob_ever <- survey_prob(timing = "ever", gen = T, hiv = T, df_to_add)
list_raw_prob_ever[[1]]$hivstat = NA
raw_prob_ever <- rbind(list_raw_prob_ever[[1]], list_raw_prob_ever[[2]])

list_raw_prob_3y <- survey_prob(timing = "3y", gen = T, hiv = T, df_to_add)
list_raw_prob_3y[[1]]$hivstat = NA
raw_prob_3y <- rbind(list_raw_prob_3y[[1]], list_raw_prob_3y[[2]])

raw_prob_ever$time_length <- "ever"
raw_prob_3y$time_length <- "3y"
raw_prob <- rbind(raw_prob_3y, raw_prob_ever)

#change regions according to gbd classifications
raw_prob <- select(raw_prob, -Region)
region_country <- surveys[,c("Country", "Region")]
raw_prob <- left_join(raw_prob, region_country)
raw_prob$Region[raw_prob$Region == "Central" | raw_prob$Region == "Western"] <- "Central/Western"

df_comp <- left_join(pred_prob, raw_prob)
df_comp <- df_comp[!(df_comp$'HIV Biomarker' == "Y" & is.na(df_comp$hivstat)),]
df_comp$hivstat[is.na(df_comp$hivstat)] <- "Either"
#df_comp <- filter(df_comp, !is.na(df_comp$hivstat))
df_comp$'Modelled Estimates' <- paste(df_comp$'50%', df_comp$'2.5%', df_comp$'97.5%', sep = ", ")
df_comp$Data <- paste(df_comp$raw_mean, df_comp$raw_lower, df_comp$raw_upper, sep = ", ")

df_comp_long <- pivot_longer(df_comp, cols = c("Modelled Estimates", "Data"), values_to = "val", names_to = "type")
df_comp_long <- df_comp_long %>%
  select(-c("50%", "2.5%", "97.5%", "raw_mean", "raw_lower", "raw_upper")) %>%
  separate(col = "val", into = c("median", "lower", "upper"), sep = ", ")

df_comp_long[,c("median", "lower", "upper")] <- lapply(df_comp_long[,c("median", "lower", "upper")], as.numeric)

#in-samp comp
df_comp$Error <- df_comp$raw_mean - df_comp$'50%'
df_comp$Abs_Error <- abs(df_comp$Error)
df_comp$Rltv_Error <- abs(df_comp$Error)/df_comp$raw_mean
df_comp$Below <- 0
df_comp$Above <- 0
df_comp$Below[df_comp$raw_mean < df_comp$'2.5%'] <- 1
df_comp$Above[df_comp$raw_mean > df_comp$'97.5%'] <- 1

error_gr_combined73 <- df_comp %>%
  group_by(time_length) %>%
  summarise(Error = round(median(Error)*100, 2), 
            Abs_Error = round(median(Abs_Error)*100, 2), 
            Rltv_Error = round(median(Rltv_Error, na.rm = T)*100, 2),
            Below = round(mean(Below)*100, 2), 
            Above = round(mean(Above)*100, 2))

total <- data.frame(time_length = "Any",
                    Error = round(median(df_comp$Error)*100, 2),
                    Abs_Error = round(median(abs(df_comp$Error))*100, 2),
                    Rltv_Error = round(median((abs(df_comp$Error)/df_comp$raw_mean)*100, na.rm = T), 2),
                    Below = round(mean(df_comp$Below)*100, 2),
                    Above = round(mean(df_comp$Above)*100, 2))

error_gr_combined73 <- rbind(error_gr_combined73, total)

error_gr_combined73$time_length <- factor(error_gr_combined73$time_length, levels = c("3y", "ever", "Any"), labels = c("Within the past 3Y", "Ever", "Any"))
names(error_gr_combined73) <- c("Screening Timing", "Error (%)", "Absolute Error (%)", "Absolute Relative Error (%)", "Below 95% CrI (%)", "Above 95% CrI (%)")

kbl(error_gr_combined73, caption = "In-Sample Comparisons: Model restricted to ages 25-49") %>%
  kable_styling() %>%
  footnote(general = "Random effects: Country, Survey, TimexRegion, HIV StatusxCountry\nFixed effects: Age")

#plots
#splitting the df into 2 because there's so many terms
# first_half <- sort(unique(df_comp_long$ISO3))[1:13]
# second_half <- sort(unique(df_comp_long$ISO3))[14:25]
# 
# df_comp_long1 <- df_comp_long[with(df_comp_long, ISO3 %in% first_half),]
# df_comp_long2 <- df_comp_long[with(df_comp_long, ISO3 %in% second_half),]
# 
#p_isc_gr_combined73_sa <- 
  ggplot(df_comp_long[df_comp_long$Region == "Southern",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.4)) +
  ylim(0,1) +
  ylab("Proportion") +
  #ggtitle("Posterior Predictive Check: Southern Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
  theme_bw() +
  guides(x=guide_axis_nested(angle = 90)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        # axis.text.y = element_text(size = 30),
        # axis.title.y = element_text(size = 30),
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, vjust = 1, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        # legend.text = element_text(size = 20),
        # plot.title = element_text(size = 30),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(angle = 90)) +
  scale_colour_nejm(labels = c("Survey Data", "Modelled Estimates"))

#p_isc_gr_combined73_ea <- 
  ggplot(df_comp_long[df_comp_long$Region == "Eastern",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.4)) +
  ylim(0,1) +
  ylab("Proportion") +
  #ggtitle("Posterior Predictive Check: Eastern Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
  theme_minimal() +
  guides(x=guide_axis_nested(angle = 90)) +
  theme(axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        # axis.text.y = element_text(size = 30),
        # axis.title.y = element_text(size = 30),
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, vjust = 1, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        # legend.text = element_text(size = 20),
        # plot.title = element_text(size=30),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(angle = 90)) +
    scale_colour_nejm(labels = c("Survey Data", "Modelled Estimates"))

#p_isc_gr_combined73_cwa <- 
  ggplot(df_comp_long[df_comp_long$Region == "Central/Western",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.4)) +
  ylim(0,1) +
  ylab("Proportion") +
  #ggtitle("Posterior Predictive Check: Central/Western Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
  theme_minimal() +
  guides(x=guide_axis_nested(angle = 90)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        # axis.text.y = element_text(size = 30),
        # axis.title.y = element_text(size = 30),
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, vjust = 1, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        # legend.text = element_text(size = 20),
        # plot.title = element_text(size=30),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(angle = 90)) +
    scale_colour_nejm(labels = c("Survey Data", "Modelled Estimates"))

p_isc_gr_combined73_sa; p_isc_gr_combined73_ea; p_isc_gr_combined73_cwa

save(df_comp, df_comp_long, 
     waic_gr_combined73, loo_gr_combined73,
     error_gr_combined73,
     p_isc_gr_combined73_sa,
     p_isc_gr_combined73_ea,
     p_isc_gr_combined73_cwa,
     #p_isc_gr_combined73_1, p_isc_gr_combined73_2,
     file = "mod_val_gr_combined73")


ggplot(df_comp_long[df_comp_long$Region == "Southern",], aes(x = Year, y = median, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.4)) +
  ylim(0,1) +
  ylab("Proportion") +
  #ggtitle("Posterior Predictive Check: Southern Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(~ hivstat + time_length + agegr + ISO3 , switch = "x", scales = "free_x", space = "free_x") +
  theme_bw() +
  guides(x=guide_axis_nested(angle = 90)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        # axis.text.y = element_text(size = 30),
        # axis.title.y = element_text(size = 30),
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, vjust = 1, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        # legend.text = element_text(size = 20),
        # plot.title = element_text(size = 30),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(angle = 90)) +
  scale_colour_nejm(labels = c("Survey Data", "Modelled Estimates"))

# ----  i) Extra Code ----
#add screening program to df_pred
# any_scr_pgm <- filter(screening_program, pgm_type == "Any")
# any_scr_pgm$screen_end[any_scr_pgm$screen_end == "Present"] <- 3000
# any_scr_pgm$screen_end <- as.numeric(any_scr_pgm$screen_end)
# df_pred$pgm <- 0
# time_lag_start <- 0
# time_lag_end <- 0
# for(c in 1:length(unique(df_pred$Country))){
#   country <- unique(df_pred$Country)[c]
#   df_pgm <- filter(any_scr_pgm, Country == country)
#   if(nrow(df_pgm) > 0){
#     for(j in 1:nrow(df_pgm)){
#       df_pred$pgm[df_pred$Country == country & (df_pred$Year + 2000) >= (df_pgm$screen_start[j] + time_lag_start) & (df_pred$Year + 2000) <= (df_pgm$screen_end[j] + time_lag_end)] <- 1
#     }
#   }
# }

#add gdp to df_stan
# gni_raw <- read.csv("gni.csv")
# gni <- pivot_longer(gni_raw, cols = -c(1:4), names_to = "Year", values_to = "gni")
# gni$Year <- as.numeric(sub(".", "", gni$Year))
# names(gni)[names(gni) == "Country.Name"] <- "Country"
# gni <- gni[,-c(2:4)]
# gni$Country[gni$Country == "Congo, Rep."] <- "Congo (Rep)"
# df_pred$Year <- df_pred$Year + 2000
# df_pred <- left_join(df_pred, gni)
# #create income group variable
# gni_thresholds <- read_excel("gni_thresholds.xlsx")
# gni_thresholds <- pivot_wider(gni_thresholds, names_from = Class, values_from = Max_Threshold)
# df_pred <- left_join(df_pred, gni_thresholds)
# df_pred$income_group[df_pred$gni < df_pred$Low_LowerMiddle] <- "Low"
# df_pred$income_group[df_pred$gni >= df_pred$Low_LowerMiddle & df_pred$gni <= df_pred$LowerMiddle_UpperMiddle] <- "LowerMiddle"
# df_pred$income_group[df_pred$gni > df_pred$LowerMiddle_UpperMiddle & df_pred$gni <= df_pred$UpperMiddle_High] <- "UpperMiddle"
# df_pred$income_group[df_pred$gni > df_pred$UpperMiddle_High] <- "High"
# df_pred$income_group_num <- as.numeric(factor(df_pred$income_group, levels = c("Low", "LowerMiddle", "UpperMiddle", "High")))
# df_pred$Year <- df_pred$Year - 2000
# 
# #create income group matrix
# G <- length(unique(df_pred$income_group_num))
# for(g in 2:G){
#   df_pred$x[df_pred$income_group_num == g] <- 1
#   df_pred$x[df_pred$income_group_num != g] <- 0
#   name <- paste0("income_group", g)
#   names(df_pred)[names(df_pred) == "x"] <- name
# }
# 
# G <- G-1 #subtract one for the model as there will be n-1 parameters
# 
# gni_pred <- as.matrix(df_pred[,c("income_group2", "income_group3")])

#15+
# pop_age <- pop_age[,c("Location",
#                       "15-19", "20-24", "25-29", "30-34",
#                       "35-39", "40-44", "45-49", "50-54",
#                       "55-59", "60-64", "65-1000")]
# pop_age$`65-1000`<- NULL
# for(i in 1:nrow(pop_age)){
#   pop_age[i,"65-1000"] <- sum(pop_age[i,"65-69"],
#                               pop_age[i,"70-74"],
#                               pop_age[i,"75-79"],
#                               pop_age[i,"80-84"],
#                               pop_age[i,"85-89"],
#                               pop_age[i,"90-94"],
#                               pop_age[i,"95-99"],
#                               pop_age[i,"100+"])
# }
#50-65+
#pop_age <- pop_age[,c("Location", "50-54", "55-59", "60-64", "65-1000")]

#change to french labels
# count_fr <- c("Bénin", "Botswana", "Burkina Faso", "Cameroun", "TChad", "Comores", "Congo (Rép)", "Cote d'Ivoire",
#               "Eswatini","Ethiopie", "Ghana", "Kenya", "Lesotho", "Malawi", "Mali", "Mauritanie", "Maurice", "Namibie",
#               "Rwanda", "Senegal", "Afrique du Sud", "Tanzanie", "Ouganda", "Zambie", "Zimbabwe")
# prob_by_count_year$Country_fr <- factor(prob_by_count_year$Country, labels = count_fr)
# prob_by_reg_year$Region_fr <- factor(prob_by_reg_year$Region, labels = c("centrale/de l'Ouest", "de l'Est", "du sud"))
# prob_by_hiv_count_year_3y$Country_fr <- factor(prob_by_count_year$Country, labels = count_fr)
# prob_by_count_year_ever$Country_fr <- factor(prob_by_count_year$Country, labels = count_fr)
# 
# raw_prob$Country_fr <- factor(raw_prob$Country, labels = c("Bénin",
#                                                                                "Botswana",
#                                                                                "Burkina Faso",
#                                                                                "Cameroun",
#                                                                                "TChad",
#                                                                                "Comores",
#                                                                                "Congo (Rép)",
#                                                                                "Cote d'Ivoire",
#                                                                                "Eswatini",
#                                                                                "Ethiopie",
#                                                                                "Ghana",
#                                                                                "Kenya",
#                                                                                "Lesotho",
#                                                                                "Malawi",
#                                                                                "Mali",
#                                                                                "Mauritanie",
#                                                                                "Maurice",
#                                                                                "Namibie",
#                                                                                "Rwanda",
#                                                                                "Senegal",
#                                                                                "Afrique du Sud",
#                                                                                "Tanzanie",
#                                                                                "Ouganda",
#                                                                                "Zambie",
#                                                                                "Zimbabwe"))
# raw_prob$Region_fr <- factor(raw_prob$Region, labels = c("centrale/de l'Ouest", "de l'Est", "du sud"))

# #p_count_fr <- 
#   ggplot(prob_by_count_year, aes(x = Year, y = mean)) +
#     geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
#     geom_line(aes(colour = time_length)) + 
#     geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
#     ylim(0, 1) +
#     ylab("Proportion") +
#     xlab("Année") +
#     ggtitle("Time trends for CC screening for for those aged 15-49 by country") +
#     labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
#          caption = "*Points represent weighted estimates from survey data") +
#     theme_minimal() +
#     theme(axis.text=element_text(size=15),
#           axis.title=element_text(size=17),
#           axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
#           strip.text = element_text(size = 12)) +
#     scale_colour_brewer(name = "Calendrier de dépistage", palette = "Set2", labels = c("3a", "à vie")) +
#     scale_fill_brewer(name = "Calendrier de dépistage", palette = "Set2", labels = c("3a", "à vie")) +
#     facet_wrap(~Country_fr)
#   
# #p_reg_fr <- 
#   ggplot(prob_by_reg_year, aes(x = Year, y = mean)) +
#     geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
#     geom_line(aes(colour = time_length)) + 
#     geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = time_length), alpha = 0.25) +
#     geom_hline(yintercept = 0.7, lty = "dashed", colour = "darkred") +
#     ylim(0, 1) +
#     ylab("Proportion") +
#     xlab("Année") +
#     ggtitle("Time trends for CC screening for for those aged 15-49 by region") +
#     labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
#          caption = "*Points represent weighted estimates from survey data") +
#     theme_minimal() +
#     theme(axis.text=element_text(size=15),
#           axis.title=element_text(size=17),
#           axis.text.x = element_text(size = 15, angle = 45, hjust=0.95, vjust=1.25),
#           strip.text = element_text(size = 12),
#           legend.position = "bottom") +
#     scale_colour_brewer(name = "Calendrier de dépistage", palette = "Set2", labels = c("3a", "à vie")) +
#     scale_fill_brewer(name = "Calendrier de dépistage", palette = "Set2", labels = c("3a", "à vie")) +
#     facet_wrap(~Region_fr)


p4 <- ggplot(data = combine_q, aes(x = Year, y = Country)) +
  facet_grid(rows = vars(Country), scales = 'free', space = "free") +
  # add point for each survey
  geom_point(aes(colour = Survey, fill = Survey, shape = time_length), size = 2.5) +
  facet_grid(Region~Type, space = "free", scales = "free_y") +
  # add small white point for hiv biomarker data
  geom_point(data = combine_q[combine_q$hivbio == "Y", ], aes(fill = "white", colour = "white", shape = time_length), size = 1.5) +
  # text
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour = "black")) +
  #ggtitle("Survey Screening Data Availability") +
  labs(x = 'Year', y = 'Country') +
  theme_light() +
  theme(strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(colour = "black", size = 12),
        strip.text.y = element_text(colour = "black", size = 7.5)) +
  scale_x_continuous(breaks=seq(2000, 2021), labels=as.character(seq(2000, 2021))) +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2)) +
  #theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_text(colour = "black"), panel.grid.minor.x = element_blank()) +
  theme(axis.text.y = element_text(colour = "black")) +
  scale_colour_nejm() +
  scale_fill_nejm() +
  scale_shape_manual("Screening Interval", values = c(24, 25, 21)); p4

p2 <- ggplot(data = all_q[(all_q$Type == "CC Screening" | all_q$Type == "CC Treatment"),], aes(x = Year, y = Country)) +
  facet_grid(rows = vars(Country), scales = 'free', space = "free") +
  # add point for each survey
  geom_point(aes(colour = Survey, shape = hivbio), size = 2.5) +
  # text
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black")) +
  labs(x = 'Year', y = 'Country') +
  #ggtitle("Survey Data Availability") +
  facet_grid(Region~Type, space = "free", scales = "free") +
  theme_light() +
  theme(strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(colour = "black", size = 12),
        strip.text.y = element_text(colour = "black", size = 7.5)) +
  scale_x_continuous(breaks=seq(2000, 2020), labels=as.character(seq(2000, 2020)), limits = c(c(2000, 2020), c(2010, 2020))) +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2)) +
  #theme(strip.text.y = element_blank()) +
  theme(axis.text.x = element_text(colour = "black"), panel.grid.minor.x = element_blank()) +
  theme(axis.text.y = element_text(colour = "black")) +
  scale_colour_nejm() +
  scale_shape_manual("HIV Biomarker Data", values = c(17, 16)); p2














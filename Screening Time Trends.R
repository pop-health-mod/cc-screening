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
load("New Cleaned Pooled Surveys")

source("predict_poststrat_func.R")  

# ---- A) Prep ----
##Create datasets that remove NAs for each outcome variable (i.e. either pap smear ever, or cc test ever)
dat_pap <- filter(pooled_surveys, !is.na(ever_had_pap_smear) & !is.na(pap_smear_past_3y))
dat_cctest <- filter(pooled_surveys, !is.na(ever_cc_tested) & !is.na(cc_test_past_3y))

##Create new dataset which has a variable for either pap smear or a general cc test (the surveys asked questions regarding pap smears only to many of the southern african countries - this approach will assume that pap smears were the only/primary form of cc testing at time of survey)
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

design <- 2.0 #median value from previously created design effect

#we create a dataset where women with the same covariate patterns are grouped together
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

#to speed computations, we only select covariate patterns that are found in the data
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

#create grouped dataset
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
  design <- 2.0
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
  #here i will apply the design effect to the denominators 
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
  
  #add values from Zimbabwe PHIA 2020 (note only values for women HIV+ are available)
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
  f_hivprev_1549$Country[f_hivprev_1549$Country == "Cabo Verde"] <- "Cape Verde"
  f_hivprev_1549$Country[f_hivprev_1549$Country == "Congo"] <- "Congo (Rep)"
  f_hivprev_1549$Country[f_hivprev_1549$Country == "Gambia"] <- "The Gambia"
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
  f_hivprev_1524$Country[f_hivprev_1524$Country == "Cabo Verde"] <- "Cape Verde"
  f_hivprev_1524$Country[f_hivprev_1524$Country == "Congo"] <- "Congo (Rep)"
  f_hivprev_1524$Country[f_hivprev_1524$Country == "Gambia"] <- "The Gambia"
  f_hivprev_1524$hivprev1524[f_hivprev_1524$hivprev1524 == "<0.1 "] <- 0.1
  f_hivprev_1524$Year <- as.numeric(f_hivprev_1524$Year)
  f_hivprev_1524$hivprev1524 <- as.numeric(f_hivprev_1524$hivprev1524)
  f_hivprev_1524$hivprev1524 <- f_hivprev_1524$hivprev1524/100
  
  pop_age <- read_excel("pop_age0020.xlsx", sheet = 2) #i should change this so i can use the age structure for every year
  colnames(pop_age) <- pop_age[1,]
  pop_age <- pop_age[-c(1:3), -c(1, 3, 5)]
  pop_age[,-c(1:2)] <- lapply(pop_age[,-c(1:2)], as.numeric)
  pop_age_long <- pivot_longer(pop_age, cols = -(1:2), names_to = "Year", values_to = "f_pop")
  pop_age_long$Year <- as.numeric(pop_age_long$Year)
  names(pop_age_long)[names(pop_age_long) == "Location"] <- "Country"
  names(pop_age_long)[names(pop_age_long) == "Age"] <- "agegr"
  #clean up data to match 
  pop_age_long$Country[pop_age_long$Country == "Cabo Verde"] <- "Cape Verde"
  pop_age_long$Country[pop_age_long$Country == "Congo"] <- "Congo (Rep)"
  pop_age_long$Country[pop_age_long$Country == "Côte d'Ivoire"] <- "Cote d'Ivoire"
  pop_age_long$Country[pop_age_long$Country == "United Republic of Tanzania"] <- "Tanzania"
  pop_age_long$Country[pop_age_long$Country == "Gambia"] <- "The Gambia"
  
  agegr_1549 <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  agegr_1524 <- c("15-19", "20-24")
  agegr_2549 <- c("25-29", "30-34", "35-39", "40-44", "45-49")
  
  pop_age_1549 <- unique(pop_age_long %>%
                           group_by(Country, Year) %>%
                           filter(agegr %in% agegr_1549) %>%
                           summarise(Country, Year, f_pop_1549 = sum(f_pop)))
  pop_age_1524 <- unique(pop_age_long %>%
                           group_by(Country, Year) %>%
                           filter(agegr %in% agegr_1524) %>%
                           summarise(Country, Year, f_pop_1524 = sum(f_pop)))
  pop_age_2549 <- unique(pop_age_long %>%
                           group_by(Country, Year) %>%
                           filter(agegr %in% agegr_2549) %>%
                           summarise(Country, Year, f_pop_2549 = sum(f_pop)))
  
  pop_age <- left_join(pop_age_1549, pop_age_1524)
  pop_age <- left_join(pop_age, pop_age_2549)
  f_hivprev <- left_join(f_hivprev_1549, f_hivprev_1524)
  f_hivprev <- left_join(pop_age, f_hivprev)
  
  #obtain values for 25-49 using 15-24 and 25-49
  f_hivprev$hivprev2549 <- (f_hivprev$hivprev1549*f_hivprev$f_pop_1549 - f_hivprev$hivprev1524*f_hivprev$f_pop_1524)/f_hivprev$f_pop_2549

  #add hiv prevalence to surveys that didn't have hiv biomarker
  dat_timetrends_gr <- left_join(dat_timetrends_gr, det_avail[,c("df", "HIV Biomarker")])

  df_imp <- dat_timetrends_gr[dat_timetrends_gr$'HIV Biomarker' == "N" | is.na(dat_timetrends_gr$hivstat),] #includes indv from surveys with biomarkers but don't have data on their biomarker
  df_imp <- (left_join(df_imp, f_hivprev[,c("Country", "Year", "hivprev1549", "hivprev2549")]))

  df_noimp <- dat_timetrends_gr[dat_timetrends_gr$'HIV Biomarker' != "N" & !is.na(dat_timetrends_gr$hivstat),]
  df_noimp$hivprev1549 <- df_noimp$hivstat
  df_noimp$hivprev2549 <- df_noimp$hivstat

  df_cmb <- rbind(df_imp, df_noimp)
  
  saveRDS(df_cmb, file = "new_hiv_dat_timetrends_gr.rds")
}

# read in the data
hiv_df_gr <- readRDS("new_hiv_dat_timetrends_gr.rds")

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
  vector[R] re_region; 
  real<lower = 0, upper = 10> sd_re_region;
  
  vector[R] rs_time;
  real<lower = 0, upper = 10> sd_rs_region;
  real slope_time; 
  
  //vector[C] rs_pgm; 
  //real<lower = 0, upper = 10> sd_rs_country_pgm;
  //real slope_pgm;
  
  row_vector[A] beta_age;
  //row_vector[G] beta_gni;
  
  real<lower = -5, upper = 5> beta_hiv;
  vector[C] rs_hiv_country;
  vector[R] rs_hiv_region;
  real<lower = 0, upper = 10> sd_rs_country_hiv;
  real<lower = 0, upper = 10> sd_rs_region_hiv;
  
  //real beta_pgm;
  
  row_vector[A + 1] beta_time_length;
}

transformed parameters{
  vector[N] logit_prd;
  vector[N] logit_prob;  
  vector[n_hiv] hiv_effect;
  matrix[N, n_hiv] prd;
  
  hiv_effect[1] = 0;
  hiv_effect[2] = 1;
  for(n in 1:N) {
    logit_prd[n] = alpha + re_region[region[n]] + re_country[country[n]] + re_survey[survey[n]] 
    + rs_time[region[n]] * time[n] 
    //+ rs_pgm[country[n]] * pgm[n]
    //+ beta_time * time[n] 
    + beta_age * age[, n]
    //+ beta_gni * gni[, n]
    //+ beta_pgm * pgm[n]
    + beta_time_length * time_length[, n];
    for (h in 1:n_hiv) {
      prd[n, h] = inv_logit(logit_prd[n] + (beta_hiv + rs_hiv_region[region[n]] + rs_hiv_country[country[n]]) * hiv_effect[h]);
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
  //beta_gni ~ normal(0, 5);
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
#Compile the model
cc_trends_gr_stan <- rstan::stan_model(model_code = cc_trends_gr)

# ----  C) Fit ----
#change regions from un classification to gbd
df_stan <- select(hiv_df_gr, -Region)
GBD = T
if(GBD == T){
  region_country_2 <- unique(surveys[,c("Country", "Region")]) #either un or gbd classification was chosen straight in the excel document
  df_stan <- left_join(df_stan, region_country_2)
}

#restrict analysis by age 25-49
df_stan <- filter(df_stan, agegr == "25-29" | agegr == "30-34" | agegr == "35-39" | agegr == "40-44" | agegr == "45-49")

#remove values with a denominator of 0
df_stan <- filter(df_stan, den != 0)

#prep variables for stan model
N <- nrow(df_stan)

S <- length(unique(df_stan$df))
survey <- as.numeric(factor(df_stan$df))

C <- length(unique(df_stan$Country))
df_stan$count_num <- as.numeric(factor(df_stan$Country))
country <- df_stan$count_num #levels will be in alphabetical order

df_stan$Region[df_stan$Region == "Central" | df_stan$Region == "Western"] <- "Central/Western"
R <- length(unique(df_stan$Region))
df_stan$reg_num <- as.numeric(factor(df_stan$Region))
region <- df_stan$reg_num 

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
A <- n_age - 1 #subtract one for the model as there will be n-1 parameters

#center year around 2000 (lowest survey year)
df_stan$year_cent <- df_stan$Year - 2000
time <- df_stan$year_cent

#create binary variable for time length, 1 = 3y, 0 = ever 
df_stan$time_length2[df_stan$time_length == "3y"] <- 1
df_stan$time_length2[df_stan$time_length == "ever"] <- 0

#creating time_length matrix if trying to have a different odds ratio for screening interval by age
n_age <- length(unique(df_stan$agegr))
for(a in 1:n_age){
  df_stan$x[df_stan$time_length2 == 1 & df_stan$agegr_num == a] <- 1
  df_stan$x[df_stan$time_length2 == 1 & df_stan$agegr_num != a] <- 0
  df_stan$x[df_stan$time_length2 == 0] <- 0
  name <- paste0("time_length_agegr", a)
  names(df_stan)[names(df_stan) == "x"] <- name
}

time_length = t(as.matrix(df_stan[,c("time_length_agegr1", "time_length_agegr2", "time_length_agegr3", "time_length_agegr4", "time_length_agegr5")]))

#group logistic reg specific variables
num <- df_stan$num
den <- df_stan$den

#variables for hiv prevalence
prv_wgt <- matrix(data = NA, nrow = N, ncol = 2)
prv_wgt[, 1] <- 1 - df_stan$hivprev2549 # hiv neg
prv_wgt[, 2] <- df_stan$hivprev2549     # hiv pos
n_hiv <- ncol(prv_wgt)

data_stan <- list(N = N,
                  S = S,
                  survey = survey,
                  C = C,
                  country = country,
                  R = R,
                  region = region,
                  A = A,
                  age = age,
                  time = time,
                  time_length = time_length,
                  n_hiv = n_hiv,
                  prv_wgt = t(prv_wgt),
                  num = num,
                  den = den)

# ----  D) Run ----
options(mc.cores = parallel::detectCores())
fit <- rstan::sampling(cc_trends_gr_stan, 
                       data = data_stan, 
                       iter = 3000, chains = 6, refresh = 50,
                       warmup = 1000, thin = 1,
                       control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15))
save(data_stan, fit, df_stan, cc_trends_gr_stan, file = "gr_combined77")
#load("gr_combined77")
rstan::rstan_options(auto_write = TRUE) 
   
#view convergence + model parameter summaries
rstan::stan_trace(fit)
rstan::summary(fit, probs = c(0.025, 0.5, 0.975))$summary

##obtain odds ratios for effects of HIV status
draws <- rstan::extract(fit)
rs_hiv_overall <- draws$beta_hiv
rs_hiv_region <- draws$rs_hiv_region
rs_hiv_country <- draws$rs_hiv_country

hiv_overall_or <- exp(quantile(rs_hiv_overall, probs = c(0.025, 0.5, 0.975)))
hiv_overall_or <- as.data.frame(t(hiv_overall_or))
hiv_overall_or$Country <- "Sub-Saharan Africa"
hiv_overall_or$Region <- "Overall"
hiv_overall_or$Type <- "Overall"

hiv_region_or <- matrix(data = NA, ncol = 3, nrow = 3)
for(i in 1:ncol(rs_hiv_region)){
  a <- rs_hiv_region[,i] + rs_hiv_overall
  hiv_region_or[i,] <- exp(quantile(a, probs = c(0.025, 0.5, 0.975)))
}
hiv_region_or <- as.data.frame(hiv_region_or)
hiv_region_or$Country <- c("Western/Central Africa", "Eastern Africa", "Southern Africa")
hiv_region_or$Region <- c("Western/Central", "Eastern", "Southern")
hiv_region_or$Type <- "Region"

df_rs_hiv <- data.frame(Country = unique(df_stan$Country), 
                        count_num = unique(df_stan$count_num),
                        t(rs_hiv_country[unique(df_stan$count_num),]))

df_rs_hiv <- left_join(df_rs_hiv, region_country)
df_rs_hiv$Region[df_rs_hiv$Region == "Central" | df_rs_hiv$Region == "Western"] <- "Western/Central"
df_rs_hiv$reg_num <- as.numeric(factor(df_rs_hiv$Region, levels = c("Western/Central", "Eastern", "Southern")))

hiv_country_or <- matrix(data = NA, ncol = 3, nrow = nrow(df_rs_hiv))
for(i in 1:nrow(df_rs_hiv)){
  a <- rs_hiv_overall + rs_hiv_region[,df_rs_hiv$reg_num[i]] + rs_hiv_country[,df_rs_hiv$count_num[i]]
  hiv_country_or[i,] <- exp(quantile(a, probs = c(0.025, 0.5, 0.975)))
}

hiv_country_or <- data.frame(hiv_country_or)
hiv_country_or$Country <- df_rs_hiv$Country
hiv_country_or$Region <- df_rs_hiv$Region
hiv_countries <- unique(df_stan[!is.na(df_stan$hivstat),]$Country)
hiv_country_or$Type <- "Country"

#create dataframe with national, regional, and overall OR
lab <- c("Lower", "Median", "Upper")
colnames(hiv_overall_or)[1:3] <- lab
colnames(hiv_region_or)[1:3] <- lab
colnames(hiv_country_or)[1:3] <- lab
hiv_country_or_1 <- hiv_country_or[hiv_country_or$Country %in% hiv_countries,]
hiv_or <- rbind(hiv_overall_or, hiv_region_or, hiv_country_or_1)
hiv_or$Country <- factor(hiv_or$Country, levels = rev(c("Sub-Saharan Africa", "Western/Central Africa", "Cameroon", "Cote d'Ivoire",
                                                    "Eastern Africa", "Ethiopia", "Kenya", "Malawi", "Rwanda", "Tanzania", "Zambia",
                                                    "Southern Africa", "Lesotho", "Namibia", "South Africa", "Zimbabwe")))
                         #labels = rev(c("SSA", "WCA", "CMR", "CIV", "EA", "ETH", "KEN", "MWI", "RWA", "TZA", "ZMB", "SA", "LSO", "NAM", "ZAF", "ZWE")))
hiv_or$Region <- factor(hiv_or$Region, levels = c("Overall", "Eastern", "Southern", "Western/Central"))
hiv_or$Type <- factor(hiv_or$Type, levels = c("Overall", "Region", "Country"))

#create forest plot for hiv odds ratios
ggplot(hiv_or[hiv_or$Type != "Overall",], aes(y = Country, x = Median, colour = Type)) +
  geom_point(shape = 18, size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, size = 1.5) +
  geom_vline(xintercept = 1, lty = "dotted") +
  ylab(" ") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
  scale_colour_nejm(name = "Estimate Type", labels = c("Regional", "National")) + 
  scale_x_continuous(name = "Odds Ratio (OR)", breaks = seq(0, 13, by = 1))

# ----  E) Predictions ----
df_stan$predicted <- "Y"
region_country$Region[region_country$Region == "Central" | region_country$Region == "Western"] <- "Central/Western"

#add all countries (even if no survey available) to data frame for imputations to be used for aggregate values
df_stan_all_count <- left_join(region_country, df_stan)
df_stan_all_count$predicted[is.na(df_stan_all_count$predicted)] <- "N"

source("predict_poststrat_func.R")

##find weights by country
df_weights_hiv <- weights(agegr = c("25-29", "30-34", "35-39", "40-44", "45-49"), hiv = T, years = c(2000:2020))

#Obtain unweighted predicted probabilities for each iteration and category. This can be done in generated quantities however I'm going to do it outside generated quantities 
#create data set to predict
draws <- rstan::extract(fit)

pred_prob_weighted_2549 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(0, 5, 10, 15, 20), agegr = 1:count(!is.na(unique(df_stan$agegr))), 
                                    hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)
pred_prob_weighted_3049 <- pred_weighted(draws, df_stan = df_stan_all_count, Country = unique(df_stan_all_count$Country), Year = c(0, 5, 10, 15, 20), agegr = 2:5, 
                                         hivstat_log = T, gni_log = F, pgm_log = F, multi_time_length = T, df_weights_hiv)

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

#combined ever and past 3y dataframes
prob_by_year <- dat_combine(prob_by_year_ever, prob_by_year_3y)
prob_by_reg_year <- dat_combine(prob_by_reg_year_ever, prob_by_reg_year_3y)
prob_by_count_year <- dat_combine(prob_by_count_year_ever, prob_by_count_year_3y)

prob_by_hiv_year <- dat_combine(prob_by_hiv_year_ever, prob_by_hiv_year_3y)
prob_by_hiv_reg_year <- dat_combine(prob_by_hiv_reg_year_ever, prob_by_hiv_reg_year_3y)
prob_by_hiv_count_year <- dat_combine(prob_by_hiv_count_year_ever, prob_by_hiv_count_year_3y)

prob_by_reg_year_age <- dat_combine(prob_by_reg_year_age_ever, prob_by_reg_year_age_3y)
prob_by_count_year_age <- dat_combine(prob_by_count_year_age_ever, prob_by_count_year_age_3y)

View(prob_by_year[,c("Year", "2.5%", "50%", "97.5%")])
hiv_countries <- c("Cameroon","Cote d'Ivoire", "Ethiopia", "Kenya", "Lesotho", "Malawi", "Namibia", "Rwanda", "South Africa", "Tanzania", "Zambia", "Zimbabwe")

hiv_count <- prob_by_hiv_count_year[prob_by_hiv_count_year$Country %in% hiv_countries,c("Country", "Year", "50%", "2.5%", "97.5%", "Region", "hivstat", "time_length")]
hiv_count$Type <- "Country"
hiv_region <- prob_by_hiv_reg_year[,c("Year","50%", "2.5%", "97.5%", "Region", "hivstat", "time_length")]
hiv_region$Country <- paste(hiv_region$Region, "Africa")
hiv_region$Type <- "Region"
hiv_overall <- prob_by_hiv_year[,c("Year","50%", "2.5%", "97.5%", "hivstat", "time_length")]
hiv_overall$Country <- "Sub-Saharan Africa"
hiv_overall$Region <- "Overall"
hiv_overall$Type <- "Overall"

hiv <- rbind(hiv_count, hiv_region, hiv_overall)
hiv$Country[hiv$Country == "Central/Western Africa"] <- "Western/Central Africa"
hiv$Region[hiv$Region == "Central/Western"] <- "Western/Central"
hiv$Country <- factor(hiv$Country, levels = c("Sub-Saharan Africa", "Western/Central Africa", "Cameroon", "Cote d'Ivoire",
                                                        "Eastern Africa", "Ethiopia", "Kenya", "Malawi", "Rwanda", "Tanzania", "Zambia",
                                                        "Southern Africa", "Lesotho", "Namibia", "South Africa", "Zimbabwe"))
                      #labels = c("SSA", "WCA", "CMR", "CIV", "EA", "ETH", "KEN", "MWI", "RWA", "TZA", "ZMB", "SA", "LSO", "NAM", "ZAF", "ZWE"))
hiv$Region <- factor(hiv$Region, levels = c("Overall", "Western/Central", "Eastern", "Southern"))
hiv$Type <- factor(hiv$Type, levels = c("Overall", "Region", "Country"))
hiv$time_length[hiv$time_length == "3y"] <- "Past 3Y"
hiv$time_length[hiv$time_length == "ever"] <- "Lifetime"

hiv_2020 <- filter(hiv, Year == 2020, time_length == "Lifetime")

#plot
ggplot(hiv_2020, aes(x = Country, y = `50%`*100, fill = hivstat)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=`2.5%`*100, ymax=`97.5%`*100), width=.5, position = position_dodge(0.9)) +
  geom_hline(yintercept = 70, lty = "dotted") +
  facet_grid(time_length~Type, space = "free", scales = "free") +
  ylim(0, 100) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(name = "HIV Status", values = c("#7876B1FF", "#E18727FF")) +
  ylab("Percentage(%)")

#create summary dataframes
year = 2020
agegr = "45-49"
hivstat = "HIV+"

a <- prob_by_year_age_ever[prob_by_year_age_ever$Year == year & prob_by_year_age_ever$agegr == agegr,]
b <- prob_by_count_year_age_ever[prob_by_count_year_age_ever$Year == year & prob_by_count_year_age_ever$Region == "Southern" & prob_by_count_year_age_ever$agegr == agegr,]
c <- prob_by_count_year_age_ever[prob_by_count_year_age_ever$Year == year & prob_by_count_year_age_ever$Region == "Central/Western" & prob_by_count_year_age_ever$agegr == agegr ,]
d <- prob_by_count_year_age_ever[prob_by_count_year_age_ever$Year == year & prob_by_count_year_age_ever$Region == "Eastern" & prob_by_count_year_age_ever$agegr == agegr,]
e <- prob_by_reg_year_age_ever[prob_by_reg_year_age_ever$Year == year & prob_by_reg_year_age_ever$agegr == agegr,]

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

##raw probabilities from survey CC Screening in the Past 3Y and Lifetime Screening survey estimates
raw_prob_ever <- survey_prob(timing = "ever", gen = T, hiv = T, df_to_add)
raw_prob_3y <- survey_prob(timing = "3y", gen = T, hiv = T, df_to_add)

raw_prob <- dat_combine(raw_prob_ever[[1]], raw_prob_3y[[1]])
hiv_raw_prob <- dat_combine(raw_prob_ever[[2]], raw_prob_3y[[2]])

agegr <- c("30-34", "35-39", "40-44", "45-49")
raw_prob_age <- raw_prob[raw_prob$agegr %in% agegr, ]
  
raw_prob <- raw_prob[raw_prob$agegr == "30-49",]
hiv_raw_prob <- hiv_raw_prob[hiv_raw_prob$agegr == "25-49",]
# ----  F) Plots ----
p_overall <- 
  ggplot(prob_by_year, aes(x = Year, y = mean*100)) +
  #geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`*100, ymax = `97.5%`*100, fill = time_length), alpha = 0.25) +
  geom_hline(yintercept = 70, lty = "dashed", colour = "darkred") +
  ylim(0, 100) +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20)) +
  scale_color_lancet(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  scale_fill_lancet(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime"));p_overall

p_hiv_overall <- 
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
prob_by_count_year$Region[prob_by_count_year$Region == "Central/Western"] <- "Western/Central"
raw_prob$Region[raw_prob$Region == "Central/Western"] <- "Western/Central"
p_count_2p <- 
  ggplot(prob_by_count_year[prob_by_count_year$Country %in% countries_2p & prob_by_count_year$Region == "Western/Central",], aes(x = Year, y = mean*100)) +
  geom_pointrange(data = raw_prob[raw_prob$Country %in% countries_2p & raw_prob$Region == "Western/Central",], aes(y = raw_mean*100, ymin = raw_lower*100, ymax = raw_upper*100, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_hline(yintercept = 70, lty = "dashed", colour = "darkred") +
  geom_ribbon(aes(ymin = `2.5%`*100, ymax = `97.5%`*100, fill = time_length), alpha = 0.25) +
  ylim(0, 100) +
  ylab("Percentage (%)") +
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
  scale_color_lancet(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  scale_fill_lancet(name = "Recall Period", labels = c("Within the Past 3Y", "Lifetime")) +
  facet_grid(Region~Country); p_count_2p

prob_by_reg_year$Region[prob_by_reg_year$Region == "Central/Western"] <- "Western/Central"
p_reg <- 
  ggplot(prob_by_reg_year, aes(x = Year, y = mean*100)) +
  #geom_pointrange(data = raw_prob, aes(y = raw_mean, ymin = raw_lower, ymax = raw_upper, colour = time_length), size = 0.2) +
  geom_line(aes(colour = time_length)) + 
  geom_ribbon(aes(ymin = `2.5%`*100, ymax = `97.5%`*100, fill = time_length), alpha = 0.25) +
  geom_hline(yintercept = 70, lty = "dashed", colour = "darkred") +
  ylim(0, 100) +
  ylab("Percentage (%)") +
  # ggtitle("Time trends for CC screening for for those aged 30-49 by region") +
  # labs(subtitle = "Random effects: intercept by country, survey; slopes for time by region, HIV status by Country\nFixed effects: Age",
  #      caption = "*Points represent weighted estimates from survey data") +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5),
        strip.text = element_text(size = 20),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20),
        legend.title=element_text(size=20)) +
  scale_colour_lancet(name = "Recall Period", labels = c("Within the past 3Y", "Lifetime")) +
  scale_fill_lancet(name = "Recall Period", labels = c("Within the past 3Y", "Lifetime")) +
  facet_wrap(~Region); p_reg
  
hiv_countries_ever <- c("Cameroon","Cote d'Ivoire", "Ethiopia", "Kenya", "Lesotho", "Malawi", "Namibia", "Rwanda", "South Africa", "Tanzania", "Zambia", "Zimbabwe")
hiv_countries_3y <- c("Cameroon", "Ethiopia", "Lesotho", "Malawi", "Rwanda", "South Africa", "Tanzania", "Zambia", "Zimbabwe")

hiv_count_year_ever <- prob_by_hiv_count_year_ever[prob_by_hiv_count_year_ever$Country %in% hiv_countries_ever,]
hiv_count_year_3y <- prob_by_hiv_count_year_3y[prob_by_hiv_count_year_3y$Country %in% hiv_countries_3y,]
hiv_count_year_ever$Region[hiv_count_year_ever$Region == "Central/Western"] <- "Western/Central"
hiv_raw_prob$Region[hiv_raw_prob$Region == "Central/Western"] <- "Western/Central"

p_hiv_count_ever <-
  ggplot(hiv_count_year_ever[hiv_count_year_ever$Region == "Eastern",], aes(x = Year, y = mean*100)) +
    geom_line(aes(colour = hivstat)) +
    geom_ribbon(aes(ymin = `2.5%`*100, ymax = `97.5%`*100, fill = hivstat), alpha = 0.25) +
    geom_hline(yintercept = 70, lty = "dashed", colour = "darkred") +
    geom_pointrange(data = hiv_raw_prob[hiv_raw_prob$time_length == "ever" & hiv_raw_prob$Region == "Eastern",], aes(y = raw_mean*100, ymin = raw_lower*100, ymax = raw_upper*100, colour = hivstat), size = 0.2) +
    ylim(0, 100) +
    ylab("Percentage (%)") +
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
    facet_grid(Region~Country);p_hiv_count_ever
  
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

# ----  G) Save ----
save(prob_by_count_year_3y, prob_by_count_year_ever, prob_by_count_year,
     prob_by_reg_year_3y, prob_by_reg_year_ever,
     prob_by_hiv_count_year_3y, prob_by_hiv_count_year_ever,
     p_overall, p_hiv_overall,
     p_count, p_count_2p, p_reg,
     p_hiv_count_3y, p_hiv_count_ever,
     pred_prob_weighted_2549,
     pred_prob_weighted_3049,
     file = "pred_prob_gr_combined77_all_ssa")

# load("pred_prob_gr_combined77_all_ssa")

# ----  H) Model validation ----
## Model validation 
#load("gr_combined77")
#WAIC
ll <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = T)
waic_gr_combined77 <- waic(ll)
#loo
ll_2 <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = F)
r_eff <- relative_eff(exp(ll_2), cores = 6) 
loo_gr_combined77 <- loo(ll_2, r_eff = r_eff, cores = 6)
print(waic_gr_combined77); print(loo_gr_combined77)

##In sample comparisons
#obtain values from model
draws <- rstan::extract(fit)
pred_prob <- data.frame(t(draws$pred_prob))
pred_prob <- row_sumstats(pred_prob, prob = c(0.5, 0.025, 0.975))
pred_prob <- cbind(df_stan, pred_prob[,c("50%", "2.5%", "97.5%")])
pred_prob$hivstat <- factor(pred_prob$hivstat, levels = c(0, 1), labels = c("HIV-", "HIV+"))

#values from survey
source("predict_poststrat_func.R")
list_raw_prob_ever <- survey_prob(timing = "ever", gen = T, hiv = T, df_to_add) #df_to_add can be obtained using the code from lines 151 to 207 
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

error_gr_combined77 <- df_comp %>%
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

error_gr_combined77 <- rbind(error_gr_combined77, total)

error_gr_combined77$time_length <- factor(error_gr_combined77$time_length, levels = c("3y", "ever", "Any"), labels = c("Within the past 3Y", "Ever", "Any"))
names(error_gr_combined77) <- c("Screening Timing", "Error (%)", "Absolute Error (%)", "Absolute Relative Error (%)", "Below 95% CrI (%)", "Above 95% CrI (%)")

kbl(error_gr_combined77, caption = "In-Sample Comparisons: Model restricted to ages 25-49") %>%
  kable_styling() %>%
  footnote(general = "Random effects: Country, Survey, TimexRegion, HIV StatusxCountry\nFixed effects: Age")

#plots
df_comp_long$Year <- factor(df_comp_long$Year, levels = unique(df_comp_long$Year), labels = paste0(" ", unique(df_comp_long$Year), " "))
df_comp_long$hivstat <- factor(df_comp_long$hivstat, levels = c("Either", "HIV-", "HIV+"), labels = c(" Either ", " HIV- ", " HIV+ "))
df_comp_long$time_length <- factor(df_comp_long$time_length, levels = c("ever", "3y"), labels = c(" Lifetime ", " Past 3Y "))
df_comp_long$agegr <- factor(df_comp_long$agegr, levels = unique(df_comp_long$agegr), labels = paste0(" ", unique(df_comp_long$agegr), " "))
p_isc_gr_combined77_sa <- 
  ggplot(df_comp_long[df_comp_long$Region == "Southern",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median*100, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100), position = position_dodge(0.4)) +
  ylim(0,100) +
  ylab("Percentage (%)") +
  #ggtitle("Posterior Predictive Check: Southern Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
  theme_bw() +
  guides(x=guide_axis_nested(angle = 90)) +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 30, angle = 0, vjust = 1, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        # legend.text = element_text(size = 20),
        # plot.title = element_text(size = 30),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(angle = 90)) +
  scale_colour_lancet(labels = c("Survey Data", "Modelled Estimates"))

p_isc_gr_combined77_ea <- 
  ggplot(df_comp_long[df_comp_long$Region == "Eastern",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median*100, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100), position = position_dodge(0.4)) +
  ylim(0,100) +
  ylab("Percentage (%)") +
  #ggtitle("Posterior Predictive Check: Eastern Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
    theme_bw() +
    guides(x=guide_axis_nested(angle = 90)) +
    theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          panel.border = element_blank(), 
          panel.spacing.x = unit(0,"line"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = 30, angle = 0, vjust = 1, hjust = 0),
          legend.position = "bottom",
          legend.title = element_blank(),
          # legend.text = element_text(size = 20),
          # plot.title = element_text(size = 30),
          axis.ticks = element_line(colour = "black"),
          ggh4x.axis.nesttext.x = element_text(angle = 90)) +
    scale_colour_lancet(labels = c("Survey Data", "Modelled Estimates"))

p_isc_gr_combined77_cwa <- 
  ggplot(df_comp_long[df_comp_long$Region == "Central/Western",], aes(x = interaction(Year, hivstat, time_length, agegr, lex.order = T), y = median*100, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100), position = position_dodge(0.4)) +
  ylim(0,100) +
  ylab("Percentage (%)") +
  #ggtitle("Posterior Predictive Check: Central/Western Africa\nRandom effects: Country, Survey, TimexRegion, HIV hivstatxCountry\nFixed effects: Age") + 
  facet_grid(.~ISO3, switch = "x", scales = "free_x", space = "free_x") +
    theme_bw() +
    guides(x=guide_axis_nested(angle = 90)) +
    theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          panel.border = element_blank(), 
          panel.spacing.x = unit(0,"line"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = 30, angle = 0, vjust = 1, hjust = 0),
          legend.position = "bottom",
          legend.title = element_blank(),
          # legend.text = element_text(size = 20),
          # plot.title = element_text(size = 30),
          axis.ticks = element_line(colour = "black"),
          ggh4x.axis.nesttext.x = element_text(angle = 90)) +
    scale_colour_lancet(labels = c("Survey Data", "Modelled Estimates"))

p_isc_gr_combined77_sa; p_isc_gr_combined77_ea; p_isc_gr_combined77_cwa

save(df_comp, df_comp_long, 
     waic_gr_combined77, loo_gr_combined77,
     error_gr_combined77,
     p_isc_gr_combined77_sa,
     p_isc_gr_combined77_ea,
     p_isc_gr_combined77_cwa,
     #p_isc_gr_combined77_1, p_isc_gr_combined77_2,
     file = "mod_val_gr_combined77")

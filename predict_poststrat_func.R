#functions for prediction and poststratification
#model prediction functions
#will need to update based on variables added 
combined_pred <- function(df_pred, draws, hivstat, gni, pgm, age_pred, time_length_age_pred, gni_pred, multi_time_length){
  iter <- nrow(draws[[1]])
  #iter <- 1500
  
  df_3y <- matrix(NA, nrow = nrow(df_pred), ncol = iter)
  df_ever <- matrix(NA, nrow = nrow(df_pred), ncol = iter)
  
  #link country to re_country and rs_hiv value - since we're adding countries that weren't originally predicted
  predicted_countries <- unique(select(filter(df_stan_all_count, predicted == "Y"), c(Country, Region)))
  pred_re_country <- data.frame(t(draws$re_country)[,1:iter], predicted_countries[order(predicted_countries$Country),])
  pred_rs_hiv <- data.frame(t(draws$rs_hiv_country)[,1:iter], predicted_countries[order(predicted_countries$Country),])
  
  #add values for countries that don't have a predicted value - adding it here so it will change for each iteration
  non_predicted_countries <- unique(select(filter(df_stan_all_count, predicted == "N"), c(Country, Region)))
  
  non_pred_draws <- non_pred_draws_hiv <- matrix(data = NA, nrow = nrow(non_predicted_countries), ncol = iter)
  for(i in 1:iter){
    re_country <- rnorm(n = nrow(non_predicted_countries), mean = 0, sd = draws$sd_re_country[i])
    rs_hiv <- rnorm(n = nrow(non_predicted_countries), mean = 0, sd = draws$sd_rs_country_hiv[i])
    non_pred_draws[,i] <- re_country
    non_pred_draws_hiv[,i] <- rs_hiv
  }
  
  non_pred_re_country <- data.frame(non_pred_draws, non_predicted_countries)
  df_pred_re_country <- rbind(pred_re_country, non_pred_re_country)
  
  non_pred_rs_hiv <- data.frame(non_pred_draws_hiv, non_predicted_countries)
  df_pred_rs_hiv <- rbind(pred_rs_hiv, non_pred_rs_hiv)
  
  #predict
    for(p in 1:nrow(df_pred)){
      print(p)
      #recent
      pred_logit_prob_3y <- as.numeric(draws$alpha[1:iter] + 
                                         draws$re_region[1:iter, df_pred$reg_num[p]] +
                                         draws$rs_time[1:iter, df_pred$reg_num[p]] * df_pred$Year[p] +
                                         #draws$beta_time[i] * df_pred$Year[p] +
                                         filter(df_pred_re_country, Country == df_pred$Country[p])[1:iter])
      
      if(multi_time_length == T){
        for(a in 1:ncol(time_length_age_pred)){
          pred_logit_prob_3y = pred_logit_prob_3y + draws$beta_time_length[1:iter,a]*time_length_age_pred[p,a]
        }
      } else{
        pred_logit_prob_3y = pred_logit_prob_3y + draws$beta_time_length[1:iter]
      }
      
      for(a in 1:ncol(age_pred)){
        pred_logit_prob_3y = pred_logit_prob_3y + draws$beta_age[1:iter,a]*age_pred[p,a]
      }
       
      #ever
      pred_logit_prob_ever <- as.numeric(draws$alpha[1:iter] + 
                                           draws$re_region[1:iter, df_pred$reg_num[p]] +
                                           draws$rs_time[1:iter, df_pred$reg_num[p]] * df_pred$Year[p] +
                                           #draws$beta_time[i] * df_pred$Year[p] +
                                           filter(df_pred_re_country, Country == df_pred$Country[p])[1:iter]) #no beta_time_length value for ever
      
      for(a in 1:ncol(age_pred)){
        pred_logit_prob_ever = pred_logit_prob_ever + draws$beta_age[1:iter,a]*age_pred[p,a]
      }
      if(hivstat == T){
        pred_logit_prob_3y = as.numeric(pred_logit_prob_3y + 
                                          (draws$beta_hiv[i:iter] + 
                                             draws$rs_hiv_region[i:iter, df_pred$reg_num[p]] + 
                                          filter(df_pred_rs_hiv, Country == df_pred$Country[p])[1:iter]) * df_pred$hivstat[p])
        
        pred_logit_prob_ever = as.numeric(pred_logit_prob_ever + 
                                            (draws$beta_hiv[i:iter] + 
                                            draws$rs_hiv_region[i:iter, df_pred$reg_num[p]] +
                                            filter(df_pred_rs_hiv, Country == df_pred$Country[p])[1:iter]) * df_pred$hivstat[p]) 
                                            #draws$beta_recall_hiv[i] * df_pred$hivstat[p] 
      }
      
      if(gni == T){
        for(g in 1:G){
          pred_logit_prob_3y = pred_logit_prob_3y + draws$beta_gni[1:iter,g]*gni_pred[p,g]
          pred_logit_prob_ever = pred_logit_prob_ever + draws$beta_gni[1:iter,g]*gni_pred[p,g]
        }
      }
      
      if(pgm == T){
        pred_logit_prob_3y = pred_logit_prob_3y + draws$beta_pgm[i]*df_pred$pgm[p]
        pred_logit_prob_ever = pred_logit_prob_ever + draws$beta_pgm[i]*df_pred$pgm[p]
        # pred_logit_prob_3y = pred_logit_prob_3y + 
        #   draws$rs_pgm[1:iter, df_pred$count_num[p]] * df_pred$pgm[p]
        # 
        # pred_logit_prob_ever = pred_logit_prob_ever + 
        #   draws$rs_pgm[1:iter, df_pred$count_num[p]] * df_pred$pgm[p]  
      }
      df_3y[p, ] = plogis(as.numeric(pred_logit_prob_3y)) 
      df_ever[p, ] = plogis(as.numeric(pred_logit_prob_ever)) 
    }
  
  df <- list(df_3y, df_ever)
  return(df)
}

#function to obtain weights
weights <- function(agegr, hiv, years){
  pop_age <- read_excel("pop_age0020.xlsx", sheet = 2) #i should change this so i can use the age structure for every year
  colnames(pop_age) <- pop_age[1,]
  pop_age <- pop_age[-c(1:3), -c(1, 3, 5)]
  pop_age[,-c(1:2)] <- lapply(pop_age[,-c(1:2)], as.numeric)
  #pop_age <- pop_age[pop_age$Age %in% agegr, ] #rmbr to use the age groups that were used for the model even if the post-model analysis uses other age groups
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
  
  #isolate to only those of a certain age and year
  pop_age_long <- pop_age_long[pop_age_long$Year %in% years & pop_age_long$agegr %in% agegr,]
  df_weights_long <- pop_age_long
  df_weights_long$poststrat_weight <- df_weights_long$f_pop
  
  # add weights by hiv prevalence
  # add hiv prevalence data - can try to obtain values for 24-49 by using 15-49 and 15-24
  # right now I'll just add the 15+ or 15-49
  if(hiv == T){
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
    
    f_hivprev <- left_join(f_hivprev_1549, f_hivprev_1524)
    f_hivprev <- left_join(pop_age, f_hivprev)
    
    f_hivprev$hivprev2549 <- (f_hivprev$hivprev1549*f_hivprev$f_pop_1549 - f_hivprev$hivprev1524*f_hivprev$f_pop_1524)/f_hivprev$f_pop_2549
    
    df_weights <- left_join(pop_age_long, f_hivprev[f_hivprev$Year %in% years,])
    df_weights$hivpos_pop <- df_weights$hivprev2549*df_weights$f_pop
    df_weights$hivneg_pop <- df_weights$f_pop - df_weights$hivprev2549*df_weights$f_pop
    df_weights_long <- pivot_longer(df_weights, cols = c("hivpos_pop", "hivneg_pop"), names_to = "hivstat_1", values_to = "hiv_pop")
    df_weights_long$hivstat <- ifelse(df_weights_long$hivstat_1 == "hivpos_pop", 1, 0)
    
    #df_weights_long <- select(df_weights_long, -"Year")
    df_weights_long$poststrat_weight <- df_weights_long$hiv_pop
  }
  
  df_weights_long$agegr <- as.numeric(factor(df_weights_long$agegr))
  df_weights_long$Year <- df_weights_long$Year - 2000
  return(df_weights_long)
}

#to create the weighted probability matrices - uses combined_pred
pred_weighted <- function(draws, df_stan, Country, Year, agegr, hivstat_log, gni_log, pgm_log, multi_time_length, df_weights_long){
  df_pred <- expand.grid(Country = Country, 
                         Year = Year,
                         hivstat = 0:1,
                         agegr = agegr)
  
  #add region to df_pred
  #region_country <- unique(df_stan[,c("Country", "count_num", "Region", "reg_num")])
  region_country$Region[region_country$Region == "Central" | region_country$Region == "Western"] <- "Central/Western"
  df_pred <- left_join(df_pred, region_country)
  df_pred$reg_num <- as.numeric(factor(df_pred$Region))
  
  #add whether Country is predicted or not
  predicted <- unique(df_stan[,c("Country", "predicted")])
  df_pred <- left_join(df_pred, predicted)
  
  n_age <- count(!is.na(unique(df_stan$agegr)))
  for(a in 1:n_age){
    df_pred$x[df_pred$agegr == a] <- 1
    df_pred$x[df_pred$agegr != a] <- 0
    name <- paste0("agegr", a)
    names(df_pred)[names(df_pred) == "x"] <- name
  }
  
  if(n_age == 5){
    age_pred <- as.matrix(df_pred[,c("agegr2", "agegr3", "agegr4", "agegr5")])
    time_length_age_pred <- as.matrix(df_pred[,c("agegr1", "agegr2", "agegr3", "agegr4", "agegr5")])
    
  } else{stop("add new age options")}
  
  A <- n_age-1
  
  #add pgm to df_pred
  if(pgm_log == T){
    any_scr_pgm <- filter(screening_program, pgm_type == "Any")
    any_scr_pgm$screen_end[any_scr_pgm$screen_end == "Present"] <- 2020
    any_scr_pgm$screen_end <- as.numeric(any_scr_pgm$screen_end)
    any_scr_pgm$screen_start <- any_scr_pgm$screen_start - 2000
    any_scr_pgm$screen_end <- any_scr_pgm$screen_end - 2000
    
    for(i in 1:nrow(df_pred)){
      c <- df_pred$Country[i]
      pgm <- filter(any_scr_pgm, Country == c)
      if(nrow(pgm) == 1){
        if(df_pred$Year[i] >= pgm$screen_start & df_pred$Year[i] <= pgm$screen_end){
          df_pred$pgm[i] <- 1
        }else{
          df_pred$pgm[i] <- 0
        }
      }else if(nrow(pgm) == 2){
        if(df_pred$Year[i] >= pgm$screen_start[1] & df_pred$Year[i] <= pgm$screen_end[1]){
          df_pred$pgm[i] <- 1
        }else if(df_pred$Year[i] >= pgm$screen_start[2] & df_pred$Year[i] <= pgm$screen_end[2]){
          df_pred$pgm[i] <- 1
        }else{
          df_pred$pgm[i] <- 0
        }
      }else{
        df_pred$pgm[i] <- 0
      }
    }
  }
  
  #add gni to df_pred
  if(gni_log == T){
    gni_raw <- read.csv("gni.csv")
    gni <- pivot_longer(gni_raw, cols = -c(1:4), names_to = "Year", values_to = "gni")
    gni$Year <- as.numeric(sub(".", "", gni$Year)) - 2000
    names(gni)[names(gni) == "Country.Name"] <- "Country"
    gni <- gni[,-c(2:4)]
    gni$Country[gni$Country == "Cabo Verde"] <- "Cape Verde"
    gni$Country[gni$Country == "Congo, Rep."] <- "Congo (Rep)"
    gni$Country[gni$Country == "Gambia, The"] <- "Gambia"
    df_pred <- left_join(df_pred, gni)
    #create income group variable
    gni_thresholds <- read_excel("gni_thresholds.xlsx")
    gni_thresholds <- pivot_wider(gni_thresholds, names_from = Class, values_from = Max_Threshold)
    gni_thresholds$Year <- gni_thresholds$Year - 2000
    df_pred <- left_join(df_pred, gni_thresholds)
    df_pred$income_group[df_pred$gni < df_pred$Low_LowerMiddle] <- "Low"
    df_pred$income_group[df_pred$gni >= df_pred$Low_LowerMiddle & df_pred$gni <= df_pred$LowerMiddle_UpperMiddle] <- "LowerMiddle"
    df_pred$income_group[df_pred$gni > df_pred$LowerMiddle_UpperMiddle & df_pred$gni <= df_pred$UpperMiddle_High] <- "UpperMiddle"
    df_pred$income_group[df_pred$gni > df_pred$UpperMiddle_High] <- "High"
    df_pred$income_group_num <- as.numeric(factor(df_pred$income_group, levels = c("Low", "LowerMiddle", "UpperMiddle", "High")))
    
    G <- length(unique(df_pred$income_group_num))
    for(g in 2:G){
      df_pred$x[df_pred$income_group_num == g] <- 1
      df_pred$x[df_pred$income_group_num != g] <- 0
      name <- paste0("income_group", g)
      names(df_pred)[names(df_pred) == "x"] <- name
    }
    
    G <- G-1 #subtract one for the model as there will be n-1 parameters
    
    gni_pred <- as.matrix(df_pred[,c("income_group2", "income_group3")])
  }
  df <- combined_pred(df_pred, draws, hivstat = hivstat_log, gni = gni_log, pgm = pgm_log, age_pred = age_pred, time_length_age_pred = time_length_age_pred, gni_pred = gni_pred, multi_time_length = multi_time_length) #from my predict_poststrat_func.R file
  pred_prob_3y <- df[[1]]
  pred_prob_ever <- df[[2]]
  
  dat_pred <- left_join(df_pred, df_weights_long)
  
  pred_prob_weighted_3y <- sweep(pred_prob_3y, 1, dat_pred$poststrat_weight, "*")
  pred_prob_weighted_ever <- sweep(pred_prob_ever, 1, dat_pred$poststrat_weight, "*")
  
  pred_prob_weighted <- list(pred_prob_weighted_3y, pred_prob_weighted_ever, dat_pred)
  
  return(pred_prob_weighted)
}

# to find the probabilities by certain variables - will need to update based on num variables used
prob_by_var <- function(pred_prob, dat_pred, numvar, var1, var2, var3, var4, count_restrict = F, countries){ #will add the if statements for the number of variables as they arise
  
  df <- cbind(dat_pred, pred_prob)
  
  if(count_restrict == T){
    df <- df[df$Country %in% countries,]
    dat_pred <- dat_pred[dat_pred$Country %in% countries,]
  }
  
  if(numvar == 0){
    df_by_var <- pred_prob
    df_by_var <- rbind(df_by_var, colSums(pred_prob)/sum(dat_pred$poststrat_weight))
    return(df_by_var)
  }
  if(numvar == 1){
    list_by_var <- split(df, list(dat_pred[,var1]))
    
    iter <- ncol(pred_prob) #number of iterations
    col <- ncol(list_by_var[[1]]) #number of columns in combined df
    pred_index <- col - iter + 1 #column number where predictions start
    
    df_by_var <- data.frame(x = names(list_by_var))
    names(df_by_var)[names(df_by_var) == "x"] <- var1
    for(i in 1:length(list_by_var)){
      df_by_var[i,2:(iter+1)] <- colSums(list_by_var[[i]][,pred_index:col])/sum(list_by_var[[i]]$poststrat_weight)
    }
    return(df_by_var)
  }
  if(numvar == 2){
    list_by_var <- split(df, list(dat_pred[,var1], dat_pred[,var2]))
    
    iter <- ncol(pred_prob) #number of iterations
    col <- ncol(list_by_var[[1]]) #number of columns in combined df
    pred_index <- col - iter + 1 #column number where predictions start
    
    df_by_var <- data.frame(x = names(list_by_var))
    df_by_var <- tidyr::extract(df_by_var, x, c(var1, var2), regex = "^([^.]+)\\.(.*)")
    for(i in 1:length(list_by_var)){
      df_by_var[i,3:(iter+2)] <- colSums(list_by_var[[i]][,pred_index:col])/sum(list_by_var[[i]]$poststrat_weight)
    }
    return(df_by_var)
  }
  if(numvar == 3){
    list_by_var <- split(df, list(dat_pred[,var1], dat_pred[,var2], dat_pred[,var3]))
    
    iter <- ncol(pred_prob) #number of iterations
    col <- ncol(list_by_var[[1]]) #number of columns in combined df
    pred_index <- col - iter + 1 #column number where predictions start
    
    df_by_var <- data.frame(x = names(list_by_var))
    df_by_var <- tidyr::extract(df_by_var, x, c(var1, var2), regex = "^([^.]+)\\.(.*)")
    df_by_var <- tidyr::extract(df_by_var, var2, c(var2, var3), regex = "^([^.]+)\\.(.*)")
    for(i in 1:length(list_by_var)){
      df_by_var[i,4:(iter+3)] <- colSums(list_by_var[[i]][,pred_index:col])/sum(list_by_var[[i]]$poststrat_weight)
    }
    return(df_by_var)
  }
  if(numvar == 4){
    list_by_var <- split(df, list(dat_pred[,var1], dat_pred[,var2], dat_pred[,var3], dat_pred[,var4]))
    
    iter <- ncol(pred_prob) #number of iterations
    col <- ncol(list_by_var[[1]]) #number of columns in combined df
    pred_index <- col - iter + 1 #column number where predictions start
    
    df_by_var <- data.frame(x = names(list_by_var))
    df_by_var <- tidyr::extract(df_by_var, x, c(var1, var2), regex = "^([^.]+)\\.(.*)")
    df_by_var <- tidyr::extract(df_by_var, var2, c(var2, var3), regex = "^([^.]+)\\.(.*)")
    df_by_var <- tidyr::extract(df_by_var, var3, c(var3, var4), regex = "^([^.]+)\\.(.*)")
    for(i in 1:length(list_by_var)){
      df_by_var[i,5:(iter+4)] <- colSums(list_by_var[[i]][,pred_index:col])/sum(list_by_var[[i]]$poststrat_weight)
    }
    return(df_by_var)
  }
}

# to find and output the mean and 95% CrI for each row
row_sumstats <- function(df, prob){
  df$mean <- apply(df, 1, mean)
  cri <- rowQuantiles(as.matrix(df), probs = c(prob[1], prob[2], prob[3]))
  df_to_return <- NULL
  if(length(cri) == length(prob)){
    df_to_return <- cbind(df, t(cri))
  } else{
    df_to_return <- cbind(df, cri)
  }
  return(df_to_return)
}

# function which groups the data by their characteristics to get the final proportions
prob_final <- function(pred_prob, dat_pred, numvar, var1, var2, var3, var4, count_restrict = F, countries, n_agegr, agegr_labels){

  prob_summary = NULL
  region_country <- read_excel("surveys.xlsx", sheet = 7)
  if(numvar == 1){
    #obtaining weighted probabilities by variable characteristics
    prob <- prob_by_var(pred_prob, dat_pred, numvar, var1, count_restrict = count_restrict, countries = countries)
    
    #finding mean and quantiles for each combination of variables
    prob_summary <- cbind(prob[,1:numvar], row_sumstats(prob[,(numvar+1):ncol(prob)], prob = c(0.5, 0.025, 0.975)))
    names(prob_summary)[1] <- var1
    
    #changing value of year to numeric (if variable for year exists)
    if(var1 == "Year"){
      prob_summary$Year <- as.numeric(prob_summary$Year) + 2000
    }
    
    #adding regions if it's not one of the variables
    if(var1 == "Country"){
      #region_country <- unique(surveys[,c("Region", "Country")])
      region_country$Region[region_country$Region == "Western" | region_country$Region == "Central"] <- "Central/Western"
      prob_summary <- left_join(prob_summary, region_country)
    }
    
    #cleaning HIV status labels if one of the variables is HIV
    if(var1 == "hivstat"){
      prob_summary$hivstat <- factor(prob_summary$hivstat, labels = c("HIV-", "HIV+"))
    }
    
    #cleaning age group labels if one of the variables is agegr
    if(var1 == "agegr"){
      prob_summary$agegr <- factor(prob_summary$agegr, levels = n_agegr,
                                          labels = agegr_labels)
    }
  } else if(numvar == 2){
    #obtaining weighted probabilities by variable characteristics
    prob <- prob_by_var(pred_prob, dat_pred, numvar, var1, var2, count_restrict = count_restrict, countries = countries)
    
    #finding mean and quantiles for each combination of variables
    prob_summary <- cbind(prob[,1:numvar], row_sumstats(prob[,(numvar+1):ncol(prob)], prob = c(0.5, 0.025, 0.975)))
    
    #changing value of year to numeric (if variable for year exists)
    if(var1 == "Year" | var2 == "Year"){
      prob_summary$Year <- as.numeric(prob_summary$Year) + 2000
    }
    
    #adding regions if it's not one of the variables
    if(var1 == "Country" | var2 == "Country"){ #I have it as == Country here and not != Region in case one of the variables isn't Country --> this is possible for numvar == 2 but not as possible for numvar == 3+
      #region_country <- unique(surveys[,c("Region", "Country")])
      region_country$Region[region_country$Region == "Western" | region_country$Region == "Central"] <- "Central/Western"
      prob_summary <- left_join(prob_summary, region_country)
    }
    
    #cleaning HIV status labels if one of the variables is HIV
    if(var1 == "hivstat" | var2 == "hivstat"){
      prob_summary$hivstat <- as.numeric(prob_summary$hivstat)
      prob_summary$hivstat <- factor(prob_summary$hivstat, labels = c("HIV-", "HIV+"))
    }
    
    #cleaning age group labels if one of the variables is agegr
    if(var1 == "agegr" | var2 == "agegr"){
      prob_summary$agegr <- factor(prob_summary$agegr, levels = n_agegr,
                                   labels = agegr_labels)
    }
  } else if(numvar == 3){
    #obtaining weighted probabilities by variable characteristics
    prob <- prob_by_var(pred_prob, dat_pred, numvar, var1, var2, var3, count_restrict = count_restrict, countries = countries)
    
    #finding mean and quantiles for each combination of variables
    prob_summary <- cbind(prob[,1:numvar], row_sumstats(prob[,(numvar+1):ncol(prob)], prob = c(0.5, 0.025, 0.975)))
    
    #changing value of year to numeric (if variable for year exists)
    if(var1 == "Year" | var2 == "Year" | var3 == "Year"){
      prob_summary$Year <- as.numeric(prob_summary$Year) + 2000
    }
    
    #adding regions if it's not one of the variables
    if(var1 != "Region" & var2 != "Region" & var3 != "Region"){
      #region_country <- unique(surveys[,c("Region", "Country")])
      region_country$Region[region_country$Region == "Western" | region_country$Region == "Central"] <- "Central/Western"
      prob_summary <- left_join(prob_summary, region_country)
    }
    
    #cleaning HIV status labels if one of the variables is HIV
    if(var1 == "hivstat" | var2 == "hivstat" | var3 == "hivstat"){
      prob_summary$hivstat <- as.numeric(prob_summary$hivstat)
      prob_summary$hivstat <- factor(prob_summary$hivstat, labels = c("HIV-", "HIV+"))
    }
    
    #cleaning age group labels if one of the variables is agegr
    if(var1 == "agegr" | var2 == "agegr" | var3 == "agegr"){
      prob_summary$agegr <- factor(prob_summary$agegr, levels = n_agegr,
                                   labels = agegr_labels)
    }
  } else if(numvar == 4){
    #obtaining weighted probabilities by variable characteristics
    prob <- prob_by_var(pred_prob, dat_pred, numvar, var1, var2, var3, var4, count_restrict = count_restrict, countries = countries)
    
    #finding mean and quantiles for each combination of variables
    prob_summary <- cbind(prob[,1:numvar], row_sumstats(prob[,(numvar+1):ncol(prob)], prob = c(0.5, 0.025, 0.975)))
    
    #changing value of year to numeric (if variable for year exists)
    if(var1 == "Year" | var2 == "Year" | var3 == "Year" | var4 == "Year"){
      prob_summary$Year <- as.numeric(prob_summary$Year) + 2000
    }
    
    #adding regions if it's not one of the variables
    if(var1 != "Region" & var2 != "Region" & var3 != "Region" & var4 != "Region"){
      #region_country <- unique(surveys[,c("Region", "Country")])
      region_country$Region[region_country$Region == "Western" | region_country$Region == "Central"] <- "Central/Western"
      prob_summary <- left_join(prob_summary, region_country)
    }
    
    #cleaning HIV status labels if one of the variables is HIV
    if(var1 == "hivstat" | var2 == "hivstat" | var3 == "hivstat" | var4 == "hivstat"){
      prob_summary$hivstat <- as.numeric(prob_summary$hivstat)
      prob_summary$hivstat <- factor(prob_summary$hivstat, labels = c("HIV-", "HIV+"))
    }
    
    #cleaning age group labels if one of the variables is agegr
    if(var1 == "agegr" | var2 == "agegr" | var3 == "agegr" | var4 == "agegr"){
      prob_summary$agegr <- factor(prob_summary$agegr, levels = n_agegr,
                                   labels = agegr_labels)
    }
  }
  return(prob_summary)
}

#survey probabilities
survey_prob <- function(timing, gen, hiv, df_to_add){
  raw_prob <- list(NA)
  n = 1
  if(timing == "3y"){
    load("CC Testing Past 3Y")
    if(gen == T){
      raw_prob_3y <- testpast3y_df[!is.na(testpast3y_df$Region) & 
                                     !is.na(testpast3y_df$Mean), c(1:5, 8, 10:12)]
      names(raw_prob_3y)[names(raw_prob_3y) == "Mean"] <- "raw_mean"
      names(raw_prob_3y)[names(raw_prob_3y) == "Lower"] <- "raw_lower"
      names(raw_prob_3y)[names(raw_prob_3y) == "Upper"] <- "raw_upper"
      raw_prob[[n]] <- raw_prob_3y
      n = n + 1
    }
    if(hiv == T){
      hiv_raw_prob_3y <- hiv_testpast3y_df[!is.na(hiv_testpast3y_df$Region), c(1:6, 9, 11:13)]
      names(hiv_raw_prob_3y)[names(hiv_raw_prob_3y) == "Mean"] <- "raw_mean"
      names(hiv_raw_prob_3y)[names(hiv_raw_prob_3y) == "Lower"] <- "raw_lower"
      names(hiv_raw_prob_3y)[names(hiv_raw_prob_3y) == "Upper"] <- "raw_upper"
      names(hiv_raw_prob_3y)[names(hiv_raw_prob_3y) == "Status"] <- "hivstat"
      raw_prob[[n]] <- hiv_raw_prob_3y
      n = n + 1
    }
  } else if(timing == "ever"){
    if(gen == T){
      load("CC Testing Ever")
      raw_prob_ever <- evertest_df[!is.na(evertest_df$Region) &
                                     !is.na(evertest_df$Mean), c(1:5, 7, 9:11)]
      
      names(raw_prob_ever)[names(raw_prob_ever) == "Mean"] <- "raw_mean"
      names(raw_prob_ever)[names(raw_prob_ever) == "Lower"] <- "raw_lower"
      names(raw_prob_ever)[names(raw_prob_ever) == "Upper"] <- "raw_upper"
      raw_df_to_add <- df_to_add[is.na(df_to_add$hivstat),]
      raw_df_to_add$raw_mean <- raw_df_to_add$num/raw_df_to_add$den
      for(i in 1:nrow(raw_df_to_add)){
        raw_df_to_add$raw_lower[i] <- prop.test(raw_df_to_add$num[i], raw_df_to_add$den[i])$conf.int[1]
        raw_df_to_add$raw_upper[i] <- prop.test(raw_df_to_add$num[i], raw_df_to_add$den[i])$conf.int[2]
      }
      region_country <- unique(surveys[,c("df", "Country", "Region", "Year", "Survey", "ISO3")])
      raw_df_to_add <- left_join(raw_df_to_add, region_country)
      col_to_add <- colnames(raw_prob_ever)
      raw_prob_ever <- rbind(raw_prob_ever, raw_df_to_add[,c(col_to_add)])
      raw_prob[[n]] <- raw_prob_ever
      n = n + 1
    }
    if(hiv == T){
      hiv_raw_prob_ever <- hiv_evertest_df[!is.na(hiv_evertest_df$Region),]
      names(hiv_raw_prob_ever)[names(hiv_raw_prob_ever) == "Mean"] <- "raw_mean"
      names(hiv_raw_prob_ever)[names(hiv_raw_prob_ever) == "Lower"] <- "raw_lower"
      names(hiv_raw_prob_ever)[names(hiv_raw_prob_ever) == "Upper"] <- "raw_upper"
      names(hiv_raw_prob_ever)[names(hiv_raw_prob_ever) == "Status"] <- "hivstat"
      hiv_raw_df_to_add <- df_to_add[!is.na(df_to_add$hivstat),]
      #for a weighted average for 25-49 from the zw20phia tabulations
      a <- hiv_raw_df_to_add[hiv_raw_df_to_add$agegr %in% c("25-29", "30-34", "35-39", "40-44", "45-49"),]
      hiv_raw_df_to_add[(nrow(hiv_raw_df_to_add) + 1),] <- hiv_raw_df_to_add[nrow(hiv_raw_df_to_add),]
      hiv_raw_df_to_add[nrow(hiv_raw_df_to_add),c("agegr", "prob", "denominator_original", "numerator_original","den", "num")] <- c("25-49", NA, sum(a$denominator_original), sum(a$numerator_original), sum(a$den), sum(a$num))
      hiv_raw_df_to_add$den <- as.numeric(hiv_raw_df_to_add$den)
      hiv_raw_df_to_add$num <- as.numeric(hiv_raw_df_to_add$num)
      
      hiv_raw_df_to_add$raw_mean <- hiv_raw_df_to_add$num/hiv_raw_df_to_add$den
      for(i in 1:nrow(hiv_raw_df_to_add)){
        hiv_raw_df_to_add$raw_lower[i] <- prop.test(hiv_raw_df_to_add$num[i], hiv_raw_df_to_add$den[i])$conf.int[1]
        hiv_raw_df_to_add$raw_upper[i] <- prop.test(hiv_raw_df_to_add$num[i], hiv_raw_df_to_add$den[i])$conf.int[2]
      }
      region_country <- unique(surveys[,c("df", "Country", "Region", "Year", "Survey", "ISO3")])
      hiv_raw_df_to_add <- left_join(hiv_raw_df_to_add, region_country)
      hiv_raw_df_to_add$hivstat <- "HIV+"
      col_to_add <- colnames(hiv_raw_prob_ever)
      hiv_raw_prob_ever <- rbind(hiv_raw_prob_ever, hiv_raw_df_to_add[,c(col_to_add)])
      
      raw_prob[[n]] <- hiv_raw_prob_ever
      n = n + 1
    }
  }
  
  for(i in 1:length(raw_prob)){
    raw_prob[[i]] <- select(raw_prob[[i]], -Region)
    surveys <- read_excel("surveys.xlsx", sheet = 1)
    region_country <- surveys[,c("Country", "Region")]
    raw_prob[[i]] <- left_join(raw_prob[[i]], region_country)
    
    raw_prob[[i]]$Region[raw_prob[[i]]$Region == "Western" | raw_prob[[i]]$Region == "Central"] <- "Central/Western"
  }
  return(raw_prob)
}

#to combine the 3y and ever dataframes
dat_combine <- function(dat_ever, dat_3y){
  dat_ever$time_length <- "ever"
  dat_3y$time_length <- "3y"
  
  dat_combine <- rbind(dat_3y, dat_ever)
  
  return(dat_combine)
}


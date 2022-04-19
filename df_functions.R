#functions to create dataframes
#for dataframes stratified by age only
df_age_fig <- function(age.groups, df, func, df_name){
  k=1
  for(i in 1:nrow(age.groups)){
    n <- nrow(func(df, age.groups[1,1], age.groups[1,2]))
    m <- ncol(func(df, age.groups[1,1], age.groups[1,2]))
    df_name[k:(k+n-1), c(2:(m+1))] <- func(df, age.groups[i,1], age.groups[i,2])
    k=k+n
  }
  return(df_name)
}

df_age_xl <- function(age.groups, df, func, df_name){
  for(i in 1:nrow(age.groups)){
    prob <- func(df, age.groups[i,1], age.groups[i,2])
    TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
    TMP <- paste0(format(prob[1],digits=3,nsmall=1)," [",TMP,"]")
    df_name[i, 2] <- TMP
  }
  return(df_name)
}

#for dataframes stratified by age and sexual activity level
df_age_sa_fig <- function(age.groups, sa.level, df, func, df_name){
  n = 1
  for(i in 1:nrow(age.groups)){
    for(k in 1:nrow(sa.level)){
      df_name[n, c(3:5)] <- func(df, age.groups[i, 1], age.groups[i, 2], sa.level[k, 1], sa.level[k, 2])
      n = n+1
    }
  }
  return(df_name)
}

df_age_sa_xl <- function(age.groups, sa.level, df, func, df_name){
  for(i in 1:nrow(age.groups)){
    for(j in 1:nrow(sa.level)){
      prob <- func(df, age.groups[i, 1], age.groups[i, 2], sa.level[j, 1], sa.level[j, 2])
      TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
      TMP <- paste0(format(prob[1],digits=3,nsmall=1)," [",TMP,"]")
      df_name[j, i+1] <- TMP
    }
  }
  return(df_name)
}

#for dataframes stratified by age, and sex level
df_age_sex_fig <- function(age.groups, df, func, df_name){
  n = 1
  for(i in 1:nrow(age.groups)){
    for(j in 1:length(df)){
        df_name[n, c(4:6)] <- func(df[[j]], age.groups[i, 1], age.groups[i, 2])
        n = n+1
    }
  }
  return(df_name)
}

df_age_sex_xl <- function(age.groups, sa.level, df, func, df_name){
  n = 2
  for(i in 1:nrow(age.groups)){
    m = 1
    for(j in 1:length(df)){
        prob <- func(df[[j]], age.groups[i, 1], age.groups[i, 2])
        TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
        TMP <- paste0(format(prob[1],digits=3,nsmall=1)," [",TMP,"]")
        df_name[m, n] <- TMP
        m = m+1
      }
      n = n+1
  } 
  return(df_name)
}

#for dataframes stratified by age, sex and sexual activity level
df_age_sex_sa_fig <- function(age.groups, sa.level, df, func, df_name){
  n = 1
  for(i in 1:nrow(age.groups)){
    for(j in 1:length(df)){
      for(k in 1:nrow(sa.level)){
        df_name[n, c(4:6)] <- func(df[[j]], age.groups[i, 1], age.groups[i, 2], sa.level[k, 1], sa.level[k, 2])
        n = n+1
      }
    }
  }
  return(df_name)
}

df_age_sex_sa_xl <- function(age.groups, sa.level, df, func, df_name){
  n = 2
  for(i in 1:nrow(age.groups)){
    for(j in 1:length(df)){
      m = 1
      for(k in 1:nrow(sa.level)){
        prob <- func(df[[j]], age.groups[i, 1], age.groups[i, 2], sa.level[k, 1], sa.level[k, 2])
        TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
        TMP <- paste0(format(prob[1],digits=3,nsmall=1)," [",TMP,"]")
        df_name[m, n] <- TMP
        m = m+1
      }
      n = n+1
    }
  } 
  return(df_name)
}

df_age_sex_sa_prior <- function(age.groups, sa.level, df, func, df_name){
  n = 3
  for(i in 1:nrow(age.groups)){
    m = 1
    for(k in 1:nrow(sa.level)){
      for(j in 1:length(df)){
        prob <- func(df[[j]], age.groups[i, 1], age.groups[i, 2], sa.level[k, 1], sa.level[k, 2])
        TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
        df_name[m, n] <- TMP
        m = m+1
      }
    }
    n = n+1
  } 
  return(df_name)
}

#for dataframes stratified by 2 variables
df_2_fig <- function(var1, var2, df, func, df_name){
  k = 1
  for(i in 1:length(var2)){
    for(j in 1:nrow(var1)){
      df_name[k, c(3:5)] <- func(df, var1[j, 1], var1[j, 2], var2[i])
      k = k+1
    }
  }
  return(df_name)
}

df_2_xl <- function(var1, var2, df, func, df_name){
  for(i in 1:length(var2)){
    for(j in 1:nrow(var1)){
      prob <- func(df, var1[j, 1], var1[j, 2], var2[i])
      TMP <- paste(format(prob[2],digits=3,nsmall=1),format(prob[3],digits=3,nsmall=1),sep=", ")
      TMP <- paste0(format(prob[1],digits=3,nsmall=1)," [",TMP,"]")
      df_name[i, j+1] <- TMP
    }
  }
  return(df_name)
}

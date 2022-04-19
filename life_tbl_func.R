library(abind) #to allow for the combination of two matrices into an array

#life table function
risk_to_rate <- function(risk, time){
  rate = -(log(1 - risk)) / time
  return(rate)
}

rates <- c(rep(0.98, 3), rep(0.98, 3))
lower_age <- 30
upper_age <- 50
dt <- 0.1
#View(life_table_ras(rates[1,], lower_age, upper_age, dt))

# life table function with a realistic age structure (ras)
life_table_ras_lag <- function(rates, 
                           lower_age, upper_age,
                           dt) {
  
  wom <- 1000
  age <- seq(lower_age, upper_age - 1, by = 1)
  
  if (length(rates) != length(age)) { stop("vector for rates should be same size as age") }
  
  start <- 0
  end <- upper_age - lower_age
  time <- seq(start, end, by = dt)
  aged_at <- seq(1, end, by = 1)
  
  #create initial matrix
  lt <- matrix(data = 0, nrow = length(time), ncol = length(age) * 2 + 2 + 1)
  id_never <- paste(age, "never", sep = "_")
  id_once <- paste(age, "once", sep = "_")
  id_never_age <- rev(paste(age, "never", sep = "_"))
  id_once_age <- rev(paste(age, "once", sep = "_"))
  colnames(lt) <- c("time", 
                    id_never,
                    id_once,
                    "twice",
                    "can_test_twice")
  
  #rename rates for code clarity
  id_rates <- paste(age, "rate", sep = "_")
  names(rates) <- id_rates
  
  # assign initial conditions
  lt[1, "30_never"] <- wom
  lt[, "time"] <- time
  
  for (i in 2:nrow(lt)) {
    
    twice_by_age_dt <- 0
    val <- 3/dt
    if(i > val){
      lt[i, "can_test_twice"] <- sum(lt[i - val, id_once]) - lt[i-1, "twice"] #this works because there's only one age group that tests once per row
      if(lt[i, "can_test_twice"] < 0){
        lt[i, "can_test_twice"] <- 0
      }
    }
    for (j in 1:length(age)) {
      id_n <- id_never[j]
      id_o <- id_once[j]
      if(lt[i-1, id_n] != 0){
        lt[i, id_n] <- lt[i - 1, id_n] + dt * (-rates[j] * lt[i - 1, id_n])
        lt[i, id_o] <- lt[i - 1, id_o] + dt * (rates[j] * lt[i - 1, id_n])
        if(i > (val + 1) & j > 3){
          twice_by_age_dt <- twice_by_age_dt + dt * rates[j] * (lt[i, "can_test_twice"])
          lt[i, id_o] <- lt[i - 1, id_o] + dt * rates[j] * lt[i - 1, id_n] - dt * rates[j] * (lt[i, "can_test_twice"])
        }
      }
    }
    lt[i, "twice"] <- lt[i - 1, "twice"] + twice_by_age_dt
    
    if (lt[i, "time"] %in% aged_at) {
      for (j in 2:length(id_once_age)) {
        lt[i, id_never_age[j - 1]] <- lt[i, id_never_age[j]]
        lt[i, id_once_age[j - 1]] <- lt[i, id_once_age[j]]
        lt[i, id_never_age[j]] <- 0
        lt[i, id_once_age[j]] <- 0
      }}
  }
  # rowSums(lt[, colnames(lt) != "time"])
  tested_twice <- lt[, "twice"] / wom
  to_sel <- seq(1, length(time), by = 1 /dt)
  tt_f <- tested_twice[to_sel]
  return(lt)
}

life_table_ras2_lag <- function(rates, 
                           lower_age, upper_age,
                           dt) {
  
  wom <- 1000
  age <- seq(lower_age, upper_age - 1, by = 1)
  
  if (ncol(rates) != length(age)) { stop("vector for rates should be same size as age") }
  
  start <- 0
  end <- upper_age - lower_age
  time <- seq(start, end, by = dt)
  aged_at <- seq(1, end, by = 1)
  
  #create initial matrix
  n_row = length(time)
  n_col = length(age) * 2 + 3
  n_tbl = nrow(rates)
  lt <- array(data = 0, dim = c(n_row, n_col, n_tbl))
  id_never <- paste(age, "never", sep = "_")
  id_once <- paste(age, "once", sep = "_")
  id_never_age <- rev(paste(age, "never", sep = "_"))
  id_once_age <- rev(paste(age, "once", sep = "_"))
  colnames(lt) <- c("time", 
                    id_never,
                    id_once,
                    "twice",
                    "can_test_twice")
  
  #rename rates for code clarity
  id_rates <- paste(age, "rate", sep = "_")
  colnames(rates) <- id_rates

  # assign initial conditions
  lt[1, "30_never",] <- wom
  lt[, "time",] <- time
  
  for (i in 2:nrow(lt)) {
    
    twice_by_age_dt <- 0
    val <- 3/dt
    if(i > val){
      for(k in 1:nrow(rates))
      lt[i, "can_test_twice", k] <- sum(lt[i - val, id_once, k]) - lt[i-1, "twice", k] #this works because there's only one age group that tests once per row
    }
    for (j in 1:length(age)) {
      id_n <- id_never[j]
      id_o <- id_once[j]
      if(lt[i-1, id_n, 1] != 0){
        lt[i, id_n,] <- lt[i - 1, id_n,] + dt * (-rates[,j] * lt[i - 1, id_n,])
        lt[i, id_o,] <- lt[i - 1, id_o,] + dt * (rates[,j] * lt[i - 1, id_n,])
        if(i > (val + 1) & j > 3){
          twice_by_age_dt <- twice_by_age_dt + dt * rates[,j] * (lt[i, "can_test_twice",])
          lt[i, id_o,] <- lt[i - 1, id_o,] + dt * rates[,j] * lt[i - 1, id_n,] - dt * rates[,j] * (lt[i, "can_test_twice",])
        }
      }
    }
    lt[i, "twice",] <- lt[i - 1, "twice",] + twice_by_age_dt
    
    if (lt[i, "time", 1] %in% aged_at) {
      for (j in 2:length(id_once_age)) {
        lt[i, id_never_age[j - 1],] <- lt[i, id_never_age[j],]
        lt[i, id_once_age[j - 1],] <- lt[i, id_once_age[j],]
        lt[i, id_never_age[j],] <- 0
        lt[i, id_once_age[j],] <- 0
      }}
  }
  # rowSums(lt[, colnames(lt) != "time"])
  tested_twice <- lt[, "twice",] / wom
  to_sel <- seq(1, length(time), by = 1 /dt)
  tt_f <- tested_twice[to_sel,]
  return(tt_f)
}

# life table function with a realistic age structure (ras)
life_table_ras <- function(rates,
                           lower_age, upper_age,
                           dt) {

  wom <- 1000
  age <- seq(lower_age, upper_age - 1, by = 1)

  if (length(rates) != length(age)) { stop("vector for rates should be same size as age") }

  start <- 0
  end <- upper_age - lower_age
  time <- seq(start, end, by = dt)
  aged_at <- seq(1, end, by = 1)

  #create initial matrix
  lt <- matrix(data = 0, nrow = length(time), ncol = length(age) * 2 + 2)
  id_never <- paste(age, "never", sep = "_")
  id_once <- paste(age, "once", sep = "_")
  id_never_age <- rev(paste(age, "never", sep = "_"))
  id_once_age <- rev(paste(age, "once", sep = "_"))
  colnames(lt) <- c("time",
                    id_never,
                    id_once,
                    "twice")

  #rename rates for code clarity
  id_rates <- paste(age, "rate", sep = "_")
  names(rates) <- id_rates

  # assign initial conditions
  lt[1, "30_never"] <- wom
  lt[, "time"] <- time

  for (i in 2:nrow(lt)) {

    twice_by_age_dt <- 0
    for (j in 1:length(age)) {
      id_n <- id_never[j]
      id_o <- id_once[j]
      lt[i, id_n] <- lt[i - 1, id_n] + dt * (-rates[j] * lt[i - 1, id_n])
      lt[i, id_o] <- lt[i - 1, id_o] + dt * (rates[j] * lt[i - 1, id_n] - rates[j] * lt[i - 1, id_o])
      twice_by_age_dt <- twice_by_age_dt + dt * (rates[j] * lt[i - 1, id_o])
    }
    lt[i, "twice"] <- lt[i - 1, "twice"] + twice_by_age_dt

    if (lt[i, "time"] %in% aged_at & i != nrow(lt)) {
      for (j in 2:length(id_once_age)) {
        lt[i, id_never_age[j - 1]] <- lt[i, id_never_age[j]]
        lt[i, id_once_age[j - 1]] <- lt[i, id_once_age[j]]
        lt[i, id_never_age[j]] <- 0
        lt[i, id_once_age[j]] <- 0
      }}
  }
  # rowSums(lt[, colnames(lt) != "time"])
  tested_once <- (lt[, "twice"] + lt[, tail(id_once, 1)]) / wom
  tested_twice <- lt[, "twice"] / wom
  to_sel <- seq(1, length(time), by = 1 /dt)
  to_f <- tested_once[to_sel]
  tt_f <- tested_twice[to_sel]
  
  f <- c(tail(to_f, 1), tail(tt_f, 1))
  return(tt_f)
}

# #function which allows the input of an dataframe or matrix for faster computations
life_table_ras2 <- function(rates, beta,
                            lower_age, upper_age,
                            dt,
                            wom_never, wom_once) {

  wom_never <- as.numeric(1000*wom_never)
  wom_once <- as.numeric(1000*wom_once)
  age <- seq(lower_age, upper_age - 1, by = 1)

  if (ncol(rates) != length(age)) { stop("vector for rates should be same size as age") }

  start <- 0
  end <- upper_age - lower_age
  time <- seq(start, end, by = dt)
  aged_at <- seq(1, end, by = 1)

  #create initial matrix
  n_row = length(time)
  n_col = length(age) * 2 + 2
  n_tbl = nrow(rates)
  lt <- array(data = 0, dim = c(n_row, n_col, n_tbl))
  id_never <- paste(age, "never", sep = "_")
  id_once <- paste(age, "once", sep = "_")
  id_never_age <- rev(paste(age, "never", sep = "_"))
  id_once_age <- rev(paste(age, "once", sep = "_"))
  colnames(lt) <- c("time",
                    id_never,
                    id_once,
                    "twice")

  #rename rates for code clarity
  id_rates <- paste(age, "rate", sep = "_")
  colnames(rates) <- id_rates

  # assign initial conditions
  lt[1, "30_never",] <- wom_never
  lt[1, "30_once",] <- wom_once
  lt[, "time",] <- time
  
  #if (table(round(lt[1, "30_never",] + lt[1, "30_once",]) != 1000)) { stop("initial values should equal 1000") }

  for (i in 2:nrow(lt)) {

    twice_by_age_dt <- 0
    for (j in 1:length(age)) {
      id_n <- id_never[j]
      id_o <- id_once[j]
      lt[i, id_n,] <- lt[i - 1, id_n,] + dt * (-rates[,j] * lt[i - 1, id_n,])
      lt[i, id_o,] <- lt[i - 1, id_o,] + dt * (rates[,j] * lt[i - 1, id_n,] - rates[,j]*beta[,j] * lt[i - 1, id_o,])
      twice_by_age_dt <- twice_by_age_dt + dt * (rates[,j]*beta[,j] * lt[i - 1, id_o,])
    }
    lt[i, "twice",] <- lt[i - 1, "twice",] + twice_by_age_dt

    if (lt[i, "time",1] %in% aged_at & i != nrow(lt)) {
      for (j in 2:length(id_once_age)) {
        lt[i, id_never_age[j - 1],] <- lt[i, id_never_age[j],]
        lt[i, id_once_age[j - 1],] <- lt[i, id_once_age[j],]
        lt[i, id_never_age[j],] <- 0
        lt[i, id_once_age[j],] <- 0
      }}
  }
  # rowSums(lt[, colnames(lt) != "time"])
  to_sel <- seq(1, length(time), by = 1 /dt)
  
  tested_once <- (lt[, "twice",] + lt[, tail(id_once, 1),])/ round(lt[1, "30_never",] + lt[1, "30_once",])
  to_f <- tested_once[to_sel,]
  
  tested_twice <- lt[, "twice",] / round(lt[1, "30_never",] + lt[1, "30_once",])
  tt_f <- tested_twice[to_sel,]
  
  f <- data.frame(rbind(tail(to_f, 1),
                        tail(tt_f, 1)))
  return(f)
}

# #function which allows the input of an array for faster computations
#a <- rates

#beta <- array(as.numeric(unlist(beta)), dim = c(nrow(b), ncol(b), 2))
life_table_ras3 <- function(rates, beta,
                            lower_age, upper_age,
                            dt,
                            wom_never, wom_once) {
  wom <- 1000
  wom_never <- as.numeric(wom*wom_never)
  wom_once <- as.numeric(wom*wom_once)
  age <- seq(lower_age, upper_age - 1, by = 1)
  
  if (ncol(rates[,,1]) != length(age)) { stop("vector for rates should be same size as age") }
  
  start <- 0
  end <- upper_age - lower_age
  time <- seq(start, end, by = dt)
  aged_at <- seq(1, end, by = 1)
  
  #create initial matrix
  n_row = length(time)
  n_col = length(age) * 2 + 2
  n_tbl = nrow(rates[,,1])
  n_beta_samp = dim(rates)[3]
  lt <- array(data = 0, dim = c(n_row, n_col, n_tbl, n_beta_samp))
  id_never <- paste(age, "never", sep = "_")
  id_once <- paste(age, "once", sep = "_")
  id_never_age <- rev(paste(age, "never", sep = "_"))
  id_once_age <- rev(paste(age, "once", sep = "_"))
  colnames(lt) <- c("time",
                    id_never,
                    id_once,
                    "twice")
  
  #rename rates for code clarity
  id_rates <- paste(age, "rate", sep = "_")
  colnames(rates) <- id_rates
  
  # assign initial conditions
  lt[1, "30_never",,] <- wom_never
  lt[1, "30_once",,] <- wom_once
  lt[, "time",,] <- time
  
  #if (table(round(lt[1, "30_never",] + lt[1, "30_once",]) != 1000)) { stop("initial values should equal 1000") }
  
  for (i in 2:nrow(lt)) {
    
    twice_by_age_dt <- 0
    for (j in 1:length(age)) {
      id_n <- id_never[j]
      id_o <- id_once[j]
      lt[i, id_n,,] <- lt[i - 1, id_n,,] + dt * (-rates[,j,] * lt[i - 1, id_n,,])
      lt[i, id_o,,] <- lt[i - 1, id_o,,] + dt * (rates[,j,] * lt[i - 1, id_n,,] - rates[,j,]*beta[,j,] * lt[i - 1, id_o,,])
      twice_by_age_dt <- twice_by_age_dt + dt * (rates[,j,]*beta[,j,] * lt[i - 1, id_o,,])
    }
    lt[i, "twice",,] <- lt[i - 1, "twice",,] + twice_by_age_dt
    
    if (lt[i, "time", 1, 1] %in% aged_at & i != nrow(lt)) {
      for (j in 2:length(id_once_age)) {
        lt[i, id_never_age[j - 1],,] <- lt[i, id_never_age[j],,]
        lt[i, id_once_age[j - 1],,] <- lt[i, id_once_age[j],,]
        lt[i, id_never_age[j],,] <- 0
        lt[i, id_once_age[j],,] <- 0
      }}
    print(i)
  }
  # rowSums(lt[, colnames(lt) != "time"])
  to_sel <- seq(1, length(time), by = 1 /dt)
  
  tested_once <- (lt[, "twice",,] + lt[, tail(id_once, 1),,])/wom
  to_f <- tested_once[to_sel,,]
  
  tested_twice <- lt[, "twice",,] /wom
  tt_f <- tested_twice[to_sel,,]
  
  f <- data.frame(rbind(tail(to_f, 1),
                        tail(tt_f, 1)))
  return(f)
}

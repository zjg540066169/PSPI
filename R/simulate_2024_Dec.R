library(mvtnorm) 
library(arm)
library(dplyr)

simulate_standard = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 7
  n_binary <- 3
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  correlation = matrix(
    c(1, 0.2, 0, 0, 0, 0, 0,
      0.2, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0.9, 0, 0, 0,
      0, 0, 0.9, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1),
    nrow = n_continuous, ncol = n_continuous
  )
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = correlation)
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  data$Z = rbinom(n_samples, 1, prop)
  
  data$outcome1 = 
    2 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1 + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2 + #log(0.5 + abs(data$cont_var2 + data$cont_var4)) +
    1 * (10 * data$cont_var1 - 5 * data$cont_var2 -4 * data$cont_var3 - 4 * data$cont_var4 + 8 * data$bin_var1 +
           3 * data$cont_var1 * data$cont_var4 + 2 * data$cont_var3^2 - 1 * data$cont_var2 * data$bin_var1 +
           0.3 * data$cont_var4^3 * data$bin_var1) #+ invlogit(data$cont_var1 * data$cont_var2))
  
  data$outcome0 = 
    0 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1  + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2
  
  if(binary){
    #data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    #data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  pi = invlogit(-4.169 + 1 * data$cont_var1 + 2 * data$cont_var2 + -1 * data$cont_var3 + 1 * data$cont_var4 + 2 * data$bin_var1)
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}


simulate_standard_linear = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 7
  n_binary <- 3
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  correlation = matrix(
    c(1, 0.2, 0, 0, 0, 0, 0,
      0.2, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0.9, 0, 0, 0,
      0, 0, 0.9, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1),
    nrow = n_continuous, ncol = n_continuous
  )
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = correlation)
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  data$Z = rbinom(n_samples, 1, prop)
  
  data$outcome1 = 
    2 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1 + 
    1 * (6 * data$cont_var1 - 5 * data$cont_var2 -4 * data$cont_var3 - 4 * data$cont_var4 + 8 * data$bin_var1) 
  
  data$outcome0 = 
    0 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1
  
  if(binary){
    data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  pi = invlogit(-4.169 + 1 * data$cont_var1 + 2 * data$cont_var2 + -1 * data$cont_var3 + 1 * data$cont_var4 + 2 * data$bin_var1)
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}

#a = sapply(1:200, function(i) sum(simulate_standard_linear(n_samples = 5000, seed = i)$data$selected))
#mean(a)




simulate_simple_quadratic = function(n_samples = 5000, seed = NULL, binary = FALSE){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 10
  n_binary <- 10
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = diag(n_continuous))
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  #data$Z = rbinom(n_samples, 1, 0.5)
  data$Z = rbinom(n_samples, 1, 0.1)
  
  
  # outcome model
  fixed_coefficients <- c(1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1.0, 2, 0, 0, 0, 0, 0) # Coefficients for Z and interactions
  
  # fixed_coefficients <- c(2.0, -1.5, 1.2, 0, 0, 0, 0, 0, 0, 0, 
  #                         2.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  #                         3.0, 0.5, 0.9, -1.4, 0.6, 0, 0) # Coefficients for Z and interactions
  # 
  
  data$outcome1 = as.matrix(data[, 1:20]) %*% fixed_coefficients[1:20] + 
    fixed_coefficients[21] * 1 + 
    fixed_coefficients[22] * 1 * data$cont_var1 + 
    fixed_coefficients[23] * 1 * data$cont_var2 + 
    fixed_coefficients[24] * 1 * data$cont_var3 + 
    #fixed_coefficients[25] * data$Z * data$cont_var4 + 
    #fixed_coefficients[26] * data$Z * data$cont_var5 + 
    fixed_coefficients[25] * 1 * data$bin_var1 + 
    fixed_coefficients[26] * 1 * data$bin_var2 + 
    fixed_coefficients[27] * 1 * data$bin_var3 + 
    1 * 1 * data$cont_var1^2
  #data[, 1:20]
  #fixed_coefficients[30] * data$Z * data$bin_var4 + 
  #fixed_coefficients[31] * data$Z * data$bin_var5 +
  
  
  
  data$outcome0 = as.matrix(data[, 1:20]) %*% fixed_coefficients[1:20] + 
    fixed_coefficients[21] * 0 + 
    fixed_coefficients[22] * 0 * data$cont_var1 + 
    fixed_coefficients[23] * 0 * data$cont_var2 + 
    fixed_coefficients[24] * 0 * data$cont_var3 + 
    #fixed_coefficients[25] * data$Z * data$cont_var4 + 
    #fixed_coefficients[26] * data$Z * data$cont_var5 + 
    fixed_coefficients[25] * 0 * data$bin_var1 + 
    fixed_coefficients[26] * 0 * data$bin_var2 + 
    fixed_coefficients[27] * 0 * data$bin_var3
  #fixed_coefficients[30] * data$Z * data$bin_var4 + 
  #fixed_coefficients[31] * data$Z * data$bin_var5 +
  
  
  if(binary){
    data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  selection_coefficients <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  pi = invlogit(as.matrix(data[,1:20]) %*% selection_coefficients - 1.95)
  
  pi = pi# / sum(pi) * selected
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  #ID = sample(1:n_samples, selected, prob = pi)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}


simulate_main = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 7
  n_binary <- 3
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  correlation = matrix(
    c(1, 0.2, 0, 0, 0, 0, 0,
      0.2, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0.9, 0, 0, 0,
      0, 0, 0.9, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1),
    nrow = n_continuous, ncol = n_continuous
  )
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = correlation)
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  data$Z = rbinom(n_samples, 1, prop)
  
  data$outcome1 = 
    2 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1 + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2
  
  data$outcome0 = 
    0 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1 + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2
  
  
  if(binary){
    #data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    #data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  pi = invlogit(-4.169 + 1 * data$cont_var1 + 2 * data$cont_var2 + -1 * data$cont_var3 + 1 * data$cont_var4 + 2 * data$bin_var1)
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}



simulate_standard_nonlinear_selection = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 7
  n_binary <- 3
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  correlation = matrix(
    c(1, 0.2, 0, 0, 0, 0, 0,
      0.2, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0.9, 0, 0, 0,
      0, 0, 0.9, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1),
    nrow = n_continuous, ncol = n_continuous
  )
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = correlation)
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  data$Z = rbinom(n_samples, 1, prop)
  
  data$outcome1 = 
    2 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1 + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2 + #log(0.5 + abs(data$cont_var2 + data$cont_var4)) +
    1 * (2 + 10 * data$cont_var1 - 5 * data$cont_var2 -4 * data$cont_var3 - 4 * data$cont_var4 + 8 * data$bin_var1) 
  #3 * data$cont_var1 * data$cont_var4 + 2 * data$cont_var3^2 - 1 * data$cont_var2 * data$bin_var1 +
  # 0.3 * data$cont_var4^3 * data$bin_var1) #+ invlogit(data$cont_var1 * data$cont_var2))
  
  data$outcome0 = 
    2 + 2 * data$cont_var1 + -1 * data$cont_var2 + 3 * data$cont_var3 + 2 * data$cont_var4 + 3 * data$bin_var1  + 
    2 * data$cont_var1 * data$cont_var3 * data$bin_var1 - 1 * data$cont_var2 * data$cont_var4 + 0.5 * data$cont_var1^2
  
  if(binary){
    #data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    #data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  pi = with(data,
            invlogit(-3.0115 + -2 * cont_var1^2 * bin_var1 + 2 * cont_var2 * data$cont_var3 + cont_var4^2)   
  )
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}

#a = sapply(1:1000, function(i) sum(simulate_standard_nonlinear_selection(n_samples = 5000, seed = i)$data$selected))
#mean(a)


simulate_overlapping_simplest = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
  # Set the seed for reproducibility
  if(!is.null(seed))
    set.seed(seed)
  
  # Number of samples
  
  
  # Number of continuous and binary variables
  n_continuous <- 7
  n_binary <- 3
  
  # Generate continuous variables (normally distributed)
  continuous_mean = rep(0, n_continuous)
  continuous_vars <- rmvnorm(n_samples, mean = continuous_mean, sigma = diag(n_continuous))
  
  
  
  # Generate binary variables (Bernoulli distributed)
  binary_mean = rep(0.5, n_binary)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  #data$Z = rbinom(n_samples, 1, 0.5)
  data$Z = rbinom(n_samples, 1, prop)
  
  
  # outcome model
  fixed_coefficients <- c(1.0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0,
                          1.0, 2) # Coefficients for Z and interactions
  
  
  data$outcome1 = as.matrix(data[, 1:10]) %*% fixed_coefficients[1:10] + 
    fixed_coefficients[11] * 1 + 
    fixed_coefficients[12] * 1 * data$cont_var1
  
  data$outcome0 = as.matrix(data[, 1:10]) %*% fixed_coefficients[1:10] 
  
  if(binary){
    data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  selection_coefficients <- c(-1, 0, 0, 0, 0, 0, 0, 
                              0, 0, 0)
  
  pi = invlogit(as.matrix(data[,1:10]) %*% selection_coefficients - 1.96)
  
  pi = pi# / sum(pi) * selected
  
  ID = which(rbinom(length(pi), 1, pi) == 1)
  #ID = sample(1:n_samples, selected, prob = pi)
  data$selected = F
  data$selected[ID] = T
  data$ps = pi
  
  #lm(outcome~. , data = data)
  EHR = as_tibble(data)
  EHR$Z = NULL
  EHR$outcome = NULL
  EHR$selected = NULL
  EHR$outcome0 = NULL
  EHR$outcome1 = NULL
  
  trials = data[ID,]
  trials$selected = NULL
  trials$outcome = ifelse(trials$Z == 1, trials$outcome1, trials$outcome0)
  trials$outcome0 = NULL
  trials$outcome1 = NULL
  
  
  return(list(data = data, trials = trials, EHR = EHR))
}

#a = sapply(1:1000, function(i) sum(simulate_overlapping_simplest(n_samples = 5000, seed = i)$data$selected))
#mean(a)

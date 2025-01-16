library(mvtnorm) 
library(arm)
library(dplyr)

simulate_simplest = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
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





simulate_simplest_200 = function(n_samples = 5000, seed = NULL, binary = FALSE, prop = 0.5){
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
  
  pi = invlogit(as.matrix(data[,1:10]) %*% selection_coefficients - 3.618)
  
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


# summary(foreach(seed = 1:1000, .combine = c) %dopar% {
#   data = simulate_simplest_200(seed = seed, binary = binary, prop = prop)
#   set.seed(seed)
#   
#   dim(data$trials)[1]
# })


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
    2 * data$cont_var1 + -1.5 * data$cont_var2 + 2.5 * data$cont_var3 + -2 * data$cont_var4 + 1 * data$bin_var1 + 
    1 * (2 + 3 * data$cont_var1 + 2 * data$cont_var2 + 4 * data$cont_var3 - 4 * data$cont_var4 + 3 * data$bin_var1) 
  
  data$outcome0 = 
    2 * data$cont_var1 + -1.5 * data$cont_var2 + 2.5 * data$cont_var3 + -2 * data$cont_var4 + 1 * data$bin_var1
  
  if(binary){
    data$outcome1 = rbinom(n_samples, 1, invlogit(data$outcome1))
    data$outcome0 = rbinom(n_samples, 1, invlogit(data$outcome0))
  }else{
    data$outcome1 =  data$outcome1 + rnorm(n_samples)
    data$outcome0 =  data$outcome0 + rnorm(n_samples)
  }
  
  # selection model
  #pi = invlogit(-1.977 + 1 * data$cont_var1 + 0.5 * data$cont_var2 + -1* data$cont_var3 + 1.5 * data$cont_var4 -0.5 * data$bin_var1)
  #pi = invlogit(-1.687 + 1 * data$cont_var1 - 1.5 * data$cont_var2 + -1* data$cont_var3 + 1.5 * data$cont_var4 - 2 * data$bin_var1)
  #summary(pi)
  pi = invlogit(-1.96 + 1 * data$cont_var1)
  
  #ID = which(pi >= runif(length(pi)))
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
# summary(foreach(seed = 1:1000, .combine = c) %dopar% {
#   data = simulate_standard_linear(seed = seed, binary = binary, prop = prop)
#   set.seed(seed)
# 
#   dim(data$trials)[1]
# })


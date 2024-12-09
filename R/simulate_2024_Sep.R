# Load necessary library
library(mvtnorm)
library(rstanarm)
library(tidyverse)

simulate_overlapping_both = function(selected = 300, n_samples = 5000, seed = NULL){
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
  data$Z = rbinom(n_samples, 1, 0.5)
  
  
  # outcome model
  fixed_coefficients <- c(2.0, -1.5, 1.2, 0, 0, 0, 0, 0, 0, 0,
                          2.8, -0.5, 3.7, 0, 0, 0, 0, 0, 0, 0,
                          3.0, 2.5, 2.9, -3.4, 2.6, -3.3, 2.9) # Coefficients for Z and interactions
  
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
    #data[, 1:20]
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples)
  
  
  data$outcome0 = as.matrix(data[, 1:20]) %*% fixed_coefficients[1:20] + 
    fixed_coefficients[21] * 0 + 
    fixed_coefficients[22] * 0 * data$cont_var1 + 
    fixed_coefficients[23] * 0 * data$cont_var2 + 
    fixed_coefficients[24] * 0 * data$cont_var3 + 
    #fixed_coefficients[25] * data$Z * data$cont_var4 + 
    #fixed_coefficients[26] * data$Z * data$cont_var5 + 
    fixed_coefficients[25] * 0 * data$bin_var1 + 
    fixed_coefficients[26] * 0 * data$bin_var2 + 
    fixed_coefficients[27] * 0 * data$bin_var3 + 
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples)
  
  
  #data$outcome = ifelse(data$Z == 1, data$outcome1, data$outcome0)
  #summary(lm(as.formula(paste("outcome~", paste(paste("cont_var", 1:10, sep = "", collapse = " + "), paste("bin_var", 1:10, sep = "", collapse = " + "), "Z", sep = "+"))), data = data ))
  
  
  
  
  
  
  
  
  # selection model
  selection_coefficients <- c(-4.8, 5.2, -3.4, 0, 0, 0, 0, 0, 0, 0, 
                              6, -4.1, 6.2, 0, 0, 0, 0, 0, 0, 0)
  
  pi = invlogit(as.matrix(data[,1:20]) %*% selection_coefficients - 3.5)
  
  pi = pi / sum(pi) * selected
  
  ID = sample(1:n_samples, selected, prob = pi)
  data$selected = F
  data$selected[ID] = T
  
  
  
  
  
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


simulate_high_outcome_low_selection = function(selected = 300, n_samples = 5000, seed = NULL){
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
  binary_mean =  rep(0.5, n_binary) #rbeta(n_binary, 5, 5)
  #print(binary_mean)
  binary_vars <- sapply(1:n_binary, function(i) rbinom(n_samples, size = 1, prob = binary_mean[i]))
  
  # Combine continuous and binary variables into one data frame
  data <- data.frame(continuous_vars, binary_vars)
  
  # Rename the columns
  colnames(data) <- c(paste0("cont_var", 1:n_continuous), paste0("bin_var", 1:n_binary))
  
  #data$Z = rbinom(n_samples, 1, 0.5)
  data$Z = rbinom(n_samples, 1, 0.5)
  
  
  # outcome model
  fixed_coefficients <- c(1.0, -1.5, 1.2, 0, 0, 0, 0, 0, 0, 0,
                          1.8, -2.5, 3.7, 0, 0, 0, 0, 0, 0, 0,
                          3.0, 2.5, 2.9, -3.4, 2.6, -3.3, 2.9) # Coefficients for Z and interactions
  
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
    #data[, 1:20]
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples)
  
  
  data$outcome0 = as.matrix(data[, 1:20]) %*% fixed_coefficients[1:20] + 
    fixed_coefficients[21] * 0 + 
    fixed_coefficients[22] * 0 * data$cont_var1 + 
    fixed_coefficients[23] * 0 * data$cont_var2 + 
    fixed_coefficients[24] * 0 * data$cont_var3 + 
    #fixed_coefficients[25] * data$Z * data$cont_var4 + 
    #fixed_coefficients[26] * data$Z * data$cont_var5 + 
    fixed_coefficients[25] * 0 * data$bin_var1 + 
    fixed_coefficients[26] * 0 * data$bin_var2 + 
    fixed_coefficients[27] * 0 * data$bin_var3 + 
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples)
  
  
  #data$outcome = ifelse(data$Z == 1, data$outcome1, data$outcome0)
  #summary(lm(as.formula(paste("outcome~", paste(paste("cont_var", 1:10, sep = "", collapse = " + "), paste("bin_var", 1:10, sep = "", collapse = " + "), "Z", sep = "+"))), data = data ))
  
  
  
  
  
  
  
  
  # selection model
  selection_coefficients <- c(-3.3, 5.6, 0, 0, 0, 0, 0, 0, 0, 0, 
                              4.5, 5, 0, 0, 0, 0, 0, 0, 0, 0)
  
  pi = invlogit(as.matrix(data[,1:20]) %*% selection_coefficients - 20)
  
  pi = pi / sum(pi) * selected
  
  ID = sample(1:n_samples, selected, prob = pi)
  data$selected = F
  data$selected[ID] = T
  
  
  
  
  
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


simulate_low_outcome_high_selection = function(selected = 300, n_samples = 5000, seed = NULL){
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
  data$Z = rbinom(n_samples, 1, 0.5)
  
  
  # outcome model
  fixed_coefficients <- c(5.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          -4.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          3.0, 4.5, 0, 0, 7.6, 0, 0) # Coefficients for Z and interactions
  
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
    #data[, 1:20]
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples) + 0
  
  
  data$outcome0 = as.matrix(data[, 1:20]) %*% fixed_coefficients[1:20] + 
    fixed_coefficients[21] * 0 + 
    fixed_coefficients[22] * 0 * data$cont_var1 + 
    fixed_coefficients[23] * 0 * data$cont_var2 + 
    fixed_coefficients[24] * 0 * data$cont_var3 + 
    #fixed_coefficients[25] * data$Z * data$cont_var4 + 
    #fixed_coefficients[26] * data$Z * data$cont_var5 + 
    fixed_coefficients[25] * 0 * data$bin_var1 + 
    fixed_coefficients[26] * 0 * data$bin_var2 + 
    fixed_coefficients[27] * 0 * data$bin_var3 + 
    #fixed_coefficients[30] * data$Z * data$bin_var4 + 
    #fixed_coefficients[31] * data$Z * data$bin_var5 +
    rnorm(n_samples)
  
  
  #data$outcome = ifelse(data$Z == 1, data$outcome1, data$outcome0)
  #summary(lm(as.formula(paste("outcome~", paste(paste("cont_var", 1:10, sep = "", collapse = " + "), paste("bin_var", 1:10, sep = "", collapse = " + "), "Z", sep = "+"))), data = data ))
  
  
  
  
  
  
  
  
  # selection model
  selection_coefficients <- c(2, 5.2, -4.4, 0, 0, 0, 0, 0, 0, 0, 
                              2, -4.1, 6.2, 0, 0, 0, 0, 0, 0, 0)
  
  pi = invlogit(as.matrix(data[,1:20]) %*% selection_coefficients - 1.5)
  #print(summary(pi))
  
  pi = pi / sum(pi) * selected
  
  ID = sample(1:n_samples, selected, prob = pi)
  data$selected = F
  data$selected[ID] = T
  
  
  
  
  
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








# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# summary(lm(outcome1 - outcome0~ 0 + ., data = data[, -24]))
# summary(glm(selected~ 1 + .*Z, data = data[,-c(22, 23)]))
# 
# cor(data$cont_var3, data$selected)
# 
# 
# 
# summary(sapply(1:100, function(i){
#   sim = simulate_low_outcome_high_selection(selected = 300, n_samples = 5000, seed = i)
#   
#   data = sim$data
#   mean(data$outcome0)
#   mean(data$outcome1)
#   mean(data$outcome0[data$selected])
#   mean(data$outcome1[data$selected])
#   mean(sim$trials$outcome[sim$trials$Z == 1]) - mean(sim$trials$outcome[sim$trials$Z == 0]) - (mean(data$outcome1) - mean(data$outcome0))
#   #mean(data$outcome1[data$selected]) - mean(data$outcome0[data$selected]) - (mean(data$outcome1) - mean(data$outcome0))
#   #mean(data$outcome1[data$selected]) - mean(data$outcome1)
# }))
# 
# 
# data %>% 
#   mutate(selected = factor(selected, levels = c(F, T))) %>% 
#   ggplot() + geom_point(aes(x = cont_var3, y = outcome1 - outcome0, color = selected))
# 
# 

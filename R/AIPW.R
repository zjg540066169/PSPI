apiw = function(Y, X, Z, ps, X_pop, ps_pop, normalize = F){
  outcome_model_treat <- SuperLearner(Y = Y[Z == 1],
                                      X = X[Z == 1, ],
                                      SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.nnet", "SL.rpart"))
  outcome_model_control <- SuperLearner(Y = Y[Z == 0],
                                      X = X[Z == 0, ],
                                      SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.nnet", "SL.rpart"))
  pred_treat = predict(outcome_model_treat, newdata = X_pop)$pred
  pred_control = predict(outcome_model_control, newdata = X_pop)$pred
  if(normalize == FALSE){
    weights1 <- 1 / (ps * mean(Z)) / length(pred_treat)
    weights0 <- 1 / (ps * (1 - mean(Z))) / length(pred_control)
  }else{
    weights1 <- 1 / (ps * mean(Z))
    weights0 <- 1 / (ps * (1 - mean(Z)))
    weights1 = weights1 / sum(weights1)
    weights0 = weights0 / sum(weights0)
  }
  
  
  outcome1 = mean(pred_treat) + sum(weights1[Z == 1] * (Y[Z == 1] -  predict(outcome_model_treat, newdata = X)$pred[Z == 1])) 
  outcome0 = mean(pred_control) + sum(weights0[Z == 0] * (Y[Z == 0] -  predict(outcome_model_control, newdata = X)$pred[Z == 0]))
  return(c(
    outcome1 = outcome1,
    outcome0 = outcome0,
    PATE = outcome1 - outcome0)
  )
}





apiw_bart = function(Y, X, Z, ps, X_pop, ps_pop, normalize = F, binary = F, nburn = 5000L, npost = 5000L){
  if(binary == F){
    invisible(capture.output(suppressMessages(outcome_model_treat <- BART::wbart(y.train = Y[Z == 1],
                                                                                 x.train = X[Z == 1, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(outcome_model_control <- BART::wbart(y.train = Y[Z == 0],
                                                                                   x.train = X[Z == 0, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- predict(outcome_model_treat, X_pop))))
    invisible(capture.output(suppressMessages(pred_control <- predict(outcome_model_control, X_pop))))
  }else{
    invisible(capture.output(suppressMessages(outcome_model_treat <- BART::pbart(y.train = Y[Z == 1],
                                                                                 x.train = X[Z == 1, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(outcome_model_control <- BART::pbart(y.train = Y[Z == 0],
                                                                                   x.train = X[Z == 0, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- predict(outcome_model_treat, X_pop)$prob.test)))
    invisible(capture.output(suppressMessages(pred_control <- predict(outcome_model_control, X_pop)$prob.test)))
  }
  
  
  
  
  
  if(normalize == FALSE){
    weights1 <- 1 / (ps * mean(Z)) / length(ps_pop)
    weights0 <- 1 / (ps * (1 - mean(Z))) / length(ps_pop)
  }else{
    weights1 <- 1 / (ps * mean(Z))
    weights0 <- 1 / (ps * (1 - mean(Z)))
    weights1 = weights1 / sum(weights1)
    weights0 = weights0 / sum(weights0)
  }
  if(binary == F){
    invisible(capture.output(suppressMessages(pred_treat_group1 <- predict(outcome_model_treat, X))))
    invisible(capture.output(suppressMessages(pred_treat_group0 <- predict(outcome_model_control, X))))
  }else{
    invisible(capture.output(suppressMessages(pred_treat_group1 <- predict(outcome_model_treat, X)$prob.test)))
    invisible(capture.output(suppressMessages(pred_treat_group0 <- predict(outcome_model_control, X)$prob.test)))
  }
  outcome1 = sapply(1:npost, function(i){
    outcome1 <- mean(pred_treat[i,]) + sum(weights1[Z == 1] * (Y[Z == 1] -  pred_treat_group1[i, Z == 1]))
    return(outcome1)
  })
  
  outcome0 = sapply(1:npost, function(i){
    outcome0 <- mean(pred_control[i,]) + sum(weights0[Z == 0] * (Y[Z == 0] -  pred_treat_group0[i, Z == 0]))
    return(outcome0)
  })
  
  
  
  
  # 
  # pred_treat = colMeans(pred_treat)
  # pred_control = colMeans(pred_control)
  # pred_treat_group1 = colMeans(pred_treat_group1)
  # pred_treat_group0 = colMeans(pred_treat_group0)
  # 
  # 
  # boot_outcome1_result = c()
  # boot_outcome0_result = c()
  # 
  # 
  # 
  # 
  # #boot = foreach(b = 1:bootstrap_n, .combine = rbind) %dopar%{
  # for (b in 1:bootstrap_n) {
  #   boot_pop_indices = sample(1:length(ps_pop), size = length(ps_pop), replace = T)
  #   boot_mcmc_indices = sample(1:npost, size = npost, replace = T)
  #   boot_Z_1 = sample(which(Z == 1), size = sum(Z == 1), replace = T)
  #   boot_Z_0 = sample(which(Z == 0), size = sum(Z == 0), replace = T)
  #   boot_tri_indices = c(boot_Z_1, boot_Z_0)
  #   
  #   boot_ps = ps[boot_tri_indices]
  #   boot_ps_pop = ps_pop[boot_pop_indices]
  #   boot_Z = Z[boot_tri_indices]
  #   boot_Y = Y[boot_tri_indices]
  #   
  #   if(normalize == FALSE){
  #     boot_weights1 <- 1 / (boot_ps * mean(boot_Z)) / length(boot_ps_pop)
  #     boot_weights0 <- 1 / (boot_ps * (1 - mean(boot_Z))) / length(boot_ps_pop)
  #   }else{
  #     boot_weights1 <- 1 / (boot_ps * mean(boot_Z))
  #     boot_weights0 <- 1 / (boot_ps * (1 - mean(boot_Z)))
  #     boot_weights1 = boot_weights1 / sum(boot_weights1)
  #     boot_weights0 = boot_weights0 / sum(boot_weights0)
  #   }
  #   
  #   boot_outcome1 = mean(pred_treat[boot_pop_indices]) + sum(boot_weights1[1:length(boot_Z_1)] * (boot_Y[1:length(boot_Z_1)] -  pred_treat_group1[boot_Z_1]))
  #   boot_outcome0 = mean(pred_control[boot_pop_indices]) + sum(boot_weights0[(length(boot_Z_1)+1):length(boot_tri_indices)] * (boot_Y[(length(boot_Z_1)+1):length(boot_tri_indices)] -  pred_treat_group0[boot_Z_0]))
  #   # boot_outcome1 = mean(sapply(boot_mcmc_indices, function(i){
  #   #   outcome1 <- mean(pred_treat[i, boot_pop_indices]) + sum(boot_weights1[1:length(boot_Z_1)] * (boot_Y[1:length(boot_Z_1)] -  pred_treat_group1[i, boot_Z_1]))
  #   #   return(outcome1)
  #   # }))
  #   # 
  #   # boot_outcome0 = mean(sapply(boot_mcmc_indices, function(i){
  #   #   outcome0 <- mean(pred_control[i, boot_pop_indices]) + sum(boot_weights0[(length(boot_Z_1)+1):length(boot_tri_indices)] * (boot_Y[(length(boot_Z_1)+1):length(boot_tri_indices)] -  pred_treat_group0[i, boot_Z_0]))
  #   #   return(outcome0)
  #   # }))
  #  # return(c(boot_outcome0, boot_outcome1))
  #   boot_outcome1_result = c(boot_outcome1_result, boot_outcome1)
  #   boot_outcome0_result = c(boot_outcome0_result, boot_outcome0)
  #   
  # }
  # 
  # 
  # 
  # 
  # 
  # 
  # c(quantile(boot_outcome1_result, 0.025), quantile(boot_outcome1_result, 0.975))
  # c(quantile(outcome1, 0.025), quantile(outcome1, 0.975))
  # 
  # c(quantile(boot_outcome0_result, 0.025), quantile(boot_outcome0_result, 0.975))
  # c(quantile(outcome0, 0.025), quantile(outcome0, 0.975))
  # 
  # c(quantile(boot_outcome1_result - boot_outcome0_result, 0.025), quantile(boot_outcome1_result - boot_outcome0_result, 0.975))
  # c(quantile(outcome1 - outcome0, 0.025), quantile(outcome1 - outcome0, 0.975))
  
  
  return(list(
    outcome1 = outcome1,
    outcome0 = outcome0,
    ATE = outcome1 - outcome0)
  )
}


# 
# aipw_function <- function(data, indices) {
#   d <- data[indices, ]  # Resample data
#   # Fit outcome models for treatment and control groups
#   outcome_model_treat <- SuperLearner(Y = d[d$Z == 1, c(23)],
#                         X = d[d$Z == 1, c(1:20)],
#                         SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.nnet", "SL.rpart"))
#   
#   
#   outcome_model_treat <- lm(outcome ~ ., data = d[d$Z == 1, c(1:20, 23)])
#   outcome_model_control <- lm(outcome ~ ., data = d[d$Z == 0, c(1:20, 23)])
#   d$pred_treat <- predict(outcome_model_treat, newdata = d[, c(1:20, 23)])
#   d$pred_control <- predict(outcome_model_control, newdata = d[, c(1:20, 23)])
#   d$weights <- 1 / d$ps
#   d$aipw <- with(d, 
#                  Z * outcome / ps - 
#                    (1 - Z) * outcome / (1 - ps) + 
#                    (1 - Z) * pred_control + 
#                    Z * pred_treat - 
#                    pred_control
#   )
#   # Return mean AIPW estimate (ATE)
#   return(mean(d$aipw))
# }
# 
# 
# aipw_function_treat_mean <- function(data, indices) {
#   d <- data[indices, ]
#   
#   # Fit outcome models
#   outcome_model_treat <- lm(outcome ~ ., data = d[d$Z == 1, c(1:20, 23)])
#   outcome_model_control <- lm(outcome ~ ., data = d[d$Z == 0, c(1:20, 23)])
#   
#   # Predict potential outcomes for everyone
#   d$pred_treat <- predict(outcome_model_treat, newdata = d[, c(1:20, 23)])
#   d$pred_control <- predict(outcome_model_control, newdata = d[, c(1:20, 23)])
#   
#   # Propensity score-based weights
#   d$weights <- 1 / d$ps
#   
#   # Estimate E[Y(1)] using AIPW
#   d$aipw_treat <- with(d,
#                        pred_treat + (Z * (outcome - pred_treat)) / ps
#   )
#   
#   return(mean(d$aipw_treat))
# }
# 
# aipw_function_control_mean <- function(data, indices) {
#   d <- data[indices, ]
#   
#   # Fit outcome models
#   outcome_model_treat <- lm(outcome ~ ., data = d[d$Z == 1, c(1:20, 23)])
#   outcome_model_control <- lm(outcome ~ ., data = d[d$Z == 0, c(1:20, 23)])
#   
#   # Predict potential outcomes for everyone
#   d$pred_treat <- predict(outcome_model_treat, newdata = d[, c(1:20, 23)])
#   d$pred_control <- predict(outcome_model_control, newdata = d[, c(1:20, 23)])
#   
#   # Propensity score-based weights
#   d$weights <- 1 / d$ps
#   
#   # Estimate E[Y(0)] using AIPW
#   d$aipw_control <- with(d,
#                          pred_control + ((1 - Z) * (outcome - pred_control)) / (1 - ps)
#   )
#   
#   return(mean(d$aipw_control))
# }

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





apiw_bart = function(Y, X, Z, ps, X_pop, ps_pop, s, normalize = F, binary = F, nburn = 5000L, npost = 5000L){
  if(binary == F){
    invisible(capture.output(suppressMessages(outcome_model_treat <- BART::wbart(y.train = Y[Z == 1],
                                                                                 x.train = X[Z == 1, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(outcome_model_control <- BART::wbart(y.train = Y[Z == 0],
                                                                                   x.train = X[Z == 0, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- (predict(outcome_model_treat, X_pop)))))
    invisible(capture.output(suppressMessages(pred_control <- (predict(outcome_model_control, X_pop)))))
  }else{
    invisible(capture.output(suppressMessages(outcome_model_treat <- BART::pbart(y.train = Y[Z == 1],
                                                                                 x.train = X[Z == 1, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(outcome_model_control <- BART::pbart(y.train = Y[Z == 0],
                                                                                   x.train = X[Z == 0, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- predict(outcome_model_treat, X_pop)$prob.test)))
    invisible(capture.output(suppressMessages(pred_control <- predict(outcome_model_control, X_pop)$prob.test)))
  }
  
  if(binary == F){
    pred_treat_group1 = pred_treat[,s]
    pred_treat_group0 = pred_control[,s]
  }else{
    invisible(capture.output(suppressMessages(pred_treat_group1 <- predict(outcome_model_treat, X)$prob.test)))
    invisible(capture.output(suppressMessages(pred_treat_group0 <- predict(outcome_model_control, X)$prob.test)))
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
  

  outcome1 = t(sapply(1:npost, function(i){
    outcome1 <- mean(pred_treat[i,]) + sum(weights1[Z == 1] * (Y[Z == 1] -  pred_treat_group1[i, Z == 1]))
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = weights1[Z == 1], outcome = Y[Z == 1] -  pred_treat_group1[i, Z == 1]))
    SE = unname(survey::SE(survey::svymean(~outcome, des)))
    return(c(outcome1, SE))
  }))
  
  outcome0 = t(sapply(1:npost, function(i){
    outcome0 <- mean(pred_control[i,]) + sum(weights0[Z == 0] * (Y[Z == 0] -  pred_treat_group0[i, Z == 0]))
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = weights0[Z == 0], outcome = Y[Z == 0] -  pred_treat_group0[i, Z == 0]))
    SE = unname(survey::SE(survey::svymean(~outcome, des)))
    return(c(outcome0, SE))
  }))
  
  se1 = sqrt(mean(outcome1[,2]^2) + var(outcome1[,1]) * (1 + 1 / dim(outcome1)[1]))
  se0 = sqrt(mean(outcome0[,2]^2) + var(outcome0[,1]) * (1 + 1 / dim(outcome0)[1]))
  
  se = sqrt(mean((sapply(1:npost, function(i){
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = ifelse(Z == 0, weights0, weights1), Z = Z, outcome = ifelse(Z == 0, Y -  pred_treat_group0[i, ],  Y -  pred_treat_group1[i, ]) ))

    ttest_result = survey::svyttest(outcome~I(Z == 1), des)
    diff_means <- ttest_result$estimate
    t_value <- ttest_result$statistic
    return(unname(abs(diff_means / t_value)))
  }))^2) + var(outcome1[,1] + outcome0[,1]) * (1 + 1 / dim(outcome1)[1]))
  
  return(list(
    outcome1 = mean(outcome1[,1]),
    outcome0 = mean(outcome0[,1]),
    ATE = mean(outcome1[,1] - outcome0[,1]),
    se = se,
    se1 = se1,
    se0 = se0)
  )
}






apiw_bart_joint = function(Y, X, Z, ps, X_pop, ps_pop, s, normalize = F, binary = F, nburn = 5000L, npost = 5000L){
  if(binary == F){
    invisible(capture.output(suppressMessages(outcome_model <- BART::wbart(y.train = Y,
                                                                                 x.train = cbind(X, Z), ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- (predict(outcome_model, cbind(X_pop, 1))))))
    invisible(capture.output(suppressMessages(pred_control <- (predict(outcome_model, cbind(X_pop, 0))))))
    
    pred_treat_group1 = pred_treat[,s]
    pred_treat_group0 = pred_control[,s]
  
  }else{
    invisible(capture.output(suppressMessages(outcome_model_treat <- BART::pbart(y.train = Y[Z == 1],
                                                                                 x.train = X[Z == 1, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(outcome_model_control <- BART::pbart(y.train = Y[Z == 0],
                                                                                   x.train = X[Z == 0, ], ndpost=npost, nskip=nburn))))
    invisible(capture.output(suppressMessages(pred_treat <- predict(outcome_model_treat, X_pop)$prob.test)))
    invisible(capture.output(suppressMessages(pred_control <- predict(outcome_model_control, X_pop)$prob.test)))
  }
  
  if(binary == F){
    
  }else{
    invisible(capture.output(suppressMessages(pred_treat_group1 <- predict(outcome_model_treat, X)$prob.test)))
    invisible(capture.output(suppressMessages(pred_treat_group0 <- predict(outcome_model_control, X)$prob.test)))
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
  
  
  outcome1 = t(sapply(1:npost, function(i){
    outcome1 <- mean(pred_treat[i,]) + sum(weights1[Z == 1] * (Y[Z == 1] -  pred_treat_group1[i, Z == 1]))
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = weights1[Z == 1], outcome = Y[Z == 1] -  pred_treat_group1[i, Z == 1]))
    SE = unname(survey::SE(survey::svymean(~outcome, des)))
    return(c(outcome1, SE))
  }))
  
  outcome0 = t(sapply(1:npost, function(i){
    outcome0 <- mean(pred_control[i,]) + sum(weights0[Z == 0] * (Y[Z == 0] -  pred_treat_group0[i, Z == 0]))
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = weights0[Z == 0], outcome = Y[Z == 0] -  pred_treat_group0[i, Z == 0]))
    SE = unname(survey::SE(survey::svymean(~outcome, des)))
    return(c(outcome0, SE))
  }))
  
  se1 = sqrt(mean(outcome1[,2]^2) + var(outcome1[,1]) * (1 + 1 / dim(outcome1)[1]))
  se0 = sqrt(mean(outcome0[,2]^2) + var(outcome0[,1]) * (1 + 1 / dim(outcome0)[1]))
  
  se = sqrt(mean((sapply(1:npost, function(i){
    des <- survey::svydesign(ids = ~1, weights = ~weight, data = data.frame(weight = ifelse(Z == 0, weights0, weights1), Z = Z, outcome = ifelse(Z == 0, Y -  pred_treat_group0[i, ],  Y -  pred_treat_group1[i, ]) ))
    
    ttest_result = survey::svyttest(outcome~I(Z == 1), des)
    diff_means <- ttest_result$estimate
    t_value <- ttest_result$statistic
    return(unname(abs(diff_means / t_value)))
  }))^2) + var(outcome1[,1] + outcome0[,1]) * (1 + 1 / dim(outcome1)[1]))
  
  return(list(
    outcome1 = mean(outcome1[,1]),
    outcome0 = mean(outcome0[,1]),
    ATE = mean(outcome1[,1] - outcome0[,1]),
    se = se,
    se1 = se1,
    se0 = se0)
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

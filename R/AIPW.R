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


apiw_bart = function(Y, X, Z, ps, X_pop, ps_pop, normalize = F){
  outcome_model_treat <- BART::wbart(y.train = Y[Z == 1],
                                      x.train = X[Z == 1, ], ndpost=5000L, nskip=5000L)
  outcome_model_control <- BART::wbart(y.train = Y[Z == 0],
                                        x.train = X[Z == 0, ], ndpost=5000L, nskip=5000L)
  pred_treat = colMeans(predict(outcome_model_treat, X_pop))
  pred_control = colMeans(predict(outcome_model_control, X_pop))
  if(normalize == FALSE){
    weights1 <- 1 / (ps * mean(Z)) / length(pred_treat)
    weights0 <- 1 / (ps * (1 - mean(Z))) / length(pred_control)
  }else{
    weights1 <- 1 / (ps * mean(Z))
    weights0 <- 1 / (ps * (1 - mean(Z)))
    weights1 = weights1 / sum(weights1)
    weights0 = weights0 / sum(weights0)
  }
  
  
  outcome1 = mean(pred_treat) + sum(weights1[Z == 1] * (Y[Z == 1] -  colMeans(predict(outcome_model_treat, X[Z == 1,]))))
  outcome0 = mean(pred_control) + sum(weights0[Z == 0] * (Y[Z == 0] -  colMeans(predict(outcome_model_control, X[Z == 0,]))))
  return(c(
    outcome1 = outcome1,
    outcome0 = outcome0,
    PATE = outcome1 - outcome0)
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

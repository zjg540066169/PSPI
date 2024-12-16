library(SuperLearner)
library(glmnet)
library(stats)

TMLE_generalizable <- function(Y, X, Z, ps, X_pop, ps_pop,
                               use_logit_link = FALSE) {
  Y_src <- Y
  X_src <- X
  Z_src <- Z
  ps_src <- ps # P(S=1|X) for source samples
  
  # Compute marginal treatment probabilities in the source population:
  pZ1 <- mean(Z_src)
  pZ0 <- 1 - pZ1
  
  # If using a logit link for bounding Y
  if(use_logit_link) {
    y_min <- min(Y_src) - 1e-10
    y_max <- max(Y_src) + 1e-10
    Y_rescaled <- (Y_src - y_min)/(y_max - y_min)
    Y_fit <- Y_rescaled
    family <- gaussian()
  } else {
    Y_fit <- Y_src
    family <- gaussian()
  }
  
  # Fit outcome model Q(A,X) using SuperLearner on the source data
  # Here, we include Z as a covariate. Another approach: Fit separate models for Z=1 and Z=0.
  Q_fit <- SuperLearner(Y = Y_fit,
                        X = data.frame(X_src, Z=Z_src),
                        family = family,
                        SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.nnet", "SL.rpart"))
  
  # Get initial predictions for Q(A,X)
  newdata_1 <- data.frame(X_src, Z=1)
  newdata_0 <- data.frame(X_src, Z=0)
  Q1_init <- predict(Q_fit, newdata=newdata_1)$pred
  Q0_init <- predict(Q_fit, newdata=newdata_0)$pred
  
  # Q(A,W) for the observed A
  QAW_init <- ifelse(Z_src == 1, Q1_init, Q0_init)
  
  # Construct clever covariates according to Schmid et al.
  # H1 = I(Z=1)/(P(Z=1)*P(S=1|X)), H0 = I(Z=0)/(P(Z=0)*P(S=1|X))
  H1 <- ifelse(Z_src == 1, 1/(pZ1 * ps_src), 0)
  H0 <- ifelse(Z_src == 0, -1/(pZ0 * ps_src), 0)
  
  # Combine into a single clever covariate H = H1 - H0
  # This is often how the fluctuation is done: one model with H as the covariate
  H <- H1 + H0
  
  if(use_logit_link) {
    # Targeting step: logistic fluctuation
    # Logit link: We find epsilon by solving a logistic regression with offset.
    # That is: logit(QAW_init) + epsilon*H ~ Y
    # We'll do a one-step approximation:
    
    # Ensure QAW_init is in (0,1)
    QAW_init_clamped <- pmin(pmax(QAW_init, 1e-10), 1 - 1e-10)
    # Offset in logistic regression:
    offset_vals <- qlogis(QAW_init_clamped)
    
    # Fit logistic regression: Y_rescaled ~ H with offset
    fit_epsilon <- glm(Y_rescaled ~ 0 + H, offset=offset_vals, family=binomial())
    epsilon <- coef(fit_epsilon)
    
    # Update Q1 and Q0
    Q1_star <- invlogit(qlogis(pmin(pmax(Q1_init,1e-10),1-1e-10)) + epsilon*H1)
    Q0_star <- invlogit(qlogis(pmin(pmax(Q0_init,1e-10),1-1e-10)) + epsilon*H0)
    
    # Rescale back to original Y scale
    Q1_update <- Q1_star*(y_max - y_min) + y_min
    Q0_update <- Q0_star*(y_max - y_min) + y_min
  } else {
    # Identity link targeting:
    # Solve for epsilon in a linear model: (Y - QAW_init) ~ H
    # This gives epsilon = cov(H, Y - QAW_init)/var(H)
    fit_epsilon <- lm(I(Y_fit - QAW_init) ~ 0 + H)
    epsilon <- coef(fit_epsilon)
    
    # Update Q
    Q1_update <- Q1_init + epsilon*H1
    Q0_update <- Q0_init + epsilon*H0
  }
  
  # The TMLE estimate in the source population
  tmle_est_source <- mean(Q1_update - Q0_update)
  
  H1_pop = 1/(pZ1 * ps_pop)
  H0_pop = -1/(pZ0 * ps_pop)
  
  newdata_1_pop <- data.frame(X_pop, Z=1)
  newdata_0_pop <- data.frame(X_pop, Z=0)
  Q1_init_pop <- predict(Q_fit, newdata=newdata_1_pop)$pred
  Q0_init_pop <- predict(Q_fit, newdata=newdata_0_pop)$pred
  
  if(use_logit_link) {
    # Update Q1 and Q0
    Q1_star_pop <- invlogit(qlogis(pmin(pmax(Q1_init_pop,1e-10),1-1e-10)) + epsilon*H1_pop)
    Q0_star_pop <- invlogit(qlogis(pmin(pmax(Q0_init_pop,1e-10),1-1e-10)) + epsilon*H0_pop)
    
    # Rescale back to original Y scale
    Q1_update_pop <- Q1_star_pop*(y_max - y_min) + y_min
    Q0_update_pop <- Q0_star_pop*(y_max - y_min) + y_min
  } else {

    # Update Q
    Q1_update_pop <- Q1_init_pop + epsilon*H1_pop
    Q0_update_pop <- Q0_init_pop + epsilon*H0_pop
  }
  tmle_est_pop = mean(Q1_update_pop - Q0_update_pop)
  
  Psi1 <- mean(Q1_update_pop)
  Psi0 <- mean(Q0_update_pop)
  
  # Compute Influence Curves for Each Mean
  # IC1 is based on Psi1 = mean(Q1_update_pop)
  # IC0 is based on Psi0 = mean(Q0_update_pop)
  
  # Influence Curve for Psi1 (Mean under treatment)
  IC1 <- H1 * (Y_src - Q1_update) + (Q1_update - mean(Q1_update))
  
  # Influence Curve for Psi0 (Mean under control)
  IC0 <- H0 * (Y_src - Q0_update) + (Q0_update - mean(Q0_update))
  
  # Influence Curve for ATE (Psi = Psi1 - Psi0)
  IC_diff <- IC1 - IC0
  
  # Variance Estimation using Sample Variance of the Influence Curve
  # As per Schmid et al. (2022), Var(TMLE) = Var(IC) / n
  n <- length(Y_src)
  
  var_Psi1 <- var(IC1) / n
  var_Psi0 <- var(IC0) / n
  var_diff <- var(IC_diff) / n
  
  # Standard Errors
  se_Psi1 <- sqrt(var_Psi1)
  se_Psi0 <- sqrt(var_Psi0)
  se_diff <- sqrt(var_diff)
  
  return(c(
    tmle_est_pop = tmle_est_pop,
    outcome1 = Psi1,
    outcome0 = Psi0,
    
    se_diff = sqrt(var_diff),
    se_Psi1 = se_Psi1,
    se_Psi0 = se_Psi0,
    epsilon = epsilon
  ))
}

# X_pop = data$data[,1:20]
# ps_pop = data$data$ps
# use_logit_link = FALSE
# bootstrap <- boot(data = trials, statistic = TMLE_function, R = 100, parallel = "multicore", ncpus = 5)
# 
# TMLE_function <- function(data, indices) {
#   use_logit_link = TRUE
#   data = data[indices, ]
#   Y_src <- data$outcome
#   X_src <- data[,1:20]
#   Z_src <- data$Z
#   ps_src <- data$ps # P(S=1|X) for source samples
#   
#   # Compute marginal treatment probabilities in the source population:
#   pZ1 <- mean(Z_src)
#   pZ0 <- 1 - pZ1
#   
#   # If using a logit link for bounding Y
#   if(use_logit_link) {
#     y_min <- min(Y_src) - 1e-10
#     y_max <- max(Y_src) + 1e-10
#     Y_rescaled <- (Y_src - y_min)/(y_max - y_min)
#     Y_fit <- Y_rescaled
#     family <- gaussian()
#   } else {
#     Y_fit <- Y_src
#     family <- gaussian()
#   }
#   
#   # Fit outcome model Q(A,X) using SuperLearner on the source data
#   # Here, we include Z as a covariate. Another approach: Fit separate models for Z=1 and Z=0.
#   Q_fit <- SuperLearner(Y = Y_fit,
#                         X = data.frame(X_src, Z=Z_src),
#                         family = family,
#                         SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.nnet", "SL.rpart"))
#   
#   # Get initial predictions for Q(A,X)
#   newdata_1 <- data.frame(X_src, Z=1)
#   newdata_0 <- data.frame(X_src, Z=0)
#   Q1_init <- predict(Q_fit, newdata=newdata_1)$pred
#   Q0_init <- predict(Q_fit, newdata=newdata_0)$pred
#   
#   # Q(A,W) for the observed A
#   QAW_init <- ifelse(Z_src == 1, Q1_init, Q0_init)
#   
#   # Construct clever covariates according to Schmid et al.
#   # H1 = I(Z=1)/(P(Z=1)*P(S=1|X)), H0 = I(Z=0)/(P(Z=0)*P(S=1|X))
#   H1 <- ifelse(Z_src == 1, 1/(pZ1 * ps_src), 0)
#   H0 <- ifelse(Z_src == 0, 1/(pZ0 * ps_src), 0)
#   
#   # Combine into a single clever covariate H = H1 - H0
#   # This is often how the fluctuation is done: one model with H as the covariate
#   H <- H1 + H0
#   
#   if(use_logit_link) {
#     # Targeting step: logistic fluctuation
#     # Logit link: We find epsilon by solving a logistic regression with offset.
#     # That is: logit(QAW_init) + epsilon*H ~ Y
#     # We'll do a one-step approximation:
#     
#     # Ensure QAW_init is in (0,1)
#     QAW_init_clamped <- pmin(pmax(QAW_init, 1e-10), 1 - 1e-10)
#     # Offset in logistic regression:
#     offset_vals <- qlogis(QAW_init_clamped)
#     
#     # Fit logistic regression: Y_rescaled ~ H with offset
#     fit_epsilon <- glm(Y_rescaled ~ 0 + H, offset=offset_vals, family=binomial())
#     epsilon <- coef(fit_epsilon)
#     
#     # Update Q1 and Q0
#     Q1_star <- invlogit(qlogis(pmin(pmax(Q1_init,1e-10),1-1e-10)) + epsilon*H1)
#     Q0_star <- invlogit(qlogis(pmin(pmax(Q0_init,1e-10),1-1e-10)) + epsilon*H0)
#     
#     # Rescale back to original Y scale
#     Q1_update <- Q1_star*(y_max - y_min) + y_min
#     Q0_update <- Q0_star*(y_max - y_min) + y_min
#   } else {
#     # Identity link targeting:
#     # Solve for epsilon in a linear model: (Y - QAW_init) ~ H
#     # This gives epsilon = cov(H, Y - QAW_init)/var(H)
#     fit_epsilon <- lm(I(Y_fit - QAW_init) ~ 0 + H)
#     epsilon <- coef(fit_epsilon)
#     
#     # Update Q
#     Q1_update <- Q1_init + epsilon*H1
#     Q0_update <- Q0_init + epsilon*H0
#   }
#   
#   # The TMLE estimate in the source population
#   tmle_est_source <- mean(Q1_update - Q0_update)
#   
#   H1_pop = 1/(pZ1 * ps_pop)
#   H0_pop = 1/(pZ0 * ps_pop)
#   
#   newdata_1_pop <- data.frame(X_pop, Z=1)
#   newdata_0_pop <- data.frame(X_pop, Z=0)
#   Q1_init_pop <- predict(Q_fit, newdata=newdata_1_pop)$pred
#   Q0_init_pop <- predict(Q_fit, newdata=newdata_0_pop)$pred
#   
#   if(use_logit_link) {
#     # Update Q1 and Q0
#     Q1_star_pop <- invlogit(qlogis(pmin(pmax(Q1_init_pop,1e-10),1-1e-10)) + epsilon*H1_pop)
#     Q0_star_pop <- invlogit(qlogis(pmin(pmax(Q0_init_pop,1e-10),1-1e-10)) + epsilon*H0_pop)
#     
#     # Rescale back to original Y scale
#     Q1_update_pop <- Q1_star_pop*(y_max - y_min) + y_min
#     Q0_update_pop <- Q0_star_pop*(y_max - y_min) + y_min
#   } else {
#     
#     # Update Q
#     Q1_update_pop <- Q1_init_pop + epsilon*H1_pop
#     Q0_update_pop <- Q0_init_pop + epsilon*H0_pop
#   }
#   tmle_est_pop = mean(Q1_update_pop - Q0_update_pop)
#   
#   Psi1 <- mean(Q1_update_pop)
#   Psi0 <- mean(Q0_update_pop)
#   return(mean(Q1_update))
#   # Compute Influence Curves for Each Mean
#   # IC1 is based on Psi1 = mean(Q1_update_pop)
#   # IC0 is based on Psi0 = mean(Q0_update_pop)
#   
#   # Influence Curve for Psi1 (Mean under treatment)
#   IC1 <- H1 * (Y_src - Q1_update) + (Q1_update - Psi1)
#   
#   # Influence Curve for Psi0 (Mean under control)
#   IC0 <- H0 * (Y_src - Q0_update) + (Q0_update - Psi0)
#   
#   # Influence Curve for ATE (Psi = Psi1 - Psi0)
#   IC_diff <- IC1 - IC0
#   
#   # Variance Estimation using Sample Variance of the Influence Curve
#   # As per Schmid et al. (2022), Var(TMLE) = Var(IC) / n
#   n <- length(Y_src)
#   
#   var_Psi1 <- var(IC1) / n
#   var_Psi0 <- var(IC0) / n
#   var_diff <- var(IC_diff) / n
#   
#   # Standard Errors
#   se_Psi1 <- sqrt(var_Psi1)
#   se_Psi0 <- sqrt(var_Psi0)
#   se_diff <- sqrt(var_diff)
#   
#   return(list(
#     tmle_est_pop = tmle_est_pop,
#     outcome1 = Psi1,
#     outcome0 = Psi0,
#     
#     se_diff = sqrt(var_diff),
#     se_Psi1 = se_Psi1,
#     se_Psi0 = se_Psi0,
#     epsilon = epsilon
#   ))
# }

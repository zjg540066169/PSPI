library(dplyr)
source("./simulate_2024_Sep_simple.R")
library(rstan)
library(rstanarm)
library(RcppArmadillo)
library(Rcpp)
library(boot)
library(SuperLearner)
library(readr)
sourceCpp("src/MCMC_BART_Causal.cpp", args = c("-fopenmp"))


path = "./simulation_nov/"
nburn = 4
npost = 4

args = commandArgs(trailingOnly = T)
seed = as.integer(args[1])
scenario = as.integer(args[2])
cat("seed:", seed, "\n")
cat("scenario:", scenario, "\n")
#dir.create(paste0(path, seed))







aipw_function <- function(data, indices) {
  d <- data[indices, ]  # Resample data
  # Fit outcome models for treatment and control groups
  
  outcome_model_treat <- lm(outcome ~ ., data = d[d$Z == 1, c(1:20, 23)])
  outcome_model_control <- lm(outcome ~ ., data = d[d$Z == 0, c(1:20, 23)])
  d$pred_treat <- predict(outcome_model_treat, newdata = d[, c(1:20, 23)])
  d$pred_control <- predict(outcome_model_control, newdata = d[, c(1:20, 23)])
  d$weights <- 1 / d$ps
  d$aipw <- with(d, 
                 Z * outcome / ps - 
                   (1 - Z) * outcome / (1 - ps) + 
                   (1 - Z) * pred_control + 
                   Z * pred_treat - 
                   pred_control
  )
  # Return mean AIPW estimate (ATE)
  return(mean(d$aipw))
}






if(scenario == 1){
  data = simulate_overlapping_both(selected = 800, seed = seed)
  cat("scenario: overlapping_both\n")
  set.seed(seed)
}
if(scenario == 2){
  data = simulate_high_outcome_low_selection(selected = 800, seed = seed)
  cat("scenario: high_outcome_low_selection\n")
  set.seed(seed)
}
if(scenario == 3){
  data = simulate_low_outcome_high_selection(selected = 800, seed = seed)
  cat("scenario: low_outcome_high_selection\n")
  set.seed(seed)
}
if(scenario == 4){
  data = simulate_overlapping_both_quadratic(selected = 800, seed = seed)
  cat("scenario: overlapping_both_quadratic\n")
  set.seed(seed)
}




trials = data$trials
EHR = as.data.frame(data$EHR)

true_TE = mean(data$data$outcome1 - data$data$outcome0)

sample_mean = mean(trials$outcome[trials$Z == 1]) - mean(trials$outcome[trials$Z == 0])

se = sqrt(var(trials$outcome[trials$Z == 1]) / length(trials$outcome[trials$Z == 1]) + var(trials$outcome[trials$Z == 0]) / length(trials$outcome[trials$Z == 0]))




true = stan_lm(outcome ~ Z * (cont_var1 + cont_var2 + cont_var3 + bin_var1 + bin_var2 + bin_var3), data = trials, chains = 1, prior=NULL, iter = nburn + npost)
true_predict1 = posterior_predict(true, cbind(EHR, Z = 1))
true_predict0 = posterior_predict(true, cbind(EHR, Z = 0))
true_posterior = rowMeans(true_predict1 - true_predict0)

sel = c(rep(1, dim(trials)[1]), rep(0, dim(EHR)[1]))
co = rbind(trials[,1:20], EHR[,1:20])

bart = stan_glm(sel ~ . , family = binomial(link = "logit"), data = cbind(sel, co), chains = 1, prior=NULL, iter = nburn + npost)

trials$e_ps = colMeans(posterior_epred(bart, trials[,1:20]))
EHR$e_ps = colMeans(posterior_epred(bart, EHR[,1:20]))
data$data$e_ps = colMeans(posterior_epred(bart, EHR[,1:20]))
X_s = EHR$e_ps


des <- survey::svydesign(ids = ~1, weights = ~ 1 / ps, data = trials)
IPW = unname(survey::svyttest(outcome~Z, des)$estimate)
IPW_2_5 = survey::svyttest(outcome~Z, des)$conf.int[1]
IPW_97_5 = survey::svyttest(outcome~Z, des)$conf.int[2]


AIPW = aipw_function(trials, 1:dim(trials)[1])
bootstrap_results <- boot(data = trials, statistic = aipw_function, R = 1000)
AIPW_CI = c(boot.ci(bootstrap_results, type = "perc")$percent[4], boot.ci(bootstrap_results, type = "perc")$percent[5])

cbart = MCMC_BART_Causal(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 2, nburn, npost)
model_1 = MCMC_BART_Causal(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 3, nburn, npost)
model_2_4 = MCMC_BART_Causal(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 4, nburn, npost)
model_2_4_spline = MCMC_BART_Causal(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 5, nburn, npost)
BART_model = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z, trials$e_ps)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
bart_pure_TE = rowMeans(predict(BART_model, cbind(as.matrix(EHR[,1:20]), 1, EHR$e_ps)) - predict(BART_model, cbind(as.matrix(EHR[,1:20]), 0, EHR$e_ps)))
BART_model_no_pi = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
bart_pure_TE_no_pi = rowMeans(predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 1)) - predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 0)))
result = c(
  replicates = seed,
  true_TE = true_TE,
  sample_mean = sample_mean,
  sample_mean_2_5 = t.test(trials$outcome[trials$Z == 1], trials$outcome[trials$Z == 0])$conf.int[1], #sample_mean - 1.96 * se,
  sample_mean_97_5 = t.test(trials$outcome[trials$Z == 1], trials$outcome[trials$Z == 0])$conf.int[2], #sample_mean + 1.96 * se
  IPW = IPW,
  IPW_2_5 = IPW_2_5,
  IPW_97_5 = IPW_97_5,
  AIPW = AIPW,
  AIPW_2_5 = AIPW_CI[1],
  AIPW_97_5 = AIPW_CI[2],
  CBART = mean(cbart$post_te),
  CBART_2_5 = quantile(rowMeans(cbart$post_te), 0.025, names = F),
  CBART_97_5 = quantile(rowMeans(cbart$post_te), 0.975, names = F),
  model_1 = mean(model_1$post_te),
  model_1_2_5 = quantile(rowMeans(model_1$post_te), 0.025, names = F),
  model_1_97_5 = quantile(rowMeans(model_1$post_te), 0.975, names = F),
  model_2_4 = mean(model_2_4$post_te),
  model_2_4_2_5 = quantile(rowMeans(model_2_4$post_te), 0.025, names = F),
  model_2_4_97_5 = quantile(rowMeans(model_2_4$post_te), 0.975, names = F),
  model_2_4_spline = mean(model_2_4_spline$post_te),
  model_2_4_spline_2_5 = quantile(rowMeans(model_2_4_spline$post_te), 0.025, names = F),
  model_2_4_spline_97_5 = quantile(rowMeans(model_2_4_spline$post_te), 0.975, names = F),
  BART_pure = mean(bart_pure_TE),
  BART_pure_2_5 = quantile(bart_pure_TE, 0.025, names = F),
  BART_pure_97_5 = quantile(bart_pure_TE, 0.975, names = F),
  BART_pure_no_pi = mean(bart_pure_TE_no_pi),
  BART_pure_no_pi_2_5 = quantile(bart_pure_TE_no_pi, 0.025, names = F),
  BART_pure_no_pi_97_5 = quantile(bart_pure_TE_no_pi, 0.975, names = F),
  true_model = mean(true_posterior),
  true_model_2_5 = quantile(true_posterior, 0.025, names = F),
  true_model_97_5 = quantile(true_posterior, 0.975, names = F)
)

write_csv(as.data.frame(t(result)), paste0(path, "/", scenario, "_", seed, ".csv"))

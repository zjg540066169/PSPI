library(dplyr)
source("./R/simulate_2024_Sep_simple.R")
source("./R/AIPW.R")
source("./R/TMLE.R")
library(rstan)
library(brms)
library(RcppArmadillo)
library(Rcpp)
library(RcppDist)
library(boot)
library(SuperLearner)
library(readr)
source("./R/MCMC_BART_Causal_R.R")

path = "./simulation_dec/"
nburn = 3000
npost = 3000

library(doParallel)
library(foreach)
library(MCMCpack)
registerDoParallel(4)
a = foreach(i = 1:4, .combine = list, .multicombine = T) %dopar%{
  MCMC_BART_Causal_R(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 5, nburn, npost, binary, F, seed + i, reverse = F, ntrees_s = 50 * i)
}

for (i in 1:dim(a[[1]]$post_te)[2]) {
  mcmc = mcmc.list(mcmc(a[[1]]$post_te[,i]), mcmc(a[[2]]$post_te[,i]), mcmc(a[[3]]$post_te[,i]), mcmc(a[[4]]$post_te[,i]))
  if(gelman.diag(mcmc, multivariate=F)$psrf[1] >= 1.01)
    print(i)
}

for (i in 1:dim(a[[1]]$post_outcome1)[2]) {
  mcmc = mcmc.list(mcmc(a[[1]]$post_outcome1[,i]), mcmc(a[[2]]$post_outcome1[,i]), mcmc(a[[3]]$post_outcome1[,i]), mcmc(a[[4]]$post_outcome1[,i]))
  if(gelman.diag(mcmc, multivariate=F)$psrf[1] >= 1.01)
    print(i)
}

for (i in 1:dim(a[[1]]$post_outcome0)[2]) {
  mcmc = mcmc.list(mcmc(a[[1]]$post_outcome0[,i]), mcmc(a[[2]]$post_outcome0[,i]), mcmc(a[[3]]$post_outcome0[,i]), mcmc(a[[4]]$post_outcome0[,i]))
  if(gelman.diag(mcmc, multivariate=F)$psrf[1] >= 1.01)
    print(i)
}


args = commandArgs(trailingOnly = T)
seed = as.integer(args[1])
scenario = as.integer(args[2])
binary = as.logical(args[3])
cat("seed:", seed, "\n")
cat("scenario:", scenario, "\n")
#dir.create(paste0(path, seed))


if(scenario == 1){
  data = simulate_overlapping_both(seed = seed, binary = binary)
  cat("scenario: overlapping_both\n")
  set.seed(seed)
}
if(scenario == 2){
  data = simulate_high_outcome_low_selection(seed = seed, binary = binary)
  cat("scenario: high_outcome_low_selection\n")
  set.seed(seed)
}
if(scenario == 3){
  data = simulate_low_outcome_high_selection(seed = seed, binary = binary)
  cat("scenario: low_outcome_high_selection\n")
  set.seed(seed)
}
if(scenario == 4){
  data = simulate_overlapping_both_quadratic(seed = seed, binary = binary)
  cat("scenario: overlapping_both_quadratic\n")
  set.seed(seed)
}



if(binary == FALSE){
  trials = data$trials
  EHR = as.data.frame(data$EHR)
  
  true_TE = mean((data$data$outcome1) - (data$data$outcome0))
  
  sample_mean = mean(trials$outcome[trials$Z == 1]) - mean(trials$outcome[trials$Z == 0])
  
  se = sqrt(var(trials$outcome[trials$Z == 1]) / length(trials$outcome[trials$Z == 1]) + var(trials$outcome[trials$Z == 0]) / length(trials$outcome[trials$Z == 0]))
  
  true = brm(
    formula = outcome ~ Z * (cont_var1 + cont_var2 + cont_var3 + bin_var1 + bin_var2 + bin_var3),
    data = trials,
    family = gaussian(),  # equivalent to a normal response distribution
    chains = 1,
    iter = nburn + npost
  )
  true_predict1 = posterior_predict(true, cbind(EHR, Z = 1))
  true_predict0 = posterior_predict(true, cbind(EHR, Z = 0))
  true_posterior = rowMeans(true_predict1 - true_predict0)
  
  sel = c(rep(1, dim(trials)[1]), rep(0, dim(EHR)[1]))
  co = rbind(trials[,1:20], EHR[,1:20])
  
  bart = brm(
    sel ~ .,
    data = data.frame(sel, co),
    family = bernoulli(link = "logit"),
    chains = 1,
    iter = nburn + npost
  )
  
  trials$e_ps = colMeans(posterior_epred(bart, trials[,1:20]))
  EHR$e_ps = colMeans(posterior_epred(bart, EHR[,1:20]))
  data$data$e_ps = colMeans(posterior_epred(bart, EHR[,1:20]))
  X_s = EHR$e_ps
  
  trials$inv_weight = 1 / (trials$e_ps * ifelse(trials$Z == 1, mean(trials$Z), 1 - mean(trials$Z)))
  
  des <- survey::svydesign(ids = ~1, weights = ~inv_weight, data = trials)
  IPW = unname(survey::svyttest(outcome~I(Z == 1), des)$estimate)
  IPW_2_5 = survey::svyttest(outcome~I(Z == 1), des)$conf.int[1]
  IPW_97_5 = survey::svyttest(outcome~I(Z == 1), des)$conf.int[2]
  
  IPW1 = unname(survey::svyttest(outcome~0, subset(des, Z == 1))$estimate)
  IPW1_2_5 = survey::svyttest(outcome~0, subset(des, Z == 1))$conf.int[1]
  IPW1_97_5 = survey::svyttest(outcome~0, subset(des, Z == 1))$conf.int[2]
  
  IPW0 = unname(survey::svyttest(outcome~0, subset(des, Z == 0))$estimate)
  IPW0_2_5 = survey::svyttest(outcome~0, subset(des, Z == 0))$conf.int[1]
  IPW0_97_5 = survey::svyttest(outcome~0, subset(des, Z == 0))$conf.int[2]
  
  AIPPW = apiw(Y = trials$outcome, Z = trials$Z, X = trials[,1:20], ps = trials$e_ps, X_pop = data$data[,1:20], ps_pop = data$data$e_ps)
  AIPPW_normalize = apiw(Y = trials$outcome, Z = trials$Z, X = trials[,1:20], ps = trials$e_ps, X_pop = data$data[,1:20], ps_pop = data$data$e_ps, normalize = TRUE)
  
  
  
  # AIPW = aipw_function(trials, 1:dim(trials)[1])
  # bootstrap_results <- boot(data = trials, statistic = aipw_function, R = 1000)
  # AIPW_CI = c(boot.ci(bootstrap_results, type = "perc")$percent[4], boot.ci(bootstrap_results, type = "perc")$percent[5])
  # 
  # AIPW1 = aipw_function_treat_mean(trials, 1:dim(trials)[1])
  # bootstrap_results1 <- boot(data = trials, statistic = aipw_function_treat_mean, R = 1000)
  # AIPW_CI1 = c(boot.ci(bootstrap_results1, type = "perc")$percent[4], boot.ci(bootstrap_results1, type = "perc")$percent[5])
  # 
  # AIPW0 = aipw_function_control_mean(trials, 1:dim(trials)[1])
  # bootstrap_results0 <- boot(data = trials, statistic = aipw_function_control_mean, R = 1000)
  # AIPW_CI0 = c(boot.ci(bootstrap_results0, type = "perc")$percent[4], boot.ci(bootstrap_results0, type = "perc")$percent[5])
  
  
  tmle_fit = TMLE_generalizable(Y = trials$outcome, Z = trials$Z, X = trials[,1:20], ps = trials$e_ps, X_pop = data$data[,1:20], ps_pop = data$data$e_ps)
  tmle_fit_logit = TMLE_generalizable(Y = trials$outcome, Z = trials$Z, X = trials[,1:20], ps = trials$e_ps, X_pop = data$data[,1:20], ps_pop = data$data$e_ps, use_logit_link = T)
  
  
  cbart = MCMC_BART_Causal_R(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 2, nburn, npost, binary, F, seed)
  model_1 = MCMC_BART_Causal_R(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 3, nburn, npost, binary, F, seed)
  model_2_4 = MCMC_BART_Causal_R(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 4, nburn, npost, binary, F, seed)
  model_2_4_spline = MCMC_BART_Causal_R(as.matrix(trials[,1:20]), trials$outcome, trials$Z, trials$e_ps, as.matrix(EHR[,1:20]), EHR$e_ps, 5, nburn, npost, binary, F, seed)
  BART_model = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z, invlogit(trials$e_ps))), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
  bart_pure_outcome1 = rowMeans(predict(BART_model, cbind(as.matrix(EHR[,1:20]), 1, invlogit(EHR$e_ps))))
  bart_pure_outcome0 = rowMeans(predict(BART_model, cbind(as.matrix(EHR[,1:20]), 0, invlogit(EHR$e_ps))))
  bart_pure_TE = bart_pure_outcome1 - bart_pure_outcome0
  rm(BART_model)
  
  BART_model_no_pi = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
  bart_pure_no_pi_outcome1 = rowMeans(predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 1)))
  bart_pure_no_pi_outcome0 = rowMeans(predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 0)))
  bart_pure_TE_no_pi = bart_pure_no_pi_outcome1 - bart_pure_no_pi_outcome0
  rm(BART_model_no_pi)
  
  
  result_ATE = data.frame(
    replicates = seed,
    type = "ATE",
    true_TE = true_TE,
    sample_mean = sample_mean,
    sample_mean_2_5 = t.test(trials$outcome[trials$Z == 1], trials$outcome[trials$Z == 0])$conf.int[1], #sample_mean - 1.96 * se,
    sample_mean_97_5 = t.test(trials$outcome[trials$Z == 1], trials$outcome[trials$Z == 0])$conf.int[2], #sample_mean + 1.96 * se
    IPW = IPW,
    IPW_2_5 = IPW_2_5,
    IPW_97_5 = IPW_97_5,
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
    true_model_97_5 = quantile(true_posterior, 0.975, names = F),
    tmle = tmle_fit[1],
    tmle_logit = tmle_fit_logit[1],
    AIPPW = AIPPW[3],
    AIPPW_normalize = AIPPW_normalize[3]
  )
  
  
  result_Y1 = data.frame(
    replicates = seed,
    type = "Y1",
    true_TE = mean(data$data$outcome1),
    sample_mean = mean(trials$outcome[trials$Z == 1]),
    sample_mean_2_5 = t.test(trials$outcome[trials$Z == 1])$conf.int[1], #sample_mean - 1.96 * se,
    sample_mean_97_5 = t.test(trials$outcome[trials$Z == 1])$conf.int[2], #sample_mean + 1.96 * se
    IPW = IPW1,
    IPW_2_5 = IPW1_2_5,
    IPW_97_5 = IPW1_97_5,
    CBART = mean(cbart$post_outcome1),
    CBART_2_5 = quantile(rowMeans(cbart$post_outcome1), 0.025, names = F),
    CBART_97_5 = quantile(rowMeans(cbart$post_outcome1), 0.975, names = F),
    model_1 = mean(model_1$post_outcome1),
    model_1_2_5 = quantile(rowMeans(model_1$post_outcome1), 0.025, names = F),
    model_1_97_5 = quantile(rowMeans(model_1$post_outcome1), 0.975, names = F),
    model_2_4 = mean(model_2_4$post_outcome1),
    model_2_4_2_5 = quantile(rowMeans(model_2_4$post_outcome1), 0.025, names = F),
    model_2_4_97_5 = quantile(rowMeans(model_2_4$post_outcome1), 0.975, names = F),
    model_2_4_spline = mean(model_2_4_spline$post_outcome1),
    model_2_4_spline_2_5 = quantile(rowMeans(model_2_4_spline$post_outcome1), 0.025, names = F),
    model_2_4_spline_97_5 = quantile(rowMeans(model_2_4_spline$post_outcome1), 0.975, names = F),
    BART_pure = mean(bart_pure_outcome1),
    BART_pure_2_5 = quantile(bart_pure_outcome1, 0.025, names = F),
    BART_pure_97_5 = quantile(bart_pure_outcome1, 0.975, names = F),
    BART_pure_no_pi = mean(bart_pure_no_pi_outcome1),
    BART_pure_no_pi_2_5 = quantile(bart_pure_no_pi_outcome1, 0.025, names = F),
    BART_pure_no_pi_97_5 = quantile(bart_pure_no_pi_outcome1, 0.975, names = F),
    true_model = mean(rowMeans(true_predict1)),
    true_model_2_5 = quantile(rowMeans(true_predict1), 0.025, names = F),
    true_model_97_5 = quantile(rowMeans(true_predict1), 0.975, names = F),
    tmle = tmle_fit[2],
    tmle_logit = tmle_fit_logit[2],
    AIPPW = AIPPW[1],
    AIPPW_normalize = AIPPW_normalize[1]
  )
  
  
  result_Y0 = data.frame(
    replicates = seed,
    type = "Y0",
    true_TE = mean(data$data$outcome0),
    sample_mean = mean(trials$outcome[trials$Z == 0]),
    sample_mean_2_5 = t.test(trials$outcome[trials$Z == 0])$conf.int[1], #sample_mean - 1.96 * se,
    sample_mean_97_5 = t.test(trials$outcome[trials$Z == 0])$conf.int[2], #sample_mean + 1.96 * se
    IPW = IPW0,
    IPW_2_5 = IPW0_2_5,
    IPW_97_5 = IPW0_97_5,
    CBART = mean(cbart$post_outcome0),
    CBART_2_5 = quantile(rowMeans(cbart$post_outcome0), 0.025, names = F),
    CBART_97_5 = quantile(rowMeans(cbart$post_outcome0), 0.975, names = F),
    model_1 = mean(model_1$post_outcome0),
    model_1_2_5 = quantile(rowMeans(model_1$post_outcome0), 0.025, names = F),
    model_1_97_5 = quantile(rowMeans(model_1$post_outcome0), 0.975, names = F),
    model_2_4 = mean(model_2_4$post_outcome0),
    model_2_4_2_5 = quantile(rowMeans(model_2_4$post_outcome0), 0.025, names = F),
    model_2_4_97_5 = quantile(rowMeans(model_2_4$post_outcome0), 0.975, names = F),
    model_2_4_spline = mean(model_2_4_spline$post_outcome0),
    model_2_4_spline_2_5 = quantile(rowMeans(model_2_4_spline$post_outcome0), 0.025, names = F),
    model_2_4_spline_97_5 = quantile(rowMeans(model_2_4_spline$post_outcome0), 0.975, names = F),
    BART_pure = mean(bart_pure_outcome0),
    BART_pure_2_5 = quantile(bart_pure_outcome0, 0.025, names = F),
    BART_pure_97_5 = quantile(bart_pure_outcome0, 0.975, names = F),
    BART_pure_no_pi = mean(bart_pure_no_pi_outcome0),
    BART_pure_no_pi_2_5 = quantile(bart_pure_no_pi_outcome0, 0.025, names = F),
    BART_pure_no_pi_97_5 = quantile(bart_pure_no_pi_outcome0, 0.975, names = F),
    true_model = mean(rowMeans(true_predict0)),
    true_model_2_5 = quantile(rowMeans(true_predict0), 0.025, names = F),
    true_model_97_5 = quantile(rowMeans(true_predict0), 0.975, names = F),
    tmle = tmle_fit[3],
    tmle_logit = tmle_fit_logit[3],
    AIPPW = AIPPW[2],
    AIPPW_normalize = AIPPW_normalize[2]
  )
  
  result = rbind(result_ATE, result_Y1, result_Y0)
  
  write_csv(result, paste0(path, "/", scenario, "_", seed, "_", binary, ".csv"))
  
}

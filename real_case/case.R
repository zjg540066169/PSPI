library(tidyverse)
library(dplyr)
library(knitr)
library(mice)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(rstanarm)
library(BART3)
library(RcppProgress)

raw_data = 
  readxl::read_xlsx("/Volumes/DOM_CBCH$/Study Folders/Moise_iHeart DepCare_2018/Data Management/Data Requests/LOCKED DATASETS/Combined Implementation Data/Locked Datasets/Mengxiao_Luan/Clean_Pre_Post_Comb_Implementation_Dataset.xlsx") |>
  janitor::clean_names() 


# baseline
baseline_data = raw_data %>% 
  filter(visit_type == 0)

analyzable_data = baseline_data |>
  select(mrn, visit_type, implementation_phase, 
         dcs_pref_tx:dcs_8c_tool_rec, pam_score, phq_score,
         dob, enrollment_date_preimp, enrollment_date_postimp,
         gender, gender_other, sex, sex_birth, sex_birth_other, ethnicity, 
         race_amer_ind_nat_amer:race_decline, race,
         education, education_years, 
         partner_status, partner_how_long_years, partner_how_long_mths,
         ipaq_30mins, phq8_score, phq_score,
         tx1_depression:tx11_old_e_txnot,
         step_name, clinic_name)


analyzable_data = analyzable_data |>
  mutate(implementation_phase = case_when(implementation_phase == 0 ~ "pre",
                                          implementation_phase == 1 ~ "post"),
         sex = case_when(sex == 1 ~ "Female", sex == 2 ~ "Male"),
         gender = case_when(gender == 0 ~ "Woman",
                            gender == 1 ~ "Man",
                            gender == 6 ~ "Don't Know/Not sure",
                            TRUE ~ "Unknown"),
         ethnicity = case_when(ethnicity == 1 ~ "Hispanic or Latino",
                               ethnicity == 0 ~ "Not Hispanic or Latino",
                               ethnicity == -1 ~ "Decline to respond",
                               TRUE ~ "Ethnicity Unknown/Not Reported"))


analyzable_data = analyzable_data |> 
  mutate(dcs_pref_tx = 
           case_when(dcs_pref_tx == 1 ~ "Antidepressants/Medications",
                     dcs_pref_tx == 2 ~ "Talk Therapy/ Counseling",
                     dcs_pref_tx == 3 ~ "Exercise Program/Cardiac rehab",
                     dcs_pref_tx == 0 ~ "None",
                     dcs_pref_tx == -1 ~ "Unsure"))

analyzable_data = analyzable_data |>
  mutate(across(c(implementation_phase, sex, gender, ethnicity, dcs_pref_tx, step_name, clinic_name), ~ as.factor(.)))

analyzable_data = analyzable_data |>
  mutate(dcs_1_option = dcs_1_option - 1,
         dcs_2_benefit = dcs_2_benefit - 1,
         dcs_3_risks = dcs_3_risks - 1,
         dcs_4_clearben = dcs_4_clearben - 1,
         dcs_5_clearrisk = dcs_5_clearrisk - 1,
         dcs_6_clearimp = dcs_6_clearimp - 1,
         dcs_7_support = dcs_7_support - 1,
         dcs_8_pressure = dcs_8_pressure - 1,
         dcs_9_advice = dcs_9_advice - 1,
         dcs_10_clearbest = dcs_10_clearbest - 1,
         dcs_11_sure = dcs_11_sure - 1,
         dcs_12_easy = dcs_12_easy - 1,
         dcs_13_inform = dcs_13_inform - 1,
         dcs_14_import = dcs_14_import - 1,
         dcs_15_expect = dcs_15_expect - 1)

analyzable_data = analyzable_data |>
  mutate(dcs_1_option_low = case_when(dcs_1_option_low == 0 ~ 4,
                                      dcs_1_option_low == 1 ~ 0,
                                      dcs_1_option_low == 2 ~ 2,
                                      dcs_1_option_low == 3 ~ 2,
                                      dcs_1_option_low == -1 ~ NA),
         dcs_2_benefit_low = case_when(dcs_2_benefit_low == 0 ~ 4,
                                       dcs_2_benefit_low == 1 ~ 0,
                                       dcs_2_benefit_low == 2 ~ 2,
                                       dcs_2_benefit_low == 3 ~ 2,
                                       dcs_2_benefit_low == -1 ~ NA),
         dcs_3_risks_low = case_when(dcs_3_risks_low == 0 ~ 4,
                                     dcs_3_risks_low == 1 ~ 0,
                                     dcs_3_risks_low == 2 ~ 2,
                                     dcs_3_risks_low == 3 ~ 2,
                                     dcs_3_risks_low == -1 ~ NA),
         dcs_4_clearben_low = case_when(dcs_4_clearben_low == 0 ~ 4,
                                        dcs_4_clearben_low == 1 ~ 0,
                                        dcs_4_clearben_low == 2 ~ 2,
                                        dcs_4_clearben_low == 3 ~ 2,
                                        dcs_4_clearben_low == -1 ~ NA),
         dcs_5_clearrisk_low = case_when(dcs_5_clearrisk_low == 0 ~ 4,
                                         dcs_5_clearrisk_low == 1 ~ 0,
                                         dcs_5_clearrisk_low == 2 ~ 2,
                                         dcs_5_clearrisk_low == 3 ~ 2,
                                         dcs_5_clearrisk_low == -1 ~ NA),
         dcs_6_support_low = case_when(dcs_6_support_low == 0 ~ 4,
                                       dcs_6_support_low == 1 ~ 0,
                                       dcs_6_support_low == 2 ~ 2,
                                       dcs_6_support_low == 3 ~ 2,
                                       dcs_6_support_low == -1 ~ NA),
         dcs_7_pressure_low = case_when(dcs_7_pressure_low == 0 ~ 4,
                                        dcs_7_pressure_low == 1 ~ 0,
                                        dcs_7_pressure_low == 2 ~ 2,
                                        dcs_7_pressure_low == 3 ~ 2,
                                        dcs_7_pressure_low == -1 ~ NA),
         dcs_8_advice_low = case_when(dcs_8_advice_low == 0 ~ 4,
                                      dcs_8_advice_low == 1 ~ 0,
                                      dcs_8_advice_low == 2 ~ 2,
                                      dcs_8_advice_low == 3 ~ 2,
                                      dcs_8_advice_low == -1 ~ NA),
         dcs_9_clearbest_low = case_when(dcs_9_clearbest_low == 0 ~ 4,
                                         dcs_9_clearbest_low == 1 ~ 0,
                                         dcs_9_clearbest_low == 2 ~ 2,
                                         dcs_9_clearbest_low == 3 ~ 2,
                                         dcs_9_clearbest_low == -1 ~ NA),
         dcs_10_sure_low = case_when(dcs_10_sure_low == 0 ~ 4,
                                     dcs_10_sure_low == 1 ~ 0,
                                     dcs_10_sure_low == 2 ~ 2,
                                     dcs_10_sure_low == 3 ~ 2,
                                     dcs_10_sure_low == -1 ~ NA),)

## race
analyzable_data = analyzable_data |> 
  rowwise() |>
  mutate(race_cnt = sum(race_amer_ind_nat_amer, race_asian, race_black, race_hawaiian_pi, race_white, race_other, race_unknown, na.rm = TRUE))

analyzable_data = analyzable_data |>
  mutate(race_fct = case_when(
    race_amer_ind_nat_amer == 1 ~ 1, 
    race_asian == 1 ~ 2, 
    race_black == 1 ~ 3, 
    race_hawaiian_pi == 1 ~ 4, 
    race_white == 1 ~ 5, 
    race_cnt > 1 | race_more_than_one == 1 ~ 6,
    race_other == 1 ~ 7,
    TRUE ~ 8
  ))

analyzable_data$race_fct = 
  factor(analyzable_data$race_fct, ordered=TRUE, levels=c(1,2,3,4,5,6,7,8))

levels(analyzable_data$race_fct) = c(
  "American Indian/Alaska Native",
  "Asian",
  "Black or African American",
  "Native Hawaiian/Pacific Islander",
  "White",
  "More than one race",
  "Non-specific",
  "Race Unknown/Not Reported"
)

## age
analyzable_data = analyzable_data |>
  mutate(enrollment_date_preimp = 
           as.Date(enrollment_date_preimp, format="%m/%d/%Y"),
         enrollment_date_postimp = 
           as.Date(enrollment_date_postimp, format="%m/%d/%Y"),
         dob = as.Date(dob, format = "%m/%d/%Y")) %>%
  mutate(enrollment_date = ifelse(implementation_phase == "pre", enrollment_date_preimp, enrollment_date_postimp)) |>
  mutate(enrollment_date = as.Date(enrollment_date, format = "%m/%d/%Y")) |>
  mutate(age = floor(as.numeric((enrollment_date - dob)/365.25)))

analyzable_data |> filter(is.na(age)) |> pull(dob)

## education, partner
analyzable_data = analyzable_data |>
  mutate(education = 
           case_when(education == 1 ~ "Less than high school",
                     education == 2 ~ "Some high school",
                     education == 3 ~ "High school diploma/ GED",
                     education == 4 ~ "Trade school/ Vocational school",
                     education == 5 ~ "Some college",
                     education == 6 ~ "College graduate",
                     education == 7 ~ "Graduate school/ Professional school",
                     education == -1 ~ "Decline to respond"),
         partner_status = 
           case_when(partner_status == 1 ~ "Single",
                     partner_status == 2 ~ "Partner / Spouse",
                     partner_status == 3 ~ "Separated",
                     partner_status == 4 ~ "Widowed",
                     partner_status == 5 ~ "Divorced",
                     partner_status == -1 ~ "Decline to respond"))

## dcs score
### both versions missing
analyzable_data = analyzable_data |>
  mutate(dcs_na = ifelse(is.na(dcs_1_option) & is.na(dcs_2_benefit) & is.na(dcs_3_risks) & is.na(dcs_4_clearben) & is.na(dcs_5_clearrisk) & is.na(dcs_6_clearimp) & is.na(dcs_7_support) & is.na(dcs_8_pressure) & is.na(dcs_9_advice) & is.na(dcs_10_clearbest) & is.na(dcs_11_sure) & is.na(dcs_12_easy) & is.na(dcs_13_inform) & is.na(dcs_14_import) & is.na(dcs_15_expect & is.na(dcs_1_option_low) & is.na(dcs_2_benefit_low) & is.na(dcs_3_risks_low) & is.na(dcs_4_clearben_low) & is.na(dcs_5_clearrisk_low) & is.na(dcs_6_support_low) & is.na(dcs_7_pressure_low) & is.na(dcs_8_advice_low) & is.na(dcs_9_clearbest_low) & is.na(dcs_10_sure_low)), 1, 0))

analyzable_data = analyzable_data |>
  mutate(dcs_version = case_when(is.na(dcs_1_option_low) & is.na(dcs_2_benefit_low) & is.na(dcs_3_risks_low) & is.na(dcs_4_clearben_low) & is.na(dcs_5_clearrisk_low) & is.na(dcs_6_support_low) & is.na(dcs_7_pressure_low) & is.na(dcs_8_advice_low) & is.na(dcs_9_clearbest_low) & is.na(dcs_10_sure_low) & is.na(dcs_1_option) & is.na(dcs_2_benefit) & is.na(dcs_3_risks) & is.na(dcs_4_clearben) & is.na(dcs_5_clearrisk) & is.na(dcs_6_clearimp) & is.na(dcs_7_support) & is.na(dcs_8_pressure) & is.na(dcs_9_advice) & is.na(dcs_10_clearbest) & is.na(dcs_11_sure) & is.na(dcs_12_easy) & is.na(dcs_13_inform) & is.na(dcs_14_import) & is.na(dcs_15_expect) ~ "missing",
                                 is.na(dcs_1_option_low) & is.na(dcs_2_benefit_low) & is.na(dcs_3_risks_low) & is.na(dcs_4_clearben_low) & is.na(dcs_5_clearrisk_low) & is.na(dcs_6_support_low) & is.na(dcs_7_pressure_low) & is.na(dcs_8_advice_low) & is.na(dcs_9_clearbest_low) & is.na(dcs_10_sure_low) ~ "high",
                                 is.na(dcs_1_option) & is.na(dcs_2_benefit) & is.na(dcs_3_risks) & is.na(dcs_4_clearben) & is.na(dcs_5_clearrisk) & is.na(dcs_6_clearimp) & is.na(dcs_7_support) & is.na(dcs_8_pressure) & is.na(dcs_9_advice) & is.na(dcs_10_clearbest) & is.na(dcs_11_sure) & is.na(dcs_12_easy) & is.na(dcs_13_inform) & is.na(dcs_14_import) & is.na(dcs_15_expect) ~ "low",
                                 TRUE ~ "both")) |> 
  mutate(dcs_score = case_when(dcs_version == "low" ~ sum(dcs_1_option_low, dcs_2_benefit_low, dcs_3_risks_low, dcs_4_clearben_low, dcs_5_clearrisk_low, dcs_6_support_low, dcs_7_pressure_low, dcs_8_advice_low, dcs_9_clearbest_low, dcs_10_sure_low) / 10 * 25,
                               dcs_version == "high" ~ sum(dcs_1_option, dcs_2_benefit, dcs_3_risks, dcs_4_clearben, dcs_5_clearrisk, dcs_7_support, dcs_8_pressure, dcs_9_advice, dcs_10_clearbest, dcs_11_sure) / 10 * 25,
                               dcs_version == "both" ~ sum(dcs_1_option_low, dcs_2_benefit_low, dcs_3_risks_low, dcs_4_clearben_low, dcs_5_clearrisk_low, dcs_6_support_low, dcs_7_pressure_low, dcs_8_advice_low, dcs_9_clearbest_low, dcs_10_sure_low) / 10 * 25,
                               dcs_version == "missing" ~ NA))

analyzable_data = analyzable_data |>
  mutate(dcs_score_l = sum(dcs_1_option_low, dcs_2_benefit_low, dcs_3_risks_low, dcs_4_clearben_low, dcs_5_clearrisk_low, dcs_6_support_low, dcs_7_pressure_low, dcs_8_advice_low, dcs_9_clearbest_low, dcs_10_sure_low) / 10 * 25,
         dcs_score_h = sum(dcs_1_option, dcs_2_benefit, dcs_3_risks, dcs_4_clearben, dcs_5_clearrisk, dcs_7_support, dcs_8_pressure, dcs_9_advice, dcs_10_clearbest, dcs_11_sure) / 10 * 25)

analyzable_data = analyzable_data |>
  mutate(dcs_1_option = dcs_1_option + 1,
         dcs_2_benefit = dcs_2_benefit + 1,
         dcs_3_risks = dcs_3_risks + 1,
         dcs_4_clearben = dcs_4_clearben + 1,
         dcs_5_clearrisk = dcs_5_clearrisk + 1,
         dcs_6_clearimp = dcs_6_clearimp + 1,
         dcs_7_support = dcs_7_support + 1,
         dcs_8_pressure = dcs_8_pressure + 1,
         dcs_9_advice = dcs_9_advice + 1,
         dcs_10_clearbest = dcs_10_clearbest + 1,
         dcs_11_sure = dcs_11_sure + 1,
         dcs_12_easy = dcs_12_easy + 1,
         dcs_13_inform = dcs_13_inform + 1,
         dcs_14_import = dcs_14_import + 1,
         dcs_15_expect = dcs_15_expect + 1,
         dcs_1_option_low = dcs_1_option_low / 2 + 1,
         dcs_2_benefit_low = dcs_2_benefit_low / 2 + 1,
         dcs_3_risks_low = dcs_3_risks_low / 2 + 1,
         dcs_4_clearben_low = dcs_4_clearben_low / 2 + 1,
         dcs_5_clearrisk_low = dcs_5_clearrisk_low / 2 + 1,
         dcs_6_support_low = dcs_6_support_low / 2 + 1,
         dcs_7_pressure_low = dcs_7_pressure_low / 2 + 1,
         dcs_8_advice_low = dcs_8_advice_low / 2 + 1,
         dcs_9_clearbest_low = dcs_9_clearbest_low / 2 + 1,
         dcs_10_sure_low = dcs_10_sure_low / 2 + 1)

impute_data = analyzable_data |>
  mutate(gender = ifelse(gender == "Don't Know/Not sure" | gender == "Unknown",
                         NA, gender),
         race_fct = ifelse(race_fct == "Race Unknown/Not Reported" | 
                             race_fct == "Non-Specific",
                           NA, race_fct),
         ethnicity = ifelse(ethnicity == "Decline to respond" | 
                              ethnicity == "Ethnicity Unknown/Not Reported",
                            NA, ethnicity)) |>
  select(mrn:dcs_pref_tx, education, education_years, partner_status, ipaq_30mins, phq_score, phq8_score, age, sex, gender, race_fct, ethnicity, step_name, clinic_name, dcs_version, pam_score) |>
  mutate(across(c(implementation_phase, education, partner_status, sex, gender, ethnicity, dcs_pref_tx, step_name, clinic_name), ~ as.factor(.)))

ini = mice(impute_data, m = 5, maxit = 0, print = FALSE, seed = 2024)
pred = ini$pred
meth = ini$meth
imp = mice(impute_data, pred = pred, meth = meth, m = 5, maxit = 5,
           print = T, seed = 2024)

plot(imp)
com_data = complete(imp, "all")

data = com_data[[1]]
data$calc_eligibility_status = baseline_data$calc_eligibility_status


data = data %>% 
  mutate(
    trial = ifelse(!is.na(calc_eligibility_status) & (calc_eligibility_status == "Eligible" | 
                     calc_eligibility_status == "Eligible PHQ8"), 1, 0)
  ) %>% 
  dplyr::select(-mrn, -visit_type)
data[data$trial == 0, "pam_score"] = NA



library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppProgress)
sourceCpp("./src/MCMC_BART_Causal.cpp")
nburn = 10000
npost = 10000


bart = gbart(dplyr::select(data, -pam_score, -trial, -implementation_phase, -calc_eligibility_status) %>% mutate(
  dcs_pref_tx = as.factor(dcs_pref_tx),
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  step_name = as.factor(step_name),
  clinic_name = as.factor(clinic_name),
  dcs_version = as.factor(dcs_version),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
),data$trial, type = "pbart")

logistic = stan_glm(trial ~ . , family = binomial(link = "logit"), data = dplyr::select(data, -pam_score, -implementation_phase, -calc_eligibility_status) %>% mutate(
  dcs_pref_tx = as.factor(dcs_pref_tx),
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  step_name = as.factor(step_name),
  clinic_name = as.factor(clinic_name),
  dcs_version = as.factor(dcs_version),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
), chains = 1, prior=NULL, iter = 1000 + 1000)

e_ps = predict(bart, dplyr::select(data, -pam_score, -trial, -implementation_phase, -calc_eligibility_status) %>% mutate(
  dcs_pref_tx = as.factor(dcs_pref_tx),
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  step_name = as.factor(step_name),
  clinic_name = as.factor(clinic_name),
  dcs_version = as.factor(dcs_version),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
))$prob.test.mean

data_X = dplyr::select(data, -calc_eligibility_status) %>% 
  #filter(trial == 1) %>% 
  mutate(
  dcs_pref_tx = as.factor(dcs_pref_tx),
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  step_name = as.factor(step_name),
  clinic_name = as.factor(clinic_name),
  dcs_version = as.factor(dcs_version),
  pam_score = sqrt(as.numeric(pam_score)),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
) %>% 
  model.matrix(~ 0 + dcs_pref_tx + education + partner_status + sex + clinic_name + race_fct + ethnicity + step_name  + dcs_version + education_years + ipaq_30mins + phq_score + phq8_score + age, data = .)
  


cbart = MCMC_BART_Causal(as.matrix(data_X[data$trial == 1, ]), sqrt(data$pam_score[data$trial == 1]), (as.integer(data$implementation_phase) - 1)[data$trial == 1], e_ps[data$trial == 1], as.matrix(data_X[data$trial == 0, ]), e_ps[data$trial == 0], 2, nburn, npost, F)
model_1 = MCMC_BART_Causal(as.matrix(data_X[data$trial == 1, ]), sqrt(data$pam_score[data$trial == 1]), (as.integer(data$implementation_phase) - 1)[data$trial == 1], e_ps[data$trial == 1], as.matrix(data_X[data$trial == 0, ]), e_ps[data$trial == 0], 3, nburn, npost, F)
model_2_4 = MCMC_BART_Causal(as.matrix(data_X[data$trial == 1, ]), sqrt(data$pam_score[data$trial == 1]), (as.integer(data$implementation_phase) - 1)[data$trial == 1], e_ps[data$trial == 1], as.matrix(data_X[data$trial == 0, ]), e_ps[data$trial == 0], 4, nburn, npost, F)
model_2_4_spline = MCMC_BART_Causal(as.matrix(data_X[data$trial == 1, ]), sqrt(data$pam_score[data$trial == 1]), (as.integer(data$implementation_phase) - 1)[data$trial == 1], e_ps[data$trial == 1], as.matrix(data_X[data$trial == 0, ]), e_ps[data$trial == 0], 5, nburn, npost, F)
#BART_model = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z, trials$e_ps)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
#bart_pure_TE = rowMeans(predict(BART_model, cbind(as.matrix(EHR[,1:20]), 1, EHR$e_ps)) - predict(BART_model, cbind(as.matrix(EHR[,1:20]), 0, EHR$e_ps)))
#BART_model_no_pi = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
#bart_pure_TE_no_pi = rowMeans(predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 1)) - predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 0)))

mean(cbart$post_outcome1^2 - cbart$post_outcome0^2)
mean(model_1$post_outcome1^2 - model_1$post_outcome0^2)
mean(model_2_4$post_outcome1^2 - model_2_4$post_outcome0^2)
mean(model_2_4_spline$post_outcome1^2 - model_2_4_spline$post_outcome0^2)


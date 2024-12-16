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
library(table1)

raw_data = 
  readxl::read_xlsx("/Volumes/DOM_CBCH$/Study Folders/Moise_iHeart DepCare_2018/Data Management/Data Requests/LOCKED DATASETS/Combined Implementation Data/Locked Datasets/Mengxiao_Luan/Clean_Pre_Post_Comb_Implementation_Dataset.xlsx") |>
  janitor::clean_names()

raw_data = raw_data %>% 
  filter(clinic_name %in% c("RKL HUD 26 Indian Rock Suffern", "CUDOC CIVT", "AIM WEST", "AIM EAST")) %>% 
  filter(mrn != "1103171844")# %>% 
  #filter(!(mrn %in% c("1000120720", "1002470551", "1002530885", "1002755608", "1003717284", "1004188270", "1004730796", "1007542328", "1007560519", "1009339185", "1009547268", "1101060289")))


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
                            gender == 6 ~ NA,
                            TRUE ~ NA),
         ethnicity = case_when(ethnicity == 1 ~ "Hispanic or Latino",
                               ethnicity == 0 ~ "Not Hispanic or Latino",
                               ethnicity == -1 ~ NA,
                               TRUE ~ NA))


analyzable_data = analyzable_data |> 
  mutate(dcs_pref_tx = 
           case_when(dcs_pref_tx == 1 ~ "Antidepressants/Medications",
                     dcs_pref_tx == 2 ~ "Talk Therapy/ Counseling",
                     dcs_pref_tx == 3 ~ "Exercise Program/Cardiac rehab",
                     dcs_pref_tx == 0 ~ NA,
                     dcs_pref_tx == -1 ~ NA))

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
    race_amer_ind_nat_amer == 1 ~ "Other", 
    race_asian == 1 ~ "Other",
    race_black == 1 ~ "Black or African American",
    race_hawaiian_pi == 1 ~ "Other", 
    race_white == 1 ~ "White",
    race_cnt > 1 | race_more_than_one == 1 ~ "More than one race",
    race_other == 1 ~ "Other",
    TRUE ~ NA
  )) 

analyzable_data$race_fct = as.factor(analyzable_data$race_fct)
  #factor(analyzable_data$race_fct, ordered=TRUE, levels=c(1,2,3,4,5,6,7,8))

# levels(analyzable_data$race_fct) = c(
#   "American Indian/Alaska Native",
#   "Asian",
#   "Black or African American",
#   "Native Hawaiian/Pacific Islander",
#   "White",
#   "More than one race",
#   "Non-specific",
#   "Race Unknown/Not Reported"
# )

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
           case_when(education == 1 ~ "High School and Below",
                     education == 2 ~ "High School and Below",
                     education == 3 ~ "High School and Below",
                     education == 4 ~ "Above High School",
                     education == 5 ~ "Above High School",
                     education == 6 ~ "Above High School",
                     education == 7 ~ "Above High School",
                     education == -1 ~ NA),
         partner_status = 
           case_when(partner_status == 1 ~ "Single",
                     partner_status == 2 ~ "Partner / Spouse",
                     partner_status == 3 ~ "Previously Married",
                     partner_status == 4 ~ "Previously Married",
                     partner_status == 5 ~ "Previously Married",
                     partner_status == -1 ~ NA))

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
  select(mrn:dcs_pref_tx, education, education_years, partner_status, phq_score, phq8_score, age, sex, race_fct, ethnicity, clinic_name, pam_score) |>
  mutate(across(c(implementation_phase, education, partner_status, sex, ethnicity, dcs_pref_tx, clinic_name), ~ as.factor(.))) %>% 
  mutate(
    phq_score = ifelse(phq_score == -1, NA, phq_score),
    phq8_score = ifelse(phq8_score == -1, NA, phq8_score),
  )





# 
# table1(~ clinic_name + age + sex + gender + as.factor(race_fct) +ethnicity + education + education_years + dcs_pref_tx + dcs_version + partner_status + ipaq_30mins + phq_score + phq8_score + step_name + pam_score | trial, 
#        data = impute_data %>% as_tibble() %>% 
#          mutate(calc_eligibility_status = baseline_data$calc_eligibility_status, 
#                 trial = ifelse(!is.na(calc_eligibility_status) & (calc_eligibility_status == "Eligible" | calc_eligibility_status == "Eligible PHQ8"), "Trial", "Out of Trial")),
#        extra.col=list(`P-value`=pvalue))

impute_data = 
  impute_data %>% as_tibble() %>% select(-dcs_pref_tx, -visit_type) %>% 
  mutate(calc_eligibility_status = baseline_data$calc_eligibility_status, 
         trial = ifelse(!is.na(calc_eligibility_status) & (calc_eligibility_status == "Eligible" | calc_eligibility_status == "Eligible PHQ8"), "Trial", "Out of Trial"))   %>% 
  select(-calc_eligibility_status)

table1(~ clinic_name + age + sex + race_fct + ethnicity + education + education_years +  partner_status + phq_score + phq8_score +  pam_score  | trial * implementation_phase, 
       data = impute_data)
                                                                                                                                                                                                                                                                                                                                       

  

ini = mice(impute_data, m = 5, maxit = 0, print = FALSE, seed = 2024)
pred = ini$pred
meth = ini$meth
imp = mice(impute_data, pred = pred, meth = meth, m = 5, maxit = 5,
           print = T, seed = 2024)

plot(imp)
com_data = complete(imp, "all")

data = com_data[[1]]
#data$calc_eligibility_status = baseline_data$calc_eligibility_status


data = data %>% 
  mutate(
    implementation_phase = ifelse(implementation_phase == "pre", 0, 1)
  ) %>% 
  dplyr::select(-mrn)






data[data$trial == "Out of Trial", "pam_score"] = NA



library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppProgress)
sourceCpp("./src/MCMC_BART_Causal.cpp")
nburn = 5000
npost = 5000


bart = gbart(dplyr::select(data, -pam_score, -trial, -implementation_phase) %>% mutate(
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  clinic_name = as.factor(clinic_name),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
), data$trial == "Trial", type = "pbart")

logistic = stan_glm(I(trial == "Trial") ~ . , family = binomial(link = "logit"), data = dplyr::select(data, -pam_score, -implementation_phase) %>% mutate(
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  clinic_name = as.factor(clinic_name),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
), chains = 1, prior=NULL, iter = 1000 + 1000)

e_ps = predict(bart, dplyr::select(data, -pam_score, -trial, -implementation_phase) %>% mutate(
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  clinic_name = as.factor(clinic_name),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
))$prob.test.mean

e_ps2 = colMeans(posterior_epred(logistic, data = dplyr::select(data, -pam_score, -implementation_phase) %>% mutate(
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  clinic_name = as.factor(clinic_name),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity))))



data_X = data %>% 
  #filter(trial == 1) %>% 
  mutate(
  education = as.factor(education),
  partner_status = as.factor(partner_status),
  sex = as.factor(sex),
  clinic_name = as.factor(clinic_name),
  pam_score = sqrt(as.numeric(pam_score)),
  race_fct = as.factor(race_fct),
  ethnicity = as.factor(ethnicity)
) %>% 
  model.matrix(~ 0 + education + partner_status + sex + clinic_name + race_fct + ethnicity + education_years + phq_score + phq8_score + age, data = .)
  


cbart = MCMC_BART_Causal_R(as.matrix(data_X[data$trial == "Trial", ]), sqrt(data$pam_score[data$trial == "Trial"]), (as.integer(data$implementation_phase))[data$trial == "Trial"], e_ps2[data$trial == "Trial"], as.matrix(data_X), e_ps2, 2, nburn, npost, F, F, 123)
model_1 = MCMC_BART_Causal_R(as.matrix(data_X[data$trial == "Trial", ]), sqrt(data$pam_score[data$trial == "Trial"]), (as.integer(data$implementation_phase))[data$trial == "Trial"], e_ps2[data$trial == "Trial"], as.matrix(data_X), e_ps2, 3, nburn, npost, F, F, 123)
model_2_4 = MCMC_BART_Causal_R(as.matrix(data_X[data$trial == "Trial", ]), sqrt(data$pam_score[data$trial == "Trial"]), (as.integer(data$implementation_phase))[data$trial == "Trial"], e_ps[data$trial == "Trial"], as.matrix(data_X), e_ps, 4, nburn, npost, F, F, 123)
model_2_4_spline = MCMC_BART_Causal_R(as.matrix(data_X[data$trial == "Trial", ]), sqrt(data$pam_score[data$trial == "Trial"]), (as.integer(data$implementation_phase))[data$trial == "Trial"], e_ps[data$trial == "Trial"], as.matrix(data_X), e_ps, 5, nburn, npost, F, F, 123)
#BART_model = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z, trials$e_ps)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
#bart_pure_TE = rowMeans(predict(BART_model, cbind(as.matrix(EHR[,1:20]), 1, EHR$e_ps)) - predict(BART_model, cbind(as.matrix(EHR[,1:20]), 0, EHR$e_ps)))
#BART_model_no_pi = BART::wbart(as.matrix(cbind(trials[,1:20], trials$Z)), as.numeric(trials$outcome), ndpost=npost, nskip = nburn, rm.const = F)
#bart_pure_TE_no_pi = rowMeans(predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 1)) - predict(BART_model_no_pi, cbind(as.matrix(EHR[,1:20]), 0)))

mean(cbart$post_outcome1^2 - cbart$post_outcome0^2)
mean(model_1$post_outcome1^2 - model_1$post_outcome0^2)
mean(model_2_4$post_outcome1^2 - model_2_4$post_outcome0^2)
mean(model_2_4_spline$post_outcome1^2 - model_2_4_spline$post_outcome0^2)


library(tidyverse)
dic = "./simulation_reverse/"
dic = "./simulation_fit0/"
# for (i in 1:200) {
#   if (file.exists(paste0(dic, "1_", i, "_no_reverse.csv")))
#     data = rbind(data, read_csv(paste0(dic, "1_", i, "_no_reverse.csv"), show_col_types = F))
#   else
#       print(i)
# }

# scn 1: generation, prop of treatment = 0.5
# scn 2: generation, prop of treatment = 0.2
# scn 3: reverse generation, prop of treatment = 0.5
# scn 4: reverse generation, prop of treatment = 0.2

# outcome = f(X) + A*g(X)


a = 0
scn = 4
data = tibble()
data_reverse = tibble()
for (i in 1:100) {
  if (file.exists(paste0(dic, scn, "_", i, "_no_reverse.csv")) & file.exists(paste0(dic, scn, "_", i, "_reverse.csv"))){
    data = rbind(data, read_csv(paste0(dic,scn, "_", i, "_no_reverse.csv"), show_col_types = F))
    data_reverse = rbind(data_reverse, read_csv(paste0(dic, scn,"_", i, "_reverse.csv"), show_col_types = F))
  }else{
    a = a + 1
  }
}


data_reverse %>% 
  as_tibble() %>%
  dplyr::select(c(1,2, 3,4,7,10,13,16,19,22,23,24,25)) %>% 
  pivot_longer(cols = c(4,5,6,7,8,9,10,11,12,13), names_to = "model", values_to = "estimates") %>%
  mutate(bias = (estimates - true_TE), mse = (estimates - true_TE)^2) %>% 
  group_by(model, type) %>%
  summarize(total_bias = mean(bias), rmse = sqrt(mean(mse))) %>% 
  filter(type == "ATE") %>% 
  cbind(
    data %>% 
      as_tibble() %>%
      dplyr::select(c(1,2, 3,4,7,10,13,16,19,22,23,24,25)) %>% 
      pivot_longer(cols = c(4,5,6,7,8,9,10,11,12,13), names_to = "model", values_to = "estimates") %>%
      mutate(bias = (estimates - true_TE), mse = (estimates - true_TE)^2) %>% 
      group_by(model, type) %>%
      summarize(total_bias = mean(bias), rmse = sqrt(mean(mse))) %>% 
      filter(type == "ATE")
  )


data_reverse %>% 
  as_tibble() %>%
  dplyr::select(c(1,2, 3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)) %>% 
  pivot_longer(cols = c(4,6,8,10,12,14,16,18,20), names_to = "model", values_to = "CI_2_5") %>%
  pivot_longer(cols = c(4,5,6,7,8,9,10,11,12), names_to = "model2", values_to = "CI_97_5") %>% 
  mutate(model = str_remove(model, "_2_5"), model2 = str_remove(model2, "_97_5")) %>% 
  filter(model == model2) %>% 
  mutate(width = CI_97_5 - CI_2_5, coverage = (true_TE <= CI_97_5 & true_TE >= CI_2_5)) %>%
  group_by(model, type) %>%
  summarize(coverage_rate = mean(coverage), average_width = mean(width)) %>% 
  filter(type == "ATE") %>% 
  cbind(data %>% 
          as_tibble() %>%
          dplyr::select(c(1,2, 3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)) %>% 
          pivot_longer(cols = c(4,6,8,10,12,14,16,18,20), names_to = "model", values_to = "CI_2_5") %>%
          pivot_longer(cols = c(4,5,6,7,8,9,10,11,12), names_to = "model2", values_to = "CI_97_5") %>% 
          mutate(model = str_remove(model, "_2_5"), model2 = str_remove(model2, "_97_5")) %>% 
          filter(model == model2) %>% 
          mutate(width = CI_97_5 - CI_2_5, coverage = (true_TE <= CI_97_5 & true_TE >= CI_2_5)) %>%
          group_by(model, type) %>%
          summarize(coverage_rate = mean(coverage), average_width = mean(width)) %>% 
          filter(type == "ATE"))
 
  

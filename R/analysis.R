library(tidyverse)
dic = "./simulation_dec/simulation_dec/"
data = tibble()
for (i in 1:200) {
  if (file.exists(paste0(dic, "4_", i, ".csv")))
    data = rbind(data, read_csv(paste0(dic, "4_", i, ".csv"), show_col_types = F))
  else
      print(i)
}
data %>% 
  as_tibble() %>%
  pivot_longer(cols = c(3,6,9, 12, 15, 18, 21, 24, 27, 30), names_to = "model", values_to = "estimates") %>%
  mutate(
    replicates = unlist(replicates),
    estimates = unlist(estimates)
  ) %>%
  pivot_longer(cols = c(3,5, 7, 9, 11, 13, 15, 17, 19, 21), names_to = "model_2", values_to = "CI_2_5") %>%
  pivot_longer(cols = c(3,4, 5, 6, 7, 8, 9, 10, 11, 12), names_to = "model_3", values_to = "CI_97_5") %>%
  filter(str_detect(model_2, model) & str_detect(model_3, model)) %>%
  mutate(bias = (estimates - true_TE), mse = (estimates - true_TE)^2, width = CI_97_5 - CI_2_5, coverage = (true_TE <= CI_97_5 & true_TE >= CI_2_5)) %>%
  group_by(model) %>%
  summarize(total_bias = mean(bias), rmse = sqrt(mean(mse)), coverage_rate = mean(coverage), average_width = mean(width))

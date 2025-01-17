library(tidyverse)
dic = "./R/Jan/results_simplest_200//"
r = tibble()
files = list.files(dic,full.names = T)
for (i in files) {
  #if(str_detect(i, "standard_linear.csv"))
  r = rbind(r, read_csv(i, show_col_types = F))
  #if (file.exists(paste0("./R/test_standard/simulation_standard_linear2/", i, "_0.5_standard_linear.csv")))
  #result = rbind(result, read_csv(paste0(dic, i, "_0.5_conti_simplest.csv"), show_col_types = F))
  #else
  #  print(i)
}


r %>% as_tibble() %>% 
  #filter(prop == 0.2) %>% 
  pivot_longer(
    cols = ends_with("_2_5"),
    names_to = "model",
    values_to = "CI_25",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_2_5"
  ) %>% 
  pivot_longer(
    cols = ends_with("_97_5"),
    names_to = "model2",
    values_to = "CI_975",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_97_5"
  ) %>% 
  pivot_longer(
    #cols = 4:17,
    cols = ends_with("_est"),
    names_to = "model3",
    values_to = "est",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_est"
  ) %>% 
  group_by(replicates, type, model3) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(type, model3) %>% 
  summarise(
    bias = mean(est - true_TE),
    rmse = sqrt(mean((est - true_TE)^2))
  ) %>% 
  left_join(
    r %>% 
      #filter(prop == 0.2) %>% 
      pivot_longer(
        cols = ends_with("_2_5"),
        names_to = "model",
        values_to = "CI_25",
        # This regex says: capture everything before "_data"
        names_pattern = "(.*)_2_5"
      ) %>% 
      pivot_longer(
        cols = ends_with("_97_5"),
        names_to = "model2",
        values_to = "CI_975",
        # This regex says: capture everything before "_data"
        names_pattern = "(.*)_97_5"
      ) %>% 
      pivot_longer(
        cols = ends_with("_est"),
        names_to = "model3",
        values_to = "est",
        # This regex says: capture everything before "_data"
        names_pattern = "(.*)_est"
      ) %>% 
      filter(model == model2 & model == model3) %>% 
      mutate(width = CI_975 - CI_25, coverage = (true_TE <= CI_975 & true_TE >= CI_25)) %>% 
      group_by(type, model3) %>% 
      summarise(
        CP = mean(coverage),
        AW = mean(width)
      ) 
  )  %>%  View()


Y ~ f(X, Z)
f = p1*f1 + p2*f2 + ...
Y - f(X, Z) ~ W * (-1)^Z * g(ps)

Y_hat_1 = f(X, 1) - g(ps)
Y_hat_0 = f(X, 0) + g(ps)



data %>% 
  pivot_longer(
    cols = ends_with("_2_5"),
    names_to = "model",
    values_to = "CI_25",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_2_5"
  ) %>% 
  pivot_longer(
    cols = ends_with("_97_5"),
    names_to = "model2",
    values_to = "CI_975",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_97_5"
  ) %>% 
  pivot_longer(
    #cols = 4:17,
    cols = ends_with("_est"),
    names_to = "model3",
    values_to = "est",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_est"
  ) %>% 
  group_by(replicates, type, model3) %>% 
  slice(1) %>% 
  ungroup() %>% 
  #group_by(type, model3) %>% 
  mutate(
    bias = ((est - true_TE)),
    rmse = ((est - true_TE)^2)
  ) %>% 
  filter(type == "ATE" & model3 == "model_2_4_spline") %>% 
  summary(.$rmse)

  













data %>% 
  pivot_longer(
    cols = ends_with("_2_5"),
    names_to = "model",
    values_to = "CI_25",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_2_5"
  ) %>% 
  pivot_longer(
    cols = ends_with("_97_5"),
    names_to = "model2",
    values_to = "CI_975",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_97_5"
  ) %>% 
  pivot_longer(
    cols = ends_with("_est"),
    names_to = "model3",
    values_to = "est",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_est"
  ) %>% 
  filter(model == model2 & model == model3) %>% 
  filter(type == "ATE") %>% 
  filter(model3 == "IPW")

















r %>% 
  pivot_longer(
    cols = ends_with("_2_5"),
    names_to = "model",
    values_to = "CI_25",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_2_5"
  ) %>% 
  pivot_longer(
    cols = ends_with("_97_5"),
    names_to = "model2",
    values_to = "CI_975",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_97_5"
  ) %>% 
  pivot_longer(
    #cols = 4:17,
    cols = ends_with("_est"),
    names_to = "model3",
    values_to = "est",
    # This regex says: capture everything before "_data"
    names_pattern = "(.*)_est"
  ) %>% 
  group_by(replicates, type, model3) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(type, model3) %>% 
  mutate(
    bias = (est - true_TE),
    rmse = (est - true_TE)^2
  ) %>% 
  filter(model3 != "IPW" & model3!= "sample_mean") %>% 
  #filter(model3 == "model_2_4_spline_logit") %>% 
  ggplot() + geom_violin(aes(model3, bias)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

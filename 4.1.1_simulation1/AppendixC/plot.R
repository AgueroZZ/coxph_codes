library(tidyverse)
TEXT_SIZE <- 25

load("result_diff_cens_low.rda")
load("result_diff_cens_med.rda")
load("result_diff_cens_high.rda")

result2 <- as_tibble(result2)
colnames(result2) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
agg_means_diff_low <- apply(result2, 2, mean)

as.data.frame(result2) %>% filter(frailty_cov_aghq >= 0) %>%
  mutate(Proposed = beta_mse_aghq, INLA = beta_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-10,0))

as.data.frame(result2) %>% filter(frailty_cov_aghq >= 0) %>% mutate(Proposed = frailty_mse_aghq,INLA = frailty_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-1.5,0.5))


result2_med <- as_tibble(result2_med)
colnames(result2_med) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
agg_means_diff_med <- apply(result2_med, 2, mean)

as.data.frame(result2_med) %>% filter(frailty_cov_aghq >= 0) %>%
  mutate(Proposed = beta_mse_aghq, INLA = beta_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-10,0))

as.data.frame(result2_med) %>% filter(frailty_cov_aghq >= 0) %>% mutate(Proposed = frailty_mse_aghq,INLA = frailty_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-1.5,0.5))



result2_high <- as_tibble(result2_high)
colnames(result2_high) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
agg_means_diff_high <- apply(result2_high, 2, mean)

as.data.frame(result2_high) %>% filter(frailty_cov_aghq >= 0) %>%
  mutate(Proposed = beta_mse_aghq, INLA = beta_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-10,0))

as.data.frame(result2_high) %>% filter(frailty_cov_aghq >= 0) %>% mutate(Proposed = frailty_mse_aghq,INLA = frailty_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-1.5,0.5))

agg_means_diff_low <- apply(result2, 2, mean)
agg_means_diff_med <- apply(result2_med, 2, mean)
agg_means_diff_high <- apply(result2_high, 2, mean)














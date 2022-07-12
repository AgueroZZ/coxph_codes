### Load:
library(tidyverse)
TEXT_SIZE = 25
load("aggresult5000.Rda")

### change names:
result <- as_tibble(aggresult5000)
result$m <- c(1,2,3,4,5,10)
names(result)


### point plot for beta MSE:
betaMSE <- tibble(MSE = c(result$beta_mse_aghq, result$beta_mse_inla), 
                  m = c(result$m, result$m),
                  Type = c(rep("Proposed", nrow(result)), rep("INLA", nrow(result))))

betaMSE %>% ggplot(aes(x = m, y = MSE, color = Type)) + geom_point() + 
  geom_line() +
  theme_classic(base_size = TEXT_SIZE) +
  scale_x_continuous(breaks=c(1:10)) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) + xlab("m (group size)")



### point plot for beta Cov:
betaCov <- tibble(Cov = c(result$beta_cov_aghq, result$beta_cov_inla), 
                  m = c(result$m, result$m),
                  Type = c(rep("Proposed", nrow(result)), rep("INLA", nrow(result))))

betaCov %>% ggplot(aes(x = as.factor(m), y = Cov, fill = Type)) + 
  geom_col(position = "dodge") + coord_cartesian(ylim=c(0.45,1)) +
  theme_classic(base_size = TEXT_SIZE) +
  geom_hline(yintercept = 0.95, color = "red") + 
  theme(legend.position = "none") + ylab("Coverage Prob") + 
  xlab("m (group size)")



### point plot for frailty MSE:
fraMSE <- tibble(MSE = c(result$frailty_mse_aghq, result$frailty_mse_inla), 
                  m = c(result$m, result$m),
                  Type = c(rep("Proposed", nrow(result)), rep("INLA", nrow(result))))

fraMSE %>% ggplot(aes(x = m, y = MSE, color = Type)) + geom_point() + 
  geom_line() +
  theme_classic(base_size = TEXT_SIZE) +
  scale_x_continuous(breaks=c(1:10)) +
  theme(legend.position = "none") +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) + xlab("m (group size)")



### point plot for frailty Cov:
fraCov <- tibble(Cov = c(result$frailty_cov_aghq, result$frailty_cov_inla), 
                  m = c(result$m, result$m),
                  Type = c(rep("Proposed", nrow(result)), rep("INLA", nrow(result))))

fraCov %>% ggplot(aes(x = as.factor(m), y = Cov, fill = Type)) + 
  geom_col(position = "dodge") + coord_cartesian(ylim=c(0.45,1)) +
  theme_classic(base_size = TEXT_SIZE) +
  geom_hline(yintercept = 0.95, color = "red") + 
  theme(legend.position = "none") + ylab("Coverage Prob") + xlab("m (group size)")








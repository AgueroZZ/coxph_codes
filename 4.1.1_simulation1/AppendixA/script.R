library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(rstan)
library(TMB)
library(INLA)
library(tmbstan)
library(foreach)
library(doMC)
library(parallel)
library(foreach)
library(abcoxph)
library(mvQuad)
library(doParallel)

TEXT_SIZE <- 25


compile("03_coxph_frailty.cpp")
dyn.load(dynlib("03_coxph_frailty"))


##### The new case with different group size:
beta = 0.2
sd = 1
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.5)))
Simulate_baseline3 <- function(timelim = 300, breaks = 0.001, cut = 5){
  timelim <- timelim
  tdom <- seq(0, timelim, by = breaks)
  haz <- rep(0, length(tdom))
  cut <- cut
  for (i in 1:cut) {
    low <- as.numeric(quantile(tdom,(i-1)/cut))
    high <- as.numeric(quantile(tdom,(i)/cut))
    if(i %% 2 == 1){
      a <- runif(1,0,1)
      if(a > 0.3) haz[tdom<=high & tdom > low] <- 0.1
      else {
        c <- tdom[tdom<=high & tdom > low]
        haz[tdom<=high & tdom > low] <-0.01
      }
    }
    if(i %% 2 == 0){
      a <- runif(1,0,1)
      if(a > 0.8){
        c <- tdom[tdom<=high & tdom > low]
        haz[tdom<=high & tdom > low] <- 0.25
      }
      else{
        haz[tdom<=high & tdom > low] <- sample(c(0.05,0.15),size = 1,prob = c(0.5,0.5))
      }
    }
  }
  baseline <- data.frame(time = tdom, hazard = haz, timelim = timelim)
}
### Try this baseline on the smoothing example as well
set.seed(1234)
baseline <- Simulate_baseline3()


Simulate_grouped_data_diff_size <- function(group_size_vec, bas = "piecewiseconstant", beta = 0.2, sdtheta = 0.8){
    n <- sum(group_size_vec) ### Total samples
    K <- length(group_size_vec) ### Number of groups
    if(bas == "piecewiseconstant") {
      timelim <- baseline$timelim[1]
      tdom <- baseline$time
      haz <- baseline$hazard
      true <- data.frame(time = tdom, hazard = haz)
      u <- rnorm(K, sd = sdtheta)
      u <- rep(u, times = group_size_vec)
      x <- rnorm(n, sd = 1)
      eta <- u + beta*x
      failtimes <- c()
      r <- runif(n)
      for (i in 1:n) {
        hazz <- haz * exp(eta[i])
        cumhaz <- cumsum(hazz*0.001)
        Surv <- exp(-cumhaz)
        Surv[1] <- 1
        failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
      }
      data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
      for (i in 1:length(data$censoring)) {
        if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
      }
      data$group <- rep(1:K, times = group_size_vec)
      data$true <- u
    }
    else if (bas == "regular") {
      timelim <- 200
      tdom <- seq(0, timelim, by = 0.001)
      haz <- rep(0, length(tdom))
      haz <- 0.2 * cos(0.15*tdom) + 0.3
      true <- data.frame(time = tdom, hazard = haz)
      u <- rnorm(K, sd = sdtheta)
      u <- rep(u, times = group_size_vec)
      x <- rnorm(n, sd = 3)
      eta <- u + beta*x
      failtimes <- c()
      r <- runif(n)
      for (i in 1:n) {
        hazz <- haz * exp(eta[i])
        cumhaz <- cumsum(hazz*0.001)
        Surv <- exp(-cumhaz)
        Surv[1] <- 1
        failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
      }
      data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
      for (i in 1:length(data$censoring)) {
        if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
      }
      data$group <- rep(1:K, times = group_size_vec)
      data$true <- u
    }
    else if (bas == "constant") {
      timelim <- 200
      tdom <- seq(0, timelim, by = 0.001)
      haz <- rep(0.05, length(tdom))
      true <- data.frame(time = tdom, hazard = haz)
      u <- rnorm(K, sd = sdtheta)
      u <- rep(u, times = group_size_vec)
      x <- rnorm(n, sd = 3)
      eta <- u + beta*x
      failtimes <- c()
      r <- runif(n)
      for (i in 1:n) {
        hazz <- haz * exp(eta[i])
        cumhaz <- cumsum(hazz*0.001)
        Surv <- exp(-cumhaz)
        Surv[1] <- 1
        failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
      }
      data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
      for (i in 1:length(data$censoring)) {
        if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
      }
      data$group <- rep(1:K, times = group_size_vec)
      data$true <- u
    }
    else{
      timelim <- 200
      tdom <- seq(0, timelim, by = 0.001)
      haz <- rep(0, length(tdom))
      cut <- 40
      for (i in 1:cut) {
        low <- as.numeric(quantile(tdom,(i-1)/cut))
        high <- as.numeric(quantile(tdom,(i)/cut))
        if(i %% 2 == 1){
          haz[tdom<=high & tdom > low] <- 0.01
        }
        else if(i %% 2 == 0){
          haz[tdom<=high & tdom > low] <- 0.25
        }
      }
      true <- data.frame(time = tdom, hazard = haz)
      u <- rnorm(K, sd = sdtheta)
      u <- rep(u, times = group_size_vec)
      x <- rnorm(n, sd = 3)
      eta <- u + beta*x
      failtimes <- c()
      r <- runif(n)
      for (i in 1:n) {
        hazz <- haz * exp(eta[i])
        cumhaz <- cumsum(hazz*0.001)
        Surv <- exp(-cumhaz)
        Surv[1] <- 1
        failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))] 
      }
      data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
      for (i in 1:length(data$censoring)) {
        if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
      }
      data$group <- rep(1:K, times = group_size_vec)
      data$true <- u
    }
    data
  }

do_once_diff_size <- function(seed, beta, group_size_vec, sdtheta, bas = "constant"){
  set.seed(seed)
  n <- sum(group_size_vec)
  K <- length(group_size_vec)
  data <- Simulate_grouped_data_diff_size(group_size_vec = group_size_vec,
                                          beta = beta, sdtheta = sdtheta, bas = "piecewiseconstant")
  data <- abcoxph:::arrange_data(data)
  dat <- tibble(x = data$x, t = data$times, cens = data$censoring, group = data$group)
  dat$ranks <- rank(dat$t, ties.method = "min")
  X <- as(as.matrix(dat$x),"dgTMatrix")
  B <- as(abcoxph:::create_blist_element(u = dat$group)$B,"dgTMatrix")
  D <- as(abcoxph:::create_diff_matrix(n), "dgTMatrix") ### n = K * N
  ### Setup TMB:
  tmbdat <- list(
    # Design matrix (random and fixed)
    B = as(B,"dgTMatrix"),
    X = as(X,"dgTMatrix"),
    # Differencing matrix
    D = as(D,"dgTMatrix"),
    # Response
    ranks = as.integer(dat$ranks),
    cens = as.integer(dat$cens),
    # Prior params
    u = 1,
    alpha = 0.5,
    betaprec = 0.001)
  
  tmbparams <- list(
    W = rep(0,ncol(B)+ncol(X)),
    theta = 0 # -2log(sigma)
  )
  ##### Fitting:
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "03_coxph_frailty",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  # AGHQ
  quad <- tryCatch(aghq::marginal_laplace_tmb(ff,15,0), error=function(err) 0)
  if(!is.list(quad)){
    return(rep(-1,8))
  }
  samps <- sample_marginal(quad,10000,interpolation = 'spline')
  beta_est <- samps$samps[(K+1),]
  post_means <- apply(samps$samps, 1, mean)
  post_sds <- apply(samps$samps, 1, sd)
  post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
  post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
  post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
  rownames(post_sum_aghq) <- c("mean", "sd", "upper", "lower")
  colnames(post_sum_aghq) <- c(c(1:K),"beta")
  beta_cov_aghq <- ifelse(beta <= post_sum_aghq[,K+1][3] & beta >= post_sum_aghq[,K+1][4], 1, 0)
  beta_mse_aghq <- (post_sum_aghq[,K+1][1] - beta)^2
  ### INLA:
  prior.prec <- list(prec = list(prior = "pc.prec",
                                 param = c(tmbdat$u, tmbdat$a)))
  formula <- inla.surv(t,cens)~ x + f(group, model = "iid", hyper = prior.prec)

  Inlaresult <- tryCatch(inla(formula = formula, control.fixed = list(prec = 0.001), data = dat, family = "coxph"), error=function(err) 0)
  if(!is.list(Inlaresult)){
    return(rep(-1,8))
  }
  
  beta_cov_inla <- ifelse(beta <= Inlaresult$summary.fixed[2,]$'0.975quant' & beta >= Inlaresult$summary.fixed[2,]$'0.025quant', 1, 0)
  beta_mse_inla <- (Inlaresult$summary.fixed[2,]$mean - beta)^2
  frailty <- data %>% select(c(group,true)) %>% arrange(group) %>% unique(by = group)
  frailty$AGHQ <- as.numeric(post_sum_aghq[1,][-(K+1)])
  frailty$INLA <- as.numeric(Inlaresult$summary.random$group$mean)
  ### Random effects: MSE
  frailty_mse_aghq <- mean((frailty$true - frailty$AGHQ)^2)
  frailty_mse_inla <- mean((frailty$true - frailty$INLA)^2)
  ### Random effects: coverage
  frailty$AGHQ_up <- as.numeric(post_sum_aghq[3,][-(K+1)])
  frailty$INLA_up <- as.numeric(Inlaresult$summary.random$group$`0.975quant`)
  frailty$AGHQ_lo <- as.numeric(post_sum_aghq[4,][-(K+1)])
  frailty$INLA_lo <- as.numeric(Inlaresult$summary.random$group$`0.025quant`)
  ### AGHQ:
  frailty_cov_aghq <- mean((frailty$true >= frailty$AGHQ_lo & frailty$true <= frailty$AGHQ_up))
  ### INLA:
  frailty_cov_inla <- mean((frailty$true >= frailty$INLA_lo & frailty$true <= frailty$INLA_up))
  result <- c(beta_cov_aghq,beta_cov_inla,beta_mse_aghq,beta_mse_inla,frailty_cov_aghq,frailty_cov_inla,frailty_mse_aghq,frailty_mse_inla)
  result
}

######### Implementation:
group_size_vec <- c(rep(1,30),rep(2,25),rep(3,5))

M <- 5000
time_begin <- Sys.time()
resultdiff <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once_diff_size(seed = i, beta = beta, group_size_vec = group_size_vec, sdtheta = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
colnames(resultdiff) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
save(resultdiff, file = "resultdiff5000_case1.Rda")

agg_means_diff <- apply(resultdiff, 2, mean)

as.data.frame(resultdiff) %>% filter(frailty_cov_aghq >= 0) %>%
  mutate(Proposed = beta_mse_aghq, INLA = beta_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-10,0))

ggsave("beta_mse_diff_size.png", width = 5, height = 5)

as.data.frame(resultdiff) %>% filter(frailty_cov_aghq >= 0) %>% mutate(Proposed = frailty_mse_aghq,INLA = frailty_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-1,0.5))

ggsave("fra_mse_diff_size.png", width = 5, height = 5)


### Also a sparse case:
group_size_vec <- c(rep(1,30),rep(2,28),rep(6,2))
time_begin <- Sys.time()
resultdiff2 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once_diff_size(seed = i, beta = beta, group_size_vec = group_size_vec, sdtheta = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
colnames(resultdiff2) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
save(resultdiff2, file = "resultdiff5000_case2.Rda")



agg_means_diff2 <- apply(resultdiff2, 2, mean)

as.data.frame(resultdiff2) %>% filter(frailty_cov_aghq >= 0) %>%
  mutate(Proposed = beta_mse_aghq, INLA = beta_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-10,0))
ggsave("beta_mse_diff_size2.png", width = 5, height = 5)

as.data.frame(resultdiff2) %>% filter(frailty_cov_aghq >= 0) %>% mutate(Proposed = frailty_mse_aghq,INLA = frailty_mse_inla) %>% 
  pivot_longer(cols = Proposed:INLA, names_to = "method", values_to = "MSE") %>% 
  ggplot() + geom_boxplot(aes(x = method, y = log(MSE), fill = method), outlier.shape = NA) +
  theme_classic(base_size = TEXT_SIZE) + 
  theme(legend.position = "none") +
  xlab("") + ylim(c(-1.5,0.5))
ggsave("fra_mse_diff_size2.png", width = 5, height = 5)


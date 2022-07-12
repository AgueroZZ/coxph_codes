## Coxph regression with frailty ###
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


precompile()
TEXT_SIZE = 25
ncores = 2
registerDoMC(ncores)
# workers = makeCluster(ncores, type="SOCK")
# registerDoParallel(workers)

### Simulating function:
K = 60
beta = 0.2
sd = 1
M <- 5000
N <- 2

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

Simulate_grouped_data <- function(N = 2, bas = "piecewiseconstant", K = 50, beta = 0.2, sdtheta = 0.8, cens = 0.1){
  n <- N*K ### Total samples
  if(bas == "piecewiseconstant") {
    timelim <- baseline$timelim[1]
    tdom <- baseline$time
    haz <- baseline$hazard
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
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
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = (1-cens))}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  else if (bas == "regular") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz <- 0.2 * cos(0.15*tdom) + 0.3
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
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
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = (1-cens))}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  else if (bas == "constant") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0.05, length(tdom))
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
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
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = (1-cens))}
    }
    data$group <- rep(1:K, each = N)
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
    u <- rep(u, each = N)
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
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = (1-cens))}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  data
}

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.5)))



# TMB function template
compile("03_coxph_frailty.cpp")
dyn.load(dynlib("03_coxph_frailty"))



######### Speeding up the simulation function:
do_once <- function(seed,beta, N, K, sd, bas = "constant", cens){
  set.seed(seed)
  n <- K*N
  data <- Simulate_grouped_data(N = N, bas = bas, K = K, beta = beta, sdtheta = sd, cens = cens)
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
  ### Check if the model works:
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
  ### Check if the model works:
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


###### Step Hazard function:

### Sigma = 1
sd = 1

### N = 2 (cens rate of 20 percent)
time_begin <- Sys.time()
result2_med <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 2, K = K, sd = sd, bas = "piecewiseconstant", cens = 0.2)
time_end <- Sys.time()
time_end - time_begin
agg_means2_med <- apply(result2_med, 2, mean)
save(file = "result_diff_cens_med.rda", result2_med)

### N = 2 (cens rate of 40 percent)
time_begin <- Sys.time()
result2_high <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 2, K = K, sd = sd, bas = "piecewiseconstant", cens = 0.4)
time_end <- Sys.time()
time_end - time_begin
agg_means2_high <- apply(result2_high, 2, mean)
save(file = "result_diff_cens_high.rda", result2_high)






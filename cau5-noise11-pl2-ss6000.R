rm(list = ls())

library(BWMR)
library(TwoSampleMR)
source("gsmr_md.R")
library(ggplot2)

library(MASS)  # Needed when sample from multivariate normal distribution


sim_null_4GM_select_sd <- function(P, N_1, N_2, Pi_11, Pi_10, Pi_01, Pi_00, beta, rho, noise, threshold){
  
  # genotype data
  # X_1: Individual-SNP matrix for exposure
  X_1 <- matrix(1, N_1, P)
  X_1f <- matrix(runif(N_1*P), N_1, P)
  # X_2: Individual-SNP matrix for outcome
  X_2 <- matrix(1, N_2, P)
  X_2f <- matrix(runif(N_2*P), N_2, P)
  # X_1: Individual-SNP matrix for exposure (MR)
  X_3 <- matrix(1, N_1, P)
  X_3f <- matrix(runif(N_1*P), N_1, P)
  # f: minor allele frequencies of P SNPs
  f <- runif(P, 0.05, 0.5)
  f2_1 <- matrix(rep(f^2, each = N_1), N_1, P)
  f2_2 <- matrix(rep(f^2, each = N_2), N_2, P)
  f0_1 <- matrix(rep((1 - f)^2, each = N_1), N_1, P)
  f0_2 <- matrix(rep((1 - f)^2, each = N_2), N_2, P)
  # assign 2
  X_1[X_1f < f2_1]  <- 2
  X_2[X_2f < f2_2]  <- 2
  X_2[X_3f < f2_1]  <- 2
  # assign 0
  X_1[(1 - X_1f) < f0_1] <- 0
  X_2[(1 - X_2f) < f0_2] <- 0
  X_3[(1 - X_3f) < f0_1] <- 0
  
  # Sigma: covariance matrix for (gamma, Gamma), under Pi_11
  s1 <- 1 / (P*(Pi_11+Pi_10))
  s2 <- 1 / (P*(Pi_11+Pi_01))
  cov <- rho*sqrt(s1*s2)
  Sigma <- matrix(c(s1, cov, cov, s2), 2, 2)
  
  # W: SNP-(gamma, Gamma) matrix (ground truth) 
  W <- matrix(0, P, 2)
  W[1:(P*Pi_11), ] <- mvrnorm(P*Pi_11, rep(0, 2), Sigma)
  W[(P*Pi_11+1):(P*(Pi_11+Pi_10)), 1] <- rnorm(P*Pi_10, 0, sqrt(1/(P*(Pi_11+Pi_10))))
  W[(P*(Pi_11+Pi_10)+1):(P*(Pi_11+Pi_10+Pi_01)), 2] <- rnorm(P*Pi_01, 0, sqrt(1/(P*(Pi_11+Pi_01))))
  W[, 2] <- beta*W[, 1] + W[, 2]
  
  # Exposure data
  y0 <- X_1%*%W[, 1]
  y <- y0 + noise*sqrt(drop(var(y0)))*rnorm(N_1)
  # Oucome data
  z0 <- X_2%*%W[, 2]
  z <- z0 + noise*sqrt(drop(var(z0)))*rnorm(N_2)
  # Exposure data (MR)
  y03 <- X_3%*%W[, 1]
  y3 <- y03 + noise*sqrt(drop(var(y03)))*rnorm(N_1)
  
  # b.exposure, se.exposure
  b.exposure <- rep(0, P)
  se.exposure <- rep(0, P)
  pval.exposure <- rep(0, P)
  for(j in 1:P){
    fit_lm <- lm(y~., data = data.frame(y, X_1[, j]))
    b.exposure[j] <- summary(fit_lm)$coefficients[2, 1]
    se.exposure[j] <- summary(fit_lm)$coefficients[2, 2]
    pval.exposure[j] <- summary(fit_lm)$coefficients[2, 4]
  }
  # b.outcome, se.outcome
  b.outcome <- rep(0, P)
  se.outcome <- rep(0, P)
  pval.outcome <- rep(0, P)
  for(j in 1:P){
    fit_lm <- lm(z~., data = data.frame(z, X_2[, j]))
    b.outcome[j] <- summary(fit_lm)$coefficients[2, 1]
    se.outcome[j] <- summary(fit_lm)$coefficients[2, 2]
    pval.outcome[j] <- summary(fit_lm)$coefficients[2, 4]
  }
  # b.exposure, se.exposure (MR)
  b.exposure.mr <- rep(0, P)
  se.exposure.mr <- rep(0, P)
  pval.exposure.mr <- rep(0, P)
  for(j in 1:P){
    fit_lm <- lm(y3~., data = data.frame(y3, X_3[, j]))
    b.exposure.mr[j] <- summary(fit_lm)$coefficients[2, 1]
    se.exposure.mr[j] <- summary(fit_lm)$coefficients[2, 2]
    pval.exposure.mr[j] <- summary(fit_lm)$coefficients[2, 4]
  }
  
  res_bias <- list()
  res_unbias <- list()
  for (j in 1:length(threshold)) {
    indx <- which(pval.exposure < threshold[j])
    res_bias[[j]] <- data.frame(b.exposure = b.exposure[indx], b.outcome = b.outcome[indx],
                                se.exposure = se.exposure[indx], se.outcome = se.outcome[indx])
    res_unbias[[j]] <- data.frame(b.exposure = b.exposure.mr[indx], b.outcome = b.outcome[indx],
                                  se.exposure = se.exposure.mr[indx], se.outcome = se.outcome[indx])
  }
  # res_all <- data.frame(b.exposure=b.exposure.mr, b.outcome=b.outcome,
  #                       se.exposure=se.exposure.mr, se.outcome=se.outcome)
  res <- list('unbias'=res_unbias, 'bias'=res_bias, 'index'=indx)
}


Rp <- 100
L <- 3

pval.bwmr.unbias <- matrix(nrow = L, ncol = Rp)
pval.raps.unbias <- matrix(nrow = L, ncol = Rp)
pval.egger.unbias <- matrix(nrow = L, ncol = Rp)
pval.gsmr.unbias <- matrix(nrow = L, ncol = Rp)
pval.gsmr.F.unbias <- matrix(nrow = L, ncol = Rp)

b.bwmr.unbias <- matrix(nrow = L, ncol = Rp) 
b.raps.unbias <- matrix(nrow = L, ncol = Rp)
b.egger.unbias <- matrix(nrow = L, ncol = Rp)
b.gsmr.unbias <- matrix(nrow = L, ncol = Rp)
b.gsmr.F.unbias <- matrix(nrow = L, ncol = Rp)

pval.bwmr.bias <- matrix(nrow = L, ncol = Rp)
pval.raps.bias <- matrix(nrow = L, ncol = Rp)
pval.egger.bias <- matrix(nrow = L, ncol = Rp)
pval.gsmr.bias <- matrix(nrow = L, ncol = Rp)
pval.gsmr.F.bias <- matrix(nrow = L, ncol = Rp)

b.bwmr.bias <- matrix(nrow = L, ncol = Rp) 
b.raps.bias <- matrix(nrow = L, ncol = Rp)
b.egger.bias <- matrix(nrow = L, ncol = Rp)
b.gsmr.bias <- matrix(nrow = L, ncol = Rp)
b.gsmr.F.bias <- matrix(nrow = L, ncol = Rp)


for (i in 1:Rp) {
  message(i)
  set.seed(i)
  
  res <- sim_null_4GM_select_sd(P = 10000, N_1 = 6000, N_2 = 6000,
                                Pi_11 = 0.02, Pi_10 = 0.08, Pi_01 = 0.08, Pi_00 = 0.82,
                                beta = 0.5, rho = 0, noise = sqrt(1/1), threshold = c(1e-5, 5*1e-5, 1e-4))
  
  for (k in 1:3){ # k thresholds
    ## unbias
    dat <- res$unbias[[k]]
    N <- dim(dat)[1]
    message("N=", N)
    
    # BWMR
    res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                         sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
    if (inherits(res_bwmr, "try-error")) {next}
    b.bwmr.unbias[k, i] <- res_bwmr$beta
    pval.bwmr.unbias[k, i] <- res_bwmr$P_value
    
    # GSMR
    res_gsmr <-  try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(N), 
                          heidi_outlier_flag = T,
                          heidi_outlier_thresh=0.01,
                          nsnps_thresh = 1))
    if (inherits(res_gsmr, "try-error")) {next}
    b.gsmr.unbias[k, i] <- drop(res_gsmr$bxy)
    pval.gsmr.unbias[k, i] <- drop(res_gsmr$bxy_pval)
    
    # GSMR HEIDI False
    res_gsmr <-  try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(N), 
                          heidi_outlier_flag = F,
                          heidi_outlier_thresh=0.01,
                          nsnps_thresh = 1))
    if (inherits(res_gsmr, "try-error")) {next}
    b.gsmr.F.unbias[k, i] <- drop(res_gsmr$bxy)
    pval.gsmr.F.unbias[k, i] <- drop(res_gsmr$bxy_pval)
    
    ## Two sample MR
    dat$id.exposure <- rep("uVXQCX", nrow(dat))
    dat$id.outcome <- rep("uVXQCX", nrow(dat))
    dat$mr_keep <- rep(TRUE, nrow(dat))
    dat$exposure <- rep('exposure', nrow(dat))
    dat$outcome <- rep('outcome', nrow(dat))
    colnames(dat) <-  c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
    result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
    if (inherits(result, "try-error")) {next}
    b.raps.unbias[k, i] <- result$b[1]
    pval.raps.unbias[k, i] <- result$pval[1]
    b.egger.unbias[k, i] <- result$b[2]
    pval.egger.unbias[k, i] <- result$pval[2]
    
    
    ## bias
    dat <- res$bias[[k]]
    N <- dim(dat)[1]
    # BWMR
    res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                         sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
    if (inherits(res_bwmr, "try-error")) {next}
    b.bwmr.bias[k, i] <- res_bwmr$beta
    pval.bwmr.bias[k, i] <- res_bwmr$P_value
    
    # GSMR
    res_gsmr <-  try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(N), 
                          heidi_outlier_flag = T,
                          heidi_outlier_thresh=0.01,
                          nsnps_thresh = 1))
    if (inherits(res_gsmr, "try-error")) {next}
    b.gsmr.bias[k, i] <- drop(res_gsmr$bxy)
    pval.gsmr.bias[k, i] <- drop(res_gsmr$bxy_pval)
    
    # GSMR HEIDI False
    res_gsmr <- try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                         bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                         ldrho = diag(N), 
                         heidi_outlier_flag = F,
                         heidi_outlier_thresh=0.01,
                         nsnps_thresh = 1))
    if (inherits(res_gsmr, "try-error")) {next}
    b.gsmr.F.bias[k, i] <- drop(res_gsmr$bxy)
    pval.gsmr.F.bias[k, i] <- drop(res_gsmr$bxy_pval)
    
    ## Two sample MR
    dat$id.exposure <- rep("uVXQCX", nrow(dat))
    dat$id.outcome <- rep("uVXQCX", nrow(dat))
    dat$mr_keep <- rep(TRUE, nrow(dat))
    dat$exposure <- rep('exposure', nrow(dat))
    dat$outcome <- rep('outcome', nrow(dat))
    colnames(dat) <-  c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
    result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
    if (inherits(result, "try-error")) {next}
    b.raps.bias[k, i] <- result$b[1]
    pval.raps.bias[k, i] <- result$pval[1]
    b.egger.bias[k, i] <- result$b[2]
    pval.egger.bias[k, i] <- result$pval[2]
  }
}

save.image("cau5-noise11-pl2-ss6000.RData")

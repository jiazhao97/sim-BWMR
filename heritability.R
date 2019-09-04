rm(list = ls())

set.seed(1)

source("simu_sum_func.R")
source("rnormalmix_func.R")
library(BWMR)
library(ggplot2)
library(TwoSampleMR)
source("gsmr_md.R")

Rp <- 200
beta_choice <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
l <- length(beta_choice)

bset.bwmr <- matrix(NA, nrow = l, ncol = Rp)
pvalset.bwmr <- matrix(NA, nrow = l, ncol = Rp)
bset.raps <- matrix(NA, nrow = l, ncol = Rp)
pvalset.raps <- matrix(NA, nrow = l, ncol = Rp)
bset.gsmr <- matrix(NA, nrow = l, ncol = Rp)
pvalset.gsmr <- matrix(NA, nrow = l, ncol = Rp)
bset.egger <- matrix(NA, nrow = l, ncol = Rp)
pvalset.egger <- matrix(NA, nrow = l, ncol = Rp)

for (i in 1:l) {
  beta <- beta_choice[i]
  for (rp in 1:Rp) {
    dat <- simu_sum(p = 100000, 
                    n1 = 50000, 
                    n2 = 50000, 
                    h1 = 0.2, # 0.3, 0.4
                    h2 = 0.2, # 0.3, 0.4
                    neffect1 = 5000, 
                    neffect2 = 5000,
                    beta = beta,
                    pval_thresh = 5e-7)
    
    # BWMR
    res <- BWMR(gammahat = dat$beta.exposure, Gammahat = dat$beta.outcome,
                sigmaX = dat$se.exposure, sigmaY = dat$se.outcome)
    message(length(dat$beta.exposure))
    bset.bwmr[i, rp] <- res$beta
    pvalset.bwmr[i, rp] <- res$P_value
    
    # GSMR
    res_gsmr <-  try(gsmr(bzx = dat$beta.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$beta.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(length(dat$beta.exposure)), 
                          heidi_outlier_flag = T,
                          heidi_outlier_thresh=0.01,
                          nsnps_thresh = 1))
    if (inherits(res_gsmr, "try-error")) {next}
    bset.gsmr[i, rp] <- drop(res_gsmr$bxy)
    pvalset.gsmr[i, rp] <- drop(res_gsmr$bxy_pval)
    
    ## Two sample MR
    dat$id.exposure <- rep("uVXQCX", nrow(dat))
    dat$id.outcome <- rep("uVXQCX", nrow(dat))
    dat$mr_keep <- rep(TRUE, nrow(dat))
    dat$exposure <- rep('exposure', nrow(dat))
    dat$outcome <- rep('outcome', nrow(dat))
    colnames(dat) <-  c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
    result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
    if (inherits(result, "try-error")) {next}
    bset.raps[i, rp] <- result$b[1]
    pvalset.raps[i, rp] <- result$pval[1]
    bset.egger[i, rp] <- result$b[2]
    pvalset.egger[i, rp] <- result$pval[2]
  }
}

save.image("heritability-02.RData")
# save.image("heritability-03.RData")
# save.image("heritability-04.RData")

# Rp <- length(pvalset.bwmr[i, ])
# type.I.error <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
#                   "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
# type.I.error

##### Case-1 in summary-level simulations (estimates) #####
# Figures S3 and S4 in Supplementary Document


rm(list = ls())

library(BWMR)
library(TwoSampleMR)
source("gsmr_md.R")
library(ggplot2)

sim_without_outlier <- function(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true){
  gamma_true <- rnorm(N, 0, sqrt(sqsigma_true))
  gammahat <- rnorm(N, gamma_true, sigmaX_true)
  Gammahat <- rnorm(N, beta_true*gamma_true, sqrt(sigmaY_true^2+sqtau_true))
  df_data <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX_true,
    sigmaY = sigmaY_true
  )
}

Rp <- 50 # Repeat 50 times
tau_choice <- c(0.1, 0.2, 0.3)
sigma_choice <- c(0.6, 0.8, 1.0)
beta_choice <- c(0.0, 0.2, 0.5)
L <- length(beta_choice)
L1 <- length(tau_choice)
L2 <- length(sigma_choice)

sigmaX_min <- 0.3
sigmaX_max <- 0.5
sigmaY_min <- 0.3
sigmaY_max <- 0.5
N <- 50 # N=50 SNPs


### change tau
sqsigma_true <- 0.8^2 # fix sigma=0.8
# declare sets
b.bwmr.tau <- matrix(nrow = L, ncol = Rp*L1)
b.raps.tau <- matrix(nrow = L, ncol = Rp*L1)
b.gsmr.tau <- matrix(nrow = L, ncol = Rp*L1)
b.egger.tau <- matrix(nrow = L, ncol = Rp*L1)
for (k in 1:L) {
  beta_true <- beta_choice[k]
  i <- 0
  for (j in 1:L1) {
    sqtau_true <- tau_choice[j]^2
    for (iter in 1:Rp) {
      i <- i+1
      set.seed(iter)
      sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
      sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
      dat <- sim_without_outlier(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true)
      colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
      # BWMR
      res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                           sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
      b.bwmr.tau[k, i] <- res_bwmr$beta
      # GSMR
      res_gsmr = try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(N), 
                          heidi_outlier_flag = T, 
                          nsnps_thresh = 1))
      b.gsmr.tau[k, i] <- drop(res_gsmr$bxy)
      ## Two sample MR
      dat$id.exposure <- rep("uVXQCX", nrow(dat))
      dat$id.outcome <- rep("uVXQCX", nrow(dat))
      dat$mr_keep <- rep(TRUE, nrow(dat))
      dat$exposure <- rep('exposure', nrow(dat))
      dat$outcome <- rep('outcome', nrow(dat))
      colnames(dat) = c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
      result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
      b.raps.tau[k, i] <- result$b[1]
      b.egger.tau[k, i] <- result$b[2]
    }
  }
}


### change sigma
sqtau_true <- 0.3^2 # fix tau=0.3
# declare sets
b.bwmr.sigma <- matrix(nrow = L, ncol = Rp*L2)
b.raps.sigma <- matrix(nrow = L, ncol = Rp*L2)
b.gsmr.sigma <- matrix(nrow = L, ncol = Rp*L2)
b.egger.sigma <- matrix(nrow = L, ncol = Rp*L2)
for (k in 1:L) {
  beta_true <- beta_choice[k]
  i <- 0
  for (j in 1:L2) {
    sqsigma_true <- sigma_choice[j]^2
    for (iter in 1:Rp) {
      i <- i+1
      set.seed(iter)
      sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
      sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
      dat <- sim_without_outlier(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true)
      colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
      # BWMR
      res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                           sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
      b.bwmr.sigma[k, i] <- res_bwmr$beta
      # GSMR
      res_gsmr = try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                          ldrho = diag(N), 
                          heidi_outlier_flag = T, 
                          nsnps_thresh = 1))
      b.gsmr.sigma[k, i] <- drop(res_gsmr$bxy)
      ## Two sample MR
      dat$id.exposure <- rep("uVXQCX", nrow(dat))
      dat$id.outcome <- rep("uVXQCX", nrow(dat))
      dat$mr_keep <- rep(TRUE, nrow(dat))
      dat$exposure <- rep('exposure', nrow(dat))
      dat$outcome <- rep('outcome', nrow(dat))
      colnames(dat) = c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
      result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
      b.raps.sigma[k, i] <- result$b[1]
      b.egger.sigma[k, i] <- result$b[2]
    }
  }
}


a <- ls()
rm(list=a[which(a!='b.bwmr.tau' & a !='b.egger.tau' & a !='b.gsmr.tau' & a !='b.raps.tau' & a!='b.bwmr.sigma' & a !='b.egger.sigma' & a !='b.gsmr.sigma' & a !='b.raps.sigma')])
save.image("case-1-estimate.RData")


#-----------------------------------------------------------------------------------------------------------
### Plot
rm(list = ls())
library(ggplot2)
library(reshape2)

method_0 <- c("BWMR", "RAPS", "GSMR", "Egger")
method_level <- c("BWMR", "Egger", "GSMR", "RAPS")
beta_0 <- c("beta=0.0", "beta=0.2", "beta=0.5")
beta_num <- c(0.0, 0.2, 0.5)
tau_num <- c(0.1, 0.2, 0.3)
sig_num <- c(0.6, 0.8, 1.0)

# load data
load("case-1-estimate.RData")

## Boxplot (chang tau)
# df
df_boxplot_est <- data.frame(
  beta.num = rep(beta_num, rep(150, 3)),
  beta.title = rep(beta_0, rep(150, 3)),
  tau.num = factor(rep(rep(tau_num, rep(50, 3)), 3)),
  BWMR = c(b.bwmr.tau[1, ], b.bwmr.tau[2, ], b.bwmr.tau[3, ]),
  RAPS = c(b.raps.tau[1, ], b.raps.tau[2, ], b.raps.tau[3, ]),
  GSMR = c(b.gsmr.tau[1, ], b.gsmr.tau[2, ], b.gsmr.tau[3, ]),
  Egger = c(b.egger.tau[1, ], b.egger.tau[2, ], b.egger.tau[3, ])
)
df_boxplot_est <- melt(df_boxplot_est, id.vars = c("beta.num", "beta.title", "tau.num"),
                       variable.name = 'Method', value.name = 'beta.est')
df_boxplot_est$Method <- factor(df_boxplot_est$Method, levels = method_level)
# plt
plt_boxplot_est <- ggplot(df_boxplot_est, aes(x = tau.num, y = beta.est)) + 
  geom_boxplot(aes(fill = Method)) + 
  geom_hline(aes(yintercept = beta.num), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "True value of parameter tau", title = "Estimation") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plt_boxplot_est <- plt_boxplot_est + facet_grid(~ beta.title) +
  theme(strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))
plt_boxplot_est 

## Boxplot (chang sigma)
df_boxplot_est <- data.frame(
  beta.num = rep(beta_num, rep(150, 3)),
  beta.title = rep(beta_0, rep(150, 3)),
  sig.num = factor(rep(rep(sig_num, rep(50, 3)), 3)),
  BWMR = c(b.bwmr.sigma[1, ], b.bwmr.sigma[2, ], b.bwmr.sigma[3, ]),
  RAPS = c(b.raps.sigma[1, ], b.raps.sigma[2, ], b.raps.sigma[3, ]),
  GSMR = c(b.gsmr.sigma[1, ], b.gsmr.sigma[2, ], b.gsmr.sigma[3, ]),
  Egger = c(b.egger.sigma[1, ], b.egger.sigma[2, ], b.egger.sigma[3, ])
)
df_boxplot_est <- melt(df_boxplot_est, id.vars = c("beta.num", "beta.title", "sig.num"),
                       variable.name = 'Method', value.name = 'beta.est')
df_boxplot_est$Method <- factor(df_boxplot_est$Method, levels = method_level)
plt_boxplot_est <- ggplot(df_boxplot_est, aes(x = sig.num, y = beta.est)) + 
  geom_boxplot(aes(fill = Method)) + 
  geom_hline(aes(yintercept = beta.num), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "True value of parameter sigma", y = "Estimation", title = "Estimation") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plt_boxplot_est <- plt_boxplot_est + facet_grid(~ beta.title) +
  theme(strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))
plt_boxplot_est 


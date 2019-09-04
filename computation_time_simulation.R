rm(list = ls())

setwd("/Users/jiazhao/Documents/HKUST/201809 BWMR/BWMR-bioinformatics/major_revision/code-RData/computation_time")
library(BWMR)
library(TwoSampleMR)
library(MRPRESSO)
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

Rp <- 20 # Repeat 20 times
sigmaX_min <- 0.1
sigmaX_max <- 0.2
sigmaY_min <- 0.1
sigmaY_max <- 0.2
N_choice <- c(10, 20, 30, 40, 50) # N=50 SNPs

sqsigma_true <- 0.8^2 # fix sigma=0.8
sqtau_true <- 0.2^2
beta_true <- 0.0

recordtime_mat_bwmr <- matrix(nrow = length(N_choice), ncol = Rp)
recordtime_mat_raps <- matrix(nrow = length(N_choice), ncol = Rp)
recordtime_mat_egger <- matrix(nrow = length(N_choice), ncol = Rp)
recordtime_mat_gsmr <- matrix(nrow = length(N_choice), ncol = Rp)
recordtime_mat_presso <- matrix(nrow = length(N_choice), ncol = Rp)

for (i in 1:length(N_choice)) {
  N <- N_choice[i]
  for (iter in 1:Rp) {
    set.seed(iter)
    sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
    sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
    dat <- sim_without_outlier(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true)
    if (sum(is.na(dat)) > 0) {
      print("data error!")
      break
    }
    colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
    
    # BWMR
    recordtime_mat_bwmr[i, iter] <- system.time(try(BWMR_quick(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                                                             sigmaX = dat$se.exposure, sigmaY = dat$se.outcome)))[3]
    # GSMR
    recordtime_mat_gsmr[i, iter] <-  system.time(try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                                                          bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                                                          ldrho = diag(N), 
                                                          heidi_outlier_flag = T, 
                                                          nsnps_thresh = 1)))[3]
    ## Two sample MR
    dat$id.exposure <- rep("uVXQCX", nrow(dat))
    dat$id.outcome <- rep("uVXQCX", nrow(dat))
    dat$mr_keep <- rep(TRUE, nrow(dat))
    dat$exposure <- rep('exposure', nrow(dat))
    dat$outcome <- rep('outcome', nrow(dat))
    colnames(dat) = c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
    # RAPS
    recordtime_mat_raps[i, iter] <- system.time(try(mr(dat, method_list=c("mr_raps"))))[3]
    # Egger
    recordtime_mat_egger[i, iter] <- system.time(try(mr(dat, method_list=c("mr_egger_regression"))))[3]
    
    # MRPRESSO
    # Run MR-PRESSO global method
    recordtime_mat_presso[i, iter] <- system.time(try(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)))[3]
  }
}


a <- ls()
rm(list=a[which(a!='recordtime_mat_bwmr' & a !='recordtime_mat_raps' & a !='recordtime_mat_gsmr' & a !='recordtime_mat_egger' & a !='recordtime_mat_presso')])
save.image("computation_time_N.RData")


#------------------------------------------
rm(list = ls())
load("computation_time_N.RData")
library(ggplot2)
library(reshape2)
df_plot <- data.frame(
  BWMR = rowMeans(recordtime_mat_bwmr),
  RAPS = rowMeans(recordtime_mat_raps),
  Egger = rowMeans(recordtime_mat_egger),
  GSMR = rowMeans(recordtime_mat_gsmr),
  Nsnps = seq(10, 50, 10)
)
df_plot <- melt(df_plot, id=c('Nsnps'))
colnames(df_plot) <- c("Nsnps", "Method", "time")
method_level <- c("BWMR", "Egger", "GSMR", "RAPS")
df_plot$Method <- factor(df_plot$Method, levels = method_level)
plt <- ggplot(df_plot, aes(x = Nsnps, y = time, color = Method)) +
  geom_line() + 
  geom_point(size = 3, shape = 20) +
  labs(x = "Number of SNPs", y = "Time (seconds)", title = "Computational time") +
  theme(legend.position = "top") +
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plt


#------------------------------------------
rm(list = ls())
load("computation_time_N.RData")
library(ggplot2)
library(reshape2)
df_plot <- data.frame(
  BWMR = rowMeans(recordtime_mat_bwmr),
  RAPS = rowMeans(recordtime_mat_raps),
  Egger = rowMeans(recordtime_mat_egger),
  GSMR = rowMeans(recordtime_mat_gsmr),
  PRESSO = rowMeans(recordtime_mat_presso),
  Nsnps = seq(10, 50, 10)
)
df_plot <- melt(df_plot, id=c('Nsnps'))
colnames(df_plot) <- c("Nsnps", "Method", "time")
method_level <- c("BWMR", "Egger", "GSMR", "RAPS", "PRESSO")
df_plot$Method <- factor(df_plot$Method, levels = method_level)
df_plot$Nsnps <- factor(df_plot$Nsnps)
plt <- ggplot(df_plot, aes(x = Nsnps, y = time, fill = Method)) +
  geom_bar(position=position_dodge(), stat="identity") +
  labs(x = "Number of SNPs", y = "Time (seconds)", title = "Computational time") +
  theme(legend.position = "top") +
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plt


##### Case-5 in summary-level simulations (qq-plot, type I error rate, power) ####
# Figures S19 in Supplementary Document


rm(list = ls())

library(BWMR)
library(TwoSampleMR)
source("gsmr_md.R")
library(ggplot2)

library(distr)
D <- DExp(rate = 1) 

sim_without_outlier_laplace <- function(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true){
  gamma_true <- rnorm(N, 0, sqrt(sqsigma_true))
  gammahat <- rnorm(N, gamma_true, sigmaX_true)
  Gammahat <- rnorm(N, beta_true*gamma_true, sigmaY_true) + sqrt(sqtau_true)*r(D)(N)
  
  df_data <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX_true,
    sigmaY = sigmaY_true
  )
}

Rp <- 500
beta_0 <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
sqtau_true <- 0.3^2
sqsigma_true <- 0.8^2
sigmaX_min <- 0.3
sigmaX_max <- 0.5
sigmaY_min <- 0.3
sigmaY_max <- 0.5
N <- 50

L <- length(beta_0)

b.bwmr <- matrix(nrow = L, ncol = Rp)
b.raps <- matrix(nrow = L, ncol = Rp)
b.gsmr <- matrix(nrow = L, ncol = Rp)
b.egger <- matrix(nrow = L, ncol = Rp)
pval.bwmr <- matrix(nrow = L, ncol = Rp)
pval.raps <- matrix(nrow = L, ncol = Rp)
pval.gsmr <- matrix(nrow = L, ncol = Rp)
pval.egger <- matrix(nrow = L, ncol = Rp)

for (k in 1:L) {
  beta_true <- beta_0[k]
  for (iter in 1:Rp) {
    set.seed(iter)
    sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
    sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
    dat <- sim_without_outlier_laplace(beta_true, sqtau_true, sqsigma_true, N, sigmaX_true, sigmaY_true)
    colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
    # BWMR
    res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                         sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
    b.bwmr[k, iter] <- res_bwmr$beta
    pval.bwmr[k, iter] <- res_bwmr$P_value
    # GSMR
    res_gsmr = try(gsmr(bzx = dat$b.exposure, bzx_se = dat$se.exposure,  
                        bzy = dat$b.outcome, bzy_se = dat$se.outcome, 
                        ldrho = diag(N), 
                        heidi_outlier_flag = T, 
                        nsnps_thresh = 1))
    b.gsmr[k, iter] <- drop(res_gsmr$bxy)
    pval.gsmr[k, iter] <- drop(res_gsmr$bxy_pval)
    ## Two sample MR
    dat$id.exposure <- rep("uVXQCX", nrow(dat))
    dat$id.outcome <- rep("uVXQCX", nrow(dat))
    dat$mr_keep <- rep(TRUE, nrow(dat))
    dat$exposure <- rep('exposure', nrow(dat))
    dat$outcome <- rep('outcome', nrow(dat))
    colnames(dat) = c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "id.exposure", "id.outcome", "mr_keep", "exposure", "outcome")
    result <- try(mr(dat, method_list=c("mr_raps", "mr_egger_regression")))
    b.raps[k, iter] <- result$b[1]
    pval.raps[k, iter] <- result$pval[1]
    b.egger[k, iter] <- result$b[2]
    pval.egger[k, iter] <- result$pval[2]
  }
}

save.image("case-5-typeIerror-power.RData")


#################################################################################################################
### Plot
rm(list = ls())

library(ggplot2)
library(reshape2)
load("case-5-typeIerror-power.RData")

Rp <- 500
method_level <- c("BWMR", "Egger", "GSMR", "RAPS")
type.I.error <- c("BWMR" = sum(pval.bwmr[1, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[1, ] < 0.05)/Rp,
                  "GSMR" = sum(pval.gsmr[1, ] < 0.05)/Rp, "Egger" = sum(pval.egger[1, ] < 0.05)/Rp)
se <- c("BWMR" = sqrt(type.I.error["BWMR"]*(1 - type.I.error["BWMR"])/Rp),
        "RAPS" = sqrt(type.I.error["RAPS"]*(1 - type.I.error["RAPS"])/Rp),
        "GSMR" = sqrt(type.I.error["GSMR"]*(1 - type.I.error["GSMR"])/Rp),
        "Egger" = sqrt(type.I.error["Egger"]*(1 - type.I.error["Egger"])/Rp))
df_typeIerror <- data.frame(
  typeIerror = type.I.error,
  se = se,
  Method = names(type.I.error)
)
df_typeIerror$Method <- factor(df_typeIerror$Method, levels = method_level)
plt_typeIerror <- ggplot(df_typeIerror, aes(x = Method, y = typeIerror, fill = Method)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = typeIerror - se, ymax = typeIerror + se),
                width=.2,
                position=position_dodge(.9)) +
  geom_hline(aes(yintercept = 0.05), colour = "tomato1", linetype = "dashed", size = 1) +
  labs(x = "Method", y = "Type I error rate", title = "Type I error rate") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 25, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 30),
        axis.text = element_text(size = 25))
plt_typeIerror

## Power
power_1 <- c("BWMR" = sum(pval.bwmr[2, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[2, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr[2, ] < 0.05)/Rp, "Egger" = sum(pval.egger[2, ] < 0.05)/Rp)
power_2 <- c("BWMR" = sum(pval.bwmr[3, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[3, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr[3, ] < 0.05)/Rp, "Egger" = sum(pval.egger[3, ] < 0.05)/Rp)
power_3 <- c("BWMR" = sum(pval.bwmr[4, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[4, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr[4, ] < 0.05)/Rp, "Egger" = sum(pval.egger[4, ] < 0.05)/Rp)
power_4 <- c("BWMR" = sum(pval.bwmr[5, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[5, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr[5, ] < 0.05)/Rp, "Egger" = sum(pval.egger[5, ] < 0.05)/Rp)
power_5 <- c("BWMR" = sum(pval.bwmr[6, ] < 0.05)/Rp, "RAPS" = sum(pval.raps[6, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr[6, ] < 0.05)/Rp, "Egger" = sum(pval.egger[6, ] < 0.05)/Rp)
df_power <- data.frame(
  beta_true = rep(c(0.1, 0.2, 0.3, 0.4, 0.5), rep(4, 5)),
  power = c(power_1, power_2, power_3, power_4, power_5),
  Method = rep(names(power_5), 5)
)
df_power$Method <- factor(df_power$Method, levels = method_level)
plt_power <- ggplot(df_power, aes(x = beta_true, y = power, fill = Method)) +
  geom_bar(position=position_dodge(), stat="identity") +
  ylim(0, 1) +
  scale_x_continuous(breaks=seq(0.1, 0.5, 0.1)) +
  labs(x = "Beta", y = "Power", title = "Power") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 25, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 30),
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"))
plt_power

## qqplot
source('qqunifplot.R')
my.pvalue.list <- list("GSMR" = pval.gsmr[1, ][pval.gsmr[1, ] > 1e-8], "RAPS" = pval.raps[1, ][pval.raps[1, ] > 1e-8],
                       "Egger" = pval.egger[1, ][pval.egger[1, ] > 1e-8], "BWMR" = pval.bwmr[1, ][pval.bwmr[1, ] > 1e-8])
plt_qq <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
plt_qq



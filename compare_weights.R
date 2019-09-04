# Figure S6 in Supplementary Document

rm(list = ls())

source('func_bayesian_ivw_linearregressionxerror.R')
source('func_no_weighting.R')
set.seed(1234)

N <- 50
b_true <- 0.5
sqsigma_true <- 0.8^2
sqtau_true <- 0.3^2
sigma_X_set <- c(0.001, 0.01, 0.1, 0.2, 0.3)

Rp <- 500
L <- length(sigma_X_set)
b.bwmr <- matrix(NA, nrow = L, ncol = Rp)
b.bivw <- matrix(NA, nrow = L, ncol = Rp)
b.now <- matrix(NA, nrow = L, ncol = Rp)

for (i in 1:L) {
  sqsigma_X <- rep(sigma_X_set[i]^2, N)
  for (j in 1:Rp) {
    # simu
    x <- rnorm(N, 0, sqrt(sqsigma_true))
    xhat <- rnorm(N, x, sqrt(sqsigma_X))
    yhat <- rnorm(N, b_true*x, sqrt(sqtau_true))
    # BWMR
    res <- BWMR::BWMR(gammahat = xhat, Gammahat = yhat, sigmaX = sqrt(sqsigma_X), sigmaY = rep(0, N))
    b.bwmr[i,j] <- res$beta
    # Bayesian ivw
    b.bivw[i,j] <- bayesian_ivw(xhat = xhat, yhat = yhat, sqsigma_X = sqsigma_X)
    # no weighting
    b.now[i,j] <- no_weighting(xhat = xhat, yhat = yhat, sqsigma_X = sqsigma_X)
  }
}

save.image('compare_weights.RData')

#---------------------------------------------------------------------------------
rm(list = ls())

load('compare_weights.RData')
library(ggplot2)
library(reshape2)

df_plt <- data.frame(
  BWMR = c(b.bwmr[1,], b.bwmr[2,], b.bwmr[3,], b.bwmr[4,], b.bwmr[5,]),
  Bayesian_ivw = c(b.bivw[1,], b.bivw[2,], b.bivw[3,], b.bivw[4,], b.bivw[5,]),
  No_weighting = c(b.now[1,], b.now[2,], b.now[3,], b.now[4,], b.now[5,]),
  sigma_of_epsilonx = as.factor(rep(sigma_X_set, rep(Rp, L)))
)
df_plt <- melt(df_plt, id.vars = c("sigma_of_epsilonx"), variable.name = 'Method', value.name = 'beta.est')
plt <- ggplot(df_plt, aes(x = sigma_of_epsilonx, y = beta.est)) + 
  geom_boxplot(aes(fill = Method)) + 
  geom_hline(aes(yintercept = 0.5), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "sigma_x", title = "beta estimation") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size=25),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plt

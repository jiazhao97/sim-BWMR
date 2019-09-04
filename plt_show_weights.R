# Figure S7 in supplementary document

rm(list = ls())

source('func_bayesian_ivw_linearregressionxerror.R')

set.seed(1234)

N <- 50
b_true <- 0.5
sqsigma_true <- 0.8^2
sqtau_true <- 0.3^2
sqsigma_X <- 0.3^2
x <- rnorm(N, 0, sqrt(sqsigma_true))
xhat <- rnorm(N, x, sqrt(sqsigma_X))
yhat <- rnorm(N, b_true*x, sqrt(sqtau_true))
# BWMR
res <- BWMR::BWMR(gammahat = xhat, Gammahat = yhat, sigmaX = sqrt(sqsigma_X), sigmaY = rep(0, N))
w.bwmr <- res$weights
# Bayesian ivw
res <- bayesian_ivw(xhat = xhat, yhat = yhat, sqsigma_X = sqsigma_X)
w.bivw <- res$weights


df_plt <- data.frame(
  idx = seq(1,N),
  BWMR = w.bwmr,
  bivw = w.bivw
)
plt <- ggplot(df_plt, aes(x=idx, y=BWMR)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_hline(aes(yintercept = 1.0), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "observation No.", title = "weights (BWMR)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 20, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 28),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))
plt
plt <- ggplot(df_plt, aes(x=idx, y=bivw)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_hline(aes(yintercept = 1.0), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "observation No.", title = "weights (Bayesian ivw)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 20, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 28),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))
plt


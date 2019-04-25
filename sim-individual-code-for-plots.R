##### Code for plots of individual-level simulation results #####
# Vary "pl2", "pl5", "pl8" to get 
# Top panel of Figure 1 in Main text,
# Figure S20 in Supplementary Document,
# Bottom panel of Figure 1 in Main text,
# respectively.


### Plot
rm(list = ls())

library(ggplot2)
library(reshape2)

Rp <- 100
method_level <- c("BWMR", "Egger", "GSMR", "RAPS")

## Type I error
load("null-pl2.RData")
i <- 1
Rp <- 100
type.I.error <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
                  "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
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

### QQplot
load("null-pl2.RData")
source('qqunifplot.R')
my.pvalue.list <- list("GSMR" = pval.gsmr.unbias[1, ][pval.gsmr.unbias[1, ] > 1e-30], "RAPS" = pval.raps.unbias[1, ][pval.raps.unbias[1, ] > 1e-30],
                       "Egger" = pval.egger.unbias[1, ][pval.egger.unbias[1, ] > 1e-30], "BWMR" = pval.bwmr.unbias[1, ][pval.bwmr.unbias[1, ] > 1e-30])
plt_qq <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
plt_qq


## Power
load("cau1-pl2.RData")
i <- 1
Rp <- 100
power_1 <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
load("cau2-pl2.RData")
i <- 1
Rp <- 100
power_2 <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
load("cau3-pl2.RData")
i <- 1
Rp <- 100
power_3 <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
load("cau4-pl2.RData")
i <- 1
Rp <- 100
power_4 <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
load("cau5-pl2.RData")
i <- 1
Rp <- 100
power_5 <- c("BWMR" = sum(pval.bwmr.unbias[i, ] < 0.05)/Rp, "RAPS" = sum(pval.raps.unbias[i, ] < 0.05)/Rp,
             "GSMR" = sum(pval.gsmr.unbias[i, ] < 0.05)/Rp, "Egger" = sum(pval.egger.unbias[i, ] < 0.05)/Rp)
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


### Estimation (boxplot)
method_0 <- c("BWMR", "RAPS", "Egger", "GSMR")
beta_0 <- c("beta=0.0", "beta=0.1", "beta=0.2", "beta=0.3", "beta=0.4", "beta=0.5")
beta_value_0 <- seq(0.0, 0.5, 0.1)
# store estimate
Rp <- 100
b.mat <- matrix(nrow = length(beta_0), ncol = length(method_0)*Rp)
b.mat.b <- matrix(nrow = length(beta_0), ncol = length(method_0)*Rp)
rownames(b.mat) <- beta_0
load("null-pl2.RData")
i <- 1
Rp <- 100
b.mat[1, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[1, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])
load("cau1-pl2.RData")
i <- 1
Rp <- 100
b.mat[2, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[2, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])
load("cau2-pl2.RData")
i <- 1
Rp <- 100
b.mat[3, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[3, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])
load("cau3-pl2.RData")
i <- 1
Rp <- 100
b.mat[4, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[4, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])
load("cau4-pl2.RData")
i <- 1
Rp <- 100
b.mat[5, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[5, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])
load("cau5-pl2.RData")
i <- 1
Rp <- 100
b.mat[6, ] <- c(b.bwmr.unbias[i, ], b.raps.unbias[i, ], b.egger.unbias[i, ], b.gsmr.unbias[i, ])
b.mat.b[6, ] <- c(b.bwmr.bias[i, ], b.raps.bias[i, ], b.egger.bias[i, ], b.gsmr.bias[i, ])

beta <- rep(beta_0, rep(length(method_0)*Rp, length(beta_0)))
beta_value <- rep(beta_value_0, rep(length(method_0)*Rp, length(beta_0)))
method <- rep(rep(method_0, rep(Rp, length(method_0))), length(beta_0))

est_df <- data.frame(
  beta = beta,
  beta_value = beta_value,
  Method = method,
  beta_est = c(b.mat[1, ], b.mat[2, ], b.mat[3, ], b.mat[4, ], b.mat[5, ], b.mat[6, ]),
  beta_est_bias = c(b.mat.b[1, ], b.mat.b[2, ], b.mat.b[3, ], b.mat.b[4, ], b.mat.b[5, ], b.mat.b[6, ])
)

est_plt <- ggplot(est_df, aes(x = Method, y = beta_est)) +
  geom_boxplot(aes(fill = Method)) +
  geom_hline(aes(yintercept = beta_value), colour="blue", linetype="dashed", size=1) +
  labs(x = "Method", title="Estimation") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 25, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 30),
        axis.text = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
est_plt <- est_plt + facet_grid( ~ beta) +
  theme(strip.text.x = element_text(size = 20))
est_plt



# Figure S25 in supplementary document

rm(list = ls())

library(ggplot2)
library(reshape2)

method_level <- c("BWMR", "Egger", "GSMR", "RAPS")

## Type I error
load("heritability-02.RData")
#load("heritability-03.RData")
#load("heritability-04.RData")

i <- 1
Rp <- length(pvalset.bwmr[i, ])
type.I.error <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
                  "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
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
  theme(axis.title.x = element_text(size = 25, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text = element_text(size = 25))
plt_typeIerror

### QQplot
source('qqunifplot.R')
my.pvalue.list <- list("GSMR" = pvalset.gsmr[1, ][pvalset.gsmr[1, ] > 1e-30], "RAPS" = pvalset.raps[1, ][pvalset.raps[1, ] > 1e-30],
                       "Egger" = pvalset.egger[1, ][pvalset.egger[1, ] > 1e-30], "BWMR" = pvalset.bwmr[1, ][pvalset.bwmr[1, ] > 1e-30])
plt_qq <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
plt_qq


## Power
i <- 2
power_1 <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
             "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
i <- 3
power_2 <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
             "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
i <- 4
power_3 <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
             "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
i <- 5
power_4 <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
             "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
i <- 6
power_5 <- c("BWMR" = sum(pvalset.bwmr[i, ] < 0.05)/Rp, "RAPS" = sum(pvalset.raps[i, ] < 0.05)/Rp,
             "GSMR" = sum(pvalset.gsmr[i, ] < 0.05)/Rp, "Egger" = sum(pvalset.egger[i, ] < 0.05)/Rp)
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
  theme(axis.title.x = element_text(size = 25, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"))
plt_power


### Estimation (boxplot)
method_0 <- c("BWMR", "RAPS", "Egger", "GSMR")
beta_0 <- c("beta=0.0", "beta=0.1", "beta=0.2", "beta=0.3", "beta=0.4", "beta=0.5")
beta_value_0 <- seq(0.0, 0.5, 0.1)
# store estimate
b.mat <- matrix(nrow = length(beta_0), ncol = length(method_0)*Rp)
b.mat.b <- matrix(nrow = length(beta_0), ncol = length(method_0)*Rp)
rownames(b.mat) <- beta_0
for (i in 1:6) {
  b.mat[i, ] <- c(bset.bwmr[i, ], bset.raps[i, ], bset.egger[i, ], bset.gsmr[i, ])
}
beta <- rep(beta_0, rep(length(method_0)*Rp, length(beta_0)))
beta_value <- rep(beta_value_0, rep(length(method_0)*Rp, length(beta_0)))
method <- rep(rep(method_0, rep(Rp, length(method_0))), length(beta_0))

est_df <- data.frame(
  beta = beta,
  beta_value = beta_value,
  Method = method,
  beta_est = c(b.mat[1, ], b.mat[2, ], b.mat[3, ], b.mat[4, ], b.mat[5, ], b.mat[6, ])
)

est_plt <- ggplot(est_df, aes(x = Method, y = beta_est)) +
  geom_boxplot(aes(fill = Method)) +
  geom_hline(aes(yintercept = beta_value), colour="blue", linetype="dashed", size=1) +
  labs(x = "Method", title="Estimation") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 25, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
est_plt <- est_plt + facet_grid( ~ beta) +
  theme(strip.text.x = element_text(size = 20))
est_plt

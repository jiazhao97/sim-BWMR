# Figure S24 in supplementary document

rm(list = ls())

library(ggplot2)
library(reshape2)

Rp <- 100
method_level <- c("BWMR", "Egger", "GSMR", "RAPS")

## Type I error
load("pl8-seed1-gpu-U-qu01-bu1.RData")
#load("pl8-seed1-gpu-U-qu01-bu05.RData")
#load("pl8-seed1-gpu-U-qu01-bu2.RData")
#load("pl8-seed1-gpu-U-qu02-bu1.RData")
#load("pl8-seed1-gpu-U-qu04-bu1.RData")

i <- 1
Rp <- length(pval.bwmr.bias[1, ])
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
  labs(x = "Method", y = "Type I error rate", title = "Type I error rate (qu = 0.4, bu = 1)") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 25, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 30),
        axis.text = element_text(size = 25))
plt_typeIerror

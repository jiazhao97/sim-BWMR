### Plot
rm(list = ls())

library(ggplot2)
library(reshape2)

Rp <- 100
method_level <- c("BWMR.u", "BWMR.b", "Egger.u", "Egger.b", 
                  "GSMR.u", "GSMR.b", "RAPS.u", "RAPS.b")

### Estimation (boxplot)
method_0 <- c("BWMR.u", "BWMR.b", "RAPS.u", "RAPS.b", "Egger.u", "Egger.b", "GSMR.u", "GSMR.b")
# samplesize_0 <- c("samplesize=4000", "samplesize=5000", "samplesize=6000", "samplesize=8000")
samplesize_0 <- c("sample size=4000", "sample size=5000", "sample size=6000", "sample size=8000")
beta_value_0 <- seq(0.0, 0.5, 0.1)
# store estimate
b.mat <- matrix(nrow = length(samplesize_0), ncol = length(method_0)*Rp)
rownames(b.mat) <- samplesize_0
# load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/new-small-sderror/0405/cau5-noise11-pl2-ss4000.RData")
# i <- 1
# b.mat[1, ] <- c(b.bwmr.unbias[i, ], b.bwmr.bias[i, ], b.raps.unbias[i, ], b.raps.bias[i, ], 
#                 b.egger.unbias[i, ], b.egger.bias[i, ], b.gsmr.unbias[i, ], b.gsmr.bias[i, ])
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss4000.RData")
i <- 1
b.mat[1, ] <- c(b.bwmr.unbias[i, ], b.bwmr.bias[i, ], b.raps.unbias[i, ], b.raps.bias[i, ], 
                b.egger.unbias[i, ], b.egger.bias[i, ], b.gsmr.unbias[i, ], b.gsmr.bias[i, ])
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss5000.RData")
i <- 1
b.mat[2, ] <- c(b.bwmr.unbias[i, ], b.bwmr.bias[i, ], b.raps.unbias[i, ], b.raps.bias[i, ], 
                b.egger.unbias[i, ], b.egger.bias[i, ], b.gsmr.unbias[i, ], b.gsmr.bias[i, ])
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss6000.RData")
i <- 1
b.mat[3, ] <- c(b.bwmr.unbias[i, ], b.bwmr.bias[i, ], b.raps.unbias[i, ], b.raps.bias[i, ], 
                b.egger.unbias[i, ], b.egger.bias[i, ], b.gsmr.unbias[i, ], b.gsmr.bias[i, ])
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss8000.RData")
i <- 1
b.mat[4, ] <- c(b.bwmr.unbias[i, ], b.bwmr.bias[i, ], b.raps.unbias[i, ], b.raps.bias[i, ], 
                b.egger.unbias[i, ], b.egger.bias[i, ], b.gsmr.unbias[i, ], b.gsmr.bias[i, ])

samplesize <- rep(samplesize_0, rep(length(method_0)*Rp, length(samplesize_0)))
# beta_value <- rep(beta_value_0, rep(length(method_0)*Rp, length(beta_0)))
method <- rep(rep(method_0, rep(Rp, length(method_0))), length(samplesize_0))

est_df <- data.frame(
  samplesize = samplesize,
  # beta_value = beta_value,
  Method = method,
  # beta_est = c(b.mat[1, ], b.mat[2, ], b.mat[3, ], b.mat[4, ])
  beta_est = c(b.mat[1, ], b.mat[2, ], b.mat[3, ], b.mat[4, ])
)
est_df$Method <- factor(est_df$Method, levels = method_level)
est_plt <- ggplot(est_df, aes(x = Method, y = beta_est)) +
  geom_boxplot(aes(fill = Method)) +
  geom_hline(aes(yintercept = 0.5), colour="blue", linetype="dashed", size=1) +
  labs(x = "Method", title="Estimation") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 25, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"))   
est_plt <- est_plt + facet_grid( ~ samplesize) +
  theme(strip.text.x = element_text(size = 20))
est_plt

#---------------------------------------------------------------------------------------
## same method; different sample size
samplesize_0 <- c("4000", "5000", "6000", "8000")
# store estimate
b.mat <- matrix(nrow = length(samplesize_0), ncol = Rp)
rownames(b.mat) <- samplesize_0
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss4000.RData")
i <- 1
b.mat[1, ] <- b.egger.bias[i, ]
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss5000.RData")
i <- 1
b.mat[2, ] <- b.egger.bias[i, ]
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss6000.RData")
i <- 1
b.mat[3, ] <- b.egger.bias[i, ]
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/sim-BWMR/individual-level/genotype012-new/selectionbias/cau5-noise11-pl2-ss8000.RData")
i <- 1
b.mat[4, ] <- b.egger.bias[i, ]

est_df <- data.frame(
  samplesize = rep(samplesize_0, rep(Rp, length(samplesize_0))),
  beta_est = c(b.mat[1, ], b.mat[2, ], b.mat[3, ], b.mat[4, ])
)

est_plt <- ggplot(est_df, aes(x = samplesize, y = beta_est)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0.5), colour="blue", linetype="dashed", size=1) +
  labs(x = "Sample size", title="Estimation with selection bias (Egger)") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 25, vjust = -0.5),
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 30),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "top") +
  # theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold"))   
est_plt


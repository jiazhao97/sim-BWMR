rm(list = ls())

library(BWMR)
library(ggplot2)

sim_beta_corrupt <- function(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true){
  gamma_true <- rnorm(N, 0, sqrt(sqsigma_true))
  gammahat <- rnorm(N, gamma_true, sigmaX_true)
  
  # the first N*(1-C) data of Gammahat are reliable
  # the last N*C data of Gammahat are corrupted
  h <- as.integer(N*(1-C))
  Gammahat <- rnorm(h, beta_true*gamma_true[1:h], sqrt((sigmaY_true^2 + sqtau_true)[1:h]))
  Gammahat <- c(Gammahat, rnorm((N-h), beta_corrupted*gamma_true[(h+1):N], sqrt((sigmaY_true^2 + sqtau_true)[(h+1):N])))
  
  df_data <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX_true,
    sigmaY = sigmaY_true
  )
}

Rp <- 200
N <- 20

sigmaX_min <- 0.3
sigmaX_max <- 0.5
sigmaY_min <- 0.3
sigmaY_max <- 0.5

beta_corrupted <- 3
sqsigma_true <- 1.0^2
sqtau_true <- 0.2^2
beta_true <- 0.0
C <- 0.1

# declare sets
w_mat <- matrix(nrow = Rp, ncol = N)
b_set <- numeric(length = Rp)

for (iter in 1:Rp) {
  set.seed(iter)
  sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
  sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
  dat <- sim_beta_corrupt(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true)
  colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
  # BWMR
  res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                       sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
  w_mat[iter, ] <- res_bwmr$weights
  b_set[iter] <- res_bwmr$beta
}

a <- ls()
rm(list=a[which(a!='w_mat' & a !='b_set')])
save.image("weights-N20-C01-corbeta3.RData")


#-------------------------
rm(list = ls())

library(BWMR)
library(ggplot2)

sim_beta_corrupt <- function(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true){
  gamma_true <- rnorm(N, 0, sqrt(sqsigma_true))
  gammahat <- rnorm(N, gamma_true, sigmaX_true)
  
  # the first N*(1-C) data of Gammahat are reliable
  # the last N*C data of Gammahat are corrupted
  h <- as.integer(N*(1-C))
  Gammahat <- rnorm(h, beta_true*gamma_true[1:h], sqrt((sigmaY_true^2 + sqtau_true)[1:h]))
  Gammahat <- c(Gammahat, rnorm((N-h), beta_corrupted*gamma_true[(h+1):N], sqrt((sigmaY_true^2 + sqtau_true)[(h+1):N])))
  
  df_data <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX_true,
    sigmaY = sigmaY_true
  )
}

Rp <- 200
N <- 20

sigmaX_min <- 0.3
sigmaX_max <- 0.5
sigmaY_min <- 0.3
sigmaY_max <- 0.5

beta_corrupted <- 4
sqsigma_true <- 1.0^2
sqtau_true <- 0.2^2
beta_true <- 0.0
C <- 0.1

# declare sets
w_mat <- matrix(nrow = Rp, ncol = N)
b_set <- numeric(length = Rp)

for (iter in 1:Rp) {
  set.seed(iter)
  sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
  sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
  dat <- sim_beta_corrupt(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true)
  colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
  # BWMR
  res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                       sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
  w_mat[iter, ] <- res_bwmr$weights
  b_set[iter] <- res_bwmr$beta
}

a <- ls()
rm(list=a[which(a!='w_mat' & a !='b_set')])
save.image("weights-N20-C01-corbeta4.RData")


#-------------------------
rm(list = ls())

library(BWMR)
library(ggplot2)

sim_beta_corrupt <- function(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true){
  gamma_true <- rnorm(N, 0, sqrt(sqsigma_true))
  gammahat <- rnorm(N, gamma_true, sigmaX_true)

  # the first N*(1-C) data of Gammahat are reliable
  # the last N*C data of Gammahat are corrupted
  h <- as.integer(N*(1-C))
  Gammahat <- rnorm(h, beta_true*gamma_true[1:h], sqrt((sigmaY_true^2 + sqtau_true)[1:h]))
  Gammahat <- c(Gammahat, rnorm((N-h), beta_corrupted*gamma_true[(h+1):N], sqrt((sigmaY_true^2 + sqtau_true)[(h+1):N])))

  df_data <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX_true,
    sigmaY = sigmaY_true
  )
}

Rp <- 200
N <- 20

sigmaX_min <- 0.3
sigmaX_max <- 0.5
sigmaY_min <- 0.3
sigmaY_max <- 0.5

beta_corrupted <- 2
sqsigma_true <- 1.0^2
sqtau_true <- 0.2^2
beta_true <- 0.0
C <- 0.1

# declare sets
w_mat <- matrix(nrow = Rp, ncol = N)
b_set <- numeric(length = Rp)

for (iter in 1:Rp) {
  set.seed(iter)
  sigmaX_true <- runif(N, sigmaX_min, sigmaX_max)
  sigmaY_true <- runif(N, sigmaY_min, sigmaY_max)
  dat <- sim_beta_corrupt(beta_true, sqtau_true, sqsigma_true, N, C, beta_corrupted, sigmaX_true, sigmaY_true)
  colnames(dat) <- c("b.exposure", "b.outcome", "se.exposure", "se.outcome")
  # BWMR
  res_bwmr <- try(BWMR(gammahat = dat$b.exposure, Gammahat = dat$b.outcome,
                       sigmaX = dat$se.exposure, sigmaY = dat$se.outcome))
  w_mat[iter, ] <- res_bwmr$weights
  b_set[iter] <- res_bwmr$beta
}

a <- ls()
rm(list=a[which(a!='w_mat' & a !='b_set')])
save.image("weights-N20-C01-corbeta2.RData")




#------
rm(list = ls())
library(ggplot2)
library(reshape2)
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/BWMR-bioinformatics/major_revision/code-RData/weights-reviewer1_Major_comment/weights-N20-C01-corbeta3.RData")

median_set <- numeric(ncol(w_mat))
for (i in 1:ncol(w_mat)) {median_set[i] <- summary(w_mat[, i])['Median']}
df_plot <- data.frame(
  indx = factor(seq(1, ncol(w_mat), 1)),
  Mean = colMeans(w_mat),
  Median = median_set
)
df_plot <- melt(df_plot, id=c('indx'))
colnames(df_plot)[2] <- 'Statistics'
plt <- ggplot(df_plot, aes(x=indx, y=value, fill=Statistics)) +
  scale_fill_brewer(palette="Set3") +
  geom_bar(position=position_dodge(), stat="identity") +
  labs(x = "Observation No.", title = "beta_corrupted=2") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 20, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 28),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))
plt


#------
rm(list = ls())
library(ggplot2)
library(reshape2)
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/BWMR-bioinformatics/major_revision/code-RData/weights-reviewer1_Major_comment/weights-N20-C01-corbeta2.RData")
b_mat <- matrix(nrow = nrow(w_mat), ncol = 3)
b_mat[,1] <- b_set
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/BWMR-bioinformatics/major_revision/code-RData/weights-reviewer1_Major_comment/weights-N20-C01-corbeta3.RData")
b_mat[,2] <- b_set
load("/Users/jiazhao/Documents/HKUST/201809 BWMR/BWMR-bioinformatics/major_revision/code-RData/weights-reviewer1_Major_comment/weights-N20-C01-corbeta4.RData")
b_mat[,3] <- b_set

colnames(b_mat) <- seq(2, 4, 1)
df_plot <- as.data.frame(b_mat)
df_plot <- melt(df_plot, id=c())

plt <- ggplot(df_plot, aes(x=variable, y=value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0.0), colour = "blue", linetype = "dashed", size = 1) +
  labs(x = "beta_corrupted", title = "beta estimates") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 20, vjust = -1.5),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 28),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))
plt


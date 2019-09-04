###### function for generating data in summary-level simulation ######
##
### INPUT
## p:                                number of snps;
## n1:                               sample size in database for the exposure trait;
## n2:                               sample size in database for the outcome trait;
## h1:                               heritability of the exposure trait;
## h2:                               heritability of the outcome trait;
## neffect1:                         average number of the exposure trait effect snps;
## neffect2:                         average number of the outcome trait direct effect snps (horizonal pleiotropy);
## beta:                             effect of the exposure trait on the outcome trait;
## pval_thresh:                      p-value threshold for IVs selection.
##
## * Note that INPUT should statisfy the condition:  beta^2*h1 <= h2
##
### OUTPUT
## dataframe containing beta.exp, beta.out, se.exp, se.out.


simu_sum <- function(p, n1, n2, h1, h2, neffect1, neffect2, beta, pval_thresh) {
  # check the condition for INPUT
  if (beta^2*h1 > h2) {stop("Please check the condition for the input data.")}
  
  # exposure trait
  p1 <- neffect1/p
  sigma_1 <- sqrt(h1/neffect1)                               # note: sigma_1 corresponds to the parameter sigma 
                                                             #       (sqrt(var(gamma))) in the BWMR model.
  b1 <- rnormalmix(n = p, pi = c(p1, 1-p1), mean = c(0, 0), sd = c(sigma_1, 0))   
                                                             # note: b1 corresponds to gamma in the BWMR model.
  b1.se <- sqrt(1/n1)
  b1.hat1 <- rnorm(n = p, mean = b1, sd = b1.se)             # note: we generate two datasets for the exposure traits
  b1.hat2 <- rnorm(n = p, mean = b1, sd = b1.se)             #       to correct the selection bias.
  
  # outcome trait
  p2 <- neffect2/p
  sigma_2 <- sqrt((h2-beta^2*h1)/neffect2)                   # note: sigma_2 corresponds to the parameter tau
                                                             #       (sqrt(var(Gamma-beta*gamma))) in the BWMR model.
  b2 <- beta*b1 + rnormalmix(n = p, pi = c(p2, 1-p2), mean = c(0, 0), sd = c(sigma_2, 0))  
                                                             # note: b2 corresponds to Gamma in the BWMR model.
  b2.se <- sqrt(1/n2)
  b2.hat <- rnorm(n = p, mean = b2, sd = b2.se)
  
  # IVs selection
  b1.pval1 <- pchisq((b1.hat1/b1.se)^2, 1, lower.tail = FALSE)
  indx <- which(b1.pval1 < pval_thresh)
  nsnps <- length(indx)
  return(data.frame(beta.exposure = b1.hat2[indx],
                    beta.outcome  = b2.hat[indx],
                    se.exposure   = rep(b1.se, nsnps),
                    se.outcome    = rep(b2.se, nsnps)))
}

no_weighting <- function(xhat, yhat, sqsigma_X) {
  N <- length(xhat)
  # vem Bayesian_ivw
  # fixed 
  a0 <- 1
  b0 <- 1
  sqsigma0 <- (1e+6)^2
  a_w <- a0+0.5
  total_iter <- 1e+3
  # init
  sqtau <- 1
  sqsigma <- 1
  mu_x <- xhat
  sqsigma_x <- rep(1, N)
  b_w <- rep(a0+0.5, N)
  mu_w <- a_w/b_w
  elbo_set <- numeric(total_iter)
  # vem
  for (iter in 1:total_iter) {
    # b
    sqsigma_b <- 1 / (1/sqsigma0 + sum(mu_w*(mu_x^2+sqsigma_x))/sqtau)
    mu_b <- sqsigma_b * sum(mu_w*yhat*mu_x) / sqtau
    # x
    sqsigma_x <- 1 / (1/sqsigma + 1/sqsigma_X + mu_w*(mu_b^2+sqsigma_b)/sqtau)
    mu_x <- sqsigma_x * (xhat/sqsigma_X + mu_w*yhat*mu_b/sqtau)
    # w
    b_w <- b0 + 0.5/sqtau*((mu_b^2+sqsigma_b)*(mu_x^2+sqsigma_x) - 2*mu_b*mu_x*yhat + yhat^2)
    mu_w <- a_w/b_w
    mu_logw <- digamma(a_w) - log(b_w)
    mu_w <- rep(1, N)
    # sqsigma
    sqsigma <- sum(mu_x^2+sqsigma_x)/N
    # sqtau
    sqtau <- sum(mu_w*((mu_b^2+sqsigma_b)*(mu_x^2+sqsigma_x) - 2*mu_b*mu_x*yhat + yhat^2)) / N
    # elbo
    elbo <- -0.5*N*log(sqsigma) - 0.5/sqsigma*sum(mu_x^2+sqsigma_x) -
      0.5*sum(((xhat-mu_x)^2+sqsigma_x) / sqsigma_X) -
      0.5*N*log(sqtau) + sum(0.5*mu_logw - 0.5/sqtau*mu_w*((mu_b^2+sqsigma_b)*(mu_x^2+sqsigma_x) - 2*mu_b*mu_x*yhat + yhat^2)) -
      0.5/sqsigma0*(mu_b^2+sqsigma_b) +
      sum((a0-1)*mu_logw - b0*mu_w) +
      0.5*log(sqsigma_b) +
      sum(0.5*log(sqsigma_x)) -
      sum(a_w*log(b_w) + (a_w-1)*mu_logw - b_w*mu_w - lgamma(a_w))
    elbo_set[iter] <- elbo
    if (iter > 1) {
      if (elbo_set[iter] < elbo_set[iter-1]) {
        message('elbo decreases!')
        break
      }
      if (abs((elbo_set[iter]-elbo_set[iter-1])/elbo_set[iter-1]) < 1e-6) {
        message('elbo converges...')
        elbo_set <- elbo_set[1:iter]
        message('iter=', iter, ', b=', mu_b, ', sigma=', sqrt(sqsigma), ', tau=', sqrt(sqtau), '.')
        break
      }
    }
  }
  return(mu_b)
}
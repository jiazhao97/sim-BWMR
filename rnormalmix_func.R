###### function for generating data from mixed normal distribution ######
##
### INPUT
## n:                               number of observations to simulate;
## pi:                              proportion of each component in mixed Gaussian;
## mean:                            mean of each component in mixed Gaussian;
## sd:                              standard error of each component in mixed Gaussian;
##
### OUTPUT
## output                           samples from mixed normal distribution (n-vector).


rnormalmix <- function(n, pi, mean, sd) {
  # check the validity of inputs
  if (sum(pi) != 1) {stop("Sum of the proportions pi is not equivalent to one.")}
  n_comp <- length(pi)              # number of components
  if (length(mean) != n_comp | length(sd) != n_comp) {stop("Please check the validity of inputs.")}
  z <- sample(1:n_comp, size = n, replace = TRUE, prob = pi)
  output <- rep(NA, n)
  for (k in 1:n_comp) {
    nk <- sum(z==k)
    if (nk == 0) {next}
    output[z==k] <- rnorm(n = nk, mean = mean[k], sd = sd[k])
  }
  return(output)
}

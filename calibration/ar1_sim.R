##==============================================================================
## ar1_sim.R
##
## Original code for BRICK written by Tony Wong and Alexander Bakker
##   (https://github.com/scrim-network/BRICK, DOI:10.5194/gmd-10-2741-2017)
## simulate stationary AR(1) process (approximate - faster, better convergence, and
## results not sensitive to use of this as opposed to exact AR1)
##==============================================================================
ar1.sim = function(N,rho1,sigma) {
  x = rep(NA,N)
  if(length(sigma) > 1) {
    x[1] = rnorm(n=1,sd=sigma[1]/sqrt(1-rho1^2))
    for (i in 2:N) {
      x[i] = rho1*x[i-1] + rnorm(1,sd=sigma[i])
    }
  } else {
    x[1] = rnorm(n=1,sd=sigma/sqrt(1-rho1^2))
    for (i in 2:N) {
      x[i] = rho1*x[i-1] + rnorm(1,sd=sigma)
    }
  }
  return(x)
}

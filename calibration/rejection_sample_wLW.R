##==============================================================================
## rejection_sample.R
## Perform rejection sampling on a set of model slr time series,
## accounting for landwater contribution to total sea-level rise
##
## Original version for BRICK by Tony Wong
## Modified slightly for Hector by Ben Vega-Westhoff
##
## Output: model member indices that pass the rejection sampling
##==============================================================================

rejection_sample_wLW = function( 
                    obs.sl,         # observed total SL
                    obs.sl.time,    # corresponding times (years)
                    obs.sl.err,     # corresponding errors
                    slr.norm.stat ) # 2-d array containing the normalized model total sl, w/
                                    # added statistical noise. 1st index: ensemble member #,
                                    # second index: time.
{
  n.ensemble = length( slr.norm.stat[,1] )
  # 1901-1990: -0.11 [-0.16 to -0.06] (5-95% range)
  lw.time.1900 <- 1900:1989
  i1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
  lw.1900 <- (-0.11/1000)*(lw.time.1900 - 1900)
  lw.err.1900 <- (0.25*(-0.06--0.16)/1000)*sqrt(lw.time.1900 - lw.time.1900[1])
  # 1971-2010: 0.12 [0.03 to 0.22]
  lw.time.1970 <- 1970:2005 #only ran to 2005, changing for now (BRICK calib ran thru 2009)
  i1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
  lw.1970 <- (0.12/1000)*(lw.time.1970 - lw.time.1970[1])
  lw.err.1970 <- (0.25*(0.2-0.03)/1000)*sqrt(lw.time.1970 - lw.time.1970[1])
  # 1993-2010: 0.38 [0.26 to 0.49]
  lw.time.1992 <- 1992:2005 #only ran thru 2005, changing for now (BRICK calib ran thru 2009)
  i1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
  lw.1992 <- (0.38/1000)*(lw.time.1992 - lw.time.1992[1])
  lw.err.1992 <- (0.25*(0.49-0.26)/1000)*sqrt(lw.time.1992 - lw.time.1992[1])

  # normalize, subtract and add error in quadrature
  obs.sl.lw.1900 <- obs.sl[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])] - obs.sl[which(obs.sl.time==lw.time.1900[1])]
  obs.sl.lw.1970 <- obs.sl[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])] - obs.sl[which(obs.sl.time==lw.time.1970[1])]
  obs.sl.lw.1992 <- obs.sl[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])] - obs.sl[which(obs.sl.time==lw.time.1992[1])]

  obs.sl.lw.1900 <- obs.sl.lw.1900 - lw.1900
  obs.sl.lw.1970 <- obs.sl.lw.1970 - lw.1970
  obs.sl.lw.1992 <- obs.sl.lw.1992 - lw.1992

  obs.sl.lw.err.1900 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]^2 + lw.err.1900^2)
  obs.sl.lw.err.1970 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]^2 + lw.err.1970^2)
  obs.sl.lw.err.1992 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]^2 + lw.err.1992^2)

  # calculate likelihood as the product of the three independent likelihoods
  resid.1900 <- obs.sl.lw.1900 - obs.sl.lw.1900
  llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
  resid.1970 <- obs.sl.lw.1970 - obs.sl.lw.1970
  llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
  resid.1992 <- obs.sl.lw.1992 - obs.sl.lw.1992
  llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
  lik.max <- (llik.1900 + llik.1970 + llik.1992)/n.ensemble

  imod.1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
  imod.1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
  imod.1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
  
  survive = rep(0, n.ensemble)
  uni.rnd = log(runif(n.ensemble))
  for (i in 1:n.ensemble) {
    resid.1900 <- obs.sl.lw.1900 - (slr.norm.stat[i,imod.1900]-slr.norm.stat[i,imod.1900[1]])
    resid.1970 <- obs.sl.lw.1970 - (slr.norm.stat[i,imod.1970]-slr.norm.stat[i,imod.1970[1]])
    resid.1992 <- obs.sl.lw.1992 - (slr.norm.stat[i,imod.1992]-slr.norm.stat[i,imod.1992[1]])
    llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
    llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
    llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
    lik.mem <- llik.1900 + llik.1970 + llik.1992
    if( uni.rnd[i] <= lik.mem-lik.max) {survive[i]=1}
  }
  ind.survive = which( as.logical(survive))
  print(paste('Calibration to sea level data by rejection sampling leaves ',length(ind.survive),' full calibrated ensemble members',sep=''))
    
  return(ind.survive) 
}

##==============================================================================
## assim_chains.R
## 
## Original code for BRICK written by Tony Wong and Alexander Bakker
##   (https://github.com/scrim-network/BRICK, DOI:10.5194/gmd-10-2741-2017)
## Modified for Hector by Ben Vega-Westhoff
##==============================================================================
## Assimilate all the parallel chains from their respective .RData files
## Use Gelman-Rubin diagnostic to determine burn-in
## Save posterior parameter sets (after burn-in) to csv file
##==============================================================================

rm(list=ls())           # Clear all previous variables

## Required packages/libraries
library(optparse)       # Option parsing
library(adaptMCMC)      # Use robust adaptive Metropolis MCMC method

## Inputs
option_list = list(
    make_option( c("-f", "--folder"), type="character", default=NULL       ,
                 help = "folder where chains were saved", metavar="character") )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$folder)){
    print_help(opt_parser)
    stop("Calibration folder must be supplied.",
    call.=FALSE)
}

calib.folder = opt$folder #Where were the chains saved

if( strsplit( calib.folder,"/" )[[1]][1] != "" ) {
    calib.folder = paste0( getwd(), "/", calib.folder ) # Convert to absolute path
}

load( paste0( calib.folder, 'hector_calib_MCMC_ready.RData' ) )

##==============================================================================

## Determine when (in increments of ??,000 iterations, using Gelman and Rubin
## diagnostic) the two parallel chains converge.
## Initialize the testing of the Gelman and Rubin diagnostics
gr.step = floor( niter.mcmc/100 )
niter.test = seq( from = max(gr.step, 100), to = niter.mcmc, by = gr.step )
gr.test = rep( NA, length(niter.test) )
gr.max = 1.1

## Calculate the statistic at a few spots throughout the chain. Once it is
## close to 1 (people often use GR<1.1 or 1.05), the between-chain variability
## is indistinguishable from the within-chain variability, and they are
## converged. It is only after the chains are converged that you should use the
## parameter values as posterior draws, for analysis.
string.mcmc.list <- 'mcmc1'
for ( m in 2:nparallel.mcmc ) {
    string.mcmc.list <- paste0( string.mcmc.list, ', mcmc', m )
}

amcmc.par = list()
for( m in 1:nparallel.mcmc ) {
  dat = paste0( calib.folder, "chain_", m, '.RData' )
  load(dat)
  amcmc.par[[m]] = amcmc.out
}

for ( i in 1:length(niter.test) ) {
    for( m in 1:nparallel.mcmc ) {
        eval(parse(text=paste0('mcmc',m,' <- as.mcmc(amcmc.par[[m]]$samples[1:niter.test[i],])')))
    }
    eval(parse(text=paste0('mcmc_chain_list = mcmc.list(list(',string.mcmc.list,'))')))
    gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
}

## Plot GR statistics as a function of iteration
pdf(paste0( calib.folder, 'gr_stat_hector.pdf' ) )
plot( niter.test, gr.test )
abline(h=gr.max,lty=2)
dev.off()
##==============================================================================
# Chop off burn-in
#===============================================================================
#

# Note: here, we are only using the Gelman and Rubin diagnostic. But this is
# only after looking at the quantile stability as iterations increase, as well
# as the Heidelberger and Welch diagnostics, which suggest the chains are okay.
# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models.
# save a separate ifirst for each experiment
ifirst <- NULL
lgr <- rep(NA, length(niter.test))
for (i in 1:length(niter.test)) {lgr[i] <- gr.test[i] < gr.max}
for (i in seq(from=length(niter.test), to=1, by=-1)) {
  if( all(lgr[i:length(lgr)]) ) {ifirst <- niter.test[i]}
}
if ( is.null(ifirst) ) {
  print( "Unable to find convergence w/ GR diagnostic." )
  chains_file = paste0(calib.folder, "/unconverged_chain_list.RData")
  unconverged_chain_list = mcmc_chain_list
  save(unconverged_chain_list, file = chains_file)
  print(paste0("MCMC list saved to ", chains_file))
}

chains_burned <- vector( 'list', nparallel.mcmc )
for( m in 1:nparallel.mcmc) {
  chains_burned[[m]] <-  amcmc.par[[m]]$samples[(ifirst+1):niter.mcmc,]
}

# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior.
parameters.posterior <- chains_burned[[1]]
if( length( chains_burned ) > 1 ) {
  for (m in 2:length ( chains_burned ) ) {
    parameters.posterior <- rbind(parameters.posterior, chains_burned[[m]])
  }
}

# Only saving the transition covariance matrix for one of the chains
covjump.posterior <- amcmc.par[[1]]$cov.jump

n.parameters = ncol(parameters.posterior)

## Histograms
pdf( paste0( calib.folder, "/posterior_hist_hector.pdf" ) )
par(mfrow=c(3,3))
for ( pp in 1:length(parnames) ) {
        hist( parameters.posterior[,pp], xlab = params[pp], main='')
}
dev.off()

## Fit PDFs to the parameter distributions
pdf.all=vector('list',n.parameters)
n.node=200
for (pp in 1:n.parameters){
  tmp = density(parameters.posterior[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=params[pp]
}

##==============================================================================
## Plot the PDFs
pdf( paste0( calib.folder, "/posterior_pdf_hector.pdf" ) )
par(mfrow=c(3,3))
for ( pp in 1:n.parameters ){
  plot( pdf.all[[pp]]$x, pdf.all[[pp]]$y, type='l', xlab=params[pp], ylab="Density" );
}
dev.off()
##==============================================================================

## Write a CSV file with the successful parameter combinations
## Structure of the CSV file is as follows. n.parameters columns, 
##   (# of chains) x (# of burnt-in parameter sets / chain ) rows.
##     First row: Parameter names.
##     Rows 2-??: The calibrated parameter values.

to.file = parameters.posterior
rownames(to.file) = NULL
colnames(to.file) = params
filename = paste0( calib.folder, "/hector_calibrated_MCMC_parameters.csv")
write.table( to.file, file = filename, sep = ",", qmethod = "double", row.names = FALSE )

print( paste0( "Calibration completed. ", length(parameters.posterior[,1]),
               " posterior parameter sets saved to ", filename ) )
post_medians = rep( NA, n.parameters )
for (pp in 1:n.parameters){
  post_medians[pp] = median( parameters.posterior[,pp] )
}
print( "Hector parameter posterior medians: " )
print( params[in.hector] )
print( post_medians[in.hector] )
##==============================================================================
## End
##============================================================================== 
  



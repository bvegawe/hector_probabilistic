##==============================================================================
## hector_chain.R
## 
## Original code for BRICK written by Tony Wong and Alexander Bakker
##   (https://github.com/scrim-network/BRICK, DOI:10.5194/gmd-10-2741-2017)
## Modified for Hector by Ben Vega-Westhoff
##==============================================================================
## Because we read text file outputs from Hector, have to pass separate filename
## arguments to each parallel chain of our MCMC, thus the need for this code. 
##
## Implementation of MCMC (RAM) calibration of Hector simple climate model
## (https://github.com/JGCRI/hector, DOI:10.5194/gmd-8-939-2015)
##  
## Can add parameters to calibration in hector_calib_params.R.
## Can add/switch out observational constraints in obs_readData.R. 
##
## Robust adaptive MCMC (RAM) calibration for model parameters (Vihola, 2001)
## Differential evolution (DE) optimization to find suitable initial parameter guesses
## (Storn and Price, 1997; Price et al, 2005)
##
## Typical default use: $ Rscript hector_calib_driver.R -f default_calib
## See the input options: $ Rscript hector_calib_driver.R -h
##==============================================================================

rm(list=ls())	        # Clear all previous variables

library(optparse)
library(DEoptim)        # Use differential evolution for initial parameter guesses
library(adaptMCMC)      # Use robust adaptive Metropolis MCMC method

## Inputs
option_list = list(
    make_option( c("--chain")     , type="integer", default=1  ,
                 help = "Which chain is this in our parallel MCMC? (default = %default)",
		 metavar="integer" ), 
    make_option( c("--init")    , type="character", default=NULL       ,
		 help = "Which .RData file contains all the MCMC arguments and stuff? (default = %default)",
		 metavar="character" ) )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$init)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (path to the RData file with the MCMC args).",
    call.=FALSE)
}

## Set input variables
chain.mcmc = opt$chain	# RData file containing MCMC chain list to continue (we will set p0 from here)
init       = opt$init   # RData file containing the entire R state needed for these MCMC chains

load(init)
set.seed( 1234 + chain.mcmc )

source("../subparam.R")
source("../hectorwrapper.R")
source("../forcing_total.R")
source("../convertVars.R")
source("../hector_assimLikelihood.R")

##==============================================================================

## MCMC calibration
## + Cite Metropolis et al (1953) and Hasting (1973) for any Metropolis-Hastings
## + Also cite Vihola (2012) if you use "adaptMCMC" library's "MCMC" with adapt=TRUE
## + log.post is in the 'hector_assimLikelihood.R' module, and defines
##   the statistical model for calibration


print(paste0("Starting MCMC chain ",chain.mcmc))
t.beg=proc.time()
amcmc.out = MCMC( log.post                     , niter.mcmc                  , 
		     p0                           , scale = step.mcmc           , 
		     adapt = TRUE                 , acc.rate = accept.mcmc      ,
  		     gamma = gamma.mcmc           , list = TRUE                 , 
		     n.start = round(0.01*niter.mcmc),
		     parnames.in = parnames       , in.hector.in = in.hector    , 
		     sections.in = sections       , chain.str = chain.mcmc      , 
		     calib.folder = calib.folder  , forcing.in = forcing        ,
		     ini.template = ini.template  , output.vars = output.vars   ,
		     output.components = output.components,
		     mod.time = mod.time          , l.project = l.project       ,
		     bound.lower.in = bound.lower , bound.upper.in = bound.upper,
	             trends = trends              , oidx = oidx, midx = midx    , 
		     obs = obs, obs.err = obs.err , ind.norm.data = ind.norm.data)
t.end=proc.time()											# save timing
print(paste0(niter.mcmc," runs took: ",(t.end[3]-t.beg[3])/60.," min"))

## Save the MCMC object from this chain for later
save( amcmc.out, file = paste0( calib.folder, "/chain_", chain.mcmc, ".RData" ) )

##==============================================================================
## End
##==============================================================================

##==============================================================================
## hector_calib_parallel_init.R
## 
## Original code for BRICK written by Tony Wong and Alexander Bakker
##   (https://github.com/scrim-network/BRICK, DOI:10.5194/gmd-10-2741-2017)
## Modified for Hector by Ben Vega-Westhoff
##==============================================================================
## Set all relevant paramaters and determine initial parameter vals for chains.
##
## These are all passed to run the parallel chains in an:
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
## Typical default use: $ Rscript hector_calib_parallel_init.R -f default_calib
## See the input options: $ Rscript hector_calib_parallel_init.R -h
##==============================================================================

rm(list=ls())	        # Clear all previous variables

## Required packages/libraries
library(optparse)       # Option parsing
library(DEoptim)        # Use differential evolution for initial parameter guesses
library(adaptMCMC)	# Use robust adaptive Metropolis MCMC method

## Inputs
option_list = list(
    make_option( c("-f", "--folder"), type="character", default=NULL       ,
		 help = "folder for this calibration", metavar="character"),
    make_option( c("-n", "--niter") , type="integer"  , default=5E4        ,
		 help = "length of each chain (default= %default)",
	 	 metavar= "integer"),
    make_option( c("-g", "--gr")    , type="logical"  , default=TRUE       ,
		 help = "Use Gelman Rubin burn-in diagnostic? (default= %default)",
		 metavar= "logical"),
    make_option( c("--forcing")     , type="logical"  , default=FALSE      ,
		 help = "Do forcing-based runs, largely ignoring Hector's carbon cycle? (default= %default)",
		 metavar= "logical"), 
    make_option( c("--model_set")   , type="character", default="all_model",
		 help = "Which set of free model parameters for calibration? (default = %default, see hector_calib_params.R for options)",
                 metavar= "character"), 
    make_option( c("--obs_set")     , type="character", default="all_obs"  ,
                 help = "Which set of observation data sets against which to calibrate? (default = %default, see hector_calib_params.R for options)",
		 metavar= "character"), 
    make_option( c("--continue")    , type="character", default=NULL       ,
		 help = "MCMC list to continue, saved as a .RData file (default = %default)",
		 metavar= "character"), 
    make_option( c("--nparallel")   , type="integer"  , default=4          ,
                 help = "How many parallel MCMC chains? (default = %default)",
                 metavar= "integer"), 
    make_option( c("--name")        , type="character", default=NULL       ,
                 help = "Name for this MCMC calibration, in case of multiples running at once (default = %default)?",
                 metavar= "character") )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$folder)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (calibration folder).",
    call.=FALSE)
}

## Set input variables
calib.folder       = opt$folder 	# Folder for temporary files, diagnostic plots, and posterior parameters
forcing.based      = opt$forcing 	# If FALSE, use emissions-based Hector runs     
niter.mcmc         = opt$niter 		# Iterations per MCMC chain, 50k has been ample for DOECLIM parameters (based on GR diagnostic)
gr.mcmc            = opt$gr 		# If TRUE, perform Gelman Rubin burn-in diagnostic (a second chain is run for this if we only have one). If FALSE, assume the second half of the chain has converged (not necessarily a great assumption!).
model.set          = opt$model_set 	# Set of free parameters for the calibration. Set names are found in hector_calib_params.R.
obs.set            = opt$obs_set 	# Set of observation data sets against which to calibrate. Set names are found in hector_calib_params.R.
continue.mcmc      = opt$continue	# RData file containing MCMC chain list to continue (we will set p0 from here)
nparallel.mcmc     = opt$nparallel      # How many parallel chains?
name               = opt$name           # Name for this MCMC calibration? Needed if running more than one at a time.

today = Sys.Date(); today=format(today,format="%d%b%Y")
name = paste0(name,"_",today)

## Set the seed (for reproducibility)
set.seed(1234)

## Make projections? Only needed if calibrating against other projections (for emulation)
l.project=FALSE

## Create calibration folder
if( strsplit( calib.folder,"/" )[[1]][1] != "" ) {
    calib.folder = paste0( getwd(), "/", calib.folder ) # Convert to absolute path
}
if( !file.exists( calib.folder ) ) { system( paste0( "mkdir ", calib.folder ) ) }

## Get the forcing data (before possible aerosol scaling)
if( forcing.based ) {
    forcing.root = "obs_constraints/forcing/"
    file = ifelse( !l.project, "forcing_hindcast.csv",
           ifelse( scenario == 2.6, "forcing_rcp26.csv",
  	   ifelse( scenario == 4.5, "forcing_rcp45.csv",
  	   ifelse( scenario == 6.0, "forcing_rcp6.csv" ,
  	   ifelse( scenario == 8.5, "forcing_rcp85.csv")))))
    forcing.file = paste0( getwd(), "/", forcing.root, file )
    forcing = read.csv( forcing.file, header=TRUE )
    begyear = forcing$year[1]
    endyear = forcing$year[length(forcing$year)]
    print( "Doing calibration of forcing-based Hector" ) 
} else { forcing = NULL; begyear = endyear = NULL }

## Get emissions file
emissions.root = "obs_constraints/emissions/"
file = ifelse( !l.project, "hindcast_emissions.csv",
       ifelse( scenario == 2.6, "RCP26_emissions.csv",
       ifelse( scenario == 4.5, "RCP45_emissions.csv",
       ifelse( scenario == 6.0, "RCP6_emissions.csv" ,
       ifelse( scenario == 8.5, "RCP85_emissions.csv")))))
emissions.file = paste0( getwd(), "/", emissions.root, file )
volcanic.file = paste0( getwd(), "/", emissions.root, "volcanic_RF.csv" )
emissions = read.csv( emissions.file, header=TRUE, skip=3 )

## Set mod.time to include the forcing and/or emissions data needed to run Hector
if(is.null(begyear)){ begyear = emissions$Date[1] 
} else{ begyear = max(forcing$year[1],emissions$Date[1]) }
if(is.null(endyear)){ endyear = emissions$Date[length(emissions$Date)] 
} else{ endyear = min(forcing$year[length(forcing$year)],emissions$Date[length(emissions$Date)])}
mod.time = begyear:endyear
print( paste0( "Calibration will include Hector runs from ", begyear, " to ", endyear ) )

## Write a template .ini file for this calibration. 
## Individual parameters are updated in hectorwrapper().
## Here we just change things that are the same for the entire calibration.
source("../subparam.R") # Useful function to replace param values in .ini file
#if(forcing.based){ template.default = "input_templates/default_forcing_nowrite.ini"
#} else { template.default = "input_templates/default_nowrite.ini" }
if(forcing.based){ template.default = "input_templates/brick_forcing_nowrite.ini"
} else { template.default = "input_templates/brick_nowrite.ini" }
lines = readLines(template.default)
lines = gsub("obs_constraints/emissions/RCP45_emissions.csv", emissions.file, lines) #Replace default emissions file location
lines = gsub("obs_constraints/emissions/volcanic_RF.csv",volcanic.file, lines) #Replace default volcanic emissions file location
lines = subparam( lines, "core", "startDate", begyear - 1 ) #Hector's first calcs are for startDate+1
lines = subparam( lines, "core", "endDate", endyear )
lines = subparam( lines, "forcing", "baseyear", begyear )
ini.template = paste0( calib.folder, "/template.ini" )
cat( lines, file=ini.template, sep="\n" )

## Get the details for all of the free parameters in this calibration
source('../hector_calib_params.R')
param.list   = param_calib_details( model_set = model.set, obs_set = obs.set )
params       = param.list$params      ; parnames      = param.list$parnames   
sections     = param.list$sections    ; in.hector     = param.list$in.hector   
in.DEoptim   = param.list$in.DEoptim  ;
p0           = param.list$p0          ; bound.lower   = param.list$bound.lower         
bound.upper  = param.list$bound.upper ; step.mcmc     = param.list$step.mcmc
obs.ts       = param.list$obs.ts      ; output.vars   = param.list$compare.var 
output.components = param.list$compare.component
index.DEoptim = param.list$index.DEoptim

print( "Calibrating these parameters: " )
for( sec in unique(sections) ){ 
    print( paste( params[sections == sec], collapse = ", " ) ) 
}
print( "Using these observational constraints: " )
print( paste( obs.ts, collapse = ", " ) )

## Read in useful functions
source('../hectorwrapper.R')	# Runs Hector
source('../compute_indices.R')	# Computes the overlap indices of two time series
source('../forcing_total.R')	# Adjusts total forcing to account for any aerosol scaling
source('../convertVars.R')         # Function to convert model output to match observational output variables/units

## Read in all data sets against which we will calibrate
source('../obs_readData.R')
## Then gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions.
## Also gather actual observation/error values.
trends = list(); midx = list(); oidx = list(); obs = list(); obs.err = list()
norm.lower = list(); norm.upper = list()
for ( i in 1:length(obs.ts) ) {
    ts = obs.ts[[i]]
    trends[[ts]]     = obs.all[[ts]]$trends
    midx[[ts]]       = obs.all[[ts]]$midx
    oidx[[ts]]       = obs.all[[ts]]$oidx
    obs[[ts]]        = obs.all[[ts]]$obs
    obs.err[[ts]]    = obs.all[[ts]]$obs.err
    norm.lower[[ts]] = obs.all[[ts]]$norm.lower
    norm.upper[[ts]] = obs.all[[ts]]$norm.upper
}

## Which model indices should be used to normalize in same way as data?
l.idx = c(); u.idx = c()
for( i in 1:length(obs.ts) ) {
    l = norm.lower[[obs.ts[[i]]]]
    u = norm.upper[[obs.ts[[i]]]]
    if( !is.na(l) ) { l.idx = append( l.idx, which( mod.time == l ) ) 
    } else { l.idx = append( l.idx, NA ) }
    if( !is.na(u) ) { u.idx = append( u.idx, which( mod.time == u ) )
    } else { u.idx = append( u.idx, NA ) }
}
ind.norm.data = data.frame( obs.ts, l.idx, u.idx) 		

##==============================================================================

## Use 'DEoptim' (differential evolution optimization) to find better initial parameters
##   p0 initial parameter guesses based on Urban and Keller (2010)
##   These are okay, and work, but can improve using differential evolution optimization
##   (as long as you use a large enough vector population (at least 10*[# parameters]))

#Can skip to go straight to mcmc, instead using the best values from a previous DEoptim run (see commented example below)
if(is.null(continue.mcmc)){
    print( "Starting DEoptim (differential evolution optimization) to find initial parameters for the MCMC chain" )
    source('../hector_DEoptim.R')

    niter.deoptim=200          # number of iterations for DE optimization
    NP.deoptim=10*length(p0[index.DEoptim])   # population size for DEoptim (do at least 10*[N parameters])
    F.deoptim=0.8              # as suggested by Storn et al (2006)
    CR.deoptim=0.9             # as suggested by Storn et al (2006)

    t.beg=proc.time()          # save timing
    outDEoptim <- DEoptim( minimize_residuals_hector            , 
			   bound.lower[index.DEoptim]           , bound.upper[index.DEoptim],
   		           DEoptim.control(NP = NP.deoptim, itermax = niter.deoptim, F = F.deoptim, 
		                           CR = CR.deoptim, trace = FALSE)                               , 
			   parnames.in = parnames[index.DEoptim], in.hector.in = in.hector[index.DEoptim],
			   sections.in = sections[index.DEoptim], calib.folder = calib.folder            ,
			   forcing.in = forcing                 , ini.template = ini.template            ,
			   output.vars = output.vars            , output.components = output.components  ,
			   mod.time = mod.time                  , l.project = l.project                  , 
		           trends = trends                      , oidx = oidx, midx = midx               , 
		           obs = obs, obs.err = obs.err         , ind.norm.data = ind.norm.data          )

    p0[index.DEoptim] = outDEoptim$optim$bestmem
    print("Finished DEoptim")
    t.end=proc.time()          # save timing
    print(paste0("DEoptim took: ",(t.end[3]-t.beg[3])/60.," min"))
    print(paste0(params," ",p0))
}

if(!is.null(continue.mcmc)){

    # 12/27/17 Forcing-based Hector. Took 170 minutes. DOECLIM DEoptim output.
#   p0[index.DEoptim] = c( 1.91804058912437, 1.53993386594774, 0.520910696383523, 
#                          -0.00949539681064582, -20.1757086003454 )

    # 12/28/17 Emissions-based Hector. Took 170 minutes. DOECLIM DEoptim output.
#    DEoptim.output.str = "Emissions-based Hector, DOECLIM parameters and observation constraints"
#    p0[index.DEoptim] = c( 2.01823906164744, 2.99881682965184, 0.50315986159485, 
#		         0.0289778671684981, -27.8713241773148 )

    print(paste0("Skipping DEoptim, continuing first chain from ",continue.mcmc) )
    load(continue.mcmc)
    p0 = unconverged_chain_list[[1]][length(unconverged_chain_list[[1]][,1]),]
    print(paste0(params," ",p0))
}
##==============================================================================

## Getting ready for MCMC calibration
## + Cite Metropolis et al (1953) and Hasting (1973) for any Metropolis-Hastings
## + Also cite Vihola (2012) if you use "adaptMCMC" library's "MCMC" with adapt=TRUE
## + log.post is in the 'hector_assimLikelihood.R' module, and defines
##   the statistical model for calibration

## Source the statistical model
source('../hector_assimLikelihood.R')

accept.mcmc = 0.234   		     # Optimal acceptance rate as # parameters->infinity (Gelman et al, 1996; Roberts et al, 1997)
gamma.mcmc = 0.5	             # rate of adaptation (between 0.5 and 1, lower is faster adaptation)
stopadapt.mcmc = round(niter.mcmc*1) # stop adapting after how long (if ever)?

#Load all of this stuff for every parallel chain of the MCMC
filename.mcmc_init <- paste0(calib.folder,"/hector_calib_MCMC_ready.RData")
save.image(file = filename.mcmc_init)
##==============================================================================
## Write batch scripts to run each of the chains 
## (this is specific to my UIUC department cluster and would need 
##  tweaking for another cluster)

batch_template = "hector_chain_template.batch"
for (i in 1:nparallel.mcmc){
    lines = readLines(batch_template)
    lines = gsub( "hcal_c0", paste0( "hcal_c", i ), lines )
    lines = gsub( "filename.mcmc_init", filename.mcmc_init, lines )
    lines = gsub( "--chain 0", paste0( " --chain ", i ), lines )
    chain.file = paste0(name,"_",i,"chain.batch")
    cat( lines, file=chain.file, sep="\n" )
}    
##==============================================================================
## End
##==============================================================================

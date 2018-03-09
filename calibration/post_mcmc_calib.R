##==============================================================================
## post_mcmc_calib.R
##
## Original code for BRICK written by Tony Wong and Alexander Bakker
##   (https://github.com/scrim-network/BRICK, DOI:10.5194/gmd-10-2741-2017)
## Modified for Hector by Ben Vega-Westhoff
##==============================================================================
##
## Pipeline for processing DAIS calibration results and Hector rest-of-model
## calibration results.
##
## 1. Data are read, hindcasts are set up. Parameters for the DAIS paleoclimatic
##              calibration and the rest-of-model modern calibration are drawn. The
##              number of parameter combinations (initial ensemble members) is specified
##              in the section of the script for the user to modify. (n.ensemble)
## 2. Hector model parameter sets run the full model to obtain hindcasts
##              of global mean sea level. The parameters are calibrated to global mean sea
##              level data (Church and White, 2011). This is done using rejection sampling.
## 3. These fully calibrated parameter sets are written to a netCDF file whose
##              name is given by [filename.parameters], specified by the user.
##
##  Typical default use: $ Rscript post_mcmc_calib.R -f *calib_folder*
##==============================================================================

rm(list=ls())

## Required packages/libraries
library(optparse)       # Option parsing
library(ncdf4)

## Inputs
option_list = list(
    make_option( c("-f", "--folder")  , type="character", 
                 default=NULL                                                      ,
                 help="existing calibration folder location"                       , 
                 metavar="character") ,
    make_option( c("-d", "--dais")    , type="character", 
                 default="DAIS_calibratedParameters_18Oct2017.nc"                  ,
                 help ="DAIS calibrated params file (default= %default)"           , 
                 metavar="character") ,
    make_option( c("-w", "--working") , type="character", 
                 default="../calibration"                                          ,
                 help="working drectory w/ all calib functions (default= %default)",
                 metavar="character") , 
    make_opion( c("-n", "--nensemble"), type="integer"  , 
                 default=20000                                                     ,
                 help="# of draws from our posteriors (default= %default)"         ,
                 metavar="integer")   )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$folder)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (calibration folder).",
    call.=FALSE)
}

## Set input variables
calib.folder       = opt$folder         # Where are our calibration results and where will we save?
                                        # Required files in this folder: 
					#    hector_calib_MCMC_ready.RData
                                        #    hector_calibrated_MCMC_parameters.csv
dais.file          = opt$dais           # Calibrated DAIS parameters file, in netCDF format
working.folder     = opt$working        # Where are all of the calibration functions?
n.ensemble         = opt$nensemble      # How many draws will we make from our posteriors?

if( strsplit( calib.folder,"/" )[[1]][1] != "" ) {
    calib.folder = paste0( getwd(), "/", calib.folder ) # Convert to absolute path
}

load( paste0 ( calib.folder, "/hector_calib_MCMC_ready.RData" ) #All the parameter set-up is in here

## Read calibrated Hector parameter sets
parameters.hector = read.csv( paste0( calib.folder, "hector_calibrated_MCMC_parameters.csv" ), header=TRUE )
print( paste0( 'read ', length( parameters.hector[[1]] ), ' calibrated Hector model parameter sets' ) )

## Read calibrated DAIS parameter sets
ncdata <- nc_open( dais.file )
parameters.dais = ncvar_get( ncdata, 'DAIS_parameters' )
parnames.dais = ncvar_get( ncdata, 'parnames' )
nc_close(ncdata)
parameters.dais = t(parameters.dais)
colnames(parameters.dais) = parnames.dais
print( paste0( 'read ', nrow( parameters.dais ), ' calibrated DAIS model parameter sets') )

## Make sure we aren't drawing more samples than we have
n.ensemble.in = n.ensemble
n.ensemble = min( c( n.ensemble, nrow( parameters.dais ), length( parameters.hector[[1]] ) ) )
if ( n.ensemble < n.ensemble.in ) {
    print( 'Reduced sample size to match the posterior chain length' )
    print( paste0( 'Input: ', n.ensemble.in, ', actual: ', n.ensemble ) )
}

## Draw parameters for calibrated DAIS model
parameters.ensemble.dais = mat.or.vec(n.ensemble , ncol(parameters.dais) )
ind.dais=sample( seq(1,nrow(parameters.dais)), size=n.ensemble, replace=FALSE)
for( p in 1:ncol(parameters.dais) ) {
    for( i in 1:n.ensemble ) {
        parameters.ensemble.dais[i,p] = parameters.dais[ind.dais[i],p]
    }
    
    sections.dais = rep( 'slr_brick', ncol(parameters.dais) )
    in.hector.dais = rep( TRUE, ncol(parameters.dais) )
}


## Draw parameters for calibrated Hector model
print(paste('Creating possible parameter combinations for calibrated models...'))
parameters.ensemble.hector = mat.or.vec(n.ensemble , ncol(parameters.hector) )
ind.ensemble = sample( seq(1,nrow(parameters.hector)), size=n.ensemble, replace=FALSE)
for( p in 1:ncol(parameters.hector) ) {
    for( i in 1:n.ensemble ) {
        parameters.ensemble.hector[i,p] = parameters.hector[ind.ensemble[i],p]
    }
}
print(paste(' ... done creating parameter sets'))

## Add DAIS parameter values/details to those from the Hector MCMC
parameters = cbind(parameters.ensemble.hector, parameters.ensemble.dais)
parnames = c(parnames, parnames.dais)
sections = c(sections, sections.dais)
in.hector = c(in.hector, in.hector.dais)
rownames(parameters)=NULL

## Prepare for hindcasts
original.wd = getwd()
setwd(working.folder)
source('subparam.R')
source('hectorwrapper.R')
source('forcing_total.R')
source('convertVars.R')

## Initialize matrix to store model ensemble output
hector.out = vector("list", n.ensemble)

## Initialize flag for possibly bad runs (DAIS+Hector parameters could go wrong,
## because the other model components were calibrated without DAIS, and vice
## versa)
badruns = rep(0, n.ensemble)

## Run the sample, with a progress bar!
print( paste0( 'Starting ', n.ensemble, ' model hindcasts' ) )
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for( i in 1:n.ensemble ) {
    parvals = as.numeric( parameters[i,] )

    ## Scale forcing if running in forcing mode
    if(!is.null(forcing)){
        alpha = parvals[match("alphatemperature",paste0(parnames,sections))]
        if(is.na(alpha)){ alpha = 1 }
        forcing.total = forcing_total(  forcing  =forcing.in, alpha.doeclim =alpha      ,
                                        l.project=l.project , begyear       =mod.time[1],
                                        endyear  =mod.time[length(mod.time)])
        df_out = data.frame(year=mod.time,Ftot_constrain=forcing.total)
        forcing.file = paste0( calib.folder,"/DEoptim_temp_RF.csv" )
        write.table(df_out,file=forcing.file, quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
        parvals[match("alphatemperature",paste0(parnames,sections))] = 1 #Forcing already scaled, no need to do it in Hector
    }
    else {
        forcing.file = NULL
    }
    
        model.out = hectorwrapper( parnames[in.hector],
                                   parvals[in.hector],
                                   sections[in.hector],
                                   chain.str=chain.str,
                                   working.dir=calib.folder,
                                   forcing.file=forcing.file,
                                   ini.template=ini.template,
                                   output.vars=output.vars,
                                   output.components=output.components,
                                   mod.time )












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
                 default="separate_calib_output/DAIS_calibratedParameters_18Oct2017.nc"                  ,
                 help ="DAIS calibrated params file (default= %default)"           , 
                 metavar="character") ,
    make_option( c("-w", "--working") , type="character", 
                 default="../calibration"                                          ,
                 help="working directory w/ all calib functions (default= %default)",
                 metavar="character") , 
    make_option( c("-n", "--nensemble"), type="integer"  , 
                 default=200000                                                     ,
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
					#    hector_calib_after_mcmcs.RData
                                        #    hector_calibrated_MCMC_parameters.csv
dais.file          = opt$dais           # Calibrated DAIS parameters file, in netCDF format
working.folder     = opt$working        # Where are all of the calibration functions?
n.ensemble         = opt$nensemble      # How many draws will we make from our posteriors?

if( strsplit( calib.folder,"/" )[[1]][1] != "" ) {
    calib.folder = paste0( getwd(), "/", calib.folder ) # Convert to absolute path
}

filename.saveprogress.post = paste0( calib.folder, "/hector_calib_after_rs.RData" ) 
                        #Will contain images as we run thru this code

load( paste0 ( calib.folder, "/hector_calib_after_mcmcs.RData" ) ) #All the parameter set-up is in here

## Read calibrated Hector parameter sets
parameters.hector = read.csv( paste0( calib.folder, "/hector_calibrated_MCMC_parameters.csv" ), header=TRUE )
print( paste0( 'read ', length( parameters.hector[[1]] ), ' calibrated Hector model parameter sets' ) )

## Read calibrated DAIS parameter sets
ncdata <- nc_open( dais.file )
parameters.dais = ncvar_get( ncdata, 'DAIS_parameters' )
parnames.dais = ncvar_get( ncdata, 'parnames' )
nc_close(ncdata)
parameters.dais = t(parameters.dais)
parameters.dais = parameters.dais[,1:(ncol(parameters.dais)-1)] #drop the var.dais parameter that is only for likelihood calculations
colnames(parameters.dais) = parnames.dais[1:(length(parnames.dais)-1)]
print( paste0( 'read ', nrow( parameters.dais ), ' calibrated DAIS model parameter sets') )

## Make sure we aren't drawing more samples than we have
n.ensemble.in = n.ensemble
n.ensemble = min( c( n.ensemble, nrow( parameters.dais ), length( parameters.hector[[1]] ) ) )
if ( n.ensemble < n.ensemble.in ) {
    print( 'Reduced sample size to match the posterior chain length' )
    print( paste0( 'Input: ', n.ensemble.in, ', actual: ', n.ensemble ) )
}

## Draw parameters for calibrated DAIS model
ind.dais=sample( seq(1,nrow(parameters.dais)), size=n.ensemble, replace=FALSE)
parameters.ensemble.dais = data.matrix(parameters.dais[ind.dais,])
    
sections.dais = rep( 'slr_brick', ncol(parameters.dais) )
in.hector.dais = rep( TRUE, ncol(parameters.dais) ) 

## It's ugly...but adding these ranges in by brute force. They are specified in the DAIS calibration file:
## https://github.com/scrim-network/BRICK/blob/master/calibration/DAIS_calib_driver.R
bound.lower.dais = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 735.5, 47.5, 740 , 0.00045, 0.005  , -20 )
bound.upper.dais = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 ,2206.5,142.5, 820 , 0.00075, 0.015  , -10 )
parnames.dais = c( "a_anto", "b_anto", "gamma_dais", "alpha_dais", "mu_dais", "nu_dais", 
                       "P0_dais", "kappa_dais", "f0_dais", "h0_dais", "c_dais",
                       "b0_dais", "slope_dais", "lambda_dais", "Tcrit_dais" )



## Draw parameters for calibrated Hector model
print(paste('Creating possible parameter combinations for calibrated models...'))
ind.ensemble = sample( seq(1,nrow(parameters.hector)), size=n.ensemble, replace=FALSE)
parameters.ensemble.hector = data.matrix(parameters.hector[ind.ensemble,])
print(paste(' ... done creating parameter sets'))

## Add DAIS parameter values/details to those from the Hector MCMC
parameters = cbind(parameters.ensemble.hector, parameters.ensemble.dais)
parnames = c(parnames, parnames.dais)
params = c(params, paste0(parnames.dais,".slr_brick"))
sections = c(sections, sections.dais)
in.hector = c(in.hector, in.hector.dais)
bound.lower = c( bound.lower, bound.lower.dais )
bound.upper = c( bound.upper, bound.upper.dais )
rownames(parameters)=NULL
colnames(parameters)=params

## Add necessary slr components if not part of the original calibration
output.vars.orig = output.vars
output.components.orig = output.components
slr_vars = c("slr_gis","slr_gsic","slr_te","slr_ais","slr")
slr_components = c("slr_brick","slr_brick","slr_brick","slr_brick","slr_brick")
for( var in slr_vars){
  if( !(var %in% output.vars) ){
    output.vars = c(output.vars,var)
  }
}
for( comp in slr_components){
  if( !(comp %in% output.components) ){
    output.components = c(output.components,comp)
  }
}

## Prepare for hindcasts
original.wd = getwd()
setwd(working.folder)
source('subparam.R')
source('hectorwrapper.R')
source('forcing_total.R')
source('convertVars.R')
source('ar1_sim.R')

ntime = length(mod.time)
slr_tot    = mat.or.vec(n.ensemble,length(mod.time))

## Run the sample, with a progress bar!
print( paste0( 'Starting ', n.ensemble, ' model hindcasts' ) )
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for( i in 1:n.ensemble ) {
    parvals = as.numeric( parameters[i,] )

    ## Scale forcing if running in forcing mode
    if(!is.null(forcing)){
        alpha = parvals[match("alphatemperature",paste0(parnames,sections))]
        if(is.na(alpha)){ alpha = 1 }
        forcing.total = forcing_total(  forcing  =forcing, alpha.doeclim =alpha      ,
                                        l.project=l.project , begyear       =mod.time[1],
                                        endyear  =mod.time[length(mod.time)])
        df_out = data.frame(year=mod.time,Ftot_constrain=forcing.total)
        forcing.file = paste0( calib.folder,"/rejsample_temp_RF.csv" )
        write.table(df_out,file=forcing.file, quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
        parvals[match("alphatemperature",paste0(parnames,sections))] = 1 #Forcing already scaled, no need to do it in Hector
    } else {
        forcing.file = NULL
    }
    hector.out = hectorwrapper( parnames[in.hector],
                                parvals[in.hector],
                                sections[in.hector],
                                chain.str="rs",
                                working.dir=calib.folder,
                                forcing.file=forcing.file,
                                ini.template=ini.template,
                                output.vars=output.vars,
                                output.components=output.components,
                                mod.time )

    ## Before post-calibration, need to add the modeled statistical noise back in.
    ## Only using sea-level rise data, so only need to modify GSIC, GIS.
    ## using the statistical parameters for AR1, AR1 and Gaussian noise, respectively
    ## Do not do for AIS, because var.dais was fit to paleo data-model mismatch, not
    ## representative of the current era.

    #GSIC - note, obs.err for gsic is much larger than its year-to-year variability,
    #       therefore, not including it in the statistical noise, for now.
    sigma.gsic = parameters[i,match("sigmaslr_gsic_obs",paste0(parnames,sections))]
    rho.gsic   = parameters[i,match("rhoslr_gsic_obs"  ,paste0(parnames,sections))]
    err.gsic = rep(sigma.gsic,ntime)
    err.gsic[midx[["slr_gsic"]]] = sqrt( sigma.gsic^2 + obs.err[["slr_gsic"]]^2 )
    gsic = hector.out$slr_gsic + ar1.sim( ntime, rho.gsic, sigma.gsic )

    #GIS
    sigma.gis  = parameters[i,match("sigmaslr_gis_obs",paste0(parnames,sections))]
    rho.gis    = parameters[i,match("rhoslr_gis_obs"  ,paste0(parnames,sections))]
    err.gis    = rep(sigma.gis,ntime)
    err.gis[midx[["slr_gis"]]] = sqrt( sigma.gis^2 + obs.err[["slr_gis"]]^2 )
    gis = hector.out$slr_gis + ar1.sim( ntime, rho.gis, err.gis )

    #Other components of the total SLR don't have these associated errors (limited obs)
    slr_tot[i,] = gsic + gis + hector.out$slr_te + hector.out$slr_ais 
    
    #Normalize to match Church and White obs
    norm.lower = church_slr_obs$norm.lower; norm.upper = church_slr_obs$norm.upper
    slr_tot[i,] = slr_tot[i,] - mean(slr_tot[i,which(mod.time == norm.lower):which(mod.time == norm.upper)]) 

    setTxtProgressBar(pb, i)  
}
close(pb)
print( "Done running hindcats and adding up and adding noise to hindcast total slr." )
save.image( file = filename.saveprogress.post )

## Perform rejection sampling
source("rejection_sample_wLW.R")
if(!is.null(perfFile)){
  ind.survive = rejection_sample_wLW( perf.all$slr$obs,
                                      perf.all$slr$obs.time,
                                      perf.all$slr$obs.err,
                                      slr_tot  )
} else {
  ind.survive = rejection_sample_wLW( church_slr_obs$obs, 
                                      church_slr_obs$obs.time, 
                                      church_slr_obs$obs.err,
                                      slr_tot )
}         
parameters.good = parameters[ind.survive,]
colnames(parameters.good) = parnames

save.image( file = filename.saveprogress.post )

## Fit PDFs to the parameter distributions
n.parameters = length(parnames)
pdf.rs.all = vector( 'list', n.parameters )
n.node = 200
for (pp in 1:n.parameters) {
  tmp = density(parameters.good[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.rs.all[[pp]] = tmp; names(pdf.rs.all)[pp]=params[pp]
}

## Also for the mcmcs ones? Should be in the .RData file we read - try deleting this later?
pdf.all=vector('list',n.parameters)
n.node=200
for (pp in 1:ncol(parameters.hector)){
  tmp = density(parameters.hector[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=params[pp]
}
## Also for the dais mcmc
for(pp in 1:(ncol(parameters.dais)-1)){
  tmp = density(parameters.dais[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[ncol(parameters.hector)+pp]] = tmp
  names(pdf.all)[ncol(parameters.hector)+pp] = parnames.dais[pp]
}

## Plot the PDFs, comparing before vs after rejection sampling
pdf( paste0( calib.folder, "/pdfs_rejectionSampling.pdf" ) )
par(mfrow=c(3,3))
for (pp in 1:n.parameters ) {
  plot( pdf.rs.all[[pp]]$x, pdf.rs.all[[pp]]$y, type='l', xlab=params[pp], ylab="Density" )
  lines( pdf.all[[pp]]$x, pdf.all[[pp]]$y, lty=2 )
  legend( "topright", c("After RS", "Before"), lty=c(1,2) ) 
}

## Write post-calibrated parameters to output csv file
filename = paste0( calib.folder, "hector_postcalibrated_parameters.csv" )
to.file = parameters.good
rownames(to.file) = NULL
colnames(to.file) = params

write.table( to.file, file = filename, sep = ",", qmethod = "double", row.names = FALSE )
print( paste0( "Rejection sampling completed. ", length(ind.survive),
               " final parameter sets saved to ", filename ) )
save.image( file = filename.saveprogress.post )
##==============================================================================
## End
##==============================================================================

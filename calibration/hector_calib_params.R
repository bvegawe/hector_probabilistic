##==============================================================================
## hector_calib_params.R
##
## Ben Vega-Westhoff
##==============================================================================
## Specify parameters to be calibrated (using hector_calib_driver.R).
## Additional uncertain parameters can be added, just make sure they are also
##  added to the appropriate parameter lists.
##
## Required details for each calibrated parameter:
## param       - A unique parameter name
## parname     - The parameter's name in Hector
## component   - In which Hector component?
## in.hector   - Is this actually a parameter used in Hector? 
##               (parameters used to match the observations are not) ? 
## in.DEoptim  - Is this a free parameter in the differential evolution 
##		 optimization (pretty much everything except for the 
##		 autoregression observation-match parameters) ?
## p0	       - Best guess for the value (most are updated in DEoptim before MCMC)
## bound.lower - Value lower bound
## bound.upper - Value upper bound
## step.mcmc   - MCMC proposal step sizes
##
## Note: Each observational time series constraint should have 
##       three associated free parameters:
##
##	 offset - The offset at the first time step, after any normalization
##	 sigma  - The year-to-year internal variability of the time series
##	 rho    - The lag-1 autocorrelation of the time series	 	 
##
## param_calib_details:
## Input:
##  model_params     - A list of model parameters to calibrate
##  model_components - If no model_params, a list of model components to 
##			 calibrate
##  model_set        - If no model_params or model_components, a predefined
##			 set of parameters to calibrate
##  obs_ts	     - A list of observation time series against which to 
##			 calibrate
##  obs_set	     - If no obs_ts, a predefined set of observations against 
##			 which to calibrate
## 
## Output:
##  A list of the details of all of these parameters
##
##==============================================================================
## Hector uncertain parameters
##==============================================================================
## Uncertain parameters in the temperature component
S.temperature     = list( param = "S.temperature"    , parname = "S"    , component = "temperature", in.hector = TRUE, in.DEoptim = TRUE, p0 = 3.1, bound.lower = 0.1, bound.upper = 10, step.mcmc = 0.16  )
diff.temperature  = list( param = "diff.temperature" , parname = "diff" , component = "temperature", in.hector = TRUE, in.DEoptim = TRUE, p0 = 3.5, bound.lower = 0.1, bound.upper = 4 , step.mcmc = 0.17  )
alpha.temperature = list( param = "alpha.temperature", parname = "alpha", component = "temperature", in.hector = TRUE, in.DEoptim = TRUE, p0 = 1.1, bound.lower = 0  , bound.upper = 2 , step.mcmc = 0.025 )
parlist.temperature = list( S.temperature, diff.temperature, alpha.temperature )

## Uncertain parameters in the slr_brick component
## Glaciers and small ice caps subcomponent (GSIC-MAGICC)
beta0_gsic.slr_brick = list( param = "beta0_gsic.slr_brick", parname = "beta0_gsic", component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 0.00058, bound.lower = 0      , bound.upper = 0.01 , step.mcmc = 0.01 )
V0_gsic.slr_brick    = list( param = "V0_gsic.slr_brick"   , parname = "V0_gsic"   , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 0.41   , bound.lower = 0.3    , bound.upper = 0.5   , step.mcmc = 0.01 )
n_gsic.slr_brick     = list( param = "n_gsic.slr_brick"    , parname = "n_gsic"    , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 0.82   , bound.lower = 0.55   , bound.upper = 1.0   , step.mcmc = 0.1  )
Gs0_gsic.slr_brick   = list( param = "Gs0_gsic.slr_brick"  , parname = "Gs0_gsic"  , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 0.0    , bound.lower = -0.0041, bound.upper = 0.0041, step.mcmc = 0.01 )
## Explicit thermal expansion subcomponent (TEE)
a_tee.slr_brick      = list( param = "a_tee.slr_brick"     , parname = "a_tee"     , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 0.16   , bound.lower = 0.05   , bound.upper = 0.3   , step.mcmc = 0.05 * ( 0.3 - 0.05 ) )
## Greenland ice sheet subcomponent (GIS/SIMPLE)
a_simple.slr_brick     = list( param = "a_simple.slr_brick"     , parname = "a_simple"    , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = -0.825 , bound.lower = -4  , bound.upper = -1e-3, step.mcmc = 0.2  )
b_simple.slr_brick     = list( param = "b_simple.slr_brick"     , parname = "b_simple"    , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 7.36   , bound.lower = 0   , bound.upper = 8.832, step.mcmc = 0.05 )
alpha_simple.slr_brick = list( param = "alpha_simple.slr_brick" , parname = "alpha_simple", component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 1.63e-5, bound.lower = 0   , bound.upper = 1e-3 , step.mcmc = 1e-5 )
beta_simple.slr_brick  = list( param = "beta_simple.slr_brick"  , parname = "beta_simple" , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 2.85e-5, bound.lower = 0   , bound.upper = 1e-3 , step.mcmc = 1e-5 )
V0_simple.slr_brick    = list( param = "V0_simple.slr_brick"    , parname = "V0_simple"   , component = "slr_brick", in.hector = TRUE, in.DEoptim = TRUE, p0 = 7.36   , bound.lower = 7.16, bound.upper = 7.56 , step.mcmc = 0.05 )
parlist.slr_brick = list( beta0_gsic.slr_brick, V0_gsic.slr_brick, n_gsic.slr_brick, Gs0_gsic.slr_brick,
                          a_tee.slr_brick,
                          a_simple.slr_brick, b_simple.slr_brick, alpha_simple.slr_brick, beta_simple.slr_brick, V0_simple.slr_brick )

## Add uncertain parameters in other components here

## !!!Always update parlist.all_model to include all of the free parameters in Hector!!!
parlist.all_model = c( parlist.temperature, parlist.slr_brick )

## Additional predefined model parameter sets, called with "model_set" in param_calib_details()
## an example: just DOECLIM parameters: S.temperature, diff.temperature, alpha.temperature
parlist.doeclim_model = list( S.temperature, diff.temperature, alpha.temperature )
parlist.noTE_model = list( S.temperature, diff.temperature, alpha.temperature,
                           beta0_gsic.slr_brick, V0_gsic.slr_brick, n_gsic.slr_brick,
                           Gs0_gsic.slr_brick, a_simple.slr_brick, b_simple.slr_brick,
                           alpha_simple.slr_brick, beta_simple.slr_brick, V0_simple.slr_brick )
parlist.onlyTE_model = list( S.temperature, diff.temperature, alpha.temperature,
                             a_tee.slr_brick ) 
parlist.onlyT_model = list( S.temperature, diff.temperature, alpha.temperature)
## Add more sets as desired


##==============================================================================
## Parameters used to match the observation time series
##==============================================================================
## Free parameters of the globally averaged surface temperature time series
offset.Tgav_obs    = list( param = "offset.Tgav_obs", parname = "offset", component = "Tgav_obs", in.hector = FALSE, in.DEoptim = TRUE , p0 = -0.06, bound.lower = -0.3, bound.upper = 0.3 , step.mcmc = 0.003 )
sigma.Tgav_obs     = list( param = "sigma.Tgav_obs" , parname = "sigma" , component = "Tgav_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.1  , bound.lower = 0.05, bound.upper = 5   , step.mcmc = 5e-4  )
rho.Tgav_obs       = list( param = "rho.Tgav_obs"   , parname = "rho"   , component = "Tgav_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.55 , bound.lower = 0   , bound.upper = 0.99, step.mcmc = 0.007 )
parlist.Tgav_obs   = list( offset.Tgav_obs   , sigma.Tgav_obs   , rho.Tgav_obs )
compare.Tgav_obs   = list( var = "Tgav", component = "temperature", obsvar = "Tgav" ) # Which Hector variable are we using for comparing?

## Free parameters of the ocean heat time series
# Note: because we don't normalize ocean heat content, its offset relative to the observations 
# (which start in mid-20th century) will depend on Hector's start-date. 
# That means p0 and bounds may need to be changed. 
# These initial values were originally for runs starting in 1850.
offset.ocheat_obs  = list( param = "offset.ocheat_obs", parname = "offset", component = "ocheat_obs", in.hector = FALSE, in.DEoptim = TRUE , p0 = -33, bound.lower = -50, bound.upper = 0   , step.mcmc = 0.9   )
sigma.ocheat_obs   = list( param = "sigma.ocheat_obs" , parname = "sigma" , component = "ocheat_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 2  , bound.lower = 0.1, bound.upper = 10  , step.mcmc = 0.025 )
rho.ocheat_obs     = list( param = "rho.ocheat_obs"   , parname = "rho"   , component = "ocheat_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.9, bound.lower = 0  , bound.upper = 0.99, step.mcmc = 0.006 )
parlist.ocheat_obs = list( offset.ocheat_obs , sigma.ocheat_obs , rho.ocheat_obs )
compare.ocheat_obs = list( var = "heatflux", component = "temperature", obsvar = "ocheat" ) # Which Hector variable are we using for comparison?

## Free parameters of the glaciers and small ice caps SLR contribution time series
# Note: Observations are normalized to the value in 1960 - no need for an offset parameter
sigma.slr_gsic_obs = list( param = "sigma.slr_gsic_obs", parname = "sigma", component = "slr_gsic_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.00045, bound.lower = 0.0   , bound.upper = 0.00150, step.mcmc = 0.0001 )
rho.slr_gsic_obs   = list( param = "rho.slr_gsic_obs"  , parname = "rho"  , component = "slr_gsic_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.5    , bound.lower = -0.999, bound.upper = 0.999  , step.mcmc = 0.01   )
parlist.slr_gsic_obs = list( sigma.slr_gsic_obs, rho.slr_gsic_obs )
compare.slr_gsic_obs = list( var = "slr_gsic", component = "slr_brick",  obsvar = "slr_gsic" ) # Which Hector variable are we using for comparing? 
## No free parameters of the thermal expansion time series currently, we just compare trends in two time periods
compare.slr_te_obs = list( var = "slr_te", component = "slr_brick", obsvar = "slr_te" )
## Free parameters of the Greenland ice sheet SLR contribution time series
# Note: I guess V0.simple takes any offset into account. Obs are normalized to 1961-1990
sigma.slr_gis_obs = list( param = "sigma.slr_gis_obs", parname = "sigma", component = "slr_gis_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 5e-4, bound.lower = 0     , bound.upper = 0.002, step.mcmc = 0.0001 )
rho.slr_gis_obs   = list( param = "rho.slr_gis_obs"  , parname = "rho"  , component = "slr_gis_obs", in.hector = FALSE, in.DEoptim = FALSE, p0 = 0.5 , bound.lower = -0.999, bound.upper = 0.999, step.mcmc = 0.1    )
parlist.slr_gis_obs = list( sigma.slr_gis_obs, rho.slr_gis_obs )
compare.slr_gis_obs = list( var = "slr_gis", component = "slr_brick", obsvar = "slr_gis" ) # Which Hector variable are we using for comparison?

## Add free parameters used to match additional observation time series here. component should end with "_obs"
parlist.all_obs = c( parlist.Tgav_obs, parlist.ocheat_obs, 
                     parlist.slr_gsic_obs, parlist.slr_gis_obs )
## !!!Always update compare.all_obs to include all of the Hector outputs needed for comparison!!!
compare.all_obs  = list( compare.Tgav_obs, compare.ocheat_obs, 
                         compare.slr_gsic_obs, compare.slr_te_obs, compare.slr_gis_obs )

#Everything but ocean heat observational constraints. 
#(this avoids double counting of ocean heat info from thermal expansion and ocean heat)
#  !!!Always update this as well because it's the default
parlist.noOcheat_obs = c( parlist.Tgav_obs, parlist.slr_gsic_obs,
                          parlist.slr_gis_obs )
compare.noOcheat_obs = list( compare.Tgav_obs, compare.slr_gsic_obs,
                             compare.slr_te_obs, compare.slr_gis_obs )

## Additional predefined observation parameter sets, called with "obs_set" in param_calib_details()
## An example: DOECLIM-related obs: just temperature and ocean heat
parlist.doeclim_obs = c( parlist.Tgav_obs, parlist.ocheat_obs )
compare.doeclim_obs = list( compare.Tgav_obs, compare.ocheat_obs )

parlist.noTE_obs = c( parlist.Tgav_obs, parlist.ocheat_obs,
                      parlist.slr_gsic_obs, parlist.slr_gis_obs )
compare.noTE_obs = list( compare.Tgav_obs, compare.ocheat_obs,
                         compare.slr_gsic_obs, compare.slr_gis_obs )

parlist.onlyTE_obs = c( parlist.Tgav_obs, parlist.ocheat_obs )
compare.onlyTE_obs = list( compare.Tgav_obs, compare.ocheat_obs,
                           compare.slr_te_obs )

#Only temperature and thermosteric slr
parlist.doeclimTE_obs = c( parlist.Tgav_obs )
compare.doeclimTE_obs = list( compare.Tgav_obs, compare.slr_te_obs )

#Only temperature
parlist.onlyT_obs = c( parlist.Tgav_obs )
compare.onlyT_obs = list( compare.Tgav_obs )
## Add more sets as desired

##==============================================================================
## Inputs: -A list of model parameters, components, or a predefined set of these, to calibrate
##         -A list of observation time series (or a predefined set of them), against which to calibrate
##
## Outputs: -Details for all of the free parameters in this calibration
param_calib_details = function( model_params = NULL,     # Model parameters to calibrate
                                model_components = NULL, # Model components to calibrate
				model_set = "all_model", # Predefined set of model parameters to calibrate (defined above)
                                obs_ts = NULL,           # Observation time series against which to calibrate
				obs_set = "noOcheat_obs" )    # Predefined set of observations against which to calibrate (defined above)
{
    # Retrieve model parameter details
    if ( !is.null( model_params ) ) {
        nmodel = length(model_params)
        model_param_list = vector( "list", nmodel )
        for ( i in 1:nmodel ) {
            for ( j in 1:length(parlist.all_model) ) {
                if ( parlist.all_model[[j]]$param == model_params[i] ) {
                    model_param_list[[i]] <- parlist.all_model[[j]]
                }
            }
        }
    } else if ( is.null( model_params ) & !is.null( model_components ) ) {
        model_param_list = list()
        for ( i in 1:length(model_components) ) {
	    model_param_list = c( model_param_list, eval( parse( text = paste0( "parlist.", model_components[i] ) ) ) )
        }
    } else {
        model_param_list <- eval( parse( text = paste0( "parlist.", model_set ) ) )
    }

    # Retrieve observation-related parameter details
    if( !is.null( obs_ts ) ) { 
        obs_param_list = list()
	compare_list   = list()
        for ( i in 1:length( obs_ts ) ) {
            obs_param_list = c( obs_param_list, eval( parse( text = paste0 ("parlist.", obs_ts[i] ) ) ) )
	    compare_list   = append( compare_list  , eval( parse( test = paste0 ("compare.", obs_ts[i] ) ) ) ) 
        }
    } else {
        obs_param_list <- eval( parse( text = paste0( "parlist.", obs_set ) ) )
	compare_list   <- eval( parse( text = paste0( "compare.", obs_set ) ) )
    }
    
    # Combine model and observation parameters into a complete list of free parameters
    free_param_list = c( model_param_list, obs_param_list )
 
    params      = character(length(free_param_list))
    parnames    = character(length(free_param_list))
    sections    = character(length(free_param_list))
    in.hector   = logical(length(free_param_list))
    in.DEoptim  = logical(length(free_param_list))
    p0          = double(length(free_param_list))
    bound.lower = double(length(free_param_list))
    bound.upper = double(length(free_param_list))
    step.mcmc   = double(length(free_param_list))

    for ( i in 1:length(free_param_list) ) {
        params[i]      <- free_param_list[[i]]$param
        parnames[i]    <- free_param_list[[i]]$parname
        sections[i]    <- free_param_list[[i]]$component
        in.hector[i]   <- free_param_list[[i]]$in.hector
        in.DEoptim[i]  <- free_param_list[[i]]$in.DEoptim
        p0[i]          <- free_param_list[[i]]$p0
        bound.lower[i] <- free_param_list[[i]]$bound.lower
        bound.upper[i] <- free_param_list[[i]]$bound.upper
        step.mcmc[i]   <- free_param_list[[i]]$step.mcmc
    }
    index.DEoptim = which( in.DEoptim )
    
    compare.var       = character(length(compare_list))
    compare.component = character(length(compare_list))
    obs.ts            = character(length(compare_list))
    for ( i in 1:length(compare_list) ) {
        compare.var[i]       <- compare_list[[i]]$var
        compare.component[i] <- compare_list[[i]]$component
        obs.ts[i]            <- compare_list[[i]]$obsvar
    }
    return( list( params = params            , parnames = parnames                  , 
		  sections = sections        , in.hector = in.hector                , 
		  in.DEoptim = in.DEoptim    , p0 = p0                              ,
                  bound.lower = bound.lower  , bound.upper = bound.upper            ,  
		  step.mcmc = step.mcmc      , obs.ts = obs.ts                      ,
                  compare.var = compare.var  , compare.component = compare.component, 
                  index.DEoptim = index.DEoptim        ) )
}  

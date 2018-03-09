##==============================================================================
## obs_readData.R
##
## Original GSIC version by Kelsey Ruckert
## Modified for DOECLIM by Tony Wong
## Modified for Hector by Ben Vega-Westhoff
##
## Load observational time series and their associated uncertainty time series,
## for use in Hector calibration
##
## Currently includes:
##   -Surface temperature, HadCRUT4 global means (1850-2016), Morice et al (2012)
##   -Ocean heat, 0-3000 m (1953-1996), Gouretski and Koltermann (2007)
##   -SLR (global means, 1880-2013), Church and White (2011)
##   -Glaciers and small ice caps sea level contribution, 1961-2003, NSIDC (2005)
##   -Thermal expansion trends, 1971-2009 and 1993-2009, IPCC AR5 Ch. 13
##   -Greenland ice sheet sea level contribution, 1958-2009, Sasgen et al (2012)
##==============================================================================

##==============================================================================
## Temperature
##==============================================================================
# HADCRUT4 annual global mean surface temperature
# Note: all ensemble member files have same time series of uncertainties, so just
# grabbing the first one.
dat = read.table("obs_constraints/temperature/HadCRUT.4.4.0.0.annual_ns_avg.txt")
obs.temp = dat[,2]
obs.temp.time = dat[,1]
dat = read.table("obs_constraints/temperature/HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt")
obs.temp.err = dat[,3]

# Normalize temperature anomaly so 1850-1870 mean is 0
norm.lower = 1850; norm.upper = 1870
ibeg=which(obs.temp.time==norm.lower)
iend=which(obs.temp.time==norm.upper)
obs.temp = obs.temp - mean(obs.temp[ibeg:iend])

idx = compute_indices(obs.time=obs.temp.time, mod.time=mod.time)
oidx.temp = idx$oidx; midx.temp = idx$midx

hadcrut_sst_obs = list( midx = midx.temp, oidx = oidx.temp,
		 obs  = obs.temp,  obs.time = obs.temp.time,
		 obs.err = obs.temp.err, norm.lower = norm.lower,
		 norm.upper = norm.upper )

##==============================================================================
## Ocean heat content
##==============================================================================
# Gouretski annual global ocean heat content (0-3000 m)
dat = read.table("obs_constraints/ocheat/gouretski_ocean_heat_3000m.txt",skip=1)
obs.ocheat = dat[,2]
obs.ocheat.time = dat[,1]
obs.ocheat.err = dat[,3]

# Don't normalize these ocheat anomalies
norm.lower = NA; norm.upper = NA

idx = compute_indices(obs.time=obs.ocheat.time, mod.time=mod.time)
oidx.ocheat = idx$oidx; midx.ocheat = idx$midx

gouretski_ocheat_obs = list( midx = midx.ocheat, oidx = oidx.ocheat,
		   obs = obs.ocheat, obs.time = obs.ocheat.time,
                   obs.err = obs.ocheat.err, norm.lower = norm.lower,
		   norm.upper = norm.upper )  

##==============================================================================
## Global mean sea level
##==============================================================================
# Church and White
dat = read.table("obs_constraints/slr/GMSL_ChurchWhite2011_yr_2015.txt")
obs.slr      = dat[,2]/1000              # data are in mm
obs.slr.time = dat[,1]-0.5 # they give half-year so it is unambiguous as to which year they mean
obs.slr.err  = dat[,3]/1000

## Normalize SLR so 1986-2005 mean is 0
norm.lower = 1986; norm.upper = 2005
ibeg=which(obs.slr.time==norm.lower)
iend=which(obs.slr.time==norm.upper)
obs.slr = obs.slr - mean(obs.slr[ibeg:iend])

idx = compute_indices(obs.time=obs.slr.time, mod.time=mod.time)
oidx.slr = idx$oidx; midx.slr = idx$midx

church_slr_obs = list( midx = midx.slr, oidx = oidx.slr,
			obs = obs.slr, obs.time = obs.slr.time,
			obs.err = obs.slr.err, norm.lower = norm.lower,
			norm.upper = norm.upper )

##==============================================================================
## Glaciers and small ice caps
##==============================================================================
# NSIDC
# Historical global mean sea level contribution from Glacier and Small Ice Cap melt
dat = read.csv("obs_constraints/slr/GSICobservations_UPDATED.csv", skip = 1)
obs.gsic.time = dat[1:43, 1] #1961-2003
obs.gsic = dat[1:43, 5]/1000 # m of melt contribution to sea level rise (Note -- data are in mm)
obs.gsic.err = dat[1:43,9]/1000 # m (standard error) (Note -- data are in mm)

## Already normalized to 1960
norm.lower = 1960; norm.upper = 1960

idx = compute_indices(obs.time=obs.gsic.time, mod.time=mod.time)
oidx.gsic = idx$oidx; midx.gsic = idx$midx

nsidc_gsic_obs = list( midx = midx.gsic, oidx = oidx.gsic,
			obs = obs.gsic, obs.time = obs.gsic.time,
			obs.err = obs.gsic.err, norm.lower = norm.lower,
			norm.upper = norm.upper )

##==============================================================================
## Thermal expansion
##==============================================================================
# IPCC Thermal expansion
# Trends in thermal expansion from IPCC AR5 Ch. 13, m/yr

## trends.te = [Column 1=trend ; Column 2-3=90% window ; Column 4-5=beginning/ending year ;
##              Column 6-7=beginning/ending model indices for trend period]
## Note: IPCC trends are through 2010, but for the hindcast calibration, the
## forcing data go through 2009.
trends.te = mat.or.vec( 2 , 7)

trends.te[1,1:5] = c( 0.8/1000. , 0.5/1000. , 1.1/1000. , 1971 , 2009 )
trends.te[2,1:5] = c( 1.1/1000. , 0.8/1000. , 1.4/1000. , 1993 , 2009 )

## Compute the indices for the model to compare trends
for (i in 1:nrow(trends.te)) {
  idx = compute_indices(obs.time=(trends.te[i,4]:trends.te[i,5]) , mod.time=mod.time)
  trends.te[i,6:7] = c(idx$midx[1],(idx$midx[length(idx$midx)]))
}

# Don't normalize these (it's just trends)
norm.lower = NA; norm.upper = NA

ipcc_te_obs = list(trends = trends.te, norm.lower = norm.lower, norm.upper = norm.upper)

##==============================================================================
## Greenland ice sheet SLR contribution
##==============================================================================
# Sasgen et al (2012)
dat = read.csv("obs_constraints/slr/Greenland_OBS_MAR_InSAR_updated.csv")
obs.gis.time = dat[1:52,9]
obs.gis = dat[1:52,12] # Annual mass balance in meters sea level equivalence
obs.gis.err = rep( dat[1,15], length(obs.gis) )

idx = compute_indices(obs.time=obs.gis.time, mod.time=mod.time)
oidx.gis = idx$oidx; midx.gis = idx$midx

# These GIS values are already relative to 1961-1990
norm.lower = 1961; norm.upper = 1990

sasgen_gis_obs = list( midx = midx.gis, oidx = oidx.gis,
                       obs = obs.gis, obs.time = obs.gis.time,
                       obs.err = obs.gis.err, norm.lower = norm.lower,
                       norm.upper = norm.upper )

##==============================================================================

# !!!Add additional observational data sets to calibrate against

##==============================================================================
# !!!Add/switch in additional observational data sets to the obs.all list
obs.all = list( hadcrut_sst_obs, gouretski_ocheat_obs, nsidc_gsic_obs, ipcc_te_obs, sasgen_gis_obs )
names(obs.all) = c( "Tgav", "ocheat", "slr_gsic", "slr_te", "slr_gis" ) #Should match names in hector_calib_params.R
##==============================================================================
## End
##==============================================================================

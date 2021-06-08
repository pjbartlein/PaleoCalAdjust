## Download netCDF files ##

Example CMIP/PMIP-formatted netCDF source files (and resulting month-length adjusted files), along with the source data and files used to create the calendar-effects figures can be downloaded using the following URLs:

**Dropbox:  
[[https://www.dropbox.com/sh/hromb9qagzqk9kl/AAAcItqBRJzTnjW5VhjkSMOOa?dl=0]](https://www.dropbox.com/sh/hromb9qagzqk9kl/AAAcItqBRJzTnjW5VhjkSMOOa?dl=0)**

**Globus:  
[[https://app.globus.org/file-manager?origin_id=eb8ed828-b1c2-11eb-866b-d16fa0cfc9e7&origin_path=%2F]](https://app.globus.org/file-manager?origin_id=eb8ed828-b1c2-11eb-866b-d16fa0cfc9e7&origin_path=%2F)**

### Example source and adjusted CMIP5/PMIP3 and CMIP6/PMIP4 netCDF files ###

There are three sets of folders that contain model output for the *midHolocene* experiment that can be used to demonstrate `cal_adjust.f90`:  `/data/nc_files/test1/`, `/data/nc_files/test2/`, and `/data/nc_files/test3/` each of which contains a `source` and `adjusted` folder.  In addition, the folder `/data/nc_files/ctrl_nc_files/` contains *piControl* simulations that can be used to illustrate the impact of the calendar effect on long-term mean differences.

- `/data/nc_files/test1/` contains some typical CMIP5/PMIP3 netCDF files, including long-term means (`Aclim`) files, monthly time series (`Amon` files), and daily time series (`day`) files):
	
		tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		tas_Amon_IPSL-CM6A-LR_midHolocene_r1i1p1f1_gr_185001-204912.nc
		tas_Amon_IPSL-CM6A-LR_lig127k_r1i1p1f1_gr_185001-194912.nc
		tas_Amon_CNRM-CM5_midHolocene_r1i1p1_205001-214912.nc
		tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim.nc
		tas_Aclim_CNRM-CM5_midHolocene_r1i1p1_195001-214912-clim.nc
		tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		tas_day_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231.nc
		tas_day_CNRM-CM5_midHolocene_r1i1p1_19550101-19591231.nc
		tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc
		pr_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		pr_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc

- `/data/nc_files/test2/` contains contains a variety of CMIP5/PMIP3 and CMIP6/PMIP4 3-D and 4-D files, including some on rotated-pole ocean grids:

		tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim.nc
		tas_Aclim_CFSR_reanalysis_ana4mips_198101-201012-clim.nc
		tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		ta_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		ta_Amon_CCSM4_midHolocene_r1i1p1_120001-130012.nc
		tas_day_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231.nc
		tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc
		sos_Omon_CESM2_midHolocene_r1i1p1f1_gr_000101-005012.nc
		sos_Omon_CESM2_midHolocene_r1i1p1f1_gn_000101-005012.nc
		sic_OImon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		sic_OImon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		msftmyz_Oclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim.nc
		msftmyz_Oclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc

- `/data/nc_files/test3/` contains additional files that show the impact of the calendar effect on a range of other variables:

		clt_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		hfls_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		hfss_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		psl_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		rlnet_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		rnet_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		rsnet_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		sic_OIclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		ta_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		tos_Oclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		ua_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		va_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		wap_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		zg_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc

- `/data/nc_files/ctrl_nc_files/` contains *piControl* files:

		pr_Aclim_CCSM4_piControl_r1i1p1_025001-130012-clim.nc
		pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc
		pr_Amon_CCSM4_piControl_r1i1p1_100001-130012.nc
		pr_Amon_MPI-ESM-P_piControl_r1i1p1_185001-194912.nc
		ta_Amon_CCSM4_piControl_r1i1p1_120001-130012.nc
		tas_Aclim_CCSM4_piControl_r1i1p1_025001-130012-clim.nc
		tas_Aclim_CNRM-CM5_piControl_r1i1p1_185001-269912-clim.nc
		tas_Aclim_MPI-ESM-P_piControl_r1i1p1_185001-300512-clim.nc
		tas_Amon_CCSM4_piControl_r1i1p1_100001-130012.nc
		tas_Amon_CNRM-CM5_piControl_r1i1p1_260001-269912.nc
		tas_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-194912.nc
		tas_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-204912.nc
		tas_Amon_MPI-ESM-P_piControl_r1i1p1_185001-194912.nc

Note that the daily output files, e.g. `tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc`, are renamed by `cal_adjust.f90` to e.g. `tas_Amon2_CCSM4_midHolocene_r1i1p1_10000101-10491231_cal_adj.nc` after aggregating to monthly time-step data, using the file-name element `Amon2` to distinguish such files from standard monthly time-step files with the file-name element of `Amon`.

### netCDF files for calendar-effects figures ###

Data files used to create Figs. 10-12 and 16:

- `data/nc_files/cal_effects_nc_files/`

		CFSR_xm_gdd5_txn_cal_effects_006ka_04.nc
		CFSR_xm_gdd5_txn_cal_effects_097ka_04.nc
		CFSR_xm_gdd5_txn_cal_effects_116ka_04.nc
		CFSR_xm_gdd5_txn_cal_effects_127ka_04.nc
		CMAP_xm_cal_effects_006ka_04.nc
		CMAP_xm_cal_effects_097ka_04.nc
		CMAP_xm_cal_effects_116ka_04.nc
		CMAP_xm_cal_effects_127ka_04.nc

Source data for calendar-effects figures:

- `/data/nc_files/CFSR_and_CMAP_data/`
	
		enhanced.precip.mon.ltm.nc
		tas_Aclim_CFSR_reanalysis_ana4mips_198101-201012-clim.nc
		ana4mips_URL.txt
		cdo_ltm.txt
		ESRL-PSD_CMAP_URL.txt
		tas_Aclim_CFSR_reanalysis_ana4mips_198101-201012-clim_cal_adj.nc
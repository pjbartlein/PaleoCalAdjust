## Download netCDF files ##

Example CMIP5/PMIP3-formatted netCDF source files (and resulting month-length adjusted files), along with the source data and files used to create the calendar-effects figures can be downloaded using the following URLs:

**Dropbox:  
[[https://www.dropbox.com/sh/hromb9qagzqk9kl/AAAcItqBRJzTnjW5VhjkSMOOa?dl=0]](https://www.dropbox.com/sh/hromb9qagzqk9kl/AAAcItqBRJzTnjW5VhjkSMOOa?dl=0)**

**Globus:  
[[https://app.globus.org/file-manager?origin_id=3c131f8a-6df3-11ea-af52-0201714f6eab&origin_path=%2F]](https://app.globus.org/file-manager?origin_id=3c131f8a-6df3-11ea-af52-0201714f6eab&origin_path=%2F)**

### Example source and adjusted CMIP5/PMIP3 and CMIP6/PMIP4 netCDF files ###

There are four folders:  

- `/data/nc_files/PMIP3_source/` (which contains some typical CMIP5/PMIP3 netCDF files, including long-term means (`Aclim`) files, monthly time series (`Amon` files), and daily time series (`day`) files):
	
		tas_day_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231.nc
		tas_day_CNRM-CM5_midHolocene_r1i1p1_19550101-19591231.nc
		tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc
		tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		tas_Amon_CNRM-CM5_midHolocene_r1i1p1_205001-214912.nc
		tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim.nc
		tas_Aclim_CNRM-CM5_midHolocene_r1i1p1_195001-214912-clim.nc
		tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc
		pr_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc
		pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
		pr_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc

- `/data/nc_files/PMIP3_adjusted/` (which contains month-length-adjusted files created by `cal_adjust_PMIP.f90`):

		tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim_cal_adj.nc
		tas_Aclim_CNRM-CM5_midHolocene_r1i1p1_195001-214912-clim_cal_adj.nc
		tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim_cal_adj.nc
		tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012_cal_adj.nc
		tas_Amon_CNRM-CM5_midHolocene_r1i1p1_205001-214912_cal_adj.nc
		tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912_cal_adj.nc
		tas_Amon2_CCSM4_midHolocene_r1i1p1_10000101-10491231_cal_adj.nc
		tas_Amon2_CNRM-CM5_midHolocene_r1i1p1_19550101-19591231_cal_adj.nc
		tas_Amon2_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231_cal_adj.nc
		pr_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim_cal_adj.nc
		pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012_cal_adj.nc
		pr_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912_cal_adj.nc

- ` /data/nc_files/PMIP4_source/` (which contains two example CMIP6/PMIP4 files):

		tas_Amon_IPSL-CM6A-LR_lig127k_r1i1p1f1_gr_185001-194912.nc
		tas_Amon_IPSL-CM6A-LR_midHolocene_r1i1p1f1_gr_185001-204912.nc	

- ` /data/nc_files/PMIP4_adjusted/` (which contains month-length-adjusted CMIP6/PMIP4 files):

		tas_Amon_IPSL-CM6A-LR_lig127k_r1i1p1f1_gr_185001-194912_cal_adj.nc
		tas_Amon_IPSL-CM6A-LR_midHolocene_r1i1p1f1_gr_185001-204912_cal_adj.nc

Note that daily output files, e.g. `tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc`, are renamed by `cal_adjust_PMIP.f90` to e.g. `tas_Amon2_CCSM4_midHolocene_r1i1p1_10000101-10491231_cal_adj.nc` after aggregating to monthly time-step data, using the file-name element `Amon2` to distinguish such files from standard monthly time-step files with the file-name element of `Amon`.

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

	This folder also contains an adjusted netCDF files of CFSR *tas* data created by `cal_adjust_PMIP.f90` using the info file `cal_adj_info_reanalysis.csv`:

		tas_Aclim_CFSR_reanalysis_ana4mips_198101-201012-clim_cal_adj.nc
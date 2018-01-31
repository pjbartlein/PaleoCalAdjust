## netCDF files on Dropbox ##

Example CMIP5/PMIP3-formatted netCDF source files and resulting month-length adjusted files can be downloaded using the following links:

### Dropbox:  [[https://www.dropbox.com/sh/waerdamq9jgwtav/AAAQPgeJV-npNiTidUGyfvHJa?dl=0]](https://www.dropbox.com/sh/waerdamq9jgwtav/AAAQPgeJV-npNiTidUGyfvHJa?dl=0) ###


### Globus:  [[https://www.globus.org/app/transfer?origin_id=d74454ce-e6b4-11e7-80f5-0a208f818180&origin_path=%2FPaleoCalendarAdjust%2F]](https://www.globus.org/app/transfer?origin_id=d74454ce-e6b4-11e7-80f5-0a208f818180&origin_path=%2FPaleoCalendarAdjust%2F) ###


There are two folders:  
- `/source`, which contains some typical CMIP5/PMIP3 netCDF files, including long-term means (`Aclim`) files, monthly time series (`Amon` files), and daily time series (`day`) files:
	
	`/data/nc_files/source/`
		`tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc`
		`tas_Aclim_CNRM-CM5_midHolocene_r1i1p1_195001-214912-clim.nc`
		`tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim.nc`
		`tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc`
		`tas_Amon_CNRM-CM5_midHolocene_r1i1p1_205001-214912.nc`
		`tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc`
		`tas_day_CCSM4_midHolocene_r1i1p1_10000101-10491231.nc`
		`tas_day_CNRM-CM5_midHolocene_r1i1p1_19550101-19591231.nc`
		`tas_day_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231.nc`
		`pr_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim.nc`
		`pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc`
		`pr_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc`

- `/adjusted`, which contains monthlength-adjusted files created by `cal_adjust_PMIP.f90`:

	`/data/nc_files/adjusted/`
		`tas_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim_cal_adj.nc`
		`tas_Aclim_CNRM-CM5_midHolocene_r1i1p1_195001-214912-clim_cal_adj.nc`
		`tas_Aclim_MPI-ESM-P_midHolocene_r1i1p1_185001-194912-clim_cal_adj.nc`
		`tas_Amon_CCSM4_midHolocene_r1i1p1_100001-130012_cal_adj.nc`
		`tas_Amon_CNRM-CM5_midHolocene_r1i1p1_205001-214912_cal_adj.nc`
		`tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912_cal_adj.nc`
		`tas_Amon2_CCSM4_midHolocene_r1i1p1_10000101-10491231_cal_adj.nc`
		`tas_Amon2_CNRM-CM5_midHolocene_r1i1p1_19550101-19591231_cal_adj.nc`
		`tas_Amon2_MPI-ESM-P_midHolocene_r1i1p1_18500101-18591231_cal_adj.nc`
		`pr_Aclim_CCSM4_midHolocene_r1i1p1_100001-130012-clim_cal_adj.nc`
		`pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012_cal_adj.nc`
		`pr_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912_cal_adj.nc`



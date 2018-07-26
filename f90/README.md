## Programs ##

Main programs, including `month_lenth.f90` and `cal_adjust_PMIP3.f90`, plus additional driver and demonstration programs are in the folder `/main_programs`:

	month_length.f90			! month-length tables
	cal_adjust_PMIP3.f90			! paleo calendar adjustment
	GISS_orbpar_driver.f90			! orbital parameters (eccentricity, climatic precession)
	GISS_srevents_driver.f90		! equinox, solstice, perihelion dates
	demo_01_pseudo_daily_interp.f90		! demo of pseudo-daily interpolation
	demo_02_adjust_1yr.f90			! demo of paleo calendar adjustment
	demo_03_adjust_TraCE_ts.f90		! demo of transient simulation adjustment

The `/modules` folder contains the following:

	calendar_effects_subs.f90
	CMIP5_netCDF_subs.f90
	GISS_orbpar_subs.f90
	GISS_srevents_subs.f90
	month_length_subs.f90
	pseudo_daily_interp_subs.f90

The programs are used as follows:

- `GISS_orbar_driver.f90` and `GISS_srevents.f90` write orbital-parameter output to the folder `/GISS_orbital`, using specific parameter values set in the programs;
- `month_length.f90` reads the info file `month_length_info.csv` in the folder `/info_files` and writes month-length tables to the folder `/month_lengths`;
- `cal_adjust_PMIP3.f90` reads the info file `cal_adjust_info.csv` in the folder `/info_files`, and source netCDF files from the folder `/nc_files/source` and writes paleo calendar-adjusted netCDF files to the folder `/nc_files/adjusted`;
- `demo_01_pseudo_daily_interp.f90` and `demo_02_adjust_1yr.f90` are stand-alone programs, writing only to the console;
- `demo_03_adjust_TraCE_ts.f90` reads, for example, the file `TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv` in the folder `TraCE_example` and writes a paleo calendar-adjusted output file into the same folder.

All programs write progress information to the console.

With the exception of directory paths, the same code compiles and runs on the following systems:

- Windows 10: Intel Parallel Studio XE 2018 Update 3 Composer Edition for Fortran Windows, with netCDF version 4.1.3, using the Visual Studio 2015 IDE; 
- MacOS: gfortran version 8.1.0, with netCDF version 4.6.1 (from Homebrew), using the Eclipse Oxygen.3a Release (4.7.3a) IDE, and Eclipse for Parallel Application Developers, Version: 9.1.4.201802282054.
## Programs (PaleoCalAdjust v1.1) ##

Main programs, including `month_length.f90` and `cal_adjust_PMIP.f90`, plus additional driver and demonstration programs are in the folder `/main_programs`:

	month_length.f90				! month-length tables
	cal_adjust.f90					! paleo calendar adjustment
	GISS_orbpar_driver.f90			! orbital parameters (eccentricity, climatic precession)
	GISS_srevents_driver.f90		! equinox, solstice, perihelion dates
	demo_01_pseudo_daily_interp.f90	! demo of pseudo-daily interpolation
	demo_02_adjust_1yr.f90			! demo of paleo calendar adjustment
	demo_03_adjust_TraCE_ts.f90		! demo of transient simulation adjustment

The `/modules` folder contains the following:

	CMIP_netCDF_subs.f90			! subroutines for reading and writing netCDF files
	GISS_orbpar_subs.f90			! GISS orbital parameter calculations
	GISS_srevents_subs.f90			! GISS subroutines for finding days of equinoxes, solstices, aphelion and perihelion
	month_length_subs.f90			! subroutines for month-length calculations
	pseudo_daily_interp_subs.f90	! subroutines for managing pseudo-daily interpolation
	mp_interp_epstein_subs.f90		! Epstein (1991) harmonic-regression interpolation
	mp_interp_harzallah_subs.f90	! Harzallah (1995) iterative-spline interpolation
	spline_subs.f90					! Burkhardt spline-interpolation subroutines 

The module `spline_subs.f90` contains several subroutines and functions from John Burkhardt's library of Fortran90 spline-fitting subroutines, and is used by `mp_interp_harzallah_subs.f90` [[https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)]](https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html))

The `/projects` folder contains a set of subfolders, one for each main program, containing example GNU Make makefiles for the individual main programs.  The makefiles must be localized for a particular operating system or file structure (to correctly point to the compiler and source code).  The examples are for macOS.

The programs can be executed from the command line (as well as in Visual Studio or Eclipse), with the path and filename of an appropriate infofile appended to the program name.  For example, to build and run the program `month_length.f90` on macOS, after editing the makefile, open a terminal window in the `../PaleoCalAdjust/f90/makefiles/month_length/` folder and type, for example:

	make
	./month_length ../../../../PaleoCalAdjust/data/info_files/month_length_info_Mac.csv

The programs are used as follows:

- `GISS_orpar_driver.f90` and `GISS_srevents_driver.f90` write orbital-parameter output to the folder `/GISS_orbital`, using specific parameter values set in the programs;
- `month_length.f90` reads the info file `month_length_info.csv` in the folder `/info_files` and writes month-length tables to the folder `/month_lengths`;
- `cal_adjust_PMIP.f90` reads the info file `cal_adj_info.csv` in the folder `/info_files`, and source netCDF files from the folders `/nc_files/PMIP3_source` and `/nc_files/PMIP4_source` and writes paleo calendar-adjusted netCDF files to the folder `/nc_files/PMIP3_adjusted` and `/nc_files/PMIP4_adjusted`;
- `demo_01_pseudo_daily_interp.f90`, and `demo_02_adjust_1yr.f90`, are stand-alone programs, writing only to the console;
- `demo_03_adjust_TraCE_ts.f90` reads, for example, the file `TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv` in the folder `/TraCE_example` and writes a paleo calendar-adjusted output file into the same folder.

Most of the programs write progress information to the console.

The same code compiles and runs on the following systems:

- Windows 10: Intel Parallel Studio XE 2019 Update 4 Composer Edition for Fortran, with netCDF version 4.1.3, using the Visual Studio 2019 IDE; 
- MacOS: gfortran version 10.2.0 (Homebrew GCC 10.2.0_4), with netCDF version 4.8.0 (from Homebrew), using the Eclipse IDE for Scientific Computing, Version 4.19.0 (2021-03), with Eclipse for Parallel Application Developers, Version 9.4.0.
- Linux: (Ubuntu 18.04.5 LTS (bionic)):  gfortran version 7.5.0, with netcdf version 4.6.0.
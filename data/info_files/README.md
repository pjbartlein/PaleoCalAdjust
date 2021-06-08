## Info files ##

"Infofiles" are `.csv` files used to specify the specific characteristics of the simulations that are used to generate month-length tables that illustrate the effects of the changing orbit using `month_length.f90` or to specify the files that are to be adjusted using `cal_adjust.f90`.

Note:  In both verions of the info files, "simulation ages" control the calculation of the orbital parameters, while "simulation years" determine the occurrence of leap years and the date of occurrence of the vernal equinox.  Most paleoclimatic simulations use what might be regarded as a "fake" span of years for the simulations. For example, the *midHolocene* simulation by the `MPI-ESM-P` model in the file `tas_Amon_MPI-ESM-P_midHolocene_r1i1p1_185001-194912.nc` spans the 100-yr interval from January 1850 CE through December 1949 CE.  Because this model uses a proleptic Gregorian calendar, leap years are explicitly considered, and occur on an appropriate schedule.  Inspection of the `time_bnds` variable in the file that months 2 and 14 (February 1850 CE and 1851 CE respectively are each 28 days long, while month 26 (February 1852) is 29 days long as it should be because 1852 CE is a leap year.  See Bartein & Shafer (2019) for further discussion of the difference between simulation ages and simulation years.  Ignoring the timing of the leap years, or assuming that every fourth years is a leap year leads to artifacts that can be seen in time series of the the calendar effects. 

The info file for `month_length.f90` contains the following   

	! prefix        :: (string) a short name for the output month-length table
	! calendar_type :: (string) CF calendar type (e.g. "noleap", "proleptic_gregorian", etc.)
	! begageBP      :: (integer) beginning age of the table (negative, e.g. 10 ka = -10000 BP)
	! endageBP      :: (integer) ending age of the table (e.g. 0 = 0ka)
	! agestep       :: (integer) interval between age calculations
	! begyrCE       :: (integer) beginning calendar (actual or pseudo) year of simulation
	! nsimyrs       :: (integer) number of simulation ages 
	! outpath       :: (string)  path to and existing outpout folder, enclosed in single or double quote marks

And that for `cal_adjust.f90` contains the following

	! header        :: (string) header for individual files (ignored), then for each file:
	! activity      :: (string) PMIP3, PMIP4, etc. (used simply to label lines in the info file)
	! variable      :: (string) CMIP/PMIP variable name (e.g. "tas", "pr")
	! time_freq     :: (string) CMIP/PMIP output time-frequency type (e.g. "Amon", "Aclim", "day", "Oclim" ...)
	! calendar_type :: (string) calendar type (e.g. "noleap", "proleptic_gregorian", etc.)
	! begageBP      :: (integer) beginning simulation age (year BP) (e.g. 21000 (= 21 ka)
	! endageBP      :: (integer) ending simulation age (year BP) 
	! agestep       :: (integer) interval between age calculations
	! begyrCE       :: (integer) beginning calendar year of simulation for multi-year simulations at each age
	! nsimyrs       :: (integer) number of calendar years
	! interp_method :: (character) interpolation method ("Epstein" or "Harzallah")
	! no_negatives  :: (logical) true or false (constrrain interpolated values to greater or equal to 0.0)
	! match_mean    :: (logical) true or false (constrain the monthly mean of interpolated values to match input)
	! tol           :: (real(8)) tolerance value for enforce_mean() subroutine 
	! source_path   :: (string) path to (input) source file (enclosed in single or double quotation marks)
	! source_file   :: (string) source file name (enclosed in single or double quotation marks)
	! adjusted_path :: (string) path to (output) adjusted files (enclosed in single or double quotation marks)
	! adjusted_file :: (string) adjusted file name (enclosed in single or double quotation marks)

To facilitate creation of the files, example `.xlsx` files are provided that contain "Custom" formatting of the file and path columns that will enclose them in single quotations when saved as `.csv` files.

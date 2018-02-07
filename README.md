# PaleoCalendarAdjust #

#### Work-in-progress versions of programs for calculating orbitally determined paleo month lengths, and using these to adjust CMIP5/PMIP3-formatted paleo simulation files that were created by summarizing data using a fixed modern calendar. ####

**Overview**

The main program here **cal_adjust_PMIP3.f90** reads a CMIP5/PMIP3-formatted file, and replicates it by copying dimensions and attributes and replacing the main variable values with those that have been adjusted to reflect the appropriate calendar for the time of the simulation.

**Individual programs**

**cal_adjust_PMIP3.f90** uses modules:

- month_length.f90
- pseudo_daily_interp.f90
- etc.

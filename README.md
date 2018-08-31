PaleoCalendarAdjust
===================

The programs here implement an approach for adjusting for the "paleo calendar effect" -- the impact that the changes in the length of months or seasons over time, related to the changes in the eccentricity of Earth's orbit and to precession, have on summarization of paleoclimatic model output. The key program is `cal_adjust_PMIP3.f90` (in the folder `/f90`), which applies the adjustment to CMIP5/PMIP3-formatted files, and a related program, `month_length.f90`, that can be used to produce tables of the changing length of months over time. Figures illustrating the paleo calendar effect are in the folder `/figures`, and relevant data sets for exercising the programs are in the folder `/data`.

A draft of the full documentation of the paleo calendar effect in PMIP4 time-slice and transient experiments, including a discussion of the impact of the effect and the application of the programs, can be obtained by contacting Pat Bartlein (<bartlein@uoregon.edu>).

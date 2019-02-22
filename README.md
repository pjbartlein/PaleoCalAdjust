PaleoCalAdjust
===================

With the exception of several changes noted below, this is the repository that accompanies the paper:

Bartlein, P. J. and Shafer, S. L.: Paleo calendar-effect adjustments in time-slice and transient climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis, *Geosci. Model Dev. Discuss.*, 2018, 1-36, [https://doi.org/10.5194/gmd-2018-283](https://doi.org/10.5194/gmd-2018-283), 2018.

The programs implement an approach for adjusting for the "paleo calendar effect" -- the impact that the changes in the length of months or seasons over time (related to the changes in the eccentricity of Earth's orbit and to precession), have on the summarization of paleoclimatic model output. The current version of the key program is `cal_adjust_PMIP.f90` (in the folder `/f90`), which applies the adjustment to CMIP5/PMIP3- or CMIP6/PMIP4-formatted files, and a related program, `month_length.f90`, that can be used to produce tables of the changing length of months over time. (The original main program `cal_adjust_PMIP3.f90` was modified to accommodate CMIP6/PMIP4-formatted files.)  Figures illustrating the paleo calendar effect are in the folder `/figures`, and relevant data sets for exercising the programs are in the folder `/data`.  

Several minor modifications to the main program `cal_adjust_PMIP.f90` and its modules were made to accommodate CMIP6/PMIP4 files, the filenames of which contain an additional "grid_label" field not present in CMIP5/PMIP3 filenames.  Addditionally, following a referee's suggestion, we replaced the approach for calculating month lengths using the approximation of Kutzbach and Gallimore (1988, J. Geophys. Res. 93(D1):803-821), with a direct approach based on Kepler's equation. This substitution of approaches had no practical significance.  Several other code modifications were made in the interests of transparency.

The main changes from the original submission therefore include:

- the main program `cal_adjust_PMIP3.f90` and module `CMIP5_netCDF_subs.f90` were renamed as `cal_adjust_PMIP.f90` and `CMIP_netCDF_subs.f90` respectivly, because they are now generic.  No changes were made to `CMIP_netCDF_subs.f90`, and only minor changes were made to `cal_adjust_PMIP.f90` to allow reading and writing of both CMIP5/PMIP3 and CMIP6/PMIP4 formatted files;
- the info file read by `cal_adjust_PMIP.f90` (e.g. `cal_adj_info.csv`) was modified to include an "activity" field (either "PMIP3" or "PMIP4") and a "grid_label" field (blank for PMIP3 files);
- the paths to the "source" and "adjusted" netCDF files are now explicitly given in the info file, as opposed to being set in the main program;
- the subroutine `monlen(...)` in the module file `month_length_subs.f90` now computes month length, beginning, middle and ending days using an approach based on Kepler's equation.  To view the code changes and their impact on the results, the code and figures (plus the data used to plot the figures) in this release can be compared with those in previous releases.  
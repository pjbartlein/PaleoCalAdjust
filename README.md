PaleoCalAdjust
===================

This is the repository that accompanies the paper:

Bartlein, P. J. and Shafer, S. L.: Paleo calendar-effect adjustments in time-slice and transient climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis, *Geosci. Model Dev.*,  [https://doi.org/10.5194/gmd-12-3889-2019](https://doi.org/10.5194/gmd-12-3889-2019), 2019.

## Abstract ##

The “paleo calendar effect” is a common expression for the impact that the changes in the length of months or seasons over time, related to changes in the eccentricity of Earth’s orbit and precession, have on the analysis or summarization of climate-model output. This effect can have significant implications for paleoclimate analyses. In particular, using a “fixed-length” definition of months (i.e. defined by a fixed number of days), as opposed to a “fixed-angular” definition (i.e. defined by a fixed number of degrees of the Earth’s orbit), leads to comparisons of data from different positions along the Earth’s orbit when comparing paleo with modern simulations. This effect can impart characteristic spatial patterns or signals in comparisons of time-slice simulations that otherwise might be interpreted in terms of specific paleoclimatic mechanisms, and we provide examples for 6, 97, 116, and 127 ka. The calendar effect is exacerbated in transient climate simulations, where, in addition to spatial or map-pattern effects, it can influence the apparent timing of extrema in individual time series and the characterization of phase relationships among series. We outline an approach for adjusting paleo simulations that have been summarized using a modern fixed-length definition of months and that can also be used for summarizing and comparing data archived as daily data. We describe the implementation of this approach in a set of Fortran 90 programs and modules (PaleoCalAdjust v1.0 and v1.1).

## Programs ##

The current version of the key program is `cal_adjust.f90` (in the folder `/f90/main_programs`), which applies the paleo calendar-effect adjustment to CMIP/PMIP-type netCDF files, but should also work in the case of files that are "CF-compliant" (or nearly so), and they are also known to work with netCDF files of transient paleoclimatic simulations.  There is a related program, `month_length.f90`, that can be used to produce tables of the changing length of months over time that are used by `cal_adjust.f90`.   Figures illustrating the paleo calendar effect are in the folder `/figures`, and relevant data sets for exercising the programs are in the folder `/data`.  

Several minor modifications to the main program and its modules were made since the original *Geoscientific Model Development Discussions (GMDD)* manuscript submission to accommodate the adjustment of CMIP6-PMIP4 files.  Additionally, following a referee's suggestion, we replaced the approach used in the initial submission of the paper for calculating month lengths (i.e., the approximation of Kutzbach and Gallimore (1988, *J. Geophys. Res.* 93(D1):803-821)), with a direct approach based on Kepler's equation. This substitution of approaches had no practical significance.  Several other code modifications were made over time in the interests of transparency.

The current version is v1.1.  Relative to previous versions, this version includes:

- updated versions of `month_length.f90` and `cal_adjust.f90` (renamed from `cal_adjust_PMIP.f90` because the updated version is  no longer tied to the specific CMIP/PMIP-file naming structure);
- a new format for the "infofile" for each of these two programs that allows for the explicit specification of file paths and names, replacing CMIP6/PMIP4-style names assembled from infofile data;
- specification of the infofile path and name on the command line, so that once built locally, the programs can be run in a terminal window;
- addition of the adjusted month lengths, and beginning, middle, and ending dates to the output file;
- a choice of two mean-preserving interpolation methods, including the Epstein (1991) approach implemented in v1.0, as well as the Harzallah (1995) iterated-spline approach;
- the inclusion of a subroutine, `enforce_mean()` that requires the pseudo-daily interpolated values to have the same monthly mean as the input monthly values.

## Interpolation methods ##

The Epstein (1991) interpolation approach is intrinsically periodic, meaning that when applied to interpolate pseudo-daily values from monthly input values, the interpolated daily values at the end of the year will be consistent with those at the beginning, which is a desirable feature.  However, when iteratively applied to multi-annual time series of monthly data, small discontinuities will arise between years.  In the v1.0 implementation, this discontinuity was removed by smoothing the interpolated daily values at the end and beginning of the year.  The Harzallah (1995) approach, which involves iteratively fitting splines to the input data (and to the residuals from the original fit) is intrinsically not periodic, meaning the interpolated daily values at the end of the year will not be consistent with those at the beginning of a single year.  However, because this approach involves local as opposed to global fitting (as in the Epstein approach), the input data can be padded, either cyclically in the case of a single year's data, or with data from adjacent years in the case of time series.  This effectively eliminates the discontinuity between years.  The Epstein (1991) approach is recommended for adjusting "climatology" data sets (e.g. CMIP/PMIP "Aclim" time-frequency data), while the Harzallah (1995) approach is better suited for adjusting time-series data (e.g. "Amon"-type data sets), or transient-simulation data.

Despite the name, "mean-preserving" interpolation methods do not necessarily yield interpolated data that exactly reproduce the input data.  This can be addressed by setting a "tolerance" (`tol`) value for reproduction of the input values (typically 0.01 or 0.001 times the mean value of the data), that when exceeded, causes the discrepancy to be redistributed among the interpolated values.  Values of `tol` that are too large may lead to differences between the input monthly means and the means of the pseudo-daily interpolated values.  Values of `tol` that are too low may result in anomalously large adjusted values.  This will be easily seen in maps of the calendar effect, the difference between the adjusted and input monthly values.

Further discussion of mean-preserving interpolation, and comparisons among several practical approaches for its application can be found in the GitHub repository at [[https://github.com/pjbartlein/mp-interp]](https://github.com/pjbartlein/mp-interp).  Contact Pat Bartlein (bartlein@uoregon.edu) for further information.

Animations used in a presentation at the Fall 2019 AGU Meeting can be found in the `/animations` folder.  

References:  
Epstein, E), On obtaining daily climatological values from monthly means, *J. Climate* 4:365-368.  
Harzallah, A. (1995) The interpolation of data series using a constrained iterating technique *Monthly Weather Review* 123:2251-2254.  
Bartlein, P.J. and S.L. Shafer, 2019, Paleo calendar effects on radiation, atmospheric circulation, and surface temperature, moisture, and energy-balance variables can produce interpretable but spurious large-scale patterns and trends in analyses of paleoclimatic simulations. PP31A-08, AGU 2019 Fall Meeting.  [[https://agu.confex.com/agu/fm19/meetingapp.cgi/Paper/525140]](https://agu.confex.com/agu/fm19/meetingapp.cgi/Paper/525140)

## Version history ##

### v1.0 ###

Original release.

### v1.0a ###

Minor modifications for consistency with the Bartlein and Shafer (2018, *Geoscientific Model Development Discussions*) submission.

### v1.0b ###

Several minor modifications to the main program `cal_adjust_PMIP.f90` and its modules were made to accommodate CMIP6/PMIP4 files, the filenames of which contain an additional "grid_label" field not present in CMIP5/PMIP3 filenames.  

- The main program `cal_adjust_PMIP3.f90` and the module `CMIP5_netCDF_subs.f90` were renamed as `cal_adjust_PMIP.f90` and `CMIP_netCDF_subs.f90` respectively, because they now accommodate CMIP5/PMIP3- and CMIP6/PMIP4-formatted files.  No changes were made to `CMIP_netCDF_subs.f90`, and only minor changes were made to `cal_adjust_PMIP.f90` to allow reading and writing of both CMIP5/PMIP3- and CMIP6/PMIP4-formatted files;

- The info file read by `cal_adjust_PMIP.f90` (e.g. `cal_adj_info.csv`) was modified to include an "activity" field (either "PMIP3" or "PMIP4") and a "grid_label" field (blank for PMIP3 files), and the paths to the "source" and "adjusted" netCDF files are now explicitly given in the info file, as opposed to being set in the main program.

### v1.0c ###

Following a referee's suggestion, we replaced the approach for calculating month lengths using the approximation of Kutzbach and Gallimore (1988, *J. Geophys. Res.* 93(D1):803-821), with a direct approach based on Kepler's equation.  This substitution of approaches had no practical significance.  Several other code modifications were made in the interests of transparency.

- The subroutine `monlen(...)` in the module file `month_length_subs.f90` now computes month length, beginning, middle and ending days using an approach based on Kepler's equation.  To view the code changes and their negligible impact on the results, the code and figures (plus the data used to plot the figures) in this release can be compared with those in previous releases (e.g. v1.0b).    

### v1.0d ###

The main changes from the original submission include:

- The main program `cal_adjust_PMIP.f90`, was modified to handle rotated-pole and 4-dimensional (e.g. time x level x latitude x longitude) files.  This change required development of an additional subroutine, `get_var_diminfo(...)` in the module `CMIP_netCDF_subs.f90`. This release (v1.0d) is consistent with the revised and accepted version of the paper in the GMDD discussion, but has been modified to also handle 4-D files (e.g. *lon* x *lat* x *level* x *time*), and rotated-pole files (e.g. *Oclim* or *OIclim* files).
- the main program `cal_adjust_PMIP3.f90` and module `CMIP5_netCDF_subs.f90` were renamed as `cal_adjust_PMIP.f90` and `CMIP_netCDF_subs.f90` respectivly, because they are now generic.  An additional subroutine was added to the module `CMIP_netCDF_subs.f90` (see below), and only minor changes were made to `cal_adjust_PMIP.f90` to allow reading and writing of both CMIP5-PMIP3 and CMIP6-PMIP4 formatted files;
- the info file read by `cal_adjust_PMIP.f90` (e.g. `cal_adj_info.csv`) was modified to include an "activity" field (either "PMIP3" or "PMIP4") and a "grid_label" field (blank for PMIP3 files);
- the subroutine `monlen(...)` in the module file `month_length_subs.f90` now computes month length, beginning, middle and ending days using an approach based on Kepler's equation.  To view the code changes and their impact on the results, the code and figures (plus the data used to plot the figures) in this release can be compared with those in previous releases (e.g. v1.0b);
- The main program `cal_adjust_PMIP.f90` was modified to allow the adjustment of 4-D data sets and generalized to work with "rotated-pole" data sets (by eliminating explicit references to longitude and latitude).  This change required the addition of a subroutine `get_var_diminfo(...)` to `CMIP_netCDF_subs.f90` for getting dimension-variable information to `CMIP_netCDF_subs.f90`     



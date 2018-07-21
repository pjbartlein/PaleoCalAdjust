PaleoCalendarAdjust
===================

**Work-in-progress versions of programs for calculating orbitally determined paleo month lengths, and using these to adjust CMIP5/PMIP3-formatted paleo simulation files that were created by summarizing data using a fixed modern calendar.**

The Paleo Calendar Effect
-------------------------

The material below describes the "paleo calendar effect" \-- the common expression for the impact that the changes in the length of months or seasons over time, related to the changes in the eccentricity of Earth's orbit and to precession, has on summarization of paleoclimatic model output. This is followed by a description of an approach implemented in Fortran programs (and modules), that can be used to determine the changing length of months (on various calendars) over time, and to adjust existing data sets in the "CMIP5" format to an appropriate paleo calendar.

### Introduction

The calendar effect arises from a consequence of Kepler's second law of planetary motion: Earth moves faster along its elliptical orbit near perihelion, and slower near aphelion. Because the time of year of perihelion and aphelion vary over time, the length of time that it takes the Earth to traverse one-quarter or one-twelfth of its orbit (a nominal season or month) also varies, so that months are shorter near perihelion and longer near aphelion. When examining present day and paleo simulations, summarizing data using a fixed-length definition of a particular month (e.g. January as the first 31 days of a 365-day year), as opposed to a "fixed-angular" definition (e.g. (31/365.2425)•360 degrees of orbit), will therefore result in comparing conditions that prevailed over different portions of the orbit.

The impact of the calendar effect on the analysis of paleoclimatic simulations and their comparison with present-day or "control" simulations is well known and not trivial (e.g. Kutzbach and Gallimore, 1988; Joussaume and Braconnot, 1997), and various approaches have been proposed for incorporating the effect or "adjusting" monthly values in analyses of present and paleo simulations (e.g. Pollard and Reusch, 2002; Timm et al., 2008; Chen et al., 2011). Despite this work, the calendar effect is generally ignored or dismissed, and so our motivation here is to provide an adjustment method that is relatively simple and can be applied generally to "CMIP5-formatted" files like those that PMIP simulations are distributed as. With a suitable "wrapper," the approach can also be applied to transient simulations (e.g. Liu et al., 2009; Ivanovic et al., 2016).

The impact of the calendar effect is large and spatially variable, and can produce apparent map patterns that might otherwise be interpreted as evidence of, for example, latitudinal amplification or damping of temperature changes, development of continental/marine temperature contrasts, interhemispheric contrasts (the "bipolar see-saw"), changes in the latitude of the intertropical convergence zone (ITCZ), variations in strength of the global monsoon, and others. In transient climate-model simulations, time series of data aggregated using a fixed modern calendar as opposed to an appropriately changing one can differ not only in the overall shape of long-term trends in the series, but also in variations in the timing of, for example, Holocene "thermal maxima" which, depending on the time of year, can be on the order of several thousand years. The impact arises not only from the orbitally controlled changes in insolation amount and the length of months or seasons, but also from the advancement or delay in the starting and ending days of months or seasons relative to the solstices. Even if daily data are available, the calendar effect must still be considered when summarizing those data by months or seasons, or when calculating climatic indices such as the temperature of the warmest or coldest month---values that are often used for comparisons with paleoclimatic observations.

A number of approaches for adjusting monthly data that were averaged using present-day calendar definitions to a "paleo calendar" are described in the above references. A simple one (broadly similar to Pollard and Reusch, 2002) that we implement here involves 1) determining the appropriate fixed-angular month lengths for a paleo experiment (e.g., Kutzbach and Gallimore, 1988), 2) interpolating the data to a daily time step using a mean-preserving interpolation method (e.g., Epstein, 1991), and then 3) averaging or accumulating the interpolated daily data using the appropriate (paleo) month starting and ending days, thereby explicitly incorporating the changing month lengths. In cases where daily data are available (e.g. in CMIP5/PMIP3 "day" files), only the third step is necessary.

In the following discussion we describe a) the changing lengths of months and their beginning, middle and ending days over the past 150 kyr; b) the spatial patterns of the calendar effects on temperature and precipitation rate for several key times, and c) the methods and programs that can be used to calculate month lengths and to adjust monthly and daily paleo model output to an appropriate paleo calendar.

### Month-length variations

The length of months as they vary over time can be calculated using the algorithm in Appendix A of Kutzbach and Gallimore (1988). This algorithm yields the length of time (in days, or fractions thereof) required to traverse given number of degrees of celestial (as opposed to geographical) longitude starting from the vernal equinox, the common "origin" for orbital calculations (see Jousaumme and Braconnot, 1997, for discussion). Although developed for a 360-day year, the algorithm can easily be adapted for other calendars and for cases where the present-day months are not equal in length. The beginnings and ends of each month in a 365-day "noleap" calendar (<http://cfconventions.org>) are shown in Fig. 1, with the month-length "anomalies" (or long-term differences between paleo and present) shown in color, with (paleo) months that are shorter than those at present (i.e. 1950 CE) in green shades, and months that are longer in blue. Not only do the lengths of months vary, but so do their middle, beginning and ending days (Fig. 2), with mid-month days that are closer to the June solstice indicated in red and those that are farther in blue. The patterns of variations in month length (Fig. 1) obviously track the changing time of year of perihelion, while the beginning and ending days reflect the climatic precession parameter (Fig. 2).

The calendar effect is illustrated below for four times: 6 and 127 ka are the target times for the planned warm-interval *midHolocene* and *lig127k* CMIP6/PMIP4 simulations (Otto-Bliesner et al., 2017) and illustrate the calendar effects when perihelion occurs in the boreal summer or autumn (Fig. S1); 116 ka is the time of a proposed sensitivity experiment for the onset of glaciation, and illustrates the effect when perihelion occurs in boreal winter; and 97 ka was chosen to illustrate an orbital configuration not represented by the other times (i.e. one with boreal spring months occurring closer to the June solstice).

At 6ka, perihelion occurred in September (Fig. S1), and the months from May through October were shorter than today (Fig. 1), with the greatest differences in August (1.65 days shorter than present). This contraction of month lengths moved the middle all of the months from April through December closer to the June solstice (Fig. 2), with the greatest difference in November (4.5 days closer to the June solstice, and so 4.5 days farther from the December solstice). At 127 ka, perihelion was in late June, and the months May through September were shorter than today, with the greatest difference in July (3.23 days shorter than present). As at 6 ka, the shorter boreal summer months move the middle of the months between July and December closer to the June solstice, with the greatest difference in September and October (12.4 and 12.3 days closer, respectively). At both 6 and 127 ka, the longer boreal winter months begin and end earlier in the year, placing the middle of January 3.3 (6 ka) and 4.3 (127 ka) farther away from the June solstice than at present. As can be noted on Figs. 1 and 2, 127 ka does not represent a simple amplification of 6 ka conditions. Although broadly similar in having shorter late boreal summer and autumn months, that begin earlier in the year (and hence closer to the June solstice), the two times are only broadly similar in the relative differences from present in month length and beginning and ending days.

At 116 ka perihelion was in late December, and consequently the months from October through March were shorter than present. This has the main effect of moving the middle of the months July through December farther from the June solstice (with a maximum in September of 5.6 days), somewhat opposite to the pattern at 6 and 127 ka. At 97 ka, perihelion occurred in mid November, in between its occurrence in September at 6ka and December at 116 ka. The impact on month length and mid-month timing is complicated, with both the mid-month days of January and July occurring closer to the June solstice.

The first-order impact of the calendar effect can be gauged by comparing (at a particular latitude) daily insolation values for mid-month days determined using the appropriate paleo calendar (which assumes fixed-angular definitions of months) with insolation values for mid-month days using the present-day calendar (which assumes fixed-length definitions of months). At 6 ka, at 45ºN, the shorter (than present), but earlier (relative to the June solstice) months of September through November had insolation values over 10 Wm^-2^ greater for mid-month days defined using the paleo calendar, in comparison with values determined using the present-day calendar (Fig. 3), and at 127 ka, the differences exceeded 40 Wm-2 for the months of August through October. These positive differences were accompanied by negative differences during the winter months. At 116 ka, the longer but later occurring months of September and October had negative differences in mid-month insolation that exceeded 10 Wm^-2^. For regions where surface temperatures are strongly tied to insolation with little lag, such as the interiors of the northern continents, these calendar effects on insolation will directly be reflected by the calendar effects on temperatures. By moving the beginning, middle and end of individual months (and seasons) closer to or farther from the solstices, the "apparent temperature" of those intervals will be affected (i.e. months or seasons that start or end closer to the summer solstice will be warmer). The calendar effect on insolation varies strongly with latitude, with the sign of the difference broadly reversing in the southern hemisphere (Figs. S2 and S3).

### Impact of the calendar effect

The impact on paleoclimatic simulations of the calendar effect can be assessed by assuming that the long-term mean difference in climate (also referred to as the experiment minus control "anomaly") is zero everywhere. Pseudo-daily interpolated values (or actual daily output, if available) of present-day monthly data can then simply reaggregated using the appropriate paleo calendar and compared with the present-day data. (The pseudo-daily values were obtained here using a mean-preserving algorithm described below.)

The "pure" calendar effect is demonstrated here using monthly long-term mean (1981-2010) values of near-surface air temperature (*tas*) from the Climate Forecast System Reanalysis (CFSR; <https://esgf.nccs.nasa.gov/projects/ana4mips/>), and monthly precipitation rate (*pr*) from the CPC Merged Analysis of Precipitation (CMAP; ([https://www.esrl.noaa.gov/psd/data/
gridded/data.cmap.html](https://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html)) (Fig. S4). (These data were chosen because they are global in extent, and of reasonably high resolution.) If it is assumed that there is no long-term mean difference between a present-day and paleo simulation (by adopting the present-day data as the simulated paleo data), then the unadjusted data can be compared with data adjusted to the appropriate paleo month lengths. The adjusted minus unadjusted differences will therefore reveal the inverse of the built-in "signal" in the unadjusted data that might easily be interpreted in terms of some specific paleoclimatic mechanisms while being instead a data analytical artifact. Positive values on the maps indicate, for example, where temperatures would be higher or precipitation greater if a fixed-angular calendar were used to summarize the paleo data.

#### Monthly temperature

The impacts of using the appropriate calendar to summarize the data (as opposed to not) are large, often exceeding 1ºC in absolute value (Fig. 4). The effects are spatially variable, and not simple functions of latitude as might be initially expected, because the effect increases with the amplitude of the annual cycle (which has a substantial longitudinal component) for temperature regimes that are in phase with the annual cycle of insolation. For temperature regimes that are out of phase with insolation, the effect would be negative, and largest when the temperature variations were exactly out of phase. (If there were no annual cycle, i.e. if

a climate variable remained constant over the course of a year, the calendar effect would be zero.) The interaction between the annual cycle and the direct calendar effect on insolation produces patterns of the overall calendar effect that happen to resemble some of the large-scale responses that are frequently found in climate simulations, both past and future, such as high-latitude amplification or damping, continental-ocean contrasts, interhemispheric contrasts and changes in seasonality (c.f. Izumi et al., 2013). Because the month-length calculations use the N.H. vernal equinox as a fixed origin for the location of Earth along its orbit, the effects seem to be small during the month surrounding the equinox (i.e. February through April), and indeed the selection of a different origin would produce different apparent effects (see Jousaumme and Braconnot, 1997, sec. 2.1). However, the selection of a different origin would not change the relative length of time it would take Earth to transit any particular angular segment of its orbit.

At 6 ka, the largest calendar effects on temperature can be observed over the Northern Hemisphere continents for the months from September through December (Fig. 4), consistent with the earlier beginning of these months (Fig. 2) and the direct calendar effect on insolation at 45ºN (Fig. 3). For example, in the interior of the northern continents, as well as North Africa, temperature is in phase with insolation, and so the calendar effect on insolation (Fig. 3) which is strongly positive from August through November is reflected by the calendar effect on temperature. Over the northern oceans, temperature is broadly in phase with insolation, but with a lag, which reduces the magnitude of the effect, and gives rise to an apparent land-ocean contrast which otherwise might be interpreted as some particular paleoclimatic mechanism. The calendar effect on temperature from January through March is negative in the northern continental interiors (Fig. 4), which is also consistent with the calendar effect on insolation. In the Southern Hemisphere at 6 ka, the calendar effects on temperature are negative, which is consistent with the calendar effects on mid-month insolation at 45ºS (Fig. S3), which are generally negative throughout the year, particularly during the months of August through November. Like the continent -- ocean contrast in the Northern Hemisphere, the Northern Hemisphere -- Southern Hemisphere contrast in the calendar effect on temperature also could be interpreted in terms of one or another of the mechanisms thought to be responsible for interhemispheric temperature contrasts.

At 127 ka, the calendar effect on temperature is broadly similar to that at 6 ka over the months from September through March, but differs in sign from May through July, and in magnitude in August (Fig. 4). This is also consistent with the direct calendar effects on insolation. At 127 ka, the calendar effect on insolation becomes positive in the Northern Hemisphere earlier in the northern summer than at 6 ka (Fig. 3), while at 45ºS the calendar effects on insolation become strongly negative in July and persist that way through November (Fig. S3). At 116 ka, perihelion occurs in late December, in comparison late June at 127 ka (Figs. 1 and S1), and not surprisingly the calendar effect on temperature is nearly the inverse of that at 127 ka. This is alarming, because in addition to all of the changes in forcing and paleoclimatic responses accompanying the transition out of the last interglacial, the possibility that some of the apparent change in simulations between 127 and 116 ka may be an artifact of data analysis procedures cannot be discounted.

At 97 ka, a time selected to illustrate a different orbital configuration than the similar (6 ka and 127 ka) or contrasting (127 and 116 ka) configurations, the calendar effect on temperature in the Northern Hemisphere (Fig. 4) shows a switch from positive in the early boreal summer (May and June) to negative in the late summer (August and September). This is again consistent with the direct calendar effect on insolation (Fig. 3). Like the other times, the spatial variations in the calendar effect could easily be interpreted in terms of one kind of paleoclimatic mechanism or another.

#### Mean temperature of the warmest and coldest months

Although the calendar effects on monthly temperature show some sub-continental scale variability, the overall patterns are of relatively large spatial scales, and are interpretable in terms of the direct orbital effects on month lengths and insolation. The calendar effects on the mean temperature of the warmest (MTWA) and coldest (MTCO) months (and their differences) are much more spatially variable (Fig. 5). This variability arises in large part because of the way these variables are usually defined (e.g. as the temperature of the warmest or coldest conventionally defined month, as opposed to the temperature of the warmest or coldest 30-day interval), but also because the calendar adjustment can result in a change in the specific month that is warmest or coldest. These effects are compounded when calculating seasonality (as MTWA minus MTCO). In the particular set of example times chosen here, the magnitudes of the calendar effects are also smaller than those of individual months because, as it happens, the calendar effects in January, February, July and August are not large. There are also some surprising patterns. The inverse relationship between the calendar effects at 116 ka and 127 ka that might be expected from inspection of the monthly effects (Fig. 4) are not present, while the calendar effects on MTCO and MTWA at 97 ka and 116 ka tend to resemble one another. Across the four example times, there is indistinct, but still noticeable pattern in reduced seasonality (MTWA minus MTCO) between the adjusted and unadjusted values, which like the other patterns described above could tempt interpretation in terms of specific mechanisms.

#### Monthly precipitation

In contrast to the large spatial-scale patterns of the calendar effect on temperature, the patterns of the calendar effect on precipitation rate are much more complex, showing both continental-scale patterns (like those for temperature), but also features that are apparently related to precipitation associated with the ITCZ and regional and global monsoons (Fig. 6). The first pattern is evident in the calendar effects at 6 and 127 ka, particularly in the months from September through November (Fig. 6), where it also can be noted (especially over the mid-latitude continents in both hemispheres) that there is a positive association with the calendar effect on temperature. This association is related simply to similarities in the shapes of the annual cycles of those variables, and not to some kind of thermodynamic process. At 116 ka, as for temperature, the large-scale calendar-effect patterns appear to be nearly the inverse of those at 127 ka. The second pattern is well illustrated at 127 ka in the region from the tropical North Atlantic, sub-Saharan Africa and south Asia. There, a negative calendar effect can be noted for June through August, giving way to a positive effect from September through November, and the same transition appears inversely at 116 ka. Another example can be found in the South Pacific Convergence Zone in spring and early summer (September through November) at 6 and 127 ka. This second kind of pattern, most evident in the subtropics, is not mirrored by the calendar effects on temperature.

#### Summary

Overall, the magnitude and spatial patterns of the calendar effects on temperature and precipitation (Figs. 4-6) resemble those in the paleoclimatic simulations and observations that we attempt to explain in mechanistic terms. Depending on the sign of the effect, neglecting to account for the calendar effects could spuriously amplify some "signals" in long-term mean differences between experiment and control simulations, while damping others.

### Calendar effects and transient experiments

Calendar effects must also be considered in the analysis of transient climate-model simulations (even if those data are available on the daily time step). This can be illustrated for a variety of variables and regions using data from the TraCE-21k transient simulations (Liu, et al., 2009; <https://www.earthsystemgrid.org/project/trace.html>). The series plotted in Fig. 7 are area-averages for individual months on a yearly time step, with a 100-yr (half-window width) locally weighted regression curves added to emphasize century-timescale variations. The original yearly time-step data were aggregated using a perpetual "no leap" (365-day) calendar (using the present-day month lengths for all years). The gray and black curves on Fig. 7 show these unadjusted "original" values, while the colored curves show month-length adjusted values (i.e. daily interpolated values, reaggregated using the appropriate paleo calendar). Area averages were calculated for ice-free land points.

Figure 7a shows area-weighted averages for 2-m air temperature for a region that spans 15 to 75ºN and ‑170 to 60 ºE, the region used by Marsicek et al. (2018) to discuss Holocene temperature trends in simulations and reconstructions. The largest differences between month-length adjusted values and unadjusted values (for the four months shown) occur in October between 14 and 6 ka, when perihelion occurred during the northern summer months. October month lengths during this interval were within one day of that at present (Fig. 1), but the generally shorter months of April through September resulted in Octobers beginning up to ten days earlier in the calendar than at present, i.e. closer in time to the boreal summer solstice (Fig 2). The calendar-effect adjusted values therefore average up to 4ºC higher than the unadjusted values during this interval, consistent with the direct calendar effects on insolation at 45ºN (Fig. 3). The calendar effect also changes the shape of the temporal trends in the data, particularly during the Holocene. October temperatures in the unadjusted data showed generally increasing trend over the Holocene (i.e. since 11.7 ka), reaching a maximum around 3 ka, comparable with present-day values, while the adjusted data reached levels consistently above present-day values by 7.5 ka. The unadjusted October temperature data could be described as reaching a "Holocene thermal maximum" only in the late Holocene (after 4 ka), while the adjusted data displays more of a mid-Holocene maximum. As is the case with the mapped assessments of the "pure" calendar effect, the differences between unadjusted and adjusted time series are of the kind that could be interpreted in terms of various hypothesized mechanisms. For example, the calendar-effect adjustment delays the time of occurrence of a Holocene thermal maximum in October by about 3 kyr.

As in North America and Europe, the adjusted temperature trends in Australia (10 to 50ºS and 110 to 160ºE) (Fig. 7b) are consistent with the direct calendar effects on insolation (i.e. for 45ºS, Fig. S3). The difference between adjusted and unadjusted values are again largest in October, but the difference is the inverse of that for the North America and Europe region, because annual cycle of temperature is inversely related to the annual cycle of the insolation anomalies themselves (Fig. S7) and so to the direct calendar effects on insolation (Fig. S3). and again, the shapes of the Holocene trends in the adjusted and unadjusted data are noticeably different. In these two examples, relatively large areas are being averaged, and the calendar effect becomes more apparent as the size of the area decreases. Notably, the effect does not completely disappear at the largest scales, i.e. for area-weighted averages for the globe (for ice-free land grid cells) (Fig. 7c). The differences are small, but still discernable.

In the Northern Hemisphere (African-Asian) Monsoon region (0 to 30ºN and ‑30 to 120ºE), the calendar effects on precipitation rate are similar to those on temperature in the midlatitudes because the annual cycle of precipitation is roughly in phase with that of insolation (Fig. S5). There is little effect in the winter and spring, but a substantial effect in summer and autumn over the interval from 17 ka to about 3 ka (Fig. 7d). The calendar effect reverse sign between July and August (when the month-length adjusted precipitation rate values are less than the unadjusted ones) and September and October (when the adjusted values are greater than the unadjusted ones). In July, the timing of relative maxima and minima in the two data sets is similar, while in October, in particular, the Holocene precipitation maximum is several thousand years earlier in the adjusted data than in the unadjusted.

The time-series expression of the latitudinally reversing calendar effect on precipitation rate evident in Fig. 6 can be illustrated by comparing precipitation or precipitation minus evaporation ($P\ –\ E$) for the North African (sub-Saharan) Monsoon region (5 to 17ºN and ‑5 to 30ºE) with the Mediterranean region (31 to 43ºN and -5 to 30ºE) (Figs. 7e and 7f). The differences between the adjusted and unadjusted data in the North African region (Fig. 7e) parallels that of the larger monsoon region (Fig. 7d). The Mediterranean region, which is characteristically moister in winter and drier in summer shows the reverse pattern: when the adjusted minus unadjusted difference is positive in the monsoon region, it is negative in the Mediterranean. Dipoles are frequently observed in climatic data, both present-day and paleo, and usually interpreted in terms of broad-scale circulation changes in the atmosphere or ocean. This example illustrates that they could also be artifacts of the calendar effect. The impacts of the calendar effect on temporal trends in transient simulations (Fig. 7), when compounded by the spatial effects (Figs. 4-6), makes it even more likely to infer spurious climatic mechanisms in analyzing transient simulations that in the simpler time-slice simulations.

Methods
-------

The approach we describe here for adjusting model output reported either as monthly data (using fixed-length definitions of months) or as daily data to reflect the calendar effect (i.e. to make month-length corrections) has two fundamental steps: 1) pseudo-daily interpolation of the data on a fixed-month-length calendar (which when actual daily data are available is not necessary), followed by 2) aggregation of those data to fixed-angular months defined for the particular time of the simulations. The second step obviously requires the calculation of the beginning and ending days of each month as they vary over ("geological") time, which in turn depends on the orbital parameters. The definition of the beginning and ending days of a month in a "leap-year," "Gregorian," or "proleptic Gregorian" calendar additionally depends on timing of the (northern) vernal equinox, which varies from year to year. To calculate the orbital parameters using the Berger (1978) solution, and the timing of the (northern) vernal equinox (as well as insolation itself) we adapted a set of programs provided by National Aeronautics and Space Administration, Goddard Institute for Space Studies (<https://data.giss.nasa.gov/ar5/SOLAR/>). For clarity, the month-length calculations will be described first, followed by the a discussion of pseudo-daily interpolation methods. Then the calendar-program will be described, along with a few demonstration programs that exercise some of the individual procedures.

### Month-length calculation


What does code look like?

	! Fortran code
	stop "dave"

In-line code `month_lengths_subs.f90` is a subroutine.

### Pseudo-daily interpolation

It turns out that the most common way of producing pseudo-daily values, linear interpolation be-tween monthly means, is not mean preserving; the monthly means of the interpolated daily values will generally not match the original values. Here we use the mean-preserving "harmonic" interpolation method of Epstein (1991, J. Climate 4:365-368) to gauge this mismatch, and compare it with the calendar-effect adjustment (which itself is sensitive to the daily interpolation method).

Two approaches are compared here using Climate Forecast System Reanalysis (CFSR) near-surface air temperature and CPC Merged Analysis of Precipitation (CMAP) 1981-2101 long-term mean values (left-hand column of monthly maps far right above, and black curves below):

1\) linear interpolation (blue, interpolation errors in the middle column of maps), and

2\) mean-preserving interpolation (Epstein, 1991) (red, interpolation errors in the right-hand column of maps; note the different scales for interpolation errors).

Like other harmonic-based approaches, the Epstein approach can create interpolated curves that are wavy (see Pollard and Reusch (2002) for discussion), but these effects are small enough to not be practically important). The pathological case for precipitation is shown in the lower-right panel below.

The linear interpolation errors are quite large, with absolute values exceeding 1ºC and 1 mm d-1 and have distinct seasonal and spatial patterns: underpredictions of temperature in summers (and overpredictions in winter, and underpredictions of precipitation in the wet season and locations (and overpredictions in the dry season and locations). The magnitude and patterns of these effects rival those we attempt to infer or interpret in the paleo record. The mean-preserving interpolation errors for temperature are very small, and show only vague spatial patterns (note the differing scales). The errors for precipitation area also quite small, but can be locally larger, as in the pathological case illustrated below

### Paleo calendar adjustments

### Further examples

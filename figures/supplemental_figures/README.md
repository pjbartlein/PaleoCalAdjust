Supplemental Figures
-----------------

The figures in this folder illustrate the calendar effects on a selection of simulations used to exercise the calendar-effect adjustment programs (`cal_adjust_PMIP.f90` (version 1) and `call_adjust.f90` (version 1.1).  Each figure has four panels:  

- the long-term mean differences (experiment minus control, sometimes referred to as "anomalies") for each simulation calculated using adjusted paleo data (upper left); 
- the long-term mean differences (experiment minus control) calculated using unadjusted data (upper right);
- the calendar effect, or the adjusted long-term mean differences (upper left) minus the unadjusted long-term mean differences (upper right). For example, positive differences on this panel for near-surface air temperature (*tas*) indicate where the adjusted temperatures would be higher than unadjusted ones, and similarly for precipitation rate (*pr*) where adjusted data would be greater than unadjusted data. (lower left); and
- long-term mean difference sign reversals (which indicate regions where the magnitude of the calendar effect is large enough to reverse the sign of the long-term mean difference) (bottom right). 

For *AClim* files, with time dimension lengths of 12, the calendar-effect maps are simple differences.  For *Amon* files, with varying time dimensions (often 1200, i.e. 100 years of monthly data), the maps show differences in the long-term means (for each month of the year, i.e. the "climatologies") of experiment minus control monthly values.  In other words, the *Aclim* show data that were averaged first, and then differenced, while the *Amon* maps are based on data that were differenced first and then averaged.

## Discussion

The figures here show "real" calendar effects on paleoclimatic simulations, in contrast to the "pure" effects shown in Figs. 11-13 of Bartlein and Shafer (2019).  There, the calendar effects were shown using present-day data, in effect assuming that the experiment-minus-control long-term mean difference was zero.  In practice, the calendar effect is influenced both by the nature of the annual cycle (with greater effects on variables with higher-amplitude annual cycles) and by the long-term mean differences.

In the maps shown here, the range of values plotted are the same for the long-term mean differences as for the calendar effects.  This allows one to evaluate the magnitude of the calendar effect relative to the climate changes that are the focus of paleo simulations.

#### Calendar effects on temperature ####

	tas_Aclim_CCSM4_midHolocene.pdf
	tas_Aclim_CNRM-CM5_midHolocene.pdf
	tas_Aclim_MPI-ESM-P_midHolocene.pdf
		
	tas_Amon_CCSM4_midHolocene.pdf
	tas_Amon_CNRM-CM5_midHolocene.pdf
	tas_Amon_MPI-ESM-P_midHolocene.pdf

	tas_Amon_IPSL-CM6-LR_midHolocene.pdf
	tas_Amon_IPSL-CM6-LR_lig127k.pdf

There are broad similarities across individual models in the calendar effects (lower left), and resemble the "pure" calendar effects shown in Fig. 11 and 12 of Bartlein and Shafer (2019), including those for the *midHolocene* simulations as well as a single *lig127k* simulation of near-surface air temperature (and that one is also consistent with the *midHolocene* and *lig127k* temperature-effect difference on Fig. 11 of Bartlein and Shafer).  The larger calendar effects can be seen to rival in size the paleo minus *piControl* long-term mean differences (upper left and upper right).  In many cases the calendar effects simply reinforce the long-term mean differences, but in some cases they are large enough to *change the sign* of the difference (see the long-term mean difference sign reversal panels, lower right).  The spatial patterns of the sign-change maps for temperature have two components:  1) large, or continental-scale patches that indicate where and when during the year the calendar effect could completely reverse the sense of change between a paleo and control experiment, and 2) long looping, or hollow cell-like features that indicate where the location of zero isopleths of the long-term mean differences are sensitive to the calendar effect, i.e. "blobular" shaped differences expanding or contracting.

The three *midHolocene* simulations of *tas* with *Aclim* files (CCSM4, CNRM-CM5, and MPI-ESM-P) can be compared with those with *Amon* files.  In general, the calendar-effect patterns are quite similar, with perhaps slightly stronger map patterns expressed by data from the *Amon* files, which is probably related to the order of differencing and averaging.  (This is more evident in the two *pr* simulations with both *Aclim* and *Amon* files, see below.)  However, there is little practical difference in the patterns that emerge from the two kinds of files.  The two PMIP4 simulations with IPSL-CM6A-LR, *midHolocene* and *lig12ka* show calendar-effect differences that strongly resemble those on Fig. 11 of Bartlein and Shafer despite large differences in the amplitude and sign of the long-term mean differences, particularly in May and June. 

Comparisons among the different *tas* simulations suggests that the calendar effects vary little among models with larger amplitude "experiment minus control" long-term mean differences and those with smaller amplitude differences. This implies that adjusting for the calendar effect does not suppress inter-model (or inter "warm-climate" experiment) differences (although this is a very small sample, both in terms of number of models and variables).

Overall, the calendar effect for temperature shows very strong interhemispheric contrasts, as well as noticeable continent-ocean contrasts, particularly from September through January.  The sense of the effect on temperature during these months, positive in the Northern Hemisphere, negative in the Southern, could, if not adjusted for, give the reader a sense of reduced interhemispheric contrast in these particular simulations. In particular the impression of higher positive long-term mean differences at high Southern Hemisphere latitudes in the austral spring and summer is reduced in the adjusted data, while the continent-ocean contrast is amplified in the adjusted data in the Northern Hemisphere.

#### Calendar effects on energy-balance variables ####

	rsnet_Aclim_CCSM4_midHolocene.pdf
	rlnet_Aclim_CCSM4_midHolocene.pdf
	rnet_Aclim_CCSM4_midHolocene.pdf
	hfls_Aclim_CCSM4_midHolocene.pdf
	hfss_Aclim_CCSM4_midHolocene.pdf

#### Calendar effects on atmospheric circulation-related variables ####

	zg_Aclim_CCSM4_midHolocene.pdf
	psl_Aclim_CCSM4_midHolocene.pdf
	wap_Aclim_CCSM4_midHolocene.pdf
	clt_Aclim_CCSM4_midHolocene.pdf

#### Calendar effects on precipitation ####

	pr_Aclim_CCSM4_midHolocene.pdf
	pr_Amon_CCSM4_midHolocene.pdf
	pr_Amon_MPI-ESM-P_midHolocene.pdf

The three *midHolocene* precipitation-rate simulations generally show smaller spatial-scale long-term mean difference sign-reversal patterns than do those for temperature, reflecting both the smaller spatial scale of the map-patterns of precipitation, and the tendency for the long-term mean differences to show a combination of amplification and damping of broad-scale features (like the monsoons) in place, as well as some shifts in the location of features.  For example, the calendar-effect maps show a consistent pattern of decreases in N. African monsoon precipitation in the adjusted data in July and August, and increases in September, October and November.  In contrast to temperature, these changes do not produce sign changes in the long-term mean differences that are of continental scale and consistent from month to month.  The sign-change patterns are still large enough to influence data-model comparisons at regional scales.

#### Calendar effects on SSTs, sea ice and salinity variables ####


	tos_Aclim_CCSM4_midHolocene.pdf
	sic_Aclim_CCSM4_midHolocene.pdf
	sic_Aclim_CCSM4_midHolocene_polar.pdf
	sos_Aclim_CCSM4_midHolocene.pdf






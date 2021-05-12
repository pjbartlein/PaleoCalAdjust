program demo_02_adjust_1yr
! demonstrates (paleo) calendar adjustment for one year (of monthly values)
! 6 ka example using a 365-day "noleaps" calendar
    
use pseudo_daily_interp_subs    ! pseudo-daily interpolation module
use mp_interp_epstein_subs      ! Epstein pseudo-daily interpolation subroutines

implicit none

integer, parameter  :: nm=12, nd=365, ny=1    ! number of months in year, number of days in year, number of years
integer, parameter  :: nt=nm * ny               ! total number of months

! month length, middle, beginning and ending days (from month_length.f90), present (1950 CE) and 6 ka
real(8)         :: rmonlen_00(nm), rmonmid_00(nm)
real(8)         :: rmonlen_06(nm), rmonmid_06(nm), rmonbeg_06(nm), rmonend_06(nm)

! number of days in year
integer(4)      :: ndays(ny)

! data, "raw" (i.e. on 06ka calendar) and adjusted to 6 ka calendar
real(8)         :: ym_06ka_cal(nm), ym_adj_cal(nm)

! weighted annual averages
real(8)         :: yann_06ka_cal, yann_adj_cal

integer(4)      :: max_nctrl_in         ! maximum number of inner intervals (months) in an outer interval (years)
integer(4)      :: max_ntargs_in        ! maximum number of subintervals (days) in an outer interval (years)

! pseudo-daily interpolated data and monthly means of interpolated data
real(8)         :: ydh(nd), ym_int(nm), yfill

! smoothing parameters for multi-year pseudo-daily interpolation
logical                 :: smooth=.false.
logical                 :: no_negatives = .false. ! restrict pseudo-daily interpolated values to positive values?    
logical                 :: match_mean = .true.    ! force monthly means of interpolated values to equal input
real(8)                 :: tol = 0.01             ! tolerance value for enforce_mean()

! month-length data, 0 ka
data rmonlen_00 / 31.0d0,  28.0d0,  31.0d0,  30.0d0,  31.0d0,  30.0d0,  31.0d0,  31.0d0,  30.0d0,  31.0d0,  30.0d0,  31.0d0 /
data rmonmid_00 / 15.5d0,  45.0d0,  74.5d0, 105.0d0, 135.5d0, 166.0d0, 196.5d0, 227.5d0, 258.0d0, 288.5d0, 319.0d0, 349.5d0 /

! paleo month-length data (e.g. 6 ka)

data rmonlen_06 / 32.48411855d0,  29.53995521d0,  32.47169095d0,  30.80519959d0,  30.97070708d0,  29.16304027d0, & 
    29.54741870d0,  29.34677894d0,  28.62575603d0,  30.18400184d0,  30.01008258d0,  31.85125027d0 / 
data rmonmid_06 / 12.11927457d0,  43.16312612d0,  74.21920341d0, 105.9122942d0,  136.8461822d0,  166.89607396d0, &
    196.1953744d0,  225.5599778d0,  254.4893412d0,  283.869959d0,  313.9840274d0,  344.9283972 /
data rmonbeg_06 / -4.05417526d0,  28.40375842d0,  57.92559617d0,  90.40319079d0, 121.23431371d0, 152.22542876d0, & 
    181.38183671d0, 210.90212219d0, 240.22885383d0, 268.86102212d0, 299.07135244d0, 329.10071578d0 / 
data rmonend_06 / 28.40375842d0,  57.92559617d0,  90.40319079d0, 121.23431371d0, 152.22542876d0, 181.38183671d0, & 
    210.90212219d0, 240.22885383d0, 268.86102212d0, 299.07135244d0, 329.10071578d0, 360.94582474d0 /

! data to adjust -- example tas values c65r36 from TraCE-21ka data
data ym_06ka_cal /-7.68829d0, -2.08282d0, -2.47379d0, 0.732910d0, 3.59814d0, 10.4207d0, &
    15.6256d0, 13.9701d0, 10.3976d0, 7.65222d0, -1.86514d0, -4.61826d0 / 
yfill = 1.0e32
ndays(1) = nd
max_nctrl_in = 12 
max_ntargs_in = 366 

write (out_unit,'(a)') "Paleo calendar adjustment of monthly means:"

write (out_unit,'(/a)') "Input data (monthly mean tas on 0 ka 365-day noleaps calendar):"
write (out_unit,'("  ym_06ka_cal: ",12f9.3)') ym_06ka_cal

! weighted annual mean of monthly input
call ann_wmean(nm, ym_06ka_cal, rmonlen_00, yann_06ka_cal)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') yann_06ka_cal

! interpolate monthly data to pseudo-daily values

call mp_interp_epstein(ny, nm, nt, ym_06ka_cal, yfill, rmonmid_00, int(rmonlen_00), &
    smooth, no_negatives, match_mean, tol, nd, max_nctrl_in, max_ntargs_in, &
    ydh, ym_int)

! reaggregate pseudo-daily values to new calendar
call day_to_rmon_ts(ny, ndays, rmonbeg_06, rmonend_06, nd, ydh, yfill, ym_adj_cal)

write (out_unit,'(/a)') "Paleo month-length adjusted values on 6ka calendar:"
write (out_unit,'("   ym_adj_cal: ",12f9.3)') ym_adj_cal

! weighted annual mean of adjusted data
call ann_wmean(nm, ym_adj_cal, rmonlen_06, yann_adj_cal)
write (out_unit,'("Weighted annual mean of adjusted data: ", f12.6)') yann_adj_cal

write (out_unit,'(/a)') "Calendar effect::"
write (out_unit,'("   ym_adj_cal: ",12f9.3)') ym_adj_cal - ym_06ka_cal

end program
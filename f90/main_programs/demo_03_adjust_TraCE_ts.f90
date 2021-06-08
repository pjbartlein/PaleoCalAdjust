program demo_03_adjust_TraCE_ts
! adjust a single TraCE-21ka grid-point time series
! interpolate to pseudo-daily values, and reaggregate using appropriate paleo month lengths
! this verion uses the Harzallah iterative-spline interpolation method
    
use pseudo_daily_interp_subs
use mp_interp_harzallah_subs

implicit none

integer(4), parameter   :: nm=12, ny=22040              ! number of months in year, number of days in year, number of years
integer(4), parameter   :: nt=nm*ny                     ! total number of months

! month lengths, present (1950 CE) (implied by TraCE-21ka noleap calendar)
integer(4)              :: imonlen_00(nm)
real(8)                 :: rmonmid_00(nm)
integer(4), allocatable :: imonlen(:)

! input monthly time series
real(8), allocatable    :: ym(:)        ! (nt = ny*nm) input monthly data for a single grid cell
real(8), allocatable    :: yrbp(:)      ! (ny) time (-ka)
real(8)                 :: yfill        ! fill value

! month- and year-length variables
integer(4)              :: iyrbp        ! integer YearBP
integer(4)              :: iyrce        ! integer YearCE
real(8), allocatable    :: rmonlen(:,:),rmonmid(:,:)    ! (ny,nm) real-value month lengths and mid days (from month_length.f90)
real(8), allocatable    :: rmonbeg(:,:),rmonend(:,:)    ! (ny,nm) real-value month beginning and ending days (from month_length.f90)
real(8), allocatable    :: VE_day(:)                    ! (ny) vernal equinox day
real(8), allocatable    :: SS_day(:)                    ! (ny) summer solstice day
integer(4), allocatable :: ndays(:)                     ! (ny) number of days in year

! pseudo-daily values
integer(4)              :: ndtot                ! total number of days
real(8), allocatable    :: x_ctrl(:)            ! (nt) real-value control points for pseudo-daily interpolation
real(8), allocatable    :: x_targ(:)            ! (ndtot) real_value target points for pseudo-daily interpolation
real(8)                 :: yrdy                 ! decimal day time (BP)
real(8), allocatable    :: ydh(:)               ! pseudo-daily values (ndtot)

! output calendar-adjusted monthly averages
real(8), allocatable    :: ym_adj(:)            ! calendar-adjusted monthly means
real(8), allocatable    :: ym_int(:)            ! monthly means of calendar-adjusted pseudo-daily values

!  parameters for Harzallah pseudo-daily interpolation
logical                 :: no_negatives = .false.   ! restrict pseudo-daily interpolated values to positive values?
logical                 :: match_mean = .true.      ! force monthly means of interpolated values to equal input
real(8)                 :: tol = 0.01               ! tolerance value for enforce_mean()

integer(4)              :: spline_case = 2      ! spline used for Harzallah interpolation (1 = cubic, 2 = PCHIP)
integer(4)              :: npad = 2             ! padding with annual chunks of data, Harzallah interpolation
integer(4)              :: max_nctrl_in         ! maximum number of inner intervals (months) in an outer interval (years)
integer(4)              :: max_ntargs_in        ! maximum number of subintervals (days) in an outer interval (years)


integer(4)              :: n, m, i, nd

character(2048)         :: datapath, infile, dailyfile, outfile, monlenpath, monfile
character(128)          :: header_out
character(32)           :: voutname='tas'
character(3)            :: monname_JanDec(nm)           ! month names (for header)
character(1)            :: header_in

data imonlen_00 /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
data rmonmid_00 / 15.5d0,  45.0d0,  74.5d0, 105.0d0, 135.5d0, 166.0d0, 196.5d0, 227.5d0, 258.0d0, 288.5d0, 319.0d0, 349.5d0 /
data monname_JanDec/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

datapath = "/Projects/Calendar/PaleoCalAdjust/data/TraCE_example/"
infile = "TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv"
!dailyfile = "TraCE_c30r40_tas_land_monlen0ka_Jan-Dec_ts_daily.csv" ! uncomment to see daily data
outfile = "TraCE_c30r40_tas_land_monlenadj_Jan-Dec.csv"

monlenpath = "/Projects/Calendar/PaleoCalAdjust/data/month_lengths/"
monfile = "tr21_cal_noleap_rmonlen.csv"

! open output files and write headers
open (3, file=trim(datapath)//trim(outfile))
header_out="age "
do m=1,nm
    header_out=trim(header_out)//", "//trim(voutname)//"_"//trim(monname_JanDec(m))
end do
write (3,'(a)') trim(header_out)
open (10, file=trim(datapath)//"debug.dat")

! allocate large arrays
allocate (yrbp(ny),ym(nt))
allocate (rmonlen(ny,nm),rmonmid(ny,nm),rmonbeg(ny,nm),rmonend(ny,nm),imonlen(nt),VE_day(ny),SS_day(ny),ndays(ny))
allocate (ym_adj(ny * nm), ym_int(ny * nm))

! read the input data
write (*,*) "Reading input data..."
open (1, file=trim(datapath)//trim(infile))
read (1,*) header_in
do n=1,ny
    read (1,*) yrbp(n),(ym((n-1)*nm+m),m=1,nm)
end do
close (1)
yfill = 1.0e32

! replicate 0 ka month lengths over the input time series -- all years have the same noleap month lengths
do n=1,ny
    do m=1,nm
        imonlen((n-1)*nm + m) = imonlen_00(m)
    end do
end do

! read the month-length variables and get total number of days
open (1, file=trim(monlenpath)//trim(monfile))
read (1,'(a)') header_in
ndtot = 0
do n=1,ny
    read (1,*) iyrbp,iyrce,rmonlen(n,1:nm),rmonmid(n,1:nm),rmonbeg(n,1:nm),rmonend(n,1:nm),VE_day(n),SS_day(n),ndays(n)
    ndtot = ndtot + ndays(n)
end do
close (1)

! generate control and target points for interpolation
allocate (x_ctrl(nt), x_targ(ndtot))
    
! control points
x_ctrl(1:nm) = rmonmid_00(1:nm)
nd = 0
do n = 2,ny
    nd = nd + ndays(n - 1)
    do m = 1,nm
        i = (n-1)*nm + m
        x_ctrl(i) = dble(nd) + rmonmid_00(m)
    end do
end do
!write (10, '(12f12.2)') x_ctrl

! target points
do i = 1, ndtot
    x_targ(i) = dble(i)
end do
!write (10, '(10f12.2)') x_targ
max_nctrl_in = 12 + 2 * npad        ! include padding
max_ntargs_in = 366 + 2 * npad * 31  ! include padding

write (*,'("ny, nt, ndtot: ",3i8)') ny, nt, ndtot

! allocate array of pseudo-daily values
allocate (ydh(ndtot))

! do the pseudo-daily interpolation
write (*,*) "Interpolating..."

call mp_interp_harzallah(ny, nm, nt, ym, yfill, x_ctrl, imonlen, &
    spline_case, npad, no_negatives, match_mean, tol, ndtot, x_targ, max_nctrl_in, max_ntargs_in, &
    ydh, ym_int)

!! write out the daily data
!write (*,*) "Writing interpolated pseudo-daily data..."
!open (2, file=trim(datapath)//trim(dailyfile))
!write (2,'(a)') "YearBP, day, YrDy, tas"
!do n=1,ny
!    do i=1,nd
!        iyrbp=int(yrbp(n)*1000.0d0)
!        yrdy=dble(iyrbp)+dble(i-1)/dble(nd)
!        write (2,'(i6,", ",i3,",",f16.8,", ",g14.6)') iyrbp,i,yrdy,yd((n-1)*nd + i)
!    end do
!end do
!close (2)

! calendar adjustment
write (*,*) "Calendar adjustment..."

call day_to_rmon_ts(ny, ndays, rmonbeg, rmonend, ndtot, ydh, yfill, ym_adj)

do n=1,ny
    write (3,'(f8.3,12(", ",g14.6))') yrbp(n),(ym_adj((n-1)*nm + m), m=1,nm)
end do

write (*,'(a)') "Done (demo_03_adjust_TraCE_ts)"

end program

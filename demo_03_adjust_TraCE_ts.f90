program demo_03_adjust_TraCE_ts
! adjust a single TraCE-21ka grid-point time series
! interpolate to pseudo-daily values, and reaggregate using appropriate paleo month lengths
    
use calendar_effects

implicit none

integer(4), parameter   :: nm=12, nd=365, ny=22040      ! number of months in year, number of days in year, number of years
integer(4), parameter   :: nt=nm*ny                     ! total number of months

! integer-value month lengths, present (1950 CE) (implied by TraCE-21ka noleap calendar)
integer(4)         :: imonlen_00(nm), imonlen(nt)

! input monthly time series
real(8), allocatable    :: xm(:)        ! input monthly data for a single grid cell (nt = ny*nm)
real(8), allocatable    :: yrbp(:)      ! time (-ka) (ny)
real(8)                 :: vfill        ! fill value

! month- and year-length variables
integer(4)              :: iyrbp        ! integer YearBP
integer(4)              :: iyrce        ! integer YearCE
real(8), allocatable    :: rmonlen(:,:),rmonmid(:,:)    ! ny*nm real-value month lengths and mid days (from month_length.f90)
real(8), allocatable    :: rmonbeg(:,:),rmonend(:,:)    ! ny*nm real-value month beginning and ending days (from month_length.f90)
real(8), allocatable    :: VE_day(:)                    ! ny vernal equinox day
integer(4), allocatable :: ndays(:)                     ! ny number of days in year

! pseudo-daily values
integer(4)              :: ndtot                        ! total number of days
real(8)                 :: yrdy         ! decimal day time (BP)
real(8), allocatable    :: xd(:)                        ! pseudo-daily values (ndtot)

! output calendar-adjusted monthly averages
real(8), allocatable    :: xm_adj(:)                  ! calendar-adjusted monthly means

! smoothing parameters for multi-year pseudo-daily interpolation
integer(4)              :: nw_tmp=21, nsw_tmp=20    ! smoothing parameters
logical                 :: smooth=.true., restore=.true.
logical                 :: no_negatives = .false. ! restrict pseudo-daily interpolated values to positive values?

integer(4)              :: n, m, i

character(2048)         :: datapath, monfile, infile, dailyfile, outfile, debugfile
character(128)          :: header_out
character(32)           :: voutname='tas'
character(3)            :: monname_JanDec(nm)           ! month names (for header)
character(1)            :: header_in

data imonlen_00 /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
data monname_JanDec/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

datapath = "/Projects/Calendar/data/work02/"
monfile = "tr21_cal_noleap_rmonlen.csv"
infile = "TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv"
dailyfile = "TraCE_c30r40_tas_land_monlen0ka_Jan-Dec_ts_daily_v2.csv"
outfile = "TraCE_c30r40_tas_land_monlenadj_Jan-Dec_v2.csv"
debugfile="debug_demo_03_c30r40_v2.dat"

! open output files and write headers
open (3, file=trim(datapath)//trim(outfile))
header_out="age "
do m=1,nm
    header_out=trim(header_out)//", "//trim(voutname)//"_"//trim(monname_JanDec(m))
end do
write (3,'(a)') trim(header_out)
open (10, file=trim(datapath)//trim(debugfile))

! allocate large arrays
allocate (yrbp(ny),xm(nt))
allocate (rmonlen(ny,nm),rmonmid(ny,nm),rmonbeg(ny,nm),rmonend(ny,nm),VE_day(ny),ndays(ny))
allocate (xm_adj(ny * nm))

! read the input data
write (*,*) "Reading input data..."
open (1, file=trim(datapath)//trim(infile))
read (1,*) header_in
do n=1,ny
    read (1,*) yrbp(n),(xm((n-1)*nm+m),m=1,nm)
    !write (10,'(f8.3, 12f8.2)') yrbp(n),(xm((n-1)*nm+m),m=1,nm)
end do
close (1)
vfill = 1.0e32

! replicate 0 ka month lengths over the input time series -- all years have the same noleap month lengths
do n=1,ny
    do m=1,nm
        imonlen((n-1)*nm + m) = imonlen_00(m)
        !write (10,*) n,m,(n-1)*nm + m,imonlen((n-1)*nm + m)
    end do
end do

! read the month-length variables and get total number of days
open (1, file=trim(datapath)//trim(monfile))
read (1,'(a)') header_in
ndtot = 0
do n=1,ny
    read (1,*) iyrbp,iyrce,rmonlen(n,1:nm),rmonmid(n,1:nm),rmonbeg(n,1:nm),rmonend(n,1:nm),VE_day(n),ndays(n)
    ndtot = ndtot + ndays(n)
    !write (10,'(i8,", ",i8,48(", ",f12.8),", ",f10.6,", ",i4,", ",i9)') &
    !    iyrbp,iyrce,rmonlen(n,1:nm),rmonmid(n,1:nm),rmonbeg(n,1:nm),rmonend(n,1:nm),VE_day(n),ndays(n),ndtot
end do
close (1)

! allocate array of pseudo-daily values
allocate (xd(ndtot))

! do the pseudo-daily interpolation
! this version makes one pass over the input data, and lightly smooths the year-to-year discontinuities
write (*,*) "Interpolating..."
call mon_to_day_ts(nt, imonlen, xm, vfill, no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xd)

!! write out the daily data
!write (*,*) "Writing interpolated pseudo-daily data..."
!open (2, file=trim(datapath)//trim(dailyfile))
!write (2,'(a)') "YearBP, day, YrDy, tas"
!do n=1,ny
!    do i=1,nd
!        iyrbp=int(yrbp(n)*1000.0d0)
!        yrdy=dble(iyrbp)+dble(i-1)/dble(nd)
!        write (2,'(i6,", ",i3,",",f16.8,", ",g14.6)') iyrbp,i,yrdy,xd((n-1)*nd + i)
!    end do
!end do
!close (2)

! calendar adjustment
write (*,*) "Calendar adjustment..."

call day_to_mon_ts(ny, ndays, rmonbeg, rmonend, ndtot, xd, vfill, xm_adj)

do n=1,ny
    write (3,'(f8.3,12(", ",g14.6))') yrbp(n),(xm_adj((n-1)*nm + m), m=1,nm)
end do

end program

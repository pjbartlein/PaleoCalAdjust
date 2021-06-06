program cal_adjust
    
! Performs calendar (month-length) adjustments of data in a CMIP/PMIP-formatted netCDF file,
! creating a new netCDF file by copying dimension variables and global attributes from the input file.
! This version supports 3-D (e.g. longitude, latitude, time) and 4-D (e.g. longitude, latitude, level, time),
! long-term mean (AClim, Oclim), monthly (Amon, Omon, OImon), or daily input files.

! The program requires the modules: CMIP_netCDF_subs.f90, GISS_orbpar_subs.f90,
! GISS_srevents_subs.f90, month_length_subs.f90 and pseudo_daily_interp_subs.f90,
! mp_interp_epstein_subs, and mp_interp_harzallah_subs
    
! The program must be compiled with local netCDF support, and it will use OpenMP if available

! An info .csv file (with an appropriate header) containing the following information is read:
    
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

! (See below for the distinction between "simulation age" and "simulation year".)

! Author: Patrick J. Bartlein, Univ. of Oregon (bartlein@uoregon.edu), with contributions by S.L. Shafer (sshafer@usgs.gov)
! This program and related subroutines are part of PaleoCalAdjust:
!
! Version: 1.0 (original release)
!   - see P.J. Bartlein & S.L. Shafer (2019) Paleo calendar-effect adjustments in time-slice and transient
!         climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis,
!         Geoscientific Model Development, 12:3889–3913, doi: 10.5194/gmd-12-3889-2019
!   - available from GitHub:  https://github.com/pjbartlein/PaleoCalAdjust or Zenodo:  https://doi.org/10.5281/zenodo.1478824
!
! Version: 1.1 (this version:
!
!   - renamed main program to cal_adjust.f90 (from cal_adjust_PMIP.f90)
!   - modified the info file
!   - infofile name and path are now read as a command-line argument
!   - added Harzallah spline pseudo-daily interpolation approach for adjustment of "Amon"-type (timeseries) files
!   - added enforce_mean() subroutine (to require that the mean of the interpolated values should match that of the input data)
!   - added interpolation method attributes to adjusted files
!   - combined the calendar_effects_subs.f90 and pseudo_daily_interp_subs.f90 modules (as pseudo_daily_interp_subs.f90)
!
! Please cite Bartlein & Shafer (2019) if you use this code.

!    
! Last update: 2021-05-11

use CMIP_netCDF_subs
use month_length_subs
use pseudo_daily_interp_subs
use mp_interp_epstein_subs
use mp_interp_harzallah_subs
use omp_lib

implicit none

! Variables defined in modules
! integer(4)              :: nm      !number of months in year
! integer(4)              :: ny      !number of years

! Note:  "simulation ages" control the orbital parameters, "simulation years" determine the occurrence of leap years and 
! the date of occurrence of the vernal equinox.
! See Bartein & Shafer (2019) for further discussion of the difference between simulation ages and simulation years.

! simulation age-related variables (controls of orbital parameters)
! past simulation ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0
integer(4)              :: begageBP             ! beginning simulation age (year BP) (negative, e.g. 10 ka = -10000 BP)
integer(4)              :: endageBP             ! ending simulation age (year BP)
integer(4)              :: agestep              ! interval/step size between age calculations
integer(4)              :: nages                ! number of calendar years of simulation
integer(4), allocatable :: iageBP(:)            ! (ny) year BP 1950 (negative, e.g. 1900 CE = -50 years BP 1950)

! simulation year-related variables (controls of vernal equinox day and leap-year status)
! simulation years are specified as years CE
integer(4)              :: begyrCE              ! beginning calendar year of simulation for multi-year simulations at each age
integer(4)              :: nsimyrs              ! number of years of simulation
integer(4), allocatable :: iyearCE(:)           ! (ny) yearCE simulation year (e.g. 1850CE, 850CE, etc.)

! month-length variables
integer(4), allocatable :: imonlen_0ka(:,:)     ! (ny,nm) integer-value month lengths (0ka)
integer(4), allocatable :: imonmid_0ka(:,:)     ! (ny,nm) integer-value mid days (0ka)
integer(4), allocatable :: imonbeg_0ka(:,:)     ! (ny,nm) integer-value month beginning days (0ka)
integer(4), allocatable :: imonend_0ka(:,:)     ! (ny,nm) integer-value month ending days (0ka)
integer(4), allocatable :: imonlen(:,:)         ! (ny,nm) integer-value month lengths days (paleo)
integer(4), allocatable :: imonmid(:,:)         ! (ny,nm) integer-value mid days (paleo)
integer(4), allocatable :: imonbeg(:,:)         ! (ny,nm) integer-value month beginning days (paleo)
integer(4), allocatable :: imonend(:,:)         ! (ny,nm) integer-value month ending days (paleo)
real(8), allocatable    :: rmonlen_0ka(:,:)     ! (ny,nm) real-value month lengths days (0ka)
real(8), allocatable    :: rmonmid_0ka(:,:)     ! (ny,nm) real-value mid days (0ka)
real(8), allocatable    :: rmonbeg_0ka(:,:)     ! (ny,nm) real-value month beginning days (0ka)
real(8), allocatable    :: rmonend_0ka(:,:)     ! (ny,nm) real-value month ending days days (0ka)
real(8), allocatable    :: rmonlen(:,:)         ! (ny,nm) real-value month lengths days (paleo)
real(8), allocatable    :: rmonmid(:,:)         ! (ny,nm) real-value mid days (paleo)
real(8), allocatable    :: rmonbeg(:,:)         ! (ny,nm) real-value month beginning days (paleo)
real(8), allocatable    :: rmonend(:,:)         ! (ny,nm) real-value month ending days (paleo)
real(8), allocatable    :: VE_day(:)            ! (ny) vernal equinox day in each simulation year
real(8), allocatable    :: SS_day(:)            ! (ny) (northern) summer solstice day in each simulation year
integer(4), allocatable :: ndays(:)             ! (ny) number of days in in each simulation year
real(8), allocatable    :: VEtoSS(:),SStoAE(:)  ! (nages*nsimyrs) number of days VE to SS, SS to AE
real(8), allocatable    :: AEtoWS(:),WStoVE(:)  ! (nages*nsimyrs) number of days AE to WS, WS to VE

integer(4), allocatable :: imonlen_0ka_ts(:)    ! (nt) integer-value month lengths at present (0 ka) as time series
real(8), allocatable    :: rmonmid_0ka_ts(:)    ! (nt) real-value mid days at present (0 ka) as time series
real(8), allocatable    :: rmonlen_ts(:)        ! (nt) real-value paleo month lengths as time series
real(8), allocatable    :: rmonmid_ts(:)        ! (nt) real-value paleo month mid days as time series
real(8), allocatable    :: rmonbeg_ts(:)        ! (nt) real-value paleo month beginning days as time series
real(8), allocatable    :: rmonend_ts(:)        ! (nt) real-value paleo month ending days as time series
integer(4), allocatable :: imonmid_ts(:)        ! (nt) integer-value paleo month mid days as time series
integer(4), allocatable :: imonbeg_ts(:)        ! (nt) integer-value paleo month beginning days as time series
integer(4), allocatable :: imonend_ts(:)        ! (nt) integer-value paleo month ending as time series
integer(4), allocatable :: ndays_ts(:)          ! (nt) integer-value times series of paleo year lengths
real(8), allocatable    :: x_ctrl(:)            ! (nt) real-value control points for pseudo-daily interpolation
real(8), allocatable    :: x_targ(:)            ! (ndtot) real_value target points for pseudo-daily interpolation
integer(8), allocatable :: idaynum(:)           ! (integer) day number
character(256)          :: time_comment         ! source of new (adjusted) monthly time values
real(8), allocatable    :: mon_time(:)          ! (nt) new (adjusterd) monthly time values for daily-input files
real(8), allocatable    :: mon_time_bnds(:,:)   ! (2,nt) new (adjusted) monthly time-bounds values for daily input files

! other names
character(32)           :: calendar_type        ! CF calendar type (e.g. "noleap", "proleptic_gregorian", etc.)

! strings
character(12)            :: activity             ! simple label (e.g. "PMIP3" or "PMIP4"), or some other label (e.g. "transient")
character(64)           :: variable             ! variable name
character(8)            :: time_freq            ! CMIP/PMIP-style output time-frequency type (e.g. "Amon", "Aclim", "day", "Oclim" ...)
character(14)           :: tol_string           ! enforce_mean() tolerance value as string
character(14)           :: buffer               ! string-conversion buffer

! data
! nlon = number of longitudes, nlat = number of latitudes, nlev = number of levels, nt = number of observations in input)
integer(4)              :: ny, nt               ! number of years and total nmber of months (nt = ny*nm)
integer(4)              :: ndtot                ! total number of days and years in simulation
real(4), allocatable    :: var3d_in(:,:,:)      ! 3-D (e.g. nlon,nlat,nt or nlon,nlat,ndtot) input data (either daily or monthly)
real(4), allocatable    :: var3d_out(:,:,:)     ! 3-D (e.g. nlon,nlat,nt) output adjusted data (monthly)
real(8), allocatable    :: ydh(:,:)             ! 2-D (e.g. ivar_dimlen(1),ndtot) pseudo- or actual daily data
real(8), allocatable    :: ym_int(:,:)          ! 2-D (e.g. ivar_dimlen(1),nt) monthly means of interpolated data
real(8), allocatable    :: var_adj(:,:)         ! 2-D (e.g. ivar_dimlen(1),nt) adjusted data 
real(8), allocatable    :: var4d_in(:,:,:,:)    ! 4-D (e.g. nlon,nlat,nlev,nt) input data
real(4), allocatable    :: var4d_out(:,:,:,:)   ! 4-D (e.g. nlon,nlat,nlev,nt) output adjusted data
real(8)                 :: vfill                ! fill value

! interpolation method information
character(9)            :: interp_method        ! interpolation method ("Epstein" or "Harzallah")
character(5)            :: no_negatives_in      ! character true or false
character(5)            :: match_mean_in        ! character true or false
real(8)                 :: tol                  ! tolerance value for enforce_mean() subroutine

logical                 :: smooth=.true.        ! smooth year-to-year discontinuities, Epstein method
logical                 :: no_negatives=.false. ! constrain pseudo-daily interpolated values to positive values?
logical                 :: match_mean = .true.  ! constrain the monthly mean of interpolated values to match input (via enforce_mean())

integer(4)              :: spline_case = 2      ! spline used for Harzallah interpolation (1 = cubic, 2 = PCHIP)
integer(4)              :: npad = 2             ! padding with annual chunks of data, Harzallah interpolation
integer(4)              :: max_nctrl_in         ! maximum number of inner intervals (months) in an outer interval (years)
integer(4)              :: max_ntargs_in        ! maximum number of subintervals (days) in an outer interval (years)

! indices, etc.
integer(4)              :: n,m,j,k,l,i,nd       ! indices
integer(4)              :: max_threads          ! if OpenMP enabled
integer(4)              :: iostatus             ! IOSTAT value
integer(4)              :: infofile_len         ! length of the info file
integer(4)              :: infofile_status      ! info file read status

! file paths and names
character(2048)         :: source_path, source_file, ncfile_in, adjusted_path, adjusted_file, ncfile_out
character(2048)         :: infofile     ! infofile name
character(2048)         :: info_header

! timers
real(4)                 :: overall_secs, monlen_secs, define_secs, read_secs, allocate_secs
real(4)                 :: openmp_secs, write_secs, total_secs

! progress bar
integer(4)              :: ntotalpts, numberdone, nprogress, nextra
real(4)                 :: onepercent, nextonepercent

! if OpenMP enabled
max_threads = omp_get_max_threads()
write (*,'("OMP max_threads: ",i4)') max_threads
max_threads = max_threads !- 2 ! subtract 2 to free up threads for other computational processes if necessary
call omp_set_num_threads(max_threads)

! get the info-file name:
call get_command_argument (1, infofile, infofile_len, infofile_status)
if (infofile_status .ne. 0) then
    write (*,*) "Getting command argument failed with status = ", infofile_status
    stop "infofile"
end if
write (*,'("Info file: ", a)') trim(infofile)

!open(10, file="\Projects\Calendar\PaleoCalAdjust\data\debug.dat") ! open debug file if needed

! open the info file, and read source-file and adjusted-file paths
open (3, file=trim(infofile))
read (3,'(a)') info_header
write(*, '("Info_header: ", a80)') trim(info_header)

! main loop (over individual model-output files)
iostatus = 1
overall_secs = secnds(0.0)
do
    total_secs = secnds(0.0); monlen_secs = secnds(0.0)
    ! Step 1:  Read info file, allocate month-length arrays, etc.
    
    ! read a line from the info file
    write (*,'(120("*"))')
    read (3,*,iostat=iostatus) activity, variable, time_freq, calendar_type, &
        begageBP, endageBP, agestep, begyrCE, nsimyrs, interp_method, no_negatives_in, &
        match_mean_in, tol, source_path, source_file, adjusted_path, adjusted_file
    write (*,'("iostatus: ", i4)') iostatus
    if (iostatus.lt.0) then
        write (*,'(a)') "*** Done ***"
        exit
    end if
    
    write (*,'("     activity: ", a)') trim(activity)
    write (*,'("     variable: ", a)') trim(variable)
    write (*,'("    time_freq: ", a)') trim(time_freq)
    write (*,'("calendar_type: ", a)') trim(calendar_type)
    write (*,'("interp_method: ", a)') trim(interp_method)
    write (*,'("  source_path: ", a)') trim(source_path)
    write (*,'("  source_file: ", a)') trim(source_file)
    write (*,'("adjusted_path: ", a)') trim(adjusted_path)
    write (*,'("adjusted_file: ", a)') trim(adjusted_file)
    
    ncfile_in = trim(source_path)//trim(source_file)
    ncfile_out = trim(adjusted_path)//trim(adjusted_file)

    varinname = trim(variable); varoutname = varinname

    write (*,'("    ncfile_in: ",a)') trim(ncfile_in)
    write (*,'("   ncfile_out: ",a)') trim(ncfile_out)
    write (*,'(a, 1x, a, 1x, 5i7)') trim(variable),trim(calendar_type), begageBP, endageBP, agestep, begyrCE, nsimyrs
    
    ! set interpolation method and max number of control and target points per year
    select case (interp_method)
    case ("epstein", "Epstein")
        interp_method = "Epstein"
        max_nctrl_in = 12 
        max_ntargs_in = 366 
    case ("harzallah", "Harzallah")
        interp_method = "Harzallah"
        max_nctrl_in = 12 + 2 * npad            ! include padding
        max_ntargs_in = 366 + 2 * npad * 31         ! include padding
    case ("none", "None")
        interp_method = "none"
    end select
    
    ! set flags
    select case (no_negatives_in)
    case ("T", "true", "True", "TRUE")
        no_negatives = .true.
    case default
        no_negatives = .false.
    end select
    
    select case (match_mean_in)
    case ("T", "true", "True", "TRUE")
        match_mean = .true.
    case default
        match_mean = .false.
    end select
    
    ! convert tol to string
    write (buffer, '(g14.6)') tol
    read (buffer, *) tol_string
     
    write (*,'("no_negatives: ", l1)') no_negatives
    write (*,'("match_mean: ", l1)') match_mean
    write (*,'("tolerance: ", g14.6, 2x, a)') tol, trim(tol_string)

    ! allocate month-length arrays
    nages = (endageBP - begageBP)/agestep + 1
    ny = nages * nsimyrs
    nt = ny * nm
    write (*,'("nages, nsimyrs, nt: ",3i8)') nages, nsimyrs, nt
    allocate (iageBP(ny), iyearCE(ny))
    allocate (imonlen_0ka(ny,nm),imonmid_0ka(ny,nm),imonbeg_0ka(ny,nm),imonend_0ka(ny,nm))
    allocate (imonlen(ny,nm), imonmid(ny,nm), imonbeg(ny,nm), imonend(ny,nm))
    allocate (rmonlen_0ka(ny,nm), rmonmid_0ka(ny,nm), rmonbeg_0ka(ny,nm), rmonend_0ka(ny,nm))
    allocate (rmonlen(ny,nm), rmonmid(ny,nm), rmonbeg(ny,nm), rmonend(ny,nm))
    allocate (VE_day(ny), SS_day(ny), ndays(ny), imonlen_0ka_ts(nt), rmonmid_0ka_ts(nt), ndays_ts(nt))  
    allocate (VEtoSS(nages*nsimyrs),SStoAE(nages*nsimyrs),AEtoWS(nages*nsimyrs),WStoVE(nages*nsimyrs))
    allocate (imonmid_ts(nt), imonbeg_ts(nt),imonend_ts(nt), rmonlen_ts(nt), rmonmid_ts(nt), rmonbeg_ts(nt), rmonend_ts(nt))
    allocate (mon_time(nt), mon_time_bnds(2,nt))

    ! Step 2:  get month lengths and control and target values for interpolation

    ! 0 ka month lengths (used in pseudo-daily interpolation)
    if (trim(time_freq) .ne. 'day') then
        write (*,'(a)') "0 ka month lengths for pseudo-daily interpolation..."
        call get_month_lengths(calendar_type, 0, 0, agestep, nages, begyrCE, nsimyrs, &
            iageBP, iyearCE, imonlen_0ka, imonmid_0ka, imonbeg_0ka, imonend_0ka, & 
            rmonlen_0ka, rmonmid_0ka, rmonbeg_0ka, rmonend_0ka, VE_day, SS_day, ndays, &
            VEtoSS, SStoAE, AEtoWS, WStoVE)
    end if

    ! paleo month lengths (used for calendar adjustment)
    write (*,'(a)') "Paleo month-lengths for aggregation of daily data..."
    call get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, &
        iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays, &
        VEtoSS, SStoAE, AEtoWS, WStoVE)

    ! reshape month lengths into time series, and get cumulative number of days 
    ndtot = 0
    do n=1,ny
        ndtot = ndtot + ndays(n)
        do m=1,nm
            i = (n-1)*nm + m
            imonlen_0ka_ts(i) = imonlen_0ka(n,m)
            rmonmid_0ka_ts(i) = rmonmid_0ka(n,m)
            rmonlen_ts(i) = rmonlen(n,m)
            rmonmid_ts(i) = rmonmid(n,m); rmonbeg_ts(i) = rmonbeg(n,m); rmonend_ts(i) = rmonend(n,m)
            imonmid_ts(i) = imonmid(n,m); imonbeg_ts(i) = imonbeg(n,m); imonend_ts(i) = imonend(n,m)
            ndays_ts(i) = ndtot
        end do
    end do
    
    ! generate control and target points for interpolation
    allocate (x_ctrl(nt), x_targ(ndtot), idaynum(ndtot))
    
    ! control points
    x_ctrl(1:nm) = rmonmid_0ka_ts(1:nm)
    nd = 0
    do n = 2,ny
        nd = nd + ndays(n - 1)
        do m = 1,nm
            i = (n-1)*nm + m
            x_ctrl(i) = dble(nd) + rmonmid_0ka_ts(i)
        end do
    end do

    ! target points
    do i = 1, ndtot
        x_targ(i) = dble(i)
        idaynum(i) = i
    end do

    write (*,'("ny, nt, ndtot: ",3i8)') ny, nt, ndtot
    write (*,'(a,f7.2)') "Month-length time: ", secnds(monlen_secs)
    
    ! Step 3:  Open input and create output netCDF files 

    ! input netCDF file
    define_secs = secnds(0.0)
    print '(" nc_fname (in) = ",a)', trim(ncfile_in)
    call check( nf90_open(ncfile_in, nf90_nowrite, ncid_in) )
    if (print_out) print '("  ncid_in = ",i8)', ncid_in

    ! output netCDF file
    print '(" nc_fname (out) = ",a)', trim(ncfile_out)
    call check( nf90_create(ncfile_out, 0, ncid_out) )
    if (print_out) print '("  ncid_out = ",i8)', ncid_out
    
    ! Step 4:  Redefine the time-coordinate variable 

    ! redefine the time coordinate
    if (trim(time_freq) .eq. 'day') then
        time_comment = "paleo monthly mid, beginning and ending dates from original daily values"
        call new_time_day(ncid_in, ny, nm, nt, ndtot, &
            imonmid_ts, imonbeg_ts, imonend_ts, ndays_ts, mon_time, mon_time_bnds)
    else
        time_comment = "paleo monthly values from get_month_lengths()"
        call new_time_month(calendar_type, ncid_in, ny, nm, nt, &
            rmonmid_ts, rmonbeg_ts, rmonend_ts, ndays_ts, mon_time, mon_time_bnds)
    end if
    
    ! Step 5:  Define the new netCDF file
    
    ! get number of input variable dimensions, dimension names and lengths
    call get_var_diminfo(ncid_in, varinname, invar_ndim, invar_dimid, invar_dimlen, invar_dimname)

    ! define the new netCDF file, and copy dimension variables and global attributes
    call current_time(current)
    addattname = "paleo_calendar_adjustment"
    addatt = trim(current)//" paleo calendar adjustment by cal_adjust_PMIP.f90"
    call copy_dims_and_glatts(ncid_in, ncid_out, addattname, addatt, nt, &
        mon_time, mon_time_bnds, time_comment, varid_out)

    ! define the output (adjusted) variable, and copy attributes
    write (*,'(a)') "Defining (new) adjusted variable..."
    addvarattname = "paleo_calendar_adjustment"
    select case (trim(time_freq))
    case ('Aclim', 'Oclim', 'OIclim', 'LIclim')
        addvaratt = "long-term mean values adjusted for the appropriate paleo calendar, "
        addvaratt = trim(addvaratt)//" "//trim(interp_method)//" pseudo-daily interpolation"
        if (match_mean) addvaratt = trim(addvaratt)//", called enforce_mean() tolerance: "//trim(tol_string)
    case ('Amon', 'Omon', 'OImon', 'LImon')
        addvaratt = "monthly values adjusted for the appropriate paleo calendar, "
        addvaratt = trim(addvaratt)//" "//trim(interp_method)//" pseudo-daily interpolation"
        if (match_mean) addvaratt = trim(addvaratt)//", called enforce_mean() tolerance: "//trim(tol_string)
    case ('day')
        addvaratt = "daily values aggregated to monthly values using the appropriate paleo month lengths"
        interp_method = "none"
    case default 
        stop "defining new variable"
    end select

    ! define the output variable
    call define_outvar(ncid_in, ncid_out, varinname, varid_out, varoutname, addvarattname, addvaratt, varid_in) 
    write (*,'(a,f7.2)') "Define time: ", secnds(define_secs)

    ! Step 6:  Get the input variable to be adjusted
    
    ! allocate variables
    allocate_secs = secnds(0.0)
    write (*,'(a)') "Allocating arrays"
    select case(invar_ndim)
    case (3)
        write (*,*) invar_dimlen(1),invar_dimlen(2),nt,invar_dimlen(1)*invar_dimlen(2)*nt
        if (trim(time_freq) .eq. 'day') then
            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),ndtot))
        else
            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),nt))
        end if
        allocate(var3d_out(invar_dimlen(1),invar_dimlen(2),nt))
    case (4)
        write (*,*) invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt,invar_dimlen(1)*invar_dimlen(2)*invar_dimlen(3)*nt
        if (trim(time_freq) .eq. 'day') then
            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),ndtot))
        else
            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
        end if
        allocate(var4d_out(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
    case default
        stop "allocating variables"
    end select
    allocate(ydh(invar_dimlen(1),ndtot), ym_int(invar_dimlen(1),nt), var_adj(invar_dimlen(1),nt))
   write (*,'(a,f7.2)') "Allocate time: ", secnds(allocate_secs)
    
    ! get input data
    write (*,'(a)') "Reading input data..."
    read_secs = secnds(0.0)
    select case (invar_ndim)
    case (3)
        call check( nf90_get_var(ncid_in, varid_in, var3d_in) )
    case (4)
        call check( nf90_get_var(ncid_in, varid_in, var4d_in) )
    case default
        stop "input data"
    end select

    ! get _FillValue
    call check( nf90_get_att(ncid_in, varid_in, '_FillValue', vfill) )
    write (*,'("_FillValue:", g14.6)') vfill
    write (*,'(a,f7.2)') "Read time: ", secnds(read_secs)
    
    ! Step 7:  Get calendar-adjusted values
    
    openmp_secs = 0.0
    if (trim(time_freq) .eq.'day') then
        write (*,'(a)') "Aggregating daily data ..."
    else
        write (*,'(a)') "Interpolating and aggregating ..."
    end if
    
    openmp_secs = secnds(0.0)
    select case(invar_ndim)
    case (3) 
        ntotalpts = invar_dimlen(1) * invar_dimlen(2)
        write (*,'("Number of grid points = ",i6)') ntotalpts
    case (4)
        ntotalpts = invar_dimlen(1) * invar_dimlen(2) * invar_dimlen(3)
        write (*,'("Number of grid points x number of levels = ",i6)') ntotalpts
    end select
        
    numberdone = 0; onepercent = (ntotalpts/100); nextonepercent = onepercent; nprogress = 0
    write (*,'("|",9("---------+"),"---------|")'); write (*,'(a1,$)') " "
    
    ! Note: OpenMP seems to work best here if only innermost loop is parallelized
    select case (invar_ndim)
    case (3)
        do k=1,invar_dimlen(2)      
            ! if daily data, copy into ydh
            !if (trim(time_freq) .eq.'day') ydh(:,:) = var3d_in(:,k,:)
        
            !$omp parallel do
            do j=1,invar_dimlen(1)
                !write (*,'(/2i5)') j,k
                if (trim(time_freq) .eq.'day') ydh(j,:) = var3d_in(j,k,:)
                ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
                if (trim(time_freq) .ne. 'day') then
                    ! interpolate
                    select case (trim(interp_method))
                    case ('Epstein')
                        call mp_interp_epstein(ny, nm, nt, dble(var3d_in(j,k,:)), dble(vfill), x_ctrl, imonlen_0ka_ts, &
                            smooth, no_negatives, match_mean, tol, ndtot, max_nctrl_in, max_ntargs_in, &
                            ydh(j,:), ym_int(j,:))
                    case ('Harzallah')
                        call mp_interp_harzallah(ny, nm, nt, dble(var3d_in(j,k,:)), dble(vfill), x_ctrl, imonlen_0ka_ts, & 
                            spline_case, npad, no_negatives, match_mean, tol, ndtot, x_targ, max_nctrl_in, max_ntargs_in, &
                            ydh(j,:), ym_int(j,:))
                    case default
                        stop "interp_method"
                    end select
                    ! reaggregate daily data using correct calendar
                    call day_to_rmon_ts(ny,ndays,imonbeg,imonend,rmonbeg,rmonend,ndtot,idaynum,ydh(j,:), &
                        dble(var3d_in(j,k,:)),dble(vfill),var_adj(j,:))
                else
                    ! input data are already daily, so just reaggregate using correct calendar
                    call day_to_rmon_ts(ny,ndays,imonbeg,imonend,rmonbeg,rmonend,ndtot,idaynum,ydh(j,:), &
                        dble(var3d_in(j,k,:)),dble(vfill),var_adj(j,:))
                end if
            
                var3d_out(j,k,:)=sngl(var_adj(j,:))

            end do
            !$omp end parallel do
          
            ! update progress bar
            numberdone = numberdone + invar_dimlen(1)
            if (numberdone .ge. nextonepercent) then
                if (nprogress .lt. 99) then
                    write (*,'(a1,$)') "="
                    nextonepercent = nextonepercent + onepercent
                    nprogress = nprogress + 1
                    !write (*,*) numberdone, nextonepercent
                end if
            end if

        end do
        
    case(4)
        do l=1,invar_dimlen(3)
            do k=1,invar_dimlen(2)      
                ! if daily data, copy into ydh
                !if (trim(time_freq) .eq.'day') ydh(:,:) = var4d_in(:,k,l,:)
        
                !$omp parallel do
                do j=1,invar_dimlen(1)
                    !write (*,'(/3i5)') j,k,l
                    if (trim(time_freq) .eq.'day') ydh(j,:) = var4d_in(j,k,l,:)
                    ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
                    if (trim(time_freq) .ne. 'day') then
                        ! interpolate
                        select case (trim(interp_method))
                        case ('Epstein')
                            call mp_interp_epstein(ny, nm, nt, dble(var4d_in(j,k,l,:)), dble(vfill), x_ctrl, imonlen_0ka_ts, &
                                smooth, no_negatives, match_mean, tol, ndtot, max_nctrl_in, max_ntargs_in, &
                                ydh(j,:), ym_int(j,:))
                        case ('Harzallah')
                            call mp_interp_harzallah(ny, nm, nt, dble(var4d_in(j,k,l,:)), dble(vfill), x_ctrl, imonlen_0ka_ts, &
                                spline_case, npad, no_negatives, match_mean, tol, ndtot, x_targ, max_nctrl_in, max_ntargs_in, & 
                                ydh(j,:), ym_int(j,:))
                        case default
                            stop "interp_method"
                        end select
                        ! reaggregate daily data using correct calendar
                        call day_to_rmon_ts(ny,ndays,imonbeg,imonend,rmonbeg,rmonend,ndtot,idaynum,ydh(j,:), &
                            dble(var4d_in(j,k,l,:)),dble(vfill),var_adj(j,:))
                    else
                        ! input data are already daily, so just reaggregate using correct calendar
                        call day_to_rmon_ts(ny,ndays,imonbeg,imonend,rmonbeg,rmonend,ndtot,idaynum,ydh(j,:), &
                            dble(var4d_in(j,k,l,:)),dble(vfill),var_adj(j,:))
                    end if
            
                    var4d_out(j,k,l,:)=sngl(var_adj(j,:))

                end do
                !$omp end parallel do
          
                ! update progress bar
                numberdone = numberdone + invar_dimlen(1)
                if (numberdone .ge. nextonepercent) then
                    if (nprogress .lt. 99) then
                        write (*,'(a1,$)') "="
                        nextonepercent = nextonepercent + onepercent
                        nprogress = nprogress + 1
                        !write (*,*) numberdone, nextonepercent
                    end if
                end if

            end do
        end do
    case default
        stop "adjusting"
    end select
    
    ! finish progress bar if short
    if (nprogress .lt. 99) then
        do nextra=1,(99-nprogress)
            write (*,'(a1,$)') "="
        end do
    end if
    write (*,'(a)') " "
    
    write (*,'(a,f7.2)') "OpenMP time: ", secnds(openmp_secs)

    ! Step 8:  Write out the adjusted data, and close the output netCDF file.
    
    ! write out adjusted data
    write_secs = secnds(0.0)
    write (*,'(/a)') "Writing adjusted data..."
    select case(invar_ndim)
    case (3)
        call check( nf90_put_var(ncid_out, varid_out, var3d_out) )
    case (4)
        call check( nf90_put_var(ncid_out, varid_out, var4d_out) )
    case default
        stop "adjusted output"
    end select
    
    ! write out month-length variables
    call check( nf90_put_var(ncid_out, rmonlenid, rmonlen_ts) )
    call check( nf90_put_var(ncid_out, rmonmidid, rmonmid_ts) )
    call check( nf90_put_var(ncid_out, rmonbegid, rmonbeg_ts) )
    call check( nf90_put_var(ncid_out, rmonendid, rmonend_ts) )
    
    ! close the output file
    call check( nf90_close(ncid_out) )
    write (*,'(a,f7.2)') "write time: ", secnds(write_secs)

    deallocate (iageBP, iyearCE)
    deallocate (imonlen_0ka,imonmid_0ka,imonbeg_0ka,imonend_0ka)
    deallocate (imonlen, imonmid, imonbeg, imonend)
    deallocate (rmonlen, rmonmid, rmonbeg, rmonend)
    deallocate (VE_day, SS_day, ndays, imonlen_0ka_ts, rmonmid_0ka_ts, ndays_ts) 
    deallocate (VEtoSS,SStoAE,AEtoWS,WStoVE)
    deallocate (imonmid_ts, imonbeg_ts,imonend_ts, rmonlen_ts, rmonmid_ts, rmonbeg_ts, rmonend_ts)
    deallocate (rmonlen_0ka, rmonmid_0ka, rmonbeg_0ka, rmonend_0ka)
    deallocate (x_ctrl, x_targ, idaynum)
    deallocate (mon_time, mon_time_bnds)
    select case (invar_ndim)
    case (3)
        deallocate (var3d_in, var3d_out)
    case (4)
        deallocate (var4d_in, var4d_out)
    case default
        stop "deallocating"
    end select
    deallocate (ydh, ym_int, var_adj)

    write (*,'(a)') " "
    write (*,'(a,f7.2)') "Total time: ", secnds(Total_secs)
    
end do

    write (*,'(a,f7.2)') "Overall time: ", secnds(overall_secs)
    write (*,'(a)') " "

end program

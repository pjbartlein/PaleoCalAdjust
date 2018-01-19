program cal_adjust_PMIP3
! adjusts data in a CMIP5/PMIP3 netCDF file
! creates a new netCDF file by copying dimension variables and global attributes from input file

! this version uses an info file

use calendar_effects
use pseudo_daily_interp
use month_length_subs
use CMIP5_netCDF
use netcdf
use typesizes
use omp_lib

implicit none

! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0
! simulation age-related variables (controls orbital parameters)
integer(4)              :: begageBP             ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
integer(4)              :: endageBP             ! ending year (BP) 
integer(4)              :: agestep              ! age step size
integer(4)              :: nages                ! number of simulation ages 
integer(4), allocatable :: iageBP(:)            ! year BP 1950 (negative, e.g. 1900 CE = -50.0d0 BP)

! simulation year-related variables (controls VE_day and leap-year status)
integer(4)              :: begyrCE              ! beginning (pseudo-) year of individual model simulation
integer(4)              :: nsimyrs              ! number of years of simulation
integer(4), allocatable :: iyearCE(:)           ! yearCE simulation year (e.g. 1850CE, 850CE, etc.)

! month-length variables
integer(4), allocatable :: imonlen_0ka(:,:),imonmid_0ka(:,:)    ! integer-value month lengths -- 0ka
integer(4), allocatable :: imonbeg_0ka(:,:),imonend_0ka(:,:)    ! integer-value month beginning and ending days -- 0ka
integer(4), allocatable :: imonlen(:,:),imonmid(:,:)    ! integer-value month lengths (paleo)
integer(4), allocatable :: imonbeg(:,:),imonend(:,:)    ! integer-value month beginning and ending days (paleo)
real(8), allocatable    :: rmonlen(:,:),rmonmid(:,:)    ! real-value month lengths and mid days (paleo) 
real(8), allocatable    :: rmonbeg(:,:),rmonend(:,:)    ! real-value month beginning and ending days (paleo) 
real(8), allocatable    :: VE_day(:)                    ! vernal equinox day
integer(4), allocatable :: ndays(:)                     ! number of days in year

integer(4), allocatable :: imonlen_0ka_ts(:)            ! integer-value month lengths at present as time series
real(8), allocatable    :: rmonmid_ts(:)                ! real-value paleo month mid days as time series
real(8), allocatable    :: rmonbeg_ts(:)                ! real-value paleo month beginning days as time series
real(8), allocatable    :: rmonend_ts(:)                ! real-value paleo month ending as time series
integer(4), allocatable :: imonmid_ts(:)                ! integer-value paleo month mid days as time series
integer(4), allocatable :: imonbeg_ts(:)                ! integer-value paleo month beginning days as time series
integer(4), allocatable :: imonend_ts(:)                ! integer-value paleo month ending as time series
integer(4), allocatable :: ndays_ts(:)                  ! integer-value times series of paleo year lengths
character(256)          :: time_comment                 ! source fo new monthly time values
real(8), allocatable    :: mon_time(:)                  ! new monthly time values for daily-input files
real(8), allocatable    :: mon_time_bnds(:,:)           ! new monthly time-bounds values for daily input files

! other names
character(32)           :: calendar_type                ! calendar type

! components of file names
character(64)           :: variable                     ! variable name
character(8)            :: time_freq                    ! type of CMIP5/PMIP3 time frequency (e.g. Aclim, Amon, day, etc.)
character(8)            :: time_freq_output             ! time_freq output label
character(64)           :: model                        ! model name
character(64)           :: experiment                   ! experiment name   
character(16)           :: ensemble                     ! ensemble designator
character(8)            :: begdate, enddate             ! string beginning and ending dates of simulation   
character(32)           :: suffix                       ! file name suffix (e.g. "-clim")
character(32)           :: adj_name                     ! adjustment name (e.g. "_adj")

! data
integer(4)              :: nlon, nlat, ny, nt           ! number of longitudes, latitudes, years, obs nt = ny*nm
integer(4)              :: ndtot,ndtot_0ka,ndyr         ! total number of days
real(4), allocatable    :: var3d_in(:,:,:)              ! input data (nlon x nlat x nt)
real(4), allocatable    :: var3d_out(:,:,:)             ! output adjusted data (nlon x nlat x nt)
real(8), allocatable    :: xdh(:,:)                     ! pseudo- or actual daily data (lat x ndtot)
real(8), allocatable    :: var3d_adj(:,:)               ! adjusted data (nlat x nt)
real(8)                 :: vfill                        ! fill value

! smoothing parameters for multi-year pseudo-daily interpolation
integer(4)              :: nw_tmp=21, nsw_tmp=20        ! smoothing parameters
logical                 :: smooth=.true., restore=.true.
logical                 :: no_negatives = .false.       ! restrict pseudo-daily interpolated values to positive values?

integer(4)              :: n,m,j,k,i                    ! indices
integer(4)              :: max_threads                  ! if OpenMP enabled
integer(4)              :: iostatus                     ! IOSTAT value

! file paths and names
character(2048)         :: nc_path, ncfile_in, ncfile_out, nc_fname, infopath, debugpath, debugfile
character(64)           :: infofile
character(1)            :: csvheader                            ! info .csv file header

! if OpenMP enabled
max_threads = omp_get_max_threads()
write (*,'("OMP max_threads: ",i4)') max_threads
max_threads = max_threads - 4 ! to be able to do other things
call omp_set_num_threads(max_threads)

! path to netCDF files
nc_path = "../../data/nc_files/"

! debugging output files
debugpath="../../debug_files"
debugfile="debug_cal_adjust.dat"
open (10, file=trim(debugpath)//trim(debugfile))

! info files
infopath = "../../PaleoCalendarAdjust/info_files/"
infofile = "cal_adj_info.csv"

! open the info file, and loop over specified calendar tables
open (3,file=trim(infopath)//trim(infofile))
read (3,'(a)') csvheader

iostatus = 1
do
    ! read a line from the info file, and construct netCDF file names
    write (*,'(125("="))') 
    suffix = ""
    read (3,*,iostat=iostatus) variable, time_freq, model, experiment, ensemble, begdate, enddate, suffix, adj_name, &
        calendar_type, begageBP, endageBP, agestep, begyrCE, nsimyrs
    !write (*,'("iostatus = ",i2)') iostatus
    if (iostatus.lt.0) exit
    
    varinname = trim(variable); varoutname = varinname
    time_freq_output = trim(time_freq)
    if (trim(time_freq_output) .eq. 'day') time_freq_output = "Amon2"

    ncfile_in = trim(variable)//"_"//trim(time_freq)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
        trim(ensemble)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//".nc"
    ncfile_out = trim(variable)//"_"//trim(time_freq_output)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
        trim(ensemble)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//trim(adj_name)//".nc"
    
    write (*,'(" ncfile_in: ",a)') trim(ncfile_in)
    write (*,'("ncfile_out: ",a)') trim(ncfile_out)    
    write (*,'(a, 1x, a, 1x, 5i7)') trim(variable),trim(calendar_type), begageBP, endageBP, agestep, begyrCE, nsimyrs

    ! set no_negatives flag
    no_negatives = .false.
    select case (trim(variable))
    case ('pr','clt')
        no_negatives = .true.
    case default    
        continue
    end select
    write (*,'("no_negatives: ", l1)') no_negatives

    ! allocate month-length arrays
    nages = (endageBP - begageBP)/agestep + 1
    ny = nages * nsimyrs
    nt = ny * nm
    write (*,'("nages, nsimyrs, nt: ",3i8)') nages, nsimyrs, nt
    allocate (iageBP(ny), iyearCE(ny))
    allocate (imonlen_0ka(ny,nm),imonmid_0ka(ny,nm),imonbeg_0ka(ny,nm),imonend_0ka(ny,nm))
    allocate (imonlen(ny,nm), imonmid(ny,nm), imonbeg(ny,nm), imonend(ny,nm))
    allocate (rmonlen(ny,nm), rmonmid(ny,nm), rmonbeg(ny,nm), rmonend(ny,nm))
    allocate (VE_day(ny), ndays(ny), imonlen_0ka_ts(nt),ndays_ts(nt))
    allocate (imonmid_ts(nt), imonbeg_ts(nt),imonend_ts(nt), rmonmid_ts(nt), rmonbeg_ts(nt), rmonend_ts(nt))
    allocate (mon_time(nt), mon_time_bnds(2,nt))

    ! get month lengths
    
    ! 0 ka month lengths (used in pseudo-daily interpolation)
    if (trim(time_freq) .ne. 'day') then
        write (*,'(a)') "0 ka month lengths for pseudo-daily interpolation..."
        call get_month_lengths(calendar_type, 0, 0, agestep, nages, begyrCE, nsimyrs, & 
            iageBP, iyearCE, imonlen_0ka, imonmid_0ka, imonbeg_0ka, imonend_0ka, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, ndays)
    end if

    ! paleo month lengths
    write (*,'(a)') "Paleo month-lengths for aggregation of daily data..."
    call get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, & 
        iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, ndays)

    ! reshape month lengths into time series, and get cumulative number of days in years
    ndtot_0ka = 0; ndyr = 0
    do n=1,ny
        ndtot_0ka=ndtot_0ka + ndays(n)   
        ndyr = ndyr + ndays(n)
        !write (10,'(2i6,f12.6,2i8)') iageBP(n),iyearCE(n),VE_day(n),ndays(n),ndtot_0ka
        do m=1,nm
            i = (n-1)*nm + m
            imonlen_0ka_ts(i) = imonlen_0ka(n,m)
            rmonmid_ts(i) = rmonmid(n,m); rmonbeg_ts(i) = rmonbeg(n,m); rmonend_ts(i) = rmonend(n,m)
            imonmid_ts(i) = imonmid(n,m); imonbeg_ts(i) = imonbeg(n,m); imonend_ts(i) = imonend(n,m)
            ndays_ts(i) = ndyr
            write (10,'(3i6,i4,f12.6,i9,3i4)') & 
                n,m,i,imonlen_0ka_ts(i),rmonmid_ts(i),ndays_ts(i), imonmid_ts(i),imonbeg_ts(i),imonend_ts(i)
        end do
    end do

    ! total number of days in simulation
    ndtot = 0
    do n=1,ny
        ndtot = ndtot + ndays(n)
        !write (10,'(2i6,f12.6,2i8)') iageBP(n),iyearCE(n),VE_day(n),ndays(n),ndtot
    end do
    write (*,'("ny, nt, ndtot: ",4i8)') ny, nt, ndtot, ndtot_0ka

    ! input netCDF file
    nc_fname = trim(nc_path)//trim(ncfile_in)
    print '(" nc_fname (in) = ",a)', trim(nc_fname)
    call check( nf90_open(nc_fname, nf90_nowrite, ncid_in) )
    if (nc_print) print '("  ncid_in = ",i8)', ncid_in

    ! output netCDF file
    nc_fname = trim(nc_path)//trim(ncfile_out)
    print '(" nc_fname (out) = ",a)', trim(nc_fname)
    call check( nf90_create(nc_fname, 0, ncid_out) )
    if (nc_print) print '("  ncid_put = ",i8)', ncid_out

    ! redefine the time variables
    if (trim(time_freq) .eq. 'day') then
        time_comment = "paleo monthly mid, beginning and ending dates from original daily values"
        call new_time_day(ncid_in, ny, nm, nt, ndtot, &
            imonmid_ts, imonbeg_ts, imonend_ts, ndays_ts, mon_time, mon_time_bnds)
    else
        time_comment = "paleo monthly values from get_month_lengths()"
        call new_time_mon(calendar_type, ncid_in, ny, nm, nt, &
            rmonmid_ts, rmonbeg_ts, rmonend_ts, ndays_ts, mon_time, mon_time_bnds)
    end if

    ! create a new netCDF file, and copy dimension variables and global attributes
    call current_time(current)
    addglattname = "paleo_calendar_adjustment"
    addglatt = trim(current)//" paleo calendar adjustment by cal_adjust_PMIP3.f90"
    call copy_dims_and_glatts_redef_time(ncid_in, ncid_out, addglattname, addglatt, nt, & 
        mon_time, mon_time_bnds, time_comment, varid_out)

    ! define the output (adjusted) variable, and copy attributes
    write (*,'(a)') "Defining (new) adjusted variable..."
    addvarattname = "paleo_calendar_adjustment"
    addvaratt = "long-term mean values adjusted for the appropriate paleo calendar"
    call define_outvar(ncid_in, ncid_out, varinname, varid_out, varoutname, addvarattname, addvaratt, varid_in, nlon, nlat, nt)

    ! allocate variables
    if (trim(time_freq) .eq. 'day') then
        allocate(var3d_in(nlon,nlat,ndtot))
    else
        allocate(var3d_in(nlon,nlat,nt))
    end if
    allocate(xdh(nlat,ndtot), var3d_adj(nlat,nt), var3d_out(nlon,nlat,nt))

    ! get input data to be adjusted
    write (*,'(a)') "Reading input data..."
    if (trim(time_freq) .eq. 'day') then
        call check( nf90_get_var(ncid_in, varid_in, var3d_in) )
    else
        call check( nf90_get_var(ncid_in, varid_in, var3d_in) )
    end if
    
    ! get _FillValue
    call check( nf90_get_att(ncid_in, varid_in, '_FillValue', vfill) )
    write (*,'("_FillValue:", g14.6)') vfill
    !var3d_in(:,:,1:60)=vfill
    !var3d_in(:,:,131:150)=vfill
    !write (10,'(12g14.6)') var3d_in(40,80,:)

    ! loop over lons and lats
    write (*,'(a)') "Interpolating (if necessary) and aggregating..."
    !!$omp parallel do
    do j=1,nlon ! 
    !do j=40,40
        write(*,'(\,i5)') j; if (mod(j,25).eq.0) write (*,'(" ")')
        if (trim(time_freq) .eq.'day') xdh(:,:) = var3d_in(j,:,:)
        !$omp parallel do
        do k=1,nlat !
        !do k=80,80
            ! unless the input data is daily, read the monthly input data
            if (trim(time_freq) .ne. 'day') then
                call mon_to_day_ts(nt, imonlen_0ka_ts, dble(var3d_in(j,k,:)), dble(vfill), & 
                    no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xdh(k,:))
                call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(k,:),dble(vfill),var3d_adj(k,:))
            else
                call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(k,:),dble(vfill),var3d_adj(k,:))
            end if
        
            var3d_out(j,k,:)=sngl(var3d_adj(k,:))
        end do
        !$omp end parallel do
    end do
    !!$omp end parallel do

    where (var3d_in .eq. vfill) var3d_out = vfill
    write (10,'(a)') " "
    write (10,'(12g14.6)') var3d_out(40,80,:)


    ! write out adjusted data
    write (*,'(/a)') "Writing adjusted data..."
    call check( nf90_put_var(ncid_out, varid_out, var3d_out) )

    ! close the output file
    call check( nf90_close(ncid_out) )
    
    deallocate (iageBP, iyearCE)
    deallocate (imonlen_0ka,imonmid_0ka,imonbeg_0ka,imonend_0ka)
    deallocate (imonlen, imonmid, imonbeg, imonend)
    deallocate (rmonlen, rmonmid, rmonbeg, rmonend)
    deallocate (VE_day, ndays, imonlen_0ka_ts,ndays_ts)
    deallocate (imonmid_ts, imonbeg_ts,imonend_ts, rmonmid_ts, rmonbeg_ts, rmonend_ts)
    deallocate (mon_time, mon_time_bnds)
    deallocate (var3d_in, xdh, var3d_adj, var3d_out)
    
    write (*,'(a)') " "

end do

end program 
    

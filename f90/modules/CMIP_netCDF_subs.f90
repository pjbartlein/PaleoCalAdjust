module CMIP_netCDF_subs
! This module contains several subroutines to facilitate reading an existing netCDF file of daily or monthly data
!   and creating a new file (of adjusted data), amending and copying the attributes of the original file.
! get_var_diminfo() gets information on dimension variables in the netCDF file (including the number of dimensions).
! copy_dims_and_glatts() copies dimension, dimension variables and attributes from the old file to the new.
! new_time_day() and new_time_month() redefine the time-dimension variable to reflect the adjusted month lengths.
! define_outvar() defines the output (adjusted) variable.
! check() and current_time() are utility subroutines.

    use netcdf
    use typesizes

    implicit none

    integer, parameter          :: maxdims = 5, maxvars = 10, maxatts = 100, maxdimsize = 1440

    ! dimensions
    integer(4)                  :: ncid_in, ncid_out, ndim, nvar, natt, nglatt, unlimid, ncformat
    integer(4)                  :: dimid(maxdims), dimlen(maxdims)
    character(nf90_max_name)    :: dimname(maxdims)

    ! variables
    integer(4)                  :: varid_in, varid_out
    integer(4)                  :: xtype(maxvars), nvardims(maxvars), vardimids(maxvars,nf90_max_var_dims), nvaratts(maxvars)
    character(nf90_max_name)    :: varname(maxvars)
    integer(4)                  :: dataid(maxvars), timedimid, rmonlenid, rmonmidid, rmonbegid, rmonendid
    
    ! input variable dimensions
    integer(4)                  :: invar_ndim
    integer(4)                  :: invar_dimid(maxdims), invar_dimlen(maxdims)
    character(nf90_max_name)    :: invar_dimname(maxdims)

    ! attributes
    character(nf90_max_name)    :: attname, addattname
    character(2048)             :: addatt
    character(19)               :: current      ! current time
    integer(4)                  :: atttype
    
    ! attribute definitions
    byte                        :: att_byte       ! xtype 0, 1        nf90_byte
    character(2048)             :: att_char       ! xtype 2           nf90_char
    integer(2)                  :: att_short      ! xtype 3, 4        nf90_short
    integer(4)                  :: att_int        ! xtype 5, 6        nf90_int
    real(4)                     :: att_float      ! xtype 7, 8, 9     nf90_float
    real(8)                     :: att_double     ! xtype 10, 11      nf90_double 

    ! input and output variable names
    character(nf90_max_name)    :: varinname, varoutname
    character(nf90_max_name)    :: addvarattname
    character(2048)             :: addvaratt

    logical                     :: print_out = .true.

    contains   
    
subroutine get_var_diminfo(ncid_in, invar_name, invar_ndim, invar_dimid, invar_dimlen, invar_dimname)
! gets the dimension characteristics of a variable

    implicit none

    ! netCDF file and variable name
    integer(4), intent(in)      :: ncid_in
    character(*), intent(in)    :: invar_name

    ! dimensions
    integer(4), intent(out)     :: invar_ndim
    integer(4), intent(out)     :: invar_dimid(maxdims), invar_dimlen(maxdims)
    character(nf90_max_name), intent(out)    :: invar_dimname(maxdims)

    integer(4)                  :: invarid, i

    if (print_out) print '("In get_var_diminfo...")'
    if (print_out) print '("Variable dimensions")'
    call check( nf90_inq_varid(ncid_in, invar_name, invarid) )
    call check( nf90_inquire_variable(ncid_in, invarid, ndims=invar_ndim) )
    if (print_out) print '("invar_name, invarid, invar_ndim: ",a,1x,i5,i3)', trim(invar_name), invarid, invar_ndim

    call check( nf90_inquire_variable(ncid_in, invarid, dimids = invar_dimid(:invar_ndim)) )

    do i = 1, invar_ndim
        call check( nf90_inquire_dimension(ncid_in, invar_dimid(i), invar_dimname(i), invar_dimlen(i)) )
        if (print_out) print '("i, invar_dimid, invar_dimlen, invar_dimname: ",i3,i5,i6,1x,a)', &
        i, invar_dimid(i), invar_dimlen(i), trim(invar_dimname(i))
    end do

end subroutine get_var_diminfo

subroutine copy_dims_and_glatts(ncid_in, ncid_out, addattname, addatt, nt, time, time_bnds, comment, varid_out)
! copies dimensions and global attributes from an existing to a new netCDF file
! while replacing the old time values

    implicit none

    integer(4), intent(in)      :: ncid_in, ncid_out
    character(*), intent(in)    :: addattname, addatt
    integer(4), intent(in)      :: nt
    real(8), intent(in)         :: time(nt), time_bnds(2,nt)
    character(*), intent(in)    :: comment
    integer(4), intent(out)     :: varid_out

    real(8), allocatable        :: var1d(:), var2d(:,:), var3d(:,:,:) ! see select case() statements below for dimensions
    real(8)                     :: var0d
    integer(4)                  :: i, ii
    character(256)              :: attname

    if (print_out) print '("In copy_dims_and_glatts")'

    ! structure of input file
    if (print_out) print '("Define and copy dimensions, dimension variables, and global attributes")'
    call check( nf90_inquire(ncid_in, ndim, nvar, natt, unlimid, ncformat) )
    if (print_out) print '(" ndim, nvar, nglatt, unlimid: ",4i6)', ndim, nvar, natt, unlimid

    ! get dimensions
    if (print_out) print '("Define dimensions:")'
    write (*,*) nt
    do i = 1,ndim
        dimid(i) = i
        call check( nf90_inquire_dimension(ncid_in, dimid(i), dimname(i), dimlen(i)) )
        if (print_out) print '("   in:  dimid, dimlen, dimname = " ,2i7,1x,a)', i, dimlen(i), trim(dimname(i))

        if (i .eq. unlimid) then
            call check( nf90_def_dim(ncid_out, dimname(i), nf90_unlimited, dimid(i)) )
        else
            call check( nf90_def_dim(ncid_out, dimname(i), dimlen(i), dimid(i)) )
        end if

        ! if dimension is time, replace length with new length (nt), and save timedimid
        if (dimname(i) .eq. 'time') then
            dimlen(i) = nt
            timedimid = dimid(i)
        end if

    end do

    ! define dimension variables
    if (print_out) print '("Define dimension variables:")'
    varid_out = 0
    do i = 1,nvar

        varid_in = i
        call check( nf90_inquire_variable(ncid_in, varid_in, varname(i), xtype(i), ndims=nvardims(i), natts=nvaratts(i)) )
        if (print_out) print '(" i, xtype, nvardims, nvaratts = ", 4i6, 1x, a)', i,xtype(i),nvardims(i),nvaratts(i),trim(varname(i))
        call check( nf90_inquire_variable(ncid_in, varid_in, dimids=vardimids(i,:nvardims(i))) )
        if (print_out) print '("    vardimids = ", 6i6)', vardimids(i,:nvardims(i))

        ! define dimension variables only
        select case (varname(i))
        case ('lon', 'lon_bnds', 'lat', 'lat_bnds', 'time', 'time_bnds', 'climatology_bnds', 'height', &
                'plev', 'j', 'i', 'vertices', 'lon_vertices', 'lat_vertices', 'lev', 'lev_bnds')

            ! define variable
            varid_out=varid_out+1
            call check( nf90_def_var(ncid_out, varname(i), xtype(i), vardimids(i,:nvardims(i)), dataid(i)) )
            if (print_out) print '("     defining varid_out : ",i4,1x,a)', varid_out, trim(varname(i))

            ! copy attributes
            do ii=1,nvaratts(i)
                call check( nf90_inq_attname(ncid_in, varid_in, ii, attname) )
                call check( nf90_copy_att(ncid_in, varid_in, attname, ncid_out, varid_out) )
                if (print_out) print '("  ii, attname: ",i4,1x,a)', ii, trim(attname)
            end do

            ! dimension variable is time, add comment
            if (varname(i) .eq. 'time') then
                call check ( nf90_put_att(ncid_out, varid_out, 'comment', comment) )
                if (print_out) print '("  ii, attname: ",i4,1x,a)', ii, "comment"
            end if

        case default
            continue
        end select
    end do

    ! copy global attributes
    if (print_out) print '("Copy global attributes")'
    do i = 1, natt
        call check ( nf90_inq_attname(ncid_in, nf90_global, i, attname) )
        call check ( nf90_inquire_attribute(ncid_in, nf90_global, trim(attname), atttype) )

        if (print_out) print '(" global attribute: ", i4, 1x, i12, 1x, a)', i, atttype, trim(attname)

        select case(atttype)
        case (0, 1)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_byte) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_byte) )
        case (2)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_char) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_char) )
        case (3, 4)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_short) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_short) )
        case (5, 6)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_int) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_int) )
        case (7, 8, 9)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_float) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_float) )
        case (10, 11)
            call check ( nf90_get_att(ncid_in, nf90_global, attname, att_double) )
            call check ( nf90_put_att(ncid_out, nf90_global, attname, att_double) )
        case default
            stop "copying global attributes"
        end select
    end do

    ! add new global attribute
    call check ( nf90_put_att(ncid_out, nf90_global, addattname, addatt) )

    ! end definition of new netCDF file (temporarily)
    call check( nf90_enddef(ncid_out) )

    ! copy dimension variable values, replacing time and time_bnds
    if (print_out) print '("Copy dimension variables")'
    varid_out = 0
    do i = 1,nvar

        varid_in = i
        select case (varname(i))
        case ('lon', 'lon_bnds', 'lat', 'lat_bnds', 'time', 'time_bnds', 'climatology_bnds', 'height', &
                'plev', 'j', 'i', 'vertices', 'lon_vertices', 'lat_vertices', 'lev', 'lev_bnds')

            varid_out=varid_out + 1
            if (print_out) print '(" varname: ",a)', trim(varname(i))
            ! how many dimensions?
            select case(nvardims(i))
            case (0)
                if (print_out) print '("  i: nvardims(i): ", 2i6, 1x, a)', &
                    i,nvardims(i),trim(varname(i))
                call check( nf90_get_var(ncid_in, varid_in, var0d) )

                ! if variable is time, replace the existing values with new ones
                if (varname(i) .eq. 'time') var1d=time
                call check( nf90_put_var(ncid_out, varid_out, var0d) )
            case (1)
                if (print_out) print '("  i, nvardims(i), vardimids, dimlen: ", 4i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), dimlen(vardimids(i,1)),trim(varname(i))
                allocate(var1d(dimlen(vardimids(i,1))))
                call check( nf90_get_var(ncid_in, varid_in, var1d) )

                ! if variable is time, replace the existing values with new ones
                if (varname(i) .eq. 'time') var1d=time
                call check( nf90_put_var(ncid_out, varid_out, var1d) )
                deallocate(var1d)
            case (2)
                if (print_out) print '("  i, nvardims(i), vardimids, dimlen: ", 6i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), vardimids(i,2), dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), trim(varname(i))
                allocate(var2d(dimlen(vardimids(i,1)), dimlen(vardimids(i,2))))
                call check( nf90_get_var(ncid_in, varid_in, var2d) )

                ! if variable is time_bnds or climatology_bnds, replace the existing values with new ones
                if (varname(i) .eq. 'time_bnds') var2d=time_bnds
                if (varname(i) .eq. 'climatology_bnds') var2d=time_bnds
                call check( nf90_put_var(ncid_out, varid_out, var2d) )
                deallocate(var2d)
            case (3)
                if (print_out) print '("  i, nvardims(i), vardimids, dimlen: ", 8i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), vardimids(i,2), vardimids(i,3), &
                    dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), dimlen(vardimids(i,3)), trim(varname(i))
                allocate(var3d(dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), dimlen(vardimids(i,3))))
                call check( nf90_get_var(ncid_in, varid_in, var3d) )

                ! if variable is time_bnds or climatology_bnds, replace the existing values with new ones
                if (varname(i) .eq. 'time_bnds') var2d=time_bnds
                if (varname(i) .eq. 'climatology_bnds') var2d=time_bnds
                call check( nf90_put_var(ncid_out, varid_out, var3d) )
                deallocate(var3d)

            case default
                continue
            end select

        case default
            if (print_out) print '(a,1x,"was not a recognizable dimension variable")', trim(varname(i))
        end select
    end do

    varid_out = varid_out + 1   ! id of next variable to be written

end subroutine copy_dims_and_glatts

subroutine new_time_day(ncid_in, ny, nm, nt, ndtot, &
    imonmid_ts, imonbeg_ts, imonend_ts, ndays_ts, time, time_bnds)

    implicit none

    integer(4), intent(in)      :: ncid_in
    integer(4), intent(in)      :: ny, nm, nt           ! number of years, months and total number of months
    integer(4), intent(in)      :: ndtot                ! total number of days
    integer(4), intent(in)      :: imonmid_ts(nt)       ! month mid-days as time series
    integer(4), intent(in)      :: imonbeg_ts(nt)       ! month beginning days as time series
    integer(4), intent(in)      :: imonend_ts(nt)       ! month ending days as time series
    integer(4), intent(in)      :: ndays_ts(nt)         ! number of days in year

    real(8), intent(out)        :: time(nt), time_bnds(2,nt)                ! new time variables

    integer(4)                  :: timeid               ! variable id
    real(8)                     :: day_time(ndtot)      ! (old) time value for each day
    integer(4)                  :: ndyr                 ! number of days in previous years

    integer(4)                  :: imid, ibeg, iend     ! indices
    integer(4)                  :: n

    write (*,'("new time values: ny, nm, nt, ndtot: ",4i8)') ny,nm,nt,ndtot

    ! get the existing daily time values
    call check ( nf90_inq_varid(ncid_in, 'time', timeid) )
    call check ( nf90_get_var(ncid_in, timeid, day_time) )

    ! new time variables -- copy appropriate existing daily time values
    imid = imonmid_ts(1); ibeg = imonbeg_ts(1); iend = imonend_ts(1)
    if (ibeg .lt. 1) ibeg = 1
    time(1) = day_time(imid)
    time_bnds(1,1) = day_time(ibeg)
    time_bnds(2,1) = day_time(iend)
    do n = 2,nt
        ndyr = ndays_ts(n) - ndays_ts(1)
        imid = imonmid_ts(n); ibeg = imonbeg_ts(n); iend = imonend_ts(n)
        time(n) = day_time(imid + ndyr)
        time_bnds(1,n) = time_bnds(2,n-1)
        time_bnds(2,n) = day_time(iend + ndyr)
    end do

end subroutine new_time_day

subroutine new_time_month(calendar_type, ncid_in, ny, nm, nt, &
     rmonmid_ts, rmonbeg_ts, rmonend_ts, ndays_ts, time, time_bnds)

    ! redefines monthly time variables

    implicit none

    character(*), intent(in)    :: calendar_type
    integer(4), intent(in)      :: ncid_in
    integer(4), intent(in)      :: ny, nm, nt           ! number of years, months and total number of months
    real(8), intent(in)         :: rmonmid_ts(nt)       ! month mid-days as time series
    real(8), intent(in)         :: rmonbeg_ts(nt)       ! month beginning days as time series
    real(8), intent(in)         :: rmonend_ts(nt)       ! month ending days as time series
    integer(4), intent(in)      :: ndays_ts(nt)         ! number of days in year

    real(8), intent(out)        :: time(nt), time_bnds(2,nt)                ! new time variables

    integer(4)                  :: timeid               ! variable id
    real(8)                     :: mon_time(nt)         ! (old) time value for each month
    real(8)                     :: ref_time             ! reference time (e.g. Jan "0" of the first year)
    integer(4)                  :: ndyr                 ! number of days in previous years

    integer(4)                  :: n

    write (*,'("new_time: ny, nm, nt: ",3i8)') ny,nm,nt

    ! get the existing monthly time values
    call check ( nf90_inq_varid(ncid_in, 'time', timeid) )
    call check ( nf90_get_var(ncid_in, timeid, mon_time) )

    ! new time variables -- calculate appropriate monthly values
    if (trim(calendar_type) .eq. '360_day') then
        ref_time = mon_time(1) - 15.0
    else
        ref_time = mon_time(1) - 15.5
    end if

    time(1) = rmonmid_ts(1) + ref_time
    time_bnds(1,1) = rmonbeg_ts(1) + ref_time
    time_bnds(2,1) = rmonend_ts(1) + ref_time

    do n = 2,nt
        ndyr = ndays_ts(n) - ndays_ts(1)
        time(n) = rmonmid_ts(n) + dble(ndyr) + ref_time
        time_bnds(1,n) = time_bnds(2,n-1)
        time_bnds(2,n) = rmonend_ts(n) + dble(ndyr) + ref_time
    end do

end subroutine new_time_month


subroutine define_outvar(ncid_in, ncid_out, varinname, varid_out, varoutname, addvarattname, addvaratt, varid_in)

    implicit none

    integer(4), intent(in)      :: ncid_in, ncid_out
    character(*), intent(in)    :: varinname, varoutname
    integer(4), intent(in)      :: varid_out
    character(*), intent(in)    :: addvarattname, addvaratt
    integer(4), intent(out)     :: varid_in

    integer(4)                  :: i, ii

    character(256)              :: attname

    ! find input variable number
    do i=1,nvar
        if (trim(varname(i)) .eq. trim(varinname)) exit
    end do
    varid_in = i
    if (print_out) print '(" i, varid, nvardims(i), varname(i), varinname: ",3i4,1x,a,1x,a)', &
        i, varid_in, nvardims(i), trim(varname(i)), trim(varinname)

    ! enter define mode again
    call check( nf90_redef(ncid_out) )

    ! define variable
    call check( nf90_def_var(ncid_out, varoutname, xtype(i), vardimids(i,:nvardims(i)), dataid(i)) )
    if (print_out) print '(" defining output variable : ",i4,1x,a)', varid_out, trim(varname(i))

    ! copy attributes
    do ii=1,nvaratts(i)
        call check( nf90_inq_attname(ncid_in, varid_in, ii, attname) )
        call check( nf90_copy_att(ncid_in, varid_in, attname, ncid_out, varid_out) )
        if (print_out) print '("  ii, attname: ",i4,1x,a)', ii, trim(attname)
    end do

    ! add new variable attribute
    call check ( nf90_put_att(ncid_out, varid_out, addvarattname, addvaratt) )

    ! define month-length time series variables

    if (print_out) print '("Defining month-length variables...")'

    call check ( nf90_def_var(ncid_out, 'rmonlen', nf90_double, timedimid, rmonlenid) )
    call check ( nf90_put_att(ncid_out, rmonlenid, 'long_name', 'month length') )
    call check ( nf90_put_att(ncid_out, rmonlenid, 'units', 'd') )
    call check ( nf90_put_att(ncid_out, rmonlenid, 'comment', 'real-value month length (days)') )

    call check ( nf90_def_var(ncid_out, 'rmonmid', nf90_double, timedimid, rmonmidid) )
    call check ( nf90_put_att(ncid_out, rmonmidid, 'long_name', 'mid-month day in year') )
    call check ( nf90_put_att(ncid_out, rmonmidid, 'units', 'd') )
    call check ( nf90_put_att(ncid_out, rmonmidid, 'comment', 'real-value mid-month day in year') )

    call check ( nf90_def_var(ncid_out, 'rmonbeg', nf90_double, timedimid, rmonbegid) )
    call check ( nf90_put_att(ncid_out, rmonbegid, 'long_name', 'beginning day in year') )
    call check ( nf90_put_att(ncid_out, rmonbegid, 'units', 'd') )
    call check ( nf90_put_att(ncid_out, rmonbegid, 'comment', 'real-value month beginning day in year') )

    call check ( nf90_def_var(ncid_out, 'rmonend', nf90_double, timedimid, rmonendid) )
    call check ( nf90_put_att(ncid_out, rmonendid, 'long_name', 'ending day in year') )
    call check ( nf90_put_att(ncid_out, rmonendid, 'units', 'd') )
    call check ( nf90_put_att(ncid_out, rmonendid,'comment', 'real-value month ending day in year') )

    ! end definition of new netCDF file
    call check( nf90_enddef(ncid_out) )
    if (print_out) print '("Done defining variables")'

end subroutine define_outvar

subroutine check(status)
! netCDF error message handler

    use netcdf
    use typesizes

    implicit none

    integer, intent(in) :: status
    if (status.ne.nf90_noerr) then
        print *, status,trim(nf90_strerror(status))
        stop 'stopped on error'
    end if

end subroutine check

subroutine current_time(current)
! gets the current time

    character(19), intent(out) :: current
    character(8)        :: cdate
    character(10)       :: ctime

    call date_and_time(cdate,ctime)

    current = cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)//" "// &
        ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)

end subroutine current_time

end module CMIP_netCDF_subs

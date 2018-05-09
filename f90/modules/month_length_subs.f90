module month_length_subs

    implicit none

    integer(4), parameter   :: nm = 12              ! number of months in the year
    integer(4), parameter   :: nd_360 = 360         ! number of days in a 360-day year
    integer(4), parameter   :: daysinmonth360 = 30  ! number of days in a month of a 360-day year
    integer(4), parameter   :: nd_365 = 365         ! number of days in a 365-day "noleaps" year
    integer(4), parameter   :: nd_366 = 366         ! number of days in a 366-day leap year

    ! other calendar-related variables
    real(8)     :: veqday_360 = 80.0d0              ! fixed vernal equinox day, 360-day year
    real(8)     :: veqday_365 = 80.5d0              ! fixed vernal equinox day, 365-day year
    real(8)     :: veqday_366 = 81.5d0              ! fixed vernal equinox day, 366-day year
    real(8)     :: tropical_year = 365.24219876d0   ! length of a tropical year (days)
    real(8)     :: progreg_year = 365.2425d0        ! length of a Gregorian year (days)

    integer(4)  :: nd_progreg                       ! number of days in a 365 or 366-day year proleptic_gregorian calendar
    real(8)     :: veqday_progreg                   ! vernal equinox day in a 365 or 366-day year proleptic_gregorian calendar
    real(8)     :: midMarch                         ! mid-March day
    real(8)     :: perihelion                       ! perihelion day (in a proleptic Gregorian calendar)
    real(8)     :: aphelion                         ! aphelion day (in a proleptic Gregorian calendar)

    ! month-length definitions
    real(8)     :: present_mon_360(nm) = 30.0d0     ! present-day month lengths in 360-day year
    real(8)     :: present_mon_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
        (/ 31.0d0, 28.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
        (/ 31.0d0, 29.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_365_trop(nm) = &     ! present-day month lengths in a tropical year (note Feb.)
        (/ 31.0d0, 28.24219876d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_365_progreg(nm) = &  ! present-day month lengths in a Gregorian year (note Feb.)
        (/ 31.0d0, 28.2425d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)

contains

subroutine get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, &
    iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, ndays)

    use GISS_orbital_subs
    use GISS_srevents_subs

    implicit none

    ! calendar type
    character(23), intent(in)   :: calendar_type

    ! simulation age-related variables (controls orbital elements)
    integer(4), intent(in)  :: begageBP                     ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
    integer(4), intent(in)  :: endageBP                     ! ending year (BP)
    integer(4), intent(in)  :: agestep                      ! age step size
    integer(4), intent(in)  :: nages                        ! number of simulation ages

    ! individual model simulation year-related variables (controls VE_day and leap-year status)
    integer(4), intent(in)  :: begyrCE                      ! beginning (pseudo-) year of individual model simulation
    integer(4), intent(in)  :: nsimyrs                      ! number of years of simulation

    ! (output) month-length variables
    integer(4), intent(out) :: iageBP(nages*nsimyrs)        ! year BP 1950 (negative, e.g. 1900 CE = -50.0d0 BP)
    integer(4), intent(out) :: iyearCE(nages*nsimyrs)       ! yearCE simulation year
    integer(4), intent(out) :: imonlen(nages*nsimyrs,nm)    ! integer-value month lengths
    integer(4), intent(out) :: imonmid(nages*nsimyrs,nm)    ! integer-value mid-month days
    integer(4), intent(out) :: imonbeg(nages*nsimyrs,nm)    ! integer-value beginning days
    integer(4), intent(out) :: imonend(nages*nsimyrs,nm)    ! integer-value ending days
    real(8), intent(out)    :: rmonlen(nages*nsimyrs,nm)    ! real-value month lengths
    real(8), intent(out)    :: rmonmid(nages*nsimyrs,nm)    ! real-value mid-month days
    real(8), intent(out)    :: rmonbeg(nages*nsimyrs,nm)    ! real-value month beginning day
    real(8), intent(out)    :: rmonend(nages*nsimyrs,nm)    ! real-value month ending days
    real(8), intent(out)    :: VE_day(nages*nsimyrs)        ! real-value vernal equinox day
    integer(4), intent(out) :: ndays(nages*nsimyrs)         ! integer number of days in year

    ! subroutine GISS_orbpars() and GISS_srevents() input and output arguments
    character(2)            :: year_type = 'BP'             ! AD (AD/BC), CE (CE/BCE), BP (bp 1950)
    real(8)                 :: AgeBP                        ! age (BP 1950) (input)
    real(8)                 :: eccen                        ! eccentricity of orbital ellipse
    real(8)                 :: obliq_deg                    ! obliquity (degrees)
    real(8)                 :: perih_deg                    ! longitude of perihelion (degrees)
    real(8)                 :: precc                        ! climatological precession parameter = eccen * sin(omegvp)
    real(8)                 :: veqday                       ! real day of vernal equinox

    ! kg_monlen() subroutine arguments
    real(8)                 :: yrlen                        ! real number of days in year (year length)
    integer(4)              :: ndyr                         ! integer number of days in year

    ! present-day month-length variables
    real(8)                 :: rmonlen_0ka(nm)              ! real calculated month lengths at 0 ka
    real(8)                 :: ryeartot_0ka                 ! real calculated total number of days at 0 ka
    real(8)                 :: rmonlen_0ka_leap(nm)         ! real calculated month lengths at 0 ka in a leap year
    real(8)                 :: ryeartot_0ka_leap            ! real calculated total number of days at 0 ka in a leap year
    real(8)                 :: present_monlen(nm)           ! "present day" month lengths

    ! arrays for calculating various month-length statistics
    real(8)                 :: rmonlen_rel(nm)              ! difference between real month lengths and present
    real(8)                 :: rmonlen_ratio(nm)            ! ratio of real month lengths and present
    real(8)                 :: ryeartot(nages*nsimyrs)      ! real total number of days in year
    integer(4)              :: iyeartot(nages*nsimyrs)      ! integer total number of days in year

    ! indices
    integer(4)              :: n        ! simulation-age index
    integer(4)              :: i        ! simulation-year index
    integer(4)              :: ii       ! age and year index
    integer(4)              :: m        ! month index

    logical                 :: debug_write = .true.  ! write additional debugging info
    logical                 :: debug_kgmonlen = .true.

    ! check for supported calendar types
    select case (trim(calendar_type))
    case ('360_day','noleap','365_day','365.2425','proleptic_gregorian','progreg','gregorian','standard')
        continue
    case default
        stop "Calendar type not supported"
    end select

    ! debugging output
    if (debug_write) then
        open (22,file="/Projects/Calendar/data/work01/"//trim(calendar_type)//"_monlen_debug.dat")
        open (11,file="/Projects/Calendar/data/work01/"//trim(calendar_type)//"_cal_rmonlen_raw.dat")
        open (12,file="/Projects/Calendar/data/work01/"//trim(calendar_type)//"_cal_rmonlen_rel.dat")
        open (13,file="/Projects/Calendar/data/work01/"//trim(calendar_type)//"_cal_rmonlen_ratio.dat")
    end if
    if (debug_kgmonlen) open (23,file="/Projects/Calendar/data/work01/"//trim(calendar_type)//"_debug_kg.dat")

    ! generate target years -- experiment ages (in yrs BP) and simulation years (in yrs CE)
    write (*,'(a)') "Calculating month lengths..."
    write (*,'("begageBP, endageBP, agestep, nages, begyrCE, nsimyrs: ",6i7)') &
        begageBP, endageBP, agestep, nages, begyrCE, nsimyrs
    if (debug_write) write (22,'("begageBP, endageBP, agestep, nages, begyrCE, nsimyrs: ",6i7)') &
            begageBP, endageBP, agestep, nages, begyrCE, nsimyrs

    ii=0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii+1
            iageBP(ii) = begageBP + (n - 1) * agestep
            iyearCE(ii) = begyrCE + (i - 1)
            if (debug_write) write (22,'("n,i,ii,iageBP,iyearCE ", 5i8)') n,i,ii,iageBP(ii),iyearCE(ii)
        end do
    end do

    ! initialize arrays
    imonlen = 0; imonmid = 0; imonbeg = 0; imonend = 0
    rmonlen = 0.0d0; rmonmid = 0.0d0; rmonbeg = 0.0d0; rmonend = 0.0d0
    VE_day = 0.0d0; ndays = 0
    ryeartot = 0.0d0; iyeartot = 0

    ! ===================================================================================================================

    ! get 0 ka (1950CE) calculated month lengths, which will be used to restore all other calculated month lengths
    ! to nominal present-day values

    ! orbital elements for 0 ka
    write (*,'(a)') "0 ka orbital elements..."
    ageBP = 0.0d0
    call GISS_orbpars('BP', ageBP, eccen, obliq_deg, perih_deg, precc)
    if (debug_write) write (22,'("ageBP, eccen, obliq_deg, perih_deg, precc: ",f10.1,4f17.12)') &
        ageBP, eccen, obliq_deg, perih_deg, precc

    ! 0 ka calculated month lengths, following Kutzbach and Gallimore (1988), also set subroutine arguments
    write (*,'(a)') "0 ka month lengths..."
    select case (trim(calendar_type))
    case ('360_day')
        yrlen = dble(nd_360); ndyr = nd_360; veqday = veqday_360; present_monlen = present_mon_360
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, ryeartot_0ka)
    case ('noleap', '365_day')
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, ryeartot_0ka)
    case ('366_day')
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, ryeartot_0ka)
    case ('365.2425')
        yrlen = dble(progreg_year); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_365_progreg
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, ryeartot_0ka)
    case ('proleptic_gregorian','progreg','gregorian','standard')
        ! 365-day year
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, ryeartot_0ka)
        ! 366-day year
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka_leap, ryeartot_0ka_leap)
    case default
        stop "calendar type not supported"
    end select

    if (debug_write) write (22,'("rmonlen_0ka: ",21x,13f12.6)') rmonlen_0ka(1:nm), ryeartot_0ka

    ! loop over simulation ages and years

    ii = 0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii + 1
            if (mod(n,1000) .eq. 0) write(*,'(\,i8)') n

            ! orbital elements for simulation age (e.g. 6 ka)
            call GISS_orbpars('BP', dble(iageBP(n)), eccen, obliq_deg, perih_deg, precc)
            !if (debug_write) write (22,'("ageBP, eccen, obliq_deg, perih_deg, precc: ",f10.1,4f17.12)') &
            !    ageBP, eccen, obliq_deg, perih_deg, precc

            ! proleptic_gregorian-like calendars
            select case (trim(calendar_type))
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! check for leap year
                nd_progreg= yearlen_CE(iyearCE(ii))
                if (debug_write) write (22,'("ii, iyearCE, nd_progreg: ",3i6,i4)') ii, iageBP(ii), iyearCE(ii), nd_progreg
                if (nd_progreg .eq. 365) then ! (non leap year)
                    yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
                    ! get real-valued month lengths
                    call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen(ii,:), ryeartot(ii))
                    ! adjust values such so that 0 ka (1950 CE) will have nominal present-day month lengths
                    call adjust_to_reference(rmonlen(ii,:), rmonlen_0ka, present_monlen)
                else ! ndprogreg = 366 (leap year)
                    yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
                    call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen(ii,:), ryeartot(ii))
                    ! adjust values such so that 0 ka (1950 CE) will have nominal present-day month lengths
                    call adjust_to_reference(rmonlen(ii,:), rmonlen_0ka_leap, present_monlen)
                end if
            case default ! other calendars ('360_day','noleap','365_day','365.2425')
                ! get real-value month lengths
                call kg_monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen(ii,:), ryeartot(ii))
                ! adjust values such so that 0 ka (1950 CE) will have nominal present-day month lengths
                call adjust_to_reference(rmonlen(ii,:), rmonlen_0ka, present_monlen)
            end select

            if (debug_write) write (22,'("ageBP, rmonlength:        ",i8,13f12.6)') iageBP(ii),rmonlen(ii,1:nm),ryeartot(ii)
            if (debug_write) write (11,'(2i8,13f12.6)') iageBP(ii),iyearCE(ii),rmonlen(ii,1:nm),ryeartot(ii) ! raw month lengths

            ! require the sum of month lengths each year to equal the year length
            call adjust_to_yeartot(rmonlen(ii,:), yrlen, ryeartot(ii))

            ! integer month lengths
            call integer_monlen(rmonlen(ii,:), ndyr, imonlen(ii,:), iyeartot(ii))

            ! various month-length statistics
            do m=1,nm
                rmonlen_rel(m) = rmonlen(ii,m)-present_monlen(m)
                rmonlen_ratio(m) = rmonlen(ii,m)/present_monlen(m)
            end do

            ! get mid-March day
            select case (trim(calendar_type))
            case ('360_day')
                midMarch = veqday_360 - (5.0d0/30.0d0) * rmonlen(ii,3)
            case ('noleap', '365_day')
                midMarch = veqday_365 - (6.0d0/31.0d0) * rmonlen(ii,3)
            case ('366_day')
                midMarch = veqday_366 - (6.0d0/31.0d0) * rmonlen(ii,3)
            case ('365.2425')
                midMarch = veqday_365 - (6.0d0/31.0d0) * rmonlen(ii,3)
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! get vernal equinox day for model simulation year (not age)
                year_type = 'CE'
                call GISS_srevents(year_type, iyearCE(ii), progreg_year, veqday, perihelion, aphelion, ndyr)
                midMarch = veqday - (6.0d0/31.0d0) * rmonlen(ii,3)
            case default
                stop "paleo calendar type not defined"
            end select

            ! real-value mid-month, beginning and ending days
            call rmon_midbegend(rmonlen(ii,:), midMarch, rmonmid(ii,:), rmonbeg(ii,:), rmonend(ii,:))

            ! integer mid-month, beginning and ending days
            call imon_midbegend(imonlen(ii,:), int(midMarch), imonmid(ii,:), imonbeg(ii,:), imonend(ii,:))

            ! save VE day and number of days in year
            VE_day(ii) = veqday; ndays(ii) = ndyr

            if (debug_write) then
                write (12,'(2i8,13f12.8)') iageBP(ii),iyearCE(ii),rmonlen_rel
                write (13,'(2i8,13f12.8)') iageBP(ii),iyearCE(ii),rmonlen_ratio
            end if

        end do
    end do

    if (debug_write) close(22); close(11); close(12); close(13)
    if (debug_kgmonlen) close(23)

end subroutine get_month_lengths

subroutine kg_monlen_360(eccen, perih, rmonlen, ryeartot)

    implicit none

    real(8), intent(in)     :: eccen, perih         ! orbital parameters
    real(8), intent(out)    :: rmonlen(nm)          ! real 360-day month length
    real(8), intent(out)    :: ryeartot             ! real total number of days in year

    real(8)     :: verneq_angle(0:nd_360)           ! (angular) days-since-vernal equinox

    real(8)     :: phi                              ! phase angle
    real(8)     :: phip                             ! phase-angle phip (such that sin(perih-phip) = -1)
    real(8)     :: diff, mindiff                    ! phase-angle difference and minimum difference
    real(4)     :: phi_inc                          ! phase-ange estimation increment
    integer(4)  :: nphi=360001                      ! number of phase-angle estimation steps
    integer(4)  :: imin                             ! phase-angle minimum index

    real(8)     :: t(0:nd_360)                      ! traverse time since vernal equinox (K&G Eqn A2)
    real(8)     :: day_length(0:nd_360)             ! relative length of each day

    integer(4)  :: day(0:nd_360)                    ! integer day index
    integer(4)  :: i,ii,m                           ! indices

    verneq_angle=0.0d0; rmonlen=0.0d0

    ! set day numbers (day(0) is the day before Jan 1) and angular difference from vernal equinox
    do i=0,nd_360
        day(i) = i
        verneq_angle(i) = (-1.0d0*veqday_360)+dble(i)
        write (23,'("i,day,verneq_angle: ",i4,i4,f10.3)') i,day(i),verneq_angle(i)
    end do

    ! find phip (such that sin(perih-phip) = -1)
    phi_inc = 360.0d0 / dble(nphi - 1)
    mindiff = 99999.0d0; phip = 0.0d0; imin = 99999
    do i=0,nphi+1
        phi = dble(i) * phi_inc
        diff = -1.0d0 - dsin(radians(perih-phi))
        if (dabs(diff).le.mindiff) then
            mindiff = dabs(diff)
            phip = phi
            imin = i
        end if
        ! write (23,'("i,phi,perih,diff,phip,mindiff,imin: ",i6,5f11.6,i6)') i,phi,perih,diff,phip,mindiff,imin
    end do
    write (23,'("perih,phip,mindiff,imin: ",3g14.6,i6)') perih,phip,mindiff,imin

    ! calculate t (traverse time since vernal equinox) (K&G Eqn A2)
    do i=0,nd_360
        t(i) = verneq_angle(i)-2.0d0*eccen*degrees(dcos(radians(verneq_angle(i)-phip)) - dcos(radians(phip)))
    end do

    ! relative length of each day
    day_length(0) = t(0) - (t(nd_360-1) - 360.0d0)
    ryeartot = 0.0d0
    write (23, '("i, t(0), day_length(i), ryeartot:   0",3f12.6)') t(0), day_length(0), ryeartot
    do i=1,nd_360
        day_length(i) = t(i) - t(i-1)
        ryeartot = ryeartot + day_length(i)
        write (23, '("i, t(i), day_length(i), ryeartot: ",i3,3f12.6)') i, t(i), day_length(i), ryeartot
    end do

    ! length of each month -- total daylength in 1/12 = daysinmonth360 / 360 proportion of year
    i=0; ryeartot = 0.0d0
    do m=1,nm
        do ii=1,daysinmonth360
            i=i+1
            rmonlen(m) = rmonlen(m) + day_length(i)
            write (23,'("m, i, ii, day_length(i), rmonlen(m):", 3i4,2f12.6)') m,i,ii,day_length(i),rmonlen(m)
        end do
        ryeartot = ryeartot + rmonlen(m)
    end do
    write (23,'("rmonlen: ",13f12.6)') rmonlen(1:nm), ryeartot

end subroutine kg_monlen_360

subroutine kg_monlen(yrlen, ndyr, veqday, ipresent_monlen, eccen, perih, rmonlen, ryeartot)

    implicit none

    real(8), intent(in)     :: yrlen                ! total length of year (days)
    integer(4), intent(in)  :: ndyr                 ! number of days in year
    real(8), intent(in)     :: veqday               ! vernal equinox day
    integer(4), intent(in)  :: ipresent_monlen(nm)  ! number of days in month at present
    real(8), intent(in)     :: eccen, perih         ! orbital parameters

    real(8), intent(out)    :: rmonlen(nm)          ! real month lengths
    real(8), intent(out)    :: ryeartot             ! real total number of days in year (as check)

    real(8)     :: verneq_angle(0:ndyr)             ! (angular) days-since-vernal equinox

    real(8)     :: phi                              ! phase angle (degrees)
    real(8)     :: phip                             ! phase-angle phip (such that sin(perih-phip) = -1)
    real(8)     :: diff, mindiff                    ! phase-angle difference and minimum difference
    real(4)     :: phi_inc                          ! phase-ange estimation increment
    integer(4)  :: nphi=360001                      ! number of phase-angle estimation steps
    integer(4)  :: imin                             ! phase-angle minimum index

    real(8)     :: pi                               ! pi
    real(8)     :: t(0:ndyr)                        ! traverse time since vernal equinox (K&G Eqn A2)
    real(8)     :: day_length(0:ndyr)               ! relative length of each day

    integer(4)  :: day(0:ndyr)                      ! integer day index
    integer(4)  :: i,ii,m                           ! indices

    pi=4.0d0*datan(1.0d0)

    verneq_angle=0.0d0; rmonlen=0.0d0

    ! set day numbers (day(0) is the day before Jan 1) and angular differences (in degrees) from vernal equinox
    do i=0,ndyr
        day(i) = i
        verneq_angle(i) = ((-1.0d0*veqday)+dble(i))*(360.0d0/dble(ndyr))
        write (23,'("i,day,verneq_angle: ",i4,i4,f12.6)') i,day(i),verneq_angle(i)
    end do

    ! find phip (in degrees) such that sin((2 * pi/yrlen)*(perih-phip)) = -1
    phi_inc = 360.0d0 / dble(nphi - 1)
    mindiff = 99999.0d0; phip = 0.0d0; imin = 99999
    do i=0,nphi+1
        phi = dble(i) * phi_inc
        diff = -1.0d0 - dsin(radians(perih-phi))
        if (dabs(diff).le.mindiff) then
            mindiff = dabs(diff)
            phip = phi
            imin = i
        end if
        ! write (23,'("i,phi,perih,diff,phip,mindiff,imin: ",i6,5f11.6,i6)') i,phi,perih,diff,phip,mindiff,imin
    end do
    write (23,'("perih,phip,mindiff,imin: ",3g14.6,i6)') perih,phip,mindiff,imin

    ! calculate t (traverse time) from vernal equinox for each day (K&G Eqn A2)
    do i=0,ndyr
        t(i) = (verneq_angle(i)-2.0d0*eccen*degrees(dcos(radians(verneq_angle(i)-phip)) &
            - dcos(radians(phip)))) * (dble(yrlen)/360.0d0)
    end do

    ! relative length of each day
    day_length(0) = (t(ndyr) - t(ndyr-1))
    ryeartot = 0.0d0
    write (23, '("i, t(0), day_length(i), ryeartot:   0",3f12.6)') t(0), day_length(0), ryeartot
    do i=1,ndyr
        day_length(i) = (t(i) - t(i-1))
        ryeartot = ryeartot + day_length(i)
        write (23, '("i, t(i), day_length(i), ryeartot: ",i3,3f12.6)') i, t(i), day_length(i), ryeartot
    end do

    ! length of each month -- total relative day length in ipresent_monlen(m)/ndyr proportion of year
    i=0; ryeartot = 0.0d0
    do m=1,nm
        do ii=1,ipresent_monlen(m)
            i=i+1
            rmonlen(m) = rmonlen(m) + day_length(i)
            !write (23,'("m, i, ii, day_length(i), rmonlen(m):", 3i4,2f12.6)') m,i,ii,day_length(i),rmonlen(m)
        end do
        ryeartot = ryeartot + rmonlen(m)
    end do
    write (23,'("rmonlen: ",13f12.6)') rmonlen(1:nm), ryeartot

end subroutine kg_monlen

subroutine adjust_to_reference(rmonlen, rmonlenref, rmonlentarg)

    implicit none

    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: rmonlenref(nm)   ! reference month lengths (usually calculated 0 ka)
    real(8), intent(in)     :: rmonlentarg(nm)  ! target month lengths (usually conventional 0 ka)

    integer(4)              :: m

    ! adjust rmonlen to reference year
    do m=1,nm
        rmonlen(m) = rmonlen(m) * (rmonlentarg(m)/rmonlenref(m))
    end do

end subroutine adjust_to_reference

subroutine adjust_to_yeartot(rmonlen, ryeartottarg, ryeartot)

    implicit none

    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: ryeartottarg     ! total annual month lengths target
    real(8), intent(out)    :: ryeartot         ! total annual month lengths

    integer(4)              :: m

    ! get ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

    ! adjust rmonlen to ryeartot target
    do m=1,nm
        rmonlen(m) = rmonlen(m) * (ryeartottarg/ryeartot)
    end do

    ! recalc ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

end subroutine adjust_to_yeartot

subroutine integer_monlen(rmonlen, ndtarg, imonlen, iyeartot)

    implicit none

    real(8), intent(in)     :: rmonlen(nm)      ! real month lengths
    integer(4), intent(in)  :: ndtarg           ! total annual integer number of days target
    integer(4), intent(out) :: imonlen(nm)      ! integer month lengths
    integer(4), intent(out) :: iyeartot

    integer(4)              :: diff_sign
    real(8)                 :: ryeartot
    real(8)                 :: inc
    integer(4)              :: i,m

    ! integer month lengths
    iyeartot=0; ryeartot=0.0d0
    do  m=1,nm
        imonlen(m)=idint(dnint(rmonlen(m)))
        iyeartot=iyeartot+imonlen(m)
        ryeartot=ryeartot+rmonlen(m)
    end do
    !write (23,'(13(i4,7x))') imonlen(1:nm),iyeartot

    ! force integer month lengths to sum to ndtarg
    if (iyeartot.ne.ndtarg) then

        write (22,'(13f12.6)') rmonlen(1:nm),ryeartot
        write (22,'(13(i4,7x))') imonlen(1:nm),iyeartot

        ! add (diff_sign = 1) or subtract (diff_sign = -1) a small increment to each rmonlen value
        ! iterate until iyeartot = ndtarg, incrementing small increment via i (i*inc)
        diff_sign=1
        if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1
        !write (20,*) diff_sign
        inc=0.000001; i=0
        do
            if (iyeartot.eq.ndtarg) exit
            i=i+1
            iyeartot=0
            do m=1,nm
                imonlen(m)=idint(dnint(rmonlen(m)+i*inc*diff_sign))
                iyeartot=iyeartot+imonlen(m)
            end do
            !write (20,'(i4,i3,f10.6,i4)') i,diff_sign,i*inc*diff_sign,iyeartot

            diff_sign=1
            if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1

        end do

        write (22,'(13(i4,7x))') imonlen(1:nm),iyeartot
        write (22,'(a)')
    end if

end subroutine integer_monlen

subroutine rmon_midbegend(rmonlen, midMarch, rmonmid, rmonbeg, rmonend)

    implicit none

    real(8), intent(in)  :: rmonlen(nm)  ! real month length
    real(8), intent(in)  :: midMarch     ! real vernal equinox day
    real(8), intent(out) :: rmonmid(nm)  ! real mid-month day
    real(8), intent(out) :: rmonbeg(nm)  ! real beginning day of month
    real(8), intent(out) :: rmonend(nm)  ! real ending day of month

    integer(4)              :: m

    ! month middle, beginning and end, relative to day veday (e.g. Mar 21 12:00 - 6 days = Mar 15 12:00)
    rmonmid(3) = midMarch
    rmonbeg(3) = rmonmid(3) - rmonlen(3)/2.0d0
    rmonend(3) = rmonmid(3) + rmonlen(3)/2.0d0
    do m=2,1,-1
        rmonmid(m)=rmonmid(m+1)-(rmonlen(m+1)/2.0d0)-(rmonlen(m)/2.0d0) ! Jan-Feb midmonth days
        rmonbeg(m)=rmonmid(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmonmid(m)+(rmonlen(m)/2.0d0)
    end do
    do m=4,12
        rmonmid(m)=rmonmid(m-1)+(rmonlen(m-1)/2.0d0)+(rmonlen(m)/2.0d0) ! Apr-Dec midmonth days
        rmonbeg(m)=rmonmid(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmonmid(m)+(rmonlen(m)/2.0d0)
    end do

end subroutine rmon_midbegend

subroutine imon_midbegend(imonlen, imidMarch, imonmid, imonbeg, imonend)

    implicit none

    integer(4), intent(in)  :: imonlen(nm)  ! integer month length
    integer(4), intent(in)  :: imidMarch    ! integer midMarch day
    integer(4), intent(out) :: imonmid(nm)  ! integer mid-month day
    integer(4), intent(out) :: imonbeg(nm)  ! integer beginning day of month
    integer(4), intent(out) :: imonend(nm)  ! integer ending day of month

    integer(4)              :: m

    ! month middle, beginning and end, relative to day iveday-6 (e.g. Mar 21 - 6 = Mar 15)
    imonmid(3) = imidMarch
    imonbeg(3) = imonmid(3) - floor(dble(imonlen(3))/2.0d0) + 1
    imonend(3) = imonmid(3) + imonlen(3)-floor(dble(imonlen(3))/2.0d0)
    m=3

    do m=2,1,-1
        imonbeg(m)=imonbeg(m+1) - imonlen(m)
        imonend(m)=imonbeg(m+1) - 1
        imonmid(m)=imonbeg(m) + floor(dble(imonlen(3))/2.0d0)
    end do
    do m=4,12
        imonbeg(m)=imonend(m-1) + 1
        imonend(m)=imonbeg(m) + imonlen(m) - 1
        imonmid(m)=imonbeg(m) + floor(dble(imonlen(3))/2.0d0)
    end do

end subroutine imon_midbegend

subroutine compare_monthdefs(nyrs, imondef, rmondef)
! used for debugging

    implicit none

    integer(4), intent(in)  :: nyrs
    integer(4), intent(in)  :: imondef(nyrs,nm)    ! integer month definition
    real(8), intent(in)     :: rmondef(nyrs,nm)    ! real month definition

    real(8)                 :: diff(nm),ssq(nm),mse(nm)
    integer(4)              :: n,m

    diff=0.0d0; ssq=0.0d0; mse=0.0d0

    do m=1,nm
        do n=1,nyrs
            diff(m) = diff(m)+dble(imondef(n,m))-rmondef(n,m)
            ssq(m) = ssq(m)+(dble(imondef(n,m))-rmondef(n,m))*(dble(imondef(n,m))-rmondef(n,m))
            !write (24,'("m,n,imondef,rmondef,diff: ", 2i4,3f12.6)') m,n,dble(imondef(n,m)),rmondef(n,m),diff(m)
        end do
    end do
    mse(:) = dsqrt(ssq(:)/dble(nyrs))

    write (24,'("diff: ",12f12.6)') diff(1:nm)
    write (24,'(" ssq: ",12f12.6)') ssq(1:nm)
    write (24,'(" mse: ",12f12.6)') mse(1:nm)

end subroutine compare_monthdefs

integer(4) function yearlen_BP(ageBP)

! gets number of days in a BP year -- no year zero or Gregorian-Julian adjustment

    implicit none

    integer(4), intent(in)  :: ageBP   ! Year BP (negative, e.g. 1900 CE = -50 BP)
    integer(4)              :: yearCE

    yearCE = ageBP + 1950 ! no Year 0 BCE/CE adjustment

    yearlen_BP = 365
    if (mod(yearCE, 4) .eq. 0) yearlen_BP = 366
    if (mod(yearCE, 100) .eq. 0) yearlen_BP = 365
    if (mod(yearCE, 400) .eq. 0) yearlen_BP = 366

end function yearlen_BP

integer(4) function yearlen_CE(yearCE)

! gets number of days in a CE year -- no year zero or Gregorian-Julian adjustment

    implicit none

    integer(4), intent(in)  :: yearCE   ! YearCE

    yearlen_CE = 365
    if (mod(yearCE, 4) .eq. 0) yearlen_CE = 366
    if (mod(yearCE, 100) .eq. 0) yearlen_CE = 365
    if (mod(yearCE, 400) .eq. 0) yearlen_CE = 366

end function yearlen_CE

real(8) function radians(d)
! decimal degrees to radians

    implicit none

    real(8) d, pi
    pi=4.0d0*datan(1.0d0)
    radians=d*((2.0d0*pi)/360.0d0)

end function radians

real(8) function degrees(r)
! radians to decimal degrees

    implicit none

    real(8) r, pi
    pi=4.0d0*datan(1.0d0)
    degrees=r/((2.0d0*pi)/360.0d0)

end function degrees

end module month_length_subs

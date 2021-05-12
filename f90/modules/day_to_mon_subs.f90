module day_to_mon_subs
! subroutines for aggregating daily dato to monthly data
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_day = .false.  
    
    
contains
    
subroutine day_to_rmon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xd,xfill,xm_adj)
! aggregation of pseudo- or actual daily data to real-length months using a paleo calendar

    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    real(8), intent(in)         :: rmonbeg(ny,nm), rmonend(ny,nm)   ! beginning and ending days of each month
    real(8), intent(in)         :: xd(ndtot)        ! daily values
    real(8), intent(in)         :: xfill            ! _FillValue
    real(8), intent(out)        :: xm_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ibegday, iendday             ! beginning day and ending day of each year
    integer(4)              :: ibeg(nm), iend(nm)           ! beginning and ending (integer) day of each month
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    real(8)                 :: xdx(-29:nd+30)               ! daily data for current year, padded by data from adjacent years
    real(8)                 :: wgt(-29:nd+30), wsum         ! weights (for interpolating over fractional days)
    integer(4)              :: nfill                        ! number of days with fill values
    
    integer(4)              :: n, m, i, nn
    
    logical                 :: debug_write_cal_effects = .false.

    
    !write (*,'(a)') "in day_to_mon_ts"
    if (debug_day) write (debug_unit,'(a)') "day_to_mon"
    if (debug_day) write (debug_unit,*) ny, ndtot
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    iendday = 0; nn = 0
    xm_adj=0.0d0
    do n=1,ny
        if (debug_day) write (debug_unit,'("n, ndays:", 2i5)') n,ndays(n)
        ibegday = iendday + 1
        iendday = ibegday + ndays(n) - 1
        if (debug_day) write (debug_unit,*) n,ibegday,iendday
        
        if (ny .eq. 1) then       ! single-year Aclim data  
            ! wrap the input daily data
            xdx(-29:0)=xd(ndays(n)-30+1:ndays(n))
            xdx(1:ndays(n))=xd(1:ndays(n))
            xdx(ndays(n)+1:ndays(n)+30)=xd(1:30)
        else 
            ! copy current year into xdx
            xdx(1:ndays(n)) = xd(ibegday:iendday)
            ! pad beginning and end of xdx
            if (n .eq. 1) then
                xdx(-29:0) = xd(ndays(n)-30+1:ndays(n))
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            elseif (n .eq. ny) then
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(ibegday+1:ibegday+30)
            else
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            end if
        end if
    
        ! integer beginning and end of each month, and number of days in each month
        ! ndays_in_month should be equal to the integer month length + 1
        ibeg=ceiling(rmonbeg(n,:)); iend=ceiling(rmonend(n,:)); ndays_in_month=(iend-ibeg+1)
        if (debug_day) write (debug_unit,'("rmonbeg:  ",12f14.8)') rmonbeg(n,:)
        if (debug_day) write (debug_unit,'("rmonend:  ",12f14.8)') rmonend(n,:)
        if (debug_day) write (debug_unit,'("ibeg:  ",12i4)') ibeg
        if (debug_day) write (debug_unit,'("iend:  ",12i4)') iend
        if (debug_day) write (debug_unit,'("ndays: ",13i4)') ndays_in_month, sum(ndays_in_month)
 
        ! monthly means
        do m=1,nm
            nn = nn + 1
            nfill = 0; wgt=1.0d0; wsum=0.0d0
            wgt(ibeg(m))=abs(rmonbeg(n,m)-dble(ibeg(m)))
            wgt(iend(m))=abs(rmonend(n,m)-dble(iend(m)-1))
            do i=ibeg(m),iend(m)
                if (xdx(i) .ne. xfill) then
                    xm_adj(nn)=xm_adj(nn)+xdx(i)*wgt(i)
                    wsum=wsum+wgt(i)
                    if (debug_day) &
                    write (debug_unit,'("m, i, xd(i), xm_adj(nn), wgt(i), wsum: ",2i4,2f12.6,2f12.6)') &
                        m,i,xdx(i),xm_adj(nn),wgt(i),wsum
                else
                    nfill = nfill + 1
                    if (debug_day) write (debug_unit, *) m, i, xdx(i), nfill
                end if
            end do
            if (debug_day) write (debug_unit, *) wsum, nfill, xm_adj(nn)
            if (wsum .ne. 0.0d0 .and. nfill .eq. 0) then
                xm_adj(nn)=xm_adj(nn)/wsum
            else
                xm_adj(nn)=xfill
            end if
            !if (debug_day) write(10,'("m, xm_adj(nn): ",i4,2f14.6)') &
            !    m,xm_adj(nn),sngl(xm_adj(nn))
            if (debug_day) write (debug_unit, *) m,xm_adj(nn),sngl(xm_adj(nn))
        end do
    end do 
    if (debug_day) write (debug_unit,'(12g14.6)') xm_adj

end subroutine day_to_rmon_ts

subroutine day_to_imon_ts(ny,ndays,imonbeg,imonend,ndtot,xd,xfill,xm_adj)
! aggregation of pseudo- or actual daily data to months using a paleo calendar

    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    integer(4), intent(in)      :: imonbeg(ny,nm), imonend(ny,nm)   ! beginning and ending days of each month
    real(8), intent(in)         :: xd(ndtot)        ! daily values
    real(8), intent(in)         :: xfill            ! _FillValue
    real(8), intent(out)        :: xm_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    !real(8)                 :: xdx(-29:nd+31)               ! daily data for current year, padded by data from adjacent years
    integer(4)              :: nfill                        ! number of days with fill value
    
    integer(4)              :: n, m, i, nn, npts

    !write (*,'(a)') "in day_to_mon_ts"
    if (debug_day) write (debug_unit,'(a)') "day_to_mon"
    if (debug_day) write (debug_unit,*) ny, ndtot
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    nn = 0
    xm_adj=0.0d0
    do n=1, 2 !ny
        ! monthly means
        do m=1,nm
            nn = nn + 1
            nfill = 0; npts = 0
            do i=imonbeg(n,m),imonend(n,m)
                if (xd(i) .ne. xfill) then
                    xm_adj(nn)=xm_adj(nn)+xd(i)
                    npts=npts+1
                    if (debug_day) &
                    write (debug_unit,'("m, i, xd(i), xm_adj(nn), npts: ",2i8,2f12.6,i5)') &
                        m,i,xd(i),xm_adj(nn),npts
                else
                    nfill = nfill + 1
                    if (debug_day) write (debug_unit, *) m, i, xd(i), nfill
                end if
            end do
            if (debug_day) write (debug_unit, *) npts, nfill, xm_adj(nn)
            if (npts .ne. 0.0d0 .and. nfill .eq. 0) then
                xm_adj(nn)=xm_adj(nn)/dble(npts)
            else
                xm_adj(nn)=xfill
            end if
            !if (debug_day) write(10,'("m, xm_adj(nn): ",i4,2f14.6)') &
            !    m,xm_adj(nn),sngl(xm_adj(nn))
            if (debug_day) write (debug_unit, *) m,xm_adj(nn),sngl(xm_adj(nn))
        end do
    end do 
    !if (debug_day) write (debug_unit,'(12g14.6)') xm_adj

end subroutine day_to_imon_ts

end module day_to_mon_subs

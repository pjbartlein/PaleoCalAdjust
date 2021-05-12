module pseudo_daily_interp_subs
! subroutines for mean-preserving interpolation
    
    implicit none

    integer             :: debug_unit = 10
    integer             :: out_unit = 6
    logical             :: debug_em = .false.
    logical             :: debug_day = .false.  

contains
    
subroutine enforce_mean(nctrl, ntargs, nsubint, tol, ym_targ, yd, ymiss)
! adjustes daily values to match target mean

    implicit none
    
    integer(4), parameter   :: maxtargs = 7305
    integer(4), intent(in)  :: nctrl
    integer(4), intent(in)  :: ntargs
    integer(4), intent(in)  :: nsubint(nctrl)
    real(8), intent(in)     :: tol
    real(8), intent(in)     :: ym_targ(nctrl)
    real(8), intent(inout)  :: yd(ntargs)
    real(8), intent(in)     :: ymiss
    !
    integer(4), parameter   :: maxwgts = 31, maxiter=100
    real(8)                 :: ydm(nctrl)
    real(8)                 :: nonzero_sum(nctrl), yd_adjust_sum(nctrl), yd_next_sum(nctrl)
    real(8)                 :: yd_next_mean(nctrl), new_diff(nctrl)
    real(8)                 :: diff(nctrl), yd_adjust
    
    integer(4)              :: ii, ib, ie, j, jj, k, niter(nctrl), ntotiter
    
    if (debug_em) write (debug_unit,'(a)') "In enforce_mean"
    if (debug_em) write (debug_unit,'("tol: ", f8.5)') tol
    if (debug_em) write (debug_unit,'("nctrl, ntargs: ", 2i8/)') nctrl, ntargs
        
    ! get means of current daily values
    
    call interval_mean(nctrl, nsubint, ntargs, yd, ymiss, ydm)
    
    ii = 0
    ntotiter = 1
    yd_adjust_sum = 0.0d0;  yd_next_sum = 0.0d0
    do k = 1, nctrl
        niter(k) = 0
        if (debug_em) write (debug_unit, '(/"Interval ", i8)') k 
        ib = ii + 1; ie = ib + nsubint(k) - 1
        
        if (debug_em) write (debug_unit, '("ym_targ,ydm,ym_targ-ydm: ",3g16.9)') ym_targ(k), ydm(k), ym_targ(k) - ydm(k)
        
        ! if ym(k) = 0.0, set all days to 0.0
        if (ym_targ(k) .eq. 0.0d0) then
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays ym_targ"
            if (debug_em) write(debug_unit, '(4i8, g16.9)')  k, ib, ie, ie - ib + 1, ym_targ(k)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to 0.0"
            yd(ib:ie) = 0.0d0
            if (debug_em) write (debug_unit, '(10g14.6)') yd(ib:ie)
            
        else if (ym_targ(k) .eq. ymiss .or. ydm(k) .eq. ymiss) then
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays ym_targ"
            if (debug_em) write(debug_unit, '(4i8, g16.9)')  k, ib, ie, ie - ib + 1, ym_targ(k)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to missing"
            yd(ib:ie) = ymiss
            if (debug_em) write (debug_unit, '(10g16.9)') yd(ib:ie)
        else    
            
            ! check for differences in actual mean and mean of interpolated values
            diff(k) = ym_targ(k) - ydm(k)
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays ym_targ"// &
                "         ydm            diff         tol  dabs(diff)"
            if (debug_em) write(debug_unit, '( 4i8, 5g16.9)')  & 
                k, ib, ie, ie - ib + 1, ym_targ(k), ydm(k), diff(k), tol, dabs(diff(k))
            
            ! if absolute value of difference is greater than tol, adjust subinterval values            
            if (dabs(diff(k)) .gt. tol) then
                
                if (debug_em) write (debug_unit, '(a)') "Adjusting subinterval values..."
                
                ! get sum of nonzero values
                if (debug_em) write (debug_unit, '(a)') "       j      jj          yd nonzero_sum"
                nonzero_sum = 0.0d0     
                jj = 0
                do j = ib, ie
                    jj = jj + 1
                    if (yd(j) .ne. 0.0d0)  nonzero_sum(k) = nonzero_sum(k) + yd(j)
                    if (debug_em) write (debug_unit, '(2i8, 2g16.9)') j, jj, yd(j), nonzero_sum(k)
                end do 
                    
                ! apply porportions of differences
                if (debug_em) write (debug_unit, '(a)') "       j      jj          yd   yd_adjust"
                jj = 0
                do j = ib, ie
                    jj = jj + 1
                    if (yd(j) .ne. 0.0d0) then
                        yd_adjust = ((yd(j) / nonzero_sum(k)) * diff(k)) * nsubint(k)
                        yd(j) = yd(j) + yd_adjust
                        yd_next_sum(k) = yd_next_sum(k) + yd(j)
                        yd_adjust_sum(k) = yd_adjust_sum(k) + yd_adjust
                        if (debug_em) write (debug_unit, '(2i8, 3f12.6)') j, jj, yd(j), yd_adjust
                    else
                        if (debug_em) write (debug_unit, '(2i8, 3f12.6)') j, jj, yd(j)
                    end if   
                end do
                yd_next_mean(k) = yd_next_sum(k) / nsubint(k)
                new_diff(k) = ym_targ(k) - yd_next_mean(k)
                if (debug_em) write (debug_unit, '(a)') "yd_adjust_sum yd_next_sum yd_next_mean  new_diff"
                if (debug_em) write (debug_unit, '(4f12.6)') yd_adjust_sum(k), yd_next_sum(k), yd_next_mean(k), new_diff(k)
            else  
            
                if (debug_em) write (debug_unit, '(a)') "Adopting subinterval values..."
            
            end if
            
        end if
        ii = ii + nsubint(k)

    end do
    
end subroutine enforce_mean

subroutine day_to_rmon_ts(ny,ndays,rmonbeg,rmonend,ndtot,yd,yfill,ym_adj)
! aggregation of pseudo- or actual daily data to real-length months using a paleo calendar

    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    real(8), intent(in)         :: rmonbeg(ny,nm), rmonend(ny,nm)   ! beginning and ending days of each month
    real(8), intent(in)         :: yd(ndtot)        ! daily values
    real(8), intent(in)         :: yfill            ! _FillValue
    real(8), intent(out)        :: ym_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ibegday, iendday             ! beginning day and ending day of each year
    integer(4)              :: ibeg(nm), iend(nm)           ! beginning and ending (integer) day of each month
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    real(8)                 :: ydx(-29:nd+30)               ! daily data for current year, padded by data from adjacent years
    real(8)                 :: wgt(-29:nd+30), wsum         ! weights (for interpolating over fractional days)
    integer(4)              :: nfill                        ! number of days with fill values
    
    integer(4)              :: n, m, i, nn
    
    logical                 :: debug_day = .false.

    if (debug_day) write (debug_unit,'(a)') "In day_to_rmon..."
    if (debug_day) write (debug_unit,*) ny, ndtot
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    iendday = 0; nn = 0
    ym_adj=0.0d0
    do n=1,ny
        if (debug_day) write (debug_unit,'("n, ndays:", 2i5)') n,ndays(n)
        ibegday = iendday + 1
        iendday = ibegday + ndays(n) - 1
        if (debug_day) write (debug_unit,*) n,ibegday,iendday
        
        if (ny .eq. 1) then       ! single-year Aclim data  
            ! wrap the input daily data
            ydx(-29:0)=yd(ndays(n)-30+1:ndays(n))
            ydx(1:ndays(n))=yd(1:ndays(n))
            ydx(ndays(n)+1:ndays(n)+30)=yd(1:30)
        else 
            ! copy current year into ydx
            ydx(1:ndays(n)) = yd(ibegday:iendday)
            ! pad beginning and end of ydx
            if (n .eq. 1) then
                ydx(-29:0) = yd(ndays(n)-30+1:ndays(n))
                ydx(ndays(n)+1:ndays(n)+30) = yd(iendday+1:iendday+30)
            elseif (n .eq. ny) then
                ydx(-29:0) = yd(ibegday-30:ibegday-1)
                ydx(ndays(n)+1:ndays(n)+30) = yd(ibegday+1:ibegday+30)
            else
                ydx(-29:0) = yd(ibegday-30:ibegday-1)
                ydx(ndays(n)+1:ndays(n)+30) = yd(iendday+1:iendday+30)
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
                if (ydx(i) .ne. yfill) then
                    ym_adj(nn)=ym_adj(nn)+ydx(i)*wgt(i)
                    wsum=wsum+wgt(i)
                    if (debug_day) &
                    write (debug_unit,'("m, i, yd(i), ym_adj(nn), wgt(i), wsum: ",2i4,2f12.6,2f12.6)') &
                        m,i,ydx(i),ym_adj(nn),wgt(i),wsum
                else
                    nfill = nfill + 1
                    if (debug_day) write (debug_unit, *) m, i, ydx(i), nfill
                end if
            end do
            if (debug_day) write (debug_unit, *) wsum, nfill, ym_adj(nn)
            if (wsum .ne. 0.0d0 .and. nfill .eq. 0) then
                ym_adj(nn)=ym_adj(nn)/wsum
            else
                ym_adj(nn)=yfill
            end if

            if (debug_day) write (debug_unit, *) m,ym_adj(nn),sngl(ym_adj(nn))
        end do
    end do 
    
    if (debug_day) write (debug_unit,'(12g14.6)') ym_adj

end subroutine day_to_rmon_ts

subroutine dayinterp(nm,nd,monlen,ym,yd)
! Interpolate pseudo-daily values of monthly data.  Not mean-preserving.

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: ym(nm)
    real(8), intent(out)    :: yd(nd)
    
    integer(4)              :: nm1
    !integer(4)              :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: ym2(0:nm+1)
    integer                 :: i,m
    
    !call midmonth_int(nm,monlen,midmon)
    call midmonth_real(nm,dble(monlen),midmon)
    
    ! pad data at beginning (m=0) and end (m=13)
    nm1=nm+1
    ym2(1:nm)=ym
    ym2(0)=ym(nm)
    ym2(nm1)=ym(1)
    midmon2(1:nm)=midmon
    midmon2(0)=1-(nd-midmon(nm))-2
    midmon2(nm1)=nd+midmon(1)   
    
    ! linear pseudo-daily interpolation
    do i=1,nd
        ! find month day i lies in
        do m=1,nm+1
            if (i.gt.midmon2(m-1) .and. i.le.midmon2(m)) exit
        end do   
        yd(i)=(dble(i-midmon2(m-1))/dble(midmon2(m)-midmon2(m-1)))*(ym2(m)-ym2(m-1))+ym2(m-1)
    end do
    
end subroutine dayinterp

subroutine dayspread(nm,nd,monlen,ym,yd)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: ym(nm)
    real(8), intent(out)    :: yd(nd)
    
    integer                 :: i,j,m
    
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            yd(i)=ym(m)
        end do
    end do
    
end subroutine dayspread

subroutine midmonth_int(nm,monlen,midmon)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    integer(4), intent(in)  :: monlen(nm)
    integer(4), intent(out) :: midmon(nm)
    
    integer(4)              :: m,endday(nm)

    ! midmonth day numbers
    m=1
    midmon(m)=ceiling(dble(monlen(m))/2.0d0)
    endday(m)=monlen(m)
    do m=2,nm
        midmon(m)=ceiling(dble(monlen(m))/2.0d0)+endday(m-1)
        endday(m)=endday(m-1)+monlen(m)
    end do
    
end subroutine midmonth_int

subroutine midmonth_real(nm,monlen,midmon)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    real(8), intent(in)     :: monlen(nm)
    real(8), intent(out)    :: midmon(nm)
    
    real(8)              :: endday(nm)
    integer(4)           :: m

    ! midmonth day numbers
    m=1
    midmon(m)=monlen(m)/2.0d0
    endday(m)=monlen(m)
    do m=2,nm
        midmon(m)=(monlen(m)/2.0d0)+endday(m-1)
        endday(m)=endday(m-1)+monlen(m)
    end do
    
end subroutine midmonth_real

subroutine mon_mean(nm,nd,monlen,yd,ym)
! gets monthly means of interpolated daily data

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: yd(nd)
    real(8), intent(out)    :: ym(nm)
    integer                 :: i,j,m
    
    ym=0.0d0
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            ym(m)=ym(m)+yd(i)
        end do
        ym(m)=ym(m)/dble(monlen(m))
    end do

end subroutine mon_mean

subroutine ann_wmean(n,x,w,ym)
! gets a weighted (by month length) annual mean value from monthly data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n),w(n)
    real(8), intent(out)    :: ym
    real(8)                 :: wsum
    integer                 :: i
    
    ym=0.0d0; wsum=0.0d0
    do i=1,n
        ym=ym+x(i)*w(i)
        wsum=wsum+w(i)
    end do

    ym=ym/wsum
    
end subroutine ann_wmean
    
subroutine ann_mean(n,x,ym)
! gets an unweighted annual mean value from monthly or daily data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n)
    real(8), intent(out)    :: ym
    integer                 :: i
    
    ym=0.0d0
    do i=1,n
        ym=ym+x(i)
    end do

    ym=ym/dble(n)
    
end subroutine ann_mean   

subroutine dzero(nm,nd,monlen,ym,yd0)
! enforces 0.0 values of interpolated daily data when the monthly mean is 0.0

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: ym(nm)
    real(8), intent(inout)  :: yd0(nd)
    
    integer(4)              :: i,m,j

    ! zero daily values in months were ym=0.0
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            if (ym(m).eq.0.0) yd0(i)=0.0
        end do
    end do
    
    ! zero all other daily values
    do i=1,nd
        if (yd0(i).le.0.0) yd0(i)=0.0
    end do

end subroutine dzero

subroutine interval_mean(ninterval, nsubint, n, y_int, ymiss, ym_int)
    
    implicit none
        
    integer(4), intent(in)          :: ninterval, n
    integer(4), intent(in)          :: nsubint(ninterval)          ! interval lengths
    real(8), intent(in)             :: y_int(n)                    ! input values
    real(8), intent(in)             :: ymiss                       ! missing/fill value
    real(8), intent(out)            :: ym_int(ninterval)           ! output means
    
    integer(4)                      :: i, j, k, npts
    
    ym_int = 0.0d0
    i = 0
    do k = 1, ninterval
        npts = 0
        do j = 1, nsubint(k)
            i = i + 1 
            if (y_int(i) .ne. ymiss) then
                ym_int(k) = ym_int(k) + y_int(i)
                npts = npts + 1 
            end if
        end do

        if (npts .gt. 0) then     
            ym_int(k) = ym_int(k) / dble(npts)
        else 
            ym_int(k) = ymiss
        end if
    end do
    
end subroutine interval_mean

subroutine interp_stat(nctrl, ym, ym_int, rmse)

    implicit none
    
    integer(4), intent(in)      :: nctrl
    real(8), intent(in)         :: ym(nctrl), ym_int(nctrl)
    real(8), intent(out)        :: rmse
    
    integer(4)                  :: k
    
    rmse = 0.0d0

    do k = 1, nctrl
        rmse = rmse + (ym(k) - ym_int(k))*(ym(k) - ym_int(k))
    end do
    rmse = sqrt(rmse / (nctrl - 1))
    
end subroutine interp_stat

subroutine step_plot(nctrl, ntarg, nsubint, ym, x, yint)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nctrl, ntarg
    integer(4), intent(in)  :: nsubint(nctrl)
    real(8), intent(in)     :: ym(nctrl)
    real(8), intent(out)    :: x(ntarg)   ! adjusted abcissas to make nice step plot
    real(8), intent(out)    :: yint(ntarg)
    
    real(8)                 :: halfstep
    
    integer                 :: i, j, k
    
    i = 0
    do k = 1, nctrl
        do j = 1, nsubint(k)
            i = i + 1
            x(i) = dble(i)
            yint(i) = ym(k)
        end do
    end do
    
    ! adjust abcissas
    
    halfstep = (x(2) - x(1)) / 2.0d0
    i = 0
    do k = 1, nctrl - 1
        i = i + nsubint(k)
        x(i) = x(i) + halfstep
        x(i + 1) = x(i)
    end do
   
end subroutine step_plot

end module pseudo_daily_interp_subs

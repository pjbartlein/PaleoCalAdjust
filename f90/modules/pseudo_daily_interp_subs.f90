module pseudo_daily_interp_subs
! subroutines for mean-preserving interpolation
    
    implicit none

    integer             :: debug_unit = 10
    integer             :: out_unit = 6
    logical             :: debug_em = .false.
    logical             :: debug_day = .false.  

contains
    
subroutine enforce_mean(nctrl, ntargs, nsubint, tol, ym, yd, ymiss)
! adjustes daily values to match target mean

    implicit none
    
    integer(4), intent(in)  :: nctrl
    integer(4), intent(in)  :: ntargs
    integer(4), intent(in)  :: nsubint(nctrl)
    real(8), intent(in)     :: tol
    real(8), intent(in)     :: ym(nctrl)
    real(8), intent(inout)  :: yd(ntargs)
    real(8), intent(in)     :: ymiss
    
    real(8)                 :: ydm(nctrl)
    real(8)                 :: nonzero_sum, nonzero_mean, yd_adjust_sum, yd_adjust_mean, yd_out_sum, yd_out_mean
    real(8)                 :: diff, yd_old, yd_adjust, new_diff
    
    integer(4)              :: i, ii, ib, ie, j, jj, n, nn, nzero 
    
    if (debug_em) write (debug_unit,'(a)') "In enforce_mean"
    if (debug_em) write (debug_unit,'("tol: ", f8.5)') tol
    if (debug_em) write (debug_unit,'("nctrl, ntargs: ", 2i8/)') nctrl, ntargs
        
    ! get means of current daily values    
    call interval_mean(nctrl, nsubint, ntargs, yd, ymiss, ydm)
    
    nn = 0
    do n = 1, nctrl
        if (debug_em) write (debug_unit, '(/"Interval ", 2i8)') n, nn
        
        ! beginning and ending subintervals
        ib = nn + 1; ie = ib + nsubint(n) - 1    
        ! update for next iteration
        nn = nn + nsubint(n)
        
        if (debug_em) write (debug_unit, '("ym,ydm,ym-ydm: ",3f12.6)') ym(n), ydm(n), ym(n) - ydm(n)
        if (debug_em) write (debug_unit, *) ym(n)
        
        ! if ym(n) = 0.0, set all days to 0.0
        if (ym(n) .eq. 0.0d0) then
            if (debug_em) write(debug_unit, '(a)') "ym(n) = 0.0d0"
            if (debug_em) write(debug_unit, '(a)') "       n      ib      ie   ndays          ym"
            !if (debug_em) write(debug_unit, '(4i8, f12.6)')  n, ib, ie, ie - ib + 1, ym(n)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to 0.0"
            yd(ib:ie) = 0.0d0
            if (debug_em) write (debug_unit, '(10f12.6)') yd(ib:ie)
            cycle
        end if
        
        ! or set them to missing
        if (ym(n) .eq. ymiss .or. ydm(n) .eq. ymiss) then
            if (debug_em) write(debug_unit, '(a)') "       n      ib      ie   ndays          ym"
            if (debug_em) write(debug_unit, '(4i8, f12.6)')  n, ib, ie, ie - ib + 1, ym(n)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to missing"
            yd(ib:ie) = ymiss
            if (debug_em) write (debug_unit, '(10g16.9)') yd(ib:ie)
            cycle
        end if
        
        ! otherwise, adjust or adopt the daily values
            
        ! check for differences in actual mean and mean of interpolated values
        diff = ym(n) - ydm(n)
        if (debug_em) write(debug_unit, '(a)') "       n      ib      ie   ndays          ym"// &
            "         ydm        diff         tol  dabs(diff)"
        if (debug_em) write(debug_unit, '( 4i8, 5f12.6)')  & 
            n, ib, ie, ie - ib + 1, ym(n), ydm(n), diff, tol, dabs(diff)
            
        ! if absolute value of difference is greater than tol, adjust subinterval values            
        if (dabs(diff) .gt. tol) then
                
            if (debug_em) write (debug_unit, '(a)') "Adjusting subinterval values..."
                
            ! get sum, total number, and mean of nonzero values
            if (debug_em) write (debug_unit, '(a)') "       i      ii          yd nonzero_sum, nzero"
            nzero = 0; nonzero_sum = 0.0d0     
            ii = 0
            do i = ib, ie
                ii = ii + 1
                if (yd(i) .ne. 0.0d0)  then
                    nonzero_sum = nonzero_sum + yd(i)
                    nzero = nzero + 1
                end if
                if (debug_em) write (debug_unit, '(2i8, 2f12.6, i4)') i, ii, yd(i), nonzero_sum, nzero
            end do 
                
            if (nzero .gt. 0) nonzero_mean = nonzero_sum / dble(nzero)
            if (debug_em) write (debug_unit, '("nzero_sum, nzero, nzero_mean, diff: ", f12.6, i4, 2f12.6)') &
                nonzero_sum, nzero, nonzero_mean, diff
                    
            ! apply porportions of difference to each nonzero subinterval
            if (debug_em) write (debug_unit, '(a)') "       i      ii      yd_old   yd_adjust   yd(i)"
            ii = 0; yd_adjust_sum = 0.0d0; yd_out_sum = 0.0d0; yd_adjust_mean = 0.0d0; yd_out_mean = 0.0d0
            do i = ib, ie
                yd_old = yd(i)
                ii = ii + 1
                if (yd(i) .ne. 0.0d0) then
                    yd_adjust = yd_old * ((dble(nsubint(n)) * diff) / nonzero_sum)
                    yd(i) = yd(i) + yd_adjust
                    yd_out_sum = yd_out_sum + yd(i)
                    yd_adjust_sum = yd_adjust_sum + yd_adjust
                    if (debug_em) write (debug_unit, '(2i8, 3f12.6)') i, ii, yd_old, yd_adjust, yd(i)
                else
                    if (debug_em) write (debug_unit, '(2i8, 3f12.6)') i, ii, yd_old, 0.0d0, yd(i)
                end if   
            end do
            yd_out_mean = yd_out_sum / nsubint(n)
            yd_adjust_mean = yd_adjust_sum / nsubint(n)
            new_diff = ym(n) - yd_out_mean
            if (debug_em) write (debug_unit, '(a)') "yd_adjust_sum, yd_out_sum, yd_out_mean,    ym(n),   new_diff"
            if (debug_em) write (debug_unit, '(5f12.6)') yd_adjust_sum, yd_out_sum, yd_out_mean, ym(n), new_diff
        else  
            
            if (debug_em) write (debug_unit, '(a)') "Adopting subinterval values..."
            
        end if

    end do
    
end subroutine enforce_mean

subroutine day_to_rmon_ts(ny, ndays, imonbeg, imonend, rmonbeg, rmonend, ndtot, idaynum, yd, ym, yfill, ym_adj)
! aggregation of pseudo- or actual daily data to real-length months using a paleo calendar

    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    integer(4), intent(in)      :: imonbeg(ny,nm), imonend(ny,nm)   ! integer beginning and ending days of each maonth
    real(8), intent(in)         :: rmonbeg(ny,nm), rmonend(ny,nm)   ! real beginning and ending days of each month
    integer(8), intent(in)      :: idaynum(ndtot)   ! integer day number
    real(8), intent(in)         :: yd(ndtot)        ! daily values
    real(8), intent(in)         :: ym(ny*nm)        ! monthly means
    real(8), intent(in)         :: yfill            ! _FillValue
    real(8), intent(out)        :: ym_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ibegday, iendday             ! beginning day and ending day of each year
    integer(4)              :: ibeg(nm), iend(nm)           ! beginning and ending (integer) day of each month
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    integer(8)              :: idaynumx(-29:nd+30)          ! integer day number for current year, padded 
    real(8)                 :: ydx(-29:nd+30)               ! daily data for current year, padded by data from adjacent years
    real(8)                 :: wgt(-29:nd+30), wsum         ! weights (for interpolating over fractional days)
    integer(4)              :: nfill                        ! number of days with fill values
    
    integer(4)              :: n, m, i, nn, ii

    if (debug_day) write (debug_unit,'(a)') "In day_to_rmon..."
    if (debug_day) write (debug_unit,*) ny, ndtot
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    iendday = 0; nn = 0
    ym_adj=0.0d0
    do n=1,ny
        if (debug_day) write (debug_unit, '(/"n, ndays:", i6, i4)') n,ndays(n)
        if (debug_day) write (debug_unit, '("imonbeg(n,1), imonbeg(n,nm), rmonbeg(n,1), rmonend(n,nm): ", 4f12.6)') &
            imonbeg(n,1), imonbeg(n,nm), rmonbeg(n,1), rmonend(n,nm)
        ibegday = iendday + 1
        iendday = ibegday + ndays(n) - 1
        if (debug_day) write (debug_unit, '("ibegday, iendday: ", 2i8)') ibegday,iendday
        
        if (ny .eq. 1) then       ! single-year Aclim data  
            ! wrap the input daily data
            ydx(-29:0)=yd(ndays(n)-30+1:ndays(n))
            ydx(1:ndays(n))=yd(1:ndays(n))
            ydx(ndays(n)+1:ndays(n)+30)=yd(1:30)
            idaynumx(-29:0)=yd(ndays(n)-30+1:ndays(n))
            idaynumx(1:ndays(n))=yd(1:ndays(n))
            idaynumx(ndays(n)+1:ndays(n)+30)=yd(1:30)
        else 
            ! copy current year into ydx
            ydx(1:ndays(n)) = yd(ibegday:iendday)
            idaynumx(1:ndays(n)) = idaynum(ibegday:iendday)
            ! pad beginning and end of ydx
            if (n .eq. 1) then
                ydx(-29:0) = yd(ndays(n)-30+1:ndays(n))
                ydx(ndays(n)+1:ndays(n)+30) = yd(iendday+1:iendday+30)
                idaynumx(-29:0) = idaynum(ndays(n)-30+1:ndays(n))
                idaynumx(ndays(n)+1:ndays(n)+30) = idaynum(iendday+1:iendday+30)
            elseif (n .eq. ny) then
                ydx(-29:0) = yd(ibegday-30:ibegday-1)
                ydx(ndays(n)+1:ndays(n)+30) = yd(ibegday+1:ibegday+30)
                idaynumx(-29:0) = idaynum(ibegday-30:ibegday-1)
                idaynumx(ndays(n)+1:ndays(n)+30) = idaynum(ibegday+1:ibegday+30)
            else
                ydx(-29:0) = yd(ibegday-30:ibegday-1)
                ydx(ndays(n)+1:ndays(n)+30) = yd(iendday+1:iendday+30)
                idaynumx(-29:0) = idaynum(ibegday-30:ibegday-1)
                idaynumx(ndays(n)+1:ndays(n)+30) = idaynum(iendday+1:iendday+30)
            end if
        end if
    
        ! integer beginning and end of each month, and number of days in each month
        ! ndays_in_month should be equal to the integer month length + 1
        ibeg=ceiling(rmonbeg(n,:)); iend=ceiling(rmonend(n,:)); ndays_in_month=(iend-ibeg+1)

        if (debug_day) write (debug_unit,'("rmonbeg:  ",12f12.6)') rmonbeg(n,:)
        if (debug_day) write (debug_unit,'("rmonend:  ",12f12.6)') rmonend(n,:)
        if (debug_day) write (debug_unit,'("1monbeg:  ",12i12)') imonbeg(n,:)
        if (debug_day) write (debug_unit,'("1monend:  ",12i12)') imonend(n,:)
        if (debug_day) write (debug_unit,'("ibeg:  ",12i12)') ibeg
        if (debug_day) write (debug_unit,'("iend:  ",12i12)') iend
        if (debug_day) write (debug_unit,'("ndays: ",13i12)') ndays_in_month, sum(ndays_in_month)
 
        ! monthly means
        do m=1,nm
            nn = nn + 1
            nfill = 0; wgt=1.0d0; wsum=0.0d0
            wgt(ibeg(m)) = abs(dble(ibeg(m)) - rmonbeg(n,m))
            wgt(iend(m)) = 1.0d0 - abs(rmonend(n,m) - dble(iend(m)))
            if (debug_day) & 
                write (debug_unit,'(/, "m,rmonbeg,ibeg,rmonend,iend,wgt(ibeg(m)),wgt(iend(mm)): ", i3, 2(f14.7,i4), 2f12.6)') &
                    m,rmonbeg(n,m),ibeg(m),rmonend(n,m),iend(m),wgt(ibeg(m)),wgt(iend(m))
            ii = 0
            do i=ibeg(m),iend(m)
                ii = ii + 1
                if (ydx(i) .ne. yfill) then
                    ym_adj(nn)=ym_adj(nn)+ydx(i)*wgt(i)
                    wsum=wsum+wgt(i)
                    if (debug_day) &
                    write (debug_unit,'("n,m,nn,i,idaynum,ii,yd,ym_adj,wgt,wsum: ", i6, i3, 2i6, i8, i3, 4f12.6)') &
                        n, m, nn, i, idaynumx(i), ii, ydx(i), ym_adj(nn), wgt(i), wsum
                else
                    nfill = nfill + 1
                    if (debug_day) & 
                    write (debug_unit, '("n,m,nn,i,idaynum,ii,yd,ym_adj: ", i6, i3, 2i6, i8, i3, 2f12.6)'), &
                        n, m, nn, i, idaynumx(i), ii, ydx(i), nfill
                end if
            end do
            if (wsum .ne. 0.0d0 .and. nfill .eq. 0) then
                ym_adj(nn)=ym_adj(nn)/wsum
            else
                ym_adj(nn)=yfill
            end if

            if (debug_day) write (debug_unit, '("n, m, nn, ym(nn), ym_adj(nn), ym_adj-ym, : ", i6, i3, i6, 4f12.5/)') &
                n, m, nn, ym(nn), ym_adj(nn), ym_adj(nn) - ym(nn) !, sngl(ym_adj(nn))
        end do
    end do 
    
    !if (debug_day) write (debug_unit,'(12g14.6)') ym_adj

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

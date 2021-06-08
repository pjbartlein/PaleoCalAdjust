module pseudo_daily_interp_subs
    
implicit none

    logical             :: debug=.false.
    integer             :: debug_unit = 10
    integer             :: out_unit = 6

contains

subroutine hdaily(nm,nd,xm,monlen,no_negatives,xdh)

    implicit none
    
<<<<<<< HEAD
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nd,nm
    real(8), intent(in)     :: xm(nm)
    integer(4), intent(in)  :: monlen(nm)
    logical, intent(in)     :: no_negatives
    real(8), intent(out)    :: xdh(nd)

    real(8)             :: a(0:nh),b(0:nh)
    integer(4)          :: m
    
    if (debug) then
        do m=1,nm
            write (debug_unit,*) m,xm(m),monlen(m)
        end do
    end if

    ! interpolate daily values
    call harmonic_coeffs(nm,xm,a,b)
    call xdhat(nm,nd,monlen,a,b,xdh)
    if (no_negatives) call dzero(nm,nd,monlen,xm,xdh)

end subroutine hdaily

subroutine harmonic_coeffs(nm,y,a,b)
! calculates a's and b's of an "adjusted" harmonic fit to monthly values of a variable,
! which preserves the monthly (and annual) mean values by interpolatred daily values
! adapted from Epstein, E.S. (1991) On obtaining daily climatological values from monthly means
! J. Climate 4:365-368
=======
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
            if (debug_em) write (debug_unit, '("nzero_sum, nzero, nzero_mean: ", f12.6, i4, f12.6)') &
                nonzero_sum, nzero, nonzero_mean
                    
            ! apply porportions of difference to each nonzero subinterval
            if (debug_em) write (debug_unit, '("diff, nzero, yd_adjust: ", f12.6, i3, f12.6)') diff, nzero, yd_adjust
            if (debug_em) write (debug_unit, '(a)') "       i      ii      yd_old   yd_adjust   yd(i)"
            ii = 0; yd_adjust_sum = 0.0d0; yd_out_sum = 0.0d0; yd_adjust_mean = 0.0d0; yd_out_mean = 0.0d0
            do i = ib, ie
                yd_old = yd(i)
                ii = ii + 1
                if (yd(i) .ne. 0.0d0) then
                    yd_adjust = ((yd(i) / nonzero_sum) * diff) * dble(nsubint(n))
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
>>>>>>> parent of 8f0ee60 (Update pseudo_daily_interp_subs.f90)

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm
    real(8), intent(in)     :: y(nm)
    real(8), intent(out)    :: a(0:nh),b(0:nh)
    
    real(8)                 :: pi
    real(8)                 :: asum,bsum,c0,c1,c2,c3,c4
    integer(4)              :: j,t
    
    a=0.0d0; b=0.0d0
    pi=4.0d0*datan(1.0d0)
    
    ! a0
    do t=1,nm
        a(0)=a(0)+y(t)
    end do
    a(0)=a(0)/dble(nm) 
    
    ! a's and b's
    do j=1,nh-1
        if (debug) write (debug_unit,'(a)') " "
        asum=0.0d0; bsum=0.0d0
        c1=pi*(dble(j)/dble(nm))
        do t=1,nm
            c0=dble(t)/dble(nm)
            c2=(2.0d0*pi*dble(j)*dble(t))/dble(nm) 
            asum=asum+y(t)*dcos(c2)/(dble(nh))
            bsum=bsum+y(t)*dsin(c2)/(dble(nh))
            if (debug) write (debug_unit,'("j,t,c0,c2,asum,bsum: ",2i3,5f10.4)') j,t,y(t),c0,c2,asum,bsum
        end do
        a(j)=(c1/dsin(c1))*asum
        b(j)=(c1/dsin(c1))*bsum
        if (debug) write (debug_unit,'("j,c1,a,b:",i3,3f10.4)') j,c1,a(j),b(j)
    end do
    
    if (debug) write (debug_unit,'(a)') " "
    asum=0.0d0
    do t=1,nm
        c3=cos(pi*dble(t))/dble(nm) 
        asum=asum+(y(t)*c3)
        if (debug) write (debug_unit,'("t,y,c3,asum: ",i3,5f10.4)') t,y(t),c3,asum,cos(pi*dble(t)),pi*dble(t)
    end do
    c4=((pi/2.0d0)/sin(pi/2.0d0))
    a(nh)=c4*asum !((pi/2.0)/sin(pi/2.0))*asum
    b(nh)=0.0d0
    if (debug) write (debug_unit,'(a)') " "
    if (debug) write (debug_unit,'("c4,asum: ",2f10.4)') c4,asum
    
    if (debug) then
        write (debug_unit,'(a)') " "
        do j=0,nh
                write (debug_unit,'("j,a,b", i3,2f10.4)') j,a(j),b(j)
        end do
        write (debug_unit,'(a)') " "
    end if 
    
end subroutine harmonic_coeffs

subroutine xdhat(nm,nd,monlen,a,b,yhat)
! calculates/interpolates pseudo-daily values of variable using the a's and b's from harmonic_coeffs()
! adapted from Epstein, E.S. (1991) On obtaining daily climatological values from monthly means
! J. Climate 4:365-368

    implicit none
    
<<<<<<< HEAD
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: a(0:nh),b(0:nh)
    real(8), intent(out)    :: yhat(nd)
    
    integer(4)              :: i,j,m,ii
    real(8)                 :: t,pi
    real(8)                 :: c2
    
    pi=4.0d0*datan(1.0d0)
    yhat=0.0
    ii=0
    do i=1,nm
        do m=1,monlen(i)
            ii=ii+1
            t=(dble(i)-0.5d0)+(dble(m)-0.5d0)/dble(monlen(i))
            do j=0,nh
                c2=((2.0d0*pi*dble(j)*t)/dble(nm))
                yhat(ii)=yhat(ii)+a(j)*dcos(c2)+b(j)*dsin(c2)
=======
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
    
    !integer(4)              :: debug_unit = 14
    !logical                 :: debug_day = .true.

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
>>>>>>> parent of 8f0ee60 (Update pseudo_daily_interp_subs.f90)
            end do
            if (debug) write (debug_unit,'("i,m,ii,t,c2,yhat: ",3i4,3f10.4)') i,m,ii,t,c2,yhat(ii)
        end do
    end do
    
end subroutine xdhat

subroutine dayinterp(nm,nd,monlen,zm,zd)
! interpolate pseudo-daily values of monthly data.  Not mean-preserving.

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer(4)              :: nm1
    integer(4)              :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: zm2(0:nm+1)
    integer                 :: i,m
    
    call midmonth_int(nm,monlen,midmon)
    
    ! pad data at beginning (m=0) and end (m=13)
    nm1=nm+1
    zm2(1:nm)=zm
    zm2(0)=zm(nm)
    zm2(nm1)=zm(1)
    midmon2(1:nm)=midmon
    midmon2(0)=1-(nd-midmon(nm))-2
    midmon2(nm1)=nd+midmon(1)   
    
    ! linear pseudo-daily interpolation
    do i=1,nd
        ! find month day i lies in
        do m=1,nm+1
            if (i.gt.midmon2(m-1) .and. i.le.midmon2(m)) exit
        end do   
        zd(i)=(dble(i-midmon2(m-1))/dble(midmon2(m)-midmon2(m-1)))*(zm2(m)-zm2(m-1))+zm2(m-1)       
    end do
    
end subroutine dayinterp

subroutine dayspread(nm,nd,monlen,zm,zd)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    real(4), intent(in)     :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer                 :: i,j,m
    integer(4)              :: imonlen(nm)

    imonlen = int(monlen)
    
    i=0
    do m=1,nm
        do j=1,imonlen(m)
            i=i+1
            zd(i)=zm(m)
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

subroutine midmonth_real(nm,veq_mon,veq_midmon_day,rmonlen,rmidmon,rmonbeg,rmonend)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    integer(4), intent(in)  :: veq_mon
    real(8), intent(in)     :: veq_midmon_day
    real(8), intent(in)     :: rmonlen(nm)
    real(8), intent(out)    :: rmidmon(nm),rmonbeg(nm),rmonend(nm)
    
    integer(4)              :: m, debug_unit = 10
    
    logical                 :: debug_write = .false.
    
    ! midmonth target values (NOTE:  first month is Jan)
    ! first month is Jan: 0ka equinox:  31.0+28.0+21.5 = 80.5; 0ka March midmonth day:  31.0+28.0+15.5 = 74.5
    ! first month is Dec: 0ka equinox:  31.0+31.0+28.0+21.5 = 111.5; 0ka March midmonth day:  31.0+31.0+28.0+15.5 = 105.5
    
    rmidmon(veq_mon)=veq_midmon_day ! fixed March midmonth relative to noleap VE
    rmonbeg(veq_mon)=rmidmon(veq_mon)-(rmonlen(veq_mon)/2.0d0)
    rmonend(veq_mon)=rmidmon(veq_mon)+(rmonlen(veq_mon)/2.0d0)
    do m=veq_mon-1,1,-1
        rmidmon(m)=rmidmon(m+1)-(rmonlen(m+1)/2.0d0)-(rmonlen(m)/2.0d0) 
        rmonbeg(m)=rmidmon(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmidmon(m)+(rmonlen(m)/2.0d0)
    end do
    do m=veq_mon+1,nm
        rmidmon(m)=rmidmon(m-1)+(rmonlen(m-1)/2.0d0)+(rmonlen(m)/2.0d0) 
        rmonbeg(m)=rmidmon(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmidmon(m)+(rmonlen(m)/2.0d0)
    end do
    
    if (debug_write) then
        write (debug_unit,'("veq_mon, veq_midmon_day: ",i3,f11.6)') veq_mon,veq_midmon_day
        write (debug_unit,'("rmonlen   ",12f11.6)') rmonlen
        write (debug_unit,'("rmidmon   ",12f11.6)') rmidmon
        write (debug_unit,'("rmonbeg   ",12f11.6)') rmonbeg
        write (debug_unit,'("rmonend   ",12f11.6)') rmonend
    end if
                
end subroutine midmonth_real


subroutine dzero(nm,nd,monlen,xm,xd0)
! enforces 0.0 values of interpolated daily data when the monthly mean is 0.0

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xm(nm)
    real(8), intent(inout)  :: xd0(nd)
    
    real(8)                 :: xdm(nm),diff,totaldiff
    integer(4)              :: i,m,j,nonzero(nm),l,maxiter=30
    integer(4)              :: imonlen(nm)
    
    imonlen = int(monlen)
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
    do l=1,maxiter  
      ! zero daily values in months where xm=0.0
      i=0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              if (xm(m).eq.0.0) xd0(i)=0.0          
          end do
      end do
  
      i=0
      xdm=0.0d0; nonzero=0; totaldiff=0.0d0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              xdm(m)=xdm(m)+xd0(i)
              if (xdm(m).gt.0.0d0) nonzero(m)=nonzero(m)+1
          end do
          xdm(m)=xdm(m)/dble(monlen(m))
          diff=dabs((xm(m)-xdm(m)))
          totaldiff=totaldiff+diff
      end do
      if (totaldiff.le.0.0001) exit
    
      i=0
      do m=1,nm
          if (nonzero(m).ne.0) then
              do j=1,imonlen(m)
                  i=i+1
                  if (xd0(i).gt.0d0) then
                      xd0(i)=xd0(i)+(xm(m)-xdm(m)) !/nonzero(m)
                  end if
              end do
          end if
      end do    
    end do    
  
  ! zero daily values in months where xm=0.0
    i=0
    do m=1,nm
        do j=1,imonlen(m)
            i=i+1
            if (xm(m).eq.0.0) xd0(i)=0.0          
        end do
    end do
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
end subroutine dzero

subroutine monmean(nm,nd,monlen,xd,xm)
! gets monthly means of interpolated daily data

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xd(nd)
    real(8), intent(out)    :: xm(nm)
    integer                 :: i,j,m
    
    xm=0.0d0
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            xm(m)=xm(m)+xd(i)
        end do
        xm(m)=xm(m)/dble(monlen(m))
    end do

end subroutine monmean

subroutine ann_wmean(n,x,w,xm)
! gets a weighted (by month length) annual mean value from monthly data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n),w(n)
    real(8), intent(out)    :: xm
    real(8)                 :: wsum
    integer                 :: i
    
    xm=0.0d0; wsum=0.0d0
    do i=1,n
        xm=xm+x(i)*w(i)
        wsum=wsum+w(i)
    end do

    xm=xm/wsum
    
end subroutine ann_wmean
    
subroutine ann_mean(n,x,xm)
! gets an unweighted annual mean value from monthly or daily data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n)
    real(8), intent(out)    :: xm
    integer                 :: i
    
    xm=0.0d0
    do i=1,n
        xm=xm+x(i)
    end do

    xm=xm/dble(n)
    
end subroutine ann_mean   

end module pseudo_daily_interp_subs
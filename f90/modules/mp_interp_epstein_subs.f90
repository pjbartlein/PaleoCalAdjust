module mp_interp_epstein_subs
! Subroutines for implementing the Epstein (1991) "harmonic" mean-preserving interpolation
    
! Epstein, E.S. (1991), On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368

! These subroutines are part of PaleoCalAdjust v1.1:
!   - see P.J. Bartlein & S.L. Shafer (2019) Paleo calendar-effect adjustments in time-slice and transient
!         climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis,
!         Geoscientific Model Development, 12:3889–3913, doi: 10.5194/gmd-12-3889-2019
!   - available from GitHub:  https://github.com/pjbartlein/PaleoCalAdjust or Zenodo:  https://doi.org/10.5281/zenodo.1478824
!   - see also:  https://github.com/pjbartlein/mp-interp (mean-preserving interpolation)
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .false.  
    
contains
    
subroutine mp_interp_epstein(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint, &
    smooth, no_negatives, match_mean, tol, ntargs, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use pseudo_daily_interp_subs

    implicit none
    
    integer(4), intent(in)      :: nctrl, ntargs                    ! number of control and target points
    integer(4), intent(in)      :: n_outer, n_inner                 ! number of outer and inner intervals (e.g. months, years)
    real(8), intent(in)         :: ym(nctrl)                        ! input
    real(8), intent(in)         :: yfill                            ! fill value
    real(8), intent(in)         :: x_ctrl(nctrl)                    ! x_ctrl
    integer(4), intent(in)      :: nsubint(nctrl)                   ! number of subintervals
    logical, intent(in)         :: smooth, no_negatives, match_mean ! logical control variables
    real(8), intent(in)         :: tol                              ! tolerance for enforce_mean()
    integer(4), intent(in)      :: max_nctrl_in                     ! max number of "inner" intervals in an "outer" interval
    integer(4), intent(in)      :: max_ntargs_in                    ! max number of subintervals in an "outer" interval 
    
    real(8), intent(inout)      :: y_int(ntargs)                    ! interpolated values
    real(8), intent(inout)      :: ym_int(nctrl)                    ! mean of interpolated values

    ! local variables

    ! smoothing window variables
    integer(4), parameter       :: nw = 21
    real(8)                     :: pi, x, wgt(nw), wsum
    real(8)                     :: y_int_temp(ntargs)
    
    ! local variables passed to hm_int() subroutine
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_nctrl_in)
    real(8)             :: rmse
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    ! indices
    integer(4)          :: mm, n, ns, i, ii, iii
    integer(4)          :: nn, beg_inner, end_inner, beg_subint, end_subint, n_inner_in, ntargs_out
    integer(4)          :: js, j, jj, jjj
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_epstein()"
    
    pi=4.0d0*datan(1.0d0)

    if (debug_write) write (debug_unit,'("tol:   ", g12.4)') tol

    ! loop over number of years
    nn = 1
    do n = 1, n_outer
    
        beg_inner = (n - 1) * n_inner + 1
        end_inner = beg_inner + (n_inner - 1)
        ntargs_out = sum(nsubint(beg_inner:end_inner))
        beg_subint = nn
        end_subint = beg_subint - 1 + sum(nsubint(beg_inner:end_inner))
        nn = end_subint + 1
        ns = (end_subint - beg_subint + 1)
    
        n_inner_in = end_inner - beg_inner + 1

        if (debug_write) write (debug_unit, '("n,beg_inner,end_inner,: ", 7i5)') & 
            n, beg_inner, end_inner
        if (debug_write) write (debug_unit, '("n_inner_in, ntargs_out, beg_subint, end_subint: ", 2i5, 2i8)') &
            n_inner_in, ntargs_out, beg_subint, end_subint
    
        ! load temporary variables
    
        ym_in = 0.0d0; x_ctrl_in = 0.0d0; nsubint_in = 0; y_int_out = 0.0d0
        ym_in(1:n_inner_in) = ym(beg_inner:end_inner)
        x_ctrl_in(1:n_inner_in) = x_ctrl(beg_inner:end_inner)
        nsubint_in(1:n_inner_in) = nsubint(beg_inner:end_inner)
    
        if (debug_write) write (debug_unit, '(16f8.2)') ym_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16f8.2)') x_ctrl_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16i8)') nsubint_in(1:n_inner_in)
        
        ! if any missing values, set output to yfill, otherwise call hm_int()
        if (any(ym_in(1:n_inner_in) .eq. yfill)) then
        
            if (debug_write) write (debug_unit, '("missing values, interval ", i8)') n
            y_int(beg_subint:end_subint) = yfill
            
        else

            if (debug_write) write (debug_unit, '("call hm_int(), interval ", i8)') n

            call hm_int(n_inner_in, ym_in, yfill, nsubint_in, no_negatives, ntargs_out, y_int_out, ym_int_out)

            if (debug_write) write (debug_unit, '(a)') "back from hm_int()"

            if (debug_write) write (debug_unit,'(10f10.6)') y_int_out
            if (debug_write) write (debug_unit,'(12f10.6)') ym_in(1:n_inner_in)
            if (debug_write) write (debug_unit,'(12f10.6)') ym_int_out(1:n_inner_in)
            if (debug_write) write (debug_unit,'(12f10.6)') ym_in(1:n_inner_in) - ym_int_out(1:n_inner_in)
            if (debug_write) then
                call interp_stat(n_inner_in, ym_in, ym_int_out, rmse)
                write (debug_unit,'("rmse: ", g12.4)') rmse
            end if

            y_int(beg_subint:end_subint) = y_int_out(1:ns) 
        
        end if

    end do
    
    ! smooth across outer intervals
    if (smooth) then

            if (debug_write) then
                iii = 0
                do i = 1, nctrl
                    do ii = 1, nsubint(i)
                        iii = iii + 1
                        write (debug_unit, '("i,ii,iii,y_int(iii): ", 3i6,f12.6)') i,ii,iii,y_int(iii)
                    end do
                end do
            end if
            
            ! generate smoothing weights
            do j=1,nw
                jj=j-((nw-1)/2)-1
                x=(dble(jj)/((nw-1)/2.0d0))*4.0d0
                wgt(j)=(1.0d0/2.0d0*pi)*(exp(-0.5d0*x**2))
                if (debug_write) write (debug_unit,'(2i6,2f12.6)') j,jj,x,wgt(j)
            end do
    
            ! save y_int
            y_int_temp(1:ntargs) = y_int(1:ntargs)
        
            ! smooth across outer intervals
            n=nsubint(1)+nsubint(2)+nsubint(3)+nsubint(4)+nsubint(5)+nsubint(6) &
                +nsubint(7)+nsubint(8)+nsubint(9)+nsubint(10)+nsubint(11)+nsubint(12)
            mm=12
            if (debug_write) write (debug_unit,'("smooth:  n,mm: ",2i6)') n,mm
            do nn=1,n_outer-1
                jjj=n-floor(real(nw)/2.0)-1
                if (debug_write) write (debug_unit,'("nn,jjj: ", 2i6)') nn,jj 
                do js=1,nw
                    jjj=jjj+1
                    wsum=0.0d0; y_int(jjj)=0.0d0
                    jj=jjj-((nw-1)/2)-1
                    if (debug_write) write (debug_unit,'("js,jjj,jj: ", 3i6)') js,jjj,jj
                    do j=1,nw
                        jj=jj+1
                        if (y_int_temp(jj) .ne. yfill) then
                            y_int(jjj)=y_int(jjj)+y_int_temp(jj)*wgt(j)
                            wsum=wsum+wgt(j)
                        if (debug_write)  &
                            write (debug_unit,'("j,jj,jjj,y_int(jjj),y_int_temp(jj),wgt,y_int_temp*wgt,wsum: ", 3i6, 5f12.6)') &
                            j,jj,jjj,y_int(jjj),y_int_temp(jj),wgt(j),y_int_temp(jj)*wgt(j),wsum
                        end if
                    end do
                    if (wsum.ne.0.0d0) then
                        y_int(jjj)=y_int(jjj)/wsum
                    else
                        y_int(jjj)=yfill
                    end if
                    if (debug_write) write (debug_unit,'("y_int_temp(jjj),y_int(jjj),: ", 2f12.6)') y_int_temp(jjj), y_int(jjj)
    
                end do
    
                n=n+nsubint(mm+1)+nsubint(mm+2)+nsubint(mm+3)+nsubint(mm+4)+nsubint(mm+5)+nsubint(mm+6) &
                    +nsubint(mm+7)+nsubint(mm+8)+nsubint(mm+9)+nsubint(mm+10)+nsubint(mm+11)+nsubint(mm+12)
                mm=mm+12
                if (debug_write) write (debug_unit,'("n,mm: ", 2i6)') n,mm
            end do
        
    end if
    
    ! no negatives (also check for NaN's)
    if (no_negatives) then
        do n = 1, ntargs
            if (y_int(n) .lt. 0.0d0) y_int(n) = 0.0d0
            if (isnan(y_int(n))) y_int(n) = 0.0d0
        end do
    end if
        
    if (debug_write) write (debug_unit,'(a)') "y_int_before_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
        
    if (debug_write) write (debug_unit,*) nctrl, ntargs
    if (match_mean) then
        if (debug_write) write (debug_unit, '(a)') "call enforce_mean() "
        call enforce_mean(nctrl, ntargs, nsubint, tol, no_negatives, ym, y_int, yfill)
    end if
    
    call interval_mean(nctrl, nsubint, ntargs, y_int, yfill, ym_int)
    
    if (debug_write) write (debug_unit,'(a)') "y_int_after_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
    if (debug_write) write (debug_unit,'(12f8.3)') ym_int
    
end subroutine mp_interp_epstein   
    
subroutine hm_int(nctrl, ym, ymiss, nsubint, no_negatives, ntargs, y_int, ym_int)
! Epstein 1991 harmonic mean-preserving interpolation
! Note:  this approach is implicitly periodic

    use pseudo_daily_interp_subs
    
    implicit none

    integer(4), intent(in)              :: nctrl, ntargs
    real(8), intent(in)                 :: ym(nctrl)                    ! input 
    real(8), intent(in)                 :: ymiss                        ! missing / fill value
    integer(4), intent(in)              :: nsubint(nctrl)               ! number of subintervals, each main interval
    logical, intent(in)                 :: no_negatives                 ! contrain interpolation to positive
    real(8), intent(out)                :: y_int(ntargs)                ! interpolated values MN[]
    real(8), intent(out)                :: ym_int(nctrl)                ! mean of interpolated values
        
    if (debug_write) write (debug_unit, '(a)') "In hm_int"
    if (debug_write) write (debug_unit, *) nctrl, ntargs
    if (debug_write) write (debug_unit, '(12f9.2)') ym
    if (debug_write) write (debug_unit, '(12i9)') nsubint
    
    call hsubint (nctrl, ntargs, ym, nsubint, no_negatives, y_int)
    
    call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss,  ym_int)
    
    if (debug_write) write (debug_unit,'(12f10.2)') ym
    if (debug_write) write (debug_unit,'(12f10.2)') ym_int
    if (debug_write) write (debug_unit,'(12f10.5)') ym_int - ym 
    
end subroutine hm_int
    
subroutine hsubint(n_inner, ns, ym, nsubint, no_negatives, ydh)

    use pseudo_daily_interp_subs

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: ns,n_inner
    real(8), intent(in)     :: ym(n_inner)
    integer(4), intent(in)  :: nsubint(n_inner)
    logical, intent(in)     :: no_negatives
    real(8), intent(out)    :: ydh(ns)

    real(8)             :: a(0:nh),b(0:nh)
    integer(4)          :: m
    
    if (debug_write) then
        do m=1,n_inner
            write (debug_unit,'(i6,f12.6,i6)') m,ym(m),nsubint(m)
        end do
    end if

    ! interpolate daily values
    call harmonic_coeffs(n_inner,ym,a,b)
    call ydhat(n_inner,ns,nsubint,a,b,ydh)
    if (no_negatives) call dzero(n_inner,ns,nsubint,ym,ydh)

end subroutine hsubint

subroutine harmonic_coeffs(n_inner,y,a,b)

! Calculates a's and b's of an "adjusted" harmonic fit to monthly values of a variable,
! which preserves e.g. the monthly (and annual) mean values by interpolated e.g., daily values.
! Adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: n_inner
    real(8), intent(in)     :: y(n_inner)
    real(8), intent(out)    :: a(0:nh),b(0:nh)
    
    real(8)                 :: pi
    real(8)                 :: asum,bsum,c0,c1,c2,c3,c4
    integer(4)              :: j,t
    
    a=0.0d0; b=0.0d0
    pi=4.0d0*datan(1.0d0)
    
    ! a0
    do t=1,n_inner
        a(0)=a(0)+y(t)
    end do
    a(0)=a(0)/dble(n_inner) 
    
    ! a's and b's
    do j=1,nh-1
        if (debug_write) write (debug_unit,'(a)') " "
        asum=0.0d0; bsum=0.0d0
        c1=pi*(dble(j)/dble(n_inner))
        do t=1,n_inner
            c0=dble(t)/dble(n_inner)
            c2=(2.0d0*pi*dble(j)*dble(t))/dble(n_inner) 
            asum=asum+y(t)*dcos(c2)/(dble(nh))
            bsum=bsum+y(t)*dsin(c2)/(dble(nh))
            if (debug_write) write (debug_unit,'("j,t,c0,c2,asum,bsum: ",2i3,5f10.4)') j,t,y(t),c0,c2,asum,bsum
        end do
        a(j)=(c1/dsin(c1))*asum
        b(j)=(c1/dsin(c1))*bsum
        if (debug_write) write (debug_unit,'("j,c1,a,b:",i3,3f10.4)') j,c1,a(j),b(j)
    end do
    
    if (debug_write) write (debug_unit,'(a)') " "
    asum=0.0d0
    do t=1,n_inner
        c3=cos(pi*dble(t))/dble(n_inner) 
        asum=asum+(y(t)*c3)
        if (debug_write) write (debug_unit,'("t,y,c3,asum: ",i3,5f10.4)') t,y(t),c3,asum,cos(pi*dble(t)),pi*dble(t)
    end do
    c4=((pi/2.0d0)/sin(pi/2.0d0))
    a(nh)=c4*asum !((pi/2.0)/sin(pi/2.0))*asum
    b(nh)=0.0d0
    if (debug_write) write (debug_unit,'(a)') " "
    if (debug_write) write (debug_unit,'("c4,asum: ",2f10.4)') c4,asum
    
    if (debug_write) then
        write (debug_unit,'(a)') " "
        do j=0,nh
                write (debug_unit,'("j,a,b", i3,2f10.4)') j,a(j),b(j)
        end do
        write (debug_unit,'(a)') " "
    end if 
    
end subroutine harmonic_coeffs


subroutine ydhat(n_inner,ns,nsubint,a,b,yhat)

! Calculates a's and b's of an "adjusted" harmonic fit to monthly values of a variable,
! which preserves e.g. the monthly (and annual) mean values by interpolated e.g., daily values.
! Adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: n_inner,ns
    integer(4), intent(in)  :: nsubint(n_inner)
    real(8), intent(in)     :: a(0:nh),b(0:nh)
    real(8), intent(out)    :: yhat(ns)
    
    integer(4)              :: i,j,m,ii
    real(8)                 :: t,pi
    real(8)                 :: c2
    
    pi=4.0d0*datan(1.0d0)
    yhat=0.0
    ii=0
    do i=1,n_inner
        do m=1,nsubint(i)
            ii=ii+1
            t=(dble(i)-0.5d0)+(dble(m)-0.5d0)/dble(nsubint(i))
            do j=0,nh
                c2=((2.0d0*pi*dble(j)*t)/dble(n_inner))
                yhat(ii)=yhat(ii)+a(j)*dcos(c2)+b(j)*dsin(c2)
            end do
            if (debug_write) write (debug_unit,'("i,m,ii,t,c2,yhat: ",3i4,3f12.6)') i,m,ii,t,c2,yhat(ii)
        end do
    end do
    
end subroutine ydhat

end module mp_interp_epstein_subs

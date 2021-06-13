module mp_interp_harzallah_subs
! Subroutines for implementing the Harzallah (1995) "iterative spline" mean-preserving interpolation
    
! Harzallah, A. (1995) The interpolation of data series using a constrained iterating technique
! Monthly Weather Review 123:2251-2254.

! These subroutines are part of PaleoCalAdjust v1.1:
!   - see P.J. Bartlein & S.L. Shafer (2019) Paleo calendar-effect adjustments in time-slice and transient
!         climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis,
!         Geoscientific Model Development, 12:3889–3913, doi: 10.5194/gmd-12-3889-2019
!   - available from GitHub:  https://github.com/pjbartlein/PaleoCalAdjust or Zenodo:  https://doi.org/10.5281/zenodo.1478824
!   - see also:  https://github.com/pjbartlein/mp-interp (mean-preserving interpolation)

! This implementation provides two versions of spline fitting:
! 1 = J. Burkhardt implementation of cubic spline (in spline.f90,
!   https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html), last accessed 8 Dec 2020)
! 2 = J. Burkhardt implementation of piecewise cubic Hermite spline (in spline.f90,
!   https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html), last accessed 8 Dec 2020)
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .false.  
    
contains

subroutine mp_interp_harzallah(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint, &
    spline_case, npad, no_negatives, match_mean, tol, ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use pseudo_daily_interp_subs

    implicit none
    
    integer(4), intent(in)      :: nctrl, ntargs                    ! number of control and target points
    integer(4), intent(in)      :: n_outer, n_inner                 ! number of outer and inner intervals (e.g. nyears, nmonths)
    real(8), intent(in)         :: ym(nctrl)                        ! input
    real(8), intent(in)         :: yfill                            ! fill value
    real(8), intent(in)         :: x_ctrl(nctrl)                    ! x_ctrl 
    integer(4), intent(in)      :: nsubint(nctrl)                   ! number of subintervals
    integer(4), intent(in)      :: spline_case                      ! selects spline-fitting procedure
    integer(4), intent(in)      :: npad                             ! number of padding each end
    logical, intent(in)         :: no_negatives, match_mean         ! logical control variables
    real(8), intent(in)         :: tol                              ! tolerance for enforce_mean()
    real(8), intent(in)         :: x_targ(ntargs)                   ! x_targ
    integer(4), intent(in)      :: max_nctrl_in                     ! max number of "inner" intervals in an "outer" interval
    integer(4), intent(in)      :: max_ntargs_in                    ! max number of subintervals in an "outer" interval 
    
    real(8), intent(inout)        :: y_int(ntargs)                  ! interpolated values
    real(8), intent(inout)        :: ym_int(nctrl)                  ! mean of interpolated values
                                                                    ! inner interval = e.g. months, outer interval = e.g. year
    ! local variables

    ! variables passed to hz_int() subroutine
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: x_targ_in(max_ntargs_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_ntargs_in)
    real(8)             :: rmse
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    ! indices
    integer(4)          :: n, nn
    integer(4)          :: beg_inner, end_inner, beg_ctrl, end_ctrl, nctrl_in, beg_targ, end_targ, ntarg_in, beg_int, end_int, nint_out
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_harzallah()"
    if (debug_write) write (debug_unit, *) n_outer, n_inner, nctrl, ntargs, max_nctrl_in, max_ntargs_in

    ! loop over number of years
    nn = 1 ! starting pseudo-daily value
    do n = 1, n_outer
        
        ! range of inner intervals (e.g. months)
        beg_inner = (n - 1) * n_inner + 1 
        end_inner = n * n_inner 
        
        ! range of control values (padded)
        beg_ctrl = beg_inner - npad
        end_ctrl = end_inner + npad

        if (n .eq. 1) beg_ctrl = 1
        if (n .eq. n_outer) end_ctrl = n * n_inner
    
        nctrl_in = end_ctrl - beg_ctrl + 1
        
        ! range of target values, including padding
        beg_targ = sum(nsubint(1:(beg_ctrl - 1))) + 1
        end_targ = beg_targ + sum(nsubint(beg_ctrl:end_ctrl)) - 1
        
        ntarg_in = end_targ - beg_targ + 1
        
        ! range of interpolated values to save
        beg_int = nn 
        end_int = nn + sum(nsubint(beg_inner:end_inner)) - 1
        
        nint_out = end_int - beg_int + 1
        
        ! update nn for next year
        nn = end_int + 1
        !write (*, '(5i6)') n, beg_int, end_int, nint_out, nn
        
        if (debug_write) write (debug_unit, '(a)') " "
        if (debug_write) write (debug_unit, '("n,beg_inner,end_inner,n_inner): ",4i9)') &
            n,beg_inner, end_inner, end_inner - beg_inner + 1
        if (debug_write) write (debug_unit, '("  beg_ctrl,end_ctrl,nctrl_in,x_ctrl(beg_ctrl),x_ctrl(end_ctrl): ",3i9,2g14.6)') &
                beg_ctrl,end_ctrl,nctrl_in,x_ctrl(beg_ctrl),x_ctrl(end_ctrl)
        if (debug_write) write (debug_unit, '("  beg_targ,end_targ,ntarg_in,x_targ(beg_targ),x_targ(end_targ): ",3i9,2g14.6)') &
                beg_targ,end_targ,ntarg_in,x_targ(beg_targ),x_targ(end_targ)
        if (debug_write) write (debug_unit, '("  beg_int,end_int,nint_in,x_targ(beg_int),x_targ(end_int): ",3i9,2g14.6)') &
                beg_int,end_int,nint_out,x_targ(beg_int),x_targ(end_int)
    
        ! load temporary variables
        
        ym_in = 0.0d0; x_ctrl_in = 0.0d0; nsubint_in = 0; y_int_out = 0.0

        ym_in(1:nctrl_in) = ym(beg_ctrl:end_ctrl)
        x_ctrl_in(1:nctrl_in) = x_ctrl(beg_ctrl:end_ctrl)
        nsubint_in(1:nctrl_in) = nsubint(beg_ctrl:end_ctrl)   
        x_targ_in(1:ntarg_in) = x_targ(beg_targ:end_targ)
        
        if (debug_write) write (debug_unit, '(a)') " ym_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16g14.6)') ym_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " x_ctrl_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16g14.6)') x_ctrl_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " nsubint_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16i16)') nsubint_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " x_targ_in(1), x_targ_in(ntarg_in)"
        if (debug_write) write (debug_unit, '(16g14.6)') x_targ_in(1), x_targ_in(ntarg_in)
        
        ! if any missing values, set output to yfill, otherwise call hz_int()
        if (any(ym_in(1:nctrl_in) .eq. yfill)) then
        
            if (debug_write) write (debug_unit, '("missing values, interval ", i8)') n
            y_int(beg_targ:end_targ) = yfill
            
        else
        
            if (debug_write) write (debug_unit, '("call hz_int(), interval ", i8)') n
            
            call hz_int(spline_case, nctrl_in, ym_in, yfill, x_ctrl_in, nsubint_in, ntarg_in, x_targ_in, y_int_out, ym_int_out)
            
            if (debug_write) write (debug_unit, '(a)') "back from hz_int()"
        
            if (debug_write) write (debug_unit,'(10g14.6)') y_int_out
            if (debug_write) write (debug_unit,'(16g14.6)') ym_in(1:nctrl_in)
            if (debug_write) write (debug_unit,'(16g14.6)') ym_int_out(1:nctrl_in)
            if (debug_write) write (debug_unit,'(16g14.6)') ym_in(1:nctrl_in) - ym_int_out(1:nctrl_in)
            if (debug_write) then
                call interp_stat(nctrl_in, ym_in, ym_int_out, rmse)
                write (debug_unit,'("rmse: ", g12.4)') rmse
            end if

            y_int(beg_int:end_int) = y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))

            if (debug_write) then
                write (debug_unit, '(a)') "n,beg_targ, beg_int, end_int, nint_out, beg_int - beg_targ + 1, beg_int - beg_targ + nint_out"
                write (debug_unit, '(8i9)') n,beg_targ, beg_int, end_int, nint_out, beg_int - beg_targ + 1, beg_int - beg_targ + nint_out
                write (debug_unit, '(a)') "y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))"
                write (debug_unit,'(10g14.6)') y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))
                write (debug_unit, '(a)') "y_int(beg_int:end_int)"
                write (debug_unit,'(10g14.6)') y_int(beg_int:end_int)
            end if
                
        end if

    end do
    
    ! no negatives (also check for NaN's)
    if (no_negatives) then
        do n = 1, ntargs
            if (y_int(n) .lt. 0.0d0) y_int(n) = 0.0d0
            if (isnan(y_int(n))) y_int(n) = 0.0d0
        end do
    end if
        
    if (debug_write) write (debug_unit,'(a)') "y_int_before_enforce_mean"
    if (debug_write) write (debug_unit,'(10g14.6)') y_int
        
    if (debug_write) write (debug_unit,*) nctrl, ntargs
    if (match_mean) then
        if (debug_write) write (debug_unit, '(a)') "call enforce_mean() "
        call enforce_mean(nctrl, ntargs, nsubint, tol, no_negatives, ym, y_int, yfill)
    end if
    
    if (debug_write) write (debug_unit,'(a)') "y_int_after_enforce_mean"
    if (debug_write) write (debug_unit,'(10g14.6)') y_int
    if (debug_write) write (debug_unit,'(12g14.6)') ym_int

    call interval_mean(nctrl, nsubint, ntargs, y_int, yfill, ym_int)   
     
end subroutine mp_interp_harzallah
    
subroutine hz_int(spline_case, nctrl, ym, ymiss, x_ctrl, nsubint, ntargs, x_targ, y_int, ym_int)
! manages Harzalla spline pseudo-daily interpolation

! This implementation provides two versions of spline fitting:
! 1 = J. Burkhardt implementation of cubic spline (in spline.f90,
!   https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html), last accessed 8 Dec 2020)
! 2 = J. Burkhardt implementation of piecewise cubic Hermite spline (in spline.f90,
!   https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html), last accessed 8 Dec 2020)

    use pseudo_daily_interp_subs
    use spline_subs
    
    implicit none
    
    integer(4), intent(in)          :: nctrl, ntargs, spline_case
    integer(4), intent(in)          :: nsubint(nctrl)                   ! interval lengths
    real(8), intent(in)             :: ym(nctrl)                        ! input y
    real(8), intent(in)             :: ymiss                            ! missing / fill value
    real(8), intent(in)             :: x_ctrl(nctrl)                    ! input x
    real(8), intent(in)             :: x_targ(ntargs)                   ! target values
    real(8), intent(out)            :: y_int(ntargs)                    ! interpolated values
    real(8), intent(out)            :: ym_int(nctrl)                    ! means of interpolated values
    
    ! local variables
    integer(4)                      :: ibcbeg = 0, ibcend = 0           ! boundary condition flags
    real(8)                         :: ybcbeg = 0.0d0, ybcend = 0.0d0   ! boundary conditions
    real(8)                         :: ypp(nctrl)                       ! derivatives
    real(8)                         :: ypval, yppval                    ! first and second derivatives
    
    integer(4), parameter           :: max_iter = 16
    real(8)                         :: tol = 0.01                       ! tolerance
    real(8)                         :: resid(nctrl)
    real(8)                         :: ym_mat(max_iter + 1, nctrl)       ! matrix
    real(8)                         :: y_int_mat(max_iter, ntargs)
    
    integer(4)                      :: i, k

    if (debug_write) write (debug_unit, '(a)') "In hz_interp ============="
    if (debug_write) write (debug_unit, '(a)') "nctrl, ntargs, spline_case"
    if (debug_write) write (debug_unit, *) nctrl, ntargs, spline_case
    if (debug_write) write (debug_unit, '(a)') "ym"
    if (debug_write) write (debug_unit, '(16g14.6)') ym
    if (debug_write) write (debug_unit, '(a)') "x_ctrl"
    if (debug_write) write (debug_unit, '(16g14.6)') x_ctrl
    if (debug_write) write (debug_unit, '(a)') "x_targ(1), x_targ(ntargs)"
    if (debug_write) write (debug_unit, '(16g14.6)') x_targ(1), x_targ(ntargs)
    
    ym_mat = 0.0d0; y_int_mat = 0.0d0; ypp = 0.0d0
    ym_mat(1, :) = ym(:)
    
    do k = 1, max_iter

        if (debug_write) write (debug_unit, '("ym_mat", i4)') k
        if (debug_write) write (debug_unit,'(16g14.6)') ym_mat(k, :)
        if (debug_write) write (debug_unit,*) spline_case
        
        select case(spline_case) ! J. Burkhardt spline interpolation methods:
            case (1)
                call spline_cubic_set ( nctrl, x_ctrl, ym_mat(k, :), ibcbeg, ybcbeg, ibcend, ybcend, ypp )
                if (debug_write) write (debug_unit,'(16g14.6)') ypp
                do i = 1, ntargs
                    call spline_cubic_val ( nctrl, x_ctrl, ym_mat(k, :), ypp, x_targ(i), y_int_mat(k,i), ypval, yppval )
                    !if (debug_write) write (debug_unit,'(i4, 4g14.6)') i, x_targ(i), y_int_mat(k,i) !, ypval, yppval
                end do     
            case (2)
                call spline_pchip_set ( nctrl, x_ctrl, ym_mat(k, :), ypp )
                if (debug_write) write (debug_unit,'(a)') "ypp"
                if (debug_write) write (debug_unit,'(16g14.6)') ypp      
                call spline_pchip_val ( nctrl, x_ctrl, ym_mat(k, :), ypp, ntargs, x_targ, y_int_mat(k, :) )
            case default
                stop "spline_case"
        end select
    
        
        ! interval mean values
        if (debug_write) write (debug_unit, '("calling interval_mean() from hz_int(), k, nctrl, ntargs: ", 3i8)') k, nctrl, ntargs
        call interval_mean(nctrl, nsubint, ntargs, y_int_mat(k, :), ymiss, ym_int)
        if (debug_write) write (debug_unit, '("ym_int", i4)') k
        if (debug_write) write (debug_unit,'(16g14.6)') ym_int
     
        ! residuals
        resid = ym_mat(k, :) - ym_int
        if (debug_write) write (debug_unit, '("residuals", i4)') k
        if (debug_write) write (debug_unit,'(16g14.6)') resid
              
        if (debug_write) write (debug_unit,'("max resid: ", g16.9)') maxval(resid)
        if (maxval(resid) .lt. tol) exit
        
        ym_mat(k + 1, :) = resid(:)
     
    end do
    
    if (debug_write) write (debug_unit,'("k (out): ", i4)') k
    
    ! sum over iterations
    y_int = 0.0d0
    do k = 1, max_iter
        do i = 1, ntargs
            y_int(i) = y_int(i) + y_int_mat(k, i)
        end do
    end do
    
    if (debug_write) then
        do i = 1, ntargs
            write (debug_unit, '(i4, 18g14.6)') i, x_targ(i), y_int(i), (y_int_mat(k, i), k = 1, max_iter)
        end do
    end if
    
    if (debug_write) write (debug_unit, '("calling interval_mean() from hz_int() end, nctrl, ntargs: ", 3i8)') nctrl, ntargs    
    call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_int)

    if (debug_write) write (debug_unit,'(16g14.6)') ym
    if (debug_write) write (debug_unit,'(16g14.6)') ym_int
    if (debug_write) write (debug_unit,'(16g14.6)') ym_int - ym
    
end subroutine hz_int

end module mp_interp_harzallah_subs

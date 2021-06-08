module spline_subs
! This module contains several suroutines and functions from John Burkhardt's library of Fortran90 software
! at https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html  
! The subroutines below were extracted from the spline.f90 source code file at:
! https://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html), last accessed 8 Dec 2020)
    
! This implementation provides two versions of spline fitting:
! 1 = cubic spline (in spline.f90, and
! 2 = J piecewise cubic Hermite spline 
    
implicit none

contains

subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )^2
!             + D(IVAL) * ( T - T(IVAL) )^3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )^2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))^2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!    0: the spline should be a quadratic over the first interval;
!    1: the first derivative at the left endpoint should be YBCBEG;
!    2: the second derivative at the left endpoint should be YBCBEG;
!    3: Not-a-knot: the third derivative is continuous at T(2).
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!    0: the spline should be a quadratic over the last interval;
!    1: the first derivative at the right endpoint should be YBCEND;
!    2: the second derivative at the right endpoint should be YBCEND;
!    3: Not-a-knot: the third derivative is continuous at T(N-1).
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of
!    the cubic spline.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ybcbeg
  real ( kind = 8 ) ybcend
  real ( kind = 8 ) ypp(n)
!
!  Check.
!
  !write (*,'(12f10.5)') t
  !write (*,'(12f10.5)') y
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n - 1
    if ( t(i+1) <= t(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop 1
    end if
  end do
!
!  Zero out the matrix.
!
  a1(1:n) = 0.0D+00
  a2(1:n) = 0.0D+00
  a3(1:n) = 0.0D+00
  a4(1:n) = 0.0D+00
  a5(1:n) = 0.0D+00
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    b(1) = 0.0D+00
    a3(1) =  1.0D+00
    a4(1) = -1.0D+00
  else if ( ibcbeg == 1 ) then
    b(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    a3(1) = ( t(2) - t(1) ) / 3.0D+00
    a4(1) = ( t(2) - t(1) ) / 6.0D+00
  else if ( ibcbeg == 2 ) then
    b(1) = ybcbeg
    a3(1) = 1.0D+00
    a4(1) = 0.0D+00
  else if ( ibcbeg == 3 ) then
    b(1) = 0.0D+00
    a3(1) = - ( t(3) - t(2) )
    a4(1) =   ( t(3)        - t(1) )
    a5(1) = - (        t(2) - t(1) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
    stop 1
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n - 1
    b(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
         - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    a2(i) = ( t(i+1) - t(i)   ) / 6.0D+00
    a3(i) = ( t(i+1) - t(i-1) ) / 3.0D+00
    a4(i) = ( t(i)   - t(i-1) ) / 6.0D+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    b(n) = 0.0D+00
    a2(n) = -1.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 1 ) then
    b(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    a2(n) = ( t(n) - t(n-1) ) / 6.0D+00
    a3(n) = ( t(n) - t(n-1) ) / 3.0D+00
  else if ( ibcend == 2 ) then
    b(n) = ybcend
    a2(n) = 0.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 3 ) then
    b(n) = 0.0D+00
    a1(n) = - ( t(n) - t(n-1) )
    a2(n) =   ( t(n)          - t(n-2) )
    a3(n) = - (        t(n-1) - t(n-2) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
    stop 1
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0D+00
    ypp(2) = 0.0D+00
!
!  Solve the linear system.
!
  else

    call penta ( n, a1, a2, a3, a4, a5, b, ypp )

  end if

  return
end

subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )^2
!             + D * ( T - T(IVAL) )^3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) h
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tval
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ypp(n)
  real ( kind = 8 ) yppval
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call r8vec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( 0.5D+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end

subroutine spline_pchip_set ( n, x, f, d )

!*****************************************************************************80
!
!! SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    This routine computes what would normally be called a Hermite
!    interpolant.  However, the user is only required to supply function
!    values, not derivative values as well.  This routine computes
!    "suitable" derivative values, so that the resulting Hermite interpolant
!    has desirable shape and monotonicity properties.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by SPLINE_PCHIP_VAL.
!
!    This routine was originally named "PCHIM".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), dependent variable values to be
!    interpolated.  F(I) is the value corresponding to X(I).
!    This routine is designed for monotonic data, but it will work for any
!    F array.  It will force extrema at points where monotonicity switches
!    direction.
!
!    Output, real ( kind = 8 ) D(N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dmin
  real ( kind = 8 ) drat1
  real ( kind = 8 ) drat2
  real ( kind = 8 ) dsave
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) hsum
  real ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  !real ( kind = 8 ) pchst
  real ( kind = 8 ) temp
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    stop 1
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      stop 1
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(2) - f(1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1) = del1
    d(n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(3) - f(2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1), del1 ) <= 0.0D+00 ) then

    d(1) = 0.0D+00
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then

     dmax = 3.0D+00 * del1

     if ( abs ( dmax ) < abs ( d(1) ) ) then
       d(1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(i+1) - f(i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(i) = 0.0D+00

    temp = pchst ( del1, del2 )

    if ( temp < 0.0D+00 ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0D+00 ) then

      if ( del2 /= 0.0D+00 ) then
        if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0D+00 * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(n) = w1 * del1 + w2 * del2

  if ( pchst ( d(n), del2 ) <= 0.0D+00 ) then
    d(n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del2

    if ( abs ( dmax ) < abs ( d(n) ) ) then
      d(n) = dmax
    end if

  end if

  return
end

subroutine spline_pchip_val ( n, x, f, d, ne, xe, fe )

!*****************************************************************************80
!
!! SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.
!
!  Description:
!
!    This routine may be used by itself for Hermite interpolation, or as an
!    evaluator for SPLINE_PCHIP_SET.
!
!    This routine evaluates the cubic Hermite function at the points XE.
!
!    Most of the coding between the call to CHFEV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFEV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!    This routine was originally called "PCHFE".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.
!
!    Input, real ( kind = 8 ) D(N), the derivative values.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    stop 1
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      stop 1
    end if
  end do

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
!
!  Loop over intervals.
!  The interval index is IL = IR-1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of the loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in the interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfev ( x(ir-1), x(ir), f(ir-1), f(ir), d(ir-1), d(ir), &
        nj, xe(j_first:j-1), fe(j_first:j-1), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFEV.'
        stop 1
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          stop 1
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
        else

          j_new = -1

          do i = j_first, j - 1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            stop 1
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end

subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr )

!*****************************************************************************80
!
!! CHFEV evaluates a cubic polynomial given in Hermite form.
!
!  Discussion:
!
!    This routine evaluates a cubic polynomial given in Hermite form at an
!    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!    The cubic polynomial is determined by function values
!    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2008
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at X1 and
!    X2, respectively.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the function is to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in NEXT.
!
!    Output, real ( kind = 8 ) FE(NE), the value of the cubic function
!    at the points XE.
!
!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of
!    extrapolation points:
!    NEXT(1) = number of evaluation points to the left of interval.
!    NEXT(2) = number of evaluation points to the right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none

  integer ( kind = 4 ) ne

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) delta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fe(ne)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xe(ne)
  real ( kind = 8 ) xma
  real ( kind = 8 ) xmi

  if ( ne < 1 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points is less than 1.'
    write ( *, '(a,i8)' ) '  NE = ', ne
    stop 1
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
    stop 1
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h
  c2 = -( del1 + del1 + del2 )
  c3 = ( del1 + del2 ) / h
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
!
!  Count the extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end

subroutine penta ( n, a1, a2, a3, a4, a5, b, x )

!*****************************************************************************80
!
!! PENTA solves a pentadiagonal system of linear equations.
!
!  Discussion:
!
!    The matrix A is pentadiagonal.  It is entirely zero, except for
!    the main diagaonal, and the two immediate sub- and super-diagonals.
!
!    The entries of Row I are stored as:
!
!      A(I,I-2) -> A1(I)
!      A(I,I-1) -> A2(I)
!      A(I,I)   -> A3(I)
!      A(I,I+1) -> A4(I)
!      A(I,I-2) -> A5(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cheney, Kincaid,
!    Numerical Mathematics and Computing,
!    1985, pages 233-236.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), A4(N), A5(N), the nonzero
!    elements of the matrix.  Note that the data in A2, A3 and A4
!    is overwritten by this routine during the solution process.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult

  do i = 2, n - 1
    xmult = a2(i) / a3(i-1)
    a3(i) = a3(i) - xmult * a4(i-1)
    a4(i) = a4(i) - xmult * a5(i-1)
    b(i) = b(i) - xmult * b(i-1)
    xmult = a1(i+1) / a3(i-1)
    a2(i+1) = a2(i+1) - xmult * a4(i-1)
    a3(i+1) = a3(i+1) - xmult * a5(i-1)
    b(i+1) = b(i+1) - xmult * b(i-1)
  end do

  xmult = a2(n) / a3(n-1)
  a3(n) = a3(n) - xmult * a4(n-1)
  x(n) = ( b(n) - xmult * b(n-1) ) / a3(n)
  x(n-1) = ( b(n-1) - a4(n-1) * x(n) ) / a3(n-1)
  do i = n - 2, 1, -1
    x(i) = ( b(i) - a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
  end do

  return
end

function pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2008
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG1, ARG2, two values to check.
!
!    Output, real ( kind = 8 ) PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) pchst

  pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

  if ( arg1 == 0.0D+00 .or. arg2 == 0.0D+00 ) then
    pchst = 0.0D+00
  end if

  return
end

subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end

end module spline_subs
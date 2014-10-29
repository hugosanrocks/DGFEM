subroutine JacobiP( m, n, alpha, beta, x, p )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! connect this subroutine to the one for DGFEM: more arguments as x,p should be provided
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!*****************************************************************************80
!
!! J_POLYNOMIAL evaluates the Jacobi polynomials J(n,a,b,x).
!
!  Differential equation:
!
!    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
!
!  Recursion:
!
!    P(0,ALPHA,BETA,X) = 1,
!
!    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
!
!    P(N,ALPHA,BETA,X)  =
!      (
!        (2*N+ALPHA+BETA-1)
!        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
!        * P(N-1,ALPHA,BETA,X)
!        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
!      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
!
!  Restrictions:
!
!    -1 < ALPHA
!    -1 < BETA
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
!      * P(N,ALPHA,BETA,X)^2 dX
!    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
!      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
!
!  Special values:
!
!    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) ALPHA, one of the parameters defining the Jacobi
!    polynomials, ALPHA must be greater than -1.
!
!    Input, real ( kind = 8 ) BETA, the second parameter defining the Jacobi
!    polynomials, BETA must be greater than -1.
!
!    Input, real ( kind = 8 ) X(M), the point at which the polynomials are
!    to be evaluated.
!
!    Output, real ( kind = 8 ) P(M), the values of the first N Jacobi
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: alpha
  real ( kind = 8 ) :: beta
  real ( kind = 8 ),dimension(:,:), allocatable :: cx
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c4
  integer ( kind = 4 ) i
  real ( kind = 8 ) r_i
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) p(m)

  real ( kind = 8 ) gamma0, gamma1
  integer ( kind=4 ) id,jd

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JacobiP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JacobiP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    write ( *, '(a)' ) '  But BETA must be greater than -1.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

id=n+1
allocate(cx(1:m,0:id))   ! easier to keep track of zero order to n order ... avoid the +1 shift

  !Polynomials P0
  gamma0=((2.00+00**(alpha+beta+1.00+00))/(alpha+beta+1.00+00))*gamma(alpha+1.00+00)*gamma(beta+1.00+00)/gamma(alpha+beta+1.00+00)
  
  cx(1:m,0) = 1.0D+00/sqrt(gamma0)
  if ( n == 0 ) then
    p(:)=cx(:,n)
    return
  end if

  !Polynomials P1
  gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
  cx(1:m,1) = ( 0.5D+00 * ( alpha + beta + 2.00+00 ) ) * x(1:m) &
    + 0.5D+00 * ( alpha - beta )
  cx(1:m,1) = cx(1:m,1)/sqrt(gamma1)

   c4 = (2.0 / (2.0 + alpha +beta))*sqrt((alpha + 1.00+00) * (beta + 1.00+00) / (alpha + beta + 3.0))

   id=n-1
   do i=1,id
   r_i = real ( i, kind = 8 )

   c1 = 2.0 * r_i + alpha + beta
   c2 = 2.0 / (c1 + 2.0) * sqrt( (r_i + 1.00+00) &
   * (r_i + 1.00+00 + alpha + beta) * (r_i + 1.00+00 + alpha) &
   * (r_i + 1.00+00 + beta) / (c1 + 1.00+00) / (c1 + 3.0))
   c3 = (alpha**2.0 - beta**2.0 ) / c1 / (-c1 - 2.0)
   cx(:,i+1) = 1.0 / c2 *( -1.0 * c4 * cx(:,i-1) + (x(:) - c3) * cx(:,i))
   c4 = c2
   enddo
  p(:)=cx(:,n)    ! only save the last one : N order ...

  return
end subroutine JacobiP

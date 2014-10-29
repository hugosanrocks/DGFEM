subroutine JacobiGQ( n, alpha, beta, x, w )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! connect this subroutine to the one for DGFEM: more arguments as x,w should be provided
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!*****************************************************************************80
!
!! J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
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
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ), N, the order.
!
!    Input, real ( kind = 8 ), ALPHA, BETA, the parameters.
!    -1 < ALPHA, BETA.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) abi
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
  ab = alpha + beta
  abi = 2.0D+00 + ab
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00 ** ( ab + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
    * gamma ( beta + 1.0D+00 ) / gamma ( abi )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / abi
  x(2:n) = 0.0D+00
  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( abi + 1.0D+00 ) * abi * abi )
  bj(2:n) = 0.0D+00

  a2b2 = beta * beta - alpha * alpha

  do i = 2, n
    i_r8 = real ( i, kind = 8 )
    abi = 2.0D+00 * i_r8 + ab
    x(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
    abi = abi ** 2
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + ab ) / ( ( abi - 1.0D+00 ) * abi )
  end do
  bj(1:n) =  sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end subroutine JacobiGQ

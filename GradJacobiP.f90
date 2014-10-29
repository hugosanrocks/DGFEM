subroutine GradJacobiP( m, n, alpha, beta, r, dp )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! connect this subroutine to the one for DGFEM: more arguments as x,dp should be provided
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!*****************************************************************************
! purpose: Evaluate teh derivative of the Jacobi polynomial of type 
!  (alpha,beta)>-1, at points r for order N and returns dP[1:length(r))]
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
!    Input, real ( kind = 8 ) r(M), the point at which the polynomials are
!    to be evaluated.
!
!    Output, real ( kind = 8 ) dP(M), the values of the derivative at the point r.
!
  implicit none

  integer ( kind = 4 ) :: m, i, nm1
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: alpha, alphaP1
  real ( kind = 8 ) :: beta, betaP1

  real ( kind = 8 ) :: r(m), pol
  real ( kind = 8 ) :: dp(m)


  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GradJacobiP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GradJacobiP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    write ( *, '(a)' ) '  But BETA must be greater than -1.'
    stop
  end if

  if ( n < 0 ) then
    write(*,*) ' GradJacobiP n should be positive',n
    return
  end if

  if ( n == 0 ) then
    dp(:)=0.d00
    return
  else
    nm1=n-1
    alphaP1=alpha+1.d00
    betaP1=beta+1.d00

    do i=1,m
    call JacobiP(1,nm1,alphaP1,betaP1,r(i),pol)
    dp(i)=sqrt(n*(n+alpha+beta+1))*pol/2.d00
    enddo

  end if

  return
end subroutine GradJacobiP

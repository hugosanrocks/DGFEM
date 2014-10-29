subroutine GradVandermonde1D(m, n, r, dVr )

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! connect this subroutine to the one for DGFEM: more arguments as r,dVr should be provided
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
  !    Input, real ( kind = 8 ) r(M), the point at which the polynomials are
  !    to be evaluated.
  !
  !    Output, real ( kind = 8 ) dP(M), the values of the derivative
  !    at the point r.
  !
  implicit none

  integer ( kind = 4 ) :: m
  integer ( kind = 4 ) :: n, i, j

  real ( kind = 8 ) :: r(m)
  real ( kind = 8 ) :: dVr(m,n+1)

  real ( kind = 8 ) :: dP(m)

  do i=0,n
     call GradJacobiP( m, i, 0.d00, 0.d00, r, dP )
                                    ! on doit faire plus rapide
     do j=1,m
     dVr(j,i+1)=dp(j)
     enddo
  enddo

  return
end subroutine GradVandermonde1D

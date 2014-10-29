subroutine Dmatrix1D(m, n, r, V, Dr )

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! connect this subroutine to the one for DGFEM: more arguments as r,dVr should be provided
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !*****************************************************************************
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
  !
  !    Input, real ( kind = 8 ) r(M), the point at which the polynomials are
  !    to be evaluated.
  !
  !    Input, real ( kind = 8 ) V(M,N+1), Vandermonde values 
  !
  !    Output, real ( kind = 8 ) dr(M,N+1), the values of the derivative
  !    at the point r.
  !
  implicit none

  integer ( kind = 4 ) :: m, i, j, k
  integer ( kind = 4 ) :: n

  real ( kind = 8 ) :: r(m), c(m,n+1), v1d(m,n+1)
  real ( kind = 8 ) :: V(m,n+1), Vr(m,n+1), Dr(m,n+1)

  call GradVandermonde1D(m, n, r, Vr)
  v1d=v

  call inverse(v1d,c,n+1)

  DO i=1,m
   DO j=1,n+1
   Dr(i,j) = 0.0
    DO k=1,n+1
    Dr(i,j) = Dr(i,j) + vr(i,k)*c(k,j)
    ENDDO
   ENDDO
  ENDDO

  return
end subroutine Dmatrix1D

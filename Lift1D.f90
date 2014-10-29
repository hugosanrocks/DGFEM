subroutine Lift1D(m, n, nfp, V, L )

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! connect this subroutine to the one for DGFEM: more arguments as r,dVr should be provided
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !*****************************************************************************
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the number of nodes per face.
  !
  !    Input, real ( kind = 8 ) V(M,N+1), Vandermonde matrix NodesperElement X N_order+1.
  !
  !    Input, real ( kind = 8 ) L(M,N+1), Lift1D matrix.
  !
  implicit none

  integer ( kind = 4 ) :: m, i, j, k
  integer ( kind = 4 ) :: n, nfp

  real ( kind = 8 ) :: VT(n+1,m), L1(m,nfp)
  real ( kind = 8 ) :: V(m,n+1), L(m,nfp), Emat(m,nfp)

  VT=transpose(V)

  do i=1,m
  Emat(i,:)=0.d00
  enddo
  
  j=1
  Emat(j,j)=1.d00
  Emat(m,nfp)=1.d00

  DO i=1,m
   DO j=1,nfp
   L1(i,j) = 0.0
    DO k=1,n+1
    L1(i,j) = L1(i,j) + VT(i,k)*Emat(k,j)
    ENDDO
   ENDDO
  ENDDO

  DO i=1,m
   DO j=1,nfp
   L(i,j) = 0.0
    DO k=1,n+1
    L(i,j) = L(i,j) + V(i,k)*L1(k,j)
    ENDDO
   ENDDO
  ENDDO


  return
end subroutine Lift1D


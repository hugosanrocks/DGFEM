subroutine Vandermonde1D( m, n, r, v1D )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! connect this subroutine to the one for DGFEM: vector x should be provided
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Parameters:
!
!    Input, integer ( kind = 4 ) m the number of zeroes
!    Input, integer ( kind = 4 ) n the order
!
!    Output, real ( kind = 8 ) V1D(m,n+1) Vandermonde values
!
! Normalement la matrice doit etre carre m=n+1 ? Jean
!=======================================================================================
  implicit none

  integer ( kind = 4 ) :: n, m, cont, cont2, im
  real ( kind = 8 ) :: r(m), alpha, beta, pol
  real ( kind = 8 ) :: v1D(m,n+1)
  
  alpha=0.0+00
  beta=0.0+00

    do cont=1,m
     im=0
     do cont2=1,n+1
     call JacobiP(1,im,alpha,beta,r(cont),pol)
     v1d(cont,cont2)=pol
     im=im+1
     enddo
    enddo
 return
 end subroutine Vandermonde1D

subroutine JacobiGL( n, alpha, beta, x )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! connect this subroutine to the one for DGFEM: vector x should be provided
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  implicit none

!  include 'grid.h'
  integer ( kind = 4 ) :: n,nm2,i

  real ( kind = 8 ),dimension(:),allocatable :: w_GQ,x_GQ
  real ( kind = 8 ) :: x(n+1)
  real ( kind = 8 ) :: alpha,alphaP1
  real ( kind = 8 ) :: beta,betaP1
n=n+1
if(n == 1) then
  x(1)=-1.0; x(2)=1.0;return  
endif

alphaP1=alpha+1.d00
betaP1=beta+1.d00
nm2=n-2

allocate(w_GQ(n-2))
allocate(x_GQ(n-2))
call JacobiGQ( nm2, alphaP1, betaP1, x_GQ, w_GQ )

x(1)=-1.d00
do i=1,n-2
x(i+1)=x_GQ(i)
enddo
x(n)=1.d00

deallocate(w_GQ)
deallocate(x_GQ)
return
end subroutine JacobiGL


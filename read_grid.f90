! Subroutine to read the grid information
!
subroutine read_grid( m, nfp, v_x, etov )

! External shared variables
implicit none
include 'grid.h'

! INternal variables
integer ( kind = 4 ) :: i, m, nfp
integer ( kind = 4 ) :: etov(m,nfp)
real ( kind = 8 ) :: v_x(m+1)


open(10,file='grid.dat',status='unknown')
! Read vertix
read(10,*) V_x
do i=1,m
read(10,*) EToV(i,:)
enddo
close(10)

return
end subroutine read_grid

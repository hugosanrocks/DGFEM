! Purpose, find the corresponding wave speed at each node tacking into account the wave speed model
! (model_input.dat)

subroutine localvelo(advec_mesh)

! External variables
implicit none
include 'grid.h'
TYPE (mesh) :: advec_mesh

! Internal variables
integer ix,idx,nx,nv,iunit,totalnodes,i,j
real*8, dimension(:), allocatable :: vel(:), coor(:), xr(:), speed(:)

! Read wave speed model header
iunit=8
open(iunit,file='vel_header.in',status='unknown')
read(iunit,*) nv
write(*,*) ' Points conforming velocity field (header_velo.in) ', nv

! Read velocity and coordinates
allocate(vel(nv))
allocate(coor(nv))

open(11,file='model_vel.in',status='unknown')
open(12,file='model_coor.in',status='unknown')
read(11,*) vel(:)
read(12,*) coor(:)
close(11)
close(12)

! Points inside elements
totalnodes=advec_mesh%N_node*advec_mesh%N_element
! Reshape node coordinates from matrix to vector
allocate(xr(totalnodes))
allocate(speed(totalnodes))

ix=1
do j=1,advec_mesh%N_element
do i=1,advec_mesh%N_node
xr(ix)=advec_mesh%r_grid(i,j)
ix=ix+1
enddo
enddo
!=====================================

! Assign corresponding velocity to each node
do idx=2,nv
do ix=1,totalnodes
if ( (xr(ix) .lt. coor(idx)) .and. (xr(ix) .ge.coor(idx-1)) )then
speed(ix)=vel(idx-1)
endif
enddo
enddo

! Reshape local (element) velocity into a matrix Nodes X Elements
ix=1
do j=1,advec_mesh%N_element
do i=1,advec_mesh%N_node
advec_mesh%speedv(i,j)=speed(ix)
ix=ix+1
enddo
enddo

! Free memory
deallocate(speed)
deallocate(xr)
deallocate(vel)
deallocate(coor)

end subroutine localvelo

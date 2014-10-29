!==================================
!Subrutine to create the velocity model in 1D medium and its coordinates
!=================================
subroutine model(advec_mesh)


! External variables
implicit none
INCLUDE 'grid.h'
TYPE (mesh) :: advec_mesh


! Internal variables
real*8 :: vel_start,coor_start,vel_i,coor_i,vel_d
integer :: nx,ix,idx,iunit,iunit1,iunit2
real*8, allocatable,dimension(:) :: vel,coor

!============================= file units
iunit=7
open(iunit,file='model_vel.in',status='unknown')
iunit1=iunit+1
open(iunit1,file='model_coor.in',status='unknown')

open(14,file='model_data.in',status='unknown')

write(*,*) ' Creating velocity model (model_data.in) '
read(14,*) nx
write(*,*) ' Entrer the number of points: ', nx
read(14,*) vel_start,coor_start
write(*,*) ' Entrer the starting velocity and its coordinate:', vel_start, coor_start
read(14,*) vel_i,coor_i
write(*,*) ' Entrer the velocity and coordinate increment to be used:', vel_i, coor_i


allocate(vel(nx));allocate(coor(nx))
vel(1)=vel_start
coor(1)=coor_start

do ix=2,nx
vel(ix)=vel(ix-1)+vel_i
coor(ix)=coor(ix-1)+coor_i
enddo

! Write parameters used to build the wave speed model
iunit2=iunit+2
open(iunit2,file='vel_header.in',status='unknown')
write(iunit2,*) nx



!=========================== 
1000 continue

write(*,*) ' starting index for changing velocity increment and its coordinate '
write(*,*) ' < 0 to continue '
read(*,*) idx
if(idx < 0) goto 2000
write(*,*) '  Enter velocity increment: '
read(*,*) vel_d
do ix=idx+1,nx
vel(ix)=vel(ix-1)+vel_d
enddo
write(iunit2,*) idx, vel_d, vel(idx+1)
goto 1000
2000 continue

!======================== output values
write(iunit,*) vel(:)
write(iunit1,*) coor(:)
close(iunit)
close(iunit1)
close(iunit2)


!======================== free memory
deallocate(vel);deallocate(coor)
return
end subroutine model

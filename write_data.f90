subroutine write_data(advec_mesh)

implicit none
INCLUDE 'grid.h'
TYPE (mesh) :: advec_mesh

integer*4 i, j, iunit

iunit=11

open(iunit,file='data.out',status='unknown',action='write')

write(iunit,*) ' Program of advection in 1D '

write(iunit,*) ' End of the reading '
write(iunit,*) ' Number of elements ',advec_mesh%N_element
write(iunit,*) ' N order of polynomials ',advec_mesh%N_order
write(iunit,*) ' Nodes inside 1 element ',advec_mesh%N_node
write (iunit,*) ' GL coordinates ',advec_mesh%x_phy
write(iunit,*) ' Vandermonde matrix '
do i=1,advec_mesh%N_node
write(iunit,*) advec_mesh%Vand(i,:)
enddo
write(iunit,*) ' Dmatrix1d '
do i=1,advec_mesh%N_node
write(iunit,*) advec_mesh%D_r(i,:)
enddo
write(iunit,*) ' Lift1D '
do i=1,advec_mesh%N_node
write(iunit,*) advec_mesh%lift(i,:)
enddo
write(iunit,*) ' Coordinates of nodes '
do i=1,advec_mesh%N_node
write(iunit,*) advec_mesh%r_grid(i,:)
enddo
write(iunit,*) ' Jacobian of coordinates '
do i=1,advec_mesh%n_node
write(iunit,*) advec_mesh%vjac(i,:)
enddo
write(iunit,*) ' Coordinates of masks '
do i=1,advec_mesh%N_face
write(iunit,*) advec_mesh%F_x(i,:)
enddo
write(iunit,*) ' F_mask ', advec_mesh%F_mask

write(iunit,*) ' Normal vectors '
do i=1,advec_mesh%N_face
write(iunit,*) advec_mesh%normal_x(i,:)
enddo
write(iunit,*) ' F scale '
do i=1,advec_mesh%N_face
write(iunit,*) advec_mesh%f_scale(i,:)
enddo
write(iunit,*) ' EToE '
do i=1,advec_mesh%N_element
write(iunit,*) advec_mesh%EToE(i,:)
enddo
write(iunit,*) ' EToF '
do i=1,advec_mesh%N_element
write(iunit,*) advec_mesh%EToF(i,:)
enddo
write(iunit,*) ' vmapM '
write(iunit,*) advec_mesh%vmapM
write(iunit,*) ' vmapP '
write(iunit,*) advec_mesh%vmapP
write(iunit,*) ' vmapB '
write(iunit,*) advec_mesh%vmapB
write(iunit,*) ' mapB '
write(iunit,*) advec_mesh%mapB

open(13,file='ufinal.out',status='unknown')

do j=1,advec_mesh%N_element
do i=1,advec_mesh%N_node
write(13,*) advec_mesh%r_grid(i,j), advec_mesh%u(i,j)
enddo
enddo

close(13)
close(iunit)











end subroutine write_data

! Main driver to run 1D advection program using a DG-FEM.
! Hugo S. Sanchez Reyes 03/07/14
! Universite de Grenoble 1 "Joseph Fourier"
! Supervisors: J. Tago and J. Virieux

program advection

!----------------------Definition of variables--------------------------------------

implicit none

!Define all the variables ontained inside grid.h, variables shared by main programs
INCLUDE 'grid.h'
TYPE (mesh) :: advec_mesh

!Variables needed only by this program
integer*4 i,j,k,l
integer*8 m1,m2,m3 !time mesurement
!------------------------------------------------------------------------------------



write(*,*) ' Program of advection in 1D '

!N_node = number of nodes inside element, N_face number of faces at each element
advec_mesh%N_nodef=1
advec_mesh%N_face=2                          ! 1D two faces per element
!

! Read input data = number of nodes, order of polynomials (see more advec_input.dat)
call read_input(advec_mesh)

allocate(advec_mesh%EToV(advec_mesh%N_element,advec_mesh%N_face))
allocate(advec_mesh%V_x(advec_mesh%N_element+1))
call read_grid(advec_mesh%N_element,advec_mesh%N_face,advec_mesh%V_x,advec_mesh%EToV)

write(*,*) ' End of the reading '
write(*,*) ' Number of elements ',advec_mesh%N_element
write(*,*) ' N order of polynomials ',advec_mesh%N_order

advec_mesh%N_node=advec_mesh%N_order+1
write(*,*) ' Nodes inside 1 element ',advec_mesh%N_node, ' Output file data.out '

! Compute basic Legendre Gauss Lobatto grid
allocate(advec_mesh%x_phy(advec_mesh%N_node))
call JacobiGL(advec_mesh%N_order, 0.d00, 0.d00, advec_mesh%x_phy )

!adjustment due to last call
advec_mesh%N_order=advec_mesh%N_order-1

! Evaluate Vandermonde matrix
allocate(advec_mesh%Vand(advec_mesh%N_node,advec_mesh%N_order+1))
call Vandermonde1D( advec_mesh%N_node, advec_mesh%N_order, advec_mesh%x_phy, advec_mesh%Vand )

! Calculate Derivative matrix
allocate(advec_mesh%D_r(advec_mesh%N_node,advec_mesh%N_order+1))
call Dmatrix1D( advec_mesh%N_node, advec_mesh%N_order, advec_mesh%x_phy, advec_mesh%Vand, advec_mesh%D_r)

! Lift1D
allocate(advec_mesh%lift(advec_mesh%N_node,advec_mesh%N_face))
call Lift1D(advec_mesh%N_node,advec_mesh%N_order,advec_mesh%N_face,advec_mesh%Vand,advec_mesh%Lift)

! Obtain physical coordinates of nodes
allocate(advec_mesh%r_grid(advec_mesh%N_node,advec_mesh%N_element))
allocate(advec_mesh%r_grid1(advec_mesh%N_node,advec_mesh%N_element))
allocate(advec_mesh%r_grid2(advec_mesh%N_node,advec_mesh%N_element))
call coor_nodes( advec_mesh%N_element, advec_mesh%N_node, advec_mesh%N_face, advec_mesh%V_x, advec_mesh%EToV, &
     advec_mesh%x_phy, advec_mesh%r_grid)

!Staggered grids
	advec_mesh%r_grid1(:,:)=advec_mesh%r_grid(:,:)-0.5d00
	advec_mesh%r_grid2(:,:)=advec_mesh%r_grid(:,:)+0.5d00

! Compute Jacobians of transformation
allocate(advec_mesh%r_x(advec_mesh%N_node,advec_mesh%N_element))
allocate(advec_mesh%vJac(advec_mesh%N_node,advec_mesh%N_element))
call GeometricFactors1D(advec_mesh%N_element,advec_mesh%N_node,advec_mesh%D_r,advec_mesh%r_grid, &
     advec_mesh%r_x,advec_mesh%vJac)

! Find masks at element faces
allocate(advec_mesh%F_x(advec_mesh%N_face,advec_mesh%N_element))
allocate(advec_mesh%F_mask(advec_mesh%N_face))
call Masks1D(advec_mesh%N_element,advec_mesh%N_node,advec_mesh%N_face, &
     advec_mesh%x_phy,advec_mesh%r_grid,advec_mesh%F_x,advec_mesh%f_mask)

! Normal vectors at faces
Allocate(advec_mesh%normal_x(advec_mesh%N_face,advec_mesh%N_element))
allocate(advec_mesh%F_scale(advec_mesh%N_face,advec_mesh%N_element))
call Normals1D(advec_mesh)

! Connection between elements
allocate(advec_mesh%EToE(advec_mesh%N_element,advec_mesh%N_face))
allocate(advec_mesh%EToF(advec_mesh%N_element,advec_mesh%N_face))
call Connect1D(advec_mesh)

! Build map of face nodes
allocate(advec_mesh%vmapM(advec_mesh%N_face*advec_mesh%N_element))
allocate(advec_mesh%vmapP(advec_mesh%N_face*advec_mesh%N_element))
allocate(advec_mesh%vmapB(advec_mesh%N_face))
allocate(advec_mesh%mapB(advec_mesh%N_face))
call BuildMaps1D(advec_mesh)



! Start initial conditions
allocate(advec_mesh%u(advec_mesh%N_node,advec_mesh%N_element))
allocate(advec_mesh%utvd(advec_mesh%N_node,advec_mesh%N_element))
allocate(advec_mesh%rhsu(advec_mesh%N_node,advec_mesh%N_element))
! SET INITIAL CONDITION
i=1
j=1
do while (advec_mesh%r_grid(1,i) .lt. 0.5d00)
! Values = 1 before the discontinuity
advec_mesh%u(:,i)=1.d00
i=i+1
enddo
k=i



if (advec_mesh%r_grid(i-1,j-1) .lt. 0.d00)then
do i=6,advec_mesh%N_element
do j=1,advec_mesh%N_node
!advec_mesh%u(j,i)=0.0d00     ! 0 After the wavefront
advec_mesh%u(j,i)=(advec_mesh%r_grid(j,i)+abs(advec_mesh%r_grid(k,l)))/10.d00 ! Slope after wavefront
enddo
enddo
else
do i=k,advec_mesh%N_element
do j=1,advec_mesh%N_node
!advec_mesh%u(j,i)=0.0d00    ! 0 after the wavefront
advec_mesh%u(j,i)=(advec_mesh%r_grid(j,i)-abs(advec_mesh%r_grid(1,k)))/10.d00 ! Slope after the wavefront
enddo
enddo
endif

! Writing down initial conditions
open(13,file='uini.out',status='unknown')

! END OF SETTING INITIAL CONDITION



do j=1,advec_mesh%N_element
do i=1,advec_mesh%N_node
write(13,*) advec_mesh%r_grid(i,j), advec_mesh%u(i,j)
enddo
enddo
close(13)


! Set homogeneous or heterogeneous wave speed advec_mesh%speedv vector of wave speed
  allocate(advec_mesh%speedv(advec_mesh%N_node,advec_mesh%N_element))

  call model(advec_mesh)
  call localvelo(advec_mesh)

  allocate(advec_mesh%dru(advec_mesh%N_node,advec_mesh%N_element))
  allocate(advec_mesh%drur(advec_mesh%N_node*advec_mesh%N_element))

advec_mesh%it=0
m1=time8()

call Advec1D(advec_mesh)
m2=time8()

m3=m2-m1

print *,' Computational time: ', m3


! Writing output
call write_data(advec_mesh)


! Close memory requested
deallocate(advec_mesh%utvd)
deallocate(advec_mesh%speedv)
deallocate(advec_mesh%drur)
deallocate(advec_mesh%dru)
deallocate(advec_mesh%rhsu)
deallocate(advec_mesh%u)
deallocate(advec_mesh%vmapB)
deallocate(advec_mesh%mapB)
deallocate(advec_mesh%vmapM)
deallocate(advec_mesh%vmapP)
deallocate(advec_mesh%EToF)
deallocate(advec_mesh%EToE)
deallocate(advec_mesh%F_scale)
deallocate(advec_mesh%normal_x)
deallocate(advec_mesh%F_mask)
deallocate(advec_mesh%F_x)
deallocate(advec_mesh%r_x)
deallocate(advec_mesh%vJac)
deallocate(advec_mesh%r_grid)
deallocate(advec_mesh%r_grid1)
deallocate(advec_mesh%r_grid2)
deallocate(advec_mesh%V_x)
deallocate(advec_mesh%EToV)
deallocate(advec_mesh%lift)
deallocate(advec_mesh%D_r)
deallocate(advec_mesh%Vand)
deallocate(advec_mesh%x_phy)

end program advection

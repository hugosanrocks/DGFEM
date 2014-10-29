!
! Purpose Connectivity and boundary tables for nodes given in the N_element #
!         of elements, each with N_order+1 degrees of freedom
!
subroutine BuildMaps1D(advec_mesh)

  implicit none
  include 'grid.h'

  type(mesh) :: advec_mesh

  integer :: ielement, iface, i, conections, j
  integer :: TotalNodes, k2, k1, f1, f2
  integer,dimension(:,:),allocatable :: nodeids
  real*8 :: x1,x2,D
  integer,dimension(:,:),allocatable :: vm, vp
  integer :: vidM,vidP
  real (kind=8),dimension(:),allocatable :: res

  TotalNodes=advec_mesh%N_element*advec_mesh%N_node
  allocate(nodeids(advec_mesh%N_node,advec_mesh%N_element))
  allocate(res(totalnodes))

  i=1
  do ielement=1,advec_mesh%N_element
  do iface=1,advec_mesh%N_node
  nodeids(iface,ielement)=i
  i=i+1
  enddo
  enddo

  i=1
  do ielement=1,advec_mesh%N_element
  do iface=1,advec_mesh%N_node
!  PRINT *, I, ADVEC_MESH%R_GRID(IFACE,IELEMENT), RES(I)
  res(i) = advec_mesh%r_grid(iface,ielement)
!  PRINT *, I, ADVEC_MESH%R_GRID(IFACE,IELEMENT), RES(I)
  i=i+1
  enddo
  enddo

  allocate(vm(advec_mesh%N_face,advec_mesh%N_element))
  allocate(vp(advec_mesh%N_face,advec_mesh%N_element))
  do iface=1,advec_mesh%N_face
  vm(iface,:)=nodeids(advec_mesh%F_mask(iface),:)
  enddo

  do ielement=1,advec_mesh%N_element
  do iface=1,advec_mesh%N_face
  k2 = advec_mesh%EToE(ielement,iface)
  f2 = advec_mesh%EToF(ielement,iface)
  vidM = vm(iface,ielement)
  vidP = vm(f2,k2)
  x1= res(vidM)
  x2= res(vidP)
!% Compute distance matrix
 D = (x1 -x2 )**2;
 if ( D .lt. 10e-10 )then
 vp(iface,ielement)=vidP
 endif
 enddo
 enddo

 conections=advec_mesh%N_element*advec_mesh%N_face
 i=1
 do ielement=1,advec_mesh%N_element
 do iface=1,advec_mesh%N_face
! print *, advec_mesh%vmapm(i),vm(iface,ielement)
 advec_mesh%vmapM(i)=vm(iface,ielement)
 advec_mesh%vmapP(i)=vp(iface,ielement)
! print *, advec_mesh%vmapM(i), vm(iface,ielement)
 i=i+1
enddo 
enddo


j=1
do i=1,conections
if (advec_mesh%vmapM(i) .eq. advec_mesh%vmapP(i))then
advec_mesh%mapB(j)=i
advec_mesh%vmapB(j)=advec_mesh%vmapM(i)
j=j+1
endif
enddo
advec_mesh%mapI=1
advec_mesh%mapO=advec_mesh%N_element*advec_mesh%N_face


  deallocate(res)
  deallocate(nodeids)
  deallocate(vm)
  deallocate(vp)
  return
end subroutine BuildMaps1D

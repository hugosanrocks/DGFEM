! Purpose: Estimate RHS of the 1D advection scheme.


subroutine AdvecRHS1D(advec_mesh)
!
!
! External variables
  implicit none
  include 'grid.h'
  include 'rk.h'
  type(mesh) :: advec_mesh

! INternal variables used
  integer :: ielement, iface, inode, i, j, k
  integer :: TotalNodes
  real (kind=8) :: speed, alpha_f, uin
  real (kind=8), dimension(:,:), allocatable :: snx, du, Dru, fsn, lfsn, jump, snxhet
  real (kind=8), dimension(:), allocatable :: ur

  CHARACTER (len=5) :: caracter


  !Variables to estimate differences between elements
  allocate(snx(advec_mesh%N_face,advec_mesh%N_element))
  allocate(du(advec_mesh%N_face,advec_mesh%N_element))
  allocate(ur(advec_mesh%N_node*advec_mesh%N_element))

  allocate(snxhet(advec_mesh%N_face,advec_mesh%N_element))

  ! Constant 'alpha=0'= upwind flux 
  alpha_f=0.d00
  ! Nornal vectors * local speed advection
  ! snx = homogeneous speed = advec_mesh%speed
!  snx=advec_mesh%speed*advec_mesh%Normal_x(:,:)
  totalnodes=advec_mesh%N_node*advec_mesh%N_element


  !Smooth variation of speed
  !snxhet heterogeneous speed 
  do i=1,advec_mesh%N_element
  snxhet(:,i)=advec_mesh%speedv(1,i)*advec_mesh%Normal_x(:,i)
  enddo


  !Field difference at faces

  !Re-shape advec_mesh%u in a row vector ur
  i=1
  do ielement=1,advec_mesh%N_element
  do inode=1,advec_mesh%N_node
  ur(i)=advec_mesh%u(inode,ielement)
  i=i+1
  enddo
  enddo

!---DERIVATIVES AND WRITING IS DONE IN Advec1D.f90------

! Jump between elements phi_j+1/2 - phi_j-1/2
  allocate(jump(advec_mesh%N_face,advec_mesh%N_element))
  i=1
!  write(11,*) 'jump'
  do ielement=1,advec_mesh%N_element
  do iface=1,advec_mesh%N_face
  jump(iface,ielement)=ur(advec_mesh%vmapM(i))-ur(advec_mesh%vmapP(i))
  i=i+1
  enddo
  enddo  

!-----------------------------------

  i=1 
  do ielement=1,advec_mesh%N_element
  do iface=1,advec_mesh%N_face
!  Constant speed
!  du(iface,ielement)=(ur(advec_mesh%vmapM(i))-ur(advec_mesh%vmapP(i)))*(snx(iface,ielement)- &
!  ((1.d00-alpha_f)*abs(snx(iface,ielement))))/2.d00
!  Modified speed
  du(iface,ielement)=(ur(advec_mesh%vmapM(i))-ur(advec_mesh%vmapP(i)))*(snxhet(iface,ielement)- &
  ((1.d00-alpha_f)*abs(snxhet(iface,ielement))))/2.d00

  i=i+1
  enddo
  enddo

!Set boundary conditions
uin=1.d00
!Beginning
du(advec_mesh%vmapB(1),1)=(ur(advec_mesh%vmapB(1))-uin)*(snxhet(1,1)-((1.d00-alpha_f)* &
abs(snxhet(1,1))))/2.d00
!End
du(advec_mesh%N_face,advec_mesh%N_element)=0.d00


  !Estimated local derivative (non-scaled)
  allocate(Dru(advec_mesh%N_node,advec_mesh%N_element))

! First part of first element of RHS term  
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
   Dru(i,j) = 0.d00
    DO k=1,advec_mesh%N_node
    Dru(i,j) = Dru(i,j) + advec_mesh%D_r(i,k)*advec_mesh%u(k,j)
    ENDDO
   ENDDO
  ENDDO

! Second part of first element of RHS term  
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
!    Constant speed
!    Dru(i,j) = -1.d00*advec_mesh%speed*Dru(i,j)*advec_mesh%r_x(i,j)
!    Modified speed
    Dru(i,j) = -1.d00*advec_mesh%speedv(1,j)*Dru(i,j)*advec_mesh%r_x(i,j)
   ENDDO
  ENDDO

! First part of second element of RHS term  
  allocate(fsn(advec_mesh%N_face,advec_mesh%N_element))
  DO i=1,advec_mesh%N_face
   DO j=1,advec_mesh%N_element
    fsn(i,j) = advec_mesh%F_scale(i,j)*du(i,j)
   ENDDO
  ENDDO

! Second part of second element of RHS term  
  allocate(lfsn(advec_mesh%N_node,advec_mesh%N_element))
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
   lfsn(i,j) = 0.d00
    DO k=1,advec_mesh%N_face
    lfsn(i,j) = lfsn(i,j) + advec_mesh%LIFT(i,k)*fsn(k,j)
    ENDDO
   ENDDO
  ENDDO

! Addition of both elements of RHS term
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
   advec_mesh%rhsu(i,j)=Dru(i,j)+lfsn(i,j)
   ENDDO
  ENDDO

! Viscocity term additon
  DO i=1,advec_mesh%N_element
   DO j=1,advec_mesh%N_face
   advec_mesh%rhsu(advec_mesh%F_mask(1),i)=advec_mesh%rhsu(advec_mesh%F_mask(1),i)+(jump(j,i))
   advec_mesh%rhsu(advec_mesh%F_mask(2),i)=advec_mesh%rhsu(advec_mesh%F_mask(2),i)+(jump(j,i))
   ENDDO
  ENDDO

! Free memory
  deallocate(snxhet)
  deallocate(jump)
  deallocate(lfsn)
  deallocate(fsn)
  deallocate(Dru)
  deallocate(du)
  deallocate(ur)

end subroutine AdvecRHS1D

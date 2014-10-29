subroutine Connect1D(advec_mesh)

  implicit none
  include 'grid.h'

  ! local quantities

  integer :: TotalFaces
  integer :: vertex_n(2)

  integer :: iface,ielement,kelement,kface,sk,TotalFaces_2
  integer, dimension(:,:), allocatable :: FToV,FtoF,speye,FTOVT
  integer, dimension(:), allocatable :: faces1,faces2
  integer, dimension(:), allocatable :: face1,face2
  integer, dimension(:), allocatable :: element1,element2

  type(mesh) :: advec_mesh

  ! find number of elements and vertices

  if(advec_mesh%N_face /= 2) then
     write(*,*) ' Connect1D - error of initialisation in 1D, one should have the value 2:  ', advec_mesh%N_face
     stop
  endif

  TotalFaces = advec_mesh%N_face*advec_mesh%N_element
  advec_mesh%N_vertex = advec_mesh%N_element+1                ! number of vertex deduced from number of elements

  ! list of local face to local vertex connections

  vertex_n=(/ 1,2 /)

  ! build global face to node sparse array

  allocate(FToV(TotalFaces,advec_mesh%N_vertex))
  allocate(FToVt(advec_mesh%N_vertex,totalfaces))
  allocate(FtoF(TotalFaces,TotalFaces))

  ! speye identity matrix of TotalFaces X TotalFaces
  allocate(speye(TotalFaces,TotalFaces))
  !  print *, totalfaces,advec_mesh%n_vertex
  ! Initializing matrix ftov 
  do sk=1,totalfaces
  do iface=1,advec_mesh%n_vertex
  ftov(sk,iface)=0
  enddo
!  print *, ftov(sk,:)
  enddo

  sk=1
  do ielement=1,advec_mesh%N_element
     do iface=1,advec_mesh%N_face
        ! make the link between faces of each element (global index) and vertice index 
        FToV(sk,advec_mesh%EToV(ielement,vertex_n(iface))) = 1 
        sk=sk+1
     enddo
  enddo
 ! do iface=1,totalfaces
 ! print *, ftov(iface,:)
 ! enddo
  ! Build global face to global face sparse array

  ! Build sparse identity matrix
  do sk=1,TotalFaces
  do iface=1,TotalFaces
  if (sk .eq. iface) then
  speye(sk,iface)=1
  endif
  if (sk .ne. iface) then
  speye(sk,iface)=0
  endif
  enddo
  enddo
! do iface=1,totalfaces
! print *, speye(iface,:)
! enddo
  !Transpose FToV matrix
  do iface=1,advec_mesh%N_vertex
  do sk=1,TotalFaces
  FToVT(IFACE,SK)=FTOV(SK,IFACE)
  enddo
  enddo

  DO sk=1,TotalFaces
   DO iface=1,TotalFaces
   FToF(sk,iface) = 0.0
    DO ielement=1,advec_mesh%N_vertex
    FToF(sk,iface) = FToF(sk,iface) + FToV(sk,ielement)*FToVT(ielement,iface)
    ENDDO
   ENDDO
  ENDDO

  FToF=FToF-speye
!  PRINT *, 'FTOF'
!  print *, FtoF
!  print *, ''

  ! Find complete face to face connections and store related indexes faces1 and faces2
  ! Faces1 and Faces2 have always TotalFaces-2 length
  allocate(faces1(TotalFaces))
  allocate(faces2(TotalFaces))
  
  totalfaces_2=totalfaces-2
  sk=1
  do kface=1,TotalFaces
     do iface=1,kface-1          ! lower triangle
        if(FToF(kface,iface) == 1 ) then
           faces2(sk)=kface
           faces1(sk)=iface
           sk=sk+1
           if(sk > TotalFaces) then
              write(*,*) ' error of connected structure Connect1D '
              stop
           endif
        endif
     enddo    ! iface
     do iface=kface+1,TotalFaces  ! upper triangle
        if(FToF(kface,iface) == 1 ) then
           faces2(sk)=kface
           faces1(sk)=iface
           sk=sk+1
           if(sk > TotalFaces) then
              write(*,*) ' error of connected structure Connect1D '
              stop
           endif
        endif
     enddo    ! iface
  enddo
!print *, faces1
!print *, faces2
  ! Convert face global number to element and face numbers
  TotalFaces_2=TotalFaces-2

  allocate(element1(TotalFaces_2))
  allocate(element2(TotalFaces_2))
  allocate(face1(TotalFaces_2))
  allocate(face2(TotalFaces_2))

  do iface=1,TotalFaces_2

     element1(iface) =   (faces1(iface)-1)/advec_mesh%N_face  + 1
     face1(iface) =      mod( (faces1(iface)-1),advec_mesh%N_face ) + 1
     element2(iface) =   (faces2(iface)-1)/advec_mesh%N_face + 1
     face2(iface) =      mod( (faces2(iface)-1),advec_mesh%N_face ) + 1
 ! print *, element1(iface), face1(iface), element2(iface), face2(iface)
  enddo

  advec_mesh%EToE(:,:)=1
  advec_mesh%EToF(:,:)=1
  sk=1
  do iface=1,advec_mesh%N_face
    do ielement=1,advec_mesh%N_element
    advec_mesh%EToE(IELEMENT,IFACE)=IELEMENT
    advec_mesh%EToF(ielement,iface)=iface 
    sk=sk+1
    enddo
   enddo
!  print *, advec_mesh%N_element,advec_mesh%N_face
!  do iface=1,advec_mesh%n_element
!  print *, advec_mesh%etoe(iface,:)
!  enddo



  do sk=1,TotalFaces_2
  advec_mesh%etoe(element1(sk),face1(sk))=element2(sk)
  advec_mesh%etof(element1(sk),face1(sk))=face2(sk)
  enddo



  deallocate(FToV)
  deallocate(FTOVT)
  deallocate(FtoF)
  deallocate(speye)
  deallocate(faces1)
  deallocate(faces2)
  deallocate(element1)
  deallocate(element2)
  deallocate(face1)
  deallocate(face2)

  return
end subroutine Connect1D

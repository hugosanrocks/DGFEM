!
! Purpose compute outward pointing normals at elements faces
!
subroutine Normals1D(advec_mesh)

implicit none
include 'grid.h'

type(mesh) :: advec_mesh
integer ( kind = 4 ) :: i, j

advec_mesh%normal_x(1,:)= -1.d00
advec_mesh%normal_x(2,:)=  1.d00

do i=1,advec_mesh%N_face
advec_mesh%F_scale(i,:)=1.d00/advec_mesh%vJac(advec_mesh%F_mask(i),:)
enddo

return
end subroutine Normals1D

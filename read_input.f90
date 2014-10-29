! Program used only to read input data.
!
!
subroutine read_input(advec_mesh)

! External variable definition

implicit none
include 'grid.h'
type(mesh) :: advec_mesh

! INternaal variables
integer ( kind = 4 ) i


open(10,file='advec_input.dat')
read(10,*) advec_mesh%N_element,advec_mesh%N_order   ! nbre of elements and the order level in each element
read(10,*) advec_mesh%FinalTime                      ! Final time of propagation
read(10,*) advec_mesh%re_set                         ! How many step the field is re-set
read(10,*) advec_mesh%imprime                        ! How many step the field is print

! (See advec_input.dat)

return
end subroutine read_input

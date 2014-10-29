!
!
! Purpose compute the metric elements for the local mappings of the 1D elements
!
subroutine GeometricFactors1D(m,n,D_r,r_grid,r_x,vJac)

integer ( kind = 4) :: m, n, i, j, k
real ( kind= 8) D_r(n,n), r_grid(n,m), r_x(n,m), vJac(n,m)

  DO i=1,n
   DO j=1,m
   vJac(i,j) = 0.d00
    DO k=1,n
    vJac(i,j) = vJac(i,j) + D_r(i,k)*r_grid(k,j)
   ENDDO
   ENDDO
  ENDDO

do i=1,n
do j=1,m
r_x(i,j)=1.d00/vJac(i,j)
enddo
enddo

return
end subroutine GeometricFactors1D

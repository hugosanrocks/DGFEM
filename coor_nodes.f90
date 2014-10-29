subroutine coor_nodes( m, np, nfp, V_x, EToV, r, x)

integer ( kind = 4) :: EToV(m,nfp), va(m), vb(m), i, j
real ( kind = 8) :: V_x(m+1), x(np,m), one(np,1), r(np), k

do i=1,m
va(i) = EToV(i,1)
vb(i) = EToV(i,2)
enddo
  DO i=1,m                  !elements
   DO j=1,np                !nodes
!    k = real ( va(i), kind = 8 )
    k=V_x(va(i))
    x(j,i)=k+(0.5*(r(j)+1.d00)*(V_x(vb(i))-V_x(va(i))))
   ENDDO
  ENDDO
return
end subroutine coor_nodes

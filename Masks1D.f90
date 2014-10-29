subroutine Masks1D(m,n,nfp,x_phy,r_grid,Fx,Fmask)

integer ( kind = 4) :: m, n, nfp, i, j
real ( kind= 8) x_phy(n), NODETOL, r_grid(n,m), Fx(nfp,m), fabs1, fabs2
integer ( kind = 4) Fmask(nfp), Fmask1, Fmask2

NODETOL=10E-10

do i=1,n
fabs1=abs(x_phy(i)+1.d00)
if (fabs1 .lt. NODETOL) then
fmask1 = i
endif
fabs2=abs(x_phy(i)-1.d00)
if (fabs2 .lt. NODETOL) then
fmask2 = i
endif
enddo
Fmask(1) = Fmask1
Fmask(2) = Fmask2

do i=1,nfp
do j=1,m
Fx(i,j) = r_grid(Fmask(i), j)
enddo
enddo

return
end subroutine Masks1D

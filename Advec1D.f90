! Purpose : Integrate 1D advection until Final time (using 3th order TVD Runge Kutta)

subroutine Advec1D(advec_mesh)

! External shared variables
implicit none
INCLUDE 'grid.h'
INCLUDE 'rk.h'
TYPE (mesh) :: advec_mesh


!Internal variables
integer*4 i,j,k,l,tstep,intrk,nstep
REAL (KIND=8) :: time, dt, CFL, dx, xmin, Nsteps, a, FINALTIME, timelocal, m, n, o, p
real*8, dimension(:,:), allocatable :: resu, resu1, resu2, resu3
real*8, dimension(:), allocatable :: xmin1
integer (kind=4) :: j8, j9, writ, writ2   ! j8 x,y pairs to fit in linear regresion, j9 pairs not used
real*8  x(4), y(4), slope, const, corr, wavef, ind, ind2, ind3, ind4 ! linear regresion variables

CHARACTER (LEN=7) :: caracter

!m, n, o, p constants to apply TVD integration (see cheng & shu 2007)
m=3.d00/4.d00
n=1.d00/4.d00
o=1.d00/3.d00
p=2.d00/3.d00


print *, ' Final time =', advec_mesh%FINALTIME
time = 0.d00;

! Runge-Kutta residual storage
allocate(resu(advec_mesh%N_node,advec_mesh%N_element))
allocate(resu1(advec_mesh%N_node,advec_mesh%N_element))
allocate(resu2(advec_mesh%N_node,advec_mesh%N_element))

!set residual arrays to zero
do i=1,advec_mesh%N_node
do j=1,advec_mesh%N_element
resu(i,j)=0.d00
resu1(i,j)=0.d00
resu2(i,j)=0.d00
enddo
enddo

! Compute time step size, Xmin and Vmax needed
! Xmin, x minimum between two sequential nodes
allocate(xmin1(advec_mesh%N_element))
do i=1,advec_mesh%N_element
 xmin1(i) = abs(advec_mesh%r_grid(1,i)-advec_mesh%r_grid(2,i))
enddo
xmin=xmin1(1)

! Esitmation of time steps
!CFL=0.75  BOOK HEASTHAVEN
CFL=0.1d00   !CHENG
print *, ' CFL', cfl
dt= (CFL/(advec_mesh%speedv(1,advec_mesh%N_element)))*xmin
dt = 0.5 * dt
Nsteps = ceiling(advec_mesh%FinalTime/dt)
dt = advec_mesh%FinalTime/NstepS
nstep = int (Nsteps,kind=4)
write(6,*) ' dt', dt, 'xmin', xmin, 'steps', nstep

PRINT *, ' Speed at first node: ', ADVEC_MESH%SPEEDV(1,1), 'Speed at last node: ', ADVEC_MESH%SPEEDV(1,advec_mesh%N_element)

! Set iterations to zero
advec_mesh%it=0



!====================================================================================================
! Loop of TVD integration
do tstep=1,Nstep

   advec_mesh%utvd(:,:)=advec_mesh%u(:,:)
   call AdvecRHS1D(advec_mesh)
    resu1(:,:) = advec_mesh%u(:,:) + (dt*advec_mesh%rhsu(:,:))
    advec_mesh%u(:,:) = resu1(:,:)
     call AdvecRHS1D(advec_mesh)
     resu2(:,:) = m*advec_mesh%utvd(:,:) + (n*resu1(:,:)) + (n*(dt*advec_mesh%rhsu(:,:)))
      advec_mesh%u(:,:) = resu2(:,:)
      call AdvecRHS1D(advec_mesh)

   advec_mesh%u(:,:) = (o*advec_mesh%utvd(:,:)) + (p*resu2(:,:)) + (p*(dt*advec_mesh%rhsu(:,:)))

!-------------Derivative---------------------------------------
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
   advec_mesh%dru(i,j) = 0.d00
    DO k=1,advec_mesh%N_node
    advec_mesh%dru(i,j) = advec_mesh%dru(i,j) + advec_mesh%D_r(i,k)*advec_mesh%u(k,j)
    ENDDO
   ENDDO
  ENDDO
! advec_mesh%r_x is the scalling factor to estimate the true derivative value
! due to the mapping in the standar space [-1,1]
  DO i=1,advec_mesh%N_node
   DO j=1,advec_mesh%N_element
    advec_mesh%dru(i,j) = advec_mesh%r_x(i,j)*advec_mesh%dru(i,j)
   ENDDO
  ENDDO

!       Increment of iteration
        advec_mesh%it=advec_mesh%it+1

!----------------TRACK WAVEFRONT----------------------------------
  ! Identify where the instabilities finish
  ! Cycle over elements to evaluate the condition to identify at each element
 do i=2,advec_mesh%N_element-1
   ! Condition to identify where instabilities finish
   ! IND substraction of derivative values at faces of two sequential elements
   ! if there is no oscilations the substraction must be very small.
   ind=abs(advec_mesh%dru(1,i)-advec_mesh%dru(1+2,i))
   ind4=abs(advec_mesh%dru(1,i+1)-advec_mesh%dru(1+2,i+1))
   ! IND2 substraction of the field at i+1 minus at i, the value should be
   ! positive if the slope has no noisy oscilations.
   ind2=advec_mesh%u(1+2,i)-advec_mesh%u(1,i)
   ind3=advec_mesh%u(1+2,i+1)-advec_mesh%u(1,i+1)

!  write(22,*) advec_mesh%It, IND,"LT E-06 AND", IND2,">9E-04 <0.1" !Just to verify values

   ! Condition to !identify where there are no oscilations
!   STRONG CRITERIA. very far away from the wavefront
!  if ((ind .Lt. 1e-06) .and. ((ind2 .gt. 9e-04) .and. (ind2 .lt. 0.1d00)) )then

!   WEAK CRITERIA, closer to the wavefront
    if (((advec_mesh%dru(1,i) .lt. 0.105d00) .and. (advec_mesh%dru(1,i) .gt. 0.097d00)) .and. &
       (((advec_mesh%u(1,i+1) .gt. advec_mesh%u(1,i)) .and. (advec_mesh%u(1,i) .gt. advec_mesh%u(1,i-1) )) &
       .and. (advec_mesh%u(1,i-1) .gt. 0.d00)) ) then

        ! Least square Linear Regresion to estimate where is the wavefront
        j8=3   !pairs (x,y) to estimate linear regresion
        j9=0   !zero pairs to remove from pairs to be used
        do j=1,j8
         x(j)=advec_mesh%r_grid(1+(j),i)
         y(j)=advec_mesh%u(1+(j),i)
        enddo
        ! Linear Regresion subroutine
        call linreg(j8,j9,const,slope,corr,x,y)
        wavef=((-1.d00*const)/slope)
        ! Write the information related to the linear regresion
        open(98,file='linreg.out',status='unknown')
        write(98,*) x(1:3)
        write(98,*) y(1:3)
        write(98,*) slope, const, 0.d00
        write(98,*) advec_mesh%it, wavef, time

        ! Write where the wavefront is located
        open(99,file='track.out',status='unknown')
        write(99,*), advec_mesh%it, wavef, time

!-------RE-SET each 'advec_mesh%re_set' tsteps. This is read from advec_input.dat
     writ=mod(tstep,advec_mesh%re_set)
     if(writ .eq. 0)then
       do k=1,advec_mesh%N_element
       do l=1,advec_mesh%N_node
        !Values = 1, before wavefront location.
        if ((wavef .gt. advec_mesh%r_grid(l,k)))then
        advec_mesh%u(l,k)=1.d00
        else
       !Values = slope, after wavefront location. 
       advec_mesh%u(l,k)=(advec_mesh%r_grid(l,k)-abs(wavef))/10.d00
        endif
       enddo
       enddo
     endif
!-----End of RE-SET--------------------------------------

       GOTO 2000  ! Jump to continue with next time step.

   !End of condition to identify where there are no oscilations
   endif

  !Cycle over elements to identify instabilities
 enddo

2000  CONTINUE

!----Instructions to save field and derivative each 'imprime' iterations.
!----'imprime' is read from advec_input.dat.

     ! Print field and derivative only each 'imprime' steps
     writ2=mod(tstep,advec_mesh%imprime)
     if(writ2 .eq. 0)then
     write(caracter,'(I7.7)') advec_mesh%it!/advec_mesh%imprime
     open(11,file='uder.'//caracter,status='unknown')
       do k=1,advec_mesh%N_element
       do l=1,advec_mesh%N_node
        write(11,*) advec_mesh%r_grid(l,k), advec_mesh%u(l,k), advec_mesh%dru(l,k)
       enddo
       enddo
     endif

!----End of writing instructions



time=time+dt
!End of integration cycle
enddo

close(99)

!Free memory
deallocate(resu1)
deallocate(resu2)
deallocate(xmin1)
deallocate(resu)



return
end subroutine Advec1D

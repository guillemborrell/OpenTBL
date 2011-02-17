!----------------------------------------------------------------------*
! here, all boundary conditions are imposed except the outflow 
! boundary conditions which will be calculated in rhsp
!----------------------------------------------------------------------* 


subroutine boun(ut,vt,wt)

  use alloc_dns,only:inby,ccon,y,idx,x,re,rkc,rkd,idy,dx,dy,delx,vmagic
  use ctesp
  use point
  use omp_lib

  implicit none

  ! -------------------- I/O --------------------------------------------!
  real*8,dimension(0:2*nz2+1,ny+1,ib:ie)::ut,wt
  real*8,dimension(0:2*nz2+1,ny,ib:ie)  ::vt

  ! -------------------- work arrays ------------------------------------!
  integer i,j,kk,k2,k
  !----------------------------------------------------------      
  !Condicion de contorno para Entrada turbulenta
  !----------------------------------------------------------

  ! ================================================
  !        wall and upper boundary 
  ! i=1 is changed in genflu? fix wall, just in case
  ! ================================================


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,kk,k,k2) SCHEDULE(STATIC)
  do kk=0,2*nz2+1,blockl
     k2=min(2*nz2+1,kk+blockl-1)
     do i=ib,ie
        do k =kk,k2
           vt(k,1,i)   =0d0
           vt(k,ny,i)  =0d0   
        enddo
        ! MAGIC NUMBER 
        vt(0,ny,i)=vmagic(i)  

        do k =kk,k2
           ut(k,1,i) = -1d0/inby(2,1)*(inby(2,2)*ut(k,2,i)+       &
                &          inby(2,3)*ut(k,3,i)+inby(2,4)*ut(k,4,i))
           ut(k,ny+1,i) = -1d0/ccon(1,1)*(ccon(1,2)*ut(k,ny,i)+   &
                &          ccon(1,3)*ut(k,ny-1,i)+ccon(1,4)*ut(k,ny-2,i))

           wt(k,1,i) = -1d0/inby(2,1)*(inby(2,2)*wt(k,2,i) +      &
                &          inby(2,3)*wt(k,3,i)+inby(2,4)*wt(k,4,i))   
           wt(k,ny+1,i) = -1d0/ccon(1,1)*(ccon(1,2)*wt(k,ny,i)+   &
                &          ccon(1,3)*wt(k,ny-1,i)+ccon(1,4)*wt(k,ny-2,i))
        enddo
     enddo
  enddo 
end subroutine boun


! ==============================================================
! B.C at the outflow  (assumes at least two planes in last node)   ???????????
! ==============================================================

subroutine outflow_correction(ut,vt,rhsupat,communicator)
  use alloc_dns,only:inby,ccon,y,idx,x,re,rkc,rkd,idy,dx,dy,delx
  use ctesp
  use point
  use omp_lib
  implicit none
  include 'mpif.h'
  integer,intent(in)::communicator
  ! -------------------- I/O --------------------------------------------!
  real*8,dimension(0:2*nz2+1,ny+1,ib:ie)::ut,rhsupat
  real*8,dimension(0:2*nz2+1,ny  ,ib:ie)::vt
  integer ierr

  ! -------------------- work arrays ------------------------------------!
  integer i,j,kk,k2,k
  real*8 sumax, alpha

  sumax=sum(vt(0,ny,ib0:ie))*dx(1)
  if (ib==1)  sumax=sumax-sum(ut(0,2:ny,1)*dy(1:ny-1))
  if (ie==nx) sumax=sumax+sum(ut(0,2:ny,nx)*dy(1:ny-1))

  call MPI_ALLREDUCE(sumax,alpha,1,MPI_real8,MPI_sum,communicator,ierr)      

  alpha= -alpha/(y(ny)-y(1))
  if (ie==nx) then
     !$OMP PARALLEL WORKSHARE
     ut(0,2:ny,nx)   = ut(0,2:ny,nx) + alpha    
     rhsupat(0,2:ny,nx) =rhsupat(0,2:ny,nx)+alpha*idx       
     !$OMP END PARALLEL WORKSHARE
  endif
end subroutine outflow_correction




!-----------------------------------------
!---------------AUXILIARY SUBROUTINE------
!-----Generate the V magic number --------
! MAGIC NUMBER from nagib, chauhan & monkewitz (PTRSA, vol 365:755-779, 2007)
subroutine magicnumber(mpiid)
use alloc_dns,only: re,pi,ax,vmagic,v0
use genmod,only:rthin
use ctesp
implicit none

real*8,dimension(nx)::x,h
real*8 e1,c1,ei1,kap,ckap,cprim,reth,uep,rd
real*8 rx,dx
integer i,mpiid

!see matlab file: blevol.m
 e1=0.8659d0; c1=0.01277d0;   !!reth=c1*rex^e1
 ei1=1/e1

 kap=0.384d0; ckap=4.127d0;   !! uep= log(reth)/kap+ckap 
 cprim=7.135d0              
 dx=ax*pi/(nx-1)

if(mpiid.eq.1) then
write(*,*) '------------------------------------------------------'
write(*,*) 'ax',ax
write(*,*) 'pi',pi
write(*,*) 'rthin',rthin
write(*,*) 'dx',dx
write(*,*) '------------------------------------------------------'

endif

do i=1,nx
 x(i)=dx*(i-1)
 reth=(rthin**ei1+c1**ei1*re*x(i))**e1
 uep= log(reth)/kap+ckap  
 h(i)=reth/((1-cprim/uep)*re)
enddo
!!!!!   vinf/uinf = d(delta*)/dx
vmagic(1)     =(h(2)-h(1))/(x(2)-x(1))
vmagic(2:nx-1)=(h(3:nx)-h(1:nx-2))/(x(3:nx)-x(1:nx-2))
vmagic(nx)    =(h(nx)-h(nx-1))/(x(nx)-x(nx-1))
!vmagic=vmagic-(vmagic(1)-v0(ny))

endsubroutine magicnumber




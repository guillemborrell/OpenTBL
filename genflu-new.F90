subroutine genflu(ut,vt,wt,y,re,dt,tiempo,mpiid,m)
  !------------------------------------------------------------------------*
  !     fluctuating velocity profile from a BL
  !------------------------------------------------------------------------*
  use genmod
  use temporal
  use point
  use alloc_dns,only:kazr,idx,inbx,dy,idy,v0
  use main_val, only:wkp
  use ctesp
  use omp_lib
  implicit none
  include "mpif.h"

  !------------------------------- I/O -------------------------------!
  real*8  y(0:ny+1)
  real*8  tiempo,re,dt
  real*8, dimension(0:2*nz2+1,ny+1,ib:ie):: ut,wt
  real*8, dimension(0:2*nz2+1,ny,ib:ie)  :: vt

  ! ----------------------- Workspace ------------------------------!      
  integer m,k,l,jtop,j,mpiid,kk,i,kll,khh,ii,ntop(2),kl(ny,2)
  real*8  tloc,den,www1,u99,uinf
  real*8  fsca(ny,2),ysca(ny),wlwa(ny)

  ! --------------------- MPI workspaces -----------------------------!
  integer istat(MPI_STATUS_SIZE),ierr,comm,countu,countv,tipo,id
  ! ------------------------ Program -------------------------------------
  jtop=0;ntop=0;kl=0;kll=0;khh=0;fsca=0d0;
  if(mpiid2.eq.0) tm1 = MPI_WTIME()

  countu=2*(nz2+1)*(ny+1)
  countv=2*(nz2+1)*ny
  comm=MPI_COMM_WORLD
  tipo=MPI_REAL8
  tloc=tiempo-timeinit

  if(mpiid2.eq.0) tm3 = MPI_WTIME()

  !  ------  bring reference ref. planes --------
  if (mpiid==mpiout) then    
     call MPI_SEND(ut(0,1,xout),countu,tipo,0,1,comm,istat,ierr)
     call MPI_SEND(vt(0,1,xout),countv,tipo,0,2,comm,istat,ierr)
     call MPI_SEND(wt(0,1,xout),countu,tipo,0,3,comm,istat,ierr)
  endif

  if (mpiid==0) then   
     call MPI_RECV(ut(0,1,1),countu,tipo,mpiout,1,comm,istat,ierr)
     call MPI_RECV(vt(0,1,1),countv,tipo,mpiout,2,comm,istat,ierr)
     call MPI_RECV(wt(0,1,1),countu,tipo,mpiout,3,comm,istat,ierr)

     if(mpiid2.eq.0) then
        tm4 = MPI_WTIME()
        tmp24=tmp24+(tm4-tm3)
     endif
     ! ------   update mean reference profiles -------
     if (tloc < 20 ) then ! ~ \delta/U_\inf
        timec  =20
     elseif (tloc<100) then ! ~ 50 * \delta/U_\inf
        timec   = 100
     else
        timec = tloc
     endif
     um= (dt/timec)*ut(0,:,1)+um*(1d0-dt/timec)
     utauout = sqrt(abs(2d0*(um(2))/(y(2)-y(1))/re))

     !  --- compute momentum thickness -----
     uinf=sum(um(ny+1-20:ny+1-5))/16d0 
     u99 = 0.99*uinf   

     do j=1,ny-2
        if (um(j).ge.u99) exit
        jtop = j+1
     enddo
     dout = y(jtop)       
     rthout = sum((uinf-um(2:jtop+3))*um(2:jtop+3)*dy(1:jtop+2))*re/uinf

     !  --------  scaling factors and interp. limits --------------
     gamma(2)  = rthin/rthout   
     gamma(1)  = gamma(2)**(0.125)

     ! --------- new grids for rescaling
     do k=1,2
        if (k==1) then 
           ysca = 5d-1*(y(1:ny)+y(2:ny+1))/din  ! --- u
        else
           ysca = y(1:ny)/din                   ! --- v
        endif

        call dobleu(ysca,wlwa,www1,ny) 
        wlwa = ysca/gamma(1)+(1d0/gamma(2)-1d0/gamma(1))/www1*wlwa

        kll=1
        fsca(:,k) = 0d0
        do i=1,ny+1       
           if (wlwa(i).le.ysca(ny)) then
              ntop(k) = i
              
              khh = ny
              do while (khh-kll.gt.1)
                 kk=(khh+kll)/2                     
                 if (ysca(kk).ge.wlwa(i)) then
                    khh=kk
                 else 
                    kll=kk
                 endif
              enddo
              kl(i,k) = kll             
              fsca(i,k) = (wlwa(i)-ysca(kll))/(ysca(khh)-ysca(kll))
           endif
        enddo
     enddo
     write(*,*) 'VALORES DE NTOP=',ntop(1),ntop(2)
     tm3 = MPI_WTIME()
     tmp25=tmp25+abs(tm4-tm3)


     tm3 = MPI_WTIME()
     ! -------  the zero mode ---------
     ut(0,:,1) = u0
     ut(1,:,1) = 0d0
     vt(0,:,1) = v0
     vt(1,:,1) = 0d0
     wt(0:1,:,1) = 0d0

#ifndef NOINFOGENFLU    
     write(36) tiempo,rthout,rthin,utauout,gamma(1),gamma(2),din,dout,jtop,&
     &         ntop(1:2),kl(1:ny,1:2),wlwa(1:ny),fsca(1:ny,1:2)                        
     write(36) ut(2:5,1:ny+1,1)                                     
#endif
      call profile(ut,wkp,ntop(1),kl(1,1),fsca(1,1),ny+1)
#ifndef NOINFOGENFLU
      write(36) ut(2:5,1:ny+1,1)     
      write(36) wt(2:5,1:ny+1,1)
#endif   
      call profile(wt,wkp,ntop(1),kl(1,1),fsca(1,1),ny+1)     
#ifndef NOINFOGENFLU
      write(36) wt(2:5,1:ny+1,1)                 
      write(36) vt(2:5,1:ny  ,1)
#endif
      call profile(vt,wkp,ntop(2),kl(1,2),fsca(1,2),ny  )
#ifndef NOINFOGENFLU
      write(36) vt(2:5,1:ny  ,1)
#endif

  endif

#ifndef NOINFOGENFLU
!   if(m.eq.3) then     
!      den=4/(y(2)*re)
!      do j=1,9
!         id=(j-1)*1023        
!         if(mpiid==id) then            
!            cfinfo(j)=ut(0,2,ib)                       
!         endif
!         call MPI_BCAST(cfinfo(j),1,MPI_REAL8,id,MPI_COMM_WORLD,ierr)
!      enddo
!      cfinfo=cfinfo*den     
!      if(mpiid==0) then        
!         write(36,'(14(d22.14))') tiempo,rthout,utauout,gamma(1),dout
!         write(36)     
!         call flush(36)
!      endif
!   endif
#endif

  if(mpiid2.eq.0) then
     tm4 = MPI_WTIME()
     tmp26=tmp26+abs(tm4-tm3)             
     tmp21 = tmp21 + abs(tm4-tm1)               
  endif
end subroutine genflu

! ---------------------------------------------------------------------- !
! ------------------------ Auxiliary functions ------------------------- !
! ---------------------------------------------------------------------- !

subroutine dobleu(eta,www,www1,n)
  implicit none
  real*8 ag,tag,bg,bg2,xxx,yyy,www1,w1default,www(n),eta(n)
  integer i,n,im1
  parameter(ag=4d0, bg=0.2d0, w1default=0.7895)

  bg2=1d0-2d0*bg
  tag=1d0/tanh(ag)
  www1 = w1default    ! --- just in case

  www(1)=0d0
  im1=1
  do i=2,n
     xxx=5d-1*(eta(i)+eta(i-1))
     if (xxx.ge.1d0) then
        yyy=1d0
     else
        yyy=5d-1*(1d0+tag*tanh(ag*(xxx-bg)/(bg2*xxx+bg)))
     endif
     www(i)=www(i-1)+yyy*(eta(i)-eta(i-1))

     if (eta(i).ge.1d0 .and. eta(i-1).lt.1d0) then
        www1=www(i-1)+(1d0-eta(i-1))*yyy        
     endif

  enddo
end subroutine dobleu
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------

subroutine profile(ut,wk,ntop,kl,fsca,n)
  use ctesp
  use genmod  
  implicit none

  real*8, dimension(0:nz1-1,n):: wk,ut
  real*8  fsca(ny),ggg
  integer ntop,kl(ny),ii,i,k,n,kll

  ggg = 1d0/gamma(1) 
  
  if (n==ny+1) then
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
     do i=1,ntop
        kll = kl(i)    
        wk(2:nz1-1,i+1)= ut(2:nz1-1,kll+2)*fsca(i)+ut(2:nz1-1,kll+1)*(1-fsca(i))
     enddo
     wk(2:nz1-1,1) = -wk(2:nz1-1,2)
  else
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
     do i=1,ntop
        kll = kl(i)    
        wk(2:nz1-1,i)= ut(2:nz1-1,kll+1)*fsca(i)+ut(2:nz1-1,kll)*(1-fsca(i))
     enddo
  endif
  
!---------------------- beyond interpolation range ---
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=ntop-2,n     
     wk(2:nz1-1,i)=wk(2:nz1-1,ntop-3) 
  enddo
 

!== COPY BACK====    
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=1,n
     ut(2:nz1-1,i) = wk(2:nz1-1,i)*ggg
  enddo

#ifndef NOINFOGENFLU 
  write(36) wk(3:6,1:n)
#endif
end subroutine profile


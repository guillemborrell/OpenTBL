subroutine genflu(ut,vt,wt,y,re,dt,tiempo,mpiid,m)
  !------------------------------------------------------------------------*
  !     fluctuating velocity profile from a BL
  !------------------------------------------------------------------------*
  use genmod
  use temporal
  use point
  use alloc_dns,only:kazr,idx,inbx,dy,idy,v0
  use main_val, only:wkp,wkpo
  use ctesp
  use omp_lib
  implicit none
  include "mpif.h"

  !------------------------------- I/O -------------------------------!

  real*8  y(0:ny+1),den
  real*8 tiempo,re,dt,uinf,wsca(ny+1)

  real*8, dimension(0:2*nz2+1,ny+1,ib:ie):: ut,wt
  real*8, dimension(0:2*nz2+1,ny,ib:ie)  :: vt

  ! ----------------------- Workspace ------------------------------!      
  integer m,k,l,jtop,j,mpiid,kk,i,kll,khh,kl2,ii,ntop(2),kl(ny+1,2)
  real*8 wlws(0:ny+1),tloc,fsca(1:ny+1,2),u99
  ! --------------------- MPI workspaces -----------------------------!
  integer istat(MPI_STATUS_SIZE),ierr,comm,countu,countv,tipo,id
  ! ------------------------ Program -------------------------------------
  jtop=0;ntop=0;kl=0;kll=0;khh=0;kl2=0;wlws=0d0;fsca=0d0;
  if(mpiid2.eq.0) tm1 = MPI_WTIME()

  countu=2*(nz2+1)*(ny+1)
  countv=2*(nz2+1)*ny
  comm=MPI_COMM_WORLD
  tipo=MPI_REAL8
  tloc=tiempo-timeinit

  if(mpiid2.eq.0) tm3 = MPI_WTIME()

  !computing the Cf...
  !   counter=counter+1
  !   den=1/(y(2)-y(1))/re
  !   cf2(ib:ie)=cf2(ib:ie)+2d0*ut(0,2,ib:ie)*den

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
     gamma(2) = rthin/rthout   
     gamma(1)  = gamma(2)**(0.125)

     do ii=1,2
        kll=1
        fsca(:,ii) = 0d0
        wsca(1:ny+1)=y(0:ny)*gamma(ii)  !inner
        do i=1,ny+1       
           if (y(i-1).le.wsca(ny+1)) then
              ntop(ii) = i

              khh = ny+1
              do while (khh-kll.gt.1);kk=(khh+kll)/2                     
                 if (wsca(kk).gt.y(i-1)) then
                    khh=kk
                 else 
                    kll=kk
                 endif
              enddo
              kl(i,ii) = kll             
              fsca(i,ii) = (y(i-1)-wsca(kll))/(wsca(khh)-wsca(kll))
           endif
        enddo
     enddo

     tm3 = MPI_WTIME()
     tmp25=tmp25+abs(tm4-tm3)


     wlws=y/din
     call dobleu(wlws,wlws,ny+1) ! weigth function of Lundt wu and squires
     tm3 = MPI_WTIME()
     ! -------  the zero mode ---------
     ut(0,:,1) = u0
     ut(1,:,1) = 0d0
     vt(0,:,1) = v0
     vt(1,:,1) = 0d0
     wt(0:1,:,1) = 0d0
     ! ----------------  ACHTUNG ---- what happens to v,w!!!!!
#ifndef NOINFOGENFLU    
     write(36) tiempo,rthout,rthin,utauout,gamma(1),gamma(2),din,dout,jtop,&
     &                          ntop(1:2),kl,wlws,fsca                        
      write(36) ut(2:5,1:ny+1,1)                                     
#endif
      call profile(ut,wkp,wkpo,ntop,kl,fsca,wlws,ny+1)
#ifndef NOINFOGENFLU
      write(36) ut(2:5,1:ny+1,1)     
      write(36) wt(2:5,1:ny+1,1)
#endif   
      call profile(wt,wkp,wkpo,ntop,kl,fsca,wlws,ny+1)     
#ifndef NOINFOGENFLU
      write(36) wt(2:5,1:ny+1,1)                 
      write(36) vt(2:5,1:ny  ,1)
#endif
      call profile(vt,wkp,wkpo,ntop,kl,fsca,wlws,ny  )
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

subroutine dobleu(eta,www,n)
  implicit none
  real*8 ag,tag,bg,bg2,www(0:n),eta(0:n)
  integer i,n
  parameter(ag=4d0, bg=0.2d0)

  bg2=1d0-2d0*bg
  tag=1d0/tanh(ag)
  do i=0,n
     if (eta(i).ge.1d0) then
        www(i)=1d0
     else
        www(i)=5d-1*(1d0+tag*tanh(ag*(eta(i)-bg)/(bg2*eta(i)+bg)))
     endif
  enddo
end subroutine dobleu
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------

subroutine profile(ut,wk,wko,ntop,kl,fsca,wlws,n)
  use ctesp
  use genmod  
  implicit none

  real*8, dimension(nz1,n):: wk,wko,ut
  real*8 wlws(0:ny+1),fsca(1:ny+1,2),ggg
  integer ntop(2),kl(ny+1,2),ii,i,k,n,kll

  ggg = 1d0/gamma(1) 
  
!---------------------- beyond interpolation range ---
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=ntop(1)+1,n     
     wko(3:nz1,i)=ut(3:nz1,n-1)  !n-1 cause @n->ut==0
  enddo
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=ntop(2)+1,n     
     wk(3:nz1,i)=ut(3:nz1,n-1) 
  enddo
 
!---------------------- in interpolation range ---
!==INNER REGION====          
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
  do i=1,min(ntop(1),n)     
     kll = min(kl(i,1),n-1)    
         wko(3:nz1,i)= ut(3:nz1,kll+1)*fsca(i,1)+ut(3:nz1,kll)*(1-fsca(i,1))       
  enddo

!==OUTER REGION====   
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
  do i=1,min(ntop(2),n)
     kll = min(kl(i,2),n-1)    
     wk(3:nz1,i)= ut(3:nz1,kll+1)*fsca(i,2)+ut(3:nz1,kll)*(1-fsca(i,2))
  enddo
  
!==COMPOUND PROFILE====    
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=1,n
     ut(3:nz1,i) = (wko(3:nz1,i)*(1d0-wlws(i))+wk(3:nz1,i)*wlws(i))*ggg
  enddo
#ifndef NOINFOGENFLU 
  write(36) wk(3:6,1:n),wko(3:6,1:n)     
#endif
end subroutine profile


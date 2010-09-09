#ifdef NEWGENFLU

subroutine genflu(ut,vt,wt,y,re,dt,tiempo,mpiid,m,communicator)
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
  integer,intent(in)::communicator
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
  comm=communicator
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
     utauout = sqrt(abs((um(2)-um(1))/(y(2)-y(1))/re))

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
        do i=1,ny       
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
     !      write(*,*) 'Ntop==============',ntop
     tm3 = MPI_WTIME()
     tmp25=tmp25+abs(tm4-tm3)


     tm3 = MPI_WTIME()
     ! -------  the zero mode ---------
     ut(0,:,1) = u0
     ut(1,:,1) = 0d0
     vt(0,:,1) = v0
     vt(1,:,1) = 0d0
     wt(0:1,:,1) = 0d0

     call profile(ut,wkp,ntop(1),kl(1,1),fsca(1,1),ny+1)
     call profile(wt,wkp,ntop(1),kl(1,1),fsca(1,1),ny+1)     
     call profile(vt,wkp,ntop(2),kl(1,2),fsca(1,2),ny  )
  endif

#ifndef NOINFOSTEP
  !Writting info per time substep at different X locations
  call info_per_step(ut,y,tiempo,dt,re,ntop,jtop,mpiid) 
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
end subroutine profile

!------------------------------------------------------
!-------------TIME PER SUBSTEP INFO---------------------
!------------------------------------------------------
subroutine info_per_step(ut,y,tiempo,dt,re,ntop,jtop1,mpiid,communicator) 
  use genmod
  use point
  use ctesp
  use omp_lib
  use alloc_dns,only:dy
  implicit none
  include "mpif.h"
  integer,intent(in)::communicator
  real*8, dimension(0:2*nz2+1,ny+1,ib:ie):: ut
  real*8::  y(0:ny+1)
  integer:: mpiid,nodo,i,j,k,jtop,jtop1,ntop(2)
  real*8:: r_thout(8),u_tau(8),u_inf(8),d_out(8),x_plane(8),buffer(5)
  real*8:: tiempo,dt,uinf,re,u99,utau
  integer istat(MPI_STATUS_SIZE),ierr,comm

!Computing Re_theta at different X position on the BL
!Equispaced nodes==equispaced X planes  
do i=1,8
   nodo=((nummpi-1)*i)/8 !Node from which collect the info 
   if(mpiid.eq.nodo) then
     um=ut(0,:,ib) !here we dont weight <U> with exp(-tau/T)     
     uinf=sum(um(ny+1-20:ny+1-5))/16d0 !instantaneous profile
     u99 = 0.99*uinf   
     utau= sqrt(abs((um(2)-um(1))/(y(2)-y(1))/re))

     do j=1,ny-2
        if (um(j).ge.u99) exit
        jtop = j+1
     enddo
     dout = y(jtop)       
     rthout = sum((uinf-um(2:jtop+3))*um(2:jtop+3)*dy(1:jtop+2))*re/uinf
     !Composing the data buffer to send: R*8     
     buffer(1)=uinf;buffer(2)=utau;buffer(3)=dout;buffer(4)=rthout;buffer(5)=ib*1d0
     call MPI_SEND(buffer,5,MPI_REAL8,0,i,communicator,istat,ierr) 
   endif
enddo

if(mpiid.eq.0) then
        do i=1,8
           nodo=((nummpi-1)*i)/8
           call MPI_RECV(buffer,5,MPI_REAL8,nodo,i,communicator,istat,ierr) 
           u_inf(i)=buffer(1);u_tau(i)=buffer(2);d_out(i)=buffer(3);r_thout(i)=buffer(4);x_plane(i)=buffer(5)                             
        enddo                                
        write(36,'(49(d22.14))') tiempo,rthout,utauout,gamma(1:2),dout,jtop1*1d0,ntop(1:2)*1d0,x_plane,u_inf,u_tau,d_out,r_thout     
        call flush(36)
!         write(*,*) 'tiempo,rthout,utauout,gamma(1:2),dout,jtop1*1d0,ntop(1:2)*1d0,ntop----------------------------'
!         write(*,*) tiempo,rthout,utauout,gamma(1:2),dout,jtop1*1d0,ntop(1:2)*1d0,ntop
!         write(*,*) 'u_inf...u_tau...d_out...r_thout en filas'
!         write(*,*) x_plane
!         write(*,*) u_inf
!         write(*,*) u_tau
!         write(*,*) d_out
!         write(*,*) r_thout
!         write(*,*) '----------------------------------------------------'
     endif
endsubroutine info_per_step


#endif
















!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
! OLD GENFLU..........................

#ifdef OLDGENFLU
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
  integer m,k,l,jtop,j,mpiid,kk,i,kll,khh,kl2,ii,ntop(2),kh(ny+1,2),kl(ny+1,2)
  real*8 wlws(0:ny+1),tloc,fsca(1:ny+1,2)
  ! --------------------- MPI workspaces -----------------------------!
  integer istat(MPI_STATUS_SIZE),ierr,comm,countu,countv,tipo,id
  ! ------------------------ Program -------------------------------------
  jtop=0
  if(mpiid2.eq.0) tm1 = MPI_WTIME()
  if(mpiid2.eq.0) write(*,*) 'USANDO OLD-GENFLU'
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
     din = 0.99*uinf   

     do j=1,ny-2
        if (um(j).ge.din) exit
        jtop = j+1
     enddo

     dout = y(jtop)
     rthout  = 0d0       
     !      uinf=um(ny+1-5)     !!!!!!!!  fix this !!!!!!!
     rthout = sum((uinf-um(2:jtop+3))*um(2:jtop+3)*dy(1:jtop+2))*re/uinf
     !  --------  scaling factors and interp. limits --------------
     gamma(2) = rthin/rthout   
     din    = dout*gamma(2)
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
              kh(i,ii) = khh
              fsca(i,ii) = (y(i-1)-wsca(kll))/(wsca(khh)-wsca(kll))
           endif
        enddo
     enddo

     tm3 = MPI_WTIME()
     tmp25=tmp25+abs(tm4-tm3)


     wlws=y/din
     call dobleu(wlws,wlws,ny+2) ! weigth function of Lundt wu and squires
     tm3 = MPI_WTIME()
     ! -------  the zero mode ---------
     ut(0,:,1) = u0
     ut(1,:,1) = 0d0
     vt(0,:,1) = v0
     vt(1,:,1) = 0d0
     wt(0:1,:,1) = 0d0
     ! ----------------  ACHTUNG ---- what happens to v,w!!!!!
     call profile(ut,wkp,wkpo,ntop,kl,fsca,wlws,ny+1)
     call profile(vt,wkp,wkpo,ntop,kl,fsca,wlws,ny  )
     call profile(wt,wkp,wkpo,ntop,kl,fsca,wlws,ny+1)     
  endif

#ifndef NOINFO
  if(m.eq.3) then     
     den=4/(y(2)*re)
     do j=1,9
        id=(j-1)*1023        
        if(mpiid==id) then            
           cfinfo(j)=ut(0,2,ib)                       
        endif
        call MPI_BCAST(cfinfo(j),1,MPI_REAL8,id,MPI_COMM_WORLD,ierr)
     enddo
     cfinfo=cfinfo*den     
     if(mpiid==0) then        
        write(36,'(14(d22.14))') tiempo,rthout,utauout,gamma(1),dout,cfinfo    
        call flush(36)
     endif
  endif
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
  real*8 ag,tag,bg,bg2,www(*),eta(*)
  integer i,n
  parameter(ag=4d0, bg=0.2d0)

  bg2=1d0-2d0*bg
  tag=1d0/tanh(ag)
  do i=1,n
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
  real*8 wlws(0:ny+1),fsca(1:ny+1,2)
  integer ntop(2),kl(ny+1,2),ii,i,k,n,kll,ggg

  ggg = 1d0/gamma(1) 
  
  if(n.eq.ny) then
   do i=1,n
     if(kl(i,1).eq.n) kl(i,1)=n-1
     if(kl(i,2).eq.n) kl(i,2)=n-1
   enddo
  endif

  !   ---- beyond interpolation range ---
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=ntop(1)+1,n     
     wko(3:nz1,i)=ut(3:nz1,n) 
  enddo
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=ntop(2)+1,n     
     wk(3:nz1,i)=ut(3:nz1,n) 
  enddo
 
  !   ---- in interpolation range ---          
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
  do i=1,min(ntop(1),n)     
     kll = kl(i,1)    
         wko(3:nz1,i)= ut(3:nz1,kll+1)*fsca(i,1)+ut(3:nz1,kll)*(1-fsca(i,1))       
  enddo
 
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,kll) SCHEDULE(STATIC)
  do i=1,min(ntop(2),n)
     kll = kl(i,2)
     wk(3:nz1,i)= ut(3:nz1,kll+1)*fsca(i,2)+ut(3:nz1,kll)*(1-fsca(i,2))
  enddo
  
  ! ump(:,1) -> U^inner_inlet/gamma(1), ump(:,2) -> U^outer_inlet/gamma(1)   
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=1,n
     ut(3:nz1,i) = (wko(3:nz1,i)*(1d0-wlws(i))+wk(3:nz1,i)*wlws(i))*ggg
  enddo
 
end subroutine profile

#endif
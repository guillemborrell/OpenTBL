!------------------------------------------------------------------
!                                                                   */
!  Module used in DNS with definitions of variables to be allocated */
!                                                                   */
!                                         Mark Simens               */
!  idea by Javier Crespo                                            */
!********************************************************************/

!Number of nodes for each BL (BL1 & BL2)
module num_nodes
integer,parameter:: numnodes_1=12,numnodes_2=4
endmodule num_nodes

module ctesp 
  !       CTES CODE FOR 2D STAGGERED GRID FINITE DIFERENCES
  !--------------------------------------------------------------------
  !       declaracion de variables
  !--------------------------------------------------------------------
  integer nx,ny,nz,nz1,nz2,ngz,npos,mcycl,nx1,ny1,nplanes
  integer nummpi,pnodes,mpiid2
  integer xin,xout 
  real*8 totmem

! Parameters for genflu and getstart!
#ifdef BLGRANDE
!====================================================
 parameter ( nx =8193,   ny =711, nz=4096)
 parameter ( xin = 1 , xout =3671*2/2) !50d99-Re2500
!====================================================
#else
!====================================================
 parameter ( nx =1297,   ny =331, nz=768)
 parameter ( xin = 1 , xout =1024) !50d99-Re2500
!====================================================
#endif

  parameter ( nz1 = 2*(nz/3), nz2=nz1/2-1,ngz=nz/2,nx1=nx-1,ny1=ny-1 )
  parameter ( nplanes = nz/3)
  parameter ( mcycl = 40, npos=6)

#ifdef LINES
 !%%%%%%%%%%%%%%%%%%Testing LInes%%%%%%%%%%%
  integer,parameter:: nlineas=500*3
  real*8:: perfil(1:nx,4,nlineas)    !for bigger testing should be ib:ie
  !%%%%%%%%%%%%%%%%%%
#endif

    real*8:: tracking(ny,28,30)

  !blocking for OMP
  integer,parameter:: blockl=50
  ! buffers for real Fourier transform !
  integer nmax
  parameter (nmax=nz+2)  
  !For the outflow BC:
  real*8:: Uinfinity
  parameter (Uinfinity=1d0)

  ! Points for spectra and correlations
  integer nspec,ltot,lxp,ncorr,tots,lxcorr
  integer xci,xco 
  parameter(nspec=20,xci=1,xco=nx)
  parameter(lxp=3,lxcorr=lxp,ncorr=nspec)
  parameter(tots=ncorr*lxcorr*4*nx)
  integer  xpoint(lxp),nxp(lxp),xcorpoint(lxcorr)
 

#ifdef BLGRANDE
  data xpoint /6912,10542,14171/ !Location=50d99;(50d99+XRtau2000)/2;XRtau2000
  data xcorpoint /6912,10542,14171/
#else
  data xpoint /200,600,1200/ !Location=50d99;(50d99+XRtau2000)/2;XRtau2000
  data xcorpoint /200,600,1200/
#endif
  data nxp /62,72,82/    !delta/4 at each X. Averaged Spectra (X-d/4)<X<(X+d/4)

#ifdef PLANESPECTRA 
   !for debugging purposes!! 
   integer  frequency(7)
!   data frequency /1,2,5,10,20,50,100/  !how often we sample
   data frequency /1,2,3,4,5,6,7/  !how often we sample
#endif  

#ifdef PLANESPECTRA2  
   integer::  ss
#endif 

 
#ifdef CREATEPROFILES  
  integer,parameter:: num_planes=300   !Number of planes to composite U0 & V0
  real*8,dimension(:,:),allocatable:: u0c,v0c,u0c_once,v0c_once,w0c_once  !Here are storaged the profiles    
  real*8:: pdiv(ny,nx)
  integer:: paso
#endif
 
end module ctesp
module ctesp2
integer,parameter:: nx=34
endmodule ctesp2 
!---------------------------------------------------------------------
! Variables for genflu 
! --------------------------------------------------------------------!
module genmod
 use ctesp,only: nx,ny,nz,nz2
 implicit none
  real*8:: um (ny+1),u0(ny+1),cf(1:nx),cf2(1:nx),cfinfo(9)
  real*8:: gamma(2),utauout,rthin,rthout,utauin,dout,din
  real*8:: timeinit,timec
  integer:: flaggen,counter
end module genmod

!---------------------------------------------------------------------
module  alloc_dns
  use ctesp,only: nx,ny,nz,nz2
  implicit none
  complex*16:: kaz(0:nz2)

  real*8:: kaz2(0:nz2),kazr(0:nz2),kvis(0:nz2)
  real*8:: delx(1:nx)

  real*8::inxu(3,nx),inyu(3,ny),inxv(3,nx),inyv(3,ny),inbx(4,4),inby(4,4)
  real*8::cofivx(2,nx),cofiux(2,nx),cofivy(2,ny),cofiuy(2,ny)

  real*8::dcxu   (3,nx)
  real*8::dcyu   (3,ny+1)
  real*8::dcxv   (3,nx)
  real*8::dcyv   (3,ny)
  real*8::dcby   (4,4)
  real*8::dcbx   (5,4)
  real*8::cofcyu (2,ny)
  real*8::cofcyv (2,ny)
  real*8::cofcxu (2,nx)
  real*8::cofcxv (2,nx)
  real*8::vmagic(nx)

  real*8::vixu  (3,nx)
  real*8::viyu  (3,ny+1)
  real*8::vixv  (3,nx)
  real*8::viyv  (3,ny)
  real*8::vyui  (3,ny+1)
  real*8::vyvi  (3,ny)      
  real*8::cofvby(4,4)
  real*8::cofvbx(4,4)
  real*8::cofvyu(3,ny)
  real*8::cofvyv(3,ny)
  real*8::cofvxu(3,nx)
  real*8::cofvxv(3,nx)
  real*8::ccon  (2,4)
  real*8::rkc   (3)
  real*8::rkd   (3)
  real*8::rkcv  (3)
  real*8::rkdv  (3)
  real*8::signo (3)
  real*8::x     (0:nx+1)
  real*8::xr    (0:nx+1)
  real*8::y     (0:ny+1)

  real*8::dy    (0:ny)
  real*8::dx    (0:nx)

  real*8::idy(0:ny),idyy(1:ny)


  real*8,dimension(:)  ,allocatable:: vtop,vtop1,dvtdx,nyy
  real*8,dimension(:),allocatable::eta

  real*8::axp(3,nx)
  real*8::ayp(3,ny)
  real*8::phix(nx-1)
  real*8::phiy(ny-1)
  real*8::ayphi(3,ny-1)
  real*8::axphi(3,nx-1) 
  real*8::kmod(nx-1)	!real8
  real*8::dymin(ny)	!real8
  real*8:: v0(ny)

  real*8 tiempo,uav,ax,ay,az,re,rex
  integer nescr,nesta,nsubstps,nhist,nsteps,istart,stest
  integer len,reav,avcrt,flag,iflagr,zero,iflag,stats
  real*8 r,pi,times,tvcrt,Wz0,alp,bet,dtret,dxmin,dzmin

  ! clf
  integer cyclx, cycly
  real*8 cfl,idx,idxx
  logical*4 setstep
#ifdef PLANESPECTRA       
  logical*4 dostat(7) 
#else  
  logical*4 dostat
#endif
  ! ------------------------ aux rhsp ----------------------------------------!  
  real*8:: wki1r(nz+2),vm(ny),vmtmp(ny),wki2r(nz+2,0:7) !8=Max num of threads
end module alloc_dns

! -------------  things for rft -------------------------------!
module fourthings 
  complex*16, dimension(:), allocatable:: bdum,fdum
  integer*8  planf,planb
  real*8     dnf,dnb
  integer*4  nf,nb
end module fourthings
! ---------------------------------------------------------------!

module shared_mem
  implicit none
  type nodedata
     integer:: startpl, endpl, planes
     integer:: startpen, endpen, pencils
     integer:: size
     integer:: bound
     integer, dimension(:), pointer:: scount, sdisp, rcount, rdisp
  end type nodedata

  type domaindata
     integer:: NX, NY, NZ
     integer:: pencils, planes
     integer:: fatnodes, widenodes
     integer:: planesize, pencilsize
     integer, dimension(:), pointer:: pb, pe, ib, ie
  end type domaindata
contains
  subroutine domaindealloc(domain)
    implicit none
    type(domaindata), intent(inout):: domain
    deallocate(domain%ib)
    deallocate(domain%ie)
    deallocate(domain%pb)
    deallocate(domain%pe)
  end subroutine domaindealloc

  subroutine nodedealloc(node)
    implicit none
    type(nodedata), intent(inout):: node
    deallocate(node%scount)
    deallocate(node%sdisp)
    deallocate(node%rcount)
    deallocate(node%rdisp)
  end subroutine nodedealloc
end module shared_mem

! ---------------------------------------------------------------!
module point
  use omp_lib
  use shared_mem
  integer,dimension(:),allocatable :: ibeg,iend,pcibeg,pciend,pcibeg2,pciend2
  integer:: ib,ie,ib0,mmx,mpu,mpv,mpiout,ntotb,ntotv,ntot_corr,ntot_corr2
  integer:: mp_corr,mp_corr2,pcib,pcie,pcib2,pcie2
  ! ---- derived types for collective change.  ----
  type(nodedata):: nodeu,nodev,node_corr,node_corr2
  type(domaindata):: domainu,domainv,domain_corr,domain_corr2
 !--------Stuff for doing the OMP FFTw & Cosinus transform
  integer ompid,nthreads,chunk1,chunk2,chunk3
  !$OMP THREADPRIVATE(ompid,nthreads,chunk1,chunk2,chunk3)
  integer jbf1,jef1,jbf2,jef2,mpvb,mpve
  !$OMP THREADPRIVATE(jbf1,jef1,jbf2,jef2,mpvb,mpve)
end module point

! ---------------------------------------------------------------!
module statistics
  use ctesp,only: nx,ny,nz,nz2
  implicit none

  real*8,dimension(:,:),allocatable:: us,vs,ws,ua,va,wa,uv,vortx,vorty,vortz,&
       &vortxa,vortya,vortza,oxp,oyp,ozp,&
       &u3,v3,u2v,v2u,w2v,w2u,dispu,dispv,dispw,dispuv,&
       &pp,pm,pup,pvp,pdudx,pdvdy,pdudy,pdvdx,pdwdz,&
       &dudx0,dudy0,dudz0,dvdx0,dvdy0,dvdz0,dwdx0,dwdy0,dwdz0

#ifdef INFOINTER 
  real*8,dimension(:,:),allocatable::v_0,u_x0,u_xy0,w_0,w_y0,dwdx_0,dudz_x0,v_y0,dudx_0
#endif       
  real*8 ener(15)
  real*8,allocatable::i1c(:),i2c(:)
  real*8:: hy(0:ny)  
  integer, allocatable::flags(:),jspecy(:,:)!,jspecor(:) ! flags for velocity and vor spectra  
  ! Espectros de Velocidad y vorticidad !
  real*8,dimension(:,:,:),allocatable:: ensu,ensv,ensw,ensuv
  real*8,dimension(:,:,:),allocatable:: ensomz,ensomx,ensomy,pesp
  ! Correlaciones (x-x',j,z) 
  real*8,dimension(:,:),allocatable::coru ,corv ,corw ,coruv
  real*8,dimension(:,:),allocatable::corox,coroy,coroz,corp
#ifdef PLANESPECTRA
  real*8,dimension(:,:,:),allocatable::plane_specu,plane_specv,plane_specw
  integer totalcal(7)
#else
  integer totalcal
#endif

#ifdef PLANESPECTRA2
  complex*16,dimension(:,:,:),allocatable::planesv 
#endif

end module statistics
! ---------------------------------------------------------------!

module names
  implicit none
  character*100 stfile,etfile,vetfile,hfile,chfile,chinit,chinfo,chinfoext,corfile,budfile,spectraplane,spectraplane2  
  integer indst,ifile
end module names

!---------------------------------------------------------------------!
! In this module we define some of the arrays that will change their
! shape during the run
!
! SHC 20/06/06
!---------------------------------------------------------------------!
module main_val
  implicit none
  real*8,dimension(:),allocatable::u,v,w,p,rhsupa,rhsvpa,rhswpa,rhsu,rhsv,&
       &    rhsw,res,wki1,wki2,wki3,resv,resw,wkp,wkpo,bufuphy,buf_corr
end module main_val


!---------------------------------------------------------------------!
!Time parameters...
module temporal 
  real*8 tc1,tc2,iratec,tm1,tm2,iratem,tr1,tr2,irater
  real*8  ttotc,ttotm,tmp1,tmp2,tmp3,ttotr,ttri,ttrv,tci,tcv
  real*8  ttotinty,ttotvx,tmp4,tmp5,tmp6,ttotintx,ttotvy,ttotdy
  real*8  ttotfftc,ttotfft,ttotaux,ttotpois,ttotrhs,tmp19,tred
  real*8  tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,ttotim,ttotvdx,tmp13,tmp14
  real*8 tp1,tp2,th1,th2,tmp15,tmp16,tmpois,tmrhs,tmp17,tmp18,ttotgen,ttotbou
  real*8 tmp20,tmp21,tmp22,tmp23,tmp24,ttot1,ttot2,ttot3,ttot4,ttot5,ttot6,ttot7,ttot8
  real*8 ttot9,ttot10
  real*8 tmp25,tmp26,tmp27,tmp28,tmp29,tmp30
  real*8 tm3,tm4,tm5,tm6
end module temporal


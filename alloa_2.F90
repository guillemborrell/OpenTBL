!********************************************************************!
!							             !
!   This routine allocates the array dimensions                      !
!							             !
!                                         Mark Simens 24-5-2005      !
!  Based on routine by Javier Crespo                                 !
!********************************************************************!
subroutine alloa_2(mpiid)
  use statistics_2
  use alloc_dns_2
  use point_2
  use genmod_2
  use ctesp_2
  use temporal_2
  implicit none
  include 'mpif.h'
  real*8 memoryt1,mem_ux,mem_corr,mem_us,mem_ensu,mem_ensp,mem_sta
  integer mpiid,ierr
  if (mpiid.eq.0) then
     write(*,*) '=========== MASTER MPI ===================='
     write(*,'(a10,i5,a10,i5,a10,i5,a10,i5)') 'nx=',nx,'ny=',ny,'nz=',nz,'nz2=',nz2 
     write(*,*) '--------------------------------------------------------------------'
     write(*,'(a10,i5,a10,i5/a10,i5,a10,i5)') 'ie=',ie,'ib=',ib
     write(*,'(a10,i5,a10,i5/a10,i5,a10,i5)') 'Pen. mpu=',mpu,'Pen. mpv=',mpv
     write(*,'(a10,i5,a10,i5,a10,i5/a10,i5,a10,i5)') & 
          &  'nspec=',nspec,'ncorr=',ncorr,'lxcorr=',lxcorr
     write(*,*) '==============================================='

     memoryt1=8d0*(4*(nz2+1)+46*(nx+ny)+(ie-ib+1)*20*ny+ &
          & 8*(nz2+1)*nspec*lxp+(10*ntot_corr2*lxcorr)+&
          & 46*(ny+1)*(ie-ib+1))/1024**2
   
     totmem=totmem+memoryt1
     write(*,*) '-------------------------------------------------------'
     write(*,'(a35,f10.2,a4)') 'Allocated Memory Alloa:', memoryt1, ' MB'
     write(*,*) 
     write(*,*) '-------------------------------------------------------'
  endif

  kaz    =0d0;   kaz2=0d0;    kazr=0d0;   kvis=0d0
  inxu   =0d0;   inyu=0d0;    inxv=0d0;   inyv=0d0;  inbx=0d0;  inby=0d0
  cofivx =0d0;cofiux =0d0;cofivy =0d0; cofiuy=0d0 
  dcxu   =0d0;   dcyu=0d0;    dcxv=0d0;   dcyv=0d0;  dcby=0d0;  dcbx=0d0
  cofcyu =0d0; cofcyv=0d0;  cofcxu=0d0; cofcxv=0d0
  vixu   =0d0;    viyu=0d0;    vixv=0d0;   viyv=0d0;  vyui=0d0;  vyvi=0d0
  cofvby =0d0;  cofvbx=0d0;  cofvyu=0d0; cofvyv=0d0;cofvxu=0d0;cofvxv=0d0;ccon=0d0
  rkc    =0d0;    rkd=0d0;    rkcv=0d0;   rkdv=0d0 
  signo  =0d0
        x=0d0;      xr=0d0;       y=0d0;     dy=0d0;    dx=0d0;    hy=0d0 
      axp=0d0;     ayp=0d0;    phix=0d0;   phiy=0d0; ayphi=0d0; axphi=0d0
  vmagic =0d0;      v0=0d0;

  !TIMERS
  ttotc=0d0; ttotinty=0d0; ttotvy=0d0;ttot1=0d0;ttot2=0d0;ttot3=0d0;ttotpois=0d0;ttotrhs=0d0
  ttotm=0d0; ttotvx=0d0; ttotdy=0d0;ttotfft=0d0;ttotfftc=0d0;tmpois=0d0;tmrhs=0d0
  ttotr=0d0; ttotintx=0d0; ttotvdx=0d0;ttotim=0d0; ttotaux=0d0;totmem=0d0;ttotbou=0d0;ttotgen=0d0
  tred=0d0;
  ttot4=0d0;ttot5=0d0;ttot6=0d0;ttot7=0d0;ttot8=0d0;ttot9=0d0;ttot10=0d0
  tmp1  =0d0; tmp4=0d0; tmp7=0d0; tmp10=0d0; tmp13=0d0;tmp16=0d0;tmp19=0d0;tmp22=0d0;tmp25=0d0;tmp28=0d0
  tmp2  =0d0; tmp5=0d0; tmp8=0d0; tmp11=0d0; tmp14=0d0;tmp17=0d0;tmp20=0d0;tmp23=0d0;tmp26=0d0;tmp29=0d0
  tmp3  =0d0; tmp6=0d0; tmp9=0d0; tmp12=0d0; tmp15=0d0;tmp18=0d0;tmp21=0d0;tmp24=0d0;tmp27=0d0;tmp30=0d0

!--Reynolds Stress Tensor ====================================
  allocate(us(ny+1,ib:ie))	!real8
  allocate(vs(ny  ,ib:ie))
  allocate(ws(ny+1,ib:ie))
  allocate(ua(ny+1,ib:ie))
  allocate(va(ny  ,ib:ie))
  allocate(wa(ny+1,ib:ie))
  allocate(uv(ny+1,ib:ie)) 
  allocate(vw(ny+1,ib:ie)) 
  allocate(uw(ny+1,ib:ie))  
  allocate(pp(ny,ib:ie))
  allocate(pm(ny,ib:ie))
  us =0d0;vs =0d0;ws =0d0;ua =0d0;wa =0d0;va =0d0;uw=0d0;vw=0d0;
  uv=0d0;pp=0d0;pm=0d0
!--Vorticities: mean and fluctuations---------------------------  
  allocate(vortx(ny+1,ib:ie))
  allocate(vorty(ny+1,ib:ie))
  allocate(vortz(ny+1,ib:ie))
  allocate(vortxa(ny+1,ib:ie))
  allocate(vortya(ny+1,ib:ie))
  allocate(vortza(ny+1,ib:ie))
  vortx=0d0;vorty=0d0;
  vortz=0d0;vortxa=0d0;vortya=0d0;vortza=0d0
!--TRIPLE PRODUCTS-----------------------------------------------
  allocate(u3(ny+1 ,ib:ie))
  allocate(v3(ny   ,ib:ie))
  allocate(u2v(ny  ,ib:ie))
  allocate(v2u(ny  ,ib:ie))
  allocate(w2v(ny  ,ib:ie))
  allocate(w2u(ny+1,ib:ie))
  u3=0d0;v3=0d0;u2v=0d0;v2u=0d0;w2v=0d0;w2u=0d0 
  
   

!*********STATISTIC IN FOURIER SPACE ********************************
!velocity spectra in z (Reynolds stresses) 
  allocate(ensu (0:nz2,nspec,lxp))
  allocate(ensv (0:nz2,nspec,lxp))
  allocate(ensw (0:nz2,nspec,lxp))
  allocate(ensuv(0:nz2,nspec,lxp))
  ensu=0d0;ensv=0d0;ensw=0d0;ensuv=0d0
 !vorticities spectra  
  allocate(ensomz (0:nz2,nspec,lxp))
  allocate(ensomx (0:nz2,nspec,lxp))
  allocate(ensomy (0:nz2,nspec,lxp))
  allocate(pesp   (0:nz2,nspec,lxp))
  ensomz =0d0;ensomx =0d0;ensomy =0d0;pesp=0d0
  
allocate(flags(nx)) !integer

!*********STATISTIC IN REAL SPACE *****
!vorticities ===================================
allocate(oxp(ny,ib:ie)) !real8
allocate(oyp(ny,ib:ie))
allocate(ozp(ny,ib:ie))
oxp=0d0;oyp=0d0;ozp=0d0

!dissipation ===================================
allocate(dispu (ny,ib:ie)) !real8
allocate(dispv (ny,ib:ie))
allocate(dispw (ny,ib:ie))
allocate(dispuv(ny,ib:ie))
dispu=0d0;dispv=0d0;dispw=0d0;dispuv=0d0;

!pressure statistics: p*div(v), p*grad(v)=======
allocate(pup(ny,ib:ie))
allocate(pvp(ny,ib:ie))
allocate(pdudx(ny,ib:ie))
allocate(pdvdy(ny,ib:ie))
allocate(pdudy(ny,ib:ie))
allocate(pdvdx(ny,ib:ie))
allocate(pdwdz(ny,ib:ie))	
pvp=0d0;pup=0d0;pdudx=0d0;pdudy=0d0;pdvdx=0d0
pdvdy=0d0;pdwdz=0d0

!K=0 for derivatives:
allocate(dudx0(ny,ib:ie))
allocate(dudy0(ny,ib:ie))
allocate(dudz0(ny,ib:ie))
allocate(dvdx0(ny,ib:ie))
allocate(dvdy0(ny,ib:ie))
allocate(dvdz0(ny,ib:ie))
allocate(dwdx0(ny,ib:ie))
allocate(dwdy0(ny,ib:ie))
allocate(dwdz0(ny,ib:ie))
dudx0=0d0;dudy0=0d0;dudz0=0d0;
dvdx0=0d0;dvdy0=0d0;dvdz0=0d0;
dwdx0=0d0;dwdy0=0d0;dwdz0=0d0;

#ifndef NOCORR 
!Correlaciones ===================================
!               1365   30     3      3   = 2.81 Mb
 allocate(coru (ntot_corr2,lxcorr))	!real8  
 allocate(corv (ntot_corr2,lxcorr))
 allocate(corw (ntot_corr2,lxcorr))
 allocate(coruv(ntot_corr2,lxcorr))
 allocate(coruw(ntot_corr2,lxcorr))
 allocate(corvw(ntot_corr2,lxcorr))
 coru=0d0;corv=0d0;corw=0d0;coruv=0d0;coruw=0d0;corvw=0d0;

 allocate(corox(ntot_corr2,lxcorr))
 allocate(coroy(ntot_corr2,lxcorr))  
 allocate(coroz(ntot_corr2,lxcorr))
 allocate(corp (ntot_corr2,lxcorr))	!real8 
 corox=0d0;coroy=0d0;coroz=0d0;corp=0d0;
#endif

#ifdef INFOINTER 
!  &,v_0,u_x0,u_xy0,w_0,w_y0,dwdx_0,dudz_x0,v_y0,dudx_0
allocate(v_0    (ny,ib:ie))
allocate(u_x0   (ny,ib:ie))
allocate(u_xy0  (ny,ib:ie))
allocate(w_0    (ny,ib:ie))
allocate(w_y0   (ny,ib:ie))
allocate(dwdx_0 (ny,ib:ie))
allocate(dudz_x0(ny,ib:ie))
allocate(v_y0   (ny,ib:ie))
allocate(dudx_0 (ny,ib:ie))
v_0=0d0;u_x0=0d0;u_xy0=0d0;w_0=0d0;w_y0=0d0;
dwdx_0=0d0;dudz_x0=0d0;v_y0=0d0;dudx_0=0d0
#endif


! genflu  ===================================
 cf=0d0;cf2=0d0;cfinfo=0d0
 counter=0
end subroutine alloa_2
!---------------------------------------------------------------------!
!
! In this module we define some of the arrays that will change their
! shape during the run.
!
!  do first only things read from file, to allow for reading buffers
!
! SHC 20/06/06,     JJS 01/10
!                                                                     !  
!---------------------------------------------------------------------!

subroutine alloavar_2(ntotb,ntotv,mpiid)
  use main_val_2
  use ctesp_2

  implicit none
  include 'mpif.h'

  integer:: ntotb,ntotv,mpiid,ierr
  real*8 memoryt2,mem_u,mem_v

  !Computing memory requirements:-------------------------------

  if (mpiid.eq.0) then
     write(*,*) '==============================================='
     write(*,'(a25,f10.2,a25,f10.2)') 'ntotb (Mpoints)',real(ntotb/1e6),'ntotv (Mpoints)',real(ntotv/1e6)
     write(*,*) '==============================================='

     memoryt2=(2d0*ntotb+2d0*ntotv)*8/1024**2
     mem_u=ntotb*8d0/1024**2
     mem_v=ntotv*8d0/1024**2

     totmem=totmem+memoryt2
     write(*,*) '-------------------------------------------------------'
     write(*,'(a35,f10.2,a4)') 'Allocated Memory Alloavar=', memoryt2,' MB'
     write(*,*)
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '2 copies: Memoria u(ntotb)=', mem_u,' MB',' total=',2*mem_u,' MB'
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '2 copies: Memoria v(ntotv)=', mem_v,' MB',' total=',2*mem_v,' MB'
     write(*,*) '-------------------------------------------------------'
  endif
  !-----------------------------------------------------------
  
  allocate(u(ntotb)) !real8
  allocate(w(ntotb))
  allocate(v(ntotv)) !real8
  allocate(p(ntotv))

  u=0d0;w=0d0;v=0d0;p=0d0
end subroutine alloavar_2


!---------------------------------------------------------------------!
!
! In this module we define some of the arrays that will change their
! shape during the run
!
!                        buffers
! 
! SHC 20/06/06,    JJs 01/10
!                                                                     !  
!---------------------------------------------------------------------!

subroutine alloabuff_2(ntotb,ntotv,ntot_corr,mpiid)
  use main_val_2
  use ctesp_2

  implicit none
  include 'mpif.h'

  integer:: ntotb,ntotv,ntot_corr,mpiid,ierr
  real*8 memoryt3,mem_u,mem_v,mem_2,mem_3,mem_c

  if (mpiid.eq.0) then   

     memoryt3=(7*ntot_corr+8d0*ntotb+4d0*ntotv+3*(nz+2)*(ny+1)+2*(nz2+1)*(ny+1))*8/1024**2
     mem_u=ntotb*8d0/1024**2
     mem_v=ntotv*8d0/1024**2
     mem_c=7*lxcorr*ntot_corr*8d0/1024**2 !7 buffers of size ntot_corr*lxcorr
     mem_2=(nz+2)*(ny+1)*8d0/1024**2
     mem_3=2*(nz2+1)*(ny+1)*8d0/1024**2

     totmem=totmem+memoryt3

     write(*,*) '-------------------------------------------------------'
     write(*,'(a35,f10.2,a4)') 'Allocated Memory AlloaBuff=', memoryt3,' MB'
     write(*,*)
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '8 copies: Memoria u(ntotb)=', mem_u,' MB',' total=',8*mem_u,' MB'
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '4 copies: Memoria v(ntotv)=', mem_v,' MB',' total=',4*mem_v,' MB'
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '7 copies: Mem. Corr=', mem_c/7d0,' MB',' total=',mem_c,' MB'
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '3 copies: ZY Physic Plane=', mem_2,' MB',' total=',3*mem_2,' MB'
     write(*,'(a35,f10.2,a4,a10,f10.2,a4)') '2 copies: ZY Complx Plane=', mem_3,' MB',' total=',2*mem_3,' MB'
     write(*,*) '-------------------------------------------------------'
  endif

  !-----------------------------------------------------------
  allocate(rhsupa(ntotb)) !rhsu paso anterior 'pa'
  allocate(rhswpa(ntotb))
  allocate(rhsu  (ntotb))
  allocate(rhsw  (ntotb))
  allocate(res   (ntotb))
  allocate(resw  (ntotb))
  allocate(wki1  (ntotb))
  allocate(wki3  (ntotb)) !real8
  allocate(wkp   ((nz+2)*(ny+1)))
  allocate(wkpo  ((nz+2)*(ny+1)))		  
  !CORRELATIONS:
  allocate(buf_corr(7*lxcorr*ntot_corr))  !ntot_corr=max[(nz1*ncorr)*planes,pencils*nx ]
 
!!!!Extra buffer for triple products
  allocate(bufuphy   ((nz+2)*(ny+1)))    
  allocate(rhsvpa(ntotv))
  allocate(rhsv  (ntotv))
  allocate(resv  (ntotv))
  allocate(wki2  (ntotv)) !real8
  buf_corr=0d0

  rhsupa  =0d0
  rhswpa  =0d0 
  rhsu    =0d0
  rhsw    =0d0 
  res     =0d0
  resw    =0d0
  wki1    =0d0 
  wki3    =0d0
  wkp     =0d0
  wkpo    =0d0
  bufuphy =0d0
  rhsvpa  =0d0 
  rhsv    =0d0
  resv    =0d0
  wki2    =0d0 
end subroutine alloabuff_2


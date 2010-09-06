!----------------------------------------------------------------------------*
! This code is based on compact finite differences in x and y on non-uniform
! grids. The grid is staggered in x and y but at the moment in z it is collo
! cated. The fractional step method used is based on the method presented in
! Perot JCP 108 Analysis of the Fractional Step Method.
! Time integration is a third order Runge-Kutta scheme based on a paper by
! Spalart,Rogers&Moser JCP 96. The method is explicit for ALL terms in x and
! z but the viscous terms in the y-direction are implicit. 
! The pressure poisson will be solved using multigrid.
! This routine is an interface to viscop and pois:
!          viscop is the interface to the routines which calculate all
!                                    derivatives
!          pois   contains the routines which assures mass conservation.     
!          est    write statistics
!          escp   write restart files
!
! Mark Simens 22-4-2005
! Improvements by G.Hauet, J Jimenez, and myself.
! New Version, Sergio Hoyas.
!--------------------------------------------------------------------------*

! ==========================================================================
!                              AAACCCHHHTTTUUUNNNGGG
!   To correct a serious miskeying of p with the rest of the variables
!     that was now becoming very inconvenient,  p HAS BEEN CHANGED
!     p(i=nx) had no meaning, while p(i=1) was meaninful,
!     that was the oposite to everything else, for which
!     i=1 is the inflow condition, which is not integrated 
!        SAME with j
!     To have corresponding indices, I have shifted p by one
!         p(2:ny,2:nx)= p(1:ny,1:nx-1)
!         p(i=1) = 0
!         p(j=1) = 0
!     This version reads new I/O files with p in the new indexing
!     and written in cross (zy) planes. It also uses changes 
!     that are completely independent of xy planes.
!               JJS/01/10
! ==========================================================================

program capalimite

  use alloc_dns
  use main_val
  use names
  use temporal
  use point
  use statistics,only: ener
  use genmod,only: rthin,din,u0
  use ctesp  
  use omp_lib

  implicit none
  include "mpif.h"


  real*8 dt,vardt
  integer isubstp,istep,ical,ierr,mpiid,numprocs
  logical:: vcontrol

  !----------------------------------------------------------------------*
  !       el proceso maestro llama a las rut. de inic. del codigo
  !----------------------------------------------------------------------*
  vcontrol=.false. !just checking correct time step
  
  ! Medimos tiempo ! 
  tiempo=0d0
  flag = 1
  iflagr =1
  zero=0
  times=0d0
  pi=4d0*atan(1d0)
  tvcrt=2d0*pi/(16*160)

  !       /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  !$ CALL OMP_SET_DYNAMIC(.FALSE.)
  mpiid2=mpiid

  !Then number of processors used when compilling is assigned to th program:
  nummpi=numprocs; pnodes=numprocs

  if (mpiid.eq.0) ical=0   
  ! reads hrem.dat.
  ! creates coordinates for parallel pros's.
  if (mpiid.eq.0) write(*,*) '======== PROGRAM BEGINS ==============='
  call  iniciap(mpiid)
  

  ! reads/generates inlet conditions initiates arrays (just in case)
  ! warning!! getstart reads y from file !!!

  if (mpiid .eq. 0) write(*,*) '======== CALLING GETSTART ==============='
  if (mpiid2.eq.0) th2 = MPI_WTIME()                
  call getstartzy(u,v,w,p,dt,mpiid)
  if (mpiid2 .eq. 0) then  
     th1 = MPI_WTIME()      
     tmp29 =tmp29+abs(th2-th1)
     write(*,*) 
     write(*,'(a40,f15.3,a3)') '===== FIELD (u,v,w,p) read in: =======',tmp29,'Sg' 
     write(*,*)    
  endif
   
  call coef(mpiid)
  call inlet_retheta(u0,rthin,din,mpiid) !compute Rth_inlet and d99_inlet
  call magicnumber(mpiid)

#ifdef CREATEPROFILES
  paso=-1 !Substeps counter
  if(mpiid.eq.0) open(17,file='extraprofiles',form='unformatted',status='unknown')
  open(27,file='plano_div.dat',form='unformatted',status='unknown')
  open(28,file='zero-mode.dat',form='unformatted',status='unknown')  
  call create_profiles(u,v,w,rthin,mpiid)
  call impose_profiles(u,v,w,mpiid) 
  close(17)  
#endif

  call alloabuff(ntotb,ntotv,ntot_corr,mpiid)

  rhsupa = 0d0
  rhswpa = 0d0
  rhsvpa = 0d0

#ifdef CREATEPROFILES
   if(mpiid.eq.0) write(*,*) 'Checking Divergence of the Field:'
   call check_divergence(u,v,w,rhsupa,mpiid)
#endif

  call mpi_barrier(mpi_comm_world,ierr)

#ifdef FLOPS
  CALL HPM_INIT()
  CALL HPM_START('WORK1')
#endif

#ifndef NOINFOSTEP
!Genflu info:
if(mpiid.eq.0) open(36,file=chinfoext,form='formatted',status='unknown')
#endif

  do istep = 1,nsteps
     if(mpiid.eq.0) then
        tc1 = MPI_WTIME()
        write(*,'(a60,i6)') '....................................................istep=',istep
     endif

     if (.TRUE.) setstep=.TRUE.
     if (mod(istep,avcrt)==0.and.istep>stest) dostat=.TRUE.

     !.............................................  
     do isubstp = 1,3
        if (mpiid2 .eq. 0) th1 = MPI_WTIME()               
        call mpi_barrier(mpi_comm_world,ierr)   !************************

        call rhsp(u,v,w,p,rhsupa,rhsvpa,rhswpa,       &
             &    res,res,resv,resv,resw,resw,        & 
             &    rhsu,rhsv,rhsw,                     &
             &    wki1,wki1,wki2,wki2,wki3,wki3,      &
             &    wkp,wkp,wkpo,wkpo,bufuphy,buf_corr, & 
             &    dt,isubstp,ical,istep,mpiid)

        call mpi_barrier(mpi_comm_world,ierr)    !********************
        if (mpiid2 .eq. 0) then  
           th2 = MPI_WTIME()      
           tmp13 =tmp13+abs(th2-th1)
        endif

        call outflow_correction(u,v,rhsupa)
        
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp18 =tmp18+abs(th2-th1)
        endif
        call mpi_barrier(mpi_comm_world,ierr)   !************************

        vardt = 5d-1/dt/rkdv(isubstp)   

      
        call pois(u,v,w,p,res,res,resw,vardt,mpiid)

#ifdef CREATEPROFILES      
call check_divergence(u,v,w,res,mpiid)
#endif 
        if (mpiid2 .eq. 0) then  
           th2 = MPI_WTIME()      
           tmp14 =tmp14+abs(th2-th1)
        endif
     enddo
     !.............................................
     call mpi_barrier(mpi_comm_world,ierr)
     tiempo = tiempo+dt
     if (times .ge. tvcrt) then
        times=times-tvcrt
     endif
     times=times+dt

     ! I/O Operations --------------------------------
     if (mod(istep,stats).eq.0) then
        if (mpiid2.eq.0) th2 = MPI_WTIME()                                                     
        call escrst(ax,ay,az,cfl,tiempo,re,x,y,mpiid,ical)
        ical=0
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp27 =tmp27+abs(th2-th1)
        endif
     endif

     if (mod(istep,reav).eq.0) then
        if (mpiid2.eq.0) th2 = MPI_WTIME() 
        call escribezy(u,v,w,p,dt,mpiid)   !Vul & BG     
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp28 =tmp28+abs(th2-th1)
        endif
     endif

     if (mpiid .eq. 0) then  ! Tiempos 
        tc2 = MPI_WTIME()      
        tmp1 = tc2-tc1     
        call summary1(istep,dt,vcontrol)             
     endif
#ifdef CHECKTIME
     call MPI_BCAST(vcontrol,1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
     IF(vcontrol) stop
#endif
  enddo


#ifdef FLOPS
  CALL HPM_STOP('WORK1')
  call hpm_print()
  call hpm_print_flops()
  call hpm_print_flops_agg()
#endif 


#ifndef NOINFOSTEP
close(36)
#endif

#ifdef CREATEPROFILES 
close(27)
close(28)
#endif


  if (mpiid.eq.0) then
     call summary2()
  end if

  call mpi_barrier(mpi_comm_world,ierr)
  if (mpiid.eq.0) then
     write(*,*)
     write (*,*) '=============================================='
     write (*,*) 'FINALIZING ALL THE PROCESSES: PROGRAM DONE'
     write (*,*) '=============================================='
  endif
  
  call MPI_FINALIZE(ierr)
  stop
end program capalimite


! ---------------------------------------------------------------------!
! ----------------------- Subroutines ---------------------------------!
! ---------------------------------------------------------------------!

!----------------------------------------------------------------------*
!        lee/genera condiciones iniciales y de contorno
!        inicializa variables
!----------------------------------------------------------------------*

subroutine iniciap(mpiid)
  use alloc_dns
  use names
  use point
  use shared_mem
  use genmod
  use statistics
  use ctesp
  use omp_lib
  implicit none
  include "mpif.h"

  integer mpiid,k,blocks,elem,stride,iproc,fcr(1:nz1),i,ii,j
  integer posx,posy(1:npos),ierr,idat(20),i3

  real(8)       dat(6)
  character(99) text
  character(3) ext,ext2


  !----------------------------------------------------------------------!
  !      el  maestro lee parametros iniciales y los reparte
  !----------------------------------------------------------------------!

  if (mpiid .eq. 0) then
     open(19,file='hrem.dat',status='old')   
     call avanza(text)
     !                  Re    ,ax    ,ay    ,az    ,r  
     read(text,*) dat(1),dat(2),dat(3),dat(4),dat(5) 
     call avanza(text)
     !                 cfl   ,istart
     read(text,*) dat(6),idat(5)   
     call avanza(text)
     write(*,*) text
     !                 nsteps,   reav,    stats   avcrt,  stest
     read(text,*) idat(1), idat(13), stats, idat(14),idat(15)
     call avanza(text)
     !           ifile1
     read(text,*) idat(7)
     write(*,*) 'idat',idat(7)
     call avanza(text)
     read(text,*) idat(8),idat(9),idat(10),idat(11),idat(12)
     ! Genflu flags 
     call avanza(text)
     read(text,*) flaggen,utauout,rthout,dout
     call avanza(text)
     read(text,'(a100)') chinit
     call avanza(text)
     read(text,'(a100)') chfile
     call avanza(text)
     read(text,'(a100)') chinfo
     close (19)
  endif

  call MPI_BCAST(dat    ,6 ,MPI_REAL8    ,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(idat   ,19,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stats  ,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(flaggen,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(chfile,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  !for the Parallel IO Writting
  call MPI_BCAST(chinit,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  !for the Parallel IO Reading
  re     = dat(1)
  ax     = dat(2)
  ay     = dat(3)
  az     = dat(4)
  r      = dat(5)
  cfl    = dat(6)

  nsteps = idat(1)
  istart = idat(5)       
  ifile  = idat(7)
  posx   = idat(8)
  posy(1:npos) = idat(9:(npos+8))
  reav = idat(13)
  avcrt= idat(14)
  stest= idat(15)

  write(ext,'(i3.3)') ifile-1
  write(ext2,'(i3.3)') ifile
  if (mpiid .eq. 0) then
     !stfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.st'
     etfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.esp'
     vetfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.vesp'
     chinfoext=chinfo(1:index(chinfo,' ')-1)//'.'//ext2//'.dat'  !File for info at each step
     !      hfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.cf'
     !      open(31,file=hfile)
     !      hfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.sc'  
     !      open(102,file=hfile)
  endif

  tiempo = 0d0
  flag   = 1
  iflagr = 1
  zero   = 0
  times  = 0d0
  pi     = 4d0*atan(1d0)
  tvcrt  = 2d0*pi/(16d0*160d0)
  totalcal=0

  !----------------------------------------------------------------------!
  !           inicializa los punteros
  !----------------------------------------------------------------------!
  allocate(ibeg   (0:nummpi-1),iend   (0:nummpi-1))
  allocate(pcibeg (0:nummpi-1),pciend (0:nummpi-1))
  allocate(pcibeg2(0:nummpi-1),pciend2(0:nummpi-1))

  call pointers_p2p(mpiid)   !! pointers et al. for changep2p & total sizes
  call pointersfft(mpiid) !pointers for thread_private indexes (used in FFTw)

  ib0 = max(2,ib)  ! -- i=1 is inflow, not integrated
  call alloa(mpiid)

! --  found out where is the reference plane: 
  mpiout=-100  
  do i=0,nummpi-1
     if (ibeg(i).le.xout .and. iend(i).ge.xout) then
        mpiout=i
        exit
     endif     
  enddo
  if(mpiout<0) then
     if(mpiid.eq.0) write(*,*) 'xout ', xout,' out of bounds, stopping'
     stop
  endif

  if (mpiid==0) write(*,*) 'ALOCATANDO...................................'
  call alloavar(ntotb,ntotv,mpiid)
  call rfti(nz)
  call cosfti(nx-1)

    
#ifndef NOSPECTRA
!======================SPECTRA
  flags=0
  do i=1,lxp    
     do ii=-nxp(i),nxp(i)
        flags(xpoint(i)+ii)=i
     enddo
  enddo

  ! jspecy data
  allocate(jspecy(nspec,lxp))
  !Readint table of spectra from file:
  if(mpiid.eq.0) then
     open(100,file='tablespectra',form='unformatted',status='unknown')
     do j=1,nspec
        read(100) jspecy(j,1:lxp)
     enddo
     write(*,*) 'valores de jspecy'
     do j=1,nspec
      write(*,'(10i6)') jspecy(j,1:lxp)
     enddo
     close(100)
  endif
  call MPI_BCAST(jspecy,size(jspecy),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif


  ! -------------  write header for output -------------
  if (mpiid==0) then
     write(*,*)
     write(*,*)
     write(*,'(a7,f10.3,a8,f8.3,a8,f8.3,a8,f8.3)') &
          &       'Re =',re, '  lon =',ax,' alt =',ay,'wid=',az
     write(*,*)
     write(*,'(a8,i5,a8,i5,a8,i5)')&
          &           'nx =',nx,'ny =',ny,'nz=',nz
     write(*,*)
     write(*,'(a9,f10.4)')'vel =',r
     write(*,*)
     write(*,'(a8,f5.2)') 'CFL =',cfl
     write(*,*)
     write(*,'(a10,i6,a9,i6,a9,i5)') &
          &           'nstep =',nsteps,'nescr =',reav,'nhist =',avcrt
     write(*,*) 
     write(*,'(a10,i2/a9,a41/a12,a41)') &
          &     'istart =',istart,' from:',chinit,'write in :',chfile
     write(*,*)
  endif
  ! ---  computes wavenumbers and fourier indices and interpolation factor
  do k=0,nz2
     fcr(2*k+1) = k
     fcr(2*k+2) = k
     kaz(k)  = cmplx(0d0 ,k/az)
     kazr(k) = k/az
  enddo

  do k=0,nz2 ! this is positive but will be made negative for viscous terms
     kaz2(k) = -kaz(k)**2
  enddo
endsubroutine iniciap

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------
subroutine inlet_retheta(u0,rthin,din,mpiid)
use ctesp,only:ny
use alloc_dns,only:re,dy,y
implicit none
include "mpif.h"
real*8:: u0(ny+1),rthin,din,uinf,u99
integer:: j,mpiid,jtop,ierr
if(mpiid.eq.0)then
  uinf=sum(u0(ny+1-20:ny+1-5))/16d0
  u99 = 0.99*uinf
     do j=1,ny-2
        if (u0(j).ge.u99) exit
        jtop = j+1
     enddo
     din = y(jtop)
     rthin = sum((uinf-u0(2:jtop+3))*u0(2:jtop+3)*dy(1:jtop+2))*re/uinf
     write(*,*) '============================================='
     write(*,'(a20,f10.4,a10,f10.4)') 'Re_theta_inlet=',rthin,'d99_inlet=',din
     write(*,*) '============================================='
endif
CALL MPI_BCAST(rthin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(din,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

endsubroutine inlet_retheta


!  =======================================
!    auxiliary to read input param. 
!  =======================================
subroutine avanza(text)
  implicit none
  character*99 text
  do
     read(19,'(a)') text
     !   write(*,*) text
     if(text(1:2)/='CC') exit
  enddo
endsubroutine avanza

!----------------------------------------------
!----------------------------------------------
!----------------------------------------------

subroutine summary1(istep,dt,vcontrol)
  use temporal
  use genmod
  use alloc_dns
  implicit none
  logical:: vcontrol
  real*8::dt
  integer::istep,ierr
  ttotm = ttotm+tmp1
  ttotc = ttotc+tmp2
  ttotr = ttotr+tmp3
  ttotinty = ttotinty+ tmp4
  ttotvdx  = ttotvdx+   tmp5
  ttotintx = ttotintx+ tmp6
  ttotvy   = ttotvy+   tmp7
  ttotdy   = ttotdy+   tmp8
  ttotim   = ttotim+   tmp9
  ttotfft   = ttotfft+   tmp11
  ttotfftc   = ttotfftc+   tmp10
  ttotaux = ttotaux+tmp12
  ttotpois=ttotpois+tmp14
  ttotrhs=ttotrhs+tmp13
  tmpois=tmpois+tmp15
  tmrhs=tmrhs+tmp16
  ttotbou=ttotbou+tmp18
  ttotgen=ttotgen+tmp17
  !                 tred   =tred+tmp19
  ttot1=ttot1+tmp20
  ttot2=ttot2+tmp21
  ttot3=ttot3+tmp22
  ttot4=ttot4+tmp23
  ttot5=ttot5+tmp24
  ttot6=ttot6+tmp25
  ttot7=ttot7+tmp26
  ttot8=ttot8+tmp27
  ttot9=ttot9+tmp28

  write(*,*)
  write(*,'(a10,i8,2e15.6,a30,2f10.4)')  'step',istep, dt,tiempo,'Re_theta: in/out-ref',rthin,rthout
  write(*,'(a35,3f10.4)')  'tiempos: Trans,Comm,Total Step',tmp3,tmp2,tmp1-tmp27-tmp28
  write(*,'(a35,2f10.4)')  'ffts: fft, cos', tmp11,tmp10
#ifdef CHECKTIME
!!!!!!!!!!!!!!!! ONLY FOR 1 PLANE PER NODE !!!!!!!!!!!! CHECKING THE CORRECT TIMING
  if(tmp1-tmp27-tmp28.gt.20) then
  WRITE(*,*) '======== EXCESIVE TIME ============== NOW stopping'
  vcontrol=.true.
  endif
#endif
  write(*,*)
  !         write(31,'(20d22.14)') tiempo,dt,ener(1:15),tmp2,tmp1
  !         call flush(31)
  tmp1  =0d0; tmp4=0d0; tmp7=0d0; tmp10=0d0; tmp13=0d0;tmp16=0d0;tmp19=0d0;tmp22=0d0;tmp25=0d0;tmp28=0d0
  tmp2  =0d0; tmp5=0d0; tmp8=0d0; tmp11=0d0; tmp14=0d0;tmp17=0d0;tmp20=0d0;tmp23=0d0;tmp26=0d0;!tmp29=0d0
  tmp3  =0d0; tmp6=0d0; tmp9=0d0; tmp12=0d0; tmp15=0d0;tmp18=0d0;tmp21=0d0;tmp24=0d0;tmp27=0d0;tmp30=0d0     
  ! END I/O Operations --------------------------------
end subroutine summary1

subroutine summary2()
  use temporal
  use ctesp
  use point
  implicit none
  ttotm=ttotm+tmp29 !tmp29: reading time
  
  write(*,*) '==============================================================='
  write(*,'(a30,1f10.0,a2)') 'TOTAL MEM ALLOCATED (node):',totmem,'Mb'
  write(*,*) 'Size, Nodes, Threads',nx,ny,nz,nummpi,nthreads
  write(*,*) '---------------------------------------------------------------'
  write(*,'(a30,f15.4)') 'tiempo total: ', ttotm
  write(*,*) '---------------------------------------------------------------'
  write(*,'(a30,2f15.4,a2)') 'Reading Fields (u,v,w,p)', tmp29, tmp29/ttotm*100,'%'  
  write(*,'(a30,2f15.4,a2)') 'Writing Fields (u,v,w,p)', ttot9, ttot9/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Writing Statistics & Spectra', ttot8, ttot8/ttotm*100,'%'
  write(*,*) '---------------------------------------------------------------'     
  write(*,'(a30,2f15.4,a2)') 'Poisson Subroutine', ttotpois, ttotpois/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'RHS Subroutine'    , ttotrhs, ttotrhs/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Outflow Subroutine'  , ttotbou, ttotbou/ttotm*100,'%'

  !      write(*,'(a30,2f15.4,a2)') 'Genflu Subroutine' , ttotgen, ttotgen/ttotm*100,'%'
  write(*,*) '---------------------------------------------------------------'
  write(*,'(a30,2f15.4,a2)') 'transpuestas  '     , ttotr, ttotr/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'comm Changes       ', ttotc, ttotc/ttotm*100,'%'     
  write(*,'(a30,2f15.4,a2)') 'Poisson Com 1st S/R' , tmpois, tmpois/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Poisson Com 2nd S/R' , ttot4, ttot4/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'RHS Com'     , tmrhs, tmrhs/ttotm*100,'%'
  !      write(*,'(a30,2f15.4,a2)') 'Time Inside Bound', tred, tred/ttotm*100,'%'
  write(*,*) '----------------------------------------------------------------'
  write(*,'(a30,2f15.4,a2)') 'Inside genflu Master', ttot2, ttot2/ttotm*100,'%'
  !     write(*,'(a30,2f15.4,a2)') 'Genflu inside max', ttot3, ttot3/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Inside genflu Send/Rec', ttot5, ttot5/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Inside genflu Inter Loop', ttot6, ttot6/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'Inside genflu K     Loop', ttot7, ttot7/ttotm*100,'%'

  write(*,'(a30,2f15.4,a2)') 'Time Step Allreduce', ttot1, ttot1/ttotm*100,'%'
  write(*,*) '----------------------------------------------------------------'
  write(*,'(a30,2f15.4,a2)') 'interpy       ', ttotinty, ttotinty/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'interpx       ', ttotintx, ttotintx/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'T.Implici     ', ttotim, ttotim/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'diff_y        ', ttotdy, ttotdy/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'vis_yy        ', ttotvy, ttotvy/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'vis_xx+dif_x  ', ttotvdx,ttotvdx/ttotm*100,'%'
  write(*,*) '----------------------------------------------------------------'
  write(*,'(a30,2f15.4,a2)') 'cos transform ', ttotfftc,ttotfftc/ttotm*100,'%'
  !write(*,'(a30,2f15.4,a2)') 'cos transform2 ', ttotaux,ttotaux/ttotm*100,'%'
  write(*,'(a30,2f15.4,a2)') 'fft           ', ttotfft, ttotfft/ttotm*100,'%'
  write(*,*) '================================================================='

  !      open(31,file='tiempo.dat',status='unknown')
  !      write(31,*) numprocs,nthreads,blockl,nx,ny,nz,nsteps,real(ttotm),&
  !                  &real(ttotc),real(ttotr),real(ttotinty),real(ttotintx),&
  !                  &real(ttotim),real(ttotdy),real(ttotvy),real(ttotvdx),&
  !                  &real(ttotfft),real(ttotfftc)
end subroutine summary2

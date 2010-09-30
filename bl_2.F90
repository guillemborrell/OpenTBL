

!=========================================================================
!=========================================================================
!=======================     SECOND BL  ==================================


subroutine bl_2(mpiid,mpiid_global,comm_global,comm_local)
  use alloc_dns_2
  use main_val_2
  use names_2
  use temporal_2
  use point_2
  use statistics_2,only: ener
  use genmod_2,only: rthin,din,u0
  use ctesp_2
  use omp_lib
  use num_nodes
  implicit none
  include "mpif.h"

  integer,intent(in):: mpiid,mpiid_global,comm_global,comm_local
  real*8 dt,vardt
  integer isubstp,istep,ical,ierr
  logical:: vcontrol
  character*60:: newfile
  vcontrol=.false. !just checking correct time step
  
  ! Medimos tiempo ! 
  tiempo=0d0
  flag = 1
  iflagr =1
  zero=0
  times=0d0
  pi=4d0*atan(1d0)
  tvcrt=2d0*pi/(16*160)

  !$ CALL OMP_SET_DYNAMIC(.FALSE.)
  mpiid2=mpiid

  !Then number of processors used when compilling is assigned to th program:
  nummpi=numnodes_2; pnodes=numnodes_2

  if (mpiid.eq.0) ical=0   
  ! reads hrem.dat.
  ! creates coordinates for parallel pros's.
  if (mpiid.eq.0) write(*,*) '======== PROGRAM 2 BEGINS ==============='
  call  iniciap_2(mpiid,comm_local)
  call MPI_BCAST(mpi_inlet,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  if(mpiid.eq.0) write(*,*) 'BL2**********: GLOBAL MPIID mpiid1',mpiid_1,'mpi_inlet',mpi_inlet
  

  ! reads/generates inlet conditions initiates arrays (just in case)
  ! warning!! getstart reads y from file !!!

  if (mpiid .eq. 0) write(*,*) '======== CALLING GETSTART 2 ==============='
  if (mpiid2.eq.0) th2 = MPI_WTIME()                
  call getstartzy_2(u,v,w,p,dt,mpiid,comm_local)
  if (mpiid2 .eq. 0) then  
     th1 = MPI_WTIME()      
     tmp29 =tmp29+abs(th2-th1)
     write(*,*) 
     write(*,'(a40,f15.3,a3)') '===== FIELD 2 (u,v,w,p) read in: =======',tmp29,'Sg' 
     write(*,*)    
  endif
   
  call coef_2(mpiid)
  call inlet_retheta_2(u0,rthin,din,mpiid,comm_local) !compute Rth_inlet and d99_inlet
  call magicnumber_2(mpiid)

#ifdef CREATEPROFILES
  paso=-1 !Substeps counter
  if(mpiid.eq.0) open(170,file='extraprofiles',form='unformatted',status='unknown',convert='BIG_ENDIAN')
  open(270,file='plano_div.dat',form='unformatted',status='unknown',convert='BIG_ENDIAN')
  open(280,file='zero-mode.dat',form='unformatted',status='unknown',convert='BIG_ENDIAN')  
  call create_profiles_2(u,v,w,rthin,mpiid,comm_local)
  call impose_profiles_2(u,v,w,mpiid,comm_local) 
  close(170)  
#endif

  call alloabuff_2(ntotb,ntotv,ntot_corr,mpiid)

  rhsupa = 0d0
  rhswpa = 0d0
  rhsvpa = 0d0

#ifdef CREATEPROFILES
   if(mpiid.eq.0) write(*,*) 'Checking Divergence of the Field:'
   call check_divergence_2(u,v,w,rhsupa,mpiid)
#endif

  call mpi_barrier(comm_local,ierr)

#ifdef FLOPS
  CALL HPM_INIT()
  CALL HPM_START('WORK1')
#endif

#ifndef NOINFOSTEP
!Genflu info:
if(mpiid.eq.0) open(360,file=chinfoext,form='formatted',status='unknown',convert='BIG_ENDIAN')
#endif

#ifdef TRACE
write(newfile,'("Uflow.",i4.4)') mpiid
open(988,file=trim(chfile)//newfile,form='unformatted')
write(988) ib,ie,ny
#endif

  do istep = 1,nsteps
     if(mpiid.eq.0) then
        tc1 = MPI_WTIME()
!         write(*,'(a60,i6)') 'Second BL....................................................istep=',istep
     endif

     if (.TRUE.) setstep=.TRUE.
     if (mod(istep,avcrt)==0.and.istep>stest) dostat=.TRUE.

     !.............................................  
     do isubstp = 1,3
        if (mpiid2 .eq. 0) th1 = MPI_WTIME()               
        call mpi_barrier(comm_local,ierr)
!         IF(MPIID.EQ.0) WRITE(*,*) 'COMUNICADOR LOCAL============================',comm_local,'substep',isubstp
        call rhsp_2(u,v,w,p,rhsupa,rhsvpa,rhswpa,       &
             &    res,res,resv,resv,resw,resw,        & 
             &    rhsu,rhsv,rhsw,                     &
             &    wki1,wki1,wki2,wki2,wki3,wki3,      &
             &    wkp,wkp,wkpo,wkpo,bufuphy,buf_corr, & 
             &    dt,isubstp,ical,istep,mpiid,comm_local)

        call mpi_barrier(comm_local,ierr)    !********************
        if (mpiid2 .eq. 0) then  
           th2 = MPI_WTIME()      
           tmp13 =tmp13+abs(th2-th1)
        endif

        call outflow_correction_2(u,v,rhsupa,comm_local)
        
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp18 =tmp18+abs(th2-th1)
        endif
        call mpi_barrier(comm_local,ierr)   !************************

        vardt = 5d-1/dt/rkdv(isubstp)   

        call pois_2(u,v,w,p,res,res,resw,vardt,mpiid,comm_local)

        if (mpiid2 .eq. 0) then  
           th2 = MPI_WTIME()      
           tmp14 =tmp14+abs(th2-th1)
        endif
     enddo
     !.............................................
     call mpi_barrier(comm_local,ierr)
     tiempo = tiempo+dt
     if (times .ge. tvcrt) then
        times=times-tvcrt
     endif
     times=times+dt

#ifdef TRACE
     !* write timetrace
     call timetrace_2(u,v,w,p,tiempo,mpiid)
#endif

     ! I/O Operations --------------------------------
     if (mod(istep,stats).eq.0) then
        if (mpiid2.eq.0) th2 = MPI_WTIME()                                                     
        call escrst_2(ax,ay,az,cfl,tiempo,re,x,y,mpiid,ical,comm_local)
        ical=0
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp27 =tmp27+abs(th2-th1)
        endif
     endif

     if (mod(istep,reav).eq.0) then
        if (mpiid2.eq.0) th2 = MPI_WTIME() 
         call escribezy_2(u,v,w,p,dt,mpiid,comm_local)   !Vul & BG     
        if (mpiid2 .eq. 0) then  
           th1 = MPI_WTIME()      
           tmp28 =tmp28+abs(th2-th1)
        endif
     endif

     if (mpiid .eq. 0) then  ! Tiempos 
        tc2 = MPI_WTIME()      
        tmp1 = tc2-tc1     
        call summary1_2(istep,dt,vcontrol)             
     endif
#ifdef CHECKTIME
     call MPI_BCAST(vcontrol,1 ,MPI_LOGICAL,0,comm_local,ierr)
     IF(vcontrol) stop
#endif
  enddo


#ifndef NOINFOSTEP
close(360)
#endif
#ifdef CREATEPROFILES 
close(270)
close(280)
#endif
  if (mpiid.eq.0) call summary2_2()
endsubroutine bl_2




! 
! ---------------------------------------------------------------------!
! ----------------------- Subroutines ---------------------------------!
! ---------------------------------------------------------------------!
! 
! ----------------------------------------------------------------------*
!        lee/genera condiciones iniciales y de contorno
!        inicializa variables
! ----------------------------------------------------------------------*
! 


subroutine iniciap_2(mpiid,communicator)
  use alloc_dns_2
  use names_2
  use point_2
  use shared_mem_2
  use genmod_2
  use statistics_2
  use ctesp_2
  use omp_lib
  implicit none
  include "mpif.h"
  integer,intent(in):: communicator
  integer mpiid,k,blocks,elem,stride,iproc,fcr(1:nz1),i,ii,j
  integer posx,posy(1:npos),ierr,idat(20),i3

  real(8)       dat(6)
  character(99) text
  character(3) ext,ext2


  !----------------------------------------------------------------------!
  !      el  maestro lee parametros iniciales y los reparte
  !----------------------------------------------------------------------!

  if (mpiid .eq. 0) then
     write(*,*) 'ABRIENDO ARCHIVO HREM_2.DAT---------------------------------------'
     open(23,file='hrem_2.dat',status='old')   
     call avanza(text,23)
     !                  Re    ,ax    ,ay    ,az    ,r  
     read(text,*) dat(1),dat(2),dat(3),dat(4),dat(5) 
     call avanza(text,23)
     !                 cfl   ,istart
     read(text,*) dat(6),idat(5)   
     call avanza(text,23)
     write(*,*) text
     !                 nsteps,   reav,    stats   avcrt,  stest
     read(text,*) idat(1), idat(13), stats, idat(14),idat(15)
     call avanza(text,23)
     !           ifile1
     read(text,*) idat(7)
     write(*,*) 'idat',idat(7)
     call avanza(text,23)
     read(text,*) idat(8),idat(9),idat(10),idat(11),idat(12)
     ! Genflu flags 
     call avanza(text,23)
     read(text,*) flaggen,utauout,rthout,dout
     call avanza(text,23)
     read(text,'(a100)') chinit
     call avanza(text,23)
     read(text,'(a100)') chfile
     call avanza(text,23)
     read(text,'(a100)') chinfo
     close (23)
  endif

  call MPI_BCAST(dat    ,6 ,MPI_REAL8    ,0,communicator,ierr)
  call MPI_BCAST(idat   ,19,MPI_INTEGER  ,0,communicator,ierr)
  call MPI_BCAST(stats  ,1 ,MPI_INTEGER  ,0,communicator,ierr)
  call MPI_BCAST(flaggen,1 ,MPI_INTEGER  ,0,communicator,ierr)
  call MPI_BCAST(chfile,100,MPI_CHARACTER,0,communicator,ierr)  !for the Parallel IO Writting
  call MPI_BCAST(chinit,100,MPI_CHARACTER,0,communicator,ierr)  !for the Parallel IO Reading
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
     etfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.esp'
     vetfile=chfile(1:index(chfile,' ')-1)//'.'//ext//'.vesp'
     chinfoext=chinfo(1:index(chinfo,' ')-1)//'.'//ext2//'.dat'  !File for info at each step 
  endif

  tiempo = 0d0
  flag   = 1
  iflagr = 1
  zero   = 0
  times  = 0d0
  pi     = 4d0*atan(1d0)
  tvcrt  = 2d0*pi/(16d0*160d0)
  totalcal=0

!   ----------------------------------------------------------------------!
!             inicializa los punteros
!   ----------------------------------------------------------------------!
  allocate(ibeg   (0:nummpi-1),iend   (0:nummpi-1))
  allocate(pcibeg (0:nummpi-1),pciend (0:nummpi-1))
  allocate(pcibeg2(0:nummpi-1),pciend2(0:nummpi-1))

  call pointers_p2p_2(mpiid)   !! pointers et al. for changep2p & total sizes
  call pointersfft_2(mpiid) !pointers for thread_private indexes (used in FFTw)

  ib0 = max(2,ib)  ! -- i=1 is inflow, not integrated
  call alloa_2(mpiid)

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
  if(mpiid.eq.0) write(*,*) 'xout ', xout,'in Node #',mpiout
  if (mpiid==0) write(*,*) 'ALOCATANDO...................................'
  call alloavar_2(ntotb,ntotv,mpiid)
  call rfti_2(nz)
  call cosfti_2(nx-1)

    
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
     open(108,file='tablespectra',form='unformatted',status='unknown',convert='BIG_ENDIAN')
     do j=1,nspec
        read(108) jspecy(j,1:lxp)
     enddo
     write(*,*) 'valores de jspecy'
     do j=1,nspec
      write(*,'(10i6)') jspecy(j,1:lxp)
     enddo
     close(108)
  endif
  call MPI_BCAST(jspecy,size(jspecy),MPI_INTEGER,0,communicator,ierr)
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
endsubroutine iniciap_2

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------
subroutine inlet_retheta_2(u0,rthin,din,mpiid,communicator)
use ctesp_2,only:ny
use alloc_dns_2,only:re,dy,y
implicit none
include "mpif.h"
integer:: communicator
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
CALL MPI_BCAST(rthin,1,MPI_REAL8,0,communicator,ierr)
CALL MPI_BCAST(din,1,MPI_REAL8,0,communicator,ierr)

endsubroutine inlet_retheta_2

subroutine timetrace_2(u,v,w,p,time,mpiid)
  use point
  use names
  use genmod
  use alloc_dns
  use ctesp
  implicit none
  include 'mpif.h'
  real*8,dimension(nz1,ny+1,ib:ie)::u,w    
  real*8,dimension(nz1,ny  ,ib:ie)::p,v
  integer i,j,k,l,irec,mpiid
  real*8 time

  write(988) time,((u(1,k,i),k=1,ny+1),i=ib,ie), &
&                 ((v(1,k,i),k=1,ny  ),i=ib,ie), &
&                 ((w(1,k,i),k=1,ny+1),i=ib,ie)
! if(mpiid.eq.0)write(*,*)time,(u(1,2,i),i=ib,ie)
  
  flush(988)

end subroutine timetrace_2

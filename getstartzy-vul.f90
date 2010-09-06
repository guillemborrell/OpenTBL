!=====================================================================
!       Genera o lee las condiciones iniciales (Campos u,v,w,p)
!       Inicializa las arrays para estadisticas
!	Lee la malla en la direccion Y
!
! Lee todos los planos YX donde estan los campos y la presion en R*4 y
! los copia en R*8 en los buffers (u,v,w,p). 
!
!  El master se copia su campo correspondiente y luego envia como R*4 los trozos a
!  los otros cores que se lo copian como R*8
!
!       ACHTUNG!!    allocates one input plane, check there is space
!==========================================================================
subroutine getstartzy(u,v,w,p,dt,mpiid)
  use alloc_dns
  use statistics
  use names
  use point
  use genmod
  use ctesp
  implicit none
  include "mpif.h"

  ! ---------------------------- I/O -----------------------------------!
  real(8),dimension(nz1,ny,ib:ie)  :: v,p
  real(8),dimension(nz1,ny+1,ib:ie):: u,w
  real*4, dimension(:,:,:),allocatable:: resu 
  real*4, dimension(:),allocatable:: ulast 

  ! -------------------------- Work ----------------------------------------!
  integer status(MPI_STATUS_SIZE),ierr,mpiid
  integer nxr,nyr,nzr,nz1r,nzz,j,i,k,l,dot,lim2,rsize,irec,rsize1,rsize2,ji
  real(8) jk,dt,dum(20)
  character text*99, uchar*1,fil1*80,fil2*80,fil3*80,fil4*80

  ! --------------------------  Programa  ----------------------------------!
  !       lee el fichero 
  if (mpiid.eq.0) then
     write(*,*) 'Leyendo del fichero'
     write(*,*) chinit(1:index(chinit,' ')-1)//'.'//'u'
     
! 
!      !!!!!!!!!!!!!!! V0
!      !Zero mode for V (from torroja statistics Re=1550)
!      open (90,file='v0-vmagic-1550',form='unformatted',status='old',convert='BIG_ENDIAN')
!      read( 90) (v0(i), i=1,ny+1)    
!      close(90)     
! 
!      !!!!!!!!!!!!!!! U0 
!      open (40,file='uminlet-yBG',form='unformatted',status='old',convert='BIG_ENDIAN')
!      read(40) (u0(i),i=1,ny+1)
!      close(40)

     !rsize = 15  !in BG change this size.... rsize*4
     
     rsize=500
     fil1=chinit(1:index(chinit,' ')-1)//'.'//'u'
     open (20,file=fil1,status='old',form='unformatted',access='direct',recl=rsize,convert='BIG_ENDIAN')
     write(*,*) 'VULCANO: file open, reading tiempo and dimensions'
     read(20,rec=1) uchar,tiempo,jk,jk,jk,jk,jk,nxr,nyr,nzr
     write(*,*) '==========================================='
     write(*,*) 'in file    ', uchar,tiempo,nxr,nyr,nzr
     write(*,*) 'in ctes    ',nx,ny,nz2
     close(20)          
  endif

  call MPI_BCAST(nxr,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nzr,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nyr,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
  if (ny.ne.nyr) then
     if (mpiid==0) write(*,*) 'changing the y grid has to be done separately'
     if (mpiid==0) write(*,*) 'ny=',ny,'nyr',nyr
     stop
  endif

  u = 0d0
  v = 0d0
  w = 0d0
  p = 0d0

  nz1r=2*(nzr+1) 
  rsize = nz1r*(ny+1)*4
  rsize1 = nz1r*(ny+1)
  rsize2 = nz1r*(ny  )
  nzz = min(nz1r,nz1)
  allocate (resu(nz1r,ny+1,4),ulast(ny+1))

  if (mpiid.eq.0) then
     fil1=chinit(1:index(chinit,' ')-1)//'.'//'u'
     open (10,file=fil1,status='unknown', &
          & form='unformatted',access='direct',recl=rsize1,convert='BIG_ENDIAN')

     !!  store end mean U profile, to extend the layer, if needed
     irec=nxr+1
     read(10,rec=irec) ulast 
     read(10,rec=1) uchar,tiempo,jk,jk,jk,jk,jk,nxr,nyr,nzr,ji,timeinit,dt, &
          & dum, (x(i), i=0,nxr+1), (y(i), i=0,nyr+1), (um(i), i=1,nyr+1)
     write(*,*) 'in file    ', uchar,tiempo,nxr,nyr,nzr
     write(*,*) '              x                   y                  um'
     write(*,*) '--------------------------------------------------------------------'      
     do i=1,10
      write(*,'(3f20.6)') x(i),y(i),um(i)
     enddo
     !opening rest of the files: 
     fil2=chinit(1:index(chinit,' ')-1)//'.'//'v'
     open (11,file=fil2,status='unknown', &
          & form='unformatted',access='direct',recl=rsize2,convert='BIG_ENDIAN')

     fil3=chinit(1:index(chinit,' ')-1)//'.'//'w'
     open (12,file=fil3,status='unknown', &
          & form='unformatted',access='direct',recl=rsize1,convert='BIG_ENDIAN')

     fil4=chinit(1:index(chinit,' ')-1)//'.'//'p'
     open (13,file=fil4,status='unknown', &
          & form='unformatted',access='direct',recl=rsize2,convert='BIG_ENDIAN')

     write(*,*) '----- files OPEN ------'
     write(*,*) fil1
     write(*,*) fil2
     write(*,*) fil3
     write(*,*) fil4
  endif

  call MPI_BCAST(ulast,ny+1,mpi_real4,0,MPI_COMM_WORLD,ierr)
  
  if (mpiid.eq.0) then
     irec=1
     !!  start reading the flow field  
     do i=ib,ie       
        if (i.le.nxr) then        
           irec = irec+1
           read(10,rec=irec) resu(1:nzz,1:ny+1,1)
           read(11,rec=irec) resu(1:nzz,1:ny,2)
           read(12,rec=irec) resu(1:nzz,1:ny+1,3)
           read(13,rec=irec) resu(1:nzz,1:ny,4)

           u(1:nzz,1:ny+1,i) = resu(1:nzz,1:ny+1,1)
           v(1:nzz,1:ny,i)   = resu(1:nzz,1:ny,2)
           w(1:nzz,1:ny+1,i) = resu(1:nzz,1:ny+1,3)
           p(1:nzz,1:ny,i)   = resu(1:nzz,1:ny,4)

           if (i==1) then
            u0=resu(1,:,1)
            v0=resu(1,:,2)
           endif
        else
           u(1,:,i) = ulast    !!! extend downstream with mean profile
        endif
     enddo
   
     do dot = 1,nummpi-1   ! -- read for the other nodes 
        do i= ibeg(dot),iend(dot)
        if (mod(i,75).eq.0) write(*,*) 'Read & Send up to:',i
           if (i.le.nxr) then
              irec = irec+1
              read(10,rec=irec) resu(1:nzz,1:ny+1,1)
              read(11,rec=irec) resu(1:nzz,1:ny,2)
              read(12,rec=irec) resu(1:nzz,1:ny+1,3)
              read(13,rec=irec) resu(1:nzz,1:ny,4)
              call MPI_SEND(resu,rsize,MPI_real4,dot,1,MPI_COMM_WORLD,ierr)
           endif
        enddo
     enddo
      close(10);close(11);close(12);close(13)
     call flush(6)
  else  !   --- the other nodes receive the information  
     do i=ib,ie
        if (i.le.nxr) then
           call MPI_RECV(resu,rsize,MPI_real4,0,1,MPI_COMM_WORLD,status,ierr)
           u(1:nzz,1:ny+1,i) = resu(1:nzz,1:ny+1,1)
           v(1:nzz,1:ny,i)   = resu(1:nzz,1:ny,2)
           w(1:nzz,1:ny+1,i) = resu(1:nzz,1:ny+1,3)
           p(1:nzz,1:ny,i)   = resu(1:nzz,1:ny,4)
        else
           u(1,:,i) = ulast
        endif
     enddo
  endif

  call MPI_BCAST(tiempo,1,mpi_real8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(y,ny+2,mpi_real8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,1,mpi_real8,0,MPI_COMM_WORLD,ierr)

  deallocate(resu,ulast)
end subroutine getstartzy

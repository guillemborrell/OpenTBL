!=====================================================================
!       Genera o lee las condiciones iniciales (Campos u,v,w,p)
!        Inicializa las arrays para estadisticas
!	Lee la malla en la direccion Y
!
! Lee todos los planos YX donde estan los campos y la presion en R*4 y
! los copia en R*8 en los buffers (u,v,w,p). 
!
!   El master se copia su campo correspondiente y luego envia como R*4 los trozos a
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
  integer nxr,nyr,nzr,nz1r,nzz,j,i,k,l,dot,lim2,rsize,irec,nyr2,nxr2,nzr2
  real(8) jk,dt,dum(20)
  character*99 text

  ! --------------------------  Programa  ----------------------------------!

  !       lee el fichero 

  if (mpiid.eq.0) then
     write(*,*) 'Leyendo del fichero'
     write(*,*) chinit

     write(*,*) 'file open, reading tiempo and dimensions'
     
     irec = 1
     
     write(*,*) 'RECL, 8 BYTES #OF ELEMENTS: =================================='
     rsize = 15*4
     open (100,file=chinit,form='unformatted',status='old',access='direct',recl=rsize)
     
     read(100,rec=irec) tiempo,jk,jk,jk,jk,jk,nxr,nyr,nzr
     close(100)
     write(*,*) '--first file read------'

     nyr2=nyr;nxr2=nxr;nzr2=nzr;nxr=nx;nyr=ny;nzr=nz2;

     write(*,*) 'in file---------------------',tiempo,nxr2,nyr2,nzr2
     write(*,*) 'in ctes',nx,ny,nz2
     
  endif

  call MPI_BCAST(nx,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nzr,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nyr,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
!  if (ny.ne.nyr) then
!     if (mpiid==0) write(*,*) 'changing the y grid has to be done separately'
!     if (mpiid==0) write(*,*) 'ny=',ny,'nyr',nyr
!     stop
!  endif



  u = 0d0
  v = 0d0
  w = 0d0
  p = 0d0

  nz1r=2*(nzr+1) 
  rsize = (nz1r*(ny+1)*4)*4
  nzz = min(nz1r,nz1)
  allocate (resu(nz1r,ny+1,4),ulast(ny+1))
  if (mpiid.eq.0) then
      write(*,*) '--second file open for read------'
     open (100,file=chinit,form='unformatted',status='old',access='direct',recl=rsize)
      
     !!  store end mean U profile, to extend the layer, if needed
     irec=nxr+1
    
     read(100,rec=1) tiempo,jk,jk,jk,jk,jk,nxr2,nyr2,nzr2,timeinit,dt, &
          & dum, (dum(1), i=0,nxr2+1), (y(i), i=0,nyr2+1), (um(i), i=1,nyr2+1)
     write(*,*) '--second file read------'
  endif

um(nyr2+2:ny+1)=1d0
do i=nyr2+2,ny+1
 y(i)=y(nyr2+1)+0.05d0*i
enddo

ulast=um

 call MPI_BCAST(ulast,ny+1,mpi_real4,0,MPI_COMM_WORLD,ierr)
  
 
     do i=ib,ie

        if (i.le.nxr) then        
         do j=1,ny
          do k=1,nzz
           u(k,j,i) =(k+j+i)/(2*k*i+j) 
           v(k,j,i) =(k+5*j+i)/(2*k+i+j) 
           w(k,j,i) =(k-j+i)/(2*k+i+j) 
           p(k,j,i)   =(k+j*3+i)/(2*k+i+j)
          enddo
         enddo
         do k=1,nzz
           u(k,ny+1,i) =(k+ny+1+i)/(2*k+i+j) 
           w(k,ny+1,i) = (k-ny-1+i)/(2*k+i+j)
         enddo 
           if (i==1) u0=um
        else
           u(1,:,i) = um    !!! extend downstream with mean profile
        endif

     enddo

    

  call MPI_BCAST(tiempo,1,mpi_real8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(y,ny+2,mpi_real8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,1,mpi_real8,0,MPI_COMM_WORLD,ierr)

  deallocate(resu,ulast)
end subroutine getstartzy

! Subroutines needed for the global transpose.

! Dependencies

!  * shared_mem: defines the derived types nodedata and domaindata
!  * point: needed by chp2x and chx2p routines that perform the transposes
!           in the bl code.

! The transpose algorithm is straightforward looking at, for instance, chp2xu.

! 0. comm_setup routine is called to build the derived types with all the 
!    necessary information to transpose the data. This routine is completely
!    atomic and does not rely on bl data at all.  It can be reused for any
!    box-like geometry that needs to be transposed.

! 1. Groups of pencils are reorganized to be contiguous and form a single
!    message.  At the same time the numbers are converted to single precision.
!    This is the transpose1 when changing from planes to pencils and
!    transpose4 when changing from pencils to planes.

! 2. MPI_ALLTOALLV is called.  The count and disp vectors are contained
!    in the derived types and, in most of times, declared and allocated in
!    a module.

! 3. The second transpose is a plain transpose to change the contiguous index.
!    This is transpose2 when changing from planes to pencils and transpose3
!    when changing from pencils to planes.

subroutine comm_setup_2(NX,NY,NZ,rank,size,node,domain,pblock)
  ! Previous configuration for the parallel transpose package
  ! IMPORTANT: Relies on the shared_mem module where nodedata and
  !            domaindata derived types are declared.
  ! Input arguments:
  !  * NX, NY, NZ: Size of the actual data box in double precision.
  !  * rank: MPI rank
  !  * size: MPI size
  !  * node: local nodedata derived type
  !  * domain: local domaindata derived type.
  use shared_mem_2
  implicit none

  integer, intent(in):: NX,NY,NZ,rank,size,pblock
  type(nodedata), intent(out):: node
  type(domaindata), intent(out):: domain

  integer:: penpn, plpn, penpn2
  integer:: i

  !Allocate the structs with size
  allocate(node%scount(size))
  allocate(node%sdisp(size))
  allocate(node%rcount(size))
  allocate(node%rdisp(size))

  allocate(domain%ib(size))
  allocate(domain%ie(size))
  allocate(domain%pb(size))
  allocate(domain%pe(size))

  !! check that pencils are multiples of pblock
  if (mod(NZ*NY,pblock).ne.0) then
     write(*,*) 'non-integer pencil blocks:', NZ*NY, pblock, 'stopping'
     stop
  endif

  !! Pencils
  penpn2 = NZ*NY/size/pblock
  domain%widenodes = NZ*NY/pblock-penpn2*size

  penpn  = penpn2*pblock

  domain%pb(1) = 1
  if(domain%widenodes > 0) then
     do i = 1,domain%widenodes
        domain%pe(i) = domain%pb(i) + penpn+pblock-1
        domain%pb(i+1) = domain%pe(i) + 1
     end do
  end if
  do i = domain%widenodes+1, size-1
     domain%pe(i) = domain%pb(i) + penpn - 1
     domain%pb(i+1) = domain%pe(i) + 1
  end do
  domain%pe(size) = domain%pb(size) + penpn - 1

  !! Planes

  plpn = NX/size
  domain%fatnodes = NX-plpn*size

  domain%ib(1) = 1
  if (domain%fatnodes > 0) then
     do i = 1,domain%fatnodes
        domain%ie(i) = domain%ib(i) + plpn
        domain%ib(i+1) = domain%ie(i) + 1
     end do
  end if
  do i=domain%fatnodes+1,size-1
     domain%ie(i) = domain%ib(i) + plpn - 1
     domain%ib(i+1) = domain%ie(i) + 1
  end do

  domain%ie(size) = domain%ib(size) + plpn - 1

  node%startpl = domain%ib(rank+1)
  node%endpl = domain%ie(rank+1)
  node%planes = node%endpl - node%startpl + 1

  node%startpen = domain%pb(rank+1)
  node%endpen = domain%pe(rank+1)
  node%pencils = node%endpen - node%startpen + 1 

  domain%NX = NX
  domain%NY = NY
  domain%NZ = NZ

  domain%pencils = NZ*NY
  domain%planes = NX
  domain%planesize = NZ*NY
  domain%pencilsize = NX

  node%size = max(node%pencils*domain%pencilsize,&
       & node%planes*domain%planesize)

  !! Counts and displacements for collective communications

  node%scount(1) = (domain%ie(rank+1)-domain%ib(rank+1)+1)*&
       &(domain%pe(1)-domain%pb(1)+1)
  node%sdisp(1) = 0

  do i = 2, size
     node%scount(i) = (domain%ie(rank+1)-domain%ib(rank+1)+1)*&
          &(domain%pe(i)-domain%pb(i)+1)
     node%sdisp(i) = node%sdisp(i-1)+node%scount(i-1)
  end do

  node%rcount(1) = (domain%ie(1)-domain%ib(1)+1)*&
       &(domain%pe(rank+1)-domain%pb(rank+1)+1)
  node%rdisp(1) = node%scount(size)+node%sdisp(size)

  do i = 2, size
     node%rcount(i) = (domain%ie(i)-domain%ib(i)+1)*&
          &(domain%pe(rank+1)-domain%pb(rank+1)+1)
     node%rdisp(i) = node%rdisp(i-1)+node%rcount(i-1)
  end do

  node%bound = node%rdisp(1)

end subroutine comm_setup_2

subroutine pointers_p2p_2(rank)
  ! Driver routine and wrapper to the comm_setup to adapt the global
  ! transpose to the bl code.  This routine is mostly junk to keep the changes
  ! to the bl code minimal.

  ! Note that the variable size of the comm_setup routine is read from
  ! the ctesp module increasing the number of dependencies needed.

  ! Tne node and domain derived types are stored in the point module.
  ! I need a two couples, one for u-like variables and other for v-like
  ! variables.
  ! 
  ! ADDED: Variables type correlation
  use point_2
  use ctesp_2
  use shared_mem_2
  implicit none
  integer:: rank,i
  write(*,*) '						NUMMPI',nummpi,'rank',rank
  !========Setting Structure Data Types:
  call comm_setup_2(nx,ny,nz1,rank,nummpi,nodev,domainv,1)  
  call comm_setup_2(nx,ny+1,nz1,rank,nummpi,nodeu,domainu,1)
  !For the correlations: 
  call comm_setup_2(nx,ncorr,nz1  ,rank,nummpi,node_corr ,domain_corr,2) !Complex*16 Plane = nz1 R8 elements (each node a even number of pencils, "2")
  call comm_setup_2(nx,ncorr,nz1/2,rank,nummpi,node_corr2,domain_corr2,1) !Real*8 plane = nz1/2 R8 elements
  !==========================================
  

  !Global List of indexes  
  ibeg = domainu%ib
  iend = domainu%ie
  pcibeg = domain_corr%pb !List of the beginning and end of the pencils
  pciend = domain_corr%pe
  pcibeg2 = domain_corr2%pb
  pciend2 = domain_corr2%pe
  !Size to allocate arrays
  ntotb = nodeu%size
  ntotv = nodev%size
  ntot_corr =node_corr%size
  ntot_corr2=node_corr2%size
  !Local indexes for each node
  ib = nodeu%startpl
  ie = nodeu%endpl   
  pcib=node_corr%startpen !Pencils start and end -for correlations-
  pcie=node_corr%endpen
  pcib2=node_corr2%startpen !Pencils start and end -for correlations-
  pcie2=node_corr2%endpen
  !Total number of pencils per node:
  mpu = nodeu%pencils
  mpv = nodev%pencils  
  mp_corr =node_corr%pencils
  mp_corr2=node_corr2%pencils !Half number of pencils (k x k*)
  !Number of planes per node
  mmx = ie-ib+1

	  
  if(rank.eq.0) then
    write(*,'(a40,i10,a4,i5,a4,i5)') 'Pencils Buffer Correlations:mp_corr',mp_corr,'pb',pcib,'pe',pcie
    write(*,'(a40,i10,a4,i5,a4,i5)') 'Pencils Correlations: mp_corr2',mp_corr2,'pb',pcib2,'pe',pcie2
    write(*,'(a40,i10)') 'Pencils u',mpu
    write(*,'(a40,i10)') 'Pencils v',mpv
    write(*,*) 'NX',NX,'NCORR',NCORR,'NZ1',NZ1,'NZ2+1',NZ2+1
    write(*,'(a40,i15)') 'ntot_corr',ntot_corr		
    write(*,'(a40,i15)') 'ntot_corr2',ntot_corr2		
    write(*,*) '-----------------------------------------------------------'
!     write(*,*) 'PENCIL INDEXES:'
!     do i=0,nummpi-1
!     write(*,'(a30,2i10,a15,i10,a15,i10,a8,i6)') 'Pencils Buffer Correlations',pcibeg(i),pciend(i),'1st P.Corr=',2*pcibeg2(i)-1, 'Last=',2*pciend2(i),'NODE',i  
!     enddo
  endif
  
end subroutine pointers_p2p_2

subroutine chp2x_2(pen,plan,buf,rank,jee,communicator)
  ! Driver routine to keep the changes to bl minimal. Look at
  ! the chp2xu and chp2xv routines to understand the whole thing.

  ! The only thing that can be misunderstood is that this routine
  ! gets the MPI rank as an argument while chp2xv wants the MPI size.
  ! The reason is that in the new transpose algorithm the MPI rank
  ! is implicitly stored in the node and domain derived types while
  ! the MPI size is still needed as an argument.

  ! OTOH, the old change reads the size from a module and passes the
  ! MPI rank to prevent anyone to overwrite the variable, something that
  ! may cause serious bugs difficult to catch.
  ! 
  ! "control" variable intented for which change for correlations must be used
  use shared_mem_2
  use point_2
  use ctesp_2
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodeu%size), intent(in):: plan
  real(kind = 8), dimension(nodeu%size), intent(out):: pen
  real(kind = 8), dimension(nodeu%size), intent(inout):: buf
  integer, intent(in):: rank,jee

  if (jee == domainv%NY) then
     call chp2xv_2(pen,plan,buf,nummpi,communicator) !here does not matter the value of control
  elseif (jee == domainu%NY) then
     call chp2xu_2(pen,plan,buf,nummpi,communicator)
  elseif (jee == domain_corr%NY) then  !Complex*16 Plane = nz1 R8 elements
     call chp2xc_2(pen,plan,buf,nummpi,communicator)    
  end if
end subroutine chp2x_2

subroutine chx2p_2(pen,plan,buf,rank,jee,communicator)
  use shared_mem_2
  use point_2
  use ctesp_2
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodeu%size), intent(in):: plan
  real(kind = 8), dimension(nodeu%size), intent(out):: pen
  real(kind = 8), dimension(nodeu%size), intent(inout):: buf
  integer, intent(in):: rank,jee

  if (jee == domainv%NY) then
     call chx2pv_2(pen,plan,buf,nummpi,communicator)
  elseif (jee == domainu%NY) then
     call chx2pu_2(pen,plan,buf,nummpi,communicator)
  elseif (jee == domain_corr%NY) then  !Complex*16 Plane = nz1 R8 elements
     call chx2pc_2(pen,plan,buf,nummpi,communicator) 
  end if

end subroutine chx2p_2



!!========= change subroutines for V kind structures:
subroutine chp2xv_2(pen,plan,buf,size,communicator)
  ! It is a good idea to write this routine to understand how this transpose
  ! works.

  ! First, the transpose1 is used to change locally the alignment and the
  ! precision of the array.  Then the global communicator is called using
  ! the appropiate node and domain structures.  The last step is transposing
  ! the received array to align the data properly.

  ! Timers are not critical and can be removed.
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodev%size), intent(in):: plan
  real(kind = 8), dimension(nodev%size), intent(out):: pen
  real(kind = 8), dimension(nodev%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr

  if (mpiid2.eq.0) tm2 = MPI_WTIME()

  call transpose1_2(plan,buf,nodev,domainv,size)

  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

  
  call MPI_ALLTOALLV(buf, nodev%scount, nodev%sdisp, MPI_REAL4,&
       &buf, nodev%rcount, nodev%rdisp, MPI_REAL4,&
       &communicator, ierr)

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif

  call transpose2_2(buf(nodev%bound/2+1:nodev%size),pen,nodev,domainv)

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

end subroutine chp2xv_2

subroutine chx2pv_2(plan,pen,buf,size,communicator)
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodev%size), intent(in):: pen
  real(kind = 8), dimension(nodev%size), intent(out):: plan
  real(kind = 8), dimension(nodev%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr


  if (mpiid2.eq.0) tm2 = MPI_WTIME()

  call transpose3_2(pen,buf(nodev%bound/2+1:nodev%size),nodev,domainv)

  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

  call MPI_ALLTOALLV(buf, nodev%rcount, nodev%rdisp, MPI_REAL4,&
       &buf, nodev%scount, nodev%sdisp, MPI_REAL4,&
       &communicator, ierr)

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif

  call transpose4_2(buf,plan,nodev,domainv,size)

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif
end subroutine chx2pv_2

!!========= change subroutines for U kind structures:
subroutine chp2xu_2(pen,plan,buf,size,communicator)
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodeu%size), intent(in):: plan
  real(kind = 8), dimension(nodeu%size), intent(out):: pen
  real(kind = 8), dimension(nodeu%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr
  if (mpiid2.eq.0) tm2 = MPI_WTIME()

 

  call transpose1_2(plan,buf,nodeu,domainu,size)
  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

  call MPI_ALLTOALLV(buf, nodeu%scount, nodeu%sdisp, MPI_REAL4,&
       &buf, nodeu%rcount, nodeu%rdisp, MPI_REAL4,&
       &communicator,ierr)


  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif

  call transpose2_2(buf(nodeu%bound/2+1:nodeu%size),pen,nodeu,domainu)
  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

end subroutine chp2xu_2

subroutine chx2pu_2(plan,pen,buf,size,communicator)
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(nodeu%size), intent(in):: pen
  real(kind = 8), dimension(nodeu%size), intent(out):: plan
  real(kind = 8), dimension(nodeu%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr
  if (mpiid2.eq.0) tm2 = MPI_WTIME()
  call transpose3_2(pen,buf(nodeu%bound/2+1:nodeu%size),nodeu,domainu)
  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif
  call MPI_ALLTOALLV(buf, nodeu%rcount, nodeu%rdisp, MPI_REAL4,&
       &buf, nodeu%scount, nodeu%sdisp, MPI_REAL4,&
       &communicator, ierr)
  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif
  call transpose4_2(buf,plan,nodeu,domainu,size)
  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

end subroutine chx2pu_2

!!========= change subroutines for Correlation kind structures. C*16 Planes:
subroutine chp2xc_2(pen,plan,buf,size,communicator)
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(node_corr%size), intent(in):: plan
  real(kind = 8), dimension(node_corr%size), intent(out):: pen
  real(kind = 8), dimension(node_corr%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr

  if (mpiid2.eq.0) tm2 = MPI_WTIME()


  call transpose1_2(plan,buf,node_corr,domain_corr,size)

  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

  call MPI_ALLTOALLV(buf, node_corr%scount, node_corr%sdisp, MPI_REAL4,&
       &buf, node_corr%rcount, node_corr%rdisp, MPI_REAL4,&
       &communicator, ierr)

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif


  call transpose2_2(buf(node_corr%bound/2+1:node_corr%size),pen,node_corr,domain_corr)

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

end subroutine chp2xc_2



subroutine chx2pc_2(plan,pen,buf,size,communicator)
  use shared_mem_2
  use point_2
  use temporal_2
  use ctesp_2,only: mpiid2
  include "mpif.h"
  integer,intent(in)::communicator
  real(kind = 8), dimension(node_corr%size), intent(in):: pen
  real(kind = 8), dimension(node_corr%size), intent(out):: plan
  real(kind = 8), dimension(node_corr%size), intent(inout):: buf
  integer, intent(in):: size
  integer:: ierr

  if (mpiid2.eq.0) tm2 = MPI_WTIME()

  call transpose3_2(pen,buf(node_corr%bound/2+1:node_corr%size),node_corr,domain_corr)

  if (mpiid2.eq.0) then
     tm1  = MPI_WTIME()
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

  call MPI_ALLTOALLV(buf, node_corr%rcount, node_corr%rdisp, MPI_REAL4,&
       &buf, node_corr%scount, node_corr%sdisp, MPI_REAL4,&
       &communicator, ierr)

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp2= tmp2 + abs((tm2-tm1))
  endif

  call transpose4_2(buf,plan,node_corr,domain_corr,size)

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()     
     tmp3 = tmp3 + abs((tm2-tm1))
  endif

end subroutine chx2pc_2




!================================================================
!================================================================
!======================       TRANSPOSES   ======================

subroutine transpose1_2(var,buf,node,domain,size)
  ! First transpose to change from stream normal planes to x aligned
  ! pencils.

  ! The module shared_mem is needed because the node and domain
  ! arguments are declared there.

  ! Previously the comm_setup must be called to allocate and build the
  ! node and domain derived types.

  ! Input arguments:
  ! * var, buf: Input and output buffers. DOES NOT WORK INPLACE!!!
  ! * node, domain: Node and domain derived types with the collective 
  !                 comm data
  ! * size: MPI size.

  ! This routine is doubly vectorised.  OMP is used to share the
  ! transpose workload amongst all the available cores and the data is
  ! moved in blocks of size the number of pencils.  There is no other
  ! way to implement blocking with higher performance.

  ! transpose4 is the reverse routine.

  use shared_mem_2
  implicit none
  type(domaindata), intent(in):: domain
  type(nodedata), intent(in):: node
  real(kind = 8), dimension(domain%planesize,node%planes), intent(in):: var
  integer, intent(in):: size
  real(kind = 4), dimension(domain%planesize*node%planes), intent(out):: buf

  integer:: rank,plan
  integer:: pencilsi,pbeg,pend
  integer:: zbeg,zend

  !$OMP PARALLEL DO PRIVATE(plan,pencilsi,pbeg,pend,zbeg,zend)
  do rank = 0,size-1
     pencilsi = domain%pe(rank+1) - domain%pb(rank+1) + 1
     do plan = 0,node%planes-1
        pbeg = ((domain%pe(rank+1)-pencilsi)*node%planes) + (pencilsi*plan) + 1
        pend = ((domain%pe(rank+1)-pencilsi)*node%planes) + (pencilsi*(plan + 1))
        zbeg = domain%pb(rank+1)
        zend = domain%pe(rank+1)
        buf(pbeg:pend) = real(var(zbeg:zend,plan+1),kind = 4)
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine transpose1_2


subroutine transpose4_2(buf,var,node,domain,size)
  use shared_mem_2
  implicit none
  type(domaindata), intent(in):: domain
  type(nodedata), intent(in):: node
  real(kind = 4), dimension(domain%planesize*node%planes), intent(in):: buf
  integer, intent(in):: size
  real(kind = 8), dimension(domain%planesize,node%planes), intent(out):: var

  integer:: rank,plan
  integer:: pencilsi,pbeg,pend
  integer:: zbeg,zend

  !$OMP PARALLEL DO PRIVATE(plan,pencilsi,pbeg,pend,zbeg,zend)
  do rank = 0,size-1
     do plan = 0,node%planes-1
        pencilsi = domain%pe(rank+1) - domain%pb(rank+1) + 1
        pbeg = ((domain%pe(rank+1)-pencilsi)*node%planes) + (pencilsi*plan) + 1
        pend = ((domain%pe(rank+1)-pencilsi)*node%planes) + (pencilsi*(plan + 1))
        zbeg = domain%pb(rank+1)
        zend = domain%pe(rank+1)
        var(zbeg:zend,plan+1) = real(buf(pbeg:pend),kind = 8)
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine transpose4_2

subroutine transpose2_2(buf,var,node,domain)
  ! TODO: Find a way to transpose inplace. It would greatly simplify
  ! the future versions of boundary layer code.

  ! Relies on the compiler's ability to set an OMP workshare
  ! environment so it is definetely NOT portable.  It seems to work
  ! fine with the IBM XL compiler OMP implementation.

  use shared_mem_2
  implicit none

  type(domaindata), intent(in):: domain
  type(nodedata), intent(in):: node
  real(kind = 4), dimension(node%pencils,domain%pencilsize), intent(in):: buf
  real(kind = 8), dimension(domain%pencilsize,node%pencils), intent(out):: var

  integer:: i,bs

  if (domain%pencilsize > 4*node%pencils) then
     ! blocking with square blocks
     bs = domain%pencilsize/node%pencils
     !$OMP PARALLEL DO
     do i = 0,bs-2
        var(i*node%pencils+1:(i+1)*node%pencils,1:node%pencils) = real(&
             & transpose(&
             & buf(1:node%pencils,i*node%pencils+1:(i+1)*node%pencils)),&
             & kind = 8)
     end do
     !$OMP END PARALLEL DO
     ! Last block is serial
     var((bs-1)*node%pencils+1:domain%pencilsize,1:node%pencils) = real(&
          & transpose(&
          & buf(1:node%pencils,(bs-1)*node%pencils+1:domain%pencilsize)),&
          & kind = 8)
  else
     ! no blocking
     !$OMP PARALLEL WORKSHARE
     var = real(transpose(buf),kind = 8)
     !$OMP END PARALLEL WORKSHARE
  end if

end subroutine transpose2_2

subroutine transpose3_2(var,buf,node,domain)
  use shared_mem_2
  implicit none
  type(domaindata), intent(in):: domain
  type(nodedata), intent(in):: node
  real(kind = 8), dimension(domain%pencilsize,node%pencils), intent(in):: var
  real(kind = 4), dimension(node%pencils,domain%pencilsize), intent(out):: buf

  integer:: i,bs

  if (domain%pencilsize > 4*node%pencils) then
     ! blocking with square blocks
     bs = domain%pencilsize/node%pencils
     !$OMP PARALLEL DO
     do i = 0,bs-2
        buf(1:node%pencils,i*node%pencils+1:(i+1)*node%pencils) = real(&
             & transpose(&
             & var(i*node%pencils+1:(i+1)*node%pencils,1:node%pencils)),&
             & kind = 4)
     end do
     !$OMP END PARALLEL DO
     ! Last block is serial
     buf(1:node%pencils,(bs-1)*node%pencils+1:domain%pencilsize) = real(&
          & transpose(&
          & var((bs-1)*node%pencils+1:domain%pencilsize,1:node%pencils)),&
          & kind = 4)
  else
     !$OMP PARALLEL WORKSHARE
     buf = real(transpose(var),kind = 4)
     !$OMP END PARALLEL WORKSHARE
  end if

end subroutine transpose3_2

!----------------------------------------------------------------------!
!   FAST TRANSFORMS INTERFACE PACKAGE (discrete & cosinus)             !
!                                                                      !
! Use the particular routines below only when you transform in         !
! Fourier in one direction (the 'z' direction).                        !
! Adapted from general fou3D by M.Simens 22-12-03                      !
!                                                                      !
!----------------------------------------------------------------------!

!/********************************************************************/
  !/*                                                                  */
  !/*    makes a 2-dimensional real to complex fourier transform       */
  !/*      from F-F-Physical   (array fou)                             */
  !/*  to/ from Phys-Phys-Phys (array phys)                            */
  !/*                                                                  */
  !/*       iopt >=0    ===>  inverse transforms( fou ---> fis)        */
  !/*       iopt < 0    ===>  direct  transforms( fis ---> fou)        */
  !/*                                                                  */
  !/*       does the expanding and compacting itself                   */
  !/*                                                                  */
  !/*  NOTE:                                                           */
  !/*    fis ---> fou : supposes high frecuency modes are not          */
  !/*                   needed (dealiasing) and throws them away       */
  !/*                                                                  */
  !/********************************************************************/

subroutine fourxz(fou,phys,iopt,jee,jbf,jef)   
  use ctesp
  use temporal
  use point
  use fourthings
  
  implicit none
  include "mpif.h"

  integer iopt,jee,jbf,jef,mpiid
  complex*16 fou (0:nz2,jbf:jef),zero        !zy plane in fourier R8
  complex*16 phys(0:ngz,jbf:jef)              !zy plane in Phys R8
  integer j,k
  real*8 tmpp1,tmpp2
  
  zero=dcmplx(0d0,0d0)
 
  	if (mpiid2.eq.0) then
	     tmpp1 = MPI_WTIME()
	endif

  if (iopt<0) then     
     do j=jbf,jef         
        call dfftw_execute_dft_r2c(planf,phys(0,j),phys(0,j))    
        fou(0:nz2,j)=dnf*phys(0:nz2,j)   
     enddo
  else      
     do j=jbf,jef           
        phys(0:nz2,j)=fou(0:nz2,j)               	
        phys(nz2+1:ngz,j)=zero      !dealiasing          
        call dfftw_execute_dft_c2r(planb,phys(0,j),phys(0,j)) 
     enddo
  endif

 if (mpiid2.eq.0) then
	     tmpp2 = MPI_WTIME()
	     !$OMP CRITICAL
	     tmp11 = tmp11 + abs(tmpp2-tmpp1)/nthreads
	     !$OMP END CRITICAL
	  endif
end subroutine fourxz

! ==============================================================
!              initialize rft_w
!     uses fftw_3.0, supposed to work as the old rftsingle 
!              tested only in aeolos 
!     compile with flags: 
!     
!     ifort -O3 -tpp7 -xW -align dcommons -c rftw3.f \
!         -lfftw3f -lm -I/usr/local/include
!
! ==============================================================
subroutine rfti(n)
  use ctesp
  use fourthings
 
  implicit none
  include "fftw3.f"

  integer*4  n

  allocate (fdum(0:n/2))
  nf=n
  nb=n
  dnf=1d0/n
  dnb=1d0

  call dfftw_plan_dft_r2c_1d(planf,n,fdum,fdum,FFTW_MEASURE )
  call dfftw_plan_dft_c2r_1d(planb,n,fdum,fdum,FFTW_MEASURE )
 
  deallocate(fdum)
end subroutine rfti



! ==============================================================
!              initialize cosft_w
!     uses fftw_3.0, supposed to work as the old cosftsingle 
!              tested only in aeolos 
!     compile with flags: 
!     
!     ifort -O3 -tpp7 -xW -align dcommons -c cosftw3.f \
!         -lfftw3f -lm -I/usr/local/include
!
! ==============================================================

module cosfthings
use omp_lib

  real*8, dimension(:), allocatable:: fdum,bdum
  real*8    dnf,dnb
  integer*4 nf,nb
  integer*8 planf,planb

end module cosfthings


subroutine cosfti(n)
  use cosfthings
  implicit none
  include "fftw3.f"

  integer*4  n
  nb=n
  nf=n
  dnf=1d0/n
  dnb=5d-1
  allocate (fdum(n+2))
 
  call dfftw_plan_r2r_1d(planf,nf,fdum,fdum,FFTW_ReDFT10,FFTW_MEASURE )
  call dfftw_plan_r2r_1d(planb,nb,fdum,fdum,FFTW_REDFT01,FFTW_MEASURE )

  deallocate(fdum)
end subroutine cosfti


! ================================================================
!             the real  cosftsingle (good orientation) 
! ================================================================

subroutine cosftx(c,isa,mb,me,iopt)

  use cosfthings
  use ctesp
  use temporal
  implicit none
  include "mpif.h"

  integer isa,m,iopt,i,j,mb,me
  real*8 c(isa,mb:me),tmpp1,tmpp2
  include "fftw3.f"

	if (mpiid2.eq.0) then
	     tmpp1 = MPI_WTIME()
	endif
	
  if (iopt<0) then
     do j=mb,me
        call dfftw_execute_r2r(planf,c(1,j),c(1,j))
          c(1:nf,j) = dnf*c(1:nf,j)
     enddo
  else   
     do j=mb,me
        call dfftw_execute_r2r(planb,c(1,j),c(1,j))
          c(1:nb,j) = dnb*c(1:nb,j)
     enddo
  endif

if (mpiid2.eq.0) then
	     tmpp2 = MPI_WTIME()
	     !$OMP MASTER
	     tmp12 = tmp12 + abs(tmpp2-tmpp1)
	     !$OMP END MASTER
	  endif
end subroutine cosftx


!  ==========================================
!    OMP indexes for FFT & cosinus transform 
!  ==========================================

subroutine pointersfft(mpiid)
use point
use ctesp
use omp_lib
implicit none
integer mpiid
!$OMP PARALLEL
nthreads=OMP_GET_NUM_THREADS()
ompid=OMP_GET_THREAD_NUM()

 chunk1=((ny+1)+nthreads-1)/nthreads
 jbf1=ompid*chunk1+1
 jef1=min((ompid+1)*chunk1,ny+1)

 chunk2=(ny+nthreads-1)/nthreads
 jbf2=ompid*chunk2+1
 jef2=min((ompid+1)*chunk2,ny)

 chunk3=(mpv+nthreads-1)/nthreads
 mpvb=ompid*chunk3+1
 mpve=min((ompid+1)*chunk3,mpv)

  if (mpiid.eq.0) write(*,'(a40,7i10)') 'jbf1,jef1,jbf2,jef2,mpvb,mpve,ompid',&
     & jbf1,jef1,jbf2,jef2,mpvb,mpve,ompid
!$OMP END PARALLEL
end subroutine pointersfft
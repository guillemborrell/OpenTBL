module adjunto
integer ompid,nthreads,chunk1,chunk2,chunk3
!$OMP THREADPRIVATE(ompid,nthreads,chunk1,chunk2,chunk3)
integer jbf1,jef1,jbf2,jef2,mpvb,mpve
!$OMP THREADPRIVATE(jbf1,jef1,jbf2,jef2,mpvb,mpve)
endmodule adjunto


program pointersfft
use adjunto
use omp_lib
implicit none
integer, parameter:: ny=646,mpv=12212

! call omp_set_num_threads(4)

!$OMP PARALLEL
nthreads=OMP_GET_NUM_THREADS()
ompid=OMP_GET_THREAD_NUM()
if(ompid.eq.0) write(*,*) 'nthreads=',nthreads
 chunk1=((ny+1)+nthreads-1)/nthreads
 jbf1=ompid*chunk1+1
 jef1=min((ompid+1)*chunk1,ny+1)

 chunk2=(ny+nthreads-1)/nthreads
 jbf2=ompid*chunk2+1
 jef2=min((ompid+1)*chunk2,ny)

 chunk3=(mpv+nthreads-1)/nthreads
 mpvb=ompid*chunk3+1
 mpve=min((ompid+1)*chunk3,mpv)


 write(*,'(a40,7i10)') 'jbf1,jef1,jbf2,jef2,mpvb,mpve,ompid',&
     & jbf1,jef1,jbf2,jef2,mpvb,mpve,ompid
!$OMP END PARALLEL
endprogram  pointersfft
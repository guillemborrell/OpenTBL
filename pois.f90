!=============================================================
!   Solves the Poisson eq. for the pressure correction
!   using a cosine transform in x (with wavenumbers of 2nd ord. FD)  
!   and Fourier in z. Tridiagonal FD in y
!
!   Everything works in (zy) planes, except for the cosines
!
!  old version from MP Simens, storage changed for BGP by
!                     JJS, Dec 24/2009
!=============================================================
subroutine pois(ut,vt,wt,pt,res,rest,rt,varstep,mpiid,communicator)

  use point
  use alloc_dns,only:idx,idy,idxx,idyy,phiy,dy,y,kaz,kaz2,kmod,ayp
  use ctesp
  use omp_lib
  use temporal
  implicit none
  include 'mpif.h'
  integer,intent(in):: communicator
  ! ---------------------- I/O -------------------------------------!
  integer mpiid
  real*8 dt,varstep

  real*8,     dimension(nx,mpv) :: res    !<<<<<<<<<<<<
  complex*16, dimension(0:nz2,ny+1,ib:ie):: wt,rt,ut
  complex*16, dimension(0:nz2,ny,ib:ie)  :: pt,vt,rest
  ! -------------------------- Work Arrays -------------------------!
  real*8  aypr(3,ny-1),dpdyH,dpdyH1,a
  integer i,j,l,k,kk,k2
  ! --------------------- MPI workspaces -----------------------------!
  integer istat(MPI_STATUS_SIZE),ierr,comm,countu,countv,tipo
  ! ----------------------------------------------------------------!

  countu=(nz2+1)*(ny+1)
  countv=(nz2+1)*ny
  comm = communicator
  tipo=MPI_COMPLEX16

  ! --- compute the divergence, we are in (zy) 
  do i=ib0,ie
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        rest(:,j,i) = wt(:,j+1,i)*kaz +(vt(:,j+1,i)-vt(:,j,i))*idy(j)
     enddo
  enddo

  ! --- add du/dx -------------
  do i=ib+1,ie
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        rest(:,j,i)=rest(:,j,i)+idx*(ut(:,j+1,i)-ut(:,j+1,i-1))
     enddo
  enddo


  if (mpiid2.eq.0) tm1 = MPI_WTIME()

  if (mpiid.eq.0) then
     call MPI_SEND(ut(0,1,ie),countu,tipo,mpiid+1,0,comm,istat,ierr)
  elseif (mpiid.eq.pnodes-1) then
     call MPI_RECV(rt(0,1,ib),countu,tipo,mpiid-1,0,comm,istat,ierr)
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        rest(:,j,ib)=rest(:,j,ib)+idx*(ut(:,j+1,ib)-rt(:,j+1,ib))
     enddo
  else
     call MPI_SENDRECV(ut(0,1,ie),countu,tipo,mpiid+1,0,  &
          &                      rt(0,1,ib),countu,tipo,mpiid-1,0,  comm,istat,ierr)
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        rest(:,j,ib)=rest(:,j,ib)+idx*(ut(:,j+1,ib)-rt(:,j+1,ib))
     enddo
  endif

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp15 = tmp15 + abs(tm2-tm1)
  endif

  !  ----  go to lines, transform, and go back to planes ----
  call chp2x(res,rest,rt,mpiid,ny,comm)
  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()
  endif

!=============================================================
!=============================================================
  !$OMP PARALLEL DEFAULT(SHARED)
  call cosftx(res(2,mpvb),nx,mpvb,mpve,-1) 
  !$OMP END PARALLEL
!=============================================================
!=============================================================


  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp10 = tmp10 + abs(tm2-tm1)
  endif
  call chx2p(res,rest,rt,mpiid,ny,comm)

  rest(0,:,ib:ie)=real(rest(0,:,ib:ie),kind=8) !Ensuring the 0th mode is Real



  ! --------  solve poisson in (zy) planes  (using equiv. wavenumbers)
  do i = ib0,ie 
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,aypr) SCHEDULE(STATIC)        
     do k = 0,nz2
        aypr = ayp(:,1:ny-1)
        aypr(2,1:ny-1) = ayp(2,1:ny-1)-kaz2(k)-kmod(i-1)
        if (i==2 .and. k==0) then
            aypr(2,ny-1) = 2d0*ayp(2,ny-1)
            rest(0,ny-1,2)=0d0
        endif 

        call soltriy(rest(0,1,i),rest(0,1,i),k,aypr)
     enddo
     
     
!      if (i .eq.2) then
!          dpdyH  = -rest(0,ny-1,1)/(dy(ny)+dy(ny-1))
!          dpdyH1 = -phiy(ny-1)/(dy(ny)+dy(ny-1))
!          a      = -dpdyH/dpdyH1
!          write(*,*) 'masas',a,dpdyH,dpdyH1
!      endif
  enddo

 

  !  ----  go to lines, back-transform, and go back to planes ----
  call chp2x(res,rest,rt,mpiid,ny,comm)

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()
  endif
  !$OMP PARALLEL DEFAULT(SHARED)
  call cosftx(res(2,mpvb),nx,mpvb,mpve,1) 
  !$OMP END PARALLEL
  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp10 = tmp10 + abs(tm2-tm1)
  endif
  call chx2p(res,rest,rt,mpiid,ny,comm)

  
  do i =ib0,ie
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j = 2,ny
        pt(:,j,i) = varstep*rest(:,j-1,i)+pt(:,j,i)    ! ------- update pressure
        wt(:,j,i) = wt(:,j,i)-kaz*rest(:,j-1,i)        ! ------- update velocities 
     enddo
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j = 2,ny-1
        vt(:,j,i) = vt(:,j,i)-idyy(j)*(rest(:,j,i)-rest(:,j-1,i))
     enddo
  enddo
  
  ! ------- update u+dp/dx ---
  do i=ib0,ie-1
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        ut(:,j+1,i)=ut(:,j+1,i)-idxx*(rest(:,j,i+1)-rest(:,j,i))
     enddo
  enddo

  if (mpiid2.eq.0) tm1 = MPI_WTIME()

  if (mpiid.eq.pnodes-1) then
     call MPI_SEND(rest,countv,tipo,mpiid-1,1,comm,istat,ierr)
  elseif (mpiid.eq.0) then
     call MPI_RECV(rt,countv,tipo,mpiid+1,1,comm,istat,ierr)
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        ut(:,j+1,ie)=ut(:,j+1,ie)-idxx*(rt(:,j,ib)-rest(:,j,ie))
     enddo
  else
     call MPI_SENDRECV(rest,countv,tipo,mpiid-1,1,  &
          &                        rt,countv,tipo,mpiid+1,1,  comm,istat,ierr)
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) SCHEDULE(STATIC)
     do j=1,ny-1
        ut(:,j+1,ie)=ut(:,j+1,ie)-idxx*(rt(:,j,ib)-rest(:,j,ie))
     enddo
  endif

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp23 = tmp23 + abs(tm2-tm1)-(tp2-tp1)
  endif

! 
!  if(paso.eq.1) then
!       if(mpiid.eq.0) write(*,*) 'WRITING THE K=0 XY PLANE TO A FILE FOR U,V & W IN POISON'
!       do i=ib,ie   
!          pdiv(1:ny,i)=real(ut(0,1:ny,i),kind=8) !Each node copy a piece of the array
!       enddo
!       call MPI_ALLREDUCE(MPI_IN_PLACE,pdiv,ny*nx,MPI_real8,MPI_SUM,comm,ierr)
!       if(mpiid.eq.0) write(28) pdiv(1:ny,1:nx) 
!       pdiv=0d0
! 
!       do i=ib,ie   
!          pdiv(1:ny,i)=real(vt(0,1:ny,i),kind=8) !Each node copy a piece of the array
!       enddo
!       call MPI_ALLREDUCE(MPI_IN_PLACE,pdiv,ny*nx,MPI_real8,MPI_SUM,comm,ierr)
!       if(mpiid.eq.0) write(28) pdiv(1:ny,1:nx)
!       pdiv=0d0
! 
!       do i=ib,ie   
!          pdiv(1:ny,i)=real(wt(0,1:ny,i),kind=8) !Each node copy a piece of the array
!       enddo
!       call MPI_ALLREDUCE(MPI_IN_PLACE,pdiv,ny*nx,MPI_real8,MPI_SUM,comm,ierr)      
!       if(mpiid.eq.0) write(28) pdiv(1:ny,1:nx) 
!       pdiv=0d0
! 
!       if(mpiid.eq.0) write(*,*) 'WRITING THE K=0 XY PLANE TO A FILE FOR U,V & W IN POISON.............DONE'
!  endif


end subroutine pois

!/*********************************************************************/
!/*                                                                   */
!/*               Resolvedor de un sistema tridiagonal                */
!/*              resuelve arrays 2D por columnas ...                  */
!/*                                                                   */
!/*  entrada: a1 a2 a3          us    u                               */
!/*              a1 a2 a3       us =  u                               */
!/*                 a1 a2 a3    us    u                               */
!/*                                                                   */
!/*           nstr      Maxima primera dimension de u y us            */
!/*           nlen      Numero de sistemas a resolver                 */
!/*           neq       Dimension del sistema a resolver              */
!/*           iflag=0   Primera llamada (factorizacion y resolucion)  */
!/*                =1   Solo resolucion                               */
!/*                                                                   */
!/*   salida: us        Solucion                                      */
!/*           iflag=1                                                 */
!/*                                                                   */
!/*                                                                   */
!/*********************************************************************/
subroutine soltriy(u,us,ii,a) 	!ii=column to be solved
  use ctesp
  implicit none

  integer ii,j
  complex*16 u(0:nz2,ny1),us(0:nz2,ny1)
  real*8     a(3,ny1),d

  !       lu-decomposition
  do  j = 2,ny1
     d=1d0/a(2,j-1)
     a(2,j-1) = d
     a(1,j  ) =        - a(1,j)*d
     a(2,j  ) = a(2,j) + a(1,j)*a(3,j-1)
  enddo
  a(2,ny1)=1d0/a(2,ny1)

  !       backsubstitution

  us(ii,1) = u(ii,1)
  do j = 2,ny1
     us(ii,j) = u(ii,j)+a(1,j)*us(ii,j-1)
  enddo

  us(ii,ny1) = us(ii,ny1)*a(2,ny1)
  do j = ny1-1,1,-1
     us(ii,j) = (us(ii,j)-a(3,j)*us(ii,j+1))*a(2,j)
  enddo
end subroutine soltriy


! ====================================================
!             upen_int -> result, works in x(lines) 
!       can work in place,         JJS 12/09
! ====================================================
! ===================================================================
!   ADDS explicit viscous term in x ---  orientation in x(lines)
!               works in place
!   Called by:  rhsp.f90,          JJS 12/09
! ===================================================================
!//////////////////// JSS MAY 2010 ///////
!interpyy() now works in place
!differzy() now works in place
!////////////////////////////////////////

!JSS: Also, apart of computing the viscous term the x convective term is computed.
!  upen_int: for computing the convective terms
! upencil: for computing the viscous terms
!result: Result

! icon,econ   --> convective    [Ac]X' =[Bc]X  == [icon]X'=[econ]X    
! iconv,econv --> viscous       [Av]X''=[Bv]X  == [iconv]X'=[econv]X   
! 
! dcbx  = Boundary coefficients for convective terms
! cofvbx= Boundary coefficients for viscous terms
! nlin=number of pencils
!d/dx: chv depending of 'v,w' or 'u' ------> For convective terms (chv for u=1, for v,w=2)
!d/dx: chv depending of (uv),(upen_int),(uw)
! vel=0: 'v,w'; vel.ne.0 u      ------> For viscous terms
! rex=-1/Reynolds

!ie:
! ===========d/dX==========================
! differx(geni,result ,icon,dcbx,econ      ,chv,jee)
! call differx(u_x ,result,dcxu,dcbx,cofcxu,1,ny+1)
! call differx(v_x ,result,dcxv,dcbx,cofcxv,2,ny)
! call differx(w_x ,result,dcxv,dcbx,cofcxv,2,ny+1)
! call differx   (v,result,dcxu,dcbx,cofcxu,1,ny)
! call differx   (w,result,dcxu,dcbx,cofcxu,1,ny+1)
! 
! ===========d/dY==========================
! differzy(geni,result,icon,dcby,econ        ,chv,npter1,fwb1,fwb2  ,jee)
! call differzy(u  ,resultu,dcyv,dcby,cofcyv,  1,   1,     0,   1,  ny+1)
! call differzy(w  ,resultu,dcyv,dcby,cofcyv,  1,   1,     0,   1,  ny+1)
! call differzy(u_y,rhsu,dcyu,dcby,cofcyu,  3,   1,    -1,   0,  ny+1)
! call differzy(v_y,rhsv,dcyv,dcby,cofcyv,  1,   1,     0,   1,  ny  )
! call differzy(w_y,rhsw,dcyu,dcby,cofcyu,  3,   1,    -1,   0,  ny+1)


!i.e: call difvisxx(wki1,resultu,wki1,dcxu,vixu,dcbx,cofcxu,cofvxu,cofvbx,1,mpu,1,rex)


subroutine difvisxx_2(upen_int,upencil,result,icon,iconv,dcbx,econ,econv,cofvbx,chv,nlin,vel,rex)
  use ctesp_2
  use temporal_2
  use alloc_dns_2,only:idx
  implicit none
  include "mpif.h"

  integer npter1,vel
  real*8 icon(3,nx),econ(2,nx),iconv(3,1:nx),econv(3,1:nx)
  real*8 dcbx(5,4),cofvbx(4,4),rex
  real*8 result(nx,nlin),upen_int(nx,nlin),upencil(nx,nlin),bb(nx)
  integer chv1,chv,i,j,k,l,nlin,ic     

  ! For the convective terms 
  chv1 = chv+2

  ! For the viscous terms 
  if (vel==0) then
     ic=1
  else
     ic=3
  end if

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()
  endif

!===========Start computing the convective terms and copy it to result()
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,i,bb) 
  !$OMP DO SCHEDULE(STATIC)
  do k =1,nlin
     !---------[A]du=[B]u Computing [B]u=[Bu]
     do i = 2,nx-1
        bb(i)=econ(1,i)*upen_int(i-1,k)+econ(2,i)*upen_int(i,k)
     enddo

     bb(1) = dcbx(chv,1)*upen_int(1,k)+ dcbx(chv,2)*upen_int(2,k)+      & 
          &  dcbx(chv,3)*upen_int(3,k)+ dcbx(chv,4)*upen_int(4,k) 
     bb(nx) = dcbx(chv1,1)*upen_int(nx,k)+ dcbx(chv1,2)*upen_int(nx-1,k)+  & 
          &   dcbx(chv1,3)*upen_int(nx-2,k)+ dcbx(chv1,4)*upen_int(nx-3,k)

     !---------SOLVING [A]du=[Bu]----------------------
     do j = 2,nx
        bb(j) = bb(j) + icon(1,j)*bb(j-1)
     enddo

     bb(nx) = bb(nx)*icon(2,nx)
     do j = nx-1,1,-1
        bb(j) = (bb(j)-icon(3,j)*bb(j+1))*icon(2,j)
     enddo
     !-------------Pade of the convective terms DONE!

     !IMPOSE outflow BC Uinf*du/dx for the last plane:
     bb(nx)=Uinfinity*(upencil(nx,k)-upencil(nx-1,k))*idx

     result(:,k)=bb !convective term saved in result() buffer

!===========Start computing the viscous terms and add it to result()  
     !---------[A]ddu=[B]u Computing [B]u=[Bu]
     do i = 2,nx-1
        bb(i) = econv(1,i)*upencil(i+1,k)+ econv(2,i)*upencil(i,k)+  &
             &     econv(3,i)*upencil(i-1,k)               
     enddo

     bb(1)  = cofvbx(ic,1)*upencil(1,k) + cofvbx(ic,2)*upencil(2,k) +       &   
          &      cofvbx(ic,3)*upencil(3,k) + cofvbx(ic,4)*upencil(4,k)

     bb(nx) = cofvbx(ic+1,1)*upencil(nx,k) + cofvbx(ic+1,2)*upencil(nx-1,k)+  &   
          &      cofvbx(ic+1,3)*upencil(nx-2,k) + cofvbx(ic+1,4)*upencil(nx-3,k)

     !---------SOLVING [A]ddu=[Bu]----------------------
     do i = 2,nx
        bb(i) = bb(i) + iconv(1,i)*bb(i-1)
     enddo

     bb(nx) = bb(nx)*iconv(2,nx)
     do i = nx-1,1,-1
        bb(i) = (bb(i)-iconv(3,i)*bb(i+1))*iconv(2,i)
     enddo

     bb(nx)=0d0  !! IMPOSE inviscid output BC 
     result(:,k)=result(:,k)+rex*bb  !convective + viscous terms
  enddo
  !$OMP END PARALLEL 

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp5 = tmp5 + abs(tm2-tm1)
  endif

end subroutine difvisxx_2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!             used in the (zy) part of rhsp 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!JSS May 20th 2010
!differzy now works in place
!
subroutine differzy_2(geni,result,icon,dcby,econ,chv,npter1,fwb1,fwb2,jee)
  use ctesp_2
  use point_2
  use temporal_2
  implicit none
  include "mpif.h"

  integer jee,npter1,kk,k2
  real*8 icon(3,npter1:jee),econ(2,1:ny)
  real*8 dcby(4,4),tm1p,tm2p

  real*8 result(0:2*nz2+1,ny+1),geni(0:2*nz2+1,jee)
  real*8:: buf_top(0:2*nz2+1),buf_bottom(0:2*nz2+1) !in order to work in place
  integer chv,fwb1,fwb2,i,j,k,l,cha

  cha = chv + 1

  if (mpiid2.eq.0) then
     tm1p = MPI_WTIME()
  endif

  !$OMP DO SCHEDULE(STATIC)
  do kk=0,2*nz2+1,blockl
     k2=min(2*nz2+1,kk+blockl-1)

     !at the boundaries
     do k =kk,k2 
        buf_bottom(k) =  (dcby(chv,1)*geni(k,1) +      & 
             &                     dcby(chv,2)*geni(k,2) +      & 
             &                     dcby(chv,3)*geni(k,3) +      & 
             &                     dcby(chv,4)*geni(k,4))         

        buf_top(k)    = (dcby(cha,1)*geni(k,ny)   + & 
             &                      dcby(cha,2)*geni(k,ny-1) + & 
             &                      dcby(cha,3)*geni(k,ny-2) + & 
             &                      dcby(cha,4)*geni(k,ny-3))
     enddo

     if(fwb2.eq.1) then
        do j = 2,jee-1 
           do k =kk,k2  
              result(k,j) = econ(2,j)*geni(k,j+fwb2)+ econ(1,j)*geni(k,j+fwb1)   
           enddo
        enddo
     elseif(fwb2.eq.0) then
        do j = jee-1,2,-1 
           do k =kk,k2  
              result(k,j) = econ(2,j)*geni(k,j+fwb2)+ econ(1,j)*geni(k,j+fwb1)   
           enddo
        enddo
     endif

     result(kk:k2,1)  =buf_bottom(kk:k2)
     result(kk:k2,jee)=buf_top(kk:k2)

     ! resolution of tridiagonal 
     do j = 2,jee
        do k =kk,k2 
           result(k,j) = result(k,j) + icon(1,j)*result(k,j-1)  
        enddo
     enddo

     do k =kk,k2 
        result(k,jee) = result(k,jee)*icon(2,jee) 
     enddo

     do j = jee-1,1,-1
        do k =kk,k2 
           result(k,j) = (result(k,j)-icon(3,j)*result(k,j+1))*icon(2,j) 
        enddo
     enddo
  enddo

  if (mpiid2.eq.0) then
     tm2p = MPI_WTIME()
     !$OMP CRITICAL
     tmp8 = tmp8 + abs(tm2p-tm1p)/nthreads
     !$OMP  END CRITICAL
  endif

end subroutine differzy_2

! ******************************************************************
!                                                                  
!    Interpolates the velocities u,v in x using a Pade scheme.     
!    Thus Au=b is solved where u is the interpolated velocity.              
!    gen   ==> input velocity to be interpolated                   
!    geni  ==> output interpolated velocity.                       
!    icon  ==> coefficients of the lefthand side (above A)      
!                  of the scheme                                   
!    econ     ==> coefficient to calculate b in the above            
!    inbx    ==> boundary scheme coefficients                       
!    loguv = 1:  'u'. loguv=0: 'v,w'
!                                                                  
!     works in (x-lines).      M.Simens 03,   JJS 12/09
! *******************************************************************
subroutine interpxx_2(gen,geni,icon,econ,inbx,chi,nlin,loguv)
  use ctesp_2
  use temporal_2
  use omp_lib
  implicit none
  include "mpif.h"

  real*8  icon(3,nx),econ(2,nx),inbx(4,4)
  real*8  gen (nx,nlin),geni(nx,nlin), bb(nx)
  integer i,j,k,l,loguv,chi,chi2,nlin

  if (loguv .eq. 1) then
     chi2 = chi+1
  else
     chi2 = chi+3
  endif

  if (mpiid2.eq.0) then
     tm1 = MPI_WTIME()
  endif

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(k,i,bb)
  do k =1,nlin

     do i = 2,nx-1
        bb(i)=econ(1,i)*gen(i,k)+econ(2,i)*gen(i+1,k)
     enddo

     !boundary scheme.
     bb(1) = inbx(chi2,1)*gen(1,k) + inbx(chi2,2)*gen(2,k) +   &  
          &  inbx(chi2,3)*gen(3,k) + inbx(chi2,4)*gen(4,k)            

     bb(nx)= inbx(chi,1)*gen(nx,k) + inbx(chi,2)*gen(nx-1,k) + &  
          &  inbx(chi,3)*gen(nx-2,k)+ inbx(chi,4)*gen(nx-3,k)

     ! solving tridiagonal
     do i = 2,nx
        bb(i) = bb(i) + icon(1,i)*bb(i-1)
     enddo

     bb(nx) = bb(nx)*icon(2,nx)
     do i = nx-1,1,-1
        bb(i) = (bb(i)-icon(3,i)*bb(i+1))*icon(2,i)
     enddo
     geni(:,k)=bb

  enddo

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp6 = tmp6 + abs(tm2-tm1)
  endif
end subroutine interpxx_2

! ====================================================
!             u -> result, works in x(lines)
!       can work in place,         JJS 12/09
! ====================================================
subroutine differxx_2(u,result,icon,dcbx,econ,chv,nlin)

  use ctesp_2
  implicit none

  integer npter1
  real*8 icon(3,nx),econ(2,nx)
  real*8 dcbx(5,4)
  real*8 result(nx,nlin),u(nx,nlin),bb(nx)
  integer chv1,chv,i,j,k,l,nlin

  chv1 = chv+2
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,k,i,bb)
  do k =1,nlin
     do i = 2,nx-1
        bb(i)=econ(1,i)*u(i-1,k)+econ(2,i)*u(i,k)
     enddo

     bb(1) = dcbx(chv,1)*u(1,k)+ dcbx(chv,2)*u(2,k)+      &
          &  dcbx(chv,3)*u(3,k)+ dcbx(chv,4)*u(4,k)
     bb(nx) = dcbx(chv1,1)*u(nx,k)  + dcbx(chv1,2)*u(nx-1,k)+  &
          &   dcbx(chv1,3)*u(nx-2,k)+ dcbx(chv1,4)*u(nx-3,k)

     do j = 2,nx
        bb(j) = bb(j) + icon(1,j)*bb(j-1)
     enddo

     bb(nx) = bb(nx)*icon(2,nx)
     do j = nx-1,1,-1
        bb(j) = (bb(j)-icon(3,j)*bb(j+1))*icon(2,j)
     enddo
     result(:,k)=bb
  enddo
end subroutine differxx_2

! *******************************************************************
!                                                                  
!    Interpolates the velocities u,v,w in y using a Pade scheme.     
!    Thus Au=b is solved where u is the final interpolated         
!               velocity.                                                       
!    gen   ==> input velocity to be interpolated                   
!    geni  ==> output interpolated velocity.                       
!    icon  ==> coefficients of the lefthand side (above A)      
!                  of the scheme                                   
!    econ   ==> coefficient to calculate b in the above            
!    inby   ==> boundary scheme coefficients                       
!     
!    loguv = 0:  'u,w'. loguv=1: 'v'                                                             
!    M.Simens 03                                                   
! *******************************************************************
!
!JSS May 20th 2010
!interpyy now works in place
subroutine interpyy_2(gen,geni,icon,econ,inby,jee,loguv,dim,m)     
  use ctesp_2
  use point_2
  use temporal_2
  implicit none
  include "mpif.h"

  integer jee,i,j,k,l,loguv,m,dim,ii,i2    
  real*8 icon(3,1:ny),inby(4,4),econ(2,1:ny),gen(dim,1:jee),         &
       &      geni(dim,1:ny+1),tm1p,tm2p

  real*8:: buffer(dim) !in order to work in place

  if (mpiid2.eq.0) then
     tm1p = MPI_WTIME()    
  endif

  if (loguv .eq. 0) then    ! ---  u & w
     !$OMP DO SCHEDULE(STATIC)
     do ii=1,m,blockl
        i2=min(m,ii+blockl-1)
        do k =ii,i2                ! boundary schemes.
           geni(k,1)  = 0d0  ! for no slip boundary conditions
           buffer(k)  = inby(1,1)*gen(k,ny+1) + inby(1,2)*gen(k,ny)  + &
                &         inby(1,3)*gen(k,ny-1) + inby(1,4)*gen(k,ny-2)
        enddo
        do j = 2,ny-1 
           do k =ii,i2         ! calculating b
              geni(k,j)  = econ(1,j)*gen(k,j) + econ(2,j)*gen(k,j+1)
           enddo
        enddo
        geni(ii:i2,ny)=buffer(ii:i2)
     enddo
     !$OMP END DO NOWAIT
  else                     ! --- v
     !$OMP DO SCHEDULE(STATIC)
     do ii=1,m,blockl
        i2=min(m,ii+blockl-1)
        do k =ii,i2                ! boundary schemes.
           buffer(k)  = inby(3,1)*gen(k,1) + inby(3,2)*gen(k,2)  +  &
                &         inby(3,3)*gen(k,3) + inby(3,4)*gen(k,4)
           geni(k,ny) = inby(4,1)*gen(k,ny) + inby(4,2)*gen(k,ny-1)+ &
                &         inby(4,3)*gen(k,ny-2)+inby(4,4)*gen(k,ny-3)
        enddo
        do j = ny-1,2,-1 !inverse loop order in order to work in place
           do k =ii,i2         ! calculating b
              geni(k,j) = econ(2,j)*gen(k,j)+econ(1,j)*gen(k,j-1)
           enddo
        enddo
        geni(ii:i2,1 )=buffer(ii:i2)
     enddo
     !$OMP END DO NOWAIT
  endif
  !  --- solve the tridiagonal ----
  !$OMP DO SCHEDULE(STATIC)
  do ii=1,m,blockl

     i2=min(m,ii+blockl-1)
     do j = 2,ny
        do k =ii,i2     
           geni(k,j) = geni(k,j) + icon(1,j)*geni(k,j-1)
        enddo
     enddo

     do k =ii,i2     
        geni(k,ny) = geni(k,ny)*icon(2,ny)
     enddo

     do j = ny-1,1,-1
        do k =ii,i2     
           geni(k,j) = (geni(k,j)-icon(3,j)*geni(k,j+1))*icon(2,j)
        enddo
     enddo
  enddo

  if (mpiid2.eq.0) then
     tm2p = MPI_WTIME()
     !$OMP CRITICAL
     tmp4 = tmp4 + abs(tm2p-tm1p)/nthreads
     !$OMP  END CRITICAL
  endif
end subroutine interpyy_2


! ===================================================================
! compute explicit viscous term in y ---  orientation in (zy)
!   Called by:  rhsp.f90,          JJS 12/09
! ===================================================================
subroutine vistyy_2(u,result,icon,econ,cofvby,jee)
  use ctesp_2
  use point_2
  use temporal_2
  implicit none
  include "mpif.h"

  integer jee,ii,i2
  real*8 icon(3,jee),econ(3,ny),cofvby(4,4)
  real*8 u(nz1,jee),result(nz1,ny+1),tm1p,tm2p
  integer i,j,k,ic

  if (jee==ny) then
     ic=1
  else
     ic=3
  end if

  if (mpiid2.eq.0) then
     tm1p = MPI_WTIME()
  endif

  !----- right hand side ---------
  !$OMP DO SCHEDULE(STATIC)
  do ii=1,nz1,blockl

     i2=min(nz1,ii+blockl-1)
     do k =ii,i2  
        result(k,1)=cofvby(ic,1)*u(k,1)+cofvby(ic,2)*u(k,2) + cofvby(ic,3)*u(k,3)+cofvby(ic,4)*u(k,4)
     enddo

     do j = 2,jee-1
        do k =ii,i2  
           result(k,j) = econ(1,j)*u(k,j+1)+econ(2,j)*u(k,j)+ econ(3,j)*u(k,j-1)
        enddo
     enddo

     do k =ii,i2  
        result(k,jee)=cofvby(ic+1,1)*u(k,jee)+cofvby(ic+1,2)*u(k,jee-1)+ & 
             &  cofvby(ic+1,3)*u(k,jee-2)+cofvby(ic+1,4)*u(k,jee-3)
     enddo

     !----- backsubstitution in place
     do j = 2,jee
        do k =ii,i2  
           result(k,j) = result(k,j) + icon(1,j)*result(k,j-1)
        enddo
     enddo

     do k =ii,i2  
        result(k,jee) = result(k,jee)*icon(2,jee)
     enddo

     do j = jee-1,1,-1
        do k =ii,i2  
           result(k,j) = (result(k,j)-icon(3,j)*result(k,j+1))*icon(2,j)
        enddo
     enddo
  enddo

  if (mpiid2.eq.0) then
     tm2p = MPI_WTIME()
     !$OMP CRITICAL
     tmp7 = tmp7 + abs(tm2p-tm1p)/nthreads
     !$OMP  END CRITICAL
  endif
end subroutine vistyy_2


! ===========================================================
!               new version, works in zy 
!    icon  ==> coefficients of the lefthand side (above A)      
!                  of the scheme                                   
!    econ     ==> coefficient to calculate b in the above 
!
!res=input	u=output
!
! ./rhsp.f90:284:call impl(u,wki1,vyui,cofvyu,kb,ke,ny+1,rex,dt,m,rkdv,1)
! ./rhsp.f90:285:call impl(w,wki3,vyui,cofvyu,kb,ke,ny+1,rex,dt,m,rkdv,3)
! ./rhsp.f90:286:call impl(v,wki2,vyvi,cofvyv,kb,ke,ny  ,rex,dt,m,rkdv,0)
! ===========================================================
subroutine implzy_2(u,res,icon,econ,jee,var1)
  use ctesp_2
  use temporal_2
  implicit none
  include "mpif.h"

  integer jee,kk,k2
  real*8, dimension(0:2*nz2+1,jee):: u,res
  real*8  icon(3,jee),econ(3,ny)
  real*8  a1(ny+1),a2(ny+1),a3(ny+1),var1,d
  integer i,j,m,k,l,vel

!!!!!!!!!  hopefully, u already has the boundary conditions in j=1 and jee 
  if (mpiid2.eq.0) tm1 = MPI_WTIME()

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,kk,k,k2)
  !$OMP DO SCHEDULE(STATIC)
  do j = 2,jee-1
     u(:,j)=icon(1,j)*res(:,j-1)+icon(2,j)*res(:,j)+icon(3,j)*res(:,j+1)

     a1(j) = icon(1,j) + var1*econ(3,j)
     a2(j) = icon(2,j) + var1*econ(2,j)
     a3(j) = icon(3,j) + var1*econ(1,j)
  enddo

  !$OMP WORKSHARE
  u(:,2)     = u(:,2) - a1(2)*u(:,1)
  u(:,jee-1) = u(:,jee-1) - a3(jee-1)*u(:,jee)
  a1(2)     = 0d0
  a3(jee-1) = 0d0
  !$OMP END WORKSHARE

  !$OMP SINGLE
  do  j=2+1,jee-1  
     d      =1d0/a2(j-1)
     a2(j-1)=d
     a1(j  )=     -a1(j)*d
     a2(j  )=a2(j)+a1(j)*a3(j-1)
  enddo
  a2(jee-1)=1d0/a2(jee-1)
  !$OMP END SINGLE

  ! backsubstitution
  !$OMP DO SCHEDULE(STATIC)
  do kk=0,2*nz2+1,blockl
     k2=min(2*nz2+1,kk+blockl-1)

     do j = 3,jee-1
        do k=kk,k2
           u(k,j) = u(k,j)+a1(j)*u(k,j-1)
        enddo
     enddo

     do k=kk,k2
        u(k,jee-1) = u(k,jee-1)*a2(jee-1)
     enddo

     do  j = jee-2,2,-1
        do k=kk,k2
           u(k,j) = (u(k,j)-a3(j)*u(k,j+1))*a2(j)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL

  if (mpiid2.eq.0) then
     tm2 = MPI_WTIME()
     tmp9 = tmp9 + abs(tm2-tm1)
  endif
end subroutine implzy_2






!!4th order interpolator for the pressure
!JSS January 2011
!ACHTUNG!!! DOES NOT WORK IN PLACE
!dp/dn=nu*d2V/dy2 @th wall
subroutine interp_pressure_2(f,fi,n,v)
  use ctesp_2
  use alloc_dns_2, only: y,l_weight,ldyy,re
  use point_2 
  implicit none
  include "mpif.h"
  integer:: n,j,k,kk,k2
  real*8, dimension(nz1,n):: f,fi,v
  real*8:: dvdy2

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,kk,k,k2,dvdy2)  
  !$OMP DO SCHEDULE(STATIC)
  do kk=1,nz1,blockl
     k2=min(nz1,kk+blockl-1)
     do k =kk,k2
        dvdy2=dot_product(ldyy,v(k,1:5))/re
        f(k,1)=(dvdy2-dot_product(l_weight(2:5,0),f(k,2:5)))/l_weight(1,0); 
     enddo
     do j=1,2
	do k=kk,k2
           fi(k,j)=dot_product(l_weight(1:5,j),f(k,j:j+4))
	enddo
     enddo
     do j=3,n-2
	do k=kk,k2
           fi(k,j)=dot_product(l_weight(1:5,j),f(k,j-2:j+2))
	enddo
     enddo
     do j=n-1,n
	do k=kk,k2
           fi(k,j)=dot_product(l_weight(1:5,j),f(k,j-4:j))
	enddo
     enddo
  enddo
  !$OMP END PARALLEL
end subroutine interp_pressure_2







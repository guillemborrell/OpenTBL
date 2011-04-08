! ****************************************************************************
! in this routine all necessary constants (constant arrays) are calculated
! inflow velocity is a constant found in routine inflowp.f
!
! Interpolation coefficients are derived.
! The Boundary Conditions are computed using LU decomposition.

! ****************************************************************************

subroutine coef(mpiid)
  use alloc_dns
  use point
  use statistics,only:hy
  use ctesp
  use num_nodes

  implicit none      

  include 'mpif.h'

  integer i,j,ii,n,np,flagv,mpiid,im,n1,n2,k

  real*8 bcc(5),a,b,c,d,am,am1,invre,dtrex,dtrez,dtrey,fact,point0,daux
  real*8 fren(nx),aux,dxi(0:nx),xxi(0:nx),dyi(1:ny),yyi(0:ny),ayy,gamy
  !MPI workspaces
  integer istat(MPI_STATUS_SIZE),ierr,stencil
  real*8:: x1,x2,x3,x4,x5,xd,nm_aux,dn(5,ny+1),nm(5,ny+1),yp(1:ny+1),xl(5)
  real*8:: h(ny+1),hi(ny+1),e,f,g

  pi = 4d0*atan(1d0)
  np = 5

  do j=1,ny-1
     hy(j) = (y(j+1)-y(j-1))/2.5d0/2d0
  enddo

  hy(0)  = (y(1   )-y(0 ))/2.5d0
  hy(ny) = (y(ny)-y(ny-1))/2.5d0



  !===================================================================
  ! COMPACT FINITE DIFFERENCES SCHEMES. JSS MARCH 011
  !
  !! As in: "COMPACT FINITE DIFFERENCE SCHEMES ON NON-UNIFORM MESHES".
  !! Gamet, Ducros, Nicoud, Poinsot)
  !! (I am doing high order 6th order CFD...as costly as 4th)
  !! Close to the wall, the order is improved using a bigger stencil
  !! (Mark: 4 points stencils, Me: 4,5 point stencils)
  !===================================================================
  
  
  !!====================================================
  !!       D/DY using an non-uniform mesh (in place)
  !!====================================================

  fd_dvdy=0d0
  h(1:ny+1)=y(1:ny+1)-y(0:ny)

  !BOUNDARY SCHEME FOR I=1
  a=h(2);b=h(2)+h(3);c=h(2)+h(3)+h(4)
  fd_dvdy(1,1) =-(a*b+a*c+2*b*c)/b/a/c
  fd_dvdy(3,1) =-a**2*c/(-c+b)/b/(-b+a)**2
  fd_dvdy(4,1) =a**2*b/c/(a-c)**2/(-c+b)
  fd_dvdy(2,1) =-(fd_dvdy(1,1)+fd_dvdy(3,1)+fd_dvdy(4,1))
  !Tridg:
  fd_dvdy(7,1)=1d0
  fd_dvdy(8,1) =b*c/(a-c)/(-b+a) !beta

  !BOUNDARY SCHEME FOR I=2
  fd_dvdy(1,2) =(-b+a)**2*(a-c)*(a*b+2*a*c+2*b*c)/a/b**3/c**2
  fd_dvdy(2,2) =1/a*(-3*a*b+2*b*c+5*a**2-4*a*c)/(a-c)/(-b+a)
  fd_dvdy(3,2) =a**2*(a-c)*(-2*a*c+3*a*b+4*b*c-5*b**2)/b**3/(-c+b)**2/(-b+a)
  fd_dvdy(4,2) =-sum(fd_dvdy(1:3,2))
  !Tridg:
  fd_dvdy(6,2) =-(-b+a)**2*(a-c)/c/b**2 !alpha
  fd_dvdy(7,2) =1d0
  fd_dvdy(8,2) =a**2*(a-c)/(-c+b)/b**2 !beta


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,a,b,c,d,e,f)
  do i=3,ny-2
     a=h(i-1); b=h(i); c=h(i+1); d=h(i+2); f=c+d; e=a+b

     fd_dvdy(1,i)=-c**2*b**2*f/(f+e)/e/(c**2+2*c*e+e**2)/(e**2-2*e*b+b**2)

     fd_dvdy(2,i)=f/b*e*c**2*(3*c*b*f-2*c*f*e&
          &+4*b**2*c-3*b*e*c+5*b**2*f-4*b*e*f&
          &+6*b**3-5*e*b**2)/(f**2+2*b*f+b**2)/&
          &(b**3+3*b*c**2+c**3+3*b**2*c)/(e**2-2&
          &*e*b+b**2)

     fd_dvdy(3,i)=(c*b*f-b*e*c-2*b*e*f+2*c*f*e)/f/b/e/c

     fd_dvdy(4,i)=-(-3*c*b*f-2*b*e*f+3*b*e*c+&
          &4*b*c**2-5*c**2*f-4*c*f*e+6*c**3+5*&
          &c**2*e)*b**2/c*e*f/(f**2-2*f*c+c**2)&
          &/(3*e**2*b*c**2+6*e*b**2*c**2+6*e*b*&
          &c**3+e**2*b**3+c**5+2*c**4*e+c**3*e*&
          &*2+3*b**2*c**3+b**3*c**2+3*b*c**4+2*&
          &e*b**3*c+3*e**2*b**2*c)

     fd_dvdy(5,i)=-sum(fd_dvdy(1:4,i))

     !Tridiagonal Coefficients:
     fd_dvdy(6,i)=-e*c**2*f/(f+b)/(b**2+2*b*c+c**2)/(b-e)  !this is alpha
     fd_dvdy(7,i)=1d0
     fd_dvdy(8,i)=-e*b**2*f/(-f+c)/(b**2+2*b*c+c**2)/(c+e) !this is beta
  enddo
  !$OMP END PARALLEL DO

<<<<<<< HEAD:coeft.f90
=======



>>>>>>> 2bls-new-stat:coeft.f90
  !BOUNDARY SCHEME FOR I=ny-1  ==============================
  a=h(ny-2); b=h(ny-1); c=h(ny); f=a+b+c; e=b+c
  fd_dvdy(1,ny-1) = c**2*(c-e)**2/f**2/(f-e)**2/(-f+c) 
  fd_dvdy(3,ny-1) =-(5*c**2-4*f*c-3*c*e+2*f*e)/c/(-f+c)/(c-e)
  fd_dvdy(4,ny-1) =-(c-e)**2*(-f+c)*(2*f*c+c*e+2*f*e)/c/f**2/e**3
  fd_dvdy(2,ny-1) = -(fd_dvdy(1,ny-1)+fd_dvdy(3,ny-1)+fd_dvdy(4,ny-1))
  !Tridg:
  fd_dvdy(6,ny-1) =-(-f+c)*c**2/e**2/(f-e)
  fd_dvdy(7,ny-1) =1d0
  fd_dvdy(8,ny-1) =-(c-e)**2*(-f+c)/f/e**2

  !BOUNDARY SCHEME FOR I=ny
  fd_dvdy(1,ny) = c**2*e/(-f+c)**2/(f-e)/f
  fd_dvdy(2,ny) = -c**2*f/e/(c-e)**2/(f-e)
  fd_dvdy(4,ny) =(f*c+c*e+2*f*e)/c/f/e
  fd_dvdy(3,ny) = -(fd_dvdy(1,ny)+fd_dvdy(2,ny)+fd_dvdy(4,ny))
  !Tridg:
  fd_dvdy(6,ny) =f*e/(-f+c)/(c-e)
  fd_dvdy(7,ny) =1d0

  !Reorganicing the coeffs to save flops in the trid-solver (LU)
  do j=2,ny
     fd_dvdy(6,j)=fd_dvdy(6,j)/fd_dvdy(7,j-1)
     fd_dvdy(7,j)=fd_dvdy(7,j)-fd_dvdy(6,j)*fd_dvdy(8,j-1)
  enddo
  fd_dvdy(7,:)=1d0/fd_dvdy(7,:)

  
  !!====================================================
  !!       D/DX using an uniform mesh (in place)
  !!====================================================

  fd_dx=0d0;
  !BOUNDARY SCHEME FOR I=1       !BOUNDARY SCHEME FOR I=2  
  fd_dx(1,1) =-37d0/12d0;        fd_dx(1,2) = -43d0/96d0;  
  fd_dx(2,1) = 2d0/3d0;          fd_dx(2,2) = -5d0/6d0;    
  fd_dx(3,1) = 3d0;              fd_dx(3,2) = 9d0/8d0;     
  fd_dx(4,1) = -2d0/3d0;         fd_dx(4,2) = 1d0/6d0;     
  fd_dx(5,1) = 1d0/12d0;         fd_dx(5,2) = -1d0/96d0;   
  !Tridg:                        !Tridg:                   
  fd_dx(7,1) =1d0;               fd_dx(6,2) =1d0/8d0;      
  fd_dx(8,1) =4d0;               fd_dx(7,2) =1d0;          
                                 fd_dx(8,2) =3d0/4d0;      
  !$OMP WORKSHARE
  fd_dx(1,3:nx-2)=-1d0/36d0;
  fd_dx(2,3:nx-2)=-7d0/9d0;
  fd_dx(3,3:nx-2)=0;
  fd_dx(4,3:nx-2)= 7d0/9d0;
  fd_dx(5,3:nx-2)= 1d0/36d0;
  !Tridg:
  fd_dx(6,3:nx-2)=1d0/3d0;
  fd_dx(7,3:nx-2)=1d0;
  fd_dx(8,3:nx-2)=1d0/3d0;
  !$OMP END WORKSHARE

  !BOUNDARY SCHEME FOR I=nx-1      !BOUNDARY SCHEME FOR I=nx    
  fd_dx(1,nx-1) =1d0/96d0;            fd_dx(1,nx) =-1d0/12d0        
  fd_dx(2,nx-1) =-1d0/6d0;            fd_dx(2,nx) =2d0/3d0          
  fd_dx(3,nx-1) =-9d0/8d0;            fd_dx(3,nx) =-3d0             
  fd_dx(4,nx-1) = 5d0/6d0;            fd_dx(4,nx) =-2d0/3d0         
  fd_dx(5,nx-1) = 43d0/96d0;          fd_dx(5,nx) =37d0/12d0        
  !Tridg:                          !Tridg:                      
  fd_dx(6,nx-1) =3d0/4d0;             fd_dx(6,nx) =4d0              
  fd_dx(7,nx-1) =1d0;                 fd_dx(7,nx) =1d0              
  fd_dx(8,nx-1) =1d0/8d0;
<<<<<<< HEAD:coeft.f90


  point0=ax*pi/(nx+1)
  fd_dx(1:5,:)=fd_dx(1:5,:)/point0

  !Reorganicing the coeffs to save flops in the trid-solver (LU)
  do j=2,nx
     fd_dx(6,j)=fd_dx(6,j)/fd_dx(7,j-1)
     fd_dx(7,j)=fd_dx(7,j)-fd_dx(6,j)*fd_dx(8,j-1)
  enddo
  fd_dx(7,:)=1d0/fd_dx(7,:)

=======


  point0=ax*pi/(nx+1)
  fd_dx(1:5,:)=fd_dx(1:5,:)/point0

  !Reorganicing the coeffs to save flops in the trid-solver (LU)
  do j=2,nx
     fd_dx(6,j)=fd_dx(6,j)/fd_dx(7,j-1)
     fd_dx(7,j)=fd_dx(7,j)-fd_dx(6,j)*fd_dx(8,j-1)
  enddo
  fd_dx(7,:)=1d0/fd_dx(7,:)



>>>>>>> 2bls-new-stat:coeft.f90

  !!====================================================
  !!    Interp X using an uniform mesh (mid-point)
  !!====================================================

  !BOUNDARY SCHEME FOR I=1  !BOUNDARY SCHEME FOR I=nx-1
  fd_ix(1,1) =5d0/24d0;     fd_ix(1,nx-1) =-1d0/84d0;
  fd_ix(2,1) =15d0/8d0;     fd_ix(2,nx-1) =1d0/4d0;
  fd_ix(3,1) =5d0/8d0;      fd_ix(3,nx-1) =5d0/4d0;
  fd_ix(4,1) =-1d0/24d0;    fd_ix(4,nx-1) =5d0/12d0;
  !Tridg:                   !Tridg:     
  fd_ix(7,1) =1d0;          fd_ix(6,nx-1) =5d0/6d0;
  fd_ix(8,1) =5d0/3d0;      fd_ix(7,nx-1) =1d0;
  fd_ix(8,nx-1) =1d0/14d0;   

  !$OMP WORKSHARE              !BOUNDARY SCHEME FOR I=nx 
  fd_ix(1,2:nx-2)=1d0/20d0;    fd_ix(1,nx) =1d0/8d0;     
  fd_ix(2,2:nx-2)=3d0/4d0;     fd_ix(2,nx) =-7d0/8d0;      
  fd_ix(3,2:nx-2)=3d0/4d0;     fd_ix(3,nx) =35d0/8d0;    
  fd_ix(4,2:nx-2)=1d0/20d0;    fd_ix(4,nx) =35d0/8d0;    
  !Tridg:                      !Tridg:                   
  fd_ix(6,2:nx-2)=3d0/10d0;    fd_ix(6,nx) =7d0;         
  fd_ix(7,2:nx-2)=1d0;         fd_ix(7,nx) =1d0;         
  fd_ix(8,2:nx-2)=3d0/10d0;
  !$OMP END WORKSHARE


  !Reorganicing the coeffs to save flops in the trid-solver (LU) 
  do j=2,nx
     fd_ix(6,j)=fd_ix(6,j)/fd_ix(7,j-1)
     fd_ix(7,j)=fd_ix(7,j)-fd_ix(6,j)*fd_ix(8,j-1)
  enddo
  fd_ix(7,:)=1d0/fd_ix(7,:)

<<<<<<< HEAD:coeft.f90
  
=======

>>>>>>> 2bls-new-stat:coeft.f90
  
  !!====================================================
  !!    Interp Y using non-uniform mesh (mid-point): u,w 
  !!           (intented for u,w & p positions only!)
  !!====================================================

  yp=(y(0:ny)+y(1:ny+1))*0.5d0    !@half cell position (u,w,p position)
  hi(1:ny+1)=y (1:ny+1)-y (0:ny); !dh for the interpolated grid
  h (2:ny+1)=yp(2:ny+1)-yp(1:ny); !dh for the original grid 
 
<<<<<<< HEAD:coeft.f90
  !BOUNDARY SCHEME FOR I=1 (for pressure)
=======
  !BOUNDARY SCHEME FOR j=1 (for pressure)
>>>>>>> 2bls-new-stat:coeft.f90
  j=1;a=hi(j+1)/2;b=a+h(j+2);c=b+h(j+3);d=c+h(j+4); 
  fd_iy(1,1)=-c*b*d/(-d+a)/(-c+a)/(-b+a);
  fd_iy(2,1)=c*a*d/(-d+b)/(-c+b)/(-b+a);
  fd_iy(3,1)=-b*a*d/(-d+c)/(c**2-a*c-b*c+b*a);
  fd_iy(4,1)=1d0-sum(fd_iy(1:3,1));

<<<<<<< HEAD:coeft.f90
  !BOUNDARY SCHEME FOR I=2
=======
  !BOUNDARY SCHEME FOR j=2
>>>>>>> 2bls-new-stat:coeft.f90
  j=2;a=hi(j+1)/2;b=a+h(j+2);c=b+h(j+3);d=c+h(j+4);e=hi(j)/2; 
  fd_iy(1,2)=2*c*b*a**2*d/(e+d)/(e+c)/(b+e)/(2*a**2+e**2+3*e*a);
  fd_iy(2,2)=-2*c*b*e*d/(a-d)/(a-c)/(a-b)/(a+e);
  fd_iy(3,2)=2*c*e*a**2*d/(-d+b)/(-c+b)/(2*a**2*e+2*a**2*b-3*a*b**2-3*a*e*b+b**3+e*b**2);
  fd_iy(4,2)=-2*a**2*b*e*d/(-d+c)/(3*e*a*c**2-2*e*a**2*c-3*b*e*a*c-2*a**2*c**2+2*a**2*e*b+3&
       &*a*c**3+b*c**3-e*c**3+2*b*a**2*c-3*b*a*c**2+b*e*c**2-c**4);
  !Tridg:
  fd_iy(7,2)=1d0;
  fd_iy(8,2)=-c*b*e*d/(2*a+e)/(2*a-d)/(2*a-c)/(2*a-b);
  fd_iy(5,2)=sum(fd_iy(7:8,2))-sum(fd_iy(1:4,2));   !because of consistency

 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,a,b,c,d,e,f)
  do j=3,ny-1
     a=hi(j+1)/2;b=a+h(j+2);e=hi(j)/2;f=e+h(j);
     fd_iy(1,j)=4*e**2*a**2*b/(b+f)/(2*a**2*f**2+f**4+4*e**2*a**2-6*a**2*e*f&
          &-9*e*a*f**2+6*e**2*a*f+2*e**2*f**2+3*a*f**3-3*e*f**3);
     fd_iy(2,j)=-4*f*a**2*b/(b+e)/(a+e)/(e-f)/(2*a+e);
     fd_iy(3,j)=-4*f*e**2*b/(a-b)/(a+f)/(2*e**2+a**2+3*e*a);
     !Tridg:   
     fd_iy(6,j)=-f*b*a**2/(a+e)/(-f+2*e)/(b+2*e)/(a+2*e);
     fd_iy(7,j)=1d0;
     fd_iy(8,j)=-f*e**2*b/(2*a-b)/(2*a+f)/(2*a**2+e**2+3*e*a);
     fd_iy(4,j)=sum(fd_iy(6:8,j))-sum(fd_iy(1:3,j)); !because of consistency
  enddo
  !$OMP END PARALLEL DO

<<<<<<< HEAD:coeft.f90
  !BOUNDARY SCHEME FOR I=ny
=======
  !BOUNDARY SCHEME FOR j=ny
>>>>>>> 2bls-new-stat:coeft.f90
  j=ny; a=hi(j+1)/2;e=hi(j)/2;f=e+h(j);g=f+h(j-1);    
  fd_iy(1,ny) =2*(4*e**2+2*e*a-a*f-2*f*e)*e**2/(-g+e)/(f-g)/g/(g+a);
  fd_iy(2,ny) =-2*(2*e*a+4*e**2-a*g-2*g*e)*e**2/(e-f)/f/(f**2+a*f-g*f-a*g);
  fd_iy(3,ny) =2*(8*e**3-4*g*e**2-2*f*a*e-2*g*e*a-4*f*e**2+g*a*f+4*a*e**2+2*f*g*e)&
       &/(e**3+a*e**2-f*e**2-g*e**2+g*a*f-g*e*a+f*g*e-f*a*e);
  fd_iy(4,ny) =2*(-2*f*e-2*g*e+g*f+4*e**2)*e**2/(g+a)/(a+f)/(a+e)/a;
  !Tridg:     
  fd_iy(6,ny) =1;
  fd_iy(7,ny) =(8*e**3-4*g*e**2-2*f*a*e-2*g*e*a-4*f*e**2+g*a*f+4*a*e**2+2*f*g*e)/g/a/f;
<<<<<<< HEAD:coeft.f90
  
=======

  !For pressure (it has only ny points the stencil)
  fd_iyp=fd_iy;
  !BOUNDARY SCHEME FOR j=ny-1
  j=ny-1;a=hi(j+1)/2;e=hi(j)/2;f=e+h(j);
  fd_iyp(1,ny-1) =4*e**2*a**2/(2*a+f)/(a+f)/(2*e**2+f**2-3*e*f);
  fd_iyp(2,ny-1) =-4*f*a**2/(3*e*a+e**2+2*a**2)/(e-f);
  fd_iyp(3,ny-1) =4*e**2*f/(a**3+2*e**2*a+3*e*a**2+2*e**2*f+f*a**2+3*f*e*a);
  fd_iyp(6,ny-1) =-f*a**2/(a+e)/(2*e-f)/(2*e+a);
  fd_iyp(7,ny-1) =1d0;
  fd_iyp(8,ny-1) =e**2*f/(4*a**3+2*e**2*a+6*e*a**2+3*f*e*a+e**2*f+2*f*a**2);
    
  !BOUNDARY SCHEME FOR j=ny
  j=ny; e=hi(j)/2;f=e+h(j);g=f+h(j-1);    
  fd_iyp(1,ny) =2*e**2*(2*e-f)/(e-g)/(f-g)/g;
  fd_iyp(2,ny) =-2*(-g+2*e)*e**2/(e-f)/f/(f-g);
  fd_iyp(3,ny) =2*(4*e**2+g*f-2*e*f-2*g*e)/(e**2-e*f+g*f-g*e);
  fd_iyp(6,ny) =1;
  fd_iyp(7,ny) =(4*e**2+g*f-2*e*f-2*g*e)/g/f;
  
  !Reorganicing the coeffs: ACHTUNG!! Start from j=2 
  !(j=1; for u,w=0d0; for p, taylor)
  do j=3,ny
     fd_iyp(6,j)=fd_iyp(6,j)/fd_iyp(7,j-1)
     fd_iyp(7,j)=fd_iyp(7,j)-fd_iyp(6,j)*fd_iyp(8,j-1)
  enddo
  fd_iyp(7,2:ny)=1d0/fd_iyp(7,2:ny)


>>>>>>> 2bls-new-stat:coeft.f90

  !Reorganicing the coeffs: ACHTUNG!! Start from j=2 
  !(j=1; for u,w=0d0; for p, taylor)
  do j=3,ny
     fd_iy(6,j)=fd_iy(6,j)/fd_iy(7,j-1)
     fd_iy(7,j)=fd_iy(7,j)-fd_iy(6,j)*fd_iy(8,j-1)
  enddo
  fd_iy(7,2:ny)=1d0/fd_iy(7,2:ny)

<<<<<<< HEAD:coeft.f90
=======


>>>>>>> 2bls-new-stat:coeft.f90
!=====================================================================================
!=====================================================================================
!=====================================================================================



 !********* Weights for Pressure Interpolation:
  !For pressure: y_int=y (pressure interpolated is at wall position)

  stencil=5
  dn=1d0;nm=1d0; !Initialize to 1

!Lagrange Interpolation:
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,xd,xl,k,i)
  !$OMP DO SCHEDULE(STATIC)
  do j=1,ny+1
     if (j.le.2) then
        xl(1:5)=yp(j:j+4);
     else if (j.ge.ny) then
        xl(1:5)=yp(j-4:j);
     else if (j.gt.2 .and. j.lt.ny) then
        xl(1:5)=yp(j-2:j+2);
     endif
     xd=y(j)
     do i=1,stencil
        do k=1,stencil
           if (k.ne.i) then
              nm(i,j)=nm(i,j)*(xd-xl(k))
              dn(i,j)=dn(i,j)*(xl(i)-xl(k))
           endif
        enddo
     enddo
     l_weight(:,j)=nm(:,j)/dn(:,j)
  enddo
  !$OMP END PARALLEL

  !SCHEME FOR 1ST DERIVATIVE:
  !Fix first ghost point so dp/dn=nu*d2V/dyy at the wall
  x1=yp(1);x2=yp(2);x3=yp(3);x4=yp(4);x5=yp(5);
  nm_aux=-x3*x4*x5-x2*x4*x5-x2*x3*x5-x2*x3*x4;  
  l_weight(1,0)=nm_aux/dn(1,1);
  nm_aux=-x3*x4*x5-x1*x4*x5-x1*x3*x5-x1*x3*x4;  
  l_weight(2,0)=nm_aux/dn(2,1);
  nm_aux=-x1*x4*x5-x2*x4*x5-x2*x1*x5-x2*x1*x4;  
  l_weight(3,0)=nm_aux/dn(3,1);
  nm_aux=-x3*x1*x5-x2*x1*x5-x2*x3*x5-x2*x3*x1;  
  l_weight(4,0)=nm_aux/dn(4,1);
  nm_aux=-x3*x4*x1-x2*x4*x1-x2*x3*x1-x2*x3*x4;
  l_weight(5,0)=nm_aux/dn(5,1);


!Weight for Second Derivative dV/dyy @th wall:
j=1;dn=1d0 !Initialize to 1
xl(1:stencil)=y(j:j+4)
do i=1,stencil
   do k=1,stencil
      if (k.ne.i) dn(i,j)=dn(i,j)*(xl(i)-xl(k)) 
   enddo
enddo
nm(1,1)=2*(xl(5)*(xl(4)+xl(2))+xl(3)*(xl(2)+xl(5))+xl(4)*(xl(2)+xl(3)));
nm(2,1)=2*(xl(5)*(xl(4)+xl(1))+xl(3)*(xl(1)+xl(5))+xl(4)*(xl(1)+xl(3)));
nm(3,1)=2*(xl(5)*(xl(4)+xl(2))+xl(1)*(xl(2)+xl(5))+xl(4)*(xl(2)+xl(1)));
nm(4,1)=2*(xl(5)*(xl(1)+xl(2))+xl(3)*(xl(2)+xl(5))+xl(1)*(xl(2)+xl(3)));
nm(5,1)=2*(xl(1)*(xl(4)+xl(2))+xl(3)*(xl(2)+xl(1))+xl(4)*(xl(2)+xl(3)));
ldyy(:)=nm(:,1)/dn(:,j)

!=====================================================================================
!=====================================================================================
!=====================================================================================





  !********* DISTANCE BETWEEN NODES AND FACES *************
  yyi(0:ny) = y(0:ny)+0.5d0*(y(1:ny+1)-y(0:ny))
  dy(0:ny) = y(1:ny+1)-y(0:ny)
  dyi(1:ny) = 0.5d0*(dy(1:ny)+dy(0:ny-1))

  point0=ax*pi/(nx+1)

  do i = 0,nx+1
     x(i) = dfloat(i)*point0
  enddo

  xxi(0:nx) = x(0:nx)+0.5d0*(x(1:nx+1)-x(0:nx))


  daux= 1000
  do i=0,nx
     dx(i) = x(i+1)-x(i)
     daux = min(daux,dx(i))
  enddo
  dxmin=daux
  do i=1,nx
     dxi(i) = 0.5d0*(dx(i)+dx(i-1))
  enddo

  idx  = 1d0/dx(1)
  idy  = 1d0/dy

  idxx = 1d0/dx(1)
  idyy(1:ny) = 2d0/(dy(0:ny-1)+dy(1:ny))
  !--------------------------------------------------------------


  ! this are the viscous timestep limitations, which we hope
  ! are higher than the cfl limitation


  daux = 1000
  do i =1,ny
     dymin(i) =dy(i)*cfl/dsqrt(3d0)
     daux = min(daux,dy(i))
  enddo

  dzmin = az*2*pi/nz

  dtrex = re*(cfl*dxmin**2)/(6.0d0)
  dtrez = re*(cfl*dzmin**2)/(pi**2)
  dtrey = re*(cfl*daux**2)/(6.0d0)
  dtret = min(dtrex,dtrez)

  ! the sqrt(3) is the maximum wavenumber of the fourth order
  ! compact scheme.
  ! and pi is for the fourier spectral method in z.


  dymin(1:ny) = dy(1:ny)*cfl/dsqrt(3d0)


  dxmin = dxmin*cfl/dsqrt(3d0)
  dzmin = dzmin*cfl/pi

  invre    = 1d0/re
  signo(1) = -invre
  signo(2) = 1d0
  signo(3) = 2d0

  fact=550d0
  alp=0.97d0
  bet=0.05d0

  do i=1,nx
     fren(i)= dtanh(fact*(1d0+0.5d0*(x(i)+x(i-1))-0.5d0*(x(1)+x(0))))
  enddo

  do i=1,nx
     fren(i)= 0.5d0*fren(i)*(1+dtanh(fact*(alp*x(nx) &
          &        -0.5d0*(x(i)+x(i-1))-0.5d0*(x(1)+x(0)))))

  enddo

  do i=1,nx
     fren(i)=fren(i)+bet*(1+dtanh(-fact*(alp*x(nx)  &
          &        -0.5d0*(x(i)+x(i-1))-0.5d0*(x(1)+x(0)))))
  enddo

  rex=-1d0/re

  ! coefficient for viscosity en z***************

  do i=0,nz2
     kvis(i)=rex* kaz2(i)
  enddo

  !      write(*,*) 'kvis'
  !*********************************************
  !***************************************************
  ! Runge-Kutta coefficients**************************


  rkc(1) = 8d0/15d0
  rkc(2) = 5d0/12d0
  rkc(3) = 3d0/4d0 

  rkd(1) = 0d0
  rkd(2) = -17d0/60d0
  rkd(3) = -5d0/12d0

  rkcv(1) = 4d0/15d0
  rkcv(2) = 1d0/15d0
  rkcv(3) = 1d0/6d0

  rkdv(1) = 4d0/15d0
  rkdv(2) = 1d0/15d0
  rkdv(3) = 1d0/6d0

  ! coefficientes para interpolar ********************************

  ! parte implicito **********************************************
  ! u,x
  do i = 2,nx-1
     c  = dxi(i+1)
     a  = dxi(i)
     b  = dx(i)
     inxu(1,i) =  c*b**2/((2*a+b)*(-b+2*a)*(c+a))
     inxu(2,i) = 1d0
     inxu(3,i) = -a*b**2/((b-2*c)*(b+2*c)*(c+a))
  enddo
  ! u,y 
  do i = 2,ny-1
     a=dy(i)
     b=dy(i-1)
     inyu(1,i) = a**2/((a+b)*(a+2*b))
     inyu(2,i) = 1d0
     inyu(3,i) = b**2/(b**2+2*a**2+3*a*b)
  enddo
  ! v,x
  do i = 2,nx-1
     a = dx(i)
     b = dx(i-1)
     inxv(1,i)=  a**2/((b+a)*(a+2d0*b)) 
     inxv(2,i) = 1d0
     inxv(3,i) = b**2/(2d0*a**2+b**2+3d0*a*b)
  enddo
  ! v,y
  do i = 2,ny-1
     a=dyi(i)
     b=dy(i-1)
     c=dyi(i-1)
     inyv(1,i) = c*b**2/((2*a+b)*(-b+2*a)*(c+a))
     inyv(2,i) = 1d0
     inyv(3,i) = -a*b**2/((b-2*c)*(b+2*c)*(c+a))
  enddo



  inxu(3,1)  = 0d0
  inxu(2,1)  = 1d0 
  inxu(1,1)  = 0d0

  inxu(1,nx) = 0d0 
  inxu(2,nx) = 1d0 
  inxu(3,nx) = 0d0 

  inxv(1,1)  = 0d0
  inxv(2,1)  = 1d0
  inxv(3,1)  = 0d0

  inxv(1,nx) = 0d0
  inxv(2,nx) = 1d0
  inxv(3,nx) = 0d0 

  inyu(1,1)  = 0d0
  inyu(2,1)  = 1d0
  inyu(3,1)  = 0d0

  inyu(1,ny) = 0d0
  inyu(2,ny) = 1d0
  inyu(3,ny) = 0d0

  inyv(1,1)  = 0d0
  inyv(2,1)  = 1d0
  inyv(3,1)  = 0d0

  inyv(1,ny) = 0d0
  inyv(2,ny) = 1d0
  inyv(3,ny) = 0d0
  ! los coefficientes de los condiciones de contorno en x para u,v ***
  ! primero v despues u
  n = 3
  flagv = 0
  a = dx(nx-1)
  b = dxi(nx-1)
  c = dxi(nx-2)
  d = dxi(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! la velocidad v en la salida
  inbx(1,1)  = bcc(1)
  inbx(1,2)  = bcc(2)
  inbx(1,3)  = bcc(3)
  inbx(1,4)  = 0d0
  a = dx(nx)
  b = dx(nx-1)
  c = dx(nx-2)
  d = dx(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv)  ! la velocidad u en la salida
  inbx(2,1)  = bcc(1)
  inbx(2,2)  = bcc(2)
  inbx(2,3)  = bcc(3)
  inbx(2,4)  = 0d0
  n = 4
  flagv = 4
  a = dx(1)
  b = dx(1)
  c = dx(2)
  d = dx(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)  ! la velocidad u en la entrada
  inbx(3,1)  = bcc(1)
  inbx(3,2)  = bcc(2)
  inbx(3,3)  = bcc(3)
  inbx(3,4)  = bcc(4)
  n = 4
  flagv = 4
  a = dx(0)
  b = dx(1)
  c = dxi(2)
  d = dxi(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)  ! la velocidad v en la entrada
  inbx(4,1)  = bcc(1)
  inbx(4,2)  = bcc(2)
  inbx(4,3)  = bcc(3)
  inbx(4,4)  = bcc(4)

  !      write(*,*) (inbx(4,i),i=1,4)
  ! end boundary coefficients in x **********************
  ! los coefficientes de los condiciones de contorno en y para u,v ***
  ! primero v despues u

  n = 4
  flagv = 4
  a = dy(ny)
  b = dy(ny-1)
  c = dyi(ny-1)
  d = dyi(ny-2)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! u en n=ny
  inby(1,1) = bcc(1)
  inby(1,2) = bcc(2)
  inby(1,3) = bcc(3)
  inby(1,4) = bcc(4)

  n = 3
  flagv = 0
  a = -dy(0)
  b = -dy(1)
  c = -dy(2)
  d =  dy(4)  ! Atencion. Esto con el termino en res es dy=-1
  call coefb(n,np,a,b,c,d,bcc,flagv) ! v en n=1
  inby(3,1) = bcc(1)
  inby(3,2) = bcc(2)
  inby(3,3) = bcc(3)
  inby(3,4) = 0d0

  n = 4
  flagv = 4
  a = dy(0)
  b = dy(1)
  c = dyi(2)
  d = dyi(3)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! u en n=1 para el termino viscoso
  inby(2,1)   = bcc(1)
  inby(2,2)   = bcc(2)
  inby(2,3)   = bcc(3)
  inby(2,4)   = bcc(4)

  a = -dy(ny-1)
  b = -dy(ny-1)
  c = -dy(ny-2)
  d = -dy(ny-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! v en n=ny para el termino viscoso
  inby(4,1)   = bcc(1)
  inby(4,2)   = bcc(2)
  inby(4,3)   = bcc(3)
  inby(4,4)   = bcc(4)

  ! end boundary coefficients in y **********************

  ! calculate inverse of coefficients
  call coefinv(inxu,nx,1)
  call coefinv(inxv,nx,1)
  call coefinv(inyu,ny,1)
  call coefinv(inyv,ny,1)
  ! parte explicito x,y *******************************************
  do i = 2,ny-1
     a = dyi(i)
     b = dy(i-1)
     c = dyi(i-1)

     cofivy(2,i) = (-2d0*c*a)/(2*a*b+b**2-4*c*a-2*c*b)
     cofivy(1,i) = (2d0*c*a)/(-b**2+2*a*b-2*c*b+4*c*a)

     a = dy(i)
     b = dy(i-1)

     cofiuy(2,i) = 4*b**2/(2*b**2+a**2+3*a*b)
     cofiuy(1,i) = 4*a**2/((a+b)*(2*a+b))
  enddo

  do i = 2,nx-1
     ! v,x
     a = dx(i)
     b = dx(i-1)

     cofivx(2,i) = (4d0*b**2)/(2d0*b**2+a**2d0+3d0*a*b)
     cofivx(1,i) =  (4d0*a**2)/((b+a)*(2d0*a+b))
     ! u,x
     c  = dxi(i+1)
     a  = dxi(i)
     b  = dx(i)
     cofiux(2,i) =    (-2d0*c*a)/(2*a*b+b**2-4*c*a-2*c*b)
     cofiux(1,i) =  ( 2d0*c*a)/(-b**2+2*a*b-2*c*b+4*c*a)
  enddo

  !***************************************************************
  ! condiciones de contorno **************************************
  ! u,x
  ! v,x
  ! coefficientes para terminino convectivos *********************
  ! parte implicito **********************************************
  ! u,x
  do i = 2,nx-1
     a=dx(i)
     b=dx(i-1)
     dcxu(1,i) =  -a*(-2d0*a*b-b**2+2d0*a**2)/   &
          &                  (2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)
     dcxu(2,i) = 1d0
     dcxu(3,i) =  ((a**2+2d0*a*b-2d0*b**2)*b)/    &
          &                (2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)
  enddo
  dcxu(2,1)   = 1d0
  dcxu(1,1)   = 0d0
  dcxu(3,1)   = 0d0
  dcxu(2,nx)  = 1d0
  dcxu(1,nx)  = 0d0
  dcxu(3,nx)  = 0d0

  ! 0d0
  ! u,y 
  do i = 2,ny
     a = dyi(i)
     b = dyi(i-1)
     c = dy(i-1)
     dcyu(1,i) = (b*c**2)/(-c**2*a-c**2*b+12d0*b**2*a+12d0*b*a**2)
     dcyu(2,i) = 1d0
     dcyu(3,i) = (c**2*a)/(-c**2*a-c**2*b+12d0*b**2*a+12d0*b*a**2)
  enddo
  dcyu(2,1)  = 1d0
  dcyu(3,1)  = 0d0
  dcyu(1,1)  = 0d0
  dcyu(2,ny+1) = 1d0
  dcyu(3,ny+1) = 0d0
  dcyu(1,ny+1) = 0d0
  ! v,x
  do i = 2,nx-1
     a = dxi(i)
     b = dxi(i-1)
     c = dx(i-1)
     dcxv(1,i) = (b*c**2)/(-c**2*a-c**2*b+12d0*b**2*a+12d0*b*a**2)
     dcxv(2,i) = 1d0
     dcxv(3,i) = (c**2*a)/(-c**2*a-c**2*b+12d0*b**2*a+12d0*b*a**2)
  enddo
  dcxv(2,1)  = 1d0
  dcxv(3,1)  = 0d0
  dcxv(1,1)  = 0d0
  dcxv(2,nx) = 1d0
  dcxv(3,nx) = 0d0
  dcxv(1,nx) = 0d0

  ! v,y
  do i = 2,ny-1
     a=dy(i)
     b=dy(i-1)
     aux= (2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)    
     dcyv(1,i) = -a*(-2d0*a*b-b**2+2d0*a**2)/aux
     dcyv(2,i) = 1d0
     dcyv(3,i) = ((a**2+2d0*a*b-2d0*b**2)*b)/aux
  enddo
  dcyv(2,1)  = 1d0
  dcyv(3,1)  = 0d0
  dcyv(1,1)  = 0d0
  dcyv(2,ny) = 1d0
  dcyv(3,ny) = 0d0
  dcyv(1,ny) = 0d0
  ! los coefficnxntes de los condiciones de contorno en x para u,v ***
  ! primero u y'uv' en N=1 despues para N=nx
  n = 4
  flagv = 3

  a = dx(1)
  b = dxi(2)
  c = dxi(3)
  d = dxi(4)

  call coefb(n,np,a,b,c,d,bcc,flagv) ! la velocidad u en la entrada 
  dcbx(1,1) = bcc(1)
  dcbx(1,2) = bcc(2)
  dcbx(1,3) = bcc(3)
  dcbx(1,4) = bcc(4)

  a = dx(0)
  b = dx(1)
  c = dx(2)
  d = dx(3)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! la velocidad uv en la entrada
  dcbx(2,1) = bcc(1)
  dcbx(2,2) = bcc(2)
  dcbx(2,3) = bcc(3)
  dcbx(2,4) = bcc(4)

  flagv = 5

  a = dx(nx)
  b = dx(nx-1)
  c = dxi(nx-1)
  d = dxi(nx-2)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! la velocidad u en la salida 
  dcbx(3,1)  = bcc(1)
  dcbx(3,2)  = bcc(2)
  dcbx(3,3)  = bcc(3)
  dcbx(3,4)  = bcc(4)

  a = dx(nx-1)
  b = dx(nx-1)
  c = dx(nx-2)
  d = dx(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv)  ! la velocidad uv en la salida
  dcbx(4,1)  = bcc(1)
  dcbx(4,2)  = bcc(2)
  dcbx(4,3)  = bcc(3)
  dcbx(4,4)  = bcc(4)

  a = -dx(nx-1)
  b = -dxi(nx-1)
  c = -dxi(nx-2)
  d = -dxi(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) ! la velocidad u en la entrada 
  dcbx(5,1) = bcc(1)
  dcbx(5,2) = bcc(2)
  dcbx(5,3) = bcc(3)
  dcbx(5,4) = bcc(4)
  ! end boundary coefficients in x **********************

  ! los coefficientes de los condiciones de contorno en y para u,v ***
  ! primero v despues u
  n = 4
  flagv = 5
  a = -dy(0)
  b = -dy(1)
  c = -dyi(2)
  d = -dyi(3)
  call coefb(n,np,a,b,c,d,bcc,flagv) !la velocidad v en y=0
  dcby(1,1)  = bcc(1)
  dcby(1,2)  = bcc(2)
  dcby(1,3)  = bcc(3)
  dcby(1,4)  = bcc(4)

  flagv = 3
  a = -dy(ny-1)
  b = -dyi(ny-1)
  c = -dyi(ny-2)
  d = -dyi(ny-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) !la velocidad v en y=ngy
  dcby(2,1)   = bcc(1)
  dcby(2,2)   = bcc(2)
  dcby(2,3)   = bcc(3)
  dcby(2,4)   = bcc(4)

  a = dy(0)
  b = dy(1)
  c = dy(2)
  d = dy(3)
  call coefb(n,np,a,b,c,d,bcc,flagv) !la velocidad uv en y=1
  dcby(3,1)   = bcc(1)
  dcby(3,2)   = bcc(2)
  dcby(3,3)   = bcc(3)
  dcby(3,4)   = bcc(4)

  a = -dy(ny)
  b = -dy(ny-1)
  c = -dy(ny-2)
  d = -dy(ny-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) !la velocidad uv en y=ngy
  dcby(4,1)   = bcc(1)
  dcby(4,2)   = bcc(2)
  dcby(4,3)   = bcc(3)
  dcby(4,4)   = bcc(4)

  flagv = 3
  a = -dy(ny-1)
  b = -dyi(ny-1)
  c = -dyi(ny-2)
  d = -dyi(ny-3)
  call coefb(n,np,a,b,c,d,bcc,flagv) !la velocidad v en y=ngy
  dcby(2,1)   = bcc(1)
  dcby(2,2)   = bcc(2)
  dcby(2,3)   = bcc(3)
  dcby(2,4)   = bcc(4)
  !      write(*,*) (dcby(2,i),i=1,4)

  n=4
  flagv = 7
  a = dy(ny)
  b = dy(ny-1)
  c = dyi(ny-1)
  d = dyi(ny-2)
  call coefb(n,np,a,b,c,d,bcc,flagv) !zero vorticity b.c.
  ccon(1,1)   = bcc(1)
  ccon(1,2)   = bcc(2)
  ccon(1,3)   = bcc(3)
  ccon(1,4)   = bcc(4)
  !      write(*,*) (ccon(1,i),i=1,4)

  flagv = 6
  a = -dy(0)
  b = -dyi(1)
  c = -dyi(2)
  d = 0d0
  call coefb(n,np,a,b,c,d,bcc,flagv) !zero vorticity b.c.
  ccon(2,1)   = bcc(1)
  ccon(2,2)   = bcc(2)
  ccon(2,3)   = bcc(3)
  ccon(2,4)   = bcc(4)
  !      write(*,*) (ccon(2,i),i=1,4)
  ! end boundary coefficients in y **********************

  ! calculate inverse of coefficients
  call coefinv(dcxu,nx  ,1)
  call coefinv(dcxv,nx  ,1)
  call coefinv(dcyu,ny+1,1)
  call coefinv(dcyv,ny  ,1)
  ! parte explicito **********************************************
  do i=2,ny-1
     a=dy(i)
     b=dy(i-1)
     cofcyv(1,i) =-24*a*b/(2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)
     cofcyv(2,i) = 24*a*b/(2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)

     a = dyi(i)
     b = dyi(i-1)
     c = dy(i-1)
     cofcyu(1,i) =-12d0*a*b/((12d0*a*b-c**2)*c)
     cofcyu(2,i) = 12d0*a*b/((12d0*a*b-c**2)*c)
  enddo
  cofcyu(1,i) = 0d0
  cofcyu(2,i) = 0d0
  i = ny
  a = dyi(i)
  b = dyi(i-1)
  c = dy(i-1)
  cofcyu(1,i) =-12d0*a*b/((12d0*a*b-c**2)*c)
  cofcyu(2,i) = 12d0*a*b/((12d0*a*b-c**2)*c)

  do i=2,nx-1
     ! v,x
     a = dxi(i)
     b = dxi(i-1)
     c = dx(i-1)
     cofcxv(1,i) = -12d0*a*b/((12d0*a*b-c**2)*c)
     cofcxv(2,i) =  12d0*a*b/((12d0*a*b-c**2)*c)
     a=dx(i)
     b=dx(i-1)
     cofcxu(1,i) = -24*a*b/(2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)
     cofcxu(2,i) =  24*a*b/(2d0*b**3+2d0*a**3+9d0*b**2*a+9d0*b*a**2)
  enddo
  ! end parte convective terms ******************************************
  ! coefficientes para terminino viscoso ********************************

  ! parte implicito *****************************************************
  ! u,x
  do i = 2,nx-1
     a = dx(i)
     b = dx(i-1)
     vixu(1,i) = b*(a**2+a*b-b**2)/(b**3+a**3+4*a**2*b+4*a*b**2) 
     vixu(2,i) = 1d0
     vixu(3,i) = -(a**2-a*b-b**2)*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo
  vixu(1,1)   = 0d0
  vixu(2,1)   = 1d0
  vixu(3,nx) = 0d0
  vixu(2,nx) = 1d0

  ! u,y 
  do i = 2,ny
     a = dyi(i)
     b = dyi(i-1)
     viyu(1,i) = b*(a**2+a*b-b**2)/(b**3+a**3+4*a**2*b+4*a*b**2) 
     viyu(2,i) = 1d0
     viyu(3,i) = -(a**2-a*b-b**2)*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo
  viyu(1,1)    = 0d0
  viyu(2,1)    = 1d0
  viyu(3,ny+1) = 0d0
  viyu(2,ny+1) = 1d0
  ! v,x
  do i = 2,nx-1
     a = dxi(i)
     b = dxi(i-1)
     vixv(1,i) = b*(a**2+a*b-b**2)/(b**3+a**3+4*a**2*b+4*a*b**2) 
     vixv(2,i) = 1d0
     vixv(3,i) = -(a**2-a*b-b**2)*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo
  vixv(1,1)  = 0d0
  vixv(2,1)  = 1d0
  vixv(3,nx) = 0d0
  vixv(2,nx) = 1d0

  ! v,y 
  do i = 2,ny-1
     a=dy(i)
     b=dy(i-1)
     viyv(1,i) = b*(a**2+a*b-b**2)/(b**3+a**3+4*a**2*b+4*a*b**2) 
     viyv(2,i) = 1d0
     viyv(3,i) = -(a**2-a*b-b**2)*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo
  viyv(1,1)  = 0d0
  viyv(2,1)  = 1d0
  viyv(3,ny) = 0d0
  viyv(2,ny) = 1d0
  ! los coefficientes de los condiciones de contorno en x para u,v ********
  ! primero v despues u

  flagv = 1
  n=5
  a = dxi(1)
  b = dxi(2)
  c = dxi(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)
  vixv(3,1)   = bcc(1)
  cofvbx(1,1) = bcc(2)
  cofvbx(1,2) = bcc(3)
  cofvbx(1,3) = bcc(4)
  cofvbx(1,4) = bcc(5)

  a = dxi(nx-1)
  b = dxi(nx-2)
  c = dxi(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv)
  vixv(1,nx) = bcc(1)
  cofvbx(2,1) = bcc(2)
  cofvbx(2,2) = bcc(3)
  cofvbx(2,3) = bcc(4)  
  cofvbx(2,4) = bcc(5)

  n=5
  a = dx(1)
  b = dx(2)
  c = dx(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)

  vixu(3,1)   = bcc(1)
  cofvbx(3,1) = bcc(2)
  cofvbx(3,2) = bcc(3)
  cofvbx(3,3) = bcc(4)
  cofvbx(3,4) = bcc(5)

  a = dx(nx-1)
  b = dx(nx-2)
  c = dx(nx-3)
  call coefb(n,np,a,b,c,d,bcc,flagv)

  vixu(1,nx) = bcc(1)
  cofvbx(4,1) = bcc(2)
  cofvbx(4,2) = bcc(3)
  cofvbx(4,3) = bcc(4)  
  cofvbx(4,4) = bcc(5)
  ! end coefficients *********************************************************
  !     los coefficientes de los condiciones de contorno en y para u,v ******** 
  n=5
  a = dy(1)
  b = dy(2)
  c = dy(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)

  viyv(3,1)   = bcc(1)
  cofvby(1,1) = bcc(2)
  cofvby(1,2) = bcc(3)
  cofvby(1,3) = bcc(4)
  cofvby(1,4) = bcc(5)

  a = dy(ny-1)
  b = dy(ny-2)
  c = dy(ny-3)

  call coefb(n,np,a,b,c,d,bcc,flagv)

  viyv(1,ny) = bcc(1)
  cofvby(2,1) = bcc(2)
  cofvby(2,2) = bcc(3)
  cofvby(2,3) = bcc(4)  
  cofvby(2,4) = bcc(5)

  n=5
  a = dyi(1)
  b = dyi(2)
  c = dyi(3)
  call coefb(n,np,a,b,c,d,bcc,flagv)
  viyu(3,1)   = bcc(1)
  cofvby(3,1) = bcc(2)
  cofvby(3,2) = bcc(3)
  cofvby(3,3) = bcc(4)
  cofvby(3,4) = bcc(5)

  a = dyi(ny)
  b = dyi(ny-1)
  c = dyi(ny-2)
  call coefb(n,np,a,b,c,d,bcc,flagv)

  viyu(1,ny+1) = bcc(1)
  cofvby(4,1) = bcc(2)
  cofvby(4,2) = bcc(3)
  cofvby(4,3) = bcc(4)  
  cofvby(4,4) = bcc(5)
  ! end coefficients *********************************************************
  ! coefficients for implicit viscous terms en y
  do i = 1,ny
     vyui(1,i) = viyu(1,i)
     vyvi(1,i) = viyv(1,i)
     vyui(2,i) = viyu(2,i)
     vyvi(2,i) = viyv(2,i)
     vyui(3,i) = viyu(3,i)
     vyvi(3,i) = viyv(3,i)
  enddo
  i=ny+1
  vyui(1,i) = viyu(1,i)
  vyui(2,i) = viyu(2,i)
  vyui(3,i) = viyu(3,i)
  ! end implicit viscous terms ***********************************************
  ! calculate inverse of coefficients

  call coefinv(vixu,nx  ,1)
  call coefinv(vixv,nx  ,1)
  call coefinv(viyu,ny+1,1)
  call coefinv(viyv,ny  ,1)
  ! parte explicito **********************************************
  do i=2,ny-1
     a=dy(i)
     b=dy(i-1)

     cofvyv(1,i) = 12*b/(b**3+a**3+4*a**2*b+4*a*b**2)
     cofvyv(2,i) =-12/(a**2+3*a*b+b**2)
     cofvyv(3,i) = 12*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo

  do i=2,ny
     a=dyi(i)
     b=dyi(i-1)

     cofvyu(1,i) = 12*b/(b**3+a**3+4*a**2*b+4*a*b**2)
     cofvyu(2,i) =-12/(a**2+3*a*b+b**2)
     cofvyu(3,i) = 12*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo

  do i=2,nx-1
     a=dxi(i)
     b=dxi(i-1)

     cofvxv(1,i) = 12*b/(b**3+a**3+4*a**2*b+4*a*b**2)
     cofvxv(2,i) =-12/(a**2+3*a*b+b**2)
     cofvxv(3,i) = 12*a/(b**3+a**3+4*a**2*b+4*a*b**2)

     a=dx(i)
     b=dx(i-1)

     cofvxu(1,i) = 12*b/(b**3+a**3+4*a**2*b+4*a*b**2)
     cofvxu(2,i) =-12/(a**2+3*a*b+b**2)
     cofvxu(3,i) = 12*a/(b**3+a**3+4*a**2*b+4*a*b**2)
  enddo
  ! condiciones de contorno ********************************************
  ! u,x
  ! u,y   
  ! v,x
  ! v,y   
  ! end gradient coefficients at boundaries
  ! end gradient coefficients ******************************************




  ! constantes para la multimalla **************************************
  do i=1,nx-1
     kmod(i) =2d0/dx(1)**2*(1-cos((pi*(i-1))/dfloat(nx-1)))
  enddo


  call mtr(ayp,cyclx,cycly,y)

  ! end constantes multimalla ********************************************




  do i = 1,nx-1
     axphi(1,i) =  1/dx(1)**2
     axphi(2,i) = -2/dx(1)**2
     axphi(3,i) =  1/dx(1)**2
  enddo
  axphi(1,1)     = 0d0
  axphi(3,nx-1)  = 0d0
  axphi(2,1)     = -1/dx(1)**2

  do i = 1,nx-1
     phix(i) = 1d0
  enddo

  iflag = 0

  call soltrix(phix,nx-1,1,nx-1,phix,1,1,axphi,iflag)

  do j = 1,ny-1
     phiy(j) = 1d0
  enddo

  do j = 1,ny-1
     ayphi(1,j) =  ayp(1,j)
     ayphi(2,j) =  ayp(2,j)
     ayphi(3,j) =  ayp(3,j)
  enddo
  ayphi(2,ny-1)   = 2d0*ayphi(2,ny-1)
  !***********************************************************************
  !  c.c.:  phi(i,2) = phi(i,2) y
  ! phi(i,ny+1) = phi(i,ny+1)-a3(ny+1)*phi(i,ny+2) pero phi(i,ny+2)=0
  !************************************************************************
  iflag = 0
  call soltrix(phiy,ny-1,1,ny-1,phiy,1,1,ayphi,iflag)

  !**********************************************************
  ! Linear interpolation to connect the two boundary layers
  !**********************************************************

 if (mpiid == mpi_inlet) then
    call MPI_SEND(y,ny+2,MPI_REAL8,mpiid_2(0),1,MPI_COMM_WORLD,istat,ierr)
 end if

<<<<<<< HEAD:coeft.f90
!  if (mpiid == mpi_inlet) then
!     call MPI_SEND(y,ny+2,MPI_REAL8,mpiid_2(0),1,MPI_COMM_WORLD,istat,ierr)
!  end if


!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
=======

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
>>>>>>> 2bls-new-stat:coeft.f90
  return 

end subroutine coef


!***************************************************************
!***************************************************************
!************ Auxiliary subroutines ****************************
!***************************************************************
!***************************************************************
!***************************************************************

subroutine mtr(ayp,cyclx,cycly,yy)
  use ctesp

  implicit none 


  integer i,iiiy,prim(4),prmy,prmx

  integer ngy2(mcycl),cyclx,cycly,lt
  real*8 bc2,b2,ac2,a3,bc1,ac1,dy1,dy2,cyclx1,cycly1            
  real*8 yy(0:ny+1),ayp(3,ny)


  ngy2(1) = ny-1

  prim(1) = 1
  prim(2) = 3
  prim(3) = 5
  prim(4) = 7

  cyclx1   = 10.5d0
  i = 1

  do while ((cyclx1-int(cyclx1)) .gt. 1d-5 .and. i .le. 4)
     cyclx1=log(float(nx-2)/prim(i))/log(dfloat(2))
     i=i+1
  enddo
  prmx = prim(i-1)

  cycly1   = 10.5d0
  i = 1
  do while ((cycly1-int(cycly1)) .gt. 1d-5 .and. i .le. 4)
     cycly1=log(float(ny-2)/prim(i))/log(dfloat(2))
     i=i+1
  enddo
  prmy = prim(i-1)

  cyclx = int(cyclx1)
  cycly = int(cycly1)

  do i = 2,cycly
     ngy2(i) = prmy*2**(cycly+1-i)+1
  enddo

  do i = cycly+1,cyclx
     ngy2(i) = prmy*2**(cycly+1-cycly)+1
  enddo


  i=1
  ac2 = 0d0
  bc2 = 0d0

  ac1 = ac2
  bc1 = bc2
  dy1 = yy(1+1)-yy(1)
  dy2 = yy(1+2)-yy(1+1)
  a3  =  1d0/dy1
  b2  = -1d0/dy1
  ac2 =  2d0/(dy1+dy2)
  bc2 = -2d0/(dy1+dy2)
  ayp(1,1) = b2*bc1
  ayp(3,1) = a3*ac2
  ayp(2,1) = a3*bc2+b2*ac1

  do iiiy = 2,ngy2(i)-1
     ac1 = ac2
     bc1 = bc2
     dy1 = yy(iiiy+1)-yy(iiiy)
     dy2 = yy(iiiy+2)-yy(iiiy+1)
     a3  =  1d0/dy1
     b2  = -1d0/dy1
     ac2 =  2d0/(dy1+dy2)
     bc2 = -2d0/(dy1+dy2)
     ayp(1,iiiy) = b2*bc1
     ayp(3,iiiy) = a3*ac2
     ayp(2,iiiy) = a3*bc2+b2*ac1

  enddo
  lt=ngy2(i)
  b2  = -1d0/dy2
  ac1  =  ac2
  bc1  =  bc2
  ayp(1,lt) = b2*bc1
  ayp(3,lt) = 0d0
  ayp(2,lt) = b2*ac1


end subroutine mtr



! ------------------------ Coefinv -------------------------------!

subroutine coefinv(a,neq,npter1)
  implicit none
  integer j,neq,npter1
  real*8 a(3,npter1:neq)
  real*8 one,d
  data one/1.0d0/

  do j = npter1+1,neq
     d       = one/a(2,j-1)
     a(2,j-1) =  d
     a(1,j  ) = -a(1,j)*d
     a(2,j  ) =  a(2,j)+a(1,j)*a(3,j-1)
  enddo
  a(2,neq) = one/a(2,neq)

end subroutine coefinv

! ****************************************************
! calculate coefficients for the boundary conditions 
! using lu decomposition
! ****************************************************

subroutine coefb(n,np,a,b,c,d,bcc,flagv)
  implicit none

  real*8 bcc(5),bc(5,5)
  real*8 indx(5)
  real*8 a,b,c,d

  integer n,np,i,j,flagv,dd

  if (flagv .eq. 1) then
     dd = 0
     do i=1,n
        bcc(i) = 0d0
        indx(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1 ) =  1d0

     bc(1,1) = -1d0
     bc(1,2) =  0d0
     bc(1,3) =  1/2d0*a**2
     bc(1,4) =  1/2d0*(a+b)**2
     bc(1,5) =  1/2d0*(a+b+c)**2

     bc(2,1) =  0d0
     bc(2,2) =  1d0
     bc(2,3) =  1d0
     bc(2,4) =  1d0
     bc(2,5) =  1d0

     bc(3,1) =  0d0
     bc(3,2) =  0d0
     bc(3,3) =  a
     bc(3,4) =  a+b
     bc(3,5) =  a+b+c

     bc(4,1) = -a
     bc(4,2) =  0d0
     bc(4,3) =  1/6d0*a**3
     bc(4,4) =  1/6d0*(a+b)**3
     bc(4,5) =  1/6d0*(a+b+c)**3

     bc(5,1) = -1/2d0*a**2
     bc(5,2) =  0d0 
     bc(5,3) =  1/24d0*a**4
     bc(5,4) =  1/24d0*(a+b)**4
     bc(5,5) =  1/24d0*(a+b+c)**4

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 0) then
     dd = 0
     do i=1,n
        bcc(i) = 0d0
        indx(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0

     bc(2,1) =  -1/2d0*a
     bc(2,2) =  -(1/2d0*a+b)
     bc(2,3) =  -(1/2d0*a+b+c)

     bc(3,1) =  1/8d0*a**2
     bc(3,2) =  1/2d0*(-(1/2d0*a+b))**2
     bc(3,3) =  1/2d0*(-(1/2d0*a+b+c))**2

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 2) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1)  =  1d0

     bc(1,1) = -1d0
     bc(1,2) =  1/2d0*a
     bc(1,3) =  (a+1/2d0*b) 
     bc(1,4) =  (a+b+1/2d0*c)
     bc(1,5) =  (a+b+c+1/2d0*d)

     bc(2,1) =  0d0
     bc(2,2) =  1d0
     bc(2,3) =  1d0
     bc(2,4) =  1d0
     bc(2,5) =  1d0

     bc(3,1) = -a
     bc(3,2) =  1/8d0*a**2
     bc(3,3) =  1/2d0*(a+1/2d0*b)**2
     bc(3,4) =  1/2d0*(a+b+1/2d0*c)**2
     bc(3,5) =  1/2d0*(a+b+c+1/2d0*d)**2

     bc(4,1) = -1/2d0*a**2
     bc(4,2) =  1/48d0*a**3
     bc(4,3) =  1/6d0*(a+1/2d0*b)**3
     bc(4,4) =  1/6d0*(a+b+1/2d0*c)**3
     bc(4,5) =  1/6d0*(a+b+c+1/2d0*d)**3

     bc(5,1) = -1/6d0*a**3
     bc(5,2) =  1/48d0*a**4
     bc(5,3) =  1/6d0*(a+1/2d0*b)**4
     bc(5,4) =  1/6d0*(a+b+1/2d0*c)**4
     bc(5,5) =  1/6d0*(a+b+c+1/2d0*d)**4

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 3) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1)  =  1d0

     bc(1,1) =  1/2d0*a
     bc(1,2) =  (1/2d0*a+b) 
     bc(1,3) =  (1/2d0*a+b+c)
     bc(1,4) =  (1/2d0*a+b+c+d)

     bc(2,1) =  1d0
     bc(2,2) =  1d0
     bc(2,3) =  1d0
     bc(2,4) =  1d0

     bc(3,1) =  1/8d0*a**2
     bc(3,2) =  1/2d0*(1/2d0*a+b)**2
     bc(3,3) =  1/2d0*(1/2d0*a+b+c)**2
     bc(3,4) =  1/2d0*(1/2d0*a+b+c+d)**2

     bc(4,1) =  1/48d0*a**3
     bc(4,2) =  1/6d0*(1/2d0*a+b)**3
     bc(4,3) =  1/6d0*(1/2d0*a+b+c)**3
     bc(4,4) =  1/6d0*(1/2d0*a+b+c+d)**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 4) then
     dd = 0
     do i=1,n
        bcc(i) = 0d0
        indx(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) = -1/2d0*a
     bc(2,2) =  1/2d0*b
     bc(2,3) =  1/2d0*b+c
     bc(2,4) =  1/2d0*b+c+d

     bc(3,1) =  1/8d0*a**2
     bc(3,2) =  1/8d0*b**2
     bc(3,3) =  1/2d0*(1/2d0*b+c)**2
     bc(3,4) =  1/2d0*(1/2d0*b+c+d)**2

     bc(4,1) = -1/48d0*a**3
     bc(4,2) =  1/48d0*b**3
     bc(4,3) =  1/6d0*(1/2d0*b+c)**3
     bc(4,4) =  1/6d0*(1/2d0*b+c+d)**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 5) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(1)  =  1d0

     bc(1,1) =  1/2d0*a
     bc(1,2) = -1/2d0*b 
     bc(1,3) = -(1/2d0*b+c)
     bc(1,4) = -(1/2d0*b+c+d)

     bc(2,1) =  1d0
     bc(2,2) =  1d0
     bc(2,3) =  1d0
     bc(2,4) =  1d0

     bc(3,1) =  1/8d0*a**2
     bc(3,2) =  1/8d0*b**2
     bc(3,3) =  1/2d0*(-(1/2d0*b+c))**2
     bc(3,4) =  1/2d0*(-(1/2d0*b+c+d))**2  

     bc(4,1) =  1/48d0*a**3
     bc(4,2) = -1/48d0*b**3
     bc(4,3) =  1/6d0*(-(1/2d0*b+c))**3
     bc(4,4) =  1/6d0*(-(1/2d0*b+c+d))**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 6) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(2)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) =  0d0
     bc(2,2) = -1/2d0*a 
     bc(2,3) = -(1/2d0*a+b)
     bc(2,4) = -(1/2d0*a+b+c)

     bc(3,1) =  0d0
     bc(3,2) =  1/8d0*a**2
     bc(3,3) =  1/2d0*(-(1/2d0*a+b))**2
     bc(3,4) =  1/2d0*(-(1/2d0*a+b+c))**2  

     bc(4,1) =  0d0
     bc(4,2) = -1/48d0*a**3
     bc(4,3) =  1/6d0*(-(1/2d0*a+b))**3
     bc(4,4) =  1/6d0*(-(1/2d0*a+b+c))**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 7) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo

     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(2)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) =  1/2d0*a 
     bc(2,2) = -1/2d0*b
     bc(2,3) = -(1/2d0*b+c)
     bc(2,4) = -(1/2d0*b+c+d)

     bc(3,1) =  1/8d0*a**2
     bc(3,2) =  1/8d0*b**2
     bc(3,3) = (1/2d0*b+c)**2  
     bc(3,4) = (1/2d0*b+c+d)**2 

     bc(4,1) =  1/48d0*a**3
     bc(4,2) = -1/48d0*b**3
     bc(4,3) = -1/6d0*(1/2d0*b+c  )**3
     bc(4,4) = -1/6d0*(1/2d0*b+c+d)**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)
  elseif (flagv .eq. 8) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo

     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(2)  =  1d0

     bc(1,1) =  0d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) = -1d0
     bc(2,2) = -1/2d0*a 
     bc(2,3) =  1/2d0*a
     bc(2,4) =  1/2d0*a+b

     bc(3,1) = -c
     bc(3,2) =  1/8d0*a**2
     bc(3,3) =  1/8d0*a**2
     bc(3,4) =  1/2d0*(1/2d0*a+b)**2 

     bc(4,1) = -1/2d0*c**2
     bc(4,2) = -1/48d0*a**3
     bc(4,3) =  1/48d0*a**3
     bc(4,4) =  1/6d0*(1/2d0*a+b)**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)

  elseif (flagv .eq. 9) then         
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo

     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(2)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) =-a
     bc(2,2) =-1/2d0*b
     bc(2,3) = 1/2d0*c
     bc(2,4) = d

     bc(3,1) = 1/2d0*a**2 
     bc(3,2) = 1/8d0*b**2
     bc(3,3) = 1/8d0*c**2
     bc(3,4) = 1/2d0*d**2

     bc(4,1) = -1/6d0 *a**3 
     bc(4,2) = -1/48d0*b**3 
     bc(4,3) =  1/48d0*c**3 
     bc(4,4) =  1/6d0 *d**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)
  elseif (flagv .eq. 10) then
     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo

     do i=1,n
        bcc(i) = 0d0
        do j=1,n
           bc(i,j) = 0d0
        enddo
     enddo
     bcc(2)  =  1d0

     bc(1,1) =  1d0
     bc(1,2) =  1d0
     bc(1,3) =  1d0
     bc(1,4) =  1d0

     bc(2,1) = -1/2d0*a
     bc(2,2) =  1/2d0*b
     bc(2,3) =  c
     bc(2,4) =  d

     bc(3,1) = 1/8d0*a**2
     bc(3,2) = 1/8d0*b**2
     bc(3,3) = 1/2d0*c**2
     bc(3,4) = 1/2d0*d**2

     bc(4,1) = -1/48d0*a**3
     bc(4,2) =  1/48d0*b**3
     bc(4,3) =  1/6d0*c**3
     bc(4,4) =  1/6d0*d**3

     call ludcmp(bc,n,np,indx,dd)
     call lubksb(bc,n,np,indx,bcc)
  endif
end subroutine coefb


!     Note: the following routine from "Numerical Recipe "
!     (Press et al. 1989)
!     uses "real*8" (or double precision) as the default
!     type for non-integers. So in your calling program
!     where you have used "implicit none", you must 
!     explicitly declare "real*8 A(n,n)", etc...

!     Given an NxN matrix A (which is stored on an
!     array of dimension Np x Np) this subroutine
!     replaces the 
!     array on output by the LU decomposition of the
!     input matrix A (L-matrix in the lower triangle,
!     U-matrix in the upper triangle of A)
!     INDX(N) is an array which you need to define
!     on the calling program, bu you won't
!     need explicitly...
!     D is just a real*8 number

!     To solve the linear eq. system A*x=b
!     you must
!      call ludcmp(a,n,np,indx,d)
!      call lubksb(a,n,np,indx,b)
!     and after those two calls, the solution x
!     is stored in the array b!


subroutine ludcmp(a,n,np,indx,d)
  implicit none

  integer nmax,i,n,np,j,imax,k,d
  real*8 tiny
  parameter (nmax=100,tiny=1.0d-20)
  real*8 a(np,np),indx(n),vv(nmax)
  real*8 sum,aamax,dum

  d=1d0
  do i=1,n
     aamax=0d0
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     enddo
     if (aamax.eq.0d0) pause 'singular matrix.'
     vv(i)=1d0/aamax
  enddo
  do j=1,n
     if (j.gt.1) then
        do i=1,j-1
           sum=a(i,j)
           if (i.gt.1)then
              do k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              enddo
              a(i,j)=sum
           endif
        enddo
     endif
     aamax=0d0
     do i=j,n
        sum=a(i,j)
        if (j.gt.1)then
           do k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           enddo
           a(i,j)=sum
        endif
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     enddo
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(j.ne.n)then
        if(a(j,j).eq.0d0)a(j,j)=tiny
        dum=1d0/a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        enddo
     endif
  enddo
  if(a(n,n).eq.0d0)a(n,n)=tiny
  return
end subroutine ludcmp


!      from numerical recipes, see comments in 
!      ludcmp.f


subroutine lubksb(a,n,np,indx,b)

  implicit none

  integer i,n,np,ll,ii,j
  real*8 a(np,np),indx(n),b(n),sum

  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        enddo
     else if (sum.ne.0d0) then
        ii=i
     endif
     b(i)=sum
  enddo
  do i=n,1,-1
     sum=b(i)
     if(i.lt.n)then
        do j=i+1,n
           sum=sum-a(i,j)*b(j)
        enddo
     endif
     b(i)=sum/a(i,i)
  enddo

  return
end subroutine lubksb

! ----------------------------------------------------------------------!

subroutine soltrix(u,nstr,nlen,neq,us,npter1,npter2,a,iflag)
  use ctesp

  implicit none



  integer i,j
  integer nstr,nlen,neq,npter1,npter2,iflag,nfin
  real*8 u(nstr,*),us(nstr,*)
  real*8 a(3,*)
  real*8 one,d

  data one/1.0d0/

  nfin = npter1 + neq - 1

  !     lu-decomposition
  if(iflag.eq.0) then
     iflag=1
     do  j=npter1+1,nfin
        d=one/a(2,j-1)
        a(2,j-1)=d
        a(1,j  )=     -a(1,j)*d
        a(2,j  )=a(2,j)+a(1,j)*a(3,j-1)
     enddo
     a(2,nfin)=one/a(2,nfin)
  endif

  !     backsubstitution
  do  j=npter2,npter2+nlen-1
     us(npter1,j) = u(npter1,j)
     do  i=npter1+1,nfin
        us(i,j)=u(i,j)+a(1,i)*us(i-1,j)
     enddo
     us(nfin,j)=us(nfin,j)*a(2,nfin)
     do  i=nfin-1,npter1,-1
        us(i,j)=(us(i,j)-a(3,i)*us(i+1,j))*a(2,i)
     enddo
  enddo

end subroutine soltrix

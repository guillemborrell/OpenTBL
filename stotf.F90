!================EVERYTHING IS THREADED==============
! interpyy,interpxx,differzy ---> should be whithin a Parallel region 
! (if not,only 1 thread works) Works in place!
! 
! ==== Notation =====
! _x: Interpolated in X direction
! _y: Interpolated in Y direction
! d()dx=Derivative in X direction
! d()dy=Derivative in Y direction
! 
! ================ One point statistics ============================
!   u_x = U interpolated to 1/2 points== u_x   INTENT IN!
!   u = U in the regular location (planes) 	INTENT IN!
!   v = V in the regular location (planes)	INTENT IN!
!   w = W in the regular location (planes) 	INTENT IN!
!   p = pressure (planes)			INTENT IN!
!   
!   dvdx = First derivative of v== d(v_x)/dx (@th half cell)
!   dwdx = First derivative of w== d(w_x)/dx (@th half cell)
!   
!   ---------BIG BUFFERS----------------------
!   upencil=u in pencils ==resu			INTENT IN!
!   dudx_zy:Big buffer of size ut (planes)     ==rhsut
!   dudx_pencil:Bif buffer of size resu (pencils) ==rhsut
!   
!   ---------SMALL BUFFERS----------------------
!   buf1 = First free buffer (FOURIER SIZE) ==wkp
!   buf2 = Second   "     "         "       ==wkpo!   
!   buf3 = Third  "     "         "         ==rhswt(0:nz2,ny+1,1)
!   buf4 = 4th  "     "         "           ==bufuphy
!   buf5 = 5th  "     "         "           ==rhsut
!   buf_cor = 6th  "     "         "         ==buf_corr (0:nz2,ncorr,ib:ie,7)
!   bphy1,2,3: 3 Buffers (REAL SIZE)        ==wkp,wkpo,bufuphy (extra buffer)
!   (bphy1==buf1 & bphy2==buf2  & bphy3==buf4 & buf3=dwdx==rhswt & buf5=bplanes=bpencils==rhsut)
!   (dvdx==rhsvt)
!=====================================================================
!NOTE:
! IN POSTPROCESSING USE CTE=1 OR 2 depending if it is k.eq.0 or k.ne.0
! For Spectra and Correlations


!IN RHSP IS CALL AS:   
!        call statsp( wki1t,ut,vt,wt,pt, &   
!      &                rhsvt,rhswt,resu,rhsut,rhsut,&
!      &                wkp,wkpo,rhswt,bufuphy,rhsut,buf_corr,buf_corr,&
!      &                wkp,wkpo,bufuphy,mpiid)



subroutine statsp(u_x,u,v,w,p, &
     &            dvdx,dwdx,upencil,dudx_zy,dudx_pencil,&
     &            buf1,buf2,buf3,buf4,buf5,buf_cor,buf_corp,&
     &            bphy1,bphy2,bphy3,mpiid,communicator)

  use alloc_dns
  use statistics
  use point
  use ctesp
  use omp_lib
  implicit none
  include 'mpif.h'
  integer,intent(in)::communicator
  !Pencils:
  real*8, dimension(nx  ,mpu),intent(in):: upencil
  real*8, dimension(nx  ,mpu):: dudx_pencil
  !Planes:
  complex*16, intent(in), dimension(0:nz2,ny+1,ib:ie):: u_x,u,w
  !-----------Big Buffers and future buffer (dwdx)
  complex*16, dimension(0:nz2,ny+1,ib:ie):: dwdx,dudx_zy
  !Correlations buffers:
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor  !in planes: To storage the 7 variables: u,v,w,p,omega_x,y,z  
  real*8,dimension(nx,mp_corr,7,lxcorr):: buf_corp	       !in pencils: mp_corr=total number of pencils for correlation buffers
  

  complex*16, intent(in), dimension(0:nz2,ny,ib:ie):: v
  complex*16, dimension(0:nz2,ny,ib:ie)::dvdx,p
  !-----------Small Buffers (FOURIER plane)
  complex*16, dimension(0:nz2,ny+1):: buf1,buf2,buf3,buf4,buf5
  !-----------Small Buffers (REAL plane:For triple products)
  real*8,dimension(nz+2,ny+1)::bphy1,bphy2,bphy3
  !--------------------------------------
  integer, intent(in):: mpiid
  real(kind = 8):: cte,aux1,aux2,aux3
  complex*16:: auxc1,auxc2,auxc3
  integer:: i,j,k,jj,i3,ii

  do i=ib,ie             
     call interp_pressure(p(0,1,i),buf5(0,1),ny,v(0,1,i)) !buf5=p(...,i) [outside parallel region!]
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,cte,jj) 
     !$OMP DO SCHEDULE(STATIC)
     do k = 0,nz2
        buf5(k,ny+1) = buf5(k,ny) !copying the last point ny-->ny+1 
     enddo
          
     !================Statistics for Reynolds stresses: uu,vv,ww and uv    
     call interpyy(u_x(0,1,i),buf1(0,1)  ,inyu,cofiuy,inby,ny+1,0,nz1,nz1)   !buf1=(u_x)_y          
     call interpyy(w(0,1,i)   ,buf2(0,1)  ,inyu,cofiuy,inby,ny+1,0,nz1,nz1)   !buf2=w_y 

#ifdef INFOINTER 
        !$OMP DO SCHEDULE(STATIC) 
        do j=1,ny
           v_0(j,i)=v_0(j,i)+v(0,j,i) 
           u_x0(j,i)= u_x0(j,i)+u_x(0,j,i)
           u_xy0(j,i)= u_xy0(j,i)+buf1(0,j)       
           w_0(j,i)= w_0(j,i)+w(0,j,i)           
           w_y0(j,i)= w_y0(j,i)+buf2(0,j)
           dwdx_0(j,i)= dwdx_0(j,i)+dwdx(0,j,i) 
        enddo
#endif      


     call interpyy(dwdx(0,1,i),dwdx(0,1,i),inyu,cofiuy,inby,ny+1,0,nz1,nz1)  !after interpyy: dwdx==[d(w_x)/dx]_y!!  

     !$OMP DO SCHEDULE(STATIC)
     do j = 1,ny
        ! First wavenumber(k=0) means ua,va,wa and Reynolds stresses       
        ua(j,i) = ua(j,i) + buf1(0,j)
        wa(j,i) = wa(j,i) + buf2(0,j)
        va(j,i) = va(j,i) + v(0,j,i)       
        do k = 0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif
           us(j,i) = us(j,i) + cte*buf1(k,j)*dconjg(buf1(k,j))
           ws(j,i) = ws(j,i) + cte*buf2(k,j)*dconjg(buf2(k,j))
           vs(j,i) = vs(j,i) + cte*v(k,j,i) *dconjg(v(k,j,i ))
           uv(j,i) = uv(j,i) + cte*dreal(buf1(k,j)*dconjg(v(k,j,i)))
           uw(j,i) = uw(j,i) + cte*dreal(buf1(k,j)*dconjg(buf2(k,j)))
           vw(j,i) = vw(j,i) + cte*dreal(v(k,j,i)*dconjg(buf2(k,j)))
        end do
     end do

#ifndef NOSPECTRA
     !================Spectra in z
     ! The spectra is computed here because buf1=u and buf2=w will be reused
     ! for pressure and vorticities.
     if(flags(i).ne.0) then
        !$OMP DO SCHEDULE(STATIC)
        do j = 1,nspec                 !Number of points in Y to sample
           jj = jspecy(j,flags(i))     !Point in Y grid for this X station=flags(i)={1,2,3...}                                                           
           do k=0,nz2
              ensu (k,j,flags(i)) = ensu (k,j,flags(i)) + buf1(k,jj)*dconjg(buf1(k,jj))          !(u_x)_y**2
              ensv (k,j,flags(i)) = ensv (k,j,flags(i)) + v(k,jj,i )*dconjg(v(k,jj,i ))          !v**2
              ensw (k,j,flags(i)) = ensw (k,j,flags(i)) + buf2(k,jj)*dconjg(buf2(k,jj))          !(w_y)**2
              ensuv(k,j,flags(i)) = ensuv(k,j,flags(i)) + dreal(buf1(k,jj)*dconjg(v (k,jj,i)))   ![(u_x)_y]*v
              pesp (k,j,flags(i)) = pesp (k,j,flags(i)) + buf5(k,jj)*dconjg(buf5(k,jj))          !p_y**2                       
           enddo
        end do
     endif
#endif
     !$OMP END PARALLEL

#ifndef NOCORR
     !SAVING BUFFERS IN ORDER TO TRANSPOSE THEM LATER
     !buf1=(u_x)_y  !buf2=w_y  !buf5=p_y
     do ii=1,lxcorr 
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,jj,k) SCHEDULE(STATIC)
        do j=1,ncorr
           jj=jspecy(j,ii)
              buf_cor(0:nz2,j,i,1,ii)=buf1(0:nz2,jj)   !u  
              buf_cor(0:nz2,j,i,2,ii)=v   (0:nz2,jj,i) !v
              buf_cor(0:nz2,j,i,3,ii)=buf2(0:nz2,jj)   !w
              buf_cor(0:nz2,j,i,4,ii)=buf5(0:nz2,jj)   !p          
        enddo
     enddo
#endif

     !-------------------Pressure and disipation statistics (part of)
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,jj,k,cte,aux1,aux2)
     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny               
        buf2(0:nz2,j)=buf2(0:nz2,j)*kaz(0:nz2)	!buf2=dwdz=d(w_y)/dz
        dwdz0(j,i)=dwdz0(j,i)+buf2(0,j)          
        do k=0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0          
           endif
           !buf1=u  !buf2=dwdz  !buf5=p IN CORRECT LOCATIONS
           !Pressure statistics
           pdwdz(j,i)=pdwdz(j,i)+cte*dreal(buf5(k,j)*dconjg(buf2(k,j)))    !p*dwdz
           pup  (j,i)=pup  (j,i)+cte*dreal(buf5(k,j)*dconjg(buf1(k,j)))    !p*u
#ifndef NODISSIPATION
           !Dissipation W (dwdz & dwdx)                  
           aux1=buf2(k,j)*dconjg(buf2(k,j))        !aux1=|dwdz|^2           
           aux2=dwdx(k,j,i)*dconjg(dwdx(k,j,i))	  !aux2=|dwdx|^2            
           dispw(j,i)=dispw(j,i)+cte*(aux1+aux2)
#endif
        enddo
     enddo

     !================Statistics for vorticity         
     call differzy(u_x(0,1,i),buf1(0,1),dcyv,dcby,cofcyv,1,1,0,1,ny+1) !buf1=dudy=d(u_x)dy
     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny+1        
        buf2(0:nz2,j) = u_x(0:nz2,j,i)*kaz(0:nz2) 			 !buf2=dudz=d(u_x)dz
     end do
     
     
#ifdef INFOINTER 
        !$OMP DO SCHEDULE(STATIC) 
        do j=1,ny          
           dudz_x0(j,i)= dudz_x0(j,i)+buf2(0,j) 
        enddo
#endif     
     
     
     call interpyy(buf2(0,1),buf2(0,1),inyu,cofiuy,inby,ny+1,0,nz1,nz1)  !buf2=dudz=[d(u_x)dz]_y                        
     !$OMP DO SCHEDULE(STATIC)
     do j = 1,ny
        !buf1=dudy  buf2=dudz  IN CORRECT LOCATIONS           
        vortza(j,i) = vortza(j,i) + (dvdx(0,j,i) - buf1(0,j))	!omega[z]=(dvdx-dudy)
        vortya(j,i) = vortya(j,i) + (buf2(0,j) - dwdx(0,j,i))	!omega[y]=(dudz-dwdx)  
        dvdx0(j,i)=dvdx0(j,i)+dvdx(0,j,i) 
        dvdz0(j,i)=dvdz0(j,i)+v(0,j,i)*kaz(0)
        dwdx0(j,i)=dwdx0(j,i)+dwdx(0,j,i) 
        dudy0(j,i)=dudy0(j,i)+buf1(0,j)  
        dudz0(j,i)=dudz0(j,i)+buf2(0,j) 
        
        do k = 0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif 
!           if(i.eq.600 .and. j.eq.15 .and. k.lt.10) write(*,*) 'value of k',k,'cte',cte,'ompid',ompid,'mpiid',mpiid         
           vortz(j,i)=vortz(j,i)+cte*(dvdx(k,j,i)-buf1(k,j))*dconjg(dvdx(k,j,i)-buf1(k,j)) !aux2=omega[z]=(dvdx-dudy)     
           vorty(j,i)=vorty(j,i)+cte*(buf2(k,j)-dwdx(k,j,i))*dconjg(buf2(k,j)-dwdx(k,j,i)) !aux3=omega[y]=(dudz-dwdx)      
#ifndef NODISSIPATION
           !Dissipation V (dvdx & dvdz)
           aux1=dvdx(k,j,i)*dconjg(dvdx(k,j,i))		  !aux1=|dvdx|^2
           aux2=v(k,j,i)*kaz(k)*dconjg(v(k,j,i)*kaz(k)) 	  !aux2=|dvdz|^2
           dispv(j,i)=dispv(j,i)+cte*(aux1+aux2)           
           !Dissipation U (dudy & dudz)
           aux1=buf1(k,j)*dconjg(buf1(k,j))		  !aux1=|dudy|^2
           aux2=buf2(k,j)*dconjg(buf2(k,j))  		  !aux2=|dudz|^2
           dispu(j,i)=dispu(j,i)+cte*(aux1+aux2)
           !Dissipation UV (duvdz)
           aux1=dreal(buf2(k,j)*dconjg(v(k,j,i)*kaz(k)))          !aux1=|dudz|*|dvdz|
           dispuv(j,i)=dispuv(j,i)+cte*aux1    
#endif                           
        end do
     end do
!----------------------------------------------------------
#ifndef NOSPECTRA     
     if(flags(i).ne.0) then
        !$OMP DO SCHEDULE(STATIC)
        do j = 1,nspec
           jj = jspecy(j,flags(i)) 
           do k=0,nz2                                      
              aux1=(buf2(k,jj)-dwdx(k,jj,i))*dconjg(buf2(k,jj)-dwdx(k,jj,i)) !aux2=omega[y]=(dudz-dwdx)   !dwdx==buf3!!!!!       
              aux2=(dvdx(k,jj,i)-buf1(k,jj))*dconjg(dvdx(k,jj,i)-buf1(k,jj)) !aux3=omega[z]=(dvdx-dudy)                     
              ensomy (k,j,flags(i)) = ensomy (k,j,flags(i)) + aux1
              ensomz (k,j,flags(i)) = ensomz (k,j,flags(i)) + aux2                                
           enddo
        end do
     endif
#endif
     !$OMP END PARALLEL         
#ifndef NOCORR   
     do ii=1,lxcorr 
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,jj) SCHEDULE(STATIC)
        do j=1,ncorr
           jj=jspecy(j,ii)
              buf_cor(0:nz2,j,i,6,ii)=buf2(0:nz2,jj)-dwdx(0:nz2,jj,i)     !omega[y]=(dudz-dwdx)
              buf_cor(0:nz2,j,i,7,ii)=dvdx(0:nz2,jj,i)-buf1(0:nz2,jj)     !omega[z]=(dvdx-dudy)          
        enddo
     enddo
#endif    
!----------------------------------------------------------


     !---------dwdx now is free!!!! Use it as buffer to compute d(u)/dx or whatever
     !---------dwdx==buf3---------   

     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,jj,k,cte,auxc1)
     call interpyy(v(0,1,i),buf3,inyv,cofivy,inby,ny,1,nz1,nz1) 	!(buf3=v_y)  
     call differzy(buf3,buf4,dcyv,dcby,cofcyv,1,1,0,1,ny)               !buf4=dvdy=d(v_y)/dy	          
#ifdef INFOINTER 
        !$OMP DO SCHEDULE(STATIC) 
        do j=1,ny          
           v_y0(j,i)= v_y0(j,i)+buf3(0,j)            
        enddo
#endif     
     call differzy(w(0,1,i),buf3(0,1),dcyv,dcby,cofcyv,1,1,0,1,ny+1)    !buf3=dwdy !!JSS:corrected the coeff May 25 (now v, and 1,1,0,1)
     !===(buf1=dudy  buf2=dudz   buf3=dwdy  buf4=dvdy)=======  
     !$OMP DO SCHEDULE(STATIC)
     do j = 1,ny        
        vortxa(j,i) = vortxa(j,i) + dreal(buf3(0,j)-v(0,j,i)*kaz(0))   !omega[x]=(dwdy-dvdz)    
        dvdy0(j,i)=dvdy0(j,i)+buf4(0,j)   
        dwdy0(j,i)=dwdy0(j,i)+buf3(0,j)               
        do k = 0,nz2 
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif
           auxc1=buf3(k,j)-v(k,j,i)*kaz(k)     !auxc1=(dwdy-dvdz)      
           vortx(j,i) = vortx(j,i)+cte*(auxc1*dconjg(auxc1))
#ifndef NODISSIPATION
           !Dissipation V (dvdx & dvdz & dvdy) ==COMPLETED
           dispv(j,i)=dispv(j,i)+cte*buf4(k,j)*dconjg(buf4(k,j))     !|dvdy|^2        
           !Dissipation W (dwdz & dwdx & dwdy)==COMPLETED     
           dispw(j,i)=dispw(j,i)+cte*buf3(k,j)*dconjg(buf3(k,j))     !|dwdy|^2      
           !Dissipation UV (duvdz & duvdy)
           dispuv(j,i)=dispuv(j,i)+cte*dreal(buf1(k,j)*dconjg(buf4(k,j))) !|dudy|*|dvdy|           
#endif
        end do
     end do

!----------------------------------------------------------
#ifndef NOSPECTRA     
     if(flags(i).ne.0) then
        !$OMP DO SCHEDULE(STATIC)
        do j = 1,nspec
           jj = jspecy(j,flags(i)) 
           do k=0,nz2         
              auxc1=buf3(k,jj)-v(k,jj,i)*kaz(k)   !omega[x]=auxc1=(dwdy-dvdz)  
              ensomx (k,j,flags(i)) = ensomx (k,j,flags(i)) + auxc1*dconjg(auxc1)          
           enddo   
        end do
     endif
#endif
     !$OMP END PARALLEL   
      
#ifndef NOCORR    
     do ii=1,lxcorr 
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,jj,k) SCHEDULE(STATIC)
        do j=1,ncorr
           jj=jspecy(j,ii)
              buf_cor(0:nz2,j,i,5,ii)=buf3(0:nz2,jj)-v(0:nz2,jj,i)*kaz(0:nz2) !omega[x]=(dwdy-dvdz)                     
        enddo
     enddo
#endif
!----------------------------------------------------------
     
     !================Statistics for pressure
     !p*u      =p*u   +p*v   +p*w    --> pup  ,pvp  ,pwp   (3)
     !p*grad(u)=p*dudx+p*dudy+p*dudz --> pdudx,pdudy,pdudz (3)
     !p*div(u) =p*dudx+p*dvdy+p*dwdz -->      ,pdvdy,pdwdz (2)
     !--------------------------------------------   
     !===(buf1=dudy  buf2=dudz   buf3=dwdy  buf4=dvdy buf5=p)=======

     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,cte)
     !$OMP DO SCHEDULE(STATIC)
     do j = 1,ny   
        pm(j,i) = pm(j,i)+buf5(0,j)        
        do k = 0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif
           pp(j,i) = pp(j,i) + cte*buf5(k,j)*dconjg(buf5(k,j))  ! (prms+p0)^2
           pvp  (j,i) = pvp  (j,i)+cte*dreal(buf5(k,j)*dconjg(v(k,j,i)))    !p*v
           pdvdy(j,i) = pdvdy(j,i)+cte*dreal(buf5(k,j)*dconjg(buf4(k,j)))   !p*dvdy
           pdudy(j,i) = pdudy(j,i)+cte*dreal(buf5(k,j)*dconjg(buf1(k,j)))   !p*dudy
           pdvdx(j,i) = pdvdx(j,i)+cte*dreal(buf5(k,j)*dconjg(dvdx(k,j,i))) !p*dvdx
        end do
     end do

#ifndef NO3PRODUCTS
     !=====================Triple products 
     call interpyy(u_x(0,1,i),buf2,inyu,cofiuy,inby,ny+1,0,nz1,nz1) 		  !buf2=[u_x]_y in Fourier (faster)         
     call fourxz(buf2(0:nz2,jbf1:jef1),bphy1(1:nz+2,jbf1:jef1),1,ny+1,jbf1,jef1)  !bphy1=[u_x]_y real              
     !$OMP BARRIER         !(Now: bphy1==buf1=wkp USED. U REAL)      
     call interpyy(w(0,1,i),buf2,inyu,cofiuy,inby,ny+1,0,nz1,nz1) 		 !w_y fourier in Fourier (faster)
     call fourxz(buf2(0:nz2,jbf1:jef1),bphy3(1:nz+2,jbf1:jef1),1,ny+1,jbf1,jef1) !bphy3=w_y real
     !$OMP BARRIER	   !(Now: bphy3==buf4=bufuphy USED. W REAL)         
     call fourxz(v(0:nz2,jbf2:jef2,i),bphy2(1:nz+2,jbf2:jef2),1,ny,jbf2,jef2) 	 !bphy2=v real 
     !$OMP BARRIER	   !(Now: bphy2==buf2=wkpo USED. V REAL) 
    
     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny+1
        do k=1,nz
           u3(j,i)=u3(j,i)+bphy1(k,j)**3             
           w2u(j,i)=w2u(j,i)+bphy1(k,j)*bphy3(k,j)**2
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny
        do k=1,nz
           v3(j,i)=v3(j,i)+bphy2(k,j)**3  
           w2v(j,i)=w2v(j,i)+bphy2(k,j)*bphy3(k,j)**2
           v2u(j,i)=v2u(j,i)+bphy1(k,j)*bphy2(k,j)**2
           u2v(j,i)=u2v(j,i)+bphy2(k,j)*bphy1(k,j)**2
        enddo
     enddo
#endif
     !$OMP END PARALLEL
  end do


  !========================== Out of the loop ib:ie (global)=============================
  !dudx_zy & dudx_pencil only used HERE!!! (==buf5;it was free to use up to here)   
  call differxx(upencil,dudx_pencil,dcxv,dcbx,cofcxv,2,mpu)  !dudx_pencil=d(u)dx pencils @th center of the cell
  call chx2p(dudx_pencil,dudx_zy,dwdx,mpiid,ny+1,communicator)            !dudx_zy =d(u)dx in planes . dwdx==buffer
  !==============================================================================
   
  do i=ib,ie
     call interp_pressure(p(0,1,i),buf1(0,1),ny,v(0,1,i)) !buf1=p(...,i) [outside parallel region!]

     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,cte)
     !$OMP DO SCHEDULE(STATIC)
     do k = 0,nz2
        buf1(k,ny+1) = buf1(k,ny) !copying the last point ny-->ny+1 
     enddo    
     !---------------------------------------------

     !Case for u: dudx_zy=dudx
#ifdef INFOINTER 
        !$OMP DO SCHEDULE(STATIC) 
        do j=1,ny          
           dudx_0(j,i)= dudx_0(j,i)+dudx_zy(0,j,i)            
        enddo
#endif                   
     call interpyy(dudx_zy(0,1,i),dudx_zy(0,1,i),inyu,cofiuy,inby,ny+1,0,nz1,nz1) !dudx_zy=dudx=(du/dx)_y           
     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny
        dudx0(j,i)=dudx0(j,i)+dudx_zy(0,j,i)                      
        do k=0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif
#ifndef NODISSIPATION 
           !Dissipation U (dudy & dudz & dudx)==COMPLETED                              
           dispu(j,i)=dispu(j,i)+cte*dudx_zy(k,j,i)*dconjg(dudx_zy(k,j,i))		!|dudx|^2  
           !Dissipation UV (duvdz & duvdy & duvdx)==COMPLETED        
           dispuv(j,i)=dispuv(j,i)+cte*dreal(dudx_zy(k,j,i)*dconjg(dvdx(k,j,i))) 	!|dudx|*|dvdx|
#endif
           pdudx(j,i)=pdudx(j,i)+cte*dreal(buf1(k,j)*dconjg(dudx_zy(k,j,i)))    !p*dudx
        enddo
     enddo
     !$OMP END PARALLEL   
  enddo

!Big buffers are free: dwdx,dvdx,rhsut
!===============CORRELATIONS==================
#ifndef NOCORR
call compute_corr(coru,corv,corw,coruv,corp,corox,coroy,coroz,buf_cor,buf_corp,dwdx,dudx_zy,mpiid,communicator)
#endif
end subroutine statsp


!==================================================================
!==================================================================
!=================================================================


subroutine compute_corr(cor1,cor2,cor3,cor4,cor5,cor6,cor7,cor8,buf_cor,buf_corp,buf_change,buf_cor2,mpiid,communicator)
  use ctesp,only: nx,ny,nz,nz2,ncorr,lxcorr,xcorpoint,nxp,xci,xco
  use point
  use omp_lib
  implicit none
  include 'mpif.h'
  integer,intent(in)::communicator
  integer:: mpiid,ii,j
   real*8,dimension(nx,mp_corr2,lxcorr)::cor1,cor2,cor3,cor4,cor5,cor6,cor7,cor8
  !Correlations buffers:
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor  !in planes: To storage the 7 variables: u,v,w,p,omega_x,y,z
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor2
  real*8,dimension(nx,mp_corr,7,lxcorr):: buf_corp	       !in pencils: mp_corr=total number of pencils for correlation buffers
  !Buffer for change
  complex*16, dimension(0:nz2,ny+1,ib:ie):: buf_change    
    
  !$OMP PARALLEL WORKSHARE
  buf_cor2=buf_cor  !Looks like the change does not work in place for 
  !$OMP END PARALLEL WORKSHARE

   !!!!!!!!!!!CHANGING CORRELATIONS BUFFERS TO PENCILS
  do ii=1,lxcorr
     do j=1,7       
        call chp2x(buf_corp(1,1,j,ii),buf_cor2(0,1,ib,j,ii),buf_change,mpiid,ncorr,communicator)       
     enddo
     !Computing correlations:     
     call c_corr(cor1 (1,1,ii),buf_corp(1,1,1,ii),buf_corp(1,1,1,ii),ii,mpiid)!cor_u   
     call c_corr(cor2 (1,1,ii),buf_corp(1,1,2,ii),buf_corp(1,1,2,ii),ii,mpiid)!cor_v   
     call c_corr(cor3 (1,1,ii),buf_corp(1,1,3,ii),buf_corp(1,1,3,ii),ii,mpiid)!cor_w   
     call c_corr(cor4 (1,1,ii),buf_corp(1,1,1,ii),buf_corp(1,1,2,ii),ii,mpiid)!cor_uv   
     call c_corr(cor5 (1,1,ii),buf_corp(1,1,4,ii),buf_corp(1,1,4,ii),ii,mpiid)!cor_p   
     call c_corr(cor6 (1,1,ii),buf_corp(1,1,5,ii),buf_corp(1,1,5,ii),ii,mpiid)!cor_ox   
     call c_corr(cor7 (1,1,ii),buf_corp(1,1,6,ii),buf_corp(1,1,6,ii),ii,mpiid)!cor_oy   
     call c_corr(cor8 (1,1,ii),buf_corp(1,1,7,ii),buf_corp(1,1,7,ii),ii,mpiid)!cor_oz   
  enddo
end subroutine compute_corr


!==================================================================
!==================================================================
!=================================================================
!!!!! SUBROUTINE C_CORR INFO:
!JAS JUNE 16TH 2010
!buffer1=buffer2 for all the correlations except Cor_uv
!buffer1 & buffer2 are Complex*16 when in planes[0:nz2,1:ncorr,ib:ie]...R8 in pencils
!correlation: It is exactly half of the buffer1 or 2, since it is one R*8 plane... R8 in pencils (half the buffer1 pencils)
!pcib,pcie:   Beginning and end of the pencils kind--> mp_corr =[2*(nz2+1)]*ncorr=nz1  *ncorr. mp_corr =total number of pencils for the node
!pcib2,pcie2: Beginning and end of the pencils kind--> mp_cor2r=(nz2+1)*ncorr    =nz1/2*ncorr. mp_corr2=total number of pencils for the node
!Note: Constant multiplying the mode to be add at postprocessing.

subroutine c_corr(cor,buffer1,buffer2,j,rank)
  use ctesp,only: nx,ny,nz,nz2,lxcorr,xcorpoint,nxp,xci,xco
  use point
  implicit none
  include 'mpif.h'
  integer:: r1,r2,jj,j,i,ii,i3,kk,liminf,limsup,ip,rank,pn,d994,xx
  real*8,dimension(nx,mp_corr):: buffer1,buffer2
  real*8,dimension(nx,mp_corr2):: cor   !EXACTLY Half pencils than buffer1 or 2  
  real*8:: auxr1,auxr2
  d994=nxp(j)
  xx=xcorpoint(j)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pn,kk,ii,auxr1,auxr2,liminf,limsup) SCHEDULE(STATIC)
  do pn=1,mp_corr2  !pencil=Each 2 buffers pencils=1 correlation pencil                     
     do kk=-d994,d994  !+/-d99/4 averaging             
        auxr1=buffer2(xx+kk,2*pn-1)
        auxr2=buffer2(xx+kk,2*pn)
        !Correlations limits
        liminf=xci+d994+kk
        limsup=xco-d994+kk
        !Shifted X units in order to save it: ip=i-i3+X=i-kk. R=i-i3=correlation index                           
        cor(xci+d994:xco-d994,pn)=cor(xci+d994:xco-d994,pn)+buffer1(liminf:limsup,2*pn-1)*auxr1+buffer1(liminf:limsup,2*pn)*auxr2                
     enddo
  enddo
end subroutine c_corr




! **********************  minimal txapuza  ************* 
subroutine ministats(uiy,viy,wiy,uner,hy,i)
  use point
  use ctesp 
  implicit none

  complex*16 uiy(0:nz2,ny+1),viy(0:nz2,ny),wiy(0:nz2,ny+1) 
  real*8 uner(15),hy(0:ny)
  integer i

  real*8 aux1(ny+1),aux2(ny),aux3(ny+1) ,kk3
  integer k,j,jj

  ! Energy
  if (i==1) then         
     aux1 = sum(uiy*conjg(uiy),1)
     aux2 = sum(viy*conjg(viy),1)
     aux3 = sum(wiy*conjg(wiy),1)
     uner(1) = sum(aux1*hy(0:ny))  
     uner(2) = sum(aux2*hy(0:ny-1 ))
     uner(3) = sum(aux3*hy(0:ny))
  endif

  if (i==(nx/2)) then
     aux1 = 2d0*sum(uiy(1:nz2,:)*conjg(uiy(1:nz2,:)),1)
     aux2 = 2d0*sum(viy(1:nz2,:)*conjg(viy(1:nz2,:)),1)
     aux3 = 2d0*sum(wiy(1:nz2,:)*conjg(wiy(1:nz2,:)),1)
     uner(4) =  real(uiy(0,1))     
     uner(5) =  real(uiy(0,2))     
     uner(6) =  aux1(1)     
     uner(7) =  aux1(2)     
     uner(8) =  real(viy(0,1))     
     uner(9) =  real(viy(0,2))     
     uner(10) = aux2(1)     
     uner(11) = aux2(2)     
     uner(12) =  real(wiy(0,1))     
     uner(13) =  real(wiy(0,2))     
     uner(14) =  aux3(1)     
     uner(15) =  aux3(2)
     aux1 = aux1 + real(uiy(0,1))**2
     aux2 = aux2 + real(viy(0,1))**2
     aux3 = aux3 + real(wiy(0,1))**2
     uner(1) = sum(aux1*hy(0:ny))
     uner(2) = sum(aux2*hy(0:ny-1 ))
     uner(3) = sum(aux3*hy(0:ny))
  endif
end subroutine ministats










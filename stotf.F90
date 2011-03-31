!================EVERYTHING IS THREADED==============
!Modified by JSS, March 2011 (new interpolations and derivatives)
!
! Interpolation and Derivative in Y direction ---> should be whithin a Parallel region 
! (if not,only 1 thread works) Does not Work in place! (not tridiagonal [B]*f)
! 
!
! ==== Notation ====================
! _x: Interpolated in X direction
! _y: Interpolated in Y direction
! d()dx=Derivative in X direction
! d()dy=Derivative in Y direction
! 
!
! ================ One point statistics ============================
!   u_x = U interpolated to 1/2 points== u_x    ==wk1t
!   v = V in the regular location (planes)	==vt
!   w = W in the regular location (planes) 	==wt
!   p = pressure (planes)			==pt
!   dudx = First derivative of u==d(u_x)dx      ==wki3t
!   dvdx = First derivative of v==dvdx          ==rhsvt
!   dwdx = First derivative of w== dwdx         ==rhswt
!
!   
!   ---------BIG BUFFERS----------------------
!   buf_big =rhsut. (Used for correlation changes)
!
!
!   
!   ---------SMALL BUFFERS----------------------
!   buf1 = First free buffer (FOURIER SIZE) ==wkp
!   buf2 = Second   "     "         "       ==wkpo!   
!   buf3 = Third  "     "         "         ==wki2t(0:nz2,ny,1)
!   buf4 = 4th  "     "         "           ==bufuphy
!   buf5 = 5th  "     "         "           ==rhsut
!   buf6 =                                  ==rhswt
!   buf7 =                                  ==wki3t
!   buf8 =                                  ==wki1t
!   buf_cor = 6th  "     "         "         ==buf_corr (0:nz2,ncorr,ib:ie,7)
!   bphy1,2,3: 3 Buffers (REAL SIZE)        ==wkp,wkpo,bufuphy (extra buffer)
!   (bphy1==buf1 & bphy2==buf2  & bphy3==buf4 & buf3=dwdx==rhswt & buf5=bplanes=bpencils==rhsut)
!   (dvdx==rhsvt)
!=====================================================================
!NOTE:
! IN POSTPROCESSING USE CTE=1 OR 2 depending if it is k.eq.0 or k.ne.0
! For Spectra and Correlations


subroutine statsp(u_x,v,w,p, &
     &            dudx,dvdx,dwdx,&
     &            buf1,buf2,buf3,buf4,&
     &            buf5,buf6,buf7,buf8,&
     &            buf_cor,buf_corp,buf_big,&
     &            bphy1,bphy2,bphy3,&
     &            mpiid,communicator)

  use alloc_dns
  use statistics
  use point
  use ctesp
  use omp_lib
  implicit none
  include 'mpif.h'
  integer,intent(in)::communicator
  !Planes:
  complex*16, dimension(0:nz2,ny+1,ib:ie):: u_x,w,dudx,dwdx,buf_big
  complex*16, dimension(0:nz2,ny,ib:ie)::dvdx
  complex*16, intent(in), dimension(0:nz2,ny,ib:ie):: v,p
  !Correlations buffers:
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor  !in planes: To storage the 7 variables: u,v,w,p,omega_x,y,z  
  real*8,dimension(nx,mp_corr,7,lxcorr):: buf_corp	       !in pencils: mp_corr=total number of pencils for correlation buffers
  !-----------Small Buffers (FOURIER plane)
  complex*16, dimension(0:nz2,ny+1):: buf1,buf2,buf6,buf4,buf5,buf7,buf8
  complex*16, dimension(0:nz2,ny  ):: buf3
  !-----------Small Buffers (REAL plane:For triple products)
  real*8,dimension(nz+2,ny+1)::bphy1,bphy2,bphy3
  !--------------------------------------
  integer, intent(in):: mpiid
  real(kind = 8):: cte,aux1,aux2,aux3
  complex*16:: auxc1,auxc2,auxc3
  integer:: i,j,k,jj,i3,ii


  do i=ib,ie             
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,jj,k,cte,aux1,aux2,aux3,auxc1)
     call interpy_new(u_x(0,1,i) ,buf1,fd_iy ,0)   !buf1=(u_x)_y          
     call interpy_new(w(0,1,i)   ,buf2,fd_iy ,0)   !buf2=w_y 
     call interpy_new(dwdx(0,1,i),buf3,fd_iy ,0)   !buf3==dwdx_y!!  
     call interpy_new(dudx(0,1,i),buf4,fd_iy ,0)   !buf4==dudx_y!!  
     call interpy_new(p(0,1,i)   ,buf5,fd_iyp,1)   !buf5=p 
     !$OMP BARRIER
     call diffy_inplace(buf1    ,buf6,fd_dvdy)  !buf6=dudy (barrier cause, buf6=dwdx)
     call diffy_inplace(buf2    ,buf7,fd_dvdy)  !buf7=dwdy (barrier cause, buf7=dudx)
     call diffy_inplace(v(0,1,i),buf8,fd_dvdy)  !buf8=dvdy (barrier cause, buf8=u_x)

     !==================================================================================
     !  At this point, everything collocated: p,v,u,w,dvdx,dwdx,dudx,dwdy,dudy,dvdy
     !  buf1=u  !buf2=w !buf3=dwdx !buf4=dudx !buf5=p !buf6=dudy !buf7=dwdy !buf8=dvdy
     !==================================================================================

     !$OMP DO SCHEDULE(STATIC)
     do j = 1,ny
        !First wavenumber(k=0) means ua,va,wa and Reynolds stresses       
        ua(j,i) = ua(j,i) + buf1(0,j)
        wa(j,i) = wa(j,i) + buf2(0,j)
        va(j,i) = va(j,i) + v(0,j,i)       
        pm(j,i) = pm(j,i) + buf5(0,j)        

        dudx0(j,i)=dudx0(j,i)+buf4(0,j)                      
        dvdx0(j,i)=dvdx0(j,i)+dvdx(0,j,i) 
        dwdx0(j,i)=dwdx0(j,i)+buf3(0,j)                       
        dudy0(j,i)=dudy0(j,i)+buf6(0,j)  
        dvdy0(j,i)=dvdy0(j,i)+buf8(0,j)   
        dwdy0(j,i)=dwdy0(j,i)+buf7(0,j)               

        vortya(j,i)=vortya(j,i)-buf3(0,j)                       !omega[y]=(dudz-dwdx) (dudz=0 for k=0)
        vortxa(j,i)=vortxa(j,i)+dreal(buf7(0,j))                !omega[x]=(dwdy-dvdz)    
        vortza(j,i)=vortza(j,i)+dreal((dvdx(0,j,i)-buf6(0,j)))	!omega[z]=(dvdx-dudy)

        do k = 0,nz2
           if(k.eq.0)then;cte=1d0;else;cte=2d0
           endif
           !================Statistics for Reynolds stresses: uu,vv,ww and uv    
           us(j,i) = us(j,i) + cte*buf1(k,j)*dconjg(buf1(k,j))
           ws(j,i) = ws(j,i) + cte*buf2(k,j)*dconjg(buf2(k,j))
           vs(j,i) = vs(j,i) + cte*v(k,j,i) *dconjg(v(k,j,i ))
           uv(j,i) = uv(j,i) + cte*dreal(buf1(k,j)*dconjg(v(k,j,i)))
           uw(j,i) = uw(j,i) + cte*dreal(buf1(k,j)*dconjg(buf2(k,j)))
           vw(j,i) = vw(j,i) + cte*dreal(v(k,j,i)*dconjg(buf2(k,j)))

           !================Statistics for pressure
           !p*u      =p*u   +p*v   +p*w    --> pup  ,pvp  ,pwp   (3)
           !p*grad(u)=p*dudx+p*dudy+p*dudz --> pdudx,pdudy,pdudz (3)
           !p*div(u) =p*dudx+p*dvdy+p*dwdz -->      ,pdvdy,pdwdz (2)
           !--------------------------------------------   
           pup  (j,i)=pup  (j,i)+cte*dreal(buf5(k,j)*dconjg(buf1(k,j)))   !p*u
           pvp  (j,i)=pvp  (j,i)+cte*dreal(buf5(k,j)*dconjg(v(k,j,i)))    !p*v
           pp   (j,i)= pp(j,i)  +cte*buf5(k,j)*dconjg(buf5(k,j))          !(prms+p0)^2

           !------------Z Derivatives-------------------------------------:
           buf2(k,j)=buf2(k,j)*kaz(k)	!Now: buf2=dwdz
           buf1(k,j)=buf1(k,j)*kaz(k)	!Now: buf1=dudz
           !buf1=dudz  !buf2=dwdz !buf3=dwdx !buf4=dudx !buf5=p !buf6=dudy !buf7=dwdy !buf8=dvdy

           pdwdz(j,i)=pdwdz(j,i)+cte*dreal(buf5(k,j)*dconjg(buf2(k,j)))   !p*dwdz
           pdudx(j,i)=pdudx(j,i)+cte*dreal(buf5(k,j)*dconjg(buf4(k,j)))   !p*dudx
           pdvdx(j,i)=pdvdx(j,i)+cte*dreal(buf5(k,j)*dconjg(dvdx(k,j,i))) !p*dvdx
           pdvdy(j,i)=pdvdy(j,i)+cte*dreal(buf8(k,j)*dconjg(buf8(k,j)))   !p*dvdy
           pdudy(j,i)=pdudy(j,i)+cte*dreal(buf6(k,j)*dconjg(buf6(k,j)))   !p*dudy

           !================Statistics for vorticity
           auxc1=buf7(k,j)-v(k,j,i)*kaz(k)     !auxc1=(dwdy-dvdz)      
           vortx(j,i) = vortx(j,i)+cte*(auxc1*dconjg(auxc1))
           vorty(j,i)=vorty(j,i)+cte*(buf1(k,j)-buf3(k,j))*dconjg(buf1(k,j)-buf3(k,j))    !o_y=(dudz-dwdx) 
           vortz(j,i)=vortz(j,i)+cte*(dvdx(k,j,i)-buf6(k,j))*dconjg(dvdx(k,j,i)-buf6(k,j)) !o_z=(dvdx-dudy)   

#ifndef NODISSIPATION
           !================Statistics for dissipation
           !Dissipation U 
           aux1=buf4(k,j)*dconjg(buf4(k,j))		 !aux1=|dudx|^2
           aux2=buf1(k,j)*dconjg(buf1(k,j))!aux2=|dudz|^2
           aux3=buf6(k,j)*dconjg(buf6(k,j))		 !aux1=|dudy|^2
           dispu(j,i)=dispu(j,i)+cte*(aux1+aux2+aux3)
           !Dissipation V 
           aux1=dvdx(k,j,i)*dconjg(dvdx(k,j,i))		 !aux1=|dvdx|^2
           aux2=v(k,j,i)*kaz(k)*dconjg(v(k,j,i)*kaz(k))  !aux2=|dvdz|^2
           aux3=buf8(k,j)*dconjg(buf8(k,j))		 !aux3=|dvdy|^2
           dispv(j,i)=dispv(j,i)+cte*(aux1+aux2+aux3)         
           !Dissipation W 
           aux1=buf2(k,j)*dconjg(buf2(k,j))              !aux1=|dwdz|^2           
           aux2=buf3(k,j)*dconjg(buf3(k,j))	         !aux2=|dwdx|^2            
           aux3=buf7(k,j)*dconjg(buf7(k,j))		 !aux3=|dwdy|^2
           dispw(j,i)=dispw(j,i)+cte*(aux1+aux2+aux3)
           !Dissipation UV 
           aux1=dreal(buf2(k,j)*dconjg(v(k,j,i)*kaz(k))) !aux1=|dudz|*|dvdz|   
           aux2=dreal(buf4(k,j)*dconjg(dvdx(k,j,i)))     !aux2=|dudx|*|dvdx|
           aux3=dreal(buf6(k,j)*dconjg(buf8(k,j)))	 !aux3=|dudy|*|dvdy|           
           dispuv(j,i)=dispuv(j,i)+cte*(aux1+aux2+aux3)   
#endif
        end do
     end do

     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny
        buf1(:,j)=buf1(:,j)/kaz(:)  !dudz->u
        buf2(:,j)=buf2(:,j)/kaz(:)  !dwdz->w
     enddo

#ifndef NOSPECTRA
     !================Spectra in z
     if(flags(i).ne.0) then
        !$OMP DO SCHEDULE(STATIC)
        do j = 1,nspec                 !Number of points in Y to sample
           jj = jspecy(j,flags(i))     !Point in Y grid for this X station=flags(i)={1,2,3...}
           do k=0,nz2
              ensu (k,j,flags(i)) =ensu (k,j,flags(i)) +buf1(k,jj)*dconjg(buf1(k,jj))          !(u_x)_y**2
              ensv (k,j,flags(i)) =ensv (k,j,flags(i)) +v(k,jj,i )*dconjg(v(k,jj,i ))          !v**2
              ensw (k,j,flags(i)) =ensw (k,j,flags(i)) +buf2(k,jj)*dconjg(buf2(k,jj))          !(w_y)**2
              ensuv(k,j,flags(i)) =ensuv(k,j,flags(i)) +dreal(buf1(k,jj)*dconjg(v (k,jj,i)))   ![(u_x)_y]*v

              auxc1=buf7(k,jj)-v(k,jj,i)*kaz(k)     !o_x: auxc1=(dwdy-dvdz)                    
              ensomx(k,j,flags(i))=ensomx(k,j,flags(i))+auxc1*dconjg(auxc1) ! o_x
              auxc1=buf1(k,jj)*kaz(k)-buf3(k,jj)    !o_y: auxc1=(dudz-dwdx) 
              ensomy(k,j,flags(i))=ensomy(k,j,flags(i))+auxc1*dconjg(auxc1) ! o_y
              auxc1=dvdx(k,jj,i)-buf6(k,jj)         !o_z: auxc1=(dvdx-dudy)   
              ensomz(k,j,flags(i))=ensomz(k,j,flags(i))+auxc1*dconjg(auxc1) ! o_z
           end do
        enddo
     endif
#endif
     !$OMP END PARALLEL

#ifndef NOCORR
     !SAVING BUFFERS IN ORDER TO TRANSPOSE THEM LATER
     !buf1=u  !buf2=w !buf3=dwdx !buf4=dudx !buf5=p !buf6=dudy !buf7=dwdy !buf8=dvdy
     do ii=1,lxcorr 
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,jj,k) SCHEDULE(STATIC)
        do j=1,ncorr
           jj=jspecy(j,ii)
           buf_cor(0:nz2,j,i,1,ii)=buf1(0:nz2,jj)   !u  
           buf_cor(0:nz2,j,i,2,ii)=v   (0:nz2,jj,i) !v
           buf_cor(0:nz2,j,i,3,ii)=buf2(0:nz2,jj)   !w
           buf_cor(0:nz2,j,i,4,ii)=buf5(0:nz2,jj)   !p          
           buf_cor(0:nz2,j,i,5,ii)=buf7(0:nz2,jj)-v(0:nz2,jj,i)*kaz(0:nz2) !omega[x]=(dwdy-dvdz)         
           buf_cor(0:nz2,j,i,6,ii)=buf1(0:nz2,jj)*kaz-buf3(0:nz2,jj)       !omega[y]=(dudz-dwdx)
           buf_cor(0:nz2,j,i,7,ii)=dvdx(0:nz2,jj,i)-buf6(0:nz2,jj)         !omega[z]=(dvdx-dudy)          
        enddo
     enddo
#endif

     !dwdx (buf3,buf6),dudx(buf4) is free,dvdx also

#ifndef NO3PRODUCTS
     !=====================Triple products Statistics
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k) 
     call fourxz(buf1(0:nz2,jbf2:jef2),bphy3(1:nz+2,jbf2:jef2),1,ny,jbf2,jef2)!bphy3=[u_x]_y real
     !$OMP BARRIER  !(bphy3==buf4=bufuphy USED. U REAL)         
     call fourxz(buf2(0:nz2,jbf2:jef2),bphy1(1:nz+2,jbf2:jef2),1,ny,jbf2,jef2)!bphy1=w_y real
     !$OMP BARRIER  !(bphy1==buf1=wkp     USED. W REAL)  	   
     call fourxz(v(0:nz2,jbf2:jef2,i),bphy2(1:nz+2,jbf2:jef2),1,ny,jbf2,jef2) !bphy2=v real 
     !$OMP BARRIER  !(bphy2==buf2=wkpo    USED. V REAL) 

     !$OMP DO SCHEDULE(STATIC)
     do j=1,ny
        do k=1,nz
           u3(j,i)=u3(j,i)+bphy3(k,j)**3             
           w2u(j,i)=w2u(j,i)+bphy3(k,j)*bphy1(k,j)**2
           v3(j,i)=v3(j,i)+bphy2(k,j)**3  
           w2v(j,i)=w2v(j,i)+bphy2(k,j)*bphy1(k,j)**2
           v2u(j,i)=v2u(j,i)+bphy3(k,j)*bphy2(k,j)**2
           u2v(j,i)=u2v(j,i)+bphy2(k,j)*bphy3(k,j)**2
        enddo
     enddo
     !$OMP END PARALLEL
#endif

  enddo



  !Big buffers are free: dwdx,dvdx,rhsut

#ifndef NOCORR
  !===============CORRELATIONS==================
  call compute_corr(coru,corv,corw,coruv,corp,corox,coroy,coroz,coruw,corvw,&
       &            buf_cor,buf_corp,dwdx,buf_big,mpiid,communicator)
#endif
end subroutine statsp


!==================================================================
!==================================================================
!=================================================================


subroutine compute_corr(cor1,cor2,cor3,cor4,cor5,cor6,cor7,cor8,cor9,cor10,&
     &buf_cor,buf_corp,buf_change,buf_cor2,mpiid,communicator)

  use ctesp,only: nx,ny,nz,nz2,ncorr,lxcorr,xcorpoint,nxp,xci,xco
  use point
  use omp_lib
  implicit none
  include 'mpif.h'
  integer,intent(in)::communicator
  integer:: mpiid,ii,j
  real*8,dimension(nx,mp_corr2,lxcorr)::cor1,cor2,cor3,cor4,cor5,cor6,cor7,cor8,cor9,cor10
  !Correlations buffers:
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor  !in planes: To storage the 7 variables: u,v,w,p,omega_x,y,z
  complex*16, dimension(0:nz2,ncorr,ib:ie,7,lxcorr):: buf_cor2
  real*8,dimension(nx,mp_corr,7,lxcorr):: buf_corp	       !in pencils: mp_corr=total number of pencils for correlation buffers
  !Buffer for change
  complex*16, dimension(0:nz2,ny+1,ib:ie):: buf_change    

  !$OMP PARALLEL WORKSHARE
  buf_cor2=buf_cor  !change does not work in place for correlations 
  !$OMP END PARALLEL WORKSHARE

!!!!!!!!!!!CHANGING CORRELATIONS BUFFERS TO PENCILS
  do ii=1,lxcorr
     do j=1,7       
        call chp2x(buf_corp(1,1,j,ii),buf_cor2(0,1,ib,j,ii),buf_change,mpiid,ncorr,communicator)       
     enddo
     !Computing correlations:     
     call c_corr(cor1 (1,1,ii),buf_corp(1,1,1,ii),buf_corp(1,1,1,ii),ii,mpiid,0)!cor_u   
     call c_corr(cor2 (1,1,ii),buf_corp(1,1,2,ii),buf_corp(1,1,2,ii),ii,mpiid,0)!cor_v   
     call c_corr(cor3 (1,1,ii),buf_corp(1,1,3,ii),buf_corp(1,1,3,ii),ii,mpiid,0)!cor_w   
     call c_corr(cor4 (1,1,ii),buf_corp(1,1,1,ii),buf_corp(1,1,2,ii),ii,mpiid,0)!cor_uv   
     call c_corr(cor5 (1,1,ii),buf_corp(1,1,4,ii),buf_corp(1,1,4,ii),ii,mpiid,0)!cor_p   
     call c_corr(cor6 (1,1,ii),buf_corp(1,1,5,ii),buf_corp(1,1,5,ii),ii,mpiid,0)!cor_ox   
     call c_corr(cor7 (1,1,ii),buf_corp(1,1,6,ii),buf_corp(1,1,6,ii),ii,mpiid,0)!cor_oy   
     call c_corr(cor8 (1,1,ii),buf_corp(1,1,7,ii),buf_corp(1,1,7,ii),ii,mpiid,0)!cor_oz   
     call c_corr(cor9 (1,1,ii),buf_corp(1,1,1,ii),buf_corp(1,1,3,ii),ii,mpiid,1)!cor_uw   
     call c_corr(cor10(1,1,ii),buf_corp(1,1,2,ii),buf_corp(1,1,3,ii),ii,mpiid,1)!cor_vw   
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
!Tipo: Specify if it has odd (1) or even (0) simmetry when computing correlation.
!Note: Constant multiplying the mode to be add at postprocessing.

subroutine c_corr(cor,buffer1,buffer2,j,rank,tipo)
  use ctesp,only: nx,ny,nz,nz2,lxcorr,xcorpoint,nxp,xci,xco
  use point
  implicit none
  include 'mpif.h'
  integer:: r1,r2,jj,j,i,ii,i3,kk,liminf,limsup,ip,rank,pn,d994,xx,tipo
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
        if(tipo.eq.0) then
           cor(xci+d994:xco-d994,pn)=cor(xci+d994:xco-d994,pn)+buffer1(liminf:limsup,2*pn-1)*auxr1+&
                &buffer1(liminf:limsup,2*pn)*auxr2        
        else
           cor(xci+d994:xco-d994,pn)=cor(xci+d994:xco-d994,pn)-buffer1(liminf:limsup,2*pn-1)*auxr2+&
                &buffer1(liminf:limsup,2*pn)*auxr1        
        endif
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










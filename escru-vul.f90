! ----------------------------------------------------------------------! 
!
! Writing Subroutines Package
!
! ----------------------------------------------------------------------! 
!----------------------------------------------------------------------*
subroutine escribezy(u,v,w,p,dt,mpiid)
  !----------------------------------------------------------------------*
  ! escribe campos de u,v,w,p
  !----------------------------------------------------------------------*
  use point
  use names
  use genmod
  use alloc_dns
  use ctesp
  implicit none
  include 'mpif.h'
  real*8,dimension(nz1,ny+1,ib:ie)::u,w
  real*8,dimension(nz1,ny  ,ib:ie)::p,v
  real*4,dimension(:,:),allocatable::resu
  integer i,j,k,l,irec
  integer status(MPI_STATUS_SIZE),ierr,size,size1,size2,dot,mpiid,lim1,lim2
  character fil1*80,fil2*80,fil3*80,fil4*80,ext1*3,uchar*1
  real*8    dt,dum(20),jk,tiempoo
  integer:: nxr3,nyr3,nzr3,xoutf
  ! ------------------------- Program ----------------------------
  ! operaciones del modo maestro
  pi=4d0*atan(1d0)
  dum=0d0
  
  allocate (resu(nz1,ny+1))
  size  = nz1*(ny+1)*4
  size1 = nz1*(ny+1)
  size2 = nz1*(ny  )

  if (mpiid.eq.0) then

     write(*,*) 'before starting to write'
     write(*,*) '              x                   y                  um'
     write(*,*) '--------------------------------------------------------------------'
      
     do i=1,10
     write(*,'(3f20.6)') x(nx-i),y(ny+1-i),um(ny+1-i)
     enddo


     write(ext1,'(i3.3)') ifile
     write(*,*) 'ifile=',ifile
     fil1=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'u'
     open (10,file=fil1,status='unknown', &
          & form='unformatted',access='direct',recl=size1,convert='BIG_ENDIAN')

     fil2=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'v'
     open (11,file=fil2,status='unknown', &
          & form='unformatted',access='direct',recl=size2,convert='BIG_ENDIAN')

     fil3=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'w'
     open (12,file=fil3,status='unknown', &
          & form='unformatted',access='direct',recl=size1,convert='BIG_ENDIAN')

     fil4=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'p'
     open (13,file=fil4,status='unknown', &
          & form='unformatted',access='direct',recl=size2,convert='BIG_ENDIAN')

     write(*,*) '----- files OPEN ------'
     write(*,*) fil1
     write(*,*) fil2
     write(*,*) fil3
     write(*,*) fil4
     write(*,*) '..... writing the headers: ...................'

     irec=1
     write(*,*) 'u'
     write(*,'(8f12.4)') tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,timeinit,dt
     write(*,'(5i10)') nx,ny,nz2,xout,nummpi
     write(*,*) '------------------------------------------------------------------'
     write(10,rec=irec) 'u',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
          & dum, (x(i), i=0,nx+1), (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

     write(11,rec=irec) 'v',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
          & dum, (x(i), i=0,nx+1), (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

     write(12,rec=irec) 'w',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
          & dum, (x(i), i=0,nx+1), (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

    write(13,rec=irec)  'p',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
          & dum, (x(i), i=0,nx+1), (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

     do i=ib,ie
        if(mod(i,100).eq.0) write (*,*) 'escribiendo hasta plano',i,' de ', nx        
        irec = irec+1
        write(10,rec=irec) real(u(:,:,i),kind=4)
        write(11,rec=irec) real(v(:,:,i),kind=4)
        write(12,rec=irec) real(w(:,:,i),kind=4)
        write(13,rec=irec) real(p(:,:,i),kind=4)
        call flush(10,11,12,13)
     enddo

     do dot = 1,nummpi-1   !recibe la informacion de cada procesador
        do i= ibeg(dot),iend(dot)
          if(mod(i,100).eq.0) write (*,*) 'escribiendo hasta plano',i,' de ', nx  
          irec = irec+1     
          call MPI_RECV(resu,size1,MPI_real4,dot,1,MPI_COMM_WORLD,status,ierr)
          write(10,rec=irec) resu(:,:)
          call MPI_RECV(resu,size2,MPI_real4,dot,2,MPI_COMM_WORLD,status,ierr)            
          write(11,rec=irec) resu(:,1:ny)
          call MPI_RECV(resu,size1,MPI_real4,dot,3,MPI_COMM_WORLD,status,ierr)         
          write(12,rec=irec) resu(:,:)
          call MPI_RECV(resu,size2,MPI_real4,dot,4,MPI_COMM_WORLD,status,ierr)  
          write(13,rec=irec) resu(:,1:ny)
          call flush(10,11,12,13)
        enddo
     enddo
     close(10);close(11);close(12);close(13)
     ifile=ifile+1
  else
     !operaciones para el resto de los nodos ********************
     do i=ib,ie     
        call MPI_SEND(real(u(:,:,i),kind=4),size1,MPI_real4,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(real(v(:,:,i),kind=4),size2,MPI_real4,0,2,MPI_COMM_WORLD,ierr)
        call MPI_SEND(real(w(:,:,i),kind=4),size1,MPI_real4,0,3,MPI_COMM_WORLD,ierr)
        call MPI_SEND(real(p(:,:,i),kind=4),size2,MPI_real4,0,4,MPI_COMM_WORLD,ierr)
     enddo
  endif
  deallocate (resu)

!   if (mpiid.eq.0) then    
!      open (10,file=fil1,status='unknown', &
!           & form='unformatted',access='direct',recl=size1,convert='BIG_ENDIAN')
!      write(*,*)
!      write(*,*) 'Checking some values of the new file:'
!      write(*,*) 'size1===',size1
!      read(10,rec=1) uchar,tiempoo,jk,jk,jk,jk,jk,nxr3,nyr3,nzr3,xoutf,jk,jk,dum,(x(i), i=0,nx+1), (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
!      write(*,*) 'in file    ', uchar,tiempoo,nxr3,nyr3,nzr3,nummpi
!      write(*,*) '              x                   y                  um'
!      write(*,*) '--------------------------------------------------------------------'
!       
!      do i=1,10
!      write(*,'(3f20.6)') x(nx-i),y(ny+1-i),um(ny+1-i)
!      enddo
!   endif
  
end subroutine escribezy

! -------------------------------------------------------------------! 
! -------------------------------------------------------------------! 
! -------------------------------------------------------------------! 

subroutine escrst(res,ax,ay,az,cfl,tiempo,re,x,y,mpiid,ical)
  use names
  use statistics
  use point
  use ctesp
  use omp_lib
  use genmod
  !The needed variables of alloc_dns are passed as arguments  
  implicit none  
  include 'mpif.h'
  
  integer status(MPI_STATUS_SIZE),ierr,i,j,k,lim1,lim2,siz,dot
  integer mpiid,ical,ii,kk,ind
  real*8 x(0:nx+1),y(0:ny+1),res(ny,nx),ax,ay,az,cfl,tiempo,re
  character ext1*2,ext*3,corfilv*99, cffile*99
  real*8,allocatable,dimension(:,:)::wkn,wknp
  real*8,allocatable,dimension(:,:,:)::spe
  !No spectra and correlations yet
  
  ! ------------------------ codigo ------------------------------------! 
  totalcal=totalcal+ical
  call MPI_ALLREDUCE(cf2,cf,nx,MPI_real8,MPI_SUM,MPI_COMM_WORLD,ierr)

  if (mpiid .eq. 0) then
     write(ext,'(i3.3)') ifile
     write(ext1,'(i2.2)') indst     
     cffile =trim(chfile)//'-'//'cf'//'.'//ext
     write(*,*) 'writting Cf in',cffile
     open (26,file=cffile,status='unknown',form='unformatted',convert='BIG_ENDIAN')
     write(26) cf(1:nx)/counter   
     close(26)

     write(ext,'(i3.3)') ifile
     write(ext1,'(i2.2)') indst
     stfile =trim(chfile)//'.'//ext1//'.'//ext//'.st'
     etfile =trim(chfile)//'.'//ext1//'.'//ext//'.esp'
     open (29,file=stfile,status='unknown',form='unformatted')
     indst=indst+1
     
  end if

   ! master only things  
  allocate(wkn(ny,11),wknp(ny+1,4))
  if (mpiid.eq.0) write(29) 0d0
  

  ! Write Reynolds stresses statistics
  if (mpiid .eq. 0) then
     write(*,*) 'writing in ',stfile
     rewind(29)
     write(*,*) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical
     write(29) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical
     !write(29) (x(i), i=0,nx+1) No longer write x
     write(29) (y(i), i=0,ny+1)
     
     do i=ib,ie 
        wknp(:,1)=us (:,i)
        wkn (:,1)=vs (:,i)
        wknp(:,2)=ws (:,i)
        wkn (:,2)=uv (:,i)
        wknp(:,3)=ua (:,i)
        wkn (:,3)=va (:,i)
        wknp(:,4)=wa (:,i)
        wkn (:,4)=pp (:,i)
        wkn (:,5)=pm (:,i)
        !Vorticities
        wkn (:,6) =vortx(1:ny,i)
        wkn (:,7) =vorty(1:ny,i)
        wkn (:,8) =vortz(1:ny,i)
        wkn (:,9) =vortxa(1:ny,i)
        wkn (:,10)=vortya(1:ny,i)
        wkn (:,11)=vortza(1:ny,i)
        
        write(29) (wknp(j,1),j=1,ny+1),&
             &    (wkn (j,1),j=1,ny  ),&
             &    (wknp(j,2),j=1,ny+1),&
             &    (wkn (j,2),j=1,ny  ),&
             &    (wknp(j,3),j=1,ny+1),&
             &    (wkn (j,3),j=1,ny  ),&
             &    (wknp(j,4),j=1,ny+1),&
             &    (wkn (j,4),j=1,ny  ),&
             &    (wkn (j,5),j=1,ny  ),&
             !Vorticities
             &    (wkn (j,6),j=1,ny  ),&
             &    (wkn (j,7),j=1,ny  ),&
             &    (wkn (j,8),j=1,ny  ),&
             &    (wkn (j,9),j=1,ny  ),&
             &    (wkn (j,10),j=1,ny  ),&
             &    (wkn (j,11),j=1,ny  )
     enddo
     
     do dot = 1,nummpi-1
        do i=ibeg(dot),iend(dot)
           
           call MPI_RECV(wknp,4*(ny+1),MPI_real8,dot,0,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(wkn ,11*ny   ,MPI_real8,dot,0,MPI_COMM_WORLD,status,ierr)
           
           write(29) (wknp(j,1),j=1,ny+1),&
                &    (wkn (j,1),j=1,ny  ),&
                &    (wknp(j,2),j=1,ny+1),&
                &    (wkn (j,2),j=1,ny  ),&
                &    (wknp(j,3),j=1,ny+1),&
                &    (wkn (j,3),j=1,ny  ),&
                &    (wknp(j,4),j=1,ny+1),&
                &    (wkn (j,4),j=1,ny  ),&
                &    (wkn (j,5),j=1,ny  ),&
                !Vorticities             
                &    (wkn (j,6),j=1,ny  ),&
                &    (wkn (j,7),j=1,ny  ),&
                &    (wkn (j,8),j=1,ny  ),&
                &    (wkn (j,9),j=1,ny  ),&
                &    (wkn (j,10),j=1,ny  ),&
                &    (wkn (j,11),j=1,ny  )
        enddo
     enddo
     call flush(29)
     close(29)
  else
     do i=ib,ie
        wknp(:,1)=us (:,i)
        wkn (:,1)=vs (:,i)
        wknp(:,2)=ws (:,i)
        wkn (:,2)=uv (:,i)
        wknp(:,3)=ua (:,i)
        wkn (:,3)=va (:,i)
        wknp(:,4)=wa (:,i)
        wkn (:,4)=pp (:,i)
        wkn (:,5)=pm (:,i)
        !Vorticities
        wkn (:,6) =vortx(1:ny,i)
        wkn (:,7) =vorty(1:ny,i)
        wkn (:,8) =vortz(1:ny,i)
        wkn (:,9) =vortxa(1:ny,i)
        wkn (:,10)=vortya(1:ny,i)
        wkn (:,11)=vortza(1:ny,i)

        call MPI_SEND(wknp,4*(ny+1),MPI_real8,0,0,MPI_COMM_WORLD,ierr)
        call MPI_SEND(wkn ,11*ny   ,MPI_real8,0,0,MPI_COMM_WORLD,ierr)
     enddo         
  endif
  ! Initialize everything to zero
  
  pm   = 0d0
  pp   = 0d0  
  us   = 0d0
  ws   = 0d0
  ua   = 0d0
  wa   = 0d0
  uv   = 0d0
  vs   = 0d0
  va   = 0d0
  vortx = 0d0 
  vorty = 0d0
  vortz = 0d0
  vortxa= 0d0
  vortya= 0d0
  vortza= 0d0
  
  deallocate(wkn,wknp)

  ! Writing the spectra in Z.
  
  allocate(spe(0:nz2,nspec,4))

  siz = (nz2+1)*nspec*4 !Size is common for all nodes and is constant.

  if (mpiid == 0) then
     open(30,file=etfile,status='unknown',form='unformatted')
     rewind(30)
     write(*,*) 'writing in ',etfile
     write(30) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical
     !write(29) (x(i), i=0,nx+1) No longer write x
     write(30) (y(i), i=0,ny+1),((jspecy(i,j),i=1,nspec),j=1,lxp)
     
     do i = ib,ie
        spe(:,:,1) = ensu (:,:,i)
        spe(:,:,2) = ensv (:,:,i)
        spe(:,:,3) = ensw (:,:,i)
        spe(:,:,4) = ensuv(:,:,i)
        write(30) (((spe(k,j,ii),k=0,nz2),j=1,nspec),ii=1,4)
     end do

     do dot = 1,nummpi-1
        call MPI_RECV(spe,siz,MPI_REAL8,dot,0,MPI_COMM_WORLD,status,ierr)
        write(30) (((spe(k,j,ii),k=0,nz2),j=1,nspec),ii=1,4)
     end do

     call flush(30)
     close(30)

  else !Rest of nodes

     do i = ib,ie
        spe(:,:,1) = ensu (:,:,i)
        spe(:,:,2) = ensv (:,:,i)
        spe(:,:,3) = ensw (:,:,i)
        spe(:,:,4) = ensuv(:,:,i)
        call MPI_SEND(spe,siz,MPI_REAL8,0,0,MPI_COMM_WORLD,ierr)
     end do

  end if

  ensu =0d0
  ensv =0d0
  ensw =0d0
  ensuv=0d0
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  deallocate(spe)

end subroutine escrst


  ! ----------------------------------------------------------------------! 
  !
  ! Writing Subroutines Package
  !
  ! ----------------------------------------------------------------------! 
  !----------------------------------------------------------------------*
  
#define MB *1024*1024
#define MAXPE 64*1024
#define MAXCHARLEN 250
  
  subroutine escribezy_2(u,v,w,p,dt,mpiid,communicator)
    !----------------------------------------------------------------------*
    ! escribe campos de u,v,w,p
    !----------------------------------------------------------------------*
    use point_2
    use names_2
    use genmod_2
    use alloc_dns_2
    use ctesp_2
    implicit none
    include 'mpif.h'
    integer,intent(in):: communicator
    real*8,dimension(nz1,ny+1,ib:ie)::u,w
    real*8,dimension(nz1,ny  ,ib:ie)::p,v
    real*4,dimension(:,:,:),allocatable::resu
    integer i,j,k,l,irec
    integer status(MPI_STATUS_SIZE),ierr,t_size,t_size1,t_size2,dot,mpiid,lim1,lim2
    character(len=MAXCHARLEN):: fil1,fil2,fil3,fil4
    character:: ext1*3,uchar*1
    real*8    dt,dum(20),jk,t0
    integer:: nxr3,nyr3,nzr3,comm,tipo,chunkfbs,nfile,sidio,mpiw1,mpiw2,mpiw3,mpiw4
    integer*8:: chunks1,chunks2,chunksM1,chunksM2

    ! ------------------------- Program ----------------------------  
    pi=4d0*atan(1d0)
    dum=0d0
    comm=communicator
    tipo=MPI_real4

    write(ext1,'(i3.3)') ifile     
    fil1=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'u'
    fil2=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'v'
    fil3=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'w'
    fil4=chfile(1:index(chfile,' ')-1)//'.'//ext1//'.'//'p'


#ifdef WPARALLEL
    !PARALLEL WRITTER ==================================================================
    !First the header and last the field
    if (mpiid.eq.0) t0=MPI_Wtime()

    nfile=1			    !Number of files for parallel IO
    chunkfbs=2*1024*1024            !File block system 2Mb
    chunks1 =nz1*(ny+1)*(ie-ib+1)*4  !Number of bytes in R4 for LocalBuffer
    chunks2 =nz1*(ny  )*(ie-ib+1)*4  !Number of bytes in R4 for LocalBuffer
    chunksM1=nz1*(ny+1)*(ie-ib+1+1)*4  !Number of bytes in R4 for the Master Node
    chunksM2=nz1*(ny  )*(ie-ib+1+1)*4  !Number of bytes in R4 for the Master Node
    if(mpiid.eq.0) then      
       write(*,*) '-------------------------CHUNK (Mb)---------------------------------------------------'
       write(*,*) '              chunks1          chunks2         chunksM1          chunksM2         chunkfbs'
       write(*,'(5F18.3)') 1.0*chunks1/1024/1024,1.0*chunks2/1024/1024,1.0*chunksM1/1024/1024,1.0*chunksM2/1024/1024,1.0*chunkfbs/1024/1024
       write(*,*) '----------------------------------------------------------------------------------'
    endif

    if(mpiid.ne.0) then       
       allocate (resu(nz1,ny+1,ie-ib+1),stat=ierr);resu=0 !R4 buffer to convert R8 variables
       if(ierr.ne.0) write(*,*) "ERROR ALLOCATING RESU"       
       !Writting u:
       resu=real(u,kind=4)       
       call blockwrite_2(fil1,comm,resu,chunks1,nfile,mpiid,sidio)     
       !Writting v:
       resu(:,1:ny,:)=real(v,kind=4)         
       call blockwrite_2(fil2,comm,resu(:,1:ny,:),chunks2,nfile,mpiid,sidio)     
       !Writting w:
       resu=real(w,kind=4)   
       call blockwrite_2 (fil3,comm,resu,chunks1,nfile,mpiid,sidio)   
       !Writting p:
       resu(:,1:ny,:)=real(p,kind=4)    
       call blockwrite_2 (fil4,comm,resu(:,1:ny,:),chunks2,nfile,mpiid,sidio)      
       deallocate (resu)  
       ifile=ifile+1 
    else      
       allocate (resu(nz1,ny+1,(ie-ib+1)+1),stat=ierr);resu=0 !R4 buffer to convert R8 variables
       if(ierr.ne.0) write(*,*) "ERROR ALLOCATING RESU"
       write(*,*)
       write(*,'(a75,f10.4,a3)') 'Size of the allocated buffer in order to write:',size(resu)*4.0/1024/1024,'Mb'              
       !Writting u:
       resu(:,:,2:)=real(u,kind=4)
       call blockwrite_2 (fil1,comm,resu,chunksM1,nfile,mpiid,sidio)       
       call writeheader_2(fil1,'u',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,xout,timeinit,dt,y,um,nummpi)
       !Writting v:
       resu(:,1:ny,2:)=real(v,kind=4)  
       call blockwrite_2 (fil2,comm,resu(:,1:ny,:),chunksM2,nfile,mpiid,sidio)      
       call writeheader_2fil2,'v',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,xout,timeinit,dt,y,um,nummpi)
       !Writting w:
       resu(:,:,2:)=real(w,kind=4)           
       call blockwrite_2 (fil3,comm,resu,chunksM1,nfile,mpiid,sidio)      
       call writeheader_2(fil3,'w',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,xout,timeinit,dt,y,um,nummpi)
       !Writting p:
       resu(:,1:ny,2:)=real(p,kind=4)           
       call blockwrite_2 (fil4,comm,resu(:,1:ny,:),chunksM2,nfile,mpiid,sidio)      
       call writeheader_2(fil4,'p',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,xout,timeinit,dt,y,um,nummpi)
       deallocate (resu)  
       ifile=ifile+1        
    endif

    call MPI_BARRIER(comm,ierr)

    if (mpiid.eq.0) then 
       t0=MPI_Wtime()-t0
       write(*,*)
       write(*,*) '=========================================================================='
       write(*,*) 'Done writting', chfile(1:index(chfile,' ')-1)//'.'//ext1,' fields'
       write(*,*) '=========================================================================='    
       write(*,'(a20,f10.3,a3)') 'TIME SPENT IN WRITING:',t0,'sc'
       write(*,*)   '--------------------------------------------------------------------------'
    endif
#endif



#ifdef WSERIAL
    !SERIAL WRITTER ==================================================================
    allocate (resu(nz1,ny+1,1))
    t_size  = nz1*(ny+1)*4
    t_size1 = nz1*(ny+1)
    t_size2 = nz1*(ny  )

    if (mpiid.eq.0) then
       write(*,*) 'before starting to write'
       write(*,*) '              x                   y                  um'
       write(*,*) '--------------------------------------------------------------------'      
       do i=1,10
          write(*,'(3f20.6)') x(nx-i),y(ny+1-i),um(ny+1-i)
       enddo

       open (40,file=fil1,status='unknown', &
            & form='unformatted',access='direct',recl=t_size1*4,convert='Big_endian')     
       open (41,file=fil2,status='unknown', &
            & form='unformatted',access='direct',recl=t_size2*4,convert='Big_endian')    
       open (42,file=fil3,status='unknown', &
            & form='unformatted',access='direct',recl=t_size1*4,convert='Big_endian')    
       open (43,file=fil4,status='unknown', &
            & form='unformatted',access='direct',recl=t_size2*4,convert='Big_endian')

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

       write(40,rec=irec) 'u',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       write(41,rec=irec) 'v',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       write(42,rec=irec) 'w',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       write(43,rec=irec)  'p',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

       do i=ib,ie
          if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx
          irec = i+1
          write(40,rec=irec) real(u(:,:,i),kind=4)
          write(41,rec=irec) real(v(:,:,i),kind=4)
          write(42,rec=irec) real(w(:,:,i),kind=4)
          write(43,rec=irec) real(p(:,:,i),kind=4)
          call flush(40,41,42,43)
       enddo

       do dot = 1,nummpi-1   !recibe la informacion de cada procesador
          do i= ibeg(dot),iend(dot)
             if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx
             irec = i+1     
             call MPI_RECV(resu,t_size1,tipo,dot,1,comm,status,ierr)
             write(40,rec=irec) resu(:,:,1)
             call MPI_RECV(resu,t_size2,tipo,dot,2,comm,status,ierr)
             write(41,rec=irec) resu(:,1:ny,1)
             call MPI_RECV(resu,t_size1,tipo,dot,3,comm,status,ierr)         
             write(42,rec=irec) resu(:,:,1)
             call MPI_RECV(resu,t_size2,tipo,dot,4,comm,status,ierr)  
             write(43,rec=irec) resu(:,1:ny,1)
             call flush(40,41,42,43)
          enddo
       enddo
       close(40);close(41);close(42);close(43)
       ifile=ifile+1
    else
       !operaciones para el resto de los nodos ********************
       do i=ib,ie     
          call MPI_SEND(real(u(:,:,i),kind=4),t_size1,tipo,0,1,comm,ierr)
          call MPI_SEND(real(v(:,:,i),kind=4),t_size2,tipo,0,2,comm,ierr)
          call MPI_SEND(real(w(:,:,i),kind=4),t_size1,tipo,0,3,comm,ierr)
          call MPI_SEND(real(p(:,:,i),kind=4),t_size2,tipo,0,4,comm,ierr)
       enddo
    endif
    deallocate (resu)  
    !END SERIAL WRITTER ===============================================================
#endif


#ifdef WSERIAL4
    !SERIAL WRITTER WITH 4 NODES==========================================================
    !MPI Writers nodes:
    mpiw1=0
    mpiw2=nummpi/4
    mpiw3=nummpi/2
    mpiw4=3*nummpi/4
    if (mpiid.eq.mpiw1) t0=MPI_Wtime()
    allocate (resu(nz1,ny+1,1))
    t_size=nz1*(ny+1)*4
    t_size1=nz1*(ny+1)
    t_size2=nz1*(ny  )

    if (mpiid.eq.mpiw1) then

       write(*,*) 'before starting to write'
       write(*,*) '              x                   y                  um'
       write(*,*) '--------------------------------------------------------------------'      
       do i=1,10
          write(*,'(3f20.6)') x(nx-i),y(ny+1-i),um(ny+1-i)
       enddo

       open (40,file=fil1,status='unknown', &
            & form='unformatted',access='direct',recl=t_size1*4,convert='Big_endian') 
       write(*,*) '----- files OPEN ------'
       write(*,*) fil1
       irec=1
       write(*,*) 'u'
       write(*,'(8f12.4)') tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,timeinit,dt
       write(*,'(5i10)') nx,ny,nz2,xout,nummpi
       write(*,*) '------------------------------------------------------------------'

       write(40,rec=irec) 'u',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi

       do i=ib,ie             
          irec = i+1
          write(40,rec=irec) real(u(:,:,i),kind=4)      
          call flush(40)
       enddo

       do dot = 1,nummpi-1   !recibe la informacion de cada procesador
          do i= ibeg(dot),iend(dot)
             if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx, '....FILE U'         
             irec = i+1     
             call MPI_RECV(resu,t_size1,tipo,dot,1,comm,status,ierr)
             write(40,rec=irec) resu(:,:,1)         
             call flush(40)
          enddo
       enddo
       close(40)
       ifile=ifile+1
    else
       !operaciones para el resto de los nodos ********************
       do i=ib,ie     
          call MPI_SEND(real(u(:,:,i),kind=4),t_size1,tipo,mpiw1,1,comm,ierr)      
       enddo
    endif

    !-----------------------------------------------
    if (mpiid.eq.mpiw2) then
       open (11,file=fil2,status='unknown', &
            & form='unformatted',access='direct',recl=t_size2*4,convert='Big_endian')    
       write(*,*) fil2
       irec=1
       write(41,rec=irec) 'v',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       do i=ib,ie               
          irec = i+1      
          write(41,rec=irec) real(v(:,:,i),kind=4)      
          call flush(41)
       enddo
       do dot = 0,nummpi-1   !recibe la informacion de cada procesador
          if(dot.ne.mpiw2) then
             do i= ibeg(dot),iend(dot)
                if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx, '....FILE V'         
                irec = i+1
                call MPI_RECV(resu,t_size2,tipo,dot,2,comm,status,ierr)
                write(41,rec=irec) resu(:,1:ny,1)       
                call flush(41)
             enddo
          endif
       enddo
       close(41)
       ifile=ifile+1
    else
       !operaciones para el resto de los nodos ********************
       do i=ib,ie        
          call MPI_SEND(real(v(:,:,i),kind=4),t_size2,tipo,mpiw2,2,comm,ierr)      
       enddo
    endif

    !-----------------------------------------------
    if (mpiid.eq.mpiw3) then
       open (12,file=fil3,status='unknown', &
            & form='unformatted',access='direct',recl=t_size1*4,convert='Big_endian')  
       write(*,*) fil3  
       irec=1
       write(42,rec=irec) 'w',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       do i=ib,ie
          irec = i+1       
          write(42,rec=irec) real(w(:,:,i),kind=4)       
          call flush(12)
       enddo
       do dot = 0,nummpi-1   !recibe la informacion de cada procesador
          if(dot.ne.mpiw3) then
             do i= ibeg(dot),iend(dot)
                if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx, '....FILE W'         
                irec = i+1             
                call MPI_RECV(resu,t_size1,tipo,dot,3,comm,status,ierr)
                write(42,rec=irec) resu(:,:,1)
                call flush(42)
             enddo
          endif
       enddo
       close(42)
       ifile=ifile+1
    else
       !operaciones para el resto de los nodos ********************
       do i=ib,ie            
          call MPI_SEND(real(w(:,:,i),kind=4),t_size1,tipo,mpiw3,3,comm,ierr)       
       enddo
    endif

    !-----------------------------------------------
    if (mpiid.eq.mpiw4) then
       open (43,file=fil4,status='unknown', &
            & form='unformatted',access='direct',recl=t_size2*4,convert='Big_endian')
       write(*,*) fil4      
       irec=1
       write(43,rec=irec)  'p',tiempo,cfl,Re,ax*pi,ay*pi,2*pi*az,nx,ny,nz2,xout,timeinit,dt, &
            & (y(i), i=0,ny+1), (um(i), i=1,ny+1),nummpi
       do i=ib,ie
          irec = i+1      
          write(43,rec=irec) real(p(:,:,i),kind=4)
          call flush(43)
       enddo
       do dot = 0,nummpi-1   !recibe la informacion de cada procesador
          if(dot.ne.mpiw4) then
             do i= ibeg(dot),iend(dot)
                if(mod(i,500).eq.0) write (*,*) 'Writting plane #',i,' of ', nx, '....FILE P'         
                irec = i+1            
                call MPI_RECV(resu,t_size2,tipo,dot,4,comm,status,ierr)  
                write(43,rec=irec) resu(:,1:ny,1)
                call flush(43)
             enddo
          endif
       enddo
       close(43)
       ifile=ifile+1
    else
       !operaciones para el resto de los nodos ********************
       do i=ib,ie           
          call MPI_SEND(real(p(:,:,i),kind=4),t_size2,tipo,mpiw4,4,comm,ierr)
       enddo
    endif

    !-----------------------------------------------
    call MPI_BARRIER(comm,ierr) 
    if (mpiid.eq.mpiw1) then 
       t0=MPI_Wtime()-t0
       write(*,*)
       write(*,*) '=========================================================================='
       write(*,*) 'Done writing', chfile(1:index(chfile,' ')-1)//'.'//ext1,' fields'
       write(*,*) '=========================================================================='
       write(*,'(a20,f10.3,a3)') 'TIME SPENT IN WRITING:',t0,'sc'
       write(*,*)   '--------------------------------------------------------------------------'
    endif
    deallocate (resu)  
#endif
  end subroutine escribezy_2



  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 

  subroutine escrst_2(ax,ay,az,cfl,tiempo,re,x,y,mpiid,ical,communicator)
    use names_2
    use statistics_2
    use point_2
    use ctesp_2
    use omp_lib
    use genmod_2
    !The needed variables of alloc_dns are passed as arguments  
    implicit none  
    include 'mpif.h'
    integer,intent(in):: communicator

    integer status(MPI_STATUS_SIZE),ierr,i,j,k,lim1,lim2,t_size,siz,dot
    integer mpiid,ii,kk,ind,comm,tipo,elements_spec,nbud
    real*8 x(0:nx+1),y(0:ny+1),ax,ay,az,cfl,tiempo,re,xxx    
    character ext1*2,ext*3,corfilv*99, cffile*99
    real*8,allocatable,dimension(:,:)    ::wkn,wknp
    real*8,allocatable,dimension(:,:,:)  ::buf_spe  
    real*8,allocatable,dimension(:,:  )  ::buf_bud,buf_bud2
#ifdef PLANESPECTRA
    integer:: ical(7) 
#else
    integer:: ical
#endif  
    comm=communicator
    tipo=MPI_real8
    ! ------------------------ codigo ------------------------------------! 
    totalcal=totalcal+ical
    if (mpiid.eq.0) then
       write(ext,'(i3.3)') ifile
       write(ext1,'(i2.2)') indst     
      
       stfile =trim(chfile)//'.'//ext1//'.'//ext//'.st'
       etfile =trim(chfile)//'.'//ext1//'.'//ext//'.esp'
       corfile=trim(chfile)//'.'//ext1//'.'//ext//'.cor'
       budfile=trim(chfile)//'.'//ext1//'.'//ext//'.budget'
       spectraplane=trim(chfile)//'.'//ext1//'.'//ext//'.extraesp'

       open(39,file=stfile,status='unknown',form='unformatted',convert='Big_endian');rewind(39)
       indst=indst!+1 !WE WRITE STATISTICS WHEN RECORD IMAGES ONLY!!
       write(*,*) 'writing in==============  ',stfile       
       write(*,'(6f12.4,4i10)') tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical
       
       write(39) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical       
       write(39) (y(i), i=0,ny+1)
    end if
    ! master only things  
    allocate(wkn(ny,11),wknp(ny+1,4))
    if (mpiid.eq.0) write(39) 0d0

    ! things in ib:ie  
    if (mpiid .eq. 0) then
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

          write(39) (wknp(j,1),j=1,ny+1),&
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
             call MPI_RECV(wknp,4*(ny+1),tipo,dot,0,comm,status,ierr)
             call MPI_RECV(wkn ,11*ny   ,tipo,dot,0,comm,status,ierr)

             write(39) (wknp(j,1),j=1,ny+1),&
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
       call flush(39)
       close(39)
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

          call MPI_SEND(wknp,4*(ny+1),tipo,0,0,comm,ierr)
          call MPI_SEND(wkn ,11*ny   ,tipo,0,0,comm,ierr)
       enddo
    endif

    ! Initialize everything to zero
    pp=0d0  !pm,ua,va must be initialized later...I'll used it for budgets
    us=0d0;ws=0d0;wa=0d0;uv=0d0;vs=0d0
    vortx=0d0;vorty=0d0;vortz=0d0;vortxa= 0d0;vortya= 0d0;vortza= 0d0

    deallocate(wkn,wknp) 

#ifndef NOSPECTRA
    elements_spec=lxp*(nz2+1)*nspec
    !=================== SPECTRA =======================
    !elements_spec=lxp*(nz2+1)*nspec. ens=ens*2*nxp (divide it when writing to disk)
    !VELOCITIES:
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensu  ,elements_spec,tipo,MPI_SUM,comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensv  ,elements_spec,tipo,MPI_SUM,comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensw  ,elements_spec,tipo,MPI_SUM,comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensuv ,elements_spec,tipo,MPI_SUM,comm,ierr)
    !VORTICITIES:
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensomz,elements_spec,tipo,MPI_SUM,comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensomx,elements_spec,tipo,MPI_SUM,comm,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ensomy,elements_spec,tipo,MPI_SUM,comm,ierr)
    !PRESSURE:
    call MPI_ALLREDUCE(MPI_IN_PLACE,pesp  ,elements_spec,tipo,MPI_SUM,comm,ierr)
   


    ! Writing the spectra in Z.  
    allocate(buf_spe(0:nz2,nspec,8))    
    if (mpiid.eq.0) then
       open(83,file=etfile,status='unknown',form='unformatted',convert='Big_endian');rewind(83)
       write(*,*) 'writing in==============  ',etfile
       write(83) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical,nspec,lxp,nxp(1:lxp),xcorpoint(1:lxcorr)            
       write(83) (y(i), i=0,ny+1),((jspecy(i,j),i=1,nspec),j=1,lxp)

       do i=1,lxp             
          buf_spe(:,:,1)=ensu  (:,:,i)  
          buf_spe(:,:,2)=ensv  (:,:,i)
          buf_spe(:,:,3)=ensw  (:,:,i)
          buf_spe(:,:,4)=ensuv (:,:,i)
          buf_spe(:,:,5)=ensomz(:,:,i)
          buf_spe(:,:,6)=ensomx(:,:,i)
          buf_spe(:,:,7)=ensomy(:,:,i)
          buf_spe(:,:,8)=pesp  (:,:,i)
          buf_spe=buf_spe/(2d0*nxp(i)+1) !Averaging the buf_spectra in nxp
          write(83) buf_spe(0:nz2,1:nspec,1:8)
       end do       
       close(83)    
    end if
    ensu =0d0;ensv =0d0;ensw =0d0;ensuv=0d0;ensomz =0d0;ensomx =0d0;ensomy =0d0;pesp=0d0   
    deallocate(buf_spe)  
#endif

#ifndef NODISSIPATION
    !=================== BUDGETS =======================
    if (mpiid.eq.0) then
       open(55,file=budfile,status='unknown',form='unformatted',convert='Big_endian');rewind(55)
       write(*,*) 'writing in==============',budfile
       !HEADER
       write(55) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical    
       write(55) (y(i), i=0,ny+1)
    end if

#ifdef INFOINTER 
   nbud=27+9
#else    
   nbud=27
#endif

    allocate(buf_bud(ny,nbud),buf_bud2(ny+1,2))
    if (mpiid.eq.0) then
       write(*,*) 'Writing a total of',nbud,'budgets'
       do i=ib,ie   
          !NY size Buffers            
          buf_bud(:,1) =dispu (1:ny,i)
          buf_bud(:,2) =dispv (1:ny,i)
          buf_bud(:,3) =dispw (1:ny,i)
          buf_bud(:,4) =dispuv(1:ny,i)
          buf_bud(:,5) =pup   (1:ny,i)
          buf_bud(:,6) =pvp   (1:ny,i)
          buf_bud(:,7) =pdudx (1:ny,i)
          buf_bud(:,8) =pdvdy (1:ny,i)
          buf_bud(:,9) =pdudy (1:ny,i)
          buf_bud(:,10)=pdvdx (1:ny,i)
          buf_bud(:,11)=pdwdz (1:ny,i)
          buf_bud(:,12)=v3    (1:ny,i)
          buf_bud(:,13)=u2v   (1:ny,i)
          buf_bud(:,14)=v2u   (1:ny,i)
          buf_bud(:,15)=w2v   (1:ny,i)
          !K=0 MODE VARIABLES
          buf_bud(:,16)=dudx0 (1:ny,i)
          buf_bud(:,17)=dudy0 (1:ny,i)
          buf_bud(:,18)=dudz0 (1:ny,i)
          buf_bud(:,19)=dvdx0 (1:ny,i)
          buf_bud(:,20)=dvdy0 (1:ny,i)
          buf_bud(:,21)=dvdz0 (1:ny,i)
          buf_bud(:,22)=dwdx0 (1:ny,i)
          buf_bud(:,23)=dwdy0 (1:ny,i)
          buf_bud(:,24)=dwdz0 (1:ny,i)
          buf_bud(:,25)=pm    (1:ny,i)
          buf_bud(:,26)=ua    (1:ny,i)
          buf_bud(:,27)=va    (1:ny,i)
#ifdef INFOINTER 
         buf_bud(:,28)=v_0    (1:ny,i)
         buf_bud(:,29)=u_x0   (1:ny,i)
         buf_bud(:,30)=u_xy0  (1:ny,i)
         buf_bud(:,31)=w_0    (1:ny,i)
         buf_bud(:,32)=w_y0   (1:ny,i)
         buf_bud(:,33)=dwdx_0 (1:ny,i)
         buf_bud(:,34)=dudz_x0(1:ny,i) 
         buf_bud(:,35)=v_y0   (1:ny,i)         
         buf_bud(:,36)=dudx_0 (1:ny,i)            
#endif          
          !NY+1 size Buffers
          buf_bud2(:,1)=u3    (1:ny+1,i)
          buf_bud2(:,2)=w2u   (1:ny+1,i)                                                           
          write(55) buf_bud(1:ny,1:nbud),buf_bud2(1:ny+1,1:2)  
       enddo
       do dot = 1,nummpi-1
          do i=ibeg(dot),iend(dot)          
             call MPI_RECV(buf_bud, size(buf_bud), tipo,dot,0,comm,status,ierr)
             call MPI_RECV(buf_bud2,size(buf_bud2),tipo,dot,1,comm,status,ierr)
             if(mod(i,500).eq.0) write(*,*) 'writing budget',i,'of',nx
             write(55) buf_bud(1:ny,1:nbud),buf_bud2(1:ny+1,1:2)                           
          enddo
       enddo
       call flush(55)
       close(55)   
    else
       do i=ib,ie
          !NY size Buffers            
          buf_bud(:,1) =dispu (1:ny,i)
          buf_bud(:,2) =dispv (1:ny,i)
          buf_bud(:,3) =dispw (1:ny,i)
          buf_bud(:,4) =dispuv(1:ny,i)
          buf_bud(:,5) =pup   (1:ny,i)
          buf_bud(:,6) =pvp   (1:ny,i)
          buf_bud(:,7) =pdudx (1:ny,i)
          buf_bud(:,8) =pdvdy (1:ny,i)
          buf_bud(:,9) =pdudy (1:ny,i)
          buf_bud(:,10)=pdvdx (1:ny,i)
          buf_bud(:,11)=pdwdz (1:ny,i)
          buf_bud(:,12)=v3    (1:ny,i)
          buf_bud(:,13)=u2v   (1:ny,i)
          buf_bud(:,14)=v2u   (1:ny,i)
          buf_bud(:,15)=w2v   (1:ny,i)
          !K=0 MODE VARIABLES
          buf_bud(:,16)=dudx0 (1:ny,i)
          buf_bud(:,17)=dudy0 (1:ny,i)
          buf_bud(:,18)=dudz0 (1:ny,i)
          buf_bud(:,19)=dvdx0 (1:ny,i)
          buf_bud(:,20)=dvdy0 (1:ny,i)
          buf_bud(:,21)=dvdz0 (1:ny,i)
          buf_bud(:,22)=dwdx0 (1:ny,i)
          buf_bud(:,23)=dwdy0 (1:ny,i)
          buf_bud(:,24)=dwdz0 (1:ny,i)
          buf_bud(:,25)=pm    (1:ny,i)
          buf_bud(:,26)=ua    (1:ny,i)
          buf_bud(:,27)=va    (1:ny,i)
#ifdef INFOINTER 
         buf_bud(:,28)=v_0    (1:ny,i)
         buf_bud(:,29)=u_x0   (1:ny,i)
         buf_bud(:,30)=u_xy0  (1:ny,i)
         buf_bud(:,31)=w_0    (1:ny,i)
         buf_bud(:,32)=w_y0   (1:ny,i)
         buf_bud(:,33)=dwdx_0 (1:ny,i)
         buf_bud(:,34)=dudz_x0(1:ny,i) 
         buf_bud(:,35)=v_y0   (1:ny,i)         
         buf_bud(:,36)=dudx_0 (1:ny,i)            
#endif           
          !NY+1 size Buffers
          buf_bud2(:,1)=u3    (1:ny+1,i)
          buf_bud2(:,2)=w2u   (1:ny+1,i)                                                          
          call MPI_SEND(buf_bud ,size(buf_bud) ,tipo,0,0,comm,ierr)
          call MPI_SEND(buf_bud2,size(buf_bud2),tipo,0,1,comm,ierr)
       enddo
    endif
    ! Initialize everything to zero
    pm=0d0;ua=0d0;va=0d0 !must be intialized here!!
    dispu=0d0;dispv=0d0;dispw=0d0;dispuv=0d0;
    pvp=0d0;pup=0d0;pdudx=0d0;pdudy=0d0;pdvdx=0d0
    pdvdy=0d0;pdwdz=0d0
    u3=0d0;v3=0d0;u2v=0d0;v2u=0d0;w2v=0d0;w2u=0d0
    dudx0=0d0;dudy0=0d0;dudz0=0d0;
    dvdx0=0d0;dvdy0=0d0;dvdz0=0d0;
    dwdx0=0d0;dwdy0=0d0;dwdz0=0d0; 

#ifdef INFOINTER 
    v_0=0d0;u_x0=0d0;u_xy0=0d0;w_0=0d0;w_y0=0d0;
    dwdx_0=0d0;dudz_x0=0d0;v_y0=0d0;dudx_0=0d0
#endif 
    deallocate(buf_bud,buf_bud2)
#endif

#ifndef NOCORR
    !===================CORRELATIONS=======================    
    if (mpiid.eq.0) then
       open(47,file=corfile,status='unknown',form='unformatted',convert='Big_endian');rewind(47)
       write(*,*) 'writing in==============  ',corfile
       !HEADER
       write(47) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical,ncorr,lxcorr,nxp(1:lxcorr),xcorpoint(1:lxcorr)      
       write(47) y(0:ny+1),jspecy(1:ncorr,1:lxcorr)
       call flush(47)     
    end if 
  
    call escr_corr_2(coru,corv,coruv,corw,corp,corox,coroy,coroz,mpiid,communicator)  

    ! Initialize everything to zero
    coru=0d0;corv=0d0;corw=0d0;coruv=0d0;
    corox=0d0;coroy=0d0;coroz=0d0;corp=0d0;    
    close(47) 
#endif

#ifdef PLANESPECTRA
     if (mpiid.eq.0) then
       open(49,file=spectraplane,status='unknown',form='unformatted',convert='Big_endian');rewind(49)
       write(*,*) 'writing in==============  ',spectraplane
       write(49) ical(1:7)
     endif

     if (mpiid.eq.0) then
        do i=ib,ie
           write(49) plane_specu(0:nz2,1:7,i),plane_specv(0:nz2,1:7,i),plane_specw(0:nz2,1:7,i)
        enddo
        
        do dot = 1,nummpi-1
          do i=ibeg(dot),iend(dot)    
             call MPI_RECV(plane_specu,(nz2+1)*7,tipo,dot,0,comm,status,ierr)
             call MPI_RECV(plane_specv,(nz2+1)*7,tipo,dot,1,comm,status,ierr)
             call MPI_RECV(plane_specw,(nz2+1)*7,tipo,dot,2,comm,status,ierr)
             write(49) plane_specu(0:nz2,1:7,1),plane_specv(0:nz2,1:7,1),plane_specw(0:nz2,1:7,1)
          enddo
       enddo
       call flush(49)
       close(49)
     else
       do i=ib,ie
       call MPI_SEND(plane_specu(0,1,i),(nz2+1)*7,tipo,0,0,comm,ierr)
       call MPI_SEND(plane_specv(0,1,i),(nz2+1)*7,tipo,0,1,comm,ierr)
       call MPI_SEND(plane_specw(0,1,i),(nz2+1)*7,tipo,0,2,comm,ierr)
       enddo          
     endif
     plane_specu=0d0;plane_specv=0d0;plane_specw=0d0;
#endif  



#ifdef PLANESPECTRA2
tipo=MPI_DOUBLE_COMPLEX

     if (mpiid.eq.0) then
       open(44,file=spectraplane,status='unknown',form='unformatted',convert='Big_endian');rewind(44)
       write(*,*) 'writing in (PLANESPECTRA2)==============  ',spectraplane
       write(44) ss
     endif

     if (mpiid.eq.0) then
        do i=ib,ie
           write(44) planesv(0:nz2,1:500,i)
        enddo
        
        do dot = 1,nummpi-1
          do i=ibeg(dot),iend(dot)    
             call MPI_RECV(planesv,(nz2+1)*500,tipo,dot,0,comm,status,ierr)            
             write(44) planesv(0:nz2,1:500,1)
          enddo
       enddo
       call flush(44)
       close(44)
     else
       do i=ib,ie
          call MPI_SEND(planesv(0,1,i),(nz2+1)*500,tipo,0,0,comm,ierr)    
       enddo          
     endif
     planesv=0d0;
#endif  









  end subroutine escrst_2

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------


! call escr_corr(coru,corv,coruv,corw,corp,corox,coroy,coroz,mpiid)
subroutine escr_corr_2(c1,c2,c3,c4,c5,c6,c7,c8,rank,communicator)
use ctesp_2,only:nx,nz2,ncorr,lxcorr,nummpi,nxp,ny
use point_2
implicit none
include 'mpif.h'
integer,intent(in):: communicator
real*8,allocatable,dimension(:,:)::buf_cor 
real*8,dimension(nx,pcib2:pcie2,lxcorr)::c1,c2,c3,c4,c5,c6,c7,c8
integer:: rank,ierr,comm,status(MPI_STATUS_SIZE),i,j,dot,tipo,k,l

 comm = communicator
 tipo = MPI_real8
  
allocate(buf_cor(1:nx,8)) !8 Correlations

  do j=1,lxcorr
   if (rank.eq.0) write(*,*) 'writing correlations',j,'of',lxcorr
    if (rank.eq.0) then
       do i=pcib2,pcie2
          buf_cor(:,1)=c1(1:nx,i,j)       
          buf_cor(:,2)=c2(1:nx,i,j)
          buf_cor(:,3)=c3(1:nx,i,j)
          buf_cor(:,4)=c4(1:nx,i,j)
          buf_cor(:,5)=c5(1:nx,i,j)               
          buf_cor(:,6)=c6(1:nx,i,j)
          buf_cor(:,7)=c7(1:nx,i,j)
          buf_cor(:,8)=c8(1:nx,i,j)          
          buf_cor=buf_cor/(2d0*nxp(j)+1d0) !averaging         
          write(40) ((buf_cor(k,l),k=1,nx),l=1,8)  
       enddo

       do dot = 1,nummpi-1
          do i=pcibeg2(dot),pciend2(dot)          
             call MPI_RECV(buf_cor,8*nx,tipo,dot,j,comm,status,ierr)
             write(47) ((buf_cor(k,l),k=1,nx),l=1,8)                                 
          enddo
          call flush(47)        
       enddo          
    else   
       do i=pcib2,pcie2
          buf_cor(:,1)=c1(1:nx,i,j)       
          buf_cor(:,2)=c2(1:nx,i,j)
          buf_cor(:,3)=c3(1:nx,i,j)
          buf_cor(:,4)=c4(1:nx,i,j)
          buf_cor(:,5)=c5(1:nx,i,j)               
          buf_cor(:,6)=c6(1:nx,i,j)
          buf_cor(:,7)=c7(1:nx,i,j)
          buf_cor(:,8)=c8(1:nx,i,j)                
          buf_cor=buf_cor/(2d0*nxp(j)+1d0) !averaging        
          call MPI_SEND(buf_cor,8*nx,tipo,0,j,comm,ierr)
       enddo    
    endif
 enddo 
 deallocate(buf_cor)
endsubroutine escr_corr_2

#ifdef WPARALLEL
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! ------------ PARALLEL WRITTING SUBROUTINES ------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  subroutine blockwrite_2(filename,comm,localbuffer,chunksize,&
       & nfiles,rank,sid)

    !  Write a buffer localbuffer to a single file concurrently using all
    !  the MPI processes
    ! Input arguments:
    !  
    !  filename: String. Name of the file
    !  comm: MPI communicator
    !  localbuffer: Buffer to be written
    !  chunksize: Amount of **bytes** written from localbuffer.  Please, do not
    !             play weird games and use the same size of localbuffer
    !  fsblksize: File system block size.  GPFS is 2 MB
    !  nfiles: Put a variable that contains 1 here.  Not the literal.  Writing
    !          files is not supported yet.
    !  rank: Rank of the MPI process
    !  sid: File id, different from the OS file and obtained from the parallel
    !       opening process

#ifndef BG
    use mpi
#endif  
    implicit none
#ifdef BG
    include 'mpif.h'
#endif
    character(len=MAXCHARLEN),intent(in):: filename
    integer,intent(in):: comm
    integer*8,intent(in):: chunksize
    character,dimension(chunksize),intent(in):: localbuffer
    integer,intent(in):: nfiles
    integer,intent(in):: rank
    integer,intent(out):: sid
    integer,parameter:: fbsize=2*1024*1024
    character(len=MAXCHARLEN) :: newfname='newfile'
    integer:: lcomm,ierr,rankl

    integer*8:: checksum_fp,left,bsumwrote,chunkcnt
    integer*8:: bwrite,bwrote,sumsize

#ifdef TIMER
    real*8:: starttime,gstarttime,opentime
    real*8:: writetime,gwritetime,closetime
    real*8:: barr1time,barr2time,barr3time
#endif

#ifdef TIMER
    starttime=MPI_Wtime()
#endif
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'OPENING FSION_PARAOPEN-----------------------------------',rank   
    call fsion_paropen_mpi(trim(filename),'bw',nfiles, comm,&
         & lcomm,chunksize,fbsize,rank,newfname,sid)
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'OPENING FSION_PARAOPEN	done',rank       
#ifdef TIMER
    opentime = MPI_Wtime() - starttime
#endif
    call MPI_COMM_RANK(lcomm, rankl, ierr)
#ifdef TIMER
    starttime = MPI_Wtime()
    call barrier_after_open(lcomm)
    barr1time = MPI_Wtime()-starttime
#endif
    checksum_fp = 0
    left = chunksize
    bsumwrote = 0
    chunkcnt = 0
#ifdef TIMER
    starttime = MPI_Wtime()
    gstarttime = starttime
#endif
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'bucle-----------------------------------',rank   
    ! Fortran 90 specific!
    do while(left > 0)
       bwrite = chunksize
       if (bwrite > left) bwrite = left

       call fsion_ensure_free_space(sid,bwrite,ierr)
       call fsion_write(localbuffer, 1, bwrite, sid, bwrote)

#ifdef CHECKSUM
       do i=1,bwrote
          checksum_fp = checksum_fp + real(IACHAR(localbuffer(i)))
       end do
#endif
       left = left - bwrote
       bsumwrote = bsumwrote + bwrote
       chunkcnt = chunkcnt + 1
    end do
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'bucle		done',rank  
#ifdef TIMER
    writetime = MPI_Wtime() - starttime
    starttime = MPI_Wtime()
    call barrier_after_write(lcomm)
    barr2time = MPI_Wtime() - starttime
    gwritetime = MPI_Wtime() - gstarttime
    starttime = MPI_Wtime()
#endif
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'fsion_parclose-----------------------------------',rank  
    call fsion_parclose_mpi(sid,ierr)
    ! if(rank.eq.0.or.rank.eq.10) write(*,*) 'fsion_parclose		done',rank 
#ifdef TIMER
    call barrier_after_close(lcomm)
    barr3time = MPI_Wtime()-starttime
    closetime = MPI_Wtime() - starttime
    starttime = MPI_Wtime()


    if (writetime == 0) writetime = -1
#endif
    call MPI_REDUCE(bsumwrote, sumsize, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
    call MPI_BARRIER(comm,ierr)
#ifdef TIMER
    if (rank == 0) then       
       write(*,'(A)') "-----------------------------------------------------------------------"
       write(*,*) 'File written:',trim(filename)
       write(*,'(a20,f10.4,a3)') 'File Size:',1.0*sumsize/1024/1024/1024,'Gb'
       write(*,'(a20,f10.4)') 'T.Time Master Node:',gwritetime
       write(*,'(a20,f10.4,a7)') 'BandWidth:',1.0*sumsize/1024/1024/1024/gwritetime,'Gb/sec'                 
    end if
#endif
  end subroutine blockwrite_2


  subroutine writeheader_2(filename,field,tiempo,cfl,re,lx,ly,lz,nx,ny,nz2,&
       & xout,timein,dt,y,um,procs)
    implicit none
    character(len = MAXCHARLEN), intent(in):: filename
    character(len=1),intent(in):: field
    real(kind = 8), intent(in):: tiempo,cfl,re  !after 8k nods..make it R8 
    real(kind = 8), intent(in):: lx,ly,lz  
    integer, intent(in):: nx,ny,nz2
    integer, intent(in):: xout
    integer,parameter:: fbsize=2*1024*1024
    real(kind = 8),intent(in):: timein,dt
    real(kind = 8),dimension(ny+1),intent(in):: um
    real(kind = 8),dimension(0:ny+1),intent(in)::y
    integer,intent(in):: procs
    integer*8:: cursor,i

    open(unit = 91,file=trim(filename), status = "unknown", access="stream")
    cursor = fbsize+1
    write(91,pos=cursor) field 
    write(91) tiempo
    write(91) cfl
    write(91) re
    write(91) lx
    write(91) ly
    write(91) lz
    write(91) nx
    write(91) ny
    write(91) nz2
    write(91) xout
    write(91) timein
    write(91) dt
    write(91) (y(i), i =0,ny+1)
    write(91) (um(i), i = 1,ny+1)
    write(91) procs
    close(91)
    !     write(*,*) 'VALORES ESCRITOS HEADER:====================='
    !     WRITE(*,*) field,tiempo,cfl,re,lx,ly,lz,nx,ny,nz2,xout,timein,dt,procs
  end subroutine writeheader_2

  !===========BARRIER SUBROUTINES=============
  subroutine barrier_after_start(comm)
    integer, intent(in) :: comm
    integer ierr
    call MPI_BARRIER(comm,ierr)
  end subroutine barrier_after_start

  subroutine barrier_after_malloc(comm)
    integer, intent(in) :: comm
    integer ierr

    call MPI_BARRIER(comm,ierr)
  end subroutine barrier_after_malloc

  subroutine barrier_after_open(comm)
    integer, intent(in) :: comm
    integer ierr
    call MPI_BARRIER(comm,ierr)
    return
  end subroutine barrier_after_open

  subroutine barrier_after_write(comm)
    integer, intent(in) :: comm
    integer:: ierr
    call MPI_BARRIER(comm,ierr)
  end subroutine barrier_after_write

  subroutine barrier_after_read(comm)
    integer, intent(in) :: comm
    integer:: ierr
    call MPI_BARRIER(comm,ierr)
  end subroutine barrier_after_read

  subroutine barrier_after_close(comm)
    integer, intent(in) :: comm
    integer:: ierr
    call MPI_BARRIER(comm,ierr)
  end subroutine barrier_after_close
#endif

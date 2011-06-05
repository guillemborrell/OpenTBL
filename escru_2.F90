 ! ----------------------------------------------------------------------! 
  !
  ! Writing Subroutines Package
  !
  ! ----------------------------------------------------------------------! 
  !----------------------------------------------------------------------*
  
  subroutine escribezy_2(u,v,w,p,dt,mpiid,communicator)
    !----------------------------------------------------------------------*
    ! escribe campos de u,v,w,p
    !----------------------------------------------------------------------*
    use point_2
    use names_2
    use genmod_2
    use alloc_dns_2
    use ctesp_2


#ifdef WPARALLEL
    use hdf5
#endif

    implicit none
    include 'mpif.h'
    integer,intent(in):: communicator
    real*8,dimension(nz1,ny+1,ib:ie)::u,w
    real*8,dimension(nz1,ny  ,ib:ie)::p,v
    real*4,dimension(:,:,:),allocatable::resu
    integer i,j,k,l,irec
    integer status(MPI_STATUS_SIZE),ierr,t_size,t_size1,t_size2,dot,mpiid,lim1,lim2
    character(len=256):: fil1,fil2,fil3,fil4
    character:: ext1*3,uchar*1
    real*8    dt,dum(20),jk,t0
    integer:: nxr3,nyr3,nzr3,comm,tipo,chunkfbs,nfile,sidio,mpiw1,mpiw2,mpiw3,mpiw4
    integer*8:: chunks1,chunks2,chunksM1,chunksM2

#ifdef WPARALLEL
    ! ------------------------- HDF5 -------------------------------
    
    integer:: info
    integer(hid_t):: fid,pid
    integer:: h5err
    integer(hsize_t),dimension(3)::dims
#endif

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


    call MPI_INFO_CREATE(info,ierr)

#ifdef WPARALLEL
    !PARALLEL WRITER ==================================================================
    !First the header and last the field
    if (mpiid.eq.0) then 
       write(*,*) 'Escribiendo el fichero'
       write(*,*) fil1   
       t0=MPI_Wtime()  
    end if
        
    !U and w
    dims =(/ nz1, ny+1, ie-ib+1 /)
    allocate (resu(nz1,ny+1,ie-ib+1),stat=ierr)
    resu=0.0 !R4 buffer to convert R8 variables
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
    call h5pset_fapl_mpio_f(pid,comm,info,h5err)
    call h5pset_sieve_buf_size_f(pid, 4*1024*1024, h5err)
    call h5fcreate_f(trim(fil1)//".h5",H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
    call h5pclose_f(pid,h5err)

    !Dump the data to the allocated array and save to the disk
    resu=real(u(1:nz1,1:ny+1,ib:ie),kind=4)
    call MPI_BARRIER(comm,ierr)
    call h5dump_parallel(fid,"value",3,dims,mpiid,nummpi,comm,info,resu,h5err)
    call h5fclose_f(fid,h5err)

    if (mpiid == 0) write(*,*) "File for U successfully closed"

    resu=0.0
    call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
    call h5pset_fapl_mpio_f(pid,comm,info,h5err)
    call h5pset_sieve_buf_size_f(pid, 4*1024*1024, h5err)
    call h5fcreate_f(trim(fil3)//".h5",H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
    call h5pclose_f(pid,h5err)

    resu=real(w(1:nz1,1:ny+1,ib:ie),kind=4)
    call MPI_BARRIER(comm,ierr)
    call h5dump_parallel(fid,"value",3,dims,mpiid,nummpi,comm,info,resu,h5err)
    call h5fclose_f(fid,h5err)
    deallocate (resu)

    !Now the v and p
    dims =(/ nz1, ny, ie-ib+1 /)
    allocate (resu(nz1,ny,ie-ib+1),stat=ierr)
    resu=0.0 !R4 buffer to convert R8 variables
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
    call h5pset_fapl_mpio_f(pid,comm,info,h5err)
    call h5pset_sieve_buf_size_f(pid, 4*1024*1024, h5err)
    call h5fcreate_f(trim(fil2)//".h5",H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
    call h5pclose_f(pid,h5err)

    resu=real(v(1:nz1,1:ny,ib:ie),kind=4)
    call MPI_BARRIER(comm,ierr)
    call h5dump_parallel(fid,"value",3,dims,mpiid,nummpi,comm,info,resu,h5err)
    call h5fclose_f(fid,h5err)

    resu=0.0

    call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
    call h5pset_fapl_mpio_f(pid,comm,info,h5err)
    call h5pset_sieve_buf_size_f(pid, 4*1024*1024, h5err)
    call h5fcreate_f(trim(fil4)//".h5",H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
    call h5pclose_f(pid,h5err)

    resu = real(p(1:nz1,1:ny,ib:ie),kind=4)
    call MPI_BARRIER(comm,ierr)
    call h5dump_parallel(fid,"value",3,dims,mpiid,nummpi,comm,info,resu,h5err)
    call h5fclose_f(fid,h5err)

    deallocate (resu)

    !! Writing the Headers
    if (mpiid == 0) then
       write(*,*) "Writing the headers"
       write(*,*) "Data dimensions ", dims(1), dims(2), dims(3)
       write(*,*) "                ", nz1, ny+1, ie-ib+1
       write(*,*) "Write ",dims(1)*dims(2)*dims(3)*4/1024/1024*nummpi," Mbytes per field"
       call writeheader(fil1,'u',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,&
            & xout,timeinit,dt,y,um,nummpi)
       call writeheader(fil2,'v',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,&
            & xout,timeinit,dt,y,um,nummpi)
       call writeheader(fil3,'w',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,&
            & xout,timeinit,dt,y,um,nummpi)
       call writeheader(fil4,'p',tiempo,cfl,re,ax*pi,ay*pi,az*2*pi,nx,ny,nz2,&
            & xout,timeinit,dt,y,um,nummpi)
    end if
  
    ifile=ifile+1        

    call MPI_BARRIER(comm,ierr)
 
    if (mpiid.eq.0) then 
       t0=MPI_Wtime()-t0
       write(*,*)
       write(*,*) '===================================================================='
       write(*,*) 'Done writting fields'
       write(*,*) '===================================================================='    
       write(*,'(a20,f10.3,a3)') 'TIME SPENT IN WRITING:',t0,'sc'
       write(*,*)   '------------------------------------------------------------------'
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

       open(39,file=stfile,status='unknown',form='unformatted');rewind(39)
       indst=indst!+1 !WE WRITE STATISTICS WHEN RECORD IMAGES ONLY!!
       write(*,*) 'writing in==============  ',stfile       
       write(*,'(6f12.4,4i10)') tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical
       
       write(39) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical       
       write(39) (y(i), i=0,ny+1)
    end if

    call MPI_BCAST(corfile,100,MPI_CHARACTER,0,comm,ierr)  !for the Parallel IO Writting
     ! master only things  
    allocate(wkn(ny,13),wknp(ny+1,4))
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
          wkn (:,12)=uw (1:ny,i)
          wkn (:,13)=vw (1:ny,i)

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
               &    (wkn (j,11),j=1,ny  ),&
               &    (wkn (j,12),j=1,ny  ),&
               &    (wkn (j,13),j=1,ny  )
       enddo

       do dot = 1,nummpi-1
          do i=ibeg(dot),iend(dot)          
                   call MPI_RECV(wknp,size(wknp),tipo,dot,0,comm,status,ierr)
                   call MPI_RECV(wkn ,size(wkn) ,tipo,dot,0,comm,status,ierr)

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
                  &    (wkn (j,11),j=1,ny  ),&
                  &    (wkn (j,12),j=1,ny  ),&
                  &    (wkn (j,13),j=1,ny  )
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
          wkn (:,12)=uw (1:ny,i)
          wkn (:,13)=vw (1:ny,i)

          call MPI_SEND(wknp,size(wknp),tipo,0,0,comm,ierr)
          call MPI_SEND(wkn ,size(wkn),tipo,0,0,comm,ierr)
       enddo
    endif
    call mpi_barrier(comm,ierr)
    ! Initialize everything to zero
    pp=0d0;pm=0d0;ua=0d0;va=0d0 
    us=0d0;ws=0d0;wa=0d0;uv=0d0;vs=0d0;uw=0d0;vw=0d0;
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
       open(55,file=budfile,status='unknown',form='unformatted');rewind(55)
       write(*,*) 'writing in==============',budfile
       !HEADER
       write(55) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical    
       write(55) (y(i), i=0,ny+1)
    end if

#ifdef INFOINTER 
   nbud=21+9
#else    
   nbud=21
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
          buf_bud(:,18)=dvdx0 (1:ny,i)
          buf_bud(:,19)=dvdy0 (1:ny,i)          
          buf_bud(:,20)=dwdx0 (1:ny,i)
          buf_bud(:,21)=dwdy0 (1:ny,i)                    
#ifdef INFOINTER 
         buf_bud(:,22)=v_0    (1:ny,i)
         buf_bud(:,23)=u_x0   (1:ny,i)
         buf_bud(:,24)=u_xy0  (1:ny,i)
         buf_bud(:,25)=w_0    (1:ny,i)
         buf_bud(:,26)=w_y0   (1:ny,i)
         buf_bud(:,27)=dwdx_0 (1:ny,i)
         buf_bud(:,28)=dudz_x0(1:ny,i) 
         buf_bud(:,29)=v_y0   (1:ny,i)         
         buf_bud(:,30)=dudx_0 (1:ny,i)            
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
          buf_bud(:,18)=dvdx0 (1:ny,i)
          buf_bud(:,19)=dvdy0 (1:ny,i)          
          buf_bud(:,20)=dwdx0 (1:ny,i)
          buf_bud(:,21)=dwdy0 (1:ny,i)                    
#ifdef INFOINTER 
         buf_bud(:,22)=v_0    (1:ny,i)
         buf_bud(:,23)=u_x0   (1:ny,i)
         buf_bud(:,24)=u_xy0  (1:ny,i)
         buf_bud(:,25)=w_0    (1:ny,i)
         buf_bud(:,26)=w_y0   (1:ny,i)
         buf_bud(:,27)=dwdx_0 (1:ny,i)
         buf_bud(:,28)=dudz_x0(1:ny,i) 
         buf_bud(:,29)=v_y0   (1:ny,i)         
         buf_bud(:,30)=dudx_0 (1:ny,i)            
#endif          
          !NY+1 size Buffers
          buf_bud2(:,1)=u3    (1:ny+1,i)
          buf_bud2(:,2)=w2u   (1:ny+1,i)                                                          
          call MPI_SEND(buf_bud ,size(buf_bud) ,tipo,0,0,comm,ierr)
          call MPI_SEND(buf_bud2,size(buf_bud2),tipo,0,1,comm,ierr)
       enddo
    endif
    ! Initialize everything to zero
    call mpi_barrier(comm,ierr)

    dispu=0d0;dispv=0d0;dispw=0d0;dispuv=0d0;
    pvp=0d0;pup=0d0;pdudx=0d0;pdudy=0d0;pdvdx=0d0
    pdvdy=0d0;pdwdz=0d0
    u3=0d0;v3=0d0;u2v=0d0;v2u=0d0;w2v=0d0;w2u=0d0
    dudx0=0d0;dudy0=0d0;
    dvdx0=0d0;dvdy0=0d0;
    dwdx0=0d0;dwdy0=0d0; 

#ifdef INFOINTER 
    v_0=0d0;u_x0=0d0;u_xy0=0d0;w_0=0d0;w_y0=0d0;
    dwdx_0=0d0;dudz_x0=0d0;v_y0=0d0;dudx_0=0d0
#endif 
    deallocate(buf_bud,buf_bud2)
#endif


#ifndef NOCORR

#ifdef WSERIAL
    !===================CORRELATIONS=======================    
    if (mpiid.eq.0) then
       open(47,file=corfile,status='unknown',form='unformatted',convert='Big_endian');rewind(47)
       write(*,*) 'writing in==============  ',corfile
       !HEADER
       write(47) tiempo,cfl,Re,ax,ay,az,nx,ny,nz2,ical,ncorr,lxcorr,nxp(1:lxcorr),xcorpoint(1:lxcorr)      
       write(47) y(0:ny+1),jspecy(1:ncorr,1:lxcorr)
       call flush(47)     
    end if 
  
    call escr_corr_2(coru,corv,coruv,corw,corp,corox,coroy,coroz,coruw,corvw,mpiid,communicator)  
    call mpi_barrier(comm,ierr)
    close(47) 
#endif

#ifdef WPARALLEL
    if(mpiid.eq.0) write(*,*) 'Correlation file name before calling escr_corr:', corfile
    call escr_corr_2(corfile,ical,coru,corv,coruv,corw,corp,corox,coroy,coroz,&
         & coruw,corvw,mpiid,nummpi,comm)
#endif


    ! Initialize everything to zero
    coru=0d0;corv=0d0;corw=0d0;coruv=0d0;coruw=0d0;corvw=0d0;
    corox=0d0;coroy=0d0;coroz=0d0;corp=0d0;    
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
     call mpi_barrier(comm,ierr)
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

! call escr_corr(coru,corv,coruv,corw,corp,corox,coroy,coroz,coruw,corvw,mpiid)
#ifdef WSERIAL
subroutine escr_corr_2(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,rank,communicator)
use ctesp_2,only:nx,nz2,ncorr,lxcorr,nummpi,nxp,ny
use point_2
implicit none
include 'mpif.h'
integer,intent(in):: communicator
real*8,allocatable,dimension(:,:)::buf_cor 
real*8,dimension(nx,pcib2:pcie2,lxcorr)::c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
integer:: rank,ierr,comm,status(MPI_STATUS_SIZE),i,j,dot,tipo,k,l

 comm = communicator
 tipo = MPI_real8
  
allocate(buf_cor(1:nx,10)) !10 Correlations

  do j=1,lxcorr
   if (rank.eq.0) write(*,*) 'writing correlations',j,'of',lxcorr
    if (rank.eq.0) then
       do i=pcib2,pcie2
          buf_cor(:,1) =c1 (1:nx,i,j)       
          buf_cor(:,2) =c2 (1:nx,i,j)
          buf_cor(:,3) =c3 (1:nx,i,j)
          buf_cor(:,4) =c4 (1:nx,i,j)
          buf_cor(:,5) =c5 (1:nx,i,j)               
          buf_cor(:,6) =c6 (1:nx,i,j)
          buf_cor(:,7) =c7 (1:nx,i,j)
          buf_cor(:,8) =c8 (1:nx,i,j)          
          buf_cor(:,9) =c9 (1:nx,i,j)          
          buf_cor(:,10)=c10(1:nx,i,j)          
          buf_cor=buf_cor/(2d0*nxp(j)+1d0) !averaging         
          write(47) ((buf_cor(k,l),k=1,nx),l=1,10)  
       enddo

       do dot = 1,nummpi-1
          do i=pcibeg2(dot),pciend2(dot)          
             call MPI_RECV(buf_cor,10*nx,tipo,dot,j,comm,status,ierr)
             write(47) ((buf_cor(k,l),k=1,nx),l=1,10)                                 
          enddo
          call flush(47)        
       enddo          
    else   
       do i=pcib2,pcie2
          buf_cor(:,1) =c1 (1:nx,i,j)       
          buf_cor(:,2) =c2 (1:nx,i,j)
          buf_cor(:,3) =c3 (1:nx,i,j)
          buf_cor(:,4) =c4 (1:nx,i,j)
          buf_cor(:,5) =c5 (1:nx,i,j)               
          buf_cor(:,6) =c6 (1:nx,i,j)
          buf_cor(:,7) =c7 (1:nx,i,j)
          buf_cor(:,8) =c8 (1:nx,i,j)          
          buf_cor(:,9) =c9 (1:nx,i,j)          
          buf_cor(:,10)=c10(1:nx,i,j)          
          buf_cor=buf_cor/(2d0*nxp(j)+1d0) !averaging        
          call MPI_SEND(buf_cor,10*nx,tipo,0,j,comm,ierr)
       enddo    
    endif
 enddo 
 deallocate(buf_cor)
endsubroutine escr_corr_2
#endif


#ifdef WPARALLEL
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! ------------ PARALLEL WRITTING SUBROUTINES ------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 
  ! -------------------------------------------------------------------! 


subroutine escr_corr_2(fname,ical,coru,corv,coruv,corw,corp,corox,coroy,coroz,&
     &coruw,corvw,mpiid,nummpi,comm)

use ctesp_2, only:nx,nz2,ncorr,lxcorr,nxp,&
     & nx,ny,nz2,xcorpoint,nspec
use alloc_dns_2, only:y,Re,ax,ay,az,tiempo,cfl
use statistics_2, only:jspecy
use point_2, only: pcib2,pcie2,mp_corr2
use hdf5
use h5lt

implicit none
include 'mpif.h'

integer(hid_t):: fid,pid
character(len=*), intent(in):: fname
integer, intent(in):: mpiid,nummpi,comm,ical
real(8), dimension(nx,pcib2:pcie2,lxcorr)::coru,corv,coruv,corw,corp
real(8), dimension(nx,pcib2:pcie2,lxcorr)::corox,coroy,coroz,coruw,corvw
real(8), dimension(:,:,:,:), allocatable:: buf_cor
real(8), dimension(:,:), allocatable:: aux_buf_cor
integer(hsize_t), dimension(1):: hdims = (/ 1 /)
integer(hsize_t), dimension(2):: hdims2 = (/ 0, 0 /) 
integer,dimension(nummpi):: npencils
integer:: h5err,mpierr
integer:: i,j,k,l,npen,id
integer(hid_t)::dset,dspace,mspace
real(8), parameter:: pi = 3.14159265358979
real(8):: timer

npen = pcie2-pcib2+1

if(mpiid == 0) timer = MPI_WTIME()

allocate(buf_cor(nx,npen,lxcorr,10))
allocate(aux_buf_cor(lxcorr*10,npen))

buf_cor(:,:,:,1)  = coru
buf_cor(:,:,:,2)  = corv
buf_cor(:,:,:,3)  = coruv
buf_cor(:,:,:,4)  = corw
buf_cor(:,:,:,5)  = corp
buf_cor(:,:,:,6)  = corox
buf_cor(:,:,:,7)  = coroy
buf_cor(:,:,:,8)  = coroz
buf_cor(:,:,:,9)  = coruw
buf_cor(:,:,:,10) = corvw

!! Scale or whatever Juan did.
do j = 1,lxcorr
   buf_cor(:,:,j,:) = buf_cor(:,:,j,:)/(2d0*nxp(j)+1d0)
end do

call mpi_barrier(comm,mpierr)

call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
call h5pset_fapl_mpio_f(pid,comm,MPI_INFO_NULL,h5err)
call h5pset_sieve_buf_size_f(pid, 4*1024*1024, h5err)
call h5fcreate_f(trim(fname)//".h5",H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
call h5pclose_f(pid,h5err) !! Close property access list

!!Write the correlations concurrently to disk
call dump_corr(fid,"corr",nx,lxcorr,pcib2,pcie2,(nz2+1)*nspec,mpiid,nummpi,comm,&
     & buf_cor,h5err)

call mpi_barrier(comm,mpierr)

if(mpiid == 0) then
   write(*,*) "time of everything"
   write(*,*) MPI_WTIME()-timer
end if
if (mpiid.eq.0) write(*,*) 'closing.....'
call h5fclose_f(fid,h5err)
if (mpiid.eq.0) write(*,*) 'closed.....'
deallocate(buf_cor)
deallocate(aux_buf_cor)

!Write the header

!Get the number of pencils written by every node
call MPI_GATHER(pcie2-pcib2+1,1,MPI_INTEGER,&
     & npencils,1,MPI_INTEGER,0,comm,mpierr)

if (mpiid == 0) then
   write(*,*) 'CORRELATION FILE TO BE WRITED:', trim(fname)
   call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,fid,h5err)
   
   call h5ltmake_dataset_double_f(fid,"tiempo",1,hdims,(/tiempo/),h5err)
   call h5ltmake_dataset_double_f(fid,"cfl"   ,1,hdims,(/cfl/),h5err)
   call h5ltmake_dataset_double_f(fid,"Re"    ,1,hdims,(/re/),h5err)
   call h5ltmake_dataset_double_f(fid,"lx"    ,1,hdims,(/ax*pi/),h5err)
   call h5ltmake_dataset_double_f(fid,"ly"    ,1,hdims,(/ay*pi/),h5err)
   call h5ltmake_dataset_double_f(fid,"lz"    ,1,hdims,(/2*az*pi/),h5err)
   call h5ltmake_dataset_int_f(fid,"nx",    1,hdims,(/nx/),    h5err)
   call h5ltmake_dataset_int_f(fid,"ny",    1,hdims,(/ny/),    h5err)
   call h5ltmake_dataset_int_f(fid,"nz2",   1,hdims,(/nz2/),   h5err)
   call h5ltmake_dataset_int_f(fid,"ical",  1,hdims,(/ical/),  h5err)
   call h5ltmake_dataset_int_f(fid,"ncorr", 1,hdims,(/ncorr/), h5err)
   call h5ltmake_dataset_int_f(fid,"lxcorr",1,hdims,(/lxcorr/),h5err)
   call h5ltmake_dataset_int_f(fid,"mpisize",1,hdims,(/nummpi/),h5err)
   hdims = (/ lxcorr /)
   call h5ltmake_dataset_int_f(fid,"nxp",1,hdims,nxp,h5err)
   call h5ltmake_dataset_int_f(fid,"xcorpoint",1,hdims,xcorpoint,h5err)
   hdims = (/ ny+2 /)
   call h5ltmake_dataset_double_f(fid,"y",1,hdims,y,h5err)
   hdims = (/ nummpi /)
   call h5ltmake_dataset_int_f(fid,"npencils",1,hdims,npencils,h5err)
   hdims2 = (/ ncorr,lxcorr /)
   call h5ltmake_dataset_int_f(fid,"jspecy",2,hdims2,jspecy,h5err)

   call h5fclose_f(fid,h5err)

   write(*,*) MPI_WTIME()-timer

end if

!End write the header


end subroutine escr_corr_2

#endif

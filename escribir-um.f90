program escribir
integer ny,nz2,nxr,nyr,nzr,nz1r,nzz,j,i,k,l,dot,lim2,rsize,irec,rsize1,rsize2,ji
real(8):: jk,dt,dum(20),a,b,tiempo,timeinit
real(8):: um(646+1),dummy(6145+2),y(646+2)
character:: chinit*100,uchar*1

 chinit='/intrepid-fs0/users/sillero/scratch/images/bg-6145x-646y.009.u'
 write(*,*) 'Leyendo del fichero',chinit
 ny=646;nz2=511 
    
     rsize=2*nz2*(ny+1)*4
     
     open (10,file=chinit,status='old',form='unformatted',access='direct',recl=rsize)!,convert='BIG_ENDIAN')
     write(*,*) 'BG: file open, reading tiempo and dimensions'
     read(10,rec=1) uchar,tiempo,jk,jk,jk,jk,jk,nxr,nyr,nzr,ji,timeinit,dt, &
          & dum, (dummy(i), i=0,nxr+1), (y(i), i=0,nyr+1), (um(i), i=1,nyr+1)     
     write(*,*) '==========================================='
     write(*,*) 'in file    ', uchar,tiempo,nxr,nyr,nzr
     close(10)
     
     open (11,file='um-extractedBG',status='unknown',form='unformatted')!,convert='BIG_ENDIAN')
     write(11) um(1:nyr+1)
     write(11) y(0:nyr+1)
     write(11) dummy(0:nxr+1)
     close(11)
     write(*,*) 'done writting um,y,x'


end program escribir


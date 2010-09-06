!===================================================== 
!!!!! THIS PROGRAM JUST COMPUTE THE RE_THETA BASED ON:
!!!!! The profile U
!!!!! The grid Y
!===================================================== 

program compute
implicit none
integer,parameter:: ny=646
integer:: jtop,j,i
real*8:: y(0:647),um(647),uinf,din,utauout,rthout,dout,dy(0:646)
real*4:: um2(647)
real*8,parameter:: re=3.5d+03


 open (103,file='v0-nagib',form='unformatted',status='old')
        read(103) um(1:ny)
        close(103)
        !um=um2
        write(*,*) '-----------------------------------------------------um leido del fichero um50d99'
        do j=1,10
         write(*,*) um(j)
        enddo
open (110,file='ybg713',form='unformatted',status='old')
        read(110) (y(i), i=0,ny+1)
        close(110) 



!um= (dt/timec)*ut(0,:,1)+um*(1d0-dt/timec)
     utauout = sqrt(abs(2d0*(um(2))/(y(2)-y(1))/re))

     !  --- compute momentum thickness -----
     uinf=sum(um(ny+1-20:ny+1-5))/16d0 
     write(*,*) 'Uinfinito=',uinf
     din = 0.99*uinf   
     write(*,*) 'd99=      ',din
     do j=1,ny-2
        if (um(j).ge.din) exit
        jtop = j+1
     enddo
     
     dout = y(jtop)
      write(*,*) 'jtop=    ',jtop
      write(*,*) 'altura=  ',y(jtop)
     rthout  = 0d0       
!      uinf=um(ny+1-5)     !!!!!!!!  fix this !!!!!!!
     dy(0:ny) = y(1:ny+1)-y(0:ny)
     rthout = sum((uinf-um(2:jtop+3))*um(2:jtop+3)*dy(1:jtop+2))*re/uinf
      write(*,*) 'Re_theta=',rthout

end program compute

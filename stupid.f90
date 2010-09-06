program stupid

integer:: nx,ny,nz,xin,xout,nz1,nz2,ngz,nx1,ny1
!====================================================
 parameter ( nx =8193,   ny =711, nz=4096)
 parameter ( xin = 1 , xout =3671*2/2) !50d99-Re2500
!====================================================

  parameter ( nz1 = 2*(nz/3), nz2=nz1/2-1,ngz=nz/2,nx1=nx-1,ny1=ny-1 )


write(*,*) 'nx,ny,nz,nz1,nz2',nx,ny,nz,nz1,nz2

endprogram stupid 

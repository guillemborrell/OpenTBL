program main_writer

  implicit none

  include "mpif.h"

  character(len = 99) :: filename = "test.dat"
  integer:: ierr, mympi, nummpi
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mympi,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nummpi,ierr)

  call posix_open(trim(filename),mympi)
  call posix_close()

  call MPI_FINALIZE(ierr)

end program main_writer

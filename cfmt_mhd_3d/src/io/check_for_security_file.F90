module check_security_file_mod

  implicit none

  contains

  function check_for_security_file() result(do_fast_exit)

  use precision
  use mo_mpi
  implicit none

  integer(kind=ik) :: fast_exit, ier
  logical          :: do_fast_exit
  logical          :: file_exists
  

  fast_exit = 0

  if (myproc .eq. 0) then
     
     inquire(file='fastexit',exist=file_exists)
     if (file_exists) then
        fast_exit = 1
     endif
     
  endif
         
  ! communicate this to all tasks
               
  call MPI_BCAST(fast_exit, 1, MPI_INTEGER, 0, &
                 MPI_COMM_WORLD, ier) 

  if (fast_exit .eq. 1) then
    do_fast_exit = .true.
  else
    do_fast_exit = .false.
 endif
 
  end function check_for_security_file

end module check_security_file_mod

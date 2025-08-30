!>
!> \details
!> This file contains "dummy" MPI subroutines, that can be called if
!> the compilation was donw without the MPI libs. All subroutines
!> mimmic a MPI behaviour of _ONE_ MPI-Task
!>
!> \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 845 $
!>   $Date: 2010-03-04 21:42:15 +0100 (Thu, 04 Mar 2010) $
!>   
!> \endverbatim
!>

module mpi_stubs
  use precision
  use abort

  implicit none

  save

  private
 
! MPI overloaded function in alphabetical order 
  public  mpi_abort,       &
          mpi_allgatherv,  &
          mpi_allreduce,   &
          mpi_barrier,     &
          mpi_bcast,       &
          mpi_cart_create, &
          mpi_cart_shift,  &
          mpi_comm_rank,   &
          mpi_comm_size,   &
          mpi_finalize,    &
          mpi_init,        &
          mpi_irecv,       &
          mpi_send,        &
          mpi_sendrecv,    &
          mpi_waitall,     &
          mpi_wtime

 
  public &  ! MPI stubs
         mpi_comm_world,        &
         MPI_DOUBLE_PRECISION,  &
         MPI_2DOUBLE_PRECISION, &
         MPI_INTEGER,           &
         MPI_MAX,               &
         MPI_MAXLOC,            &
         MPI_MAX_PROCESSOR_NAME,&
         MPI_MIN,               &
         MPI_REAL8,             &
         MPI_STATUS_SIZE,       &
         MPI_SUCCESS,           &
         MPI_SUM

  integer(kind=ik) :: MPI_SUCCESS, mpi_comm_world, MPI_SUM, MPI_DOUBLE_PRECISION
  integer(kind=ik) :: MPI_MAX, MPI_INTEGER, MPI_2DOUBLE_PRECISION, MPI_MIN
  integer(kind=ik) :: MPI_MAXLOC, MPI_REAL8
  integer(kind=ik) :: MPI_STATUS_SIZE = 1
  integer(kind=ik), parameter :: MPI_MAX_PROCESSOR_NAME=40



  interface MPI_Allgatherv
     module procedure MPI_Allgatherv_r8_3d
  end interface

  interface mpi_allreduce
     module procedure mpi_allreduce_i4
     module procedure mpi_allreduce_i4_1d
     module procedure mpi_allreduce_r8
     module procedure mpi_allreduce_r8_1d
     module procedure mpi_allreduce_r8_2d
  end interface

  interface mpi_irecv
     module procedure mpi_irecv_r8
     module procedure mpi_irecv_r8_1d
     module procedure mpi_irecv_r8_2d
  end interface

  interface mpi_send
     module procedure mpi_send_r8
     module procedure mpi_send_r8_1d
     module procedure mpi_send_r8_2d
  end interface

  interface mpi_sendrecv
     module procedure mpi_sendrecv_r8_1d
     module procedure mpi_sendrecv_r8_2d
     module procedure mpi_sendrecv_r8_3d
     module procedure mpi_sendrecv_r8_4d
     module procedure mpi_sendrecv_r8_5d
  end interface

contains

subroutine mpi_abort(mpi_comm, errorcode, ierr)

  use precision
  implicit none
  
  integer(kind=ik) :: mpi_comm
  integer(kind=ik) :: errorcode
  integer(kind=ik) :: ierr

  raise_abort("Calling MPI_ABORT: stopping now!")

  return
end subroutine mpi_abort

subroutine mpi_allgatherv_r8_3d(sbuf, sendcount, datatype1, rbuf, rcvcount, &
                                displs, datatype2, comm, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:,:), rbuf(:,:,:)
  integer(kind=ik) :: sendcount, datatype1, rcvcount(:), datatype2, displs(:), comm
  integer(kind=ik) :: ierr


  write (*,*) "This is the dummy version of MPI_allgather_3D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_allgatherv_r8_3d

subroutine mpi_allreduce_i4(source, dest, size_mpi, type, method, comm, err)

  use precision

  implicit none

  integer(kind=ik) :: size_mpi
  
  integer(kind=ik) :: source
  integer(kind=ik) :: dest
  
  integer(kind=ik) :: type
  integer(kind=ik) :: method
  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_ALLREDUCE_r8_1d"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")



  return
end subroutine mpi_allreduce_i4

subroutine mpi_allreduce_i4_1d(source, dest, size_mpi, type, method, comm, err)

  use precision

  implicit none

  integer(kind=ik) :: size_mpi
  
  integer(kind=ik) :: source(:)
  integer(kind=ik) :: dest(:)
  
  integer(kind=ik) :: type
  integer(kind=ik) :: method
  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_ALLREDUCE_r8_1d"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")



  return
end subroutine mpi_allreduce_i4_1d

subroutine mpi_allreduce_r8(source, dest, size_mpi, type, method, comm, err)

  use precision

  implicit none

  integer(kind=ik) :: size_mpi
  
  real(kind=rk) :: source
  real(kind=rk) :: dest
  
  integer(kind=ik) :: type
  integer(kind=ik) :: method
  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_ALLREDUCE_r8_1d"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")



  return
end subroutine mpi_allreduce_r8

subroutine mpi_allreduce_r8_1d(source, dest, size_mpi, type, method, comm, err)

  use precision

  implicit none

  integer(kind=ik) :: size_mpi
  
  real(kind=rk) :: source(:)
  real(kind=rk) :: dest(:)
  
  integer(kind=ik) :: type
  integer(kind=ik) :: method
  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_ALLREDUCE_r8_1d"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")



  return
end subroutine mpi_allreduce_r8_1d


subroutine mpi_allreduce_r8_2d(source, dest, size_mpi, type, method, comm, err)

  use precision

  implicit none

  integer(kind=ik) :: size_mpi
  
  real(kind=rk) :: source(:,:)
  real(kind=rk) :: dest(:,:)
  
  integer(kind=ik) :: type
  integer(kind=ik) :: method
  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_ALLREDUCE_r8_2d"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")



  return
end subroutine mpi_allreduce_r8_2d



subroutine mpi_barrier(mpi_comm, mpierr)

  use precision
  implicit none

  integer(kind=ik) :: mpi_comm
  integer(kind=ik) :: mpierr

  mpierr = 0

  return

end subroutine mpi_barrier

subroutine mpi_bcast(sen_buf, nums, type, rank, MPI_COMM_WORLD, ier)

  use precision
  implicit none

  integer(kind=ik) :: sen_buf, nums, type, rank, MPI_COMM_WORLD, ier
  integer(kind=ik) :: mpierr

  ier = 0

  return

end subroutine mpi_bcast


subroutine mpi_cart_create(MPI_COMM_WORLD,ndim,cart_dims,cart_periods, &
                           cart_reorder,cart_comm,ier)

  use precision
  implicit none

  integer(kind=ik) :: MPI_COMM_WORLD, ier, cart_comm, ndim 
  logical          :: cart_periods(0:1),cart_reorder 
  integer(kind=ik) :: cart_dims(0:1)


  write (*,*)  " This subroutine is a overloaded procedure of"
  write (*,*)  " mpi_cart_create. You are not allowed to see"
  write (*,*)  " it! "

  raise_abort("An MPI stub routine was called, this must not happen")
  return
end subroutine mpi_cart_create

subroutine mpi_cart_shift( cart, in, out, src, dest, err)
  use precision

  implicit none  

  integer(kind=ik) :: cart
  integer(kind=ik) :: in
  integer(kind=ik) :: out
  integer(kind=ik) :: src
  integer(kind=ik) :: dest
!  integer(kind=ik) :: comm
  integer(kind=ik) :: err

  write (*,*) "This is the dummy version of MPI_CARTSHIFT"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")
end subroutine mpi_cart_shift


subroutine mpi_comm_rank(mpi_comm_world, myid, mpierr)
  use precision

  implicit none

  integer(kind=ik) :: mpi_comm_world
  integer(kind=ik) :: mpierr
  integer(kind=ik) :: myid

  myid = 0
  mpierr = 0

  return
end subroutine mpi_comm_rank

subroutine mpi_comm_size(mpi_comm_world, ntasks, mpierr)

  use precision
  
  implicit none

  integer(kind=ik) :: mpi_comm_world
  integer(kind=ik) :: ntasks
  integer(kind=ik) :: mpierr

  ntasks = 1
  mpierr = 0

  return

end subroutine mpi_comm_size


subroutine mpi_init(mpierr)
  use precision
  
  implicit none

  integer(kind=ik) :: mpierr

  mpierr = -1

  return
end subroutine mpi_init

subroutine mpi_irecv_r8(rbuf, rcvcount, datatype1, src, dest, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: rbuf
  integer(kind=ik) :: rcvcount, datatype1, dest, src, comm
  integer(kind=ik) :: mpistat, ierr


  write (*,*) "This is the dummy version of MPI_irecv_r8"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_irecv_r8

subroutine mpi_irecv_r8_1d(rbuf, rcvcount, datatype1, src, dest, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: rbuf(:)
  integer(kind=ik) :: rcvcount, datatype1, dest, src, comm
  integer(kind=ik) :: mpistat(:), ierr


  write (*,*) "This is the dummy version of MPI_irecv_1D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_irecv_r8_1d


subroutine mpi_irecv_r8_2d(rbuf, rcvcount, datatype1, src, dest, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: rbuf(:,:)
  integer(kind=ik) :: rcvcount, datatype1, dest, src, comm
  integer(kind=ik) :: mpistat(:), ierr


  write (*,*) "This is the dummy version of MPI_irecv_2D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_irecv_r8_2d



subroutine mpi_finalize(mpierr)
  use precision
  implicit none
  
  integer(kind=ik) :: mpierr

  mpierr = 0

  return
end subroutine mpi_finalize

subroutine mpi_send_r8(sbuf, sendcount, datatype1, src, dest, comm, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf
  integer(kind=ik) :: sendcount, datatype1, src, dest, comm
  integer(kind=ik) :: ierr


  write (*,*) "This is the dummy version of MPI_SEND_r8"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_send_r8

subroutine mpi_send_r8_1d(sbuf, sendcount, datatype1, src, dest, comm, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:)
  integer(kind=ik) :: sendcount, datatype1, src, dest, comm
  integer(kind=ik) :: ierr


  write (*,*) "This is the dummy version of MPI_SEND_r8"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_send_r8_1d

subroutine mpi_send_r8_2d(sbuf, sendcount, datatype1, src, dest, comm, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:)
  integer(kind=ik) :: sendcount, datatype1, src, dest, comm
  integer(kind=ik) :: ierr


  write (*,*) "This is the dummy version of MPI_SEND_2D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_send_r8_2d



subroutine mpi_sendrecv_r8_1d(sbuf, sendcount, datatype1, dest, tag1, rbuf, &
                              recvcount, datatype2, src, tag2, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:), rbuf(:)
  integer(kind=ik) :: sendcount, datatype1, dest, tag1, recvcount, src, comm
  integer(kind=ik) :: datatype2, tag2
  integer(kind=ik) :: mpistat(:), ierr


  write (*,*) "This is the dummy version of MPI_SENDRECV_1D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_sendrecv_r8_1d


subroutine mpi_sendrecv_r8_2d(sbuf, sendcount, datatype1, dest, tag1, rbuf, &
                              recvcount, datatype2, src, tag2, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:), rbuf(:,:)
  integer(kind=ik) :: sendcount, datatype1, dest, tag1, recvcount, src, comm
  integer(kind=ik) :: datatype2, tag2
  integer(kind=ik) :: mpistat(:), ierr

  write (*,*) "This is the dummy version of MPI_SENDRECV_2D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_sendrecv_r8_2d

subroutine mpi_sendrecv_r8_3d(sbuf, sendcount, datatype1, dest, tag1, rbuf, &
                              recvcount, datatype2, src, tag2, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:,:), rbuf(:,:,:)
  integer(kind=ik) :: sendcount, datatype1, dest, tag1, recvcount, src, comm
  integer(kind=ik) :: datatype2, tag2
  integer(kind=ik) :: mpistat(:), ierr

  write (*,*) "This is the dummy version of MPI_SENDRECV_3D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_sendrecv_r8_3d


subroutine mpi_sendrecv_r8_4d(sbuf, sendcount, datatype1, dest, tag1, rbuf, &
                              recvcount, datatype2, src, tag2, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:,:,:), rbuf(:,:,:,:)
  integer(kind=ik) :: sendcount, datatype1, dest, tag1, recvcount, src, comm
  integer(kind=ik) :: datatype2, tag2
  integer(kind=ik) :: mpistat(:), ierr

  write (*,*) "This is the dummy version of MPI_SENDRECV_3D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_sendrecv_r8_4d

subroutine mpi_sendrecv_r8_5d(sbuf, sendcount, datatype1, dest, tag1, rbuf, &
                              recvcount, datatype2, src, tag2, comm, mpistat, ierr)

  use precision

  implicit none

  real(kind=rk) :: sbuf(:,:,:,:,:), rbuf(:,:,:,:,:,:)
  integer(kind=ik) :: sendcount, datatype1, dest, tag1, recvcount, src, comm
  integer(kind=ik) :: datatype2, tag2
  integer(kind=ik) :: mpistat(:), ierr

  write (*,*) "This is the dummy version of MPI_SENDRECV_3D"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_sendrecv_r8_5d

subroutine mpi_waitall(nr, ireq, istat, ier)

  use precision

  implicit none

  integer(kind=ik) :: nr, ireq(:), istat(:,:), ier


  write (*,*) "This is the dummy version of MPI_waitall"
  write (*,*) "you should never see this message"
  raise_abort("An MPI stub routine was called, this must not happen")

end subroutine mpi_waitall

function mpi_wtime() result(time)
  use precision
  use abort

  implicit none
  real(kind=rk) :: time

  time = 0._rk

  raise_abort("do not use mpi_wtime without using MPI")

end function mpi_wtime


end module mpi_stubs

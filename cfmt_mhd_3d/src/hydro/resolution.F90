module resolution_checker
  implicit none

  private
  public check_resolution

  contains

  function check_resolution() result(res)
    use precision
    use abort
    use totare_hy, only : dentot, xzntot, ishtot
    use mo_mpi

    implicit none

    logical          :: res
    integer(kind=ik) :: i, j, k
    integer          :: mpi_err
    real(kind=rk)    :: max_resolution_criterium, recv_buf

    res = .false.

#ifdef MPI_HYDRO
    raise_abort("Not yet tested with MPI")
#endif
   
    ! Just 1D at the moment
    max_resolution_criterium = 0
    do j = qy_s, qy_e
      do k = qz_s, qz_e
        write(*,*) "Task ",myproc, j, k, ": running up to i = ", maxloc(ishtot(1:,j,k), dim=1) - 10
        do i = 1, maxloc(ishtot(1:,j,k), dim=1) - 10
          max_resolution_criterium = max(max_resolution_criterium,  &
                                         abs(log10(dentot(i + 1, 1, 1) / dentot(i,1,1))) * 10.0_rk)
          if (ishtot(i,j,k) .eq. 1_ik) then
            exit
          endif
        enddo
      enddo
    enddo

    if (use_mpi) then
      call MPI_Allreduce(max_resolution_criterium, recv_buf, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)
      abort_if(mpi_err .ne. MPI_SUCCESS)
      max_resolution_criterium = recv_buf
    endif

    if(max_resolution_criterium .gt. 0.6) then
       res = .true.
    endif
  end function check_resolution

end module resolution_checker

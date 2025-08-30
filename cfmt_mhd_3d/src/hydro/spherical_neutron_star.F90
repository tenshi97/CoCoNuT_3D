module spherical_neutron_star_mod

  implicit none

  contains

  subroutine spherical_neutron_star
  use mo_mpi
  use precision
 
  use totare_hy, only : dentot
  use mapare_proc, only : mapare
 
  use hydro_areas_mod
  use configure
  use print_stdout_mod
  implicit none

  
  integer(kind=ik) :: i, j, k, nsb, nsb_rcv, ierr

         nsb = areas%ix_are(1,2)

         do i = 1,config%qx
            do j = qy_s,qy_e
              do k = qz_s,qz_e
                if (dentot(i,j,k) .gt. config%spherical_neutron_star) then
                  nsb = i
                endif
              enddo
            enddo
         enddo

         if (use_mpi) then
            nsb_rcv = 0
            call MPI_AllReduce(nsb, nsb_rcv, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
            nsb = nsb_rcv
         endif

         write(*,*) "Task ",myproc,' nsb ',nsb,areas%ix_are(1,2)

         if (nsb .gt. areas%ix_are(1,2)+2 .or. nsb .lt. areas%ix_are(1,2)-2) then
            areas%ix_are(1,2) = nsb
            areas%ix_are(2,1) = nsb+1
            call mapare
            call printit_taskX(0,"ix_are changed")
         endif

  end subroutine spherical_neutron_star
end module spherical_neutron_star_mod

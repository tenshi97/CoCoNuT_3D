module setinf_mod

implicit none

contains

!> \par sets up boundary conditions at the interfaces of the different calculation areas
!>
!> \author W. Keil and F. Hanke
!>
!> \detail 
!> This routine sets up boundary conditions at the interfaces of
!> the different calculation areas. \n
!> These areas are used in multi-D-runs to avoid a too restrictive
!> CFL-condition in the innermost core due to very small sizes of the
!> grid zones in the angular directions. Therefore two areas are used, a
!> spherical symmetric (i.e. 1D) core at the first (e.g. 6) zones, which
!> contains only a negligible amount of the total mass, and the actual
!> multi-D-zones. To make use of this tool set up the following values in
!> ppm.par:
!>
!> - are_nu 2
!> - nxa  6
!> - ioya ny
!> - ioza nz
!> - nxb  nx
!> - ioyb 1
!> - iozb 1
!>
!> The next values do not to be considered for 2 areas. If ioya/ioza is set
!> to the number of y- and z-zones ny/nz, a 1D calculation is done in this
!> area; a multi-D calculation is performed, if ioyb/iozb equals 1. "nxa"
!> determines the last zone of the first area (recommend value: 6) and "nxb"
!> the last zone of the second area, i.e. nx for 2 areas. \n
!> The current version supports only 2 areas.
!> As stated at the beginning this subroutine handles the boundary conditions
!> at the interfaces of two areas. \n
!> Therefore the "left/inner" boundary and the "right/outer" boundary have to
!> be treated for all possible interfaces. Two areas have only one interface,
!> of course. The innermost 1D area has a outer boundary and the multi-D area
!> an inner boundary. In VERTEX 8 "ghost-zones", four on each side, are used
!> to handle the boundary. At the inner boundary of the second area (that is
!> the first case, where ia greater 1 and ioy and ioz are equal 1) the four ghost
!> zones, which are next to the boundary, but already in the first area, have
!> to be filled with the appropiate values of the first area. Since this is a
!> 1D area, the 1D properties can just be copied into the multi-D 
!> "ghost zones", for example:
!> \f[
!>  \mbox{den\_bi}(\mbox{ibi},j,k) = \mbox{dentot}(\mbox{ito},1,1)
!>  \;\; \mbox{for all} \; j,k
!> \f]
!> "ibi" determines the number of the ghost zones and "ito" the according
!> zone number at the adjacent area. \n
!> At the outer boundary of the first area
!> (that is the second case, where ia is lower than are_nu and ny and nz are equal
!> 1) the "ghost zones" of the 1D area have to be filled with the appropriate
!> mass and volume weighted values of the multi-D area, respectively, since 
!> there is one 1D zone and ny/nz corresponding angular zoens. For
!> example, the velocity is mass weighted:
!> \f[ \mbox{vex\_bi}(\mbox{ibi},1,1) = \frac{\sum_{j,k} \mbox{vextot}
!> (\mbox{ito},j,k) \mbox{dentot}(\mbox{ito},j,k) \mbox{dvytot}(j) 
!> \mbox{dvztot}(k)}{\sum_{j,k}\mbox{dentot}(\mbox{ito},j,k) 
!> \mbox{dvytot}(j) \mbox{dvztot}(k)}
!> \f]
!> and den_bi is volume weighted:
!> \f[ \mbox{den\_bi}(\mbox{ibi},1,1) = \frac{\sum_{j,k} \mbox{dentot}
!> (\mbox{ito},j,k) \mbox{dvytot}(j) \mbox{dvztot}(k)}
!> {(\sum_{j}\mbox{dvytot}(j)) (\sum_{k}\mbox{dvztot}(k))}
!> \f]
!> \verbatim
!>   SVN - Information
!>   $Revision:$
!>   $Date:$
!> \endverbatim
subroutine setinf(selftime, childrentime)

  use precision
  use abort
      
  use bndinf_hy
  use totare_hy
  ! use arecon_hy
  use intgrs_hy

  use mo_mpi
#ifndef DEBUG_TIMINGS
  use cputim
#endif

  use hydro_areas_mod
  use configure
  use print_stdout_mod
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2)

  integer(kind=ik)           :: ixi, ixf, iox, iyi, iyf, ioy, izi, &
                                izf, ioz, ia, i, j, k, ibi, ito, n

  real(kind=rk)              :: buf(4), buf_rcv(4)
  integer(kind=ik)           :: ierr

  real(kind=rk) :: time_communication_start(2), time_communication_end(2)
!-----------------------------------------------------------------------
!     search for all interfaces:
!-----------------------------------------------------------------------
  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  do ia = 1, areas%are_nu

     ixi = areas%ix_are(ia, 1)
     ixf = areas%ix_are(ia, 2)
     iox = areas%ix_are(ia, 3)
     iyi = areas%ix_are(ia, 4)
     iyf = areas%ix_are(ia, 5)
     ioy = areas%ix_are(ia, 6)
     izi = areas%ix_are(ia, 7)
     izf = areas%ix_are(ia, 8)
     ioz = areas%ix_are(ia, 9)
         
!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------

     areas%nz = izf - izi
     areas%ny = iyf - iyi
     areas%nx = ixf - ixi

     areas%nz = areas%nz/ioz + 1
     areas%ny = areas%ny/ioy + 1
     areas%nx = areas%nx/iox + 1

     if(config%itstp .ne. 0) then
        if(mod(nstep,config%itstp) .eq. 0) then
           call printit_taskX(0,' setinf> ia = ',ia,' ixi = ',ixi,' ixf = ',ixf)
        endif
     endif

!-----------------------------------------------------------------------
!     set up "left" (= "inner") boundary:
!-----------------------------------------------------------------------

         if(ia .gt. 1 .and. ioy .eq. 1 .and. ioz .eq. 1) then
            areas%ix_are(ia,12) = 5   ! -> set bndmnx to "special boundary"

              do i=1,4
                ito = ixi - 5 + i
                ibi = (ia - 1) * 8 + i
                do k=qz_s,qz_e
                  do j=qy_s,qy_e
                    vex_bi(ibi, j, k) = vextot(ito,qy_s,qz_s)
                    vey_bi(ibi, j, k) = veytot(ito,qy_s,qz_s)
                    vez_bi(ibi, j, k) = veztot(ito,qy_s,qz_s)
                    den_bi(ibi, j, k) = dentot(ito,qy_s,qz_s)
                    ene_bi(ibi, j, k) = enetot(ito,qy_s,qz_s)
                    gra_bi(ibi, j, k) = gpotot(ito,qy_s,qz_s)
                    gac_bi(ibi, j, k) = gactot(ito,qy_s,qz_s)
                    gae_bi(ibi, j, k) = gaetot(ito,qy_s,qz_s)
                    pre_bi(ibi, j, k) = pretot(ito,qy_s,qz_s)
                    do n = 1, config%qn
                      xnu_bi(ibi, j, k,n) = xnutot(ito,qy_s,qz_s,n)
                    end do
                  end do
                end do
                dxx_bi(ibi)       = xzrtot(ito) - xzltot(ito)
                ugr_bi(ibi)       = ugrtot(ito)
              end do

        endif

!-----------------------------------------------------------------------
!     set up "right" (= "outer") boundary:
!-----------------------------------------------------------------------

         if(ia .lt. areas%are_nu .and. areas%ny .eq.1 .and. areas%nz .eq.1) then
            areas%ix_are(ia,13) = 5   ! -> set bndmxx to "special boundary"

              do i = 1, 4
                ito = ixf + i
                ibi = ia * 8 + i + 4
                vex_bi(ibi,qy_s,qz_s)=0.0_rk
                vey_bi(ibi,qy_s,qz_s)=0.0_rk
                vez_bi(ibi,qy_s,qz_s)=0.0_rk
                den_bi(ibi,qy_s,qz_s)=0.0_rk
                ene_bi(ibi,qy_s,qz_s)=0.0_rk
                gra_bi(ibi,qy_s,qz_s)=0.0_rk
                gac_bi(ibi,qy_s,qz_s)=0.0_rk
                gae_bi(ibi,qy_s,qz_s)=0.0_rk
                pre_bi(ibi,qy_s,qz_s)=0.0_rk
                xnu_bi(ibi,qy_s,qz_s,1:config%qn)=0.0_rk

                do k = qz_s,qz_e
                  do j = qy_s, qy_e
                    vex_bi(ibi, qy_s, qz_s) = vex_bi(ibi, qy_s, qz_s) + &
                      vextot(ito,j,k)*dentot(ito,j,k)*dvytot(j)*dvztot(k)
                    vey_bi(ibi, qy_s, qz_s) = vey_bi(ibi, qy_s, qz_s) + &
                      veytot(ito,j,k)*dentot(ito,j,k)*dvytot(j)*dvztot(k)
                    vez_bi(ibi, qy_s, qz_s) = vez_bi(ibi, qy_s, qz_s) + &
                      veztot(ito,j,k)*dentot(ito,j,k)*dvytot(j)*dvztot(k)
                    den_bi(ibi, qy_s, qz_s) = den_bi(ibi, qy_s, qz_s) + &
                      dentot(ito,j,k)                *dvytot(j)*dvztot(k)
                    ene_bi(ibi, qy_s, qz_s) = ene_bi(ibi, qy_s, qz_s) + &
                      enetot(ito,j,k)*dentot(ito,j,k)*dvytot(j)*dvztot(k)
                    gra_bi(ibi, qy_s, qz_s) = gra_bi(ibi, qy_s, qz_s) + &
                      0.5_rk*(gpotot(ito,j,k)+gpotot(ito,j-1,k))* &
                      dentot(ito,j,k)                *dvytot(j)*dvztot(k)
                    gac_bi(ibi, qy_s, qz_s) = gac_bi(ibi, qy_s, qz_s) + &
                      gactot(ito,j,k)                *dvytot(j)*dvztot(k)
                    gae_bi(ibi, qy_s, qz_s) = gae_bi(ibi, qy_s, qz_s) + &
                      gaetot(ito,j,k)                *dvytot(j)*dvztot(k)
                    pre_bi(ibi, qy_s, qz_s) = pre_bi(ibi, qy_s, qz_s) + &
                      pretot(ito,j,k)                *dvytot(j)*dvztot(k)

                    do n = 1, config%qn
                      xnu_bi(ibi, qy_s, qz_s,n) = xnu_bi(ibi, qy_s, qz_s,n) + &
                        xnutot(ito,j,k,n)*dentot(ito,j,k)*dvytot(j)*dvztot(k)
                    end do
                  enddo
                enddo

                dxx_bi(ibi)       = xzrtot(ito) - xzltot(ito)
                ugr_bi(ibi)       = ugrtot(ito)

              enddo


!     MPI-AllReduce (Addition) fuer vex_bi,vey_bi,...,pre_bi

              if (use_mpi) then
#ifndef DEBUG_TIMINGS
                 call second_v(time_communication_start)
#endif

                 buf(:) = vex_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                 vex_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = vey_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                 
                 vey_bi(13:16,qy_s,qz_s) = buf_rcv(:)


                 buf(:) = vez_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4,  MPI_DOUBLE_PRECISION, &
                                   MPI_SUM, MPI_COMM_WORLD, ierr)


                 
                 vez_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = den_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

         
                 den_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = ene_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

              
                 ene_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = gra_bi(13:16,qy_s,qz_s)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                 buf(:)=buf_rcv(:)
                 gra_bi(13:16,qy_s,qz_s) = buf(:)

                 buf(:) = gac_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                
                 gac_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = gae_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                 
                 gae_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 buf(:) = pre_bi(13:16,qy_s,qz_s)
                 buf_rcv(:)=buf(:)
                 call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, MPI_COMM_WORLD, ierr)

                
                 pre_bi(13:16,qy_s,qz_s) = buf_rcv(:)

                 do n = 1, config%qn
                    buf(:) = xnu_bi(13:16,qy_s,qz_s,n)
                    buf_rcv(:)=buf(:)
                    call MPI_allreduce(buf, buf_rcv, 4, MPI_DOUBLE_PRECISION, &
                                       MPI_SUM, MPI_COMM_WORLD, ierr)

                   
                    xnu_bi(13:16,qy_s,qz_s,n) = buf_rcv(:)
                 end do
#ifndef DEBUG_TIMINGS
                 call second_v(time_communication_end)
                 timer%hydro_comm = timer%hydro_comm + (time_communication_end - time_communication_start)
                 childrentime = childrentime + (time_communication_end - time_communication_start)
#endif
              endif ! use_mpi


              do i=1,4
                ito = ixf + i
                ibi = ia * 8 + i + 4
                vex_bi(ibi, qy_s, qz_s) = vex_bi(ibi, qy_s, qz_s)/ &
                  den_bi(ibi,qy_s,qz_s)
                vey_bi(ibi, qy_s, qz_s) = vey_bi(ibi, qy_s, qz_s)/ &
                  den_bi(ibi,qy_s,qz_s)
                vez_bi(ibi, qy_s, qz_s) = vez_bi(ibi, qy_s, qz_s)/ &
                  den_bi(ibi,qy_s,qz_s)
                ene_bi(ibi, qy_s, qz_s) = ene_bi(ibi, qy_s, qz_s)/ &
                  den_bi(ibi,qy_s,qz_s)
                gra_bi(ibi, qy_s, qz_s) = gra_bi(ibi, qy_s, qz_s)/ &
                          den_bi(ibi,qy_s,qz_s)

                do n=1,config%qn
                  xnu_bi(ibi, qy_s, qz_s,n) = xnu_bi(ibi, qy_s, qz_s,n)/ &
                                                 den_bi(ibi,qy_s,qz_s)
                end do

                den_bi(ibi, qy_s, qz_s) = den_bi(ibi, qy_s, qz_s)/ &
                  (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
                gac_bi(ibi, qy_s, qz_s) = gac_bi(ibi, qy_s, qz_s)/ &
                  (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
                gae_bi(ibi, qy_s, qz_s) = gae_bi(ibi, qy_s, qz_s)/ &
                  (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
                pre_bi(ibi, qy_s, qz_s) = pre_bi(ibi, qy_s, qz_s)/ &
                  (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
              end do

         end if

      enddo   !end of ia-loop

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return

    end subroutine setinf

end module setinf_mod

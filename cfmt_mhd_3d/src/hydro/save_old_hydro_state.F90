module save_old_hydro_state_mod

  public

contains

  subroutine save_old_hydro_state(dehtot, tehtot, yehtot, enhtot, &
                                      selftime, childrentime)

    use precision
    use cputim
#ifdef CFC_TRANSPORT
    use size_cfc
    use hydro_primitives_cfc
    use phycon, only: pc_geog
#endif
    use mo_mpi
    use configure
    use totare_hy
    implicit none

    integer(kind=ik)            :: k, j, jk
    real(kind=rk), intent(out)  :: selftime(2), childrentime(2)
    real(kind=rk), dimension(2) :: selftime_start(2)
    real(kind=rk), dimension(2) :: paralleltime_start, paralleltime_end

    real(kind=rk), intent(out), &
    dimension(config%qx,qy_s:qy_e,qz_s:qz_e) :: tehtot, yehtot, enhtot, dehtot

    selftime          = 0._rk
    childrentime      = 0._rk
    paralleltime_start = 0._rk
    paralleltime_end   = 0._rk
#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif



    ! save density temperature, electron fraction and internal energy:
#ifdef CFC_TRANSPORT
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
    call second_v(paralleltime_start)
#endif
!$OMP PARALLEL DO &
!$OMP SCHEDULE(static) &
!$OMP PRIVATE(j,k,jk) &
!$OMP SHARED(n_loc,o_loc,n_s,o_s,dehtot,yehtot,tehtot,enhtot,rho,xnnu,t,eps)
#endif
    do jk = 1, n_loc * o_loc
       k = int((jk + n_loc - 1) / n_loc)
       j = (n_s - 1) + (jk - (k - 1) * n_loc)
       k = k + o_s - 1

       dehtot(1:config%qx,j,k) = rho (1:config%qx,j,k)*pc_geog
       yehtot(1:config%qx,j,k) = xnnu(1:config%qx,j,k,config%qn)
       tehtot(1:config%qx,j,k) = t   (1:config%qx,j,k)
       enhtot(1:config%qx,j,k) = eps (1:config%qx,j,k)
    end do
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
    call second_v(paralleltime_end)
#endif
#endif
    timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start

#else /* CFC_TRANSPORT */
    dehtot(:,:,:) = dentot(:,:,:)
    yehtot(:,:,:) = xnutot(:,:,:,config%qn)
    tehtot(:,:,:) = temtot(:,:,:)
    enhtot(:,:,:) = enetot(:,:,:)-0.5_rk*(vextot(:,:,:)**2+veytot(:,:,:)**2+ &
                    veztot(:,:,:)**2)
#endif /* CFC_TRANSPORT */

#ifndef DEBUG_TIMINGS
    call second_v(selftime)
    selftime = selftime - selftime_start
#endif
  end subroutine save_old_hydro_state

end module save_old_hydro_state_mod

module mod_enclosed_mass
  use precision

  implicit none

  contains

  !> Calculate the enclosed baryonic mass on zone centers
  function enclosed_baryonic_mass_c() result(m_enc)
    use abort
    use gfloat_hy, only : vlfrac
    use totare_hy, only : xzntot, yzntot, zzntot, &
                          xzltot, yzltot, zzltot, &
                          xzrtot, yzrtot, zzrtot, &
                          dvytot, dvztot, dentot
    use configure
    use phycon
    use mo_mpi

    implicit none

    integer(kind=ik) :: i, j, k
    real(kind=rk) :: dsurf
    real(kind=rk) :: m_enc(1:config%nx)
    real(kind=rk), dimension(config%nx) :: dvx_l, dvx_r

    abort_if((config%ny .ne. 1) .or. (config%nz .ne. 1), "Only tested in 1D!")
    if (use_mpi) then
      raise_abort("MPI not implemented here yet")
    endif

    m_enc(:) = 0.0_rk

    dvx_l = (xzntot(1:config%nx)**3 - xzltot(1:config%nx)**3) / 3._rk
    dvx_r = (xzrtot(1:config%nx)**3 - xzntot(1:config%nx)**3) / 3._rk

    m_enc(1) = config%pmass
    do k = qz_s, qz_e
      do j = qy_s, qy_e
        dsurf = vlfrac * dvztot(k) * dvytot(j)
        m_enc(1) = m_enc(1) + dsurf * dentot(1,j,k) * dvx_l(1)
      end do
    end do
    do i = 2, config%nx
      m_enc(i) = m_enc(i-1)
      do k = qz_s, qz_e
        do j = qy_s, qy_e
          dsurf = vlfrac * dvztot(k) * dvytot(j)

          m_enc(i) = m_enc(i) + &
                     dsurf * (dentot(i  ,j,k) * dvx_l(i) + &
                              dentot(i-1,j,k) * dvx_r(i-1))
        enddo
      enddo
    enddo

    m_enc = m_enc / pc_ms

  end function

  !> Calculate the enclosed baryonic mass on zone boundaries
  function enclosed_baryonic_mass_b() result(m_enc)
    use abort
    use gfloat_hy, only : vlfrac
    use totare_hy, only : xzntot, yzntot, zzntot, &
                          xzltot, yzltot, zzltot, &
                          xzrtot, yzrtot, zzrtot, &
                          dvytot, dvztot, dentot
    use configure
    use phycon
    use mo_mpi

    implicit none

    integer(kind=ik) :: i, j, k
    real(kind=rk) :: dsurf
    real(kind=rk) :: m_enc(0:config%nx)
    real(kind=rk), dimension(config%nx) :: dvx

    abort_if((config%ny .ne. 1) .or. (config%nz .ne. 1), "Only tested in 1D!")
    if (use_mpi) then
      raise_abort("MPI not implemented here yet")
    endif

    m_enc(:) = 0.0_rk

    dvx = (xzrtot(1:config%nx)**3 - xzltot(1:config%nx)**3) / 3._rk

    m_enc(0) = config%pmass
    do i = 1, config%nx
      m_enc(i) = m_enc(i-1)
      do k = qz_s, qz_e
        do j = qy_s, qy_e
          dsurf = vlfrac * dvztot(k) * dvytot(j)
          m_enc(i) = m_enc(i) + dsurf * dvx(i) * dentot(i,j,k)
        enddo
      enddo
    enddo

    m_enc = m_enc / pc_ms

  end function
end module

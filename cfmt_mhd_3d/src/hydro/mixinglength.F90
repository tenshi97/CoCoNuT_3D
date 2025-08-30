#if defined(CONVECTION_1D) && !defined(CFC_TRANSPORT)
#ifndef ENECONS
! needs grav which is only recalculated when ENECONS is active
#error "You have to specify -DENECONS along with CONEVECTION_1D"
#endif

#undef DEBUG

module mixinglength
    implicit none

    private
    public mixinglengthconvection

    contains

    ! finite difference of df/dx, centered, first order
    function diff(f, x) result(df_dx)
      use precision
      real(kind=rk), intent(in)                     :: f(:)
      real(kind=rk), intent(in), dimension(size(f)) :: x
      real(kind=rk), dimension(size(f) - 1) :: df_dx
      integer(kind=ik) :: n

      n = size(f)
      df_dx = (f(2:n) - f(1:n-1)) / (x(2:n) - x(1:n-1))
    end function

    !> Adds fluxes due to mixing-length approximated
    !> convection treatment in 1D simulations
    !>
    !> \param alpha     O(1) parameter relating the "mixing length" to the
    !>                  pressure scale height P/(dP/dr)
    !>
    !> \param r         radius of zone centers [cm]
    !> \param p         pressure at zone centers [erg/cm^3]
    !> \param p_b       pressure at zone interfaces, size(p_b) == size(p_c) - 1, [erg/cm^3]
    !> \param rho       density at zone centers [g/cm^3]
    !> \param rho_b     density at zone interfaces [g/cm^3]
    !> \param gamma     adiabatic index at zone centers [1]
    !> \param gamma_b   adiabatic index at zone interfaces [1]
    !> \param phi       gravitational potential at zone centers [erg/g]
    !> \param v         x-velocity at zone centers [cm/s]
    !>
    !> \param x         mass fractions at zone centers [1]
    !> \param xflx      partial density flux at zone interfaces [g/cm^2/s]
    !> \param e         mass-specific energy density at zone centers [erg/g]
    !> \param eflx      energy flux at zone interfaces [erg/cm^2/s]

    subroutine mixinglengthconvection(alpha, r, p, p_b, rho, rho_b, gamma, gamma_b, phi, v, x, xflx, e, eflx)
      use precision
      implicit none

      real(kind=rk), intent(in)    :: alpha
      real(kind=rk), intent(in)    :: r(:)
      real(kind=rk), intent(in)    :: p(:), phi(:), rho(:), e(:), gamma(:), v(:)
      real(kind=rk), intent(in)    :: p_b(:), rho_b(:), gamma_b(:)
      real(kind=rk), intent(in)    :: x(:, :)
      real(kind=rk), intent(inout) :: eflx(:), xflx(:, :)

      real(kind=rk), dimension(size(r)-1) :: lmix, vmix, dp_dr, dp_drho_s, vsound
      integer(kind=rk) :: i, i_sonic
#ifdef DEBUG
      integer(kind=rk) :: j
      real(kind=rk), dimension(size(eflx)) :: dflux
#endif

#ifdef DEBUG

      debug(size(r))
      debug(size(p))
      debug(size(p_b))
      debug(size(rho))
      debug(size(rho_b))
      debug(size(gamma))
      debug(size(gamma_b))
      debug(size(phi))
      debug(shape(x))
      debug(shape(xflx))
      debug(size(e))
      debug(size(eflx))

      write(*,*) "Mixing input [nr, r, p, p_b, rho, rho_b, e, gamma, gamma_b, phi]"
      do i = 1, size(r)
        write(*,'(i4,x,9(x,1pe12.5))') i, r(i), p(i), p_b(i), rho(i), rho_b(i), e(i), gamma(i), gamma_b(i), phi(i)
      enddo
      write(*,*)
#endif

      dp_dr = diff(p, r)
      lmix = abs(alpha * p_b / dp_dr)
      dp_drho_s = gamma_b * p_b / rho_b
      vmix = sqrt(max(diff(phi, r) * lmix**2 / (2 * rho_b) * max(diff(rho, r) - dp_dr / dp_drho_s, 0.0_rk), 0.0_rk))

      ! Limit the fluxes according to CFL, the "right" way would probably
      ! be to limit the timestep according to min(vsound, vmix * lmix / delta_r)
      vsound = sqrt(gamma_b * p_b / rho_b)
      vmix(:) = min(vmix(:), vsound * (r(2:size(r)) - r(1:size(r)-1)) / lmix)

#ifdef DEBUG
      write(*,*) "Mixing: [i:i+1, lmix, vmix, dp_dr, dp_drho_s, vsound]"
      do i = 1, size(r) - 1
        write(*,'(i4,'':'',i4,5(x,1pe12.5))') i, i+1, lmix(i), vmix(i), dp_dr(i), dp_drho_s(i), vsound(i)
      enddo
      write(*,*)
#endif

      ! do it only within sub-sonic region, note that
      ! v and vsound have different shape and are defined
      ! at different locations
      i_sonic = size(vsound) + 1
      do i = 1, size(vsound)
        if (abs(0.5 * (v(i) + v(i+1))) .gt. vsound(i)) then
          i_sonic = i
          exit
        endif
      enddo


#ifdef DEBUG
      debug(i_sonic)
#endif

      if (i_sonic .le. size(vsound)) then
        vmix(i_sonic-1:) = 0.0_rk
      endif

#ifdef DEBUG
      write(*,*) "Mixing: [i:i+1, lmix, vmix, dp_dr, dp_drho_s, vsound]"
      do i = 1, size(r) - 1
        write(*,'(i4,'':'',i4,5(x,1pe12.5))') i, i+1, lmix(i), vmix(i), dp_dr(i), dp_drho_s(i), vsound(i)
      enddo
      write(*,*)
#endif

      ! update fluxes
      do i = 1, size(xflx, dim=2)
        xflx(:,i) = xflx(:, i) - 0.5_rk * vmix * rho_b * lmix * diff(x(:,i), r)
#ifdef DEBUG
        dflux =  - 0.5_rk * vmix * rho_b * lmix * diff(x(:,i), r)
        write(*,*) "xflx(:,", i, ")"
        do j = 1, size(xflx, dim=1)
          write(*,'(i4,'':'',i4,1(x,1pe12.5))') j, j + 1, dflux(j)
        enddo
        write(*,*)
#endif
      enddo
      eflx = eflx - 0.5_rk * vmix * rho_b * lmix * (diff(e - 0.5 * v**2, r) + p_b * diff(1.0_rk/rho, r))
#ifdef DEBUG
      dflux =  - 0.5_rk * vmix * rho_b * lmix * (diff(e - 0.5 * v**2, r) + p_b * diff(1.0_rk/rho, r))
      write(*,*) "eflx"
      do i = 1, size(eflx)
        write(*,'(i4,'':'',i4,1(x,1pe12.5))') i, i + 1, dflux(i)
      enddo
      write(*,*)
#endif
    end subroutine

  end module
#endif /* CONVECTION_1D */

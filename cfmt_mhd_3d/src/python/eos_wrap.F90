program eos_python
end program

module eos_wrapper
  implicit none

  contains

  subroutine init(nspecies, low_den_nse_eos, tj_ls_rho, tj_wolff_rho, wolff_ls_rho)
    use precision
    use eos_sn2, only: init_eos
    use configure
    use abort

    use configure

    implicit none
    integer(kind=ik), intent(in) :: nspecies, low_den_nse_eos
    real(kind=rk), intent(in), optional :: tj_ls_rho, tj_wolff_rho, wolff_ls_rho

    config%qn = nspecies
    config%low_den_nse_eos = low_den_nse_eos

    if (tj_ls_rho == 0.0) then
      config%tj_ls_rho = -1.0
    else
      config%tj_ls_rho = tj_ls_rho
    endif

    if (wolff_ls_rho == 0.0) then
      config%wolff_ls_rho = -1.0
    else
      config%wolff_ls_rho = wolff_ls_rho
    endif

    if (tj_wolff_rho == 0.0) then
      config%tj_wolff_rho = -1.0
    else
      config%tj_wolff_rho = tj_wolff_rho
    endif

    config%eos_sw = 0

    call init_eos

  end subroutine

  subroutine eos(n, nspecies, rho, rho_in, tem, tem_in, xi, xi_in, xhrep, za,   ede, ede_in, p, gam,  s, ccu, cce, ccn, ccp, mode,nsemode,error)
    use precision
    use eos_sn2 , only: eos_actual => eos
    use configure
    implicit none

    integer :: n, nspecies
    integer(kind=ik), intent(in)    :: mode,nsemode
    real(kind=rk),    intent(in)    :: rho_in(n), tem_in(n), xi_in(n,nspecies), ede_in(n)
    real(kind=rk),    intent(out)   :: rho(n), tem(n), ede(n), xi(n,nspecies), &
                                       p(n), gam(n), s(n), za(n,2), &
                                       xhrep(n), ccu(n), cce(n),    &
                                       ccn(n), ccp(n)
    logical, intent(out) :: error
    real(kind=rk) :: selftime(2), childtime(2)

    ede = ede_in
    tem = tem_in
    rho = rho_in
    xi  = xi_in

    if ((size(xi_in, 1) .ne. size(rho_in)) .or. (size(xi_in, 2) .ne. config%qn)) then
      write(0,*) "Error: shape of xi has to be (", size(rho_in), ", ", config%qn, ")"
      error = .false.
    else
      call eos_actual(rho, tem, xi, xhrep, za, ede, p, gam, s, ccu, cce, ccn, ccp, selftime, childtime, mode, nsemode, error)
    endif

    if (error) then
      write(0,*) "Error evaluating EoS"
      p = 0.0
      gam = 0.0
      s = 0.0
      za = 0.0
      xhrep = 0.0
      ccu = 0.0
      cce = 0.0
      ccn = 0.0
      ccp = 0.0
    endif

  end subroutine

end module

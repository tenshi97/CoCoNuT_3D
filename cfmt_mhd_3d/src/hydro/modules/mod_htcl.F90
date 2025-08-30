!> defines the parameters for the FT model
module htcl_hy

    use precision
    use phycon

    !> luminosity parameter
    real ( kind = rk ) :: htcl_lum

    !> opacities
    real(kind=rk) :: k_abs = 1._rk

    !> proton and neutron fraction and optical depth
    real(kind=rk), allocatable :: y_n(:),y_p(:),f_abs(:,:,:),f_abs_a(:,:,:)

    contains

    subroutine init_htcl
        !! set the parameters of the model
        htcl_lum    = 2.8_rk ! in 1e52 erg/s
    end subroutine init_htcl

    subroutine allocate_htcl_hy(mem_global)

      use precision
      use abort
      use mo_mpi
      use print_stdout_mod, only : print_memory_alloc

      use configure
      implicit none

      integer(kind=ik) :: istat
      real(kind=rk)    :: mem_global, mem_local

      allocate(y_n(config%q),y_p(config%q), stat=istat)

      allocate(f_abs(config%qx,qy_s:qy_e,qz_s:qz_e),    &
               f_abs_a(config%q,qy_s:qy_e,qz_s:qz_e), stat=istat)

      if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 1 failed")
      end if

      mem_local = (qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk * (config%qx + config%q) + config%q*8._rk

      mem_local = mem_local/1024._rk/1024._rk

      mem_global = mem_global + mem_local

      call print_memory_alloc(mem_local, mem_global, "htcl_hy")

    end subroutine allocate_htcl_hy

    subroutine deallocate_htcl_hy

      use precision
      use abort
      use mo_mpi

      implicit none

      integer(kind=ik) :: istat

      deallocate(y_n,y_p, stat=istat)

      deallocate(f_abs,f_abs_a, stat=istat)

      if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
      end if

    end subroutine deallocate_htcl_hy

  end module htcl_hy

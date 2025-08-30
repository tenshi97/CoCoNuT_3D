module setup_c
  use precision
  use, intrinsic :: iso_c_binding

  implicit none
  public

  interface
    subroutine version_short() bind(C, name="version_short")
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine
  end interface

  interface
    subroutine version_make_inc_config() bind(C, name="version_make_inc_config")
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine
  end interface
#ifdef WRITE_BINARY_OUTPUT
#ifdef WITH_HDF5
  interface
    subroutine version_write_hdf5(loc) bind(C, name="version_write_hdf5")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: loc
    end subroutine
  end interface
#endif
#endif

  interface
    subroutine unlimit_stack() bind(C, name="unlimit_stack")
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine
  end interface

  interface
    integer(kind=C_INT) function create_directories_c() bind(C, name="create_directories")
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface

  interface
    subroutine redirect_stdout_c(myproc) bind(C, name="redirect_stdout")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(in) :: myproc
    end subroutine
  end interface

  interface
    subroutine get_process_id_c(process_id, pprocess_id) bind(C, name="get_process_id")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(out) :: process_id, pprocess_id
    end subroutine
  end interface

#ifdef CHECK_THREAD_AFFINITY
  interface
    subroutine get_thread_affinity_c(cpu_id) bind(C, name="get_thread_affinity")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(out) :: cpu_id
    end subroutine
  end interface
  
  interface
    subroutine get_process_affinity_c(cpu_id) bind(C, name="get_process_affinity")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(out) :: cpu_id
    end subroutine
  end interface
#endif

  interface
    subroutine ms_since_epoch_c(ms) bind(C, name="ms_since_epoch")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT64_T), intent(out) :: ms
    end subroutine
  end interface

  interface
    subroutine exit(status) bind(C, name="exit")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(in), value :: status
    end subroutine
  end interface

  contains

    function create_directories() result(res)
      implicit none
      integer(kind=ik) :: res
      res = int(create_directories_c(), kind=ik)
    end function

    subroutine redirect_stdout(myproc)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=ik), intent(in) :: myproc
      call redirect_stdout_c(int(myproc, kind=C_INT))
    end subroutine
    
#ifdef CHECK_THREAD_AFFINITY
    subroutine get_thread_affinity(cpu_id)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=ik), intent(out) :: cpu_id
      integer(kind=C_INT) :: cpu_id_c
      call get_thread_affinity_c(cpu_id_c)
      cpu_id = int(cpu_id_c, kind=ik)
    end subroutine

    subroutine get_process_affinity(cpu_id)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=ik), intent(out) :: cpu_id
      integer(kind=C_INT) :: cpu_id_c
      call get_process_affinity_c(cpu_id_c)
      cpu_id = int(cpu_id_c, kind=ik)
    end subroutine
#endif

    subroutine get_process_id(process_id, pprocess_id)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=ik), intent(out) :: process_id, pprocess_id
      integer(kind=C_INT) :: process_id_c, pprocess_id_c

      call get_process_id_c(process_id_c, pprocess_id_c)

      process_id  = int(process_id_c,  kind=ik)
      pprocess_id = int(pprocess_id_c, kind=ik)
    end subroutine
    
    function ms_since_epoch() result(ms)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=8) :: ms
      integer(kind=C_INT64_T) :: ms_c
      call ms_since_epoch_c(ms_c)
      ms = int(ms_c, kind=8)
    end function

end module setup_c

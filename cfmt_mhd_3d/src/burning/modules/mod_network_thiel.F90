module netw_thiel

  use precision
  ! use abort


  use configure

  implicit none
  private
  public ye_netw, xi_dummy_netw, ye_save, xnuc_dmy_save, &
         l_renorm_xnuc, xnuc_save, allocate_network

  real(kind=rk)               :: ye_netw
  real(kind=rk) , allocatable :: ye_save(:,:,:) ! save ye of the  zones
  real(kind=rk) , allocatable :: xnuc_dmy_save(:,:,:,:) ! save xnuc of the 2 dummy nuclei
 
  real(kind=rk), allocatable  :: xi_dummy_netw(:)

  integer(kind=ik) :: l_renorm_xnuc


  real(kind=rk), allocatable  :: xnuc_save(:,:,:,:)
   
contains

  subroutine allocate_network( mem_global)
    use precision
    use abort
    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc

    
    use configure

    implicit none 
    
    integer(kind=ik)        :: istat
    real(kind=rk)           :: mem_global, mem_local

    allocate(ye_save(config%qx,config%qy,config%qz) ,           &
             xnuc_dmy_save(config%qx,config%qy,config%qz,2),    &
             xi_dummy_netw(2),                                  &
             xnuc_save(config%qx,config%qy,config%qz,config%qn), stat=istat)

 
    if (istat .ne. 0) then
        raise_abort("allocate_network(): allocation of module netw_thiel")
    end if

    mem_local = config%qx*config%qy*config%qz*3*8._rk +         &
                xi_dummy_netw(2)*8._rk                +         &
                config%qx*config%qy*config%qz*config%qn*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "totare_hy")

  end subroutine allocate_network
end module netw_thiel

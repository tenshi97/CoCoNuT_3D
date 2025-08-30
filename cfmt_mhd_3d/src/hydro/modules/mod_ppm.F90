!>
!> \par This module provides some floating-point constants
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module fpconst_hy 
  use precision

  implicit none
! LOCAL variables that are not in modules
  SAVE

  real(kind=rk), parameter :: zero  = 0.0_rk
  real(kind=rk), parameter :: half  = 0.5_rk 
  real(kind=rk), parameter :: qrtr  = 0.25_rk  
  real(kind=rk), parameter :: one   = 1.0_rk 
  real(kind=rk), parameter :: two   = 2.0_rk   
  real(kind=rk), parameter :: three = 3.0_rk  
  real(kind=rk), parameter :: four  = 4.0_rk   
  real(kind=rk), parameter :: six   = 6.0_rk   
      
  real(kind=rk), parameter :: c1by3 = 1.0_rk/3.0_rk   
  real(kind=rk), parameter :: c2by3 = 2.0_rk/3.0_rk   
  real(kind=rk), parameter :: c4by3 = 4.0_rk/3.0_rk   
  
end module fpconst_hy

!>
!> \par This module provides the new (in a time-step) hydro variables 
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module vnew_hy 
  use precision
  
      
  implicit none
  SAVE
! LOCAL variables that are not in modules

      

  real(kind=rk), dimension(:,:,:), allocatable ::            &
                 velx, vely, velz, densty, energy, press,    &
                 gammae, gammac, stot, temp

  real(kind=rk), dimension(:,:,:,:), allocatable :: cpo

  real(kind=rk), allocatable :: ugridx(:), xnuc(:,:,:,:),gpot(:,:,:)

  integer(kind=ik) , allocatable :: ishck(:,:,:), ishock_pos(:,:,:)


!     Arrays for ghost zones:
  real(kind=rk), dimension(:,:,:), allocatable :: densty_ub, vely_ub, velx_ub, &
         velz_ub, vxold_ub, vzold_ub, energy_ub, press_ub, temp_ub,     &
         gammae_ub, gammac_ub, gpocntr_ub,densty_lb, vely_lb, velx_lb, &
         velz_lb, vxold_lb, vzold_lb, energy_lb, press_lb, temp_lb,     &
         gammae_lb, gammac_lb, gpocntr_lb

  real(kind=rk), dimension(:,:,:,:), allocatable :: xnuc_ub,xnuc_lb

  real(kind=rk), dimension(:,:,:), allocatable :: densty_kb, vely_kb, velx_kb, &
         velz_kb, vxold_kb, vyold_kb, energy_kb, press_kb, temp_kb,     &
         gammae_kb, gammac_kb, gpocntr_kb,densty_pb, vely_pb, velx_pb, &
         velz_pb, vxold_pb, vyold_pb, energy_pb, press_pb, temp_pb,     &
         gammae_pb, gammac_pb, gpocntr_pb

  real(kind=rk), dimension(:,:,:,:), allocatable :: xnuc_kb,xnuc_pb



contains
!>
!> \par This subroutine allocates the arrays from module vnew_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>   
  subroutine allocate_vnew_hy(mem_global)
      
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure, only : config

    
    implicit none

    
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(velx(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             vely(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             velz(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             densty(config%qx,qy_s:qy_e,qz_s:qz_e),            &
             energy(config%qx,qy_s:qy_e,qz_s:qz_e),            &
             press(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             gammae(config%qx,qy_s:qy_e,qz_s:qz_e),            &
             gammac(config%qx,qy_s:qy_e,qz_s:qz_e),            &
             stot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             temp(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             xnuc(config%qx,qy_s:qy_e,qz_s:qz_e,config%qn),    &  
             gpot(0:config%qx,qy_s-1:qy_e+1,qz_s:qz_e),        &
             ishck(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1), &
             ishock_pos(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 1 failed")
    end if

    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*     8._rk*10._rk + &
                config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%qn * 8._rk +&
                (config%qx+1)*(qy_e-qy_s+3)*(qz_e-qz_s+1)* 8._rk        + &
                (config%qx+2)*(qy_e-qy_s+3)*(qz_e-qz_s+1)* 4._rk*2._rk  
              


    allocate(cpo(config%qx,qy_s:qy_e,qz_s:qz_e,4), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 1 failed")
    end if

    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*4 * 8._rk

    allocate(ugridx(1:config%qx),stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 2 failed")
    end if

    mem_local = mem_local + config%qx* 8._rk

    if (use_mpi) then

       allocate(densty_ub(config%qx,4,qz_s:qz_e),          &
                vely_ub(config%qx,4,qz_s:qz_e),            &
                velx_ub(config%qx,4,qz_s:qz_e),            &
                velz_ub(config%qx,4,qz_s:qz_e),            &
                vxold_ub(config%qx,4,qz_s:qz_e),           &
                vzold_ub(config%qx,4,qz_s-1:qz_e+1),       &
                energy_ub(config%qx,4,qz_s:qz_e),          &
                press_ub(config%qx,4,qz_s:qz_e),           &
                temp_ub(config%qx,4,qz_s:qz_e),            &
                gammae_ub(config%qx,4,qz_s:qz_e),          &
                gammac_ub(config%qx,4,qz_s:qz_e),          &
                gpocntr_ub(config%qx,4,qz_s:qz_e),         &
                densty_lb(config%qx,4,qz_s:qz_e),          &
                vely_lb(config%qx,4,qz_s:qz_e),            &
                velx_lb(config%qx,4,qz_s:qz_e),            &
                velz_lb(config%qx,4,qz_s:qz_e),            &
                vxold_lb(config%qx,4,qz_s:qz_e),           &
                vzold_lb(config%qx,4,qz_s-1:qz_e+1),       &
                energy_lb(config%qx,4,qz_s:qz_e),          &
                press_lb(config%qx,4,qz_s:qz_e),           &
                temp_lb(config% qx,4,qz_s:qz_e),           &
                gammae_lb(config%qx,4,qz_s:qz_e),          &
                gammac_lb(config%qx,4,qz_s:qz_e),          &
                gpocntr_lb(config%qx,4,qz_s:qz_e),stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 2 failed")
       end if

       mem_local = mem_local + config%qx*(qz_e-qz_s+1)*4*8._rk*22._rk
       mem_local = mem_local + config%qx*(qz_e-qz_s+3)*4*8._rk*2._rk

       allocate(xnuc_ub(config%qx,4,qz_s:qz_e,config%qn),    &
                xnuc_lb(config%qx,4,qz_s:qz_e,config%qn),stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 2 failed")
       end if

       mem_local = mem_local + config%qx*(qz_e-qz_s+1)*4*config%qn*8._rk*2._rk


#if !(defined(PROGRAM_remap))

       if (config%nsdim .eq. 3) then
          allocate(densty_kb(config%qx,qy_s:qy_e,4),        &
                   vely_kb(config%qx,qy_s:qy_e,4),          &
                   velx_kb(config%qx,qy_s:qy_e,4),          &
                   velz_kb(config%qx,qy_s:qy_e,4),          &
                   vxold_kb(config%qx,qy_s:qy_e,4),         &
                   vyold_kb(config%qx,qy_s-1:qy_e+1,4),     &
                   energy_kb(config%qx,qy_s:qy_e,4),        &
                   press_kb(config%qx,qy_s:qy_e,4),         &
                   temp_kb(config%qx,qy_s:qy_e,4),          &
                   gammae_kb(config%qx,qy_s:qy_e,4),        &
                   gammac_kb(config%qx,qy_s:qy_e,4),        &
                   gpocntr_kb(config%qx,qy_s:qy_e,4),       &
                   densty_pb(config%qx,qy_s:qy_e,4),        &
                   vely_pb(config%qx,qy_s:qy_e,4),          &
                   velx_pb(config%qx,qy_s:qy_e,4),          &
                   velz_pb(config%qx,qy_s:qy_e,4),          &
                   vxold_pb(config%qx,qy_s:qy_e,4),         &
                   vyold_pb(config%qx,qy_s-1:qy_e+1,4),     &
                   energy_pb(config%qx,qy_s:qy_e,4),        &
                   press_pb(config%qx,qy_s:qy_e,4),         &
                   temp_pb(config%qx,qy_s:qy_e,4),          &
                   gammae_pb(config%qx,qy_s:qy_e,4),        &
                   gammac_pb(config%qx,qy_s:qy_e,4),        &
                   gpocntr_pb(config%qx,qy_s:qy_e,4),stat=istat)
          if (istat .ne. 0) then
             raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 2 failed")
          end if
       
          mem_local=mem_local + config%qx*(qy_e-qy_s+1)*4*8._rk*22._rk
          mem_local = mem_local + config%qx*(qy_e-qy_s+3)*4*8._rk*2._rk

          allocate(xnuc_kb(config%qx,qy_s:qy_e,4,config%qn),   &
                   xnuc_pb(config%qx,qy_s:qy_e,4,config%qn),stat=istat)
          if (istat .ne. 0) then
             raise_abort("allocate_hydro_arrays(): allocation of module vnew_hy 2 failed")
          end if
          mem_local = mem_local + config%qx*(qy_e-qy_s+1)*4*config%qn*8._rk*2._rk
       endif ! config%nsdim
#endif /* PROGRAM_REMAP */

    endif ! use_mpi

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local
    
    call print_memory_alloc(mem_local, mem_global, "vnew_hy")

  end subroutine allocate_vnew_hy
!>
!> \par This subroutine deallocates the arrays from module vnew_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_vnew_hy
      
    use precision
    use abort
    use mo_mpi
    use configure

    implicit none
    
    integer(kind=ik) :: istat

    deallocate(velx, vely, velz, densty, energy, press,       &
               gammae, gammac, stot, temp, xnuc, gpot,        &
               ishck, ishock_pos, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
    end if

    deallocate(cpo, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): dallocation of module vnew_hy 2 failed")
    end if

    deallocate(ugridx,stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 3 failed")
    end if
    
    if (use_mpi) then

       deallocate(densty_ub, vely_ub, velx_ub, &
                  velz_ub, vxold_ub, vzold_ub, energy_ub, press_ub, temp_ub,     &
                  gammae_ub, gammac_ub, gpocntr_ub,densty_lb, vely_lb, velx_lb, &
                  velz_lb, vxold_lb, vzold_lb, energy_lb, press_lb, temp_lb,     &
                  gammae_lb, gammac_lb, gpocntr_lb, stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
       end if
       
       deallocate(xnuc_ub,xnuc_lb, stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
       end if
       
#if !(defined(PROGRAM_remap))
       
       if (config%nsdim .eq. 3) then
          deallocate(densty_kb, vely_kb, velx_kb, velz_kb, vxold_kb,    &
                     vyold_kb, energy_kb, press_kb, temp_kb, gammae_kb, &
                     gammac_kb, gpocntr_kb,densty_pb, vely_pb, velx_pb, &
                     velz_pb, vxold_pb, vyold_pb, energy_pb, press_pb,  &
                     temp_pb, gammae_pb, gammac_pb, gpocntr_pb, stat=istat)
          if (istat .ne. 0) then
             raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
          end if
       
          deallocate(xnuc_kb,xnuc_pb, stat=istat)
          if (istat .ne. 0) then
             raise_abort("allocate_hydro_arrays(): deallocation of module vnew_hy 1 failed")
          end if
       
       endif ! config%nsdim
#endif /* PROGRAM_REMAP */
    endif ! use_mpi
  end subroutine deallocate_vnew_hy

    
end module vnew_hy

! ---------------------------------------------------------------------

!>
!> \par This module provides the old (in a time-step) hydro variables 
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module vold_hy 
  use precision


  implicit none
  SAVE
! LOCAL variables that are not in modules

  real(kind=rk), allocatable, dimension(:,:,:) ::  vxold, vyold, &
                                                   vzold, gold


contains

!>
!> \par This subroutine allocates the arrays from module vold_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_vold_hy(mem_global)
    
    use precision
    use abort

    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc

    use configure 
    implicit none

      
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(vxold(config%qx,qy_s:qy_e,qz_s:qz_e),     &
             vyold(config%qx,qy_s-1:qy_e+1,qz_s:qz_e), &
             vzold(config%qx,qy_s:qy_e,qz_s-1:qz_e+1), &
             gold (0:config%qx,qy_s-1:qy_e,qz_s:qz_e), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module vold_hy 1 failed")
    end if
    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk
    mem_local = mem_local + config%qx*(qy_e-qy_s+3)*(qz_e-qz_s+1)*8._rk
    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+3)*8._rk
    mem_local = mem_local + (config%qx+1)*(qy_e-qy_s+2)*(qz_e-qz_s+1)*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

     call print_memory_alloc(mem_local, mem_global, "vold_hy")

  end subroutine allocate_vold_hy
!>
!> \par This subroutine deallocates the arrays from module vold_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_vold_hy
      
    use precision
    use abort

    use mo_mpi
    
    implicit none
    integer(kind=ik) :: istat

    deallocate(vxold, vyold, vzold, gold, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module vold_hy 1 failed")
    end if


  end subroutine deallocate_vold_hy

end module vold_hy
 
! ---------------------------------------------------------------------

!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module mesh_hy 
  use precision
  
  implicit none
      
! LOCAL variables that are not in modules

  save

  real(kind=rk), allocatable, dimension(:) :: xznl, xzn,  xznr,  dvx, &
                                              yznl, yzn,  yznr,  dvy, &
                                              zznl, zzn,  zznr,  dvz 

contains
!>
!> \par This subroutine allocates the arrays from module mesh_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_mesh_hy(mem_global)
    
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
    
    allocate(xznl(config%q_nqx), xzn(config%q_nqx),     &
             xznr(config%q_nqx),  dvx(config%q_nqx),    &
             yznl(config%q_nqy), yzn(config%q_nqy),     &
             yznr(config%q_nqy),  dvy(config%q_nqy),    &
             zznl(config%q_nqz), zzn(config%q_nqz),     &
             zznr(config%q_nqz),  dvz(config%q_nqz), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module mesh_hy 1 failed")
    end if

    mem_local = config%q_nqx * 8._rk * 4._rk + config%q_nqy * 8._rk * 4._rk + config%q_nqz * 8._rk * 4._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "mesh_hy")

  end subroutine allocate_mesh_hy
!>
!> \par This subroutine deallocates the arrays from module mesh_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_mesh_hy
    
    use precision
    use abort

    implicit none

      
    integer(kind=ik) :: istat

    deallocate(xznl, xzn,  xznr,  dvx, yznl, yzn,  yznr,  dvy, &
               zznl,  zzn,  zznr,  dvz, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module mesh_hy 1 failed")
    end if

  end subroutine deallocate_mesh_hy

end module mesh_hy

! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro related floats
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module gfloat_hy
  use precision
  
  implicit none
! LOCAL variables that are not in modules
  save

  real(kind=rk) :: gamma, time,    &
                        rob, alph0,   &
                       rhoin, uin, utin, uttin, pin, ein, gamein, gamcin,&
                       gravin, tin, srfint, vlfrac, &
                        trst, tout1
!, fpg, g, pmass trstrt tout cfl dtini dtmin dtmax cvisc small smlrho
! tma
! smallp, smalle smallu smallx rib gridlx gridly gridlz, dt
  real(kind=rk), allocatable :: xnin(:)


contains

!>
!> \par This subroutine allocates the arrays from module gfloat_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_gfloat_hy(mem_global)
    
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none

      
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(xnin(config%qn), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module mesh_hy 1 failed")
    end if

    mem_local = config%qn * 8._rk 

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "gfloat_hy")

  end subroutine allocate_gfloat_hy
!>
!> \par This subroutine allocates the arrays from module gfloat_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_gfloat_hy
    
    use precision
    use abort
    
    implicit none

      
    integer(kind=ik) :: istat

    deallocate(xnin, stat=istat)

    if (istat .ne. 0) then
       raise_abort("deallocate_hydro_arrays(): deallocation of module mesh_hy 1 failed")
    end if

  end subroutine deallocate_gfloat_hy
end module gfloat_hy



!>
!> \par This module provides some hydro related floats
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module lfloat_hy

  use precision
  implicit none
! LOCAL variables that are not in modules
  save
  real(kind=rk) :: xbot, xtop, ybot, ytop, ylft, yrgt, zlft, zrgt

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate(xbot, xtop, ybot, ytop, ylft, yrgt, zlft, zrgt)
#endif
#endif


end module lfloat_hy

! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module grd_hy
  use precision
  
  implicit none
! LOCAL variables that are not in modules
  save

 real(kind=rk), allocatable, dimension(:) :: x, xl, xr, dx, dvol, dtdx, areal, &
                                             area, arear

#if !(defined(PROGRAM_remap))   
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate(x, xl  , xr  , dx  ,       &
!$omp               dvol, dtdx, areal, area, arear)
#endif
#endif

contains

!>
!> \par This subroutine allocates the arrays from module grid_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_grd_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif


    allocate(x(config%q), xl(config%q), xr(config%q),     &
             dx(config%q), dvol(config%q), dtdx(config%q),&
             areal(config%q), area(config%q), arear(config%q), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module grid_hy failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = config%q*8._rk*9._rk * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "grd_hy")

  end subroutine allocate_grd_hy

!>
!> \par This subroutine deallocates the arrays from module grid_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_grd_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(x, xl, xr, dx, dvol, dtdx, areal, area, arear, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module grid_hy failed")
    end if
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif

  end subroutine deallocate_grd_hy

end module grd_hy

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module hydro_hy

  use precision

  implicit none
      ! LOCAL variables that are not in modules

  save

  real(kind=rk), allocatable, dimension(:) :: rho, v, u, ut, utt, p, e, &
                                              ek, ei, rhonu, unu, utnu, &
                                              uttnu, enu, pnu, rhoav,   &
                                              uav, utav, uttav, pav,    &
                                              rhoflx, uflx, utflx,      &
                                              uttflx, eflx, we, w,      &
                                              ugrid, urell, ugrdl, rhol,& 
                                              rhor, ul, ur, utl, utr,   &
                                              uttl, uttr, pl, pr, vl,   &
                                              vr, gravl, gravr, gamcl,  &
                                              gamcr, cflno, c, ce,      &
                                              urel, game, gamc, gamel,  &
                                              gamer, gameav, tmp, tmpnu,&
                                              edt, avis, uttp, utbt,    &
                                              utlt, utrt

  real(kind=rk), allocatable, dimension(:,:) :: xn, xnnu, xnav, xnflx, xnl, xnr

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate(rho   , v    , u    , ut    , utt   , &
!$omp               p     , e    , ek   , ei    , rhonu , &
!$omp               unu   , utnu , uttnu, enu   , pnu   , &
!$omp               rhoav , uav  , utav , uttav , pav   , &
!$omp               rhoflx, uflx , utflx, uttflx, eflx  , &
!$omp               we    , w    , ugrid, urell , ugrdl , &
!$omp               rhol  , rhor , ul   , ur   ,          &
!$omp               utl   , utr  , uttl , uttr ,          &
!$omp               pl    , pr   , vl   , vr   ,          &
!$omp               gravl , gravr, gamcl, gamcr,          &
!$omp               cflno , c    , ce   , urel ,          &
!$omp               game  , gamc , gamel, gamer, gameav,  &
!$omp               tmp   , tmpnu, edt  , avis ,          &
!$omp               uttp  , utbt , utlt , utrt ,          &
!$omp               xn    , xnnu, xnav,                   &
!$omp               xnflx , xnl , xnr)
#endif
#endif

contains
!>
!> \par This subroutine allocates the arrays from module hydro_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_hydro_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif

    allocate(rho(config%q), v(config%q), u(config%q),        &
              ut(config%q), utt(config%q), p(config%q),      &
              e(config%q), ek(config%q), ei(config%q),       &
             rhonu (config%q), unu(config%q), utnu(config%q),&
             uttnu(config%q), enu(config%q), pnu(config%q),  &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_hy 1 failed")
    end if

    allocate(rhoav(config%q), uav(config%q), utav(config%q),   &
             uttav(config%q), pav(config%q), rhoflx(config%q), & 
             uflx(config%q), utflx(config%q), uttflx(config%q),&
             eflx(config%q), we(config%q), w(config%q),        &
             ugrid(config%q), urell(config%q),                 &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_hy 2 failed")
    endif

    allocate(ugrdl(config%q), rhol(config%q), rhor(config%q),  &
             ul(config%q), ur(config%q), utl(config%q),        &
             utr(config%q), uttl(config%q), uttr(config%q),    &
             pl(config%q), pr(config%q), vl(config%q),         &
             vr(config%q), gravl(config%q),                    &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_hy 3 failed")
    endif

    allocate(gravr(config%q), gamcl(config%q), gamcr(config%q),  &
             cflno (config%q), c(config%q), ce(config%q),        &
             urel(config%q), game(config%q), gamc(config%q),     &
             gamel(config%q), gamer(config%q), gameav(config%q), & 
             tmp(config%q),  stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_hy 4 failed")
    endif

    allocate(tmpnu(config%q), edt(config%q), avis(config%q),     &
             uttp(config%q), utbt(config%q), utlt(config%q),     &
              utrt(config%q), xn(config%q,config%qn),            &
             xnnu(config%q,config%qn), xnav(config%q,config%qn), &
             xnflx(config%q,config%qn), xnl(config%q,config%qn), &
             xnr (config%q,config%qn), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_hy 5 failed")
    endif

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
  thread_num=omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = config%q * 8._rk * 63._rk + config%q*config%qn*8._rk*6._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "hydro_hy")

  end subroutine allocate_hydro_hy

!>
!> \par This subroutine deallocates the arrays from module hydro_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_hydro_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))  
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(rho, v, u, ut, utt, p, e, ek, ei, &
             rhonu, unu, utnu, uttnu, enu, pnu,        &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_hy 1 failed")
    end if

    deallocate(rhoav, uav, utav, uttav, pav, rhoflx, uflx, &
             utflx, uttflx, eflx, we, w, ugrid, urell,   &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_hy 2 failed")
    endif

    deallocate(ugrdl, rhol, rhor, ul, ur, utl, utr,  &
             uttl, uttr, pl, pr, vl, vr, gravl,    &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_hy 3 failed")
    endif

    deallocate(gravr, gamcl, gamcr, cflno, c, ce, urel, &
             game, gamc, gamel, gamer, gameav, tmp,       &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_hy 4 failed")
    endif

    deallocate(tmpnu, edt, avis, uttp, utbt, utlt, utrt, &
             xn, xnnu, xnav, xnflx, xnl,      &
             xnr, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_hy 5 failed")
    endif
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_hydro_hy

end module hydro_hy

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module physcs_hy
  use precision

  implicit none
! LOCAL variables that are not in modules
  save
!  real(kind=rk) :: grav(q), dphidx(q), dedt(q), s(q), fict(q) 

  real(kind=rk), allocatable, dimension(:) :: grav, dphidx, dedt, s, fict

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate (grav, dphidx, dedt, s, fict   )
#endif
#endif

  contains

!>
!> \par This subroutine allocates the arrays from module physcs_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
    subroutine allocate_physcs_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none

    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif

    allocate(grav(config%q), dphidx(config%q), dedt(config%q),   &
             s(config%q), fict(config%q), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module physcs_hy 1 failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local =  config%q*8._rk*5._rk * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "physcs_hy")

  end subroutine allocate_physcs_hy
!>
!> \par This subroutine deallocates the arrays from module physcs_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
    subroutine deallocate_physcs_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(grav, dphidx, dedt, s, fict, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module physcs_hy 1 failed")
    end if
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP end parallel
#endif
#endif
    end subroutine deallocate_physcs_hy


end module physcs_hy

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!>
!> \par This module provides some Riemann-solver related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module reman_hy
  use precision

  implicit none
! LOCAL variables that are not in modules

  save

  real(kind=rk), allocatable, dimension(:) :: &
                       ppl   , pml   , p0l   ,          &
                      upl   , uml   , u0l   ,           &
                      utpl  , utml  , ut0l  , utt0l ,&
                      rhopl , rhoml , rho0l ,           &
                      clft  , plft  , uttlft,           &
                      ulft  , vlft  , utlft , wlft  ,&
                      crght , prght , vrght ,           &
                      urght , utrght, uttrgt, wrght ,&
                      pstar , ustar , vstar , cestar,&
                      rhostr, westar, ps    , us    , &
                      uts   , utts  , vs    , rhos  ,&
                      ces   , ws    , wes   , gmstar,&
                      game0l, gmelft, gmergt, gamc0l,&
                      gmclft, gmcrgt, games , gamcs ,&
                      gmemin, gmemax

  real(kind=rk), allocatable, dimension(:,:) :: &
                      xnpl, xnml, xn0l, xnlft, xnrght, xnstar,   &
                      xns
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate (ppl   , pml   , p0l   ,         &
!$omp                upl   , uml   , u0l   ,         &
!$omp                utpl  , utml  , ut0l  , utt0l , &
!$omp                rhopl , rhoml , rho0l ,         &
!$omp                clft  , plft  , uttlft,         &
!$omp                ulft  , vlft  , utlft , wlft  , &
!$omp                crght , prght , vrght ,         &
!$omp                urght , utrght, uttrgt, wrght , &
!$omp                pstar , ustar , vstar , cestar, &
!$omp                rhostr, westar, ps    , us    , &
!$omp                uts   , utts  , vs    , rhos  , &
!$omp                ces   , ws    , wes   , gmstar, &
!$omp                game0l, gmelft, gmergt, gamc0l, &
!$omp                gmclft, gmcrgt, games , gamcs , &
!$omp                gmemin, gmemax,                 &
!$omp                xnpl  , xnml , xn0l  ,          &
!$omp                xnlft, xnrght, xnstar,          &
!$omp                xns)
#endif
#endif

contains

!>
!> \par This subroutine allocates the arrays from module reman_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  

  subroutine allocate_reman_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none

    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    allocate(ppl(config%q), pml(config%q), p0l(config%q),       &
             upl(config%q), uml(config%q), u0l(config%q),       &
             utpl(config%q), utml(config%q), ut0l(config%q),    &
             utt0l(config%q), rhopl(config%q), rhoml(config%q), & 
             rho0l(config%q), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 1 failed")
    end if


    allocate(clft(config%q), plft(config%q), uttlft(config%q),  &
             ulft(config%q), vlft(config%q), utlft(config%q),   &
             wlft  (config%q), crght(config%q), prght(config%q),&
             vrght(config%q), urght(config%q), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 2 failed")
    end if

    allocate(utrght(config%q), uttrgt(config%q), wrght(config%q),   &
             pstar(config%q), ustar(config%q), vstar(config%q),     &
             cestar(config%q), rhostr(config%q), westar(config%q),  &
             ps(config%q), us(config%q), uts(config%q), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 3 failed")
    end if

    allocate(utts(config%q), vs(config%q), rhos(config%q),        &
             ces(config%q), ws(config%q), wes(config%q),          &
             gmstar(config%q), game0l(config%q), gmelft(config%q),& 
             gmergt(config%q), gamc0l(config%q), gmclft(config%q),&
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 4 failed")
    end if

    allocate(gmcrgt(config%q), games (config%q), gamcs (config%q), &
             gmemin(config%q), gmemax(config%q),                   &
             xnpl (config%q,config%qn), xnml  (config%q,config%qn),&
             xn0l  (config%q,config%qn), xnlft(config%q,config%qn),&
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 5 failed")
    end if

    allocate(xnrght(config%q,config%qn), xnstar(config%q,config%qn), &
             xns  (config%q,config%qn), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module reman_hy 6 failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = config%q*8._rk*43._rk + config%q*config%qn*8._rk*7._rk 
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "reman_hy")

  end subroutine allocate_reman_hy

!>
!> \par This subroutine deallocates the arrays from module reman_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_reman_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(ppl, pml, p0l, upl, uml, u0l, utpl, &
             utml, ut0l, utt0l, rhopl, rhoml, rho0l,&
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 1 failed")
    end if

    deallocate(clft, plft, uttlft, ulft, vlft, utlft, &
             wlft  , crght, prght, vrght, urght,       &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 2 failed")
    end if

    deallocate(utrght, uttrgt, wrght, pstar, ustar, vstar, &
             cestar, rhostr, westar, ps, us, uts,        &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 3 failed")
    end if

    deallocate(utts, vs, rhos, ces, ws, wes, gmstar, &
             game0l, gmelft, gmergt, gamc0l, gmclft,     &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 4 failed")
    end if

    deallocate(gmcrgt, games , gamcs , gmemin, gmemax,    &
             xnpl , xnml  , xn0l , xnlft,     &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 5 failed")
    end if

    deallocate(xnrght, xnstar, xns , stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module reman_hy 6 failed")
    end if
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_reman_hy


end module reman_hy

! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module spez_hy
  use precision
  
  
  implicit none
  SAVE
! LOCAL variables that are not in modules

  real(kind=rk), allocatable, dimension(:,:,:) :: &
           vxvold, vyvold, vzvold, denold

contains 
!>
!> \par This subroutine allocates the arrays from module spez_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine allocate_spez_hy(mem_global)

    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local


    allocate(vxvold(config%qx,qy_s:qy_e,qz_s:qz_e),    &
             vyvold(config%qx,qy_s:qy_e,qz_s:qz_e),    &
             vzvold(config%qx,qy_s:qy_e,qz_s:qz_e),    &
             denold(config%qx,qy_s:qy_e,qz_s:qz_e), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module spez_hy 1 failed")
    end if

    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*4._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "spez_hy")

  end subroutine allocate_spez_hy

!>
!> \par This subroutine deallocates the arrays from module spez_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_spez_hy

    use precision
    use abort

    use mo_mpi
    implicit none
    integer(kind=ik) :: istat

    deallocate(vxvold, vyvold, vzvold, denold, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module spez_hy 1 failed")
    end if

  end subroutine deallocate_spez_hy

end module spez_hy

! ---------------------------------------------------------------------
!>
!> \par This module provides some hydro-grid related variables
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module intrp_hy
  use precision

  implicit none
! LOCAL variables that are not in modules
  save
!  real(kind=rk) :: coeff1(q), coeff2(q), coeff3(q), coeff4(q), &
!                       coeff5(q), dela  (q), flatn (q), flatn1(q), &
!                       de    (q), shockd(q),                       &
!                       dp    (q), du    (q), dut   (q), dutt  (q), &
!                       drho  (q), dgame (q), dgamc (q), dgrav (q), &
!                       p6    (q), u6    (q), ut6   (q), utt6  (q), &
!                       rho6  (q), game6 (q), gamc6 (q), grav6 (q), &
!                       dxn(q,qn), xn6(q,qn) 

  real(kind=rk), allocatable, dimension(:) :: & 
                 coeff1, coeff2, coeff3, coeff4, &
                 coeff5, dela  , flatn , flatn1, &
                 de    , shockd,                       &
                 dp    , du    , dut   , dutt  , &
                 drho  , dgame , dgamc , dgrav , &
                 p6    , u6    , ut6   , utt6  , &
                 rho6  , game6 , gamc6 , grav6 

  real(kind=rk), allocatable, dimension(:,:) ::dxn, xn6 

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate (coeff1, coeff2, coeff3, coeff4, &
!$omp                coeff5, dela  , flatn , flatn1, &
!$omp                de    , shockd,                 &
!$omp                dp    , du    , dut   , dutt  , &
!$omp                drho  , dgame , dgamc , dgrav , &
!$omp                p6    , u6    , ut6   , utt6  , &
!$omp                rho6  , game6 , gamc6 , grav6 , &
!$omp                dxn, xn6)
#endif
#endif

contains
  
!>
!> \par This subroutine allocates the arrays from module intrp_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>

  subroutine allocate_intrp_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local


#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    allocate(coeff1(config%q), coeff2(config%q),         &
             coeff3(config%q), coeff4(config%q),         &
             coeff5(config%q), dela  (config%q),         &
             flatn (config%q), flatn1(config%q),         &
             de    (config%q), shockd(config%q),         &
             dp    (config%q), du    (config%q),         &
             dut   (config%q), dutt  (config%q),         &
             drho  (config%q), dgame (config%q),         &
             dgamc (config%q), dgrav (config%q),         &
             p6    (config%q), u6    (config%q),         &
             ut6   (config%q), utt6  (config%q),         &
             rho6  (config%q), game6 (config%q),         &
             gamc6 (config%q), grav6 (config%q),         &
             dxn(config%q,config%qn),                    &
             xn6(config%q,config%qn), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module intrp_hy 1 failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = config%q*8._rk*26._rk + config%q*config%qn*8._rk*2._rk 
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "intrp_hy")

  end subroutine allocate_intrp_hy
!>
!> \par This subroutine deallocates the arrays from module intrp_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  

  subroutine deallocate_intrp_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(coeff1, coeff2, coeff3, coeff4, &
             coeff5, dela  , flatn , flatn1, &
             de    , shockd,                       &
             dp    , du    , dut   , dutt  , &
             drho  , dgame , dgamc , dgrav , &
             p6    , u6    , ut6   , utt6  , &
             rho6  , game6 , gamc6 , grav6 , &
             dxn, xn6, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module intrp_hy 1 failed")
    end if
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif

  end subroutine deallocate_intrp_hy
end module intrp_hy



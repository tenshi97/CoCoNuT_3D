!>
!> \par This module provides the arrays that are needed for the viscosity in the transport code
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
module arvisco_rt
  use precision
  
  implicit none
  save
!c LOCAL variables that are not in modules
  real(kind=rk) ,parameter :: avisc_0=1.0_rk, tau_0=0.5_rk

  real(kind=rk), allocatable :: avisc(:,:,:), avis_j(:,:,:,:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(avisc,avis_j)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif

contains

!>
!> \par This subroutine allocates the arrays from module arvisco_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine allocate_arvisco_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure

    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    allocate(avisc(-1:config%imaxp +1,config%iemax,config%isma),   &
            avis_j(0:config%imaxp +1,config%iemax,config%isma,2), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module arvisco_rt failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+3)*config%iemax*config%isma*8._rk + (config%imaxp +1+2)*config%iemax*config%isma*2*8._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "arvisco_rt")

  end subroutine allocate_arvisco_rt
!>
!> \par This subroutine deallocates the arrays from module arvisco_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_arvisco_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(avisc, avis_j, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module arvisco_rt failed")
    end if
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif

  end subroutine deallocate_arvisco_rt
end module arvisco_rt

!c ------------------------------------------------------------

!>
!> \par This module provides the arrays that are needed for the neutrino matter interactions
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>
module matkoeff_rt
  use precision

  use configure
  implicit none
  save
!c LOCAL varibales that are not in modules

  real(kind=rk), allocatable :: bplaq (:,:,:,:), &
                                bpla  (:,:,:,:), &
                                aabsq (:,:,:,:), &  ! index 1 = sum of all rates
                                                    !       2 = LMS
                                                    !       3 = bruenn rates
                                                    !       4 = FFN rates
                                                    !       5 = KJT rates
                                aabs  (:,:,:), &
                                sikq  (:,:,:), &
                                sik   (:,:,:), &
                                sukq  (:,:,:), &
                                suk   (:,:,:), &
                                babsq (:,:,:)

! derivatives

  real(kind=rk), allocatable :: dbqdy(:,:,:,:), &
                                dbqde(:,:,:,:), &
                                dkqdy(:,:,:,:), &
                                dkqde(:,:,:,:), &
                                dktqdy(:,:,:),  &
                                dktqde(:,:,:),  &
                                dkthlp(:,:,:)

! for neutrino-electron/positron-scattering and alike

  real(kind=rk), allocatable :: sknesq(:,:,:,:,:), &
                                sknes (:,:,:,:,:), &
                                skeprq(:,:,:,:,:), &
                                skepr (:,:,:,:,:) 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(bplaq,bpla,aabsq,aabs,sikq,sik,sukq, &
!$omp               suk,babsq,dbqdy,dbqde,dkqdy,dkqde,   &
!$omp               dktqdy,dktqde,dkthlp, &
!$omp               sknesq,sknes,skeprq,skepr)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif
contains
!>
!> \par This subroutine allocates the arrays from module matkoeff_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine allocate_matkoeff_rt(mem_global)
   use precision
   use abort
   use print_stdout_mod, only : print_memory_alloc
   
   use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif

    allocate(bplaq( 0:config%imaxp +1,config%iemax,config%isma,2),        &
              bpla(-1:config%imaxp +1,config%iemax,config%isma,2),        &
             aabsq( 0:config%imaxp +1,config%iemax,config%isma,5),        &
            aabs  (-1:config%imaxp +1,config%iemax,config%isma),          &
             sikq  ( 0:config%imaxp +1,config%iemax,config%isma),         &
             sik   (-1:config%imaxp +1,config%iemax,config%isma),         &
             sukq  ( 0:config%imaxp +1,config%iemax,config%isma),         &
             suk   (-1:config%imaxp +1,config%iemax,config%isma),         &
             babsq ( 0:config%imaxp +1,config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkoeff_rt failed")
    end if

    allocate(dbqdy(0:config%imaxp +1,config%iemax,config%isma,2),         &
             dbqde(0:config%imaxp +1,config%iemax,config%isma,2),         &
             dkqdy(0:config%imaxp +1,config%iemax,config%isma,5),         &
             dkqde(0:config%imaxp +1,config%iemax,config%isma,5),         &
             dktqdy(0:config%imaxp +1,config%iemax,config%isma),          &
             dktqde(0:config%imaxp +1,config%iemax,config%isma),          &
             dkthlp(0:config%imaxp +1,config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkoeff_rt 2 failed")
    end if

    allocate(sknesq(4, 0:config%imaxp +1,config%iemax,config%isma,2), &
             sknes (4,-1:config%imaxp +1,config%iemax,config%isma,2), &
             skeprq(6, 0:config%imaxp +1,config%iemax,config%isma,2), &
             skepr (6,-1:config%imaxp +1,config%iemax,config%isma,2), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkoeff_rt 3 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+1)*config%iemax*config%isma*2*8._rk *3._rk +  &
                (config%imaxp +1+3)*config%iemax*config%isma*2*8._rk  +        &
                (config%imaxp +1+1)*config%iemax*config%isma*5*8._rk *3._rk +  &
                (config%imaxp +1+3)*config%iemax*config%isma*8._rk*3._rk     + &
                (config%imaxp +1+1)*config%iemax*config%isma*8._rk*6._rk    +  &
                4*(config%imaxp +1+1)*config%iemax*config%isma*2*8._rk      +  &
                4*(config%imaxp +1+3)*config%iemax*config%isma*2*8._rk      +  &
                6*(config%imaxp +1+1)*config%iemax*config%isma*2*8._rk      +  &
                6*(config%imaxp +1+3)*config%iemax*config%isma*2*8._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local
    
    call print_memory_alloc(mem_local, mem_global, "matkoeff_rt")

  end subroutine allocate_matkoeff_rt

!>
!> \par This subroutine allocates the arrays from module matkoeff_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_matkoeff_rt
   use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(bplaq, bpla, aabsq, aabs, sikq, sik, sukq, suk, babsq, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module matkoeff_rt failed")
    end if

    deallocate(dbqdy, dbqde, dkqdy, dkqde, dktqdy, dktqde, dkthlp, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module matkoeff_rt 2 failed")
    end if

    deallocate(sknesq, sknes, skeprq, skepr, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module matkoeff_rt 3 failed")
    end if
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_matkoeff_rt


end module matkoeff_rt

!c ------------------------------------------------------------
module neshelp_rt
  use precision

  implicit none 
!c LOCAL variables that are not in modules

  SAVE
  real(kind=rk), allocatable :: al10(:,:,:), ah10(:,:,:), &
                                al20(:,:,:), ah20(:,:,:) 

  real(kind=rk), allocatable :: al11(:,:,:), ah11(:,:,:), &
                                al21(:,:,:), ah21(:,:,:)


! weights for the Gauss-Legendre and gauss-Laguerre integration 
  integer(kind=ik) , parameter :: nla=10,nle=30


  real(kind=rk), allocatable ::  xle(:), wle(:), xla(:), wla(:)

contains

!>
!> \par This subroutine allocates the arrays from module neshelp_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine allocate_neshelp_rt(mem_global)
   use precision
   use abort
   use print_stdout_mod, only : print_memory_alloc
   use configure
   implicit none
   integer(kind=ik) :: istat
   real(kind=rk)    :: mem_global, mem_local

    allocate(al10(0:7,config%iemax,config%iemax),      &
             ah10(0:2,config%iemax,config%iemax),      &
             al20(0:7,config%iemax,config%iemax),      &
             ah20(0:2,config%iemax,config%iemax), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module neshelp_rt failed")
    end if    

    mem_local = 8*config%iemax*config%iemax*8._rk*2._rk + 3*config%iemax*config%iemax*8._rk*2._rk

    allocate(al11(0:7,config%iemax,config%iemax),      &
             ah11(0:2,config%iemax,config%iemax),      &
             al21(0:7,config%iemax,config%iemax),      &
             ah21(0:2,config%iemax,config%iemax), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module neshelp_rt 2 failed")
    end if

    mem_local = mem_local + 8*config%iemax*config%iemax*8._rk*2._rk + 3*config%iemax*config%iemax*8._rk*2._rk

    allocate(xle(nle),wle(nle),xla(nla),wla(nla), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module neshelp_rt 3 failed")
    end if

    mem_local = mem_local + nle*8._rk*2._rk + nla*8._rk * 2._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "neshelp_rt")

  end subroutine allocate_neshelp_rt

!>
!> \par This subroutine deallocates the arrays from module neshelp_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_neshelp_rt
   use precision
   use abort
   implicit none
   integer(kind=ik) :: istat

    deallocate(al10, ah10, al20, ah20, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module neshelp_rt failed")
    end if    

    deallocate(al11, ah11, al21, ah21, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module neshelp_rt 2 failed")
    end if

    deallocate(xle, wle, xla, wla, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module neshelp_rt 3 failed")
    end if


  end subroutine deallocate_neshelp_rt

end module neshelp_rt

!-----------------------------------------------------------------------
module neskernel_rt
  use precision
  
  implicit none
  save
! LOCAL variabls that are not in modules

  real(kind=rk), allocatable :: p0scat(:,:,:,:), p1scat(:,:,:,:) 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(p0scat,p1scat)
#endif  /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif
contains
!>
!> \par This subroutine allocates the arrays from module neskernel_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine allocate_neskernel_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif

    allocate(p0scat(0:config%imaxp +1,config%iemax,config%iemax,config%isma), &
             p1scat(0:config%imaxp +1,config%iemax,config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module neskernel_rt failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+1)*config%iemax*config%iemax*config%isma*8._rk*2._rk * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "nes_kernel_rt")

  end subroutine allocate_neskernel_rt
!>
!> \par This subroutine deallocates the arrays from module neskernel_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_neskernel_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(p0scat, p1scat, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module neskernel_rt failed")
    end if
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_neskernel_rt
end module neskernel_rt

!-----------------------------------------------------------------------
module pairkernel_rt
  use precision

  implicit none
  save
! LOCAL variables that are not in modules

  real(kind=rk), allocatable :: p0pabr(:,:,:,:), &
                                p1pabr(:,:,:,:), &
                                p2pabr(:,:,:,:) 
  real(kind=rk), allocatable :: zeta(:,:,:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(p0pabr,p1pabr,p2pabr,zeta)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif
contains
!>
!> \par This subroutine allocates the arrays from module pairkernel_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine allocate_pairkernel_rt(mem_global)
    use precision
    use print_stdout_mod, only : print_memory_alloc
    use abort
    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif

    allocate( p0pabr(0:config%imaxp +1,config%iemax,config%iemax,config%isma), &
              p1pabr(0:config%imaxp +1,config%iemax,config%iemax,config%isma), &
              p2pabr(0:config%imaxp +1,config%iemax,config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module pairernel_rt failed")
    end if

    allocate( zeta(0:config%imaxp +1,config%iemax,config%iemax), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module pairernel_rt 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+1)*config%iemax*config%iemax*config%isma*8._rk * 3._rk + &
                (config%imaxp +1+1)*config%iemax*config%iemax*8._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "pairkernel_rt")

  end subroutine allocate_pairkernel_rt

!>
!> \par This subroutine deallocates the arrays from module pairkernel_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine deallocate_pairkernel_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate( p0pabr, p1pabr, p2pabr, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module pairernel_rt failed")
    end if

    deallocate( zeta, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module pairernel_rt 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_pairkernel_rt
end module pairkernel_rt

!c-----------------------------------------------------------------------

module matkhelp_rt
  use precision

  implicit none
  save
! LOCAL varibales that are not in modules

  real(kind=rk), allocatable :: adxidv(:,:,:)  , &
                                asumdv(:,:,:)  , &
                                aabdv (:,:,:,:)

  real(kind=rk), allocatable :: redj(:,:,:), redh(:,:,:)
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(adxidv,asumdv,aabdv,redj,redh)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif
contains
!>
!> \par This subroutine allocates the arrays from module matkhelp_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine allocate_matkhelp_rt(mem_global)
    use precision
    use print_stdout_mod, only : print_memory_alloc
    use abort
    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    allocate(adxidv(0:config%imaxp +1,config%iemax,config%isma)  , &
             asumdv(0:config%imaxp +1,config%iemax,config%isma)  , &
             aabdv (0:config%imaxp +1,config%iemax,config%isma,2),  stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkhelp_rt 1 failed")
    end if 

    allocate(redj(0:config%imaxp +1,0:config%iemax+1,config%isma), &
             redh(0:config%imaxp +1,0:config%iemax+1,config%isma) , stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkhelp_rt 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+1)*config%iemax*config%isma*8._rk * 4._rk + &
                (config%imaxp +1+1)*(config%iemax+2)*config%isma*8._rk*2._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "matkhelp_rt")

  end subroutine allocate_matkhelp_rt
!>
!> \par This subroutine deallocates the arrays from module matkhelp_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_matkhelp_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(adxidv, asumdv ,  aabdv ,  stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module matkhelp_rt 2 failed")
    end if

    deallocate(redj, redh , stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module matkhelp_rt 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_matkhelp_rt
end module matkhelp_rt

!-----------------------------------------------------------------------

module boundary_rt
  use precision

! rip_ib is used to specify Intensity at inner boundary
!      on the "core-rays"

  implicit none
  SAVE

!!   use bo(config%cmin:config%imaxp +1,iemax,isma) for more general OBC
  real(kind=rk), allocatable :: rip_ib(:,:,:)
  real(kind=rk) :: bo

contains

!>
!> \par This subroutine allocates the arrays from module boundary_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine allocate_boundary_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(rip_ib(config%cmin:-1,config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module boundary_rf failed")
    end if 
    mem_local = (-1 - config%cmin +1)*config%iemax*config%isma*8._rk
    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "boundary_rt")

  end subroutine allocate_boundary_rt

!>
!> \par This subroutine allocates the arrays from module boundary_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine deallocate_boundary_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(rip_ib, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module boundary_rf failed")
    end if 
  end subroutine deallocate_boundary_rt
end module boundary_rt
!-----------------------------------------------------------------------

module backquants_rt
  use precision

  implicit none
!c LOCAL variabls that are not in modules

  save

  real(kind=rk), allocatable :: be(:)  , beq(:),    &
                                bey(:) , beyq(:),   &
                                beqs(:) , besq(:),  &
                                sqbe(:) , sbe(:),   &
                                dbedr(:), dbedrq(:),&
                                bedr(:), bedrq(:),  &
                                bea(:)  , beaq(:)

  real(kind=rk), allocatable :: ex(:)    , exq(:), &
                                alp(:), alpq(:),    &
                                gam(:)   , gamq(:), &
                                acc(:)   , acq(:),  &
                                dtlga(:) , dtlgaq(:), &
                                sqdex(:) , sdexq(:),  &
                                dxi1ex(:), dxi1exq(:)

  real(kind=rk), allocatable :: beythq(:)    , beyqth(:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(be,beq,bey,beyq,beqs,besq,sqbe,sbe, &
!$omp     dbedr,dbedrq,bedr,bedrq,bea,beaq,ex,exq,alp,alpq,gam,gamq, &
!$omp     acc,acq,dtlga,dtlgaq,sqdex,sdexq,dxi1ex,dxi1exq, &
!$omp     beythq,beyqth)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif

 real(kind=rk), allocatable ::  wl     (:), wlq     (:),  &
                                detg   (:), detgq   (:),  &
                                phigr  (:), phigrq  (:)

#ifdef CFC_TRANSPORT
 real(kind=rk), allocatable ::  phi2   (:), phi2q   (:),  &
                                wlalt  (:), wlaltq  (:),  &
                                detgalt(:), detgaltq(:),  &
                                wdetg  (:), wdetgq  (:),  &
                                w_phi2 (:), w_phi2q (:),  &
                                a_phi2 (:), a_phi2q (:),  &
                                bsh    (:), bshq    (:),  &
                                bshth  (:), bshthq  (:),  &
                                veff   (:), veffq   (:),  &
                                dwdt   (:), dwdtq   (:),  &
                                bealt  (:), beqalt  (:) 

 real(kind=rk), allocatable ::  daphi2dr (:), daphi2drq(:), &
                                dlgphidr (:), dlgphidrq(:), &
                                dvdtau   (:), dvdtauq  (:), &
                                bwco     (:), bwcoq    (:), &
                                bsco     (:), bscoq    (:), &
                                dwaphi2co(:), dwaphi2coq(:),&
                                waphi2co (:), waphi2coq(:), &
                                vdaphi2  (:), vdaphi2q (:), &
                                dvdr     (:), dvdrq    (:), &
                                dbphi2   (:), dbphi2q  (:), &
                                dlgphidt (:), dlgphidtq(:), &
                                dlgawdr  (:), dlgawdrq (:), &
                                veffdr   (:), veffdrq  (:), &
                                latrs    (:), latrsq   (:) 
#endif /* CFC_TRANSPORT */

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(wl, wlq, detg, detgq, phigr,    &
#ifdef CFC_TRANSPORT
!$omp     phi2, phi2q, wlalt, wlaltq, detgalt, detgaltq,    &
!$omp     wdetg, wdetgq, w_phi2, w_phi2q, a_phi2, a_phi2q,  &
!$omp     bsh, bshq, bshth, bshthq, veff, veffq, dwdt,      &
!$omp     dwdtq, bealt, beqalt, daphi2dr, daphi2drq,        &
!$omp     dlgphidr, dlgphidrq, dvdtau, dvdtauq, bwco,       &
!$omp     bwcoq, bsco, bscoq, dwaphi2co, dwaphi2coq,        &
!$omp     waphi2co, waphi2coq, vdaphi2, vdaphi2q, dvdr,     &
!$omp     dvdrq, dbphi2, dbphi2q, dlgphidt, dlgphidtq,      &
!$omp     dlgawdr, dlgawdrq, veffdr, veffdrq, latrs, latrsq, &
#endif /* CFC_TRANSPORT */
!$omp     phigrq)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif /* !PROGRAM_remap */

#ifdef CFC_TRANSPORT
 real(kind=rk), allocatable :: wllag   (:,:,:), wlqlag(:,:,:),    &
                               wllags (:,:,:), wlqlags  (:,:,:),  &
                               detglag (:,:,:), detgqlag (:,:,:), &
                               detglags(:,:,:), detgqlags(:,:,:), &
                               philag  (:,:,:), phiqlag  (:,:,:), &
                               beold   (:,:,:), beqold   (:,:,:), &
                               philags (:,:,:), phiqlags (:,:,:), &
                               beolds  (:,:,:), beqolds  (:,:,:)
#endif /* CFC_TRANSPORT */

contains
!>
!> \par This subroutine allocates the arrays from module backquants_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine allocate_backquants_rt(mem_global)
    use precision
    use print_stdout_mod, only : print_memory_alloc
    use abort
     use configure

    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

    mem_local = 0._rk

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif /* !PROGRAM_remap */

    allocate(be(-1:config%imaxp +1)  , beq(0:config%imaxp +1), bey(-1:config%imaxp +1) , beyq(0:config%imaxp +1), &
             beqs(0:config%imaxp +1), besq(0:config%imaxp +1), sqbe(0:config%imaxp +1) , sbe(0:config%imaxp +1),  &
             dbedr(0:config%imaxp +1), dbedrq(-1:config%imaxp +1), bedr(-1:config%imaxp +1),             &
             bedrq(0:config%imaxp +1), bea(0:config%imaxp +1), beaq(0:config%imaxp +1), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rt 1 failed")
    end if 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical
    mem_local = mem_local + (config%imaxp +1+1)*8._rk*10._rk + (config%imaxp +1+2)*8._rk*4._rk
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

    allocate(ex(-1:config%imaxp +1), exq(0:config%imaxp +1), alp(-1:config%imaxp +1), alpq(0:config%imaxp +1),    &
             gam(-1:config%imaxp +1), gamq(0:config%imaxp +1), acc(-1:config%imaxp +1), acq(0:config%imaxp +1),   &
             dtlga(-1:config%imaxp +1) , dtlgaq(0:config%imaxp +1), sqdex(-1:config%imaxp +1) ,          &
             sdexq(0:config%imaxp +1), dxi1ex(-1:config%imaxp +1), dxi1exq(0:config%imaxp +1), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rf 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical
    mem_local = mem_local + (config%imaxp +1+1)*8._rk*7._rk + (config%imaxp +1+2)*8._rk*7._rk
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

    allocate(beythq(-1:config%imaxp +1)    , beyqth(0:config%imaxp +1), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rt 3 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical
    mem_local = mem_local + (config%imaxp +1+1)*8._rk + (config%imaxp +1+2)*8._rK
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

    allocate(wl(-1:config%imaxp +1), wlq(0:config%imaxp +1), detg(-1:config%imaxp +1),         &
             detgq(0:config%imaxp +1),  phigr(-1:config%imaxp +1), phigrq(0:config%imaxp +1),  &
             stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rt 3a failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical    
    mem_local = mem_local + (config%imaxp +1+1)*8._rk*3._rk + (config%imaxp +1+2)*8._rk*3._rk
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

#ifdef CFC_TRANSPORT
    allocate(phi2   (-1:config%imaxp +1), phi2q   (0:config%imaxp +1),  &
             wlalt  (-1:config%imaxp +1), wlaltq  (0:config%imaxp +1),  &
             detgalt(-1:config%imaxp +1), detgaltq(0:config%imaxp +1),  &
             wdetg  (-1:config%imaxp +1), wdetgq  (0:config%imaxp +1),  &
             w_phi2 (-1:config%imaxp +1), w_phi2q (0:config%imaxp +1),  &
             a_phi2 (-1:config%imaxp +1), a_phi2q (0:config%imaxp +1),  &
             bsh    (-1:config%imaxp +1), bshq    (0:config%imaxp +1),  &
             bshth  (-1:config%imaxp +1), bshthq  (0:config%imaxp +1),  &
             veff   (-1:config%imaxp +1), veffq   (0:config%imaxp +1),  &
             dwdt   (-1:config%imaxp +1), dwdtq   (0:config%imaxp +1),  &
             bealt  (-1:config%imaxp +1), beqalt  (0:config%imaxp +1), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rf 4 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical 
    mem_local = mem_local + (config%imaxp +1+1)*8._rk*11._rk + &
                (config%imaxp +1+2)*8._rk*11._rk
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

    allocate(daphi2dr (-1:config%imaxp +1), daphi2drq(0:config%imaxp +1), &
             dlgphidr (-1:config%imaxp +1), dlgphidrq(0:config%imaxp +1), &
             dvdtau   (-1:config%imaxp +1), dvdtauq  (0:config%imaxp +1), &
             bwco     (-1:config%imaxp +1), bwcoq    (0:config%imaxp +1), &
             bsco     (-1:config%imaxp +1), bscoq    (0:config%imaxp +1), &
             dwaphi2co(-1:config%imaxp +1), dwaphi2coq(0:config%imaxp +1),&
             waphi2co (-1:config%imaxp +1), waphi2coq(0:config%imaxp +1), &
             vdaphi2  (-1:config%imaxp +1), vdaphi2q (0:config%imaxp +1), &
             dvdr     (-1:config%imaxp +1), dvdrq    (0:config%imaxp +1), &
             dbphi2   (-1:config%imaxp +1), dbphi2q  (0:config%imaxp +1), &
             dlgphidt (-1:config%imaxp +1), dlgphidtq(0:config%imaxp +1), &
             dlgawdr  (-1:config%imaxp +1), dlgawdrq (0:config%imaxp +1), &
             veffdr   (-1:config%imaxp +1), veffdrq  (0:config%imaxp +1), &
             latrs    (-1:config%imaxp +1), latrsq   (0:config%imaxp +1) , stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rf 5 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp critical 
    mem_local = (config%imaxp +1+1)*8._rk * 14_rk + (config%imaxp +1+2)*8._rk*14._rk
!$omp end critical
#endif
#endif /* !PROGRAM_remap */

#endif /* CFC_TRANSPORT */

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif /* !PROGRAM_remap */
    mem_local = mem_local * thread_num

#ifdef CFC_TRANSPORT
    allocate(wllag   (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             wlqlag   (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             wllags ( -1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             wlqlags  (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             detglag (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             detgqlag (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             detglags(-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             detgqlags(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             philag  (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             phiqlag  (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             beold   (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             beqold   (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             philags (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             phiqlags (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             beolds  (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
             beqolds  (0:config%imaxp +1,config%nystrt:config%nymom,config%nztra), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module backquants_rf 6 failed")
    end if

    mem_local = mem_local + (config%imaxp +1+2)*(config%nymom-config%nystrt+1)*config%nztra*8._rk *8._rk + &
               (config%imaxp +1+1)*(config%nymom-config%nystrt+1)*config%nztra*8._rk * 8._rk
#endif /* CFC_TRANSPORT */

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "backquants_rt")

  end subroutine allocate_backquants_rt

!>
!> \par This subroutine deallocates the arrays from module backquants_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_backquants_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif

    deallocate(be, beq, bey, beyq, beqs, besq, sqbe, sbe, dbedr, dbedrq, &
               bedr, bedrq, bea, beaq, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rt 1 failed")
    end if 

    deallocate(ex, exq, alp, alpq, gam, gamq, acc, acq, dtlga, dtlgaq, &
               sqdex, sdexq, dxi1ex, dxi1exq, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rt 2 failed")
    end if

    deallocate(beythq , beyqth, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rf 3 failed")
    end if

    deallocate(wl, wlq, detg, detgq, phigr, phigrq, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rt 3a failed")
    end if

#ifdef CFC_TRANSPORT
    deallocate(phi2, phi2q, wlalt, wlaltq, detgalt, detgaltq, wdetg,       &
               wdetgq, w_phi2, w_phi2q, a_phi2, a_phi2q, bsh, bshq, bshth, &
               bshthq, veff, veffq, dwdt, dwdtq, bealt, beqalt, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rf 4 failed")
    end if

    deallocate(daphi2dr,dlgphidr, dvdtau, bwco, bsco, dwaphi2co, waphi2co, &
               vdaphi2, dvdr, dbphi2, dlgphidt, dlgawdr, veffdr, latrs ,   &
               stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rf 5 failed")
    end if

#endif /* CFC_TRANSPORT */

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif

#ifdef CFC_TRANSPORT
    deallocate(wllag, wlqlag, wllags, wlqlags, detglag, detgqlag,  &
              detglags, detgqlags, philag, phiqlag, beold, beqold,&
              philags, phiqlags , beolds  , beqolds  , stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module backquants_rf 6 failed")
    end if

#endif /* CFC_TRANSPORT */

  end subroutine deallocate_backquants_rt

end module backquants_rt

!-----------------------------------------------------------------------

module lotdqua_rt
  use precision
      
  implicit none
  save
! LOCAL varibales that are not in modules

! --
! local thermodynamical quantities of the TRANSP-grid
! defined on cell centers rq
! "overlap" config%imaxp +1+2 .... novl(j) contains the corresponding
! values of the hydro-grid
!--
!  real(kind=rk) :: rhq    (0:config%imaxp +1), tmq    (0:config%imaxp +1), &
!                   eiq    (0:config%imaxp +1), ceq    (0:config%imaxp +1), &
!                   cnq    (0:config%imaxp +1), cpq    (0:config%imaxp +1), &
!                   cuq    (0:config%imaxp +1), xnq    (0:config%imaxp +1,qn), &
!                   xhrepq (0:config%imaxp +1), zah    (0:config%imaxp +1,2)

  real(kind=rk), allocatable :: rhq(:), tmq(:), eiq(:), ceq(:), cnq(:), &
                                cpq(:), cuq(:), xnq(:,:), xhrepq(:),    &
                                zah(:,:)

  integer(kind=ik), allocatable :: ish_ra(:)

!  real (kind=rk) :: eiqold(0:config%imaxp +1),yeqold(0:config%imaxp +1), &
!                    eiqstar(0:config%imaxp +1),yeqstar(0:config%imaxp +1)

  real(kind=rk), allocatable :: eiqold(:),yeqold(:), &
                                eiqstar(:),yeqstar(:) 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(rhq,tmq,eiq,ceq,cnq,cpq,cuq, &
!$omp               xnq,xhrepq,zah,eiqold,yeqold,eiqstar,yeqstar)
!$omp threadprivate(ish_ra)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif
contains

!>
!> \par This subroutine allocates the arrays from module lotdqua_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine allocate_lotdqua_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    allocate(rhq    (0:config%imaxp +1), tmq    (0:config%imaxp +1), &
             eiq    (0:config%imaxp +1), ceq    (0:config%imaxp +1), &
             cnq    (0:config%imaxp +1), cpq    (0:config%imaxp +1), &
             cuq    (0:config%imaxp +1), xnq    (0:config%imaxp +1,config%qn), &
             xhrepq (0:config%imaxp +1), zah    (0:config%imaxp +1,2), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module lotdqua_rf 1 failed")
    end if 

    allocate(eiqold(0:config%imaxp +1),yeqold(0:config%imaxp +1), &
             eiqstar(0:config%imaxp +1),yeqstar(0:config%imaxp +1), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module lotdqua_rf 2 failed")
    end if

    if (config%jvisc .eq. 3) then
       allocate(ish_ra(0:config%imaxp +1), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lotdqua_rf 3 failed")
       end if
    endif

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    if (config%jvisc .eq. 3) then
       mem_local = (config%imaxp +1+1)*8._rk *12._rk +      &
                   (config%imaxp +1+1)*2*8._rk     +        &
                   (config%imaxp +1+1)*config%qn*8._rk    + &
                   (config%imaxp +1+1)*4._rk
    else
       mem_local = (config%imaxp +1+1)*8._rk *12._rk + &
                   (config%imaxp +1+1)*2*8._rk      +  &
                   (config%imaxp +1+1)*config%qn*8._rk
    endif

    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "lotdqua_rt")

  end subroutine allocate_lotdqua_rt
!>
!> \par This subroutine allocates the arrays from module lotdqua_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine deallocate_lotdqua_rt
    use precision
    use abort
    use configure
    implicit none
    integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(rhq, tmq, eiq, ceq, cnq, cpq, cuq, xnq, &
               xhrepq, zah, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module lotdqua_rf 1 failed")
    end if 

    deallocate(eiqold, yeqold, eiqstar, yeqstar, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module lotdqua_rf 2 failed")
    end if

    if (config%jvisc .eq. 3) then
       deallocate(ish_ra, stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): deallocation of module lotdqua_rf 3 failed")
       end if
    endif

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_lotdqua_rt
end module lotdqua_rt

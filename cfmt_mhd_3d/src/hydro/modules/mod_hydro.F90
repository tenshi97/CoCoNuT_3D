!>
!> \par This module provides the tot-arrays of the hydro-code with the CFC or the PROMETHEUS version
!>
!> \author  M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module totare_hy
  use precision
 
  implicit none
! LOCAL variables that are not in modules

  save
!-----------------------------------------------------------------------
!     quantities of the total calculation area:
!-----------------------------------------------------------------------

  
#ifndef PROGRAM_remap
!  integer(kind=ik) :: nutot = NQN
#else
!  integer(kind=ik) :: nutot = 0
#endif /* PROGRAM_remap */



  real(kind=rk), allocatable, dimension(:,:,:) :: &
                 vextot, veytot, veztot, vexold, veyold, vezold,  &
                 dnsold, acxtot, qgrtot, dentot, enetot, pretot,  &
                 gaetot, gactot, stotot, temtot, wltot, gpotot, gpoold 

  real(kind=rk), allocatable, dimension(:,:,:,:) :: xnutot

  integer(kind=ik), allocatable :: ishtot(:,:,:)

  real(kind=rk), allocatable :: tmatot(:), ugrtot(:)
  real(kind=rk), allocatable :: dmdtiotot(:,:,:)


  real(kind=rk),allocatable ::    xzlnew(:), xznnew(:), xzrnew(:),  &
                                  yzlnew(:), yznnew(:), yzrnew(:),  &
                                  zzlnew(:), zznnew(:), zzrnew(:),  &
                                  xzrold(:)


  real(kind=rk), allocatable, dimension(:) ::  &
                        xzltot, xzntot, xzrtot, dvxtot,  &
                        yzltot, yzntot, yzrtot, dvytot,  &
                        zzltot, zzntot, zzrtot, dvztot
contains 
!>
!> \par This subroutine allocates the arrays from module totare_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_totare_hy( mem_global)
    use precision
    use abort
    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc

    
    use configure

    implicit none

    integer(kind=ik)        :: istat
    real(kind=rk)           :: mem_global, mem_local


    allocate(xzrnew(0:config%q_nqx), xzrold(0:config%q_nqx), stat=istat)
    
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy failed")
    end if

    mem_local = (config%q_nqx+1)*8._rk*2._rk

#ifdef PROGRAM_remap

    ! with remaper: use still the value q = max(qx,qy,qz) + 20
    ! for definition of arrays, and just write later in the correct
    ! dimension
    allocate(xzlnew(0:config%q), xznnew(0:config%q), &
             yzlnew(0:config%q), yznnew(0:config%q), &
             yzrnew(0:config%q), zzlnew(0:config%q), &
             zznnew(0:config%q), zzrnew(0:config%q), &
             stat=istat)

    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy failed")
    end if
#endif /* PROGRAM_remap */

    allocate(xzltot(config%q_nqx), xzntot(config%q_nqx),  &
             xzrtot(config%q_nqx), dvxtot(config%q_nqx),  &
             yzltot(config%q_nqy), yzntot(config%q_nqy),  &
             yzrtot(config%q_nqy), dvytot(config%q_nqy),  &
             zzltot(config%q_nqz), zzntot(config%q_nqz),  &
             zzrtot(config%q_nqz), dvztot(config%q_nqz), stat=istat)

    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module 1a totare_hy failed")
    end if

    mem_local = mem_local + (config%q_nqx)*8._rk*4._rk + &
                             config%q_nqy*8._rk*4._rk  + &
                             config%q_nqz*8._rk*4._rk

    allocate(tmatot(0:config%qx), ugrtot(config%qx),     &
             dmdtiotot(6,config%qy,config%qz), stat=istat)
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy 2 failed")
    end if

     mem_local = mem_local + (config%qx+1)*8._rk +   &
                              config%qx * 8._rk  +   &
                              6*config%qy*config%qz*8._rk   

    allocate(ishtot(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),  &
             vextot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             veytot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             veztot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             dentot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             enetot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             pretot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             gaetot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             gactot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             stotot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             temtot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             vexold(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             veyold(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             vezold(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             dnsold(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             acxtot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             qgrtot(config%qx,qy_s:qy_e,qz_s:qz_e),              &
             xnutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%qn),    &
             gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e),          &
             gpoold(0:config%qx,qy_s-1:qy_e,qz_s:qz_e),          &
             wltot(config%qx,qy_s:qy_e,qz_s:qz_e), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module totare_hy 1 failed")
    end if

     mem_local = mem_local + (config%qx+2)*(qy_e-qy_s+3)*        &
                             (qz_e-qz_s+3)*4._rk +               &
                             config%qx*(qy_e-qy_s+1)*            &
                             (qz_e-qz_s+1)*8._rk * 17._rk +      &
                             config%qx*(qy_e-qy_s+1)*            &
                             (qz_e-qz_s+1)*config%qn*8._rk     + &
                             (config%qx+1)*(qy_e-qy_s+2)*        &
                             (qz_e-qz_s+1)*8._rk *2._rk 

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "totare_hy")

  end subroutine allocate_totare_hy
!>
!> \par This subroutine deallocates the arrays from module totare_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_totare_hy
    use precision
      use abort
      
    implicit none
    integer(kind=ik) :: istat

    deallocate(xzrnew, xzrold, stat=istat)
    
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy failed")
    end if

#ifdef PROGRAM_remap
    deallocate(xzlnew, xznnew, yzlnew, yznnew,  &
               yzrnew, zzlnew, zznnew, zzrnew,  &
               stat=istat)
     
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy failed")
    end if
#endif

    deallocate(xzltot, xzntot, xzrtot, dvxtot ,  &
               yzltot, yzntot, yzrtot, dvytot ,  &
               zzltot, zzntot, zzrtot, dvztot, stat=istat)


    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module 1a totare_hy failed")
    end if

    deallocate(tmatot, ugrtot, dmdtiotot, stat=istat)
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy 2 failed")
    end if


  
    deallocate(ishtot, vextot, veytot, veztot, dentot, &
               enetot, pretot, stotot, temtot, vexold, &
               veyold, vezold, dnsold, acxtot, qgrtot, &
               xnutot, gpotot, gpoold, stat=istat)
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module totare_hy 1 failed")
    end if


  end subroutine deallocate_totare_hy

end module totare_hy

!> \verbatim
!> this module provides necessary quantities of the neutrino transport:
!>
!>     isma     = number of neutrino species
!>     qc     = number of chemical potentials
!>
!>     qyetot = local source term for lepton number (g/cm^3/s)
!>     qentot = local source term for energy   (erg/cm^3/s)
!>     qmotot = local source term for momentum (g/s/cm^2/s)
!>
!>     sumqye = total lepton fraction (\rho Y_e) loss-rate (g/s)
!>     sumqen = total energy                     loss-rate (erg/s)
!>     sumqmo = total momentum                   loss-rate (g cm/s)
!>
!>     enutot = local neutrino energies (MeV)
!>     fnutot = local neutrino momentum (MeV)
!>     dnutot = local neutrino density  (1/cm^3)
!>     gnutot = local neutrino flux  (?)
!>
!>     cpotot(*,*,*,n): chemical potentials in MeV
!>     n = 1 -> neutrinos
!>     n = 2 -> electrons
!>     n = 3 -> neutrons
!>     n = 4 -> protons
!>
!>     parameters of the neutrino transport:
!>
!>     el0obs = initial neutrino luminosities (erg/s) at infinity
!>     eltobs = current neutrino luminosities (erg/s) at infinity
!>     etnobs = current neutrino energy (MeV) at infinity
!>     tmp_in = initial neutrino temperature (MeV)
!>     tmp_fi = final neutrino temperature (MeV)
!>     dgp_in = initial degeneracy parameter
!>     del_ye = average reduction of lepton fraction in the inner core
!>     del_en = total energy loss of the inner core (erg)
!>     del_tl = characteristic time scale (s) of the luminosity evolution
!>
!>     r_nusp = radius (cm) of the neutrino sphere
!>     i_nusp = index of the neutrino sphere
!>
!>
!>     sumdenu = (\int neutrino energy sourceterm dV) * dt  [erg/ccm]
!>     sumdynu = (\int neutrino lepton sourceterm dV) * dt  [?/ccm]
!>                         , summed over all timesteps between two outputs
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>

module nutrio_hy
  use precision
  
!#ifndef NOTRA

!#endif

  implicit none
! LOCAL varibales that are not in modules

  SAVE
!-----------------------------------------------------------------------
!     "nutrio" is used in: (PNS-Version)
!     qterms, ppm_nutra, pltout
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------


  integer(kind=ik) ::     i_nusp , i_nusc
  integer(kind=ik), parameter :: qc=4

  real(kind=rk), allocatable :: qyetot(:,:,:,:), qentot(:,:,:),   & 
                                qmotot(:,:,:),   qmytot(:,:,:),   &
                                enutot(:,:,:,:), fnutot(:,:,:,:), &
                                pnutot(:,:,:,:), dnutot(:,:,:,:), &
                                gnutot(:,:,:,:), cpotot(:,:,:,:)

  real(kind=rk), allocatable :: qyeinp(:,:,:), qeninp(:,:,:),     & 
                                qmoinp(:,:,:),                    &
                                enuinp(:,:,:,:), fnuinp(:,:,:,:), &
                                pnuinp(:,:,:,:), dnuinp(:,:,:,:), &
                                gnuinp(:,:,:,:)

  real(kind=rk), allocatable :: eltobh(:,:,:), etnobh(:,:,:),      &
                                rnusph(:,:)


#ifdef FCNC_CALC
  real(kind=rk), allocatable, dimension(:,:,:) :: qyetot_n, qentot_n, qmotot_n
#endif

  real(kind=rk) :: sumdenu, sumdynu
  real(kind=rk) :: sumqen, sumqmo, sumqye

#ifdef FCNC_CALC
  real(kind=rk) :: sumdenu_n, sumdynu_n
#endif

  real(kind=rk) ::  del_en, del_ye, del_tl, r_nusp    


  real(kind=rk),allocatable  :: el0obs(:), eltobs(:), etnobs(:), tmp_in(:), &
                                tmp_fi(:), dgp_in(:), eltobc(:)

contains

!>
!> \par This subroutine allocates the arrays from module nutrio_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine allocate_nutrio_hy(mem_global)
    use precision
    use abort

    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc

    use configure

    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
 
    allocate(el0obs(config%isma), eltobs(config%isma),                     &
             etnobs(config%isma), tmp_in(config%isma),                     &
             tmp_fi(config%isma), dgp_in(config%isma),                     &
             eltobc(config%isma), eltobh(config%isma,qy_s:qy_e,qz_s:qz_e), &
             etnobh(config%isma,qy_s:qy_e,qz_s:qz_e),                      &
             rnusph(qy_s:qy_e,qz_s:qz_e), stat=istat)

    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    endif

    mem_local = config%isma * 8._rk * 7._rk + config%isma*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*2._rk + &
               (qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk 
    
    allocate(qyetot(config%qx,qy_s:qy_e,qz_s:qz_e,5),           &
             qentot(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             qmotot(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             qmytot(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             enutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             fnutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), & 
             dnutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             gnutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             cpotot(config%qx,qy_s:qy_e,qz_s:qz_e,qc),          &
             stat=istat)
    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if 
    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*8._rk + &
                config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%isma*8._rk*4._rk          + &
                config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*qc*8._rk

#ifndef NOTRA
    allocate(qyeinp(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             qeninp(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             qmoinp(config%qx,qy_s:qy_e,qz_s:qz_e),             &
             enuinp(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             fnuinp(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), & 
             dnuinp(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             gnuinp(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), &
             stat=istat)
    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if 
    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*3._rk + &
                config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%isma*8._rk*4._rk

#endif
    
    if (use_mpi) then
       allocate(pnutot(config%qx,qy_s-1:qy_e+1,qz_s:qz_e,config%isma), stat=istat)
    else
       allocate(pnutot(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), stat=istat)
    endif
       
   if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
   end if

   mem_local = mem_local + config%qx*(qy_e-qy_s+3)*(qz_e-qz_s+1)*config%isma*8._rk
 
#ifndef NOTRA
   allocate(pnuinp(config%qx,qy_s:qy_e,qz_s:qz_e,config%isma), stat=istat)
   if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
   end if
   mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%isma*8._rk
#endif


#ifdef FCNC_CALC    
    allocate(qyetot_n(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             qentot_n(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             qmotot_n(config%qx,qy_s:qy_e,qz_s:qz_e), stat=istat)

    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if
   mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1) * 8._rk*3._rk
#endif /* FCNC_CALC */

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

     call print_memory_alloc(mem_local, mem_global, "nutrio_hy")

  end subroutine allocate_nutrio_hy
!>
!> \par This subroutine deallocates the arrays from module nutrio_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine deallocate_nutrio_hy
    use precision
      use abort


    implicit none
    integer(kind=ik) :: istat

    deallocate(el0obs, eltobs, etnobs, tmp_in, &
               tmp_fi, dgp_in, eltobc, eltobh, etnobh,      &
               rnusph, stat=istat)

    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): deallocation of module nutrio_hy 1  failed")
    end if

     
    deallocate(qyetot, qentot, qmotot, qmytot, enutot, fnutot, & 
               pnutot, dnutot, gnutot, cpotot, stat=istat)
      
    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if

#ifndef NOTRA
    deallocate(qyeinp, qeninp, qmoinp, enuinp, fnuinp, & 
               pnuinp, dnuinp, gnuinp, stat=istat)

    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if
#endif

#ifdef FCNC_CALC    
    deallocate(qyetot_n, qentot_n, qmotot_n, stat=istat)

    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module nutrio_hy 1  failed")
    end if

#endif /* FCNC_CALC */



  end subroutine deallocate_nutrio_hy


end module nutrio_hy

!> \verbatim
!> this module provides the general relativistic quantities of the 
!> total calculation area:
!>
!>     defined at zone interfaces:
!>
!>    tgmtot(0:qx)     =  total enclosed  grav. mass
!>     gamtot(0:qx)     =  radial metric component at the end of a hyst
!>     gamold(0:qx)     =  radial metric component at the beginning
!>     dgatot(0:qx)     =  d(log gamma)/dt
!>     ephtot(0:qx)     =  lapse function
!>
!>     defined at zone center:
!>
!>     grptot(qx)       =  grad\Phi at the end of a hydro time step
!>     grpold(qx)       =  grad\Phi at the beginning of a hydro time step
!>
!>     i_grv            = 0 gen. relativity off
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module totgrq_hy

  use precision
  
  implicit none
! LOCAL variables that are not in modules

  SAVE

!-----------------------------------------------------------------------
!     "totgrq" is used in:
!     ppm_grv,
!-----------------------------------------------------------------------


!  integer(kind=ik) :: i_grv
  
  real(kind=rk), allocatable ::  tgmtot(:) ,    gamtot(:),  gamold(:), &    
                                 ephtot(:) ,    grptot(:),    grpold(:)
 
  real(kind=rk), allocatable :: delegr(:,:,:)

  real(kind=rk) :: sumdegr
  real(kind=rk), allocatable :: qgrv(:,:,:)


contains 
!>
!> \par This subroutine allocates the arrays from module totgrq_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_totgrq_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(tgmtot(0:config%qx), gamtot(0:config%qx),     &
             gamold(0:config%qx), ephtot(0:config%qx) ,    &
             grptot(config%qx), grpold(config%qx), stat=istat)
    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): allocation of module totgrq_hy2 failed")
    end if
    
    mem_local = (config%qx+1) * 8._rk * 4._rk + config%qx*8._rk*2._rk
    allocate(qgrv(config%qx,qy_s:qy_e,qz_s:qz_e),stat =istat)

    if (istat .ne. 0) then
      raise_abort("allocate_hydro_arrays(): allocation of module totgrq_hy failed")
    end if

    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

     call print_memory_alloc(mem_local, mem_global, "totgrq_hy")

  end subroutine allocate_totgrq_hy

!>
!> \par This subroutine deallocates the arrays from module totgrq_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_totgrq_hy
    use precision

    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(tgmtot,    gamtot,  gamold, &    
               ephtot ,    grptot   ,    grpold, stat=istat)
    if (istat .ne. 0) then
          raise_abort("allocate_hydro_arrays(): deallocation of module totgrq_hy2 failed")
    end if

    deallocate(qgrv, stat=istat)

  if (istat .ne. 0) then
      raise_abort("allocate_hydro_arrays(): deallocation of module totgrq_hy failed")
  end if



end subroutine deallocate_totgrq_hy

end module totgrq_hy

!-----------------------------------------------------------------------


!> \par This module provides integers use in the hydro code
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
module intgrs_hy
!-----------------------------------------------------------------------
!     "intgrs" is used in:
!     grid, input, pltout, restrt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Attention: 1. PPM_V2 - version -> without variable "iter"!!!
!                2. use intgrs.cb only in context with dimen.pb
!                   because "xyzswp" must be declared as integer!!!
!-----------------------------------------------------------------------

  use precision
  implicit none
! LOCAL varibales that are not in modules

  SAVE
      

  integer(kind=ik) ::   nzn, nzn1, nzn2, nzn3, nzn4, nzn5, &
                          nzn6, nzn7, nzn8, xyzswp, nstep,   &
                          igeom,    &
                          igrav,   nout1, nrst,      &
                           nresum      

                          ! nendd irstrt nrstrt nout itstp nriemm
                          ! igodu nsdim nx ny nz igeomx igeomy igeomz
                          ! isym
#ifndef PROGRAM_remap
!> \todo Should be removed once coconut is fully allocatable
#ifdef CFC_TRANSPORT
!  integer(kind=ik), parameter :: nuc=NQN, nuc_mf=NQN-1
#endif
#else
!  integer(kind=ik)           :: nuc, nuc_mf
#endif /* PROGRAM_remap */


! new hydro integers; originally these variables were defined as real
! in module gfloat_hy. This, however, does not makse sence

!  integer(kind=ik) :: bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,  &
!                      , bndbot, bndtop, bndlft, bndrgt, bndmin

  integer(kind=ik) ::  bndbot, bndtop, bndlft, bndrgt, bndmax, bndmin
end module intgrs_hy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!> \verbatim
!>
!>     arrays to control the different calculation areas:
!>
!>     are_nu               = number of used areas, must be less
!>                            or equal to are_di
!>     are_di               = dimesion of dt_are and ix_are
!>     dt_are(*)            = time step of this area
!>     dt_cfl(*)            = cfl time step of this area
!>     ti_are(*)            = time at the end of a time step
!>     ix_are(*, 1 - 9)     = ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz
!>           (*,10)         = isd = sweep direction
!>           (*,11)         = ndt = number of time steps
!>           (*,12 -13)     = bndmnx, bndmxx
!>     ndtmax               = max. number of small time steps
!>                            during one large time step
!>     dtmaxx               = max. size of one large time step
!>     nhystp               = number of hydro time steps
!>
!>     dt_rad               = (number of small time steps)**(-1)*
!>                            (time step for Neutrino-Transport)
!>
!>  \author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module hydro_areas_mod
  use precision

  implicit none
! LOCAL varibales that are not in modules

  save
  private
  public :: hydro_areas, allocate_hydro_areas, deallocate_hydro_areas, &
            areas, init_areas

  type hydro_areas
     integer(kind=ik)              :: are_nu, nhystp

     integer(kind=ik)              :: are_di
!     real(kind=rk)                 :: dt_rad, dtmaxx
     real(kind=rk), allocatable    :: dt_are(:), dt_cfl(:), ti_are(:)
     integer(kind=ik), allocatable :: ix_are(:,:)

     integer(kind=ik)              :: nx
     integer(kind=ik)              :: ny
     integer(kind=ik)              :: nz
  end type hydro_areas


  type(hydro_areas)                :: areas

contains

  subroutine init_areas

    use precision
    use configure
    use abort
    use state
    implicit none

    integer(kind=ik) :: i, j
    integer(kind=ik) :: nxa, nxb, nxc


#ifndef PROGRAM_remap
    areas%are_nu = config%are_nu
#else
    ! in case of PROGRAM_remap are_nu is set later in writeback
    areas%are_nu = 1
#endif
    areas%nx     = config%nx
    areas%ny     = config%ny
    areas%nz     = config%nz

    areas%are_di = 3

  if (areas%are_nu .lt. 1  .or.  areas%are_nu .gt. areas%are_di) then
     print *,areas%are_nu,areas%are_di
     raise_abort("input(): are_nu")
  endif

#ifndef PROGRAM_remap
  do i = 1, areas%are_di
     areas%dt_are(i) = config%dtini
     areas%dt_cfl(i) = config%dtini
     do j = 1, 13
        areas%ix_are(i,j) = -100
     enddo
  enddo
#else
  ! the arrays dt_are and dt_cfl are allocated later in remap case
#endif

  nxa = config%nxa
  nxb = config%nxb
  nxc = config%nxc
  
  if(areas%are_nu .eq. 1) nxa = areas%nx   
  if(areas%are_nu .eq. 2) nxb = areas%nx   
  if(areas%are_nu .eq. 3) nxc = areas%nx

#ifndef PROGRAM_remap
  i = 1
  areas%ix_are(i, 1) = 1
  areas%ix_are(i, 2) = nxa
  areas%ix_are(i, 3) = 1
  areas%ix_are(i, 4) = 1
  areas%ix_are(i, 5) = areas%ny
  areas%ix_are(i, 6) = config%ioya
  areas%ix_are(i, 7) = 1
  areas%ix_are(i, 8) = areas%nz
  areas%ix_are(i, 9) = config%ioza
  areas%ix_are(i,10) = 1
  areas%ix_are(i,11) = 1
  areas%ix_are(i,12) = config%bndmnx
  areas%ix_are(i,13) = config%bndmxx ! will be reset if are_nu > 1
  
  i = 2
  areas%ix_are(i, 1) = nxa + 1
  areas%ix_are(i, 2) = nxb
  areas%ix_are(i, 3) = 1
  areas%ix_are(i, 4) = 1
  areas%ix_are(i, 5) = areas%ny
  areas%ix_are(i, 6) = config%ioyb
  areas%ix_are(i, 7) = 1
  areas%ix_are(i, 8) = areas%nz
  areas%ix_are(i, 9) = config%iozb
  areas%ix_are(i,10) = 1
  areas%ix_are(i,11) = 1
  areas%ix_are(i,12) = -100        ! will be set in "setinf"
  areas%ix_are(i,13) = config%bndmxx ! will be reset if are_nu > 2
  
  i = 3
  areas%ix_are(i, 1) = nxb + 1
  areas%ix_are(i, 2) = nxc
  areas%ix_are(i, 3) = 1
  areas%ix_are(i, 4) = 1
  areas%ix_are(i, 5) = areas%ny
  areas%ix_are(i, 6) = config%ioyc
  areas%ix_are(i, 7) = 1
  areas%ix_are(i, 8) = areas%nz
  areas%ix_are(i, 9) = config%iozc
  areas%ix_are(i,10) = 1
  areas%ix_are(i,11) = 1
  areas%ix_are(i,12) = -100        ! will be set in "setinf"
  areas%ix_are(i,13) = config%bndmxx
  
    do i = 2, areas%are_nu
     if(areas%ix_are(i,2) .le. areas%ix_are(i-1,2)) then
        raise_abort("input(): order of nxabc")
     endif
  enddo
#else
  ! the array ix_are is allocated and initialised later in the remap
  ! case
#endif

  end subroutine init_areas

!>
!> \par This subroutine allocates the arrays from module arecon_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_hydro_areas(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
    
    allocate(areas%dt_are(3), &
             areas%dt_cfl(3), &
             areas%ti_are(3), stat=istat)
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module hydro_areas_mod failed")
    end if
    
    mem_local = 3*8._rk*3._rk

    allocate(areas%ix_are(3,13), stat=istat)
    if (istat .ne. 0) then
        raise_abort("allocate_hydro_arrays(): allocation of module hydro_areas_mod failed")
    end if
   
    mem_local = mem_local + 3*8._rk*13._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "hydro_areas")

    end subroutine allocate_hydro_areas

!>
!> \par This subroutine deallocates the arrays from module arecon_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_hydro_areas
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat
    deallocate(areas%dt_are, areas%dt_cfl, &
               areas%ti_are, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module hydro_areas failed")
    end if

    deallocate(areas%ix_are, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module hydro_areas failed")
    end if
  end subroutine deallocate_hydro_areas
      
end module hydro_areas_mod

!-----------------------------------------------------------------------

!>     quantities of marker particles:
!>
!>       npar      maximum number of marker particles
!>       nop       current number of particles
!>
!>       ppx(n)    x-coordinate of particle n                           
!>       ppy(n)    y-coordinate of particle n                           
!>       ppz(n)    z-coordinate of particle n                           
!>                                                                      
!>       pvx(n)    x-velocity of at the position of particle n
!>       pvy(n)    y-velocity of at the position of particle n
!>       pvz(n)    z-velocity of at the position of particle n
!>                 always at the last time step
!>
!>       pma(n)    mass of the particle n
!>
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module marker_hy


  use precision
  
  implicit none
! LOCAL varibales that are not in modules

  SAVE
  integer(kind=ik) ::     nop
  integer(kind=ik), parameter ::npar = 1                            

  real(kind=rk), allocatable ::     ppx(:), ppy(:), ppz(:), pvx(:), pvy(:), &
                       pvz(:), pma(:)

contains 
!>
!> \par This subroutine allocates the arrays from module marker_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_marker_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
    

    allocate(ppx(npar), ppy(npar), ppz(npar), pvx(npar), pvy(npar), &
                       pvz(npar), pma(npar), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module marker_hy failed")
    end if

    mem_local = npar * 7._rk * 8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "marker_hy")

  end subroutine allocate_marker_hy

!>
!> \par This subroutine deallocates the arrays from module marker_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_marker_hy
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat

    deallocate(ppx, ppy, ppz, pvx, pvy, pvz, pma, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module marker_hy failed")
    end if

  end subroutine deallocate_marker_hy

end module marker_hy

!-----------------------------------------------------------------------


!>
!> \par This module provides the fluxes use in the hydro code
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
module massio_hy
  use precision

  implicit none
! LOCAL varibales that are not in modules

  SAVE

!-----------------------------------------------------------------------
!
!     dmdtio = rhoflx (1/ config%qx+1) (ohne 4 Pi) [g]
!
!
!-----------------------------------------------------------------------


  real(kind=rk), allocatable :: dmdtio(:,:,:)

  real(kind=rk), allocatable, dimension(:,:) ::  dflxl,dflxr,eflxl,eflxr, &
                                         vxflxl,vxflxr,vyflxl,    &
                                         vyflxr,vzflxr,vzflxl


  real(kind=rk), allocatable, dimension(:,:,:) :: xnflxl,xnflxr


  !real(kind=rk), dimension(2,qy,qz_s:qz_e,are_di) :: dflxtot,eflxtot,vxflxtot, &
  !                                                  vyflxtot,vzflxtot

  real(kind=rk), allocatable, dimension(:,:,:,:) :: dflxtot,    &
                                  eflxtot, vxflxtot, vyflxtot, vzflxtot

  real(kind=rk), allocatable :: xnflxtot(:,:,:,:,:)



contains 
!>
!> \par This subroutine allocates the arrays from module massio_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_massio_hy(mem_global)
    use precision

    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use hydro_areas_mod

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
    
    allocate(dflxl(qy_s:qy_e,qz_s:qz_e), dflxr(qy_s:qy_e,qz_s:qz_e),&
           eflxl(qy_s:qy_e,qz_s:qz_e), eflxr(qy_s:qy_e,qz_s:qz_e), &
           vxflxl(qy_s:qy_e,qz_s:qz_e), vxflxr(qy_s:qy_e,qz_s:qz_e),&
           vyflxl(qy_s:qy_e,qz_s:qz_e), vyflxr(qy_s:qy_e,qz_s:qz_e), &
           vzflxr(qy_s:qy_e,qz_s:qz_e), vzflxl(qy_s:qy_e,qz_s:qz_e),&
           xnflxl(qy_s:qy_e,qz_s:qz_e,config%qn), xnflxr(qy_s:qy_e,qz_s:qz_e,config%qn), &
           dflxtot(2,qy_s:qy_e,qz_s:qz_e,3),         &
           eflxtot(2,qy_s:qy_e,qz_s:qz_e,3),         &
           vxflxtot(2,qy_s:qy_e,qz_s:qz_e,3), &
           vyflxtot(2,qy_s:qy_e,qz_s:qz_e,3), &
           vzflxtot(2,qy_s:qy_e,qz_s:qz_e,3),&
           xnflxtot(2,qy_s:qy_e,qz_s:qz_e,config%qn,3), dmdtio(6,qy_s:qy_e,qz_s:qz_e), &
           stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module massio_hy4 failed")
    end if

    mem_local = (qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk * 10._rk
    mem_local = mem_local +  (qy_e-qy_s+1)*(qz_e-qz_s+1)*config%qn*8._rk*2._rk
    mem_local = mem_local +  2*(qy_e-qy_s+1)*(qz_e-qz_s+1)*3*8._rk * 6._rk
    mem_local = mem_local +  6*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "massio_hy")

  end subroutine allocate_massio_hy
!>
!> \par This subroutine deallocates the arrays from module massio_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_massio_hy
    use precision
    use abort

    implicit none
    integer(kind=ik) :: istat



    deallocate(dflxl, dflxr, eflxl, eflxr, vxflxl, vxflxr, vyflxl, &
               vyflxr, vzflxr, vzflxl, xnflxl, xnflxr, dflxtot,    &
               eflxtot, vxflxtot, vyflxtot, vzflxtot, xnflxtot,    &
               dmdtio, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module massio_hy4 failed")
    end if


  end subroutine deallocate_massio_hy


end module massio_hy

 
!-----------------------------------------------------------------------

!>
!> \par This module provides the quantities needen if a reverse shock is present
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module revsho_hy
  use precision


  implicit none
! LOCAL varibales that are not in modules

  SAVE

  integer(kind=ik) ::    nw, n_a

  real(kind=rk) :: rcut, rho0_w, rho1_w, t0_w, t1_w, v0_w, v1_w, rmw, &
                       tim1_w, tim2_w, s0_w,  timsf, t_max,   &
                       del_t, del_w, ome_w,       &
                       dt_o, v_5
                       
#ifdef PROGRAM_remap
  real(kind=rk), allocatable, dimension(:) :: v_o,  d_o, p_o, &
                                              e_o, gc_o, ge_o, &
                                              d_1, d_2, d_3, d_4, d_5
#endif
  real(kind=rk), allocatable :: tsave(:,:)

contains 
!>
!> \par This subroutine allocates the arrays from module revsho_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_revsho_hy(mem_global)
    use precision
    use print_stdout_mod, only : print_memory_alloc
    use abort

    use configure
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

#ifdef PROGRAM_remap
    allocate( v_o(config%q), d_o(config%q), p_o(config%q),     &
              e_o(config%q), gc_o(config%q), ge_o(config%q),   &
              d_1(config%q), d_2(config%q), d_3(config%q),     &
              d_4(config%q), d_5(config%q), stat=istat)
#endif

    allocate(tsave(4,config%qy), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module revsho_hy 1 failed")
    end if

!    mem_local = q*8._rk*11._rk + 4*qy*8._rk
    mem_local = 4*config%qy*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "revsho_hy")

  end subroutine allocate_revsho_hy
!>
!> \par This subroutine deallocates the arrays from module revsho_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_revsho_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

!    deallocate( v_o, d_o, p_o, e_o, gc_o, ge_o, &
!                d_1, d_2, d_3, d_4, d_5, tsave, stat=istat)
    deallocate(tsave, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module revsho_hy 1 failed")
    end if

  end subroutine deallocate_revsho_hy

end module revsho_hy


!-----------------------------------------------------------------------



!>     "bndinf" is used in:  bndry
!>
!>       contains all quantities which are needed for the boundaries
!>       between different PPM-calculation areas in x-direction
!>
!>       a(x,j,k):  x = 1....4 -> inner boundary of area 1
!>                 (x = 5....8 -> outer boundary of area 0) -> dummy
!>                  x = 9...12 -> inner boundary of area 2
!>                  x = 13..16 -> outer boundary of area 1
!>                  x = 17..20 -> inner boundary of area 3
!>                  x = 20..24 -> outer boundary of area 2
!>                 (x = 25..28 -> inner boundary of area 4) -> dummy
!>                  x = 20..24 -> outer boundary of area 2
!>
!>       are_id:    ID of the calculation area
!>                  index of inner boundary of area "are_id"
!>                  ->   (are_id - 1)*8 + i     (i = 1..4)
!>                  index of outer boundary of area "are_id"
!>                  ->    are_id*8  + 4 + i     (i = 1..4)
!>
!>
!> \author M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module bndinf_hy
  use precision

  implicit none
! LOCAL varibales that are not in modules
  SAVE

  integer(kind=ik) :: are_id
  integer(kind=ik) ,parameter::   qbi = 32

  real(kind=rk), allocatable, dimension(:,:,:) :: &
                 vex_bi, vey_bi, vez_bi, den_bi, ene_bi, pre_bi, gae_bi, &
                 gac_bi, gra_bi

  real(kind=rk), allocatable ::dxx_bi(:), ugr_bi(:), xnu_bi(:,:,:,:)

contains
!>
!> \par This subroutine allocates the arrays from module bndinf_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_bndinf_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(dxx_bi(qbi), ugr_bi(qbi), stat=istat)
    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module revsho_hy 1 failed")
    end if
    
    mem_local = qbi*8._rk*2._rk


    allocate( vex_bi(qbi,qy_s:qy_e,qz_s:qz_e), vey_bi(qbi,qy_s:qy_e,qz_s:qz_e), &
              vez_bi(qbi,qy_s:qy_e,qz_s:qz_e), den_bi(qbi,qy_s:qy_e,qz_s:qz_e),      &
              ene_bi(qbi,qy_s:qy_e,qz_s:qz_e), pre_bi(qbi,qy_s:qy_e,qz_s:qz_e),      &
              gae_bi(qbi,qy_s:qy_e,qz_s:qz_e), gac_bi(qbi,qy_s:qy_e,qz_s:qz_e),      &
              gra_bi(qbi,qy_s:qy_e,qz_s:qz_e), xnu_bi(qbi,qy_s:qy_e,qz_s:qz_e,config%qn),   &
              stat=istat)
    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module revsho_hy 1 failed")
    end if

    mem_local = mem_local + qbi*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*9._rk + &
                            qbi*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%qn*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "bndinf_hy")

  end subroutine allocate_bndinf_hy

!>
!> \par This subroutine deallocates the arrays from module bndinf_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_bndinf_hy
    use precision
    use abort

    implicit none

    integer(kind=ik) :: istat

    deallocate(dxx_bi, ugr_bi, stat=istat)
    if (istat .ne. 0) then
       raise_abort("deallocate_hydro_arrays(): deallocation of module revsho_hy 1 failed")
    end if


    deallocate( vex_bi, vey_bi, vez_bi, den_bi, ene_bi, pre_bi, gae_bi, & 
                gac_bi, gra_bi, xnu_bi, stat=istat)
    if (istat .ne. 0) then
       raise_abort("deallocate_hydro_arrays(): deallocation of module revsho_hy 1 failed")
    end if



  end subroutine deallocate_bndinf_hy


end module bndinf_hy
!-----------------------------------------------------------------------

!> \verbatim
!>    quantities of the total calculation area at very very old timelevel
!>
!>
!>  \author M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module ancient_hy
  use precision
  
!  use arecon_hy

  use mo_mpi

  use hydro_areas_mod
  implicit none
! LOCAL varibales that are not in modules

  SAVE

!-----------------------------------------------------------------------
!     "ancient" is used in:
!     puthyold, gethyold
!-----------------------------------------------------------------------


  real(kind=rk), allocatable, dimension(:,:,:)::  acxanc, vexanc, &
                                                  veyanc, vezanc, &
                                                  denanc, temanc, &
                                                  gpoanc, qenanc, &
                                                  qmoanc
  real(kind=rk), allocatable, dimension(:,:,:,:) :: qyeanc
 
  real(kind=rk), allocatable :: gamanc(:), xnuanc(:,:,:,:)


  real(kind=rk),allocatable, dimension(:) :: xzlanc,  xzranc,yzlanc,  yzranc , &
                                             zzlanc,  zzranc 

  integer(kind=ik) :: nhystpanc
  integer(kind=ik), allocatable :: ix_areanc(:,:)
!>
!> \par This subroutine allocates the arrays from module ancient_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
contains

  subroutine allocate_ancient_hy(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    use hydro_areas_mod
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
   
    allocate(acxanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             vexanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             veyanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             vezanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             denanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             temanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             gpoanc(0:config%qx,qy_s-1:qy_e,qz_s:qz_e),    &
             gamanc(0:config%qx),                          &
             xnuanc(config%qx,qy_s:qy_e,qz_s:qz_e,config%qn), stat=istat)

    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module ancient_hy 1 failed")
    end if

    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk * 6._rk + (config%qx+1)*(qy_e-qy_s+2)*(qz_e-qz_s+1)*8._rk + &
                (config%qx+1)*8._rk + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*config%qn*8._rk


    allocate(qyeanc(config%qx,qy_s:qy_e,qz_s:qz_e,5),      &
             qenanc(config%qx,qy_s:qy_e,qz_s:qz_e),        &
             qmoanc(config%qx,qy_s:qy_e,qz_s:qz_e),  stat=istat)

    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module ancient_hy 2 failed")
    end if

    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk * 7._rk 


    allocate(xzlanc(config%q_nqx), xzranc(config%q_nqx),        &
             yzlanc(config%q_nqy), yzranc(config%q_nqy),        &
             zzlanc(config%q_nqz), zzranc(config%q_nqz),        &
             ix_areanc(3, 13),  &
             stat=istat)


    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module ancient_hy 3 failed")
    end if
    mem_local = mem_local + config%q_nqx*8._rk*2._rk +config%q_nqy*8._rk*2._rk+ &
                            config%q_nqz*8._rk*2._rk+ 3*13*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "ancient_hy")

  end subroutine allocate_ancient_hy

!>
!> \par This subroutine deallocates the arrays from module ancient_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_ancient_hy
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat

    deallocate(acxanc, vexanc, veyanc, vezanc, denanc, temanc, &
               gpoanc, gamanc, xnuanc, stat=istat)

    if (istat .ne. 0) then
       raise_abort("deallocate_hydro_arrays(): deallocation of module ancient_hy 1 failed")
    end if

    deallocate(qyeanc, qenanc,  qmoanc,  stat=istat)

    if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module ancient_hy 2 failed")
    end if


    deallocate(xzlanc, xzranc, yzlanc, yzranc, zzlanc, zzranc, ix_areanc, &
               stat=istat)


    if (istat .ne. 0) then
         raise_abort("deallocate_hydro_arrays(): deallocation of module ancient_hy 3 failed")
    end if


  end subroutine deallocate_ancient_hy
end module ancient_hy
!-----------------------------------------------------------------------
!> \par fluxes over outer boundary
!>
!> \author M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module ioflx

  use precision

  implicit none
! LOCAL varibales that are not in modules

  SAVE
  
  real(kind=rk) :: sumioe,sumion,sumiom,sumioe_n,sumion_n

end module ioflx

!> \par changes of some quantities due to transport
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module dquants
  use precision

  implicit none
  
  SAVE

  real(kind=rk) :: delyem_o,delden_o,deltem_o,deleni_o,sigma_o

end module dquants

module cnsv
  use precision

  implicit none
! LOCAL varibales that are not in modules

  SAVE

  real(kind=rk) :: sloss(2),dtot(2),slep(2)
  real(kind=rk) :: eges_old,eges_old2,nges_old,nges_old2,tten_old,elto_old

end module cnsv

!> \par perturb model
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module perturbation
  use precision

  implicit none
! LOCAL varibales that are not in modules

  save
!  integer(kind=ik) :: noise
!  real(kind=rk)    :: ampl
end module perturbation

!> \par moving grid
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module movgrid_hy
  use precision

  implicit none
! LOCAL varibales that are not in modules

  save
  integer(kind=ik) :: nstr,nstrt,nxt
  real(kind=rk)    :: rstr

  contains

  subroutine init_movegrid

  use precision
  use hydro_areas_mod
  implicit none

  integer(kind=ik) ::  idltn

    ! --- set constants for moving grid:
  nstr  = areas%nx / 2
  ! --- dx(1) ~ (idltn+1) * dx(2)
  !      idltn  = 2
  idltn  = 0
  nstrt = nstr    + idltn
  nxt   = areas%nx   + idltn

  end subroutine init_movegrid
  ! xmfrac
end module movgrid_hy

!-----------------------------------------------------------------------


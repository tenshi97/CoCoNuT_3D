module rady_mod

implicit none

contains

#undef DEBUG

subroutine rady(dt_cyc, ti_cyc, selftime, childrentime)

!----------------------------------------------------------------------
! Autor           : Wolfgang Keil, Markus Rampp (MPA)
! Modul           : $Id: rady.F,v 1.15 2006/05/11 15:44:04 rburas Exp $
! Version         : $Revision: 1.15 $
! Date            : $Date: 2006/05/11 15:44:04 $
! Freeze Version  : $Name:  $
!
!
!     task:      controls time step of the various calculation areas,
!                calculates all areas.
!
!     input:     ti_cyc = time at the beginning of the cycle = t_0
!
!     output:    dt_cyc = time step of the cycle
!                ti_cyc = time at the end of the cycle
!                       = t_0 + dt_cyc
!
!----------------------------------------------------------------------

  use precision
  use abort
  use error

!  use arecon_hy
  use massio_hy
  use intgrs_hy
  use nutrio_hy


  use param_rt
#ifndef NOTRA
  use multigrid_rt
  use neutrinotypes
#endif

  use cputim


  use ioflx
  use specfun
  use totare_hy


  use sum_fluxes
  use savare_overload
  use hydro_interface

#ifdef CFC_TRANSPORT
  use size_cfc
  use hydro_primitives_cfc
  use phycon, only: pc_geog
#endif

  use save_old_hydro_state_mod

  use timestep_areas, only : minimum_timestep_of_areas

  use mo_mpi

#ifdef CHECK_MEMORY
  use meminfo
#endif


  use grid_mod
  use grids
  use print_stdout_mod

  use cputim
  use configure
  use state
  use hydro_areas_mod
  use state

  use neutrino_mod

  use transp_timestep, only : transport_timestep_control
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(inout) :: dt_cyc, ti_cyc

  real(kind=rk), intent(out)  :: selftime(2), childrentime(2)
  real(kind=rk), dimension(2) :: selftime_start(2)
  real(kind=rk), dimension(2) :: &
                                 nusource_self, nusource_children,     &
                                 neutrino_self, neutrino_children,     &
                                 save_hydro_self, save_hydro_children, &
                                 minimum_self, minimum_children,       &
                                 savare_self, savare_children,         &
                                 timestep_self, timestep_children,     &
                                 do_hydro_self, do_hydro_children

  real(kind=rk), dimension(2) :: tim1, tim2, tim4

  real(kind=rk), save         :: dt_hyd

  integer(kind=ik) :: i,j,k,jj,kk,jk
  logical ::         labb,  loutp, err_transp

  real(kind=rk), dimension(config%qx,qy_s:qy_e,qz_s:qz_e) :: tehtot, yehtot, enhtot, dehtot

#ifdef RAND_DT_RAD
  real(kind=rk)    :: rand
#endif
  integer(kind=ik) :: is

  integer(kind=ik) :: nray(2)

  integer(kind=ik) :: ierr

  integer(kind=ik) :: idtmin, isdmin, ndt, itnum, nenum, nhyrep
  integer(kind=ik) :: ii
  real(kind=rk)    :: sigma, dt_min, dt_nex

!c test only
!      real(kind=rk), parameter :: dyemx=7.0e-3,dtemx=8.0e-3
!      real(kind=rk), parameter :: drhmx=5.0e-3,deimx=5.0e-3



  selftime          = 0._rk
  childrentime      = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  nhyrep = 0
  labb = .FALSE.


  if(config%itstp .ne. 0) then
     if(mod(nstep,config%itstp) .eq. 0) then
        call printit_taskX(0,"================================================================================")
     end if
  end if

  call save_old_hydro_state(dehtot, tehtot, yehtot, enhtot, save_hydro_self, &
                                                            save_hydro_children)
  childrentime = childrentime + save_hydro_self


!-----------------------------------------------------------------------
!     search area with the shortest time step:
!-----------------------------------------------------------------------

  call minimum_timestep_of_areas(idtmin, isdmin, dt_nex, dt_cyc, dt_min, &
                                 dt_hyd, ndt, minimum_self, minimum_children)
  childrentime = childrentime + minimum_self

  call savare(1, savare_self, savare_children)
  childrentime = childrentime + savare_self

10 continue                  ! backjump if TRANSP fails

  areas%ix_are(1:areas%are_nu,11) = ndt
  areas%dt_are(1:areas%are_nu)    = min(dt_min,transport%dt)

!-----------------------------------------------------------------------
!     calculate ndt time steps; all areas with dt = dt_min:
!-----------------------------------------------------------------------

  dt_nex     = 0._rk     ! reset dt_nex to calculate the real dt_nex
!      ti_are(ia) = 0.     ! if the boundaries need time information

#ifdef CHECK_MEMORY
  if (meminfo_flag) call meminfo_start("HYDRO total")
#endif

  call do_hydro_steps(dt_min, dt_nex, ti_cyc, idtmin, do_hydro_self, do_hydro_children)
  childrentime = childrentime + do_hydro_self


#ifdef CHECK_MEMORY
  if (meminfo_flag) call meminfo_stop("HYDRO total")
#endif

!  neutrino transport can be switched on and off (p_ntr)
!  neutrino backreaction can be switched off     (p_nbk)
!      (idea for p_ntr=1,p_nbk=0:
!           neutrino-sourceterms are calculated and output
!           but not applied to the hydrodynamics)
  if (config%p_ntr .gt. 0) then

     call neutrino(dt_nex,nhyrep,err_transp, itnum, nenum, sigma, nstep, ti_cyc, labb, loutp, nray(:), ndt, dt_min, &
                   neutrino_self, neutrino_children)
     timer%rady_neutrino = timer%rady_neutrino + neutrino_self
     timer%rady_neutrino_children = timer%rady_neutrino_children + neutrino_children
     childrentime = childrentime + neutrino_self
     timer%transp_tot = timer%rady_neutrino

     if (err_transp) then
        if (nhyrep .le. 6) then
           nhyrep = nhyrep + 1
           goto 10
        else
           raise_error("rady(): too many repetitions of transport ")
!               call stopit('too many repetitions of transport ',0)
        endif
     endif
  else
     sigma = 0.0_rk
  endif

  areas%nhystp = areas%nhystp + ndt     ! hydro time step #


#ifdef CFC_TRANSPORT
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL DO &
!$OMP SCHEDULE(static) &
!$OMP PRIVATE(j,k,jk) &
!$OMP SHARED(n_loc,o_loc,n_s,o_s,dentot,xnutot,temtot,enetot,rho,xnnu,t,eps)
#endif
  do jk = 1, n_loc * o_loc
     k = int((jk + n_loc - 1) / n_loc )
     j = (n_s - 1) + (jk - (k - 1) * n_loc)
     k = k + o_s - 1

     dentot(1:config%qx,j,k)   =rho (1:config%qx,j,k)*pc_geog
     temtot(1:config%qx,j,k)   =t   (1:config%qx,j,k)
     xnutot(1:config%qx,j,k,config%qn)=xnnu(1:config%qx,j,k,config%qn)
     enetot(1:config%qx,j,k)   =eps (1:config%qx,j,k)
  end do
#endif
  !------TIMESTEP CONTROL FOR THE MULTI-TIMESTEPPING
  call transport_timestep_control(yehtot, dehtot, tehtot, enhtot, dt_hyd, &
                                  nhyrep, dt_nex,dt_cyc, ti_cyc, itnum,   &
                                  nenum, sigma, nray(:), dt_min, ndt,     &
                                  timestep_self, timestep_children)
!  childrentime = childrentime + timestep_self

#ifndef DEBUG_TIMINGS
 call second_v(selftime)

 tim2 = selftime
 selftime = selftime - selftime_start
#endif


 return

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

end subroutine rady

end module rady_mod

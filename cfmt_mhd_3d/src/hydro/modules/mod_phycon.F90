!>
!> \verbatim
!> This module provides the physcial constants used in VERTEX
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
module phycon

!-----------------------------------------------------------------------
!     PHYSICAL CONSTANTS
!
!     reference: Review of Particle Properties, Phys. Rev. D45, 1992
!-----------------------------------------------------------------------
  use precision

  implicit none
! LOCAL variables thta are not in modules

  SAVE 

  real(kind=rk), parameter :: pc_cl = 2.99792458e+10_rk  ! speed of light [cm/s]
  
  real(kind=rk), parameter :: pc_gc = 6.67259e-08_rk 
                                    ! gravitational constant [cc/g/s**2]

  real(kind=rk), parameter :: pc_kb = 1.380658e-16_rk 
                                    ! Boltzmann constant [erg/K]

  real(kind=rk), parameter :: pc_mb = 1.66e-24_rk        ! baryon mass [g]

  real(kind=rk), parameter :: pc_mn = 1.6749286e-24_rk ! neutron mass [g]

  real(kind=rk), parameter :: pc_mp = 1.6726231e-24_rk ! proton mass [g]

  real(kind=rk), parameter :: pc_mpi = 2.405911e-25_rk        ! pion mass [g]

  real(kind=rk), parameter :: pc_me = 9.1093897e-28_rk  ! electron mass [g]
  
  real(kind=rk), parameter :: pc_ms = 1.98892e+33_rk ! solar mass [g]

  real(kind=rk), parameter :: pc_pi = 3.141592653589793_rk ! \pi

  real(kind=rk), parameter :: pc_meverg = 1.602e-6_rk
                                     ! transformation from MeV to erg

  real(kind=rk), parameter :: pc_ergmev = 1.0_rk/pc_meverg   !(= 6.242197e5)
                                     ! transformation from erg to MeV 

  real(kind=rk), parameter :: pc_mevk = pc_meverg/pc_kb   !(= 1.16032e10)
                                     ! transformation from MeV to K

  real(kind=rk), parameter :: pc_kmev = 1.0_rk/pc_mevk      !(= 8.61829e-11)
                                     ! transformation from K to MeV

  real(kind=rk), parameter :: pc_ggeo = pc_gc/pc_cl**2
                                    ! mass: g -> geometrical units
  real(kind=rk), parameter :: pc_geog = 1.0_rk/pc_ggeo
  ! mass: geometrical units -> g
  real(kind=rk), parameter :: pc_egeo = pc_gc/pc_cl**4
                                     ! energy: erg -> geometrical units
  real(kind=rk), parameter :: pc_geoe = 1.0_rk/pc_egeo
                                     ! energy: geometrical units -> erg
  real(kind=rk), parameter :: pc_bgeo = sqrt (pc_egeo/4.0_rk/pc_pi)
                                     ! B-field: G -> geometrical units
  real(kind=rk), parameter :: pc_geob = 1.0_rk/pc_bgeo
                                     ! B-fieldy: geometrical units -> G
!-----------------------------------------------------------------------
! weak interaction constants

  real(kind=rk), parameter :: wc_ga2 = 1.254_rk**2   !g_a**2

  real(kind=rk), parameter :: wc_s2w = 0.23_rk       !sin^2(theta_W)  

  real(kind=rk), parameter :: wc_ca = 0.5_rk
  real(kind=rk), parameter :: wc_cv=0.5_rk+2._rk*wc_s2w

  real(kind=rk), parameter :: wc_s0 = 1.761e-44_rk   !sigma_0 [cm^2]

  real(kind=rk), parameter :: wc_hc = 197.327e-13_rk !hquer*c [MeV*cm]

  real(kind=rk), parameter :: wc_me = pc_me*pc_cl*pc_cl*pc_ergmev
                                             !electron rest energy [MeV]
  real(kind=rk), parameter :: wc_mp = pc_mp*pc_cl*pc_cl*pc_ergmev
                                             !proton rest energy [MeV]
  real(kind=rk), parameter :: wc_mn = pc_mn*pc_cl*pc_cl*pc_ergmev
                                             !neuton rest energy [MeV]
  real(kind=rk), parameter :: wc_mpi = pc_mpi*pc_cl*pc_cl*pc_ergmev
                                             !pion rest energy [MeV]
  real(kind=rk), parameter :: wc_mb = pc_mb*pc_cl*pc_cl*pc_ergmev
                                             !baryon rest energy [MeV]
  real(kind=rk), parameter :: wc_mq = wc_mn-wc_mp
                                       !rest-energy difference n-p [MeV]

!-----------------------------------------------------------------------
! coupling constants

!-- coupling constants for (anti/)neutrino-electron/positron scattering

!    (see e.g. Janka, Th. PHD p.69)

! scattering of (nu_e,nu_ebar,nu_mu/tau,nu_mu/taubar) off electrons
  real(kind=rk), parameter, dimension(4) :: cc1_nes=                &
                                       (/ (1._rk+2._rk*wc_s2w)**2,  &
                                          (    2._rk*wc_s2w)**2,    &
                                          (2._rk*wc_s2w-1._rk)**2,  &
                                       (    2._rk*wc_s2w)**2 /)
  real(kind=rk), parameter, dimension(4) :: cc2_nes=                &
                                       (/ (    2._rk*wc_s2w)**2,    &
                                         (1._rk+2._rk*wc_s2w)**2,   &
                                          (    2._rk*wc_s2w)**2,    &
                                         (2._rk*wc_s2w-1._rk)**2 /)
! scattering of (nu_e,nu_ebar,nu_mu/tau,nu_mu/taubar) off positrons
  real(kind=rk), parameter, dimension(4) :: cc1_nps=                &
                                        (/ (    2._rk*wc_s2w)**2,   &
                                          (1._rk+2._rk*wc_s2w)**2,  &
                                           (    2._rk*wc_s2w)**2,   &
                                         (2._rk*wc_s2w-1._rk)**2 /)
  real(kind=rk), parameter, dimension(4) :: cc2_nps=                &
                                       (/ (1._rk+2._rk*wc_s2w)**2,  &
                                           (    2._rk*wc_s2w)**2,   &
                                          (2._rk*wc_s2w-1._rk)**2,  &
                                          (    2._rk*wc_s2w)**2 /)

!-- coupling constants for (anti/)neutrino- (anti/)neutrino scattering
!    (see e.g. Buras et al 2002, Table 1)

! scattering of (nu_e,nu_ebar,nu_mu/tau,nu_mu/taubar) off nu_e
  real(kind=rk), parameter, dimension(4) :: cc1_nns=                 &
                                       (/ 0._rk,0._rk,1._rk,0._rk /)
  real(kind=rk), parameter, dimension(4) :: cc2_nns=                 &
                                       (/ 0._rk,0._rk,0._rk,1._rk /)

! scattering of (nu_e,nu_ebar,nu_mu/tau,nu_mu/taubar) off nu_ebar
  real(kind=rk), parameter, dimension(4) :: cc1_nas=                 &
                                       (/ 0._rk,0._rk,0._rk,1._rk /)
  real(kind=rk), parameter, dimension(4) :: cc2_nas=                 &
                                       (/ 0._rk,0._rk,1._rk,0._rk /)


end module phycon

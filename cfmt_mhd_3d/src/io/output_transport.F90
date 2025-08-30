#ifdef WRITE_BINARY_OUTPUT
module write_transport

  private
  public :: write_transport_output

  contains


!>
!> This subroutine writes the transport output files
!> 
!> \author M. Rampp, MPA
!>
!> \verbatim
!>
!>     xout(evens,:,:,:,:)      = log10(J)                ME
!>     xout(odds ,:,:,:,:)      = H/J   (versetzt !!!)    ME
!>     wout(evens,:,:,:,:)      = log10(N)                ME
!>     wout(odds ,:,:,:,:)      = F/N   (versetzt !!!)    ME
!>     elag(evens,:,1:isma,:,:) = log10(J)                BTE
!>     elag(odds ,:,1:isma,:,:) = H/J   (versetzt !!!)    BTE
!>     elag(evens,:,isma:*,:,:) = K/J   (zone centers)    BTE
!>     elag(odds ,:,isma:*,:,:) = L/J   (zone interfaces) BTE
!>
!>   1D-case :             config%nystrt =1 , nytra =1  , config%nymom =1
!>   2D-case (1D-Edd-Fak): config%nystrt =0,  nytra =0  , config%nymom =config%qy
!>           (2D-Edd-Fak): config%nystrt =1,  nytra =config%qy, config%nymom =config%qy
!>
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 668 $
!>   $Date: 2010-01-15 17:27:39 +0100 (Fri, 15 Jan 2010) $
!> \endverbatim
!>
subroutine write_transport_output
  use precision

  use abort

  use intgrs_hy
 ! use arecon_hy
  use radial_grid_rt
  use multigrid_rt

  use mo_mpi
  use lagradq_rt, only : xlag, elag, wlag
  use lagback_rt, only : wwlag


  use param_rt
  use time_rt
  use specfun

  use dataformat_vertex
  use energy_grid_rt, only : ener, emid
  use totare_hy, only: yzntot, zzntot

  use configure
  use state
  implicit none

  call write_timestep_marker(pfile)

  call writeq(pfile, "ener", ener(1:config%iemax+1), 'energy bin boundaries', 'MeV', 'boundary:'//frange(config%iemax+1))
  call writeq(pfile, "emid", emid(1:config%iemax), 'energy bin center', 'MeV', 'bin:'//frange(config%iemax))

  call writeq(pfile, "jmin", jmin(config%nystrt:config%nymom), &
    "?", &
    "1", &
    "zone:"//frange(config%nystrt, config%nymom))

  call writeq(pfile, "jmax", jmax(config%nystrt:config%nymom), &
    "?", &
    "1", &
    "zone:"//frange(config%nystrt, config%nymom))

  call writeq(pfile, "kmin", kmin(1:config%nztra), &
    "?", &
    "1", &
    "zone:"//frange(1, config%nztra))

  call writeq(pfile, "kmax", kmax(1:config%nztra), &
    "?", &
    "1", &
    "zone:"//frange(1, config%nztra))

  if (config%use_multid_collapse) then
     call writeq(pfile, "ralag", ralag_3d(:,config%nystrt,1), &
    'radius of zone boundaries of the transport grid. Attention! Only one slice is used', 'cm', &
    'zone:'//frange(-1,config%imaxp + 1))
     call writeq(pfile, "rqlag", (0.5_rk * (ralag_3d(-1:config%imaxp,config%nystrt,1)**3 + &
                                         ralag_3d( 0:config%imaxp +1,config%nystrt,1)**3    ))**(1./3.), &
    'radius of zone center of the transport grid. Attention! Only one slice is used', 'cm', &
    'zone:'//frange(0,config%imaxp +1))
  else
     call writeq(pfile, "ralag", ralag_1d(:), &
    'radius of zone boundaries of the transport grid.', 'cm', &
    'zone:'//frange(-1,config%imaxp + 1))
     call writeq(pfile, "rqlag", (0.5_rk * (ralag_1d(-1:config%imaxp)**3 + &
                                         ralag_1d( 0:config%imaxp +1)**3    ))**(1./3.), &
    'radius of zone center of the transport grid. Attention! Only one slice is used', 'cm', &
    'zone:'//frange(0,config%imaxp +1))
  endif ! mutlid_colapse
    
  !> \todo Angular grid is equal to hydro grid at the moment, is this correct, always?
  call writeq(pfile, "yzn", yzntot(1:config%qy), &
    "inclination (theta) of zone center", &
    "radian", &
    "zone:"//frange(config%qy))
  call writeq(pfile, "zzn", zzntot(1:config%qz), &
    "azimuth (phi) of zone center", &
    "radian", &
    "zone:"//frange(config%qz))

  call writeq(pfile, "nstep_rad", nstep, 'transport timestep number', '1')
  
  call writeq(pfile, "dt_rad", transport%dt, 'transport timestep interval', 's')

  call write_5d_radiation_quantity("J", xlag( 0::2,:,:,nymoms:nymome,nzmoms:nzmome), &
    'J, neutrino energy density, from ME', 'MeV/cm^2/s')
    
  call write_5d_radiation_quantity("H", xlag(-1::2,:,:,nymoms:nymome,nzmoms:nzmome), &
    'H, neutrino energy flux, from ME', 'MeV/cm^2/s')

  call write_5d_radiation_quantity("Jn", wlag( 0::2,:,:,nymoms:nymome,nzmoms:nzmome), &
    'J#, neutrino number density, from ME', '1/cm^2/s')
    
  call write_5d_radiation_quantity("Hn", wlag(-1::2,:,:,nymoms:nymome,nzmoms:nzmome), &
    'H#, neutrino number flux, from ME', '1/cm^2/s')

  call write_5d_radiation_quantity("J_bte", elag( 0::2,:,1:config%isma,:,:), &
    'J, neutrino energy density, from BTE', 'MeV/cm^2/s')
    
  call write_5d_radiation_quantity("H_bte", elag(-1::2,:,1:config%isma,:,:), &
    'H, neutrino energy flux, from BTE' , 'MeV/cm^2/s')
    
  call write_5d_radiation_quantity("f_K",   elag( 0::2,:,config%isma+1:,:,:), &
    'Eddington factor f_K = K/J, from BTE', '1')
    
  call write_5d_radiation_quantity("f_L",   elag(-1::2,:,config%isma+1:,:,:), &
    'Eddington factor f_L = L/J, from BTE', '1')

  ! Material coefficients, see mod_lagqua.F90
  call write_5d_radiation_quantity("emla",     wwlag(:,:,:,:,:,0), &
    'Neutrino energy equilibrium density','log10(MeV/cm^2/s)')
    
  call write_5d_radiation_quantity("abla",     wwlag(:,:,:,:,:,1), &
    'Absorption (sum of all rates in the code)','1/cm')
    
  call write_5d_radiation_quantity("ablahlms", wwlag(:,:,:,:,:,6), &
    'Absorption on heavies LMS','1/cm')
    
  call write_5d_radiation_quantity("ablahbrn", wwlag(:,:,:,:,:,7), &
    'Absorption on nucleons Bruenn','1/cm')
    
  call write_5d_radiation_quantity("ablanffn", wwlag(:,:,:,:,:,8), &
    'Absorption on heavies FFN','1/cm')
    
  call write_5d_radiation_quantity("ablankjt", wwlag(:,:,:,:,:,9), &
    'Absorption on nucleons KJT','1/cm')
    
  call write_5d_radiation_quantity("scla",     wwlag(:,:,:,:,:,2), &
    'Scattering opacity (sum of all rates in the code)', '1/cm')
    
  call write_5d_radiation_quantity("sela",     wwlag(:,:,:,:,:,3), &
    'Inelastic scattering opacity (sum of all rates in the code)','1/cm')
    
  call write_5d_radiation_quantity("sela1",    wwlag(:,:,:,:,:,4), &
    'First moment inelastic scattering','MeV/cm^3/s')
    
  call write_5d_radiation_quantity("sela2",    wwlag(:,:,:,:,:,5), &
    'Second moment inelastic scattering','1/cm')

  ! Flush file in order to get readable data in case of a crash later
  call flush_data_file(pfile)

end subroutine write_transport_output

end module write_transport
#endif

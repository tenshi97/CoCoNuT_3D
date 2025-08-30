module configure
  use precision

  implicit none

  private
#ifdef WRITE_BINARY_OUTPUT
  public config, conf_t, initialize_config_values, write_config
#else
  public config, conf_t, initialize_config_values
#endif
  save

  type conf_t
    !
    ! Perl magic is working here!
    !
    ! The following lines are saved in the perl array @type_defs,
    ! and this is later used to construct the contents of the routine
    ! "write_config" further down.
    !
    ! The comments on each field are necessary to get a suitable
    ! "description" and "unit" field for writeq. If they are missing,
    ! this field is not written out
    !
    ! The following lines have to obey this particular format:
    !
    !    type :: fieldname  !x!> some nice description, unit: "some unit"
    !
    ! or
    !
    !    type :: fieldname  !x!> some nice description
    !
    ! it the field has no sensible physical unit (i.e. a string, or flag)
    ! There is nothing magical about this format, it is just used further
    ! down as a marker.
    !
#ifdef WRITE_BINARY_OUTPUT
!PERL-START
#endif
     character(200)   :: calculation                        !x!> calculation mode: possible options are
     real(kind=rk)    :: cpumax                             !x!> maximum time for this run, unit: "s"
     integer(kind=ik) :: nend                               !x!> number of timesteps for this job
     real(kind=rk)    :: tmax                               !x!> maximum time for this job, unit: "s"
     integer(kind=ik) :: irstrt                             !x!> start new sequence, if irstrt = 0
     character(200)   :: basenm                             !x!> base-name of data files
#if defined(MPI_HYDRO) && !(defined(PROGRAM_remap))
    character(7)      :: suffix                             !x!> current suffix of data files
#else
    character(2)      :: suffix                             !x!> current suffix of data files
#endif
#if defined(NEW_OUTPUT_LABELS) || defined(PROGRAM_remap)
    character(8)      :: suffix_nr                          !x!> current suffix of date files
#endif
    integer(kind=ik)  :: nrstrt                             !x!> store restart file after NRSTRT timesteps
    real(kind=rk)     :: trstrt                             !x!> store restart file after TRSTRT seconds, unit: "s"
    character(80)     :: rst_mode                           !x!> restart mode
    logical           :: use_temp_restart                   !x!> If true, restart via (rho, T, X_i), else (rho, e, X_i)
    logical           :: use_stdout_redirect                !x!> redirect stdout and stderr in MPI case
    logical           :: use_deactivate_output              !x!> deactivate writing of binary output/restart
    logical           :: use_deactivate_stdout              !x!> deactivate writing of stdout
    logical           :: use_deactivate_controlFiles        !x!> deactivate writing of control files
    logical           :: use_print_memAlloc                 !x!> print memory allocation information
    logical           :: use_print_mpiDomain                !x!> print information about MPI domains
    integer(kind=ik)  :: nout                               !x!> maximum number of timesteps between output
    real(kind=rk)     :: tout                               !x!> maximum physical time between output, unit: "s"
    real(kind=rk)     :: max_central_density                !x!> stop when central density is reached
    integer(kind=ik)  :: stop_bounce                        !x!> stop at bounce time
    integer(kind=ik)  :: stop_bad_resolution                !x!> stop at bad resolution
    real(kind=rk)     :: bad_resolution_grace_period        !x!> wait duration after restart before checking resolution
    integer(kind=ik)  :: collapse_time_resolve              !x!> print 10 timesteps per central density decade
    integer(kind=ik)  :: itstp                              !x!> print timestep information
    real(kind=rk)     :: cfl                                !x!> courant factor
    real(kind=rk)     :: dtini                              !x!> initial size of the timestep, unit: "s"
    real(kind=rk)     :: dtmin                              !x!> minimum size of hydro timestep, unit: "s"
    real(kind=rk)     :: dtmax                              !x!> maximum size of timestep, unit: "s"
    integer(kind=ik)  :: nriem                              !x!> number of iterations in Riemann solver
    real(kind=rk)     :: cvisc                              !x!> artificial viscosity constant
    real(kind=rk)     :: small                              !x!> cut-off value
    real(kind=rk)     :: smlrho                             !x!> cut-off value for density, unit: "g/cm^3"
    real(kind=rk)     :: smallp                             !x!> cut-off value for pressure, unit: "erg/cm^3"
    real(kind=rk)     :: smalle                             !x!> cut-off value for energy, unit: "erg/g"
    real(kind=rk)     :: smallu                             !x!> cut-off value for velocity, unit: "cm/s"
    real(kind=rk)     :: smallx                             !x!> cut-off value for abundances
    integer(kind=ik)  :: igodu                              !x!> use Godunov method
    integer(kind=ik)  :: nsdim                              !x!> dimensionality of the problem
    real(kind=rk)     :: rib                                !x!> radius of inner grid boundary, unit: "cm"
    real(kind=rk)     :: xmfrac                             !x!> concentrate grid around mass-coord, unit: "g"
    real(kind=rk)     :: gridlx                             !x!> radius of outer grid boundary, unit: "cm"
    real(kind=rk)     :: gridly                             !x!> grid length in y-direction, unit: "pi"
    real(kind=rk)     :: gridlz                             !x!> grid length in z-direction, unit: "pi"
    integer(kind=ik)  :: nx                                 !x!> number of grid points in x-direction
    integer(kind=ik)  :: qx                                 !x!> number of grid points in x-direction
    integer(kind=ik)  :: ny                                 !x!> number of grid points in y-direction
    integer(kind=ik)  :: qy                                 !x!> number of grid points in y-direction
    integer(kind=ik)  :: nz                                 !x!> number of grid points in z-direction
    integer(kind=ik)  :: qz                                 !x!> number of grid points in z-direction

    integer(kind=ik)  :: q                                  !x!> size of ppm arrays
    integer(kind=ik)  :: max_dim                            !x!> not the slightest clue
    integer(kind=ik)  :: q_nqx                              !x!> size of some ppm arrays
    integer(kind=ik)  :: q_nqy                              !x!> size of some ppm arrays
    integer(kind=ik)  :: q_nqz                              !x!> size of some ppm arrays
    integer(kind=ik)  :: q_xoff                             !x!> size of some ppm arrays
    integer(kind=ik)  :: q_yoff                             !x!> size of some ppm arrays
    integer(kind=ik)  :: q_zoff                             !x!> size of some ppm arrays

    integer(kind=ik)  :: igeomx                             !x!> geometry for x-direction
    integer(kind=ik)  :: igeomy                             !x!> geometry for y-direction
    integer(kind=ik)  :: igeomz                             !x!> geometry for z-direction
    integer(kind=ik)  :: bndmnx                             !x!> inner boundary condition in x-direction
    integer(kind=ik)  :: bndmxx                             !x!> outer boundary condition in x-direction
    integer(kind=ik)  :: bndmny                             !x!> inner boundary condition in y-direction
    integer(kind=ik)  :: bndmxy                             !x!> outer boundary condition in y-direction
    integer(kind=ik)  :: bndmnz                             !x!> inner boundary condition in z-direction
    integer(kind=ik)  :: bndmxz                             !x!> outer boundary condition in z-direction
    integer(kind=ik)  :: latybc                             !x!> some boundary condition for the transport?
    integer(kind=ik)  :: latzbc                             !x!> some boundary condition for the transport?
    integer(kind=ik)  :: qn                                 !x!> number of ALL nuclear species (incl. YE and dummy nuclei)
    integer(kind=ik)  :: qn_network                         !x!> number of species used in network
    integer(kind=ik)  :: i_grv                              !x!> type of grav. potential, 0: Newtonian, 1: GR (monopol), 2: GRapprox for 2D
    integer(kind=ik)  :: gr_approximation                   !x!> type of GR approximation for grav. pot, 1: new (Marek et al. 2006), 2: wrong order of nu-pressure terms
    integer(kind=ik)  :: enforce_spherical_potential        !x!> enforce spherical gravitational potential
    integer(kind=ik)  :: isym                               !x!> assume equatorial symmetry, if set to 1
    integer(kind=ik)  :: noise                              !x!> perturb initial model, if set to 1
    real(kind=rk)     :: ampl                               !x!> amplitude of random perturbation
    real(kind=rk)     :: spherical_neutron_star             !x!> threshold density to use for a spherical inner core, or 0 to not do a spherical core
    integer(kind=ik)  :: are_nu                             !x!> number of separate areas in x-direction
    real(kind=rk)     :: dt                                 !x!> TODO: This is not set anywhere, but used in burn.F90!
    real(kind=rk)     :: dt_rad                             !x!> TODO: This is not used anywhere
    real(kind=rk)     :: dtmaxx                             !x!> maximum timestep with subcycling, unit: "s"
    character(200)    :: progenitor_name                    !x!> file name of initial progenitor model
    character(200)    :: setup_mode                         !x!> setup mode for initial progenitor model
    integer(kind=ik)  :: nxa                                !x!> number of last zone of 1. area in x-dir
    integer(kind=ik)  :: nxb                                !x!> number of last zone of 2. area in x-dir
    integer(kind=ik)  :: nxc                                !x!> number of last zone of 3. area in x-dir
    integer(kind=ik)  :: ioya                               !x!> grid resolution in y-direction in first zone
    integer(kind=ik)  :: ioyb                               !x!> grid resolution in y-direction in second zone
    integer(kind=ik)  :: ioyc                               !x!> grid resolution in y-direction in third zone
    integer(kind=ik)  :: ioza                               !x!> grid resolution in z-direction in first zone
    integer(kind=ik)  :: iozb                               !x!> grid resolution in z-direction in second zone
    integer(kind=ik)  :: iozc                               !x!> grid resolution in z-direction in third zone
    integer(kind=ik)  :: ndtmax                             !x!> TODO: Seems not to be used anywhere with consequence

    real(kind=rk)     :: gridlx_t                           !x!> transport grid radius for a new run, use hydro grid if <= 0, unit: "cm"
    real(kind=rk)     :: pmass                              !x!> mass of the central neutron star, unit: "g"
    integer(kind=ik)  :: geoen                              !x!> type of energy grid for new run, 0: log, 1: linlog
    character(7)      :: rfile                              !x!> file with the initial radiat. field (or 'no')
    integer(kind=ik)  :: p_ntr                              !x!> use neutrino transport (0=OFF/1=ON)
    integer(kind=ik)  :: p_nbk                              !x!> use neutrino backreactions (0=off/1=on)
    integer(kind=ik)  :: i_grtr                             !x!> gen. relativistic transp. (0=off/1=on/-1=CFC)
    integer(kind=ik)  :: irst_ra                            !x!> restart mode for transport: (2=RESET TIME[step])
    integer(kind=ik)  :: isnnv(3)
    integer(kind=ik)  :: isma                               !x!> number of neutrino species
    integer(kind=ik)  :: isnn                               !x!> matrix block size factor
    integer(kind=ik)  :: iemax                              !x!> number of energy bins

    character(200)    :: index_file                         !x!> .inx filename
    character(200)    :: evolution_file                     !x!> .evo filename
    character(200)    :: neutrino_file                      !x!> .ntr filename
    character(200)    :: gw_file                            !x!> .gw3d filename
    character(200)    :: domain_file                        !x!> .domain_decomposition filename
    character(200)    :: timestep_file                      !x!> .tim filename
    character(200)    :: energy_file                        !x!> .erg filename

    integer(kind=ik)  :: imaxp                              !x!> number of grid points for transport
    integer(kind=ik)  :: cmin                               !x!> number of shells penetrating IB (-cmin)
    integer(kind=ik)  :: jvisc                              !x!> J-diffusity (0=off, 1=old, 2=new, 3=only in shock)
    integer(kind=ik)  :: nymom                              !x!> number of angular rays in theta direction (ME)
    integer(kind=ik)  :: nzmom                              !x!> number of angular rays in phi direction (ME)
    integer(kind=ik)  :: nytra                              !x!> number of anuglar rays in theta direction (BTE)
    integer(kind=ik)  :: nztra                              !x!> number of angular rays in phi direction (BTE)
    integer(kind=ik)  :: nystrt                             !x!> first theta index
    logical           :: use_spherical_eddington_factor     !x!> use spherical Eddington factor
    logical           :: use_openmp_rays                    !x!> use OpenMP parallelization over angular rays
    logical           :: use_openmp_energy_bins             !x!> use OpenMP parallelization over energy bins
    logical           :: use_openmp_matrix                  !x!> use OpenMP parallelized matrix solver

    real(kind=rk)     :: dtmin_rt                           !x!> minimum size of transport timestep, unit: "s"
    real(kind=rk)     :: sigma1                             !x!> sigma1
    real(kind=rk)     :: sigma2                             !x!> sigma2 (mean temporal changes)
    real(kind=rk)     :: sigma3                             !x!> sigma3 (maximum temporal changes)

    real(kind=rk)     :: epsit                              !x!> accuracy for me-bte-iteration
    integer(kind=ik)  :: itanz                              !x!> maximum number of me-bte-iterations
    integer(kind=ik)  :: nsiout                             !x!> skip nsiout outputs of the spec. intensit

    real(kind=rk)     :: wolff_ls_rho                       !x!> density between LS EoS and Wolfflow EoS, unit: "g/cm^3"
    real(kind=rk)     :: tj_wolff_rho                       !x!> density between Wolfflow EoS and ThJ Eos, unit: "g/cm^3"
    real(kind=rk)     :: tj_ls_rho                          !x!> density between LS EoS and ThJ Eos [g/cc], unit: "g/cm^3"
    integer(kind=ik)  :: eos_sw                             !x!> switch lsrolo at first EoS call? (0=N/1=Y)
    integer(kind=ik)  :: restmass_version                   !x!> version of no-restmass scheme
    integer(kind=ik)  :: low_den_nse_eos                    !x!> type of low density NSE EoS
    logical           :: use_network                        !x!> use nuclear burning network (0=no,1=yes)
    logical           :: use_flash_si                       !x!> do flashing of Si
    logical           :: use_flash_o                        !x!> do flashing of O
    logical           :: use_flash_c                        !x!> do flashing of C
    real(kind=rk)     :: tkok                               !x!> threshold temperature for low density NSE, unit: "MeV"
    integer(kind=ik)  :: ihvers                             !x!> version number of hydro output
    integer(kind=ik)  :: irvers                             !x!> version number of transport output

    integer(kind=ik)  :: produc                             !x!> use security file
    integer(kind=ik)  :: laghyd                             !x!> use eulerian (0) or lagrangian (1) hydro grid
    integer(kind=ik)  :: intout                             !x!> interval for summary output in .tim-file
    integer(kind=ik)  :: wcrem                              !x!> interpret as cpu (0) or wall clock (1) time
    logical           :: use_multid_collapse                !x!> multid ralag
    integer(kind=ik)  :: ieul                               !x!> 1: eulerian grid used for neutrino transport, 0: lagrangian
    real(kind=rk)     :: epsitacc                           !x!> emach/config%epsit

    logical           :: kjtrates                           !x!> use rates from KJT
    logical           :: kjttable                           !x!> use tabulated KJT rates
    logical           :: nunu                               !x!> use nu-nu rates
    logical           :: lmsrates                           !x!> use LMS rates
    logical           :: lms_old                            !x!> use old implementation of LMS rates
    logical           :: lms_heavy                          !x!> use nucleus information from LMS table
    logical           :: lms_spectra                        !x!> use LMS rates with spectra
    logical           :: use_low_den_electron_captures      !x!> use electron captures in low density eos (at the moment this is only done in NSE)
    logical           :: qfit                               !x!> use qfitting value together with spectra
    logical           :: itoh_correlation                   !x!> use ion-ion correlations according to ITOH
    logical           :: itoh_gamma                         !x!> if itoh, use assemble average
    logical           :: c_los                              !x!> take polarisation at scattering
    logical           :: isnrates                           !x!> consider ISN rates
    logical           :: isn_heavy                          !x!> use nucleus information together with ISN
    logical           :: brems                              !x!> use bremstrahlung
    logical           :: nickelrates                        !x!> use nickel rates
#ifdef WRITE_BINARY_OUTPUT
!PERL my @type_defs = @_lines[0 .. $_line_number];
#endif
  end type conf_t

  type(conf_t)  :: config

  contains

    subroutine initialize_config_values

      use precision
      implicit none

      config%isnnv(1)=1
      config%isnnv(2)=2
      config%isnnv(3)=2

      config%isnn = config%isnnv(config%isma)

    end subroutine initialize_config_values
#ifdef WRITE_BINARY_OUTPUT
    subroutine write_config(fd)
      use dataformat
      type(datafile), intent(in) :: fd

!PERL my $i = 1;
!PERL foreach my $line (@type_defs) {
!PERL   if ($line =~ /^.*::\s*(\S*)\s+!x!> (.*)$/) {
!PERL     my $name = $1;
!PERL     my $desc = $2;
!PERL     if ($desc =~ /^(.*), unit: "(.*)"\s*$/) {
!PERL       $desc = $1;
!PERL       my $unit = $2;
      call writeq(fd, &
        "@[[$name]]", &
        config%@[[$name]], &
        "@[[$desc]]", &
        "@[[$unit]]")
!PERL     } else {
      call writeq(fd, &
        "@[[$name]]", &
        config%@[[$name]], &
        "@[[$desc]]")
!PERL     }
!PERL   };
!PERL }
!PERL-END
    end subroutine write_config
#endif

end module configure

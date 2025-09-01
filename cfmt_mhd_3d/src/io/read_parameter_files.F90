module parameters
  use precision

  implicit none

  private
  public :: read_parameter_files


 interface read_line

  module procedure  read_line_line,    &
                    read_line_integer, &
                    read_line_real,    &
                    read_line_string

  end interface


  contains

 subroutine read_line_line(funit, print_it)

  use abort
  use precision
  use print_stdout_mod

  implicit none

  integer(kind=ik), intent(in) :: funit, print_it
  character(len=72)                 :: text

  integer :: stat
  character(len=100) :: msg

!  read (funit,*,iostat=stat,iomsg=msg) text
  read (funit,*,iostat=stat) text

  if (stat .ne. 0) then
    raise_abort("Error reading in parameter file, expected comment text: " // trim(msg))
  endif

  call printit_taskX(print_it,text)

 end subroutine read_line_line


 subroutine read_line_integer(funit, print_it, label, val)

  use precision
  use abort
  use print_stdout_mod
  implicit none

  integer(kind=ik), intent(in)  :: funit, print_it
  character*(*), intent(in)       :: label

  integer(kind=ik), intent(inout) :: val

  character(len=72)      :: label_rd

  character(len=44)     :: txt
  character(len=5)      :: txtxt

  integer :: stat
  character(len=100) :: msg

  data txtxt /'...  '/

!  read (funit,*,iostat=stat,iomsg=msg) txt, label_rd, val
  read (funit,*,iostat=stat) txt, label_rd, val

  if (stat .ne. 0) then
    raise_abort("Error reading in parameter '" // trim(label) //"': " // trim(msg))
  endif

  call printit_taskX(print_it,txt, txtxt,label,val)

  if (label .ne. label_rd) then
    close (1)
    raise_abort("Expected parameter '" // trim(label) //"', got '" // trim(label_rd))
  endif

 end subroutine read_line_integer

 subroutine read_line_real(funit, print_it, label, val)

  use precision
  use abort
  use print_stdout_mod
  implicit none

  integer(kind=ik), intent(in)  :: funit, print_it
  character*(*), intent(in)       :: label

  real(kind=rk), intent(inout)    :: val

  character(len=72)                   :: label_rd

  character(len=44)                  :: txt
  character(len=5)                   :: txtxt

  integer :: stat
  character(len=100) :: msg

  data txtxt /'...  '/

!  read (funit,*,iostat=stat,iomsg=msg) txt, label_rd, val
  read (funit,*,iostat=stat) txt, label_rd, val

  if (stat .ne. 0) then
    raise_abort("Error reading in parameter '" // trim(label) //"': " // trim(msg))
  endif

  call printit_taskX(print_it,txt, txtxt,label,val)

  if (label .ne. label_rd) then
    close (1)
    raise_abort("Expected parameter '" // trim(label) //"', got '" // trim(label_rd))
  endif

 end subroutine read_line_real


 subroutine read_line_string(funit, print_it, label, val)

  use precision
  use abort
  use print_stdout_mod
  implicit none

  integer(kind=ik), intent(in)  :: funit, print_it
  character*(*), intent(in)       :: label

  character*(*), intent(inout)    :: val

  character(len=72)                   :: label_rd

  character(len=44)                  :: txt
  character(len=5)                   :: txtxt

  integer :: stat
  character(len=100) :: msg

  data txtxt /'...  '/

!  read (funit,*,iostat=stat,iomsg=msg) txt, label_rd, val
  read (funit,*,iostat=stat) txt, label_rd, val

  if (stat .ne. 0) then
    raise_abort("Error reading in parameter '" // trim(label) //"': " // trim(msg))
  endif

  call printit_taskX(print_it,txt, txtxt,label,val)

  if (label .ne. label_rd) then
    close (1)
    raise_abort("Expected parameter '" // trim(label) //"', got '" // trim(label_rd) // "'")
  endif

 end subroutine read_line_string


!> \par read input parameters for PPM from file ppm.par
!>
!> \author W Keil, M. Rampp (MPA), A. Marek (MPA)
!>
!> \verbatim
!>   SVN - Information
!>   $Revision:$
!>   $Date:$
!> \endverbatim
!>
subroutine read_parameter_files(read_mode)

  use precision

!  use intgrs_hy
  use charac
!  use totare_hy
!  use totgrq_hy
!  use arecon_hy
!  use revsho_hy   ! forcheck
#ifndef DEBUG_TIMINGS
  use cputim
#endif

!  use nutrio_hy
  use param_rt
  use phycon

  use gfloat_hy

  use mo_mpi


  use abort

#ifdef WRITE_BINARY_OUTPUT
  use output_hydro, only : set_filenames
  use restart, only : read_restart_filename! , use_temp_restart
#endif

  use print_stdout_mod, only : printit_taskX, printit_alltasks

  use configure
!  use hydro_areas_mod

  use nutrio_hy
  use machcons
  use totgrq_hy, only : sumdegr
  implicit none
! LOCAL variables that are not in modules

  character*(*), intent(in) :: read_mode

  integer(kind=ik) :: n

  character(len=8)      :: label1, label

  character(len=72)     :: text
  character(len=44)     :: txt
  character(len=5)      :: txtxt
  data txtxt /'...  '/

  integer(kind=ik) :: i, j

  real(kind=rk)    :: dtmaxx_over

  integer(kind=ik) :: ierr, qn_in


  integer(kind=ik) :: nsp, ner, nrad, nradi, nradj, nradk, spherical_eddington_factor
  integer(kind=ik) :: omp_energy, omp_rays, omp_matrix
  integer(kind=ik) :: alternate
  real(kind=rk)    :: dummyreal
  integer(kind=ik) :: dummyInteger

  integer(kind=ik)  :: multid

  integer(kind=ik)  :: print_it

  integer(kind=ik)  :: kjtrates, kjttable, nunurates, lmsrates, lmsspectra, &
                       qfit, lmsheavy, itoh, itohgamma, clos, isnrates,     &
                       isnheavy, brems, nickelrates, lmsold

  config%use_temp_restart = .false.

!      integer(kind=ik) :: nstr,nstrt,nxt
!      real(kind=rk) :: xmfrac,rstr
!      common /movgrid/  nstr,nstrt,nxt,xmfrac,rstr

!----------------------------------------------------------------------
!     open parameter file and read content:
!----------------------------------------------------------------------
  ! if read_mode = grid_init then we only have to read
  ! a part of the file which specifies the hydro grid dimensions
  ! and we do not want to print everything in this mode
  ! after we read these values and we were able to do the
  ! domain decompostion, we will have to read this file again
  ! in order to get all the other values on each MPI task
  !
  ! the printing is switched of by priting to myrpoc = -1
  if (trim(read_mode) .eq. "grid_init") then
     print_it = -1
  else if (trim(read_mode) .eq. "model_init") then
     print_it =  0
  else
     raise_abort("input(): unknown read_mode")
  endif

  open (1, file='config', form='formatted', status='OLD')

   call printit_taskX(print_it," ")

  do n = 1, 3
      call read_line(1,print_it)
  enddo
  call read_line(1,print_it,"calculation",config%calculation)
  call read_line(1,print_it)
  call read_line(1,print_it,"progenitor",config%progenitor_name)
  call read_line(1,print_it,"setup_mode",config%setup_mode)
    select case(config%setup_mode)
    case("rho-t-ye")
    case("rho-p-ye")
    case("rho-s-ye")
    case default
       raise_abort("Invalid value for setup_mode: '" // trim(config%setup_mode)  // "'")
    end select

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"modelname",config%basenm)
  call read_line(1,print_it,"cpumax",config%cpumax)
  call read_line(1,print_it,"wcrem",config%wcrem)
  call read_line(1,print_it,"nend",config%nend)

#ifdef BENCHMARK
  config%nend = 2
#endif

   call read_line(1,print_it,"tmax",config%tmax)
   call read_line(1,print_it,"irstrt",config%irstrt)
   call read_line(1,print_it)

#ifndef NEW_OUTPUT_LABELS
  ! use the old labeling of output files aa,ab..zz,AA..ZZ
  ! you are restricted to 1352 files
   call read_line(1,print_it,"suffix",config%suffix(1:2))
#else  /* NEW_OUTPUT_LABELS */
  ! use a new labeling 00000000 ,00000001 ...
   call read_line(1,print_it,"suffix",config%suffix_nr(1:8))
#endif /* NEW_OUTPUT_LABELS */

#ifdef WRITE_BINARY_OUTPUT
   if ( config%irstrt .ne. 0 ) then
      call read_restart_filename
   else
      call set_filenames(0,1)
   endif
#endif

   call read_line(1,print_it,"nrstrt",config%nrstrt)
   call read_line(1,print_it,"trstrt",config%trstrt)
   call read_line(1,print_it,"nout",config%nout)
   call read_line(1,print_it,"tout",config%tout)
   call read_line(1,print_it,"max_central_density",config%max_central_density)
   call read_line(1,print_it,"stop_bounce",config%stop_bounce)
   call read_line(1,print_it,"stop_bad_resolution",config%stop_bad_resolution)
   call read_line(1,print_it,"resolution_grace_period",config%bad_resolution_grace_period)
   call read_line(1,print_it)
   call read_line(1,print_it,"collapse_time_resolve",config%collapse_time_resolve)
   call read_line(1,print_it,"itstp",config%itstp)
   call read_line(1,print_it,"rst_mode",config%rst_mode)
  !

  if (trim(config%rst_mode) .eq. "energy") then
     config%use_temp_restart = .false.
  else if (trim(config%rst_mode) .eq. "remap") then
     config%use_temp_restart = .false.
  else if (trim(config%rst_mode) .eq. "temp") then
     config%use_temp_restart = .true.
  else
     raise_abort("input(): Invalid rst_mode set")
  endif

  call read_line(1,print_it,"produc",config%produc)
  call read_line(1,print_it,"redirect",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_stdout_redirect = .false.
  else
     config%use_stdout_redirect = .true.
  endif

  call read_line(1,print_it,"deactivate_output",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_deactivate_output = .false.
  else
     config%use_deactivate_output = .true.
  endif

  call read_line(1,print_it,"deactivate_stdout",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_deactivate_stdout = .false.
  else
     config%use_deactivate_stdout = .true.
  endif

  call read_line(1,print_it,"deactivate_control_files",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_deactivate_controlFiles = .false.
  else
     config%use_deactivate_controlFiles = .true.
  endif

  call read_line(1,print_it,"print_mem_alloc",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_print_memAlloc = .false.
  else
     config%use_print_memAlloc = .true.
  endif

  call read_line(1,print_it,"print_mpi_domain",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_print_mpiDomain = .false.
  else
     config%use_print_mpiDomain = .true.
  endif
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"nsdim",config%nsdim)

  call read_line(1,print_it,"nx",config%qx)

  if (config%qx .lt. 4) then
     raise_abort("input(): qx lower than 4")
  endif

  call read_line(1,print_it,"ny",config%qy)
  if (config%nsdim .gt. 1  .and.  (config%qy .lt. 4)) then
     raise_abort("input(): qy lower than 4")
  endif

  call read_line(1,print_it,"nz",config%qz)

  if (config%nsdim .gt. 2  .and.  (config%qz .lt. 4)) then
     raise_abort("input(): qz lower than 4")
  endif

  if (.not.use_mpi)  config%q = max(config%qx, config%qy, config%qz) + 20

  config%q_nqx   = config%qx+20
  config%q_nqy   = config%qy+20
  config%q_nqz   = config%qz+20
  config%max_dim = max(config%qx+20, config%qy+20, config%qz+20)
  config%nx      = config%qx
  config%ny      = config%qy
  config%nz      = config%qz

  call read_line(1,print_it,"rib",config%rib)
  call read_line(1,print_it,"xmfrac",config%xmfrac)

  config%xmfrac=config%xmfrac*pc_ms

  call read_line(1,print_it,"gridlx",config%gridlx)
  call read_line(1,print_it,"gridly",config%gridly)
  call read_line(1,print_it,"gridlz",config%gridlz)
  call read_line(1,print_it,"qn",qn_in)

  config%qn         = qn_in
  call read_line(1,print_it,"qn_network",qn_in)
  config%qn_network = qn_in

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"isym",config%isym)
  call read_line(1,print_it,"laghyd",config%laghyd)
  !
  !
  !-------  igeom = 0   ====>  planar      geometry
  !         igeom = 1   ====>  cylindrical geometry (radial)
  !         igeom = 2   ====>  spherical   geometry (radial)
  !         igeom = 3   ====>  cylindrical geometry (angular)
  !         igeom = 4   ====>  spherical   geometry (angular - theta)
  !         igeom = 5   ====>  spherical   geometry (angular - phi)
  !
  call read_line(1,print_it,"igeomx",config%igeomx)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  if (config%igeomx .gt. 2) then
     raise_abort("input(): igeomx")
  endif

  call read_line(1,print_it,"igeomy",config%igeomy)

  if ( (config%igeomx .eq. 0 .and. config%igeomy .ne. 0)   .or. &
       (config%igeomx .eq. 1 .and. config%igeomy .eq. 1)   .or. &
       (config%igeomx .eq. 1 .and. config%igeomy .gt. 3)   .or. &
       (config%igeomx .eq. 2 .and. config%igeomy .lt. 4)        )  then
     raise_abort("input(): igeomy")
  endif

  call read_line(1,print_it,"igeomz",config%igeomz)

  if ( (config%igeomx .eq. 0                   .and. config%igeomz .ne. 0)  .or. &
       (config%igeomx .eq. 1 .and. config%igeomy .eq. 0 .and. config%igeomz .ne. 3)  .or. &
       (config%igeomx .eq. 1 .and. config%igeomy .eq. 3 .and. config%igeomz .ne. 0)  .or. &
       (config%igeomx .eq. 2 .and. config%igeomy .eq. 4 .and. config%igeomz .ne. 5)  .or. &
       (config%igeomx .eq. 2 .and. config%igeomy .eq. 5 .and. config%igeomz .ne. 4)       )&
       then
     raise_abort("input(): igeomz")
  endif

  !
  !-------  bndm.. = 1  ====>  reflecting boundary
  !         bndm.. = 2  ====>  flow out   boundary
  !         bndm.. = 3  ====>  flow in    boundary
  !         bndm.. = 4  ====>  periodic   boundary
  !         bndm.. = 5  ====>  any other  boundary
  !

  call read_line(1,print_it,"bndmnx",config%bndmnx)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"bndmxx",config%bndmxx)
  call read_line(1,print_it,"bndmny",config%bndmny)
  call read_line(1,print_it,"bndmxy",config%bndmxy)
  call read_line(1,print_it,"bndmnz",config%bndmnz)
  call read_line(1,print_it,"bndmxz",config%bndmxz)

  if ( (config%bndmnx .eq. 4  .and.  config%bndmxx .ne. 4)  .or. &
       (config%bndmny .eq. 4  .and.  config%bndmxy .ne. 4)  .or. &
       (config%bndmnz .eq. 4  .and.  config%bndmxz .ne. 4)       ) then
     raise_abort("input(): boundary")
  endif
  !


  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"spherical_neutron_star",config%spherical_neutron_star)
  call read_line(1,print_it)
  call read_line(1,print_it,"are_nu",config%are_nu)
  call read_line(1,print_it,"nxa",config%nxa)

  if (config%nxa .lt. 4) then
     raise_abort("input(): nxa")
  endif

  call read_line(1,print_it,"ioya",config%ioya)

  if (config%ioya .lt. 1) then
     raise_abort("input(): ioya")
  endif

  call read_line(1,print_it,"ioza",config%ioza)

  if (config%ioza .lt. 1) then
     raise_abort("input(): ioza")
  endif

  call read_line(1,print_it,"nxb",config%nxb)

  if (config%nxb .lt. 0) then
     raise_abort("input(): nxb")
  endif

  call read_line(1,print_it,"ioyb",config%ioyb)

  if (config%ioyb .lt. 1) then
     raise_abort("input(): ioyb")
  endif

  call read_line(1,print_it,"iozb",config%iozb)

  if (config%iozb .lt. 1) then
     raise_abort("input(): iozb")
  endif

  call read_line(1,print_it,"nxc",config%nxc)

  if (config%nxc .lt. 0) then
     raise_abort("input(): nxc")
  endif

  call read_line(1,print_it,"ioyc",config%ioyc)

  if (config%ioyc .lt. 1) then
     raise_abort("input(): ioyc")
  endif

  call read_line(1,print_it,"iozc",config%iozc)

  if (config%iozc .lt. 1) then
     raise_abort("input(): iozc")
  endif

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"cfl",config%cfl)
  call read_line(1,print_it,"nriem",config%nriem)
  call read_line(1,print_it,"cvisc",config%cvisc)
  call read_line(1,print_it,"small",config%small)
  call read_line(1,print_it,"smlrho",config%smlrho)
  call read_line(1,print_it,"smallp",config%smallp)
  call read_line(1,print_it,"smalle",config%smalle)
  call read_line(1,print_it,"smallu",config%smallu)
  call read_line(1,print_it,"smallx",config%smallx)
  call read_line(1,print_it,"igodu",config%igodu)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"dtini",config%dtini)
  call read_line(1,print_it,"dtmin",config%dtmin)
  call read_line(1,print_it,"dtmax",config%dtmax)
  call read_line(1,print_it,"dtmaxx",config%dtmaxx)

  ! Catch overflow error
  dtmaxx_over=min(min((config%dtmaxx/config%dtmin),2000._rk), real(HUGE(1),kind=rk))
  config%ndtmax=int(dtmaxx_over)

  call read_line(1,print_it,"dtmn_rt",config%dtmin_rt)

  config%dtmin_rt=max(config%dtmin_rt,config%dtmin)

  call read_line(1,print_it,"intout",config%intout)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"noise",config%noise)
  call read_line(1,print_it,"ampl",config%ampl)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"i_grv",config%i_grv)
  call read_line(1,print_it)
  call read_line(1,print_it,"gr_approx",config%gr_approximation)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"enforce_spherical_potential",config%enforce_spherical_potential)
  call read_line(1,print_it,"i_grtr",config%i_grtr)


  if (trim(read_mode) .ne. "grid_init") then
     if (config%i_grtr .ne. 0 .and. config%i_grv .eq. 0 ) then
        raise_abort("inpmod(): i_grv <-> i_grtr")
     endif
  endif


  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"p_ntr",config%p_ntr)
  call read_line(1,print_it,"p_nbk",config%p_nbk)



  call read_line(1,print_it,"irst_ra",config%irst_ra)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)


  if (config%p_ntr .ne. 0 .or. config%irst_ra .eq. 0) then
     config%dt_rad = config%dtini
  endif

  call read_line(1,print_it,"isma",config%isma)
  call read_line(1,print_it,"ner",config%iemax)
  call read_line(1,print_it,"geoen",config%geoen)
  call read_line(1,print_it,"gridlx_t",config%gridlx_t)
  call read_line(1,print_it,"nrad",config%imaxp)
  call read_line(1,print_it,"cmin",config%cmin)
  call read_line(1,print_it,"jvisc",config%jvisc)

  config%cmin = -1*config%cmin

  call read_line(1,print_it,"multid",multid)
  call read_line(1,print_it)

  if (multid .eq. 1) then
     config%use_multid_collapse = .true.
  else
     config%use_multid_collapse = .false.
  endif

  call read_line(1,print_it,"nymom",config%nymom)
  call read_line(1,print_it,"nzmom",config%nzmom)
  call read_line(1,print_it,"spher_edd",spherical_eddington_factor)

!                nymom   nytra   nystrt  nzmom   nztra
!
!1D              1       1       1       x       1
!
!2D-1D-Edd       16      0       0       x       1
!
!2D-2D-Edd       16      16      1       x       1
!
!3D-1D-Edd       --------------------------------
!
!3D-3D-Edd       16      16      1       x       32


  if (spherical_eddington_factor .eq. 1) then
     ! check whether a 1D model is set up
     config%use_spherical_eddington_factor = .true.
     if (config%nymom .eq. 1 .and. config%nzmom .eq. 1) then
        ! in 1D nystrt = 1
        config%nystrt = 1
        config%nytra  = config%nymom
        config%nztra  = config%nzmom
     else if (config%nymom .gt. 1 .and. config%nzmom .eq. 1) then
        ! in 2D nystrt = 0
        config%nystrt = 0
        config%nytra  = 0
        config%nztra  = config%nzmom
     else if (config%nymom .eq. 1 .and. config%nzmom .gt. 1) then
        raise_abort("inpmod(): This 2D setup is not supported !")
     else if (config%nymom .gt. 1 .and. config%nzmom .gt. 1) then
        ! in 3D nystrt = 0
        config%nystrt = 0
        config%nytra  = 0
        config%nztra  = config%nzmom
     else
        raise_abort("inpmod(): This setup with spherical Eddington factor is not supported !")
     endif
  else if (spherical_eddington_factor .eq. 0) then
     config%use_spherical_eddington_factor = .false.
     if (config%nymom .eq. 1 .and. config%nzmom .eq. 1) then
        ! in 1D nystrt = 1
        config%nystrt = 1
        config%nytra  = config%nymom
        config%nztra  = config%nzmom
     else if (config%nymom .gt. 1 .and. config%nzmom .eq. 1) then
        ! in 2D nystrt = 1
        config%nystrt = 1
        config%nytra  = config%nymom
        config%nztra  = config%nzmom
     else if (config%nymom .eq. 1 .and. config%nzmom .gt. 1) then
        raise_abort("inpmod(): This 2D setup is not supported !")
     else if (config%nymom .gt. 1 .and. config%nzmom .gt. 1) then
        ! in 3D nystrt = 1
        config%nystrt = 1
        config%nytra  = config%nymom
        config%nztra  = config%nzmom
     else
        raise_abort("inpmod(): This setup is not supported !")
     endif

  else
     raise_abort("inpmod(): Wrong value for spherical Eddington factor !")
  endif

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"sigma1",config%sigma1)
  call read_line(1,print_it,"sigma2",config%sigma2)

  config%sigma1=config%sigma2/1.75_rk

  call read_line(1,print_it,"sigma3",config%sigma3)
  call read_line(1,print_it,"itanz",config%itanz)
  call read_line(1,print_it,"epsit",config%epsit)

  config%epsitacc=emach/config%epsit

  call read_line(1,print_it,"nsiout",config%nsiout)

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"kjt",kjtrates)

  if (kjtrates .eq. 0) then
     config%kjtrates = .false.
  else if (kjtrates .eq. 1) then
     config%kjtrates = .true.
  else
     raise_abort("read_parameters(): illegal option for kjt")
  endif

  call read_line(1,print_it,"kjttable",kjttable)

  if (kjttable .eq. 0) then
     config%kjttable = .false.
  else if (kjttable .eq. 1) then
     config%kjttable = .true.
  else
     raise_abort("read_parameters(): illegal option for kjttable")
  endif


  call read_line(1,print_it,"nunu",nunurates)

  if (nunurates .eq. 0) then
     config%nunu = .false.
  else if (nunurates .eq. 1) then
     config%nunu = .true.
  else
     raise_abort("read_parameters(): illegal option for nunurates")
  endif

 call read_line(1,print_it,"lmsrates",lmsrates)

  if (lmsrates .eq. 0) then
     config%lmsrates = .false.
  else if (lmsrates .eq. 1) then
     config%lmsrates = .true.
  else
     raise_abort("read_parameters(): illegal option for lmsrates")
  endif

 call read_line(1,print_it,"lmsspectra",lmsspectra)

  if (lmsspectra .eq. 0) then
     config%lms_spectra = .false.
  else if (lmsspectra .eq. 1) then
     config%lms_spectra = .true.
  else
     raise_abort("read_parameters(): illegal option for lmsspectra")
  endif

 call read_line(1,print_it,"qfit",qfit)

  if (qfit .eq. 0) then
     config%qfit = .false.
  else if (qfit .eq. 1) then
     config%qfit = .true.
  else
     raise_abort("read_parameters(): illegal option for qfit")
  endif


 call read_line(1,print_it,"lmsold",lmsold)

  if (lmsold .eq. 0) then
     config%lms_old = .false.
  else if (lmsold .eq. 1) then
     config%lms_old = .true.
  else
     raise_abort("read_parameters(): illegal option for lmsold")
  endif

 call read_line(1,print_it,"lmsheavy",lmsheavy)

  if (lmsheavy .eq. 0) then
     config%lms_heavy = .false.
  else if (lmsheavy .eq. 1) then
     config%lms_heavy = .true.
  else
     raise_abort("read_parameters(): illegal option for lmsheavy")
  endif

 call read_line(1,print_it,"ec_low_den",dummyInteger)

  if (dummyInteger .eq. 0) then
     config%use_low_den_electron_captures = .false.
  else if (dummyInteger .eq. 1) then
     config%use_low_den_electron_captures = .true.
  else
     raise_abort("read_parameters(): illegal option for lmsheavy")
  endif

 call read_line(1,print_it)

 call read_line(1,print_it,"itoh",itoh)

  if (itoh .eq. 0) then
     config%itoh_correlation = .false.
  else if (itoh .eq. 1) then
     config%itoh_correlation = .true.
  else
     raise_abort("read_parameters(): illegal option for itoh")
  endif

 call read_line(1,print_it,"itohgamma",itohgamma)

  if (itohgamma .eq. 0) then
     config%itoh_gamma = .false.
  else if (itohgamma .eq. 1) then
     config%itoh_gamma = .true.
  else
     raise_abort("read_parameters(): illegal option for itohgama")
  endif

 call read_line(1,print_it,"clos",clos)

  if (clos .eq. 0) then
     config%c_los = .false.
  else if (clos .eq. 1) then
     config%c_los = .true.
  else
     raise_abort("read_parameters(): illegal option for clos")
  endif

 call read_line(1,print_it,"isnrates",isnrates)

  if (isnrates .eq. 0) then
     config%isnrates = .false.
  else if (isnrates .eq. 1) then
     config%isnrates = .true.
  else
     raise_abort("read_parameters(): illegal option for isnrates")
  endif

 call read_line(1,print_it,"isnheavy",isnheavy)

  if (isnheavy .eq. 0) then
     config%isn_heavy = .false.
  else if (isnheavy .eq. 1) then
     config%isn_heavy = .true.
  else
     raise_abort("read_parameters(): illegal option for isn_heavy")
  endif

 call read_line(1,print_it,"brems",brems)

  if (brems .eq. 0) then
     config%brems = .false.
  else if (brems .eq. 1) then
     config%brems = .true.
  else
     raise_abort("read_parameters(): illegal option for brems")
  endif


 call read_line(1,print_it,"nickelrates",nickelrates)

  if (nickelrates .eq. 0) then
     config%nickelrates = .false.
  else if (nickelrates .eq. 1) then
     config%nickelrates = .true.
  else
     raise_abort("read_parameters(): illegal option for nickelrates")
  endif


  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"tj_ls_rho",config%tj_ls_rho)
!  call read_line(1,print_it,"wolff_ls_rho",config%wolff_ls_rho)
!  call read_line(1,print_it,"tj_wolff_rho",config%tj_wolff_rho)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it,"eos_sw",config%eos_sw)
  call read_line(1,print_it,"restmass",config%restmass_version)
  call read_line(1,print_it,"lowdennse",config%low_den_nse_eos)
  call read_line(1,print_it,"use_network",dummyInteger)
  if (dummyInteger .eq. 1) then
     config%use_network = .true.
  else
     config%use_network = .false.
  endif
  call read_line(1,print_it)
  call read_line(1,print_it,"flash_c",dummyInteger)
  if (dummyInteger .eq. 1) then
     config%use_flash_c = .true.
  else
     config%use_flash_c = .false.
  endif
  call read_line(1,print_it,"flash_o",dummyInteger)
  if (dummyInteger .eq. 1) then
     config%use_flash_o = .true.
  else
     config%use_flash_o = .false.
  endif
  call read_line(1,print_it,"flash_si",dummyInteger)
  if (dummyInteger .eq. 1) then
     config%use_flash_si = .true.
  else
     config%use_flash_si = .false.
  endif


  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"omp_energy",omp_energy)



#ifndef CFC_TRANSPORT
  if (omp_energy .eq. 1) then
     config%use_openmp_energy_bins =.true.
  endif
#else

#ifdef OPEN_MP_1D
  config%use_openmp_energy_bins =.true.

#endif

#endif /* CFC_TRANSPORT */


  call read_line(1,print_it,"omp_rays",omp_rays)

#ifndef CFC_TRANSPORT
  if (omp_rays .eq. 1) then
     config%use_openmp_rays =.true.

  endif
#else
#if defined(OPEN_MP_2D) || defined(OPEN_MP_3D)
  config%use_openmp_rays =.true.

#endif
#endif /* CFC_TRANSPORT */


  if ((config%use_openmp_rays) .and. (config%use_openmp_energy_bins)) then
     raise_abort("At the moment you can switch on only one OpenMP parallelisation mode")
  endif

!  call set_omp_parallelisation


  call read_line(1,print_it,"omp_matrix",omp_matrix)


  if (omp_matrix .eq. 1) then
     config%use_openmp_matrix = .true.
  else
     config%use_openmp_matrix = .false.
  endif

  call read_line(1,print_it)
  call read_line(1,print_it)
  call read_line(1,print_it)

  call read_line(1,print_it,"pmass",config%pmass)
  call read_line(1,print_it,"rfile",config%rfile)
  call read_line(1,print_it,"ihvers",config%ihvers)
  call read_line(1,print_it,"irvers",config%irvers)
  
  call read_line(1,print_it,"trpp",config%trpp)

  if (trim(read_mode) .eq. "grid_init") then
     ! we only had to read until here in this mode
     ! all necessary values obtained -> the rest of
     ! the parameter file is skipped and we start
     ! immediatly reading the transport parameter file
      call printit_taskX(0," ")
      call printit_taskX(0,"Dimension of hydro grid read:")
      call printit_taskX(0,"config%qx ",config%qx)
      call printit_taskX(0,"config%qy ",config%qy)
      call printit_taskX(0,"config%qz ",config%qz)
      if (.not.use_mpi) call printit_taskX(0,"config%q  ",config%q)
      call printit_taskX(0,"config%qn ",config%qn)
      call printit_taskX(0," ")

      call mpi_barrier(MPI_COMM_WORLD,ierr)

  endif





  close (1)

  config%index_file = './diagnostics/' // trim(config%basenm) // '.inx'
  config%evolution_file = './diagnostics/' // trim(config%basenm) // '.evo'
  config%neutrino_file = './diagnostics/' // trim(config%basenm) // '.ntr'
  config%gw_file = './diagnostics/' // trim(config%basenm) // '.gw3d'
  config%domain_file = './diagnostics/' // trim(config%basenm) // '.domain_decomposition'
 config%timestep_file = './diagnostics/' // trim(config%basenm) // '.tim'
 config%energy_file = './diagnostics/' // trim(config%basenm) // '.erg'


  call printit_taskX(print_it,'input> index file of stored models: ', config%index_file)
  call printit_taskX(print_it,'input> file of the time evolution : ', config%evolution_file)
  call printit_taskX(print_it,'input> file of neutrino evolution : ', config%neutrino_file)
  call printit_taskX(print_it,'input> file of domain decomposition : ', config%domain_file)
  call printit_taskX(print_it,'input> file of domain decomposition : ', config%timestep_file)
  call printit_taskX(print_it,'input> file of domain decomposition : ', config%energy_file)



  config%ieul=0
#ifdef EULGRID
  if (config%p_ntr .ne. 0)  &
       call printit_taskX(print_it,"init_nutra>: EULERIAN GRID used for neutrino transport")
  config%ieul =1
#endif
#ifdef LAGGRID
  if (config%p_ntr .ne. 0) &
       call printit_taskX(print_it,"init_nutra>: LAGRANGIAN GRID used for neutrino transport")
  config%ieul =0
#endif

  select case (config%i_grtr)
  case (-1)
     call printit_taskX(print_it,"inpmod>: GR neutrino transport, GR hydro")
  case (0)
     call printit_taskX(print_it,"inpmod>: NEWTONIAN neutrino transport")
  case (1)
     call printit_taskX(print_it,"inpmod>: GR neutrino transport")
  case (2)
     call printit_taskX(print_it,"inpmod>: APPROXIMATE GR neutrino transport (GR redshift and time dilation)")
  case (3)
     call printit_taskX(print_it,"inpmod>: MINIMAL GR neutrino transport (GR redshift only)")
  case default
     raise_abort("inpmod(): i_grtr: mode not implemented")
  end select

  if (trim(read_mode) .eq. "grid_init") then

     call printit_taskX(0," ")
     call printit_taskX(0,"Dimension of transport grid read:")
     call printit_taskX(0,"isma    ", config%isma)
     call printit_taskX(0,"iemax   ", config%iemax)
     call printit_taskX(0,"imaxp   ", config%imaxp)
     call printit_taskX(0,"nymom   ", config%nymom)
     call printit_taskX(0,"nzmom   ", config%nzmom)
     call printit_taskX(0," ")

  endif


  sumdegr = 0._rk
  sumdenu = 0._rk
  sumdynu = 0._rk
#ifdef FCNC_CALC
  sumdenu_n = 0._rk
  sumdynu_n = 0._rk
#endif

 call mpi_barrier(MPI_COMM_WORLD,ierr)


  return

end subroutine read_parameter_files

end module parameters

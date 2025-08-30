module read_progenitor_model
  use precision
  use function_pointers, only : function_arguments_t

  implicit none
  private
  public :: mass2radius, read_inimod, init, read_star, &
            read_star_models_new, deallocate_progenitor, progenitor_t !,&
!            eos_field, eos_field_arg_t, &
!            EOS_FIELD_PRESSURE, EOS_FIELD_ENTROPY

#ifdef NOMOTO_MODEL
  public :: read_onemg
#endif

  !> Datatype for the progenitor files
  type :: progenitor_t
    integer(kind=ik) :: n_zones            !< number of grid points
    real(kind=rk), allocatable :: rri(:)   !< radial positions at zone interfaces [cm]
    real(kind=rk), allocatable :: rqi(:)   !< radial positions at zone centers [cm]
    real(kind=rk), allocatable :: rhoi(:)  !< density [g/ccm]
    real(kind=rk), allocatable :: yei(:)   !< electron fraction  [1/by]
    real(kind=rk), allocatable :: vqi(:)   !< velocity [cm/s]
    real(kind=rk), allocatable :: omi(:)   !< angular velocity [1/s]
    real(kind=rk), allocatable :: xni(:,:) !< Mass-fractions [1]
    real(kind=rk), allocatable :: tmpi(:)  !< temperature  [K]

    ! optional
    real(kind=rk), allocatable :: sti(:)   !< entropy [kboltzmann]
    real(kind=rk), allocatable :: pri(:)   !< pressure [erg/cm^3]
    real(kind=rk), allocatable :: eni(:)   !< energy density
  end type

  integer, parameter :: EOS_FIELD_PRESSURE = 1
  integer, parameter :: EOS_FIELD_ENTROPY  = EOS_FIELD_PRESSURE + 1
#if 0
  type, extends(function_arguments_t) :: eos_field_arg_t
    real(kind=rk) :: rho
    real(kind=rk), pointer :: xnu(:)
    integer :: field
  end type
#endif

  contains

  function allocate_progenitor(n_zones) result(prog)
    use configure
    implicit none
    integer(kind=ik), intent(in) :: n_zones
    type(progenitor_t) :: prog

    prog%n_zones = n_zones

    allocate(prog%rri(0:n_zones))
    allocate(prog%rqi(1:n_zones))
    allocate(prog%rhoi(1:n_zones))
    allocate(prog%omi(1:n_zones))
    allocate(prog%yei(1:n_zones))
    allocate(prog%vqi(1:n_zones))
    allocate(prog%xni(1:n_zones, 1:config%qn - 1))
    allocate(prog%tmpi(1:n_zones))

  end function

  subroutine deallocate_progenitor(prog)
    implicit none
    type(progenitor_t), intent(inout) :: prog

    deallocate(prog%rri)
    deallocate(prog%rqi)
    deallocate(prog%rhoi)
    deallocate(prog%omi)
    deallocate(prog%yei)
    deallocate(prog%vqi)
    deallocate(prog%xni)
    deallocate(prog%tmpi)

    if (allocated(prog%sti)) then
      deallocate(prog%sti)
    endif

    if (allocated(prog%pri)) then
      deallocate(prog%pri)
    endif

    if (allocated(prog%eni)) then
      deallocate(prog%eni)
    endif

  end subroutine

subroutine mass2radius(rad,den,xm,rstr,imax)

  use phycon

  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik), intent(in) :: imax
  integer(kind=ik)             :: i
  real(kind=rk), intent(out)   :: rstr
                                  ! rstr is the radius which corresponds
                                  ! to the mass-shell xmfrac
  real(kind=rk), intent(in)    :: rad(0:imax),den(imax),xm
                               ! xm is xmfrac the mass coordinate
                               ! arround which the grid is concentrated
                               ! (see ppm.par)

  real(kind=rk)                :: tms(0:imax),dvxt(imax),dvol

  dvxt(1:imax)=(rad(1:imax)**3 - rad(0:imax-1)**3) / 3.0_rk

  tms(0) = 0._rk
  rstr   =-1.e99_rk
  do i = 1, imax
     dvol=dvxt(i) * 4.0_rk*pc_pi
     tms(i) = tms(i-1) + dvol * den(i)
     if (tms(i) .ge. xm) then
        rstr = rad(i)**3 + (xm-tms(i-1))*3._rk/(4._rk*pc_pi*den(i))
        rstr = rstr**(1.0_rk/3.0_rk)
        EXIT
     endif
  enddo

end subroutine mass2radius

#if 0
!> Simple and slow wrapper routine for the EoS to just return a single EoS
!> field, used in eos_inversion
!>
!> \param tem                 Temperature [K]
!> \param additional_args     A eos_field_arg_t structure
!> \param error               Error flag of EoS call, if not present(), abort on error
!>
!> \author    Janina von Groote, Lorenz Huedepohl
!>
function eos_field(tem, additional_args, error) result (val)

  use abort
  use eos_sn2
  use totare_hy
  use configure

  implicit none

  real(kind=rk), intent (in)                      :: tem
  class(function_arguments_t), intent(in), target :: additional_args
  logical, intent(out), optional                  :: error
  real(kind=rk)                                   :: val

  ! locals
  real(kind=rk)                 :: temcopy(1)
  real(kind=rk), dimension(1)   :: rhocopy, xhrep, ede, pre, gamc, sto, &
                                   ccu, cce, ccp, ccn
  real(kind=rk), dimension(1,2) :: za
  real(kind=rk), dimension(1,config%qn) :: xnucopy
  real(kind=rk) :: eos_self(2), eos_children(2)
  logical :: eos_error

  class(eos_field_arg_t), pointer :: eos_field_args

  select type(additional_args)
    type is (eos_field_arg_t)
      eos_field_args => additional_args
    class default
      raise_abort("Unsupported argument type!")
  end select

  rhocopy(1) = eos_field_args%rho
  temcopy(1) = tem
  xnucopy(1,:) = eos_field_args%xnu(:)

  ede(:) = 0.0_rk
  eos_error = .false.
  eos_self(:) = 0.0_rk
  eos_children(:) = 0.0_rk

  if (present(error)) then
    error = .false.
  endif

  call eos(rhocopy, temcopy, xnucopy, xhrep, &
           za, ede, pre, gamc, sto, &
           ccu, cce, ccp, ccn, &
           eos_self, eos_children, &
           mode = 1, nsemode = 0, ler = eos_error)

  if (eos_error) then
    if (present(error)) then
      val = 0.0
      error = eos_error
      return
    else
      raise_abort("eos() call resulted in error, aborting")
    endif
  endif

  select case(eos_field_args%field)
    case(EOS_FIELD_PRESSURE)
      val = pre(1)
    case(EOS_FIELD_ENTROPY)
      val = sto(1)
    case default
      raise_abort("Unknown eos field selector")
  end select

end function eos_field

#endif

!> Given a progenitor with an allocated pressure or entropy field, calculate
!> the temperature such that the code EoS reproduces the progenitor in this
!> quantity
!>
!> \param       star    A progenitor_t structure
!> \param       field   The quantity which shall be reproduced, one of
!>                      EOS_FIELD_PRESSURE, EOS_FIELD_ENTROPY
!>
!> In some cases, the field cannot be reproduced at all by our EoS,
!> there we still use the temperature of the progenitor and print a
!> warning
!>
!> \author      Janina von Groote, Lorenz Huedepohl
#if 0
subroutine eos_inversion(star, field)
  use precision
!  use secant_method
  use abort
  use configure


  type(progenitor_t), intent(inout) :: star
  integer, intent(in) :: field
  real(kind=rk) :: temunder, temabove, val, target_val
  type(eos_field_arg_t), target :: additional_args
  real(kind=rk), dimension(size(star%xni, 2) + 1), target :: xnu
  integer :: i
  logical :: eos_error

  additional_args%field = field

!  ----------------------------------
!  Solve f(rho,T,Y_e,X_i) = f_0 for T
!  ----------------------------------
  zone_loop: do i = 1, star%n_zones
    if (i > 2) then
      if (star%rqi(i-2) .gt. config%gridlx) then
        write(*,*) "Omitting zone", i, "(outside computational domain)"
        cycle
      endif
    endif

    select case(field)
      case(EOS_FIELD_PRESSURE)
        target_val = star%pri(i)
      case(EOS_FIELD_ENTROPY)
        target_val = star%sti(i)
      case default
        raise_abort("Unknown eos field selector")
    end select

    xnu(:) = [star%xni(i, :), star%yei(i)]
    additional_args%xnu => xnu
    additional_args%rho = star%rhoi(i)

    ! Get enclosing values
    temunder = star%tmpi(i)
    temabove = star%tmpi(i)

    ! Make sure temunder and temabove enclose the solution
    val = eos_field(star%tmpi(i), additional_args, eos_error)
    abort_if(eos_error)
    do while (val > target_val)
      temunder = 0.95_rk * temunder
      val = eos_field(temunder, additional_args, eos_error)
      if (eos_error) then
        write(*,*) "Unable to invert zone ", i, ", using temperature here!"
        cycle zone_loop
      endif
    end do

    ! same for temabove
    val = eos_field(star%tmpi(i), additional_args, eos_error)
    abort_if(eos_error)
    do while (val < target_val)
      temabove = 1.05_rk * temabove
      val = eos_field(temabove, additional_args, eos_error)
      if (eos_error) then
        write(*,*) "Unable to invert zone ", i, ", using temperature here!"
        cycle zone_loop
      endif
    end do

    ! We found enclosing values, determine solution via secant method
!    star%tmpi(i) = secant(eos_field,temunder,temabove,target_val,additional_args,1e-6_rk,100_ik)
    write(*,*) "Inverted zone", i, "with a temperature of", star%tmpi(i), &
      "(1-f/f_star) = ", (1.0_rk - eos_field(star%tmpi(i), additional_args) / target_val)

  end do zone_loop

end subroutine eos_inversion
#endif

!> read initial 1D-model
function read_inimod(filename) result(progenitor)
  use precision
  use abort

  use charac
  use movgrid_hy

  use configure

  use stringutils, only : endswith
  implicit none
  character*(*), intent(in) :: filename
  type(progenitor_t)        :: progenitor


!  progenitor = read_star(filename)

  if (endswith(filename, ".h5")) then
#ifdef WRITE_BINARY_OUTPUT
    progenitor = read_star_hdf5(filename)
#endif
  else if (endswith(filename, "new")) then
    progenitor = read_star_models_new(filename)
  else if (endswith(filename, ".std")) then
    progenitor = read_star(filename)
  else if (endswith(filename, ".inp")) then
    progenitor = PRO_WSLY(filename)
  else
    progenitor = read_star(filename)
!    raise_abort("Unknown progenitor file type!")
  end if

  select case(config%setup_mode)
    case("rho-t-ye")
    case("rho-p-ye")
!      if (.not. allocated(progenitor%pri)) then
!        raise_abort("Progenitor file does not contain pressure! Change setup_mode")
!      endif
!      call eos_inversion(progenitor, EOS_FIELD_PRESSURE)
      raise_abort("Invalid value for setup_mode: '" // trim(config%setup_mode)  // "'")
    case("rho-s-ye")
!      if (.not. allocated(progenitor%sti)) then
!        raise_abort("Progenitor file does not contain entropy! Change setup_mode")
!      endif
!      call eos_inversion(progenitor, EOS_FIELD_ENTROPY)
      raise_abort("Invalid value for setup_mode: '" // trim(config%setup_mode)  // "'")
    case default
      raise_abort("Invalid value for setup_mode: '" // trim(config%setup_mode)  // "'")
  end select

  call mass2radius(progenitor%rri,progenitor%rhoi,config%xmfrac,rstr,progenitor%n_zones)

end function read_inimod

! --------------------------------------------------------------

subroutine init(progenitor)

! --------------------------------------------------------------
! Autor        : Markus Rampp
! Modul        : $Id: coc_init.F,v 1.8 2003/02/27 17:02:54 mjr Exp $
! Version      : $Revision$
! Date         : $Date$
! Header       : $Name:  $
!
! Purpose: define initial model for Core Collapse
!
! Similar to pns_init.F written by W.Keil
! --------------------------------------------------------------
  use precision

  use intgrs_hy
  use totgrq_hy
!  use charac  ! forcheck
  use nutrio_hy
  use marker_hy
  use massio_hy
  use revsho_hy
  use gfloat_hy
  use phycon

  !      use vnew_hy
  !      use mesh_hy
  use totare_hy

  use perturbation
  use abort

  use mo_mpi

  use eos3d_routine, only : eos3d
  use interpol_eul
  use cpyare_mod
#ifdef CFC_TRANSPORT
  use gr_hyd_init
#endif
  use print_stdout_mod

  use hydro_areas_mod
  use configure
  implicit none
  ! LOCAL variables that are not in modules

  logical                    :: eos_error
  integer(kind=ik)           :: i, j, k, l
  real(kind=rk)              :: sumass, scr, scrdv, xtheta_j, gin
  integer(kind=ik)           :: ipos(config%qx)
  real(kind=rk) ,allocatable :: rand(:,:,:)  ! only used for 2D and 3D
  integer(kind=ik)           :: ierr
  real(kind=rk)              :: poisson_self(2),      &
                                poisson_children(2),  &
                                eos3d_self(2),        &
                                eos3d_children(2)

!  real(kind=rk) :: omscal(q) ! for rotation

  type(progenitor_t) :: progenitor

!-----------------------------------------------------------------------
!     set array boundaries:
!-----------------------------------------------------------------------

  nzn  = config%qx
  nzn1 = nzn + 1
  nzn2 = nzn + 2
  nzn3 = nzn + 3
  nzn4 = nzn + 4
  nzn5 = nzn + 5
  nzn6 = nzn + 6
  nzn7 = nzn + 7
  nzn8 = nzn + 8

  gamma = 4.0_rk/3.0_rk



! interpolate initial model onto hydro-grid


  do k = qz_s, qz_e
     do j = qy_s, qy_e

        dentot(:,j,k)   = 0.0_rk
        temtot(:,j,k)   = 0.0_rk
        vextot(:,j,k)   = 0.0_rk
        xnutot(:,j,k,:) = 0.0_rk


        call sort_vec(progenitor%rqi,progenitor%n_zones,1,xzntot(1:),ipos(1:),config%qx,1)

        call monintp(progenitor%rqi,progenitor%rhoi,progenitor%n_zones,1,ipos(1:), &
                     xzntot(1:),dentot(1:,j,k),config%qx,1,99,99)
        call monintp(progenitor%rqi,progenitor%tmpi,progenitor%n_zones,1,ipos(1:), &
                     xzntot(1:),temtot(1:,j,k),config%qx,1,99,99)
        call monintp(progenitor%rqi,progenitor%vqi  ,progenitor%n_zones,1,ipos(1:), &
                     xzntot(1:),vextot(1:,j,k),config%qx,1,99,99)


        do l = 1, config%qn-1
           call linintp(progenitor%rqi,progenitor%xni(1:,l),progenitor%n_zones,1,ipos(1:), &
                        xzntot(1:),xnutot(1:,j,k,l),config%qx,1)
        enddo
        l=config%qn
        call linintp(progenitor%rqi,progenitor%yei,progenitor%n_zones,1,ipos(1:), &
                     xzntot(1:),xnutot(1:,j,k,l),config%qx,1)



        veytot(:,j,k) = 0.0_rk
        veztot(:,j,k) = 0.0_rk
! rotation
!         where(xzntot .lt. 1.74703e8)
!            omscal=0.5
!         elsewhere
!            omscal=0.5*(xzntot/1.74703e8)**(-1.5)
!         endwhere
        call linintp(progenitor%rqi,progenitor%omi,progenitor%n_zones,1,ipos(1:), &
                     xzntot(1:),veztot(1:,j,k),config%qx,1)
        veztot(1:config%qx,j,k) = veztot(1:config%qx,j,k)*xzntot(1:config%qx)*sin(yzntot(j))

        acxtot(:,j,k) = 0.0_rk
        enetot(:,j,k) = -1.e33_rk ! will be set in eos3d


     enddo
  enddo

! the mapping is not conservative, so better check
!      write(*,*) vlfrac,dvztot,dvytot,dvxtot(1)
  sumass=0.0_rk

  dvxtot(1:config%qx) = (xzrtot(1:config%qx)**3 - xzltot(1:config%qx)**3)  &
         / 3._rk

  if (config%nsdim .ge. 2)  then
     do j = 1, config%qy
        dvytot(j) = cos(yzltot(j)) - cos(yzrtot(j))
     enddo
     if (config%nsdim .eq. 2)   dvztot(1) = 2.0_rk * pc_pi
  end if

  if (config%nsdim  .eq. 3)  then
     do k = 1, config%qz
        dvztot(k) = zzrtot(k) - zzltot(k)
     enddo
  end if
  
  
  if (myproc .eq. 0) then
       scr= 4.0_rk/3.0_rk*pc_pi
     do i = 1, config%qx
        scrdv=xzrtot(i)**3-xzltot(i)**3
        sumass=sumass+dentot(i,qy_s,qz_e)* scr*scrdv
     enddo
  endif


! --------- IN THE 2D CASE IMPOSE RANDOM PERTURBATIONS

! in the 2D case impose non-spherical perturbations
  if (config%qy .gt.1 .and. config%noise .ne. 0) then

     IF (config%noise .eq. 11) then
        allocate (rand(config%qx,1:config%qy,1:config%qz))
     else if (config%noise .eq. 12) then
        allocate (rand(config%qx,1:config%qy,1:config%qz*2))
     else
        allocate (rand(config%qx,qy_s:qy_e,qz_s:qz_e))
     end if

     call RANDOM_SEED
     call RANDOM_NUMBER(rand)
     select case(config%noise)
     case(1)  !randomly perturb velocity
        vextot(:,:,:)=vextot(:,:,:)*(1.0_rk+2.0_rk* &
                       (rand(:,:,:)-0.5_rk)*config%ampl)
     case(2)  !randomly perturb density
        dentot(:,:,:)=dentot(:,:,:)*(1.0_rk+2.0_rk* &
                       (rand(:,:,:)-0.5_rk)*config%ampl)
     case(11) !density perturbation from file
        open (70, file = "density_pert.i3e", form = "unformatted")
        read (70) rand
        close (70)
        dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= &
             dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)* &
             (1.0_rk+config%ampl*rand(1:config%qx,qy_s:qy_e,qz_s:qz_e))
     case(12) !velocity perturbation from file
        open (70, file = "velocity_pert.i3e", form = "unformatted")
        read (70) rand(1:config%qx,1:config%qy,1:config%qz), &
             rand(1:config%qx,1:config%qy,config%qz+1:2*config%qz)
        close (70)
        vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= &
             vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)+ &
             config%ampl*rand(1:config%qx,qy_s:qy_e,1:config%qz)
        veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= &
             veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)+ &
             config%ampl*rand(1:config%qx,qy_s:qy_e,qz_s+config%qz:qz_e+config%qz)

     case(22) !sinusoidal-mode
         k=1

         do j = qy_s, qy_e

            xtheta_j=(yzntot(j)-yzltot(1))/(yzrtot(config%qy)-yzltot(1))
            dentot(:,j,k)=dentot(:,j,k)* &
                          (1.0_rk+sin(xtheta_j*2*pc_pi)*config%ampl)
         enddo
      case(-1) !deterministically perturb velocity and density

         do k= qz_s, qz_e
            do j = qy_s, qy_e

               do i=1,config%qx
                  rand(i,j,k)=sin(real(i,kind=rk))* &
                              cos(real(j,kind=rk))* &
                              cos(real(k,kind=rk))
               enddo
            enddo
         enddo
         dentot(:,:,:)=dentot(:,:,:)*(1.0_rk+2.0_rk* &
                       (rand(:,:,:)-0.5_rk)*config%ampl)
         vextot(:,:,:)=vextot(:,:,:)*(1.0_rk+2.0_rk* &
                      (rand(:,:,:)-0.5_rk)*config%ampl)
      case default
         raise_abort("coc_init(): nocase")
      end select
      deallocate (rand)
   endif

!c Initialize parameters for grdvel
!      call massfrac(Rstr)
!      write(*,*) 'Rstr=',Rstr


   call printit_taskX(0," ")
   call printit_taskX(0," init> Interpolation of model: ")
   call printit_taskX(0," ")

   call printit_taskX(0," total Mass on grid:",sumass/pc_ms)
   call printit_taskX(0,"Task ",myproc,"   r [cm]     rho [g/cc]     T [K]      Y_e [1/by]   v_x [cm/s] ")


   j = qy_e
   k = qz_e


   if (myproc.eq.0) then

      open(17,file='ppm.grd',form='formatted')

      do i = 1, config%qx
         write(*,'(1x,i4,5(1x,1pe12.5))') &
                 i,xzntot(i),dentot(i,j,k),temtot(i,j,k), &
                 xnutot(i,j,k,config%qn),vextot(i,j,k)
         write(17,'(1x,I4,4(1x,1pe13.6))') &
                 i,xzntot(i),dentot(i,j,k), &
                 temtot(i,j,k),xnutot(i,j,k,config%qn)
      enddo

      close(17)

   end if

   call printit_taskX(0," ")

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)


! -- calculate remaining state variables:

#ifdef CFC_TRANSPORT
   call cpyare(2)
   call init_cfc_hydro(.false.)
#else
   call cpyare(2)            ! EOS cannot use ***tot-arrays
   areas%nx=config%qx
   areas%ny=config%qy
   areas%nz=config%qz

   call eos3d (1,eos_error, eos3d_self, eos3d_children)
   abort_if(eos_error)
   call cpyare(0)
#endif /* CFC_TRANSPORT */

! initialize neutrino quantities:
!     it is possible that there is an output before they are
!     calculated by the neutrino transport

   qyetot(:,:,:,:) = 0.0_rk
   qentot(:,:,:) = 0.0_rk
   qmotot(:,:,:) = 0.0_rk
   qmytot(:,:,:) = 0.0_rk

   fnutot(:,:,:,:) = 0.0_rk
   enutot(:,:,:,:) = 0.0_rk
   dnutot(:,:,:,:) = 0.0_rk
   pnutot(:,:,:,:) = 0.0_rk
   gnutot(:,:,:,:) = 0.0_rk

   eltobs(:) = 0.0_rk
   eltobc(:) = 0.0_rk
   etnobs(:) = 0.0_rk

! reset marker particles:
   ppx(:) = 0.0_rk
   ppy(:) = 0.0_rk
   ppz(:) = 0.0_rk
   pvx(:) = 0.0_rk
   pvy(:) = 0.0_rk
   pvz(:) = 0.0_rk
   pma(:) = 0.0_rk

! set inflow quantities:
   dmdtio(:,:,:) = 0.0_rk

   rhoin  = 1.0_rk
   uin    = 0.0_rk
   utin   = 0.0_rk
   uttin  = 0.0_rk
   pin    = 1.0_rk
   gin    = 0.0_rk
   tin    = 1.0_rk
   gamein = gamma
   gamcin = gamma
   ein    = pin/((gamein-1.0_rk)*rhoin)+0.5_rk*(uin**2+utin**2+uttin**2)
   xnin(:) = 1.0_rk


! define some grid related quantities:

   tsave(:,:)=0.0_rk
   s0_w=0.0_rk
   nop = 0
   igeom = config%igeomx

   ugrtot(:) = 0.0_rk

   dvxtot(1:config%qx) = (xzrtot(1:config%qx)**3 - xzltot(1:config%qx)**3) / 3._rk
   dvytot(1) = 1.0_rk
   dvztot(1) = 1.0_rk
   srfint = 1.0_rk
   vlfrac = 4.0_rk * pc_pi

   if (config%nsdim .ge. 2)  then
      call printit_taskX(0,"pns_init> j/dvy")

      if (config%igeomy .eq. 4)  then
         srfint = cos(yzltot(1)) - cos(yzrtot(config%qy))
         !-PNS            srfint = 2. - float(isym)
         vlfrac = 2.0_rk / srfint
         do j = 1, config%qy
            dvytot(j) = cos(yzltot(j)) - cos(yzrtot(j))
            if (myproc .eq. 0) then
               write(*,'('' Task '',i4,i5,1x,1pe12.5)') myproc,j,dvytot(j)
            endif
         enddo

         call printit_taskX(0," ")

         if (config%nsdim .eq. 2)   dvztot(1) = 2.0_rk * pc_pi
         if (config%nsdim .eq. 2) then
            call printit_taskX(0,"pns_init> srfint = ",srfint)
            call printit_taskX(0,"          vlfrac = ",vlfrac)
         endif

      else
         raise_abort("init(): case not implemented")
      end if
   end if


   if (config%nsdim  .eq. 3)  then
      call printit_taskX(0,"pns_init> k/dvz")

      if (config%igeomx .eq. 2  .and. &
          config%igeomy .eq. 4  .and.  config%igeomz .eq. 5 .and. &
          ((config%bndmny .eq. 4 .and. config%bndmnz .eq. 4) .or. &
          (config%bndmny .eq. 1 .and. config%bndmnz .eq. 4) )) then

             srfint  = 2.0_rk * config%gridlz * sin(0.5_rk * config%gridly)
             vlfrac  = 4.0_rk * pc_pi / srfint
             do k = 1, config%qz
                dvztot(k) = zzrtot(k) - zzltot(k)
                if (myproc .eq. 0) then
                   write(*,'('' Task '',i4,i5,1x,1pe12.5)') myproc, k,dvztot(k)
                endif

             enddo
             call printit_taskX(0," ")
             call printit_taskX(0,"pns_init> srfint = ",srfint)
             call printit_taskX(0,"          vlfrac = ",vlfrac)
          else
             raise_abort("init(): case not implemented")
          end if
       end if

! -- compute potential:

       gpotot(:,:,:)=0.0_rk

       if(config%i_grv .eq. 0) then
          call printit_taskX(0,"init> Newtonian potential")
       else
          call printit_taskX(0,"init> general relativistic potential")
       end if

       ephtot(:) = 1.0_rk
       gamtot(:) = 1.0_rk
       gamold(:) = 1.0_rk

       igrav = 1
       igrav = 0
       call printit_taskX(0,"=========================================================================")


     end subroutine init

     subroutine fix_mass_fractions(xnu, ye)
       use nucparam, only : pc_nuc, name_xnuc
       use mo_mpi
       implicit none
       real(kind=rk) :: xnu(:,:)
       real(kind=rk) :: ye(size(xnu, dim=1))
       integer :: i,n,k_max

       ! attribute the error due to imprecisions of input data
       ! to the dominant species, in order to enforce sum(x_i) == 1
       do i=1, size(xnu, dim=1)
          if (myproc .eq. 0 .and. abs(sum(xnu(i,:)) - 1.0_rk) > 1e-6) then
            write(*,'(a,i4,a,es12.5)') "Sum of mass fractions in zone ", i, &
                " of progenitor deviate from 1 by", abs(sum(xnu(i,:)) - 1.0_rk)
            do n=1, size(xnu, dim=2)
              write(*,'(a,1x,f20.18)') name_xnuc(n), xnu(i, n)
            end do
            write(*,*)
          endif
          k_max=maxloc(xnu(i,:), dim = 1)
          xnu(i,k_max)=1.0_rk - (sum(xnu(i,:))-xnu(i,k_max))
       enddo

       ! check charge neutrality
       do i=1, size(xnu, dim=1)
          if (myproc .eq. 0 .and. abs(sum(xnu(i,:) * (pc_nuc(:,1) / pc_nuc(:,2))) - ye(i)) > 1e-6_rk) then
            write(*,'(a,i4,a,es12.5)') "Charge neutrality in zone ", i, &
                " of progenitor violated by ", abs(sum(xnu(i,:) * (pc_nuc(:,1) / pc_nuc(:,2))) - ye(i))
            write(*,*)
          endif
       enddo
     end subroutine fix_mass_fractions
#ifdef WRITE_BINARY_OUTPUT
     !> read data of initial stellar model in HDF5 format
     !>
     !> \param          stellar_model   file name of progenitor
     !> \returns        p               an allocated progenitor structure,
     !>                                 the calling routine has to deallocate it!
     function read_star_hdf5(filename) result(p)
       use precision

       use phycon
       use nucparam, only : name_xnuc, n_fe56
       use print_stdout_mod

       use configure
       use abort
       use dataformat
       implicit none

       character(*), intent (in) :: filename
       type(progenitor_t) :: p

       integer(kind=ik) :: n_zones, n_species, i, n
       real(kind=rk), allocatable :: xnu(:,:)
       character(len=5), allocatable :: nuclei(:)
       integer(kind=ik), allocatable :: xni_index(:)
       type(datafile) :: f

       call open_data_file(f, "tables/" // trim(filename))

       call readq(f, "vertex/n_zones", n_zones)
       call readq(f, "vertex/n_species", n_species)

       p = allocate_progenitor(n_zones)

       call readq(f, "vertex/rri", p%rri)
       call readq(f, "vertex/rqi", p%rqi)
       call readq(f, "vertex/rhoi", p%rhoi)
       call readq(f, "vertex/tmpi", p%tmpi)

!PERL-START for my $var ("sti", "pri", "eni") {
       if (exists_in_file(f, "vertex/@[[$var]]")) then
         write(*,*) "Progenitor file contains @[[$var]]"
         allocate(p%@[[$var]](1:p%n_zones))
         call readq(f, "vertex/@[[$var]]", p%@[[$var]])
       else
         write(*,*) "Progenitor file does not contain @[[$var]]"
       endif
!PERL-END }

       call readq(f, "vertex/yei", p%yei)
       call readq(f, "vertex/vqi", p%vqi)
       call readq(f, "vertex/vai", p%omi)

       allocate(xnu(n_zones,n_species))
       allocate(nuclei(n_species))
       call readq(f, "vertex/xni", xnu)
       call readq(f, "vertex/nuclei", nuclei)

       ! sort in the composition
       allocate(xni_index(n_species))
       do i = 1, n_species
         xni_index(i) = 0
         ! search corresponding index in code arrays
         do n = 1, config%qn - 1
           if (trim(adjustl(nuclei(i))) .eq. trim(adjustl(name_xnuc(n)))) then
             xni_index(i) = n
             exit
           endif
         end do
         write(*,*) nuclei(i), xni_index(i)
       end do

       ! \todo here: temperature inversion via p, E, or s

       ! fill in the composition array, whose nuclei setup might differ
       ! from the one used in the code
       p%xni = 1e-25_rk
       do i = 1, n_zones
         do n = 1, size(nuclei)
           if (xni_index(n) == 0) then
             if (xnu(i,n) > 1e-25_rk) then
               ! Have to deal with non-standard nuclei
               select case(trim(adjustl(nuclei(n))))
                 case("Fe54")
                   p%xni(i,n_fe56) = p%xni(i,n_fe56) + xnu(i,n)
                 case("Fe")
                   p%xni(i,n_fe56) = p%xni(i,n_fe56) + xnu(i,n)
                 case default
                   p%xni(i,1) = p%xni(i,1) + xnu(i,1) / 3.0_rk
                   p%xni(i,2) = p%xni(i,2) + xnu(i,2) / 3.0_rk * 2.0_rk
               end select
             end if
           else
             p%xni(i,xni_index(n)) = p%xni(i,xni_index(n)) + xnu(i,n)
           end if
         end do
       end do

       do i = 1, n_zones
          if (sum(xnu(i,1:size(nuclei))).lt.0.5_rk) then
             p%xni (i,n_fe56) = 1.0_rk
          end if
       end do

       call close_data_file(f)

       call fix_mass_fractions(p%xni, p%yei)

       call printit_taskX(0," ")
       call printit_taskX(0,"read_star_hdf5> Max. Radius   [cm]   = ",p%rri(n_zones))
       call printit_taskX(0,"==========================================================================")
       call printit_taskX(0," ")

     end function read_star_hdf5
#endif

     !> read data of initial stellar model
     !>
     !> \param          stellar_model   file name of progenitor
     !> \returns        p               an allocated progenitor structure,
     !>                                 the calling routine has to deallocate it!
     function read_star (stellar_model) result(p)
       use precision

       use phycon
       use nucparam
       use print_stdout_mod

       use configure
       use abort
       implicit none
! LOCAL variables that are not in modules

       character(*), intent (in) :: stellar_model

       character(LEN=256)        :: string1,string2
       integer(kind=ik) :: nn,i,iz_rd,inuc,kn,nnuc,ndat
       real(kind=rk)              :: wc_mas(config%qn-1),scrz,xnucleon,xp

       real(kind=rk)              :: xm_rd,rr_rd,ve_rd,ro_rd,tt_rd,om_rd, &
                                    du_rd,ye_rd,xnn_rd,xpp_rd,xhe3_rd, &
                                    xhe4_rd,xc_rd,xn_rd,xox_rd,xne_rd, &
                                    xmg_rd,xsi_rd,xs_rd,xar_rd,xca_rd, &
                                    xti_rd,xcr_rd,xfe52_rd,            &
                                    xfe54_rd,xni_rd,xfe56_rd,xfe_rd,ve_old
       type(progenitor_t) :: p


       do nn=1,config%qn-1
          wc_mas(nn)=pc_nuc(nn,3)+pc_nuc(nn,2)*wc_mb
       enddo


       open(8,file = './tables/'//stellar_model,form = 'formatted')

       read(8,'(I5)') ndat
       read(8,'(A256)') string1
       read(8,'(A256)') string2
!       read(8,*)

       call printit_taskX(0,"read_star> Model/n: ",stellar_model)
       call printit_taskX(0,"                    ",ndat)
       call printit_taskX(0,"read_star> Comment: ",string1)
       call printit_taskX(0,"read_star> Comment: ",string2)
       call printit_taskX(0," ")

       p = allocate_progenitor(ndat)

       p%rri(0)=0.0_rk
       ve_old=0.0_rk

 2001 format(I4,31(2X,1PE23.16))
! 2002 format (I6,2X,11(E23.16,2X),25X,20(2X,E23.16))
 2002 format (I4,4X,11(E23.16,2X),25X,20(2X,E23.16))
! 2002 format (I4,4X,11(E23.16,2X),26X,20(2X,E23.16))
 2014 format (I6,12(X,E24.17),20(X,E24.17))

       do i=1,ndat
          ! read(8,2002) &
          !      iz_rd,xm_rd,rr_rd,ve_rd,ro_rd,tt_rd,            &
          !      du_rd,du_rd,du_rd,om_rd,du_rd,ye_rd,            &
          !      xnn_rd,xpp_rd,xhe3_rd,xhe4_rd,xc_rd,            &
          !      xn_rd,xox_rd,xne_rd,xmg_rd,xsi_rd,              &
          !      xs_rd,xar_rd,xca_rd,xti_rd,xcr_rd,              &
          !      xfe52_rd,xfe54_rd,xni_rd,xfe56_rd,xfe_rd
          read(8,2014) &
               iz_rd,du_rd,xm_rd,rr_rd,ve_rd,ro_rd,tt_rd,      &
               du_rd,du_rd,du_rd,om_rd,du_rd,ye_rd,            &
               xnn_rd,xpp_rd,xhe3_rd,xhe4_rd,xc_rd,                    &
               xn_rd,xox_rd,xne_rd,xmg_rd,xsi_rd,              &
               xs_rd,xar_rd,xca_rd,xti_rd,xcr_rd,              &
               xfe52_rd,xfe54_rd,xni_rd,xfe56_rd,xfe_rd
          print *,i,iz_rd,xm_rd,rr_rd
          ! sanity check
          abort_if(iz_rd /= i)

          p%rri(i)=rr_rd
          p%rqi(i)=0.5_rk*(p%rri(i)+p%rri(i-1))
          p%tmpi(i)=tt_rd
          p%rhoi(i)=ro_rd
          p%vqi(i)=0.5_rk*(ve_rd + ve_old)
          ve_old = ve_rd
          p%omi(i)=0.0!om_rd
          p%yei(i)=ye_rd

          p%xni(i,:)     = 0.0_rk
          p%xni(i,n_n)   = xnn_rd
          p%xni(i,n_p)   = xpp_rd
          if (n_he3 .ne. -99999) then
             p%xni(i,n_he3)  = xhe3_rd
             p%xni(i,n_he4)  = xhe4_rd
          else
             p%xni(i,n_he4)  = xhe3_rd+xhe4_rd
          endif

          p%xni(i,n_c12)   = xc_rd
          p%xni(i,n_o16)   = xox_rd
          p%xni(i,n_ne20)  = xne_rd
          p%xni(i,n_mg24)  = xmg_rd
          p%xni(i,n_si28)  = xsi_rd
          p%xni(i,n_s32)   = xs_rd
          p%xni(i,n_ar36)  = xar_rd
          p%xni(i,n_ca40)  = xca_rd
          p%xni(i,n_ti44)  = xti_rd
          p%xni(i,n_cr48)  = xcr_rd
          p%xni(i,n_ni56)  = xni_rd
          if (n_fe52 .ne. -99999) then
             p%xni(i,n_fe52) = xfe52_rd
             p%xni(i,n_fe56) = xfe54_rd+xfe56_rd+xfe_rd
          else
             p%xni(i,n_fe56)= xfe52_rd+xfe54_rd+xfe56_rd+xfe_rd
          endif
       enddo

       call fix_mass_fractions(p%xni, p%yei)

       call printit_taskX(0," ")
       call printit_taskX(0,"read_star> Enclosed Mass [M_sol]= ",du_rd)
       call printit_taskX(0,"read_star> Max. Radius   [cm]   = ",p%rri(ndat))
       call printit_taskX(0,"==========================================================================")
       call printit_taskX(0," ")
     end function read_star


     !> read data of initial stellar model in ".new_new" format
     !>
     !> \param          stellar_model   file name of progenitor
     !> \returns        p               an allocated progenitor structure,
     !>                                 the calling routine has to deallocate it!
     function read_star_models_new(stellar_model) result(p)
       use precision
       use phycon
       use nucparam
       use print_stdout_mod

       use configure
       use abort
       implicit none

       character(*), intent (in) :: stellar_model

       character(LEN=256)        :: string1

       integer(kind=ik) :: nn,i,iz_rd,inuc,kn,nnuc,id,ndat
       real(kind=rk)              :: wc_mas(config%qn-1),scrz,xnucleon,xp

       real(kind=rk)              :: xm_rd,rr_rd,ve_rd,ro_rd,tt_rd,   &
                                     ye_rd,xnn_rd,xpp_rd,xhe3_rd,&
                                     xhe4_rd,xc_rd,xn14_rd,xox_rd,xne_rd,&
                                     xmg_rd,xsi_rd,xs_rd,xar_rd,xca_rd,&
                                     xti_rd,xcr_rd,xfe52_rd,           &
                                     xfe54_rd,xni_rd,xfe56_rd,xfe_rd,  &
                                     xni120_rd,xzr200_rd, &
                                     xmn_rd,xni70_rd,xfe60_rd,ve_old
       type(progenitor_t) :: p

       do nn=1,config%qn-1
          wc_mas(nn)=pc_nuc(nn,3)+pc_nuc(nn,2)*wc_mb
       enddo

       open(8,file = './tables/'//stellar_model,form = 'formatted')

       read(8,*) ndat,id
       if (id .ne. 20040212) then
          call printit_taskX(0," ")
          call printit_taskX(0,"ID wrong!! Check progenitor model")
          call printit_taskX(0," ")
       endif
       read(8,'(A256)') string1
       read(8,*)
       call printit_taskX(0,"read_star> Model/n: ",stellar_model)
       call printit_taskX(0,"                    ",ndat)
       call printit_taskX(0,"read_star> Comment: ",string1)
       call printit_taskX(0," ")

       p = allocate_progenitor(ndat)

       p%rri(0)=0.0_rk
       ve_old =0.0_rk

       do i=1,ndat

          if (config%use_network) then
             call printit_taskX(0,"CHECK read in compostion with NETWORK")
             read(8,*) iz_rd,rr_rd,ve_rd,ro_rd,tt_rd,                     &
                       ye_rd,xnn_rd,xpp_rd,xhe3_rd,xhe4_rd,xc_rd,xn14_rd, &
                       xox_rd,xne_rd , xmg_rd,xsi_rd,xs_rd,xar_rd,        &
                       xca_rd,xti_rd,xcr_rd,xmn_rd ,xfe56_rd,xfe60_rd,    &
                       xni_rd,xni70_rd,xni120_rd !,xzr200_rd

          else
             read(8,*) iz_rd,rr_rd,ve_rd,ro_rd,tt_rd,                     &
                       ye_rd,xnn_rd,xpp_rd,xhe4_rd,xc_rd,                 &
                       xox_rd,xne_rd , xmg_rd,xsi_rd,xs_rd,xar_rd,        &
                       xca_rd,xti_rd,xcr_rd,xmn_rd ,xfe56_rd,xfe60_rd,    &
                       xni_rd,xni70_rd,xni120_rd,xzr200_rd

          endif

! order in prog: zone, radius, vel, den, tem,y_e,n,p,he4,c12,o16,ne20,mg24
! si28,s32,ar36,ca40,ti44,cr48,mn54,fe56,fe60,ni56,ni70,ni129,zr200

          p%rri(i)=rr_rd
          p%rqi(i)=0.5_rk*(p%rri(i)+p%rri(i-1))
          p%tmpi(i)=tt_rd
          p%rhoi(i)=ro_rd
          p%vqi(i)=0.5_rk*(ve_rd + ve_old)
          ve_old = ve_rd
          p%yei(i)=ye_rd


          p%xni(i,:)     = 0.0_rk
          p%xni(i,n_n)   = xnn_rd
          p%xni(i,n_p)   = xpp_rd

          if (config%use_network) then
             !       p%xni(i,n_he3)  = xhe3_rd
          else
             !       p%xni(i,n_he3)  = 0.0_rk
          endif
          p%xni(i,n_he4)  = xhe4_rd
          p%xni(i,n_c12)   = xc_rd
          if (config%use_network) then
             !       p%xni(i,n_n14)  = xn14_rd
          else
             !       p%xni(i,n_n14)  = 0.0_rk
          endif
          p%xni(i,n_o16)   = xox_rd
          p%xni(i,n_ne20)  = xne_rd
          p%xni(i,n_mg24)  = xmg_rd
          p%xni(i,n_si28)  = xsi_rd
          p%xni(i,n_s32)   = xs_rd
          p%xni(i,n_ar36)  = xar_rd
          p%xni(i,n_ca40)  = xca_rd
          p%xni(i,n_ti44)  = xti_rd
          p%xni(i,n_cr48)  = xcr_rd
          p%xni(i,n_mn54)  = xmn_rd
          p%xni(i,n_fe56)  = xfe56_rd
          p%xni(i,n_fe60)  = xfe60_rd
          p%xni(i,n_ni56)  = xni_rd
          p%xni(i,n_ni70)  = xni70_rd
          p%xni(i,n_ni120) = xni120_rd
          p%xni(i,n_zr200) = xzr200_rd
          !         xni(i,n_fe56)= xfe52_rd+xfe54_rd+xfe56_rd+xfe_rd
       enddo

       call fix_mass_fractions(p%xni, p%yei)

       call printit_taskX(0," ")
       call printit_taskX(0,"read_star> Max. Radius   [cm]   = ",p%rri(ndat))
       call printit_taskX(0," ")

     end function read_star_models_new

     !>   read Woosley progenitor model
     !>     input:  MODEL = file name
     function PRO_WSLY (model) result(p)

       use precision

       use phycon
       use nucparam
       use abort

       use configure
       implicit none

       integer(kind=ik)           :: i,inuc,kn,iz_rd,nn,ndat
       character(*) ,intent (in) :: model
       character(LEN=256)        :: string1
       real(kind=rk)              :: wc_mas(config%qn-1),xfe_rd,xsi_rd,xp, &
                                     scrz,xnucleon,dxfe,xmg_rd,xne_rd, &
                                     xox_rd,xnu_rd,xhe_rd,ye_rd, &
                                     xc_rd,rr_rd,tt_rd,ro_rd,ve_rd, &
                                     du_rd,ve_old
       type(progenitor_t) :: p

       do nn=1,config%qn-1
          wc_mas(nn)=pc_nuc(nn,3)+pc_nuc(nn,2)*wc_mb
       enddo

!      if (ngp .gt. 681) stop 'PRO_WSLY: too many points'

      open(8,file = './tables/'//model,form = 'formatted')
      read(8,'(I5)') ndat

      p = allocate_progenitor(ndat)

      write(*,'('' PRO_WSLY> Model  = '',A15)') model
      write(*,'('' PRO_WSLY> Shells = '',I3)') ndat

      write(*,*) '            mass        radius     temperature    density     velocity       Ye'
      write(*,*)

      p%rri(0)=0.0_rk
      ve_old=0.0_rk

      read(8,*)
      do i=1,ndat
         read(8,'(i5,6(1x,es12.5))') iz_rd,du_rd,rr_rd,tt_rd,ro_rd,ve_rd,ye_rd

         p%rri(i)=rr_rd
         p%rqi(i)=0.5_rk*(p%rri(i)+p%rri(i-1))
         p%tmpi(i)=tt_rd
         p%rhoi(i)=ro_rd
         p%vqi(i)=0.5_rk*(ve_rd + ve_old)
         ve_old = ve_rd
         p%yei(i)=ye_rd

         write(*,'(i8,6(1x,es12.5))') iz_rd,du_rd,p%rqi(i),p%tmpi(i),p%rhoi(i),p%vqi(i),p%yei(i)
      enddo
      read(8,*)
      read(8,*)
      do i=1,ndat
         read(8,'(i5,8(1x,es12.5))') iz_rd,xnu_rd,xhe_rd,xc_rd,xox_rd,xne_rd,xmg_rd, xsi_rd,xfe_rd
         p%xni(i,:)=0.0_rk
         p%xni(i,n_n)=0.0_rk
         p%xni(i,n_p)=xnu_rd
         p%xni(i,n_he4)=xhe_rd
         p%xni(i,n_c12)=xc_rd
         p%xni(i,n_o16)=xox_rd
         p%xni(i,n_ne20)=xne_rd
         p%xni(i,n_mg24)=xmg_rd
         p%xni(i,n_si28)=xsi_rd
         p%xni(i,n_fe56)=xfe_rd
      enddo

      close(8)

      call fix_mass_fractions(p%xni, p%yei)

      do i=1,ndat
         xnucleon=p%xni(i,n_p)
         scrz=0.0_rk
         do inuc=3,config%qn-1
            scrz=scrz+pc_nuc(inuc,1)*p%xni(i,inuc)/pc_nuc(inuc,2)
         enddo
!         xp=scrz-ye(i)
!         xp=max( min(xp,xnucleon) , 0.0 )
!         xi(i,1)=xnucleon-xp
!         xi(i,2)=xp
         xp=p%yei(i)-scrz
         if (xp.lt.0.0_rk) then
            p%xni(i,n_n)=xnucleon
            p%xni(i,n_p)=0.0_rk
            dxfe=xp/(pc_nuc(n_fe56,1)/pc_nuc(n_fe56,2)-pc_nuc(n_fe60,1)/pc_nuc(n_fe60,2))
            if (dxfe.lt.-p%xni(i,n_fe56)) then
               dxfe=xp/ (pc_nuc(n_fe56,1)/pc_nuc(n_fe56,2)-pc_nuc(n_ni70,1)/pc_nuc(n_ni70,2))
               if (dxfe.lt.-p%xni(i,n_fe56)) then
                  write (*,*) i,xp,p%yei(i),dxfe,p%xni(i,n_fe56)
                  raise_abort('PRO_WSLY> no charge neutrality!')
               endif
               p%xni(i,n_fe56)=p%xni(i,n_fe56)+dxfe
               p%xni(i,n_ni70)=p%xni(i,n_ni70)-dxfe
            else
               p%xni(i,n_fe56)=p%xni(i,n_fe56)+dxfe
               p%xni(i,n_fe60)=p%xni(i,n_fe60)-dxfe
            endif
         elseif (xp.gt.xnucleon) then
            p%xni(i,n_n)=0.0_rk
            p%xni(i,n_p)=xnucleon
            xp=xp-xnucleon
            dxfe=xp/ (pc_nuc(n_fe56,1)/pc_nuc(n_fe56,2)-pc_nuc(n_ni56,1)/pc_nuc(n_ni56,2))
            if (dxfe.lt.-p%xni(i,n_fe56)) then
               write (*,*) i,xp,p%yei(i),dxfe,p%xni(i,n_fe56)
               raise_abort('PRO_WSLY> no charge neutrality!')
            endif
            p%xni(i,n_fe56)=p%xni(i,n_fe56)+dxfe
            p%xni(i,n_ni56)=p%xni(i,n_ni56)-dxfe
         else
            p%xni(i,n_n)=xnucleon-xp
            p%xni(i,n_p)=xp
         endif
      enddo

      write(*,*)
      write(*,'('' PRO_WSLY> Enclosed Mass = '',1PE12.5)') du_rd
      write(*,'('' PRO_WSLY> Max. Radius   = '',1PE12.3)') p%rri(ndat)
      write(*,'(80(''=''))')

    end function PRO_WSLY



#ifdef NOMOTO_MODEL
!=======================================================================

    SUBROUTINE read_onemg

!=======================================================================
! Autor        : Markus Rampp / Bernhard Mueller
! Modul        : $Id: coc_init.F,v 1.8 2003/02/27 17:02:54 mjr Exp $
! Version      : $Revision$
! Date         : $Date$
! Header       : $Name:  $
!
! Purpose: define initial model for Core Collapse
!
! Similar to pns_init.F written by W.Keil
! --------------------------------------------------------------

      USE precision
!      USE arecon_hy
      USE intgrs_hy
      USE totgrq_hy
      USE charac
      USE nutrio_hy
      USE marker_hy
      USE massio_hy
      USE revsho_hy
      USE gfloat_hy
      USE phycon
      USE totare_hy
      use abort
#ifdef CFC_TRANSPORT
      USE size_cfc
      USE hydro_primitives_cfc
      USE gr_hyd_init
#endif
      use cpyare_mod
      use poisson_solver
      use eos3d_routine, only : eos3d

      USE eos_sn2

      USE mo_mpi
      USE error
      use print_stdout_mod

      use configure
      implicit none

      logical eos_error, err

      integer(kind=ik) :: i,j,k,ii,kn
      integer(kind=ik) :: iit
      real(kind=rk), dimension(config%qx) :: pre_tmp,xhrep_tmp,gamc_tmp, &
           ccu_tmp,cce_tmp,ccn_tmp,ccp_tmp,cx_tmp, &
           pre1,pre2,pre3,tem1,tem2,tem3,ene1,ene2,ene3, &
           delta_p,den_tmp,tsc1,tsc2
      real(kind=rk) :: ye_dum,xh_dum,cx_dum,den_dum,pre_dum,den_rand, &
           tem_dum,xi_dum,varpi,varpi0,vez1000,pre_rand,flx_rand
      real(kind=rk), dimension(config%qx,2) :: za_tmp
      real(kind=rk) :: eos_self(2), eos_children(2)

      integer(kind=ik) ::    ipos(config%qx)
      real(kind=rk), allocatable :: rand(:,:,:)  ! only used for 2D
      real(kind=rk) :: scrdv, gin,scr,sumass

      integer(kind=ik) :: ngp
      real(kind=rk) :: poisson_self(2), poisson_children(2), &
                       eos3d_self(2), eos3d_children(2)
!-----------------------------------------------------------------------
!     set array boundaries:
!-----------------------------------------------------------------------


      nzn  = config%qx
      nzn1 = nzn + 1
      nzn2 = nzn + 2
      nzn3 = nzn + 3
      nzn4 = nzn + 4
      nzn5 = nzn + 5
      nzn6 = nzn + 6
      nzn7 = nzn + 7
      nzn8 = nzn + 8

      gamma = 4.0_rk/3.0_rk



! interpolate initial model onto hydro-grid

      do k = qz_s, qz_e
      do j = qy_s, qy_e
         open(2,file='nomoto.txt',form='formatted')
         do i = 1, config%qx
            read (2,*) ii,za_tmp(i,1),dentot(i,j,k), &
                 vextot(i,j,k),temtot(i,j,k), &
                 pretot(i,j,k),xnutot(i,j,k,config%qn), &
                 xnutot(i,j,k,1:config%qn-1)
!            read (2,'(26E13.6)') za_tmp(i,1),dentot(i,j,k), &
!                 vextot(i,j,k),temtot(i,j,k), &
!                 pretot(i,j,k),xnutot(i,j,k,config%qn), &
!                 xnutot(i,j,k,1:config%qn-1)
#undef TEST
#ifdef TEST
            if (xnutot(i,j,k,config%qn).lt.0.5_rk) then
               xnutot(i,j,k,config%qn)=0.5_rk+ &
                    (xnutot(i,j,k,config%qn)-0.5_rk)*0.3
            end if
#endif
         end do
         close(2)

! -- attribute the error due to finite precision of input data (~1e-3)
!        to the dominant species, in order to enforce SUM(xi)=1.0
         do i=1, config%qx
            kn=MAXLOC(xnutot(i,j,k,1:config%qn-1), dim=1)
            xnutot(i,j,k,kn)=1.0_rk-(SUM(xnutot(i,j,k,1:config%qn-1))-xnutot(i,j,k,kn))
         enddo

         veytot(:,j,k) = 0.0_rk
         veztot(:,j,k) = 0.0_rk
! rotation

         acxtot(:,j,k) = 0.0_rk
         enetot(:,j,k) = -1.e+33_rk ! will be set in eos3d

         if (j.eq.1) then


!     -------------------------------------------------------------------
!     Solve P(rho,eps,y_e)=P_0 to determine the internal energy eps.
!     -------------------------------------------------------------------
!     Initialize secant method
            pre_tmp(:)=pretot(:,1,1)
            pre1(:)=pre_tmp(:)
            pre2(:)=pre_tmp(:)
!            where(dentot.gt.1.0d10) temtot=1.2d10
            tem1(:)=temtot(:,1,1)*0.9d0
            tem2(:)=temtot(:,1,1)*1.0d0
            !tem1(1:config%qx)=temtot(1:config%qx,1,1)*0.7d0
            !tem2(1:config%qx)=temtot(1:config%qx,1,1)*1.4d0
            tsc1(1:config%qx)=tem1(1:config%qx)
            tsc2(1:config%qx)=tem2(1:config%qx)
            do i=1,config%qx
               write(99,*) i,dentot(i,j,k),temtot(i,j,k), &
                   xnutot(i,j,k,config%qn)
            end do

            call printit_taskX(0,"eos1")
            call eos(dentot(1:config%qx,1,1),tem1(1:config%qx), &
                 xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                 za_tmp(1:config%qx,:),ene1(1:config%qx), &
                 pre1(1:config%qx),gamc_tmp(1:config%qx), &
                 stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                 cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                 ccn_tmp(1:config%qx),eos_self,eos_children,1,0,eos_error)

            call printit_taskX(0,"eos3")
            call eos(dentot(1:config%qx,1,1),tsc1(1:config%qx), &
                xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                za_tmp(1:config%qx,:),ene1(1:config%qx), &
                pre1(1:config%qx),gamc_tmp(1:config%qx), &
                stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                ccn_tmp(1:config%qx),eos_self,eos_children,3,0,eos_error)


            call printit_taskX(0,"eos2")
            call eos(dentot(1:config%qx,1,1),tem2(1:config%qx), &
                 xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                 za_tmp(1:config%qx,:),ene2(1:config%qx), &
                 pre2(1:config%qx),gamc_tmp(1:config%qx), &
                 stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                 cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                 ccn_tmp(1:config%qx),eos_self,eos_children,1,0,eos_error)

            call printit_taskX(0,"eos4")
            call eos(dentot(1:config%qx,1,1),tsc2(1:config%qx), &
                 xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                 za_tmp(1:config%qx,:),ene2(1:config%qx), &
                 pre2(1:config%qx),gamc_tmp(1:config%qx), &
                 stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                 cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                 ccn_tmp(1:config%qx),eos_self,eos_children,3,0,eos_error)

!     Secant method
            delta_p(:)=1.0d0
            call printit_taskX(0,"Bestimme Anfangsmodell...")

            do iit=1,75
           call printit_taskX(0,"Schritt ",iit)

               call eos(dentot(1:config%qx,1,1),tem2(1:config%qx), &
                    xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                    za_tmp(1:config%qx,:),ene2(1:config%qx), &
                    pre2(1:config%qx),gamc_tmp(1:config%qx), &
                    stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                    cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                   ccn_tmp(1:config%qx),eos_self,eos_children,1,0,eos_error)
               tsc2(1:config%qx)=tem2(1:config%qx)
               ene3(1:config%qx)=ene2(1:config%qx)
               tem3(1:config%qx)=tem2(1:config%qx)
               call eos(dentot(1:config%qx,1,1),tsc2(1:config%qx), &
                    xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                    za_tmp(1:config%qx,:),ene2(1:config%qx), &
                    pre2(1:config%qx),gamc_tmp(1:config%qx), &
                    stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                    cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                    ccn_tmp(1:config%qx),eos_self,eos_children,3,0,eos_error)
               if (eos_error) stop 'restrt: construction of outer envelope'


               do i=1,config%qx
!                  if(abs(pre2(i)-pre1(i))/ pre1(i).eq.0.0) delta_p(i)=0.0
                  delta_p(i)=(pre2(i)-pre_tmp(i))/pre_tmp(i)
!                  write(37,*) i,delta_p(i),pre2(i),pre_tmp(i), &
!                       tsc2(i),ene1(i),ene2(i),ene3(i),stotot(i,1,1), &
!                       dentot(i,1,1)
                  if (abs(delta_p(i)).ge.1.0d-2) then ! .and. dentot(i,1,1) .gt. 2.5e10_rk) then
                     ene1(i)=ene3(i)
                     tem2(i)=tem2(i)- &
                         (pre2(i)-pre_tmp(i))/ &
                         (pre2(i)-pre1(i))* &
                         (tem2(i)-tem1(i))

                     if (min(pre2(i),pre1(i)).gt.pre_tmp(i)) &
                          tem2(i)=0.9*min(tem3(i),tem1(i))
                     if (max(pre2(i),pre1(i)).lt.pre_tmp(i)) &
                          tem2(i)=1.1*max(tem3(i),tem1(i))

                     tem1(i)=tem3(i)
                     pre1(i)=pre2(i)
                  else
                     ene1(i)=ene3(i)
                  end if
               end do
               write (*,*) maxval(abs(delta_p(1:config%qx))), &
                    maxloc(abs(delta_p(1:config%qx)))
               if (maxval(abs(delta_p(1:config%qx))).lt.1.0d-4) exit
            end do

            call eos(dentot(1:config%qx,1,1),tem2(1:config%qx), &
                 xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                 za_tmp(1:config%qx,:),ene2(1:config%qx), &
                 pre2(1:config%qx),gamc_tmp(1:config%qx), &
                 stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                 cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                 ccn_tmp(1:config%qx),eos_self,eos_children,3,0,eos_error)
            call eos(dentot(1:config%qx,1,1),tem2(1:config%qx), &
                 xnutot(1:config%qx,1,1,1:config%qn),xhrep_tmp(1:config%qx), &
                 za_tmp(1:config%qx,:),ene2(1:config%qx), &
                 pre2(1:config%qx),gamc_tmp(1:config%qx), &
                 stotot(1:config%qx,1,1),ccu_tmp(1:config%qx), &
                 cce_tmp(1:config%qx),ccp_tmp(1:config%qx), &
                 ccn_tmp(1:config%qx),eos_self,eos_children,3,0,eos_error)

            pretot(1:config%qx,1,1)=pre2(1:config%qx)
            temtot(1:config%qx,1,1)=tem2(1:config%qx)
            enetot(1:config%qx,1,1)=ene2(1:config%qx)/ &
                 dentot(1:config%qx,1,1)+ &
                 0.5d0*vextot(1:config%qx,1,1)**2
            close(93)
            do i=1,config%qx
               write(38,*) i,pre1(i)/pre_tmp(i)-1,pretot(i,1,1), &
                    dentot(i,1,1),temtot(i,1,1),ene1(i)
            end do

            do i=1,config%qx
               write(38,*) i,pretot(i,1,1)/pre_tmp(i)-1, &
                    dentot(i,1,1),temtot(i,1,1), &
!                    xnutot(i,1,1,n_n), &
!                    xnutot(i,1,1,n_p),xnutot(i,1,1,n_he), &
                    stotot(i,1,1)
            end do
         else

           dentot(:,j,1)=dentot(:,1,1)
            temtot(:,j,1)=temtot(:,1,1)
            xnutot(:,j,1,:)=xnutot(:,1,1,:)
            enetot(:,j,1)=enetot(:,1,1)
            pretot(:,j,1)=pretot(:,1,1)

         endif

      enddo
      enddo


      sumass=0.0_rk
      do k=1,config%qz
      do j=qy_s,qy_e
         scr= 4./3.*pc_pi
         do i = 1, config%qx
            scrdv=xzrtot(i)**3-xzltot(i)**3
            sumass=sumass+dentot(i,j,k)* scr*scrdv
         enddo
      enddo
      enddo


!c Initialize parameters for grdvel
!      call massfrac(Rstr)
!      write(*,*) 'Rstr=',Rstr

      call printit_taskX(0," ")
      call printit_taskX(0,"init> Interpolation of model: ")
      call printit_taskX(0," ")
      call printit_taskX(0,"total Mass on grid: ",sumass/pc_ms)

      write(*,'(4x,1x,''   r [cm]   '', &
           &       1x,'' rho [g/cc] '', &
           &       1x,''   T [K]    '', &
           &       1x,'' Y_e [1/by] '', &
           &       1x,'' v_x [cm/s] ''  &
           &     )')

      j = qy_e
      k = qz_s

#ifdef MPI_HYDRO
      if (myproc.eq.0) then
#endif
      open(17,file='ppm.grd',form='formatted')

      do i = 1, config%qx
         write(*,'(1x,i4,5(1x,1pe12.5))') &
             i,xzntot(i),dentot(i,j,k),temtot(i,j,k), &
             xnutot(i,j,k,config%qn),vextot(i,j,k)
         write(17,'(1x,I4,4(1x,1pe13.6))') &
             i,xzntot(i),dentot(i,j,k), &
             temtot(i,j,k),xnutot(i,j,k,config%qn)
      enddo
      call printit_taskX(0," ")

      close(17)
#ifdef MPI_HYDRO
      end if
#endif

! -- calculate remaining state variables:

#ifdef CFC_TRANSPORT
      call cpyare(2)
      call init_cfc_hydro(.false.)
#else
      call cpyare(2)            ! EOS cannot use ***tot-arrays
      areas%nx=config%qx
      areas%ny=config%qy
      areas%nz=config%qz
!      nuc = nutot
      call eos3d (1,eos_error, eos3d_self, eos3d_children)
      abort_if(eos_error)
      call cpyare(0)
#endif

! initialize neutrino quantities:
!     it is possible that there is an output before they are
!     calculated by the neutrino transport

      qyetot(:,:,:,:) = 0.0
      qentot(:,:,:) = 0.0
      qmotot(:,:,:) = 0.0
      qmytot(:,:,:) = 0.0

      fnutot(:,:,:,:) = 0.0
      enutot(:,:,:,:) = 0.0
      dnutot(:,:,:,:) = 0.0
      pnutot(:,:,:,:) = 0.0
      gnutot(:,:,:,:) = 0.0

      eltobs(:) = 0.0
      eltobc(:) = 0.0
      etnobs(:) = 0.0

! reset marker particles:
      ppx(:) = 0.0
      ppy(:) = 0.0
      ppz(:) = 0.0
      pvx(:) = 0.0
      pvy(:) = 0.0
      pvz(:) = 0.0
      pma(:) = 0.0

! set inflow quantities:
      dmdtio(:,:,:) = 0.0

      rhoin  = 1.
      uin    = 0.
      utin   = 0.
      uttin  = 0.
      pin    = 1.
      gin    = 0.0
      tin    = 1.0
      gamein = gamma
      gamcin = gamma
      ein    = pin/((gamein-1.)*rhoin)+0.5*(uin**2+utin**2+uttin**2)
      xnin(:) = 1.0


! define some grid related quantities:

      tsave(:,:)=0.0
      s0_w=0.0
      nop = 0
      igeom = config%igeomx

      ugrtot(:) = 0.0

      dvytot(1) = 1.0
      dvztot(1) = 1.0
      srfint = 1.0
      vlfrac = 4.0 * pc_pi

      call printit_taskX(0,"pns_init> j/dvy")

      if (config%nsdim .ge. 2)  then
         if (config%igeomy .eq. 4)  then
            srfint = cos(yzltot(1)) - cos(yzrtot(config%qy))
!-PNS            srfint = 2. - float(isym)
            vlfrac = 2._rk / srfint
            do j = 1, config%qy
               dvytot(j) = cos(yzltot(j)) - cos(yzrtot(j))
               write(*,'(i3,1x,1pe12.5)') j,dvytot(j)
            end do
            if (config%nsdim .eq. 2)   dvztot(1) = 2.0 * pc_pi
            call printit_taskX(0,"pns_init> srfint = ",srfint)
            call printit_taskX(0,"          vlfrac = ",vlfrac)
         else
            raise_error("read_onemg(): case not implemented")
         end if
      end if

      if (config%nsdim  .eq. 3)  then
         if (config%igeomx .eq. 2  .and. &
              config%igeomy .eq. 4  .and.  config%igeomz .eq. 5  .and. &
              config%bndmny .eq. 4  .and.  config%bndmnz .eq. 4       )  then
            srfint  = 2._rk * config%gridlz * sin(0.5_rk * config%gridly)
            vlfrac  = 4._rk * pc_pi / srfint
            do k = 1, config%qz
               dvztot(k) = zzrtot(k) - zzltot(k)
            end do
         else
           raise_error("read_onemg(); case not implemented")
         end if
      end if

! -- compute potential:

      gpotot(:,:,:)=0.0_rk

      if(config%i_grv .eq. 0) then
         call printit_taskX(0,"init> Newtonian potential")
      else
         call printit_taskX(0,"init> general relativistic potential")
      end if

      ephtot(:) = 1.0_rk
      gamtot(:) = 1.0_rk
      gamold(:) = 1.0_rk

      igrav = 1
      igrav = 0


      write(*,'(80(''=''))')

      END SUBROUTINE read_onemg

!=======================================================================

#endif /* NOMOTO_MODEL */

end module read_progenitor_model

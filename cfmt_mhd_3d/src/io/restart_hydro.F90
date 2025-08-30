#undef DEBUG
module restart

  use precision

  implicit none
  private 
 
#ifndef PROGRAM_remap
 public        &
         restrt,&
         read_restart_filename !, &

#else
 public restrt

#endif

  contains

#ifndef PROGRAM_remap
!>
!> \verbatim
!> This subroutine reads the name of the restart file which is used
!> to CONTINUE a previously computed model 
!>
!>  Author: A. Marek (MPA)
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 663 $
!>   $Date: 2010-01-14 18:41:36 +0100 (Thu, 14 Jan 2010) $
!>
!> \endverbatim
!>
subroutine read_restart_filename
  use precision
  use abort
  use charac, only :  file_number
  use configure
  implicit none

  open  (4, file='./restart/'//trim(config%basenm)//'_next', form='formatted')
  read  (4, '(a)', end=11, err=11)  config%suffix_nr(1:8)
  close (4)

  ! convert suffix into a number
  read ( unit = config%suffix_nr(1:8), fmt = * ) file_number 

  return

11 continue
  write(6,*) 'WARNING:  FILE '// './restart/' // trim(config%basenm)//'_next is empty'
  raise_abort("input():  ")

end subroutine read_restart_filename

#endif /* PROGRAM_remap */

!>
!> \verbatim
!> ----------------------------------------------
!> Author(s) : Wolfgang Keil,Markus Rampp
!>
!> Purpose:      read and write restart files
!> ----------------------------------------------
!>
!>-----------------------------------------------------------------------
!>     arrays to control the different calculation areas:
!>
!>     are_nu               = number of used areas, must be less
!>                            or equal to areas%are_di
!>     areas%are_di               = dimension of dt_are and areas%ix_are
!>     dt_are(*)            = time step of this area
!>     areas%dt_cfl(*)            = cfl time step of this area
!>     ti_are(*)            = time at the end of a time step
!>     areas%ix_are(*, 1 - 9)     = ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz
!>           (*,10)         = isd = sweep direction
!>           (*,11)         = ndt = number of time steps
!>           (*,12 -13)     = config%bndmnx, config%bndmxx
!>     config%ndtmax               = max. number of small time steps
!>                            during one large time step
!>     areas%nhystp               = number of hydro time steps
!>
!>-----------------------------------------------------------------------
!>     Handling of restart files in the MPI version: 
!>
!>     Every MPI process writes its OWN hydro and restart file
!>
!>     The neutrino quantities are written into separate restart
!>     files by each process, because the quantities rip and rim are
!>     never exchanged. In order to facilitate switching the number
!>     of processes, an additional parameter file (.MPI) for MPI is used:
!>     and written in the restart directory

!>     Format: 
!>     2                                        (number of processes)
!>     1 33                                     (Processes operating on
!>     32 64                                    angular bins 1-32, 33-64)
!>     _rank_00000000                           (Suffixes of the MPI
!>     _rank_00000001                            restart files)
!>-----------------------------------------------------------------------
!>
!>+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> VERSION NUMBER SUMMARY
!> these version numbers are known to the code and the remaper
!> keep in mind that both programs have to be able to handle all of these
!> versions !>
!>
!> HYDRO
!> Version number           Explanation
!>      20050101           Standard hydro file, if no version number is
!>                         found in the restart file, than this file format
!>                         (see below) is assumed
!>
!>      20050407           almost standard hydro file, _BUT_ instead of
!>                         the temperature the energy density is input
!>
!>
!> TRANSPORT
!> Version number           Explanation 
!>     20050101           Standard (1D-Eddington factor) transport file,
!>                        if no version number is found in the restart 
!>                        file, than this file format (see below) is assumed
!>
!>     20050408           the quantites rimlag,riplag are seperated in different
!>                        records (needed for 2D-Eddington factor)
!>                        read(82) ralag,xlag
!>                        do i=nystrt,nytra
!>                        read(82) riplag
!>                        read(82) rimlag
!>                        enddo
!>
!>     20091118           Standard (1D-Eddington factor) transport file
!>                        additional quantities kmin,kmax for 3D-Hydro-coupling
!>
!>++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
!> \param ird  restart-modus
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 667 $
!>   $Date: 2010-01-15 17:19:01 +0100 (Fri, 15 Jan 2010) $
!>   
!> \endverbatim
!>
subroutine restrt (ird)

  use precision
  use abort
  use phycon
  use time_rt
  
  use charac
  use intgrs_hy
  use totgrq_hy
 ! use arecon_hy
  use nutrio_hy
  use marker_hy
  use massio_hy
  use revsho_hy
  use gfloat_hy

#ifndef NOTRA  
  use speccut_rt
  use theta_grid_rt

  use transport_files
#ifdef READ_OLD_RESTART_FILES
  use transport_files_old
#ifndef NOHYDRO
  use hydro_files_old
#endif /* NOHYDRO */
#endif
#endif /* NOTRA */
  
  !      use vnew_hy
  !      use mesh_hy
  use totare_hy
  use mo_mpi

#ifndef PROGRAM_remap
#ifndef DEBUG_TIMINGS
  use cputim
#endif
#endif

#ifdef CFC_TRANSPORT
#ifndef PROGRAM_remap
  use gr_hyd_init
#endif
#endif
  use gr_input
  use gr_output


#ifdef PROGRAM_remap
  use dimensions, only : qx_i, qy_i, qz_i, qn_i, q_i
  use oldhydro
#endif

  use hydro_areas_mod
  use configure
  use state
  use print_stdout_mod
 
 IMPLICIT NONE

! LOCAL variables that are not in modules

  integer(kind=ik), intent(in) :: ird

  character(256)               :: suffix_rank, converted_string

  character(40)                :: mpi_rstfil 

  integer(kind=ik)             :: ierr

  integer(kind=ik)             :: i, j, n, ibc_rd, im, km, nys, nyt, nym, &
                                  nzt, iem, ism, ihvers_r 
  integer(kind=ik)             :: bndhnx, bndhxx, bndhny, bndhxy, bndhnz, bndhxz

#ifdef READ_OLD_RESTART_FILES
  real(kind=rk)                :: bndhnx_r, bndhxx_r, bndhny_r, bndhxy_r, &
                                  bndhnz_r, bndhxz_r
#endif

  integer(kind=ik)             :: are_nh, ndtmah, nxhtot, isymh
  real(kind=rk)                :: dt_arh(areas%are_di), dt_cfh(areas%are_di)
  integer(kind=ik)             :: ix_arh(areas%are_di, 13)

  integer(kind=ik)             :: restmass_vers_rd, istat

  real(kind=rk)                :: tim1(2), tim2(2)
  real(kind=rk)                :: poisson_self(2), poisson_children(2)
!-----------------------------------------------------------------------

#ifdef CRAY
  call asnfile ( rstfil,'-F f77 -N ieee_dp',ierr)
  call asnfile ( rstfil_ra,'-F f77 -N ieee_dp',ierr)
#endif


!=========================
  if (ird .eq. 1) then
!=========================

#ifdef PROGRAM_remap
     ! we have to allocate the arrays to read the restart files
     ! after a successful read and storage these will be deallocated

     
  allocate(temtot(qx_i, qy_is:qy_ie, qz_is:qz_ie), enetot(qx_i, qy_is:qy_ie, qz_is:qz_ie), &
           dentot(qx_i, qy_is:qy_ie, qz_is:qz_ie), vextot(qx_i, qy_is:qy_ie, qz_is:qz_ie), &
           veytot(qx_i, qy_is:qy_ie, qz_is:qz_ie), veztot(qx_i, qy_is:qy_ie, qz_is:qz_ie), &         
           xnutot(qx_i, qy_is:qy_ie, qz_is:qz_ie, config%qn), gamtot(0:qx_i),       &
           xzltot(q_i), xzntot(q_i), xzrtot(q_i),              &
           yzltot(q_i), yzntot(q_i), yzrtot(q_i),              &
           zzltot(q_i), zzrtot(q_i), zzntot(q_i),              &
           tsave(4,qy_i),                                      &
           stat=istat )
  if (istat .ne. 0) then
     raise_abort("could not allocate arrays to read hydro restart file")
  endif

  tsave(:,:) = 0._rk

#endif /* PROGRAM_remap */

!-----------------------------------------------------------------------
!     read restart file:
!-----------------------------------------------------------------------
#ifndef NOHYDRO
#ifndef READ_OLD_RESTART_FILES
     call read_hydro(ihvers_r, nxhtot,are_nh, ndtmah, bndhnx,  &
                     bndhxx, bndhny, bndhxy, bndhnz, bndhxz,   &
                     dt_arh, dt_cfh, ix_arh, isymh, restmass_vers_rd)
#else
     call read_hydro_old(ihvers_r, nxhtot,are_nh, ndtmah, bndhnx_r,  &
                         bndhxx_r, bndhny_r, bndhxy_r, bndhnz_r, bndhxz_r,   &
                         dt_arh, dt_cfh, ix_arh, isymh)

     ! for some obscure reason, the boundary flags where stored as reals
     ! in the old restart file format
     bndhnx = int(bndhnx_r, kind=ik) 
     bndhxx = int(bndhxx_r, kind=ik) 
     bndhny = int(bndhny_r, kind=ik) 
     bndhxy = int(bndhxy_r, kind=ik) 
     bndhnz = int(bndhnz_r, kind=ik) 
     bndhxz = int(bndhxz_r, kind=ik)
#endif /*  READ_OLD_RESTART_FILES */
#endif /* NOHYDRO */


#ifdef PROGRAM_remap
     ! take care that a correct restmass version is later written to the 
     ! restart file

     ! if remaper starts from energy density, then the restmass version
     ! is _NOT_ changed, and thus no re-normalisation is necessary
     
     ! if remaper starts from temperature then we set RESTMASS_VERS = 0
     ! and enetot is normalised (by the eos) with that restmass version
     
     if (trim(config%rst_mode) .eq. "energy") then
        config%restmass_version = restmass_vers_rd
     else if (trim(config%rst_mode) .eq. "temp") then
       config% restmass_version = 0
     else
        raise_abort("wrong rst_mode during remap")
     endif
#endif /* PROGRAM_remap */


! temporaer
!         call random_number(rarr)
!         dentot(:,:,:)=dentot(:,:,:)*(1.0+2.0*(rarr(:,:,:)-0.5)*0.05)
!         vextot(:,:,:)=vextot(:,:,:)*(1.0+2.0*(rarr(:,:,:)-0.5)*0.2)
!         temtot(:,:,:)=temtot(:,:,:)*(1.0+2.0*(rarr(:,:,:)-0.5)*0.05)

!-----------------------------------------------------------------------
!     test if parameters of the restart-file are supposed to be used
!     (irstrt = 1) or not (irstrt = 2):
!-----------------------------------------------------------------------
#ifndef PROGRAM_remap

     select case (config%irstrt)

     case (1)
        call printit_taskX(0,"restrt> are_nu/ ndtmax/ dt_are/ dt_cfl/")
        call printit_taskX(0,"        ix_are/ bndmXX taken from restart-file!")

        config%qx  = nxhtot
        areas%are_nu = are_nh
        config%ndtmax = ndtmah
        config%bndmnx = bndhnx
        config%bndmxx = bndhxx
        config%bndmny = bndhny
        config%bndmxy = bndhxy
        config%bndmnz = bndhnz
        config%bndmxz = bndhxz

        config%isym = isymh

        areas%dt_are(1:areas%are_di) = dt_arh(1:areas%are_di)
        areas%dt_cfl(1:areas%are_di) = dt_cfh(1:areas%are_di)
        areas%ix_are(1:areas%are_di,1:13) = ix_arh(1:areas%are_di,1:13)

     case (2)
        write(*,'('' restrt> parameters taken from "ppm.par"!'')')

        do i = 1, areas%are_nu
           do j = 1, are_nh
              if(areas%ix_are(i,1) .ge. ix_arh(j,1) .and. &
                 areas%ix_are(i,1) .le. ix_arh(j,2)) then
                 areas%dt_are(i) = dt_arh(j)
                 areas%dt_cfl(i) = dt_cfh(j)
                 areas%ix_are(i,10) = ix_arh(1,10)
                 areas%ix_are(i,11) = ix_arh(i,11)
                 write(*,'('' restrt> old area: '',i1,'' -> '','' new area: '',i1)') j,i

                 EXIT

              endif
           enddo
        enddo
     case default
        raise_abort("restrt(): restart mode not implemented!")
     end select

     write(*,'('' restrt> are_nu = '',i2,'' ndtmax = '',i4)') areas%are_nu, config%ndtmax
     write(*,'('' restrt> ixfa: '',i4,'' ixfb: '',i4,'' ixfc: '',i4)') &
           areas%ix_are(1,2),areas%ix_are(2,2),areas%ix_are(3,2)
     write(*,'('' restrt> isda: '',i4,'' isdb: '',i4,'' isdc: '',i4)') &
           areas%ix_are(1,10),areas%ix_are(2,10),areas%ix_are(3,10)

!-----------------------------------------------------------------------
!     compute remaining state variables by calling the EOS:
!     Attention: enetot is used in the eos3d-frame, it must contain
!     values even if these are not used!
!     The same is true for pretot, gaetot, gactot, and stotot, because
!     they are used in cpyare!
!-----------------------------------------------------------------------

#ifndef NOHYDRO
     call init_hydro(ihvers_r, restmass_vers_rd)
#endif /* NOHYDRO */

!#ifdef CFC_TRANSPORT
     call init_cfc_hydro(.true.)
     call read_coconut_data
!#endif

#else /* PROGRAM_remap */

     ! copy the read arrays to the "old" remaper arrays

     tem_old = temtot
     ene_old = enetot
     den_old = dentot
     vex_old = vextot
     vey_old = veytot
     vez_old = veztot
     xnu_old = xnutot
     gam_old = gamtot
     xzl_old = xzltot
     xzn_old = xzntot
     xzr_old = xzrtot
     yzl_old = yzltot
     yzn_old = yzntot
     yzr_old = yzrtot
     zzl_old = zzltot
     zzn_old = zzntot
     zzr_old = zzrtot

     tsave_old = tsave

     config%qx  = nxhtot
     areas%are_nu = are_nh
     config%ndtmax = ndtmah
     config%bndmnx = bndhnx
     config%bndmxx = bndhxx
     config%bndmny = bndhny
     config%bndmxy = bndhxy
     config%bndmnz = bndhnz
     config%bndmxz = bndhxz


     areas%dt_are(1:areas%are_di) = dt_arh(1:areas%are_di)
     areas%dt_cfl(1:areas%are_di) = dt_cfh(1:areas%are_di)
     areas%ix_are(1:areas%are_di,1:13) = ix_arh(1:areas%are_di,1:13)

     deallocate(temtot, enetot, dentot, vextot, &
                veytot, veztot, xnutot, gamtot, &
                xzltot, xzntot, xzrtot, yzltot, &
                yzntot, yzrtot, zzltot, zzrtot, &
                zzntot, tsave, stat=istat )
  if (istat .ne. 0) then
     raise_abort("could not deallocate arrays to read hydro restart file")
  endif


#endif /* PROGRAM_remap */

!-----------------------------------------------------------------------
!     initialize neutrino quantities:
!     - the sourceterms are required for the next hydro step 
!     - it is possible that there is an output before they are
!       calculated by the neutrino transport
!-----------------------------------------------------------------------
#ifdef PROGRAM_remap
      if (config%p_ntr .ne. 0) then
         call read_transport
      endif
#else /* PROGRAM_remap */

     r_nusp = 0.0_rk

     i_nusp = 0
     i_nusc = 0

     eltobs = 0.0_rk
     eltobc = 0.0_rk
     etnobs = 0.0_rk
     dmdtio = 0.0_rk

!-----------------------------------------------------------------------
!     compute gravitational potential:     
!-----------------------------------------------------------------------

     igrav = 1
     igrav = 0

     if (use_mpi) then
#ifndef DEBUG_TIMINGS
        call second_v(tim1)
#endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
#ifndef DEBUG_TIMINGS
        call second_v(tim2)

        timer%hydro_comm =timer%hydro_comm +(tim2-tim1)
#endif
     endif

#endif /* PROGRAM_remap */

!=========================
  else ! (ird .eq. 1)
!=========================

!-----------------------------------------------------------------------
!     write restart and ini_next file:  
!-----------------------------------------------------------------------

#ifndef PROGRAM_remap

        nrst = 0
        trst = 0._rk

        if (config%nrstrt .eq. 999999)  return ! no restart files are written


#ifdef CFC_TRANSPORT
!        call write_coconut_data
#endif

#endif /* PROGRAM_remap */

        call write_coconut_data

#ifndef NOHYDRO
        call write_hydro()
#endif

#ifndef NOTRA
        if (config%p_ntr .ne. 0) then
           call write_transport 
        endif

#endif /* NOTRA */

!-----------------------------------------------------------------------

#ifndef PROGRAM_remap

        call printit_taskX(0," ")
        call printit_taskX(0,"********************************************************************************")
        call printit_taskX(0,"Hydro   restart: ",rstfil)
        call printit_taskX(0,"Transp. restart: ",rstfil_ra)
        call printit_taskX(0,"nstep:     ",nstep)
        call printit_taskX(0,"time:  [s] ",time)
        call printit_taskX(0,"dt: [s]    ",hydro%dt)
        call printit_taskX(0,"********************************************************************************")

! only proc zero writes next file and index file
           if (myproc .eq. 0) then

              open  (4, file='./restart/' // trim(config%basenm)//'_next', form='formatted')
              write (4, '(a)')  config%suffix_nr(1:8)
              close (4, status='keep')

! -- write model index file:

              open(13,file = config%index_file,position = filpos,form = 'formatted')
              write(13,'(1x,i9,2(1x,a),2(1x,1pe12.5),1x,i9)') &
                   nstep,trim(rstfil),trim(rstfil_ra),time,hydro%dt,areas%nhystp
              close(13)

           end if


           if (use_mpi) then
#ifndef DEBUG_TIMINGS
              call second_v(tim1)
              call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif
              ! this is necessary to guarantee that the hydrofile is written
              ! completely on return
#ifndef DEBUG_TIMINGS
              call second_v(tim2)
        
              timer%hydro_comm =timer%hydro_comm +(tim2-tim1)
#endif
           endif

!     caught condition config%nrstrt .eq. 99999 at beginning
!         end if 

#endif /* PROGRAM_remap */

     end if ! ird

   end subroutine restrt


#ifndef NOHYDRO

#ifndef READ_OLD_RESTART_FILES
!>
!> \verbatim
!> Read in the hydro-files
!>
!>  Author: B. Mueller
!> \endverbatim
!>
!> \param  ihvers_r      version identifier
!> \param  nxhtot      
!> \param  are_nh
!> \param  ndtmah
!> \param  bndhnx
!> \param  bndhxx
!> \param  bndhny
!> \param  bndhxy
!> \param  bndhnz
!> \param  bndhxz
!> \param  dt_arh
!> \param  dt_cfg
!> \param  ix_arh
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 667 $
!>   $Date: 2010-01-15 17:19:01 +0100 (Fri, 15 Jan 2010) $
!>   
!> \endverbatim
!>
subroutine read_hydro(ihvers_r, nxhtot,are_nh, ndtmah, bndhnx, & 
                      bndhxx, bndhny, bndhxy, bndhnz, bndhxz,  &
                      dt_arh, dt_cfh, ix_arh, isymh, restmass_vers_rd)
!----------------------------------------------------------------------
! Purpose: read hydro restart file
  use precision
  use abort
#ifdef WRITE_BINARY_OUTPUT
  use dataformat_vertex
#endif
 ! use arecon_hy 
  use totare_hy 
  use totgrq_hy, only: ephtot, gamold, gamtot, grptot

  use charac
  use gfloat_hy
  use intgrs_hy
  use marker_hy
  use param_rt
  use revsho_hy

  use mo_mpi


#ifdef PROGRAM_remap
  use dimensions
#endif

  use hydro_areas_mod

  use configure
  use state

  use print_stdout_mod
! LOCAL variables that are not in modules
  IMPLICIT NONE
      
  integer(kind=ik) :: ioerror
  integer(kind=ik), intent(out) :: are_nh, ndtmah, nxhtot, isymh
  integer(kind=ik), intent(out) :: ix_arh(areas%are_di, 13)
  integer(kind=ik), intent(out) :: ihvers_r
  real(kind=rk), intent(out)    :: dt_arh(areas%are_di), dt_cfh(areas%are_di)
  integer(kind=ik), intent(out) :: bndhnx, bndhxx, &
                                   bndhny, bndhxy, &
                                   bndhnz, bndhxz
  
  integer(kind=ik) :: i, j
  integer(kind=ik), dimension (2), parameter :: &
                    ihvers_allowed = (/ 20050101, 20050407 /)


  integer(kind=ik) :: qx_in, qy_in, qz_in, qn_in
  integer(kind=ik) :: qxe_in, qye_in, qze_in, qe_in
  integer(kind=ik) :: qxs_in, qys_in, qzs_in, qs_in
  integer(kind=ik), intent(inout) :: restmass_vers_rd

  type(datafile) :: fd
            

#ifndef PROGRAM_remap
#ifndef MPI_HYDRO
  qxs_in=1
  qxe_in=config%qx
  qys_in=1
  qye_in=config%qy
  qzs_in=1
  qze_in=config%qz
  qs_in=1
  qe_in=config%q

  qx_in = config%qx
  qy_in = config%qy
  qz_in = config%qz
#else
  qxs_in=1
  qxe_in=config%qx
  qys_in=qy_s
  qye_in=qy_e
  qzs_in=qz_s
  qze_in=qz_e
  qs_in=1
  qe_in=config%q
  qx_in = config%qx
  qy_in = config%qy
  qz_in = config%qz
#endif /* MPI_HYDRO */
#else /* PROGRAM_remap */
  ! set to "in" values of restart file to be read
  if (.not.use_mpi) then
    qxs_in=1
    qxe_in=qx_i
    qys_in=1
    qye_in=qy_i
    qzs_in=1
    qze_in=qz_i
    qs_in=1
    qe_in=q_i

    qx_in = qx_i
    qy_in = qy_i
    qz_in = qz_i
  else
    qxs_in=1
    qxe_in=qx_i !config%qx
    qys_in=qy_is
    qye_in=qy_ie
    qzs_in=qz_is
    qze_in=qz_ie
    qs_in=1
    qe_in=q_i
    qx_in = qx_i !config%qx
    qy_in = qy_i !config%qy
    qz_in = qz_i !config%qz
  endif
#endif /* PROGRAM_remap */

#ifdef DEBUG
  write(*,*) "qxs_in = ", qxs_in
  write(*,*) "qxe_in = ", qxe_in
  write(*,*) "qys_in = ", qys_in
  write(*,*) "qye_in = ", qye_in
  write(*,*) "qzs_in = ", qzs_in
  write(*,*) "qze_in = ", qze_in
  write(*,*) "qs_in = ", qs_in
  write(*,*) "qe_in = ", qe_in
  write(*,*) "qx_in = ", qx_in
  write(*,*) "qy_in = ", qy_in
  write(*,*) "qz_in = ", qz_in
#endif


#ifdef PROGRAM_remap
  call printit_taskX(0,"reading ",rstfil)
#endif


  call open_data_file(fd, trim(rstfil))

  call readq(fd, "ihvers", ihvers_r)

  if (config%use_temp_restart) then
     select case (ihvers_r)
     case (20050101)
        call printit_taskX(0,"Reading hydro file (temperature given)")
        config%ihvers = ihvers_r
     case (20050407)
        call printit_taskX(0,"Reading hydro file (temperature given)")
        call printit_taskX(0,"but restart file was saved in mode temperature given!")
        config%ihvers = 20050101
     case default
        raise_abort("read_hydro(): Invalid hydro restart file version!")
     end  select
  else
     select case (ihvers_r)
     case (20050101)
        call printit_taskX(0,"Reading hydro file (energy density given)")
        call printit_taskX(0,"but restart file was saved in mode temperature given!")
        config%ihvers = 20050407
     case (20050407)
        call printit_taskX(0,"Reading hydro file (energy density given)")
        call printit_taskX(0," ") 
        config%ihvers = ihvers_r
     case default
        raise_abort("read_hydro(): Invalid hydro restart file version!")
     end  select
  endif


  print *,'tem'
  call readq(fd, "tem", temtot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
       (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  print *,'ene'
  call readq(fd, "ene", enetot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
       (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  print *,'restmass_version'
  call readq(fd, "restmass_version", restmass_vers_rd)
     
  print *,'den'
  call readq(fd, "den", dentot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
    (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  call readq(fd, "vex", vextot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
    (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  call readq(fd, "vey", veytot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
    (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  call readq(fd, "vez", veztot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in), &
    (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in/), (/1,qx_in, 1,qy_in, 1,qz_in/))

  print *,'xnu'
  call readq(fd, "xnu", xnutot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in,1:config%qn), &
    (/qxs_in,qxe_in, qys_in,qye_in, qzs_in,qze_in,1,config%qn/), (/1,qx_in, 1,qy_in, 1,qz_in, 1,config%qn/))

  print *,'gam'
  call readq(fd, "gam", gamtot(0:qxe_in))

  call readq(fd, "xzl", xzltot(1:qx_in+20), allow_partial_read=.true.)
  call readq(fd, "xzn", xzntot(1:qx_in+20), allow_partial_read=.true.)
  call readq(fd, "xzr", xzrtot(1:qx_in+20), allow_partial_read=.true.)
  call readq(fd, "yzl", yzltot(1:qy_in+20), allow_partial_read=.true.)
  call readq(fd, "yzn", yzntot(1:qy_in+20), allow_partial_read=.true.)
  call readq(fd, "yzr", yzrtot(1:qy_in+20), allow_partial_read=.true.)
  call readq(fd, "zzl", zzltot(1:qz_in+20), allow_partial_read=.true.)
  call readq(fd, "zzn", zzntot(1:qz_in+20), allow_partial_read=.true.)
  call readq(fd, "zzr", zzrtot(1:qz_in+20), allow_partial_read=.true.)

  call readq(fd, "time", time)
  call readq(fd, "dt", hydro%dt)
  call readq(fd, "bndmnx", bndhnx)
  call readq(fd, "bndmxx", bndhxx)
  call readq(fd, "bndmny", bndhny)
  call readq(fd, "bndmxy", bndhxy)
  call readq(fd, "bndmnz", bndhnz)
  call readq(fd, "bndmxz", bndhxz)
  call readq(fd, "rhoin", rhoin)
  call readq(fd, "uin", uin)
  call readq(fd, "utin", utin)
  call readq(fd, "uttin", uttin)
  call readq(fd, "s0_w", s0_w)
  call readq(fd, "pin", pin)
  call readq(fd, "ein", ein)
  call readq(fd, "gridlx", config%gridlx)
  call readq(fd, "gridly", config%gridly)
  call readq(fd, "gridlz", config%gridlz)
  call readq(fd, "rib", config%rib)
  call readq(fd, "pmass", config%pmass)
  call readq(fd, "dt_are", dt_arh)
  call readq(fd, "dt_cfl", dt_cfh)
  call readq(fd, "ix_are", ix_arh)
  call readq(fd, "are_nu", are_nh)
  call readq(fd, "ndtmax", ndtmah)
  call readq(fd, "nhystp", areas%nhystp)
  call readq(fd, "nx", nxhtot)
  call readq(fd, "ny", config%qy)
  call readq(fd, "nz", config%qz)
  qn_in = config%qn
  call readq(fd, "nu", qn_in)
  call readq(fd, "igodu", config%igodu)
  call readq(fd, "nsdim", config%nsdim)
  call readq(fd, "isym", isymh)
  call readq(fd, "nstep", nstep)
  call readq(fd, "igeom", igeom)
  call readq(fd, "igeomx", config%igeomx)
  call readq(fd, "igeomy", config%igeomy)
  call readq(fd, "igeomz", config%igeomz)

  call printit_taskX(0,"Read restmass_version: ",restmass_vers_rd)
  call printit_taskX(0," ")
  call close_data_file(fd)

end subroutine read_hydro
#endif /* READ_OLD_RESTART_FILES */
      
#ifndef PROGRAM_remap
!>
!> \verbatim
!> Initialize hydro related quantities
!>
!>  Author: B. Mueller
!> \endverbatim
!>
!> \param  ihvers_r      version identifier
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 667 $
!>   $Date: 2010-01-15 17:19:01 +0100 (Fri, 15 Jan 2010) $
!>   
!> \endverbatim
!>
subroutine init_hydro(ihvers_r, restmass_vers_rd)
!----------------------------------------------------------------------
! Purpose: init hydro-related quantities
  use precision
  use abort
#ifndef NOTRA
#endif

  use totare_hy 
  use totgrq_hy !, only: ephtot, gamold, gamtot, grptot
  
  use nutrio_hy
  use gfloat_hy
  use intgrs_hy
  
  use mapare_proc
  use mo_mpi
  use eos3d_routine, only : eos3d
  use nucparam
  use cpyare_mod

  use configure
  use hydro_areas_mod
  use print_stdout_mod
  implicit none

! LOCAL variables that are not in modules

  integer(kind=ik), intent(in) :: ihvers_r
  logical                      :: ler
  integer(kind=ik)             :: i, j, k, ierr
  real(kind=rk)                :: totmas, totene, totele
  real(kind=rk)                :: totmas_snd, totmas_rcv, totene_rcv, &
                                  totele_rcv, offs
  real(kind=rk), allocatable :: tem_save(:,:,:)
  integer(kind=ik), intent(in) :: restmass_vers_rd
  real(kind=rk)                :: eos3d_self(2), eos3d_children(2)
! hydro mpi works automatically for arrays with new dimensions 

  if (config%restmass_version .gt. 0 .and. .not.(config%use_temp_restart)) then

     if (restmass_vers_rd .ne. config%restmass_version) then
        
     ! restmass version has been changed during writing
        ! the restart file and the current code !
        call printit_taskX(0," ")
        call printit_taskX(0,"Restmass version is changed!")
        call printit_taskX(0,"Denormalize energy density with old restmass version")
        call printit_taskX(0,"  and renormalize with the new one")    
        call printit_taskX(0," ")  
        ! Denormalize
        ! enetot [erg/g], moffs[erg/g], mbar[erg/g]

        do i=1, size(enetot,dim=1)
           do j=1, size(enetot,dim=2)
              do k=1,size(enetot, dim=3)

                 if (restmass_vers_rd .gt. 0) then

                    if (restmass_vers_rd .eq. 1) then
                       offs = moffs + SUM(xnutot(i,j,k,1:n_he4)*mbar(1:n_he4)) + (1._rk-SUM(xnutot(i,j,k,1:n_he4)))*mbar(n_rep)
                    endif

                    if (restmass_vers_rd .eq. 2) then
                       offs = moffs + SUM(xnutot(i,j,k,1:config%qn-1)*mbar(1:config%qn-1))
                    endif

                    if (restmass_vers_rd .eq. 3) then
                       offs = moffs + SUM(xnutot(i,j,k,1:config%qn)*mbar(1:config%qn))
                    endif

                    enetot(i,j,k) = enetot(i,j,k) + offs


                 endif

                 if (config%restmass_version .eq. 1) then
                    offs = moffs + SUM(xnutot(i,j,k,1:n_he4)*mbar(1:n_he4)) + (1._rk-SUM(xnutot(i,j,k,1:n_he4)))*mbar(n_rep)
                 endif

                 if (config%restmass_version .eq. 2) then
                    offs = moffs + SUM(xnutot(i,j,k,1:config%qn-1)*mbar(1:config%qn-1))
                 endif

                 if (config%restmass_version .eq. 3) then
                    offs = moffs + SUM(xnutot(i,j,k,1:config%qn)*mbar(1:config%qn))
                 endif

                 enetot(i,j,k) = enetot(i,j,k) - offs
              enddo  
           enddo
        enddo
        
     endif
  endif

  pretot(:,:,:) = -1.0e33_rk
  gaetot(:,:,:) = -1.0e33_rk
  gactot(:,:,:) = -1.0e33_rk
  stotot(:,:,:) = -1.0e33_rk
  
  gpotot(:,:,:) = -1.0e33_rk
  grptot(:) = -1.0e33_rk
  ephtot(:) = -1.0e33_rk
  gamold(:) = -1.0e33_rk
      
  acxtot(:,:,:) = 0.0_rk
  qyetot(:,:,:,:) = 0.0_rk
  qentot(:,:,:) = 0.0_rk
  qmotot(:,:,:) = 0.0_rk
  qmytot(:,:,:) = 0.0_rk
  enutot(:,:,:,:) = 0.0_rk
  fnutot(:,:,:,:) = 0.0_rk
  pnutot(:,:,:,:) = 0.0_rk
  dnutot(:,:,:,:) = 0.0_rk
  gnutot(:,:,:,:) = 0.0_rk

!-----------------------------------------------------------------------
!     define some grid related quantities:
!-----------------------------------------------------------------------

  dvytot(1) = 1.0_rk
  dvztot(1) = 1.0_rk
  srfint    = 1.0_rk
  vlfrac    = 4.0_rk * pc_pi

  if (config%igeomx .eq. 2)  then
     dvxtot(1:config%qx)= (xzrtot(1:config%qx)**3 - xzltot(1:config%qx)**3)/3._rk
  else
     raise_abort("init_hydro(): nocase ")
  end if

  if (config%nsdim .ge. 2)  then
     if (config%igeomy .eq. 4)  then
        srfint = cos(yzltot(1)) - cos(yzrtot(config%qy))
        !-PNS            srfint = 2. - float(config%isym)
        vlfrac = 2._rk / srfint
        do j = 1, config%qy
           dvytot(j) = cos(yzltot(j)) - cos(yzrtot(j))
        enddo
        if (config%nsdim .eq. 2)   dvztot(1) = 2.0_rk * pc_pi
     else
        write(6,7001)
7001    format(' case not implemented')
        raise_abort("init_hydro(): nocase")
     end if
  end if

  if (config%nsdim  .eq. 3)  then
     if (config%igeomx .eq. 2  .and. &
          config%igeomy .eq. 4  .and.  config%igeomz .eq. 5 .and. &
          ((config%bndmny .eq. 4 .and. config%bndmnz .eq. 4) .or. &
          (config%bndmny .eq. 1 .and. config%bndmnz .eq. 4) )) then
        
        srfint  = 2.0_rk * config%gridlz * sin(0.5_rk * config%gridly)
        vlfrac  = 4.0_rk * pc_pi / srfint
        do k = 1, config%qz
           dvztot(k) = zzrtot(k) - zzltot(k)
        enddo
     else
        write(6,7001)
        raise_abort("init_hydro(): nocase")
     end if
  end if
  
  do  i = 1, config%qx
     ugrtot(i) = 0.0_rk
  enddo
  
  areas%nx  = config%qx
  areas%ny  = config%qy
  areas%nz  = config%qz
  
#ifndef CFC_TRANSPORT
  call cpyare(2)       ! EOS cannot use ***tot-arrays 

  call printit_taskX(0,"eos3d restart from ",trim(config%rst_mode))

  if (config%use_temp_restart) then
     print *,'Restart from temperature'
     call eos3d(1,ler, eos3d_self, eos3d_children)
  else
     print *,'Restart from energy'
     call eos3d(3,ler, eos3d_self, eos3d_children)
  endif
  abort_if(ler, "init_hydro(): eos3d() failed")

  if (trim(config%rst_mode) .eq. "remap") then
    allocate(tem_save(size(temtot, 1), size(temtot, 2), size(temtot, 3)))
    tem_save = temtot

    call cpyare(0)

    print *,"rst_mode 'remap' set: ", count(abs(1.0_rk - temtot / tem_save) .gt. 0.1_rk), " zones have temperature outliers"

    ! restore temperature and composition in zones where
    ! the deviation in temperature is more than 10%
    temtot = merge(tem_save, temtot, abs(1.0_rk - temtot / tem_save) .gt. 0.1_rk)

    call cpyare(2)
    call eos3d(1,ler, eos3d_self, eos3d_children)
    abort_if(ler, "init_hydro(): eos3d() failed")

    deallocate(tem_save)
    !deallocate(xnu_save)
  endif

  call cpyare(0)

#ifdef FIXCORE
  if (config%irstrt .eq. 2) then
     dentot(1,:,:)=dentot(2,:,:)
     enetot(1,:,:)=enetot(2,:,:)
     xnutot(1,:,:,:)=xnutot(2,:,:,:)
     call cpyare(2)
     call eos3d(3,ler, eos3d_self, eos3d_children)
     if(ler) stop "restrt> eos2 failed"
     call cpyare(0)
  endif
#endif /* FIXCORE */

! if a new arrangement of hydro-grids is enforced by ppm.par the
!  hydro-quantities must be mapped (conservatively)
  if (config%irstrt .eq. 2 .and. config%nsdim .ge. 2) then
     call printit_taskX(0,"restrt> hydro-quantities are mapped")

     totmas=0._rk
     totene=0._rk
     totele=0._rk


     do k = qz_s,qz_e
        do j = qy_s,qy_e
          do i = 1,config%qx
            if (config%nsdim .eq. 3) then
              totmas = totmas + dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
              totene = totene + enetot(i,j,k)*dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
              totele = totele + xnutot(i,j,k,config%qn)*dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
            else if (config%nsdim .eq. 2) then
              totmas = totmas + dentot(i,j,k)*dvxtot(i)*dvytot(j)
              totene = totene + enetot(i,j,k)*dentot(i,j,k)*dvxtot(i)*dvytot(j)
              totele = totele + xnutot(i,j,k,config%qn)*dentot(i,j,k)*dvxtot(i)*dvytot(j)
            endif
          enddo
        enddo
      enddo

      if (use_mpi) then

         ! MPI reduce(SUM)
         !            call MPI_reduce(MPI_IN_PLACE, totmas, 1, 
         !     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         !            call MPI_reduce(MPI_IN_PLACE, totene, 1,  
         !     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         !            call MPI_reduce(MPI_IN_PLACE, totele, 1, 
         !     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         
         totmas_rcv = 0._rk
         call MPI_AllReduce(totmas, totmas_rcv, 1, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, ierr)

         totmas = totmas_rcv


         totene_rcv = 0._rk
         call MPI_AllReduce(totene, totene_rcv, 1,  &
                            MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, ierr)

         totene = totene_rcv


         totele_rcv = 0._rk
         call MPI_AllReduce(totele, totele_rcv, 1,       &
                            MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, ierr)

         totele = totele_rcv


      endif ! use_mpi
         call printit_taskX(0,"old mass, energy, number of charged leptons:", &
              totmas,totene,totele)

     call mapare

     totmas=0._rk
     totene=0._rk
     totele=0._rk


     do k = qz_s,qz_e
        do j = qy_s,qy_e
          do i = 1,config%qx
            if (config%nsdim .eq. 3) then
              totmas = totmas + dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
              totene = totene + enetot(i,j,k)*dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
              totele = totele + xnutot(i,j,k,config%qn)*dentot(i,j,k)*dvxtot(i)*dvytot(j)*dvztot(k)
            else if (config%nsdim .eq. 2) then
              totmas = totmas + dentot(i,j,k)*dvxtot(i)*dvytot(j)
              totene = totene + enetot(i,j,k)*dentot(i,j,k)*dvxtot(i)*dvytot(j)
              totele = totele + xnutot(i,j,k,config%qn)*dentot(i,j,k)*dvxtot(i)*dvytot(j)
            endif
          enddo
        enddo
     enddo

     if (use_mpi) then

        ! MPI reduce(SUM)
        !            call MPI_reduce(MPI_IN_PLACE, totmas, 1, 
        !     &         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        !            call MPI_reduce(MPI_IN_PLACE, totene, 1,  
        !     &         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        !            call MPI_reduce(MPI_IN_PLACE, totele, 1, 
        !     &         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)


        totmas_rcv = totmas
        call MPI_AllReduce((/totmas/), (/totmas_rcv/), 1,       &
                           MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)

        totmas = totmas_rcv


        totene_rcv = totene
        call MPI_AllReduce((/totene/), (/totene_rcv/), 1,       &
                           MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
     
        totene = totene


        totele_rcv = totele
        call MPI_AllReduce((/totele/), (/totele_rcv/), 1,       &
                           MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
        totele = totele_rcv


     endif ! use_mpi

     call printit_taskX(0,"new mass, energy, number of charged leptons:", &
                 totmas,totene,totele)

  endif

#endif /*CFC_TRANSPORT*/

  return

end subroutine init_hydro

#endif /* PROGRAM_remap */
!----------------------------------------------------------------------
!>
!> \verbatim
!> Write hydro restart files
!>
!>  Author: B. Mueller, changed by A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 667 $
!>   $Date: 2010-01-15 17:19:01 +0100 (Fri, 15 Jan 2010) $
!>   
!> \endverbatim
!>
subroutine write_hydro()
  !----------------------------------------------------------------------
  ! Purpose: write hydro restart file
  use precision
  use abort
#ifdef WRITE_BINARY_OUTPUT
  use dataformat_vertex
#endif
  use nucparam, only : name_xnuc
  ! use arecon_hy
  use totare_hy
  use totgrq_hy !, only: ephtot, gamold, gamtot, grptot

  use marker_hy
  use intgrs_hy
  use gfloat_hy
  use param_rt
  use charac
  use revsho_hy

  use mo_mpi

#ifdef PROGRAM_remap
  use dimensions
#endif

  use print_stdout_mod
  use configure
  use hydro_areas_mod
  use state

  IMPLICIT NONE
! LOCAL variables that are not in modules

  integer(kind=ik) :: qx_out, qy_out, qz_out
  integer(kind=ik) :: qxe_out, qye_out, qze_out, qe_out
  integer(kind=ik) :: qe_nqx_out, qe_nqy_out, qe_nqz_out
  integer(kind=ik) :: qxs_out, qys_out, qzs_out, qs_out

  character(len=9) :: x_slice, y_slice, z_slice

  type(datafile) :: fd

#ifndef PROGRAM_remap

  ! local MPI start and end index of radial arrays 
   qxs_out   = 1
   qxe_out   = config%qx
   ! local MPI start and end index of theta arrays
   qys_out   = qy_s
   qye_out   = qy_e
   ! local MPI start and end index of phi arrays
   qzs_out   = qz_s
   qze_out   = qz_e

   ! start and end index of "q" (sweep) arrays
   qs_out    = 1
   qe_out    = config%q

   ! end index of radial sweep (config%qx + 20) arrays
   qe_nqx_out= config%q_nqx
   ! end index of theta sweep (config%qy + 20) arrays
   qe_nqy_out= config%q_nqy
   ! end index of phi sweep (config%qz + 20) arrays
   qe_nqz_out= config%q_nqz

   ! global dimension of radial arrays
   qx_out    = config%qx
   ! global dimension of theta arrays 
   qy_out    = config%qy
   ! global dimension of phi arrays 
   qz_out    = config%qz
#else /* PROGRAM_remap */
   ! set to "out" values of restart file to be written

   ! local MPI start and end index of new radial arrays 
   qxs_out   = 1
   qxe_out   = qx_o
   ! local MPI start and end index of new theta arrays
   qys_out   = qy_os
   qye_out   = qy_oe
   ! local MPI start and end index of new phi arrays
   qzs_out   = qz_os
   qze_out   = qz_oe

   ! start and end index of new "q" sweep arrays
   qs_out=1
   qe_out=q_o

   ! end index of radial sweep (config%qx + 20) arrays
   !   qe_nqx_out= q_o_nqx
   qe_nqx_out = config%qx + 20
   ! end index of theta sweep (config%qy + 20) arrays
   !  qe_nqy_out= q_o_nqy
   qe_nqy_out = config%qy + 20
   ! end index of phi sweep (config%qz + 20) arrays
   !  qe_nqz_out= q_o_nqz
   qe_nqz_out = config%qz + 20

   ! global definition of new radial arrays
   qx_out=qx_o
   ! global definition of new theta arrays
   qy_out=qy_o
   ! global definition of new phi arrays
   qz_out=qz_o
#endif /* PROGRAM_remap */

   ! default setting, ihvers was not reset in model input parameter 
   ! file tra.par

   call printit_taskX(0," Next restart from ",trim(config%rst_mode))

     if (config%use_temp_restart) then
        config%ihvers = 20050101
     else
        config%ihvers = 20050407
     endif

  write(x_slice,"('[0:',i5,']')") qx_out
  write(y_slice,"('[0:',i5,']')") qy_out
  write(z_slice,"('[0:',i5,']')") qz_out

  call create_vertex_file(fd, rstfil, "hydro restart file")
  call writeq(fd, "ihvers", config%ihvers, "type of restart file: 20050101: with  temperature, 20050407: with energy", "enum")

  call writeq(fd, "xzl", xzltot(1:qe_nqx_out), &
    "radius at left rim in hydro grid", "cm", "zone:"//frange(qe_nqx_out))
  call writeq(fd, "xzn", xzntot(1:qe_nqx_out), &
    "radius of zone center", "cm", "zone:"//frange(qe_nqx_out))
  call writeq(fd, "xzr", xzrtot(1:qe_nqx_out), &
    "radius at right rim in hydro grid", "cm", "zone:"//frange(qe_nqx_out))

  call writeq(fd, "yzl", yzltot(1:qe_nqy_out), &
    "inclination (theta) at left rim in hydro grid", "radian", "zone:"//frange(qe_nqy_out))
  call writeq(fd, "yzn", yzntot(1:qe_nqy_out), &
    "inclination (theta) of zone center", "radian", "zone:"//frange(qe_nqy_out))
  call writeq(fd, "yzr", yzrtot(1:qe_nqy_out), &
    "inclination (theta) at right rim in hydro grid", "radian", "zone:"//frange(qe_nqy_out))

  call writeq(fd, "zzl", zzltot(1:qe_nqz_out), &
     "azimuth (phi) at left rim in hydro grid", "radian", "zone:"//frange(qe_nqz_out))
  call writeq(fd, "zzn", zzntot(1:qe_nqz_out), &
    "azimuth (phi) of zone center", "radian", "zone:"//frange(qe_nqz_out))
  call writeq(fd, "zzr", zzrtot(1:qe_nqz_out), &
     "azimuth (phi) at right rim in hydro grid", "radian", "zone:"//frange(qe_nqz_out))

  call writeq(fd, "tem", temtot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    "temperature", "K", "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "ene", enetot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    "specific energy", "erg/g", "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "restmass_version", config%restmass_version, &
                        "Energy normalization: different version for the subtraction of rest masses\n&
                        &from the energy used in PPM:\n&
                        &  0: uses energy defined as in EoS\n&
                        &  1: subtracts from EoS energy the baryon rest masses, assuming\n&
                        &     that heavy elements have the mass of fe56\n&
                        &  2: subtracts from EoS energy the baryon rest masses\n&
                        &  3: subtracts from EoS energy the baryon and unpaired electron\n&
                        &     rest masses. Caution! This version violates energy!", "1")

  call writeq(fd, "den", dentot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    'density', 'g/cm**3', "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "vex", vextot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    'velocity in radial direction', 'cm/s', "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "vey", veytot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    'velocity in theta direction', 'cm/s', "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "vez", veztot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out), &
    'velocity in phi direction', 'cm/s', "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out/), (/1,qx_out, 1,qy_out, 1,qz_out/))

  call writeq(fd, "xnu", xnutot(qxs_out:qxe_out,qys_out:qye_out,qzs_out:qze_out,1:config%qn), &
    'composition (mass fractions)', '1', "radius:xzn" // x_slice, "theta:yzn" // y_slice, "phi:zzn" // z_slice, &
    'species:' // stringlist(name_xnuc), &
    (/qxs_out,qxe_out, qys_out,qye_out, qzs_out,qze_out, 1,config%qn/), (/1,qx_out, 1,qy_out, 1,qz_out, 1,config%qn/))

  write(x_slice,"('[0:',i5,']')") qx_out + 1
#ifndef PROGRAM_remap
  call writeq(fd, "gam", gamold(0:qxe_out), "Gamma factor (grav. potential)", "1", "radius:xzn" // x_slice)
#else /* PROGRAM_remap */
  call writeq(fd, "gam", gamtot(0:qxe_out), "Gamma factor (grav. potential)", "1", "radius:xzn" // x_slice)
#endif /* PROGRAM_remap */

  call writeq(fd, "time", time, "physical time", "s")
  call writeq(fd, "dt", hydro%dt, "timestep?", "s")
  call writeq(fd, "bndmnx", config%bndmnx, "Inner boundary condition in x-direction", "enum")
  call writeq(fd, "bndmxx", config%bndmxx, "Outer boundary condition in x-direction", "enum")
  call writeq(fd, "bndmny", config%bndmny, "Inner boundary condition in y-direction", "enum")
  call writeq(fd, "bndmxy", config%bndmxy, "Outer boundary condition in y-direction", "enum")
  call writeq(fd, "bndmnz", config%bndmnz, "Inner boundary condition in z-direction", "enum")
  call writeq(fd, "bndmxz", config%bndmxz, "Outer boundary condition in z-direction", "enum")
  call writeq(fd, "rhoin", rhoin, "density value for flow in boundary", "g/cm**3")
  call writeq(fd, "uin", uin, 'velocity value for flow in boundary', 'cm/s')
  call writeq(fd, "utin", utin, 'some boundary flow', 'unknownunit')
  call writeq(fd, "uttin", uttin, 'some boundary flow', 'unknownunit')
  call writeq(fd, "s0_w", s0_w, 's0_w - what is this?', 'unknownunit')
  call writeq(fd, "pin", pin, "pressure value for flow in boundary", "erg/cm**3")
  call writeq(fd, "ein", ein, "energy value for flow in boundary", "erg/g")
  call writeq(fd, "gridlx", config%gridlx, "Radius of outer grid boundary", "cm")
  call writeq(fd, "gridly", config%gridly, "Grid length in theta direction", "pi")
  call writeq(fd, "gridlz", config%gridlz, "Grid length in phi direction", "pi")
  call writeq(fd, "rib", config%rib, "Radius of inner grid boundary", "cm")
  call writeq(fd, "pmass", config%pmass, "inner point mass", "unknownuint")
  call writeq(fd, "dt_are", areas%dt_are, "time steps of calculation areas", "s", "zone:"//frange(areas%are_di))
  call writeq(fd, "dt_cfl", areas%dt_cfl, "cfl time steps of calculation areas", "s", "zone:"//frange(areas%are_di))
  call writeq(fd, "ix_are", areas%ix_are, "calculation area structure", "enum", "zone:"//frange(areas%are_di), "zone:"//frange(13))
  call writeq(fd, "are_nu", areas%are_nu, "number of areas", "1")
  call writeq(fd, "ndtmax", config%ndtmax, "max. number of small time steps during one large time step", "1")
  call writeq(fd, "nhystp", areas%nhystp, "number of hydro timesteps", "1")
  call writeq(fd, "nx", config%qx, 'number of grid points in x-direction', '1')
  call writeq(fd, "ny", config%qy, 'number of grid points in y-direction', '1')
  call writeq(fd, "nz", config%qz, 'number of grid points in z-direction', '1')
  call writeq(fd, "nu", config%qn, 'number of nuclear species', '1')
  call writeq(fd, "igodu", config%igodu, "unknown flag", "unknownunit")
  call writeq(fd, "itstp", config%itstp, 'print timestep information, if ITSTP.ne.0', 'enum')
  call writeq(fd, "nsdim", config%nsdim, "Dimensionality of the problem", "1")
  call writeq(fd, "isym", config%isym, "Assume equatorial symmetry, if isym == 1", "enum")
  call writeq(fd, "nstep", nstep, "hydro timestep number", "1")
  call writeq(fd, "igeom", igeom, "Geometry type for first dimension in case of 1D (?)", "enum")
  call writeq(fd, "igeomx", config%igeomx, "Geometry type for first dimension", "enum")
  call writeq(fd, "igeomy", config%igeomy, "Geometry type for second dimension", "enum")
  call writeq(fd, "igeomz", config%igeomz, "Geometry type for third dimension", "enum")

  call close_data_file(fd)

end subroutine write_hydro
#endif /* NOHYDRO */
end module restart


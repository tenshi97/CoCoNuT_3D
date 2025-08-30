module hydro_files_old

  public read_hydro_old
  private

  contains

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
!> \param  r0fin
!> \param  etlos
!> \param  tmlos
!> \param  are_nh
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 667 $
!>   $Date: 2010-01-15 17:19:01 +0100 (Fri, 15 Jan 2010) $
!>   
!> \endverbatim
!>
subroutine read_hydro_old(ihvers_r, nxhtot,are_nh, ndtmah, bndhnx, & 
                          bndhxx, bndhny, bndhxy, bndhnz, bndhxz,  &
                          dt_arh, dt_cfh, ix_arh, isymh)
                        
!----------------------------------------------------------------------
! Purpose: read hydro restart file
  use precision
  use abort

  !  use arecon_hy 
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
! LOCAL variables that are not in modules
  IMPLICIT NONE
      
  integer(kind=ik) :: ioerror
  integer(kind=ik), intent(out) :: are_nh, ndtmah, nxhtot, isymh
  integer(kind=ik), intent(out) :: ix_arh(areas%are_di, 13)
  integer(kind=ik), intent(out) :: ihvers_r
  real(kind=rk)                 :: r0fin, etlos, tmlos    ! no more used
  real(kind=rk) , intent(out)   :: dt_arh(areas%are_di), &
                                   dt_cfh(areas%are_di)
  real(kind=rk) , intent(out)   :: bndhnx, bndhxx, bndhny, bndhxy, &
                                   bndhnz, bndhxz
  
  integer(kind=ik)              :: itsth, i, j
  integer(kind=ik), dimension (2), parameter :: &
                    ihvers_allowed = (/ 20050101, 20050407 /)


  integer(kind=ik)              :: qxe_in, qye_in, qze_in, qe_in
  integer(kind=ik)              :: qe_nqx_in, qe_nqy_in, qe_nqz_in
  integer(kind=ik)              :: qxs_in, qys_in, qzs_in, qs_in
            

#ifndef PROGRAM_remap
  ! local MPI start and end index of radial arrays 
  qxs_in = 1
  qxe_in = config%qx
  ! local MPI start and end index of theta arrays
  qys_in = qy_s
  qye_in = qy_e
  ! local MPI start and end index of phi arrays
  qzs_in = qz_s
  qze_in = qz_e

  ! start and end index of "q" (sweep) arrays
  qs_in  = 1
  qe_in  = config%q

  ! end index of radial sweep (qx + 20) arrays
  qe_nqx_in= config%q_nqx
  ! end index of theta sweep (qy + 20) arrays
  qe_nqy_in= config%q_nqy
  ! end index of phi sweep (qz + 20) arrays
  qe_nqz_in= config%q_nqz

#else /* PROGRAM_remap */
  ! local MPI start and end index of radial arrays 
  qxs_in = 1
  qxe_in = qx_a
  ! local MPI start and end index of theta arrays
  qys_in = 1
  qye_in = qy_a
  ! local MPI start and end index of phi arrays
  qzs_in = 1
  qze_in = qz_a

  ! start and end index of "q" (sweep) arrays
  qs_in  = 1
  qe_in  = q_a


  ! end index of radial sweep (qx + 20) arrays
  qe_nqx_in= q_a_nqx
  ! end index of theta sweep (qy + 20) arrays
  qe_nqy_in= q_a_nqy
  ! end index of phi sweep (qz + 20) arrays
  qe_nqz_in= q_a_nqz

  write (*,*) q_a_nqx, q_a_nqy, q_a_nqz
  raise_abort("restart_hydro_old")

#endif /* PROGRAM_remap */

#ifdef PROGRAM_remap
  call printit_taskX(0,"reading ",rstfil)
#endif 

  open (2,file = rstfil,form = 'unformatted')
  read (2,iostat=ioerror) ihvers_r
  if (ihvers_r .ne. 20050407 .and. ihvers_r .ne. 20050101) then
     ! old version where ihvers was not stored in the file
     ! select correct value by define-statement
     close(2)

     open (2,file = rstfil,form = 'unformatted')
#ifdef READ_TEM
     ihvers_r=20050101
#else
     ihvers_r=20050407
#endif
  end if
  
  select case (ihvers_r)
  case (20050101)
     write (*,*) 'Reading old version of hydro file (temperature given)'

     read (2, iostat=ioerror)                                                &
               dentot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               vextot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veytot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veztot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               gamtot(0:qxe_in),                                             &
               temtot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               xnutot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in,1:config%qn),       &
               xzltot(1:qe_nqx_in), xzntot(1:qe_nqx_in), xzrtot(1:qe_nqx_in),&
               yzltot(1:qe_nqy_in), yzntot(1:qe_nqy_in), yzrtot(1:qe_nqy_in),&
               zzltot(1:qe_nqz_in), zzntot(1:qe_nqz_in), zzrtot(1:qe_nqz_in),&
               tsave(1:4,1:qye_in), r0fin, etlos, tmlos,                     &
               time, transport%dt, bndhnx, bndhxx, bndhny,                   &
               bndhxy, bndhnz, bndhxz, rhoin,                                &
               uin, utin, uttin, s0_w, pin, ein,                             &
               config%gridlx, config%gridly, config%gridlz, config%rib, config%pmass,                           &
               ppx, ppy, ppz, pvx, pvy, pvz, nop,                            &
               dt_arh, dt_cfh, ix_arh, are_nh,                               &
               ndtmah,areas%nhystp,nxhtot , config%qy ,                                &
               config%qz , config%qn , config%igodu , itsth ,                               &
               config%nsdim, isymh, nstep , igeom ,                                 &
               config%igeomx, config%igeomy, config%igeomz

     if (IOERROR /= 0) then
        raise_abort("read_hydro(): Error on reading hydro file with tem!")
     endif

  case (20050407)
! read energy instead of temperature
     write (*,*)  &
          'Reading hydro file version with energy density given.'

     read (2,iostat=ioerror)                                                 &
               dentot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               vextot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veytot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veztot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &  
               gamtot(0:qxe_in),                                             &
               enetot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &  
               xnutot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in,1:config%qn),       &
               xzltot(1:qe_nqx_in), xzntot(1:qe_nqx_in), xzrtot(1:qe_nqx_in),&
               yzltot(1:qe_nqy_in), yzntot(1:qe_nqy_in), yzrtot(1:qe_nqy_in),&
               zzltot(1:qe_nqz_in), zzntot(1:qe_nqz_in), zzrtot(1:qe_nqz_in),&
               tsave(1:4,1:qye_in), r0fin, etlos, tmlos,                     &
               time, transport%dt, bndhnx, bndhxx, bndhny,                             &
               bndhxy, bndhnz, bndhxz, rhoin,                                &
               uin, utin, uttin, s0_w, pin, ein,                             &
               config%gridlx, config%gridly, config%gridlz, config%rib, config%pmass,                           &
               ppx, ppy, ppz, pvx, pvy, pvz, nop,                            &
               dt_arh, dt_cfh, ix_arh, are_nh,                               &
               ndtmah,areas%nhystp,nxhtot , config%qy ,                                &
               config%qz , config%qn , config%igodu , itsth ,                               &
               config%nsdim, isymh, nstep , igeom ,                                 &
               config%igeomx, config%igeomy, config%igeomz
     
     if (IOERROR /= 0) then
        raise_abort("read_hydro(): Error on reading hydro file with ene!")
     endif




  case default           !old restart files
     write (*,*) 'Reading old version of hydro file (temperature given).'
     ihvers_r = 20050101
     close(2)
     open (2,file = rstfil,form = 'unformatted')
     

     read (2,iostat=ioerror) &
               dentot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               vextot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veytot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               veztot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               gamtot(0:qxe_in),                                             &
               temtot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in),            &
               xnutot(qxs_in:qxe_in,qys_in:qye_in,qzs_in:qze_in,1:config%qn),       &
               xzltot(1:qe_nqx_in), xzntot(1:qe_nqx_in), xzrtot(1:qe_nqx_in),&
               yzltot(1:qe_nqy_in), yzntot(1:qe_nqy_in), yzrtot(1:qe_nqy_in),&
               zzltot(1:qe_nqz_in), zzntot(1:qe_nqz_in), zzrtot(1:qe_nqz_in),&
               tsave(1:4,1:qye_in), r0fin, etlos, tmlos,                     &
               time, transport%dt, bndhnx, bndhxx, bndhny,                             &
               bndhxy, bndhnz, bndhxz, rhoin,                                &
               uin, utin, uttin, s0_w, pin, ein,                             &
               config%gridlx, config%gridly, config%gridlz, config%rib, config%pmass,                    &
               ppx, ppy, ppz, pvx, pvy, pvz, nop,                            &
               dt_arh, dt_cfh, ix_arh, are_nh,                               &
               ndtmah,areas%nhystp,nxhtot , config%qy ,                                & 
               config%qz , config%qn , config%igodu , itsth ,                               &
               config%nsdim, isymh, nstep , igeom ,                                 &
               config%igeomx, config%igeomy, config%igeomz
     if (IOERROR /= 0) then
        raise_abort("read_hydro(): Error on reading (default) hydro file!")
     endif
     
  end select
  
  close (2)

  if (config%ihvers .eq. 0) then
     config%ihvers=ihvers_r
  else
     if (.not.ANY(ihvers_allowed(:) .eq. config%ihvers)) then
        raise_abort("read_hydro(): This ihvers is not implemented!")
     endif
  endif

end subroutine read_hydro_old

end module


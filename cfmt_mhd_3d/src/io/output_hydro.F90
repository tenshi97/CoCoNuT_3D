#ifdef WRITE_BINARY_OUTPUT
module output_hydro

  private
  public :: set_filenames
  public ::  open_files, close_files, write_output_files

  contains

!>
!> \verbatim
!> This subroutine creates the new filenames for the output and restart
!> files. In MPI-mode each task creates its own filename
!> Be carefull: the first filename after (re)-starting a model is set due
!> to historicall reasons in input.F90
!>
!>  Author: A. Marek, MPA
!> \endverbatim
!>
!> \param irstrt  restart modus
!> \param irst_ra transport on / off ?
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine set_filenames(irstrt,irst_ra)

  use precision
  use abort
  use charac
  use mo_mpi
  use configure

  implicit none
  integer(kind=ik),intent(in) :: irstrt,irst_ra
! the vallue of irstrt is only meaningful when this subroutine is called
! for the first time after starting the code!

  logical :: restart_mode=.false.

  select case(irstrt)

     case(0)
        ! starting a new model
        restart_mode=.false.
        config%suffix_nr(1:8)="00000000"

     case(1:2)
        ! simulation runs already
        restart_mode=.true.

        ! suffix / suffix_nr are set in subroutine read_restart_filename
     case(99)
        file_number=file_number+1
   end select

        write(unit=config%suffix_nr, fmt='(i8.8)') file_number
        if (use_mpi) then
           write(unit=myproc_character, fmt='(i8.8)') myproc
        endif

#ifdef CFC_TRANSPORT
        rstfil_cfc= './restart/' // trim(config%basenm) // '.' // 'q' // config%suffix_nr(1:8)
        outfil_cfc= './output/' // trim(config%basenm) // '.' // 'm' // config%suffix_nr(1:8)
#endif

        rstfil    = './restart/' // trim(config%basenm) // '.' // 'r' // config%suffix_nr(1:8)
        rstfil_ra = './restart/' // trim(config%basenm) // '.' // 's' // config%suffix_nr(1:8)

        outfil    = './output/' // trim(config%basenm) // '.' // 'o' // config%suffix_nr(1:8)
        outfil_ra = './output/' // trim(config%basenm) // '.' // 'p' // config%suffix_nr(1:8)

        if (restart_mode) then
           write(6,*) "restarting from ", rstfil
        endif

  return
end subroutine set_filenames

!>
!> \verbatim
!> This subroutine opens ieee-output files and writes headers
!>
!>  Author: M. Rampp, MPA
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine open_files
  use precision
  use dataformat_vertex, only : create_vertex_file, hfile, pfile
!  use param_rt, only : p_ntr
  use charac

! LOCAL variables that are not in modules

  use configure

  implicit none

  call create_vertex_file(hfile, outfil, "hydro output file")
  if (config%p_ntr .ne. 0) then
    call create_vertex_file(pfile, outfil_ra, "transport output file")
  endif

end subroutine open_files

!>
!> \verbatim
!> This subroutine closes ieee-output files 
!>
!>  Author: M. Rampp, MPA
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine close_files
  use precision
!  use param_rt, only : p_ntr
  use dataformat_vertex, only : close_data_file, hfile, pfile

  use configure

  implicit none

  call close_data_file(hfile)
  if (config%p_ntr .ne. 0) call close_data_file(pfile)

  return
end subroutine close_files


!>
!> \verbatim
!> This subroutine writes (some) output to standard io and all
!> of it to the output files as definied in subroutine open_filenames
!> 
!> Notes:  additional output of gpot! (wfk)
!>         additional output of neutrino-quants
!>         output sequencing: more than one
!>         timestep in a file 
!>         in MPI-mode each task writes its own output file
!>
!>  Author: W. Keil, and M. Rampp, MPA
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine write_output_files

  use precision

  use abort

  use phycon
  use intgrs_hy, only : nstep, nout1
  use charac, only : outfil, outfil_ra, filpos
  use totare_hy, only :                      xzntot, yzntot, zzntot, &
                                             xzltot, yzltot, zzltot, &
                                             xzrtot, yzrtot, zzrtot, &
#ifdef CFC_TRANSPORT2
                                             wltot, &
#endif
                                             dentot, vextot, veytot, &
                                             veztot, pretot, enetot, &
                                             temtot, stotot, gactot, &
                                             xnutot, gpotot, dvxtot, &
                                             dvytot, dvztot, ishtot
  use totgrq_hy, only : gamtot, ephtot, tgmtot
  use nutrio_hy, only : cpotot, qyetot, qentot, qmotot, dnutot, enutot, &
                        fnutot, pnutot, eltobs, etnobs 
!  use param_rt , only : p_ntr
  use gfloat_hy, only : time, tout1, vlfrac! , pmass
  use nucparam, only : name_xnuc
  use dataformat_vertex

  use mo_mpi
#ifndef NOTRA
  use write_transport, only : write_transport_output
#endif

#ifdef CFC_TRANSPORT
#ifdef CFC_MHD
  use hydro_primitives_cfc, ONLY: psi, b_1, b_2, b_3
#endif
  use mesh_coarsening
  use metric_cfc, ONLY: sqrt_gamma, alpha, phi, beta_up_1, beta_up_2, beta_up_3
#endif /* CFC_TRANSPORT */
  use grids

  use configure
  use state
  use print_stdout_mod
  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik) :: i,j,k
  real(kind=rk) :: gamhlp(config%qx), ephhlp(config%qx) 
  real(kind=rk) :: dm(config%qx), tm(config%qx), tgm(config%qx)
  real(kind=rk) :: tm_rcv(config%qx), tgm_rcv(config%qx)

  real(kind=rk) :: dsurf, dvol, dmas
#ifdef CFC_TRANSPORT
  real(kind=rk), allocatable :: xcart(:,:,:),ycart(:,:,:),zcart(:,:,:) !Cartesian coordinates
  real(kind=rk) :: bx(0:config%qx), by (0:config%qy), bz(0:config%qz)
#endif /* CFC_TRANSPORT */
  integer(kind=ik) :: ierr, istat

  ! hydro domain, start, end
  !> \todo these should be replaced when the mpi-cleanup branch is incorporated
  integer(kind=ik) :: process
  integer(kind=ik) :: ys, ye, ys0, ys_gpo
  integer(kind=ik) :: zs, ze, zs0
  real(kind=rk)     :: mapra2hyd_self(2), mapra2hyd_children(2)
#ifdef MPI_HYDRO
  process = myproc
  ys = qy_s
  ye = qy_e
  zs = qz_s
  ze = qz_e
#else /* no MPI */
  process = 0
  ys = 1
  ye = config%qy
  zs = 1
  ze = config%qz
#endif

! this is not clear to me (amarek). Check!!
!  if (myproc .gt. 0) then
!     nout1 = 0              !important: if nout1 and tout1 are not
!     tout1 = 0.             !reset, I/O attempts after the next
!     return                 !timestep will cause the code to crash
!  end if


!
#ifndef NOTRA
  if (config%p_ntr .ne. 0) call map_ra2hyd(.false., mapra2hyd_self, mapra2hyd_children)
#endif /* NOTRA */


#ifdef CFC_TRANSPORT
#if defined(COARSENING) && (DIMN==3)
  call prolongation_output
#endif  
  call init_tot_arrays(qy_s,qy_e,qz_s,qz_e,-1)
#endif /* CFC_TRANSPORT */

  dm(:)  = 0.0_rk
  tgm(:) = 0.0_rk 
  dvxtot(1:config%qx) = (xzrtot(1:config%qx)**3 - xzltot(1:config%qx)**3)  &
                        / 3._rk


  do k = qz_s, qz_e
    do j = qy_s, qy_e

           dsurf = vlfrac * dvztot(k) * dvytot(j)
           do i = 1, config%qx
!c          if(xzltot(i) .gt. rcut  .and.  nstep .le. 1) goto 12
              dvol  = dsurf * dvxtot(i)

#ifndef CFC_TRANSPORT2

              dmas  = dentot(i,j,k) * dvol

#else /* CFC_TRANSPORT2 */

              dmas  = dentot(i,j,k) * dvol &
                     * sqrt_gamma(i,j,k) * wltot(i,j,k)
#endif /* CFC_TRANSPORT2 */

              dm(i) = dm(i) + dmas

           enddo
        enddo
     enddo

     tm(1) = config%pmass + dm(1)
     do i = 2, config%qx
        tm(i) = tm(i-1) + dm(i)
     enddo
!
     do i = 1, config%qx
        tm(i)  = tm(i) / pc_ms
        tgm(i) = tgmtot(i) / pc_ms
     enddo
     !

! sum tm and tgm up
     tm(config%qx+1:config%qx)  = 0._rk
     tgm(config%qx+1:config%qx) = 0._rk


     if (use_mpi) then
! sum tm and tgm up
        call MPI_AllReduce(tm(1:config%qx), tm_rcv(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        tm(1:config%qx) = tm_rcv(1:config%qx)


        call MPI_AllReduce(tgm(1:config%qx), tgm_rcv(1:config%qx), config%qx, MPI_DOUBLE_PRECISION,&
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        tgm(1:config%qx) = tgm_rcv(1:config%qx)
     endif

!     
!----------------------------------------------------------------------
!     report output to the protocol:
!----------------------------------------------------------------------

     
     call printit_taskX(0,"********************************************************************************")
     call printit_taskX(0,"Hydro   output: ",outfil)
     call printit_taskX(0,"Transp. output: ",outfil_ra)
     call printit_taskX(0,"nstep:    ",nstep)
     call printit_taskX(0,"time: [s] ",time)
     call printit_taskX(0,"dt: [s]   ",hydro%dt)   
     call printit_taskX(0,"********************************************************************************")
  call printit_taskX(0," ")

!----------------------------------------------------------------------
!     3D-IEEE-Interface:
!----------------------------------------------------------------------

     gamhlp(1:config%qx)           = gamtot(1:config%qx) 
     gamhlp(config%qx+1:config%qx) = 0._rk
     ephhlp(1:config%qx)           = ephtot(1:config%qx) 
     ephhlp(config%qx+1:config%qx) = 0._rk

     !> \todo these should be replaced when the mpi-cleanup branch is incorporated

     call write_timestep_marker(hfile)

     call writeq(hfile, "xzn", xzntot(1:config%qx), "radius of zone center", "cm", "zone:"//frange(config%qx))
     call writeq(hfile, "yzn", yzntot(1:config%qy), "inclination (theta) of zone center", "radian", "zone:"//frange(config%qy))
     call writeq(hfile, "zzn", zzntot(1:config%qz), "azimuth (phi) of zone center", "radian", "zone:"//frange(config%qz))

     call writeq(hfile, "bx", (/xzltot(1), xzrtot(1:config%qx)/), &
       "radius of zone boundaries", &
       "cm", &
       "zone:"//frange(0,config%qx))

     call writeq(hfile, "by", (/yzltot(1), yzrtot(1:config%qy)/), &
       "inclination (theta) of zone boundaries", &
       "radian", &
       "zone:"//frange(0,config%qy))

     call writeq(hfile, "bz", (/zzltot(1), zzrtot(1:config%qz)/), &
       "azimuth (phi) of zone boundaries", &
       "radian", &
       "zone:"//frange(0,config%qz))

     call writeq(hfile, "xzl", xzltot(1:config%qx), "radius at left rim in hydro grid", "cm", &
       "zone:"//frange(config%qx))

     call writeq(hfile, "yzl", yzltot(1:config%qy), "inclination (theta) at left rim in hydro grid", "radian", &
       "zone:"//frange(config%qy))

     call writeq(hfile, "zzl", zzltot(1:config%qz), "azimuth (phi) at left rim in hydro grid", "radian", &
       "zone:"//frange(config%qz))
         
     call writeq(hfile, "xzr", xzrtot(1:config%qx), "radius at right rim in hydro grid", "cm", &
       "zone:"//frange(config%qx))

     call writeq(hfile, "yzr", yzrtot(1:config%qy), "inclination (theta) at right rim in hydro grid", "radian", &
       "zone:"//frange(config%qy))

     call writeq(hfile, "zzr", zzrtot(1:config%qz), "azimuth (phi) at right rim in hydro grid", "radian", &
       "zone:"//frange(config%qz))

#ifdef CFC_TRANSPORT
     ! Cartesian grid (for VISIT)
     if (ys .eq. 1) then
        ys0=ys-1
     else
        ys0=ys
     end if
     if (zs .eq. 1) then
        zs0=zs-1
     else
        zs0=zs
     end if

     bx(1:config%qx)=xzrtot(1:config%qx)
     by(0)=xzltot(1)
     by(1:config%qy)=yzrtot(1:config%qy)
     by(0)=yzltot(1)
     bz(1:config%qz)=zzrtot(1:config%qz)
     bz(0)=zzltot(1)

     allocate(xcart(0:config%qx,ys0:ye,zs0:ze),ycart(0:config%qx,ys0:ye,zs0:ze),zcart(0:config%qx,ys0:ye,zs0:ze),stat=istat)

     do k=zs0,ze
        do j=ys0,ye
           do i=0,config%qx
              xcart(i,j,k)=bx(i)!*sin(by(j))*cos(bz(k))
              ycart(i,j,k)=bx(i)!*sin(by(j))*sin(bz(k))
              zcart(i,j,k)=bx(i)!*cos(by(j))
           end do
        end do
     end do

     call writeq(hfile, "xcart", xcart, "X coordinate (Cartesian)", "cm", &
       "radius:bx", "theta:by", "phi:bz", (/0,config%qx, ys0,ye, zs0,ze/), (/0,config%qx,0,config%qy,0,config%qz/))
     call writeq(hfile, "ycart", ycart, "Y coordinate (Cartesian)", "cm", &
       "radius:bx", "theta:by", "phi:bz", (/0,config%qx, ys0,ye, zs0,ze/), (/0,config%qx,0,config%qy,0,config%qz/))
     call writeq(hfile, "zcart", zcart, "Z coordinate (Cartesian)", "cm", &
       "radius:bx", "theta:by", "phi:bz", (/0,config%qx, ys0,ye, zs0,ze/), (/0,config%qx,0,config%qy,0,config%qz/))

     deallocate(xcart,ycart,zcart,stat=istat)
#endif


     call writeq(hfile, "tm", tm, "Enclosed baryonic mass", "solarmass", "radius:xzn")
     call writeq(hfile, "tgm", tgm, "Enclosed gravitational mass", "solarmass", "radius:xzn")
     call writeq(hfile, "eph", ephhlp, "Lapse function (grav. potential)", "1", "radius:xzn")
     call writeq(hfile, "gam", gamhlp, "Gamma factor (grav. potential)", "1", "radius:xzn")

     call writeq(hfile, "den", dentot, "density", "g/cm**3", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "vex", vextot, "velocity in radial direction", "cm/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "vey", veytot, "velocity in theta direction", "cm/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "vez", veztot, "velocity in phi direction", "cm/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "pre", pretot, "pressure", "erg/cm**3", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "ene", enetot, "specific energy", "erg/g", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "restmass_version", int(config%restmass_version, kind=ik), &
                        "Energy normalization: different version for the subtraction of rest masses\n&
                        &from the energy used in PPM:\n&
                        &  0: uses energy defined as in EoS\n&
                        &  1: subtracts from EoS energy the baryon rest masses, assuming\n&
                        &     that heavy elements have the mass of fe56\n&
                        &  2: subtracts from EoS energy the baryon rest masses\n&
                        &  3: subtracts from EoS energy the baryon and unpaired electron\n&
                        &     rest masses. Caution! This version violates energy!", "1")

     call writeq(hfile, "tem", temtot, "temperature", "K", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "sto", stotot, "entropy per baryon", "kboltzmann", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "gac", gactot, "adiabatic index", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "xnu", xnutot, "composition (mass fractions)","1", &
       "radius:xzn", &
       "theta:yzn", &
       "phi:zzn", &
       "species:" // stringlist(name_xnuc), &
       (/1,config%qx, ys,ye, zs,ze, 1,config%qn/), (/1,config%qx,1,config%qy,1,config%qz,1,config%qn/))

     call writeq(hfile, "cpo", cpotot, "chemical potentials", "MeV", &
       "radius:xzn", "theta:yzn", "phi:zzn", 'species:["neutrino", "electron", "neutron", "proton"]', &
       (/1,config%qx, ys,ye, zs,ze, 1,4/), (/1,config%qx,1,config%qy,1,config%qz,1,4/))


     ! take care of border values
     if (ys .eq. 1) then
       ys_gpo = ys - 1
     else
       ys_gpo = ys
     endif
     call writeq(hfile, "gpo", gpotot(0:config%qx, ys_gpo:ye, zs:ze), "grav. pot. in hydro grid", "erg/g", &
       "radius:bx", "theta:by", "phi:zzn", &
       (/0,config%qx, ys_gpo,ye, zs,ze /), (/0,config%qx, 0,config%qy, 1,config%qz /))

     call writeq(hfile, "qye", qyetot, "Quell-Ye", "g/cm**3/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", 'list:["total", "part1", "part2", "part3", "part4"]', &
       (/1,config%qx, ys,ye, zs,ze, 1,5/), (/1,config%qx,1,config%qy,1,config%qz,1,5/))

#ifdef NEW_QMO
     call writeq(hfile, "qen", (qentot+qmotot*vextot), "Quell-ene", "erg/cm**3/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
#else
     call writeq(hfile, "qen", qentot, "Quell-ene", "erg/cm**3/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
#endif  /* NEW_QMO */

     call writeq(hfile, "qmo", qmotot, "Quell-mom", "erg/cm**4", & 
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))

     call writeq(hfile, "dnu", dnutot, "neutrino number density (lab)", "1/cm**3", &
       "radius:xzn", "theta:yzn", "phi:zzn", "species:"//specieslist(), &
       (/1,config%qx, ys,ye, zs,ze, 1,config%isma/), (/1,config%qx,1,config%qy,1,config%qz,1,config%isma/))

     call writeq(hfile, "enu", enutot, "neutrino energy density (lab)", "erg/cm**3", &
       "radius:xzn", "theta:yzn", "phi:zzn", "species:"//specieslist(), &
       (/1,config%qx, ys,ye, zs,ze, 1,config%isma/), (/1,config%qx,1,config%qy,1,config%qz,1,config%isma/))

     call writeq(hfile, "fnu", fnutot, "neutrino energy flux (lab)", "erg/cm^2/s", &
       "radius:xzn", "theta:yzn", "phi:zzn", "species:"//specieslist(), &
       (/1,config%qx, ys,ye, zs,ze, 1,config%isma/), (/1,config%qx,1,config%qy,1,config%qz,1,config%isma/))
     
     call writeq(hfile, "pnu", pnutot(1:config%qx,ys:ye,zs:ze,:), "neutrino pressure", "erg/cm**3", &
       "radius:xzn", "theta:yzn", "phi:zzn", "species:"//specieslist(), &
       (/1,config%qx, ys,ye, zs,ze, 1,config%isma/), (/1,config%qx,1,config%qy,1,config%qz,1,config%isma/))

     call writeq(hfile, "nstep", nstep, "hydro timestep number", "1")
     call writeq(hfile, "dt", hydro%dt, "hydro timestep interval", "s")
     call writeq(hfile, "ish", ishtot(1:config%qx, ys:ye, zs:ze), &
       "shock zones", &
       "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", &
       (/1,config%qx,ys,ye,zs,ze/), &
       (/1,config%qx,1,config%qy,1,config%qz/))


#ifdef CFC_TRANSPORT2
     call writeq(hfile, "phi", phi(1:config%qx,ys:ye,zs:ze), "conformal factor", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "alpha", alpha(1:config%qx,ys:ye,zs:ze), "lapse function", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "beta1", beta_up_1(1:config%qx,ys:ye,zs:ze), "shift vector (radial)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "beta2", beta_up_2(1:config%qx,ys:ye,zs:ze), "shift vector (meridional)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "beta3", beta_up_3(1:config%qx,ys:ye,zs:ze), "shift vector (zonal)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
#endif /* CFC_TRANSPORT2 */
#ifdef CFC_MHD
     call writeq(hfile, "b_1", pc_geob*b_1(1:config%qx,ys:ye,zs:ze), "B field (radial)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "b_2", pc_geob*b_2(1:config%qx,ys:ye,zs:ze), "B field (meridional)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "b_3", pc_geob*b_3(1:config%qx,ys:ye,zs:ze), "B field (zonal)", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
     call writeq(hfile, "psi", pc_geob*psi(1:config%qx,ys:ye,zs:ze), "psi field", "1", &
       "radius:xzn", "theta:yzn", "phi:zzn", (/1,config%qx, ys,ye, zs,ze/), (/1,config%qx,1,config%qy,1,config%qz/))
#endif /* CFC_MHD */

     ! Flush file in order to get readable data in case of a crash later
     call flush_data_file(hfile)

#ifndef NOTRA
     if (config%p_ntr .ne. 0) call write_transport_output
#endif

!----------------------------------------------------------------------
!     write model index file:
!----------------------------------------------------------------------

     if (myproc .eq. 0) then
        ! only task 0 should write the index file
        open(13,file = config%index_file,position = filpos,form = 'formatted')
        write(13,'(1x,i9,2(1x,a),2(1x,1pe12.5),3x,a)')  &
             nstep,trim(outfil),trim(outfil_ra),time,hydro%dt,"r8"
        close(13)
     endif

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!                                                                       
      nout1 = 0
      tout1 = 0.0_rk

      if(tm(config%qx) .le. 0.0_rk) raise_abort("output_hydro(): tm(config%qx) <= 0 !!!")

   return
end subroutine write_output_files

end module output_hydro
#endif

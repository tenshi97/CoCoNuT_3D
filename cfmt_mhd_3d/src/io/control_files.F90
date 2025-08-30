module timeinfo

  private
  public :: &
#ifndef NOTRA
            print_timestep_evolution, &
#endif
            write_timestep_information

  contains

#ifndef NOTRA
!>
!> \verbatim
!> This subroutine writes (some) output to standard output that
!> gives a summary of the computed timestep
!> 
!> Printed are : the energy budget
!>               the lepton number budget
!>
!>               the evolution file model_name.evo
!>               the neutrino file model_name.ntr
!>               the energy file model_name.erg

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
subroutine print_timestep_evolution

  use precision

#ifndef NOTRA
  use neutrinotypes
#endif

  use param_rt
  use phycon
  use intgrs_hy
!  use arecon_hy
  use totare_hy
  use totgrq_hy
  use nutrio_hy
  use massio_hy
  use charac
!  use revsho_hy
  use gfloat_hy
  use nucparam

  use ioflx
  use cnsv
!  use specfun
!      use vnew_hy
!      use mesh_hy

#ifdef CFC_TRANSPORT
  use size_cfc, ONLY: n_s,n_e,o_s,o_e
  use parameters_cfc, ONLY: m_rest,m_rest_ini,m_grav_1, &
                            m_grav_1_ini
  use metric_cfc, ONLY: sqrt_gamma
  use grids, only : init_tot_arrays, map_nutra2cfc
  use cons_check
#endif

  use mo_mpi
#ifndef DEBUG_TIMINGS
  use cputim
#endif
  use print_stdout_mod

  use configure
  use hydro_areas_mod
  use state
  implicit none
! LOCAL variables that are not in modules
      ! "global" variables for hydro MPI
  real(kind=rk) :: sbuf(2*config%isma+config%qx+34), rbuf(2*config%isma+config%qx+34)
  integer(kind=ik) :: ierr, offset

  integer(kind=ik) :: is,j,i,k
  real(kind=rk) :: dgrav,tnut,tnut_e,slnu,tent,edis,dmas
  real(kind=rk) :: gamfac,dmdtout,elto_new,dmint,ebindm
  real(kind=rk) :: egrav,ekinz,ebindp,ekinx,ekiny,tten
  real(kind=rk) :: etot,elto,snnu
  real(kind=rk) :: dedtout,fe,dydtout,tten_new
  real(kind=rk) :: dvol,dsurf,tnmass,phi,eint,ekin,ytot
!-----------------------------------------------------------------------
!     one dimensional arrays of the IEEE -interface:
!-----------------------------------------------------------------------

  real(kind=rk)::dm(config%qx), tm(config%qx),dmdt(config%qx), &
                 dedt(config%qx),                              & ! ph(qx) ve(qx), 
                 dydt(config%qx)

  real(kind=rk):: tens(6),tnus(6),sln(6),snn(6)

  real(kind=rk):: uave (config%qx), pave (config%qx), dave (config%qx),  &
                  tave(config%qx), vave (config%qx), wave (config%qx),   &
                  save (config%qx), yave(config%qx),                     &
                  u2ave(config%qx), v2ave(config%qx), w2ave(config%qx),  &
                  qeave(config%qx), qyave(config%qx)

  real(kind=rk)::          xp, amoz, moiz

!-----------------------------------------------------------------------


  real(kind=rk)::  eloss_h,nloss_h,eloss_n,nloss_n,eloss,nloss, &
                   nges_new,nges_new2,eges_new,eges_new2
  real(kind=rk)::  enoffs

  integer(kind=ik) ,parameter :: modout=1

  real(kind=rk) :: tim1(2), tim2(2)
!c----------------------------------------------------------------------
!c     print only every modout'th model:
!c----------------------------------------------------------------------
!
  if(nstep .gt. 2 .and. mod(nstep,modout) .ne. 0 .and. areas%ix_are(1,11) .lt. 100) return   


#ifdef CFC_TRANSPORT
  call init_tot_arrays(n_s,n_e,o_s,o_e,-1)
  call map_nutra2cfc
  call do_conservation_check
#endif

!-------------------------s---------------------------------------------
!     calculate energies and mass:
!----------------------------------------------------------------------

  tten   = 0._rk
  sln(:) = 0._rk          ! surface neutrino-luminosity for any species
  snn(:) = 0._rk          ! surface number-flux for any species
  
  tens(:) = 0._rk          ! total neutrino energy for a species
  tnus(:) = 0._rk          ! total neutrino number for a species
  
  elto   = 0._rk          ! total electron number
  sumqen = 0._rk          ! total neutrino energy source
  sumqye = 0._rk          ! total neutrino number source
  sumqmo = 0._rk          ! total neutrino momentum source
  
  etot   = 0._rk
  ekinx  = 0._rk
  ekiny  = 0._rk
  ekinz  = 0._rk
  egrav  = 0._rk
  ebindp = 0._rk
  ebindm = 0._rk
  amoz   = 0._rk                    ! z-component: angular momentum
  moiz   = 0._rk                    ! z-component: moment of inertia
  xp     = 1._rk                    ! x-projection


  do  i = 1, size(dm)
     dm (i)    = 0._rk
     dvxtot(i) = (xzrtot(i)**3 - xzltot(i)**3) / 3._rk
  enddo


  do k = qz_s, qz_e
   do j = qy_s, qy_e

     dsurf = vlfrac * dvztot(k) * dvytot(j)

     if(config%igeomy .eq. 4) then
        xp = sin(yzntot(j))
     endif
     do is=1,config%isma
        sln(is) = sln(is) + dsurf * fnutot(config%qx,j,k,is) * &
                 xzntot(config%qx)**2
        snn(is) = snn(is) + dsurf * gnutot(config%qx,j,k,is) * &
                 xzntot(config%qx)**2
! --test only
        snn(is) = snn(is)/pc_cl
     enddo

     do i = 1, size(dm)
        dvol  = dsurf * dvxtot(i)
        dmas  = &
#ifdef CFC_TRANSPORT2
                sqrt_gamma(i,j,k) * wltot(i,j,k) * &
#endif
                dentot(i,j,k) * dvol

        dm(i) = dm(i) + dmas


        if (config%restmass_version .eq. 1) then
!#if LOW_EOS_NSE_FLAG != 22
           enoffs= moffs + SUM(xnutot(i,j,k,1:n_he4)*mbar(1:n_he4)) &
                      +(1._rk-SUM(xnutot(i,j,k,1:n_he4)))*mbar(n_rep)
!#else
!         enoffs = moffs + SUM(xnutot(i,j,k,1:n_he4)*mbar(1:n_he4)) 
!     &                - SUM(xnutot(i,j,k,n_d:n_he3)*mbar(n_d:n_he3))
!     &           + (1.e0_rk-SUM(xnutot(i,j,k,1:n_he4)))*mbar(n_rep)
!     &           - (1.e0_rk-SUM(xnutot(i,j,k,n_d:n_he3)))*mbar(n_rep)
!#endif
        endif

        if (config%restmass_version .eq. 2) then 
           enoffs= moffs + SUM(xnutot(i,j,k,1:config%qn-1)*mbar(1:config%qn-1))
        endif

        if (config%restmass_version .eq. 3) then 
           enoffs= moffs + SUM(xnutot(i,j,k,1:config%qn)*mbar(1:config%qn))
        endif

        if (config%restmass_version .gt. 0) then 
           tten  = tten  +       dmas * (enetot(i,j,k) + enoffs)
        else
           tten  = tten  +       dmas * enetot(i,j,k)
        endif

        ekinx = ekinx + 0.5_rk * dmas * vextot(i,j,k)**2
        ekiny = ekiny + 0.5_rk * dmas * veytot(i,j,k)**2
        ekinz = ekinz + 0.5_rk * dmas * veztot(i,j,k)**2
        do is=1,config%isma
           tens(is)= tens(is) +       dvol * enutot(i,j,k,is)
           tnus(is)= tnus(is) +       dvol * dnutot(i,j,k,is)
        enddo

        sumqen = sumqen + dvol * qentot(i,j,k)
        sumqye = sumqye + dvol * qyetot(i,j,k,1)
        sumqmo = sumqmo + dvol * qmotot(i,j,k)

#ifdef CFC_TRANSPORT2
        dgrav = 0.0_rk
#else
        dgrav = 0.125_rk * dmas * ( &
                          + gpotot(i  ,j  ,k) &
                          + gpotot(i  ,j-1,k) &
                          + gpotot(i-1,j  ,k) &
                          + gpotot(i-1,j-1,k))
#endif /* CFC_TRANSPORT2 */

        egrav = egrav + dgrav
        edis  = enetot(i,j,k) * dmas  +  dgrav

        elto  = elto + xnutot(i,j,k,config%qn)*dmas

        if (edis .gt. 0._rk)  then
           ebindp = ebindp + edis
        else
           ebindm = ebindm + edis
        end if

        moiz = moiz + dmas * xzntot(i)**2 * xp**2
        amoz = amoz + dmas * xzntot(i) * xp * veztot(i,j,k)

     enddo
   enddo
  enddo


  if (use_mpi) then
     ! MPI reduce
     !    sln(1:config%isma),  snn(1:config%isma), dm(1:qx)
     !    tten, ekinx, ekiny, ekinz, tens, tnus, sumqen, sumqye,
     !    sumqmo, egrav
     !    elto, ebindp, ebindm, moiz, amoz,
     !      sumdegr sumdenu sumdynu sumioe sumion sumioe_n sumion_n
     do i = 1, config%isma
        sbuf(     i) = sln(i)
        sbuf(config%isma+i) = snn(i)
     end do
     do i = 1, config%qx
        sbuf(2*config%isma+i) = dm(i)
     end do
     offset = 2*config%isma+config%qx
     do i = 1, 6
        sbuf(offset+  i) = tens(i)
        sbuf(offset+6+i) = tnus(i)
     end do
     offset = offset + 12
     
     sbuf(offset+ 1) = tten
     sbuf(offset+ 2) = ekinx
     sbuf(offset+ 3) = ekiny 
     sbuf(offset+ 4) = ekinz
     sbuf(offset+ 5) = sumioe_n
     sbuf(offset+ 6) = sumion_n
     sbuf(offset+ 7) = sumqen 
     sbuf(offset+ 8) = sumqye
     sbuf(offset+ 9) = sumqmo
     sbuf(offset+10) = egrav
     sbuf(offset+11) = elto
     sbuf(offset+12) = ebindp
     sbuf(offset+13) = ebindm
     sbuf(offset+14) = moiz
     sbuf(offset+15) = amoz
     sbuf(offset+16) = sumdegr
     sbuf(offset+17) = sumdenu
     sbuf(offset+18) = sumdynu
     sbuf(offset+19) = sumioe
     sbuf(offset+20) = sumion
#ifndef DEBUG_TIMINGS  
     call second_v(tim1)
#endif
     call MPI_allreduce(sbuf, rbuf, offset+20, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
#ifndef DEBUG_TIMINGS
     call second_v(tim2)
  
     timer%transp_comm = timer%transp_comm +(tim2-tim1)
#endif
  

     do i = 1, config%isma
        sln(i) = rbuf(     i)
        snn(i) = rbuf(config%isma+i)
     end do
     do i = 1, config%qx
        dm(i) = rbuf(2*config%isma+i)
     end do
     offset = 2*config%isma+config%qx
     do i = 1, 6
        tens(i) = rbuf(offset+  i)
        tnus(i) = rbuf(offset+6+i)
     end do
     offset = offset + 12
     
     tten     = rbuf(offset+ 1)
     ekinx    = rbuf(offset+ 2)
     ekiny    = rbuf(offset+ 3)
     ekinz    = rbuf(offset+ 4)
     sumioe_n = rbuf(offset+ 5)
     sumion_n = rbuf(offset+ 6)
     sumqen   = rbuf(offset+ 7) 
     sumqye   = rbuf(offset+ 8) 
     sumqmo   = rbuf(offset+ 9)
     egrav    = rbuf(offset+10)
     elto     = rbuf(offset+11)
     ebindp   = rbuf(offset+12)
     ebindm   = rbuf(offset+13)
     moiz     = rbuf(offset+14)
     amoz     = rbuf(offset+15)
     sumdegr  = rbuf(offset+16)
     sumdenu  = rbuf(offset+17)
     sumdynu  = rbuf(offset+18)
     sumioe   = rbuf(offset+19)
     sumion   = rbuf(offset+20)

  endif ! use_mpi


     tent=SUM( tens(:) )        ! total energy of all neutrino-species
     tnut=SUM( tnus(:) )       ! total number of all neutrino-species

     tnut_e=0._rk
     do is=1,config%isma
        tnut_e=tnut_e + siglep(is)*tnus(is)
     enddo

     slnu=SUM( sln(:) )        ! surface-luminosity of neutrinos
     snnu=SUM( snn(:) )        ! surface-number-loss of neutrinos

     tm(1) = dm(1)
     
     do  i = 2, size(tm)
        tm(i) = tm(i-1) + dm(i)
     enddo

     elto = elto/pc_mb
     if (config%p_nbk .gt. 0) then
        etot   = tten + egrav + tent
        !         ytot   = rin + elto
        do is=1,config%isma
           ytot   = elto + siglep(is)*tnus(is)
        enddo
     else
        etot   = tten + egrav
        ytot   = elto
     endif
      
     ekin   = ekinx + ekiny + ekinz
     eint   = tten - ekin

     if (myproc.eq.0) then
!      etot   = etot*1.D-50
!      ekin   = ekin*1.D-50
!      ekinx = ekinx*1.D-50
!      ekiny = ekiny*1.D-50
!      ekinz = ekinz*1.D-50
!      eint   = eint*1.D-50
!      tten  = tten*1.D-50
!      egrav = egrav*1.D-50
!      ebindp= ebindp*1.D-50
!      ebindm= ebindm*1.D-50
     phi    = gpotot(config%qx,1,1)
     moiz   = moiz*1.e-45_rk
     amoz   = amoz*1.e-48_rk
!      ytot   = ytot*1.D-50
     endif ! myproc .eq. 0

     tnmass = tm(config%qx)/pc_ms

!c--  Sourceterms
!      sumdegr = sumdegr*1.D-50
!      sumdenu = sumdenu*1.D-50
!      sumdynu = sumdynu/pc_mb*1.D-50

!c-- neutrino-quantities
!      tent   = tent*1.D-50
!      tnut   = tnut*1.D-50
!      slnu   = slnu*1.D-50
!      snnu   = snnu*1.D-50
!c      spnu   = spnu*1.D-30
!      sumqen = sumqen*1.D-50
!      sumqye = sumqye/pc_mb*1.D-50
!      sumqmo = sumqmo*1.D-50

     if (myproc.eq.0) then
     do i = 1, size(uave)
        dmint    = 0._rk
        uave (i) = 0._rk
        vave (i) = 0._rk
        wave (i) = 0._rk
        u2ave(i) = 0._rk
        v2ave(i) = 0._rk
        w2ave(i) = 0._rk
        dave (i) = 0._rk
        pave (i) = 0._rk
        tave (i) = 0._rk
        save (i) = 0._rk
        yave (i) = 0._rk
        gamfac   = 2._rk/(gamtot(i-1) + gamtot(i)) ! rest mass
        dmdt(i)  = 0._rk
        dedt(i)  = 0._rk
        dydt(i)  = 0._rk
        
        qeave    = 0._rk
        qyave    = 0._rk

        uave (i) = vextot(i,1,1)
        u2ave(i) = vextot(i,1,1)**2
        dave (i) = dentot(i,1,1)
        pave (i) = pretot(i,1,1)
        tave (i) = temtot(i,1,1)
        save (i) = stotot(i,1,1)
        yave (i) = xnutot(i,1,1,config%qn)

        qeave(i) = qentot(i,1,1)
        qyave(i) = qyetot(i,1,1,1)
     enddo

!c Fluxes of mass, energy and lepton-number at outer boundary
     dmdtout   = 4._rk * pc_pi * dmdtio(2,1,1)/pc_ms
     dedtout   = 4._rk * pc_pi * dmdtio(4,1,1)*1.e-50_rk
     dydtout   = 4._rk * pc_pi * dmdtio(6,1,1)/pc_mb*1.e-50_rk
     endif ! myproc .eq. 0

     eloss=sumioe - (sumdegr+sumdenu)
     nloss=sumion/pc_mb - sumdynu/pc_mb

     eloss_h=sumioe - sumdegr
     nloss_h=sumion/pc_mb 

     eloss_n=sumioe_n
     nloss_n=sumion_n

     tten_new=tten+eloss
     elto_new=elto+nloss

#ifdef FCNC_CALC
!      eloss_h = eloss_h + sumdenu_n
!      nloss_h = nloss_h + sumdynu_n
#endif


#ifdef CFC_TRANSPORT
!----CoCoNuT version
     eges_new=0.5_rk*m_grav_1*pc_geoe  +eloss_h+eloss_n
     eges_new2=(tten+dtot(2)*pc_meverg)  +eloss_h+eloss_n
     nges_new=(elto+tnut_e)+nloss_h+nloss_n
     nges_new2=(elto+dtot(1))+nloss_h+nloss_n
     tten_new=tten+eloss
#else
!----PROMETHEUS version
     eges_new=(tten+tent)  +eloss_h+eloss_n
     eges_new2=(tten+dtot(2)*pc_meverg)  +eloss_h+eloss_n
     nges_new=(elto+tnut_e)+nloss_h+nloss_n
     nges_new2=(elto+dtot(1))+nloss_h+nloss_n
#endif /* CFC_TRANSPORT2 */

!      write (*,*) 'e',eges_new,eges_old,sumdenu_n
!      write (*,*) 'n',nges_new,nges_old,sumdynu_n
     
     call printit_taskX(0," ")
     
     call printit_taskX(0,"ECONS -->",(tten_new-tten_old)/tten_new, &
                                      (eges_new-eges_old)/eges_new/ &
                                       real((modout),kind=rk))

#ifdef FCNC_CALC

     call printit_taskX(0,"ECONS -->",(tten_new-tten_old)/tten_new, &
                          (eges_new+ sumdenu_n-eges_old)/eges_new/  &
                          real(modout,kind=rk) )
#endif

     call printit_taskX(0,"NCONS -->",(elto_new-elto_old)/elto_new, &
                                      (nges_new-nges_old)/nges_new/ &
                                       real(modout,kind=rk) )

#ifdef FCNC_CALC
     call printit_taskX(0,"NCONS -->",(elto_new-elto_old)/elto_new, &
                           (nges_new+ sumdynu_n-nges_old)/nges_new/ &
                            real(modout,kind=rk) )
#endif
     call printit_taskX(0,"ECrel -->",(tten_new-tten_old)/tten_new/hydro%dt, &
                                        (eges_new-eges_old)/eges_new/  &
                                        real(modout,kind=rk)/hydro%dt )


     call printit_taskX(0,"NCrel -->",(elto_new-elto_old)/elto_new/hydro%dt, &
                                         (nges_new-nges_old)/nges_new/ &
                                        real(modout,kind=rk)/hydro%dt)

     tten_old=tten
     elto_old=elto
     eges_old=tten+tent
     eges_old2=tten+dtot(2)*pc_meverg
     nges_old=elto+tnut_e
     nges_old2=elto+dtot(1)

#ifdef CFC_TRANSPORT
     eges_old=0.5_rk*m_grav_1*pc_geoe
     m_grav_1_ini=m_grav_1_ini-2*pc_egeo*(eloss_h+eloss_n)
     write(*,1201)  time,pc_geoe*(m_grav_1-m_grav_1_ini)
1201 format(1x, 'Accuracy of ADM energy conservation: ', 2(1X,1pe16.9), ' erg')
#endif

!----------------------------------------------------------------------
!     write output file of the time evolution: hydro quantities
!----------------------------------------------------------------------

      if (myproc.eq.0) then
         open(16,file = config%evolution_file,position = filpos,form = 'formatted')

         if(nstep .eq. 1) then
        write(16,'(''#'','' step:  '', 1x,'' time [s]   '', &
     &             1x,''E_tot[^50e]'',1x,''E_kin[^50e]'',   &
     &             1x,''E_int[^50e]'',1x,''E_ene[^50e]'',   &
     &             1x,''E_gra[^50e]'',1x,''            '',  &
     &             1x,''QE_gr[^50e/s]]'',1x,''N_ele[^50  ]'',&
     &             1x,''M_b [M_sol]'',1x,''ro_c[g/cc] '',   &
     &             1x,'' r_s [cm]   '',1x,''p_s [erg/cc]'', &
     &             1x,''ro_s [g/cc] '',1x,'' v_s [cm/s] '', &
     &             1x,'' u_gs[cm/s] '',1x,''dmdt[M_so/s]'', &
     &             1x,''dedt[^50e/s]'',1x,''dydt[    /s]'', &
     &             1x,''Phi_s[^50 e]'')')
         write(16,'( &
     &             1x,''iuma'',1x,''umax[cm/s]'',1x,''iqemi'', &
     &             1x,''1/dte [/s]'',1x,''iqymi'',1x,''1/dty [/s]'', &
     &             1x,''E_nut[^50e]'',1x,''N_enu[^50]'',             &
     &             1x,''L_nus[^50e/s]'',1x,''P_nus[^30e/s]'',        &
     &             1x,''QE_nu[^50e/s]'',1x,''QY_nu[^50 /s]'',        &
     &             1x,''QP_nu[^50e/cm/s]'')')
         write(16,'(''#'',243(''=''))')


      endif
      fe=1e-50_rk

      write(16,'(1x,i8,18(1x,1pe11.4))') &
               nstep,time,etot*fe,ekin*fe,eint*fe,tten*fe,egrav*fe,      &
               ytot*fe,tnmass,dave(1),xzrtot(config%qx),pave(config%qx), &
               dave(config%qx),uave(config%qx),ugrtot(config%qx),        &
               dmdtout,dedtout,dydtout,phi*fe                       

      write(16,'(1x,7(1x,1pe11.4))') &
               tent*fe, tnut*fe, slnu*fe, snnu*fe, sumqen*fe,       &
               sumqye*fe,sumqmo*fe

!-to check conservation of energy and particle-number
      write(16,'(1x,7(1x,1pe15.8))') &
               sumioe*fe, sumdegr*fe, sumdenu*fe,  &
               sumion/pc_mb*fe,sumdynu/pc_mb*fe,   &
               sumioe_n*fe,sumion_n/pc_mb*fe
      close(16)

   end if ! myproc

      sumdegr = 0._rk  ! total potential energy sourceterm*dt 
                    !   summed over the timesteps with no output
      sumdenu = 0._rk  ! total neutrino  energy sourceterm*dt 
                    ! summed over the timesteps with no output
      sumdynu = 0._rk  ! total neutrino  number sourceterm*dt 
                    !  summed over the timesteps with no output
#ifdef FCNC_CALC
      sumdenu_n = 0._rk  ! total neutrino  energy sourceterm*dt 
                    ! summed over the timesteps with no output
      sumdynu_n = 0._rk  ! total neutrino  number sourceterm*dt 
                    !  summed over the timesteps with no output
#endif
      sumioe = 0._rk   ! total matter-energy loss via upper boundary
                    !   summed over the timesteps with no output
      sumion = 0._rk   ! total electron-number loss via upper boundary
                    ! summed over the timesteps with no output
      sumiom = 0._rk   ! total mass loss via upper boundary
                    !  summed over the timesteps with no output
      sumioe_n = 0._rk ! total neutrino-energy loss via upper boundary
                    !   summed over the timesteps with no output
      sumion_n = 0._rk ! total (e-)neutrino-number loss via upper boundary
                    ! summed over the timesteps with no output

!----------------------------------------------------------------------
!     write output file of the time evolution: neutrino quantities
!----------------------------------------------------------------------

      if (myproc.eq.0) then
      if (config%isma .gt. 1) then
         open(17,file = config%neutrino_file,position = filpos, form = 'formatted')

         if(nstep .eq. 1) then
            write(17,'(''#'','' step:  '',1x,'' time [s]     '',   &
     &                1x,''  E_en [foe]  '',1x,''  E_ean[foe]  '', &
     &                1x,''  E_mn [foe]  '',1x,''  E_man[foe]  '', &
     &                1x,''  E_tn [foe]  '',1x,''  E_tan[foe]  '', &
     &                1x,''  N_en [fo ]  '',1x,''  N_ean[fo ]  '', &
     &                1x,''  N_mn [fo ]  '',1x,''  N_man[fo ]  '', &
     &                1x,''  N_tn [fo ]  '',1x,''  N_tan[fo ]  '', &
     &                1x,''  L_en [foe/s]'',1x,''  L_ean[foe/s]'', &
     &                1x,''  L_mn [foe/s]'',1x,''  L_man[foe/s]'', &
     &                1x,''  L_tn [foe/s]'',1x,''  L_tan[foe/s]''  &
     &           )')
            write(17,'(''#'',195(''=''))')
         endif

         write(17,'(1x,i8,19(1x,1pe12.5))') nstep,time, &
              tens(1)*1.e-51_rk,tens(2)*1.e-51_rk,tens(3)*1.e-51_rk, &
              tens(4)*1.e-51_rk,tens(5)*1.e-51_rk,tens(6)*1.e-51_rk, &
              tnus(1)*1.e-51_rk,tnus(2)*1.e-51_rk,tnus(3)*1.e-51_rk, &
              tnus(4)*1.e-51_rk,tnus(5)*1.e-51_rk,tnus(6)*1.e-51_rk, &
              sln(1)*1.e-51_rk,sln(2)*1.e-51_rk,sln(3)*1.e-51_rk,    &
              sln(4)*1.e-51_rk,sln(5)*1.e-51_rk,sln(6)*1.e-51_rk
         close(17)
      endif

!----------------------------------------------------------------------
!     detailed informations about energies:
!----------------------------------------------------------------------

      open(18,file = config%energy_file,position = filpos,form = 'formatted')

      if(nstep .eq. 1) then
         write(18,'(''#'','' step:  '',1x,'' time [s]   '', &
                &  1x,''E_tot[^50 e]'',1x,''E_kin[^50 e]'', &
                &  1x,''E_int[^50 e]'',1x,''E_ene[^50 e]'', &
                &  1x,''E_nut[^50 e]'',1x,''E_gra[^50 e]'', &
                &  1x,''E_pos[^50 e]'',1x,''E_neg[^50 e]'', &
                &  1x,''E_kix[^50 e]'',1x,''E_kiy[^50 e]'', &
                &  1x,''E_kiz[^50 e]'',1x,''L_z    [^48]'', &
                &  1x,''T_zz   [^45]'')')
         write(18,'(''#'',183(''=''))')
!         write(18,'(''#'',i8,14(1x,1pe12.5))')                          &
!              nstep,time,etot,ekin,eint,tten,tent,egrav,ebindp,ebindm,  &
!              ekinx,ekiny,ekinz,amoz,moiz
      endif
      fe=1e-50_rk
!      else
         write(18,'('' '',i8,14(1x,1pe12.5))') &
              nstep,time,etot*fe,ekin*fe,eint*fe,tten*fe,tent*fe,      &
              egrav*fe,ebindp*fe,ebindm*fe,ekinx*fe,ekiny*fe,ekinz*fe, &
              amoz,moiz
!      endif

      close(18)

      end if !myproc.eq.0

#ifndef DEBUG_TIMINGS
      call second_v(tim1)
#endif
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#ifndef DEBUG_TIMINGS
      call second_v(tim2)

      timer%transp_comm = timer%transp_comm + (tim2-tim1)
#endif

end subroutine print_timestep_evolution
#endif /* NOTRA */


!> \verbatim
!> This subroutine writes information on the choice of the actuall
!> timestep
!>
!>  Author: A. Marek, MPA, March 2009
!> \endverbatim
!>
!> \param nstep     number of timestep
!> \param ti_cyc    actuall time
!> \param dt_min    hydro timestep
!> \param dt_nex    transport timestep
!> \param delyem_o  maximum change of Ye
!> \param iyem_o    zone of maximum change of Ye
!> \param deltem_o  maximum change of temperature
!> \param item_o    zone of maximum change of temperature
!> \param deleni_o  maximum change of energy density
!> \param ieni_o    zone of maximum change of energy density
!> \param delden_o  maximum change of density
!> \param iden_o    zone of maximum change of density
!> \param sigma_o   change of neutrino quantities (J,H..)
!> \param sigma2_o  change of neutrino quantities (J,H..)
!> \param nray      ray (in multi-d simulations) where changes occured
!> \param nhystp    number of hydro sub timesteps per transport timesteps
!> \param nresum    number of repitions of timestep
!> \param ndt    
!> \param itum
!> \param nenum
!> \param dtmaxx    size of maximum allowed timestep
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine write_timestep_information(nstep,ti_cyc,dt_min,dt_nex,  &
                                      delyem_o, iyem_o, deltem_o,  &
                                      item_o, deleni_o, ieni_o,    &
                                      delden_o, iden_o, sigma_o,   &
                                      sigma2,nray, nhystp,nresum,  &
                                      ndt,itnum,nenum,dtmaxx)    

  use precision
!  use abort

  use mo_mpi
  use charac
#ifndef DEBUG_TIMINGS
  use cputim
#endif
#ifdef EXTRACT_GW
  use gw3d
#endif /*EXTRACT_GW*/

  use configure
  implicit none
! local variables that are not in modules

  integer(kind=ik), intent(in):: nstep,nhystp,nresum,ndt,itnum,nenum
  integer(kind=ik), intent(in):: iyem_o,item_o,ieni_o,iden_o,nray
  real(kind=rk), intent(in):: ti_cyc,dt_min,dt_nex,delyem_o
  real(kind=rk), intent(in):: deltem_o, deleni_o,delden_o
  real(kind=rk), intent(in):: sigma_o,sigma2,dtmaxx

  real(kind=rk) :: tim1(2), tim2(2)
  integer(kind=ik) :: ierr

  if (myproc.eq.0) then
       open(19,file = config%timestep_file,position = filpos,form = 'formatted')             

     if(nstep .eq. 1) then
        write(19,'(''#'','' step:  '',1x,'' time [s]'',1x,'' dthy [s]'', &
                 & 1x,'' dtnu [s]'',1x,''      delyem  '',               &
                 & 1x,''     deltem  '',1x,''     deleni  '',            &
                 & 1x,''     delrho  '',1x,''     sigma  '',             &
                 & 1x,''ray   '',1x,'' nhy'',1x,''srep'',1x,''ndt'',   &
                 & 1x,''itn'',1x,''nen'',1x,''max'')')
        write(19,'(''#'',141(''=''))')
     end if

     write(19,'(i7,1pe11.4,2(1x,1pe9.2),4(1x,1pe9.2,1i4),1(1x,1pe9.2),   &
              & 1x,i4,i8,i3,1x,i4,2(1x,i2),1x,1pe9.2)') nstep,ti_cyc,    &
                 dt_min,dt_nex,delyem_o, iyem_o, deltem_o, item_o,       &
                 deleni_o, ieni_o,delden_o, iden_o, sigma_o/sigma2,nray, &
                 nhystp,nresum,ndt,itnum,nenum,dtmaxx
     
     close(19)

#ifdef EXTRACT_GW
      open(300,file = config%gw_file,position = filpos,form = 'formatted')
      if (any(qijq2(:,:).ne.0.0_rk)) then
         write(300,'(1x,i8,19(1x,1pe24.12))') &
              nstep,ti_cyc,qijq1,qijq2
      end if
      close(300)
#endif /*EXTRACT_GW*/

  end if
#ifndef DEBUG_TIMINGS
  call second_v(tim1)
#endif
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#ifndef DEBUG_TIMINGS
  call second_v(tim2)
  
  timer%hydro_comm =timer%hydro_comm + (tim2-tim1)
#endif

  return
end subroutine write_timestep_information

end module timeinfo

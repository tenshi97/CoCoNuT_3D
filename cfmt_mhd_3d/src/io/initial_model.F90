#undef DEBUG_OUTPUT

module initial_model

  private
  public :: initial_model_energy_summary

  contains


!>
!> \verbatim
!> This subroutine writes some energy information for an new model
!> to standard output
!> 
!>  Author: A. Marek MPA
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine initial_model_energy_summary
  use precision
  use phycon


  use totare_hy, only :                      xzntot, yzntot, zzntot, &
                                             xzltot, yzltot, zzltot, &
                                             xzrtot, yzrtot, zzrtot, &
#ifdef CFC_TRANSPORT2
                                             wltot, &
#endif
                                             dentot, dvxtot, dvytot, &
                                             dvztot, xnutot, enetot, &
                                             vextot, veytot, veztot, &
                                             gpotot, pretot, temtot, &
                                             stotot


  use totgrq_hy, only : tgmtot

  use intgrs_hy, only :  nstep
  use gfloat_hy, only : vlfrac, srfint
 
!  use intgrs_hy, only : nsdim
#ifndef NOTRA
  use neutrinotypes, only : siglep
#endif
  use nutrio_hy, only : dnutot, enutot, sumqen, sumqye, qentot, qyetot

#ifdef CFC_TRANSPORT
  use size_cfc, only: n_s,n_e,o_s,o_e
  use metric_cfc, only: sqrt_gamma
#ifdef NOTRA
  use grids, only: init_tot_arrays
#else
  use grids, only: init_tot_arrays, map_nutra2cfc
#endif
#endif /* CFC_TRANSPORT */


!  use lagradq_rt
!  use lagback_rt

  use mo_mpi
  use print_stdout_mod
  use configure

  implicit none

  integer(kind=ik) :: i,j,k,is

  real(kind=rk) :: tm(config%qx), dm(config%qx), tgm(config%qx)
  real(kind=rk) :: dsurf, dvol, dmas, rin, elto, etnu, dgrav, edis
  real(kind=rk) :: ebindp, ebindm, egrav, etot, dmint, gamfac, dvlme
  real(kind=rk) :: etn(config%isma), tten, ekinx, ekiny, ekinz, ytot, ekin,eint
  integer(kind=ik) :: ierr


  real(kind=rk) :: dm_rcv(config%qx), etn_rcv(config%isma)


  real(kind=rk) :: sumqen_rcv
  real(kind=rk) :: ekinx_rcv
  real(kind=rk) :: ekiny_rcv
  real(kind=rk) :: ekinz_rcv
  real(kind=rk) :: ebindp_rcv
  real(kind=rk) :: ebindm_rcv
  real(kind=rk) :: egrav_rcv
  real(kind=rk) :: elto_rcv
  real(kind=rk) :: rin_rcv
  real(kind=rk) :: tten_rcv
  real(kind=rk) :: etot_rcv
  real(kind=rk) :: sumqye_rcv


#ifdef DEBUG_OUTPUT

  real(kind=rk) ::    uave (config%qx), pave (config%qx), dave (config%qx), &
                       tave(config%qx), vave (config%qx), wave (config%qx), &
                      save (config%qx), yave(config%qx),  u2ave(config%qx), &
                      v2ave(config%qx), w2ave(config%qx)
#endif

  dm(:)  = 0.0_rk
  dmas   = 0.0_rk
  rin    = 0.0_rk
  elto   = 0.0_rk
  etn(:) = 0.0_rk
  tten   = 0.0_rk


!      pmass  = tmatot(0)

!      emec2  = 511.e3_rk * 1.6021e-12_rk / 1.66e-24_rk   !  m_e * c^2 / m_u

  etot   = 0.0_rk
  ekinx  = 0.0_rk
  ekiny  = 0.0_rk
  ekinz  = 0.0_rk
  sumqen = 0.0_rk
  sumqye = 0.0_rk
!
  egrav  = 0.0_rk
  ebindp = 0.0_rk
  ebindm = 0.0_rk
  
  dvxtot(1:config%qx) = (xzrtot(1:config%qx)**3 - xzltot(1:config%qx)**3)  &
                     / 3._rk

#ifdef CFC_TRANSPORT
  call init_tot_arrays(n_s,n_e,o_s,o_e,-1)
#ifndef NOTRA
  if (config%p_ntr .ne. 0) call map_nutra2cfc
#endif
#endif /* CFC_TRANSPORT */

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

           tten  = tten   +         dmas * enetot(i,j,k)
           ekinx = ekinx + 0.5_rk * dmas * vextot(i,j,k)**2
           ekiny = ekiny + 0.5_rk * dmas * veytot(i,j,k)**2
           ekinz = ekinz + 0.5_rk * dmas * veztot(i,j,k)**2
           etn(:)= etn(:)+       dvol * enutot(i,j,k,:)

#ifndef NOTRA
           if (config%p_ntr .ne. 0) then
             do is=1,config%isma         
                rin = rin + dvol * siglep(is)*dnutot(i,j,k,is)
             enddo
          endif
#endif /* NOTRA */

           sumqen = sumqen + qentot(i,j,k) * dvol 
           !            sumqmo = sumqmo + qmotot(i,j,k)* dvol
           sumqye = sumqye + qyetot(i,j,k,1) * dvol

#ifdef CFC_TRANSPORT2
           dgrav = 0.0_rk
#else
           dgrav = 0.125_rk * dmas * (         &
                       + gpotot(i  ,j  ,k)  &
                       + gpotot(i  ,j-1,k)  &
                       + gpotot(i-1,j  ,k)  &
                       + gpotot(i-1,j-1,k))
#endif /* CFC_TRANSPORT2 */

            egrav = egrav + dgrav
            edis  = enetot(i,j,k) * dmas  +  dgrav

!-PNS        edis  = (enetot(i,j,k) - emec2*xnutot(i,j,k,config%qn)) * dmas 
!-PNS     &              + dgrav
!#else
!-PNS        edis  = (enetotg(i,j,k) - emec2*xnutotg(i,j,k,config%qn)) * dmas


            elto  = elto + xnutot(i,j,k,config%qn)*dmas

            if (edis .gt. 0._rk)  then
               ebindp = ebindp + edis
            else
               ebindm = ebindm + edis
            end if
            
         enddo ! i-loop
      enddo !j - loop
   enddo ! k-loop


   if (use_mpi) then
      dm_rcv = 0._rk
      call MPI_AllReduce(dm, dm_rcv, config%qx, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      dm = dm_rcv

      rin_rcv = 0._rk
      call MPI_AllReduce(rin, rin_rcv, 1, MPI_DOUBLE_PRECISION,&
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      rin = rin_rcv
      
      elto_rcv = 0._rk
      call MPI_AllReduce(elto, elto_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      elto = elto_rcv

      etn_rcv = 0._rk
      call MPI_AllReduce(etn, etn_rcv, config%isma, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      etn = etn_rcv
      
      tten_rcv = 0._rk
      call MPI_AllReduce(tten, tten_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      tten = tten_rcv

      etot_rcv = 0._rk
      call MPI_AllReduce(etot, etot_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      etot = etot_rcv

      ekinx_rcv = 0._rk
      call MPI_AllReduce(ekinx, ekinx_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      ekinx = ekinx_rcv

      ekiny_rcv = 0._rk
      call MPI_AllReduce(ekiny, ekiny_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      ekiny = ekiny_rcv

      ekinz_rcv = 0._rk
      call MPI_AllReduce(ekinz, ekinz_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      ekinz = ekinz_rcv

      sumqen_rcv = 0._rk
      call MPI_AllReduce(sumqen, sumqen_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      sumqen = sumqen_rcv

      sumqye_rcv = 0._rk
      call MPI_AllReduce(sumqye, sumqye_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      sumqye = sumqye_rcv


      egrav_rcv = 0._rk
      call MPI_AllReduce(egrav, egrav_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      egrav = egrav_rcv

      ebindp_rcv = 0._rk
      call MPI_AllReduce(ebindp, ebindp_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      ebindp = ebindp_rcv

      ebindm_rcv = 0._rk
      call MPI_AllReduce(ebindm, ebindm_rcv, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      ebindm = ebindm_rcv
   endif ! use_mpi


   etnu = 0.0_rk
 
   do is=1,config%isma
      etnu=etnu+etn(is)
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
         

   elto = elto/pc_mb
   if (config%p_nbk .gt. 0) then
      etot = tten + egrav + etnu
      ytot = rin + elto
   else
      etot = tten + egrav
      ytot = elto
   endif

   ekin = ekinx + ekiny + ekinz
   eint = tten - ekin
         

             

#ifdef DEBUG_OUTPUT
   if (use_mpi) then
       raise_abort("initial_model(): debug output not supported in mpi mode")
    endif

       do 100 i = 1, config%qx
!
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
      !         gamfac   = 2.0/(gamtot(i-1) + gamtot(i)) ! rest mass
          gamfac   = 1._rk
         

!----------------------------
          if (config%nsdim .eq. 3)  then
!----------------------------
!
             do k = 1, config%qz
                do j = 1, config%qy
                   dvlme    = dvztot(k) * dvytot(j)
               ! vlfrac*dvxtot(i) cancels out!


#ifndef CFC_TRANSPORT2
                   dmas     = dentot(i,j,k) * dvlme
#else /*CFC_TRANSPORT2 */
                   dmas     = dentot(i,j,k) * dvlme &
                          * sqrt_gamma(i,j,k) * wltot(i,j,k)
#endif /*CFC_TRANSPORT2 */

                   dmint    = dmint + dmas * gamfac
                   uave (i) = uave (i)  +  vextot(i,j,k)    * dvlme
                   vave (i) = vave (i)  +  veytot(i,j,k)    * dvlme
                   wave (i) = wave (i)  +  veztot(i,j,k)    * dvlme
                   u2ave(i) = u2ave(i)  +  vextot(i,j,k)**2 * dmas
                   v2ave(i) = v2ave(i)  +  veytot(i,j,k)**2 * dmas
                   w2ave(i) = w2ave(i)  +  veztot(i,j,k)**2 * dmas
                   dave (i) = dave (i)  +  dentot(i,j,k)    * dvlme
                   pave (i) = pave (i)  +  pretot(i,j,k)    * dvlme
                   tave (i) = tave (i)  +  temtot(i,j,k)    * dvlme
                   save (i) = save (i)  +  stotot(i,j,k)    * dmas
                   yave (i) = yave (i)  +  xnutot(i,j,k,config%qn)* dmas
               enddo ! j-loop
           enddo  ! k-loop

!----------------------------
         end if !config%nsdim eq 3


         if (config%nsdim .eq. 2)  then
!----------------------------
!
            k = 1
            !
            do j = 1, config%qy
               dvlme    = dvytot(j)
                          ! vlfrac*dvxtot(i)*dvztot(k) cancels out!
               dmas     = dentot(i,j,k) * dvlme
               dmint    = dmint + dmas * gamfac
               uave (i) = uave (i)  +  vextot(i,j,k)    * dvlme
               vave (i) = vave (i)  +  veytot(i,j,k)    * dvlme
               wave (i) = wave (i)  +  veztot(i,j,k)    * dvlme
!-old               uave (i) = uave (i)  +  vextot(i,j,k)
!-old               vave (i) = vave (i)  +  veytot(i,j,k)

               u2ave(i) = u2ave(i)  +  vextot(i,j,k)**2 * dmas
               v2ave(i) = v2ave(i)  +  veytot(i,j,k)**2 * dmas
               w2ave(i) = w2ave(i)  +  veztot(i,j,k)**2 * dmas

               dave (i) = dave (i)  +  dentot(i,j,k)    * dvlme
               pave (i) = pave (i)  +  pretot(i,j,k)    * dvlme
               tave (i) = tave (i)  +  temtot(i,j,k)    * dvlme

!-old               dave (i) = dave (i)  +  dentot(i,j,k)
!-old               pave (i) = pave (i)  +  pretot(i,j,k)
!-old               tave (i) = tave (i)  +  temtot(i,j,k)

               save (i) = save (i)  +  stotot(i,j,k)    * dmas
!               yave (i) = yave (i)  +  (xnutot(i,j,k, 1) *  1.0 +
!     &                                  xnutot(i,j,k, 2) *  4.0 +
!     &                                  xnutot(i,j,k, 3) * 12.0 +
!     &                                  xnutot(i,j,k, 4) * 16.0 +
!     &                                  xnutot(i,j,k, 5) * 20.0 +
!     &                                  xnutot(i,j,k, 6) * 24.0 +
!     &                                  xnutot(i,j,k, 7) * 28.0 +
!     &                                  xnutot(i,j,k, 8) * 56.0 +
!     &                                  xnutot(i,j,k, 9) * 56.0 +
!     &                                  xnutot(i,j,k,10) * 56.0 )
!     &                                  * dmas
               yave (i) = yave (i)  +  xnutot(i,j,k,config%qn)* dmas ! -> Y_e
            enddo ! j-loop

!----------------------------
      end if ! config%nsdim eq 2

      if (config%nsdim .eq. 1)  then
!----------------------------
!
         uave (i) = vextot(i,1,1)
         u2ave(i) = vextot(i,1,1)**2
         dave (i) = dentot(i,1,1)
         pave (i) = pretot(i,1,1)
         tave (i) = temtot(i,1,1)
         save (i) = stotot(i,1,1)
         yave (i) = xnutot(i,1,1,config%qn)
!
         go to 100

!----------------------------
      end if ! config%nsdim eq 1
!----------------------------

      uave (i) = uave (i) / srfint
      vave (i) = vave (i) / srfint
      wave (i) = wave (i) / srfint
!-old         uave (i) = uave (i) * vnorm
!-old         vave (i) = vave (i) * vnorm
!-old         wave (i) = wave (i) * vnorm

      u2ave(i) = u2ave(i) / dmint
      v2ave(i) = v2ave(i) / dmint
      w2ave(i) = w2ave(i) / dmint
      dave (i) = dave (i) / srfint
      pave (i) = pave (i) / srfint
      tave (i) = tave (i) / srfint
!-old         dave (i) = dave (i) * vnorm
!-old         pave (i) = pave (i) * vnorm
!-old         tave (i) = tave (i) * vnorm

      save (i) = save (i) / dmint
      yave (i) = yave (i) / dmint

100   continue

#endif /* DEBUG_OUTPUT */

      if((config%itstp .ne. 0) .or. (nstep .eq. 0) ) then

         call printit_taskX(0,"Punktmasse ............ [Msol] ",config%pmass/pc_ms)
         call printit_taskX(0,"Gesamtmasse ........... [Msol] ",tm(config%qx))
         call printit_taskX(0,"Gesamtmasse (grav.).... [Msol] ",tgm(config%qx))
         call printit_taskX(0," ")
         call printit_taskX(0,"Neutrino Number  ......        ",rin)
         call printit_taskX(0,"Elec-Posit Number .....        ",elto)
         call printit_taskX(0,"Neutrino Energie ......  [erg] ",etnu)
         call printit_taskX(0,"thermische Energie ....  [erg] ",eint)
         call printit_taskX(0,"Gravitationsenergie ...  [erg] ",egrav) 
         call printit_taskX(0,"kinetische Energie ....  [erg] ",ekin)
         call printit_taskX(0,"E_kin_x ...............  [erg] ",ekinx)
         call printit_taskX(0,"E_kin_y ...............  [erg] ",ekiny)
         call printit_taskX(0,"E_kin_z ...............  [erg] ",ekinz)
         call printit_taskX(0," ")
         call printit_taskX(0,"Gesamtenergie .........  [erg] ",etot)
         call printit_taskX(0," ")
         call printit_taskX(0,"Gesamtleptonenzahl ....  [#]   ",ytot)
         call printit_taskX(0," ")
         call printit_taskX(0,"Sourcet Neutrino Energy [erg/s]",sumqen)
         call printit_taskX(0," ")
         call printit_taskX(0,"E_bin_positiv .........  [erg] ",ebindp)
         call printit_taskX(0,"E_bin_negativ .........  [erg] ",ebindm)
         call printit_taskX(0," ")
         call printit_taskX(0," ")
         
      endif ! config%itstp .ne. 0 

end subroutine initial_model_energy_summary

end module initial_model

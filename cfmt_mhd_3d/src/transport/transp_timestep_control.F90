module transp_timestep

  private

  public :: transport_timestep_control

  contains

    subroutine transport_timestep_control(yehtot, dehtot, tehtot, enhtot,    &
                                          dt_hyd, nhyrep, dt_nex,dt_cyc,     &
                                          ti_cyc, itnum, nenum, sigma, nray, &
                                          dt_min, ndt, selftime, childrentime)

      use precision
      use configure
      use totare_hy
      use multigrid_rt
      use mo_mpi
      use cputim
      use state
#ifndef PROGRAM_remap
      use timeinfo, only : write_timestep_information
#endif
      use grid_mod
      use grids
      use print_stdout_mod

      use cputim
      use configure
      use state
      use hydro_areas_mod
      use state

      use neutrino_mod

#ifdef CHECK_MEMORY
      use meminfo
#endif
#ifdef CFC_TRANSPORT
      use size_cfc
      use hydro_primitives_cfc
      use phycon, only: pc_geog
#endif
      use ioflx
      use specfun
      use totare_hy

      use sum_fluxes
      use savare_overload
      use hydro_interface

      use abort
      use error

      !  use arecon_hy
      use massio_hy
      use intgrs_hy
      use nutrio_hy

      use param_rt
#ifndef NOTRA
      use multigrid_rt
      use neutrinotypes
#endif

      implicit none

      real(kind=rk), &
      intent(in), dimension(config%qx,qy_s:qy_e,qz_s:qz_e) :: &
                                         yehtot, dehtot, tehtot, enhtot

      real(kind=rk), intent(inout)    :: dt_hyd, ti_cyc, dt_cyc
      real(kind=rk), intent(in)       :: dt_nex, sigma, dt_min
      integer(kind=ik), intent(in)    :: nhyrep, ndt
      integer(kind=ik), intent(inout) :: itnum, nenum, nray(:)
      real(kind=rk), intent(out)      :: selftime(2), childrentime(2)


      integer(kind=ik)                :: imrad, ii, iyem(3), iden(3), &
                                         item(3), ieni(3)

      real(kind=rk)                   :: delyem, delden_drhmx, deltem, deleni
      real(kind=rk), parameter        :: dyemx=1.0e-2_rk, dtemx=4.0e-1_rk
      real(kind=rk), parameter        :: drhmx=0.0e-2_rk, deimx=3.0e-2_rk
      real(kind=rk), parameter        :: drhmx_l=3.0e-2_rk, drhmx_h=1.0e-2_rk
      real(kind=rk), parameter        :: rh_cr=3.0e10_rk, drh_stf=3._rk
      real(kind=rk), parameter        :: deimx_h=3.0e-2_rk, deimx_l=1.0e-2_rk

      real(kind=rk), dimension(2)     :: time_communication_start,           &
                                         time_communication_end
      real(kind=rk), dimension(2)     :: selftime_start(2)
      real(kind=rk)                   :: delyem_o=0._rk, delden_o=0._rk,   &
                                         deltem_o=0._rk, deleni_o=0._rk,   &
                                         sigma_o=0._rk
      integer(kind=ik)                :: iyem_o=0, iden_o=0, item_o=0, ieni_o=0


      real(kind=rk)                   :: dd
      real(kind=rk)                   :: rbuf(1,2), rbuf_buf(1,2)
      integer(kind=ik)                :: ndt_rcv
      real(kind=rk)                   :: dt_min_rcv

      integer(kind=ik)                :: ierr

      integer(kind=ik)                :: is, i

      selftime          = 0._rk
      childrentime      = 0._rk

#ifndef DEBUG_TIMINGS
      call second_v(selftime_start)
#endif

      imrad=config%qx

#ifndef NOTRA
      if (config%p_ntr .ne. 0) then
         if (config%use_multid_collapse) then
            do ii=config%qx-1,1,-1
               if (xzrtot(ii) .gt. ralag_3d(config%imaxp+1,config%nystrt,config%nztra)) imrad=ii
            enddo
         else
            do ii=config%qx-1,1,-1
               if (xzrtot(ii) .gt. ralag_1d(config%imaxp+1)) imrad=ii
            enddo
         endif ! multid_collapse
      endif
#endif /* NOTRA */

      imrad = min(imrad+3,config%qx)
      delyem = 0._rk
      delden_drhmx = 0._rk
      deltem = 0._rk
      deleni = 0._rk

      if (config%p_nbk .gt. 0) then ! do not compute if no neutrino changes allowed
         delyem=maxval( abs(yehtot(1:imrad,:,:) - xnutot(1:imrad,:,:,config%qn)) &
                       /yehtot(1:imrad,:,:)/dyemx )

         iyem  =MAXLOC( abs(yehtot(1:imrad,:,:) - xnutot(1:imrad,:,:,config%qn)) &
              /yehtot(1:imrad,:,:)/dyemx )


         !  if (p_nbk .gt. 0) then ! do not compute if no neutrino changes allowed
         !     iyem  ( abs(yehtot(1:imrad,:,:) - xnutot(1:imrad,:,:,config%qn)) &
!                       /yehtot(1:imrad,:,:)/dyemx )
         !  endif

         delden_drhmx=maxval( abs(dehtot(1:imrad,:,:) - dentot(1:imrad,:,:)) &
                          /dentot(1:imrad,:,:)/ &
                         ((drhmx_h*dentot(1:imrad,:,:)**drh_stf  &
                         + drhmx_l*rh_cr**drh_stf)/ &
                          (dentot(1:imrad,:,:)**drh_stf + rh_cr**drh_stf)))

         iden  =maxloc( abs(dehtot(1:imrad,:,:) - dentot(1:imrad,:,:)) &
                    /dentot(1:imrad,:,:)/ &
                    ((drhmx_h*dentot(1:imrad,:,:)**drh_stf  &
                    + drhmx_l*rh_cr**drh_stf)/ &
                     (dentot(1:imrad,:,:)**drh_stf + rh_cr**drh_stf)))

         deltem=maxval( abs(tehtot(1:imrad,:,:) - temtot(1:imrad,:,:)) &
                    /tehtot(1:imrad,:,:)/dtemx )

         item  =maxloc( abs(tehtot(1:imrad,:,:) - temtot(1:imrad,:,:)) &
                    /tehtot(1:imrad,:,:)/dtemx )

         deleni=maxval( abs(enhtot(1:imrad,:,:) - (enetot(1:imrad,:,:) &
                              -0.5_rk*(vextot(1:imrad,:,:)**2   &
                                      +veytot(1:imrad,:,:)**2   &
                                      +veztot(1:imrad,:,:)**2)  &
                                                             )) &
             /enhtot(1:imrad,:,:)/((deimx_h*dentot(1:imrad,:,:)**drh_stf  &
             + deimx_l*1.0e+9_rk**drh_stf)/ &
              (dentot(1:imrad,:,:)**drh_stf + 1.0e+9_rk**drh_stf)))

         ieni= maxloc( abs(enhtot(1:imrad,:,:) -  (enetot(1:imrad,:,:) &

                -0.5_rk*(vextot(1:imrad,:,:)**2 +veytot(1:imrad,:,:)**2 &
                         +veztot(1:imrad,:,:)**2) &
                                                             )) &
             /enhtot(1:imrad,:,:)/((deimx_h*dentot(1:imrad,:,:)**drh_stf  &
               + deimx_l*1.0e+9_rk**drh_stf)/ &
              (dentot(1:imrad,:,:)**drh_stf + 1.0e+9_rk**drh_stf)))


         if (use_mpi) then
            call second_v(time_communication_start)

            rbuf(1,1) = delyem
            rbuf(1,2) = real(iyem(1),kind=rk) ! Attention we only exchange the radial position

            call MPI_AllReduce(rbuf, rbuf_buf, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

            delyem = rbuf_buf(1,1)
            iyem(1)   = int(rbuf_buf(1,2),kind=ik)


            rbuf(1,1) = delden_drhmx
            rbuf(1,2) = real(iden(1),kind=rk) ! Attention we only exchange the radial position

            call MPI_AllReduce(rbuf, rbuf_buf, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

            delden_drhmx = rbuf_buf(1,1)
            iden(1)   = int(rbuf_buf(1,2),kind=ik)


            rbuf(1,1) = deltem
            rbuf(1,2) = real(item(1),kind=rk)! Attention we only exchange the radial position

            call MPI_AllReduce(rbuf, rbuf_buf, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

            deltem = rbuf_buf(1,1)
            item(1)   = int(rbuf_buf(1,2),kind=ik)

            rbuf(1,1) = deleni
            rbuf(1,2) = real(ieni(1),kind=rk)! Attention we only exchange the radial position

            call MPI_AllReduce(rbuf, rbuf_buf, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

            deleni = rbuf_buf(1,1)
            ieni(1)   = int(rbuf_buf(1,2),kind=ik)

            call second_v(time_communication_end)

            timer%transp_comm = timer%transp_comm + (time_communication_end-time_communication_start)

         endif ! use_mpi

         !#endif

      endif ! p_nbk .ne. 0

!      dd=max(delyem,deltem,deleni,delden_drhmx)
      dd=max(delyem,delden_drhmx)
      if (dd .gt. 1.0_rk) then
         dt_hyd=dt_hyd*0.75_rk
      elseif (dd.lt.0.7_rk) then
         dt_hyd=dt_hyd*min(1.2_rk,1.0_rk/max(dd,1.0e-4_rk))
         dt_hyd=min(dt_hyd,config%dtmaxx)
      endif

      !   do not use a larger timestep if transp was repeated recently
      !    (dt_nex still contains the old timestep)
      if (nhyrep .gt. 0) transport%dt=min(transport%dt,dt_nex)


!------------------
#ifndef NOTRA
      if (config%p_ntr .ne. 0) then
         !   Neutrino-energy Luminosity
         do is=1,config%isma
            sumioe_n=sumioe_n + neu_qw(is) * sumcqflux( &

                 fnutot(config%qx,qy_s:qy_e,qz_s:qz_e,is)*xzrtot(config%qx)**2,dt_nex)

         enddo
         !   Net lepton Luminosity
         do is=1,config%isma
            sumion_n=sumion_n+siglep(is) * sumcqflux( &

                 gnutot(config%qx,qy_s:qy_e,qz_s:qz_e,is)*xzrtot(config%qx)**2,dt_nex)

         enddo

      endif
#endif

      ! SUMMARY OUTPUT
      nresum = nresum + nhyrep

      if (delyem.gt.delyem_o) then
         delyem_o=delyem
         iyem_o=iyem(1)
      endif
      if (delden_drhmx.gt.delden_o) then
         delden_o=delden_drhmx
         iden_o=iden(1)
      endif
      if (deltem.gt.deltem_o) then
         deltem_o=deltem
         item_o=item(1)
      endif
      if (deleni.gt.deleni_o) then
         deleni_o=deleni
         ieni_o=ieni(1)
      endif
      !      delyem_o=max(delyem_o,delyem)
      !      delden_o=max(delden_o,delden_drhmx)
      !      deltem_o=max(deltem_o,deltem)
      !      deleni_o=max(deleni_o,deleni)
      sigma_o=max(sigma_o,sigma)


      if (config%p_ntr .eq. 0) then
         itnum   = 0
         nenum   = 0
         nray(:) = 0
      endif


      if(nstep .lt. 10 .or. mod(nstep,config%intout) .eq. 0) then
         if (.not.(config%use_deactivate_controlFiles)) then
            call write_timestep_information(nstep,ti_cyc,dt_min,dt_nex, &
                 delyem_o, iyem_o, deltem_o, item_o, deleni_o, ieni_o, &
                 delden_o, iden_o, sigma_o,config%sigma2,nray(1), &
                 areas%nhystp,nresum,ndt,itnum,nenum,config%dtmaxx)
         endif

         delden_o=0._rk
         delyem_o=0._rk
         deltem_o=0._rk
         deleni_o=0._rk
         iden_o=0
         iyem_o=0
         item_o=0
         ieni_o=0
         sigma_o=0._rk


         nresum = 0

      endif


      dt_cyc = dt_nex           ! return dt of this cycle
      ti_cyc = ti_cyc + dt_cyc  ! return time at the end of the cycle

      if (dt_nex .lt. config%dtmin_rt) then
         call printit_taskx(0," W A R N I N G :   dt = ",dt_nex)
         call printit_taskx(0,"             dtmin_rt = ",config%dtmin_rt)
         call raise_error_x("src/transport/rady.F90:664", "rady(): => arecon")
         !         call stopit('=> arecon',0)
      endif

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
    end subroutine transport_timestep_control


end module transp_timestep

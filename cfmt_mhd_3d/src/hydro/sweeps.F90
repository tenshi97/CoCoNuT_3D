#undef KEIL_EPHI
#undef CONTROL_OUTPUT
#undef CONTROL_OUTPUT_2

#ifdef IBM_COMPILER
! this is necessary on the POWER6 since the xlf90 version 10.xx
! causes problems
@PROCESS HOT
#endif

module sweeps_mod

implicit none

contains

!=======================================================================

!>
!> \verbatim
!> Driver for hydro sweeps
!> \endverbatim
!>
!> \author W. Keil
!> \param isweep sweep (x,y,z) discriminator
!> \param hydro 
!> \param accel
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine sweep( isweep ,hydro_calculation, accel, caller, selftime, childrentime)

  use precision
  use abort

  use intgrs_hy

  use vnew_hy
  use vold_hy
  use mesh_hy
  use gfloat_hy
  use grd_hy

  use cputim

  use mo_mpi
  use hlpare_hy

  use mpi_comm_routines
  use eos3d_routine, only : eos3d
  use burn_mod

#ifdef HTCL
  use htcl_mod
#endif

  use hydro_areas_mod
  use configure
  use state
  implicit none
! LOCAL varibales that are not in modules

  integer(kind=ik), intent(in) :: caller
  integer(kind=ik)             :: sweep_mode
  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk)                :: selftime_start(2), hydrow_self(2),  &
                                  hydrow_children(2), accelx_self(2), &
                                  accelx_children(2), accely_self(2), &
                                  accely_children(2), eos3d_self(2),  &
                                  eos3d_children(2), burn_self(2),    &
                                  burn_children(2)
  real(kind=rk)                :: paralleltime_start(2), paralleltime_end(2)
  real(kind=rk)                :: hydrow_hlp1(2,2,3,2), hydrow_hlp2(2,2,3,2)

  integer(kind=ik)             :: n1, n2, ij, j, i, i4, j4, k, k4, ii, jj

  integer(kind=ik), intent(in) :: isweep 
  logical, intent(in)          :: hydro_calculation, accel
  logical                      :: ler
  real(kind=rk), dimension(2)  :: time_communication_start,    &
                                  time_communication_end, tim1,tim2

  real(kind=rk)                ::  sbuf2(config%qx,qz_s:qz_e), &
                                   rbuf2(config%qx,qz_s:qz_e), &
                                   sbuf3(config%qx,qy_s:qy_e), &
                                   rbuf3(config%qx,qy_s:qy_e)
  integer(kind=ik)             :: dest, src, req1, ierr,       &
                                  mpistat(MPI_STATUS_SIZE), sendcount

  integer(kind=ik)             :: j_end,k_end

  hydrow_hlp1  = 0._rk
  hydrow_hlp2  = 0._rk
  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  xyzswp = isweep

  if (use_mpi) then
     if (ioy.eq.1 .and. ioz.eq.1) then
        j_end=qy_e
        k_end=qz_e
     else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
        j_end=qy_s
        k_end=qz_s
     end if
  else
     j_end = areas%ny
     k_end = areas%nz
  endif ! use _mpi

  if     ( xyzswp .eq. 1 ) then
     if (use_mpi) then
        n1=k_end-qz_s+1
        n2=j_end-qy_s+1
     else
        n1 = areas%nz
        n2 = areas%ny
     endif

     bndmin = config%bndmnx
     bndmax = config%bndmxx
     bndbot = config%bndmny
     bndtop = config%bndmxy
     bndlft = config%bndmnz
     bndrgt = config%bndmxz

     if (use_mpi) then
        if (areas%ny .gt. 1) then
           if (qy_s .ne.  1 .or. config%bndmny .eq. 4) bndbot = 6
           if (qy_e .ne. areas%ny .or. config%bndmxy .eq. 4) bndtop = 6
        endif
        if (areas%nz .gt. 1) then
           if (qz_s .ne.  1 .or. config%bndmnz .eq. 4) bndlft = 7
           if (qz_e .ne. areas%nz .or. config%bndmxz .eq. 4) bndrgt = 7
        endif
     endif

     igeom = config%igeomx

     nzn  = areas%nx
     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8

  else if( xyzswp .eq. 2 ) then

     if (use_mpi) then
        n1=k_end-qz_s+1
     else
        n1 = areas%nz
     endif

     n2 = areas%nx

     bndmin = config%bndmny
     bndmax = config%bndmxy
     bndbot = config%bndmnx
     bndtop = config%bndmxx
     bndlft = config%bndmnz
     bndrgt = config%bndmxz

     if (use_mpi) then
        if (areas%ny .gt. 1) then
           if (qy_s .ne.  1       .or. config%bndmny .eq. 4) bndmin = 6
           if (qy_e .ne. areas%ny .or. config%bndmxy .eq. 4) bndmax = 6
           if (use_2neighbour_comm) then
              ! two zones reflective, two zones communication
              if (qy_s .eq. 3          .and. config%bndmny .eq. 1) bndmin = 9
              if (qy_e .eq. areas%ny-2 .and. config%bndmny .eq. 1) bndmax = 9
           endif
           if (use_4neighbour_comm) then
              ! one zone reflective, three communciation
              if (qy_s .eq. 4           .and. config%bndmny .eq. 1) bndmin = 8
              if (qy_e .eq. areas%ny-4  .and. config%bndmny .eq. 1) bndmax = 8

              ! two zones reflective, two communciation
              if (qy_s .eq. 3           .and. config%bndmny .eq. 1) bndmin = 9
              if (qy_e .eq. areas%ny-3  .and. config%bndmny .eq. 1) bndmax = 9

              ! three zones reflective, one communciation
              if (qy_s .eq. 2           .and. config%bndmny .eq. 1) bndmin = 10
              if (qy_e .eq. areas%ny-1  .and. config%bndmny .eq. 1) bndmax = 10

           endif

        end if
        if (areas%nz .gt. 1) then   
           if (qz_s .ne.  1 .or. config%bndmnz .eq. 4) bndlft = 7
           if (qz_e .ne. areas%nz .or. config%bndmxz .eq. 4) bndrgt = 7
        end if
     endif

     igeom = config%igeomy

     if (use_mpi) then
        nzn = qy_proc
     else
        nzn  = areas%ny 
     endif


     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8

  else if( xyzswp .eq. 3 ) then

     if (use_mpi) then
        n1 = j_end-qy_s+1
     else
        n1 = areas%ny
     endif

     n2 = areas%nx

     bndmin = config%bndmnz
     bndmax = config%bndmxz
     bndbot = config%bndmnx
     bndtop = config%bndmxx
     bndlft = config%bndmny
     bndrgt = config%bndmxy

     if (use_mpi) then
        if (areas%ny .gt. 1) then
           if (qy_s .ne.  1 .or. config%bndmny .eq. 4) bndlft = 6
           if (qy_e .ne. areas%ny .or. config%bndmxy .eq. 4) bndrgt = 6
        endif
        if (areas%nz .gt. 1) then
           if (qz_s .ne.  1 .or. config%bndmnz .eq. 4) bndmin = 7
           if (qz_e .ne. areas%nz .or. config%bndmxz .eq. 4) bndmax = 7
        endif
     endif

     igeom = config%igeomz

     if (use_mpi) then
        nzn  = qz_proc
     else
        nzn  = areas%nz
     endif

     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8

  else
     write(*,*) "Task ",myproc
     write(*,*) 'SWEEP: direction ',xyzswp,' is not supported'
     raise_abort("sweep():  ")

  end if

#ifdef CONTROL_OUTPUT
!-----------------------------------------------------------------------
!     control output:
!-----------------------------------------------------------------------
!
  write(*,*) "Task ",myproc
  write(*,*) 'bndmin/ bndmax/ bndbot/ bndtop/ bndlft/ bndrgt'
  write(*,'(6(1x,1pe12.5))') &
        bndmin, bndmax, bndbot, bndtop, bndlft, bndrgt

  k = 1
  write(*,*) "Task ",myproc
  write(*,*) 'ABC: r/densty/ vx / vy / T / Ye /gpot'

  do j = qy_s, j_end

     write(*,'('' j: '',i3,'' Y = '',1pe12.5)') j,yzn(j)
     do i = 1, areas%nx
              write(*,'(1x,i3,7(1x,1pe12.5))') i,xzn(i),   &
                    densty(i,j,k),velx(i,j,k),vely(i,j,k), &
                    temp(i,j,k),xnuc(i,j,k,config%qn),gpot(i,j,k)
      enddo
   enddo
   call printit_taskX(0," ")
!
!-----------------------------------------------------------------------
#endif /* CONTROL_OUTPUT */


   vxold(1:areas%nx,qy_s:j_end,qz_s:k_end) = velx(1:areas%nx,qy_s:j_end,qz_s:k_end)
   vyold(1:areas%nx,qy_s:j_end,qz_s:k_end) = vely(1:areas%nx,qy_s:j_end,qz_s:k_end)
   vzold(1:areas%nx,qy_s:j_end,qz_s:k_end) = velz(1:areas%nx,qy_s:j_end,qz_s:k_end)

   if (hydro_calculation) then


      if (use_mpi) then

         if (xyzswp .eq. 1 .and. nprocs .gt. 1) then
            ! do MPI next neigbor communication because of access in getrwx
            ! boundary conditions for config%qy=1, n are set in bndry()

            ! exchange vyold with neighbors
            ! send left in theta direction
#ifndef DEBUG_TIMINGS
            call second_v(time_communication_start)
#endif
            sbuf2(1:config%qx,qz_s:qz_e) = vyold(1:config%qx,qy_s,qz_s:qz_e)
            call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
            sendcount = config%qx*qz_proc
            call MPI_SendRecv(sbuf2, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep1, rbuf2, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vyold(1:config%qx,qy_e+1,qz_s:qz_e) = rbuf2(1:config%qx,qz_s:qz_e)

            ! send right in theta direction
            sbuf2(1:config%qx,qz_s:qz_e) = vyold(1:config%qx,qy_e,qz_s:qz_e)
            call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
            sendcount = config%qx*qz_proc
            call MPI_SendRecv(sbuf2, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep2, rbuf2, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep2,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vyold(1:config%qx,qy_s-1,qz_s:qz_e) = rbuf2(1:config%qx,qz_s:qz_e)

            ! exchange vzold with neighbors
            ! send down in phi direction
            sbuf3(1:config%qx,qy_s:qy_e) = vzold(1:config%qx,qy_s:qy_e,qz_s)
            call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
            sendcount = config%qx*qy_proc
            call MPI_SendRecv(sbuf3, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep1, rbuf3, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vzold(1:config%qx,qy_s:qy_e,qz_e+1) = rbuf3(1:config%qx,qy_s:qy_e)

            ! send up in phi direction
            sbuf3(1:config%qx,qy_s:qy_e) = vzold(1:config%qx,qy_s:qy_e,qz_e)
            call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
            sendcount = config%qx*qy_proc
            call MPI_SendRecv(sbuf3, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep2, rbuf3, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep2,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vzold(1:config%qx,qy_s:qy_e,qz_s-1) = rbuf3(1:config%qx,qy_s:qy_e)

#ifndef DEBUG_TIMINGS
            call second_v(time_communication_end)
            childrentime = childrentime + time_communication_end - time_communication_start
            timer%hydro_comm = timer%hydro_comm + time_communication_end - time_communication_start
#endif
         end if ! (xyzswp .eq. 1 .and. nprocs .gt. 1)

         if (xyzswp.eq.2 .and. nprocs .gt. 1) then
#ifndef DEBUG_TIMINGS
            call second_v(time_communication_start)
#endif
            ! exchange vzold with neighbors
            ! send down in phi direction


            sbuf3(1:config%qx,qy_s:qy_e) = vzold(1:config%qx,qy_s:qy_e,qz_s)
            call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
            sendcount = config%qx*qy_proc
            call MPI_SendRecv(sbuf3, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep1, rbuf3, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vzold(1:config%qx,qy_s:qy_e,qz_e+1) = rbuf3(1:config%qx,qy_s:qy_e)

            ! send up in phi direction
            sbuf3(1:config%qx,qy_s:qy_e) = vzold(1:config%qx,qy_s:qy_e,qz_e)
            call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
            sendcount = config%qx*qy_proc
            call MPI_SendRecv(sbuf3, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep2, rbuf3, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep2,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vzold(1:config%qx,qy_s:qy_e,qz_s-1) = rbuf3(1:config%qx,qy_s:qy_e)

            if (use_1neighbour_comm) call mpi_comm_lb_ub_1neighbour
            if (use_2neighbour_comm) call mpi_comm_lb_ub_2neighbour
            if (use_4neighbour_comm) call mpi_comm_lb_ub_4neighbour
#ifndef DEBUG_TIMINGS
            call second_v(time_communication_end)
            childrentime = childrentime + time_communication_end - time_communication_start
            timer%hydro_comm = timer%hydro_comm + time_communication_end - time_communication_start
#endif


         end if ! (xyzsweep .eq. 2)

         if (config%nsdim .eq. 3) then
         if (xyzswp .eq. 3 .and. nprocs .gt. 1 ) then
            ! do MPI next neigbor communication because of access in getrwz
            ! boundary conditions for config%qy=1, n are set in bndry()

            ! exchange vyold with neighbors
            ! send left in theta direction
#ifndef DEBUG_TIMINGS
            call second_v(time_communication_start)
#endif
            sbuf2(1:config%qx,qz_s:qz_e) = vyold(1:config%qx,qy_s,qz_s:qz_e)
            call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
            sendcount = config%qx*qz_proc
            call MPI_SendRecv(sbuf2, sendcount, MPI_DOUBLE_PRECISION, &
                              dest, tag_sweep1, rbuf2, sendcount,     &
                              MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vyold(1:config%qx,qy_e+1,qz_s:qz_e) = rbuf2(1:config%qx,qz_s:qz_e)

            ! send right in theta direction
            sbuf2(1:config%qx,qz_s:qz_e) = vyold(1:config%qx,qy_e,qz_s:qz_e)
            call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
            sendcount = config%qx*qz_proc
            call MPI_SendRecv(sbuf2, sendcount, MPI_DOUBLE_PRECISION,  &
                              dest, tag_sweep2, rbuf2, sendcount,      &
                              MPI_DOUBLE_PRECISION, src, tag_sweep2,  &
                              MPI_COMM_WORLD, mpistat, ierr)
            vyold(1:config%qx,qy_s-1,qz_s:qz_e) = rbuf2(1:config%qx,qz_s:qz_e)

            if (use_1neighbour_comm) call mpi_comm_pb_kb_1neighbour
            if (use_2neighbour_comm) call mpi_comm_pb_kb_2neighbour
            if (use_4neighbour_comm) call mpi_comm_pb_kb_4neighbour
           
#ifndef DEBUG_TIMINGS
            call second_v(time_communication_end)
            childrentime = childrentime + time_communication_end - time_communication_start
            timer%hydro_comm = timer%hydro_comm + time_communication_end - time_communication_start
#endif

         end if ! (xyzswp .eq. 3 .and. nprocs .gt. 1 )
         endif ! config%nsdim

      endif ! use_mpi

      sweep_mode = timer%hydro_sweep_mode

      if ( areas%ny .ge. 4 ) then
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_start)
#endif
!$omp parallel do                                           &
!$omp private ( ij,i,j,hydrow_self, hydrow_children )       &
!$omp shared  ( n1,n2,xyzswp ) &
!$omp schedule (static)
#endif
         do ij = 1, n1*n2
            j = int( ( ij+n2-1 ) / n2 )
            i = ij - (j-1) * n2
            if     ( xyzswp .eq. 1 ) then

               call getrwx (i+qy_s-1,j+qz_s-1,1)
               call hydrow (i+qy_s-1,j+qz_s-1, hydrow_self, hydrow_children)
               call putrwx (i+qy_s-1,j+qz_s-1)

            else if ( xyzswp .eq. 2 ) then

               call getrwy (i,j+qz_s-1,1)
               call hydrow (i,j+qz_s-1, hydrow_self, hydrow_children)
               call putrwy (i,j+qz_s-1)

            else
               call getrwz (i,j+qy_s-1,1)
               call hydrow (i,j+qy_s-1, hydrow_self, hydrow_children)
               call putrwz (i,j+qy_s-1)

            end if

            hydrow_hlp1(sweep_mode,:,isweep,caller) = hydrow_hlp1(sweep_mode,:,isweep,caller) + hydrow_self
            hydrow_hlp2(sweep_mode,:,isweep,caller) = hydrow_hlp2(sweep_mode,:,isweep,caller) + hydrow_children

            childrentime = childrentime + hydrow_self
         end do
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
         call second_v(paralleltime_end)
#endif
#endif
         timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start
         timer%sweep_hydrow(sweep_mode,:,isweep,caller) = hydrow_hlp1(sweep_mode,:,isweep,caller)
         timer%sweep_hydrow_children(sweep_mode,:,isweep,caller) = hydrow_hlp2(sweep_mode,:,isweep,caller)
      else
!     areas%ny < 4 cannot happen when -DMPI_HYDRO is used

         i = qy_s
         j = qz_s

         if     ( xyzswp .eq. 1 ) then
            call getrwx (i,j,1)
            call hydrow (i,j, hydrow_self, hydrow_children)
            call putrwx (i,j)
!     the following branch alternatives are irrelevant,
!     since areas%ny<4 implies areas%ny=1 in PROMETHEUS
!         else if ( xyzswp .eq. 2 ) then
!            call getrwy (i,j,1)
!            call hydrow (i,j)
!            call putrwy (i,j)
!         else
!            call getrwz (i,j,1)
!            call hydrow (i,j)
!            call putrwz (i,j)
            timer%sweep_hydrow(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_hydrow(timer%hydro_sweep_mode,:,isweep,caller)+hydrow_self
            timer%sweep_hydrow_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_hydrow_children(timer%hydro_sweep_mode,:,isweep,caller)+hydrow_children
            childrentime = childrentime + hydrow_self
         end if
      end if


!
!-------  update grid (only required in case of moving grid)
!
      if     ( xyzswp .eq. 1 ) then
         do i = 1, areas%nx
            i4 = i + 4
            xznl(i) = xl(i4)
            xzn (i) = x (i4)
            xznr(i) = xr(i4)
            dvx (i) = dvol(i4)
         end do

      endif

      if( xyzswp .eq. 2 .and. .not.(use_mpi)) then
         do j = 1, areas%ny
            j4 = j + 4
            yznl(j) = xl(j4)
            yzn (j) = x (j4)
            yznr(j) = xr(j4)
         end do
      endif

      if (xyzswp .ne. 1 .and.  xyzswp .ne. 2 .and. .not.(use_mpi)) then
         do k = 1, areas%nz
            k4 = k + 4
            zznl(k) = xl(k4)
            zzn (k) = x (k4)
            zznr(k) = xr(k4)
         end do
      end if


#ifdef CONTROL_OUTPUT
!-----------------------------------------------------------------------
!     control output:
!-----------------------------------------------------------------------

      k = 1
      write(*,*)"Task ",myproc
      write(*,*)'BBB: r/densty/ vx / vy / T / Ye /gpot'



        do j = qy_s, qy_e

           write(*,'('' j: '',i3,'' Y = '',1pe12.5)') j,yzn(j)
           do i = 1, areas%nx
              write(*,'(1x,i3,7(1x,1pe12.5))') i,xzn(i),  &
                    densty(i,j,k),velx(i,j,k),vely(i,j,k),&
                    temp(i,j,k),xnuc(i,j,k,config%qn),gpot(i,j,k)
           enddo
        enddo


!
!-----------------------------------------------------------------------
#endif /* CONTROL_OUTPUT */

!-------  compute gravitational potential and the modification of the
!         velocity and energy due to the action of gravity
!
!-----------------------------------------------------------------------
!         NOTE: 
!
!         POISON needs to be called only once for 3D problems, ACCELY
!         needs only to be called for 2-D problems, and ACCELZ needs 
!         not to be called at all, because the gravitational potential
!         is either assumed to be spherically symmetric (in 1D and 3D 
!         problems) or axially symmetric (in 2D problems).
!
!         WARNING: Program modifications are necessary, if Cartesian
!                  or cylindrical coordinates are used!
!-----------------------------------------------------------------------
      !
      !  in PNS poisson is called by promet
      !
!      if (xyzswp.eq.1  .or.  (xyzswp.eq.2 .and. config%nsdim.eq.2)) then
!c         i_grv = 1
!
!         if(i_grv .eq. 0) then
!c            call poisson
!         else
!c           call ppm_grv
!         end if
!      end if
!-----------------------------------------------------------------------
!      igrav=2
!      call poisson
!      igrav=0

!1D      if (xyzswp .eq. 1) then
!1D         call accelx (1, 1)
!1D      end if


     endif ! hydro_calculation

     if (accel) then

        if ( areas%ny.ge.4 ) then
           do ij = 1, n1*n2
              jj = int( ( ij+n2-1 ) / n2 )
              ii = ij - (jj-1) * n2
              if (xyzswp .eq. 1) then
#ifdef HTCL

                 call htcl(ii+qy_s-1,jj+qz_s-1)
#endif


                 call accelx (ii+qy_s-1, jj+qz_s-1, accelx_self, accelx_children)
                 timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)+accelx_self
                 timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)+accelx_children
                 childrentime = childrentime + accelx_self

                 if (config%use_flash_c .eqv. .true. .or. config%use_flash_o .eqv. .true. &
                     .or. config%use_flash_si .eqv. .true.) then

                    call burn(ii+qy_s-1,jj+qz_s-1, burn_self, burn_children)
                    timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)+burn_self
                    timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)+burn_children
                    childrentime = childrentime + burn_self

                 endif

#if  defined(BURN_SS) || defined(BURN_NETWORK) || defined(BURN_NETW_NOA)
                 call burn(ii,jj, burn_self, burn_children)
                 timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)+burn_self
                 timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)+burn_children
                 childrentime = childrentime + burn_self
#endif

#ifdef MIX
                 call mix(ii,jj)
#endif
#ifdef NSMIX
                 call nsmix(ii,jj)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! if multi-D transport or gravity will be implemtend accely must be called here in 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              else if (xyzswp .eq. 2  .and.  config%nsdim .eq. 2) then
                 call accely (ii, jj, accely_self, accely_children)
                 timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)+accely_self
                 timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                        timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)+accely_children
                 childrentime = childrentime + accely_self
!           else
!c-off         call accelz (ii, jj)
              end if
           end do
        else
           if (xyzswp .eq. 1) then

              ii = qy_s
              jj = qz_s


#ifdef HTCL
              call htcl(ii,jj)
#endif

              call accelx (ii, jj, accelx_self, accelx_children)
              timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_accel(timer%hydro_sweep_mode,:,isweep,caller)+accelx_self
              timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_accel_children(timer%hydro_sweep_mode,:,isweep,caller)+accelx_children
              childrentime = childrentime + accelx_self

              if (config%use_flash_c .eqv. .true. .or. config%use_flash_o .eqv. .true. &
                  .or. config%use_flash_si .eqv. .true.) then
                 call burn(ii,jj, burn_self, burn_children)
                 timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)= &
                    timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)+burn_self
                 timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                     timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)+burn_children
                 childrentime = childrentime + burn_self
              endif


#if defined(BURN_SS) || defined(BURN_NETWORK)|| defined(BURN_NETW_NOA)

              call burn(ii,jj, burn_self, burn_children)
              timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_burn(timer%hydro_sweep_mode,:,isweep,caller)+burn_self
              timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_burn_children(timer%hydro_sweep_mode,:,isweep,caller)+burn_children
              childrentime = childrentime + burn_self
#endif

#ifdef MIX
              call mix(ii,jj)
#endif
#ifdef NSMIX
            call nsmix(ii,jj)
#endif
         end if
      end if

#ifdef CONTROL_OUTPUT
!-----------------------------------------------------------------------
!     control output:
!-----------------------------------------------------------------------

      k = 1
      write(*,*) "Task ",myproc
      write(*,*) 'CCC: r/densty/ vx / vy / T / Ye /gpot'



        do j = qy_s,qy_e

           write(*,'('' j: '',i3,'' Y = '',1pe12.5)') j,yzn(j)
           do i = 1, 1!nx
              write(*,'(1x,i3,7(1x,1pe12.5))') i,xzn(i),   &
                    densty(i,j,k),velx(i,j,k),vely(i,j,k), &
                    temp(i,j,k),xnuc(i,j,k,config%qn),gpot(i,j,k) 
           enddo
        enddo
        call printit_taskX(0," ")
        !-----------------------------------------------------------------------
#endif

!      if (xyzswp .eq. 1) then
!
!-------  call equation of state (driver)     
        call eos3d (2,ler, eos3d_self, eos3d_children)
        timer%sweep_eos3d(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_eos3d(timer%hydro_sweep_mode,:,isweep,caller)+eos3d_self
        timer%sweep_eos3d_children(timer%hydro_sweep_mode,:,isweep,caller)= &
                timer%sweep_eos3d_children(timer%hydro_sweep_mode,:,isweep,caller)+eos3d_children
        childrentime = childrentime + eos3d_self

        if (ler) then
           write(*,*) "Task ",myproc,' sweep> isweep= ',isweep
           !         call stopit('sweep> eos failed',0) 
           raise_abort("sweeps(): eos failed")

        endif

!      endif

     endif !accel

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
     return
   end subroutine sweep


!=======================================================================
!>
!> \verbatim
!> copy 3D arrays in 1D array
!> \endverbatim
!>
!> \author W. Keil
!> \param j
!> \param k
!> \param isw
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
   subroutine getrwx ( j,k,isw )

     use precision

     use totgrq_hy
     use intgrs_hy

     use vnew_hy
     use vold_hy
     use mesh_hy
     use gfloat_hy
     use lfloat_hy
     use hydro_hy
     use grd_hy
     use physcs_hy
     use vnuw_hy

#ifdef KEIL_EPHI
     use ephhlp_hy
#endif
     use bndry_mod
     use ppm

     use mo_mpi

     use hydro_areas_mod
     use configure
     use state

     implicit none
! LOCAL variables that are not in modules

     integer(kind=ik) :: j1,j,ISW,i,n,j2,k1,k2,i4,k
     j1 = max0 (j - 1, 1)
     k1 = max0 (k - 1, 1)

     j2 = min0 (j + 1, max(qy_s,areas%ny))
     k2 = min0 (k + 1, max(qz_s,areas%nz))

!-------  transverse velocities uttp and utbt are given temporary
!         values in zones outside of grid (j=1, j=areas%ny, k=1, k=areas%ny).
!         they will be reset in subroutine bndry.
!
     do i = 1, areas%nx
        i4 = i + 4
        rho  (i4) = densty(i, j, k)
        u    (i4) = velx  (i, j, k)
        ut   (i4) = vely  (i, j, k)
        utt  (i4) = velz  (i, j, k)
        uttp (i4) = vyold (i,j2, k)
        utbt (i4) = vyold (i,j1, k)
        utrt (i4) = vzold (i, j,k2)
        utlt (i4) = vzold (i, j,k1)
        e    (i4) = energy(i, j, k)
        p    (i4) = press (i, j, k)
        tmp  (i4) = temp  (i, j, k)
        grav (i4) = 0.5_rk * (gpot(i,j,k) + gpot(i,j-1,k))
        game (i4) = gammae(i, j, k)
        gamc (i4) = gammac(i, j, k)
        ugrid(i4) = ugridx(i)
        xl   (i4) = xznl  (i)
        x    (i4) = xzn   (i)
        xr   (i4) = xznr  (i)
        dx   (i4) = xr(i4) - xl(i4)
#ifdef KEIL_EPHI
        ephi (i4) = ephtot(i)
#endif
     end do

#ifdef KEIL_EPHI
     ephi(4) = ephtot(0)
     ephi(3) = ephtot(0)
     ephi(2) = ephtot(0)
     ephi(1) = ephtot(0)
     ephi(0) = ephtot(0)

     ephi(areas%nx + 5) = ephtot(areas%nx)
     ephi(areas%nx + 6) = ephtot(areas%nx)
     ephi(areas%nx + 7) = ephtot(areas%nx)
     ephi(areas%nx + 8) = ephtot(areas%nx)
#endif

     do n = 1, config%qn
        do i = 1, areas%nx
           i4 = i + 4
           xn(i4,n) = xnuc(i,j,k,n)
        end do
     end do

     ybot = yzn(j1)
     ytop = yzn(j2)
     zlft = zzn(k1)
     zrgt = zzn(k2)

! ... set potential in the center
!
     grav(4) = 0.5_rk * (gpot(0,j,k) + gpot(0,j-1,k))
     call bndry (j, k)
     if ( isw.ne.0 ) then
        call geom  (j, k)
        call force (j, k)
     end if

     return
   end subroutine getrwx

!=======================================================================
!>
!> \verbatim
!> copy 3D arrays in 1D array
!> \endverbatim
!>
!> \author W. Keil
!> \param i x-index
!> \param k z-index
!> \param isw
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
   subroutine getrwy ( i,k,isw )

     use precision
     use totgrq_hy
     use intgrs_hy

     use vnew_hy
     use vold_hy
     use mesh_hy
     use gfloat_hy
     use lfloat_hy
     use hydro_hy
     use grd_hy
     use physcs_hy

#ifdef KEIL_EPHI
     use ephhlp_hy
#endif

     use bndry_mod
     use ppm

     use mo_mpi

     use hydro_areas_mod
     use configure
     use state
     implicit none
! LOCAL variables that are not in modules

     integer(kind=ik) :: i1,i,i2,k1,k,k2,j,j4,joff,n,isw, nloop
     i1 = max0 (i - 1, 1)
     i2 = min0 (i + 1, areas%nx)
     k1 = max0 (k - 1, 1)

     k2 = min0 (k + 1, max(qz_s,areas%nz))


!-------  transverse velocities uttp and utbt are given temporary
!         values in zones outside of grid (i=1, i=areas%nx, k=1, k=areas%ny).
!         they will be reset in subroutine BNDRY.
!

     if (use_mpi) then
        nloop = qy_proc
     else
        nloop = areas%ny
     endif

     do j = 1, nloop
        j4 = j + 4
   ! MPI_HYDRO:  j instead of joff for array access
        joff = qy_s + j - 1
        rho  (j4) = densty( i,joff, k)
        u    (j4) = vely  ( i,joff, k)
        ut   (j4) = velx  ( i,joff, k)
        utt  (j4) = velz  ( i,joff, k)
        uttp (j4) = vxold (i2,joff, k)
        utbt (j4) = vxold (i1,joff, k)
        utrt (j4) = vzold ( i,joff,k2)
        utlt (j4) = vzold ( i,joff,k1)
        e    (j4) = energy( i,joff, k)
        p    (j4) = press ( i,joff, k)
        tmp  (j4) = temp  ( i,joff, k)
        grav (j4) = 0.5_rk * (gpot(i,joff,k) + gpot(i-1,joff,k))
        game (j4) = gammae( i,joff, k)
        gamc (j4) = gammac( i,joff, k)

! MPI-Kommunikation in sweeps.F      
        ugrid(j4) = 0._rk    

        xl   (j4) = yznl  (joff)
        x    (j4) = yzn   (joff)
        xr   (j4) = yznr  (joff)


        dx   (j4) = xr(j4) - xl(j4)
#ifdef KEIL_EPHI
        ephi (j4) = ephtot(i)
#endif
     end do

#ifdef KEIL_EPHI
     ephi(4) = ephtot(i)
     ephi(3) = ephtot(i)
     ephi(2) = ephtot(i)
     ephi(1) = ephtot(i)
     ephi(0) = ephtot(i)

     ephi(areas%ny + 5) = ephtot(i)
     ephi(areas%ny + 6) = ephtot(i)
     ephi(areas%ny + 7) = ephtot(i)
     ephi(areas%ny + 8) = ephtot(i)
#endif

     do n = 1, config%qn
        do j = 1, nloop

           j4 = j + 4
           joff = qy_s + j - 1
           xn(j4,n) = xnuc(i,joff,k,n)
        end do
     end do

     xbot = xzn(i1)
     xtop = xzn(i2)
     zlft = zzn(k1)
     zrgt = zzn(k2)

! ... set potential in the center

     grav(4) = 0.5_rk * (gpot(i,qy_s-1,k) + gpot(i-1,qy_s-1,k))

     call bndry (i, k)

     if ( isw.ne.0 ) then
        call geom  (i, k)
        call force (i, k)
     end if

     return
  end subroutine getrwy

!=======================================================================
!>
!> \verbatim
!> copy 3D arrays in 1D array
!> \endverbatim
!>
!> \author W. Keil
!> \param i x-index
!> \param j y-index
!> \param isw
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
  subroutine getrwz ( i,j,isw )
!=======================================================================
    use precision
    use totgrq_hy
    use intgrs_hy

    use vnew_hy
    use vold_hy
    use mesh_hy
    use gfloat_hy
    use lfloat_hy
    use hydro_hy
    use grd_hy
    use physcs_hy
    use bndry_mod
    use ppm

    use mo_mpi

    use hydro_areas_mod
    use configure
    use state

    implicit none
! LOCAL variables that are not in modules

    integer(kind=ik) :: i,j,k,i1,i2,j1,j2,k4,n,isw,koff

    i1 = max0 (i - 1, 1)
    i2 = min0 (i + 1, areas%nx)

    j1 = max0 (j - 1, 1)

    j2 = min0 (j + 1, max(qy_s,areas%ny))


!-------  transverse velocities uttp and utbt are given temporary
!         values in zones outside of grid (i=1, i=areas%nx, k=1, k=areas%ny).
!         they will be reset in subroutine BNDRY.
    !

    do k = 1, qz_proc

       k4 = k + 4

       koff = qz_s + k - 1
       rho  (k4) = densty( i,j, koff)
       u    (k4) = velz  ( i,j, koff)
       ut   (k4) = velx  ( i,j, koff)
       utt  (k4) = vely  ( i,j, koff)
       uttp (k4) = vxold (i2,j, koff)
       utbt (k4) = vxold (i1,j, koff)
       utrt (k4) = vyold ( i,j2,koff)
       utlt (k4) = vyold ( i,j1,koff)
       e    (k4) = energy( i,j, koff)
       p    (k4) = press ( i,j, koff)
       tmp  (k4) = temp  ( i,j, koff)
       grav (k4) = 0.5_rk * (gpot(i,j,  koff) + gpot(i-1,j,  koff) )
       game (k4) = gammae( i,j, koff)
       gamc (k4) = gammac( i,j, koff)

       ugrid(k4) = 0._rk        

       xl   (k4) = zznl  (koff)
       x    (k4) = zzn   (koff)
       xr   (k4) = zznr  (koff)

       dx   (k4) = xr(k4) - xl(k4)
    end do

    do n = 1, config%qn

       do k = 1, qz_proc

          koff = qz_s + k - 1
          k4 = k + 4

          xn(k4,n) = xnuc(i,j,koff,n)

       end do
    end do

    xbot = xzn(i1)
    xtop = xzn(i2)
    ylft = yzn(j1)
    yrgt = yzn(j2)
!
! ... set potential in the center
!

    grav(4) = 0.5_rk * (gpot(i,j,  qz_s) + gpot(i-1,j,  qz_s) )


    call bndry (i, j)

    if ( isw.ne.0 ) then 
       call geom  (i, j)
       call force (i, j)
    end if

    return
  end subroutine getrwz

#ifdef IBM_COMPILER
! this is necessary on the POWER6 since the xlf90 version 10.xx
! causes problems
@PROCESS NOHOT
#endif
!=======================================================================
!>
!> \verbatim
!> perform hydro step on one row of zones
!> \endverbatim
!>
!> \author W. Keil
!> \param j
!> \param k
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine hydrow (j, k, selftime, childrentime)

  use precision

  use intgrs_hy

  use reman_hy
  use vnew_hy
  use vold_hy
  use mesh_hy
  use gfloat_hy
  use lfloat_hy
  use hydro_hy
  use grd_hy
  use physcs_hy
  use ppm

#ifdef ENECONS
  use param_rt
#endif
#if defined(CONVECTION_1D) && !defined(CFC_TRANSPORT)
  use mixinglength
#endif


  use mo_mpi
  use bndry_mod
  use specfun, only : fastsqrt

  use hydro_areas_mod
  use configure
  use cputim
  use state
  implicit none


! LOCAL variables that are not in modules
  real(kind=rk), intent(out)  :: selftime(2), childrentime(2)
  real(kind=rk)               :: selftime_start(2)
  integer(kind=ik)            :: i, j, nuc0, n, NUCBEG, K, i4
  real(KIND=RK)               :: ugmax, aux1, aux2, hxnnu, hunu, eint, &
                                 ekin, henu
  REAL(KIND=RK)               :: huttnu, hutnu
  real(KIND=RK)               :: xlold(config%q), dvold(config%q), &
                                 aold(config%q), scrch(config%q),  &
                                 fict2(config%q), wrksm(config%q)

#ifdef ENECONS
  integer(KIND=IK)            :: iene_sw
#endif

  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

#ifdef ENECONS
  iene_sw=1-config%laghyd
#endif

#ifdef CONTROL_OUTPUT_2
  write(*,*) "Task ",myproc
  write(*,*) 'hydrow 1: x/ dx/ rho/ u/ ce/ dtdx/ game/ e/ ephi/ grav'
#endif

  do i = 1, nzn8
     c(i) = gamc(i) * p(i) * rho(i)
  enddo

  call fastsqrt(c,nzn8)


  do i = 1, nzn8
     v    (i) = 1._rk / rho(i)
     ! c    (i) = sqrt (gamc(i) * p(i) * rho(i))
     ce   (i) = c(i) / rho(i)
     !        dtdx (i) = hydro%dt / dx(i) * 0.5 * (ephi(i) + ephi(i-1))
     dtdx (i) = hydro%dt / dx(i) 
     fict2(i) = 0._rk

#ifdef CONTROL_OUTPUT_2
     write(*,*) "Task ",myproc
     write(*,'(1x,i3,9(1x,1pe12.5))') &
          i,x(i),dx(i),rho(i),u(i),ce(i),dtdx(i),game(i),e(i), &
          grav(i)
#endif
  enddo

  do i = 2, nzn8
     ugrdl(i) = 0.5_rk * (ugrid(i) + ugrid(i-1))
  enddo
  ugrdl(1) = ugrdl(2)

!-----------------------------------------------------------------------
!     take into account non Cartesian coordinates:
!     if igeom = 3, 4, or 5, the x coordinate will be a radial
!     coordinate, ie. the index j will always refer to the
!     x-direction.
!-----------------------------------------------------------------------

  if (igeom .ge. 3)   then
     do i = 1, nzn8
        dtdx(i) = dtdx(i) / xzn(j)                                  
     enddo
  endif

#if 0
  if(config%nsdim .eq. 2 .and. config%igeomz .eq. 5) then
     if(xyzswp .eq. 1) then           ! radial sweep -> u = u_rad
        do i = 1, nzn8
           fict2(i) = - u(i)*utt(i)/x(i)
        enddo
     else                             ! theta sweep  -> u = u_theta
        do i = 1, nzn8
           !               fict2(i) = - u(i)*utt(i)*cot(x(i))/xzn(j)
           fict2(i) = - u(i)*utt(i)/(tan(x(i))*xzn(j))
        enddo
     endif
  endif
#endif


!-----------------------------------------------------------------------
!     Riemann-solver:
!-----------------------------------------------------------------------

  call intrfc
  call states  ( j )
  call riemann 

!-----------------------------------------------------------------------
#ifdef ENECONS
  if (xyzswp .eq. 1) then
     do i = 1, nzn
        i4 = i + 4
        grav (i4) = 0.5_rk * (gpot(i,j,k) + gpot(i,j-1,k))
     enddo

     grav (4) = 0.5 * (gpot(0,j,k) + gpot(0,j-1,k))

  else if (xyzswp.eq.2) then
     do i = 1, nzn
        i4 = i + 4

        grav (i4) = 0.5 * (gpot(j,i+qy_s-1,k) + gpot(j-1,i+qy_s-1,k))

     end do
     if (use_mpi) then
        grav (4) = 0.5 * (gpot(j-1,i,k) + gpot(j,i,k))
     else

     endif


  end if
#endif /* ENECONS */
!
!      do 200  i = 5, nzn4
!         dxzn      = xznr(i-4) - xznl(i-4)
!         dphidx(i) = (grav (i) - grav (i-1)) / dxzn 
! 200  continue
!-----------------------------------------------------------------------

!
!------  compute fluxes (using the solution of the Riemann problem)
!
  do i = 5, nzn5
!        ucor  = hydro%dt * 0.5 * (ephi(i) + ephi(i-1)) *  dphidx(i) * 0.0015
!        urell (i) = urell(i) - ucor
!        uav   (i) = uav(i) - ucor
     rhoflx(i) = rhoav (i) * urell(i)
     uflx  (i) = rhoflx(i) * uav  (i)
     utflx (i) = rhoflx(i) * utav (i)
     uttflx(i) = rhoflx(i) * uttav(i)
     eflx  (i) = rhoflx(i) * ( pav(i) / (rhoav(i) * (gameav(i)-1._rk)) &
                 + 0.5_rk * (uav(i)**2 + utav(i)**2 + uttav(i)**2) )&
                 + uav(i) * pav(i)
  enddo


  nuc0 = 1
  do n = 1, config%qn-3, 4
     do i = 5, nzn5
        xnflx(i,n  ) = xnav(i,n  ) * rhoflx(i)
        xnflx(i,n+1) = xnav(i,n+1) * rhoflx(i)
        xnflx(i,n+2) = xnav(i,n+2) * rhoflx(i)
        xnflx(i,n+3) = xnav(i,n+3) * rhoflx(i)
     end do
     nuc0 = n + 4
  end do

  nucbeg = nuc0
  do n = nucbeg, config%qn-1, 2
     do i = 5, nzn5
        xnflx(i,n  ) = xnav(i,n  ) * rhoflx(i)
        xnflx(i,n+1) = xnav(i,n+1) * rhoflx(i)
     end do
     nuc0 = n + 2
  end do

  nucbeg = nuc0
  do n = nucbeg, config%qn
     do i = 5, nzn5
        xnflx(i,n  ) = xnav(i,n  ) * rhoflx(i)
     end do
  end do


  if (config%cvisc .gt. 0.0)  call avisco (j, k)

#if defined(CONVECTION_1D) && !defined(CFC_TRANSPORT)
  if(xyzswp .eq. 1) then
     gamcl(6:nzn4)=0.5_rk*(gamc(5:nzn4-1)+gamc(6:nzn4))
    call mixinglengthconvection(&
      CONV_PARAM_1, &
      x(5:nzn4), &
      p(5:nzn4), pav(6:nzn4), &
      rho(5:nzn4), rhoav(6:nzn4), &
      gamc(5:nzn4), &
      gamcl(6:nzn4), &
      grav(5:nzn4), &
      u(5:nzn4), &
      xn(5:nzn4,:), &
      xnflx(6:nzn4,:), &
      e(5:nzn4), &
      eflx(6:nzn4))
  endif
#endif

  do i = 5, nzn5
     rhoflx(i) = rhoflx(i) * areal(i)
     uflx  (i) = uflx  (i) * areal(i)
     utflx (i) = utflx (i) * areal(i)
     uttflx(i) = uttflx(i) * areal(i)
     eflx  (i) = eflx  (i) * areal(i)
  end do

  if (xyzswp .eq. 1) then
     do i = 5, nzn5
        uttflx(i) = uttflx(i)*xl(i)
     end do
  else if (xyzswp .eq. 2) then
     do i = 5, nzn5
        uttflx(i) = uttflx(i)*sin(xl(i))
     end do
  end if

  do n = 1, config%qn
     do i = 5, nzn5
        xnflx(i,n) = xnflx(i,n) * areal(i) 
     end do
  end do

!-------  save grid quantities 

  do i = 1, nzn8
     xlold(i) = xl  (i)
     dvold(i) = dvol(i)
     aold (i) = area(i)
  end do

!-------  update grid

  do i = 2, nzn8
     xl(i  ) = xlold(i) + hydro%dt * ugrdl(i)
     xr(i-1) = xl   (i)
  end do

  do i = 2, nzn7
     dx(i) = xr(i) - xl(i)
     x (i) = 0.5_rk * (xr(i) + xl(i))
  end do

  call geom (j, k)


!------- update conserved quantities

  ugmax = 0.0_rk
  do i = 1, nzn8
     ugmax     = max (ugmax, abs( ugrid(i) ), abs(ugrdl(i) ) )
     scrch(i) = 1._rk
  end do

  if (ugmax .ne. 0._rk)  then
     do i = 5, nzn4
        scrch(i) = dvold(i) / dvol(i)
     enddo
  end if


  do i = 5, nzn4
     !        dtdx (i) = hydro%dt / dvol(i) * 0.5 * (ephi(i) + ephi(i-1))
     dtdx (i) = hydro%dt / dvol(i) 
     rhonu(i) = rho(i)*scrch(i) - dtdx(i) * (rhoflx(i+1)-rhoflx(i))
     rhonu(i) = max(config%smlrho, rhonu(i))
  end do





  do n = 1, config%qn
!pp_IVL
!pp_PVL
     do i = 5, nzn4
        aux1  = - dtdx(i) * (xnflx(i+1,n) - xnflx(i,n))
        aux2  = abs (aux1)  + abs ( scrch(i) * (rhonu(i) - rho(i)) )
        hxnnu = (rho(i) * xn(i,n) * scrch(i) + aux1) / rhonu(i)
        xnnu(i,n) = hxnnu
        if ( aux2.eq.0._rk ) xnnu(i,n) = xn(i,n)
     end do
  end do

#ifdef CMA

!     ---------
!     CMA start
!     ---------

!     ----------------------------------------------------
!     Attention: Only densities of nuclei are renormalized. 
!                Ye should not be changed at all (i.e. use 
!                index config%qn-1, instead of config%qn)!!!!

  do n = 1, config%qn-1
!pp_IVL
!pp_PVL
     do i = 5, nzn4
        xnnu(i,n) = rhonu(i) * xnnu(i,n) 
     end do
  end do

!pp_IVL
  do i = 5,nzn4
     wrksm(i) = 0._rk
  end do

  do n = 1,config%qn-1
!pp_IVL
!pp_PVL
     do i = 5,nzn4
        wrksm(i) = wrksm(i) + xnnu(i,n)
     end do
  end do

!pp_IVL
  do i = 5,nzn4
     rhonu(i) = max( config%smlrho , wrksm(i) )
     wrksm(i) = 1._rk/rhonu(i)
  end do

  do n = 1,config%qn-1
!pp_IVL
!pp_PVL
     do i = 5,nzn4
        xnnu(i,n) = xnnu(i,n) * wrksm(i)
     end do
  end do

!     -------
!     CMA end
!     -------
#endif /* CMA */


  do i = 5, nzn4

!        aux1 = - dtdx(i) * (uflx(i+1) - uflx(i) +
!     &           aold(i) * (pav(i+1) - pav(i)))
!     &         + hydro%dt * 0.5 * (rhonu(i) + rho(i)) * fict(i)
!     &         * 0.5 * (ephi(i) + ephi(i-1))
     aux1 = - dtdx(i) * (uflx(i+1) - uflx(i) + &
              aold(i) * (pav(i+1) - pav(i))) &
              + hydro%dt * 0.5_rk * (rhonu(i) + rho(i)) * fict(i)
     aux2 = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
     hunu = (rho(i) * u(i) * scrch(i) + aux1)  /  rhonu(i)

     if (aux2 .eq. 0._rk)  then          
        unu(i) = u(i)
     else
        unu(i) = hunu
     end if

     aux1  = - dtdx(i) * (utflx(i+1) - utflx(i))
     aux2  = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
     hutnu = (rho(i) * ut(i) * scrch(i) + aux1) / rhonu(i)

     if (aux2 .eq. 0._rk)  then          
        utnu(i) = ut(i)
     else
        utnu(i) = hutnu
     end if

  end do

  !Special treatment for uttnu (which is the component v_varphi in
  !case of the x- and y-sweep): effectively, we solve an advection
  !equation of the angular momentum density rho*v_varphi*r*sin(theta)

  if (xyzswp .eq. 1) then
     !We're transporting r*rho*v_varphi
     do i = 5, nzn4

#if 0
        aux1   = - dtdx(i) * (uttflx(i+1) - uttflx(i)) &
              + hydro%dt * 0.5_rk * (rhonu(i) + rho(i)) * fict2(i)
#else
        aux1   = - dtdx(i) * (uttflx(i+1) - uttflx(i)) / x(i)
#endif
        aux2   = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
        huttnu = (rho(i) * utt(i) * scrch(i) + aux1) / rhonu(i)

        if (aux2 .eq. 0._rk)  then          
           uttnu(i) = utt(i)
        else
           uttnu(i) = huttnu
        end if

     end do

  else if (xyzswp .eq. 2) then
     !We're transporting sin_theta*rho*v_varphi
     do i = 5, nzn4

#if 0
        aux1   = - dtdx(i) * (uttflx(i+1) - uttflx(i)) &
              + hydro%dt * 0.5_rk * (rhonu(i) + rho(i)) * fict2(i)
#else
        aux1   = - dtdx(i) * (uttflx(i+1) - uttflx(i)) / sin(x(i))
#endif
        aux2   = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
        huttnu = (rho(i) * utt(i) * scrch(i) + aux1) / rhonu(i)

        if (aux2 .eq. 0._rk)  then          
           uttnu(i) = utt(i)
        else
           uttnu(i) = huttnu
        end if

     end do

  else if (xyzswp .eq. 3) then
     !nothing needs to be done in phi-sweep
     do i = 5, nzn4

        aux1   = - dtdx(i) * (uttflx(i+1) - uttflx(i))         
        aux2   = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
        huttnu = (rho(i) * utt(i) * scrch(i) + aux1) / rhonu(i)

        if (aux2 .eq. 0._rk)  then          
           uttnu(i) = utt(i)
        else
           uttnu(i) = huttnu
        end if

     end do

  end if


  do i = 5, nzn4

     aux1 = - dtdx(i) * (eflx(i+1) - eflx(i))

#ifdef ENECONS
     aux1 = aux1 - 0.5_rk *dtdx(i) * iene_sw * &
            (rhoflx(i+1)+rhoflx(i))*(grav(i)-grav(i-1))
#endif
     aux2 = abs (aux1)  +  abs (scrch(i) * (rhonu(i) - rho(i)))
     henu = (rho(i) * e(i) * scrch(i) + aux1) / rhonu(i)

     if (aux2 .eq. 0._rk)  then          
        enu(i) = e(i)
     else
        enu(i) = henu
     end if

     ekin   = 0.5_rk * (unu(i)**2 + utnu(i)**2 + uttnu(i)**2)
     eint   = enu(i) - ekin
     eint   = max (eint, config%smalle)
     enu(i) = eint + ekin

  end do
#ifndef DEBUG_TIMINGS
  call second_v(selftime)
  selftime = selftime - selftime_start
#endif
  return
end subroutine hydrow

#ifdef IBM_COMPILER
! this is necessary on the POWER6 since the xlf90 version 10.xx
! causes problems
@PROCESS HOT
#endif

!=======================================================================
!>
!> \verbatim
!> copy 1D arrays in 3D array
!> \endverbatim
!>
!> \author W. Keil
!> \param j
!> \param k
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine putrwx (j, k)

  use precision

  use intgrs_hy
  ! use arecon_hy
  use bndinf_hy
  use massio_hy

  use vnew_hy
  use hydro_hy
  use grd_hy
  use mesh_hy

  use hydro_areas_mod
  use state
  use configure
  IMPLICIT NONE
! LOCAL variables that are not in modules

  INTEGER(KIND=IK) :: I,I4,J,K,N
  do i = 1, areas%nx
     i4 = i + 4
     densty(i,j,k) = rhonu(i4)
     velx  (i,j,k) = unu  (i4)
     vely  (i,j,k) = utnu (i4)
     velz  (i,j,k) = uttnu(i4)
     energy(i,j,k) = enu  (i4)
     press (i,j,k) = p    (i4)
     temp  (i,j,k) = tmp  (i4)
     gammae(i,j,k) = game (i4)
     gammac(i,j,k) = gamc (i4)
  enddo

!      if (j .eq. areas%ny  .and.  k .eq. areas%nz) then
!         do i = 1, areas%nx
!            i4 = i + 4
!            xznl(i) = xl(i4)
!            xzn (i) = x (i4)
!            xznr(i) = xr(i4)
!
!            dvx (i) = dvol(i4)
!         end do
!      end if

  do n = 1, config%qn
     do i = 1, areas%nx
        i4 = i + 4
        xnuc(i,j,k,n) = xnnu(i4,n)
     enddo
  enddo

  if (are_id .eq. 1) then
     dmdtio(1,j,k) = rhoflx(5)
     dmdtio(3,j,k) = eflx(5)
     dmdtio(5,j,k) = xnflx(5,config%qn)
  end if

  if (are_id .eq. areas%are_nu) then
     dmdtio(2,j,k) = rhoflx(nzn5)
     dmdtio(4,j,k) = eflx(nzn5)
     dmdtio(6,j,k) = xnflx(nzn5,config%qn)
  end if


  dflxl  (j,k)      = rhoflx(5)
  dflxr  (j,k)      = rhoflx(nzn5)
  eflxl  (j,k)      = eflx(5)
  eflxr  (j,k)      = eflx(nzn5)
  vxflxl  (j,k)     = uflx(5)
  vxflxr  (j,k)     = uflx(nzn5)
  vyflxl (j,k)      = utflx(5)
  vyflxr (j,k)      = utflx(nzn5)
  if (are_id .eq. areas%are_nu .and. config%nsdim .gt. 1_ik) then
    vzflxl(j,k)     = uttflx(5)/xl(5)
  else
    vzflxl(j,k)     = 0._rk
  endif
  vzflxr(j,k)       = uttflx(nzn5)/xl(nzn5)
  xnflxl (j,k,1:config%qn) = xnflx(5,1:config%qn)
  xnflxr (j,k,1:config%qn) = xnflx(nzn5,1:config%qn)


  return
end subroutine putrwx

!=======================================================================
!>
!> \verbatim
!> copy 1D arrays in 3D array
!> \endverbatim
!>
!> \author W. Keil
!> \param i
!> \param k
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine putrwy (i, k)
!=======================================================================

  use precision
  use intgrs_hy
  ! use arecon_hy
  use bndinf_hy
  use massio_hy
  use vnew_hy
  use hydro_hy
  use grd_hy
  use mesh_hy

  use mo_mpi

  use hydro_areas_mod
  use state
  use configure
  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik) :: j,j4,i,k,n

  integer(kind=ik) :: joff


  do j = 1, qy_proc

     j4 = j + 4
   ! MPI_HYDRO  joff statt j
     joff = j+qy_s-1 
     densty(i,joff,k) = rhonu(j4)
     vely  (i,joff,k) = unu  (j4)
     velx  (i,joff,k) = utnu (j4) 
     velz  (i,joff,k) = uttnu(j4)
     energy(i,joff,k) = enu  (j4)
     press (i,joff,k) = p    (j4)
     temp  (i,joff,k) = tmp  (j4)
     gammae(i,joff,k) = game (j4)
     gammac(i,joff,k) = gamc (j4)


  end do

  do n = 1, config%qn

     do j = 1, qy_proc 

        j4 = j + 4

        joff = j+qy_s-1 
        xnuc(i,joff,k,n) = xnnu(j4,n)

     end do
  end do




  return
end subroutine putrwy

!======================================================================
!>
!> \verbatim
!> copy 1D arrays in 3D array
!> \endverbatim
!>
!> \author W. Keil
!> \param i
!> \param j
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine putrwz (i, j)
!======================================================================
  use precision

  use intgrs_hy
  ! use arecon_hy
  use bndinf_hy
  use massio_hy

  use vnew_hy
  use hydro_hy
  use grd_hy
  use mesh_hy

  use mo_mpi

  use hydro_areas_mod
  use state
  use configure
  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik) :: k,k4,i,j,n

  integer(kind=ik) :: koff



  do k = 1, qz_proc 

     k4 = k + 4

   ! MPI_HYDRO  koff statt k
     koff = k+qz_s-1 
     densty(i,j,koff) = rhonu(k4)
     velz  (i,j,koff) = unu  (k4)
     velx  (i,j,koff) = utnu (k4) 
     vely  (i,j,koff) = uttnu(k4)
     energy(i,j,koff) = enu  (k4)
     press (i,j,koff) = p    (k4)
     temp  (i,j,koff) = tmp  (k4)
     gammae(i,j,koff) = game (k4)
     gammac(i,j,koff) = gamc (k4)

  end do

  do n = 1, config%qn

     do k = 1, qz_proc 

        k4 = k + 4

        koff = k+qz_s-1 
        xnuc(i,j,koff,n) = xnnu(k4,n)

     end do
  end do

  return
end subroutine putrwz

!=======================================================================
!>
!> \verbatim
!> Purpose: compute change of x-velocity and energy and electron fraction
!>          due to the action of gravitational forces and neutrino transport
!>
!>
!>  Notes:   PNS-Version by wfk, Neutrino-Terms by mjr
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param jj
!> \param kk
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
subroutine accelx (jj, kk, selftime, childrentime)

!=======================================================================
  use precision

!      use bndinf_hy ! not used
  use intgrs_hy
  use totgrq_hy

  use vnew_hy
  use vold_hy
  use mesh_hy
  use gfloat_hy
  use hydro_hy
  use physcs_hy
  use spez_hy
  use vnuw_hy
  use hlpare_hy
  use nucparam

  use eos_sn2, ONLY : lsrolo
#ifdef ENECONS
  use param_rt
#endif

  use hydro_areas_mod
  use configure
  use cputim
  use state
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2)
  integer(kind=ik)           :: i, i4, jj, kk, nst
  real(KIND=RK)              :: dpn_max, dpn, dth, dxzn, &
                                dpn_min, dthe

#ifdef AVDAMP
  real(KIND=RK),save         :: velsav(400,5000), timsav(5000)
  real(KIND=RK)              :: velav(400)
  logical,save               :: firstavdamp=.true.
#endif
  real(KIND=RK)              :: dphold(config%q), &
                                grold(config%q),  &
                                dydt(config%q),   &
                                qegr(config%q)

#ifdef ENECONS
  real(KIND=RK)              :: scr
  integer(KIND=IK)           :: iene_sw

  iene_sw=1-config%laghyd
#endif

  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

      do i = 1, nzn
         i4 = i + 4
         grold(i4)    = 0.5_rk * (gold(i,jj,kk) + gold(i,jj-1,kk))
         grav (i4)    = 0.5_rk * (gpot(i,jj,kk) + gpot(i,jj-1,kk))
         u    (i4)    = vxvold(i,jj,kk)
         unu  (i4)    = velx  (i,jj,kk)
         enu  (i4)    = energy(i,jj,kk)
         rhonu(i4)    = densty(i,jj,kk)
         rho  (i4)    = denold(i,jj,kk)
         xnnu (i4,1)  = xnuc  (i,jj,kk,1)     !=X_n
         xnnu (i4,2)  = xnuc  (i,jj,kk,2)     !=X_p
         xnnu (i4,config%qn) = xnuc  (i,jj,kk,config%qn)    !=Y_e

         ! -- neutrino source terms
         dedt (i4) = qen(i,jj,kk)/densty(i,jj,kk)
         s    (i4) = qmo(i,jj,kk)/densty(i,jj,kk)
         dydt (i4) = qye(i,jj,kk,1)/densty(i,jj,kk)

      end do

      grold(4) = 0.5_rk * (gold(0,jj,kk) + gold(0,jj-1,kk))
      grav (4) = 0.5_rk * (gpot(0,jj,kk) + gpot(0,jj-1,kk))

      do i = 5, nzn4
         dxzn      = xznr(i-4) - xznl(i-4)
         dphold(i) = (grold(i) - grold(i-1)) / dxzn 
         dphidx(i) = (grav (i) - grav (i-1)) / dxzn 
      end do


!-------  compute new velocity and energy and net-electron fraction
!
!         N O T E:  
!
!         The following loop starts at i=7, if an 
!         inner boundary condition at finite radius is used.  A
!         parabolic velocity profile is assumed for UREL!
!
!         B U T:
!
!         This is not reasonable, if the inner boundary is
!         an interface between two calculation areas!
!
!-----------------------------------------------------------------------

      dth = 0.5_rk * hydro%dt

      nst = 5
!      if(are_id .eq. 1) nst = 7 ! not used


      do i = nst, nzn4
         i4 = i-4

!         dthe   = dth * 0.5 * (ephtot(i4) + ephtot(i4-1))
         dthe   = dth 
! gravi
         unu(i)     = unu(i) - dthe * ( dphidx(i) +  dphold(i))
         qegr(i)    = - (unu(i)*dphidx(i) + u(i)*dphold(i))
#ifdef ENECONS
         scr = (grav(i)+grav(i-1))-(grold(i)+grold(i-1))
         qegr(i)    = iene_sw*(rho(i)/rhonu(i)-1.0_rk)*scr+config%laghyd*qegr(i)
#endif
         enu(i)     = enu(i) + dthe * qegr(i)

! neutrino
#ifndef NEW_QMO
         unu(i)     = unu(i)     + s(i)    * hydro%dt
#endif
         enu(i)     = enu(i)     + dedt(i) * hydro%dt
         xnnu(i,config%qn) = xnnu(i,config%qn) + dydt(i) * hydro%dt

         if (config%restmass_version .eq. 3) then
            ! quick fix !!!! Should be done implicitly in neutrino transport
            enu(i) = enu(i) - mbar(config%qn) * dydt(i) * hydro%dt
         endif

         if (rhonu(i).lt.lsrolo) then
            dpn_max=min(  xnnu(i,1),1._rk-xnnu(i,2) )
            dpn_min=max( -xnnu(i,2),xnnu(i,1)-1._rk )
            dpn=min( dydt(i)*hydro%dt,dpn_max)
            dpn=max( dpn       ,dpn_min)

            if (config%restmass_version .gt. 0) then
               enu(i) = enu(i) + (mbar(n_n)-mbar(n_p))*dpn
            endif

!==================================================
! Attention: the next 2 lines were switched of in
! order to not regard the changes in composition
! when neutrinos are absorbed/emitted, because
! we cannot regard these changes for nuclei.
! so we also switch of nucleons. A.Marek
!
! This is only ok for collapse. In post-bounce it
! must be switched back on!
!==================================================

            xnnu(i,n_n) = xnnu(i,n_n) - dpn
            xnnu(i,n_p) = xnnu(i,n_p) + dpn
         endif

      end do


#ifndef ENECONS
#ifdef FIX_INNER
!c gravi & neutrino for two innermost zones
!         unu(6)     = (unu(7) - ugridx(3)) / 9.0  +  ugridx(2)
!         qegr(6)    = - (unu(6)*dphidx(6) + u(6)*dphold(6))
!         enu(6)     = enu(6) + dth * qegr(6) + dedt(6) * hydro%dt
!         xnnu(6,config%qn) = xnnu(6,config%qn) + dydt(6) * hydro%dt
!
!         unu(5)     = (unu(6) - ugridx(2)) / 9.0  +  ugridx(1)
!         qegr(5)    = - (unu(5)*dphidx(5) + u(5)*dphold(5))
!         enu(5)     = enu(5) + dth * qegr(5) + dedt(5) * hydro%dt
!         xnnu(5,config%qn) = xnnu(5,config%qn) + dydt(5) * hydro%dt
      dvdx=unu(7)/xzn(3)

      unu(6)     = xzn(2)*dvdx
      qegr(6)    = - (unu(6)*dphidx(6) + u(6)*dphold(6))
      enu(6)     = enu(6) + dth * qegr(6) + dedt(6) * hydro%dt
      xnnu(6,config%qn) = xnnu(6,config%qn) + dydt(6) * hydro%dt

      unu(5)     = xzn(1)*dvdx
      qegr(5)    = - (unu(5)*dphidx(5) + u(5)*dphold(5))
      enu(5)     = enu(5) + dth * qegr(5) + dedt(5) * hydro%dt
      xnnu(5,config%qn) = xnnu(5,config%qn) + dydt(5) * hydro%dt
#endif
#endif
! --- update GR-Potential energy sourceterm for output
      do i = 1, nzn
         i4 = i + 4
         qgrv(i,jj,kk) = qgrv(i,jj,kk) + 0.5_rk * qegr(i4)*rhonu(i4)
      end do

!-------  update velocity and energy, net-electron fraction
!            and mass fractions of free nucleons

      do i = 1, nzn
         i4 = i + 4
         velx  (i,jj,kk)    = unu(i4)
         energy(i,jj,kk)    = enu(i4)
         xnuc  (i,jj,kk, n_n) = xnnu(i4,n_n)
         xnuc  (i,jj,kk, n_p) = xnnu(i4,n_p)
         xnuc  (i,jj,kk,config%qn) = xnnu(i4,config%qn)
      end do

#ifdef DAMPNSV
! DAMP PNS OSCILLATIONS VIA UNDERRELAXATION
      idamp = 1
      do i = 1, nzn
         if (densty(i,jj,kk).gt.1e13) idamp = i
      end do
      dampdelta = 0.9_rk
      velx(1:idamp,jj,kk) = (1._rk - dampdelta) * velx(1:idamp,jj,kk) &
                          +       dampdelta  * vxvold(1:idamp,jj,kk)
#endif
#ifdef AVDAMP
      m4=min(nzn,400)
      if (firstavdamp) then
         firstavdamp=.false.
         timsav(:) = time
         do i=1,5000
            velsav(1:m4,i) = 0._rk !velx(1:m4,jj,kk)
         enddo
         goto 111
      endif
      idamp = 1
      do i = 1, m4
         if (densty(i,jj,kk).gt.1.e13_rk) idamp = i
      end do
      do i=2,5000
         if (time-timsav(i).lt.1.e-3_rk) itsav=i
      enddo
      if (timsav(itsav-1).eq.timsav(itsav)) goto 111
      velav(1:m4) = velsav(1:m4,1)*(time-timsav(1))
      do i=2,itsav
         velav(1:m4) = velav(1:m4) + &
                      velsav(1:m4,i)*(timsav(i-1)-timsav(i))
      enddo
      velav(:) = velav(:) / (time - timsav(itsav))

      damper=0.001_rk
      iti = mod(nint(time*1.e3_rk),100)

      if (iti.le.10) &
           velx(1:idamp,jj,kk) = velx(1:idamp,jj,kk) + &
           damper * ( velav(1:idamp) - velx(1:idamp,jj,kk) )

 111  velsav(1:m4,2:5000)=velsav(1:m4,1:4999)
      velsav(1:m4,1)     =velx(1:m4,jj,kk)

      timsav(2:5000)=timsav(1:4999)
      timsav(1) = time



#endif

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return
    end subroutine accelx


!=======================================================================
!>
!> \verbatim
!> Purpose: compute change of y-velocity and energy and electron fraction
!>          due to the action of gravitational forces and neutrino transport
!>
!>
!>  Notes:   PNS-Version by wfk, Neutrino-Terms by mjr (e^\Phi not active)
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param jj
!> \param kk
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
    subroutine accely (ii, kk, selftime, childrentime)

      use precision

      use intgrs_hy
      use totgrq_hy

      use vnew_hy
      use vold_hy
      use mesh_hy
      use gfloat_hy
      use hydro_hy
      use physcs_hy
      use spez_hy
      use vnuw_hy

#ifdef ENECONS
      use param_rt
#endif
      use mo_mpi

      use hydro_areas_mod
      use configure
      use cputim
      use state
      implicit none
! LOCAL variables that are not in modules

      real(kind=rk), intent(out) :: selftime(2), childrentime(2)
      real(kind=rk)              :: selftime_start(2)
      integer(kind=ik)           :: j, j4, joff, ii, kk
      real(KIND=RK)              :: dyzn, dth
      real(KIND=RK)              :: dphold(config%q),  &
                                    grold(config%q),   &
                                    qegr(config%q)


      selftime     = 0._rk
      childrentime = 0._rk

#ifndef DEBUG_TIMINGS
      call second_v(selftime_start)
#endif

!-------  WARNING:
!
!         For non Eulerian grids GROLD should be interpolated to the
!         positions where GRAV is evaluated
!
!
      do j = 1, nzn
        j4 = j + 4

! wompi j instead of j+qy_s-1
        joff = j+qy_s-1
        grold(j4)    = 0.5 * (gold(ii,joff,kk) + gold(ii-1,joff,kk))
        grav (j4)    = 0.5 * (gpot(ii,joff,kk) + gpot(ii-1,joff,kk))
        ut   (j4)    = vyvold(ii,joff,kk)
        utnu (j4)    = vely  (ii,joff,kk)
        enu  (j4)    = energy(ii,joff,kk)
        rhonu(j4)    = densty(ii,joff,kk)
        rho  (j4)    = denold(ii,joff,kk)
        xnnu (j4,config%qn) = xnuc(ii,joff,kk,config%qn)

! -- neutrino source terms

        s    (j4) = qmy(ii,joff,kk)/densty(ii,joff,kk)

      end do

      grold(4) = 0.5_rk * ( gold(ii,qy_s-1,kk) + gold(ii-1,qy_s-1,kk) )
      grav (4) = 0.5_rk * ( gpot(ii,qy_s-1,kk) + gpot(ii-1,qy_s-1,kk) )


      do j = 5, nzn4

         dyzn      = yznr(j+qy_s-5) - yznl(j+qy_s-5)
! MPI grav, grold ist OK

         dphidx(j) = (grav (j) - grav (j-1)) / (dyzn * xzn(ii))
         dphold(j) = (grold(j) - grold(j-1)) / (dyzn * xzn(ii))
      end do



!-------  compute new velocity and energy and net-electron fraction

      dth = 0.5_rk * hydro%dt
      do j = 5, nzn4
! gravi
         utnu(j) = utnu(j) - dth * (dphidx(j) + dphold(j))
         qegr(j) = - (utnu(j)*dphidx(j) + ut(j)*dphold(j))
#ifdef ENECONS
         qegr(j) = config%laghyd*qegr(j)
#endif
         enu(j)  = enu(j) + dth * qegr(j)
!         enu (j) = enu (j) - dth * (utnu(j)*dphidx(j) + ut(j)*dphold(j))

! neutrino
         utnu(j)    = utnu(j)    + hydro%dt * s(j)          
         enu (j)    = enu (j)    + hydro%dt * s(j)*utnu(j)
!         xnnu(j,config%qn) = xnnu(j,config%qn) + dt * dydt(j)/rhonu(j) 
      end do


! --- update GR-Potential energy sourceterm for output
      do j = 1, nzn
         j4 = j + 4

         qgrv(ii,j+qy_s-1,kk) = qgrv(ii,j+qy_s-1,kk) + 0.5*qegr(j4)*rhonu(j4)

      end do

!-------  update velocity and energy and net-electron fraction

      do j = 1, nzn
         j4 = j + 4

         vely  (ii,j+qy_s-1,kk)    = utnu(j4)
         energy(ii,j+qy_s-1,kk)    = enu (j4)
         xnuc  (ii,j+qy_s-1,kk,config%qn) = xnnu(j4,config%qn)

      end do

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return
    end subroutine accely



#ifdef MIX
!=======================================================================
!>
!> \verbatim
!> simple "mixing alogorithm"
!> \endverbatim
!> \author  M. Rampp
!> \param jj
!> \param kk
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
    subroutine mix (jj, kk)

      use precision
      use abort

      use intgrs_hy
      use hydro_hy
      use gfloat_hy
      use mesh_hy

      use vnew_hy
      use phycon
      use nucparam
      use physcs_hy
      use specfun
      use eos_sn2

      use hydro_areas_mod
      use configure

      use state
      implicit none
! LOCAL variables that are not in modules

      real(KIND=RK), parameter :: vmix=-1.e9_rk,time_0=0.225_rk
      real(KIND=RK), parameter :: epsnewt=2.e-14_rk,ds=1.e-7_rk
      integer(kind=ik), parameter :: maxnewt=7

      logical :: lconv, ler

      real(KIND=RK)  ccu(config%q), cce(config%q), ccn(config%q), ccp(config%q)
      real(KIND=RK)  dmy1(config%q), dmy2(q,2) !dummys

!-- find region where mixing should occur


!- index of velocity minimum
      ishk_v=isrmin_v(areas%nx,velx(1,jj,kk),1)
!- in the vicinity of shock find largest gradient of entropy
      ishk=ishk_v-1 ! default
      dsdrmn=0._rk
      do i=ishk_v-4,ishk_v+4
         dsdr=(stot(i+1,jj,kk)-stot(i-1,jj,kk))/(xzn(i+1)-xzn(i-1))
         if (dsdr .lt. dsdrmn) then
            ishk=i
            dsdrmn=dsdr
         endif

      enddo
      ishk=ishk-1

      rmin=max( xzn(ishk)+vmix*max(time-time_0,0._rk) , 5.e5_rk )
      do i = 1,ishk-1
         ilow=i
         if (xzn(i) .ge. rmin) exit
      enddo


!      smx=0.0
!      ismx=ishk-1 ! default
!      do i=ishk-1,ilow,-1
!         if (stot(i,jj,kk).ge.smx) then
!            smx=stot(i,jj,kk)
!            ismx=i
!         endif
!         if (stot(i,jj,kk).lt.1.2*stot(ishk-1,jj,kk)) exit
!      enddo
!      imin=max(ilow,ismx)

      ismx=isrmax_v(ishk-ilow+1,stot(ilow,jj,kk),1)+ilow-1
      imin=max(ilow,ismx)
      if (stot(imin-1,jj,kk).ge.0.95_rk*stot(imin,jj,kk)) imin=imin-1


      imax=ishk-1



!      do i=ishk_v-4,ishk_v+1
!         write(*,'(i4,2e10.2)') i,stot(i,jj,kk),velx(i,jj,kk)
!      enddo
!      write(*,*)

      if (imax-imin .lt. 2) then
!         write(*,*) 'no mix'
         return
      else
         write(*,*) "Task ",myproc
         write(*,'(a5,f12.7,4I4,2f12.5)') 'mix: ',time, &
               ilow,ishk,imin,imax,stot(ismx,jj,kk),stot(ishk,jj,kk)
      endif

!-- mixit

      do i = imin, imax
         rho (i) = densty(i,jj,kk)
         tmp (i) = temp  (i,jj,kk)
         ek  (i) = 0.5_rk * (velx(i,jj,kk)**2 + vely(i,jj,kk)**2 + &
                             velz(i,jj,kk)**2                  )
         ei  (i) = (energy(i,jj,kk) - ek(i)) * rho(i)
         ei  (i) = max (ei(i), config%smallp)
      enddo
      do n=1,config%qn
         do i = imin, imax
            xn(i,n) = xnuc(i,jj,kk,n)
         enddo
      enddo
      kt   = imax-imin+1

!--
      eitot=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))
      scnst=SUM(stot(imin:imax,jj,kk))/real((imax-imin+1),kind=rk)

      lconv=.false.
      do inewt=1,maxnewt


         s(imin:imax)=scnst*(1.0_rk-ds)
         call eos (rho(imin:imax), tmp(imin:imax), xn(imin:imax,:), &
                   dmy1(imin:imax),dmy2(imin:imax,:),ei(imin:imax), &
                   p(imin:imax), gamc(imin:imax), s(imin:imax),     &
                   ccu(imin:imax),cce(imin:imax),ccn(imin:imax),    &
                  ccp(imin:imax),mode=4, nsemode=0,ler=ler)
         if (ler) raise_abort("mix():  eos failed (1)")
!         if (ler) call stopit('mix> eos failed (1)',0)
         dfunc=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))-eitot


         s(imin:imax)=scnst
         call eos (rho(imin:imax), tmp(imin:imax), xn(imin:imax,:), &
                   dmy1(imin:imax),dmy2(imin:imax,:),ei(imin:imax), &
                   p(imin:imax), gamc(imin:imax), s(imin:imax),     &
                   ccu(imin:imax),cce(imin:imax),ccn(imin:imax),    &
                   ccp(imin:imax),mode=4, nsemode=0,ler=ler)
         if (ler) raise_abort("mix():  eos failed (2)")
!         if (ler) call stopit('mix> eos failed (2)',0)
         func=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))-eitot

         dfunc=(func-dfunc)/(scnst-scnst*(1._rk-ds))
         scnst=scnst-func/dfunc

         if (abs(func/eitot).le.epsnewt .or. abs(func/dfunc/scnst).le.epsnewt) then
            lconv=.true.
            exit
         endif


      enddo

      if (.not.lconv) then
         write(*,*) "Task ",myproc
         write(*,*) 'mix PANIC:',inewt,func/eitot,func/dfunc,func,scnst
         raise_abort("mix(): no convergence")
!         call stopit('mix> no convergence',0)
      endif

! --  copy back

      do i = imin, imax
         energy(i,jj,kk) = ei(i) / rho(i) + ek(i)
         temp  (i,jj,kk) = tmp (i)
         stot  (i,jj,kk) = s   (i)
      enddo

!      call stopit('mixx: test',0)

      return      
    end subroutine mix
#endif
#ifdef NSMIX
#define RETARD
!=======================================================================
!=======================================================================
!>
!> \verbatim
!> simple neutron star "mixing alogorithm"
!> \endverbatim
!>
!> \author  R. Buras
!> \param jj
!> \param kk
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
    subroutine nsmix (jj, kk)

      use precision

      use intgrs_hy
      use hydro_hy
      use gfloat_hy
      use mesh_hy
      use nutrio_hy

      use vnew_hy
      use phycon
      use nucparam
      use physcs_hy
      use specfun
      use eos_sn2
      use revsho_hy

      use hydro_areas_mod
      use configure

      use state
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik), intent(in) :: jj,kk

#ifdef RETARD
      real(KIND=RK), parameter :: alpha=1._rk,tret=20.e-3_rk
#else
      real(KIND=RK), parameter :: alpha=2._rk
#endif

      real(KIND=RK), parameter :: epsnewt=2.e-14_rk,epsnewt2=2.e-13_rk, &
                                  ds=1.e-7_rk
      integer(kind=ik), parameter :: maxnewt=14,maxnewt2=100

      logical ::lconv, lconj, ler

      real(KIND=RK), dimension(config%q) :: ccu, cce, ccn, ccp, dmy1, &
                                     ylepscr, ynu, ylep, eint, nusto
      real(KIND=RK) scr, dsto, yleptot, yetot, eitot, scnst, dfunc, func, &
                    ylepton, dyedei
      real(KIND=RK) dmy2(config%q,2)
      integer(kind=ik) ismin,ismax,imin,imax,kt
      integer(kind=ik) i,j,n,inewt,jnewt
#ifdef RETARD
      integer(kind=ik) iti,itera
      real(KIND=RK), save :: tstart
      real(KIND=RK) ylcnst
      logical, save :: firstmix = .TRUE.
#endif

#ifdef RETARD
      if (firstmix) then
         firstmix = .FALSE.
         if (igodu.eq.1) tstart = s0_w
      endif
#endif

!-- find region where mixing should occur

      do i = 1,areas%nx
         ismin=i
         if (stot(i,jj,kk) .gt. 3.) exit
      enddo
      do i = 1,areas%nx
         ismax=i
         if (stot(i,jj,kk) .gt. 7.) exit
         if (densty(i,jj,kk) .lt. 1.e12_rk) exit
      enddo

      nusto(1:areas%nx)=stot(1:areas%nx,jj,kk) +  &
           1._rk/6._rk * ( temp(1:areas%nx,jj,kk) * pc_kmev * 1.e-5_rk /  &
           wc_hc)**3 * &
           (3._rk * 7._rk*pc_pi**2/15._rk +  &
            (cpotot(1:areas%nx,jj,kk,1)/temp(1:areas%nx,jj,kk) / pc_kmev )**2) / &
           densty(1:areas%nx,jj,kk) * (pc_mb * 1.e15_rk)

      imin = ismin
      do i = ismin,ismax-1
         dsto = nusto(i+1) - nusto(i)
         imin = i+1
         if (dsto .lt. 0._rk) exit
      enddo
      imax = imin -1
      scr = 7.
      do i=imin,ismax-1
         if  ( nusto(i) .lt. scr) then
            scr = nusto(i)
            imax = i
         endif
      enddo
      kt = imax - imin + 1

!      write (*,*) imax,kt,imin,ismax,ismin
!      write (*,*) nusto(90:100)
!      write (*,*) temp(90:100,jj,kk)

!QUICK FIX
      if (imax .lt. 10) return
      if (kt .lt. 2) then
         return
      else
         write(*,*) "Task ",myproc
         write(*,'(a5,f12.7,2I4,2f12.5)') 'mix: ',time, imin,imax,nusto(imin),nusto(imax)
      endif

#ifdef RETARD
      if (igodu.eq.0) then
         igodu=1
         tstart=time
         s0_w=tstart
         call printit_taskX(0,"nsmix: first encounter of negative gradient:",time)
      endif
      if (time-tret-tstart.lt.0) return

      if (igodu.eq.1) then
         igodu=2
         itera=500
         call printit_taskX(0,"nsmix: starting massive reshuffling:",time)
      else
         itera=1
      endif

!-- mixit

      do iti=1,itera
      if (igodu.eq.2) then      
         do i = 1,areas%nx
            ismin=i
            if (stot(i,jj,kk) .gt. 3._rk) exit
         enddo
         do i = 1,areas%nx
            ismax=i
            if (stot(i,jj,kk) .gt. 7.) exit
            if (densty(i,jj,kk) .lt. 1.e12_rk) exit
         enddo

         nusto(1:areas%nx)=stot(1:areas%nx,jj,kk) + &
              1._rk/6._rk * ( temp(1:areas%nx,jj,kk) * pc_kmev * 1.e-5_rk /  &
               wc_hc)**3 * &
              (3._rk * 7._rk*pc_pi**2/15._rk +  &
              (cpotot(1:areas%nx,jj,kk,1)/temp(1:areas%nx,jj,kk) / pc_kmev )**2)/ &
              densty(1:areas%nx,jj,kk) * (pc_mb * 1.e15_rk)

         imin = ismin
         do i = ismin,ismax-1
            dsto = nusto(i+1) - nusto(i)
            imin = i+1
            if (dsto .lt. 0._rk) exit
         enddo
         imax = imin -1
         scr = 7._rk
         do i=imin,ismax-1
            if  ( nusto(i) .lt. scr) then
               scr = nusto(i)
               imax = i
            endif
         enddo
         kt = imax - imin + 1
         write(*,*) "Task ",myproc
         write(*,'(a5,f12.7,2I4,2f12.5)') 'mix: ',time, imin,imax,nusto(imin),nusto(imax)
      endif
#endif

      do i = imin,imax
         rho (i) = densty(i,jj,kk)
         tmp (i) = temp  (i,jj,kk)
         ek  (i) = 0.5_rk * (velx(i,jj,kk)**2 + vely(i,jj,kk)**2 + &
                             velz(i,jj,kk)**2                )
         ei  (i) = (energy(i,jj,kk) - ek(i)) * rho(i)
         ei  (i) = max (ei(i), config%smallp)
         ynu (i) = ( dnutot(i,jj,kk,1)-dnutot(i,jj,kk,2) ) * pc_mb / rho(i)
!         write (*,*) i,ynu(i)
         ylep(i) = xnuc(i,jj,kk,config%qn) + ynu(i)

         if (config%restmass_version .ne. 2) then
            raise_abort('nsmix and restmass version .ne. 2 does not work')
         endif

         eint(i) = ei(i)
      enddo
      do n=1,config%qn
         do i = imin, imax
            xn(i,n) = xnuc(i,jj,kk,n)
         enddo
      enddo
      kt   = imax-imin+1

!--
      eitot=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))
      yetot=SUM(xn(imin:imax,config%qn)*rho(imin:imax)*vlfrac*dvx(imin:imax))
      yleptot=SUM(ylep(imin:imax)*rho(imin:imax)*vlfrac*dvx(imin:imax))
      scnst=SUM(nusto(imin:imax))/real((imax-imin+1),kind=rk)
      ylcnst=yleptot/SUM(rho(imin:imax)*vlfrac*dvx(imin:imax))

#ifdef RETARD
      if (igodu.eq.2) then
         xn(imin:imax,config%qn) = ylcnst - ynu(imin:imax)
      endif
#endif

      lconj=.false.
      do jnewt=1,maxnewt2

      lconv=.false.
      do inewt=1,maxnewt


         s(imin:imax)=scnst*(1._rk-ds)
         call eos (rho(imin:imax), tmp(imin:imax), xn(imin:imax,:), &
                   dmy1(imin:imax),dmy2(imin:imax,:),ei(imin:imax), &
                   p(imin:imax), gamc(imin:imax), s(imin:imax),     &
                   ccu(imin:imax),cce(imin:imax),ccn(imin:imax),    &
                   ccp(imin:imax),mode=4, nsemode=0,ler=ler)
         if (ler) raise_abort("mix(): eos failed (1)")
!         if (ler) call stopit('mix> eos failed (1)',0)
         dfunc=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))-eitot


         s(imin:imax)=scnst
         call eos (rho(imin:imax), tmp(imin:imax), xn(imin:imax,:), &
                   dmy1(imin:imax),dmy2(imin:imax,:),ei(imin:imax), &
                   p(imin:imax), gamc(imin:imax), s(imin:imax),     &
                   ccu(imin:imax),cce(imin:imax),ccn(imin:imax),    &
                   ccp(imin:imax),mode=4, nsemode=0,ler=ler)
         if (ler) raise_abort("mix(): eos failed (2)")
!         if (ler) call stopit('mix> eos failed (2)',0)
         func=SUM(ei(imin:imax)*vlfrac*dvx(imin:imax))-eitot

         dfunc=(func-dfunc)/(scnst-scnst*(1._rk-ds))
         scnst=scnst-func/dfunc

         if (abs(func/eitot).le.epsnewt .or. abs(func/dfunc/scnst).le.epsnewt) then
            lconv=.true.
            exit
         endif

      enddo

      if (.not.lconv) then
         write(*,*) "Task ",myproc
         write(*,*) 'mix PANIC:',inewt,func/eitot,func/dfunc,func,scnst
         raise_abort("mix(): no convergence")
!         call stopit('mix> no convergence',0)
      endif

#ifdef RETARD
      dyedei = 0.
      if (igodu.eq.2) goto 120
#endif
      ylepscr(imin:imax) = ylep(imin:imax) * (1._rk + alpha* &
                (ei(imin:imax) - eint(imin:imax))/eint(imin:imax))
      ylepton = SUM(ylepscr(imin:imax)*rho(imin:imax)* &
                                     vlfrac*dvx(imin:imax))
      ylepscr(imin:imax) = ylepscr(imin:imax) * yleptot / ylepton

      dyedei = 0._rk
      do i=imin,imax
         dyedei = max(dyedei,abs(ylepscr(i)- ynu(i) - &
                                 xn(i,config%qn))/xn(i,config%qn))
      enddo

      xn(imin:imax,config%qn) = ylepscr(imin:imax) - ynu(imin:imax)

#ifdef RETARD
 120  continue
#endif
!         write (*,*) 'mi2',epsnewt2,dyedei

      if (dyedei.le.epsnewt2) then
         lconj=.true.
         exit
      endif

      enddo

      if (.not.lconj) then
         write(*,*) "Task ",myproc
         write(*,*) 'mix2 PANIC:',jnewt,dyedei
         raise_abort("mix(): no convergence")
!         call stopit('mix> no convergence',0)
      endif


! --  copy back

      do i = imin, imax
         energy(i,jj,kk) = ei(i) / rho(i) + ek(i)
         temp  (i,jj,kk) = tmp (i)
         xnuc(i,jj,kk,config%qn) = xn(i,config%qn)
         call eos (rho(imin:imax), tmp(imin:imax), xn(imin:imax,:), &
              dmy1(imin:imax),dmy2(imin:imax,:),ei(imin:imax),      &
              p(imin:imax), gamc(imin:imax), s(imin:imax),          &
              ccu(imin:imax),cce(imin:imax),ccn(imin:imax),         &
              ccp(imin:imax),mode=1, nsemode=0,ler=ler)
         stot  (i,jj,kk) = s   (i)
      enddo

#ifdef RETARD
      enddo !itera
      if (igodu.eq.2) igodu=3
#endif

!      call stopit('mixx: test',0)


    end subroutine nsmix
#endif /* NSMIX */

end module sweeps_mod

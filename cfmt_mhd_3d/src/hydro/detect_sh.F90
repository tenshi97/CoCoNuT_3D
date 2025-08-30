module detect_sh_mod

implicit none

contains
!>
!> \verbatim
!> 
!> pseudo multi-dimensional shock detection (algorithm operates in sweeps)
!>
!>  Author: W. Keil
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine detect_sh
!
!
!     SUBROUTINE DETECT_SH
!     --------------------
!
!     
!     
!
  use precision

  use fpconst_hy
  use intgrs_hy
      
  use gfloat_hy
!  use mesh_hy ! forcheck
  use hydro_hy
  use vnew_hy
  use vold_hy
!  use grd_hy ! forcheck
  use intrp_hy

  use abort

  use mo_mpi
  use hlpare_hy
  use mpi_comm_routines
#ifndef DEBUG_TIMINGS
  use cputim
#endif
  use sweeps_mod

  use hydro_areas_mod
  use configure
  implicit none
! LOCAL variabls that are not in modules

  integer(kind=ik)        :: ibnd,k,i,j,j4,i4
  integer(kind=ik)        :: k4

  real(kind=rk)           :: du1,ptest,utest,dutest,dptest
  real(kind=rk),parameter :: epsiln  = 1.0_rk

  real(kind=rk)           :: dp1(config%q)
  integer(kind=ik) :: dest, src, mpistat(MPI_STATUS_SIZE), sendcount
  real(kind=rk) :: sbuf(0:config%qx+1, qy_s-1:qy_e+1), rbuf(0:config%qx+1, qy_s-1:qy_e+1)
  integer(kind=ik) :: ierr, j_end, k_end

  real(kind=rk) :: tim1(2), tim2(2)

!
!
!

     if (use_mpi) then
        if (ioy .eq. 1 .and. ioz .eq. 1) then
           j_end=qy_e
           k_end=qz_e
        else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
           j_end=qy_s
           k_end=qz_s
        end if
     else
        k_end = areas%nz
        j_end = areas%ny
     endif

     do k = qz_s-1, k_end+1
        do j = qy_s-1, j_end+1
           do i = 0, areas%nx+1
              ishck(i,j,k) = 0
           end do
        end do
     end do

     vxold(1:areas%nx,qy_s:j_end,qz_s:k_end) = velx(1:areas%nx,qy_s:j_end,qz_s:k_end)
     vyold(1:areas%nx,qy_s:j_end,qz_s:k_end) = vely(1:areas%nx,qy_s:j_end,qz_s:k_end)
     vzold(1:areas%nx,qy_s:j_end,qz_s:k_end) = velz(1:areas%nx,qy_s:j_end,qz_s:k_end)

     xyzswp = 1
     bndmin = config%bndmnx
     bndmax = config%bndmxx
     bndbot = config%bndmny
     bndtop = config%bndmxy
     bndlft = config%bndmnz
     bndrgt = config%bndmxz

     nzn  = areas%nx
     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8


     do k = qz_s,k_end
           do j=qy_s,j_end

!           --------------------
!           obtain boundary data
              call getrwx (j, k, 0)
!           -------------
!           detect shocks

              do i = 2,nzn7
                 dp1(i) = p(i+1) - p(i-1)
                 du1    = u(i+1) - u(i-1)
                 ptest  = epsiln * min ( p(i+1),p(i-1) ) - abs( dp1(i) )
                 utest  = config%smallu - abs(du1)
               
                 if ( utest .lt. 0.0_rk )  then
                    dutest = du1
                 else
                    dutest = 0.0_rk
                    dptest = 0.0_rk
                 end if
               
                 if ( ptest .lt. 0.0_rk )  then
                    dptest = 1.0_rk
                 else
                    dptest = 0.0_rk
                 end if
               
                 if ( du1   .ge. 0.0_rk ) dptest = 0.0_rk
!                 if ( dutest.eq.0.0_rk ) dptest = 0.0_rk
               
                 shockd(i) = -sign( 1.0_rk,dp1(i) )*dptest
                 shockd(i) = abs( shockd(i) )
              end do

              do i = 1, nzn
                 i4 = i + 4
                 ishck (i,j,k) = max(ishck(i,j,k),nint(shockd(i4)))
              end do
           end do
        end do

!     ----------------
!     second dimension

        if ( config%nsdim .ge. 2 ) then 

           if ((use_mpi) .and. areas%ny .le. 1) then

           ! do nothing in this case

           else

              xyzswp = 2
              bndmin = config%bndmny
              bndmax = config%bndmxy
              bndbot = config%bndmnx
              bndtop = config%bndmxx
              bndlft = config%bndmnz
              bndrgt = config%bndmxz

              nzn  = areas%ny

              if (use_mpi) then

                 ! default value, might be overriden in the next lines
                 bndmin = 6          ! or other value, might change later on

                 if (qy_s .eq. 1) then 
                    bndmin = config%bndmny
                 else if (qy_s .eq. 2) then
                    if (use_4neighbour_comm) bndmin = 9 ! three zones reflective, one communication
                 else if (qy_s .eq. 3 ) then
                    if (use_2neighbour_comm) bndmin = 8 ! two zones reflective, two communication
                    if (use_4neighbour_comm) bndmin = 8 ! two zones reflective, two communication
                 else if (qy_s .eq. 4) then
                    if (use_4neighbour_comm) bndmin = 10 ! one zone reflective, three communication
                 endif

                ! default value, might be overriden in the next lines
                 bndmax = 6          ! or other value, might change later on

                 if (qy_e .eq. areas%ny) then
                    bndmax = config%bndmxy
                 else if (qy_e .eq. areas%ny-1) then
                    if (use_4neighbour_comm) bndmax = 9 ! three zones reflective, one communication
                 else if (qy_e .eq. areas%ny-2) then
                    if (use_2neighbour_comm) bndmax = 8 ! two zones reflective, two communication
                    if (use_4neighbour_comm) bndmax = 8 ! two zones reflective, two communication
                 else if (qy_e .eq. areas%ny-3) then
                    if (use_4neighbour_comm) bndmax = 10 ! one zone reflective, three communication
 
                 end if

                 nzn = qy_proc

              endif ! use_mpi

              nzn1 = nzn + 1
              nzn2 = nzn + 2
              nzn3 = nzn + 3
              nzn4 = nzn + 4
              nzn5 = nzn + 5
              nzn6 = nzn + 6
              nzn7 = nzn + 7
              nzn8 = nzn + 8

              if (use_mpi) then
#ifndef DEBUG_TIMINGS
                 call second_v(tim1)
#endif

                 
                 if (use_1neighbour_comm) call mpi_comm_lb_ub_1neighbour
                 if (use_2neighbour_comm) call mpi_comm_lb_ub_2neighbour
                 if (use_4neighbour_comm) call mpi_comm_lb_ub_4neighbour

#ifndef DEBUG_TIMINGS
                 call second_v(tim2)
                 timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif
              endif

              do k = qz_s, qz_e

                 do i = 1, areas%nx
               
                    call getrwy (i, k, 0)
            
                    do j = 2,nzn7
                       dp1(j) = p(j+1) - p(j-1)
                       du1    = u(j+1) - u(j-1)
                       ptest  = epsiln * min ( p(j+1),p(j-1) ) - abs( dp1(j) )
                       utest  = config%smallu - abs(du1)
               
                       if ( utest.lt.0.0_rk )  then
                          dutest = du1
                       else
                          dutest = 0.0_rk
                       end if
                  
                       if ( ptest.lt.0.0_rk )  then
                          dptest = 1.0_rk
                       else
                          dptest = 0.0_rk
                       end if
                  
                       if ( du1   .ge.0.0_rk ) dptest = 0.0_rk
                       if ( dutest.eq.0.0_rk ) dptest = 0.0_rk
                  
                       shockd(j) = -sign( 1.0_rk,dp1(j) )*dptest
                       shockd(j) = abs( shockd(j) )
                    end do

                    do j = 1, nzn
                       j4 = j + 4
                       ishck (i,j+qy_s-1,k) = &
                            max(ishck(i,j+qy_s-1,k),nint(shockd(j4)))
                    end do
                 end do
              end do
           endif ! use_mpi
        end if ! config%nsdim = 2


!     ---------------
!     third dimension
        if ( config%nsdim .eq. 3 ) then 

           if ((use_mpi) .and. areas%nz .le. 1) then
              ! do nothing in this case
           else

              xyzswp = 3
              bndmin = config%bndmnz
              bndmax = config%bndmxz
              bndbot = config%bndmnx
              bndtop = config%bndmxx
              bndlft = config%bndmny
              bndrgt = config%bndmxy

              nzn  = areas%nz
              
              if (use_mpi) then
                 bndmin = 7          ! or other value, might change later on
                 bndmax = 7          ! or other value, might change later on 

                 nzn = qz_proc
              endif

              nzn1 = nzn + 1
              nzn2 = nzn + 2
              nzn3 = nzn + 3
              nzn4 = nzn + 4
              nzn5 = nzn + 5
              nzn6 = nzn + 6
              nzn7 = nzn + 7
              nzn8 = nzn + 8

              if (use_mpi) then
#ifndef DEBUG_TIMINGS
                 call second_v(tim1)
#endif

                 if (use_1neighbour_comm) call mpi_comm_pb_kb_1neighbour
                 if (use_2neighbour_comm) call mpi_comm_pb_kb_2neighbour
                 if (use_4neighbour_comm) call mpi_comm_pb_kb_4neighbour
#ifndef DEBUG_TIMINGS        
                 call second_v(tim2)
              
                 timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif
              endif

                 do j = qy_s, qy_e
                    do i = 1, areas%nx
               
                       call getrwz (i, j, 0)
               
                       do k = 2,nzn7
                          dp1(k) = p(k+1) - p(k-1)
                          du1    = u(k+1) - u(k-1)
                          ptest  = epsiln * min ( p(k+1),p(k-1) )- abs( dp1(k) )
                          utest  = config%smallu - abs(du1)
                  
                          if ( utest.lt.0.0_rk )  then
                             dutest = du1
                          else
                             dutest = 0.0_rk
                          end if
                  
                          if ( ptest.lt.0.0_rk )  then
                             dptest = 1.0_rk
                          else
                             dptest = 0.0_rk
                          end if
                  
                          if ( du1   .ge.0.0_rk ) dptest = 0.0_rk
                          if ( dutest.eq.0.0_rk ) dptest = 0.0_rk
                          
                          shockd(k) = -sign( 1.0_rk,dp1(k) )*dptest
                          shockd(k) = abs( shockd(k) )
                       end do
               
                       do k = 1, nzn
                          k4 = k + 4
                          ishck (i,j,k+qz_s-1) = &
                              max(ishck(i,j,k+qz_s-1),nint(shockd(k4)))
                       end do

                    end do
                 end do
              endif ! use_mpi
           end if ! nsdim = 3

! set boundaries according to specified boundary conditions

! Inner boundary condition in x-direction
              ibnd=config%bndmnx

              select case(ibnd)
              case(1:3)
                 ishck (0,:,:)   =ishck (1,:,:)
              case(4)
                 ishck (0,:,:)   =ishck (areas%nx,:,:)

              case(5)
!        add nonstandard condition here
              case(6:7)
!        add nonstandard condition here

              case default
                 raise_abort("detect_sh(): nocase")
              end select

! Outer boundary condition in x-direction
              ibnd=config%bndmxx

              select case(ibnd)
              case(1:3)
                 ishck (areas%nx+1,:,:)=ishck (areas%nx,:,:)
              case(4)
                 ishck (areas%nx+1,:,:)=ishck (1,:,:)

              case(5)

              case(6:7)

!        add nonstandard condition here
              case default
                 raise_abort("detect_sh(): nocase")
              end select
! Inner boundary condition in y-direction
              ibnd=config%bndmny

              if (use_mpi) then
                 if (qy_s .eq. 1) then 
                    ibnd = config%bndmny
                 else 
                    ibnd = 6               ! or other value, might change later on
                 end if
              endif

              select case(ibnd)
              case(1:3)
                 ishck (:,0,:)=ishck (:,1,:)
              case(4)
 
              case(5)

                 if (.not.(use_mpi)) then
                    ishck (:,0,:)=ishck (:,areas%ny,:)
                 endif

              case(6)


!        add nonstandard condition here
              case default
                 raise_abort("detect_sh(): nocase")
              end select

! Outer boundary condition in y-direction
              ibnd=config%bndmxy

              if (use_mpi) then
                 if (qy_e .eq. areas%ny) then
                    ibnd = config%bndmxy
                 else
                    ibnd = 6             ! or other value, might change later on 
                 end if
              endif

              select case(ibnd)
              case(1:3)
                 ishck (:,areas%ny+1,:)=ishck (:,areas%ny,:)
              case(4)

              case(5)

                 if (.not.(use_mpi)) then
                    ishck (:,areas%ny+1,:)=ishck (:,1,:)
                 endif

              case(6)


!        add nonstandard condition here
              case default
                 raise_abort("detect_sh(): nocase")
              end select

! Inner boundary condition in z-direction
              ibnd=config%bndmnz

              select case(ibnd)

              case(1:3)
                 ishck (:,:,0)=ishck (:,:,1)

              case(4)

                 if (use_mpi) then
#ifndef DEBUG_TIMINGS
                    call second_v(tim1)
#endif                    
                    sbuf(:,:) = real(ishck(:,:,qz_e),kind=rk)
                    call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
                    sendcount = (config%qx+2)*(qy_proc+2)
                    call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, &
                         dest, tag_sweep1, &
                         rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                         MPI_COMM_WORLD, mpistat, ierr)
                    if (qz_s .eq. 1) ishck(:,:,qz_s-1) = int(rbuf(:,:),kind=ik)
#ifndef DEBUG_TIMINGS
                    call second_v(tim2)
                    
                    timer%hydro_comm = timer%hydro_comm +(tim2-tim1)
#endif
                    if (areas%nz .eq. 1) ishck (:,:,qz_s-1)=ishck (:,:,k_end) ! special treatment for first area
                 else

                    ishck (:,:,0)=ishck (:,:,areas%nz)
                 endif


              case(5)

              case(7)
               
                 if (use_mpi) then
#ifndef DEBUG_TIMINGS
                    call second_v(tim1)
#endif
                    
                    sbuf(:,:) = real(ishck(:,:,qz_e),kind=rk)
                    call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
                    sendcount = (config%qx+2)*(qy_proc+2)
                    call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, &
                         dest, tag_sweep1, &
                         rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                         MPI_COMM_WORLD, mpistat, ierr)
                    if (qz_s .eq. 1) ishck(:,:,qz_s-1) = int(rbuf(:,:),kind=ik)
#ifndef DEBUG_TIMINGS
                    call second_v(tim2)
                    
                    timer%hydro_comm = timer%hydro_comm +(tim2-tim1)
#endif
                    if (areas%nz .eq. 1) ishck (:,:,qz_s-1)=ishck (:,:,k_end) ! special treatment for first area
                 endif



!        add nonstandard condition here
              case default
                 raise_abort("detect_sh(): nocase")
              end select

! Outer boundary condition in z-direction
              ibnd=config%bndmxz

              select case(ibnd)
              case(1:3)
                 ishck (:,:,areas%nz+1)=ishck (:,:,areas%nz)

              case(4)

                 if (use_mpi) then
#ifndef DEBUG_TIMINGS
                    call second_v(tim1)
#endif
                    sbuf(:,:) = real(ishck(:,:,qz_s),kind=rk)
                    call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
                    sendcount = (config%qx+2)*(qy_proc+2)
                    call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, &
                         dest, tag_sweep1, &
                      rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                      MPI_COMM_WORLD, mpistat, ierr)
                    if (qz_e .eq. config%qz) ishck(:,:,qz_e+1) = int(rbuf(:,:),kind=ik)
#ifndef DEBUG_TIMINGS
                    call second_v(tim2)

                    timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif
                    if (areas%nz .eq. 1) ishck (:,:,k_end+1)=ishck (:,:,qz_s) ! special treatment for first area
                 else
                    ishck (:,:,areas%nz+1)=ishck (:,:,1)
                 endif

                 case (5)


                 case(7)

                if (use_mpi) then
#ifndef DEBUG_TIMINGS
                    call second_v(tim1)
#endif
                    sbuf(:,:) = real(ishck(:,:,qz_s),kind=rk)
                    call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
                    sendcount = (config%qx+2)*(qy_proc+2)
                    call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, &
                         dest, tag_sweep1, &
                      rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                      MPI_COMM_WORLD, mpistat, ierr)
                    if (qz_e .eq. config%qz) ishck(:,:,qz_e+1) = int(rbuf(:,:),kind=ik)
#ifndef DEBUG_TIMINGS
                    call second_v(tim2)
                    timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif
                    if (areas%nz .eq. 1) ishck (:,:,k_end+1)=ishck (:,:,qz_s) ! special treatment for first area
                 endif


!        add nonstandard condition here
              case default
                 raise_abort("detect_sh(): nocase")
              end select
!
! return
              !
              return
!
! end of DETECT_SH
!
end subroutine detect_sh

end module detect_sh_mod

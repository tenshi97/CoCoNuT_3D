module tstep_mod

implicit none

contains

!>
!> \verbatim
!> Compute the size of the hydro time step
!>
!> \endverbatim
!>
!>  \author M. Rampp
!>  \todo Check 3D case
!>  \param dtnew  size of new timestep
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date$
!>   
!> \endverbatim
!>
subroutine tstep (dtnew, selftime, childrentime)
  use precision
  use error

  use intgrs_hy

  use vnew_hy
!  use vold_hy ! forcheck
  use mesh_hy
  use gfloat_hy
  
  use specfun

  use mo_mpi

  use hydro_areas_mod
  use configure
  use state
  use print_stdout_mod
  use cputim
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2)
  integer(kind=ik) :: i, ic, jc, kc, j, k, im
  real(kind=rk)    :: dtc, olddt, dx, dyi, cs2, dxi, dzi, cspd, vabs
  real(kind=rk)    :: dtnew

  character*6      :: label

  real(kind=rk)    :: scrch(config%max_dim), scratch(config%max_dim),     &
                      ceul(config%max_dim), scratch2(config%max_dim)   
                      !, iscrch(q)    
 

  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  dtc   = 1.e35_rk
  olddt = hydro%dt
!
!
!-------  check CFL-condition
!
! attention: the nsdim=1 case will be called in multi-D runs 
! for the first 1 area (this is important for the MPI-mode)!!!
  if (config%nsdim .eq. 1 .or. areas%ny .lt. 4)  then
!
     ceul(:)    = 0._rk
     scratch(:) = 0._rk
     do i = 1, areas%nx
        ceul(i) = (gammac(i,qy_s,qz_s) * press(i,qy_s,qz_s) &
                  / densty(i,qy_s,qz_s))
     enddo

     call fastsqrt(ceul, areas%nx)

     do  i = 1, areas%nx
        dx       = xznr(i) - xznl(i)
        scrch(i) = dx / (abs (velx(i,qy_s,qz_s) ) + ceul(i))

     enddo
!
     im = ISAMIN_V (areas%nx, scrch, 1_ik)
     if (dtc .ge. scrch(im))  then
        ic  = im
        dtc = scrch(im)
     endif
!

     jc = qy_s
     kc = qz_s

!

  endif ! nsdim .eq. 1 .or. ny .lt. 4
!
!
!      if (nsdim .eq. 2)  then
!     temporaer:
 if (config%nsdim .eq. 2 .and. areas%ny .ge. 4)  then
!
    if (config%igeomy .le. 2)   then
!

       do j = qy_s, qy_e

          dyi = 1.0_rk  /  ( yznr(j) - yznl(j) )
          do i = 1, areas%nx
             cs2      = gammac(i,j,1) * press(i,j,1) / densty(i,j,1)
             dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )

             scratch(i) = cs2 * (dxi**2 + dyi**2)

             scrch(i) = abs (velx(i,j,1)) * dxi  +     &
                        abs (vely(i,j,1)) * dyi   
          enddo

          call fastsqrt(scratch, areas%nx)

          do i = 1, areas%nx
             scrch(i) = scrch(i)  + scratch(i)
             scrch(i) = 1.0_rk / scrch(i)
          enddo
!
          im = ISAMIN_V (areas%nx, scrch, 1_ik)
          if (dtc .ge. scrch(im))  then
             ic  = im
             jc  = j
             dtc = scrch(im)
          end if
       enddo
       !
    else

       do j = qy_s, qy_e

          do i = 1, areas%nx
             cs2      = gammac(i,j,1) * press(i,j,1) / densty(i,j,1)
             dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )
             dyi      = 1.0_rk  /  ( xzn(i) * ( yznr(j) - yznl(j) ) )
             scratch(i) = cs2 * (dxi**2 + dyi**2)
             scrch(i) = abs (velx(i,j,1)) * dxi  +      &
                        abs (vely(i,j,1)) * dyi 
          enddo

          call fastsqrt(scratch, areas%nx)

          do i = 1, areas%nx
             scrch(i) = scrch(i) + scratch(i)
             scrch(i) = 1.0_rk / scrch(i)
            enddo
!
!DIR$ INLINE
            im = ISAMIN_V (areas%nx, scrch, 1_ik)
            if (dtc .ge. scrch(im))  then
               ic  = im
               jc  = j
               dtc = scrch(im)
            end if
         enddo
      end if ! config%igeomy .le. 2
!
      kc = 1
       !
   end if ! nsdim .eq. 2
!
!
!      if (nsdim .eq. 2)  then
!     temporaer:
 if (config%nsdim .eq. 3 .and. areas%ny .ge. 4 .and. areas%nz .ge. 4)  then
      if (config%igeomx .eq. 0)   then
!
         scratch(:) = 0._rk
         do k = qz_s, qz_e

            dzi = 1.0_rk  /  ( zznr(k) - zznl(k) )

            do j = qy_s, qy_e

               dyi = 1.0_rk  /  ( yznr(j) - yznl(j) )

               do i = 1, areas%nx
                  cs2      = gammac(i,j,k)*press(i,j,k) / densty(i,j,k)
                  dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )
                  scratch(i) =  cs2 * (dxi**2 + dyi**2 + dzi**2)
                  scrch(i) = abs (velx(i,j,k)) * dxi  + &
                       abs (vely(i,j,k)) * dyi  +       &
                       abs (velz(i,j,k)) * dzi 
               enddo

               call fastsqrt(scratch, areas%nx)
               
               do i=1, areas%nx
                  scrch(i) = scrch(i) + scratch(i)
                  scrch(i) = 1.0_rk / scrch(i)
               enddo

               im = ISAMIN_V (areas%nx, scrch, 1_ik)
               if (dtc .ge. scrch(im))  then
                  ic  = im
                  jc  = j
                  kc  = k
                  dtc = scrch(im)
               end if
            enddo
         enddo
!
      end if ! (config%igeomx .eq. 0) 
!
      if (config%igeomx .eq. 2)   then
!
!-------- WARNING: The coordinates are assumed to be ordered (r,th,phi)
!
         scratch(:) = 0._rk

         do j = qy_s, qy_e
            scratch(j) = yzn(j)
         enddo

         call fastsin(scratch, (qy_e-qy_s + 1))
         


         do k = qz_s, qz_e
            do j = qy_s, qy_e

               do i = 1, areas%nx
                  cs2      = gammac(i,j,k) * press(i,j,k) / densty(i,j,k)
                  dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )
                  dyi      = 1.0_rk  /  ( xzn(i) * ( yznr(j) - yznl(j) ) )
                  dzi      = 1.0_rk  /  ( xzn(i) * scratch(j) &
                                    *  ( zznr(k) - zznl(k) ) )
                  scratch2(i) =  cs2 * (dxi**2 + dyi**2 +  dzi**2)

                  scrch(i) = abs (velx(i,j,k)) * dxi  + &
                             abs (vely(i,j,k)) * dyi  + &
                             abs (velz(i,j,k)) * dzi 
               enddo

               call fastsqrt(scratch2, areas%nx)

               do i=1, areas%nx
                  scrch(i) = scrch(i) + scratch2(i)
                  scrch(i) = 1.0_rk / scrch(i)
               enddo
!
               im = ISAMIN_V (areas%nx, scrch, 1_ik)
               if (dtc .ge. scrch(im))  then
                  ic  = im
                  jc  = j
                  kc  = k
                  dtc = scrch(im)
               end if
            enddo
         enddo
!
      end if !(config%igeomx .eq. 2) 
!
   end if ! nsdim .eq. 3
!
!
   dtc = config%cfl * dtc
   hydro%dt  = min (dtc, config%dtmax)
!
   if (hydro%dt .eq. dtc) then
      label = 'cfl'
   else
      label = 'max'
   end if
!
   dtnew = hydro%dt

   if (hydro%dt .gt. 1.2_rk*olddt) then
      hydro%dt    = 1.2_rk * olddt
      label = 'old'
   end if

!   print *,hydro%dt, olddt, dtnew, config%cfl,dtc
!   stop 'a'
!
!      if ((mod (nstep, config%itstp) .ne. 0).and.(nstep.ne.1))  return
!
   if (config%nsdim .eq. 1)   then
      vabs  = abs (velx(ic,jc,kc) )
   else if (config%nsdim .eq. 2) then
! is this an error in case of rotation
      vabs  = sqrt( velx(ic,jc,kc)**2 + vely(ic,jc,kc)**2 + velz(ic,jc,kc)**2 )
   else if (config%nsdim .eq. 3) then
      vabs  = sqrt( velx(ic,jc,kc)**2 + vely(ic,jc,kc)**2 + velz(ic,jc,kc)**2 )
   end if
!
   if(config%itstp .ne. 0) then
      if(mod(nstep,config%itstp) .eq. 0) then
         cspd = sqrt(gammac(ic,jc,kc) * press(ic,jc,kc) / densty(ic,jc,kc))

         write (*,1001) myproc, nstep, label, time, hydro%dt, ic, jc, kc, cspd, vabs
 1001       format("Task ",i0,1x, i8, 1x, a3, 2x, 't = ', 1pe12.5,                 &
                 ' [s]  dt = ', 1pe9.3,                                 &
                 ' [s]', '    ic = ', i3, '  jc = ', i3, '  kc = ', i3, &
                 '  cspd = ', e8.2, '  vabs = ', e8.2)
!
      endif
   endif

   if (hydro%dt .lt. config%dtmin) then
      call printit_taskX(0,"W A R N I N G :   dt = ",hydro%dt)
      call printit_taskX(0,"               dtmin = ",config%dtmin)
      raise_error("tstep(): => tstep")

   end if
!
!

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
   return
end subroutine tstep

#ifdef OPTIMAL_13_SWITCH
subroutine optimal_13_switch

   use precision
   use vnew_hy
   use hlpare_hy
   use totare_hy
!   use arecon_hy
   use mesh_hy
   use intgrs_hy
   use specfun

   use mo_mpi

   use mapare_proc
   
   use hydro_areas_mod
   use configure
   implicit none
! LOCAL variables that are not in modules

   integer(kind=ik) :: i, j, k, im, sw123, ierr, sw123_rcv
   real(kind=rk)    :: dtc_1d, dtc_2d, dtc_3d, dx, dyi, cs2, dxi, dzi

   real(kind=rk)    :: scrch(config%q), scratch(config%q), ceul(config%q) 

!> calcultate at first the 1d cfl timestep

   ceul(:) = 0._rk

   do i = 1, config%qx
      ceul(i) = gammac(i,qy_s,qz_s) * press(i,qy_s,qz_s) &
                  / densty(i,qy_s,qz_s)
   enddo
   
   call fastsqrt( ceul, config%qx)


   do  i = 1, config%qx
        dx       = xznr(i) - xznl(i)
        scrch(i) = dx / (abs (velx(i,qy_s,qz_s) ) + ceul(i))
   enddo

   im = ISAMIN_V (config%qx, scrch, 1_ik)
   dtc_1d = scrch(im)

!> set the switch to the old value
   sw123 = areas%ix_are(1,2)

!> now find the zone where the 2d/3d cfl-timestep is larger than 1d-timestep, 
!> choose this zone to switch between 1d and 2d/3d
   if (config%nsdim .eq. 2) then

      scratch(:) = 0._rk

       do i = 1, config%qx
          !> calculate 2d-cfl-timestep for the i-th radial zone
          dtc_2d   = 1.e35_rk
          do j = qy_s, qy_e
             cs2      = gammac(i,j,1) * press(i,j,1) / densty(i,j,1)
             dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )
             dyi      = 1.0_rk  /  ( xzn(i) * ( yznr(j) - yznl(j) ) )

             scratch(j) = cs2 * (dxi**2 + dyi**2)
             scrch(j) = abs (velx(i,j,1)) * dxi  +      &
                  abs (vely(i,j,1)) * dyi 
          enddo

          call fastsqrt(scratch, (qy_e-qy_s+1) )

          do j = qy_s, qy_e
             scrch(j) = scrch(j) + scratch(j)
             scrch(j) = 1.0_rk / scrch(j)
          enddo

          im = ISAMIN_V (qy_proc, scrch, 1_ik)
          dtc_2d = scrch(im)

          if (dtc_2d .le. dtc_1d) then
             sw123 = i
          endif

       enddo

   endif

   if (config%nsdim .eq. 3) then
      
      scratch(:) = 0._rk
      scratch2(:)= 0._rk
      do j = qy_s, qy_e
         scratch(j) = yzn(j)
      enddo

      call fastsin(scratch, (qy_e-qy_s+1) )

      

       do i = 1, config%qx
          !> calculate 2d-cfl-timestep for the i-th radial zone
          dtc_3d   = 1.e35_rk
          do k = qz_s, qz_e
             do j = qy_s, qy_e
                cs2      = gammac(i,j,k) * press(i,j,k) / densty(i,j,k)
                dxi      = 1.0_rk  /  ( xznr(i) - xznl(i) )
                dyi      = 1.0_rk  /  ( xzn(i) * ( yznr(j) - yznl(j) ) )
                dzi      = 1.0_rk  /  ( xzn(i) * scratch(j) &
                                  *  ( zznr(k) - zznl(k) ) )

                scratch2(j) = cs2 * (dxi**2 + dyi**2 +  dzi**2)

                scrch(j) = abs (velx(i,j,k)) * dxi  + &
                           abs (vely(i,j,k)) * dyi  + &
                           abs (velz(i,j,k)) * dzi

             enddo
             
             call fastsqrt(scratch2, (qy_e-qy_s+1))

             scrch(j) = scrch(j) = scratch2(j)
             scrch(j) = 1.0_rk / scrch(j)
             enddo
             im = ISAMIN_V (qy_proc, scrch, 1_ik)
             dtc_2d = scrch(im)

!         write(*,*) 'dtc_2d, dtc_3d ',dtc_2d,dtc_3d
             if (dtc_3d .ge. dtc_2d)  then
!         write(*,*) 'dtc_2d, dtc_3d ',im,k,dtc_2d,dtc_3d
                dtc_3d = dtc_2d
             end if
          enddo

          if (dtc_3d .lt. dtc_1d) then
!             write(*,*) 'i,dtc_3d,dtc_1d ',i,dtc_3d,dtc_1d
             sw123 = i
          endif

       enddo

   endif

!> now set the new switch zone, if this zone changes more than 2 zones
!   write(*,*) 'sw123, dtc_1d, dtc_3d ',sw123,ix_are(1,2),dtc_1d,dtc_3d,xzn(sw123),densty(sw123,1,1)
   if (use_mpi) then
      sw123_rcv = 0
       call MPI_AllReduce(sw123, sw123_rcv, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

       sw123 = sw123_rcv

   endif
   write(*,*) "Task ",myproc
   write(*,*) 'sw123, dtc_1d, dtc_3d ',sw123,ix_are(1,2),dtc_1d,dtc_3d,xzn(sw123),densty(sw123,qy_s,qz_s)

   if (sw123 .gt. ix_are(1,2)+2 .or. sw123 .lt. ix_are(1,2)-2) then
       ix_are(1,2) = sw123
       ix_are(2,1) = sw123+1
       call mapare
   endif

end subroutine optimal_13_switch

#endif /* OPTIMAL_13_SWITCH */

end module tstep_mod

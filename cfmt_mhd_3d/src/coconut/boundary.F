c     ==================================================================

      MODULE gr_boundary

c     ==================================================================

      IMPLICIT NONE

      CONTAINS




c     ==================================================================

      SUBROUTINE hydro_boundary_conditions

c     ==================================================================
      
      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE boundary_cfc
      USE conserved_cfc
      USE perm_aux_cfc

      use configure
      IMPLICIT NONE
 
      integer (kind=ik) :: i_nuc

c Note: In 3D, we should actually set boundary conditions at r=0, 
c       theta=0, and theta=pi by copying the values from zones with 
c       (r',theta',varphi') = (r,-theta,varphi+pi) at the origin
c       (r',theta',varphi') = (r, theta,varphi+pi) at the axis

      call apply_boundary_conditions (rho, 1, 1, 1)
      call apply_boundary_conditions (v_1, symm_v_r (1),
     &     symm_v_theta_lower (1), symm_v_theta_upper (1))
      call apply_boundary_conditions (v_2, symm_v_r (2),
     &     symm_v_theta_lower (2), symm_v_theta_upper (2))
      call apply_boundary_conditions (v_3, symm_v_r (3),
     &     symm_v_theta_lower (3), symm_v_theta_upper (3))
      call apply_boundary_conditions (eps, 1, 1, 1)
      do i_nuc=1,config%qn
         call apply_boundary_conditions(xnnu(:,:,:,i_nuc),1,1,1)
      end do
      call apply_boundary_conditions (t, 1, 1, 1)
      call apply_boundary_conditions (gamm, 1, 1, 1)
      call apply_boundary_conditions (p, 1, 1, 1)
      call apply_boundary_conditions (h, 1, 1, 1)
      call apply_boundary_conditions (w, 1, 1, 1)
      call apply_boundary_conditions (c_sound_squared, 1, 1, 1)

#ifdef CFC_MHD
c      call apply_boundary_conditions (b_1, 1, 1, 1)
c      call apply_boundary_conditions (b_2, 1,
c     &     symm_v_theta_lower (3), symm_v_theta_upper (3))
c      call apply_boundary_conditions (b_3, -1, 1, 1)
      call apply_boundary_conditions (b_1, symm_v_r (1),
     &     symm_v_theta_lower (1), symm_v_theta_upper (1))
      call apply_boundary_conditions (b_2, symm_v_r (2),
     &     symm_v_theta_lower (2), symm_v_theta_upper (2))
c      call apply_boundary_conditions (b_3, symm_v_r (3),
c     &     symm_v_theta_lower (3), symm_v_theta_upper (3))
      call apply_boundary_conditions (b_3, -1, 1, -1)
      call apply_boundary_conditions (psi, 1, 1, 1)
#endif /* CFC_MHD */
      
      call hydrostatic_atmosphere
      
      return

      END SUBROUTINE hydro_boundary_conditions

c     ==================================================================
  


#ifdef CFC_MHD
c     ==================================================================

      SUBROUTINE induction_boundary_conditions

c     ==================================================================
      
      USE precision
      USE abort

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE fluxes_cfc
c      USE metric_cfc
c      USE boundary_cfc
c      USE conserved_cfc
      USE perm_aux_cfc

      USE configure

      IMPLICIT NONE
 
      integer (kind=ik) :: i, j, k, jk, js, je, ks, ke

ccc   radial boundary conditions (inner boundary condition should
c     not matter since the line integral over the E-field vanishes
c     anyway)

#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     E_r
            emf_face (0, j, k, 1, 2) = - emf_face (1, j, k, 1, 2) 
c     E_varphi
            emf_face (0, j, k, 2, 2) = - emf_face (1, j, k, 2, 2)
c     E_r
            emf_face (0, j, k, 1, 3) = - emf_face (1, j, k, 1, 3)
c     E_theta
            emf_face (0, j, k, 2, 3) = - emf_face (1, j, k, 2, 3)

c     constant extrapolation at outer boundary
c     E_r
            emf_face (m+1, j, k, 1, 2) = emf_face (m, j, k, 1, 2) 
c     E_varphi
            emf_face (m+1, j, k, 2, 2) = emf_face (m, j, k, 2, 2)
c     E_r
            emf_face (m+1, j, k, 1, 3) = emf_face (m, j, k, 1, 3)
c     E_theta
            emf_face (m+1, j, k, 2, 3) = emf_face (m, j, k, 2, 3)

            if (j .eq. n_s) then
               emf_face (0, j-1, k, 1, 2) = 
     &              - emf_face (1, j-1, k, 1, 2) 
               emf_face (0, j-1, k, 2, 2) =
     &              - emf_face (1, j-1, k, 2, 2)
               emf_face (m+1, j-1, k, 1, 2) =
     &              emf_face (m, j-1, k, 1, 2) 
               emf_face (m+1, j-1, k, 2, 2) = 
     &              emf_face (m, j-1, k, 2, 2)
            end if

            if (k .eq. o_s) then
               emf_face (0, j, k-1, 1, 3) = 
     &              - emf_face (1, j, k-1, 1, 3)
               emf_face (0, j, k-1, 2, 3) = 
     &              - emf_face (1, j, k-1, 2, 3)
               emf_face (m+1, j, k-1, 1, 3) = 
     &              emf_face (m, j, k-1, 1, 3)
               emf_face (m+1, j, k-1, 2, 3) = 
     &              emf_face (m, j, k-1, 2, 3)
            end if

      enddo

      if (config%nsdim .gt. 1) then
ccc   meridional boundary conditions

#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc !loop over all theta/varphi-zones:
                                  !this ensures that all OpenMP Threads
                                  !work on LOCAL data

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            if (j .eq. 1) then
c     E_theta, E_varphi on radial interfaces (not used or
c     weighted by zero anyway)
               do i = 0, m
                  emf_face (i, 0, k, :, 1) =
     &                 -emf_face (i, 1, k, :, 1)
               enddo
c     E_r and E_theta on varphi interfaces (only E_r
c     is used later on)
               do i = 0, m+1
                  emf_face (i, 0, k, 1, 3) =
     &                 -emf_face (i, 1, k, 1, 3)
c     need E_r -> 0 since v_theta B_phi -> 0 and v_phi -> B_theta -> 0
                  emf_face (i, 0, k, 2, 3) =
     &                 -emf_face (i, 1, k, 2, 3)
               end do

               if (k .eq. o_s) then
                  do i = 0, m+1
                     emf_face (i, 0, k-1, 1, 3) =
     &                    -emf_face (i, 1, k-1, 1, 3)
c     need E_r -> 0 since v_theta B_phi -> 0 and v_phi -> B_theta -> 0
                     emf_face (i, 0, k-1, 2, 3) =
     &                    -emf_face (i, 1, k-1, 2, 3)
                  end do
               end if
            endif
            
            if (j. eq. n) then
               do i = 0, m
                  emf_face (i, n+1, k, :, 1) =
     &                 -emf_face (i, n, k, :, 1)
               enddo
c     E_r and E_theta on varphi interfaces (only E_r
c     is used later on)
               do i = 0, m+1
                  emf_face (i, n+1, k, 1, 3) =
     &                 -emf_face (i, n, k, 1, 3)
                  emf_face (i, n+1, k, 2, 3) =
     &                 -emf_face (i, n, k, 2, 3)
               end do

               if (k .eq. o_s) then
                  do i = 0, m+1
                     emf_face (i, n+1, k-1, 1, 3) =
     &                    -emf_face (i, n, k-1, 1, 3)
                     emf_face (i, n+1, k-1, 2, 3) =
     &                    -emf_face (i, n, k-1, 2, 3)
                  end do
               end if

            endif
            
         enddo

         if (config%nsdim .eq. 3) then

c     -----------------------------------------------------------------
c     3D case

            if (config%bndmxz .eq. 1 .and. config%bndmnz .eq. 1) then
               raise_abort("Check MHD BC (wedge in varphi)")

ccc   boundary conditions in rotational direction
#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc

                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
                  
                  if (k .eq. 1) then
                     do i = 0, m+1
                        emf_face (i, j, 0, 1:2, 1:2) =
     &                       -emf_face (i, j, 1, 1:2, 1:2)
                     enddo
                     if (j .eq. 1) then
                        do i = 0, m+1
                           emf_face (i, 0, 0, 1:2, 1:2) =
     &                          -emf_face (i, 0, 1, 1:2, 1:2)
                        enddo
                     end if
                     if (j .eq. n) then
                        do i = 0, m+1
                           emf_face (i, n+1, 0, 1:2, 1:2) =
     &                          -emf_face (i, n+1, 1, 1:2, 1:2)
                        enddo
                     end if
                  endif
                  
                  if (k .eq. o) then
                     do i = 0, m+1
                        emf_face (i, j, o+1, 1:2, 1:2) =
     &                       -emf_face (i, j, o, 1:2, 1:2)
                     enddo
                     if (j .eq. 1) then
                        do i = 0, m+1
                           emf_face (i, 0, o+1, 1:2, 1:2) =
     &                          -emf_face (i, 0, o, 1:2, 1:2)
                        enddo
                     end if
                     if (j .eq. n) then
                        do i = 0, m+1
                           emf_face (i, n+1, o+1, 1:2, 1:2) =
     &                          -emf_face (i, n+1, o, 1:2, 1:2)
                        enddo
                     end if
                  endif
                  
               enddo
            else if (config%bndmxz .ne. 4 .or.config%bndmnz .ne. 4) then
               raise_abort("varphi-boundary cond. not allowed")
            end if
            
         endif ! config%nsdim = 3
      endif ! config%nsdim .gt. 1

      
      return

      END SUBROUTINE induction_boundary_conditions

c     ================================================================== 
#endif /* CFC_MHD */


c     ==================================================================

      SUBROUTINE apply_boundary_conditions
     &     (q, symm_r, symm_theta_lower, symm_theta_upper)

c     ==================================================================
      
      USE precision
      USE size_cfc
      USE abort

      use configure
      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk
      
      integer (kind=ik) :: symm_r, symm_theta_upper, symm_theta_lower
      
      real (kind=rk) :: q (-3:m+4,n_s-4:n_e+4,o_s-4:o_e+4)



ccc   radial boundary conditions

#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 0, - 3, - 1
               
            q (i, j, k) = q (1 - i, j, k) * symm_r
            
         enddo

      enddo
      
      call extrapolate_out (q, 4, 0)
       
      if (config%nsdim .gt. 1) then
ccc   meridional boundary conditions

#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc !loop over all theta/varphi-zones:
                                  !this ensures that all OpenMP Threads
                                  !work on LOCAL data

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            if (j .le. 4) then
               do i = 1, m + 4
                  q (i, 1 - j, k) =
     &                 q (i, j, k) * symm_theta_upper
               enddo
            endif
            
            if (j. gt. n - 4) then
               do i = 1, m + 4
                  q (i, 2 * n + 1 - j, k) =
     &                 q (i, j, k) * symm_theta_lower
               enddo
            endif
            
         enddo

#if 1 /* does not work yet */         
         if (config%nsdim .eq. 3) then

c     -----------------------------------------------------------------
c     3D case

            if (config%bndmxz .eq. 1 .and. config%bndmnz .eq. 1) then

ccc   boundary conditions in rotational direction
#if defined(OPENMP_CFC) && defined(OPEN_MP_2D)
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc

                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
                  
                  if (k .le. 4) then
                     do i = - 3, m + 4
                        q (i, j, 1 - k) = q (i, j, 1)
                     enddo
                  endif
                  
                  if (k .ge. n - 3) then
                     do i = - 3, m + 4
                        q (i, j, 2 * o + 1 - k) = q (i, j, k)
                     enddo
                  endif
                  
               enddo
            else if (config%bndmxz .ne. 4 .or.config%bndmnz .ne. 4) then
               raise_abort("varphi-boundary cond. not allowed")
            end if
            
         endif ! config%nsdim = 3
#endif /* does not work yet */
      endif ! config%nsdim .gt. 1
      return

      END SUBROUTINE apply_boundary_conditions

c     ==================================================================




c     ==================================================================

      SUBROUTINE extrapolate_out (q, range, order)

c     ==================================================================
      
      USE precision

      USE size_cfc
      USE grid_cfc, ONLY: delta_r_if

      IMPLICIT NONE

      real (kind=rk), intent (inout) ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      integer (kind=ik), intent (in) :: range, order

      integer (kind=ik) :: i, j, k, jk
      
      if (order .eq. 0) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = m + 1, m + range
               q (i, j, k) = q (i - 1, j, k)
            enddo

         enddo

      elseif (order .eq. 1) then


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = m + 1, m + range
               q (i, j, k) = 
     &              (delta_r_if (i - 2) + delta_r_if (i - 1)) /
     &              delta_r_if (i - 2) * 
     &              q (i - 1, j, k) +
     &              
     &              delta_r_if (i - 1) / (- delta_r_if (i - 2)) *
     &              q (i - 2, j, k)
            enddo

         enddo

      elseif (order .eq. 2) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = m + 1, m + range

               q (i, j, k) = 
     &              (delta_r_if (i - 3) + delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) * 
     &              (delta_r_if (i - 2) + delta_r_if (i - 1)) /
     &              ((delta_r_if (i - 3) + delta_r_if (i - 2)) *
     &              delta_r_if (i - 2)) *
     &              q (i - 1, j, k) +
     &              
     &              (delta_r_if (i - 3) + delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) * 
     &              delta_r_if (i - 1) /
     &              (delta_r_if (i - 3) * (- delta_r_if (i - 2))) * 
     &              q (i - 2, j, k) +
     &                 
     &              (delta_r_if (i - 2) + delta_r_if (i - 1)) *
     &              delta_r_if (i - 1) / 
     &              ((- delta_r_if (i - 3)) *
     &              (- delta_r_if (i - 3) - delta_r_if (i - 2))) *
     &              q (i - 3, j, k)
            enddo

         enddo

      elseif (order .eq. 3) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = m + 1, m + range
               q (i, j, k) =
     &              (delta_r_if (i - 4) + delta_r_if (i - 3) +
     &              delta_r_if (i - 2) + 
     &              delta_r_if (i - 1)) * (delta_r_if (i - 3) +
     &              delta_r_if (i - 2) + 
     &              delta_r_if (i - 1)) * (delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) / 
     &              ((delta_r_if (i - 4) + delta_r_if (i - 3) +
     &              delta_r_if (i - 2)) * 
     &              (delta_r_if (i - 3) + delta_r_if (i - 2)) *
     &              delta_r_if (i - 2)) * 
     &              q (i - 1, j, k) +
     &                 
     &              (delta_r_if (i - 4) + delta_r_if (i - 3) +
     &              delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) * (delta_r_if (i - 3) +
     &              delta_r_if (i - 2) + 
     &              delta_r_if (i - 1)) * delta_r_if (i - 1) /
     &              ((delta_r_if (i - 4) + 
     &              delta_r_if (i - 3)) * delta_r_if (i - 3) *
     &              (- delta_r_if (i - 2))) *
     &              q (i - 2, j, k) +
     &              
     &              (delta_r_if (i - 4) + delta_r_if (i - 3) +
     &              delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) * (delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) *
     &              delta_r_if (i - 1) / (delta_r_if (i - 4) *
     &              (- delta_r_if (i - 3)) *
     &              (- delta_r_if (i - 3) - delta_r_if (i - 2))) * 
     &              q (i - 3, j, k) +
     &                 
     &              (delta_r_if (i - 3) + delta_r_if (i - 2) +
     &              delta_r_if (i - 1)) * 
     &              (delta_r_if (i - 2) + delta_r_if (i - 1)) *
     &              delta_r_if (i - 1) /
     &              ((- delta_r_if (i - 4)) * (- delta_r_if (i - 4) - 
     &              delta_r_if (i - 3)) *
     &              (- delta_r_if (i - 4) - delta_r_if (i - 3) -
     &              delta_r_if (i - 2))) *
     &              q (i - 4, j, k)
               enddo

         enddo

      endif
      
      return

      END SUBROUTINE extrapolate_out

c     ==================================================================



c     ==================================================================

      SUBROUTINE hydrostatic_atmosphere

c     ==================================================================
c
c      Imposes density profile of hydrostatic atmosphere.
c
c      Author: B.Mueller
c      Date:   31.5.2007
c
c     ------------------------------------------------------------------

      USE precision
      
      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc

      IMPLICIT NONE
      
      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = m + 1, m + 4
            v_1 (i, j, k) = v_1 (m, j, k) * (r (m) / r (i))
            rho (i, j, k) = rho (m, j, k) * (r (m) / r (i)) ** 3
            p   (i, j, k) = p   (m, j, k) * (r (m) / r (i)) ** 4
            eps (i, j, k) = eps (m, j, k) * (r (m) / r (i))
#ifdef CFC_MHD
            psi (i, j, k) = 0.0_rk!psi (m, j, k) * (r (m) / r (i))
#endif /* CFC_MHD */            
         end do

      end do

      return

      END SUBROUTINE hydrostatic_atmosphere

c     ==================================================================



#ifdef MPI_HYDRO
c     ==================================================================

      SUBROUTINE set_ghost_zones_mpi

c     ==================================================================
c
c      Author: B.Mueller
c      Date:   15.2.2010
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc      

      USE mo_mpi

      use configure
      IMPLICIT NONE

#ifdef CFC_MHD
      real (kind=rk) :: sbufy (-3:m+4, 1:4, o_s:o_e, 11+config%qn)
      real (kind=rk) :: rbufy (-3:m+4, 1:4, o_s:o_e, 11+config%qn)
      real (kind=rk) :: sbufz (-3:m+4, n_s-4:n_e+4, 1:4, 11+config%qn)
      real (kind=rk) :: rbufz (-3:m+4, n_s-4:n_e+4, 1:4, 11+config%qn)
#else /* CFC_MHD */
      real (kind=rk) :: sbufy (-3:m+4, 1:4, o_s:o_e, 7+config%qn)
      real (kind=rk) :: rbufy (-3:m+4, 1:4, o_s:o_e, 7+config%qn)
      real (kind=rk) :: sbufz (-3:m+4, n_s-4:n_e+4, 1:4, 7+config%qn)
      real (kind=rk) :: rbufz (-3:m+4, n_s-4:n_e+4, 1:4, 7+config%qn)
#endif /* CFC_MHD */

      integer (kind=ik) :: cycles, i

      integer (kind=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      integer (kind=ik), parameter :: tag_ysndl=1000_ik
      integer (kind=ik), parameter :: tag_ysndr=1001_ik
      integer (kind=ik), parameter :: tag_zsndl=1002_ik
      integer (kind=ik), parameter :: tag_zsndr=1003_ik

      integer :: nvars


c     number of variables that need to be communicated
#ifdef CFC_MHD
      nvars = 11_ik+config%qn
#else
      nvars = 7_ik+config%qn
#endif /* CFC_MHD */

c     ------------------------------------------------------------------
c     determine required number of next-neighbour-communications:
c     if any domain consists of less than 4 (4*4 in 3D) zones,
c     we may data from non-immediate neighbours. Repeating the
c     send-receive operation takes care of this problem in a
c     straightforward manner, although this may not be the most
c     efficient solution.


      if (config%nsdim .lt. 3) then
         if (nprocs*4_ik .le. n) then
            cycles=1_ik
         else if (nprocs*2_ik .le. n) then
            cycles=2_ik
         else
            cycles=4_ik
         end if
      else
         if (nprocs*16_ik .le. n*o) then
            cycles=1_ik
         else if (nprocs*4_ik .le. n*o) then
            cycles=2_ik
         else
            cycles=4_ik
         end if
      endif
c      print *,'cycles',cycles

c     ------------------------------------------------------------------
c     1st step: send to LEFT neighbour in THETA direction

      do i=1,cycles
c         print *,'cycle i',i
c     copy to buffer
         sbufy(-3:m+4,1:4,o_s:o_e,1)=
     &        rho(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,2)=
     &        v_1(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,3)=
     &        v_2(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,4)=
     &        v_3(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,5)=
     &        eps(-3:m+4,n_s:n_s+3,o_s:o_e)
         
         sbufy(-3:m+4,1:4,o_s:o_e,6:5+config%qn)=
     &        xnnu(-3:m+4,n_s:n_s+3,o_s:o_e,1:config%qn)
         
         sbufy(-3:m+4,1:4,o_s:o_e,6+config%qn)=
     &        p   (-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,7+config%qn)=
     &        gamm(-3:m+4,n_s:n_s+3,o_s:o_e)

#ifdef CFC_MHD
         sbufy(-3:m+4,1:4,o_s:o_e,8+config%qn)=
     &        b_1(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,9+config%qn)=
     &        b_2(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,10+config%qn)=
     &        b_3(-3:m+4,n_s:n_s+3,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,11+config%qn)=
     &        psi(-3:m+4,n_s:n_s+3,o_s:o_e)
#endif /* CFC_MHD */

c     determine source and destination ranks and send data
c         print *,'cart_shift 1'
         call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
         
         call MPI_Sendrecv (sbufy,(m+8_ik)*4_ik*o_loc*nvars,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy,
     &        (m+8_ik)*4_ik*o_loc*nvars,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndl,cart_comm,mpistat,ierr)

c     copy back to CoCoNuT arrays
         rho(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,1)
         v_1(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,2)
         v_2(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,3)
         v_3(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,4)
         eps(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,5)
         
         xnnu(-3:m+4,n_e+1:n_e+4,o_s:o_e,1:config%qn)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,6:5+config%qn)
         
         p   (-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,6+config%qn)
         gamm(-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,7+config%qn)

#ifdef CFC_MHD
         b_1  (-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,8+config%qn)
         b_2  (-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,9+config%qn)
         b_3  (-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,10+config%qn)
         psi  (-3:m+4,n_e+1:n_e+4,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,11+config%qn)
#endif /* CFC_MHD */

      end do

      if (src .ne. MPI_PROC_NULL) then
#ifdef CFC_TRANSPORT2
         h   (1:m,n_e+1:n_e+4,o_s:o_e)=
     &        1.0_rk + eps (1:m,n_e+1:n_e+4,o_s:o_e) +
     &        p (1:m,n_e+1:n_e+4,o_s:o_e) /
     &        rho (1:m,n_e+1:n_e+4,o_s:o_e) 
         c_sound_squared (1:m,n_e+1:n_e+4,o_s:o_e)=
     &        gamm (1:m,n_e+1:n_e+4,o_s:o_e) *
     &        p (1:m,n_e+1:n_e+4,o_s:o_e) /
     &        (rho (1:m,n_e+1:n_e+4,o_s:o_e) *
     &        h (1:m,n_e+1:n_e+4,o_s:o_e))
#else 
         c_sound_squared (1:m,n_e+1:n_e+4,o_s:o_e)=
     &        gamm (1:m,n_e+1:n_e+4,o_s:o_e) *
     &        p (1:m,n_e+1:n_e+4,o_s:o_e) /
     &        rho (1:m,n_e+1:n_e+4,o_s:o_e)
#endif /* CFC_TRANSPORT2 */
      end if


c     ------------------------------------------------------------------
c     2nd step: send to RIGHT neighbour in THETA direction

      do i=1,cycles

c     copy to buffer
         sbufy(-3:m+4,1:4,o_s:o_e,1)=
     &        rho(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,2)=
     &        v_1(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,3)=
     &        v_2(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,4)=
     &        v_3(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,5)=
     &        eps(-3:m+4,n_e-3:n_e,o_s:o_e)
         
         sbufy(-3:m+4,1:4,o_s:o_e,6:5+config%qn)=
     &        xnnu(-3:m+4,n_e-3:n_e,o_s:o_e,1:config%qn)
         
         sbufy(-3:m+4,1:4,o_s:o_e,6+config%qn)=
     &        p   (-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,7+config%qn)=
     &        gamm(-3:m+4,n_e-3:n_e,o_s:o_e)

#ifdef CFC_MHD
         sbufy(-3:m+4,1:4,o_s:o_e,8+config%qn)=
     &        b_1(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,9+config%qn)=
     &        b_2(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,10+config%qn)=
     &        b_3(-3:m+4,n_e-3:n_e,o_s:o_e)
         sbufy(-3:m+4,1:4,o_s:o_e,11+config%qn)=
     &        psi(-3:m+4,n_e-3:n_e,o_s:o_e)
#endif /* CFC_MHD */

c     determine source and destination ranks and send data
c         print *,'cart_shift 2'
         call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
         
         call MPI_Sendrecv (sbufy,(m+8_ik)*4_ik*o_loc*nvars,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndr,rbufy,
     &        (m+8_ik)*4_ik*o_loc*nvars,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndr,cart_comm,mpistat,ierr)

c     copy back to CoCoNuT arrays
         rho(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,1)
         v_1(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,2)
         v_2(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,3)
         v_3(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,4)
         eps(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,5)
         
         xnnu(-3:m+4,n_s-4:n_s-1,o_s:o_e,1:config%qn)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,6:5+config%qn)
         
         p   (-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,6+config%qn)
         gamm(-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,7+config%qn)

#ifdef CFC_MHD
         b_1  (-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,8+config%qn)
         b_2  (-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,9+config%qn)
         b_3  (-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,10+config%qn)
         psi  (-3:m+4,n_s-4:n_s-1,o_s:o_e)=
     &        rbufy(-3:m+4,1:4,o_s:o_e,11+config%qn)
#endif /* CFC_MHD */

      end do

      if (src .ne. MPI_PROC_NULL) then
#ifdef CFC_TRANSPORT2
         h   (1:m,n_s-4:n_s-1,o_s:o_e)=
     &        1.0_rk + eps (1:m,n_s-4:n_s-1,o_s:o_e) +
     &        p (1:m,n_s-4:n_s-1,o_s:o_e) /
     &        rho (1:m,n_s-4:n_s-1,o_s:o_e)
         c_sound_squared (1:m,n_s-4:n_s-1,o_s:o_e)=
     &        gamm (1:m,n_s-4:n_s-1,o_s:o_e) *
     &        p (1:m,n_s-4:n_s-1,o_s:o_e) /
     &        (rho (1:m,n_s-4:n_s-1,o_s:o_e) *
     &        h (1:m,n_s-4:n_s-1,o_s:o_e))
#else 
         c_sound_squared (1:m,n_s-4:n_s-1,o_s:o_e)=
     &        gamm (1:m,n_s-4:n_s-1,o_s:o_e) *
     &        p (1:m,n_s-4:n_s-1,o_s:o_e) /
     &        rho (1:m,n_s-4:n_s-1,o_s:o_e)
#endif /* CFC_TRANSPORT2 */
      end if

      if (config%nsdim .eq. 3) then
c     ------------------------------------------------------------------
c     3rd step: send to LEFT neighbour in VARPHI direction
c               (note that sbufz and rbufz are larger than sbufy
c                and rbufy, because they also include ghost zones
c                in theta-direction)

         do i=1,cycles
            
c     copy to buffer
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,1)=
     &           rho(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,2)=
     &           v_1(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,3)=
     &           v_2(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,4)=
     &           v_3(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,5)=
     &           eps(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,6:5+config%qn)=
     &           xnnu(-3:m+4,n_s-4:n_e+4,o_s:o_s+3,1:config%qn)
            
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,6+config%qn)=
     &           p   (-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,7+config%qn)=
     &           gamm(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)

#ifdef CFC_MHD
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,8+config%qn)=
     &           b_1(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,9+config%qn)=
     &           b_2(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,10+config%qn)=
     &           b_3(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,11+config%qn)=
     &           psi(-3:m+4,n_s-4:n_e+4,o_s:o_s+3)
#endif /* CFC_MHD */
            
c     determine source and destination ranks and send data
            call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)

c            print *,'boundary3'
            if (mod (src, 2_ik) .eq. 1_ik) then           
               call MPI_Send (sbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &              cart_comm,ierr)
               call MPI_Recv (rbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &              mpistat,ierr)
            else
               call MPI_Recv (rbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &              mpistat,ierr)
               call MPI_Send (sbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &              cart_comm,ierr)
            end if
            
c     copy back to CoCoNuT arrays
            rho(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,1)
            v_1(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,2)
            v_2(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,3)
            v_3(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,4)
            eps(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,5)
            
            xnnu(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4,1:config%qn)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,6:5+config%qn)
            
            p   (-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,6+config%qn)
            gamm(-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,7+config%qn)

#ifdef CFC_MHD
            b_1 (-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,8+config%qn)
            b_2 (-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,9+config%qn)
            b_3 (-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,10+config%qn)
            psi (-3:m+4,n_s-4:n_e+4,o_e+1:o_e+4)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,11+config%qn)
#endif /* CFC_MHD */
            
         end do
         
         if (src .ne. MPI_PROC_NULL) then
#ifdef CFC_TRANSPORT2
            h   (1:m,n_s:n_e,o_e+1:o_e+4)=
     &           1.0_rk + eps (1:m,n_s:n_e,o_e+1:o_e+4) +
     &           p (1:m,n_s:n_e,o_e+1:o_e+4) /
     &           rho (1:m,n_s:n_e,o_e+1:o_e+4)
            c_sound_squared (1:m,n_s:n_e,o_e+1:o_e+4)=
     &           gamm (1:m,n_s:n_e,o_e+1:o_e+4) *
     &           p (1:m,n_s:n_e,o_e+1:o_e+4) /
     &           (rho (1:m,n_s:n_e,o_e+1:o_e+4) *
     &           h (1:m,n_s:n_e,o_e+1:o_e+4))
#else 
            c_sound_squared (1:m,n_s:n_e,o_e+1:o_e+4)=
     &           gamm (1:m,n_s:n_e,o_e+1:o_e+4) *
     &           p (1:m,n_s:n_e,o_e+1:o_e+4) /
     &           rho (1:m,n_s:n_e,o_e+1:o_e+4)
#endif /* CFC_TRANSPORT2 */
         end if
         
         
c     ------------------------------------------------------------------
c     4th step: send to RIGHT neighbour in VARPHI direction
         
         do i=1,cycles
            
c     copy to buffer
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,1)=
     &           rho(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,2)=
     &           v_1(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,3)=
     &           v_2(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,4)=
     &           v_3(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,5)=
     &           eps(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,6:5+config%qn)=
     &           xnnu(-3:m+4,n_s-4:n_e+4,o_e-3:o_e,1:config%qn)
            
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,6+config%qn)=
     &           p   (-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,7+config%qn)=
     &           gamm(-3:m+4,n_s-4:n_e+4,o_e-3:o_e)

#ifdef CFC_MHD
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,8+config%qn)=
     &           b_1 (-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,9+config%qn)=
     &           b_2 (-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,10+config%qn)=
     &           b_3 (-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
            sbufz(-3:m+4,n_s-4:n_e+4,1:4,11+config%qn)=
     &           psi (-3:m+4,n_s-4:n_e+4,o_e-3:o_e)
#endif /* CFC_MHD */
            
c     determine source and destination ranks and send data
            call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
            
c            print *,'boundary4'
            if (mod (src, 2_ik) .eq. 1_ik) then           
               call MPI_Send (sbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &              cart_comm,ierr)
               call MPI_Recv (rbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &              mpistat,ierr)
            else
               call MPI_Recv (rbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &              mpistat,ierr)
               call MPI_Send (sbufz,(m+8_ik)*(n_loc+8_ik)*4_ik*nvars,
     &              MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &              cart_comm,ierr)
            end if
            
c     copy back to CoCoNuT arrays
            rho(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,1)
            v_1(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,2)
            v_2(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,3)
            v_3(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,4)
            eps(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,5)
            
            xnnu(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1,1:config%qn)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,6:5+config%qn)
            
            p   (-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,6+config%qn)
            gamm(-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,7+config%qn)

#ifdef CFC_MHD
            b_1 (-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,8+config%qn)
            b_2 (-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,9+config%qn)
            b_3 (-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,10+config%qn)
            psi (-3:m+4,n_s-4:n_e+4,o_s-4:o_s-1)=
     &           rbufz(-3:m+4,n_s-4:n_e+4,1:4,11+config%qn)
#endif /* CFC_MHD */

         end do
         
         if (src .ne. MPI_PROC_NULL) then
#ifdef CFC_TRANSPORT2
            h   (1:m,n_s:n_e,o_s-4:o_s-1)=
     &           1.0_rk + eps (1:m,n_s:n_e,o_s-4:o_s-1) +
     &           p (1:m,n_s:n_e,o_s-4:o_s-1) /
     &           rho (1:m,n_s:n_e,o_s-4:o_s-1)
            c_sound_squared (1:m,n_s:n_e,o_s-4:o_s-1)=
     &           gamm (1:m,n_s:n_e,o_s-4:o_s-1) *
     &           p (1:m,n_s:n_e,o_s-4:o_s-1) /
     &           (rho (1:m,n_s:n_e,o_s-4:o_s-1) *
     &           h (1:m,n_s:n_e,o_s-4:o_s-1))
#else 
            c_sound_squared (1:m,n_s:n_e,o_s-4:o_s-1)=
     &           gamm (1:m,n_s:n_e,o_s-4:o_s-1) *
     &           p (1:m,n_s:n_e,o_s-4:o_s-1) /
     &           rho (1:m,n_s:n_e,o_s-4:o_s-1)
#endif /* CFC_TRANSPORT2 */
         end if

      endif ! config%nsdim .eq. 3

      return

      END SUBROUTINE set_ghost_zones_mpi

c     ==================================================================



#ifdef CFC_MHD
c     ==================================================================

      SUBROUTINE set_ghost_zones_induction_mpi

c     ==================================================================
c
c      Author: B.Mueller
c      Date:   16.2.2017
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE fluxes_cfc

      USE mo_mpi

      use configure
      IMPLICIT NONE

      real (kind=rk) :: sbufy (0:m+1, 1, o_s-1:o_e, 1:2, 1:2)
      real (kind=rk) :: rbufy (0:m+1, 1, o_s-1:o_e, 1:2, 1:2)
      real (kind=rk) :: sbufz (0:m+1, n_s-1:n_e+1, 1, 1:2, 1:2)
      real (kind=rk) :: rbufz (0:m+1, n_s-1:n_e+1, 1, 1:2, 1:2)

      integer (kind=ik) :: cycles, i

      integer (kind=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      integer (kind=ik), parameter :: tag_ysndl=2000_ik
      integer (kind=ik), parameter :: tag_ysndr=2001_ik
      integer (kind=ik), parameter :: tag_zsndl=2002_ik
      integer (kind=ik), parameter :: tag_zsndr=2003_ik

c     ------------------------------------------------------------------
c     determine required number of next-neighbour-communications:
c     if any domain consists of less than 4 (4*4 in 3D) zones,
c     we may data from non-immediate neighbours. Repeating the
c     send-receive operation takes care of this problem in a
c     straightforward manner, although this may not be the most
c     efficient solution.
c     ------------------------------------------------------------------
c     1st step: send to LEFT neighbour in THETA direction

      sbufy (0:m+1, 1, o_s:o_e, 1:2, 1)=
     &     emf_face (0:m+1, n_s, o_s:o_e, 1:2, 1)
      sbufy (0:m+1, 1, o_s-1:o_e, 1:2, 2)=
     &     emf_face (0:m+1, n_s, o_s-1:o_e, 1:2, 3)

      call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
      
      call MPI_Sendrecv (sbufy,(m+2_ik)*(o_loc+1)*4_ik,
     &     MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy,
     &     (m+2_ik)*(o_loc+1)*4_ik,MPI_DOUBLE_PRECISION,
     &     src,tag_ysndl,cart_comm,mpistat,ierr)

c     copy back to CoCoNuT arrays
      if (src .ne. MPI_PROC_NULL) then
         emf_face (0:m+1, n_e+1, o_s:o_e, 1:2, 1)=
     &        rbufy(0:m+1, 1, o_s:o_e, 1:2, 1)
         emf_face (0:m+1, n_e+1, o_s-1:o_e, 1:2, 3)=
     &        rbufy(0:m+1, 1, o_s-1:o_e, 1:2, 2)
      end if
c     ------------------------------------------------------------------
c     2nd step: send to RIGHT neighbour in THETA direction

c     copy to buffer
      sbufy (0:m+1, 1, o_s:o_e, 1:2, 1)=
     &     emf_face (0:m+1, n_e, o_s:o_e, 1:2, 1)
      sbufy (0:m+1, 1, o_s-1:o_e, 1:2, 2)=
     &     emf_face (0:m+1, n_e, o_s-1:o_e, 1:2, 3)
      
      call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
      
      call MPI_Sendrecv (sbufy,(m+2_ik)*(o_loc+1)*4_ik,
     &     MPI_DOUBLE_PRECISION,dest,tag_ysndr,rbufy,
     &     (m+2_ik)*(o_loc+1)*4_ik,MPI_DOUBLE_PRECISION,
     &     src,tag_ysndr,cart_comm,mpistat,ierr)

c     copy back to CoCoNuT arrays
      if (src .ne. MPI_PROC_NULL) then
         emf_face (0:m+1, n_s-1, o_s:o_e, 1:2, 1)=
     &        rbufy(0:m+1, 1, o_s:o_e, 1:2, 1)
         emf_face (0:m+1, n_s-1, o_s-1:o_e, 1:2, 3)=
     &        rbufy(0:m+1, 1, o_s-1:o_e, 1:2, 2)
      end if

      if (config%nsdim .eq. 3) then
c     ------------------------------------------------------------------
c     3rd step: send to LEFT neighbour in VARPHI direction
c               (note that sbufz and rbufz are larger than sbufy
c                and rbufy, because they also include ghost zones
c                in theta-direction)

c     copy to buffer
         sbufz(0:m+1,n_s-1:n_e+1,1,1:2,1:2)=
     &        emf_face(0:m+1,n_s-1:n_e+1,o_s,1:2,1:2)
c         print *,'a'
c     determine source and destination ranks and send data
         call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
c         print *,'b',src,dest

         if (mod (src, 2_ik) .eq. 1_ik) then           
            call MPI_Send (sbufz, (m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &           cart_comm,ierr)
c            print *,'c',(m+2_ik)*(n_loc+2_ik)*4_ik,
c     &           size(rbufz)
            call MPI_Recv  (rbufz,(m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &           mpistat,ierr)
c            print *,'d'
         else
            call MPI_Recv (rbufz, (m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &           mpistat,ierr)
c            print *,'c',(m+2_ik)*(n_loc+2_ik)*4_ik,size(sbufz)
            call MPI_Send (sbufz, (m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &           cart_comm,ierr)
c            print *,'d'
         end if
            
c     copy back to CoCoNuT arrays
         if (src .ne. MPI_PROC_NULL) then
            emf_face(0:m+1,n_s-1:n_e+1,o_e+1,1:2,1:2)=
     &           sbufz(0:m+1,n_s-1:n_e+1,1,1:2,1:2)
c         print *,'e'
         end if
c     ------------------------------------------------------------------
c     4th step: send to RIGHT neighbour in VARPHI direction
         
c     copy to buffer
         sbufz(0:m+1,n_s-1:n_e+1,1,1:2,1:2)=
     &        emf_face(0:m+1,n_s-1:n_e+1,o_e,1:2,1:2)

c     determine source and destination ranks and send data
         call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
         
         if (mod (src, 2_ik) .eq. 1_ik) then           
            call MPI_Send (sbufz,(m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
            call MPI_Recv (rbufz,(m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
         else
            call MPI_Recv (rbufz,(m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
            call MPI_Send (sbufz,(m+2_ik)*(n_loc+2_ik)*4_ik,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
         end if
         
c     copy back to CoCoNuT arrays
         if (src .ne. MPI_PROC_NULL) then
            emf_face(0:m+1,n_s-1:n_e+1,o_s-1,1:2,1:2)=
     &           sbufz(0:m+1,n_s-1:n_e+1,1,1:2,1:2)
         end if
c         print *,'emf_bnd1',emf_face(90,n_s+2,o_s-1:o_e+1,1,2)
c         print *,'emf_bnd2',emf_face(90,n_s+2,o_s-1:o_e,2,3)
      endif ! config%nsdim .eq. 3

      
      return

      END SUBROUTINE set_ghost_zones_induction_mpi

c     ==================================================================
#endif /* CFC_MHD */

#endif /*MPI_HYDRO*/




      END MODULE gr_boundary

c     ==================================================================

c     ==================================================================

      MODULE poisson

c     ==================================================================

c      USE fft_mpi

      IMPLICIT NONE
      
      CONTAINS



c     ==================================================================

      SUBROUTINE comp_legendre

c     ==================================================================
c     Berechnet die Legendre-Polynome P_l(cos(theta_if)) an den 
c     Zellraendern, den Zellzentren, sowie die Integrale von P_l
c     ueber eine Zelle in theta-Richtung.
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE legendre
      USE parameters_cfc, ONLY: equatorial_symmetry

      USE configure

      IMPLICIT NONE

      integer (kind=ik) :: j,l



#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL
!$OMP& PRIVATE(j,l)
#endif
      if (config%nsdim .ge. 2) then
c     ------------------------------------------------------------------
c     Legendre polynomials ought to be availably LOCALLY on the stack

      p_leg(:,0)=1.0_rk
      p_leg(0:n,1)=cos(theta_if(0:n))
      if (equatorial_symmetry) p_leg(n,1)=0.0_rk

      p_leg_c(-1:n+2,0)=1.0_rk
      p_leg_c(-1:n+2,1)=cos(theta(-1:n+2))

      p1_leg(0:n,1)=-sin(theta_if(0:n))
      p1_leg(0:n,2)=-3*sin(theta_if(0:n))*cos(theta_if(0:n))

      p1_leg_c(-1:n+2,1)=-sin(theta(-1:n+2))
      p1_leg_c(-1:n+2,2)=-3*sin(theta(-1:n+2))*cos(theta(-1:n+2))

c     Rekursionsformel fuer die Legendre-Polynome:
c     (Jackson, Classical Electrodynamics, aeq. 3.29):
      do l=1,nleg
         p_leg  (:,l+1)=(p_leg  (:,l)*p_leg  (:,1)*(2*l+1)-
     &        l*p_leg  (:,l-1))/dble(l+1)
         p_leg_c(:,l+1)=(p_leg_c(:,l)*p_leg_c(:,1)*(2*l+1)-
     &        l*p_leg_c(:,l-1))/dble(l+1)
      end do

      do l=2,nleg
         p1_leg  (:,l+1)=(p1_leg  (:,l)*p_leg  (:,1)*(2*l+1)-
     &        (l+1)*p1_leg  (:,l-1))/dble(l)
         p1_leg_c(:,l+1)=(p1_leg_c(:,l)*p_leg_c(:,1)*(2*l+1)-
     &        (l+1)*p1_leg_c(:,l-1))/dble(l)
      end do

c     Exakte Werte an den Raendern setzen (nicht unbedingt noetig):
      do l=0,nleg+2
         if (.not.equatorial_symmetry) then
            if (mod(l,2) .eq. 0) then
               p_leg(n,l)=1.0_rk
            else
               p_leg(n,l)=-1.0_rk
            end if
         end if
         p_leg(0,l)=1.0_rk
      end do

c     Integrale ueber theta-Zelle mit Hilfe der Formel fuer die
c     Ableitungen der Legendre-Polynome (ibid., aeq.3.28)
      do j=1,n
         p_leg_int(j,0)=p_leg(j,1)-p_leg(j-1,1)
      end do
      do l=1,nleg
         do j=1,n
            p_leg_int(j,l)=
     &        ((p_leg(j  ,l+1)-p_leg(j  ,l-1))-
     &         (p_leg(j-1,l+1)-p_leg(j-1,l-1)))/
     &        dble(2*l+1)
         end do
      end do

      do j=1,n
         p1_leg_int(j,1)=
     &        -(sin(2*theta_if(j))-sin(2*theta_if(j-1)))/4.0_rk+
     &        (theta_if(j)-theta_if(j-1))/2.0_rk
         p1_leg_int(j,2)=
     &        sin(theta_if(j))**3-sin(theta_if(j-1))**3
      end do
      do l=3,nleg
         do j=1,n
            p1_leg_int(j,l)=
     &           (l-2)*l/dble((l+1)*(l-1))*
     &           p1_leg_int(j,l-2)-
     &           (2*l-1)/dble((l+1)*(l-1))*
     &           (sin(theta_if(j  ))**2*p1_leg(j  ,l-1)-
     &           sin(theta_if(j-1))**2*p1_leg(j-1,l-1))
         end do
      end do

#ifdef ENFORCE_SPHERICAL_METRIC
c     Renormalisieren:
      p1_leg_int(1:n,1)=p1_leg_int(1:n,1)*( 4.0_rk/3.0_rk)/
     &     sum(sin(theta(1:n))*p1_leg_int(1:n,1))
      p_leg_int (1:n,1)=p_leg_int (1:n,1)*(-2.0_rk/3.0_rk)/
     &     sum(cos(theta(1:n))*p_leg_int(1:n,1))
      if (equatorial_symmetry) then
         p1_leg_int(1:n,1)=p1_leg_int(1:n,1)*0.5_rk
         p_leg_int (1:n,1)=p_leg_int (1:n,1)*0.5_rk
      end if
      p1_leg_int(1:n,2:nleg+1)=0.0_rk
      p_leg_int (1:n,2:nleg+1)=0.0_rk
      p1_leg_c  (1:n,2:nleg+1)=0.0_rk
      p_leg_c   (1:n,2:nleg+1)=0.0_rk
#endif

      else                      ! config%nsidm .ge. 2 
      p_leg    (:,0)= 1.0_rk
      p_leg_c  (:,0)= 1.0_rk
      p_leg_int(:,0)=-2.0_rk
      endif 

c     Ghost zones:
      p_leg_c(  0,:)=p_leg_c(1  ,:)
      p_leg_c( -1,:)=p_leg_c(2  ,:)
      p_leg_c(n+1,:)=p_leg_c(n  ,:)
      p_leg_c(n+2,:)=p_leg_c(n-1,:)

      p1_leg_c(  0,:)=-p1_leg_c(1  ,:)
      p1_leg_c( -1,:)=-p1_leg_c(2  ,:)
      p1_leg_c(n+1,:)=-p1_leg_c(n  ,:)
      p1_leg_c(n+2,:)=-p1_leg_c(n-1,:)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#endif

      return

      END SUBROUTINE comp_legendre

c     ==================================================================




c     ==================================================================

      SUBROUTINE hilfsgroessen_r

c     ==================================================================

      USE precision

      USE grid_cfc
      USE mod_hilfsgroessen_r

      IMPLICIT NONE

      integer (kind=ik) :: i,l

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL
!$OMP& PRIVATE(i,l)
#endif
      do l=0,nleg+2
         r_hoch_l(1:m,l)=r(1:m)**l
         r_hoch_l_r(1:m+1,l)=1.0_rk/(r(1:m+1)**l)
      end do

      do l=0,nleg+1
         do i=1,m
            gew_r_in(i,l)=(r_if(i)**3-r_if(i-1)**3)/3.0_rk*r(i)**l
         end do

         do i=1,m
            gew_r_ex(i,l)=(r_if(i)**3-r_if(i-1)**3)/3.0_rk/r(i)**(l+1)
         end do
c         do i=1,m
c            gew_r_in(i,l)=(r_if(i)**(l+3)-r_if(i-1)**(l+3))/
c     &           dble(l+3)
c         end do
c
c         do i=1,m
c            gew_r_ex(i,l)=(r_if(i)**(2-l)-r_if(i-1)**(2-l))/
c     &           dble(2-l)
c         end do
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL 
#endif

      return

      END SUBROUTINE hilfsgroessen_r

c     ==================================================================



 
c     ==================================================================

      SUBROUTINE loese_poisson (rho,phi)

c     ==================================================================
c
c     Solves the Poisson equation:
c     div grad phi = 4 pi rho
c
c     Remarks:
c     - Uses multipole expansion (cp. Mueller & Steinmetz)
c     - Parallelized over multipole / angles.
c     
c     rho: source density at CELL CENTRES
c     phi: solution at CELL CENTRES
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE legendre
      USE mod_hilfsgroessen_r
      USE parameters_cfc, ONLY: equatorial_symmetry

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      USE configure

      IMPLICIT NONE

      real (kind=rk), intent (in) :: 
     &     rho (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk), intent (out) :: 
     &     phi (0:m+1,n_s-1:n_e+1,o_s-1:o_e+1)

      integer (kind=ik) :: i,j,k,jk,l,lp1,lstep

#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

      real (kind=rk) :: prefac

      if (.not.equatorial_symmetry) then
         lstep=1
         prefac=0.5_rk
      else
         lstep=2
         prefac=1.0_rk
      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=0,nleg,lstep

         lp1=l+1

c     ------------------------------------------------------------------
c     Compute multipole moments: 

         mom_rho(:,l)=0.0_rk
         do k=o_s,o_e
            do j=n_s,n_e
               do i=1,m
                  mom_rho(i,l)=mom_rho(i,l)+rho(i,j,k)*p_leg_int(j,l)
               end do
            end do
         end do

         if (config%nsdim .eq. 3) then
            mom_rho(:,l)=mom_rho(:,l)/real(o,kind=rk)
         endif

#ifdef MPI_HYDRO
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      call MPI_Allreduce (mom_rho,mom_rho2,m*(nleg+2_ik),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,cart_comm,ierr)
      mom_rho=mom_rho2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=0,nleg,lstep      
#endif /*MPI_HYDRO*/
         
c     ------------------------------------------------------------------
c     Compute inner/outer potentials: 
c     ------------------------------------------------------------------
c     Note: the following FD representation may look like an
c     asymmetric 1st-order scheme, but is actually the 2nd-order
c     trapezoid method in disguise. Phi_in and Phi_ex are computed
c     at CELL CENTRES.
c     ------------------------------------------------------------------

         lp1=l+1

         phi_in(1,l)=mom_rho(1,l)*gew_r_in(1,l)*r_hoch_l_r(1,lp1)
         do i=2,m
            phi_in(i,l)=(phi_in(i-1,l)*r_hoch_l(i-1,lp1)+
     &           mom_rho(i,l)*gew_r_in(i,l))*
     &           r_hoch_l_r(i,lp1)
         end do
         phi_in(m+1,l)=phi_in(m,l)*r_hoch_l(m,lp1)*
     &        r_hoch_l_r(m+1,lp1)
         
         
         phi_ex(m,l)=0.0_rk
         do i=m-1,1,-1
            phi_ex(i,l)=(phi_ex(i+1,l)*r_hoch_l_r(i+1,l)+
     &           mom_rho(i+1,l)*gew_r_ex(i+1,l))*
     &           r_hoch_l(i,l)
         end do

#ifdef ENFORCE_SPHERICAL_METRIC
         if (l .gt. 0) then     !only monpole survives
            phi_ex(:,l)=0.0_rk
            phi_in(:,l)=0.0_rk
         end if
#endif
      end do

c     ------------------------------------------------------------------
c     Construct the solution:

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         phi(:,j,k)=0.0_rk
         
         do l=0,nleg,lstep
            do i=1,m
               phi(i,j,k)=phi(i,j,k)+(phi_in(i,l)+phi_ex(i,l))*
     &              p_leg_c(j,l)
            end do
            phi(m+1,j,k)=phi(m+1,j,k)+phi_in(m,l)*p_leg_c(j,l)
         end do

c     ------------------------------------------------------------------
c     Add correct pre-factors:
         phi(:,j,k)=phi(:,j,k)*prefac

c     ------------------------------------------------------------------
c     Fill remaining ghost zones:     
         phi(0,j,k)=phi(1,j,k)

         if (j .eq. 1) then     !left domain and grid boundary
            phi(:,j-1,k)=phi(:,j,k)
         else if (j. eq. n_s) then !only left domain boundary
            phi(:,j-1,k)=0.0_rk
            do l=0,nleg,lstep
               do i=1,m
                  phi(i,j-1,k)=phi(i,j-1,k)+(phi_in(i,l)+phi_ex(i,l))*
     &                 p_leg_c(j-1,l)
               end do
               phi(m+1,j-1,k)=phi(m+1,j-1,k)+
     &                           phi_in(m,l)*p_leg_c(j-1,l)
            end do
            phi(:,j-1,k)=phi(:,j-1,k)*prefac
         endif

         if (j .eq. n) then     !right domain and grid boundarz
            phi(:,j+1,k)=phi(:,j,k)
         else if (j. eq. n_e) then !only right domain boundary
            phi(:,j+1,k)=0.0_rk
            do l=0,nleg,lstep
               do i=1,m
                  phi(i,j+1,k)=phi(i,j+1,k)+(phi_in(i,l)+phi_ex(i,l))*
     &                 p_leg_c(j+1,l)
               end do
               phi(m+1,j+1,k)=phi(m+1,j+1,k)+ 
     &                           phi_in(m,l)*p_leg_c(j+1,l)
            end do
            phi(:,j+1,k)=phi(:,j+1,k)*prefac
         end if

         if (k .eq. 1) then     !left domain and grid boundary
            phi(:,j,k-1)=phi(:,j,k)
         else if (k. eq. o_s) then !only left domain boundary
            phi(:,j,k-1)=0.0_rk
            do l=0,nleg,lstep
               do i=1,m
                  phi(i,j,k-1)=phi(i,j,k-1)+(phi_in(i,l)+phi_ex(i,l))*
     &                 p_leg_c(j,l)
               end do
               phi(m+1,j,k-1)=phi(m+1,j,k-1)+
     &                           phi_in(m,l)*p_leg_c(j,l)
            end do
            phi(:,j,k-1)=phi(:,j,k-1)*prefac
         endif

         if (k .eq. o) then     !right domain and grid boundarz
            phi(:,j,k+1)=phi(:,j,k)
         else if (k. eq. o_e) then !only right domain boundary
            phi(:,j,k+1)=0.0_rk
            do l=0,nleg,lstep
               do i=1,m
                  phi(i,j,k+1)=phi(i,j,k+1)+(phi_in(i,l)+phi_ex(i,l))*
     &                 p_leg_c(j,l)
               end do
               phi(m+1,j,k+1)=phi(m+1,j,k+1)+
     &                           phi_in(m,l)*p_leg_c(j,l)
            end do
            phi(:,j,k+1)=phi(:,j,k+1)*prefac
         end if


      end do

      return
      
      END SUBROUTINE loese_poisson

c     ==================================================================


c#if 1 
c#if (DIMN>1)
#if !defined(ENFORCE_SPHERICAL_METRIC) && (DIMN>1)
c     ==================================================================

      SUBROUTINE loese_vektor_poisson (s_vec,beta)

c     ==================================================================
c
c     Solves the generalized vector Poisson equation:
c     div grad beta + 1/3 grad div beta= 4 pi S
c
c     Remarks:
c     - Multipole expansion of components in global Cartesian basis
c     - Uses the Shibata scheme to treat the grad div term
c     - Parallelized over multipole / angles.
c
c     s_vec: source at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c     beta:  solution at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c
c     WARNING: currently restricted to axisymmetric vector field!
c
c     ------------------------------------------------------------------

      USE precision

 !    USE size_cfc
      USE grid_cfc
      USE legendre
      USE mod_hilfsgroessen_r
      USE parameters_cfc, ONLY: equatorial_symmetry
      USE perm_aux_cfc

      USE mod_poisson_aux, ONLY: beta_x, beta_y, beta_z,
     &     s_x, s_y, s_z, s_tmp, chi

#ifdef MPI_HYDRO
      USE mo_mpi
#endif
      USE configure

      IMPLICIT NONE

      real (kind=rk), intent (in) ::
     &     s_vec (1:m,n_s:n_e,o_s:o_e,3)
      real (kind=rk), intent (out) :: 
     &     beta (0:m+1,n_s-1:n_e+1,o_s-1:o_e+1,3)

      real (kind=rk) :: llp1_r

      integer (kind=ik) :: i,j,k,jk,l,lp1,ji,lstep

#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

      real (kind=rk) :: prefac

      if (.not.equatorial_symmetry) then
         lstep=1
         prefac=0.5_rk
      else
         lstep=2
         prefac=1.0_rk
      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i=1,m

            s_x(i,j,k)=
     &           sinus_theta(j)     *s_vec(i,j,k,1)+
     &           cosin_theta(j)*r(i)*s_vec(i,j,k,2)
            s_y(i,j,k)=
     &           sinus_theta(j)*r(i)*s_vec(i,j,k,3)
            s_z(i,j,k)=
     &           cosin_theta(j)     *s_vec(i,j,k,1)-
     &           sinus_theta(j)*r(i)*s_vec(i,j,k,2)
            s_tmp(i,j,k)=-r(i)*s_vec(i,j,k,1)

         end do
      end do


c     ------------------------------------------------------------------
c     Expansion in associated Legendre functions P_l^1 for the
c     x- and -y component

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=1,nleg,lstep  !important: l has to start from 1

         lp1=l+1

c     ------------------------------------------------------------------
c     Compute multipole moments:
         mom_s(:,l,1)=0.0_rk
         mom_s(:,l,2)=0.0_rk
         do k=o_s,o_e
            do j=n_s,n_e
               do i=1,m
                  mom_s(i,l,1)=mom_s(i,l,1)+s_x(i,j,k)*p1_leg_int(j,l)
                  mom_s(i,l,2)=mom_s(i,l,2)+s_y(i,j,k)*p1_leg_int(j,l)
               end do
            end do
         end do

         if (config%nsdim .eq. 3) then
            mom_s(:,l,1)=mom_s(:,l,1)/dble(o)
            mom_s(:,l,2)=mom_s(:,l,2)/dble(o)
         endif

#ifdef MPI_HYDRO
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      call MPI_Allreduce (mom_s(:,:,1:2),mom_s2(:,:,1:2),
     &     2_ik*m*(nleg+1_ik),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,cart_comm,ierr)
      mom_s(:,:,1:2)=mom_s2(:,:,1:2)
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=0,nleg,lstep      
#endif /*MPI_HYDRO*/

c     ------------------------------------------------------------------
c     Compute inner/outer potentials: 

         lp1=l+1

         phi_inv(1,l,1)=mom_s(1,l,1)*gew_r_in(1,l)*r_hoch_l_r(1,lp1)
         phi_inv(1,l,2)=mom_s(1,l,2)*gew_r_in(1,l)*r_hoch_l_r(1,lp1)
         do i=2,m
            phi_inv(i,l,1)=(phi_inv(i-1,l,1)*r_hoch_l(i-1,lp1)+
     &           mom_s(i,l,1)*gew_r_in(i,l))*
     &           r_hoch_l_r(i,lp1)
            phi_inv(i,l,2)=(phi_inv(i-1,l,2)*r_hoch_l(i-1,lp1)+
     &           mom_s(i,l,2)*gew_r_in(i,l))*
     &           r_hoch_l_r(i,lp1)
         end do
         phi_inv(m+1,l,1)=phi_inv(m,l,1)*r_hoch_l(m,lp1)*
     &        r_hoch_l_r(m+1,lp1)
         phi_inv(m+1,l,2)=phi_inv(m,l,2)*r_hoch_l(m,lp1)*
     &        r_hoch_l_r(m+1,lp1)


         phi_exv(m,l,1)=0.0_rk
         phi_exv(m,l,2)=0.0_rk
         do i=m-1,1,-1
            phi_exv(i,l,1)=(phi_exv(i+1,l,1)*r_hoch_l_r(i+1,l)+
     &           mom_s(i+1,l,1)*gew_r_ex(i+1,l))*
     &           r_hoch_l(i,l)
            phi_exv(i,l,2)=(phi_exv(i+1,l,2)*r_hoch_l_r(i+1,l)+
     &           mom_s(i+1,l,2)*gew_r_ex(i+1,l))*
     &           r_hoch_l(i,l)
         end do

      end do

c     ------------------------------------------------------------------
c     Expansion in Legendre polynomials P_l for the
c     z-component and auxiliary potential chi:

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=0,nleg,lstep  !here l starts form 0, as for the
                                !scalar Poisson equation
         lp1=l+1

c     ------------------------------------------------------------------
c     Compute multipole moments: 

         mom_s(:,l,3)=0.0_rk
         mom_s(:,l,4)=0.0_rk
         do k=o_s,o_e
            do j=n_s,n_e
               do i=1,m
                  mom_s(i,l,3)=mom_s(i,l,3)+s_z(i,j,k)*p_leg_int(j,l)
                  mom_s(i,l,4)=mom_s(i,l,4)+s_tmp(i,j,k)*p_leg_int(j,l)
               end do
            end do
         end do

         if (config%nsdim .eq. 3) then
            mom_s(:,l,3)=mom_s(:,l,3)/dble(o)
            mom_s(:,l,4)=mom_s(:,l,4)/dble(o)
         ENDIF

#ifdef MPI_HYDRO
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      call MPI_Allreduce (mom_s(:,:,3:4),mom_s2(:,:,3:4),
     &     2_ik*m*(nleg+1_ik),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,cart_comm,ierr)
c      call MPI_Allreduce (MPI_IN_PLACE,mom_s(:,:,4),m*(nleg+1),
c     &     MPI_DOUBLE_PRECISION,MPI_SUM,cart_comm,ierr)
      mom_s(:,:,3:4)=mom_s2(:,:,3:4)
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do l=0,nleg,lstep      
#endif /*MPI_HYDRO*/

c     ------------------------------------------------------------------
c     Compute inner/outer potentials: 

         lp1=l+1

         phi_inv(1,l,3)=mom_s(1,l,3)*gew_r_in(1,l)*r_hoch_l_r(1,lp1)
         phi_inv(1,l,4)=mom_s(1,l,4)*gew_r_in(1,l)*r_hoch_l_r(1,lp1)
         do i=2,m
            phi_inv(i,l,3)=(phi_inv(i-1,l,3)*r_hoch_l(i-1,lp1)+
     &           mom_s(i,l,3)*gew_r_in(i,l))*
     &           r_hoch_l_r(i,lp1)
            phi_inv(i,l,4)=(phi_inv(i-1,l,4)*r_hoch_l(i-1,lp1)+
     &           mom_s(i,l,4)*gew_r_in(i,l))*
     &           r_hoch_l_r(i,lp1)
         end do
         phi_inv(m+1,l,3)=phi_inv(m,l,3)*
     &        (r_hoch_l(m,lp1)*r_hoch_l_r(m+1,lp1))
         phi_inv(m+1,l,4)=phi_inv(m,l,4)*
     &        (r_hoch_l(m,lp1)*r_hoch_l_r(m+1,lp1))

         phi_exv(m,l,3)=0.0_rk
         phi_exv(m,l,4)=0.0_rk
         do i=m-1,1,-1
            phi_exv(i,l,3)=(phi_exv(i+1,l,3)*r_hoch_l_r(i+1,l)+
     &           mom_s(i+1,l,3)*gew_r_ex(i+1,l))*
     &           r_hoch_l(i,l)
            phi_exv(i,l,4)=(phi_exv(i+1,l,4)*r_hoch_l_r(i+1,l)+
     &           mom_s(i+1,l,4)*gew_r_ex(i+1,l))*
     &           r_hoch_l(i,l)
         end do

      end do

#ifdef ENFORCE_SPHERICAL_METRIC
      phi_exv(:,2:nleg,1:4)=0.0_rk
      phi_inv(:,2:nleg,1:4)=0.0_rk
      phi_exv(:,0,3)=0.0_rk
      phi_inv(:,0,3)=0.0_rk
#endif

c     ------------------------------------------------------------------
c     Construct the solution vector field:

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         beta_x(:,j,k)=0.0_rk
         beta_y(:,j,k)=0.0_rk
         beta_z(:,j,k)=0.0_rk
         chi   (:,j,k)=0.0_rk

c     ------------------------------------------------------------------
c     Sum over all multipoles: 

         do l=1,nleg,lstep
            lp1=l+1
            llp1_r=1.0_rk/dble(lp1*l)
            do i=1,m
               beta_x(i,j,k)=beta_x(i,j,k)+
     &              (phi_inv(i,l,1)+phi_exv(i,l,1))*p1_leg_c(j,l)*
     &              llp1_r
               beta_y(i,j,k)=beta_y(i,j,k)+
     &              (phi_inv(i,l,2)+phi_exv(i,l,2))*p1_leg_c(j,l)*
     &              llp1_r
            end do
            beta_x(m+1,j,k)=beta_x(m+1,j,k)+
     &           phi_inv(m+1,l,1)*p1_leg_c(j,l)*llp1_r
            beta_y(m+1,j,k)=beta_y(m+1,j,k)+
     &           phi_inv(m+1,l,2)*p1_leg_c(j,l)*llp1_r
         end do

         do l=0,nleg,lstep
            do i=1,m
               beta_z(i,j,k)=beta_z(i,j,k)+ 
     &              (phi_inv(i,l,3)+phi_exv(i,l,3))*p_leg_c(j,l)
               chi   (i,j,k)=chi   (i,j,k)+
     &              (phi_inv(i,l,4)+phi_exv(i,l,4))*p_leg_c(j,l)
            end do
            beta_z(0,j,k)=beta_z(0,j,k)+ ! reflection symmetry
     &           phi_exv(1,l,3)*p_leg_c(n+1-j,l)
            beta_z(m+1,j,k)=beta_z(m+1,j,k)+
     &           phi_inv(m+1,l,3)*p_leg_c(j,l)
            chi   (m+1,j,k)=chi   (m+1,j,k)+
     &           phi_inv(m+1,l,4)*p_leg_c(j,l)
         end do

         beta_x(0,j,k)=-beta_x(1,j,k) ! MINUS sign:
                                      ! beta_x=0 on the rotation axis
         beta_y(0,j,k)=-beta_y(1,j,k) ! MINUS sign:
                                      ! beta_y=0 on the rotation axis
         chi   (0,j,k)= chi   (1,j,k)


c     ------------------------------------------------------------------
c     Add correct pre-factors
         beta_x(:,j,k)=beta_x(:,j,k)*prefac
         beta_y(:,j,k)=beta_y(:,j,k)*prefac
         beta_z(:,j,k)=beta_z(:,j,k)*prefac
         chi   (:,j,k)=chi(:,j,k)   *prefac

c     ------------------------------------------------------------------
c     ghost zones:
         if (j .le. n_s+1) then

            ji=2*n_s-1-j

            beta_x(:,ji,k)=0.0_rk
            beta_y(:,ji,k)=0.0_rk
            beta_z(:,ji,k)=0.0_rk
            chi   (:,ji,k)=0.0_rk

            do l=1,nleg,lstep
               lp1=l+1
               llp1_r=1.0_rk/dble(lp1*l)
               do i=1,m
                  beta_x(i,ji,k)=beta_x(i,ji,k)+
     &                 (phi_inv(i,l,1)+phi_exv(i,l,1))*p1_leg_c(ji,l)*
     &                 llp1_r
                  beta_y(i,ji,k)=beta_y(i,ji,k)+
     &              (phi_inv(i,l,2)+phi_exv(i,l,2))*p1_leg_c(ji,l)*
     &                 llp1_r
               end do
               beta_x(m+1,ji,k)=beta_x(m+1,ji,k)+
     &              phi_inv(m+1,l,1)*p1_leg_c(ji,l)*llp1_r
               beta_y(m+1,ji,k)=beta_y(m+1,ji,k)+
     &              phi_inv(m+1,l,2)*p1_leg_c(ji,l)*llp1_r
            end do
            
            do l=0,nleg,lstep
               do i=1,m
                  beta_z(i,ji,k)=beta_z(i,ji,k)+ 
     &                 (phi_inv(i,l,3)+phi_exv(i,l,3))*p_leg_c(ji,l)
                  chi   (i,ji,k)=chi   (i,ji,k)+
     &                 (phi_inv(i,l,4)+phi_exv(i,l,4))*p_leg_c(ji,l)
               end do
               beta_z(0,ji,k)=beta_z(0,ji,k)+ ! reflection symmetry
     &              phi_exv(1,l,3)*p_leg_c(n+1-ji,l)
               beta_z(m+1,ji,k)=beta_z(m+1,ji,k)+
     &              phi_inv(m+1,l,3)*p_leg_c(ji,l)
               chi   (m+1,ji,k)=chi   (m+1,ji,k)+
     &              phi_inv(m+1,l,4)*p_leg_c(ji,l)
            end do
            beta_x(0,ji,k)=-beta_x(1,ji,k) ! MINUS sign:
                                         ! beta_x=0 on the rotation axis
            beta_y(0,ji,k)=-beta_y(1,ji,k) ! MINUS sign:
                                         ! beta_y=0 on the rotation axis
            chi   (0,ji,k)= chi   (1,ji,k)

            beta_x(:,ji,k)=beta_x(:,ji,k)*prefac
            beta_y(:,ji,k)=beta_y(:,ji,k)*prefac
            beta_z(:,ji,k)=beta_z(:,ji,k)*prefac
            chi   (:,ji,k)=chi(:,ji,k)   *prefac
         end if

         if (j .ge. n_e-1) then

            ji=2*n_e+1-j

            beta_x(:,ji,k)=0.0_rk
            beta_y(:,ji,k)=0.0_rk
            beta_z(:,ji,k)=0.0_rk
            chi   (:,ji,k)=0.0_rk

            do l=1,nleg,lstep
               lp1=l+1
               llp1_r=1.0_rk/dble(lp1*l)
               do i=1,m
                  beta_x(i,ji,k)=beta_x(i,ji,k)+
     &                 (phi_inv(i,l,1)+phi_exv(i,l,1))*p1_leg_c(ji,l)*
     &                 llp1_r
                  beta_y(i,ji,k)=beta_y(i,ji,k)+
     &              (phi_inv(i,l,2)+phi_exv(i,l,2))*p1_leg_c(ji,l)*
     &                 llp1_r
               end do
               beta_x(m+1,ji,k)=beta_x(m+1,ji,k)+
     &              phi_inv(m+1,l,1)*p1_leg_c(ji,l)*llp1_r
               beta_y(m+1,ji,k)=beta_y(m+1,ji,k)+
     &              phi_inv(m+1,l,2)*p1_leg_c(ji,l)*llp1_r
            end do
            
            do l=0,nleg,lstep
               do i=1,m
                  beta_z(i,ji,k)=beta_z(i,ji,k)+ 
     &                 (phi_inv(i,l,3)+phi_exv(i,l,3))*p_leg_c(ji,l)
                  chi   (i,ji,k)=chi   (i,ji,k)+
     &                 (phi_inv(i,l,4)+phi_exv(i,l,4))*p_leg_c(ji,l)
               end do
               beta_z(0,ji,k)=beta_z(0,ji,k)+ ! reflection symmetry
     &              phi_exv(1,l,3)*p_leg_c(n+1-ji,l)
               beta_z(m+1,ji,k)=beta_z(m+1,ji,k)+
     &              phi_inv(m+1,l,3)*p_leg_c(ji,l)
               chi   (m+1,ji,k)=chi   (m+1,ji,k)+
     &              phi_inv(m+1,l,4)*p_leg_c(ji,l)
            end do
            beta_x(0,ji,k)=-beta_x(1,ji,k) ! MINUS sign:
                                         ! beta_x=0 on the rotation axis
            beta_y(0,ji,k)=-beta_y(1,ji,k) ! MINUS sign:
                                         ! beta_y=0 on the rotation axis
            chi   (0,ji,k)= chi   (1,ji,k)

            beta_x(:,ji,k)=beta_x(:,ji,k)*prefac
            beta_y(:,ji,k)=beta_y(:,ji,k)*prefac
            beta_z(:,ji,k)=beta_z(:,ji,k)*prefac
            chi   (:,ji,k)=chi(:,ji,k)   *prefac
         end if

      end do


c     ------------------------------------------------------------------
c     Transform back to spherical polar coordinate basis
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     treat ghost zones first
         if (j .eq. n_s+1) then
            ji=n_s-2
            do i=0,m+1
               beta(i,j,k,1)=   !use beta(:,j,:,:) for temporary storage
     &              sinus_theta(ji)*beta_x(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,j,k,2)=
     &              (cosin_theta(ji)*beta_x(i,ji,k)-
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,j,k,3)=
     &              beta_y(i,ji,k)*r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,j,k,1)*r(i) !Shibata method
            end do
         end if
         if (j .eq. n_e-1) then
            ji=n_e+2
            do i=0,m+1
               beta(i,j,k,1)=
     &              sinus_theta(ji)*beta_x(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,j,k,2)=
     &              (cosin_theta(ji)*beta_x(i,ji,k)-
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,j,k,3)=
     &              beta_y(i,ji,k)*r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,j,k,1)*r(i) !Shibata method
            end do
         end if
         if (j .eq. n_s) then
            ji=n_s-1
            do i=0,m+1
               beta(i,ji,k,1)=
     &              sinus_theta(ji)*beta_x(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,ji,k,2)=
     &              (cosin_theta(ji)*beta_x(i,ji,k)-
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,ji,k,3)=
     &              beta_y(i,ji,k)*r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,ji,k,1)*r(i) !Shibata method
            end do
         end if
         if (j .eq. n_e) then
            ji=n_e+1
            do i=0,m+1
               beta(i,ji,k,1)=
     &              sinus_theta(ji)*beta_x(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,ji,k,2)=
     &              (cosin_theta(ji)*beta_x(i,ji,k)-
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,ji,k,3)=
     &              beta_y(i,ji,k)*r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,ji,k,1)*r(i) !Shibata method
            end do
         end if

c     domain interior:
         do i=0,m+1
            beta(i,j,k,1)=
     &            sinus_theta(j)*beta_x(i,j,k)+
     &            cosin_theta(j)*beta_z(i,j,k)
#ifndef ENFORCE_SPHERICAL_METRIC
            beta(i,j,k,2)=
     &           (cosin_theta(j)*beta_x(i,j,k)-
     &            sinus_theta(j)*beta_z(i,j,k))*r_inv(i)
            beta(i,j,k,3)=
     &           beta_y(i,j,k)*r_inv(i)*sinus_theta_inv(j)
#endif
            chi (i,j,k)=chi(i,j,k)+beta(i,j,k,1)*r(i) !Shibata method
         end do

c     ------------------------------------------------------------------
c     Fill ghost zones (auxiliary potential chi only):
         chi(0,j,k)=chi(1,j,k)
         if (j .eq. 1) then
            chi(:,j-1,k)=chi(:,j,k)
         endif
         if (j .eq. n) then
            chi(:,j+1,k)=chi(:,j,k)
         end if

      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     ------------------------------------------------------------------
c     Shibata method: beta = W + grad (chi + r * W)

         do i=1,m
            beta(i,j,k,1)=beta(i,j,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(i+1,j,k)-aamp(i)*chi(i-1,j,k)+
     &           (aamp(i)-1.0_rk)*chi(i,j,k))
         end do
c         beta(0,j,k,1)=beta(1,j,k,1)-
c     &        (beta(2,j,k,1)-beta(1,j,k,1))*
c     &        2.0_rk*r(1)/(r(2)-r(1))
         beta(0,j,k,1)=beta(0,j,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(2,j,k)-aamp(1)*chi(1-1,j,k)+
     &           (aamp(1)-1.0_rk)*chi(1,j,k))
         beta(m+1,j,k,1)=beta(m+1,j,k,1)-0.125_rk*
     &        (chi(m+1,j,k)-chi(m,j,k))/delta_r_if(m)
#ifndef ENFORCE_SPHERICAL_METRIC
         do i=1,m+1
            beta(i,j,k,2)=beta(i,j,k,2)-0.0625_rk*
     &           delta_theta_inv*rr_inv(i)*
     &           (chi(i,j+1,k)-chi(i,j-1,k))
         end do
         beta(0,j,k,2)=beta(0,j,k,2)+0.0625_rk* !dchi/dtheta=0 on the rotation axis
     &        delta_theta_inv*rr_inv(i)*
     &        (chi(1,j+1,k)-chi(1,j-1,k))
#endif

c     ------------------------------------------------------------------
c     Fill ghost zones:
         if (j .eq. n_s) then
            ji=n_s-1
            do i=1,m
               beta(i,ji,k,1)=beta(i,ji,k,1)-0.125_rk*
     &              delta_r_amp_inv(i)*
     &              (chi(i+1,ji,k)-aamp(i)*chi(i-1,ji,k)+
     &              (aamp(i)-1.0_rk)*chi(i,ji,k))
            end do
            beta(0,ji,k,1)=beta(0,ji,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(2,ji,k)-aamp(1)*chi(1-1,ji,k)+
     &           (aamp(1)-1.0_rk)*chi(1,ji,k))
            beta(m+1,ji,k,1)=beta(m+1,ji,k,1)-0.125_rk*
     &           (chi(m+1,ji,k)-chi(m,ji,k))/delta_r_if(m)
#ifndef ENFORCE_SPHERICAL_METRIC
            do i=1,m+1
               beta(i,ji,k,2)=beta(i,ji,k,2)-0.0625_rk*
     &              delta_theta_inv*rr_inv(i)*
     &              (chi(i,ji+1,k)-chi(i,ji-1,k))
            end do
            beta(0,ji,k,2)=beta(0,ji,k,2)+0.0625_rk* !dchi/dtheta=0 on the rotation axis
     &           delta_theta_inv*rr_inv(i)*
     &           (chi(1,ji+1,k)-chi(1,ji-1,k))
#endif
         endif
         if (j .eq. n_e) then
            ji=n_e+1
            do i=1,m
               beta(i,ji,k,1)=beta(i,ji,k,1)-0.125_rk*
     &              delta_r_amp_inv(i)*
     &              (chi(i+1,ji,k)-aamp(i)*chi(i-1,ji,k)+
     &              (aamp(i)-1.0_rk)*chi(i,ji,k))
            end do
            beta(0,ji,k,1)=beta(0,ji,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(2,ji,k)-aamp(1)*chi(1-1,ji,k)+
     &           (aamp(1)-1.0_rk)*chi(1,ji,k))
            beta(m+1,ji,k,1)=beta(m+1,ji,k,1)-0.125_rk*
     &           (chi(m+1,ji,k)-chi(m,ji,k))/delta_r_if(m)
#ifndef ENFORCE_SPHERICAL_METRIC
            do i=1,m+1
               beta(i,ji,k,2)=beta(i,ji,k,2)-0.0625_rk*
     &              delta_theta_inv*rr_inv(i)*
     &              (chi(i,ji+1,k)-chi(i,ji-1,k))
            end do
            beta(0,ji,k,2)=beta(0,ji,k,2)+0.0625_rk* !dchi/dtheta=0 on the rotation axis
     &           delta_theta_inv*rr_inv(i)*
     &           (chi(1,ji+1,k)-chi(1,ji-1,k))
#endif
         end if

c     Preliminary solution for 3D: assume rotational symmetry of shift
c     vector field:
         if (k .eq. o_s) beta(:,j,k-1,:)=beta(:,j,k,:)
         if (k .eq. o_e) beta(:,j,k+1,:)=beta(:,j,k,:)

      end do

      return
      
      END SUBROUTINE loese_vektor_poisson

c     ==================================================================

#else


c     ==================================================================

      SUBROUTINE loese_vektor_poisson (s_vec,beta)

c     ==================================================================
c
c     Solves the generalized vector Poisson equation:
c     div grad beta + 1/3 grad div beta = S,
c     assuming a spherically symmetric vector field.
c
c     Remarks:
c     - Uses only the 0-th multipole of the source vector field
c     - Uses the Shibata scheme to treat the grad div term
c
c     s_vec: source at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c     beta:  solution at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c
c     ------------------------------------------------------------------

      USE precision

!      USE size_cfc
      USE grid_cfc
      USE legendre
      USE mod_hilfsgroessen_r
      USE parameters_cfc, ONLY: equatorial_symmetry
      USE perm_aux_cfc

      USE mod_poisson_aux, ONLY: beta_z,s_z,s_tmp,chi


#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      USE configure

      IMPLICIT NONE

      real (kind=rk), intent (in) :: s_vec (1:m,n_s:n_e,o_s:o_e,3)
      real (kind=rk), dimension(0:m+1,n_s-1:n_e+1,o_s-1:o_e+1,3),
     &     intent (out) :: beta

      integer (kind=ik) :: i

      integer (kind=ik) :: j,k,jk

      real (kind=rk) :: prefac

#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

      prefac=0.5_rk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i=1,m
            s_z(i,j,k)=s_vec(i,j,k,1)
            s_tmp(i,j,k)=-r(i)*s_vec(i,j,k,1)
         end do
      end do
         
c     ------------------------------------------------------------------
c     Expansion in Legendre polynomials P_l for the
c     z-component and auxiliary potential chi:

c     ------------------------------------------------------------------
c     Compute multipole moments: 

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      mom_s(:,0,:)=0.0_rk
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE

!$OMP DO REDUCTION (+:mom_s)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i=1,m
            mom_s(i,0,3)=mom_s(i,0,3)+s_z  (i,j,k)*p_leg_int(j,0)
            mom_s(i,0,4)=mom_s(i,0,4)+s_tmp(i,j,k)*p_leg_int(j,0)
         end do
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

#ifdef MPI_HYDRO
      call MPI_Allreduce (mom_s(:,0,3:4),mom_s2(:,0,3:4),2_ik*m,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,cart_comm,ierr)
      mom_s(:,0,3:4)=mom_s2(:,0,3:4)
#endif

      if (config%nsdim .eq. 3) then
         mom_s(1:m,0,3)=mom_s(1:m,0,3)/real(o,kind=rk)
         mom_s(1:m,0,4)=mom_s(1:m,0,4)/real(o,kind=rk)
      endif


c     ------------------------------------------------------------------
c     Compute inner/outer potentials: 
      phi_inv(1,0,3)=mom_s(1,0,3)*gew_r_in(1,1)*r_hoch_l_r(1,2)
      phi_inv(1,0,4)=mom_s(1,0,4)*gew_r_in(1,0)*r_hoch_l_r(1,1)
      do i=2,m
         phi_inv(i,0,3)=(phi_inv(i-1,0,3)*r_hoch_l(i-1,2)+
     &        mom_s(i,0,3)*gew_r_in(i,1))*
     &        r_hoch_l_r(i,2)
         phi_inv(i,0,4)=(phi_inv(i-1,0,4)*r_hoch_l(i-1,1)+
     &        mom_s(i,0,4)*gew_r_in(i,0))*
     &        r_hoch_l_r(i,1)
      end do
      phi_inv(m+1,0,3)=phi_inv(m,0,3)*
     &     (r_hoch_l(m,2)*r_hoch_l_r(m+1,2))
      phi_inv(m+1,0,4)=phi_inv(m,0,4)*
     &     (r_hoch_l(m,1)*r_hoch_l_r(m+1,1))

      phi_exv(m,0,3)=0.0_rk
      phi_exv(m,0,4)=0.0_rk
      do i=m-1,1,-1
         phi_exv(i,0,3)=(phi_exv(i+1,0,3)*r_hoch_l_r(i+1,1)+
     &        mom_s(i+1,0,3)*gew_r_ex(i+1,1))*
     &        r_hoch_l(i,1)
         phi_exv(i,0,4)=(phi_exv(i+1,0,4)*r_hoch_l_r(i+1,0)+
     &        mom_s(i+1,0,4)*gew_r_ex(i+1,0))*
     &        r_hoch_l(i,0)
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif


c     ------------------------------------------------------------------
c     Construct the solution vector field:

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         beta(:,j,k,2)=0.0_rk
         beta(:,j,k,3)=0.0_rk

c     ------------------------------------------------------------------
c     Sum over all multipoles: 
         
         do i=1,m
            beta_z(i,j,k)=
     &           (phi_inv(i,0,3)+phi_exv(i,0,3))*p_leg_c(j,0)
            chi   (i,j,k)=
     &           (phi_inv(i,0,4)+phi_exv(i,0,4))*p_leg_c(j,0)
         end do
         beta_z(0,j,k)=
     &        phi_exv(1,0,3)*p_leg_c(j,0)
         beta_z(m+1,j,k)=
     &        phi_inv(m+1,0,3)*p_leg_c(j,0)
         chi   (m+1,j,k)=
     &        phi_inv(m+1,0,4)*p_leg_c(j,0)
         chi   (0,j,k)= chi   (1,j,k)

c     ------------------------------------------------------------------
c     Add correct pre-factors
         beta_z(:,j,k)=beta_z(:,j,k)*prefac
         chi   (:,j,k)=chi(:,j,k)   *prefac


c     ------------------------------------------------------------------
c     Transform back to spherical polar coordinate basis
         do i=0,m+1
            beta(i,j,k,1)=beta_z(i,j,k)/3.0_rk
            chi (i,j,k)=chi(i,j,k)+beta(i,j,k,1)*r(i) !Shibata method
         end do

c     ------------------------------------------------------------------
c     Fill ghost zones (auxiliary potential chi only):
         chi(0,j,k)=chi(1,j,k)
         if (j .eq. n_s) chi(:,j-1,k)=chi(:,j,k)
         if (j .eq. n_e) chi(:,j+1,k)=chi(:,j,k)

#if 0
c     ------------------------------------------------------------------
c     Shibata method: beta = W + grad (chi + r * W)

         do i=1,m
            beta(i,j,k,1)=beta(i,j,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &        (chi(i+1,j,k)-aamp(i)*chi(i-1,j,k)+
     &           (aamp(i)-1.0_rk)*chi(i,j,k))
         end do
         beta(0,j,k,1)=beta(0,j,k,1)-0.125_rk*
     &        delta_r_amp_inv(i)*
     &        (chi(2,j,k)-aamp(1)*chi(1-1,j,k)+
     &        (aamp(1)-1.0_rk)*chi(1,j,k))
         beta(m+1,j,k,1)=beta(m+1,j,k,1)-0.125_rk*
     &        (chi(m+1,j,k)-chi(m,j,k))/delta_r_if(m)
#else
         beta(:,j,k,1)=beta(:,j,k,1)*0.75_rk
#endif

c     ------------------------------------------------------------------
c     Fill ghost zones:
         if (j .eq. n_s) beta(:,j-1,k,:)=beta(:,j,k,:)
         if (j .eq. n_e) beta(:,j+1,k,:)=beta(:,j,k,:)
         if (k .eq. o_s) beta(:,j,k-1,:)=beta(:,j,k,:)
         if (k .eq. o_e) beta(:,j,k+1,:)=beta(:,j,k,:)

      end do                    !j,k-loop

      return
      
      END SUBROUTINE loese_vektor_poisson

c     ==================================================================

#endif /* DIMN>1 && !defined(ENFORCE_SPHERICAL_METRIC) */

      END MODULE poisson

c     ==================================================================

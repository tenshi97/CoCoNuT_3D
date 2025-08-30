#undef SPLIT_LOOPS

module ppm

implicit none

contains


!> \verbatim
!> Calculate artificial viscosity terms. Includes HLLE fluxes
!> to eliminate "odd-even" decoupling for grid-aligned shocks.
!> Inconsistencies of older versions
!>
!>      o 0.5 were missing from the terms in tangential direction
!>      o above should be checked for 3-D
!>      o shock/ishck defined on zone centers while avis(i)
!>        is defined on zone boundaries. Fixed using maximum of 
!>        neighbouring zones.
!>
!> \endverbatim
!>
!> \author W. Keil, K. Kifonidis, M. Rampp
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
subroutine avisco (j, k)

  use precision
  use abort
  
  use intgrs_hy
      
  use mesh_hy
  use gfloat_hy
  use lfloat_hy
!  use intrp_hy ! forcheck
  use hydro_hy
  use vnew_hy
!  use physcs_hy ! forcheck
  use grd_hy
  use mo_mpi

  use specfun, only : fastsqrt

  use configure
  use state
  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik) :: i,metcav,n,k,j
  real(kind=rk) :: dytb,xtb3,xlrdth,dx3,rvm,rssi,rrs,rls
  real(kind=rk) :: chalf,a4,a1,am,bld,uhalf,betag
  real(kind=rk) :: br,bl,eta2du2,fcd1,fcd,brs,bls,bdifi
  real(kind=rk) ::  bdif,bm, bp,brd,u1r_rel,u1l_rel
  real(kind=rk) :: u1r,u1l,hdt,f4im,f3im,f1ip
  real(kind=rk) :: f5im,f2im,f4ip,f2ip,f1im,f5ip,f3ip
  real(kind=rk) :: sh_on,rsth,dth,dfi,dxtb3,rsthdth
  real(kind=rk) :: sthdfi,rsthdfi,dxtb,dzrl,dyrl
  real(kind=rk) :: sthdth,xbot2,xtop2,sinth,sinbot,xsq
  real(kind=rk) :: sintop,dxtbi,xtb2,xsqm1,dytbi,dx2
#ifdef GRAVITY
  real(kind=rk)  ::  rsxl(config%q)
#endif
  real(kind=rk) :: scratch(config%q), scratch2(config%q),       &
                  scratch3(config%q), scratch4(config%q),       &
                   scratch5(config%q), scratch6(config%q)

#ifdef HLLE_AT_SONIC_POINT
  real(kind=rk) :: mach_number, sound_speed
#endif

!     
!
!
  do i = 5,nzn5
     avis(i) = 0.0_rk
  end do

  if ( config%nsdim.eq.1 )   then

     if ( config%igeomx .eq. 0 ) then

        !           ---------
        !           1D planar


        do i = 5,nzn5
           avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))
           avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
        end do

     else if ( config%igeomx .eq. 1 ) then

        !           --------------
        !           1D cylindrical


        do i = 5,nzn5
           dx2 = x(i)*x(i) - x(i-1)*x(i-1)
           if ( dx2 .ne. 0.0_rk ) then
              avis(i) = (x(i) * u(i) - x(i-1) * u(i-1)) * 2._rk/ dx2
              avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
           end if
        end do

     else if ( config%igeomx .eq. 2 ) then

        !           ------------
        !           1D spherical


        do i = 5,nzn5
           
           xsq     = x(i)*x(i)
           xsqm1   = x(i-1)*x(i-1)
           avis(i) = (xsq * u(i) - xsqm1 * u(i-1))*3._rk /       &
                     (xsq * x(i) - xsqm1 * x(i-1)) 
           avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))

        end do

     else
        raise_abort("AVISCO(): 1D case not implemented!")

     end if
     
  end if

  if ( config%nsdim.eq.2 )  then
     
     if ( config%igeomx .eq. 0 .and. config%igeomy .eq. 0 ) then    

        !           ---------------
        !           2D planar (x,y)

        if ( xyzswp .eq. 1 ) then
           
           dytbi = 0.5_rk/(ytop - ybot)
           

           do i = 5,nzn5
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))           &
                         +(uttp(i)+uttp(i-1) -utbt(i)-utbt(i-1))    &
                         *dytbi
              avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
           end do

        else

           dxtbi = 0.5_rk/(xtop - xbot)


           do i = 5,nzn5
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))           &
                         +(uttp(i)+uttp(i-1) -utbt(i)-utbt(i-1))    &
                         *dxtbi
              avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
           end do

        end if

     else if ( config%igeomx .eq. 1 .and. config%igeomy .eq. 0 ) then 

        !           ------------------------
        !           2D cylindrical (r_cyl,z)

        if ( xyzswp.eq.1 ) then

           dytbi = 0.5_rk/(ytop - ybot)


           do i = 5,nzn5
                  dx2 = x(i)*x(i) - x(i-1)*x(i-1)
                  if ( dx2 .ne. 0.0_rk ) &
                       avis(i) = (x(i)*u(i) - x(i-1)*u(i-1)) * 2._rk / dx2
                  avis(i) = avis(i) &
                            +(uttp(i)+uttp(i-1)-utbt(i)-utbt(i-1)) &
                            *dytbi
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else

               if ( xtop*xtop .ne. xbot*xbot ) then
                  xtb2 = 2._rk/(2._rk*(xtop*xtop - xbot*xbot))
               else
                  xtb2 = -1.0_rk
               end if

               if ( xtb2 .gt. 0.0_rk ) then
                  
                  do i = 5,nzn5
                     avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))        &
                               +( xtop * (uttp(i) + uttp(i-1))          &
                               -xbot * (utbt(i) + utbt(i-1)) )*xtb2
                     avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
                  end do
               else

                  do i = 5,nzn5
                     avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))
                     avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
                  end do
               end if

            end if

         else if ( config%igeomx .eq. 1 .and. config%igeomy .eq. 3 ) then 

            !           --------------------------------
            !           2D cylindrical polar (r_cyl,phi)

            if ( xyzswp .eq. 1 ) then

               dytbi = 0.5_rk/(ytop-ybot)

               
               do i = 5,nzn5
                  dx2 = x(i)*x(i) - x(i-1)*x(i-1)
                  if ( dx2 .ne. 0.0_rk ) &
                     avis(i) = (x(i) * u(i) - x(i-1) * u(i-1)) * 2._rk/ dx2
                  if ( xl(i) .ne. 0.0_rk ) &
                       avis(i) = avis(i)+(uttp(i)+uttp(i-1) &
                                -utbt(i)-utbt(i-1))*dytbi/xl(i)
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else

               if ( xtop*xtop .ne. xbot*xbot ) then
                  xtb2 = 2._rk/(2._rk*(xtop*xtop - xbot*xbot))
               else
                  xtb2 = -1.0_rk
               end if

               if ( xtb2 .gt. 0.0_rk ) then

                  do i = 5,nzn5
                     avis(i) = (u(i) - u(i-1))/(xzn(j)*(x(i) - x(i-1)))  &
                               +( xtop * (uttp(i) + uttp(i-1))           &
                               -xbot * (utbt(i) + utbt(i-1)) )*xtb2
                     avis(i) = - config%cvisc * avis(i)* xzn(j) * ( x(i)-x(i-1) )
                  end do
               else

                  do i = 5,nzn5
                     avis(i) = (u(i) - u(i-1))/(xzn(j)*(x(i) - x(i-1)))
                     avis(i) = - config%cvisc * avis(i)*xzn(j) * ( x(i)-x(i-1) )
                  end do
               end if

            end if
            
         else if ( config%igeomx .eq. 2 .and. config%igeomy .eq. 4 ) then 
            
            !           --------------------------
            !           2D spherical (r_sph,theta)
            
            if ( xyzswp .eq. 1 ) then

               sintop = sin(ytop)
               sinbot = sin(ybot)
               sinth  = sin(yzn(j))
               sthdth = 2._rk*sinth*(ytop - ybot)


               do i = 5,nzn5
                  dx3 = x(i)**3 - x(i-1)**3
                  if ( dx3 .ne. 0.0_rk ) &
                       avis(i) = (x(i)**2 *u(i)- x(i-1)**2 *u(i-1)) &
                                 * 3._rk/dx3
                  xlrdth = xl(i)*sthdth
                  if ( xlrdth .ne. 0.0_rk ) &
                       avis(i) = avis(i)+(sintop * (uttp(i) + uttp(i-1)) - &
                                 sinbot * (utbt(i) + utbt(i-1))) /xlrdth
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else

               xtb3 = 3._rk/(2._rk*(xtop**3 - xbot**3))


               do i = 5,nzn5
                  sthdth = sin(xl(i))*(x(i)-x(i-1))
                  if ( sthdth .ne. 0.0_rk ) &
                       avis(i) = (sin(x(i)) *u(i) - sin(x(i-1)) *u(i-1)) / &
                                 (xzn(j) * sin(xl(i)) * (x(i) - x(i-1)) )
                  avis(i) = avis(i) +(xtop**2 * (uttp(i) + uttp(i-1)) - &
                            xbot**2 * (utbt(i) + utbt(i-1))) * xtb3
                  avis(i) = - config%cvisc * avis(i) * xzn(j) &
                            * (x(i) - x(i-1))
               end do

            end if

         else if ( config%igeomx .eq. 2 .and. config%igeomy .eq. 5 ) then  

            !           ----------------------------------------------
            !           2D spherical polar (r_sph,phi) ;  theta = pi/2

            if ( xyzswp.eq.1 ) then

               dytb = 2._rk*(ytop - ybot)


               do i = 5,nzn5
                  avis(i) = (x(i)**2 *u(i) - x(i-1)**2 *u(i-1)) * &
                             3._rk / (x(i)**3 - x(i-1)**3)
                  if ( xl(i) .ne. 0.0_rk ) &
                       avis(i) = avis(i)+(uttp(i)+uttp(i-1) &
                                 -utbt(i)-utbt(i-1)) / (xl(i) * dytb)
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else

               xtop2 = xtop*xtop
               xbot2 = xbot*xbot
               xtb3 = 3._rk/(2._rk*(xtop2*xtop - xbot2*xbot))


               do i = 5,nzn5
                  if ( xzn(j) .ne. 0.0_rk ) &
                       avis(i) = (u(i) - u(i-1)) / (xzn(j)*(x(i) - x(i-1)))+ &
                                 ( xtop2 * (uttp(i) + uttp(i-1)) - &
                                   xbot2 * (utbt(i) + utbt(i-1)) ) * xtb3
                  avis(i) = - config%cvisc * avis(i) * xzn(j) * (x(i) - x(i-1))
               end do

            end if

         else 

            raise_abort("AVISCO(): 2D case not implemented!")

         end if

      end if
      
      if ( config%nsdim .eq. 3 ) then
         
         if ( config%igeomx .eq. 0 .and. config%igeomy .eq. 0   &
             .and. config%igeomz .eq. 0 ) then 
            
            !           --------------------
            !           3D cartesian (x,y,z)

            if ( xyzswp.eq.1 ) then

               dytb = 0.5_rk/(ytop - ybot)
               dzrl = 0.5_rk/(zrgt - zlft)


               do i = 5,nzn5
                  avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +        &
                            (uttp(i)+uttp(i-1) -utbt(i)-utbt(i-1))*dytb+&       
                            (utrt(i)+utrt(i-1) -utlt(i)-utlt(i-1))*dzrl
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else if ( xyzswp .eq. 2 ) then

               dxtb = 0.5_rk/(xtop - xbot)
               dzrl = 0.5_rk/(zrgt - zlft)


               do i = 5,nzn5
                  avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +  &      
                            (uttp(i)+uttp(i-1) -utbt(i)-utbt(i-1))*dxtb+ &      
                            (utrt(i)+utrt(i-1) -utlt(i)-utlt(i-1))*dzrl
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else

               dxtb = 0.5_rk/(xtop - xbot)
               dyrl = 0.5_rk/(yrgt - ylft)


               do i = 5,nzn5
                  avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +        &
                            (uttp(i)+uttp(i-1) -utbt(i)-utbt(i-1))*dxtb+&       
                            (utrt(i)+utrt(i-1) -utlt(i)-utlt(i-1))*dyrl
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            end if

         else if ( config%igeomx .eq. 2 .and. config%igeomy .eq. 4    &
                   .and. config%igeomz .eq. 5) then

            !           ------------------------------
            !           3D spherical (r_sph,theta,phi)

            if ( xyzswp .eq. 1 ) then

               sintop = sin(ytop)
               sinbot = sin(ybot)
               sinth  = sin(yzn(j))
               sthdth = 2._rk*sinth*(ytop - ybot)
               sthdfi = 2._rk*sinth*(zrgt - zlft)
               

               do i = 5,nzn5
                  avis(i) = (x(i)**2 *u(i) - x(i-1)**2 *u(i-1)) *  3._rk / &
                            (x(i)**3 - x(i-1)**3)
                  rsthdth = xl(i)*sthdth
                  if ( rsthdth .ne. 0.0_rk )            &
                       avis(i) = avis(i)+( sintop * (uttp(i) + uttp(i-1)) &
                                 -sinbot * (utbt(i) + utbt(i-1)) )/rsthdth
                  rsthdfi = xl(i)*sthdfi
                  if ( rsthdfi .ne. 0.0_rk ) &
                               avis(i) = avis(i)+(utrt(i)+utrt(i-1) &
                                        -utlt(i)-utlt(i-1))/rsthdfi
                  avis(i) = - config%cvisc * avis(i) * (x(i) - x(i-1))
               end do

            else if ( xyzswp .eq. 2 ) then

               xtop2 = xtop*xtop
               xbot2 = xbot*xbot
               dxtb3 = 2._rk*(xtop2*xtop - xbot2*xbot)
               dfi   = 2._rk*(zrgt - zlft)


               do i = 5,nzn5
                  sinth   = sin( xl(i) )
                  if ( sinth .ne. 0.0_rk ) &
                       avis(i) = (sin(x(i)) * u(i) - sin(x(i-1)) * u(i-1)) / &
                                 (xzn(j) * sinth * (x(i) - x(i-1)) )
                  avis(i) = avis(i)+( xtop2 * (uttp(i) + uttp(i-1)) - &
                            xbot2 * (utbt(i) + utbt(i-1)) ) /dxtb3
                  rsthdfi = xzn(j) * sinth * dfi
                  if ( rsthdfi .ne. 0.0_rk ) &
                       avis(i) = avis(i)+(utrt(i)+utrt(i-1) &
                                -utlt(i)-utlt(i-1))/(2._rk*rsthdfi)
                  avis(i) = - config%cvisc * avis(i) * xzn(j) &
                             * (x(i) - x(i-1))
               end do

            else

               xtop2  = xtop*xtop
               xbot2  = xbot*xbot
               dxtb3  = 2._rk*(xtop2*xtop - xbot2*xbot)
               sintop = sin(yrgt)
               sinbot = sin(ylft)
               sinth  = sin(yzn(k))
               dth    = 2._rk*(yrgt - ylft)


               do i = 5,nzn5
                  rsth = xzn(j) * sinth
                  if ( rsth .ne. 0.0_rk ) &
                       avis(i) = ( u(i) - u(i-1))/((x(i) - x(i-1))*rsth ) &
                                 +( sintop * (utrt(i) + utrt(i-1))-sinbot &
                                 * (utlt(i) + utlt(i-1)) )/( rsth*dth )  
                  avis(i) = avis(i)+( xtop2 * &
                           (uttp(i) + uttp(i-1))-xbot2 * &
                           (utbt(i) + utbt(i-1)) )/dxtb3
                  avis(i) = - config%cvisc * avis(i) * xzn(j) &
                            * sinth * (x(i) - x(i-1))
               end do

            end if
            
         else
            raise_abort("AVISCO(): 3D case not implemented!")

         end if

      end if



!     -------------------------------------
!     no viscosity if zone is outside shock
!

      if ( config%nsdim .eq. 1 ) then 
         !        ----------
         !        1D version
         do i = 1,nzn+1
            sh_on     = max(ishck(i-1,j,k),ishck(i,j,k))
            !            ishock_pos(i,j,k)=ishck(i,j,k)
            avis(i+4) = sh_on*max(0.0_rk,avis(i+4))
         end do         
         
      else if ( config%nsdim .eq. 2 ) then
         
         !     ----------------------------------------------------------
         !     2D version. Shock needs to be "seen" also in 2nd dimension 
         !     in order to choose HLLE solver also for y-sweeps to avoid 
         !     "odd-even" decoupling.
                 

         do i = 1,nzn+1
            if ( xyzswp .eq. 1 ) then
               sh_on = max(ishck(i-1,j,k),ishck(i,j,k) )
!            ishock_pos(i,j,k)=ishck(i,j,k)
            else 
               sh_on = max(ishck(j,i+qy_s-2,k),ishck(j,i+qy_s-1,k) )
            end if
            avis(i+4) = sh_on*max(0.0_rk,avis(i+4))
!               ishock_pos(j,i,k)=ishck(j,i,k)
         end do


      else if ( config%nsdim .eq. 3 ) then

         do i = 1,nzn+1
            if ( xyzswp .eq. 1 ) then
               sh_on = max(ishck(i-1,j,k),ishck(i,j,k) )
!               ishock_pos(i,j,k)=ishck(i,j,k)
            else if ( xyzswp .eq. 2 ) then

               sh_on = max(ishck(j,i+qy_s-2,k),ishck(j,i+qy_s-1,k) )
!               ishock_pos(j,i,k)=ishck(j,i,k)

            else if ( xyzswp .eq. 3 ) then
               sh_on = max(ishck(j,k,i+qz_s-2),ishck(j,k,i+qz_s-1) )
!               ishock_pos(j,k,i)=ishck(j,k,i)
            end if
            avis(i+4) = sh_on*max(0.0_rk,avis(i+4))
         end do

      end if

#ifdef HLLE_AT_SONIC_POINT
      do i = 1,nzn+1
         sound_speed = 0.5_rk * (ce(i + 3) + ce(i + 4))
         mach_number = 0.5_rk * abs(u(i + 3) + u(i + 4)) / sound_speed
         if (mach_number > 0.7_rk .and. mach_number < 1.3_rk) then
!           write(*,'(a,i0,3(a,1pe12.5))') &
!             "sh_on = 1 at i = ", i, ", x = ", 0.5_rk * (x(i+4) + x(i+3)), &
!             ", mach_number = ", mach_number, ", sound_speed = ", sound_speed
            sh_on = 1_ik
         endif
      end do
#endif

      
!     ----------------
!     add viscous flux
#ifdef USE_HLLE
      scratch(:)  = 0._rk
      scratch2(:) = 0._rk    
      scratch3(:) = 0._rk    
      scratch4(:) = 0._rk    
      scratch5(:) = 0._rk    
      scratch6(:) = 0._rk  

      do i = 5,nzn5
         scratch3(i-4) = rho(i-1)
         scratch4(i-4) = rho(i)
      enddo

      call fastsqrt(scratch3, (nzn5-5+1) )
      call fastsqrt(scratch4, (nzn5-5+1) )

      do i = 5,nzn5
         scratch (i) = scratch3(i-4)
         scratch2(i) = scratch4(i-4)
      enddo

      scratch3(:) = 0._rk    
      scratch4(:) = 0._rk   

      do i = 5,nzn5

         if ( avis(i) .ne. 0.0_rk ) then

!           --------------------------------------
!>           need HLLE (Harten-Lax-Van Leer) fluxes

!>           Einfeldt, B., Munz, C.D., Roe, P.L., Sjoegren, B.,
!>                                               1991, JCP, 92, 273
!>           Einfeldt, B., 1988, SIAM JNA, 25, 294
!>           Roe, P.L., 1981, JCP, 43, 357
!>           Einfeldt, B., 1988, in "Shock Tubes and Waves",
!>           ed. H. Groenig, Proceedings of the Sixteenth Int.Symp.
!>           on Shock Tubes and Waves, Aachen July 26-31, 1987, p. 671
!>           (equation numbers refer to the first referece)

!>           metcav
!>           0   most dissipative
!>           1   standard
!>           2   sharpest
!>           3   medium

            metcav = 2

!           ---------------------------
!           calculate Roe's eigenvalues

!           --------------
!           density (5.3a)

            rls  = scratch (i)
            rrs  = scratch2(i)

            rssi = 1.0_rk/(rls+rrs)

!           ------------------------
!           velocities (5.3b & 5.3c)

 

!           ------------------
!           sound speed (5.3e)
!           Einfeldt (5.7)

            eta2du2 = 0.5_rk*rls*rrs*rssi**2 * (u(i) - u(i-1))**2
            scratch5(i-4) = (rls*ce(i-1)**2 +rrs*ce(i)**2)*rssi+eta2du2

            betag = 0.5_rk*( game(i-1) + game(i) )
            scratch6(i-4) = ( 0.5_rk*(betag-1.0_rk)/betag )

         endif ! avis ne 0
         
      enddo

      call fastsqrt(scratch5, (nzn5-5+1))
      call fastsqrt(scratch6, (nzn5-5+1))

      do i = 5,nzn5
         scratch3(i)=scratch5(i-4)
         scratch4(i)=scratch6(i-4)
      enddo

      scratch5(:) =0._rk
      scratch6(:) =0._rk


      do i = 5,nzn5
         if ( avis(i) .ne. 0.0_rk ) then
            rls  = scratch (i)
            rrs  = scratch2(i)
            rssi = 1.0_rk/(rls+rrs)
            
            rvm = ( rls*u(i-1) + rrs*u(i) )*rssi
            am = scratch3(i)
         
            !           -----------------
            !           Roe's eigenvalues

            a1 = rvm - am
            a4 = rvm + am

            !           ---------------------------
            !           numerical signal velocities
            
            !           ----------------------
            !           most dissipative (4.7)
            
            if ( metcav .eq. 0 .or. metcav .ge. 3 ) then
            
               !              (4.7a), (4.7b)
               
               chalf = 0.5_rk*(a4-a1)
               uhalf = 0.5_rk*(a4+a1) - ugrdl(i)
               
               brd = max( abs(uhalf)             + chalf, &
                          abs(u(i-1)-ugrdl(i)) + ce(i-1), &
                          abs(u(i  )-ugrdl(i)) + ce(i)  )
               bld = -brd

            end if


            !           --------------
            !           standard (4.5)

            if ( metcav .eq. 1 .or. metcav .ge. 3 ) then

            !              (4.5a), (4.5b)
               
               bl = min( a1 , u(i-1) - ce(i-1) ) - ugrdl(i)
               br = max( a4 , u(i  ) + ce(i  ) ) - ugrdl(i)

            end if

!           ---------------
!           sharpest (4.10)

            if ( metcav .eq. 2 .or. metcav .ge. 3 ) then

!              (4.9b), (4.10a), (4.10b)

!               betag = 0.5_rk*( game(i-1) + game(i) )
!               betag = sqrt( 0.5_rk*(betag-1.0_rk)/betag )

               betag = scratch4(i)
            
               bls = min(a1,u(i-1) - betag*ce(i-1)) - ugrdl(i)
               brs = max(a4,u(i  ) + betag*ce(i  )) - ugrdl(i)

            end if


            if ( metcav .eq. 0 ) then
               bl = bld
               br = brd
            else if ( metcav .eq. 2 ) then
               bl = bls
               br = brs
            else if ( metcav .eq. 3 ) then
               fcd  = 0.5_rk
               fcd1 = 1.0_rk-fcd
               bl   = fcd1*bls + fcd*bld
               br   = fcd1*brs + fcd*brd
            end if

!           ---------------------
!           see comment to (4.4b)

            bp = max( br,0.0_rk )
            bm = min( bl,0.0_rk )
            bdif = bp - bm

!           -------
!           failure

            if ( bdif .eq. 0.0_rk ) then

               write(*,*) "Task ",myproc
               write(*,*) "Task ",myproc,' HLLE failed'
               write(*,*) "Task ",myproc,' rvm     = ',rvm
               write(*,*) "Task ",myproc,' am      = ',am
               write(*,*) "Task ",myproc,' a1      = ',a1
               write(*,*) "Task ",myproc,' a4      = ',a4
               write(*,*) "Task ",myproc,' (u-c)_L = ',u(i-1) - ce(i-1)
               write(*,*) "Task ",myproc,' (u-c)_R = ',u(i  ) - ce(i  )
               write(*,*) "Task ",myproc,' rho_L   = ',rho(i-1)
               write(*,*) "Task ",myproc,' rho_R   = ',rho(i)
               write(*,*) "Task ",myproc,' gamc_L  = ',gamc(i-1)
               write(*,*) "Task ",myproc,' gamc_R  = ',gamc(i)
               write(*,*) "Task ",myproc,' ugrd    = ',ugrdl(i)
               write(*,*) "Task ",myproc,' bl      = ',bl
               write(*,*) "Task ",myproc,' br      = ',br
               write(*,*) "Task ",myproc,' bm      = ',bm
               write(*,*) "Task ",myproc,' bp      = ',bp
               raise_abort("avisco(): HLLE failed")
            end if
            
            bdifi = 1.0_rk/bdif
            

#ifdef GRAVITY
!           ----------------------------------
!           include gravitational acceleration

            hdt = 0.5_rk * hydro%dt
            
            u1l = u(i-1) + hdt * gravr(i-1)
            u1r = u(i)   + hdt * gravl(i)
#else
            u1l = u(i-1)
            u1r = u(i)
#endif

!           -----------------------
!           account for grid motion 

            u1l_rel = u1l - ugrdl(i)
            u1r_rel = u1r - ugrdl(i) 


!           --------------
!           f fluxes (1.2)

!           -----
!           i-1/2

            f1im        =       rho(i-1)*u1l_rel
            f2im        =       f1im*u1l
            f3im        =       f1im*ut(i-1)
            f4im        =       f1im*utt(i-1)
            f5im        =       f1im*e(i-1) + u(i-1)*p(i-1)

!           -----
!           i+1/2

            f1ip        =       rho(i)*u1r_rel
            f2ip        =       f1ip*u1r
            f3ip        =       f1ip*ut(i)
            f4ip        =       f1ip*utt(i)
            f5ip        =       f1ip*e(i) + u(i)*p(i)

!           ------------------
!           HLLE fluxes (4.4b)

!           ----
!           mass

            rhoflx(i) = ( bp*f1im-bm*f1ip +bp*bm*(rho(i)-rho(i-1)) )*bdifi

!           -------
!           momenta

            uflx(i)   = ( bp*f2im-bm*f2ip +bp*bm*(rho(i)*u(i)-rho(i-1)* &
                          u(i-1)) )*bdifi

            utflx(i)  = ( bp*f3im-bm*f3ip &
                         +bp*bm*(rho(i)*ut(i)-rho(i-1)*ut(i-1)) &
                         )*bdifi

            uttflx(i) = ( bp*f4im-bm*f4ip &
                         +bp*bm*(rho(i)*utt(i)-rho(i-1)*utt(i-1)) &
                         )*bdifi

!           ------------
!           total energy

            eflx(i)   = ( bp*f5im-bm*f5ip &
                         +bp*bm*(rho(i)*e(i)-rho(i-1)*e(i-1)) &
                         )*bdifi

!           -------
!           species

            do n = 1, config%qn
               xnflx(i,n) = xnav(i,n) * rhoflx(i)
            end do
            
         end if ! avis != 0
         
      end do
#endif /* use_HLLE */

        if ( bndmin .eq. 1 ) avis(5)    = 0.0_rk
        if ( bndmax .eq. 1 ) avis(nzn5) = 0.0_rk

        if ( igeom .eq. 4 .and. bndmin .eq. 4 ) then
           avis(5)    = 0.0_rk
           avis(nzn5) = 0.0_rk
        end if


      do i = 5,nzn5
         rhoflx(i) = rhoflx(i) + avis(i) * (rho(i-1) - rho(i))
         uflx  (i) = uflx  (i) + avis(i) * &
                     (rho(i-1) * u(i-1) - rho(i) * u  (i))
         utflx (i) = utflx (i) + avis(i) * &
                     (rho(i-1) * ut(i-1) - rho(i) * ut (i))
         uttflx(i) = uttflx(i) + avis(i) *  &
                     (rho(i-1) * utt(i-1) - rho(i) * utt(i))
         eflx  (i) = eflx  (i) + avis(i) * &
                     (rho(i-1) * e(i-1) - rho(i) * e(i))
      end do

      do n = 1,config%qn

         do i = 5,nzn5
            xnflx(i,n) = xnflx(i,n)  +  avis(i) * &
                         (rho(i-1) * xn(i-1,n) - rho(i) * xn(i,n))
         end do
      end do

!
! return
!
      return
!
! end of AVISCO
!
    end subroutine avisco


!> \verbatim
!> interpolate interface values (2nd order)
!> Korevaar, P., van Leer, B., 1988, A&A, 200, 153
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param al
!> \param a
!> \param ar
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
    subroutine interpl ( al, a, ar )

      use precision

      use intgrs_hy

      use intrp_hy

      use configure
      use state
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i
      real(kind=rk) :: dabias
      real(kind=rk) :: al(config%q), a(config%q), ar(config%q), scrch1(config%q)
!
!
!

!     --------------
!     z-z difference

      do i = 2,nzn8
         scrch1(i) = a(i) - a(i-1)
      end do

!     -----
!     slope

      do i = 2,nzn7
         dabias  =0.001_rk*(abs(scrch1(i))+abs(scrch1(i+1))+0.01_rk)**2

         dela(i) = 0.5_rk*( scrch1(i+1) + scrch1(i) ) &
                   *( 2.0_rk*scrch1(i+1)*scrch1(i)  + dabias ) &
                  /( scrch1(i+1)**2+scrch1(i)**2 + dabias )

!        --------------
!        slope limiting

         dela(i) = sign(  min( abs(dela(i)),                                 &
                               2.0_rk*min( abs(scrch1(i+1)),abs(scrch1(i)) ) ), &
                               dela(i) )

         if ( scrch1(i+1)*scrch1(i).le.0.0_rk ) dela(i) = 0.0_rk

         al(i) = a(i) - 0.5_rk*dela(i)
         ar(i) = a(i) + 0.5_rk*dela(i)
      end do
!
! return
!
      return
!
! end of INTERPL
!
    end subroutine interpl
!> \verbatim
!> compute left and right states for input to riemann problem
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
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
    subroutine states ( j )

      use precision

      use bndinf_hy
      use intgrs_hy
     
      use mesh_hy
      use gfloat_hy
!      use lfloat_hy ! forcheck
      use intrp_hy
      use hydro_hy
      use physcs_hy
      use grd_hy
      use reman_hy
      use specfun, only : fastsqrt

      use configure
      use state
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i,n,n0,j
      real(kind=rk) :: beta0r,betamr,betapr,rhortdi,xcrght
      real(kind=rk) :: geosrc,eta,courno,fdt,betaml,betapl
      real(kind=rk) :: xclft,rholtdi,scr2,beta0l,utmp,qdt
      real(kind=rk) :: hdt,scr1,tanx
      logical       :: lepsth
      real(kind=rk),parameter ::cmone = -1._rk

      real(kind=rk) ::slamp  (config%q), slamm  (config%q), slam0 (config%q)
      real(kind=rk) ::pltd   (config%q), ultd   (config%q), rholtd(config%q), utltd (config%q), &
                      uttltd(config%q), gameltd(config%q), gamcltd(config%q), gamcpl(config%q), &
                      gamcml(config%q),gravpl (config%q), gravml (config%q), dloga(config%q)

      real(kind=rk) ::scrch1(config%q), scrch2(config%q)
!
!
!
      hdt = 0.5_rk * hydro%dt
      qdt = 0.25_rk * hydro%dt


      do i = 4,nzn6

!        -------------------------------------
!        for moving grid use relative velocity

         utmp      = u(i)
         u(i)      = u(i) - ugrid(i)
         urel(i)   = utmp

!        -------------------------
!        characteristic velocities

         slamp(i)  = dtdx(i)*(u(i)+ce(i))
         slamm(i)  = dtdx(i)*(u(i)-ce(i))
         slam0(i)  = dtdx(i)*u(i)
      end do

!     --------------------------
!     gravitational acceleration

#ifdef GRAVITY
         scrch1(1) = 0.0_rk

         if ( igeom .eq. 4 ) then

            do i = 2,nzn8
               scrch1(i) = - (grav(i) - grav(i-1))/( xzn(j)*dx(i) )
            end do
         else

            do i = 2,nzn8
               scrch1(i) = - (grav(i) - grav(i-1))/dx(i)
            end do
         end if
#else
         do i = 1,nzn8
            scrch1(i) = 0.0_rk
         end do
#endif


!     ---------------------
!     add fictitious forces


        do i = 1,nzn8
           grav(i) = scrch1(i) + fict(i)
        end do

#ifdef NFORCE
        do i = 1,nzn8
           grav(i) = grav(i) + s(i)
        end do
#endif

        call interpl ( gravl, grav, gravr )


        do i = 4,nzn5
           dgrav(i) = gravr(i) - gravl(i)
           grav6(i) = 0.0_rk
        end do

!     -----------
!     LEFT STATES
!     -----------


      do i = 5,nzn5
         scr1       = 0.5_rk*max( 0.0_rk,slamp(i-1) )
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         pltd(i)    = pr(i-1)    - scr1*(dp(i-1)    - scr2 * p6(i-1))
         ultd(i)    = ur(i-1)    - scr1*(du(i-1)    - scr2 * u6(i-1))
         rholtd(i)  = rhor(i-1)  - scr1*(drho(i-1)  - scr2 * rho6(i-1))
         utltd(i)   = utr(i-1)   - scr1*(dut(i-1)   - scr2 * ut6(i-1))
         uttltd(i)  = uttr(i-1)  - scr1*(dutt(i-1)  - scr2 * utt6(i-1))
         gameltd(i) = gamer(i-1) - scr1*(dgame(i-1) - scr2 * game6(i-1))
         gamcltd(i) = gamcr(i-1) - scr1*(dgamc(i-1) - scr2 * gamc6(i-1))
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scr1       = 0.5_rk*slamp(i-1)
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         ppl(i)     = pr(i-1)    - scr1*(dp(i-1)    - scr2 * p6(i-1))
         upl(i)     = ur(i-1)    - scr1*(du(i-1)    - scr2 * u6(i-1))
         rhopl(i)   = rhor(i-1)  - scr1*(drho(i-1)  - scr2 * rho6(i-1))
         gamcpl(i)  = gamcr(i-1) - scr1*(dgamc(i-1) - scr2 * gamc6(i-1))

         gravpl(i)  = gravr(i-1) - scr1*(dgrav(i-1) - scr2 * grav6(i-1))

         ppl(i)     = max( config%smallp , ppl(i) )
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scr1       = 0.5_rk*slamm(i-1)
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         pml(i)     = pr(i-1)    - scr1*(dp(i-1)    - scr2 * p6(i-1))
         uml(i)     = ur(i-1)    - scr1*(du(i-1)    - scr2 * u6(i-1))
         rhoml(i)   = rhor(i-1)  - scr1*(drho(i-1)  - scr2 * rho6(i-1))
         gamcml(i)  = gamcr(i-1) - scr1*(dgamc(i-1) - scr2 * gamc6(i-1))

         gravml(i)  = gravr(i-1) - scr1*(dgrav(i-1) - scr2 * grav6(i-1))

         pml(i)     = max( config%smallp , pml(i) )
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scrch1(i)  = 0.5_rk*slam0(i-1)
         scrch2(i)  = 1.0_rk - 4.0_rk/3.0_rk*scrch1(i)

         p0l(i)     = pr(i-1) &
                      -scrch1(i) * (dp(i-1)    - scrch2(i) * p6(i-1))
         rho0l(i)   = rhor(i-1) &
                     -scrch1(i) * (drho(i-1)  - scrch2(i) * rho6(i-1))
         ut0l(i)    = utr(i-1) &
                     -scrch1(i) * (dut(i-1)   - scrch2(i) * ut6(i-1))
         utt0l(i)   = uttr(i-1) &
                     -scrch1(i) * (dutt(i-1)  - scrch2(i) * utt6(i-1))
         game0l(i)  = gamer(i-1) &
                     -scrch1(i) * (dgame(i-1) - scrch2(i) * game6(i-1))
         gamc0l(i)  = gamcr(i-1) &
                     -scrch1(i) * (dgamc(i-1) - scrch2(i) * gamc6(i-1))

         p0l(i)     = max( p0l(i) , config%smallp )
      end do

!     -----------------------------------------
!     left states for fluids (unrolled 5 times)

      n0 = 1
      do n = 1,config%qn-4,5
         do i = 5,nzn5
            xn0l(i,n  ) = xnr(i-1,n  ) &
                 -scrch1(i) * (dxn(i-1,n  ) - scrch2(i) * xn6(i-1,n  ))
            xn0l(i,n+1) = xnr(i-1,n+1) &
                 -scrch1(i) * (dxn(i-1,n+1) - scrch2(i) * xn6(i-1,n+1))
            xn0l(i,n+2) = xnr(i-1,n+2) &
                 -scrch1(i) * (dxn(i-1,n+2) - scrch2(i) * xn6(i-1,n+2))
            xn0l(i,n+3) = xnr(i-1,n+3) &
                 -scrch1(i) * (dxn(i-1,n+3) - scrch2(i) * xn6(i-1,n+3))
            xn0l(i,n+4) = xnr(i-1,n+4) &
                 -scrch1(i) * (dxn(i-1,n+4) - scrch2(i) * xn6(i-1,n+4))
         end do
         n0 = n+5
      end do

      do n = n0,config%qn
         do i = 5,nzn5
            xn0l(i,n  ) = xnr(i-1,n  ) &
                 -scrch1(i) * (dxn(i-1,n  ) - scrch2(i) * xn6(i-1,n  ))
         end do
      end do


!     ---------------------
!     effective left states

      do i = 5,nzn5
         clft(i) = gamcltd(i) * pltd(i) * rholtd(i)
      enddo

      call fastsqrt(clft, nzn5)

      do i = 5,nzn5

         xclft     = 1.0_rk / clft(i)
         rholtdi   = 1.0_rk / rholtd(i)

         betapl    = -0.5_rk*( (ultd(i) - upl(i)) &
                            +(pltd(i) - ppl(i)) * xclft &
                           )
         betaml    = +0.5_rk*( (ultd(i) - uml(i)) &
                            -(pltd(i) - pml(i)) * xclft &
                           )
         beta0l    = (pltd(i) - p0l(i)) * xclft * xclft &
                    +rholtdi -1.0_rk/rho0l(i)

!        ------------------------------
!        fix for wrong formulae in CW84

         if ( slamp(i-1).le.0.0_rk ) then
            betapl    = 0.0_rk
            gravpl(i) = 0.0_rk
         end if
         if ( slamm(i-1).le.0.0_rk ) then 
            betaml    = 0.0_rk
            gravml(i) = 0.0_rk
         end if
         if ( slam0(i-1).le.0.0_rk ) beta0l = 0.0_rk

         if ( slamm(i-1) .gt. 0.0_rk .and. slamp(i-1) .gt. 0.0_rk ) then
            fdt = qdt
         else
            fdt = hdt
         end if

         vlft(i) = rholtdi - beta0l - ( betapl + betaml )*xclft
         plft(i) = pltd(i) + clft(i)*( betapl + betaml )
         plft(i) = max( plft(i) , config%smallp )
         ulft(i) = ultd(i) + ( betapl - betaml ) &
                   +fdt*(gravml(i)+gravpl(i))

         if ( slam0(i-1).le.0.0_rk ) then
            utlft(i)  = utltd(i)
            uttlft(i) = uttltd(i)
            gmelft(i) = gameltd(i)
            gmclft(i) = gamcltd(i)
         else
            utlft(i)  = ut0l(i)
            uttlft(i) = utt0l(i)
            gmelft(i) = game0l(i)
            gmclft(i) = gamc0l(i)
         end if
      end do

!     ---------------------------------------------------
!     effective left states for fluids (unrolled 5 times)

      n0 = 1
      do n = 1,config%qn-4,5
         do i = 5,nzn5
            xnlft(i,n  ) = xn0l(i,n  )
            xnlft(i,n+1) = xn0l(i,n+1)
            xnlft(i,n+2) = xn0l(i,n+2)
            xnlft(i,n+3) = xn0l(i,n+3)
            xnlft(i,n+4) = xn0l(i,n+4)
         end do
         n0 = n+5
      end do

      do n = n0,config%qn
         do i = 5,nzn5
            xnlft(i,n  ) = xn0l(i,n  )
         end do
      end do


!     -----------------------
!     geometrical source term

      do i = 4,nzn5
         dloga(i) = 0.0_rk
      end do

      if ( igeom.le.2 ) then

         do i = 4,nzn5
            if ( x(i) .ne. 0.0_rk ) dloga(i) = igeom/x(i)
         end do

      else if ( igeom .eq. 4 ) then

         do i = 4,nzn5
            tanx = tan(x(i))
            if ( tanx .ne. 0.0_rk ) dloga(i) = 1.0_rk/(xzn(j)*tanx)
         end do

      end if

!     ---------------------------------
!     Colella hack near R=0 singularity

      do i = 4,nzn5
         if ( dloga(i) .ne. 0.0_rk ) then
            courno   = (abs(u(i))+ce(i)) * dtdx(i)
            eta      = ( 1.0_rk - courno )/( ce(i)*hydro%dt*abs(dloga(i)) )
            eta      = min(1.0_rk , abs(eta))
            dloga(i) = eta * dloga(i)
         end if
      end do

!     ---------------------------------------
!     contributions due to diverging geometry


      do i = 5,nzn5
         geosrc  = hdt * rho(i-1) * u(i-1) * dloga(i-1)
         vlft(i) = 1.0_rk/vlft(i) - geosrc
         vlft(i) = 1.0_rk/vlft(i)
         plft(i) = plft(i) - geosrc * ce(i-1)**2
         plft(i) = max( plft(i) , config%smallp )
      end do

!     ------------
!     RIGHT STATES
!     ------------


      do i = 5,nzn5
         scr1       = 0.5_rk*max( 0.0_rk,-slamm(i) )
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         pltd(i)    = pl(i)    + scr1*(dp(i)    + scr2 * p6(i))
         rholtd(i)  = rhol(i)  + scr1*(drho(i)  + scr2 * rho6(i))
         ultd(i)    = ul(i)    + scr1*(du(i)    + scr2 * u6(i))
         utltd(i)   = utl(i)   + scr1*(dut(i)   + scr2 * ut6(i))
         uttltd(i)  = uttl(i)  + scr1*(dutt(i)  + scr2 * utt6(i))
         gameltd(i) = gamel(i) + scr1*(dgame(i) + scr2 * game6(i))
         gamcltd(i) = gamcl(i) + scr1*(dgamc(i) + scr2 * gamc6(i))

         pltd(i)    = max( pltd(i) , config%smallp )
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scr1       = -0.5_rk*slamp(i)
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         ppl(i)     = pl(i)    + scr1 * (dp(i)    + scr2 * p6(i))
         upl(i)     = ul(i)    + scr1 * (du(i)    + scr2 * u6(i))
         rhopl(i)   = rhol(i)  + scr1 * (drho(i)  + scr2 * rho6(i))
         gamcpl(i)  = gamcl(i) + scr1 * (dgamc(i) + scr2 * gamc6(i))

         gravpl(i)  = gravl(i) + scr1 * (dgrav(i) + scr2 * grav6(i))

         ppl(i)     = max( ppl(i) , config%smallp )
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scr1       = -0.5_rk*slamm(i)
         scr2       = 1.0_rk - 4.0_rk/3.0_rk*scr1

         pml(i)     = pl(i)    + scr1 * (dp(i)    + scr2 * p6(i))
         uml(i)     = ul(i)    + scr1 * (du(i)    + scr2 * u6(i))
         rhoml(i)   = rhol(i)  + scr1 * (drho(i)  + scr2 * rho6(i))
         gamcml(i)  = gamcl(i) + scr1 * (dgamc(i) + scr2 * gamc6(i))

         gravml(i)  = gravl(i) + scr1 * (dgrav(i) + scr2 * grav6(i))

         pml(i)     = max( pml(i) , config%smallp )
#ifdef SPLIT_LOOPS
      end do


      do i = 5,nzn5
#endif
         scrch1(i)  = -0.5_rk*slam0(i)
         scrch2(i)  = 1.0_rk - 4.0_rk/3.0_rk*scrch1(i)

         p0l(i)     = pl(i) &
                     +scrch1(i) * (dp(i)    + scrch2(i) * p6(i))
         rho0l(i)   = rhol(i) &
                     +scrch1(i) * (drho(i)  + scrch2(i) * rho6(i))
         ut0l(i)    = utl(i) &
                     +scrch1(i) * (dut(i)   + scrch2(i) * ut6(i))
         utt0l(i)   = uttl(i) &
                     +scrch1(i) * (dutt(i)  + scrch2(i) * utt6(i))
         game0l(i)  = gamel(i) &
                     +scrch1(i) * (dgame(i) + scrch2(i) * game6(i))
         gamc0l(i)  = gamcl(i) &
                     +scrch1(i) * (dgamc(i) + scrch2(i) * gamc6(i))

         p0l(i)     = max( p0l(i) , config%smallp )
      end do   

!     ------------------------------------------
!     right states for fluids (unrolled 5 times)

      n0 = 1
      do n = 1,config%qn-4,5
         do i = 5,nzn5
            xn0l(i,n  ) = xnl(i,n  ) &
                   +scrch1(i) * (dxn(i,n  ) + scrch2(i) * xn6(i,n  ))
            xn0l(i,n+1) = xnl(i,n+1) &
                   +scrch1(i) * (dxn(i,n+1) + scrch2(i) * xn6(i,n+1))
            xn0l(i,n+2) = xnl(i,n+2) &
                   +scrch1(i) * (dxn(i,n+2) + scrch2(i) * xn6(i,n+2))
            xn0l(i,n+3) = xnl(i,n+3) &
                   +scrch1(i) * (dxn(i,n+3) + scrch2(i) * xn6(i,n+3))
            xn0l(i,n+4) = xnl(i,n+4) &
                   +scrch1(i) * (dxn(i,n+4) + scrch2(i) * xn6(i,n+4))
         end do
         n0 = n+5
      end do

      do n = n0, config%qn
         do i = 5,nzn5
            xn0l(i,n  ) = xnl(i,n  ) &
                        + scrch1(i) * (dxn(i,n  ) + scrch2(i) * xn6(i,n  ))
        end do
      end do


!     ----------------------
!     effective right states
      do i = 5,nzn5
         crght(i)  =  gamcltd(i) * pltd(i) * rholtd(i)
      enddo

      call fastsqrt(crght, nzn5)

      do i = 5,nzn5

         xcrght    = 1.0_rk / crght(i)
         rhortdi   = 1.0_rk / rholtd(i)

         betapr    = -0.5_rk*( (ultd(i) - upl(i)) &
                            +(pltd(i) - ppl(i)) * xcrght &
                           )
         betamr    = +0.5_rk*( (ultd(i) - uml(i)) &
                            -(pltd(i) - pml(i)) * xcrght &
                           )
         beta0r    = (pltd(i) - p0l(i)) * xcrght * xcrght &
                    +rhortdi &
                    -1.0_rk/rho0l(i)

!        ------------------------------
!        fix for wrong formulae in CW84

         if ( slamp(i) .ge. 0.0_rk ) then 
            betapr    = 0.0_rk
            gravpl(i) = 0.0_rk
         end if
         if ( slamm(i) .ge. 0.0_rk ) then 
            betamr    = 0.0_rk
            gravml(i) = 0.0_rk
         end if
         if ( slam0(i) .ge. 0.0_rk ) beta0r = 0.0_rk

         if ( slamm(i) .lt. 0.0_rk .and. slamp(i) .lt. 0.0_rk ) then
            fdt = qdt
         else
            fdt = hdt
         end if

         vrght(i)  = rhortdi - beta0r - ( betapr + betamr )*xcrght
         prght(i)  = pltd(i) + crght(i)*( betapr + betamr )
         prght(i)  = max( prght(i) , config%smallp )
         urght(i)  = ultd(i) + ( betapr - betamr ) &
                     +fdt*(gravml(i)+gravpl(i))

         if ( slam0(i) .ge. 0.0_rk ) then
            utrght(i) = utltd(i)
            uttrgt(i) = uttltd(i)
            gmergt(i) = gameltd(i)
            gmcrgt(i) = gamcltd(i)
         else
            utrght(i) = ut0l(i)
            uttrgt(i) = utt0l(i)
            gmergt(i) = game0l(i)
            gmcrgt(i) = gamc0l(i)
         end if
      end do

!     ----------------------------------------------------
!     effective right states for fluids (unrolled 5 times)

      n0 = 1
      do n = 1,config%qn-4,5
         do i = 5,nzn5
            xnrght(i,n  ) = xn0l(i,n  )
            xnrght(i,n+1) = xn0l(i,n+1)
            xnrght(i,n+2) = xn0l(i,n+2)
            xnrght(i,n+3) = xn0l(i,n+3)
            xnrght(i,n+4) = xn0l(i,n+4)
         end do
         n0 = n+5
      end do

      do n = n0,config%qn
         do i = 5,nzn5
            xnrght(i,n  ) = xn0l(i,n  )
         end do
      end do

!     ---------------------------------------
!     contributions due to diverging geometry


      do i = 5,nzn5
         geosrc   = hdt * rho(i) * u(i) * dloga(i)
         vrght(i) = 1.0_rk/vrght(i) - geosrc
         vrght(i) = 1.0_rk/vrght(i)
         prght(i) = prght(i) - geosrc * ce(i)**2
         prght(i) = max( prght(i) , config%smallp )
      end do

!     -----------
!     limit Gamma


      do i = 5,nzn5
         gmemin(i) = min( game(i-2), game(i-1), game(i), game(i+1))
         gmemax(i) = max( game(i-2), game(i-1), game(i), game(i+1))
         gmelft(i) = max(gmemin(i), min(gmemax(i), gmelft(i) ))
         gmergt(i) = max(gmemin(i), min(gmemax(i), gmergt(i) ))
      end do

!     -----------------------------------------
!     force exact reflecting boundary condition

      if ( bndmin .eq. 1 .or. (bndmin .eq. 5 .and. are_id .eq. 1) ) then
         ulft  (5) = -urght (5) + 2.0_rk * ugrdl(5) ! urell_L = -urell_R
         utlft (5) =  utrght(5)
         uttlft(5) =  uttrgt(5)
         plft  (5) =  prght (5)
         vlft  (5) =  vrght (5)
         gmelft(5) =  gmergt(5)
         gmclft(5) =  gmcrgt(5)
         do n = 1,config%qn
            xnlft(5,n) = xnrght(5,n)
         end do
      end if

      if ( bndmax .eq. 1 ) then
         urght (nzn5) = -ulft  (nzn5) + 2.0_rk*ugrdl(nzn5)
         utrght(nzn5) =  utlft (nzn5)
         uttrgt(nzn5) =  uttlft(nzn5)
         prght (nzn5) =  plft  (nzn5)
         vrght (nzn5) =  vlft  (nzn5)
         gmergt(nzn5) =  gmelft(nzn5)
         gmcrgt(nzn5) =  gmclft(nzn5)
         do n = 1,config%qn
            xnrght(nzn5,n) = xnlft(nzn5,n)
         end do
      end if

!     ----------------
!     restore velocity


      do i = 4,nzn6
         u(i) = urel(i)
      end do
!
! return
!
      return

!
! end of STATES
!
    end subroutine states
!> \verbatim
!> solve riemann shock tube problem
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
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
      subroutine riemann

!
      use precision

      use intgrs_hy
      
      use gfloat_hy
      use lfloat_hy
      use hydro_hy
      use reman_hy
      use specfun, only : fastsqrt
      use state
      use configure
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i,n
      real(kind=rk) :: u_side,delu2,delu1,ustrr2,ustrl2
      real(kind=rk) :: ustrr1,aux,gamfac,smalldp,ge,gc,ustrl1
      real(kind=rk)::  pstar1(config%q), pstar2(config%q), gmstrl(config%q), gmstrr(config%q), &
                       wlft1 (config%q), wrght1(config%q)

      real(kind=rk) :: scrch1(config%q), scrch2(config%q), scrch3(config%q), scrch4(config%q)
      real(kind=rk) :: scratch(config%q), scratch2(config%q), scratch3(config%q), scratch4(config%q)
!---------------------------------------------------------------------
!
!     Note that while the left and right interface values of
!     zone i are stored in the arrays
!
!             ul(i), pl(i), etc.,  and  ur(i), pr(i), etc.
!
!     the input states for the Riemann problem at the left
!     interface of zone i are stored in the arrays
!
!            uflt (i), plft (i), etc.,  and  urght(i), prght(i), etc.
!
!     The solution of the Riemann problem is stored in the arrays
!
!            uav(i), pav(i), rhoav(i), etc.
!
!---------------------------------------------------------------------

      smalldp = config%small

      scratch(:)  = 0._rk
      scratch2(:) = 0._rk
      scratch3(:) = 0._rk
      scratch4(:) = 0._rk

      do i = 5, nzn5
         pstar1(i) = prght(i) - plft(i) - crght(i) * (urght(i) - ulft(i))
         pstar1(i) = plft(i) + pstar1(i) * (clft(i) / (clft(i) + crght(i)))
         pstar1(i) = max (config%smallp, pstar1(i))
      end do


      do i = 5, nzn5
         ge        = 0.5_rk * (gmelft(i) + gmergt(i))
         gc        = 0.5_rk * (gmclft(i) + gmcrgt(i))
         gamfac    = (1.0_rk - ge / gc) * (ge - 1.0_rk)
         gmstrl(i) = gamfac * (pstar1(i) - plft(i))
         gmstrl(i) = gmelft(i) + 2.0_rk * gmstrl(i) /  (pstar1(i) + plft(i))
         gmstrr(i) = gamfac * (pstar1(i) - prght(i))
         gmstrr(i) = gmergt(i) + 2.0_rk * gmstrr(i) /  (pstar1(i) + prght(i))         
         gmstrl(i) = max(gmemin(i), min(gmemax(i), gmstrl(i) ))
         gmstrr(i) = max(gmemin(i), min(gmemax(i), gmstrr(i) ))
      end do

      do i = 5,nzn5
         scrch1(i) = pstar1(i) - (gmstrl(i) - 1.0_rk) * plft(i) &
                   / (gmelft(i) - 1.0_rk)
         if ( scrch1(i) .eq. 0.0_rk ) scrch1(i) = config%smallp

         wlft1(i)  = pstar1(i) + 0.5_rk * (gmstrl(i) - 1.0_rk) &
                     * (pstar1(i) + plft(i))

         scratch(i-4) =    abs (pstar1(i) - plft(i)) &
                         * abs (wlft1(i)) &
                         / abs (vlft(i) * scrch1(i)) 
      enddo

      call fastsqrt(scratch, (nzn5-5+1) )
      
      do i = 5,nzn5
         wlft1(i) = scratch(i-4)
      enddo

      scratch(:) = 0._rk

      do i = 5,nzn5

         scrch2(i) = pstar1(i) - (gmstrr(i) - 1.0_rk) * prght(i) &
                   / (gmergt(i) - 1.0_rk)
         if ( scrch2(i) .eq. 0.0_rk ) scrch2(i) = config%smallp
         wrght1(i) = pstar1(i) + 0.5_rk * (gmstrr(i) - 1.0_rk) &
                   * (pstar1(i) + prght(i))
         scratch(i-4) =    abs (pstar1(i) - prght(i)) &
                         * abs (wrght1(i)) &
                         / abs (vrght(i) * scrch2(i)) 

         scratch3(i-4) = 0.5_rk * (game(i) - 1.0_rk) / game(i)
      enddo
      
      call fastsqrt(scratch, (nzn5-5+1) ) 
      call fastsqrt(scratch3,(nzn5-5+1) )

      
      do i = 5,nzn5
         wrght1(i)   = scratch(i-4)
         scratch2(i) = scratch3(i-4)
      enddo

      scratch (:) = 0._rk
      scratch3(:) = 0._rk


      do i = 5,nzn5
         aux       = scratch2(i)
         if (abs (pstar1(i) /  plft(i) - 1.0_rk) .le. smalldp ) &
         wlft1(i)  = clft(i)
         wlft1(i)  = max (wlft1(i),  aux * clft(i))
         if (abs (pstar1(i) /  prght(i) - 1.0_rk) .le. smalldp ) &
         wrght1(i) = crght(i)
         wrght1(i) = max (wrght1(i), aux * crght(i))
      end do

      do i = 5, nzn5
         pstar2(i) = (prght(i) - plft(i)) - wrght1(i) &
                     * (urght(i) - ulft(i))
         pstar2(i) = plft(i) + pstar2(i) * (wlft1(i) &
                     / (wlft1(i) + wrght1(i)))
         pstar2(i) = max (config%smallp, pstar2(i))
      end do

!     -------------------------------
!     begin of Riemann iteration loop
!     -------------------------------

      do n = 1, config%nriem

         do i = 5, nzn5
            ge        = 0.5_rk * (gmelft(i) + gmergt(i))
            gc        = 0.5_rk * (gmclft(i) + gmcrgt(i))
            gamfac    = (1.0_rk - ge / gc) * (ge - 1.0_rk)
            gmstrl(i) = gamfac * (pstar2(i) - plft(i))
            gmstrl(i) = gmelft(i) + 2.0_rk * gmstrl(i) / &
                        (pstar2(i) + plft(i))
            gmstrr(i) = gamfac * (pstar2(i) - prght(i))
            gmstrr(i) = gmergt(i) + 2.0_rk * gmstrr(i) /  &
                        (pstar2(i) + prght(i))            
            gmstrl(i) = max(gmemin(i), min(gmemax(i), gmstrl(i) ))
            gmstrr(i) = max(gmemin(i), min(gmemax(i), gmstrr(i) ))
#ifdef SPLIT_LOOPS
         end do


         do i = 5,nzn5
#endif
            scrch1(i) = pstar2(i) - (gmstrl(i) - 1.0_rk) * plft(i) &
                      / (gmelft(i) - 1.0_rk)
            if ( scrch1(i) .eq. 0.0_rk ) scrch1(i) = config%smallp
            wlft  (i) = pstar2(i) + 0.5_rk * (gmstrl(i) - 1.0_rk) &
                      * (pstar2(i) + plft(i))
            scratch  (i-4) =   abs (pstar2(i) - plft(i)) &
                             * abs (wlft(i)) &
                             / abs (vlft(i) * scrch1(i))
         end do

         call fastsqrt(scratch, (nzn5-5+1) )

         do i = 5,nzn5
            wlft(i) = scratch(i-4)
         enddo

         scratch (:) = 0._rk
         scratch2(:) = 0._rk

         do i = 5,nzn5
            scrch2(i) = pstar2(i) - (gmstrr(i) - 1.0_rk) * prght(i) &
                      / (gmergt(i) - 1.0_rk)

            if ( scrch2(i) .eq. 0.0_rk ) scrch2(i) = config%smallp
            wrght(i)  = pstar2(i) + 0.5_rk * (gmstrr(i) - 1.0_rk) &
                      * (pstar2(i) + prght(i))

            scratch(i-4)  =   abs (pstar2(i) - prght(i))   &
                            * abs (wrght(i)) / abs (vrght(i) * scrch2(i)) 
            scratch3(i-4) = 0.5_rk * (game(i) - 1.0_rk) / game(i)
         end do

         call fastsqrt(scratch, (nzn5-5+1) )
         call fastsqrt(scratch3,(nzn5-5+1) )

         do i = 5,nzn5
            wrght(i) = scratch(i-4)
            scratch2(i) = scratch3(i-4)
         enddo

         scratch (:) = 0._rk
         scratch3(:) = 0._rk

         do i = 5,nzn5
            aux     = scratch2(i)
            if ( abs (pstar2(i) / plft(i) - 1.0_rk) .le. config%small ) &
                 wlft(i) = clft(i)
            wlft(i) = max (wlft(i), aux * clft(i))
            if ( abs (pstar2(i) / prght(i) - 1.0_rk) .le. config%small) &
                 &      wrght(i) = crght(i)
            wrght(i) = max (wrght(i), aux * crght(i))
#ifdef SPLIT_LOOPS
         end do


         do i = 5,nzn5
#endif
            ustrl1    =  ulft(i) - (pstar1(i) -  plft(i)) /  wlft1(i)
            ustrr1    = urght(i) + (pstar1(i) - prght(i)) / wrght1(i)
            ustrl2    =  ulft(i) - (pstar2(i) -  plft(i)) /   wlft(i)
            ustrr2    = urght(i) + (pstar2(i) - prght(i)) /  wrght(i)

            delu1     = ustrl1 - ustrr1
            delu2     = ustrl2 - ustrr2
            scrch1(i) = delu2  - delu1

            if ( scrch1(i) .eq. 0.0_rk ) then
               delu2     = 0.0_rk
               scrch1(i) = config%smallu
            end if

            pstar(i)  = pstar2(i) &
                      - delu2 * (pstar2(i) - pstar1(i))/ scrch1(i)
            pstar(i)  = max (config%smallp, pstar(i))

            pstar1(i) = pstar2(i)
            pstar2(i) = pstar(i)
            wlft1(i)  = wlft(i)
            wrght1(i) = wrght(i)
         end do
      end do

      scratch2(:) = 0._rk
!     -----------------------------
!     end of Riemann iteration loop
!     -----------------------------

      do i = 5, nzn5
         scrch3(i) = ulft (i) - (pstar(i) - plft(i) ) / wlft(i)
         scrch4(i) = urght(i) + (pstar(i) - prght(i)) / wrght(i)
         ustar(i)  = 0.5_rk * (scrch3(i) + scrch4(i))

         urell(i)  = ustar(i) - ugrdl(i)
         scrch1(i) = sign (1.0_rk, urell(i))

         scrch2(i) = 0.5_rk * ( 1.0_rk + scrch1(i) )
         scrch3(i) = 0.5_rk * ( 1.0_rk - scrch1(i) )
#ifdef SPLIT_LOOPS
      end do


      do i = 5, nzn5
#endif
         ps(i)    = plft(i)   * scrch2(i) + prght(i)  * scrch3(i)
         us(i)    = ulft(i)   * scrch2(i) + urght(i)  * scrch3(i)
         uts(i)   = utlft(i)  * scrch2(i) + utrght(i) * scrch3(i)
         utts(i)  = uttlft(i) * scrch2(i) + uttrgt(i) * scrch3(i)
         vs(i)    = vlft(i)   * scrch2(i) + vrght(i)  * scrch3(i)
         games(i) = gmelft(i) * scrch2(i) + gmergt(i) * scrch3(i)
         gamcs(i) = gmclft(i) * scrch2(i) + gmcrgt(i) * scrch3(i)
#ifdef SPLIT_LOOPS
      end do


      do i = 5, nzn5
#endif
         rhos(i)  = 1.0_rk / vs(i)
         rhos(i)  = max (config%smlrho, rhos(i))
         vs(i)    = 1.0_rk / rhos(i)
         ws(i)    = wlft(i) * scrch2(i) + wrght(i) * scrch3(i)
         scratch3(i-4) = gamcs(i) * ps(i) * vs(i)

         vstar(i)  = vs(i) - (pstar(i) - ps(i)) / ws(i) / ws(i)
         rhostr(i) = 1.0_rk / vstar(i)
         rhostr(i) = max (config%smlrho, rhostr(i))
         vstar(i)  = 1.0_rk / rhostr(i)
         scratch4(i-4) = gamcs(i) * pstar(i) * vstar(i)

      end do

      call fastsqrt(scratch3, (nzn5-5+1) )
      call fastsqrt(scratch2, (nzn5-5+1) )

      do i = 5, nzn5
         scratch (i) = scratch3(i-4)
         scratch2(i) = scratch4(i-4)
      enddo

      scratch3(:) = 0._rk
      scratch4(:) = 0._rk     

      do i = 5, nzn5
         ces(i)    = scratch(i)
         cestar(i) = scratch2(i)

         wes(i)    = ces(i)        - scrch1(i) * us(i)
         westar(i) = cestar(i)     - scrch1(i) * ustar(i)
         scrch4(i) = ws(i) * vs(i) - scrch1(i) * us(i)
         if ( pstar(i) - ps(i) .ge. 0.0_rk ) then 
            wes   (i) = scrch4(i)
            westar(i) = scrch4(i)
         end if
         wes(i)    = wes(i)    + scrch1(i) * ugrdl(i)
         westar(i) = westar(i) + scrch1(i) * ugrdl(i)
#ifdef SPLIT_LOOPS
      end do


      do i = 5, nzn5
#endif
         gamfac    = (1.0_rk - games(i) / gamcs(i)) * (games(i) - 1.0_rk)
         gmstar(i) = gamfac * (pstar(i) - ps(i))
         gmstar(i) = games(i) + 2.0_rk * gmstar(i) / (pstar(i) + ps(i))
         gmstar(i) = max(gmemin(i), min(gmemax(i), gmstar(i) ))
      end do

!     -----------------------------------------
!     compute correct state for rarefaction fan


      do i = 5, nzn5
         scrch1(i) = max (wes(i)-westar(i), wes(i)+westar(i), config%smallu)
         scrch1(i) = (wes(i) + westar(i)) / scrch1(i)
         scrch1(i) = 0.5_rk * (1.0_rk + scrch1(i))
         scrch2(i) =         1.0_rk - scrch1(i)

         rhoav(i)  = scrch1(i) * rhostr(i) + scrch2(i) * rhos (i)
         uav  (i)  = scrch1(i) * ustar(i)  + scrch2(i) * us(i)
         utav (i)  = uts(i)
         uttav(i)  = utts(i)
         pav   (i) = scrch1(i) * pstar(i)  + scrch2(i) * ps(i)
         gameav(i) = scrch1(i) * gmstar(i) + scrch2(i) * games(i)
      end do


      do i = 5, nzn5
         if ( westar(i) .ge. 0.0_rk ) then
            rhoav(i)  = rhostr(i)
            uav(i)    = ustar(i)
            pav(i)    = pstar(i)
            gameav(i) = gmstar(i)
         end if
         if ( wes(i) .lt. 0.0_rk ) then          
            rhoav(i)  = rhos (i)
            uav(i)    = us   (i)
            pav(i)    = ps   (i)
            gameav(i) = games(i)
         end if
         urell(i) = uav(i) - ugrdl(i)         
      end do


!     -------------------------
!     upwind state is known now

      do i = 5, nzn5
         u_side    = sign (1.0_rk, uav(i) - ugrdl(i) )
         scrch2(i) = 0.5_rk * ( 1.0_rk + u_side )
         scrch3(i) = 0.5_rk * ( 1.0_rk - u_side )
      end do

      do n = 1, config%qn
         do i = 5, nzn5
            xnav(i,n) = xnlft(i,n) * scrch2(i) + xnrght(i,n) * scrch3(i)
         end do
      end do
#ifdef CMA

!     ---------
!     CMA start
!     ---------

      if ( (config%qn-1) .gt. 1 ) then

         do i = 5,nzn5
            scrch1(i) = 0.0_rk
         end do

         do n = 1,config%qn-1
            do i = 5,nzn5
               scrch1(i) = scrch1(i) + xnav(i,n)
            end do
         end do

         do i = 5,nzn5
            scrch1(i) = 1.0_rk/scrch1(i)
         end do

         do n = 1,config%qn-1
            do i = 5,nzn5
               xnav(i,n) = xnav(i,n)*scrch1(i)
            end do
         end do

      end if

!     -------
!     CMA end
!     -------
#endif

      if ( bndmax .eq. 1 ) then

        uav  (nzn5) = 0.0_rk
        urell(nzn5) = 0.0_rk

      end if
!
! return
!
      return
!
! end of RIEMANN
!
    end subroutine riemann

!======================================================================
!> \verbatim
!> calculate coefficients of cubic interpolation polynomial
!> cf. C&W84 Eqs (1.6) and (1.7)
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
    subroutine coeff
!
      use precision


      use intgrs_hy

      use grd_hy
      use intrp_hy
      use state
      use configure
      implicit none
! LOCAL variables that are not in modules

       real(kind=rk):: scrch1(config%q), scrch2(config%q), scrch3(config%q), scrch4(config%q)         
       integer(kind=ik) :: i
!     
! checking equality of dx to a chached value might save a lot of time
      do i = 2, nzn8
         scrch1(i) = dx(i)     + dx(i-1)
         scrch2(i) = scrch1(i) + dx(i)
         scrch3(i) = scrch1(i) + dx(i-1)
      enddo
!
      do i = 2, nzn7
         scrch4(i) = dx(i)  /  ( scrch1(i) + dx(i+1) )
         coeff1(i) = scrch4(i) * scrch3(i)   / scrch1(i+1)
         coeff2(i) = scrch4(i) * scrch2(i+1) / scrch1(i)
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! coeff1 and coeff2 correspond to the coefficients of the
! addends in C&W84 eq. (1.7)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i = 2, nzn6
         scrch4(i) = 1._rk  /  ( scrch1(i) + scrch1(i+2) )
         coeff3(i) = -scrch4(i)*dx(i)   * scrch1(i)   / scrch3(i+1)
         coeff4(i) =  scrch4(i)*dx(i+1) * scrch1(i+2) / scrch2(i+1)
         coeff5(i) = dx(i) - &
                    2._rk * (dx(i+1) * coeff3(i) + dx(i) * coeff4(i))
         coeff5(i) = coeff5(i) / scrch1(i+1)
      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! coeff5 corresponds to the coefficients of the 2nd addend and
! of the 1st addend inside the braces of the 3rd addend in C&W84 eq. (1.6)
!, while coeff3 and coeff5 correspond to the 2nd and 3rd addends inside 
!  those braces
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
    end subroutine coeff

!======================================================================

!> \verbatim
!> search for contact discontinuities in variable a and
!     steepen the zone structure if necessary
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param al
!> \param a
!> \param ar
!> \param smalla
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!> 
    subroutine detect (al, a, ar, smalla)

      use precision

      use intgrs_hy

      use gfloat_hy
      use hydro_hy
      use grd_hy
      use intrp_hy

      use configure
      use state
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i
      real(kind=rk) :: scrch2i,scrch1i,scrch3i,scrch4i
      real(kind=rk) :: ak0,epsln,eta2,eta1,smalla
         real(kind=rk):: scrch1(config%q), scrch2(config%q), scrch3(config%q)
         real(kind=rk):: al(config%q), a(config%q), ar(config%q)
!
!
!------- the following parameters are set as in Colella and
!        Woodward (JCP, 54 (1984), 174)
!
      eta1  = 20._rk
      eta2  = 0.05_rk
      epsln = 0.01_rk
      ak0   = 0.1_rk
!
      do i = 2, nzn7
        scrch1(i) = dx(i) + dx(i-1)
        scrch2(i) = scrch1(i) + dx(i+1)
        scrch1(i) = (a(i) - a(i-1)) / scrch1(i)
      enddo
!
      do i = 2, nzn6
        scrch2(i) = (scrch1(i+1) - scrch1(i)) / scrch2(i)
      enddo
!
      do i = 2, nzn8
        scrch1(i) = x(i) - x(i-1)
        scrch1(i) = scrch1(i) * scrch1(i) * scrch1(i)
      enddo
!
      do i = 3, nzn5
        scrch3(i) = (scrch2(i-1) - scrch2(i+1)) * &
                    (scrch1(i)   + scrch1(i+1))
!
        scrch4i = a(i+1) - a(i-1)
        if (scrch4i .eq. 0._rk)  scrch4i = config%small * smalla
!
        scrch3(i) = scrch3(i) / ((x(i+1) - x(i-1)) * scrch4i)
        if (scrch2(i-1) * scrch2(i+1) .ge. 0._rk)  scrch3(i) = 0._rk
      enddo
!
!     scrch2 and scrch3 now contain finite difference approximations
!     to the second and third derivativess of a.
!
      do i = 3, nzn5
        scrch3i = scrch3(i)
!
        if (epsln * min( abs( a(i+1) ), abs( a(i-1) ) )  - &
            abs( a(i+1) - a(i-1) )  .ge.  0._rk)    scrch3i = 0._rk
!
        scrch3i = max (0._rk, min (eta1 * (scrch3i - eta2), 1._rk))
!
        scrch1i = abs (p  (i+1) - p  (i-1)) / min (p  (i+1), p  (i-1))
        scrch2i = abs (rho(i+1) - rho(i-1)) / min (rho(i+1), rho(i-1))
 
        if (game(i) * ak0 * scrch2i - scrch1i .lt. 0._rk) scrch3i =0._rk
!
!     scrch3 now contains the contact steepening coefficient
!
        scrch1i = a(i-1) + 0.5_rk * dela(i-1)
        scrch2i = a(i+1) - 0.5_rk * dela(i+1)
!
        al(i)   = al(i) + (scrch1i - al(i)) * scrch3i
        ar(i)   = ar(i) + (scrch2i - ar(i)) * scrch3i
      enddo
!
      return
    end subroutine detect

!======================================================================
!> \verbatim
!> flaten zone structure in regions where shocks are too thin
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!> 
      subroutine flaten

      use precision

      use intgrs_hy

      use gfloat_hy
      use hydro_hy
      use grd_hy
      use intrp_hy


      use configure
!      use vnew_hy
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i
      real(kind=rk) :: dptest,dp2,dutest,utest
       real(kind=rk):: scrch1(config%q), scrch2(config%q), scrch3(config%q)    
!
!
!------- This version of subroutine FLATEN only uses the simplest
!        form of dissipation as described in the appendix of
!        Colella and Woodward (JCP, 54 (1984), 174). Therefore
!        the only constants required are epsiln, omg1 and omg2.
!
       real(kind=rk),save  :: epsiln, omg1, omg2
      data  epsiln /0.33_rk/ ,  omg1 /0.75_rk/ ,  omg2 /10._rk/
!
!
      do i = 1, nzn8
         flatn (i) = 0._rk
         shockd(i) = 0._rk
      enddo
!
!
      do i = 2, nzn7
         dp(i)     = p(i+1) - p(i-1)
         du(i)     = u(i+1) - u(i-1)
         scrch1(i) = epsiln * min (p(i+1), p(i-1)) - abs( dp(i) )
         utest     = config%smallu - abs (du(i))
!
         if (utest .lt. 0._rk)  then
            dutest = du(i)
         else
            dutest = 0._rk
         end if
!
         if (scrch1(i) .lt. 0._rk)  then
            scrch1(i) = 1._rk
         else
            scrch1(i) = 0._rk
         end if
!
         if (du(i)  .ge. 0._rk)  scrch1(i) = 0._rk
         if (dutest .eq. 0._rk)  scrch1(i) = 0._rk
!
         shockd(i) = scrch1(i)
      enddo
!
!
      do i = 3, nzn6
         dp2    = p(i+2) - p(i-2)
         dptest = dp(i)
         if (dp2 .eq. 0._rk)              dp2    = config%smallp
         if (dp2 - config%smallp  .eq.  0._rk)   dptest = 0._rk
         scrch2(i) = dptest / dp2 - omg1
         scrch3(i) = scrch1(i) * max (0._rk, scrch2(i) * omg2)
      enddo
!
      do i = 4, nzn5
         if (dp(i) .lt. 0._rk)  then
            scrch2(i) = scrch3(i+1)
         else
            scrch2(i) = scrch3(i-1)
         end if
      enddo
!
      do i = 4, nzn5
         flatn(i) = max (scrch3(i), scrch2(i))
         flatn(i) = max (0._rk, min (1._rk, flatn(i)))
      enddo
!
!

      do i = 1, nzn8
!         flatn (i) = flatn(i) * (1.0 - igodu) + igodu
!         flatn1(i) = 1.0 - flatn(i)
         scrch2(i) = 0._rk
      enddo
!
!
!      do i = 1,nzn+1
!         i4 = i+4
!         sh_on = max(ishck(i-1,1,1),ishck(i,1,1))
!         if ( sh_on .ne. 0.d0 ) flatn(i4) = 1.0
!      end do         
!
!      do i = 1,nzn+1
!         write(*,'(3I4)') i,ishck(i,1,1),ishck(i,0,0)
!      end do         
!      stop 'flatn'

      

      return
    end subroutine flaten

!======================================================================
!> \verbatim
!> interpolate interface values and monotonize
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param al
!> \param a
!> \param ar
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!> 
      subroutine interp (al, a, ar)

      use precision

      use intgrs_hy

      use intrp_hy

      use configure
      implicit none
! LOCAL variables that are not modules

      integer(kind=ik) :: i

         real(kind=rk)::scrch1(config%q), scrch2(config%q), scrch3(config%q), scrch4(config%q)         
         real(kind=rk):: al(config%q), a(config%q), ar(config%q)
!
!
      do i = 2, nzn8
      scrch1(i) = a(i) - a(i-1)
      scrch2(i) = abs ( scrch1(i) + scrch1(i) )
      scrch4(i) = sign (1._rk, scrch1(i))
      enddo
!
      do i = 2, nzn6
        dela(i) = coeff1(i) * scrch1(i+1) + coeff2(i) * scrch1(i)
        if (dela(i) .lt. 0._rk)  then
           scrch3(i) = -1._rk
        else
           scrch3(i) = +1._rk
        end if
        dela(i) = min(abs(dela(i)), scrch2(i), scrch2(i+1))* scrch3(i)
        if (-scrch4(i) * scrch4(i+1)  .ge.  0._rk)   dela(i) = 0._rk
      enddo
!
      do i = 2, nzn5
      ar(i)  = a (i) + coeff5(i) * scrch1(i+1) + coeff3(i) * dela(i+1)
      ar(i)  = ar(i) + coeff4(i) * dela(i)
      al(i+1)= ar(i)
      enddo
!
      return
    end subroutine interp
!======================================================================
!> \verbatim
!> calculate interface values of all variables
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!> 
      subroutine intrfc
!
!     
!
      use precision

      use intgrs_hy

      use gfloat_hy
      use hydro_hy
      use intrp_hy
      use physcs_hy

      use specfun


      use configure
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i,n,nucy,nuc0,nucbeg,index,nfound
      real(kind=rk) :: flatm1

      real(kind=rk) ::dcheck(nzn5),checkl(nzn5),checkr(nzn5),check(nzn5)
      integer(kind=ik)::  indchk(nzn5)
!
!
      call coeff
!
      call interp (rhol, rho, rhor)
      call detect (rhol, rho, rhor, config%smlrho)
!
      do n = 1, config%qn
        call interp (xnl(1:config%q,n), xn(1:config%q,n), xnr(1:config%q,n))
        call detect (xnl(1:config%q,n), xn(1:config%q,n), xnr(1:config%q,n), config%smallx)
      enddo
!
      call interp (ul,    u,    ur   )
      call interp (utl,   ut,   utr  )
      call interp (uttl,  utt,  uttr )
      call interp (pl,    p,    pr   )
      call interp (gamel, game, gamer)
      call interp (gamcl, gamc, gamcr)
!
      call flaten
!
      do i = 4, nzn5
        flatm1 = 1.e0_rk - flatn(i)
        rhol (i) = flatn(i) * rho (i) + flatm1 * rhol (i)
        rhor (i) = flatn(i) * rho (i) + flatm1 * rhor (i)
        ul   (i) = flatn(i) * u   (i) + flatm1 * ul   (i)
        ur   (i) = flatn(i) * u   (i) + flatm1 * ur   (i)
        utl  (i) = flatn(i) * ut  (i) + flatm1 * utl  (i)
        utr  (i) = flatn(i) * ut  (i) + flatm1 * utr  (i)
        uttl (i) = flatn(i) * utt (i) + flatm1 * uttl (i)
        uttr (i) = flatn(i) * utt (i) + flatm1 * uttr (i)
        pl   (i) = flatn(i) * p   (i) + flatm1 * pl   (i)
        pr   (i) = flatn(i) * p   (i) + flatm1 * pr   (i)
        gamel(i) = flatn(i) * game(i) + flatm1 * gamel(i)
        gamer(i) = flatn(i) * game(i) + flatm1 * gamer(i)
        gamcl(i) = flatn(i) * gamc(i) + flatm1 * gamcl(i)
        gamcr(i) = flatn(i) * gamc(i) + flatm1 * gamcr(i)
      enddo
!
!
      nuc0 = 1
      do  n = 1, config%qn-3, 4
        do  i = 4, nzn5
          flatm1 = 1._rk - flatn(i)
          xnl(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnl(i,n  )
          xnr(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnr(i,n  )
          xnl(i,n+1) = flatn(i) * xn(i,n+1) + flatm1 * xnl(i,n+1)
          xnr(i,n+1) = flatn(i) * xn(i,n+1) + flatm1 * xnr(i,n+1)
          xnl(i,n+2) = flatn(i) * xn(i,n+2) + flatm1 * xnl(i,n+2)
          xnr(i,n+2) = flatn(i) * xn(i,n+2) + flatm1 * xnr(i,n+2)
          xnl(i,n+3) = flatn(i) * xn(i,n+3) + flatm1 * xnl(i,n+3)
          xnr(i,n+3) = flatn(i) * xn(i,n+3) + flatm1 * xnr(i,n+3)
        end do
        nuc0 = n + 4
      end do
!
      nucbeg = nuc0
      do  n = nucbeg, config%qn-1, 2
        do  i = 4, nzn5
          flatm1 = 1._rk - flatn(i)
          xnl(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnl(i,n  )
          xnr(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnr(i,n  )
          xnl(i,n+1) = flatn(i) * xn(i,n+1) + flatm1 * xnl(i,n+1)
          xnr(i,n+1) = flatn(i) * xn(i,n+1) + flatm1 * xnr(i,n+1)
        end do
        nuc0 = n + 2
      end do
!
      nucbeg = nuc0
      do  n = nucbeg, config%qn
        do  i = 4, nzn5
          flatm1 = 1._rk - flatn(i)
          xnl(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnl(i,n  )
          xnr(i,n  ) = flatn(i) * xn(i,n  ) + flatm1 * xnr(i,n  )
        end do
      end do
!
!
      call monot (rhol,  rho,  rhor,  drho,  rho6 )
      call monot (ul,    u,    ur,    du,    u6   )
      call monot (utl,   ut,   utr,   dut,   ut6  )
      call monot (uttl,  utt,  uttr,  dutt,  utt6 )
      call monot (pl,    p,    pr,    dp,    p6   )
      call monot (gamel, game, gamer, dgame, game6)
      call monot (gamcl, gamc, gamcr, dgamc, gamc6)
!
      do  n = 1, config%qn
      call monot (xnl(1:config%q,n), xn(1:config%q,n), xnr(1:config%q,n), dxn(1:config%q,n), xn6(1:config%q,n))
      enddo
!
!
!------- if necessary flaten distributions of mass fractions to
!        guarantee mass conservation (because of non-linear advection)
!
      if (config%qn .lt. 2)  then

      do i = 4, nzn5
        vl(i) = 1._rk / rhol(i)
        v (i) = 1._rk / rho (i)
        vr(i) = 1._rk / rhor(i)
      enddo
      return
      endif
!
      do  i = 4, nzn5
        checkl(i) = xnl(i,1) - 1._rk
        checkr(i) = xnr(i,1) - 1._rk
        check (i) = xn (i,1) - 1._rk
      end do
!
      nucy = config%qn - 1             !!!   xn(., nuc) <==> Ye
!
      nuc0 = 2
      do  n = 2, nucy-3, 4
        do  i = 4, nzn5  
          checkl(i) = checkl(i) + xnl(i,n  ) + xnl(i,n+1) &
                                + xnl(i,n+2) + xnl(i,n+3)
          checkr(i) = checkr(i) + xnr(i,n  ) + xnr(i,n+1) &
                                + xnr(i,n+2) + xnr(i,n+3)
          check (i) = check (i) + xn (i,n  ) + xn (i,n+1) &
                                + xn (i,n+2) + xn (i,n+3)
        end do
        nuc0 = n + 4
      end do
!
      nucbeg = nuc0
      do  n = nucbeg, nucy-1, 2
        do  i = 4, nzn5  
          checkl(i) = checkl(i) + xnl(i,n  ) + xnl(i,n+1)
          checkr(i) = checkr(i) + xnr(i,n  ) + xnr(i,n+1)
          check (i) = check (i) + xn (i,n  ) + xn (i,n+1)
        end do
        nuc0 = n + 2
      end do
!
      do  n = nuc0, nucy
        do  i = 4, nzn5  
          checkl(i) = checkl(i) + xnl(i,n  ) 
          checkr(i) = checkr(i) + xnr(i,n  ) 
          check (i) = check (i) + xn (i,n  ) 
        end do
      end do
!
      do  i = 4, nzn5  
        dcheck(i) = max( &
                         abs (checkl(i)) - abs (check(i)), &
                         abs (checkr(i)) - abs (check(i)) )
      end do
!
!DIR$ INLINE
      call whenfgt_v( nzn5-3, dcheck(4), 1_ik, 1.e-3_rk, indchk, nfound )
!
      if (nfound .gt. 0) then
!
        do n = 1, nucy
!DIR$ IVDEP
          do i = 1, nfound
            index = indchk(i) + 3
            xnl(index,n) = xn(index,n)
            xnr(index,n) = xn(index,n)
            dxn(index,n) = 0._rk
            xn6(index,n) = 0._rk
         enddo
      enddo
!
      end if

      do i = 4, nzn5
        vl(i) = 1._rk / rhol(i)
        v (i) = 1._rk / rho (i)
        vr(i) = 1._rk / rhor(i)
      enddo
!
      return
    end subroutine intrfc

!======================================================================
#ifdef IBM_COMPILER
!     Problems with the -qhot option of xlf90 (V9.1) for monot() 
@PROCESS NOHOT
#endif

!> \verbatim
!> apply monotonicity constraint to interpolation parabola
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
!> \param al
!> \param a
!> \param ar
!> \param da
!> \param a6
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!> 
      subroutine monot (al, a, ar, da, a6)

      use precision

      use intgrs_hy

      use configure
      implicit none
! LOCAL variables that are not in modules
      integer(kind=ik) :: i
      real(kind=rk) :: disval

       real(kind=rk):: scrch1(config%q), scrch2(config%q), scrch3(config%q)
       real(kind=rk):: al(config%q), a(config%q), ar(config%q), da(config%q), a6(config%q)

!
!
      do i = 4, nzn5
         da(i) = (ar(i) - al(i))
         da(i) = sign (1._rk, da(i))
         scrch1(i) = (ar(i) - a(i)) * (al(i) - a(i))
!
         if (scrch1(i) .ge. 0._rk)  then
            ar(i) = a(i)
            al(i) = a(i)
         end if
!
         disval = (ar(i) - a(i)) * (al(i) - a(i))
!
         if (disval .ne. 0._rk)  then
            scrch2(i) = 3._rk * a(i) - 2._rk * ar(i)
            scrch3(i) = 3._rk * a(i) - 2._rk * al(i)
         else
            scrch2(i) = al(i)
            scrch3(i) = ar(i)
         end if
!
         if (da(i) * (al(i) - scrch2(i)) .lt. 0._rk)  al(i) = scrch2(i)
         if (da(i) * (scrch3(i) - ar(i)) .lt. 0._rk)  ar(i) = scrch3(i)
!
         da(i) = ar(i) - al(i)
         a6(i) = 6._rk * a(i) - 3. * (al(i) + ar(i))
      enddo
!
      return
    end subroutine monot

!======================================================================


!> \verbatim
!> 
!> compute fictitious forces (centrifugal and coriolis)
!>     for cylindrical and spherical coordinates
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
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
    subroutine force (j, k)


      use precision

      use intgrs_hy

      use mesh_hy
      use hydro_hy
      use grd_hy
      use physcs_hy
      use abort
      use specfun, only : fastsin, fastcos

      use configure
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i,j,k
      real(kind=rk):: dens(config%q), dmass(config%q), tmass(config%q), scratch(config%q), scratch2(config%q)
!
!
!
!      go to (10, 30, 50, 70, 90, 110), igeom + 1
!

      scratch (:) = 0._rk
      scratch2(:) = 0._rk

      if (igeom .eq. 0) then
!
!-----planar
!
      do i = 1, nzn8
         fict(i) = 0._rk
      enddo
      return
      endif

      if (igeom .eq. 1) then
!
!-----cylindrical (radial)
!
      if (config%igeomy .eq. 3) then
         do i = 1, nzn8
            fict(i) = ut(i) * ut(i) / x(i)
         enddo
      else if (config%igeomz .eq. 3) then
         do i = 1, nzn8
            fict(i) = utt(i) * utt(i) / x(i)
         enddo
      else
         raise_abort("force(): check geometry")
      endif
      return
      endif
!
      if (igeom .eq. 2) then
!
!-----spherical (radial)
!
      do i = 1, nzn8
         fict(i) = (ut(i)*ut(i) + utt(i)*utt(i)) / x(i)
      enddo
      return
      endif

      if (igeom .eq. 3) then
!
!-----cylindrical (angular)
!
      do i = 1, nzn8
         fict(i) = -u(i) * ut(i) / xzn(j)
      enddo
      return
      endif

      if (igeom .eq. 4) then
!
!-----spherical (angular - theta)
!
         do i = 1, nzn8
            scratch(i) = x(i)
            scratch2(i) = x(i)
         enddo

         call fastcos(scratch, nzn8)
         call fastsin(scratch2, nzn8)

         do i = 1, nzn8
            fict(i) = u(i)*ut(i) - utt(i)*utt(i) * scratch(i) / scratch2(i)
            fict(i) = -fict(i) / xzn(j)
         enddo
         return
      endif
!
      if (igeom .eq. 5) then
!
!-----spherical (angular - phi)
!
         do i = 1, nzn8
            fict(i) = 0.0_rk
         enddo
!
      endif
      return
    end subroutine force

!======================================================================

!> \verbatim
!> define arrays needed to calculate geometric source terms
!> \endverbatim
!>
!> \author W. Keil and M. Rampp
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
    subroutine geom (j, k)

      use precision


      use intgrs_hy

      use mesh_hy
      use grd_hy
      use abort
      use specfun, only : fastsin, fastcos


      use configure
      implicit none
! LOCAL variables that are not in modules

      integer(kind=ik) :: i,j,k
      real(kind=rk) :: scratch(config%q), scratch2(config%q), scratch3(config%q)

      scratch (:) = 0._rk
      scratch2(:) = 0._rk
      scratch3(:) = 0._rk

!-------  Cartesian coordinate 
!
      if (igeom .eq. 0)  then
!
         do i = 1, nzn8
           areal(i) = 1._rk
           arear(i) = 1._rk
           area (i) = 1._rk
           dvol (i) = dx(i)
         enddo
!
      end if
!
!
!-------  radial cylindrical coordinate 
!
      if (igeom .eq. 1)  then
!
         do i = 1, nzn8
           areal(i) = abs (xl(i))
           arear(i) = abs (xr(i))
           area (i) = 0.5_rk * (arear(i) + areal(i))
           dvol (i) = area(i) * dx(i)
        enddo
!
      end if
!
!
!-------  radial spherical coordinate )
!
      if (igeom .eq. 2)  then
!
         do i = 1, nzn8
           areal(i) = xl(i) * xl(i)
           arear(i) = xr(i) * xr(i)
           dvol (i) = (xr(i) * arear(i) - xl(i) * areal(i)) / 3._rk
           area (i) = dvol(i) / dx(i)
        enddo
!
      end if
!
!
!-------  angular cylindrical coordinate 
!
      if (igeom .eq. 3)  then

         do i = 1, nzn8
           areal(i) = 1._rk
           arear(i) = 1._rk
           area (i) = 1._rk
           dvol (i) = xzn(j) * dx(i) 
        enddo
!
      end if
!
!
!-------  angular spherical coordinate (theta)
!
      if (igeom .eq. 4)  then
!
         do i = 1, nzn8
            scratch(i) = xl(i)
            scratch2(i) = xr(i)
            scratch3(i) =0.5_rk * (xl(i) + xr(i))
         enddo
         
         call fastsin(scratch, nzn8)
         call fastsin(scratch2, nzn8)
         call fastsin(scratch3, nzn8)

         do i = 1, nzn8
           areal(i) = scratch(i)
           arear(i) = scratch2(i)
           area (i) = scratch3(i)
           dvol (i) = xzn(j) * area(i) * dx(i) 
        enddo
!
      end if
!
!
!-------  angular spherical coordinate (phi)
!
      if (igeom .eq. 5)  then
!
         if (config%igeomz .eq. 5) then        !  Case:  (theta, phi)
            do i = 1, nzn8
              areal(i) = 1._rk
              arear(i) = 1._rk
              area (i) = 1._rk
              dvol (i) = dx(i) * xzn(j) * sin (yzn(k))
           enddo
!
          else if (config%igeomy .eq. 5) then  !  Case:  (phi, theta)
            do i = 1, nzn8
              areal(i) = 1._rk
              arear(i) = 1._rk
              area (i) = 1._rk
              dvol (i) = dx(i) * xzn(j) * sin (zzn(k))
           enddo
!
         else
            raise_abort("geom(): check geometry")
!
         endif
!
      endif
!
!
      return
    end subroutine geom

end module ppm

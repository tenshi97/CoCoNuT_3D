module setbnd_mod

implicit none

contains

!=======================================================================

! subroutine setbnd(timare) ! removed since timare was never used
!> \verbatim
!>  task:      sets up special boundary conditions
!> \endvarbatim
!>
!>  \author W. Keil
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>   
!> \endverbatim
!>
  subroutine setbnd

!=======================================================================
! author          :  Wolfgang Keil
! Modul           : $Id: setbnd.F,v 1.3 2002/04/25 11:22:08 mjr Exp $
! Version         : $Revision: 1.3 $
! Date            : $Date: 2002/04/25 11:22:08 $
! Freeze Version  : $Name:  $
!     
!     task:      sets up special boundary conditions
!=======================================================================
  use precision
  
  use intgrs_hy 
  use bndinf_hy 
!  use revsho_hy   ! forcheck
  use totare_hy 
  use totgrq_hy 
  ! use arecon_hy

!  use vnew_hy   ! forcheck
!  use vold_hy   ! forcheck
  use mesh_hy
  use gfloat_hy
!  use lfloat_hy   ! forcheck
  use hydro_hy
  use grd_hy
  use physcs_hy
  use phycon

  !      use mesh_hy 
  !      use vnew_hy

  use mo_mpi


  use hydro_areas_mod
  use configure
  implicit none
! LOCAL variables that are not in modules

  ! integer(kind=ik) :: i1,i2  ! forcheck
  integer(kind=ik) :: i,n,i4,j,jj,k,kk,ioz,izf,ibo,izi
  INTEGER(KIND=ik) :: iyf,iyi,iox,ixi,IA,ixf,ioy,info
  real(kind=rk) :: utt2b,ut2b ! ,TIMARE
  ! real(kind=rk) :: udamp ! forcheck

  !      info = 1                   ! print control output
  info = 0                   ! no control output

!=======================================================================
!     special outer boundary conditions (radial direction):
!     neutronstar surface
!=======================================================================

  ia  = areas%are_nu
  ixi = areas%ix_are(ia, 1)
  ixf = areas%ix_are(ia, 2)
  iox = areas%ix_are(ia, 3)
  iyi = areas%ix_are(ia, 4)
  iyf = areas%ix_are(ia, 5)
  ioy = areas%ix_are(ia, 6)
  izi = areas%ix_are(ia, 7)
  izf = areas%ix_are(ia, 8)
  ioz = areas%ix_are(ia, 9)
  
  if(areas%ix_are(ia,13) .ne. 5) return

  !      ix_are(ia,13) = 5

!-----------------------------------------------------------------------
!     outflow damping:
!-----------------------------------------------------------------------

!-v03f      if(timare .lt. p_sot) then
!-v03f         udamp = min(1.,1. + (timare - p_sot)/p_sot)
!-v03f      else
!-v03f         udamp = 1.0
!-v03f      endif

!  udamp = 0.0_rk ! forcheck

!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------

  areas%nz = izf - izi
  areas%ny = iyf - iyi
  areas%nx = ixf - ixi
  
  areas%nz = areas%nz/ioz + 1
  areas%ny = areas%ny/ioy + 1
  areas%nx = areas%nx/iox + 1

  if(info .eq. 1) then
     if (myproc .eq. 0) then
        write(*,*) "Task ",myproc
        write(*,'(''setbnd-test> ia = '',i3,'' ny = '',i3,'' nz = '',i3)') ia,areas%ny,areas%nz

        write(*,'(''setbnd-test> ixi = '',i3,'' ixf = '',i3, &
             & '' iox = '',i3,'' nuc = '',i2)') ixi, ixf, iox, config%qn
        write(*,'(''setbnd-test> iyi = '',i3,'' iyf = '',i3, &
             & '' ioy = '',i3)') iyi, iyf, ioy
        write(*,'(''setbnd-test> izi = '',i3,'' izf = '',i3, &
             &    '' ioz = '',i3)') izi, izf, ioz
     endif
  end if

  nzn  = config%qx
  nzn1 = nzn + 1
  nzn2 = nzn + 2
  nzn3 = nzn + 3
  nzn4 = nzn + 4
  nzn5 = nzn + 5
  nzn6 = nzn + 6
  nzn7 = nzn + 7
  nzn8 = nzn + 8

!-----------------------------------------------------------------------
!     gradient of the density profile in the atmosphere is choosen to
!     be independent of \theta:
!-----------------------------------------------------------------------

  j  = config%qy - 1
  k  = 1

!-----------------------------------------------------------------------
!     z/y-loop:
!-----------------------------------------------------------------------


  kk = qz_s
  do k = qz_s, qz_e
     jj = qy_s
     do j = qy_s, qy_e

!-----------------------------------------------------------------------
!     x-loop:
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     set up 1D-arrays as in "getrwx":
!-----------------------------------------------------------------------

        do i = 1, config%qx
           i4 = i + 4
           rho  (i4) = dentot(i, jj, kk)
           u    (i4) = vextot(i, jj, kk)
           ut   (i4) = veytot(i, jj, kk)
           utt  (i4) = veztot(i, jj, kk)
           e    (i4) = enetot(i, jj, kk)
           p    (i4) = pretot(i, jj, kk)
           tmp  (i4) = temtot(i, jj, kk)
           grav (i4) = 0.5_rk * (gpotot(i,jj,kk)+ gpotot(i,jj-1,kk))
           game (i4) = gaetot(i, jj, kk)
           gamc (i4) = gactot(i, jj, kk)
           ugrid(i4) = ugrtot(i)
           xl   (i4) = xzltot(i)
           x    (i4) = xzntot(i)
           xr   (i4) = xzrtot(i)
           dx   (i4) = xr(i4) - xl(i4)
           dvx  (i4) = (xr(i4)**3 - xl(i4)**3)/3._rk
        enddo
            
        grav (4) = 0.5_rk * (gpotot(0,jj,kk) + gpotot(0,jj-1,kk))
        
        do n = 1, config%qn
           do i = 1, config%qx
              i4 = i + 4
              xn(i4,n) = xnutot(i,jj,kk,n)
           enddo
        enddo

        do i = config%qx+1, config%qx+4
           i4 = i + 4
           dx   (i4) =        dx(i4-1)
           xl   (i4) =        xr(i4-1)
           xr   (i4) =        xr(i4-1) + dx(i4)
           x    (i4) = 0.5_rk * (xl(i4) + xr(i4))
           dvx  (i4) = (xr(i4)**3 - xl(i4)**3)/3._rk
        enddo
            
!-----------------------------------------------------------------------
!     here comes the actual boundary condition as in "bndry":
!     e.g. flow out boundary + a hydrostatic atmosphere
!-----------------------------------------------------------------------

!        i1 = nzn4 - 4  ! forcheck
!        i2 = nzn4 - 3  ! forcheck
            

        ut2b   = ut(nzn4)  * x(nzn4)
        utt2b  = utt(nzn4) * x(nzn4)

        do i = nzn5, nzn8
           e   (i)  = e   (nzn4)
         !--------------------------------------------------------------
         !   assume that there is nearly no matter in the atmosphere:
         !--------------------------------------------------------------
           grav(i)  = grav(nzn4) + pc_gc*tgmtot(config%qx)/xr(nzn4) &
                      - pc_gc*tgmtot(config%qx)/xr(i)
           gamc(i)  = gamc(nzn4)
           game(i)  = game(nzn4)
           ugrid(i) = ugrid(nzn4)
               
           u   (i)  = u(nzn4)
           ut  (i)  = ut2b  / x(i)
           utt (i)  = utt2b / x(i)

           rho(i) = rho (nzn4)*(x(nzn4)/x(i))**2 * (u(nzn4) + config%smallu)/(u(i) + config%smallu)

           do n = 1, config%qn
              xn  (i,n) = xn (nzn4,n)
           enddo
        enddo



        do i = nzn5, nzn8
           p(i) = max(p(nzn4),config%smallp)
           e(i) = max(e(nzn4),config%smalle)
        enddo

!-----------------------------------------------------------------------
!     copy boundary values back to the interface:
!-----------------------------------------------------------------------

        do i = nzn5, nzn8
           ibo = (i - nzn4) + ia*8 + 4
           den_bi(ibo, j , k ) =        rho (i)
           vex_bi(ibo, j , k ) =        u   (i)
           vey_bi(ibo, j , k ) =        ut  (i)
           vez_bi(ibo, j , k ) =        utt (i)
           ene_bi(ibo, j , k ) =        e   (i)
           gra_bi(ibo, j , k ) =        grav(i)
           gac_bi(ibo, j , k ) =        gamc(i)
           gae_bi(ibo, j , k ) =        game(i)
           pre_bi(ibo, j , k ) =        p   (i)
           dxx_bi(ibo)         =        dx  (i)
           ugr_bi(ibo)         =        ugrid (i)
           do n  = 1, config%qn
              xnu_bi(ibo,j,k,n) =       xn(i,n)
           enddo
        enddo


!-----------------------------------------------------------------------
!     reset 1D-arrays:
!-----------------------------------------------------------------------

        do i = 1, config%q
           rho  (i) =  0._rk 
           u    (i) =  0._rk
           ut   (i) =  0._rk
           utt  (i) =  0._rk
           e    (i) =  0._rk
           p    (i) =  0._rk
           tmp  (i) =  0._rk
           grav (i) =  0._rk
           game (i) =  0._rk
           gamc (i) =  0._rk
           ugrid(i) =  0._rk
           xl   (i) =  0._rk
           x    (i) =  0._rk
           xr   (i) =  0._rk
           dx   (i) =  0._rk
        enddo
            
        do n = 1, config%qn
           do i = 1, config%q
              xn(i,n) = 0._rk
           enddo
        enddo
            

!-----------------------------------------------------------------------
!     end of y/z-loop:
!-----------------------------------------------------------------------
        jj = jj + ioy
     enddo
     kk = kk + ioz
  enddo

  return
end subroutine setbnd
!=======================================================================


end module setbnd_mod

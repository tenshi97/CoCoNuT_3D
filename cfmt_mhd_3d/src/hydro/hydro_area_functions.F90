module hydro_area_functions

implicit none

contains

#undef DEBUG

subroutine hydro_are(iarea,iswp, selftime, childrentime)
!-----------------------------------------
! Autor           : Wolfgang Keil, Markus Rampp (MPA) 
! Modul           : $Id: rady.F,v 1.15 2006/05/11 15:44:04 rburas Exp $
! Version         : $Revision: 1.15 $
! Date            : $Date: 2006/05/11 15:44:04 $
!
!     task        :  get 3D calculation area and calculate it 
!                    separately 
!
!     input:   iarea  =  area ID
!
!=======================================================================
  use precision
  use abort
  use error

  use intgrs_hy
  use totare_hy
!  use arecon_hy
  use bndinf_hy
  use massio_hy
  use cputim
  use vnew_hy
  use mesh_hy
  use gfloat_hy
  use spez_hy
  use vnuw_hy

  use hlpare_hy
  use mo_mpi
  use setbnd_mod
  use sweeps_mod

  use hydro_areas_mod
  use configure
  use state

  implicit none
! LOCAL variables that are not in modules
  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk)                :: selftime_start(2), sweep_self(2), &
                                  sweep_children(2)
  integer(kind=ik)             :: i, j, k, ii, jj, kk, isd, iii, &
                                  jji, kki, iif, jjf, kkf
  real(kind=rk)                :: dt_old, dtnew, dvolx, dvoly, dvolz
  integer(kind=ik) ,intent(in) :: iarea, iswp
  real(kind=rk), dimension(2)  :: tim1, tim2
  integer(kind=ik)             :: j_end, k_end


  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif
!-----------------------------------------------------------------------
!     get control indices:
!-----------------------------------------------------------------------

  are_id = iarea
  ixi    = areas%ix_are(iarea, 1)
  ixf    = areas%ix_are(iarea, 2)
  iox    = areas%ix_are(iarea, 3)
  iyi    = areas%ix_are(iarea, 4)
  iyf    = areas%ix_are(iarea, 5)
  ioy    = areas%ix_are(iarea, 6)
  izi    = areas%ix_are(iarea, 7)
  izf    = areas%ix_are(iarea, 8)
  ioz    = areas%ix_are(iarea, 9)
  isd    = areas%ix_are(iarea,10)
!      ix_are(iarea,11) = 1
  config%bndmnx = areas%ix_are(iarea,12)
  config%bndmxx = areas%ix_are(iarea,13)

!-----------------------------------------------------------------------
!     the current time step is the time step of the area:
!-----------------------------------------------------------------------

  hydro%dt     = areas%dt_are(iarea)
  time   = areas%ti_are(iarea)    ! Attention: qterms and grdvel
                                ! need absolut time (end of time step)! 
  dt_old = hydro%dt

!-BN2D call grdvel               ! set grid velocity -> ugrtot(1..config%qx)

! removed since time is never used in setbnd
!  if (iswp.eq.1) call setbnd(time)         ! set inner boundaries
  if (iswp .eq. 1) call setbnd                ! set inner boundaries


!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------

  areas%nz = izf - izi
  areas%ny = iyf - iyi
  areas%nx = ixf - ixi

  areas%nz = areas%nz/ioz + 1
  areas%ny = areas%ny/ioy + 1
  areas%nx = areas%nx/iox + 1


!-----------------------------------------------------------------------
!     copy a part of the total calculation area in the PPM-arrays:
!-----------------------------------------------------------------------

  if (ioy .eq. 1 .and. ioz .eq. 1) then
     j_end=qy_e
     k_end=qz_e
  else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     j_end=qy_s
     k_end=qz_s
  end if

  velx  (1:areas%nx,qy_s:j_end,qz_s:k_end) = vextot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  vely  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veytot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  velz  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veztot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  densty(1:areas%nx,qy_s:j_end,qz_s:k_end) = dentot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  energy(1:areas%nx,qy_s:j_end,qz_s:k_end) = enetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  press (1:areas%nx,qy_s:j_end,qz_s:k_end) = pretot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  gammae(1:areas%nx,qy_s:j_end,qz_s:k_end) = gaetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  gammac(1:areas%nx,qy_s:j_end,qz_s:k_end) = gactot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  stot  (1:areas%nx,qy_s:j_end,qz_s:k_end) = stotot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  temp  (1:areas%nx,qy_s:j_end,qz_s:k_end) = temtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn) = xnutot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end,1:config%qn)

  ishck (1:areas%nx,qy_s:j_end,qz_s:k_end) = ishtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)

  ishck (0,   qy_s:j_end,qz_s:k_end) = ishtot(ixi-1,qy_s:j_end,qz_s:k_end)
  ishck (areas%nx+1,qy_s:j_end,qz_s:k_end) = ishtot(ixf+1,qy_s:j_end,qz_s:k_end)
  ishck (1:areas%nx,qy_s-1,      qz_s:k_end) = ishtot(ixi:ixf:iox,qy_s-1,qz_s:k_end)
  ishck (1:areas%nx,j_end+1,   qz_s:k_end) = ishtot(ixi:ixf:iox,j_end+1,qz_s:k_end)
  ishck (1:areas%nx,qy_s:j_end,qz_s-1    ) = ishtot(ixi:ixf:iox,qy_s:j_end,qz_s-1)
  ishck (1:areas%nx,qy_s:j_end,k_end+1   ) = ishtot(ixi:ixf:iox,qy_s:j_end,k_end+1)

  ishck (0,    qy_s-1,qz_s:k_end) = ishtot(ixi-1,qy_s-1,qz_s:k_end)
  ishck (areas%nx+1, qy_s-1,qz_s:k_end) = ishtot(ixf+1,qy_s-1,qz_s:k_end)
  ishck (0,   j_end+1,qz_s:k_end) = ishtot(ixi-1,j_end+1,qz_s:k_end)
  ishck (areas%nx+1,j_end+1,qz_s:k_end) = ishtot(ixf+1,j_end+1,qz_s:k_end)

  ishck (0,   qy_s:j_end,qz_s-1) = ishtot(ixi-1,qy_s:j_end,qz_s-1)
  ishck (areas%nx+1,qy_s:j_end,qz_s-1) = ishtot(ixf+1,qy_s:j_end,qz_s-1)
  ishck (0,   qy_s:j_end,k_end+1) = ishtot(ixi-1,qy_s:j_end,k_end+1)
  ishck (areas%nx+1,qy_s:j_end,k_end+1) = ishtot(ixf+1,qy_s:j_end,k_end+1)

  ishck (1:areas%nx, qy_s-1,qz_s-1) = ishtot(ixi:ixf:iox,qy_s-1,qz_s-1)
  ishck (1:areas%nx,j_end+1,qz_s-1) = ishtot(ixi:ixf:iox,j_end+1,qz_s-1)
  ishck (1:areas%nx, qy_s-1,k_end+1) = ishtot(ixi:ixf:iox,qy_s-1,k_end+1)
  ishck (1:areas%nx,j_end+1,k_end+1) = ishtot(ixi:ixf:iox,j_end+1,k_end+1)

  ishck (0,    qy_s-1,qz_s-1) = ishtot(ixi-1,qy_s-1,qz_s-1)
  ishck (areas%nx+1, qy_s-1,qz_s-1) = ishtot(ixf+1,qy_s-1,qz_s-1)
  ishck (0,   j_end+1,qz_s-1) = ishtot(ixi-1,j_end+1,qz_s-1)
  ishck (areas%nx+1,j_end+1,qz_s-1) = ishtot(ixf+1,j_end+1,qz_s-1)


  ishck (0,    qy_s-1,k_end+1) = ishtot(ixi-1,qy_s-1,k_end+1)
  ishck (areas%nx+1, qy_s-1,k_end+1) = ishtot(ixf+1,qy_s-1,k_end+1)
  ishck (0,   j_end+1,k_end+1) = ishtot(ixi-1,j_end+1,k_end+1)
  ishck (areas%nx+1,j_end+1,k_end+1) = ishtot(ixf+1,j_end+1,k_end+1)

  do k = qz_s, k_end
     do j = qy_s-1, j_end
        ii = ixi-1
        do i = 0, areas%nx
           gpot  (i,j,k) = gpotot(i+ii,j ,k)
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
!     Attention: For the z-direction igeomy = 5 (spherical, theta, phi) is
!     assumed! Since the resolution (ioz) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if (ioz .eq. 1) then 
     zznl(1:areas%nz) = zzltot(izi:izf:ioz)
     zzn(1:areas%nz)  = zzntot(izi:izf:ioz)
     zznr(1:areas%nz) = zzrtot(izi:izf:ioz)
     dvz(1:areas%nz)  = dvztot(izi:izf:ioz)
  else
#ifdef DEBUG
     write(*,*) 'detect_sh_are:izi, nz,ioz ',izi, areas%nz,ioz
#endif
     ii =  izi
     do i = 1, areas%nz
        zznl(i) = zzltot(ii)
        ii = ii + ioz
        if (ii .gt. config%qz+4) ii=1 !quick fix: ii could otherwise 
        !easily exceed q if ioz=qz
     enddo

     do  i = 1, areas%nz
        zznr(i) = zznl(i+1)
     enddo

     do i = 1, areas%nz
        dvz(i) = zznr(i) - zznl(i) !igeomz = 5 is assumed!
        zzn(i) = 0.5_rk * (zznr(i) + zznl(i))
     enddo
  endif


!-----------------------------------------------------------------------
!     Attention: For the y-direction igeomy = 4 (spherical, theta) is
!     assumed! Since the resolution (ioy) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if(ioy .eq. 1) then
     yznl(1:areas%ny)   = yzltot(iyi:iyf:ioy)
     yzn(1:areas%ny)    = yzntot(iyi:iyf:ioy)
     yznr(1:areas%ny)   = yzrtot(iyi:iyf:ioy)
     dvy(1:areas%ny)    = dvytot(iyi:iyf:ioy)
  else
#ifdef DEBUG
     write(*,*) 'detect_sh_are:iyi, ny,ioy ',iyi, areas%ny,ioy
#endif
     ii =  iyi
     do i = 1, areas%ny + 8
        yznl(i) = yzltot(ii)
        ii = ii + ioy
        if (ii .gt. config%qy+4) ii=1 !quick fix: ii could otherwise 
        !easily exceed q if ioy=qy
     enddo

     do  i = 1, areas%ny + 7
        yznr(i) = yznl(i+1)
     enddo

     do i = 1, areas%ny + 7
        dvy(i) = cos(yznl(i)) - cos(yznr(i)) !igeomy = 4 is assumed!
        yzn(i) = 0.5_rk * (yznr(i) + yznl(i))
     enddo
  endif

!      write(*,'(''A> yznl123.. '',6(1x,1pe12.5),'' ny= '',i3)')
!     &     yznl(1),yznl(2),yznl(3),yznl(ny-2),yznl(ny-1),yznl(ny),ny
!      write(*,'(''A> yzn123.. '',6(1x,1pe12.5),'' ny= '',i3)')
!     &     yzn(1),yzn(2),yzn(3),yzn(ny-2),yzn(ny-1),yzn(ny),ny
!      write(*,'(''A> yznr123.. '',6(1x,1pe12.5),'' ny= '',i3)')
!     &     yznr(1),yznr(2),yznr(3),yznr(ny-2),yznr(ny-1),yznr(ny),ny
!      write(*,'(''B> yzltot123.. '',6(1x,1pe12.5),'' config%qy= '',i3)')
!     &     yzltot(1),yzltot(2),yzltot(3),
!     &     yzltot(config%qy-2),yzltot(config%qy-1),yzltot(config%qy),config%qy
!      write(*,'(''B> yzntot123.. '',6(1x,1pe12.5),'' config%qy= '',i3)')
!     &     yzntot(1),yzntot(2),yzntot(3),
!     &     yzntot(config%qy-2),yzntot(config%qy-1),yzntot(config%qy),config%qy
!      write(*,'(''B> yzrtot123.. '',6(1x,1pe12.5),'' config%qy= '',i3)')
!     &     yzrtot(1),yzrtot(2),yzrtot(3),
!     &     yzrtot(config%qy-2),yzrtot(config%qy-1),yzrtot(config%qy),config%qy

!-----------------------------------------------------------------------
!     Attention: For the x-direction the resolution (iox) must be 1!
!     Otherwise the center of zone and boundaries of the zone are
!     different!!!
!-----------------------------------------------------------------------

  xznl(1:areas%nx)   = xzltot(ixi:ixf:iox)
  xzn(1:areas%nx)    = xzntot(ixi:ixf:iox)
  xznr(1:areas%nx)   = xzrtot(ixi:ixf:iox)
  dvx(1:areas%nx)    = dvxtot(ixi:ixf:iox)
  ugridx(1:areas%nx) = ugrtot(ixi:ixf:iox)

!-----------------------------------------------------------------------
!     1/2 of the original PROMETHEUS cycle: 
!-----------------------------------------------------------------------

  vxvold(1:areas%nx,qy_s:j_end,qz_s:k_end) = vexold(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  vyvold(1:areas%nx,qy_s:j_end,qz_s:k_end) = veyold(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  vzvold(1:areas%nx,qy_s:j_end,qz_s:k_end) = vezold(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  denold(1:areas%nx,qy_s:j_end,qz_s:k_end) = dnsold(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)

  select case(iswp)
  case(1)
     timer%hydro_sweep_mode = 2

     call sweep (1,.true.,.false., 1, sweep_self, sweep_children)  !  x - sweep
     timer%hydroare_sweep(timer%hydro_sweep_mode,:,1)=timer%hydroare_sweep(timer%hydro_sweep_mode,:,1)+sweep_self
     childrentime = childrentime + sweep_self
     timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,1)= &
        timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,1)+sweep_children

  case(2)
     if (config%nsdim .ge. 2 .and. areas%ny .ge. 4)  then  
        timer%hydro_sweep_mode = 2 
        call sweep (2,.true.,.false., 1, sweep_self, sweep_children)  !  !  y - sweep  (only 2D and 3D)
        timer%hydroare_sweep(timer%hydro_sweep_mode,:,2)=timer%hydroare_sweep(timer%hydro_sweep_mode,:,2)+sweep_self
        childrentime = childrentime + sweep_self
        timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,2)= &
                timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,2)+sweep_children
     endif
  case(3)
     if (config%nsdim .eq. 3 .and. areas%ny .ge. 4 .and. areas%nz.ge.4)  then  
        timer%hydro_sweep_mode = 2
        call sweep (3,.true.,.false., 1, sweep_self, sweep_children)  ! z - sweep  (only 3D) 
        timer%hydroare_sweep(timer%hydro_sweep_mode,:,3)=timer%hydroare_sweep(timer%hydro_sweep_mode,:,3)+sweep_self
        childrentime = childrentime + sweep_self
        timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,3)= &
                timer%hydroare_sweep_children(timer%hydro_sweep_mode,:,3)+sweep_children
     endif
  case default
     raise_abort("hydro_are(): nocase")
  end select


!-----------------------------------------------------------------------
!     copy PPM-arrays back to the total calculation regime:
!-----------------------------------------------------------------------

  if (ioy .eq. 1 .and. ioz .eq. 1) then
     j_end=qy_e
     k_end=qz_e
  else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     j_end=qy_s
     k_end=qz_s
  end if

  vextot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = velx  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  veytot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = vely  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  veztot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = velz  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  dentot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = densty(1:areas%nx,qy_s:j_end,qz_s:k_end)
  enetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = energy(1:areas%nx,qy_s:j_end,qz_s:k_end)
  pretot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = press (1:areas%nx,qy_s:j_end,qz_s:k_end)
  gaetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = gammae(1:areas%nx,qy_s:j_end,qz_s:k_end)
  gactot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = gammac(1:areas%nx,qy_s:j_end,qz_s:k_end)
  stotot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = stot  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  temtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = temp  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  xnutot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end,1:config%qn) = xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn)

  xzltot(ixi:ixf:iox) = xznl(1:areas%nx)   
  xzntot(ixi:ixf:iox) = xzn(1:areas%nx)    
  xzrtot(ixi:ixf:iox) = xznr(1:areas%nx)   
  dvxtot(ixi:ixf:iox) = dvx(1:areas%nx)    
  ugrtot(ixi:ixf:iox) = ugridx(1:areas%nx) 

  if (iswp .eq. 1) then
     dflxtot(1,qy_s:j_end,qz_s:k_end,iarea)=dflxl(qy_s:j_end,qz_s:k_end)
     dflxtot(2,qy_s:j_end,qz_s:k_end,iarea)=dflxr(qy_s:j_end,qz_s:k_end)
     eflxtot(1,qy_s:j_end,qz_s:k_end,iarea)=eflxl(qy_s:j_end,qz_s:k_end)
     eflxtot(2,qy_s:j_end,qz_s:k_end,iarea)=eflxr(qy_s:j_end,qz_s:k_end)
     vxflxtot(1,qy_s:j_end,qz_s:k_end,iarea)=vxflxl(qy_s:j_end,qz_s:k_end)
     vxflxtot(2,qy_s:j_end,qz_s:k_end,iarea)=vxflxr(qy_s:j_end,qz_s:k_end)
     vyflxtot(1,qy_s:j_end,qz_s:k_end,iarea)=vyflxl(qy_s:j_end,qz_s:k_end)
     vyflxtot(2,qy_s:j_end,qz_s:k_end,iarea)=vyflxr(qy_s:j_end,qz_s:k_end)
     vzflxtot(1,qy_s:j_end,qz_s:k_end,iarea)=vzflxl(qy_s:j_end,qz_s:k_end)
     vzflxtot(2,qy_s:j_end,qz_s:k_end,iarea)=vzflxr(qy_s:j_end,qz_s:k_end)
     xnflxtot(1,qy_s:j_end,qz_s:k_end,:,iarea)= xnflxl(qy_s:j_end,qz_s:k_end,:)
     xnflxtot(2,qy_s:j_end,qz_s:k_end,:,iarea)= xnflxr(qy_s:j_end,qz_s:k_end,:)
  endif

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif

  return

end subroutine hydro_are

!======================== Detect shocks
subroutine detect_sh_are(iarea, selftime, childrentime)
  use precision

  use intgrs_hy
  use totare_hy
!  use arecon_hy
  use bndinf_hy
  use massio_hy
  use cputim
  use vnew_hy
  use vold_hy
  use mesh_hy
  use gfloat_hy
  use vnuw_hy
  use hlpare_hy

  use mo_mpi
  use setbnd_mod
  use detect_sh_mod

  use hydro_areas_mod
  use configure

  use state
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk)                :: selftime_start(2)

  integer(kind=ik)             :: i, j, k, ii, jj, kk, isd, iii, &
                                   jji, kki, iif, jjf, kkf
  real(kind=rk)                :: dt_old, dtnew, dvolx, dvoly, dvolz

  integer(kind=ik), intent(in) :: iarea
  integer(kind=ik)             :: j_end, k_end

  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

!-----------------------------------------------------------------------
!     get control indices:
!-----------------------------------------------------------------------

  are_id = iarea
  ixi    = areas%ix_are(iarea, 1)
  ixf    = areas%ix_are(iarea, 2)
  iox    = areas%ix_are(iarea, 3)
  iyi    = areas%ix_are(iarea, 4)
  iyf    = areas%ix_are(iarea, 5)
  ioy    = areas%ix_are(iarea, 6)
  izi    = areas%ix_are(iarea, 7)
  izf    = areas%ix_are(iarea, 8)
  ioz    = areas%ix_are(iarea, 9)
  isd    = areas%ix_are(iarea,10)
!      areas%ix_are(iarea,11) = 1
  config%bndmnx = areas%ix_are(iarea,12)
  config%bndmxx = areas%ix_are(iarea,13)

!-----------------------------------------------------------------------
!     the current time step is the time step of the area:
!-----------------------------------------------------------------------

  hydro%dt     = areas%dt_are(iarea)
  time          = areas%ti_are(iarea)    ! Attention: qterms and grdvel
                                ! need absolut time (end of time step)! 
  dt_old = hydro%dt

!-BN2D call grdvel               ! set grid velocity -> ugrtot(1..config%qx)
! removed since time is never used in setbnd
!  call setbnd(time)         ! set inner boundaries
   call setbnd               ! set inner boundaries


!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------

  areas%nz = izf - izi
  areas%ny = iyf - iyi
  areas%nx = ixf - ixi

  areas%nz = areas%nz/ioz + 1
  areas%ny = areas%ny/ioy + 1
  areas%nx = areas%nx/iox + 1


!-----------------------------------------------------------------------
!     copy a part of the total calculation area in the PPM-arrays:
!-----------------------------------------------------------------------

  if (areas%ny .eq. 1 .and. areas%nz .eq. 1) then
     j_end=qy_s
     k_end=qz_s
  else
     j_end=qy_e
     k_end=qz_e
  end if

  velx  (1:areas%nx,qy_s:j_end,qz_s:k_end) = vextot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  vely  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veytot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  velz  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veztot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  densty(1:areas%nx,qy_s:j_end,qz_s:k_end) = dentot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  energy(1:areas%nx,qy_s:j_end,qz_s:k_end) = enetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  press (1:areas%nx,qy_s:j_end,qz_s:k_end) = pretot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  gammae(1:areas%nx,qy_s:j_end,qz_s:k_end) = gaetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  gammac(1:areas%nx,qy_s:j_end,qz_s:k_end) = gactot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  stot  (1:areas%nx,qy_s:j_end,qz_s:k_end) = stotot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  temp  (1:areas%nx,qy_s:j_end,qz_s:k_end) = temtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
  xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn) = xnutot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end,1:config%qn)


  do k = qz_s, k_end
     do j = qy_s-1, j_end
        ii = ixi-1
        do i = 0, areas%nx
           gpot  (i,j,k) = gpotot(ii,j,k)
           gold  (i,j,k) = 1.e+33_rk
           ii = ii + iox
        enddo
     enddo
  enddo


!-----------------------------------------------------------------------
!     Attention: For the z-direction igeomy = 5 (spherical, theta, phi) is
!     assumed! Since the resolution (ioz) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if (ioz .eq. 1) then 
     zznl(1:areas%nz) = zzltot(izi:izf:ioz)
     zzn(1:areas%nz)  = zzntot(izi:izf:ioz)
     zznr(1:areas%nz) = zzrtot(izi:izf:ioz)
     dvz(1:areas%nz)  = dvztot(izi:izf:ioz)
  else
#ifdef DEBUG
     write(*,*) 'detect_sh_are:izi, nz,ioz ',izi, areas%nz,ioz
#endif
     ii =  izi
     do i = 1, areas%nz
        zznl(i) = zzltot(ii)
        ii = ii + ioz
        if (ii .gt. config%qz+4) ii=1 !quick fix: ii could otherwise 
        !easily exceed q if ioz=qz
     enddo

     do  i = 1, areas%nz
        zznr(i) = zznl(i+1)
     enddo

     do i = 1, areas%nz
        dvz(i) = zznr(i) - zznl(i) !igeomz = 5 is assumed!
        zzn(i) = 0.5_rk * (zznr(i) + zznl(i))
     enddo
  endif

!-----------------------------------------------------------------------
!     Attention: For the y-direction igeomy = 4 (spherical, theta) is
!     assumed! Since the resolution (ioy) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if(ioy .eq. 1) then
     yznl(1:areas%ny)   = yzltot(iyi:iyf:ioy)
     yzn(1:areas%ny)    = yzntot(iyi:iyf:ioy)
     yznr(1:areas%ny)   = yzrtot(iyi:iyf:ioy)
     dvy(1:areas%ny)    = dvytot(iyi:iyf:ioy)
  else
     ii =  iyi
     do i = 1, areas%ny + 8
        yznl(i) = yzltot(ii)
        ii = ii + ioy
        if (ii .gt. config%qy+4) ii=1 !quick fix
     enddo

     yznr(1) = 0.
     do  i = 1, areas%ny + 7
        yznr(i) = yznl(i+1)
     enddo

     do i = 1, areas%ny + 7
        dvy(i) = cos(yznl(i)) - cos(yznr(i)) !igeomy = 4 is assumed!
        yzn(i) = 0.5_rk * (yznr(i) + yznl(i))
     enddo
  endif


  xznl(1:areas%nx)   = xzltot(ixi:ixf:iox)
  xzn(1:areas%nx)    = xzntot(ixi:ixf:iox)
  xznr(1:areas%nx)   = xzrtot(ixi:ixf:iox)
  dvx(1:areas%nx)    = dvxtot(ixi:ixf:iox)
  ugridx(1:areas%nx) = ugrtot(ixi:ixf:iox)

!     -----------------------
!     multi-D shock detection
  call detect_sh

!-----------------------------------------------------------------------
!     copy ishck-array back to the total calculation regime:
!-----------------------------------------------------------------------
  ishtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end) = ishck (1:areas%nx,qy_s:j_end,qz_s:k_end)
  ishtot(ixi:ixf:iox,qy_s-1,qz_s:k_end) = ishck (1:areas%nx,qy_s-1,qz_s:k_end)
  ishtot(ixi:ixf:iox,j_end+1,qz_s:k_end) = ishck (1:areas%nx,j_end+1,qz_s:k_end)
  ishtot(ixi:ixf:iox,qy_s:j_end,qz_s-1) = ishck (1:areas%nx,qy_s:j_end,qz_s-1)
  ishtot(ixi:ixf:iox,qy_s:j_end,k_end+1) = ishck (1:areas%nx,qy_s:j_end,k_end+1)

  if (iarea .eq. 1) then
     ishtot(ixi-1,qy_s:j_end,qz_s:k_end) = ishck (0,qy_s:j_end,qz_s:k_end) 
     ishtot(ixi-1,qy_s-1,qz_s-1) = ishck (0,   qy_s-1,    qz_s-1) 
     ishtot(ixi-1,j_end+1,qz_s-1) = ishck (0,  j_end+1,    qz_s-1) 
     ishtot(ixi-1,qy_s-1,k_end+1) = ishck (0,   qy_s-1, k_end+1) 
     ishtot(ixi-1,j_end+1,k_end+1) = ishck (0,   j_end+1,k_end+1) 
  endif

  if (iarea .eq. areas%are_nu) then
     ishtot(ixf+1,qy_s:j_end,qz_s:k_end) = ishck (areas%nx+1,qy_s:j_end,qz_s:k_end) 
     ishtot(ixf+1,qy_s-1,qz_s-1)           = ishck (areas%nx+1,   qy_s-1,   qz_s-1   ) 
     ishtot(ixf+1,j_end+1,qz_s-1)          = ishck (areas%nx+1,  j_end+1, qz_s-1   ) 
     ishtot(ixf+1,qy_s-1,k_end+1)           = ishck (areas%nx+1,   qy_s-1,   k_end+1) 
     ishtot(ixf+1,j_end+1,k_end+1)          = ishck (areas%nx+1,   j_end+1,k_end+1) 
  endif

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
  return

end subroutine detect_sh_are

!======================== Calculate Hydro ===============================

subroutine accel_are(iarea, selftime, childrentime)

!-----------------------------------------
! Autor           : Markus Rampp (MPA) 
! Modul           : $Id: rady.F,v 1.14 2005/04/07 19:29:44 rburas Exp $
! Version         : $Revision: 1.14 $
! Date            : $Date: 2005/04/07 19:29:44 $
!
!     task        :  get 3D calculation area and calculate it 
!                    separately 
!
!     input:   iarea  =  area ID
!              dta = time step [s] of the area
!              isd = sweep direction -> 1 = xyz, 2 = zyx
!              i*i = initial index of area in the total area
!              i*f = final index of area in the total area
!              io* = index offset -> full resolution = 1, etc.
!
!=======================================================================
  use precision


  use intgrs_hy
  use totare_hy
  use totgrq_hy
!  use arecon_hy
  use bndinf_hy
  use nutrio_hy
  use massio_hy
  use cputim
  use vnew_hy
  use vold_hy
  use mesh_hy
  use gfloat_hy
  use spez_hy
  use vnuw_hy
  use hlpare_hy
#ifndef NOTRA
  use qave_overload
#endif

#ifdef HTCL
  use htcl_hy
#endif

  use sweeps_mod
  use mo_mpi

  use hydro_areas_mod
  use configure
  use cputim
  use state
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2), sweep_self(2), &
                                sweep_children(2)
  integer(kind=ik) :: i, j, k, ii, jj, kk, isd, iii, jji, kki, iif, jjf, kkf
  real(kind=rk)    :: dt_old, dtnew, dvolx, dvoly, dvolz

  integer(kind=ik), intent(in) :: iarea
  real(kind=rk), dimension(2)  :: tim1, tim2

  integer(kind=ik)             :: j_end,k_end

  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif
!-----------------------------------------------------------------------
!     get control indices:
!-----------------------------------------------------------------------

  are_id = iarea
  ixi    = areas%ix_are(iarea, 1)
  ixf    = areas%ix_are(iarea, 2)
  iox    = areas%ix_are(iarea, 3)
  iyi    = areas%ix_are(iarea, 4)
  iyf    = areas%ix_are(iarea, 5)
  ioy    = areas%ix_are(iarea, 6)
  izi    = areas%ix_are(iarea, 7)
  izf    = areas%ix_are(iarea, 8)
  ioz    = areas%ix_are(iarea, 9)
  isd    = areas%ix_are(iarea,10)
!      ix_are(iarea,11) = 1
  config%bndmnx = areas%ix_are(iarea,12)
  config%bndmxx = areas%ix_are(iarea,13)

!-----------------------------------------------------------------------
!     the current time step is the time step of the area:
!-----------------------------------------------------------------------

  hydro%dt     = areas%dt_are(iarea)
  time          = areas%ti_are(iarea)    ! Attention: qterms and grdvel
                                ! need absolut time (end of time step)! 
  dt_old = hydro%dt

!      call setbnd(time)         ! set inner boundaries

!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------

  areas%nz = izf - izi
  areas%ny = iyf - iyi
  areas%nx = ixf - ixi

  areas%nz = areas%nz/ioz + 1
  areas%ny = areas%ny/ioy + 1
  areas%nx = areas%nx/iox + 1

  if (ioy .eq. 1 .and. ioz .eq. 1) then
     j_end=qy_e
     k_end=qz_e
  else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     j_end=qy_s
     k_end=qz_s
  end if

#ifdef HTCL
    !> determine the optical depth for heating and cooling model
    do k = qz_s,qz_e
      do j = qy_s,qy_e
        f_abs(:,j,k) = 0._rk
        do i = config%qx-1,1,-1
          !> determine opacity according to Janka(2001) (16)
          k_abs = 1.5e-17_rk * dentot(i,j,k) * 0.5 * (xnutot(i,j,k,1) + xnutot(i,j,k,2))
          !> now the optical depth
          f_abs(i,j,k) = f_abs(i+1,j,k) + ( xzrtot(i)-xzltot(i) ) * k_abs
        end do
      end do
    end do
    f_abs_a(1:areas%nx,qy_s:j_end,qz_s:k_end) = f_abs(ixi:ixf,qy_s:j_end,qz_s:k_end)
#endif /* HTCL */

!-----------------------------------------------------------------------
!     Attention: "poisson" in "sweep" needs the correct enclosed mass
!-----------------------------------------------------------------------
!      pmass = tmatot(0)  ! enclosed mass at the inner boundary
!      pmass = tmatot(ixi - 1)  ! enclosed mass at the inner boundary


!-----------------------------------------------------------------------
!     copy a part of the total calculation area in the PPM-arrays:
!-----------------------------------------------------------------------

  velx  (1:areas%nx,qy_s:j_end,qz_s:k_end) = vextot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  vely  (1:areas%nx,qy_s:j_end,qz_s:k_end) =  veytot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  velz  (1:areas%nx,qy_s:j_end,qz_s:k_end) =  veztot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  densty(1:areas%nx,qy_s:j_end,qz_s:k_end) = dentot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  energy(1:areas%nx,qy_s:j_end,qz_s:k_end) = enetot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  press (1:areas%nx,qy_s:j_end,qz_s:k_end) =  pretot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  gammae(1:areas%nx,qy_s:j_end,qz_s:k_end) =  gaetot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  gammac(1:areas%nx,qy_s:j_end,qz_s:k_end) =  gactot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  stot  (1:areas%nx,qy_s:j_end,qz_s:k_end) =  stotot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  temp  (1:areas%nx,qy_s:j_end,qz_s:k_end) =  temtot(ixi:ixf,qy_s:j_end,qz_s:k_end)
  xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn) = xnutot(ixi:ixf,qy_s:j_end,qz_s:k_end,1:config%qn)
  cpo   (1:areas%nx,qy_s:j_end,qz_s:k_end,1:4) = cpotot(ixi:ixf,qy_s:j_end,qz_s:k_end,1:4)

  vxvold(1:areas%nx,qy_s:j_end,qz_s:k_end) =  vexold(ixi:ixf,qy_s:j_end,qz_s:k_end)
  vyvold(1:areas%nx,qy_s:j_end,qz_s:k_end) =  veyold(ixi:ixf,qy_s:j_end,qz_s:k_end) 
  vzvold(1:areas%nx,qy_s:j_end,qz_s:k_end) =  vezold(ixi:ixf,qy_s:j_end,qz_s:k_end) 
  denold(1:areas%nx,qy_s:j_end,qz_s:k_end) =  dnsold(ixi:ixf,qy_s:j_end,qz_s:k_end)

  do k = qz_s,k_end
     do j=qy_s-1,j_end
        ii = ixi-1
        do i = 0, areas%nx
           gpot  (i,j,k) = gpotot(i+ii,j,k)
           gold  (i,j,k) = gpoold(i+ii,j,k)
        enddo
     enddo
  enddo

  qgrv(1:areas%nx,qy_s:j_end,qz_s:k_end)=0._rk

!-----------------------------------------------------------------------
!     Attention: For the z-direction igeomy = 5 (spherical, theta, phi) is
!     assumed! Since the resolution (ioz) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if (ioz .eq. 1) then 
     zznl(1:areas%nz) = zzltot(izi:izf:ioz)
     zzn(1:areas%nz)  = zzntot(izi:izf:ioz)
     zznr(1:areas%nz) = zzrtot(izi:izf:ioz)
     dvz(1:areas%nz)  = dvztot(izi:izf:ioz)
  else
#ifdef DEBUG
     write(*,*) 'detect_sh_are:izi, nz,ioz ',izi, areas%nz,ioz
#endif
     ii =  izi
     do i = 1, areas%nz
        zznl(i) = zzltot(ii)
        ii = ii + ioz
        if (ii .gt. config%qz+4) ii=1 !quick fix: ii could otherwise 
        !easily exceed q if ioz=qz
     enddo

     do  i = 1, areas%nz
        zznr(i) = zznl(i+1)
     enddo

     do i = 1, areas%nz
        dvz(i) = zznr(i) - zznl(i) !igeomz = 5 is assumed!
        zzn(i) = 0.5_rk * (zznr(i) + zznl(i))
     enddo
  endif


!-----------------------------------------------------------------------
!     Attention: For the y-direction igeomy = 4 (spherical, theta) is
!     assumed! Since the resolution (ioy) may be greater than 1 
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------

  if(ioy .eq. 1) then
     yznl(1:areas%ny)   = yzltot(iyi:iyf:ioy)
     yzn(1:areas%ny)    = yzntot(iyi:iyf:ioy)
     yznr(1:areas%ny)   = yzrtot(iyi:iyf:ioy)
     dvy(1:areas%ny)    = dvytot(iyi:iyf:ioy)
  else
     ii =  iyi
     do i = 1, areas%ny + 8
        yznl(i) = yzltot(ii)
        ii = ii + ioy
        if (ii .gt. config%qy+4) ii=1 !quick fix
     enddo

     yznr(1) = 0._rk
     do  i = 1, areas%ny + 7
        yznr(i) = yznl(i+1)
     enddo

     do i = 1, areas%ny + 7
        dvy(i) = cos(yznl(i)) - cos(yznr(i)) !igeomy = 4 is assumed!
        yzn(i) = 0.5_rk * (yznr(i) + yznl(i))
     enddo
  endif


!-----------------------------------------------------------------------
!     Attention: For the x-direction the resolution (iox) must be 1!
!     Otherwise the center of zone and boundaries of the zone are
!     different!!!
!-----------------------------------------------------------------------

!      xznl(1:areas%nx)   = xzltot(ixi:ixf:iox)
!      xzn(1:areas%nx)    = xzntot(ixi:ixf:iox)
!      xznr(1:areas%nx)   = xzrtot(ixi:ixf:iox)
!      dvx(1:areas%nx)    = dvxtot(ixi:ixf:iox)
!      ugridx(1:areas%nx) = ugrtot(ixi:ixf:iox)

  xznl(1:areas%nx)   = xzltot(ixi:ixf)
  xzn(1:areas%nx)    = xzntot(ixi:ixf)
  xznr(1:areas%nx)   = xzrtot(ixi:ixf)
  dvx(1:areas%nx)    = dvxtot(ixi:ixf)
  ugridx(1:areas%nx) = ugrtot(ixi:ixf)

!-- due to nonlocal transport effects qyetot,qentot,qmotot  vary
!    within the unresolved hydro areas; therefore averages are 
!    appropriate;

! Beware here we call the overloaded subroutine qave, whcih is either
! the purely OPENMP qave_are_OPENMP or the hybrid MPI/OPENMP 
! qave_are_MPI_HYDRO version. The overloading is done in module qave_overload
#ifndef NOTRA
 if (config%p_ntr .ne. 0) call qave_are(iarea)
#endif

!-----------------------------------------------------------------------
!     call x/y/z- or z/y/x-sweeps according to "isd":
!-----------------------------------------------------------------------

  if(isd .eq. 1) then  !-------- forward cycle
     timer%hydro_sweep_mode = 1
     call sweep (1,.false.,.true., 2, sweep_self, sweep_children)  !  x - sweep
     timer%accelare_sweep(timer%hydro_sweep_mode,:,1)=timer%accelare_sweep(timer%hydro_sweep_mode,:,1)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,1)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,1)+sweep_children

     if (config%nsdim .ge. 2 .and. areas%ny .ge. 4)  then  
        !  y - sweep  (only 2D and 3D)

        timer%hydro_sweep_mode = 1

        call sweep (2,.false.,.true., 2, sweep_self, sweep_children)  
        timer%accelare_sweep(timer%hydro_sweep_mode,:,2)=timer%accelare_sweep(timer%hydro_sweep_mode,:,2)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,2)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,2)+sweep_children
  end if

     if (config%nsdim .eq. 3 .and. areas%ny .ge. 4 .and. areas%nz .ge. 4)  then  !  z - sweep  (only 3D)
        timer%hydro_sweep_mode = 1

        call sweep (3,.false.,.true., 2, sweep_self, sweep_children)  
        timer%accelare_sweep(timer%hydro_sweep_mode,:,3)=timer%accelare_sweep(timer%hydro_sweep_mode,:,3)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,3)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,3)+sweep_children
     end if

  else                 !-------- backward cycle

     if (config%nsdim .eq. 3 .and. areas%ny .ge. 4 .and. areas%nz .ge. 4)  then  !  z - sweep  (only 3D)

        timer%hydro_sweep_mode = 1

        call sweep (3,.false.,.true., 2, sweep_self, sweep_children)  
        timer%accelare_sweep(timer%hydro_sweep_mode,:,3)=timer%accelare_sweep(timer%hydro_sweep_mode,:,3)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,3)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,3)+sweep_children

     end if

     if (config%nsdim .ge. 2 .and. areas%ny .ge. 4)  then
        !  y - sweep  (only 2D and 3D)

        timer%hydro_sweep_mode = 1

        call sweep (2,.false.,.true., 2, sweep_self, sweep_children)  
        timer%accelare_sweep(timer%hydro_sweep_mode,:,2)=timer%accelare_sweep(timer%hydro_sweep_mode,:,2)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,2)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,2)+sweep_children
     end if

     timer%hydro_sweep_mode = 1

     call sweep (1,.false.,.true., 2, sweep_self, sweep_children)  
        timer%accelare_sweep(timer%hydro_sweep_mode,:,1)=timer%accelare_sweep(timer%hydro_sweep_mode,:,1)+sweep_self
     childrentime = childrentime + sweep_self
     timer%accelare_sweep_children(timer%hydro_sweep_mode,:,1)= &
        timer%accelare_sweep_children(timer%hydro_sweep_mode,:,1)+sweep_children
  end if

!-----------------------------------------------------------------------
!     copy PPM-arrays back to the total calculation regime:
!-----------------------------------------------------------------------

  vextot(ixi:ixf,qy_s:j_end,qz_s:k_end) = velx  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  veytot(ixi:ixf,qy_s:j_end,qz_s:k_end) = vely  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  veztot(ixi:ixf,qy_s:j_end,qz_s:k_end) = velz  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  dentot(ixi:ixf,qy_s:j_end,qz_s:k_end) = densty(1:areas%nx,qy_s:j_end,qz_s:k_end)
  enetot(ixi:ixf,qy_s:j_end,qz_s:k_end) = energy(1:areas%nx,qy_s:j_end,qz_s:k_end)
  pretot(ixi:ixf,qy_s:j_end,qz_s:k_end) = press (1:areas%nx,qy_s:j_end,qz_s:k_end)
  gaetot(ixi:ixf,qy_s:j_end,qz_s:k_end) = gammae(1:areas%nx,qy_s:j_end,qz_s:k_end)
  gactot(ixi:ixf,qy_s:j_end,qz_s:k_end) = gammac(1:areas%nx,qy_s:j_end,qz_s:k_end)
  stotot(ixi:ixf,qy_s:j_end,qz_s:k_end) = stot  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  temtot(ixi:ixf,qy_s:j_end,qz_s:k_end) = temp  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  xnutot(ixi:ixf,qy_s:j_end,qz_s:k_end,1:config%qn) =xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn)
  cpotot(ixi:ixf,qy_s:j_end,qz_s:k_end,1:4) = cpo  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:4)

  qgrtot(ixi:ixf,qy_s:j_end,qz_s:k_end)   = qgrv (1:areas%nx,qy_s:j_end,qz_s:k_end)

  qentot(ixi:ixf,qy_s:j_end,qz_s:k_end)   = qen  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  qyetot(ixi:ixf,qy_s:j_end,qz_s:k_end,1) = qye  (1:areas%nx,qy_s:j_end,qz_s:k_end,1)
  qyetot(ixi:ixf,qy_s:j_end,qz_s:k_end,2) = qye  (1:areas%nx,qy_s:j_end,qz_s:k_end,2)
  qyetot(ixi:ixf,qy_s:j_end,qz_s:k_end,3) = qye  (1:areas%nx,qy_s:j_end,qz_s:k_end,3)
  qyetot(ixi:ixf,qy_s:j_end,qz_s:k_end,4) = qye  (1:areas%nx,qy_s:j_end,qz_s:k_end,4)
  qyetot(ixi:ixf,qy_s:j_end,qz_s:k_end,5) = qye  (1:areas%nx,qy_s:j_end,qz_s:k_end,5)
  qmotot(ixi:ixf,qy_s:j_end,qz_s:k_end)   = qmo  (1:areas%nx,qy_s:j_end,qz_s:k_end)
  qmytot(ixi:ixf,qy_s:j_end,qz_s:k_end)   = qmy  (1:areas%nx,qy_s:j_end,qz_s:k_end)

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
end subroutine accel_are

end module hydro_area_functions

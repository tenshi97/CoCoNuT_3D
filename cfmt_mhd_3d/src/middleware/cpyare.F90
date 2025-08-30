#undef DEBUG

module cpyare_mod

implicit none

contains

!=======================================================================
subroutine cpyare(isw)
!=======================================================================
!
!     copy area arrays completely -> calare test
!
!     task:    get or put 3D calculation area
!
!     input:   isw = 0: copy velx -> vextot
!                    1: copy velx -> vextot + set velx = -1.e33
!                    2: copy vextot -> velx
!
!=======================================================================
  use precision
  
  use totare_hy
  use nutrio_hy, only : cpotot
  
  use vnew_hy
  use mesh_hy

  use mo_mpi

  
  use configure
                                                                       
!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------


!      write(*,'(''cpyare-test> isw = '',i1,'' config%qx = '',i3,
!     &          '' config%qy = '',i3,'' config%qz = '',i3,'' config%qn = '',i2)')
!     &         isw, config%qx, config%qy, config%qz, config%qn

!-----------------------------------------------------------------------
!     copy a part of the total calculation area in the PPM-arrays:
!-----------------------------------------------------------------------
  implicit none
! LOCAL variables that are not in modules

  integer(kind=ik), intent(in) :: isw

  if(isw .eq. 2) then

! MPI_HYDRO 1:config%qy <- qy_s, qy_e
     velx  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     vely  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     velz  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = veztot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     densty(1:config%qx,qy_s:qy_e,qz_s:qz_e) = dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     energy(1:config%qx,qy_s:qy_e,qz_s:qz_e) = enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     press (1:config%qx,qy_s:qy_e,qz_s:qz_e) = pretot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     gammae(1:config%qx,qy_s:qy_e,qz_s:qz_e) = gaetot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     gammac(1:config%qx,qy_s:qy_e,qz_s:qz_e) = gactot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     stot  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = stotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     temp  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = temtot(1:config%qx,qy_s:qy_e,qz_s:qz_e)

! MPI_HYDRO: 0:config%qy instead of qy_s-1:qy_e
     gpot (0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e)
     xnuc  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn) = xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn)
     cpo  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:4) = cpotot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:4)
     
     zznl(:) = zzltot(:)
     zzn(:)  = zzntot(:)
     zznr(:) = zzrtot(:)
     dvz(:)  = dvztot(:)
 
     yznl(:) = yzltot(:)
     yzn(:)  = yzntot(:)
     yznr(:) = yzrtot(:)
     dvy(:)  = dvytot(:)

     xznl(:) = xzltot(:)
     xzn(:)  = xzntot(:)
     xznr(:) = xzrtot(:)
     dvx(:)  = dvxtot(:)
 

!-----------------------------------------------------------------------
!     copy PPM-arrays back to the total calculation regime:
!-----------------------------------------------------------------------
  else



! MPI_HYDRO: 1:config%qy instead of qy_s:qy:e
     vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = velx  (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = vely  (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     veztot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = velz  (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = densty(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = energy(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     pretot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = press (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     gaetot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = gammae(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     gactot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = gammac(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     stotot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = stot  (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     temtot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = temp  (1:config%qx,qy_s:qy_e,qz_s:qz_e)
     gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = gpot  (0:config%qx,qy_s-1:qy_e,qz_s:qz_e)
     xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn) = xnuc  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn)
     cpotot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:4) = cpo  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:4)

     zzltot(:) = zznl(:)
     zzntot(:) = zzn(:) 
     zzrtot(:) = zznr(:)
     dvztot(:) = dvz(:) 

     yzltot(:) = yznl(:)
     yzntot(:) = yzn(:) 
     yzrtot(:) = yznr(:)
     dvytot(:) = dvy(:)

     xzltot(:) = xznl(:)
     xzntot(:) = xzn(:) 
     xzrtot(:) = xznr(:)
     dvxtot(:) = dvx(:) 
 

!-----------------------------------------------------------------------
!     for test overwrite all PPM-arrays:
!-----------------------------------------------------------------------

     if(isw .ne. 0) then 

           
        velx  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        vely  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        velz  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        densty(1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        energy(1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        press (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        gammae(1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        gammac(1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        stot  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        temp  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -1.0e33_rk
        gpot  (0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = -1.0e33_rk
        xnuc  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn) = -1.0e33_rk
        cpo  (1:config%qx,qy_s:qy_e,qz_s:qz_e,1:4) = -1.0e33_rk
     end if

  end if


  return
end subroutine cpyare

end module cpyare_mod

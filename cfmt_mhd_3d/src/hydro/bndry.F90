module bndry_mod

implicit none

contains

!> \par calculate boundary conditions
!>
!> \author W. Keil
!>
!> \param j type of boundary
!> \param k phi-index
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
subroutine bndry(j, k)

  use precision

  use intgrs_hy
  use bndinf_hy

  use vold_hy
  use mesh_hy
  use gfloat_hy
  use lfloat_hy
  use hydro_hy
  use grd_hy      
  use physcs_hy     

  use phycon, only : pc_pi

  use abort

  use vnew_hy
  use hlpare_hy

  use mo_mpi

  use hydro_areas_mod
  use configure
  implicit none

! LOCAL variables that are not in modules
  integer(kind=ik), intent(in) :: j, k
  integer(kind=ik)             :: ibndbt, i, ibi, nnn, ibndtp, ibndlt
  integer(kind=ik)             :: ibndrt, ibndmn, ibndmx, n, j1, j2, k1, k2
  integer(kind=ik)             :: nzn54, i1, i2, nzn76

  real(kind=rk)                :: yrbot, ylbot, xrbot, xlbot, yltop, yrtop
  real(kind=rk)                :: xltop, xrtop, zrlft, zllft, yrlft, yllft
  real(kind=rk)                :: zlrgt, zrrgt, ylrgt, yrrgt

  integer(kind=ik)             :: j_end, k_end
  
  j1 = qy_s
  j2 = qy_e
  k1 = qz_s
  k2 = qz_e

  if (ioy .eq. 1 .and. ioz .eq. 1) then
     j_end=qy_e
  else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     j_end=qy_s
  else
     raise_abort("Invalid combinations for ioy, ioz")
  end if

  if (ioy .eq. 1 .and. ioz .eq. 1) then
     k_end=qz_e
  else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     k_end=qz_s
  else
     raise_abort("Invalid combinations for ioy, ioz")
  end if

!
!     first reset values of  utbt, uttp, utlt, utrt, xbot, xtop,
!                            ybot, ytop, ylft, yrgt, zlft, zrgt      
!
  ibndbt = bndbot
  ibndtp = bndtop
!
!     bottom of grid
!
  if (xyzswp .eq. 1) then
     nnn = qy_s !theta-direction is top-bottom-direction
  else 
     nnn = 1 !r-direction is top-bottom-direction
  endif
  if (j .ne. nnn) go to 100
!
  go to (10, 20, 30, 40, 50, 60, 70), ibndbt
      !
!     reflecting boundary
!
10 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        utbt(i+4) = -vyold(i,qy_s,k)
     enddo
!
     yrbot = yznl(qy_s)
     ylbot = 2.0_rk  * yznl(qy_s) - yznr(qy_s)
     ybot  = 0.5_rk * (yrbot  + ylbot)
!
  else if (xyzswp .eq. 2) then
!
     do i = qy_s, j_end 
        utbt(i+4-config%q_yoff) = -vxold(1,i,k)
     enddo
!
     xrbot = xznl(1)
     xlbot = 2.0_rk  * xznl(1) - xznr(1)
     xbot  = 0.5_rk * (xrbot  + xlbot)
!
  else if (xyzswp .eq. 3) then
!
     do i = qz_s, k_end
        utbt(i+4-config%q_zoff) = -vxold(1,k,i)
     enddo
!
     xrbot = xznl(1)
     xlbot = 2.0_rk  * xznl(1) - xznr(1)
     xbot  = 0.5_rk * (xrbot  + xlbot)
     !
  endif
!
  go to 100
!
!     flow out boundary
!
20 continue
!
!     value of utbt already correct
!
  if (xyzswp .eq. 1) then
!
     yrbot = yznl(1)
     ylbot = 2.0_rk  * yznl(1) - yznr(1)
     ybot  = 0.5_rk * (yrbot  + ylbot)
!
  else
!
     xrbot =  xznl(1)
     xlbot = 2.0_rk  * xznl(1) - xznr(1)
     xbot  = 0.5_rk * (xrbot  + xlbot)
     !
  endif
  !
  go to 100
!
  !     flow in boundary
!
30 continue
!
  do i = 5, nzn4
     utbt(i) = utin
  enddo
!
  if (xyzswp .eq. 1) then
     !
     yrbot = yznl(1)
     ylbot = 2.0_rk  * yznl(1) - yznr(1)
     ybot  = 0.5_rk * (yrbot  + ylbot)
     !
  else
!
     xrbot =  xznl(1)
     xlbot = 2.0_rk  * xznl(1) - xznr(1)
     xbot  = 0.5_rk * (xrbot  + xlbot)
     !
  endif
  !
  go to 100
!
!     periodic boundary
!
40 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        utbt(i+4) = vyold(i,j_end,k)
     enddo
!
     yrbot = yznl(qy_s)

     ylbot = yznl(qy_s) - yznr(j_end) + yznl(j_end)
     ybot  = 0.5_rk * (yrbot + ylbot)
!
  else if (xyzswp .eq. 2) then
!

     do i = qy_s, j_end 
        utbt(i+4-config%q_yoff) = vxold(areas%nx,i,k)

     enddo
!
     xrbot =  xznl(1)
     xlbot =  xznl(1) - xznr(areas%nx) + xznl(areas%nx)
     xbot  = 0.5_rk * (xrbot + xlbot)
     !
  else if (xyzswp .eq. 3) then
!

     do i = qz_s, k_end 
        utbt(i+4-config%q_zoff) = vxold(areas%nx,k,i)
     enddo
     !
     xrbot =  xznl(1)
     xlbot =  xznl(1) - xznr(areas%nx) + xznl(areas%nx)
     xbot  = 0.5_rk * (xrbot + xlbot)
     !
  endif
!
  go to 100
!-----------------------------------------------------------------------
!     special boundary
!-----------------------------------------------------------------------
50 continue

  if (xyzswp .eq. 1) then
     
     do i = 1, areas%nx
        utbt(i+4) = -vyold(i,1,k)
     enddo

     yrbot = yznl(1)
     ylbot = 2.0_rk  * yznl(1) - yznr(1)
     ybot  = 0.5_rk * (yrbot  + ylbot)
     
  else if (xyzswp .eq. 2) then

     ibi = (are_id - 1)*8 + 4   ! take "left" adjacent fields


     do i = qy_s, j_end
        utbt(i+4-config%q_yoff) = vex_bi(ibi, i, k)
     enddo

     xrbot = xznl(1)
     xlbot = xznl(1) - dxx_bi(ibi)
     xbot  = 0.5_rk * (xrbot  + xlbot)

     
  else if (xyzswp .eq. 3) then

     ibi = (are_id - 1)*8 + 4   ! take "left" adjacent fields

     do i = qz_s, k_end 
        utbt(i+4-config%q_zoff) = vex_bi(ibi, k , i)
     end do

     xrbot = xznl(1)
     xlbot = xznl(1) - dxx_bi(ibi)
     xbot  = 0.5_rk * (xrbot  + xlbot)

  endif
!
  go to 100
!-----------------------------------------------------------------------
!     MPI domain boundary  (theta-direction)
!-----------------------------------------------------------------------
60 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        utbt(i+4) = vyold(i,qy_s-1,k)
     enddo
!
     yrbot = yznr(qy_s-1)
     ylbot = yznl(qy_s-1)
     ybot  = 0.5_rk * (yrbot + ylbot)
  else 
     raise_abort("bndry(): error: bottom bndry")
  endif
!
  go to 100
!-----------------------------------------------------------------------
!     MPI domain boundary  (phi-direction)
!-----------------------------------------------------------------------
70 continue
!
     raise_abort("bndry(): error: bottom bndry")
!
!
!
100 continue
!
!
!=======================================================================
!     top of grid
!=======================================================================
!
     
  if (xyzswp .eq. 1) then
     nnn = qy_e
  else 
     nnn = areas%nx
  endif
  if (j .ne. nnn) go to 200
!
  go to (110, 120, 130, 140, 150, 160, 170), ibndtp
!
!     reflecting boundary
!
110 continue
!
  if (xyzswp .eq. 1) then
     !
     do i = 1, areas%nx
        uttp(i+4) = -vyold(i,j_end,k)
     enddo
     !
     yltop = yznr(j_end)
     yrtop = 2.0_rk  * yznr(j_end) - yznl(j_end)
     ytop  = 0.5_rk * (yrtop   + yltop)
     !
  else if (xyzswp .eq. 2) then
!
     do i = qy_s, j_end
        uttp(i+4-config%q_yoff) = -vxold(areas%nx,i,k)
     enddo
     !
     xltop =  xznr(areas%nx)
     xrtop = 2.0_rk * xznr(areas%nx) - xznl(areas%nx)
     xtop  = 0.5_rk * (xrtop + xltop)
     !
  else if (xyzswp .eq. 3) then
!
     do i = qz_s, k_end
        uttp(i+4-config%q_zoff) = -vxold(areas%nx,k,i)
     enddo

     !
     xltop =  xznr(areas%nx)
     xrtop = 2.0_rk * xznr(areas%nx) - xznl(areas%nx)
     xtop  = 0.5_rk * (xrtop + xltop)
     !
  endif
!
  go to 200
!
!     flow out boundary
  !
120 continue
!
!     value of uttp already correct
!
      if (xyzswp .eq. 1) then
!
         yltop = yznr(j_end)
         yrtop = 2.0_rk * yznr(j_end) - yznl(j_end)
         ytop  = 0.5_rk * (yrtop + yltop)
!
      else
!
         xltop =  xznr(areas%nx)
         xrtop = 2.0_rk  * xznr(areas%nx) - xznl(areas%nx)
         xtop  = 0.5_rk * (xrtop   + xltop)
!
      endif
!
      go to 200
!
!     flow in boundary
!
130   continue
!
      do i = 5, nzn4
         uttp(i) = utin
      enddo
!
      if (xyzswp .eq. 1) then
!
         yltop = yznr(j_end)
         yrtop = 2.0_rk  * yznr(j_end) - yznl(j_end)
         ytop  = 0.5_rk * (yrtop   + yltop)
!
      else
!
         xltop =  xznr(areas%nx)
         xrtop = 2.0_rk * xznr(areas%nx) - xznl(areas%nx)
         xtop  = 0.5_rk * (xrtop + xltop)
!
      endif
!
      go to 200
!
!     periodic boundary
!
140   continue
!
      if (xyzswp .eq. 1) then
!

         do i = 1, areas%nx
            uttp(i+4) = vyold(i,1,k)
         enddo

!
         yltop = yznr(j_end)

         yrtop = yznr(j_end) + yznr(qy_s) - yznl(qy_s)

         ytop  = 0.5_rk * (yrtop + yltop)
!
      else if (xyzswp .eq. 2) then
!
         do i = qy_s, j_end
            uttp(i+4-config%q_yoff) = vxold(1,i,k)
         enddo
!
         xltop =  xznr(areas%nx)
         xrtop =  xznr(areas%nx) + xznr(1) -xznl(1)
         xtop  = 0.5_rk * (xrtop + xltop)
!
      else if (xyzswp .eq. 3) then
!
         do i = qz_s, k_end
            uttp(i+4-config%q_zoff) = vxold(1,k,i)
         enddo
!
         xltop =  xznr(areas%nx)
         xrtop =  xznr(areas%nx) + xznr(1) -xznl(1)
         xtop  = 0.5_rk * (xrtop + xltop)
!
      endif
!
      go to 200
!
150    continue
!-----------------------------------------------------------------------
!     special boundary condition
!-----------------------------------------------------------------------

      if (xyzswp .eq. 1) then

         do i = 1, areas%nx
            uttp(i+4) = -vyold(i,j_end,k)
         enddo

         yltop = yznr(j_end)
         yrtop = 2.0_rk  * yznr(j_end) - yznl(j_end)
         ytop  = 0.5_rk * (yrtop   + yltop)

      else if (xyzswp .eq. 2) then

         ibi = are_id*8 + 4 + 1

         do i = qy_s, j_end
            uttp(i+4-config%q_yoff) = vex_bi(ibi,i,k)
         enddo

         xltop =  xznr(areas%nx)
         xrtop =  xznr(areas%nx) - dxx_bi(ibi)
         xtop  = 0.5_rk * (xrtop + xltop)

      else if (xyzswp .eq. 3) then

         ibi = are_id*8 + 4 + 1

         do i = qz_s, k_end 
             uttp(i+4-config%q_zoff) = vex_bi(ibi,k,i)
         enddo

         xltop =  xznr(areas%nx)
         xrtop =  xznr(areas%nx) - dxx_bi(ibi)
         xtop  = 0.5_rk * (xrtop + xltop)

      endif
!
      go to 200
!-----------------------------------------------------------------------
!     MPI domain boundary  (theta-direction)
!-----------------------------------------------------------------------
160 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        uttp(i+4) = vyold(i,qy_e+1,k)
     enddo
!
     yrtop = yznr(qy_e+1)
     yltop = yznl(qy_e+1)
     ytop  = 0.5_rk * (yrtop + yltop)
  else 
     raise_abort("bndry(): error: top bndry")
  endif
!
  go to 200
!-----------------------------------------------------------------------
!     MPI domain boundary  (phi-direction)
!-----------------------------------------------------------------------
170 continue
!
     raise_abort("bndry(): error: top bndry")
!
!
!
200   continue
!
!     now reset values of utlt, utrt, xlft, xrgt, ylft, yrgt, zlft, zrgt
!
      ibndlt = bndlft
      ibndrt = bndrgt
!
!
!=======================================================================
!     left of grid
!=======================================================================
!
      if (xyzswp .eq. 3) then
         nnn = qy_s !theta-direction is left-right-direction
      else 
         nnn = qz_s !phi-direction is left-right-direction
      endif
      if (k .ne. nnn) go to 300
!
      go to (210, 220, 230, 240, 250, 260, 270), ibndlt
!
!     reflecting boundary
!
210   continue
!
      if (xyzswp .eq. 1) then
!

         do i = 1, areas%nx
            utlt(i+4) = -vzold(i,j,qz_s)
         end do

         zrlft = zznl(qz_s)
         zllft = 2.0_rk  * zznl(qz_s) - zznr(qz_s)

         zlft  = 0.5_rk * (zrlft  + zllft)
!
      else if (xyzswp .eq. 2) then
!
        do i = qy_s, j_end
            utlt(i+4-config%q_yoff) = -vzold(j,i,qz_s)
         enddo
!
         zrlft = zznl(qz_s)
         zllft = 2.0_rk  * zznl(qz_s) - zznr(qz_s)

         zlft  = 0.5_rk * (zrlft  + zllft)
!
      else if (xyzswp .eq. 3) then
!
         do i = qz_s, k_end
            utlt(i+4-config%q_zoff) = -vyold(j,qy_s,i)
         enddo

         yrlft = yznl(1)
         yllft = 2.0_rk  * yznl(1) - yznr(1)
         ylft  = 0.5_rk * (yrlft  + yllft)
!
      endif
!
      go to 300
!
!     flow out boundary
!
220   continue
!
!     value of utlt already correct
!
      if (xyzswp .ne. 3) then
!
         zrlft = zznl(1)
         zllft = 2.0_rk  * zznl(1) - zznr(1)
         zlft  = 0.5_rk * (zrlft  + zllft)
!
      else
!
         yrlft =  yznl(1)
         yllft = 2.0_rk  * yznl(1) - yznr(1)
         ylft  = 0.5_rk * (yrlft  + yllft)
!
      endif
!
      go to 300
!
!     flow in boundary
!
230   continue
!
      do i = 5, nzn4
         utlt(i) = uttin
      enddo
!
      if (xyzswp .ne. 3) then
!
         zrlft = zznl(1)
         zllft = 2.0_rk  * zznl(1) - zznr(1)
         zlft  = 0.5_rk * (zrlft  + zllft)
!
      else
!
         yrlft =  yznl(1)
         yllft = 2.0_rk  * yznl(1) - yznr(1)
         ylft  = 0.5_rk * (yrlft  + yllft)
!
      endif
!
      go to 300
!
!     periodic boundary
!
240   continue
!
      if (xyzswp .eq. 1) then
!
         do i = 1, areas%nx

            utlt(i+4) = vzold(i,j,k_end)

         enddo
!
         zrlft = zznl(qz_s)
         zllft = zznl(qz_s) - zznr(k_end) + zznl(k_end)
         zlft  = 0.5_rk * (zrlft + zllft)
!
      else if (xyzswp .eq. 2) then
!
         do i = qy_s, j_end 
            utlt(i+4-config%q_yoff) = vzold(j,i,k_end)
         enddo

         zrlft = zznl(qz_s)
         zllft = zznl(qz_s) - zznr(k_end) + zznl(k_end)
         zlft  = 0.5_rk * (zrlft + zllft)
!
      else if (xyzswp .eq. 3) then
!
         do i = qz_s, k_end 
            utlt(i+4-config%q_zoff) = vyold(j,j_end,i)
         enddo
!
         yrlft =  yznl(1)
         yllft =  yznl(1) - yznr(j_end) + yznl(j_end)
         ylft  = 0.5_rk * (yrlft + yllft)
!
      endif
!
      go to 300
!
!     special boundary
!
250   continue
!
!     add any nonstandard boundary condition here
!
!
      go to 300
!-----------------------------------------------------------------------
!     MPI domain boundary  (theta-direction)
!-----------------------------------------------------------------------
260 continue
!
  if (xyzswp .eq. 3) then
!
     do i = 0, k_end - qz_s
        utlt(i+5) = vyold(j,qy_s-1,qz_s+i)
     enddo
!
     yrlft = yznr(qy_s-1)
     yllft = yznl(qy_s-1)
     ylft  = 0.5_rk * (yrlft + yllft)
  else 
     raise_abort("bndry(): error: left bndry")
  endif
!
  go to 300
!-----------------------------------------------------------------------
!     MPI domain boundary  (phi-direction)
!-----------------------------------------------------------------------
270 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        utlt(i+4) = vzold(i,j,qz_s-1)
     enddo
!
     zrlft = zznl(qz_s)
     zllft = zznl(qz_s) - zznr(k_end) + zznl(k_end)
     zlft  = 0.5_rk * (zrlft + zllft)
  else if (xyzswp .eq. 2) then
!
     do i = 0, j_end - qy_s
        utlt(i+5) = vzold(j,qy_s+i,qz_s-1)
     enddo
!
     zrlft = zznl(qz_s)
     zllft = zznl(qz_s) - zznr(k_end) + zznl(k_end)
     zlft  = 0.5_rk * (zrlft + zllft)
  else
     raise_abort("bndry(): error: left bndry")
  endif
!
!
!
300   continue
!
!
!=======================================================================
!     right of grid
!=======================================================================
!
      if (xyzswp .eq. 3) then
         nnn = qy_e !theta-direction is left-right-direction
      else 
         nnn = qz_e !phi-direction is left-right-direction
      endif

      if (k .ne. nnn) go to 400
!
      go to (310, 320, 330, 340, 350, 360, 370), ibndrt
!
!     reflecting boundary
!
310   continue
!
      if (xyzswp .eq. 1) then
!
         do i = 1, areas%nx
            utrt(i+4) = -vzold(i,j,k_end)
         enddo
!
         zlrgt = zznr(k_end)
         zrrgt = 2.0_rk  * zznr(k_end) - zznl(k_end)
         zrgt  = 0.5_rk * (zrrgt   + zlrgt)
!
      else if (xyzswp .eq. 2) then
!
         do i = qy_s, j_end
            utrt(i+4-config%q_yoff) = -vzold(j,i,k_end)
         enddo
!
         zlrgt =  zznr(k_end)
         zrrgt = 2.0_rk * zznr(k_end) - zznl(k_end)
         zrgt  = 0.5_rk * (zrrgt + zlrgt)
!
      else if (xyzswp .eq. 3) then

         do i = qz_s, k_end
            utrt(i+4-config%q_zoff) = -vyold(j,j_end,i)
         enddo
!
         ylrgt =  yznr(j_end)
         yrrgt = 2.0_rk * yznr(j_end) - yznl(j_end)
         yrgt  = 0.5_rk * (yrrgt + ylrgt)
!
      endif
!
      go to 400
!
!     flow out boundary
!
320   continue
!
!     value of utrt already correct
!
      if (xyzswp .ne. 3) then
!
         zlrgt = zznr(k_end)
         zrrgt = 2.0_rk * zznr(k_end) - zznl(k_end)
         zrgt  = 0.5_rk * (zrrgt + zlrgt)
!
      else
!
         ylrgt =  yznr(j_end)
         yrrgt = 2.0_rk  * yznr(j_end) - yznl(j_end)
         yrgt  = 0.5_rk * (yrrgt   + ylrgt)
!
      endif
!
      go to 400
!
!     flow in boundary
!
330   continue
!
      do i = 5, nzn4
         utrt(i) = uttin
      enddo
!
      if (xyzswp .ne. 3) then
!
         zlrgt = zznr(k_end)
         zrrgt = 2.0_rk  * zznr(k_end) - zznl(k_end)
         zrgt  = 0.5_rk * (zrrgt   + zlrgt)
!
      else
!
         ylrgt =  yznr(j_end)
         yrrgt = 2.0_rk * yznr(j_end) - yznl(j_end)
         yrgt  = 0.5_rk * (yrrgt + ylrgt)
!
      endif
!
      go to 400
!
!     periodic boundary
!
340   continue
!
      if (xyzswp .eq. 1) then
!

         do i = 1, areas%nx
            utrt(i+4) = vzold(i,j,qz_s)
         enddo

         zlrgt = zznr(k_end)
         zrrgt = zznr(k_end) + zznr(qz_s) - zznl(qz_s)

         zrgt  = 0.5_rk * (zrrgt + zlrgt)
!
      else if (xyzswp .eq. 2) then

         do i = qy_s, j_end
            utrt(i+4-config%q_yoff) = vzold(j,i,qz_s)
         enddo

         zlrgt = zznr(k_end)
         zrrgt = zznr(k_end) + zznr(qz_s) - zznl(qz_s)
         zrgt  = 0.5_rk * (zrrgt + zlrgt)
!
      else if (xyzswp .eq. 3) then

         do i = qz_s, k_end
            utrt(i+4-config%q_zoff) = vyold(j,qy_s,i)
         enddo

         ylrgt =  yznr(j_end)
         yrrgt =  yznr(j_end) + yznr(qy_s) -yznl(qy_s)

         yrgt  = 0.5_rk * (yrrgt + ylrgt)
!
      endif
!
      go to 400
!
!     special boundary
!
350    continue
!
!     add any nonstandard boundary condition here
!
!
      go to 400
!-----------------------------------------------------------------------
!     MPI domain boundary  (theta-direction)
!-----------------------------------------------------------------------
360 continue
!
  if (xyzswp .eq. 3) then
!
     do i = 0, k_end - qz_s
        utrt(i+5) = vyold(j,qy_e+1,qz_s+i)
     enddo
!
     yrrgt = yznr(qy_e+1)
     ylrgt = yznl(qy_e+1)
     yrgt  = 0.5_rk * (yrrgt + ylrgt)
  else 
     raise_abort("bndry(): error: left bndry")
  endif
!
  go to 400
!-----------------------------------------------------------------------
!     MPI domain boundary  (phi-direction)
!-----------------------------------------------------------------------
370 continue
!
  if (xyzswp .eq. 1) then
!
     do i = 1, areas%nx
        utrt(i+5) = vzold(i,j,qz_e+1)
     enddo
!
     zlrgt = zznr(k_end)
     zrrgt = zznr(k_end) + zznr(qz_s) - zznl(qz_s)
     zrgt  = 0.5_rk * (zrrgt + zlrgt)
  else if (xyzswp .eq. 2) then
!
     do i = 0, j_end - qy_s
        utrt(i+5) = vzold(j,qy_s+i,qz_e+1)
     enddo
!
     zlrgt = zznr(k_end)
     zrrgt = zznr(k_end) + zznr(qz_s) - zznl(qz_s)
     zrgt  = 0.5_rk * (zrrgt + zlrgt)
  else
     raise_abort("bndry(): error: left bndry")
  endif
!
!
!
400   continue
!
!
!=======================================================================
!     define zones 1 through 4
!=======================================================================
!
      ibndmn = bndmin
      ibndmx = bndmax
!

      if (.not.(use_mpi)) then
         go to (410, 420, 430, 440, 450), ibndmn
      else

         ! 460/470 is not a special boundary condition, 
         ! it's because of the parallelisation handling
         go to (410, 420, 430, 440, 450, 460, 470, 480, 490, 495), ibndmn
      endif
!
!     reflecting boundary
!
410   continue
!
      do i = 1, 4
         rho (i) =  rho (9-i)
         u   (i) = -u   (9-i)
         ut  (i) =  ut  (9-i)
         utt (i) =  utt (9-i)
         utbt(i) =  utbt(9-i)
         uttp(i) =  uttp(9-i)
         utlt(i) =  utlt(9-i)
         utrt(i) =  utrt(9-i)
         e   (i) =  e   (9-i)
!
! ...  take care !!! potential is defined on the zone interfaces
!
         grav(i) =  grav(8-i)
         gamc(i) =  gamc(9-i)
         game(i) =  game(9-i)
         p   (i) =  p   (9-i)
         dx  (i) =  dx  (9-i)
         ugrid(i)= -ugrid(9-i)
         do n = 1, config%qn
           xn(i,n) = xn(9-i,n)
        enddo
      enddo
      if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
         ! switch sign of v_phi when going across
         ! the axis through the origin in spherical
         ! polar coordinates (this is not an optimal
         ! solution)
         do i = 1, 4
            utt (i) = - utt (i)
            utrt(i) = - utrt(i)
            utlt(i) = - utlt(i)
         enddo
      endif
!
      go to 500
!
!     flow out boundary
!
420   continue
!
      do i = 1, 4
         rho (i) = rho (5)
!-new    u   (i) = amin1( 0.0, u(5) )
         u   (i) = u   (5)
         ut  (i) = ut  (5)
         utt (i) = utt (5)
         utbt(i) = utbt(5)
         uttp(i) = uttp(5)
         utlt(i) = utlt(5)
         utrt(i) = utrt(5)
         e   (i) = e   (5)
         grav(i) = grav(4)
         gamc(i) = gamc(5)
         game(i) = game(5)
         p   (i) = p   (5)
         dx  (i) = dx  (5)
         ugrid (i) = ugrid(5)
         do n = 1, config%qn
           xn  (i,n) = xn(5,n)
        enddo
      enddo

      go to 500
!
!     flow in boundary
!
430   continue
!
      do i = 1, 4
         rho (i) = rhoin
         u   (i) = uin
         ut  (i) = utin
         utt (i) = uttin
         utbt(i) = utin
         uttp(i) = utin
         utlt(i) = uttin
         utrt(i) = uttin
         e   (i) = ein
         grav(i) = gravin
         p   (i) = pin
         gamc(i) = gamcin
         game(i) = gamein
         dx  (i) = dx(5)
         ugrid (i) = 0.0_rk
      do n = 1, config%qn
           xn  (i,n) = xnin(n)
        enddo
      enddo
!
      go to 500
!
!     periodic boundary
!
440   continue
!
      do i = 1, 4
         rho  (i) = rho (i+nzn)
         u    (i) = u   (i+nzn)
         ut   (i) = ut  (i+nzn)
         utt  (i) = utt (i+nzn)
         utbt (i) = utbt(i+nzn)
         uttp (i) = uttp(i+nzn)
         utlt (i) = utlt(i+nzn)
         utrt (i) = utrt(i+nzn)
         e    (i) = e   (i+nzn)
         grav (i) = grav(i+nzn)
         gamc (i) = gamc(i+nzn)
         game (i) = game(i+nzn)
         p    (i) = p   (i+nzn)
         dx   (i) = dx  (i+nzn)
         ugrid(i)   = ugrid(i+nzn)
         do n = 1, config%qn
           xn   (i,n) = xn (i+nzn,n)
        enddo
      enddo
!
      go to 500
!
!=======================================================================
!     special inner-boundary condition (x-direction):
!=======================================================================

450   continue

      if(xyzswp .ne. 1) then
         raise_abort("bndry.F(): wrong sweep direction at 450!")
      end if
!      j1 = max0 (j - 1, 1)
!      k1 = max0 (k - 1, 1)

!      j2 = min0 (j + 1, max(qy_s,ny))
!      k2 = min0 (k + 1, max(qz_s,nz))


      j1 = qy_s
      j2 = qy_s
      k1 = qz_s
      k2 = qz_s

!-----------------------------------------------------------------------
!     copy values from the boundary interface in boundary zones:
!-----------------------------------------------------------------------

      do i = 1, 4
         ibi       = (are_id - 1)*8 + i
         rho (i)   = den_bi(ibi, j , k )
         u   (i)   = vex_bi(ibi, j , k )
         ut  (i)   = vey_bi(ibi, j , k )
         utt (i)   = vez_bi(ibi, j , k )
         uttp(i)   = vey_bi(ibi, j2, k )
         utbt(i)   = vey_bi(ibi, j1, k )
         utrt(i)   = vez_bi(ibi, j , k2)
         utlt(i)   = vez_bi(ibi, j , k1)
         e   (i)   = ene_bi(ibi, j , k )
         grav(i)   = gra_bi(ibi, j , k )
         gamc(i)   = gac_bi(ibi, j , k )
         game(i)   = gae_bi(ibi, j , k )
         p   (i)   = pre_bi(ibi, j , k )
         dx  (i)   = dxx_bi(ibi)
         ugrid (i) = ugr_bi(ibi)
         do n  = 1, config%qn
           xn(i,n) = xnu_bi(ibi,j,k,n)
        enddo
      enddo

      goto 500

460      continue
         if (xyzswp .eq. 2) then
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (1:4)=densty_lb (j, 1:4,k)       
            u   (1:4)=vely_lb   (j, 1:4,k)
            ut  (1:4)=velx_lb   (j, 1:4,k)
            utt (1:4)=velz_lb   (j, 1:4,k)
            uttp(1:4)=vxold_lb  (i2,1:4,k)  
            utbt(1:4)=vxold_lb  (i1,1:4,k)  
            utrt(1:4)=vzold_lb  (j, 1:4,k2)   
            utlt(1:4)=vzold_lb  (j, 1:4,k1)   
            e   (1:4)=energy_lb (j, 1:4,k)
            p   (1:4)=press_lb  (j, 1:4,k)
            tmp (1:4)=temp_lb   (j, 1:4,k)
            grav(1:4)=gpocntr_lb(j, 1:4,k)
            game(1:4)=gammae_lb (j, 1:4,k)
            gamc(1:4)=gammac_lb (j, 1:4,k) 
            xn  (1:4,1:config%qn) = xnuc_lb(j,1:4,k,1:config%qn)

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j.eq. 1) utbt(1:4)=-utbt(1:4)

            do i=1,4
               ugrid(i)=0.0_rk
               xl   (i) = yznl  (qy_s-5+i)
               x    (i) = yzn   (qy_s-5+i)
               xr   (i) = yznr  (qy_s-5+i)
               dx   (i) = xr(i) - xl(i)
            end do
         else
            raise_abort("bndry(): error: inner bndry")
         !         call stopit_mpi("bndry: error: inner bndry")
         end if
      
         goto 500

470      continue

         if (config%nsdim .eq. 3) then 
            if (xyzswp .eq. 3) then
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qy)
            
            rho (1:4)=densty_pb (j, k,1:4)
            u   (1:4)=velz_pb   (j, k,1:4)
            ut  (1:4)=velx_pb   (j, k,1:4)
            utt (1:4)=vely_pb   (j, k,1:4)
            uttp(1:4)=vxold_pb  (i2,k,1:4)
            utbt(1:4)=vxold_pb  (i1,k,1:4)
            utrt(1:4)=vyold_pb  (j, k2,1:4)
            utlt(1:4)=vyold_pb  (j, k1,1:4)
            e   (1:4)=energy_pb (j, k,1:4)
            p   (1:4)=press_pb  (j, k,1:4)
            tmp (1:4)=temp_pb   (j, k,1:4)
            grav(1:4)=gpocntr_pb(j, k,1:4)
            game(1:4)=gammae_pb (j, k,1:4)
            gamc(1:4)=gammac_pb (j, k,1:4)
            xn  (1:4,1:config%qn) = xnuc_pb(j,k,1:4,1:config%qn)
            
            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j.eq. 1) utbt(1:4)=-utbt(1:4)

            if (use_1neighbour_comm) then ! at least 4 zones per MPI task
               do i=1,4
                  ugrid(i)=0.0_rk
                  if (qz_s .eq. 1) then
                     xl   (i) = zznl  (config%qz-4+i) -2._rk*pc_pi
                     x    (i) = zzn   (config%qz-4+i) -2._rk*pc_pi
                     xr   (i) = zznr  (config%qz-4+i) -2._rk*pc_pi
                  else
                     xl   (i) = zznl  (qz_s-5+i)
                     x    (i) = zzn   (qz_s-5+i)
                     xr   (i) = zznr  (qz_s-5+i)
                  endif
                  dx   (i) = xr(i) - xl(i)
               end do
            endif
           
            if (use_2neighbour_comm) then ! two zones per MPI task
               do i=1,4
                  ugrid(i)=0.0_rk
                  if (qz_s .eq. 1) then
                     xl   (i) = zznl  (config%qz-4+i) -2._rk*pc_pi
                     x    (i) = zzn   (config%qz-4+i) -2._rk*pc_pi
                     xr   (i) = zznr  (config%qz-4+i) -2._rk*pc_pi
                  else if (qz_s .eq. 3) then
                     if (i .eq. 1 .or. i .eq. 2) then
                        xl   (i) = zznl  (config%qz-2+i) -2._rk*pc_pi
                        x    (i) = zzn   (config%qz-2+i) -2._rk*pc_pi
                        xr   (i) = zznr  (config%qz-2+i) -2._rk*pc_pi
                     else
                        xl   (i) = zznl  (qz_s-5+i)
                        x    (i) = zzn   (qz_s-5+i)
                        xr   (i) = zznr  (qz_s-5+i)
                     endif
                  else
                     xl   (i) = zznl  (qz_s-5+i)
                     x    (i) = zzn   (qz_s-5+i)
                     xr   (i) = zznr  (qz_s-5+i)
                  endif
                  dx   (i) = xr(i) - xl(i)
               end do
            endif
            
            if (use_1neighbour_comm) then ! one zones per MPI task
               do i=1,4
                  ugrid(i)=0.0_rk
                  if (qz_s .eq. 1) then  ! all four values are obtained
                                         ! via boundary condition
                     xl   (i) = zznl  (config%qz-4+i) -2._rk*pc_pi
                     x    (i) = zzn   (config%qz-4+i) -2._rk*pc_pi
                     xr   (i) = zznr  (config%qz-4+i) -2._rk*pc_pi
                  else if (qz_s .eq. 2) then ! thirst three values via
                                             ! boundary condition, fourth
                                             ! from zone qz_s = 1
                     if (i .le. 3) then
                        xl   (i) = zznl  (config%qz-3+i) -2._rk*pc_pi
                        x    (i) = zzn   (config%qz-3+i) -2._rk*pc_pi
                        xr   (i) = zznr  (config%qz-3+i) -2._rk*pc_pi
                     else
                        xl   (i) = zznl  (qz_s-5+i)
                        x    (i) = zzn   (qz_s-5+i)
                        xr   (i) = zznr  (qz_s-5+i)
                     endif
                     
                  else if (qz_s .eq. 3) then ! first two values via
                                             ! boundary condition
                                             ! third, fourth from zones
                                             ! qz_s-2, qz_s -1
                     if (i .le. 2) then
                        xl   (i) = zznl  (config%qz-2+i) -2._rk*pc_pi
                        x    (i) = zzn   (config%qz-2+i) -2._rk*pc_pi
                        xr   (i) = zznr  (config%qz-2+i) -2._rk*pc_pi
                     else
                        xl   (i) = zznl  (qz_s-5+i)
                        x    (i) = zzn   (qz_s-5+i)
                        xr   (i) = zznr  (qz_s-5+i)
                     endif
                     
                  else if (qz_s .eq. 4) then ! first value via
                                             ! boundary condition
                                             ! second, third, fourth from zones
                                             ! qz_s-3, qz_s-2, qz_s -1
                     if (i .eq. 1) then
                        xl   (i) = zznl  (config%qz-1+i) -2._rk*pc_pi
                        x    (i) = zzn   (config%qz-1+i) -2._rk*pc_pi
                        xr   (i) = zznr  (config%qz-1+i) -2._rk*pc_pi
                     else
                        xl   (i) = zznl  (qz_s-5+i)
                        x    (i) = zzn   (qz_s-5+i)
                        xr   (i) = zznr  (qz_s-5+i)
                     endif
                     
                  else  ! all four values from zones
                     xl   (i) = zznl  (qz_s-5+i)
                     x    (i) = zzn   (qz_s-5+i)
                     xr   (i) = zznr  (qz_s-5+i)
                  endif
                  dx   (i) = xr(i) - xl(i)
               end do
            endif
         else ! xyzswp .ne. 3
           raise_abort("bndry(): error: inner bndry")
            !         call stopit_mpi("bndry: error: inner bndry")
        end if ! xyzswp = 3 ?
     endif ! config%nsdim
         
         goto 500
!
!     special case: first zone: reflecting boundary, other three: MPI communication
!
480 continue

         if (.not.(use_1neighbour_comm)) then
            raise_abort("bndry(): error: called boundary case 480 but not use_1neighbour_comm")
         endif

         if (xyzswp .eq. 2) then
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (2:4)=densty_lb (j, 2:4,k)       
            u   (2:4)=vely_lb   (j, 2:4,k)
            ut  (2:4)=velx_lb   (j, 2:4,k)
            utt (2:4)=velz_lb   (j, 2:4,k)
            uttp(2:4)=vxold_lb  (i2,2:4,k)  
            utbt(2:4)=vxold_lb  (i1,2:4,k)  
            utrt(2:4)=vzold_lb  (j, 2:4,k2)   
            utlt(2:4)=vzold_lb  (j, 2:4,k1)   
            e   (2:4)=energy_lb (j, 2:4,k)
            p   (2:4)=press_lb  (j, 2:4,k)
            tmp (2:4)=temp_lb   (j, 2:4,k)
            grav(2:4)=gpocntr_lb(j, 2:4,k)
            game(2:4)=gammae_lb (j, 2:4,k)
            gamc(2:4)=gammac_lb (j, 2:4,k) 
            xn  (2:4,1:config%qn) = xnuc_lb(j,2:4,k,1:config%qn)

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j.eq. 1) utbt(2:4)=-utbt(2:4)

            do i=2,4
               ugrid(i)=0.0_rk
               xl   (i) = yznl  (qy_s-5+i)
               x    (i) = yzn   (qy_s-5+i)
               xr   (i) = yznr  (qy_s-5+i)
               dx   (i) = xr(i) - xl(i)
            end do
            
            ! zone 1 is reflective w.r.t. zone 2
            do i = 1, 1
               rho (i) =  rho (3-i)
               u   (i) = -u   (3-i)
               ut  (i) =  ut  (3-i)
               utt (i) =  utt (3-i)
               utbt(i) =  utbt(3-i)
               uttp(i) =  uttp(3-i)
               utlt(i) =  utlt(3-i)
               utrt(i) =  utrt(3-i)
               e   (i) =  e   (3-i)
               !
! ...  take care !!! potential is defined on the zone interfaces
!
               grav(i) =  grav(3-i)
               gamc(i) =  gamc(3-i)
               game(i) =  game(3-i)
               p   (i) =  p   (3-i)
               dx  (i) =  dx  (3-i)
               ugrid(i)= -ugrid(3-i)
               do n = 1, config%qn
                  xn(i,n) = xn(3-i,n)
               enddo
            enddo
            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
              ! the axis through the origin in spherical
              ! polar coordinates (this is not an optimal
              ! solution)
              do i = 1, 1
                utt (i) = - utt (i)
                utrt(i) = - utrt(i)
                utlt(i) = - utlt(i)
              enddo
            endif
         else
            raise_abort("bndry(): error: inner bndry in case 480")
         !         call stopit_mpi("bndry: error: inner bndry")
         end if
!
         go to 500       
!
!     special case: first two zones: reflecting boundary, other two: MPI communication
!
490   continue
         if (.not.(use_1neighbour_comm) .and. .not.(use_2neighbour_comm)) then
            raise_abort("bndry(): error: called boundary case 490 but not use_1neighbour_comm or use_2neighbour_comm" )
         endif

         if (xyzswp .eq. 2) then
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (3:4)=densty_lb (j, 3:4,k)       
            u   (3:4)=vely_lb   (j, 3:4,k)
            ut  (3:4)=velx_lb   (j, 3:4,k)
            utt (3:4)=velz_lb   (j, 3:4,k)
            uttp(3:4)=vxold_lb  (i2,3:4,k)  
            utbt(3:4)=vxold_lb  (i1,3:4,k)  
            utrt(3:4)=vzold_lb  (j, 3:4,k2)   
            utlt(3:4)=vzold_lb  (j, 3:4,k1)   
            e   (3:4)=energy_lb (j, 3:4,k)
            p   (3:4)=press_lb  (j, 3:4,k)
            tmp (3:4)=temp_lb   (j, 3:4,k)
            grav(3:4)=gpocntr_lb(j, 3:4,k)
            game(3:4)=gammae_lb (j, 3:4,k)
            gamc(3:4)=gammac_lb (j, 3:4,k) 
            xn  (3:4,1:config%qn) = xnuc_lb(j,3:4,k,1:config%qn)

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j.eq. 1) utbt(3:4)=-utbt(3:4)

            do i=3,4
               ugrid(i)=0.0_rk
               xl   (i) = yznl  (qy_s-5+i)
               x    (i) = yzn   (qy_s-5+i)
               xr   (i) = yznr  (qy_s-5+i)
               dx   (i) = xr(i) - xl(i)
            end do
            
            do i = 1, 2
               rho (i) =  rho (5-i)
               u   (i) = -u   (5-i)
               ut  (i) =  ut  (5-i)
               utt (i) =  utt (5-i)
               utbt(i) =  utbt(5-i)
               uttp(i) =  uttp(5-i)
               utlt(i) =  utlt(5-i)
               utrt(i) =  utrt(5-i)
               e   (i) =  e   (5-i)
               !
! ...  take care !!! potential is defined on the zone interfaces
!
               grav(i) =  grav(5-i)
               gamc(i) =  gamc(5-i)
               game(i) =  game(5-i)
               p   (i) =  p   (5-i)
               dx  (i) =  dx  (5-i)
               ugrid(i)= -ugrid(5-i)
               do n = 1, config%qn
                  xn(i,n) = xn(5-i,n)
               enddo
            enddo
            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
              ! the axis through the origin in spherical
              ! polar coordinates (this is not an optimal
              ! solution)
               do i = 1, 2
                  utt (i) = - utt (i)
                  utrt(i) = - utrt(i)
                  utlt(i) = - utlt(i)
               enddo
            endif
            
            
         else
            raise_abort("bndry(): error: inner bndry in case 490")
         !         call stopit_mpi("bndry: error: inner bndry")
         end if
!
      go to 500


!
!     special case: first three zones: reflecting boundary, other one: MPI communication
!
!     Carefull: this case implies that the lower boundary zone
!               1 is the "reflection" of the zones qy_s+1 (i.e.
!               _right_ of the current task qy_s = 2
!               thus the values from densty_ub(:,1,:) are needed for boundary
!               zone 1. Boundary zone to and three correspond to qy_s-1 and
!               qy_s (thus to densty_lb)
495 continue
         if (.not.(use_1neighbour_comm)) then
            raise_abort("bndry(): error: called boundary case 490 but not use_1neighbour_comm " )
         endif

         if (xyzswp .eq. 2) then
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (4:4)=densty_lb (j, 4:4,k)       
            u   (4:4)=vely_lb   (j, 4:4,k)
            ut  (4:4)=velx_lb   (j, 4:4,k)
            utt (4:4)=velz_lb   (j, 4:4,k)
            uttp(4:4)=vxold_lb  (i2,4:4,k)  
            utbt(4:4)=vxold_lb  (i1,4:4,k)  
            utrt(4:4)=vzold_lb  (j, 4:4,k2)   
            utlt(4:4)=vzold_lb  (j, 4:4,k1)   
            e   (4:4)=energy_lb (j, 4:4,k)
            p   (4:4)=press_lb  (j, 4:4,k)
            tmp (4:4)=temp_lb   (j, 4:4,k)
            grav(4:4)=gpocntr_lb(j, 4:4,k)
            game(4:4)=gammae_lb (j, 4:4,k)
            gamc(4:4)=gammac_lb (j, 4:4,k) 
            xn  (4:4,1:config%qn) = xnuc_lb(j,4:4,k,1:config%qn)

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j.eq. 1) utbt(4:4)=-utbt(4:4)

            do i=4,4
               ugrid(i)=0.0_rk
               xl   (i) = yznl  (qy_s-5+i)
               x    (i) = yzn   (qy_s-5+i)
               xr   (i) = yznr  (qy_s-5+i)
               dx   (i) = xr(i) - xl(i)
            end do
            
            do i = 2, 3
               rho (i) =  rho (7-i)
               u   (i) = -u   (7-i)
               ut  (i) =  ut  (7-i)
               utt (i) =  utt (7-i)
               utbt(i) =  utbt(7-i)
               uttp(i) =  uttp(7-i)
               utlt(i) =  utlt(7-i)
               utrt(i) =  utrt(7-i)
               e   (i) =  e   (7-i)
               !
! ...  take care !!! potential is defined on the zone interfaces
!
               grav(i) =  grav(7-i)
               gamc(i) =  gamc(7-i)
               game(i) =  game(7-i)
               p   (i) =  p   (7-i)
               dx  (i) =  dx  (7-i)
               ugrid(i)= -ugrid(7-i)
               do n = 1, config%qn
                  xn(i,n) = xn(7-i,n)
               enddo
            enddo

            do i = 1, 1
               rho (i) =  densty_ub(j,i,k)
               u   (i) = -vely_ub(j,i,k)
               ut  (i) =  velx_ub(j,i,k)
               utt (i) =  velz_ub(j,i,k)
               utbt(i) =  vxold_ub(i1,i,k)
               uttp(i) =  vxold_ub(i2,i,k)
               utlt(i) =  vzold_ub(j,i,k1)
               utrt(i) =  vzold_ub(j,i,k2)
               e   (i) =  energy_ub(j,i,k)
               !
! ...  take care !!! potential is defined on the zone interfaces
!
               grav(i) =  gpocntr_ub(j,i,k)
               gamc(i) =  gammac_ub(j,i,k)
               game(i) =  gammae_ub(j,i,k)
               p   (i) =  press_ub(j,i,k)
               dx  (i) =  dx(qy_s+1)  
               ugrid(i)= -ugrid(qy_s+1)
               do n = 1, config%qn
                  xn(i,n) = xnuc_ub(j,i,k,n)
               enddo
            enddo

            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
               ! the axis through the origin in spherical
               ! polar coordinates (this is not an optimal
               ! solution)
               do i = 1, 3
                  utt (i) = - utt (i)
                  utrt(i) = - utrt(i)
                  utlt(i) = - utlt(i)
               enddo
            endif
         else
            raise_abort("bndry(): error: inner bndry in case 495")
         !         call stopit_mpi("bndry: error: inner bndry")
         end if
         go to 500   
!=======================================================================
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


500   continue

      do i = 1, 4
         xl(5-i) =        xl(6-i) - dx(5-i)
         xr(5-i) =        xr(6-i) - dx(6-i)
         x (5-i) = 0.5_rk * (xl(5-i) + xr(5-i))
      enddo
!
!
!=======================================================================
!     define zones nzn5 through nzn8
!=======================================================================
!
      if (.not.(use_mpi)) then
         go to (510, 520, 530, 540, 550), ibndmx
      else
         go to (510, 520, 530, 540, 550, 560, 570, 580, 590, 595), ibndmx
      endif
!
!     reflecting boundary
!
510   continue
!
      nzn54 = nzn5 + nzn4
!
      do i = nzn5, nzn8
         rho (i) =  rho (nzn54-i)
         u   (i) = -u   (nzn54-i)
         ut  (i) =  ut  (nzn54-i)
         utt (i) =  utt (nzn54-i)
         utbt(i) =  utbt(nzn54-i)
         uttp(i) =  uttp(nzn54-i)
         utlt(i) =  utlt(nzn54-i)
         utrt(i) =  utrt(nzn54-i)
         e   (i) =  e   (nzn54-i)
!
! ...  take care !!! potential is defined on the zone interfaces
!
         grav(i) =  grav(nzn54-i-1)
         gamc(i) =  gamc(nzn54-i)
         game(i) =  game(nzn54-i)
         p   (i) =  p   (nzn54-i)
         dx  (i) =  dx  (nzn54-i)
!-PNS         ugrid (i) =-ugrid(nzn54-i)
!-PNS1         ugrid (i) = ugrid(nzn54-i)
         ugrid (i) = ugrid(nzn4)
         do n = 1, config%qn
           xn  (i,n) = xn (nzn54-i,n)
        enddo
      enddo

      if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
         ! switch sign of v_phi when going across
         ! the axis through the origin in spherical
         ! polar coordinates (this is not an optimal
         ! solution)
         do i = nzn5, nzn8
            utt (i) = - utt (i)
            utrt(i) = - utrt(i)
            utlt(i) = - utlt(i)
         enddo
      endif

!
      go to 600
!
!     flow out boundary
!
520   continue
!
      do i = nzn5, nzn8

!-PNS    rho (i) = rho (nzn4)
!-PNS>
         rho (i) = config%smlrho
         u   (i) = max(u(nzn4),0.0_rk) ! kein Einstroemen

!-PNS>
!-PNS    u   (i) = u   (nzn4)
         ut  (i) = ut  (nzn4)
         utt (i) = utt (nzn4)
         utbt(i) = utbt(nzn4)
         uttp(i) = uttp(nzn4)
         utlt(i) = utlt(nzn4)
         utrt(i) = utrt(nzn4)
         e   (i) = e   (nzn4)
         grav(i) = grav(nzn4)
         gamc(i) = gamc(nzn4)
         game(i) = game(nzn4)
         p   (i) = p   (nzn4)
         dx  (i) = dx  (nzn4)
         ugrid (i) = ugrid(nzn4)
         do n = 1, config%qn
           xn  (i,n) = xn (nzn4,n)
        enddo
      enddo
!
      go to 600
!
!     flow in boundary
!
530   continue
!
      do i = nzn5, nzn8
         rho (i) = rhoin
         u   (i) = uin
         ut  (i) = utin
         utt (i) = uttin
         utbt(i) = utin
         uttp(i) = utin
         utlt(i) = uttin
         utrt(i) = uttin
         e   (i) = ein
         grav(i) = gravin
         p   (i) = pin
         gamc(i) = gamcin
         game(i) = gamein
         dx  (i) = dx(nzn4)
         ugrid (i) = 0.0_rk
         do n = 1, config%qn
           xn  (i,n) = xnin(n)
        enddo
      enddo
!
      go to 600
!
!     periodic boundary
!
540   continue
!
      do i = nzn5, nzn8
         rho (i) = rho (i-nzn)
         u   (i) = u   (i-nzn)
         ut  (i) = ut  (i-nzn)
         utt (i) = utt (i-nzn)
         utbt(i) = utbt(i-nzn)
         uttp(i) = uttp(i-nzn)
         utlt(i) = utlt(i-nzn)
         utrt(i) = utrt(i-nzn)
         e   (i) = e   (i-nzn)
         grav(i) = grav(i-nzn)
         p   (i) = p   (i-nzn)
         gamc(i) = gamc(i-nzn)
         game(i) = game(i-nzn)
         dx  (i) = dx  (i-nzn)
         ugrid (i) = ugrid(i-nzn)
         do n = 1, config%qn
           xn  (i,n) = xn (i-nzn,n)
        enddo
      enddo
!
      go to 600
!
!=======================================================================
!     special outer boundary condition (x-direction):
!=======================================================================

550   continue

      if(xyzswp .ne. 1) then
         raise_abort("bndry.F(): wrong sweep direction at 550!")
      end if

!      j1 = max0 (j - 1, 1)
!      k1 = max0 (k - 1, 1)

!      j2 = min0 (j + 1, max(qy_s,ny))
!      k2 = min0 (k + 1, max(qz_s,nz))
!

      j1 = qy_s
      j2 = qy_s
      k1 = qz_s
      k2 = qz_s

!-----------------------------------------------------------------------
!     copy values from the boundary interface in boundary zones:
!-----------------------------------------------------------------------

      do i = nzn5, nzn8
         ibi       = are_id*8  + 4 + i - nzn5 + 1
         rho (i)   = den_bi(ibi, j , k )
         u   (i)   = vex_bi(ibi, j , k )
         ut  (i)   = vey_bi(ibi, j , k )
         utt (i)   = vez_bi(ibi, j , k )
         uttp(i)   = vey_bi(ibi, j2, k )
         utbt(i)   = vey_bi(ibi, j1, k )
         utrt(i)   = vez_bi(ibi, j , k2)
         utlt(i)   = vez_bi(ibi, j , k1)
         e   (i)   = ene_bi(ibi, j , k )
         grav(i)   = gra_bi(ibi, j , k )
         gamc(i)   = gac_bi(ibi, j , k )
         game(i)   = gae_bi(ibi, j , k )
         p   (i)   = pre_bi(ibi, j , k )
         dx  (i)   = dxx_bi(ibi)
         ugrid (i) = ugr_bi(ibi)
         do n  = 1, config%qn
           xn(i,n) = xnu_bi(ibi,j,k,n)
        enddo
      enddo

      goto 600

560      continue 
         
         if (xyzswp .eq. 2) then
            
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (nzn5:nzn8)=densty_ub (j ,1:4,k)     
            u   (nzn5:nzn8)=vely_ub   (j ,1:4,k)
            ut  (nzn5:nzn8)=velx_ub   (j ,1:4,k)
            utt (nzn5:nzn8)=velz_ub   (j ,1:4,k)
            uttp(nzn5:nzn8)=vxold_ub  (i2,1:4,k)    
            utbt(nzn5:nzn8)=vxold_ub  (i1,1:4,k)    
            utrt(nzn5:nzn8)=vzold_ub  (j ,1:4,k2)    
            utlt(nzn5:nzn8)=vzold_ub  (j ,1:4,k1)    
            e   (nzn5:nzn8)=energy_ub (j ,1:4,k)
            p   (nzn5:nzn8)=press_ub  (j ,1:4,k)
            tmp (nzn5:nzn8)=temp_ub   (j ,1:4,k)
            grav(nzn5:nzn8)=gpocntr_ub(j ,1:4,k)
            game(nzn5:nzn8)=gammae_ub (j ,1:4,k)
            gamc(nzn5:nzn8)=gammac_ub (j ,1:4,k)
            xn  (nzn5:nzn8,1:config%qn) = xnuc_ub(j,1:4,k,1:config%qn)                  

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j .eq. 1) utbt(nzn5:nzn8)=-utbt(nzn5:nzn8)
            if (j .eq. config%qx) uttp(nzn5:nzn8)=-uttp(nzn5:nzn8)
            
            do i = 1,4
               ugrid(nzn4+i) = 0.0_rk
               xl   (nzn4+i) = yznl  (qy_e+i)
               x    (nzn4+i) = yzn   (qy_e+i)
               xr   (nzn4+i) = yznr  (qy_e+i)
               dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
            end do
         else
            raise_abort("bndry(): error: outer bndry")
            !         call stopit_mpi("bndry: error: outer bndry")
         end if


         goto 600
         
570      continue 

      if (config%nsdim .eq. 3) then
         if (xyzswp .eq. 3) then

            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qy)
            
            rho (nzn5:nzn8)=densty_kb (j ,k,1:4)
            u   (nzn5:nzn8)=velz_kb   (j ,k,1:4)
            ut  (nzn5:nzn8)=velx_kb   (j ,k,1:4)
            utt (nzn5:nzn8)=vely_kb   (j ,k,1:4)
            uttp(nzn5:nzn8)=vxold_kb  (i2,k,1:4)
            utbt(nzn5:nzn8)=vxold_kb  (i1,k,1:4)
            utrt(nzn5:nzn8)=vyold_kb  (j ,k2,1:4)
            utlt(nzn5:nzn8)=vyold_kb  (j ,k1,1:4)
            e   (nzn5:nzn8)=energy_kb (j ,k,1:4)
            p   (nzn5:nzn8)=press_kb  (j ,k,1:4)
            tmp (nzn5:nzn8)=temp_kb   (j ,k,1:4)
            grav(nzn5:nzn8)=gpocntr_kb(j ,k,1:4)
            game(nzn5:nzn8)=gammae_kb (j ,k,1:4)
            gamc(nzn5:nzn8)=gammac_kb (j ,k,1:4)
            xn  (nzn5:nzn8,1:config%qn) = xnuc_kb(j,k,1:4,1:config%qn)                  
            
            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j .eq. 1) utbt(nzn5:nzn8)=-utbt(nzn5:nzn8)
            if (j .eq.config%qx) uttp(nzn5:nzn8)=-uttp(nzn5:nzn8)
            
              
            if (qz_s+3 .le. qz_e) then
              do i = 1,4
                ugrid(nzn4+i)=0.0_rk
                if (qz_e .eq. config%qz) then
                  xl   (nzn4+i) = zznl  (i) +2._rk*pc_pi
                  x    (nzn4+i) = zzn   (i) +2._rk*pc_pi
                  xr   (nzn4+i) = zznr  (i) +2._rk*pc_pi
                else
                  xl   (nzn4+i) = zznl  (qz_e+i)
                  x    (nzn4+i) = zzn   (qz_e+i)
                  xr   (nzn4+i) = zznr  (qz_e+i)
                endif
                dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
              end do
            else
              do i = 1,4
                ugrid(nzn4+i)=0.0_rk
                if (qz_e .eq. config%qz) then
                  xl   (nzn4+i) = zznl  (i) +2._rk*pc_pi
                  x    (nzn4+i) = zzn   (i) +2._rk*pc_pi
                  xr   (nzn4+i) = zznr  (i) +2._rk*pc_pi
                else if (qz_e .eq. config%qz-2) then
                  if (i.eq.3 .or. i.eq.4) then
                    xl   (nzn4+i) = zznl  (i-2) +2._rk*pc_pi
                    x    (nzn4+i) = zzn   (i-2) +2._rk*pc_pi
                    xr   (nzn4+i) = zznr  (i-2) +2._rk*pc_pi
                  else
                    xl   (nzn4+i) = zznl  (qz_e+i)
                    x    (nzn4+i) = zzn   (qz_e+i)
                    xr   (nzn4+i) = zznr  (qz_e+i)
                  endif
                else
                  xl   (nzn4+i) = zznl  (qz_e+i)
                  x    (nzn4+i) = zzn   (qz_e+i)
                  xr   (nzn4+i) = zznr  (qz_e+i)
                endif
                dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
              end do
            endif
         else
            raise_abort("bndry(): error: outer bndry")
         end if
         
      endif ! config%nsdim

         goto 600

580      continue 
         ! special case: last zone reflective, three communication
         if (.not.(use_1neighbour_comm)) then
            raise_abort("bndry(): error: called boundary case 580 but not use_1neighbour_comm")
         endif

        if (xyzswp .eq. 2) then
            
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (nzn5:nzn7)=densty_ub (j ,1:3,k)     
            u   (nzn5:nzn7)=vely_ub   (j ,1:3,k)
            ut  (nzn5:nzn7)=velx_ub   (j ,1:3,k)
            utt (nzn5:nzn7)=velz_ub   (j ,1:3,k)
            uttp(nzn5:nzn7)=vxold_ub  (i2,1:3,k)    
            utbt(nzn5:nzn7)=vxold_ub  (i1,1:3,k)    
            utrt(nzn5:nzn7)=vzold_ub  (j ,1:3,k2)    
            utlt(nzn5:nzn7)=vzold_ub  (j ,1:3,k1)    
            e   (nzn5:nzn7)=energy_ub (j ,1:3,k)
            p   (nzn5:nzn7)=press_ub  (j ,1:3,k)
            tmp (nzn5:nzn7)=temp_ub   (j ,1:3,k)
            grav(nzn5:nzn7)=gpocntr_ub(j ,1:3,k)
            game(nzn5:nzn7)=gammae_ub (j ,1:3,k)
            gamc(nzn5:nzn7)=gammac_ub (j ,1:3,k)
            xn  (nzn5:nzn7,1:config%qn) = xnuc_ub(j,1:3,k,1:config%qn)                  

            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j .eq. 1) utbt(nzn5:nzn6)=-utbt(nzn5:nzn6)
            if (j .eq. config%qx) uttp(nzn5:nzn6)=-uttp(nzn5:nzn6)
            
            do i = 1,3
               ugrid(nzn4+i) = 0.0_rk
               xl   (nzn4+i) = yznl  (qy_e+i)
               x    (nzn4+i) = yzn   (qy_e+i)
               xr   (nzn4+i) = yznr  (qy_e+i)
               dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
            end do
            
            do i = nzn8, nzn8
              rho (i) =  rho (nzn7)
              u   (i) = -u   (nzn7)
              ut  (i) =  ut  (nzn7)
              utt (i) =  utt (nzn7)
              utbt(i) =  utbt(nzn7)
              uttp(i) =  uttp(nzn7)
              utlt(i) =  utlt(nzn7)
              utrt(i) =  utrt(nzn7)
              e   (i) =  e   (nzn7)
!
! ...  take care !!! potential is defined on the zone interfaces
!
              grav(i) =  grav(nzn7-i-1)
              gamc(i) =  gamc(nzn7-i)
              game(i) =  game(nzn7-i)
              p   (i) =  p   (nzn7-i)
              dx  (i) =  dx  (nzn7-i)
!-PNS         ugrid (i) =-ugrid(nzn76-i)
!-PNS1         ugrid (i) = ugrid(nzn76-i)
              ugrid (i) = ugrid(nzn4)
              do n = 1, config%qn
                xn  (i,n) = xn (nzn7-i,n)
              enddo
            enddo

            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
              ! the axis through the origin in spherical
              ! polar coordinates (this is not an optimal
              ! solution)
              do i = nzn8, nzn8
                utt (i) = - utt (i)
                utrt(i) = - utrt(i)
                utlt(i) = - utlt(i)
              enddo
            endif
         else
            raise_abort("bndry(): error: outer bndry in case 590")
            !         call stopit_mpi("bndry: error: outer bndry")
         end if

         goto 600

590      continue 
         ! special case: two last zones reflective, two communication
         if (xyzswp .eq. 2) then
            
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (nzn5:nzn6)=densty_ub (j ,1:2,k)     
            u   (nzn5:nzn6)=vely_ub   (j ,1:2,k)
            ut  (nzn5:nzn6)=velx_ub   (j ,1:2,k)
            utt (nzn5:nzn6)=velz_ub   (j ,1:2,k)
            uttp(nzn5:nzn6)=vxold_ub  (i2,1:2,k)    
            utbt(nzn5:nzn6)=vxold_ub  (i1,1:2,k)    
            utrt(nzn5:nzn6)=vzold_ub  (j ,1:2,k2)    
            utlt(nzn5:nzn6)=vzold_ub  (j ,1:2,k1)    
            e   (nzn5:nzn6)=energy_ub (j ,1:2,k)
            p   (nzn5:nzn6)=press_ub  (j ,1:2,k)
            tmp (nzn5:nzn6)=temp_ub   (j ,1:2,k)
            grav(nzn5:nzn6)=gpocntr_ub(j ,1:2,k)
            game(nzn5:nzn6)=gammae_ub (j ,1:2,k)
            gamc(nzn5:nzn6)=gammac_ub (j ,1:2,k)
            xn  (nzn5:nzn6,1:config%qn) = xnuc_ub(j,1:2,k,1:config%qn)                  
            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j .eq. 1) utbt(nzn5:nzn6)=-utbt(nzn5:nzn6)
            if (j .eq. config%qx) uttp(nzn5:nzn6)=-uttp(nzn5:nzn6)
            
            do i = 1,2
               ugrid(nzn4+i) = 0.0_rk
               xl   (nzn4+i) = yznl  (qy_e+i)
               x    (nzn4+i) = yzn   (qy_e+i)
               xr   (nzn4+i) = yznr  (qy_e+i)
               dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
            end do
            
            !
            nzn76 = nzn7 + nzn6
!
            do i = nzn7, nzn8
              rho (i) =  rho (nzn76-i)
              u   (i) = -u   (nzn76-i)
              ut  (i) =  ut  (nzn76-i)
              utt (i) =  utt (nzn76-i)
              utbt(i) =  utbt(nzn76-i)
              uttp(i) =  uttp(nzn76-i)
              utlt(i) =  utlt(nzn76-i)
              utrt(i) =  utrt(nzn76-i)
              e   (i) =  e   (nzn76-i)
!
! ...  take care !!! potential is defined on the zone interfaces
!
              grav(i) =  grav(nzn76-i-1)
              gamc(i) =  gamc(nzn76-i)
              game(i) =  game(nzn76-i)
              p   (i) =  p   (nzn76-i)
              dx  (i) =  dx  (nzn76-i)
!-PNS         ugrid (i) =-ugrid(nzn76-i)
!-PNS1         ugrid (i) = ugrid(nzn76-i)
              ugrid (i) = ugrid(nzn4)
              do n = 1, config%qn
                xn  (i,n) = xn (nzn76-i,n)
              enddo
            enddo

            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
              ! the axis through the origin in spherical
              ! polar coordinates (this is not an optimal
              ! solution)
              do i = nzn7, nzn8
                utt (i) = - utt (i)
                utrt(i) = - utrt(i)
                utlt(i) = - utlt(i)
              enddo
            endif
         else
            raise_abort("bndry(): error: outer bndry")
            !         call stopit_mpi("bndry: error: outer bndry")
         end if

         goto 600

595      continue 
         ! special case: last three zones reflective, one communication

         ! special case: two last zones reflective, two communication

         !
         !     Carefull: this case implies that the higher boundary zone
         !               4 is the "reflection" of the zones qy_e-1 (i.e.
         !               _left_ of the current task nyproc-1
         !               thus the values from densty_lb(:,4,:) are needed for boundary
         !               zone 4. 
         if (xyzswp .eq. 2) then
            
            i1=max(j-1,1)
            i2=min(j+1,config%qx)
            k1=max(k-1,1)
            k2=min(k+1,config%qz)

            rho (nzn5:nzn5)=densty_ub (j ,1:1,k)     
            u   (nzn5:nzn5)=vely_ub   (j ,1:1,k)
            ut  (nzn5:nzn5)=velx_ub   (j ,1:1,k)
            utt (nzn5:nzn5)=velz_ub   (j ,1:1,k)
            uttp(nzn5:nzn5)=vxold_ub  (i2,1:1,k)    
            utbt(nzn5:nzn5)=vxold_ub  (i1,1:1,k)    
            utrt(nzn5:nzn5)=vzold_ub  (j ,1:1,k2)    
            utlt(nzn5:nzn5)=vzold_ub  (j ,1:1,k1)    
            e   (nzn5:nzn5)=energy_ub (j ,1:1,k)
            p   (nzn5:nzn5)=press_ub  (j ,1:1,k)
            tmp (nzn5:nzn5)=temp_ub   (j ,1:1,k)
            grav(nzn5:nzn5)=gpocntr_ub(j ,1:1,k)
            game(nzn5:nzn5)=gammae_ub (j ,1:1,k)
            gamc(nzn5:nzn5)=gammac_ub (j ,1:1,k)
            xn  (nzn5:nzn5,1:config%qn) = xnuc_ub(j,1:1,k,1:config%qn)                  
            !     Careful: boundary conditions for utbt and uttp (=vx_vold)
            !     are supposed to be ibndbt=1 and ibndtp=5:
            if (j .eq. 1) utbt(nzn5:nzn5)=-utbt(nzn5:nzn5)
            if (j .eq. config%qx) uttp(nzn5:nzn5)=-uttp(nzn5:nzn5)
            
            do i = 1,1
               ugrid(nzn4+i) = 0.0_rk
               xl   (nzn4+i) = yznl  (qy_e+i)
               x    (nzn4+i) = yzn   (qy_e+i)
               xr   (nzn4+i) = yznr  (qy_e+i)
               dx   (nzn4+i) = xr(nzn4+i) - xl(nzn4+i)
            end do
            
            do i = nzn6, nzn6
              rho (i) =  rho (nzn5)
              u   (i) = -u   (nzn5)
              ut  (i) =  ut  (nzn5)
              utt (i) =  utt (nzn5)
              utbt(i) =  utbt(nzn5)
              uttp(i) =  uttp(nzn5)
              utlt(i) =  utlt(nzn5)
              utrt(i) =  utrt(nzn5)
              e   (i) =  e   (nzn5)
!
! ...  take care !!! potential is defined on the zone interfaces
!
              grav(i) =  grav(nzn5-i-1)
              gamc(i) =  gamc(nzn5-i)
              game(i) =  game(nzn5-i)
              p   (i) =  p   (nzn5-i)
              dx  (i) =  dx  (nzn5-i)
!-PNS         ugrid (i) =-ugrid(nzn76-i)
!-PNS1         ugrid (i) = ugrid(nzn76-i)
              ugrid (i) = ugrid(nzn4)
              do n = 1, config%qn
                xn  (i,n) = xn (nzn5-i,n)
              enddo
            enddo
            
            do i = nzn7, nzn7
              rho (i) =  rho (nzn4)
              u   (i) = -u   (nzn4)
              ut  (i) =  ut  (nzn4)
              utt (i) =  utt (nzn4)
              utbt(i) =  utbt(nzn4)
              uttp(i) =  uttp(nzn4)
              utlt(i) =  utlt(nzn4)
              utrt(i) =  utrt(nzn4)
              e   (i) =  e   (nzn4)
!
! ...  take care !!! potential is defined on the zone interfaces
!
              grav(i) =  grav(nzn4-i-1)
              gamc(i) =  gamc(nzn4-i)
              game(i) =  game(nzn4-i)
              p   (i) =  p   (nzn4-i)
              dx  (i) =  dx  (nzn4-i)
!-PNS         ugrid (i) =-ugrid(nzn76-i)
!-PNS1         ugrid (i) = ugrid(nzn76-i)
              ugrid (i) = ugrid(nzn4)
              do n = 1, config%qn
                xn  (i,n) = xn (nzn4-i,n)
              enddo
            enddo

            do i = nzn8, nzn8
              rho (i) =  densty_lb (j, 1,k)
              u   (i) = -vely_lb   (j, 1,k)
              ut  (i) =  velx_lb   (j, 1,k)
              utt (i) =  velz_lb   (j, 1,k)
              utbt(i) =  vxold_lb  (i1,1,k)
              uttp(i) =  vxold_lb  (i2,1,k)
              utlt(i) =  vzold_lb  (j, 1,k1)
              utrt(i) =  vzold_lb  (j, 1,k2)
              e   (i) =  energy_lb (j, 1,k)
!
! ...  take care !!! potential is defined on the zone interfaces
!
              grav(i) =  gpocntr_lb(j, 1,k)
              gamc(i) =  gammac_lb (j, 1,k)
              game(i) =  gammae_lb (j, 1,k)
              p   (i) =  press_lb  (j, 1,k)
              dx  (i) =  dx  (nzn4-1)
!-PNS         ugrid (i) =-ugrid(nzn76-i)
!-PNS1         ugrid (i) = ugrid(nzn76-i)
              ugrid (i) = ugrid(nzn4)
              do n = 1, config%qn
                xn  (i,n) = xnuc_lb(j,1,k,n)
              enddo
            enddo

            if (xyzswp .ne. 3 .and. config%igeomz .eq. 5) then
              ! switch sign of v_phi when going across
              ! the axis through the origin in spherical
              ! polar coordinates (this is not an optimal
              ! solution)
              do i = nzn8, nzn8
                utt (i) = - utt (i)
                utrt(i) = - utrt(i)
                utlt(i) = - utlt(i)
              enddo
            endif
         else
            raise_abort("bndry(): error: outer bndry")
            !         call stopit_mpi("bndry: error: outer bndry")
         end if

         goto 600         


      return

600   continue

      do i = nzn5, nzn8
         xl(i) =      xl(i-1) + dx(i-1)
         xr(i) =      xr(i-1) + dx(i)
         x(i)  = 0.5_rk * (xl(i) + xr(i))
      enddo

      return
    end subroutine bndry

end module bndry_mod

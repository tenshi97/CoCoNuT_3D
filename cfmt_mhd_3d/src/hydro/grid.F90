module grid_mod

contains

!>
!> \verbatim
!> Compute grid coordinates (in theta-direction only equidistant
!> grid at the moment)
!>
!> \todo non-equidistang grid in theta-direction (2D/3D)
!>
!>  Author: W. Keil and M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine grid

  use precision

  use totare_hy     
  use intgrs_hy
  use phycon, only : pc_pi
  use gfloat_hy
! use intgrs_hy
! use mesh_hy
  use movgrid_hy

  use abort
  use mpi_vertex , only : myproc
  use print_stdout_mod
! LOCAL variables that are not in modules

  use configure
  implicit none

  integer(kind=ik) :: i
  real(kind=rk)    :: Rmax,a1,anglpm,dely,c2,a2,b2,x3
  real(kind=rk)    :: x2,x1,dxstr,delz,b1,angres
  real(kind=rk)    :: xrnew(config%qx)
      
#ifdef CFC_TRANSPORT
  real(kind=rk)    :: r1,alpha,beta
  integer(kind=ik) :: nn1,inum,idi
#endif

#if defined(NOMOTO_MODEL) || defined(INITIAL_HYDRO_GRID_FILE)
  integer(kind=ik) :: idum
#endif




! zero arrays first
  xzltot(:) = 0._rk
  xzntot(:) = 0._rk
  xzrtot(:) = 0._rk
  yzltot(:) = 0._rk
  yzntot(:) = 0._rk
  yzrtot(:) = 0._rk
  zzltot(:) = 0._rk
  zzntot(:) = 0._rk
  zzrtot(:) = 0._rk

!-----------------------------------------------------------------------
!     set up x-grid:
!-----------------------------------------------------------------------
  
  nzn  = config%qx
  nzn1 = nzn + 1
  nzn2 = nzn + 2
  nzn3 = nzn + 3
  nzn4 = nzn + 4
  nzn5 = nzn + 5
  nzn6 = nzn + 6
  nzn7 = nzn + 7
  nzn8 = nzn + 8
  
  rob   = config%gridlx


#ifdef CFC_TRANSPORT

!     Setup a hydro grid that can be used both for the collapse
!     and early post-bounce phase:

!     Innermost zone (r<10km), resolution increases with radius:

  r1=1.0e+6_rk                  !10km

  nn1=int(33.0_rk*real(config%qx,kind=rk)/400.0_rk)
  xzltot(1) = 0.0_rk
  xzrtot(1) = 3.5e+4_rk*400.0_rk/real(config%qx,kind=rk)
  
  alpha=log(r1/dble(nn1*xzrtot(1)))*2.0_rk/real(nn1-1,kind=rk)
  
  do i=2,nn1
     xzrtot(i)=xzrtot(1)*i*exp(0.5_rk*real(i-1,kind=rk)*alpha)
  end do

!     Intermediate region (10km<r<400km):
  alpha=(xzrtot(nn1)-xzrtot(nn1-1))/xzrtot(nn1-1)
  beta=0.2_rk*alpha
  alpha=0.8_rk*alpha

  do i=nn1+1,config%qx
     xzrtot(i)=xzrtot(nn1)*((1.0_rk+alpha)**(i-nn1)+ &
          beta*(i-nn1)*0.985_rk**(i-nn1-1))
     if (xzrtot(i).lt.4.1e+7_rk) inum=i
  end do
  
!     Outer region (r>400km):

  alpha=2.0_rk*(xzrtot(inum)-xzrtot(inum-1))/(xzrtot(inum)+xzrtot(inum-1))
  idi=config%qx-inum
  beta=log(alpha/((rob/xzrtot(inum))**(1.0_rk/real(idi,kind=rk))-1.0_rk))/ &
       log(real(idi,kind=rk))
  
  do i=inum+1,config%qx
     xzrtot(i)=xzrtot(inum)* &
          (1.0_rk+real(i-inum,kind=rk)**(-beta)*alpha)**(i-inum)
  end do

  do i=2,config%qx
     xzltot(i)=xzrtot(i-1)
  end do
  do i=1,config%qx
     xzntot(i)=0.5_rk*(xzrtot(i)+xzltot(i))
  end do

#if defined(NOMOTO_MODEL) || defined(INITIAL_HYDRO_GRID_FILE)
  open(1,file='gitter_cfc_hydro.dat',form='formatted')
  do i=1,config%qx
     read (1,*) idum, xzltot(i),xzrtot(i)
  end do
  close(1)
  
  do i=1,config%qx
     xzntot(i)=0.5*(xzrtot(i)+xzltot(i))
  end do
#endif /*NOMOTO_MODEL || INITIAL_HYDRO_GRID_FILE*/


#else /* CFC_TRANSPORT */

! ------------ grid ala Th.Zwerger

  if (config%rib .ne. 0._rk) then
     raise_abort("grid.F90(): rib .ne 0")
  endif

!      constants are set in input
!c --- set constants for moving grid:
!      nstr  = config%qx / 2
!c --- dx(1) ~ (idltn+1) * dx(2)
!      idltn  = 4
!      nstrt = nstr    + idltn
!      nxt   = config%qx   + idltn
!      xmfrac = 0.8*pc_ms   ! estimated mass of inner core
!
! --- rstr:
!      rstr = 5.1499799e7
!      rstr = 4.34185e+07  !s15s7b2 (Newton)
!      rstr = 4.01355E+07  !s15s7b2 (GR)
!      rstr = 3.55044e+07  !0.51 Msol; s10.8 (GR???) 
!      rstr =3.19161e+07  ! 0.64 Msol, where the shock forms in
!                         ! Mezzacappa, Liebendoerfer et al. simulation
!                         ! for nomoto 13 msol
!      rstr = 4.34185e+07 ! s1b4a (CHECK!!!) 


  Rmax = rob
      !
! --- PART I:
  a1     = log(Rmax / (Rmax - rstr)) / real(nstrt,kind=rk)
  b1     = log(Rmax) + real(nxt,kind=rk)/ &
                           real(nstrt,kind=rk)*log(1._rk - rstr/Rmax) 
  dxstr  = (Rmax - rstr) * (exp(a1) - 1._rk)
  
  do i = 1, nstr
     xrnew(i) = Rmax -  exp(a1*(config%qx-i) + b1)
  enddo
!
! --- PART II:
  x1     = log(rstr/xrnew(nstr-1)) / log(Rmax/rstr) 
  x2     = x1 * (config%qx**2 - nstr**2) + 1._rk - 2._rk* &
       real(nstr,kind=rk)
  x3     = 2._rk*( 1._rk - x1*(config%qx - nstr) )
  b2     = x2 / x3
  
  x2     = log( 1._rk - dxstr/rstr )
  x3     = 1._rk - 2._rk*( nstr + b2 )
  a2     = x2 / x3
  
  c2     = log( rstr ) - a2 * (nstr + b2)**2
  
  do i = nstr+1, config%qx
     xrnew(i) = exp(a2*(i+b2)**2 + c2) 
  enddo
      
  do i = 1, config%qx
     xzrtot(i) = xrnew(i)
  enddo
  xzltot(1) = 0._rk
  do i = 2, config%qx
     xzltot(i) = xzrtot(i-1)
  enddo
 

  do i = 1, config%qx
     xzntot(i) = 0.5_rk  * (xzrtot(i) + xzltot(i))
  enddo

#endif /* CFC_TRANSPORT */


  yzltot(1) = 0._rk
  yzntot(1) = 0._rk
  yzrtot(1) = 0._rk
  zzltot(1) = 0._rk
  zzntot(1) = 0._rk
  zzrtot(1) = 0._rk


!-----------------------------------------------------------------------
!     set up y-grid:
!-----------------------------------------------------------------------

  if (config%nsdim .ge. 2)  then
     
     nzn  = config%qy
     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8

     yzltot(1) = 0._rk
     config%gridly = pc_pi * config%gridly

     if (config%igeomy .ge. 3)  then
        call printit_taskX(0," ")
        call printit_taskX(0," ****  W A R N I N G  **** :")
        call printit_taskX(0,"    Angular grid is selected in y-direction")
        call printit_taskX(0,"    ===>  Grid length is reset to  GRIDLY * PI")
        call printit_taskX(0," ")

        dely   = config%gridly / config%qy
        !
        if (config%bndmny.eq.4 .and. config%igeomy.eq.4)  then
           yzltot(1) = pc_pi / 2._rk  -  config%gridly / 2._rk
           anglpm    = 90._rk  -  yzltot(1) * 180._rk / pc_pi
           angres    = dely * 180._rk / pc_pi
           call printit_taskX(0," ")
           call printit_taskX(0," ****  N O T E  **** :")
           call printit_taskX(0,"    Periodic boundary conditions are used at")
           call printit_taskX(0,"THETA = 90.0 - ",anglpm)
           call printit_taskX(0,"   and   THETA = 90.0 + ",anglpm)
           call printit_taskX(0,"degrees")
           call printit_taskX(0,"Angular resolution is [degrees]",angres)
           call printit_taskX(0," ")
        end if

        if (config%bndmny.eq.4 .and. (config%igeomy.eq.3 .or. config%igeomy.eq.5) )  then
           yzltot(1) = - config%gridly / 2._rk
           anglpm    = - yzltot(1) * 180._rk / pc_pi
           angres    = dely * 180._rk / pc_pi
           call printit_taskX(0," ")
           call printit_taskX(0," ****  N O T E  **** :")
           call printit_taskX(0,"    Periodic boundary conditions are used at")
           call printit_taskX(0,"PHI = 90.0 - ",anglpm)
           call printit_taskX(0,"   and   PHI = 90.0 + ",anglpm)
           call printit_taskX(0,"degrees")
           call printit_taskX(0,"Angular resolution is [degrees]",angres)
           call printit_taskX(0," ")
        end if
        
        if (config%bndmny.eq.1 .and. config%igeomy.eq.4)  then
           if (config%isym .eq. 0) then
              yzltot(1) = pc_pi / 2._rk  -  config%gridly / 2._rk
           else
              yzltot(1) = pc_pi / 2._rk  -  config%gridly 
           end if
           anglpm  = 90._rk  -  yzltot(1) * 180._rk / pc_pi
           angres  = dely * 180._rk / pc_pi
           call printit_taskX(0," ")
           call printit_taskX(0," ****  N O T E  **** :")
           call printit_taskX(0,"    Reflecting boundary conditions are used at", anglpm)
           call printit_taskX(0,"    and ",anglpm)
           call printit_taskX(0,"Angular resolution is [degrees]",angres)
           call printit_taskX(0," ")
           
        end if
        
     end if

#ifndef VAR_THETA_GRID
     dely    = config%gridly / config%qy
#else /* VAR_THETA_GRID */
     dely    = 2._rk / real(config%qy,kind=rk)
#endif /* VAR_THETA_GRID */

#ifndef VAR_THETA_GRID
     do i = 2, config%qy + 1
        yzltot(i) = yzltot(1) + real((i-1),kind=rk)*dely
     enddo
#else /* VAR_THETA_GRID */
     if (config%gridly/pc_pi .ne. 1._rk) then
      raise_abort("grid.F90(): gridly .ne. 1 and VAR_THETA_GRID activated")
     endif

     do i = 1, config%qy + 1
        yzltot(i) = acos(1._rk-real(i-1,kind=rk)*dely)
     enddo
#endif /* VAR_THETA_GRID */

     do i = 1, config%qy
        yzrtot(i) = yzltot(i+1)
     enddo
     
     do i = 1, config%qy
        yzntot(i) = 0.5_rk * (yzrtot(i) + yzltot(i))
     enddo
     
     zzltot(1) = 0._rk
     zzntot(1) = 0._rk
     zzrtot(1) = 0._rk
     

     call printit_taskX(0,"pns_grid> y-grid: i/yznl/ yzn/ yznr | deg: yznl/ yzn/ yznr")

     do i = 1, config%qy
        if (myproc .eq. 0) then
           write(*,'('' Task '',i4,i3,3(1x,1pe12.5),''|'',3(1x,1pe12.5) )')  &
                myproc, i,yzltot(i),yzntot(i),yzrtot(i), &
                yzltot(i)*180./pc_pi,yzntot(i)*180./pc_pi,yzrtot(i)*180./pc_pi 
        endif
     enddo
     
     call printit_taskX(0," ")
  endif

!-----------------------------------------------------------------------
!     set up z-grid:
!-----------------------------------------------------------------------

  if (config%nsdim .ge. 3)  then
     
     nzn  = config%qz
     nzn1 = nzn + 1
     nzn2 = nzn + 2
     nzn3 = nzn + 3
     nzn4 = nzn + 4
     nzn5 = nzn + 5
     nzn6 = nzn + 6
     nzn7 = nzn + 7
     nzn8 = nzn + 8
     
     delz    = config%gridlz / config%qz
     zzltot(1) = 0._rk
     
     if (config%igeomz .ge. 3)  then
        call printit_taskX(0," ")
        call printit_taskX(0," ****  W A R N I N G  **** :")
        call printit_taskX(0,"    Angular grid is selected in z-direction")
        call printit_taskX(0,"    ===>  Grid length is reset to  GRIDLZ * PI")
        call printit_taskX(0," ")      

        config%gridlz = pc_pi * config%gridlz
        delz   = config%gridlz / config%qz
        
        if (config%bndmnz.eq.4 .and. config%igeomz.eq.4)  then
           zzltot(1) = pc_pi / 2._rk  -  config%gridlz / 2._rk
           anglpm    = 90._rk  -  zzltot(1) * 180._rk / pc_pi
           angres    = delz * 180._rk / pc_pi
           call printit_taskX(0," ")
           call printit_taskX(0," ****  N O T E  **** ")
           call printit_taskX(0,"    Periodic boundary conditions are used at",anglpm)
           call printit_taskX(0,"    and ",anglpm)
           call printit_taskX(0,"    Angular resolution: ",angres)
           call printit_taskX(0," ")
        end if
        
        if (config%bndmnz.eq.4 .and. (config%igeomz.eq.3 .or. config%igeomz.eq.5) )  then
           zzltot(1) = - config%gridlz / 2._rk
           anglpm    = - zzltot(1) * 180._rk / pc_pi
           angres    = delz * 180._rk / pc_pi
           call printit_taskX(0," ")
           call printit_taskX(0," ****  N O T E  **** ")
           call printit_taskX(0,"    Periodic boundary conditions are used at",anglpm)
           call printit_taskX(0,"    and ",anglpm)
           call printit_taskX(0,"    Angular resolution: ",angres)
           call printit_taskX(0," ")
        end if
     end if
     
     
     do i = 2, config%qz + 1
        zzltot(i) = zzltot(i-1) + delz
     enddo
     
     do i = 1, config%qz
        zzrtot(i) = zzltot(i+1)
     enddo
     
     do i = 1, config%qz
        zzntot(i) = 0.5_rk * (zzrtot(i) + zzltot(i))
     enddo
     
     call printit_taskX(0,"ns_grid> z-grid: i/zznl/ zzn/ zznr | deg: zznl/ zzn/ zznr")
     
     do i = 1, config%qz
        if (myproc .eq. 0) then
           write(*,*) "Task ",myproc
           write(*,'('' Task '',i4,i3,3(1x,1pe12.5),''|'',3(1x,1pe12.5))')  &
                     myproc,i,zzltot(i),zzntot(i),zzrtot(i), &
                     zzltot(i)*180./pc_pi,zzntot(i)*180./pc_pi,zzrtot(i)*180./pc_pi 
        endif
     enddo
     call printit_taskX(0," ")
  endif
  
  return
end subroutine grid

#ifndef PROGRAM_remap
!> Set grid velocity in radial direction for hydro grid
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
subroutine grdvel(selftime, childrentime)

  use precision

!      use intgrs_hy      ! forcheck
  use totare_hy
  use param_rt

  use mo_mpi

#ifndef DEBUG_TIMINGS
  use cputim
#endif

  use configure

  implicit none
  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2),           &
                                time_communication_start(2), &
                                time_communication_end(2)

  integer(kind=ik)           :: ii,i,k
  real(kind=rk)              :: xi
  real(kind=rk)              :: ugrtot_tmp(config%qx)
  real(kind=rk)              :: ugrtot_tmp_buf(config%qx)
  integer(kind=ik)           :: ierr


  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
     call second_v(selftime_start)
#endif

! save old grid-positions 
  xzrold(0)       = xzltot(1)
  xzrold(1:config%qx) = xzrtot(1:config%qx)


  if (config%laghyd .eq. 1) then
! ------------------------  Quasi-Langrange (1D only):

     ugrtot_tmp(:) = 0._rk
        do k = qz_s,qz_e
           do i = qy_s, qy_e
              ugrtot_tmp(:) = ugrtot_tmp(:) + vextot(:,i,k)
           end do
        end do


     if (use_mpi) then
        ! MPI_allreduce (sum), urgtot(:)
#ifndef DEBUG_TIMINGS
        call second_v(time_communication_start)
#endif
        call MPI_allreduce(ugrtot_tmp, ugrtot_tmp_buf, config%qx, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        ugrtot_tmp(:)=ugrtot_tmp_buf(:)

#ifndef DEBUG_TIMINGS        
        call second_v(time_communication_end)
        timer%hydro_comm = timer%hydro_comm + (time_communication_end-time_communication_start)
        childrentime     = childrentime + (time_communication_end-time_communication_start)
#endif
        
     endif


        ugrtot(1:config%qx) = ugrtot_tmp(1:config%qx)/ &
             (real(config%qy,kind=rk)*real(config%qz,kind=rk))
     

     ugrtot(config%qx)=0._rk
     
     ii=5
     do i=ii,1,-1
        xi=xzntot(i)/xzntot(ii)
        ugrtot(i)=xi*ugrtot(ii)
     enddo
  else
! ------------------------- Eulerian grid
     ugrtot(1:config%qx)  = 0._rk
  endif

#ifdef FIX_OBH
  vextot(config%qx,1,1)=-1.e-9_rk
#endif

#ifndef DEBUG_TIMINGS        
        call second_v(selftime)
        selftime = selftime - selftime_start
#endif
  return

end subroutine grdvel
#endif /* PROGRAM_remap */

end module grid_mod

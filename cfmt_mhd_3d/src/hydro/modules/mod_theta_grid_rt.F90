!>
!> \par This module provides some variables needed for the lateral transport grid
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module theta_grid_rt
  use precision

  implicit none

  save

! LOCAL variables that are not in modules


!  real(kind=rk),save, dimension(0:nymom+1) :: th,thl,thr,cthl,cthr, &
!                                                  sthl,sthr

  real(kind=rk), allocatable, dimension(:) :: th, thl, thr, cthl, cthr, &
                                              sthl, sthr

contains 
!>
!> \par This subroutine sets up a equidistant lateral transport grid
!>
!> \todo non-equidistant grid
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine init_theta_grid
    use precision
    use abort
#ifndef NOTRA
    use multigrid_rt
#endif
    use totare_hy

    use configure
    implicit none
! LOCAL variables that are not in modules

    integer(kind=ik) :: j
    do j=1,config%nymom
       thl(j)=yzltot(jmin(j))
       thr(j)=yzrtot(jmax(j))
    enddo

    select case(config%latybc)
    case(4)
! periodic boundary
       thr(0)=thl(1)
       thl(0)=thr(0)-(thr(config%nymom)-thl(config%nymom))

       thl(config%nymom+1)=thr(config%nymom)
       thr(config%nymom+1)=thr(config%nymom)+(thr(1)-thl(1))

    case(1)
! reflecting boundary
       thr(0)=thl(1)
       thl(0)=thr(0)-(thr(1)-thl(1))
       
       thl(config%nymom+1)=thr(config%nymom)
       thr(config%nymom+1)=thr(config%nymom)+(thr(config%nymom)-thl(config%nymom))
    case default
       raise_abort("init_theta_grid(): case not implemented")
    end select
    
    do j=0,config%nymom+1
       th(j)=0.5_rk*(thl(j)+thr(j))
       cthl(j)=cos(thl(j))
       cthr(j)=cos(thr(j))
       sthl(j)=sin(thl(j))
       sthr(j)=sin(thr(j))
    enddo

  end subroutine init_theta_grid

!>
!> \par This subroutine allocates the arrays from module theta_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine allocate_theta_grid_rt(mem_global)
      use precision
      use abort
      use print_stdout_mod, only : print_memory_alloc
      use configure
      implicit none
      integer(kind=ik) :: istat
      real(kind=rk)    :: mem_global, mem_local

      allocate(th(0:config%nymom+1), thl(0:config%nymom+1), thr(0:config%nymom+1), cthl(0:config%nymom+1), &
               cthr(0:config%nymom+1), sthl(0:config%nymom+1), sthr(0:config%nymom+1), stat=istat)
      if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): allocation of module theta_grid_rt failed")
      end if
      
      mem_local = (config%nymom+2)*8._rk * 7._rk

      mem_local = mem_local/1024._rk/1024._rk

      mem_global = mem_global + mem_local

      call print_memory_alloc(mem_local, mem_global, "theta_grid_rt")

  end subroutine allocate_theta_grid_rt
!>
!> \par This subroutine deallocates the arrays from module theta_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine deallocate_theta_grid_rt
      use precision
      use abort
      implicit none
      integer(kind=ik) :: istat
      
      deallocate(th,thl,thr,cthl,cthr,sthl,sthr, stat=istat)
      if (istat .ne. 0) then
         raise_abort("allocate_hydro_arrays(): deallocation of module theta_grid_rt failed")
      end if


  end subroutine deallocate_theta_grid_rt


end module theta_grid_rt

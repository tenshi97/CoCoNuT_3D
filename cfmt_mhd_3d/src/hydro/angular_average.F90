#undef DIAGNOSE_2D

!> \verbatim
!> This module provides the subroutines avecoeff_mass and avecoeff_vol 
!> which are overloaded with the CFC or the PROMETHEUS version
!>
!> Here is a list of subroutines which are overloaded depending on whether
!> CFC or PROMETHEUS is uses
!>
!>  interface-name    real-name               compiled when?
!>  avecoeff_mass     avecoeff_mass_CFC       CFC_TRANSPORT
!>  avecoeff_vol      avecoeff_vol_CFC        CFC_TRANSPORT
!>  avecoeff_mass     avecoeff_mass_PROM      if not CFC_TRANSPORT
!>  avecoeff_vol      avecoeff_vol_PROM       if not CFC_TRANSPORT
!>
!>  The purpose of avecoeff_mass and avecoeff_vol is to provide 
!>  coefficients for mass-weighted angular averages
!              required e.g. for v,Y_e
!>
!>  Author: A. Marek, MPA, Nov. 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module angular_average

  use precision

  public avecoeff_mass, avecoeff_vol


  interface avecoeff_mass

#ifdef CFC_TRANSPORT2
     module procedure avecoeff_mass_CFC
#else
     module procedure avecoeff_mass_PROM
#endif

  end interface


  interface avecoeff_vol

#ifdef CFC_TRANSPORT2
     module procedure avecoeff_vol_CFC
#else
     module procedure avecoeff_vol_PROM
#endif
  end interface

contains


#ifdef CFC_TRANSPORT2
!--------------------------------------------------------------------
!>
!> \verbatim
!> Purpose : provides coefficients for mass-weighted angular averages
!>              required e.g. for v,Y_e
!>
!>  Author: M. Rampp and Bernhard Mueller
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  SUBROUTINE avecoeff_mass_CFC(mode, selftime, childrentime)


    use precision
    use abort
#ifndef NOTRA

    use multigrid_rt
#endif

    use intgrs_hy
    use totare_hy, only: dvytot, dvztot
    
    use size_cfc
    use grid_cfc
    use conserved_cfc
    use metric_cfc

    use mo_mpi
    use cputim
    use configure
    IMPLICIT NONE

! LOCAL variables that are not in modules
    real(kind=rk), intent(out)    :: selftime(2), childrentime(2)
    real(kind=rk)                 :: selftime_start(2)
    real(kind=rk)                 :: paralleltime_start(2), paralleltime_end(2)
    integer(kind=ik), intent (in) :: mode
    
    real(kind=rk)                 :: sumass(config%qx),  &
                                     secmas(config%qx,config%qy,config%qz)
    
    integer(kind=ik)              :: i, j, k, jk
    integer(kind=ik)              :: ierr

    selftime              = 0._rk
    childrentime          = 0._rk
    paralleltime_start    = 0._rk
    paralleltime_end      = 0_rk
#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif
    
    if (config%nsdim .ge. 2) then
       if (config%igeomy .ne. 4) stop 'inigrid_th> igeomy .ne. 4 '
       
!-- compute integration-weights for conservatively theta-averaging the 
!     hydro- onto the TRANSP-Grid
!     general ! 

       sumass(:) = 0.0_rk

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
       call second_v(paralleltime_start)
#endif
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(sumass,secmas) &
!$OMP DEFAULT(shared)

!$OMP DO REDUCTION(+:sumass)
#endif
       do jk = 1, n_loc * o_loc
          
          k = int((jk + n_loc - 1) / n_loc )
          j = (n_s - 1) + (jk - (k - 1) * n_loc)
          k = k + o_s - 1

          do i=1,config%qx
             thphmass(i,j,k,0)=1.0_rk
             secmas(i,j,k) = d_cap_hat(i,j,k) * dvytot(j) * dvztot(k)
             sumass(i) = sumass(i) + secmas(i,j,k)
          enddo
       end do
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
       call second_v(paralleltime_end)
#endif
#endif

       timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start
#ifdef MPI_HYDRO
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
          call MPI_Allreduce (MPI_IN_PLACE,sumass,config%qx,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
#endif /* MPI_HYDRO */

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
          call second_v(paralleltime_start)
#endif
!!$OMP DO
#endif
       do jk = 1, n_loc * o_loc

          k = int((jk + n_loc - 1) / n_loc )
          j = (n_s - 1) + (jk - (k - 1) * n_loc)
          k = k + o_s - 1
       
          do i=1,config%qx
             thphmass(i,j,k,0) = secmas(i,j,k)/sumass(i)
          enddo
       enddo
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP END PARALLEL
#ifndef DEBUG_TIMINGS
       call second_v(paralleltime_end)
#endif
#endif
       timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start
    else
       thphmass(:,1,1,1)=1._rk
    endif
#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif    
    return
    
  END SUBROUTINE avecoeff_mass_CFC


!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>
!> \verbatim
!> Purpose : provides coefficients for volume-weighted angular averages
!>              required e.g. for rho,etot
!>
!>  Author: M. Rampp and Bernhard Mueller
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  SUBROUTINE avecoeff_vol_CFC(mode)

    use precision
    use abort
#ifndef NOTRA
    use multigrid_rt
#endif
    use totare_hy

    use intgrs_hy
    
    use size_cfc
    use grid_cfc
    use metric_cfc

    use mo_mpi

    use configure
    IMPLICIT NONE
! LOCAL variables that are not in modules

    integer(kind=ik), intent(in) :: mode
    
    real(kind=rk) :: divi(1:config%qx),divi1(1:config%qx),thf(1:config%qx),div(1:config%qx)
    integer(kind=ik) :: i,j,k,jm
    integer(kind=ik) :: ierr
    
    if (config%nsdim.eq.3) stop 'inigrid_th>3D not implemented'
    if (config%nsdim.eq.2 .and. config%igeomy .ne. 4)  stop 'inigrid_th> igeomy .ne. 4 '
    
    do j=n_s,n_e
       if (config%nsdim .ge. 2) then
          sthe(j) = sin(yzntot(j))
       else
          sthe(j) = 1.0_rk
       endif
    enddo

    if (config%nsdim .ge. 2) then
       
!-- compute intergration-weights for conservatively theta-averaging the 
!     hydro- onto the TRANSP-Grid
!     general ! 


       if (config%nymom .eq. config%qy) then
          thphmean(:,:,:,1)=1.0_rk
          thphmean(:,:,:,0)=0.0_rk
       else if (config%nymom .eq. 1) then
          divi=0.0_rk
          do k=o_s,o_e
             do j=n_e,n_e
                do i=1,config%qx
                   divi(i) = divi(i)+sqrt_gamma(i,j,k)*dvytot(j)*dvztot(k)
                end do
             end do
          end do

#ifdef MPI_HYDRO
          call MPI_Allreduce(divi,divi1,config%qx,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
          divi=divi1
#endif /* MPI_HYDRO */

          divi(:)=1.0_rk/divi(:)
          
          do k=o_s,o_e
             do j=n_s,n_e
                thphmean(1:config%qx,jm,k,1) = &
                     sqrt_gamma(1:config%qx,jm,k)*dvytot(j)*dvztot(k)/divi(1:config%qx)
                thphmean(1:config%qx,jm,k,0) = thphmean(1:config%qx,jm,k,1)
             enddo
          enddo
       else
          raise_abort("Only nymom=qy or nymom=1 allowed!")
       end if
    else
       thphmean(:,:,:,1)=1._rk
       !         thphmean(1,0)=1.
       jmin(1)=1
       jmin(1)=1
    endif
    
  END SUBROUTINE avecoeff_vol_CFC

#else /* CFC_TRANSPORT2 */

!>
!> \verbatim
!> Purpose : provides coefficients for mass-weighted angular averages
!>              required e.g. for v,Y_e
!>
!>  Author: M. Rampp and F. Hanke
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine avecoeff_mass_PROM(mode, selftime, childrentime)

    use precision
#ifndef NOTRA    
!      use overlap         ! forcheck
    use multigrid_rt
#endif
    use intgrs_hy
!      use mesh_hy
!      use vnew_hy
    use totare_hy
    use abort

    use mo_mpi
    use mpi_vertex

! LOCAL variables that are not in modules
    use cputim
    use configure

    implicit none
    real(kind=rk), intent(out) :: selftime(2), childrentime(2)
    real(kind=rk)              :: selftime_start(2)
    integer(kind=ik)           :: jm, km, i, j, k ,mode
    real(kind=rk)              :: fac, div
#ifdef DIAGNOSE_2D
    real(kind=rk)              :: summe, summa
#endif
    real(kind=rk)              :: sumass(config%qx,config%qz),   &
                                  secmas(config%qx,config%qy,config%qz)

    selftime     = 0._rk
    childrentime = 0._rk
#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif

    if (config%nsdim .ge. 2) then
       if (config%igeomy .ne. 4) then
          raise_abort("inigrid_th(): igeomy .ne. 4")
       endif
       
!-- compute integration-weights for conservatively theta-phi-averaging the 
!     hydro- onto the TRANSP-Grid
!     general ! 

       sumass(:,:) = 0.0_rk

       if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
          do k=nzmoms,nzmome
             do j=nymoms,nymome
                do i=1,config%qx
                   secmas(i,j,k) = 0.0_rk
                   do jm=qy_s,qy_e
                      do km=qz_s,qz_e
                         secmas(i,j,k) = secmas(i,j,k) + dentot(i,jm,km) &
                                         * dvytot(jm) * dvztot(km)
                      enddo
                   enddo
                   sumass(i,k) = sumass(i,k) + secmas(i,j,k)
                   div = 1.0_rk / secmas(i,j,k)
                   do jm=qy_s,qy_e
                      do km=qz_s,qz_e
                         thphmass(i,jm,km,1) = dentot(i,jm,km)*dvytot(jm)*dvztot(km) * div
                      enddo
                   enddo
                enddo
             enddo

             do j=nymoms,nymome

                do i=1,config%qx
                   fac = secmas(i,j,k)/sumass(i,k)
                   do jm=qy_s,qy_e
                      do km=qz_s,qz_e
                         thphmass(i,jm,km,0) = thphmass(i,jm,km,1)*fac
                      enddo
                   enddo
                enddo
             enddo
          enddo

       else ! mpi edd
          do k=nzmoms,nzmome
             do j=nymoms,nymome
                do i=1,config%qx
                   secmas(i,j,k) = 0.0_rk
                   do jm=jmin(j),jmax(j)
                      do km=kmin(k),kmax(k)
                         secmas(i,j,k) = secmas(i,j,k) + dentot(i,jm,km) &
                                         * dvytot(jm) * dvztot(km)
                      enddo
                   enddo
                   sumass(i,k) = sumass(i,k) + secmas(i,j,k)
                   div = 1.0_rk / secmas(i,j,k)
                   do jm=jmin(j),jmax(j)
                      do km=kmin(k),kmax(k)
                         thphmass(i,jm,km,1) = dentot(i,jm,km)*dvytot(jm)*dvztot(km) * div
                      enddo
                   enddo
                enddo
             enddo
             
             do j=nymoms,nymome

                do i=1,config%qx
                   fac = secmas(i,j,k)/sumass(i,k)

                   do jm=jmin(j),jmax(j)
                      do km=kmin(k),kmax(k)
                         thphmass(i,jm,km,0) = thphmass(i,jm,km,1)*fac
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif ! mpi and edd fact
    else
       thphmass(:,1,1,1)=1.0_rk
    endif ! config%nsdim ge 2

#ifdef DIAGNOSE_2D
     if (config%nsdim .ge. 2) then
        km=1
        if (myproc .eq. 0) then
           write(*,*) 'transp_grid> th-grid: i | j | jm | jmin | jmax | ', &
               'yznl | yzn | yznr | thphmean_1 | thphmean_0 | ',        &
               'thphmass_1 | thphmass_0 |'
        endif

        do i=1,config%qx,10
#ifndef MPI_HYDRO
           do j=config%nystrt,config%nymom
#else
           do j=nymoms,nymome
#endif
              if (j .eq. 0) jat=0
              summe = zero
              summa = zero
              do jm=jmin(j),jmax(j)
                 if (myproc .eq. 0) then
                    write(*,*) "Task ",myproc
                    write(*,'(12x,I4,2I3,2I5,7(1x,1pe12.5))') i,j,jm,jmin(j), &
                             jmax(j), yzltot(jm),yzntot(jm),yzrtot(jm),       &
                              thphmean(jm,km,1),thphmean(jm,km,0),thphmass(i,jm,km,1),      &
                              thphmass(i,jm,km,0)
                 endif
                 summe = summe + thphmean(jm,1)
                 summa = summa + thphmass(i,jm,km,1)
              enddo
              if (myproc .eq. 0) then
                 write(*,*) "Task ",myproc
                 write(*,'(I3,2e12.5)') j,summe,summa
              endif
           do i=1,config%qx,10
              do j=config%nystrt,config%nymom
                 if (j .eq. 0) jat=0
                 summe = zero
                 summa = zero
                 do jm=jmin(j),jmax(j)
                    write(*,'(12x,I4,2I3,2I5,7(1x,1pe12.5))') i,j,jm,jmin(j), &
                         jmax(j), yzltot(jm),yzntot(jm),yzrtot(jm),       &
                         thphmean(jm,km,1),thphmean(jm,km,0),thphmass(i,jm,km,1),      &
                         thphmass(i,jm,km,0)
                    summe = summe + thphmean(jm,1)
                    summa = summa + thphmass(i,jm,km,1)
                 enddo
                 write(*,'(I3,2e12.5)') j,summe,summa
              enddo
           enddo
     endif


#endif /* DIAGNOSE 2D */

#ifndef DEBUG_TIMINGS
     call second_v(selftime)
     selftime = selftime - selftime_start
#endif

     return
   end subroutine avecoeff_mass_PROM
!>
!> \verbatim
!> Purpose : provides coefficients for volume-weighted angular averages
!>              required e.g. for rho,etot
!>
!>  Author: M. Rampp and F. Hanke
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
   subroutine avecoeff_vol_PROM(mode)

     use precision
#ifndef NOTRA
!      use intgrs_hy  ! foecheck
!      use overlap    ! forcheck
     use multigrid_rt
#endif

!      use mesh_hy
     use intgrs_hy
     use totare_hy
     use abort
     use mo_mpi

     use configure
! LOCAL variables that are not in modules
     implicit none
     integer(kind=ik) :: jm,km,mode,j,k
     real(kind=rk) :: div,thf,divi

     if (config%nsdim.eq.3 .and. config%igeomy .ne. 4 .and. config%igeomz .ne. 5) then
        raise_abort("inigrid_th(): igeomy .ne. 4 and igeomz .ne. 5")
     endif
      
     if (config%nsdim.eq.2 .and. config%igeomy .ne. 4) then
        raise_abort("inigrid_th(): igeomy .ne. 4")
     endif

!-- compute intergration-weights for conservatively theta/phi-averaging the 
!     hydro- onto the TRANSP-Grid
!     general ! 
     if (config%nsdim .eq. 3) then

        divi = 1.0_rk / ((cos(yzltot(1)) - cos(yzrtot(config%qy)) ) &
             *(zzrtot(config%qz)-zzltot(1)))

        if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then


           do k=nzmoms,nzmome
              do j=nymoms,nymome

                 thf = ( cos( yzltot(jmin(j)) ) - cos( yzrtot(jmax(j)) ) ) * &
                       ( zzrtot(kmax(k)) - zzltot(kmin(k)) ) * divi
                 div = 1.0_rk / (( cos( yzltot(jmin(j) )) - &
                       cos( yzrtot(jmax(j)) ) ) * ( zzrtot(kmax(k)) &
                          - zzltot(kmin(k)) ))
                 do jm=qy_s,qy_e
                    do km=qz_s,qz_e
                       thphmean(jm,km,1) = dvytot(jm)*dvztot(km)*div
                       thphmean(jm,km,0) = thphmean(jm,km,1) * thf
                    enddo
                 enddo
              enddo
           enddo



        else ! either mpi is not set or not use of spherical eddington factor
           
           do k=nzmoms,nzmome
              do j=nymoms,nymome

                 thf = ( cos( yzltot(jmin(j)) ) - cos( yzrtot(jmax(j)) ) ) * &
                       ( zzrtot(kmax(k)) - zzltot(kmin(k)) ) * divi
                 div = 1.0_rk / (( cos( yzltot(jmin(j) )) - &
                       cos( yzrtot(jmax(j)) ) ) * ( zzrtot(kmax(k)) &
                          - zzltot(kmin(k)) ))
                 do jm=jmin(j),jmax(j)
                    do km=kmin(k),kmax(k)
                       thphmean(jm,km,1) = dvytot(jm)*dvztot(km)*div
                       thphmean(jm,km,0) = thphmean(jm,km,1) * thf
                    enddo
                 enddo
              enddo
           enddo
           

        endif ! mpi
        
        if (config%nymom .eq. config%qy .and. config%nztra .eq. config%qz) then
           thphmean(:,:,1)=1.0_rk
        endif

     else if (config%nsdim .eq. 2) then
        km = 1
        divi = 1.0_rk / ( cos(yzltot(1)) - cos(yzrtot(config%qy)) )

        do j=nymoms,nymome

           thf = ( cos( yzltot(jmin(j)) ) - cos( yzrtot(jmax(j)) ) ) * divi
           div = 1.0_rk / ( cos( yzltot(jmin(j)) ) - cos( yzrtot(jmax(j)) ) )
           do jm=jmin(j),jmax(j)
              thphmean(jm,km,1) = dvytot(jm)* div
              thphmean(jm,km,0) = thphmean(jm,km,1) * thf
           enddo
        enddo
        
        if (config%nymom .eq. config%qy) thphmean(:,:,1)=1.0_rk
     else
        thphmean(1,1,1)=1.0_rk
        !         thphmean(1,1,0)=1._rk
        jmin(1)=1
        jmin(1)=1
     endif

#ifdef DIAGNOSE_2D
     if (config%nsdim .ge. 2) then
        if (myproc .eq. 0) then
           write(*,*) 'transp_grid> th-grid: j | jm | jmin | jmax | ', &
                  'yznl | yzn | yznr | thphmean_1 | thphmean_0 |'
        endif

        do j=nystrt,config%nymom
           do jm=jmin(j),jmax(j)
              if (myproc .eq. 0) then
                 write(*,*) "Task ",myproc
                 write(*,'(12x,2I3,2I5,5(1x,1pe12.5))') j,jm,jmin(j), &
                        jmax(j), yzltot(jm),yzntot(jm),yzrtot(jm),    &
                        thphmean(jm,1,1),thphmean(jm,1,0)
              endif
           enddo
        enddo
     endif
#endif 
!     stop 'inigrid_th> CHECK'

   end subroutine avecoeff_vol_PROM

#endif /* CFC_TRANSPORT */


 end module angular_average

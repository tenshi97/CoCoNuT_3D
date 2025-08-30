module mapare_proc

  use precision
  use abort

  implicit none

  public mapare, allocate_mapare, deallocate_mapare
  private

  integer(kind=ik) :: istat

  real(kind=rk), allocatable :: den_ave_buf(:), vex_ave_buf(:), &
                                vey_ave_buf(:), vez_ave_buf(:), &
                                tem_ave_buf(:), ene_ave_buf(:), &
                                xnuc_ave_buf(:,:)
  
contains

  subroutine allocate_mapare(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    real(kind=rk)    :: mem_global, mem_local

    allocate(den_ave_buf(config%qx), vex_ave_buf(config%qx), &
             vey_ave_buf(config%qx), vez_ave_buf(config%qx), &
             tem_ave_buf(config%qx), ene_ave_buf(config%qx), &
             xnuc_ave_buf(config%qx,1:config%qn), stat=istat)

    if (istat .ne. 0) then
       raise_abort("mapare.F90(): error in allocating arrays")
    endif
    
    mem_local = config%qx * 8._rk * 6._rk + config%qx * (config%qn) * 8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local
    
    call print_memory_alloc(mem_local, mem_global, "mapare")
    
  end subroutine allocate_mapare


  subroutine deallocate_mapare
    use precision
    use abort

    implicit none

    deallocate(den_ave_buf, vex_ave_buf, &
             vey_ave_buf, vez_ave_buf, &
             tem_ave_buf, ene_ave_buf, &
             xnuc_ave_buf, stat=istat)

    if (istat .ne. 0) then
       raise_abort("mapare.F90(): error in deallocating arrays")
    endif

  end subroutine deallocate_mapare

  subroutine mapare
    use precision
    use abort
    use vnew_hy
    use mesh_hy
    use totare_hy

    use intgrs_hy

#ifndef PROGRAM_remap
    use filare_overload
#endif

    use mo_mpi

    use eos3d_routine, only : eos3d

#ifndef DEBUG_TIMINGS
    use cputim
#endif
    use cpyare_mod
 
    use hydro_areas_mod
    use configure
    implicit none

! LOCAL varibales that are not in modules
    logical                             :: ler
    integer(kind=ik)                    :: l, ia
    integer(kind=ik)                    :: ixi, ixf, iox, iyi, &
                                           iyf, ioy, izi, izf, &
                                           ioz 

    real(kind=rk), dimension(config%qx) :: den_ave, vex_ave,   &
                                           vey_ave, vez_ave,   &
                                           tem_ave, ene_ave 
    real(kind=rk)                       :: xnuc_ave(config%qx,config%qn)

    integer(kind=ik)                    :: i, j, k, ierr

    real(kind=rk)                       :: tim1(2), tim2(2)
    real(kind=rk)                       :: eos3d_self(2),      &
                                           eos3d_children(2)
    do ia=1,areas%are_nu
       ixi    = areas%ix_are(ia, 1)
       ixf    = areas%ix_are(ia, 2)
       iox    = areas%ix_are(ia, 3)
       iyi    = areas%ix_are(ia, 4)
       iyf    = areas%ix_are(ia, 5)
       ioy    = areas%ix_are(ia, 6)
       izi    = areas%ix_are(ia, 7)
       izf    = areas%ix_are(ia, 8)
       ioz    = areas%ix_are(ia, 9)

       areas%nz = izf - izi
       areas%ny = iyf - iyi
       areas%nx = ixf - ixi
       
       areas%nz = areas%nz/ioz + 1
       areas%ny = areas%ny/ioy + 1
       areas%nx = areas%nx/iox + 1
! wompi: delete ifdef, keep elseif

       if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
          den_ave(:) = 0._rk
          vex_ave(:) = 0._rk
          vey_ave(:) = 0._rk
          vez_ave(:) = 0._rk
          tem_ave(:) = 0._rk
          ene_ave(:) = 0._rk
          xnuc_ave(:,:)  = 0._rk
          
          do i = ixi, ixf
             do j = qy_s, qy_e
                do k = qz_s, qz_e
                   den_ave(i) = den_ave(i) + dentot(i,j,k) * dvytot(j) * dvztot(k)
                   vex_ave(i) = vex_ave(i) + vextot(i,j,k) * dentot(i,j,k) * dvytot(j) * dvztot(k)
                   vey_ave(i) = vey_ave(i) + veytot(i,j,k) * dentot(i,j,k) * dvytot(j) * dvztot(k)
                   vez_ave(i) = vez_ave(i) + veztot(i,j,k) * dentot(i,j,k) * dvytot(j) * dvztot(k)
                   tem_ave(i) = tem_ave(i) + temtot(i,j,k)                 * dvytot(j) * dvztot(k)
                   ene_ave(i) = ene_ave(i) + enetot(i,j,k) * dentot(i,j,k) * dvytot(j) * dvztot(k)
                   do l=1,config%qn
                      xnuc_ave(i,l) = xnuc_ave(i,l) + xnutot(i,j,k,l)* dentot(i,j,k) * dvytot(j) * dvztot(k)
                   enddo
                   
                end do
             end do
          end do


          if (use_mpi) then
#ifndef DEBUG_TIMINGS
             call second_v(tim1)
#endif   
             ! mpi-allreduce (addition) fuer die Felder den_ave, vex_ave, ...

             call MPI_Allreduce(den_ave(ixi:ixf), den_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             den_ave(ixi:ixf)=den_ave_buf(ixi:ixf)

             call MPI_Allreduce(vex_ave(ixi:ixf), vex_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             vex_ave(ixi:ixf)=vex_ave_buf(ixi:ixf)

             call MPI_Allreduce(vey_ave(ixi:ixf), vey_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             vey_ave(ixi:ixf)=vey_ave_buf(ixi:ixf)

             call MPI_Allreduce(vez_ave(ixi:ixf), vez_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             vez_ave(ixi:ixf)=vez_ave_buf(ixi:ixf)

             call MPI_Allreduce(tem_ave(ixi:ixf), tem_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             tem_ave(ixi:ixf)=tem_ave_buf(ixi:ixf)


             call MPI_Allreduce(ene_ave(ixi:ixf), ene_ave_buf(ixi:ixf), areas%nx, &
                                MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
             ene_ave(ixi:ixf)=ene_ave_buf(ixi:ixf)

             do l=1,config%qn
                call MPI_Allreduce(xnuc_ave(ixi:ixf,l), xnuc_ave_buf(ixi:ixf,l), &
                                   areas%nx, MPI_DOUBLE_PRECISION, MPI_SUM,      &
                                   MPI_COMM_WORLD, ierr)
             end do

             xnuc_ave(ixi:ixf,1:config%qn)=xnuc_ave_buf(ixi:ixf,1:config%qn)
#ifndef DEBUG_TIMINGS
             call second_v(tim2)

             timer%transp_comm =timer%transp_comm + (tim2-tim1)
#endif
          endif ! use_mpi




          do i = ixi, ixf
             vex_ave(i) = vex_ave(i) / den_ave(i)
             vey_ave(i) = vey_ave(i) / den_ave(i)
             vez_ave(i) = vez_ave(i) / den_ave(i)
             ene_ave(i) = ene_ave(i) / den_ave(i)
             do l=1,config%qn
                xnuc_ave(i,l) = xnuc_ave(i,l) / den_ave(i)
             enddo

             den_ave(i) = den_ave(i) / (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
             tem_ave(i) = tem_ave(i) / (sum(dvytot(1:config%qy))*sum(dvztot(1:config%qz)))
 
              ! copy back to tot arrays
             do j = qy_s, qy_e
                do k = qz_s, qz_e
                   dentot(i,j,k) = den_ave(i)
                   vextot(i,j,k) = vex_ave(i)
                   veytot(i,j,k) = vey_ave(i)
                   veztot(i,j,k) = vez_ave(i)
                   enetot(i,j,k) = ene_ave(i)
                   temtot(i,j,k) = tem_ave(i)
                   do l=1,config%qn
                      xnutot(i,j,k,l) = xnuc_ave(i,l)
                   enddo
                end do
             end do
          end do

          call cpyare(2)
          call eos3d (2,ler, eos3d_self, eos3d_children)
          if (ler) then
             raise_abort("mapare(): eos failed")
          endif

          do i = ixi, ixf
             do j = qy_s, qy_e
                do k = qz_s, qz_e
                   temtot(i,j,k) = temp(i,j,k)
                enddo
             enddo
          enddo

       else if (ioy .eq. 1) then
          return
       else
          raise_abort("mapare(): mode not supported")
!            call stopit_mpi('mapare: mode not supported')
       endif

    enddo

  end subroutine mapare

end module mapare_proc

module mpi_comm_routines

  private


  public  :: mpi_comm_lb_ub_1neighbour, &
             mpi_comm_lb_ub_2neighbour, &
             mpi_comm_lb_ub_4neighbour, &
             mpi_comm_pb_kb_1neighbour, & 
             mpi_comm_pb_kb_2neighbour, & 
             mpi_comm_pb_kb_4neighbour, & 
             exch_allgatherv, exch, &
             stopit_mpi


  contains

 subroutine mpi_comm_lb_ub_1neighbour
!     This routine should only be used when one has to communicate only
!     with one left/right neighbour (i.e. per MPI-Task more than 4 zones
!     in each direction are used)
!         
!     to next neighbor in temporary array to use when setting the special 
!     boundary conditions
!     to densty_lower_boundary(qx,4,qz)
    use precision
      
    use vnew_hy 
    use vold_hy, only: vxold, vzold 
    use mo_mpi
    use configure
    use abort
    
    IMPLICIT NONE

    integer(kind=ik):: i, j, k, n, dest, src, ierr, &
                       mpistat(MPI_STATUS_SIZE), sendcount
    real(kind=rk)   :: gpocntr_ltmp(1:config%qx,qy_s:qy_s+3,qz_s:qz_e), &
                       gpocntr_utmp(1:config%qx,qy_e-3:qy_e,qz_s:qz_e)
    real(kind=rk)   :: sbuf(config%qx,qz_s:qz_e,44+4*config%qn),        &
                       rbuf(config%qx,qz_s:qz_e,44+4*config%qn), &
                       sbufv(config%qx,qz_s-1:qz_e+1,4), rbufv(config%qx,qz_s-1:qz_e+1,4)

    ! only next neighbour communication
    ! (i.e. at least four zones per MPI-task)
    if (qy_s+3 .gt. qy_e) then
       raise_abort("mpi_comm_lb_ub_1neighbour(): less then four zones per MPI-TASK")
    endif

    ! prepare next neighbour communciation 
    do k = qz_s, qz_e
       do j = qy_e-3, qy_e ! in case of four zones this is equal 
                           ! to qy_s, qy_e. BUT in case of more zones
                           ! we only need the _last_ four zones
          do i = 1, size(gpocntr_utmp, dim=1)
             gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
    
    do k = qz_s, qz_e
       do j = qy_s, qy_s+3 ! in case of four zones this is equal
                           ! to qy_s, qy_e. BUT in case of more zones
                           ! we only need
          do i = 1, size(gpocntr_ltmp,dim=1)
             gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
    
    if (nprocs .eq. 1) then
       ! this is necessary because otherwise rbuf will be undefined 
       ! in send_left / send_right 
       do i = 1, 4 
          densty_ub (1:config%qx,i,qz_s:qz_e) = densty      (1:config%qx,qy_s+i-1,qz_s:qz_e)
          vely_ub   (1:config%qx,i,qz_s:qz_e) = vely        (1:config%qx,qy_s+i-1,qz_s:qz_e)
          velx_ub   (1:config%qx,i,qz_s:qz_e) = velx        (1:config%qx,qy_s+i-1,qz_s:qz_e)
          velz_ub   (1:config%qx,i,qz_s:qz_e) = velz        (1:config%qx,qy_s+i-1,qz_s:qz_e)
          vxold_ub  (1:config%qx,i,qz_s:qz_e) = vxold       (1:config%qx,qy_s+i-1,qz_s:qz_e)
          vzold_ub  (1:config%qx,i,qz_s-1:qz_e+1) = vzold   (1:config%qx,qy_s+i-1,qz_s-1:qz_e+1)
          energy_ub (1:config%qx,i,qz_s:qz_e) = energy      (1:config%qx,qy_s+i-1,qz_s:qz_e)
          press_ub  (1:config%qx,i,qz_s:qz_e) = press       (1:config%qx,qy_s+i-1,qz_s:qz_e)
          temp_ub   (1:config%qx,i,qz_s:qz_e) = temp        (1:config%qx,qy_s+i-1,qz_s:qz_e)
          gammae_ub (1:config%qx,i,qz_s:qz_e) = gammae      (1:config%qx,qy_s+i-1,qz_s:qz_e)
          gammac_ub (1:config%qx,i,qz_s:qz_e) = gammac      (1:config%qx,qy_s+i-1,qz_s:qz_e)
          gpocntr_ub(1:config%qx,i,qz_s:qz_e) = gpocntr_ltmp(1:config%qx,qy_s+i-1,qz_s:qz_e)
          do n = 1, size(xnuc_ub,dim=4)
             xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = xnuc(1:config%qx,qy_s+i-1,qz_s:qz_e,n)
          end do
          
          densty_lb (1:config%qx,i,qz_s:qz_e) = densty      (1:config%qx,qy_e-4+i,qz_s:qz_e)
          vely_lb   (1:config%qx,i,qz_s:qz_e) = vely        (1:config%qx,qy_e-4+i,qz_s:qz_e)
          velx_lb   (1:config%qx,i,qz_s:qz_e) = velx        (1:config%qx,qy_e-4+i,qz_s:qz_e)
          velz_lb   (1:config%qx,i,qz_s:qz_e) = velz        (1:config%qx,qy_e-4+i,qz_s:qz_e)
          vxold_lb  (1:config%qx,i,qz_s:qz_e) = vxold       (1:config%qx,qy_e-4+i,qz_s:qz_e)
          vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = vzold   (1:config%qx,qy_e-4+i,qz_s-1:qz_e+1)
          energy_lb (1:config%qx,i,qz_s:qz_e) = energy      (1:config%qx,qy_e-4+i,qz_s:qz_e)
          press_lb  (1:config%qx,i,qz_s:qz_e) = press       (1:config%qx,qy_e-4+i,qz_s:qz_e)
          temp_lb   (1:config%qx,i,qz_s:qz_e) = temp        (1:config%qx,qy_e-4+i,qz_s:qz_e)
          gammae_lb (1:config%qx,i,qz_s:qz_e) = gammae      (1:config%qx,qy_e-4+i,qz_s:qz_e)
          gammac_lb (1:config%qx,i,qz_s:qz_e) = gammac      (1:config%qx,qy_e-4+i,qz_s:qz_e)
          gpocntr_lb(1:config%qx,i,qz_s:qz_e) = gpocntr_utmp(1:config%qx,qy_e-4+i,qz_s:qz_e)
          do n = 1, size(xnuc_lb,dim=4)
             xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = xnuc(1:config%qx,qy_e-4+i,qz_s:qz_e,n)
          end do
       end do
         
    else ! nprocs .gt. 1 
         ! do MPI next neighbor communication because of access in getrwy
         ! send e.g. densty(:,qy_e-3:qy_e,:) to right neighbor who stores 
         ! data in temporary array densty_lb(:,qy_s-4:qy_s-1,:), etc.
         ! same for  vely, velx, velz, vxold, vzold, energy, press, temp, 
         ! gammae, gammac, gpotcntr_up
         ! same for lower 4 sectors: xxx(:,qy_s:qy_s+3,:) is sent left 
         ! and stored in temporary array xxx_ub(:,qy_e+1:qy_e+4,:)
         ! send gpocntr_lb 
         
       ! send left
       if (myproc .gt. 0) then ! ensure left reciever will exist
          do i = 1, 4 ! send 4 values
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = vely        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = velx        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,13+i-1) = velz        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,17+i-1) = vxold       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,21+i-1) = energy      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,25+i-1) = press       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,29+i-1) = temp        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,33+i-1) = gammae      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,37+i-1) = gammac      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,41+i-1) = gpocntr_ltmp(1:config%qx,qy_s+i-1,qz_s:qz_e)
             do n = 1, size(xnuc_lb, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,41+i-1+4*n) = xnuc(1:config%qx,qy_s+i-1,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_s+i-1,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr) ! shift left
       
       sendcount =  config%qx*(44+4*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*4*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-1) then ! ensure right sender existed
          do i = 1, 4 ! recieve four values
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-1)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-1)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-1)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,25+i-1)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,29+i-1)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,33+i-1)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,37+i-1)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,41+i-1) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,41+i-1+4*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if ! nprocs .gt. 1
       
       ! send right
       if (myproc .lt. nprocs-1) then ! ensure right reciever will exist
          do i = 1, 4  ! send four values
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = vely        (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = velx        (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,13+i-1) = velz        (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,17+i-1) = vxold       (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,21+i-1) = energy      (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,25+i-1) = press       (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,29+i-1) = temp        (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,33+i-1) = gammae      (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,37+i-1) = gammac      (1:config%qx,qy_e-4+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,41+i-1) = gpocntr_utmp(1:config%qx,qy_e-4+i,qz_s:qz_e)
             do n = 1, size(xnuc_ub, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,41+i-1+4*n) = xnuc(1:config%qx,qy_e-4+i,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_e-4+i,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr) ! shift right
       sendcount = config%qx*(44+4*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*4*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 0) then ! ensure left sender existed
          do i = 1, 4  ! recieve four values
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-1)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-1)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-1)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,25+i-1)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,29+i-1)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,33+i-1)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,37+i-1)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,41+i-1) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,41+i-1+4*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if
       
    endif ! myprocs .eq. / .gt. 1
  end subroutine mpi_comm_lb_ub_1neighbour

 subroutine mpi_comm_lb_ub_2neighbour
!     This routine should only be used when one has to communicate only
!     with two left/right neighbours (i.e. per MPI-Task exactly 2 zones
!     in each direction are used)
!         
!     to next neighbor in temporary array to use when setting the special 
!     boundary conditions
!     to densty_lower_boundary(qx,4,qz)
    use precision
      
    use vnew_hy 
    use vold_hy, only: vxold, vzold 
    use mo_mpi
    use configure
    use abort

    IMPLICIT NONE

    integer(kind=ik) :: k, j, i, n, dest, src, ierr, &
                     mpistat(MPI_STATUS_SIZE), sendcount 
    real(kind=rk) ::  gpocntr_ltmp(1:config%qx,qy_s:qy_s+1,qz_s:qz_e), &
                      gpocntr_utmp(1:config%qx,qy_e-1:qy_e,qz_s:qz_e)
    real(kind=rk)  :: sbuf(config%qx,qz_s:qz_e,22+2*config%qn), rbuf(config%qx,qz_s:qz_e,22+2*config%qn), &
                      sbufv(config%qx,qz_s-1:qz_e+1,2), rbufv(config%qx,qz_s-1:qz_e+1,2)



    ! next-to-next neighbour communication
    ! (i.e. exactly 2 zones per MPI-task)
    if (qy_s +1 .ne. qy_e) then
       raise_abort("mpi_comm_lb_ub_2neighbour(): not 2 zones per MPI-TASK")
    endif

    do k = qz_s, qz_e
       do j = qy_s, qy_e     ! exactly two zones per MPI-task
          do i = 1, size(gpocntr_utmp, dim=1)
             gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
     
    do k = qz_s, qz_e
       do j = qy_s, qy_e ! exactly two zones per MPI-task
          do i = 1, size(gpocntr_ltmp,dim=1)
             gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
    
    if (nprocs .eq. 1) then
       write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
       call stopit_mpi("mpi_communication: employing only one MPI task requires at least four azimutal zones!!!")
    else ! nprocs .gt. 1 
       
       ! send left
       if (myproc .gt. 0) then ! ensure that left reciever will exist
          do i = 1, 2 ! send to values
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 3+i-1) = vely        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = velx        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 7+i-1) = velz        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = vxold       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,11+i-1) = energy      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,13+i-1) = press       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,15+i-1) = temp        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,17+i-1) = gammae      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,19+i-1) = gammac      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,21+i-1) = gpocntr_ltmp(1:config%qx,qy_s+i-1,qz_s:qz_e)
             do n = 1, size(xnuc_lb, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,21+i-1+2*n) = xnuc(1:config%qx,qy_s+i-1,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_s+i-1,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr) ! shift left to next neighbour
       sendcount =  config%qx*(22+2*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*2*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-1) then ! ensure that right sender existed
          do i = 1, 2 ! send 2 values
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-1)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-1)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-1)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-1)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,15+i-1)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-1)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,19+i-1)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-1) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,21+i-1+2*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if
       
       call MPI_Cart_shift(cart_comm,0,-2,src,dest,ierr) ! shift left to next-to-next neighbour
       sendcount =  config%qx*(22+2*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*2*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-2) then ! ensure that right next-to-next sender existed
          do i = 3, 4 ! recieve 2 values
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-3)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-3)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-3)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-3)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-3)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-3)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-3)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,15+i-3)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-3)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,19+i-3)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-3) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,21+i-3+2*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-2)
          end do
       end if
       
       ! send right
       if (myproc .lt. nprocs-1) then ! ensure that right reciever exists
          do i = 1, 2 ! send two values
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 3+i-1) = vely        (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = velx        (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 7+i-1) = velz        (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = vxold       (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,11+i-1) = energy      (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,13+i-1) = press       (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,15+i-1) = temp        (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,17+i-1) = gammae      (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,19+i-1) = gammac      (1:config%qx,qy_e-2+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,21+i-1) = gpocntr_utmp(1:config%qx,qy_e-2+i,qz_s:qz_e)
             do n = 1, size(xnuc_ub, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,21+i-1+2*n) = xnuc(1:config%qx,qy_e-2+i,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_e-2+i,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr) ! shift to next right neighbour
       sendcount = config%qx*(22+2*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*2*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 0) then ! ensure that left sender existed
          do i = 3, 4 ! recieve two values
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-3)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-3)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-3)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-3)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-3)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-3)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-3)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,15+i-3)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-3)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,19+i-3)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-3) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,21+i-3+2*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-2)
          end do
       end if
       
       call MPI_Cart_shift(cart_comm,0,2,src,dest,ierr) ! shift to next-to next right neighbour
       sendcount = config%qx*(22+2*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*2*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 1) then  ! ensure that next-to-next sender existed
          do i = 1, 2 ! recieve two values
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-1)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-1)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-1)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,13+i-1)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,15+i-1)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,17+i-1)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,19+i-1)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,21+i-1) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,21+i-1+2*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if
       
    endif ! myprocs .eq. / .gt. 1
     
  
end subroutine mpi_comm_lb_ub_2neighbour

 subroutine mpi_comm_lb_ub_4neighbour
!     This routine should only be used when one has to communicate 
!     with 4 left/right neighbours (i.e. per MPI-Task exactly 1 zone
!     in each direction are used)
!         
!     to next neighbor in temporary array to use when setting the special 
!     boundary conditions
!     to densty_lower_boundary(qx,4,qz)
    use precision
      
    use vnew_hy 
    use vold_hy, only: vxold, vzold 
    use mo_mpi
    use configure
    use abort

    IMPLICIT NONE

    integer(kind=ik) :: i,j,k, n, dest, src, ierr, &
                     mpistat(MPI_STATUS_SIZE), sendcount
    real(kind=rk) ::  gpocntr_ltmp(1:config%qx,qy_s:qy_s+1,qz_s:qz_e), &
                      gpocntr_utmp(1:config%qx,qy_e-1:qy_e,qz_s:qz_e)
    real(kind=rk)  :: sbuf(config%qx,qz_s:qz_e,11+1*config%qn), rbuf(config%qx,qz_s:qz_e,11+1*config%qn), &
                      sbufv(config%qx,qz_s-1:qz_e+1,1), rbufv(config%qx,qz_s-1:qz_e+1,1)



    ! next-to-next neighbour communication
    ! (i.e. exactly 1 zones per MPI-task)
    if (qy_s .ne. qy_e) then
       raise_abort("mpi_comm_lb_ub_4neighbour(): not one zone per MPI-task")
    endif

    do k = qz_s, qz_e
       do j = qy_s, qy_e     ! exactly one zones per MPI-task
          do i = 1, size(gpocntr_utmp, dim=1)
             gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
     
    do k = qz_s, qz_e
       do j = qy_s, qy_e ! exactly one zones per MPI-task
          do i = 1, size(gpocntr_ltmp,dim=1)
             gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
          end do
       end do
    end do
    
    if (nprocs .eq. 1) then
       write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
       call stopit_mpi("mpi_communication: employing only one MPI task requires at least four azimutal zones!!!")
    else ! nprocs .gt. 1 
       
       ! send left
       if (myproc .gt. 0) then ! ensure that left reciever exists
          do i = 1, 1 ! send one value
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 2+i-1) = vely        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 3+i-1) = velx        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 4+i-1) = velz        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = vxold       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 6+i-1) = energy      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 7+i-1) = press       (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 8+i-1) = temp        (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = gammae      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,10+i-1) = gammac      (1:config%qx,qy_s+i-1,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,11+i-1) = gpocntr_ltmp(1:config%qx,qy_s+i-1,qz_s:qz_e)
             do n = 1, size(xnuc_lb, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,11+i-1+1*n) = xnuc(1:config%qx,qy_s+i-1,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_s+i-1,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr) ! shift one left
       sendcount =  config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-1) then ! ensure that one right sender exists
          do i = 1, 1 ! recieved 1 value
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-1)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-1)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-1)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-1)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-1)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-1)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-1)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-1) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-1+1*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if
       
       call MPI_Cart_shift(cart_comm,0,-2,src,dest,ierr) ! shift to second left neighbour
       sendcount =  config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-2) then ! enssure that next-to-next right sender exist 
          do i = 2,2 ! send one value
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-2)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-2)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-2)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-2)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-2)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-2)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-2)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-2)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-2)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-2)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-2) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-2+1*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-1)
          end do
       end if

      
       call MPI_Cart_shift(cart_comm,0,-3,src,dest,ierr) ! shift to third left neighbour
       sendcount =  config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-3) then ! enssure that third right sender exist 
          do i = 3,3 ! send one value
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-3)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-3)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-3)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-3)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-3)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-3)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-3)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-3)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-3)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-3)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-3) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-3+1*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-2)
          end do
       end if


     
       call MPI_Cart_shift(cart_comm,0,-4,src,dest,ierr) ! shift to fourth left neighbour
       sendcount =  config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .lt. nprocs-4) then ! enssure that fourth right sender exist 
          do i = 4,4 ! send one value
             densty_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-4)
             vely_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-4)
             velx_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-4)
             velz_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-4)
             vxold_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-4)
             energy_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-4)
             press_ub  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-4)
             temp_ub   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-4)
             gammae_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-4)
             gammac_ub (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-4)
             gpocntr_ub(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-4) 
             do n = 1, size(xnuc_ub, dim=4)
                xnuc_ub(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-4+1*n) 
             end do
             
             vzold_ub   (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-3)
          end do
       end if

       
       ! send right
       if (myproc .lt. nprocs-1) then ! ensure that one right neighbour exists
          do i = 1, 1 
             sbuf(1:config%qx,qz_s:qz_e, 1+i-1) = densty      (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 2+i-1) = vely        (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 3+i-1) = velx        (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 4+i-1) = velz        (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 5+i-1) = vxold       (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 6+i-1) = energy      (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 7+i-1) = press       (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 8+i-1) = temp        (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e, 9+i-1) = gammae      (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,10+i-1) = gammac      (1:config%qx,qy_e-1+i,qz_s:qz_e)
             sbuf(1:config%qx,qz_s:qz_e,11+i-1) = gpocntr_utmp(1:config%qx,qy_e-1+i,qz_s:qz_e)
             do n = 1, size(xnuc_ub, dim=4)
                sbuf(1:config%qx,qz_s:qz_e,11+i-1+1*n) = xnuc(1:config%qx,qy_e-1+i,qz_s:qz_e,n)
             end do
             
             sbufv(1:config%qx,qz_s-1:qz_e+1,i) = vzold       (1:config%qx,qy_e-1+i,qz_s-1:qz_e+1)
          end do
       endif
       
       call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
       sendcount = config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 0) then ! ensure that one left sender exists
          do i = 4, 4 
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-4)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-4)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-4)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-4)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-4)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-4)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-4)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-4)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-4)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-4)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-4) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-4+1*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-3)
          end do
       end if
       
       call MPI_Cart_shift(cart_comm,0,2,src,dest,ierr)
       sendcount = config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 1) then ! ensure that two left senders exists
          do i = 3, 3 
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-3)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-3)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-3)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-3)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-3)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-3)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-3)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-3)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-3)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-3)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-3) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-3+1*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-2)
          end do
       end if
        
       call MPI_Cart_shift(cart_comm,0,3,src,dest,ierr)
       sendcount = config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       if (myproc .gt. 2) then ! ensure that three left senders exists
          do i = 2, 2 
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-2)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-2)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-2)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-2)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-2)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-2)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-2)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-2)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-2)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-2)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-2) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-2+1*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i-1)
          end do
       end if     
       
       call MPI_Cart_shift(cart_comm,0,4,src,dest,ierr)
       sendcount = config%qx*(11+1*config%qn)*qz_proc
       call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
            rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
            MPI_COMM_WORLD, mpistat, ierr)
       sendcount = config%qx*1*(qz_proc+2)
       call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
            rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
            MPI_COMM_WORLD, mpistat, ierr)
       
       
       if (myproc .gt. 3) then ! ensure that four left senders exists
          do i = 1, 1 
             densty_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 1+i-1)
             vely_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 2+i-1)
             velx_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 3+i-1)
             velz_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 4+i-1)
             vxold_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 5+i-1)
             energy_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 6+i-1)
             press_lb  (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 7+i-1)
             temp_lb   (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 8+i-1)
             gammae_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e, 9+i-1)
             gammac_lb (1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,10+i-1)
             gpocntr_lb(1:config%qx,i,qz_s:qz_e) = rbuf(1:config%qx,qz_s:qz_e,11+i-1) 
             do n = 1, size(xnuc_lb, dim=4)
                xnuc_lb(1:config%qx,i,qz_s:qz_e,n) = rbuf(1:config%qx,qz_s:qz_e,11+i-1+1*n) 
             end do
             
             vzold_lb  (1:config%qx,i,qz_s-1:qz_e+1) = rbufv(1:config%qx,qz_s-1:qz_e+1,i)
          end do
       end if  
    endif ! myprocs .eq. / .gt. 1
     
  
end subroutine mpi_comm_lb_ub_4neighbour


subroutine mpi_comm_pb_kb_1neighbour
! this routine should only be called if one has to communicate with
! exactly one left/right neighbour (i.e. at least 4 zones per MPI-task)
!         
!     to next neighbor in temporary array to use when setting the special boundary 
!     conditions
!     to densty_lower_boundary(qx,4,qy)
  use precision
  
  use vnew_hy 
  use vold_hy, only: vxold, vyold, vzold 
  use mo_mpi
  
  use configure
  use abort
  IMPLICIT NONE
  
  integer(kind=ik):: i, j, k, n, dest, src, ierr, &
                     mpistat(MPI_STATUS_SIZE), sendcount
  real(kind=rk) ::  gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s:qz_s+3), &
                    gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-3:qz_e)
  real(kind=rk)  ::  sbuf(config%qx,qy_s:qy_e,44+4*config%qn), rbuf(config%qx,qy_s:qy_e,44+4*config%qn), &
                     sbufv(config%qx,qy_s-1:qy_e+1,4), rbufv(config%qx,qy_s-1:qy_e+1,4)


  if (qz_s+3 .gt. qz_e) then
     raise_abort("mpi_comm_bk_kb_1neighbour(): less then four zones per MPI-TASK")
  endif

  ! prepare next neighbour communciation 
  do k = qz_e-3, qz_e
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_utmp, dim=1)
           gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  do k = qz_s, qz_s+3
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_ltmp, dim=1)
           gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  
  if (nprocs .eq. 1) then
     ! this is necessary because otherwise rbuf will be undefined 
     ! in send_left / send_right 
     do i = 1, 4 
        densty_kb (1:config%qx,qy_s:qy_e,i) = densty      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = vely        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = velx        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = velz        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = vxold       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = vyold   (1:config%qx,qy_s-1:qy_e+1,qz_s+i-1)
        energy_kb (1:config%qx,qy_s:qy_e,i) = energy      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        press_kb  (1:config%qx,qy_s:qy_e,i) = press       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = temp        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = gammae      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = gammac      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s+i-1)
        do n = 1, size(xnuc_ub, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = xnuc(1:config%qx,qy_s:qy_e,qz_s+i-1,n)
        end do
        
        densty_pb (1:config%qx,qy_s:qy_e,i) = densty      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = vely        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = velx        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = velz        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = vxold       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_e-4+i)
        energy_pb (1:config%qx,qy_s:qy_e,i) = energy      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        press_pb  (1:config%qx,qy_s:qy_e,i) = press       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = temp        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = gammae      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = gammac      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-4+i)
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = xnuc(1:config%qx,qy_s:qy_e,qz_e-4+i,n)
        end do
     end do
     
  else ! nprocs .gt. 1 
     ! do MPI next neighbor communication because of access in getrwy
     ! send e.g. densty(:,qy_e-3:qy_e,:) to right neighbor who stores 
     ! data in temporary array densty_lb(:,qy_s-4:qy_s-1,:), etc.
     ! same for  vely, velx, velz, vxold, vzold, energy, press, temp, 
     ! gammae, gammac, gpotcntr_up
     ! same for lower 4 sectors: xxx(:,qy_s:qy_s+3,:) is sent left 
     ! and stored in temporary array xxx_ub(:,qy_e+qy_s:qy_e_e+4,:)
     ! send gpocntr_lb 
     
     ! send left
     do i = 1, 4 
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,13+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,17+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,21+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,25+i-1) = press       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,29+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,33+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,37+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,41+i-1) = gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s+i-1)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,41+i-1+4*n) = xnuc(1:config%qx,qy_s:qy_e,qz_s+i-1,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_s+i-1)
     end do
     
     call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(44+4*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*4*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 4
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-1)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-1)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-1)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,25+i-1)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,29+i-1)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,33+i-1)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,37+i-1)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,41+i-1) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,41+i-1+4*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do
     
     ! send right
     do i = 1, 4 
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,13+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,17+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,21+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,25+i-1) = press       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,29+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,33+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,37+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,41+i-1) = gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-4+i)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,41+i-1+4*n) = xnuc(1:config%qx,qy_s:qy_e,qz_e-4+i,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_e-4+i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(44+4*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*4*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 4 
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-1)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-1)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-1)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,25+i-1)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,29+i-1)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,33+i-1)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,37+i-1)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,41+i-1) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,41+i-1+4*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do
     
  endif ! myprocs .eq. / .gt. 1
  
end subroutine mpi_comm_pb_kb_1neighbour


subroutine mpi_comm_pb_kb_2neighbour
! this routine should only be called if one has to communicate with
! exactly two left/right neighbours (i.e. 2 Zones per MPI-task)
         
!     to next neighbor in temporary array to use when setting the special boundary 
!     conditions
!     to densty_lower_boundary(qx,4,qy)
  use precision
  
  use vnew_hy 
  use vold_hy, only: vxold, vyold, vzold 
  use mo_mpi
  
  use configure
  use abort
  IMPLICIT NONE
  
  integer(kind=ik):: i, j, k, n, dest, src, ierr, &
                     mpistat(MPI_STATUS_SIZE), sendcount
 
  real(kind=rk) ::  gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s:qz_s+1), &
                    gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-1:qz_e)
  real(kind=rk)  ::  sbuf(config%qx,qy_s:qy_e,22+2*config%qn), rbuf(config%qx,qy_s:qy_e,22+2*config%qn), &
                     sbufv(config%qx,qy_s-1:qy_e+1,2), rbufv(config%qx,qy_s-1:qy_e+1,2)



  if (qz_s+1 .ne. qz_e) then
     raise_abort("mpi_comm_bk_kb_2neighbour(): not two zones per MPI-TASK")
  endif

    ! only two zones are in one MPI task
    ! therefore communicate with the next two tasks
    ! prepare next neighbor communciation 
  do k = qz_e-1, qz_e
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_utmp, dim=1)
           gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  do k = qz_s, qz_s+1
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_ltmp, dim=1)
           gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  
  if (nprocs .eq. 1) then
     write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
     call stopit_mpi("mpi_communication: employing only one MPI task requires at least four azimutal zones!!!")
  else ! nprocs .gt. 1 
     
     ! send left
     do i = 1, 2 
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 3+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,7+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,9+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,11+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,13+i-1) = press       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,15+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,17+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,19+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,21+i-1) = gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s+i-1)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,21+i-1+2*n) = xnuc(1:config%qx,qy_s:qy_e,qz_s+i-1,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_s+i-1)
     end do
     
     call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(22+2*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*2*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 2
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-1)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-1)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-1)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-1)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,15+i-1)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-1)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,19+i-1)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-1) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,21+i-1+2*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,-2,src,dest,ierr)
     sendcount = config%qx*qy_proc*(22+2*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*2*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 3, 4
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-3)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-3)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-3)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-3)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-3)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-3)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-3)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,15+i-3)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-3)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,19+i-3)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-3) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,21+i-3+2*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-2)
     end do
     
     ! send right
     do i = 1, 2
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e, 3+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e, 7+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,11+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,13+i-1) = press       (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,15+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,17+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,19+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_e-2+i)
        sbuf(1:config%qx,qy_s:qy_e,21+i-1) = gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-2+i)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,21+i-1+2*n) = xnuc(1:config%qx,qy_s:qy_e,qz_e-2+i,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_e-2+i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(22+2*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*2*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 3, 4
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-3)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-3)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-3)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-3)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-3)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-3)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-3)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,15+i-3)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-3)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,19+i-3)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-3) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,21+i-3+2*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-2)
     end do
     
     call MPI_Cart_shift(cart_comm,1,2,src,dest,ierr)
     sendcount = config%qx*qy_proc*(22+2*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*2*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 2
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-1)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-1)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-1)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,13+i-1)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,15+i-1)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,17+i-1)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,19+i-1)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,21+i-1) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,21+i-1+2*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do
     
  endif ! myprocs .eq. / .gt. 1
  
end subroutine mpi_comm_pb_kb_2neighbour

subroutine mpi_comm_pb_kb_4neighbour
! this routine should only be called if one has to communicate with
! exactly four left/right neighbours (i.e. 1 Zones per MPI-task)
         
!     to next neighbor in temporary array to use when setting the special boundary 
!     conditions
!     to densty_lower_boundary(qx,4,qy)
  use precision
  
  use vnew_hy 
  use vold_hy, only: vxold, vyold, vzold 
  use mo_mpi
  
  use configure
  use abort
  IMPLICIT NONE
  
  integer(kind=ik):: i, j, k, n, dest, src, ierr, &
                     mpistat(MPI_STATUS_SIZE), sendcount
 
  real(kind=rk) ::  gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s:qz_s+1), &
                    gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-1:qz_e)
  real(kind=rk)  ::  sbuf(config%qx,qy_s:qy_e,11+1*config%qn), rbuf(config%qx,qy_s:qy_e,11+1*config%qn), &
                     sbufv(config%qx,qy_s-1:qy_e+1,1), rbufv(config%qx,qy_s-1:qy_e+1,1)



  if (qz_s .ne. qz_e) then
     raise_abort("mpi_comm_bk_kb_4neighbour(): not one zone per MPI-TASK")
  endif

    ! only two zones are in one MPI task
    ! therefore communicate with the next two tasks
    ! prepare next neighbor communciation 
  do k = qz_s, qz_e
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_utmp, dim=1)
           gpocntr_utmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  do k = qz_s, qz_e
     do j = qy_s, qy_e
        do i = 1, size(gpocntr_ltmp, dim=1)
           gpocntr_ltmp(i,j,k) = 0.5_rk * ( gpot(i,j,k) + gpot(i-1,j,k) )
        end do
     end do
  end do
  
  if (nprocs .eq. 1) then
     write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
     call stopit_mpi("mpi_communication: employing only one MPI task requires at least four azimutal zones!!!")
  else ! nprocs .gt. 1 
     
     ! send left
     do i = 1, 1
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 2+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 3+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 4+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 6+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 7+i-1) = press       (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 8+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,10+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_s+i-1)
        sbuf(1:config%qx,qy_s:qy_e,11+i-1) = gpocntr_ltmp(1:config%qx,qy_s:qy_e,qz_s+i-1)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,11+i-1+1*n) = xnuc(1:config%qx,qy_s:qy_e,qz_s+i-1,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_s+i-1)
     end do
     
     call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 1
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-1)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-1)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-1)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-1)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-1)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-1)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-1)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-1) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-1+1*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,-2,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)

     do i = 2, 2
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-2)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-2)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-2)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-2)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-2)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-2)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-2)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-2)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-2)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-2)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-2) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-2+1*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-1)
     end do

     
     call MPI_Cart_shift(cart_comm,1,-3,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)

     
     do i = 3, 3
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-3)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-3)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-3)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-3)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-3)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-3)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-3)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-3)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-3)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-3)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-3) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-3+1*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-2)
     end do
     

     
     call MPI_Cart_shift(cart_comm,1,-4,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep3, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep3,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep31, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep31,  &
          MPI_COMM_WORLD, mpistat, ierr)

     do i = 4, 4
        densty_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-4)
        vely_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-4)
        velx_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-4)
        velz_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-4)
        vxold_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-4)
        energy_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-4)
        press_kb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-4)
        temp_kb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-4)
        gammae_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-4)
        gammac_kb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-4)
        gpocntr_kb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-4) 
        do n = 1, size(xnuc_kb, dim=4)
           xnuc_kb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-4+1*n) 
        end do
        
        vyold_kb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-3)
     end do


     ! send right
     do i = 1, 1
        sbuf(1:config%qx,qy_s:qy_e, 1+i-1) = densty      (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 2+i-1) = vely        (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 3+i-1) = velx        (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 4+i-1) = velz        (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-1) = vxold       (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 6+i-1) = energy      (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 7+i-1) = press       (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 8+i-1) = temp        (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-1) = gammae      (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e,10+i-1) = gammac      (1:config%qx,qy_s:qy_e,qz_e-1+i)
        sbuf(1:config%qx,qy_s:qy_e,11+i-1) = gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-1+i)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,11+i-1+1*n) = xnuc(1:config%qx,qy_s:qy_e,qz_e-1+i,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_e-1+i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 4, 4
        sbuf(1:config%qx,qy_s:qy_e, 1+i-4) = densty      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 2+i-4) = vely        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 3+i-4) = velx        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 4+i-4) = velz        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 5+i-4) = vxold       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 6+i-4) = energy      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 7+i-4) = press       (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 8+i-4) = temp        (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e, 9+i-4) = gammae      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,10+i-4) = gammac      (1:config%qx,qy_s:qy_e,qz_e-4+i)
        sbuf(1:config%qx,qy_s:qy_e,11+i-4) = gpocntr_utmp(1:config%qx,qy_s:qy_e,qz_e-4+i)
        do n = 1, size(xnuc, dim=4)
           sbuf(1:config%qx,qy_s:qy_e,11+i-4+1*n) = xnuc(1:config%qx,qy_s:qy_e,qz_e-4+i,n)
        end do
        
        sbufv(1:config%qx,qy_s-1:qy_e+1,i) = vyold       (1:config%qx,qy_s-1:qy_e+1,qz_e-4+i)
     end do
     
     call MPI_Cart_shift(cart_comm,1,2,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)

     do i = 3, 3
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-3)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-3)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-3)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-3)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-3)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-3)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-3)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-3)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-3)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-3)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-3) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-3+1*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-2)
     end do
     
     call MPI_Cart_shift(cart_comm,1,3,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 2, 2
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-2)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-2)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-2)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-2)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-2)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-2)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-2)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-2)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-2)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-2)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-2) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-2+1*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i-1)
     end do
      
     call MPI_Cart_shift(cart_comm,1,4,src,dest,ierr)
     sendcount = config%qx*qy_proc*(11+1*config%qn)
     call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep4, &
          rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep4,  &
          MPI_COMM_WORLD, mpistat, ierr)
     sendcount = config%qx*1*(qy_proc+2)
     call MPI_SendRecv(sbufv, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep41, &
          rbufv, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep41,  &
          MPI_COMM_WORLD, mpistat, ierr)
     
     do i = 1, 1
        densty_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 1+i-1)
        vely_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 2+i-1)
        velx_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 3+i-1)
        velz_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 4+i-1)
        vxold_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 5+i-1)
        energy_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 6+i-1)
        press_pb  (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 7+i-1)
        temp_pb   (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 8+i-1)
        gammae_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e, 9+i-1)
        gammac_pb (1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,10+i-1)
        gpocntr_pb(1:config%qx,qy_s:qy_e,i) = rbuf(1:config%qx,qy_s:qy_e,11+i-1) 
        do n = 1, size(xnuc_pb, dim=4)
           xnuc_pb(1:config%qx,qy_s:qy_e,i,n) = rbuf(1:config%qx,qy_s:qy_e,11+i-1+1*n)
        end do
        
        vyold_pb  (1:config%qx,qy_s-1:qy_e+1,i) = rbufv(1:config%qx,qy_s-1:qy_e+1,i)
     end do    


  endif ! myprocs .eq. / .gt. 1
  
end subroutine mpi_comm_pb_kb_4neighbour      
!>
!> Subroutine to exchange data between the MPI-processes for 
!> the MPI-version of  the VERTEX-code if an error occoured
!>
!>  Author: K. Benkert and B. Mueller
!> \endverbatim
!>
!> \param data the data to exchange
!> \param num
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision: 1241 $
!>   $Date: 2010-06-28 15:47:07 +0200 (Mon, 28 Jun 2010) $
!>   
!> \endverbatim
!>
SUBROUTINE exch_allgatherv(data,num)
!-----------------------------------------
! Autor             : Katharina Benkert 


  use precision
  use mo_mpi

  use configure
  implicit none

! LOCAL variables that are not in modules
  integer(kind=ik) :: num
  real(kind=rk)    :: data(num,config%nymom, config%nztra)
  real(kind=rk)    :: data_rcv(num,config%nymom, config%nztra)
!  real(kind=rk)    :: data_rcv(num,config%nymom, config%nztra) 
  integer(kind=ik) :: n, ier 
  integer(kind=ik), dimension(0:nprocs-1) :: displs, cnts



  do n = 0, nprocs-1
     displs(n) = (nymomsg(n)-1)*num+1
     cnts(n)   = (nymomeg(n)-nymomsg(n)+1)*num
  end do

  call MPI_Allgatherv(  & !data(1,nymoms), 
                       data, cnts(myproc), &
                       MPI_DOUBLE_PRECISION, data_rcv, cnts, &
                       displs, MPI_DOUBLE_PRECISION,     &
                       MPI_COMM_WORLD, ier)

  data = data_rcv

  !      call MPI_barrier(MPI_COMM_WORLD, ier)


END subroutine exch_allgatherv


!-----------------------------------------
!>
!> Subroutine to exchange data nbetween the MPI-processes for 
!> the MPI-version of  the VERTEX-code if an error occoured
!>
!>  Author: K. Benkert and B. Mueller
!> \endverbatim
!>
!> \param data the data to exchange
!> \param num
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision: 1241 $
!>   $Date: 2010-06-28 15:47:07 +0200 (Mon, 28 Jun 2010) $
!>   
!> \endverbatim
!>
subroutine exch(data,num)

  use precision
  use mo_mpi

  use configure
  IMPLICIT NONE

! LOCAL variables that are not in modules
  integer(kind=ik) :: num
  real(kind=rk) :: data(num,config%nymom)
  integer(kind=ik) :: j, n, nr, ier
  integer(kind=ik) :: ireq(config%nymom), istat(MPI_STATUS_SIZE,config%nymom)

  ! Set up Recv for all foreign strips
  
  nr = 0
  do n=0,nprocs-1
     if(n==myproc) CYCLE
     do j=nymomsg(n),nymomeg(n)
        nr = nr+1
        CALL MPI_Irecv(data(1,j),num,MPI_REAL8,n,j,MPI_COMM_WORLD, &
                       ireq(nr),ier)
     ENDDO
  ENDDO

  ! Send our own data

  do n=0,nprocs-1
     if(n==myproc) CYCLE
     do j=nymoms, nymome
        CALL MPI_Send(data(1,j),num,MPI_REAL8,n,j,MPI_COMM_WORLD, &
                      ier)
     end do
  end do

  ! Wait
  
  CALL MPI_Waitall(nr,ireq,istat,ier)

END subroutine exch

!-----------------------------------------


!>
!> Subroutine to stop the MPI-processes for the MPI-version
!> of  the VERTEX-code if an error occoured
!>
!>  Author: K. Benkert and B. Mueller
!> \endverbatim
!>
!> \param descr description of the error
!>  
!> \verbatim
!>   SVN - Information  
!>   $Revision: 1248 $
!>   $Date: 2010-06-30 10:15:26 +0200 (Wed, 30 Jun 2010) $
!>   
!> \endverbatim
!>
SUBROUTINE stopit_mpi(descr)
      
  use precision
  use abort
  use mo_mpi
  use print_stdout_mod
  IMPLICIT NONE

! LOCAL variables that are not in modules
  
  character(*), intent(in) :: descr
  integer(kind=ik) :: ierr, error_code


  call printit_taskX(0," ",descr)

  error_code = 99_ik

  call MPI_ABORT(MPI_COMM_WORLD, error_code, ierr)

  raise_abort("MPI_ABORT() was called")

end subroutine stopit_mpi

end module mpi_comm_routines

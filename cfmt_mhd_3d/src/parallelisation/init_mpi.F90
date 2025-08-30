module init_mpi_mod

implicit none

contains

!-----------------------------------------
!>
!> Subroutine to initiallize the MPI-processes for the MPI-version
!> of  the VERTEX-code
!>
!>  Author: K. Benkert and B. Mueller
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine init_mpi

  use precision
  use mo_mpi
  use abort

  use mpi_domains,    only : optimal_decomposition
  use  mpi_comm_routines, only : stopit_mpi
  use setup_c
  use memory_bookkeeping_mod
  use parameters

  use configure

  IMPLICIT NONE

! LOCAL variables that are not in modules

  integer(kind=ik) :: n, ier, cart_dims(0:1), cart_coords(0:1)
  logical          :: cart_periods(0:1)
  logical          :: cart_reorder

  integer(kind=ik) :: n1, n2

  integer(kind=ik) :: istat, tag, i_task
  integer :: status(MPI_STATUS_SIZE)

  character(LEN=MPI_MAX_PROCESSOR_NAME) :: procname_buf
  integer(kind=ik) :: procname_length_buf

  integer(kind=ik) :: my_pid, my_ppid, pid_buf, ppid_buf

  ! variables for host identifcation
  character(LEN=MPI_MAX_PROCESSOR_NAME) :: my_procname
  integer(kind=ik) :: my_procname_length
  real(kind=rk)    :: mem_local

  mem_local = 0._rk

  call mpi_init(ier)

  if (ier .eq. -1) then
     use_mpi = .false.
  else if (ier .ne. MPI_SUCCESS) then
     raise_abort("Error in MPI initialization")
  else
     use_mpi = .true. 
  endif
  use_1neighbour_comm = .false.
  use_2neighbour_comm = .false.
  use_4neighbour_comm = .false.

  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ier)
  call mpi_comm_rank(MPI_COMM_WORLD,myproc,ier)


  ! alloctate the arrays for domain decomposition

  allocate(nymomsg(0:nprocs-1), nymomeg(0:nprocs-1), qy_sg(0:nprocs-1), &
           qy_eg(0:nprocs-1),  nzmomsg(0:nprocs-1), nzmomeg(0:nprocs-1),&
           qz_sg(0:nprocs-1), qz_eg(0:nprocs-1), stat=istat)


  if (istat .ne. 0) then
     raise_abort("Error when allocating arrays for domain decomposition")
  endif


  mem_local = nprocs*8*4._rk

  allocate(hostnames(0:nprocs-1),  procname_lengthg(0:nprocs-1), stat=istat)

  if (istat .ne. 0) then
     raise_abort("Error when allocating array for hostnames")
  endif


  mem_local = mem_local + nprocs*2*4._rk

  allocate(pidg(0:nprocs-1), ppidg(0:nprocs-1), stat=istat)

  if (istat .ne. 0) then
     raise_abort("Error when allocating array for pid's")
  endif

  ! before setting the domain decomposition we have to read
  ! the model size, which is contained in a part of the file ppm.par
  call read_parameter_files("grid_init")

  call initialize_config_values

  mem_local = mem_local + nprocs*2*4._rk

  mem_local = mem_local/1024._rk/1024._rk
  
  global_memory = global_memory + mem_local

  if (myproc .eq. 0) then
     write (*,*) "Task ",myproc," "
     write (*,*) "Task ",myproc," In init_mpi"
     write (*,'(" Task ",i4," Locally  allocated ",1pe12.5," [Mb]")')  myproc, mem_local
     write (*,'(" Task ",i4," Globally allocated ",1pe12.5," [Mb]")') myproc, global_memory
     write (*,*) "Task ",myproc," "
  endif

  ! determine the type of communication

  if (config%qy*config%qz/nprocs .ge. 4) then
     use_1neighbour_comm=.true.
     if (myproc .eq. 0) then
        print *," "
        print *," Using one left/right neighbour for communication"
        print *," "
     endif
  endif
  ! information about MPI or not

  if (use_mpi) then
     if (myproc .eq. 0) then
        print *," "
        print *," Using the MPI-version with ",nprocs," tasks"
        print *," "
     endif
  else
     if (myproc .eq. 0) then
        print *," "
        print *," Using the OpenMP-version"
        print *," "
     endif

  endif

  if (config%qy*config%qz/nprocs .eq. 2) then
     use_2neighbour_comm=.true.
     if (myproc .eq. 0) then
        print *," "
        print *," Using two left/right neighbours for communication"
        print *," "
     endif
  endif
  if (config%qy*config%qz/nprocs .eq. 1) then
     use_4neighbour_comm=.true.
     if (myproc .eq. 0) then
        print *," "
        print *," Using four left/right neighbours for communication"
        print *," "
     endif
  endif



  !     Set global start and end for every processor
  
  if (config%nsdim .eq. 3) then

     if (use_mpi) then

        !  we need to test, whether nprocs is a square number,
        ! otherwise the domain decomposition will not work
        call optimal_decomposition(config%qy,config%qz,nprocs,n1,n2)
        
        cart_dims(0)=n1
        cart_dims(1)=n2
        
        cart_periods(0)=.false. !at present, the boundary conditions need not be
        cart_periods(1)=.true.  !hardcoded, because ppm.par is read only after
                                !init_mpi is called. If possible, this should be 
                                !changed in the future
     
        cart_reorder=.false.
        
        call MPI_Cart_create(MPI_COMM_WORLD,2,cart_dims,cart_periods,cart_reorder, &
             cart_comm,ier)
        call MPI_Comm_group(cart_comm,cart_group,ier)
        call mpi_comm_rank(MPI_COMM_WORLD,myproc,ier)!rank of process might have changed

     endif ! use_mpi

     if (config%use_spherical_eddington_factor) then

        if (use_mpi) then
        
           do n=0,nprocs-1
              nymomsg(n) = 1 !(n*config%nymom)/nprocs + 1
              nymomeg(n) = 1 !((n+1)*config%nymom)/nprocs
              nzmomsg(n) = 1
              nzmomeg(n) = 1
           enddo

           do n=0,nprocs-1
              call MPI_Cart_coords(cart_comm,n,2,cart_coords,ier)
              qy_sg(n) = (cart_coords(0)*config%qy)/n1 + 1
              qy_eg(n) = ((cart_coords(0)+1)*config%qy)/n1
              qz_sg(n) = (cart_coords(1)*config%qz)/n2 + 1
              qz_eg(n) = ((cart_coords(1)+1)*config%qz)/n2
           enddo

           !     Set local boundaries
           
           nymoms = nymomsg(myproc)
           nymome = nymomeg(myproc)
           nymom_proc = 1 ! nymome-nymoms+1
           
           qy_s = qy_sg(myproc)
           qy_e = qy_eg(myproc)
           qy_proc = qy_e-qy_s+1
           
           nzmoms = nzmomsg(myproc)
           nzmome = nzmomeg(myproc)
           nzmom_proc = 1
     
           qz_s = qz_sg(myproc)
           qz_e = qz_eg(myproc)
           qz_proc = qz_e-qz_s+1
           
        else ! mpi not used


           do n=0,nprocs-1
              nymomsg(n) = 1 
              nymomeg(n) = 1 !  NTHETA_TR
              nzmomsg(n) = 1
              nzmomeg(n) = 1
           enddo
           
           do n=0,nprocs-1
              qy_sg(n) = 1
              qy_eg(n) = config%qy
              qz_sg(n) = 1
              qz_eg(n) = config%qz
           enddo

           nymoms = 1
           nymome = 1
           nymom_proc = 1
           
           qy_s = 1
           qy_e = config%qy
           qy_proc = qy_e - qy_s + 1
           
           nzmoms = 1
           nzmome = 1
           nzmom_proc = 1
           

           raise_abort("3d 1d edd faktor no mpi do not use")

           !     nymoms = 1
           !     nymome = NTHETA_TR
           !     nymom_proc = nymome-nymoms+1


           !     qz_s = 1
           !     qz_e = NQZ
           !     qz_proc = qz_e-qz_s+1

        endif ! use_mpi

     else ! spherical_eddington_factor not used


        if (use_mpi) then
        
           do n=0,nprocs-1
              call MPI_Cart_coords(cart_comm,n,2,cart_coords,ier)
              nymomsg(n) = (cart_coords(0)*config%nymom)/n1 + 1
              nymomeg(n) = ((cart_coords(0)+1)*config%nymom)/n1
              nzmomsg(n) = (cart_coords(1)*config%nztra)/n2 + 1
              nzmomeg(n) = ((cart_coords(1)+1)*config%nztra)/n2
           enddo
           
           do n=0,nprocs-1
              call MPI_Cart_coords(cart_comm,n,2,cart_coords,ier)
              qy_sg(n) = (cart_coords(0)*config%qy)/n1 + 1
              qy_eg(n) = ((cart_coords(0)+1)*config%qy)/n1
              qz_sg(n) = (cart_coords(1)*config%qz)/n2 + 1
              qz_eg(n) = ((cart_coords(1)+1)*config%qz)/n2
           enddo

           !     Set local boundaries
           
           nymoms = nymomsg(myproc)
           nymome = nymomeg(myproc)
           nymom_proc = nymome-nymoms+1
           
           qy_s = qy_sg(myproc)
           qy_e = qy_eg(myproc)
           qy_proc = qy_e-qy_s+1
           
           nzmoms = nzmomsg(myproc)
           nzmome = nzmomeg(myproc)
           nzmom_proc = nzmome-nzmoms+1
           
           qz_s = qz_sg(myproc)
           qz_e = qz_eg(myproc)
           qz_proc = qz_e-qz_s+1



        else ! mpi is not used
     
           !     Set local boundaries

           do n=0,nprocs-1
              nymomsg(n) = 1 
              nymomeg(n) = config%nymom
              nzmomsg(n) = 1
              nzmomeg(n) = config%nzmom
           enddo
           
           do n=0,nprocs-1
              qy_sg(n) = 1
              qy_eg(n) = config%qy
              qz_sg(n) = 1
              qz_eg(n) = config%qz
           enddo
           
           nymoms = 1
           nymome = config%nymom
           nymom_proc = nymome-nymoms+1
           
           qy_s = 1
           qy_e = config%qy
           qy_proc = qy_e-qy_s+1
           
           nzmoms = 1
           nzmome = config%nzmom
           nzmom_proc = nzmome-nzmoms+1
           
           qz_s = 1
           qz_e = config%qz
           qz_proc = qz_e-qz_s+1
        endif ! use_mpi
        
     endif ! spherical_eddington_factor

#if 0
     if (use_mpi) then
       if (use_1neighbour_comm) then
           if (qy_s+3 .gt. qy_e .or. qz_s+3 .gt. qz_e) then 
              write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
        if (use_2neighbour_comm) then
           if (qy_s+1 .gt. qy_e .or. qz_s+1 .gt. qz_e) then 
              write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
        if (use_4neighbour_comm) then
           if (qy_s .ne. qy_e .or. qz_s .ne. qz_e) then 
              write(6,*) "myproc", myproc, qy_s, qy_e, qz_s, qz_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
     endif
#endif
     
  else ! config%nsdim .ne. 3


     !2D-Version: it should not be necessary to reorder the processes
     !remark: this is also the 1d non mpi version!!!
     
     if (use_mpi) then
        
        cart_dims(0)=nprocs
        cart_dims(1)=1
        
        cart_periods(0)=.false. !at present, the boundary conditions need not be
        cart_periods(1)=.true.  !hardcoded, because ppm.par is read only after
                                !init_mpi is called. If possible, this should be 
                                !changed in the future

        cart_reorder=.false.

        call MPI_Cart_create(MPI_COMM_WORLD,2,cart_dims,cart_periods,cart_reorder, &
                                cart_comm,ier)
        call MPI_Comm_group(cart_comm,cart_group,ier)
	
        call mpi_comm_rank(cart_comm,myproc,ier)
        
        do n=0,nprocs-1
           !    call MPI_Cart_coords(cart_comm,n,2,cart_coords,ier)
           !    nymomsg(n) = (cart_coords(0)*config%qy)/nprocs + 1
           !    nymomeg(n) = ((cart_coords(0)+1)*config%qy)/nprocs
           nymomsg(n) = (n*config%nymom)/nprocs + 1
           nymomeg(n) = ((n+1)*config%nymom)/nprocs
        enddo
        
        do n=0,nprocs-1
           qy_sg(n) = (n*config%qy)/nprocs + 1
           qy_eg(n) = ((n+1)*config%qy)/nprocs
        enddo
        
        !     Set local boundaries
        
        nymoms = nymomsg(myproc)
        nymome = nymomeg(myproc)
        nymom_proc = nymome-nymoms+1
        
        qy_s = qy_sg(myproc)
        qy_e = qy_eg(myproc)
        qy_proc = qy_e-qy_s+1
        
        qz_s = 1
        qz_e = 1
        qz_proc = 1
        
        nzmoms = 1
        nzmome = 1
        nzmom_proc = 1
        
     else ! mpi not used
        
        if (config%use_spherical_eddington_factor) then
           !     Set local boundaries
           
           do n=0,nprocs-1
              nymomsg(n) = 0 ! is this correct ?????????????
              nymomeg(n) = config%nymom 
           enddo
           
           do n=0,nprocs-1
              qy_sg(n) = 1
              qy_eg(n) = config%qy
              qz_sg(n) = 1
              qz_eg(n) = config%qz
           enddo
           
           nymoms = 1  ! is this correct ?????????????
           nymome = config%nymom
           nymom_proc = nymome-nymoms+1
           
           qy_s = 1
           qy_e = config%qy
           qy_proc = qy_e-qy_s+1
           
           qz_s = 1
           qz_e = 1
           qz_proc = 1
           
           nzmoms = 1
           nzmome = 1
           nzmom_proc = 1
           
        else ! spherical_eddington_factor is not used
           
           !     Set local boundaries
           
           do n=0,nprocs-1
              nymomsg(n) = 1
              nymomeg(n) = config%nzmom
           enddo
           
           do n=0,nprocs-1
              qy_sg(n) = 1
              qy_eg(n) = config%qy
              qz_sg(n) = 1
              qz_eg(n) = config%qz
           enddo
           
           nymoms = 1
           nymome = config%nymom
           nymom_proc = nymome-nymoms+1
           
           qy_s = 1
           qy_e = config%qy
           qy_proc = qy_e-qy_s+1
           
           qz_s = 1
           qz_e = 1
           qz_proc = 1
           
           nzmoms = 1
           nzmome = 1
           nzmom_proc = 1
           
        endif ! spherical_eddington_factor
        
     endif ! use_mpi
     
     
     if (use_mpi) then
        if (use_1neighbour_comm) then
           if (qy_s+3 .gt. qy_e) then
              write(6,*) "myproc", myproc, qy_s, qy_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
        if (use_2neighbour_comm) then
           if (qy_s+1 .gt. qy_e) then
              write(6,*) "myproc", myproc, qy_s, qy_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
        if (use_4neighbour_comm) then
           if (qy_s .ne. qy_e) then
              write(6,*) "myproc", myproc, qy_s, qy_e
              call stopit_mpi("init_mpi: error: condition for sweeps()&
                & (calculation of gpocntr_utmp and gpocntr_ltmp) not satisfied")
           end if
        endif
     endif
     
  endif ! problem_dimensionality



  if (use_mpi) then
     ! find out on which host each MPI task is running  
     
     call MPI_GET_PROCESSOR_NAME( my_procname, my_procname_length, ier )

     if (ier.ne.MPI_SUCCESS) then
        write(6,'(1X,A)') "*  MPI_GET_PROCESSOR_NAME failed"
        write(6,'(1X,A)') "*  running on localhost"
        my_procname_length=9
        my_procname="localhost"
     end if


     ! get MPI taks's pid and its parents ppid

     call get_process_id(my_pid, my_ppid)

    
     ! on some machines .e.g. Juelich Juropa MPI_GET_PROCESSOR_NAME does
     ! not only give the host name where the MPI-Task is running, but
     ! also the pid. Since thus the total length of the string might differ
     ! on every task we need a global maximum


     call MPI_allreduce(my_procname_length, procname_length_buf, 1, &
                        MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)

     my_procname_length=procname_length_buf


     !     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     
     if (myproc .eq. 0) then
        procname_lengthg(myproc)                          = my_procname_length
        hostnames(myproc)%procnameg(1:my_procname_length) = my_procname(1:my_procname_length)
        
         pidg(myproc) = my_pid
        ppidg(myproc) = my_ppid

     endif


     ! communicate pid, ppid, procname, and procname_length to myproc = 0
     
     do i_task=1, nprocs-1

        if ((myproc .eq. i_task) .or. (myproc.eq.0)) then
           tag=i_task

           if (myproc .eq. i_task) then 
              call MPI_SEND(my_procname_length, 1, MPI_INTEGER, 0, &
                            tag, MPI_COMM_WORLD, ier)

           endif
           
           if (myproc .eq. 0) then
              call MPI_RECV(procname_length_buf, 1, MPI_INTEGER, i_task, &
                            tag, MPI_COMM_WORLD, status, ier)

              procname_lengthg(i_task)=procname_length_buf

           endif
        
           tag = 2*i_task
           
           if (myproc .eq. i_task) then
              call MPI_SEND(my_procname, my_procname_length, MPI_CHARACTER, 0, &
                            tag, MPI_COMM_WORLD, ier)
           endif

           if (myproc .eq. 0) then
              
              call MPI_RECV(procname_buf, my_procname_length, &
                            MPI_CHARACTER, i_task, tag, MPI_COMM_WORLD,  &
                            STATUS, ier)
              hostnames(i_task)%procnameg(1:procname_lengthg(i_task))=procname_buf

           endif

        
           tag = 3*i_task
           
           if (myproc .eq. i_task) then
              call MPI_SEND(my_pid, 1, MPI_INTEGER, 0, &
                            tag, MPI_COMM_WORLD, ier)
           endif

           if (myproc .eq. 0) then
              call MPI_RECV(pid_buf, 1, MPI_INTEGER, i_task, &
                            tag, MPI_COMM_WORLD, status, ier)

              pidg(i_task)=pid_buf
           endif

        
           tag = 4*i_task
           
           if (myproc .eq. i_task) then
              call MPI_SEND(my_ppid, 1, MPI_INTEGER, 0, &
                            tag, MPI_COMM_WORLD, ier)
           endif

           if (myproc .eq. 0) then
              call MPI_RECV(ppid_buf, 1, MPI_INTEGER, i_task, &
                            tag, MPI_COMM_WORLD, status, ier)

              ppidg(i_task)=ppid_buf
           endif

        endif

     enddo

     ! now communciate global arrays back to all processes  
     call MPI_BCAST(procname_lengthg, nprocs, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     call MPI_BCAST(pidg, nprocs, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
     call MPI_BCAST(ppidg, nprocs, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     do i_task=0, nprocs-1
        call MPI_BCAST(hostnames(i_task)%procnameg(1:my_procname_length), my_procname_length, MPI_CHARACTER, 0, MPI_COMM_WORLD, ier)
        
     enddo
     
     call MPI_BARRIER(MPI_COMM_WORLD, ier)


  endif ! use_mpi


  ! now set the value q for hydro arrays

  ! we have to take the value such that the biggest
  ! sweep in every direction (plus ghostzones) will fit it.
  ! In most of the cases this will be the radial sweep
  config%q = max(config%qx, qy_proc, qz_proc) + 20

  ! However if you use MPI a problem arises here:
  ! we have only qy_proc angular rays per MPI Task (e.g. 4)
  ! however the final INDEXING might be larger than q now
  !
  ! 1..4,   5..8,   9..12,  ...  1000..1004
  ! Task 1, Task 2, Task 3, ...  Task 256
  !
  ! That is why we have to introduce an offset logic
  ! which computes offsets 


  if (config%qx .eq. qy_proc .and. config%qx .eq. qz_proc) then
     ! special case: the follwoing logic does not allow
     ! this at the moment
     raise_abort("init_mpi.F90: q logic wrong")
  endif

  if ( (config%q-20) .eq. config%qx) then
     ! Normal case: we have per MPI Task more radial zones
     ! than "rays"
     config%q_xoff = 0
     config%q_yoff = qy_s - 1
     config%q_zoff = qz_s - 1
   end if

  if ( (config%q-20) .eq. qy_proc) then
     ! We have more angular rays in theta-direction (per MPI Task) than
     ! radial zones or phi-rays
     config%q_xoff = config%qx-1
     config%q_yoff = 0
     config%q_zoff = qz_s - 1
  endif


  if ( (config%q-20) .eq. qz_proc) then
     ! We have more angular rays in phi-direction (per MPI Task) than
     ! radial zones or theta-rays
     config%q_xoff = config%qx   -1
     config%q_yoff = qy_s  -1
     config%q_zoff = 0
  endif

end subroutine init_mpi

end module init_mpi_mod

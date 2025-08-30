module mpi_domains
  private 
  public :: print_mpi_setup, optimal_decomposition
  
contains

!>
!> \verbatim
!> This subroutine prints the MPI setup, i.e. which MPI task
!> covers which part of the theta grid
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>

  

subroutine print_mpi_setup

  use precision
  use mo_mpi
  use phycon, only : pc_pi
  use totare_hy, only : yzltot, yzrtot, zzltot, zzrtot

!  use charac, only : domainfile

#ifdef CHECK_THREAD_AFFINITY
  use thread_affinity
#endif

  use configure
  implicit none

  integer(kind=ik) :: i

  if (myproc .eq. 0) then
     open(97,file=config%domain_file,form= 'formatted')
  endif
  
  if (use_mpi) then

     if (myproc .eq. 0) then
        write (*,*) " "
        write (*,*) " ======================================================="
        write (*,*) "           MPI - Domain - Decomposition INFO            "
        write (*,*) " "
        write (*,*) "    In total ",nprocs," tasks "
        write (*,*) " "
        do i=0,nprocs-1
           write(*,'("   Task ",i4," with pid ",i8," (ppid ",i8,")")') i,pidg(i),ppidg(i)
           write(*,'("   running on   ",a)') trim(adjustl(hostnames(i)%procnameg(1:procname_lengthg(i))))

        enddo
        write (*,*) " "
     endif

     if (config%nsdim .eq. 3) then
        if (myproc .eq. 0) then
           do i=0,nprocs-1
              write (*,*)  "   Task ",i, "reigns rays ", qy_sg(i)," to ",qy_eg(i)," and ", &
                   qz_sg(i)," to ", qz_eg(i)
              write (*,*)  "        ", "from ", yzltot(qy_sg(i))/pc_pi*180._rk, " to ",    &
                   yzrtot(qy_eg(i))/pc_pi*180._rk, " degrees (theta direction)"
              write (*,*)  "        ", "from ", zzltot(qz_sg(i))/pc_pi*180._rk, " to ",    &
                   zzrtot(qz_eg(i))/pc_pi*180._rk, " degrees (phi direction)" 
           enddo
           
           write (*,*) " "
           write (*,*) " ======================================================="
           write (*,*) " "

        endif
        
  
        if (myproc .eq. 0) then
           write(97,*) " "
           write(97,*) " ======================================================="
           write(97,*) "           MPI - Domain - Decomposition INFO            "
           write(97,*) " "
           write(97,*) "    In total ",nprocs," tasks "
           do i=0,nprocs-1
              write(97,'("   Task ",i4," with pid ",i8," (ppid ",i8,") running on   ",a)') i,pidg(i), &
                   ppidg(i),trim(adjustl(hostnames(i)%procnameg(1:procname_lengthg(i))))
           enddo
           write(97,*) " "
        
           do i=0,nprocs-1
              write(97,*) "   Task ",i, "reigns rays ", qy_sg(i)," to ",qy_eg(i)," and ", qz_sg(i)," to ", qz_eg(i)
              write(97,*) "        ", "from ", yzltot(qy_sg(i))/pc_pi*180._rk, " to ", &
                yzrtot(qy_eg(i))/pc_pi*180._rk, " degrees (theta direction)"
              write(97,*) "        ", "from ", zzltot(qz_sg(i))/pc_pi*180._rk, " to ", &
                   zzrtot(qz_eg(i))/pc_pi*180._rk, " degrees (phi direction)" 
           enddo
        
           write(97,*) " "
           write(97,*) " ======================================================="
           write(97,*) " "
        
        endif
     endif ! config%nsdim


     if (config%nsdim .eq. 2) then

        if (myproc .eq. 0) then
           
           do i=0,nprocs-1
              write (*,*)  "   Task ",i, "reigns rays ", qy_sg(i)," to ",qy_eg(i)
              write (*,*)  "        ", "from ", yzltot(qy_sg(i))/pc_pi*180._rk, " to ", &
              yzrtot(qy_eg(i))/pc_pi*180._rk, " degrees"
              
           enddo
           
           write (*,*) " "
           write (*,*) " ======================================================="
           write (*,*) " "
           
        endif

        
        if (myproc .eq. 0) then
           write(97,*) " "
           write(97,*) " ======================================================="
           write(97,*) "           MPI - Domain - Decomposition INFO            "
           write(97,*) " "
           write(97,*) "    In total ",nprocs," tasks "
           do i=0,nprocs-1
              write(97,'("   Task ",i4," with pid ",i8," (ppid ",i8,") running on   ",a)') i,pidg(i),ppidg(i),&
                   trim(adjustl(hostnames(i)%procnameg(1:procname_lengthg(i))))
           enddo
           write(97,*) " "
        
           do i=0,nprocs-1
              write(97,*) "   Task ",i, "reigns rays ", qy_sg(i)," to ",qy_eg(i)
              write(97,*) "        ", "from ", yzltot(qy_sg(i))/pc_pi*180._rk, " to ", &
                yzrtot(qy_eg(i))/pc_pi*180._rk, " degrees (theta direction)"
           enddo
        
           write(97,*) " "
           write(97,*) " ======================================================="
           write(97,*) " "
        endif

     endif ! config%nsdim
        
  endif ! use_mpi


#ifdef CHECK_THREAD_AFFINITY
  call print_thread_affinity
#endif



  close(97)

end subroutine print_mpi_setup

subroutine optimal_decomposition (ny,nz,nprocs,n1,n2)
  
  use precision

  implicit none

  integer (kind=ik), intent (in) :: ny,nz,nprocs
  integer (kind=ik), intent (out) :: n1,n2
  integer (kind=ik) :: i

  real (kind=rk) :: r0,r1,ropt

  !ratio of ny and ny, the ratio of n1 and n2 should be as close as possible
  r0=real(max(ny,nz),kind=rk)/real(min(ny,nz),kind=rk) 

  ropt=abs(log(real(nprocs+1,kind=rk)/r0))
  do i=1,int(sqrt(real(nprocs+1,kind=rk)))
     if ((nprocs/i)*i .eq. nprocs) then !i is a factor of nprocs
        r1=real(nprocs/i,kind=rk)/real(i,kind=rk)
        r1=abs(log(r1/r0))
        if (r1 .lt. ropt) then
           ropt=r1
           n1=nprocs/i
           n2=i
        end if
     end if
  end do

  if (nz .gt. ny) then !swap n1 and n2
     n1=n2
     n2=nprocs/n1
  end if

  if (min(n1,n2) .eq. 1 .and. nprocs .gt. 8 .and. min(ny,nz) .gt. 4) then
     write (*,*) "You may have chosen a very inefficient MPI setup!"
  end if

  return

end subroutine optimal_decomposition

end module mpi_domains

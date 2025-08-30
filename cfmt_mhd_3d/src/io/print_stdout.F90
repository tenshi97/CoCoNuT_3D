module print_stdout_mod


  implicit none

  private
  public printit_taskX, printit_alltasks, print_memory_alloc



  interface printit_taskX

     module procedure printit_taskX_1msg,             &
                      printit_taskX_1msg_logical,     &
                      printit_taskX_1msg_i4,          &
                      printit_taskX_1msg_r8,          &
                      printit_taskX_1msg_array_1d_r8, &
                      printit_taskX_1msg_i4_r8,       &
                      printit_taskX_1msg_r8_r8,       &
                      printit_taskX_1msg_r8_r8_r8,    &
                      printit_taskX_1msg_r8_r8_r8_r8, &
                      printit_taskX_1msg_i4_i4,       &
                      printit_taskX_1msg_i4_i4_r8_r8, &
                      printit_taskX_2msg,             &
                      printit_taskX_2msg_r8,          &
                      printit_taskX_2msg_i4,          &
                      printit_taskX_3msg,             &
                      printit_taskX_3msg_r8,          &
                      printit_taskX_3msg_i4,          &
                      printit_taskX_3msg_logical,     &
                      printit_taskX_4msg,             &
                      printit_taskX_msg_i4_msg,       &
                      printit_taskX_msg_i4_msg_r8,    &
                      printit_taskX_msg_r8_msg_r8,    &
                      printit_taskX_msg_i4_msg_i4_msg_i4

  end interface printit_taskX

  interface printit_alltasks
     module procedure printit_alltasks_1msg

  end interface printit_alltasks


  contains

    subroutine print_memory_alloc(memloc, memglob, modname)

      use precision
      use mpi_vertex, only : myproc

      use configure
      implicit none

      real(kind=rk), intent(in) :: memloc, memglob
      character*(*)             :: modname


      if (config%use_print_memAlloc) then
         if (.not.(config%use_deactivate_stdout)) then
            if (.not.(config%use_stdout_redirect)) then
               if (myproc .eq. 0) then
                  print *," "
                  print '(a,i0,a,a)','Task ',myproc,' In module ',modname
                  write (*,'(" Task ",i0," Locally  allocated ",1pe12.5," [Mb]")') myproc,memloc
                  write (*,'(" Task ",i0," Globally allocated ",1pe12.5," [Mb]")') myproc, memglob
               endif
            else
               print '(a)'  ," "
               print '(a,a)'," In module ",modname
               write (*,'(" Locally  allocated ",1pe12.5," [Mb]")') memloc
               write (*,'(" Globally allocated ",1pe12.5," [Mb]")') memglob
            endif
         endif
      endif

    end subroutine print_memory_alloc


    subroutine printit_taskX_1msg(task,message)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none


      integer(kind=ik) :: task
      character*(*)    :: message

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i4,a,a)','Task ',myproc," ",message
               endif
            else
               print '(a)',message
            endif
         endif
      endif

    end subroutine printit_taskX_1msg

    subroutine printit_taskX_1msg_logical(task,message,log1)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      logical          :: log1

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,L1)','Task ',myproc," ",message," ",log1
               endif
            else
               print *,               message," ",log1
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_logical

    subroutine printit_taskX_1msg_i4(task, message, val1)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      integer(kind=ik) :: val1


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i0)','Task ',myproc," ",message," ",val1
               endif
            else
               print '(a,a,i0)',message," ",val1
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_i4

    subroutine printit_taskX_1msg_r8(task, message, val1)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      real(kind=rk)    :: val1


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then

               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,1pe12.5)','Task ',myproc," ",message," ",val1
               endif
            else
               print '(a,a,1pe12.5)',message," ",val1
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_r8

    subroutine printit_taskX_1msg_i4_r8(task, message, val1,val2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      integer(kind=ik) :: val1
      real(kind=rk)    :: val2

      if (.not.(config%use_deactivate_stdout)) then

         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i0,a,1pe12.5)','Task ',myproc," ",message," ",val1," ",val2
               endif
            else
               print '(a,a,i0,a,1pe12.5)',message," ",val1," ",val2
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_i4_r8

    subroutine printit_taskX_1msg_array_1d_r8(task, message, val1)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      real(kind=rk)    :: val1(:)

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print *,"Task ",myproc," ",message," ",val1(:)
               endif
            else
               print *,               message," ",val1(:)
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_array_1d_r8

    subroutine printit_taskX_1msg_i4_i4(task, message, val1, val2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      integer(kind=ik) :: val1, val2


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i0,a,i0)','Task ',myproc," ",message," ",val1," ",val2
               endif
            else
               print '(a,a,i0,a,i0)',message," ",val1," ",val2
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_i4_i4


    subroutine printit_taskX_1msg_i4_i4_r8_r8(task, message, val1, val2, &
                                              val3, val4)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      integer(kind=ik) :: val1, val2
      real(kind=rk) :: val3, val4


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i0,a,i0,a,1pe12.5,a,1pe12.5)','Task ',myproc," ",message," ",val1," ",val2," ",val3," ",val4
               endif
            else
               print '(a,a,i0,a,i0,a,1pe12.5,a,1pe12.5)',message," ",val1," ",val2," ",val3," ",val4
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_i4_i4_r8_r8

    subroutine printit_taskX_1msg_r8_r8(task, message, val1, val2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      real(kind=rk)    :: val1, val2

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,1pe12.5,a,1pe12.5)','Task ',myproc," ",message," ",val1," ",val2
               endif
            else
               print '(a,a,1pe12.5,a,1pe12.5)',message," ",val1," ",val2
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_r8_r8

    subroutine printit_taskX_1msg_r8_r8_r8(task, message, val1, val2, val3)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      real(kind=rk)    :: val1, val2, val3


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,1pe12.5,a,1pe12.5,a,1pe12.5)','Task ',myproc," ",message," ",val1," ",val2," ",val3
               endif
            else
               print '(a,a,1pe12.5,a,1pe12.5,a,1pe12.5)',message," ",val1," ",val2," ",val3
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_r8_r8_r8

    subroutine printit_taskX_1msg_r8_r8_r8_r8(task, message, val1, val2, val3, val4)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message
      real(kind=rk)    :: val1, val2, val3, val4

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,1pe12.5,a,1pe12.5,a,1pe12.5,a,1pe12.5)','Task ', &
                       myproc," ",message," ",val1," ",val2," ",val3," ",val4
               endif
            else
               print '(a,a,1pe12.5,a,1pe12.5,a,1pe12.5,a,1pe12.5)',message," ",val1," ",val2," ",val3," ",val4
            endif
         endif
      endif

    end subroutine printit_taskX_1msg_r8_r8_r8_r8

    subroutine printit_taskX_2msg(task,message1,message2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a)',"Task ",myproc," ",message1, " ",message2
               endif
            else
               print '(a,a,a)',message1, " ", message2
            endif
         endif
      endif
    end subroutine printit_taskX_2msg


    subroutine printit_taskX_2msg_r8(task,message1,message2, val)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2
      real(kind=rk)    :: val


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a,a,1pe12.5)','Task ',myproc," ",message1, " ",message2," ", val
               endif
            else
               print '(a,a,a,a,1pe12.5)',message1, " ",message2," ", val
            endif
         endif
      endif
    end subroutine printit_taskX_2msg_r8

    subroutine printit_taskX_2msg_i4(task,message1,message2, val)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2
      integer(kind=ik)    :: val

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a,a,i0)','Task ',myproc," ",message1, " ", message2," ", val
               endif
            else
               print '(a,a,a,a,i0)', message1," ", message2," ", val
            endif
         endif
      endif
    end subroutine printit_taskX_2msg_i4

    subroutine printit_taskX_3msg(task,message1,message2,message3)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2, message3

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a)','Task ',myproc," ",message1, message2, message3
               endif
            else
               print '(a,a,a,a,a)',message1," ", message2," ", message3
            endif
         endif
      endif

    end subroutine printit_taskX_3msg

    subroutine printit_taskX_3msg_logical(task,message1,message2, message3,val)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2, message3
      logical          :: val

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print *,"Task ",myproc," ",message1," ", message2," ", message3," ", val
               endif
            else
               print *,               message1," ", message2," ", message3," ", val
            endif
         endif
      endif
      
    end subroutine printit_taskX_3msg_logical


    subroutine printit_taskX_3msg_r8(task,message1,message2, message3,val)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2, message3
      real(kind=rk)    :: val

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a,a,a,a,1pe12.5)',"Task ",myproc," ",message1," ", message2," ", message3," ", val
               endif
            else
               print '(a,a,a,a,a,a,1pe12.5)',message1," ", message2," ", message3," ", val
            endif
         endif
      endif

    end subroutine printit_taskX_3msg_r8

    subroutine printit_taskX_3msg_i4(task,message1,message2, message3, val)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2, message3
      integer(kind=ik)    :: val

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a,a,a,a,i0)','Task ',myproc," ",message1," ", message2," ", message3," ", val
               endif
            else
               print '(a,a,a,a,a,a,i0)',message1," ", message2," ", message3," ", val
            endif
         endif
      endif
    end subroutine printit_taskX_3msg_i4


    subroutine printit_taskX_4msg(task,message1,message2,message3,message4)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2, message3, message4

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,a,a,a,a,a)','Task ',myproc," ",message1," ", message2," ", message3," ", message4
               endif
            else
               print '(a,a,a,a,a,a,a)',message1," ", message2," ",message3," ",message4
            endif
         endif
      endif

    end subroutine printit_taskX_4msg

    subroutine printit_taskX_msg_i4_msg(task,message1,val1,message2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2
      integer(kind=ik)    :: val1

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,i0,a)','Task ',myproc," ",message1, val1, message2
               endif
            else
               print '(a,i0,a)', message1, val1, message2
            endif
         endif
      endif

    end subroutine printit_taskX_msg_i4_msg

    subroutine printit_taskX_msg_i4_msg_r8(task,message1,val1,message2,val2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2
      integer(kind=ik) :: val1
      real(kind=rk)    :: val2

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i4,a,a,a,1pe12.5)','Task ',myproc," ",message1," ", val1," ", message2," ", val2
               endif
            else
               print '(a,a,i4,a,a,a,1pe12.5)',message1," ", val1," ", message2," ",val2
            endif
         endif
      endif

    end subroutine printit_taskX_msg_i4_msg_r8


    subroutine printit_taskX_msg_r8_msg_r8(task,message1,val1,message2,val2)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none

      integer(kind=ik) :: task
      character*(*)    :: message1, message2
      real(kind=rk)    :: val1, val2


      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,1pe12.5,a,a,a,1pe12.5)','Task ',myproc," ",message1," ", val1," ", message2," ", val2
               endif
            else
               print '(a,a,1pe12.5,a,a,a,1pe12.5)',message1," ", val1," ", message2," ",val2
            endif
         endif
      endif
      
    end subroutine printit_taskX_msg_r8_msg_r8

    subroutine printit_taskX_msg_i4_msg_i4_msg_i4(task,msg1,val1,msg2,val2,msg3,val3)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none
      integer(kind=ik) :: task, val1,val2,val3
      character*(*)    :: msg1, msg2, msg3

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            if (.not.(config%use_stdout_redirect)) then
               if (task .eq. myproc) then
                  print '(a,i0,a,a,a,i0,a,a,a,i0,a,a,a,i0)','Task ',myproc," ",msg1," ",val1," ",msg2," ",val2," ",msg3," ",val3
               endif
            else
               print '(a,a,i0,a,a,a,i0,a,a,a,i0)',msg1," ",val1," ",msg2," ",val2," ",msg3," ",val3
            endif
         endif
      endif
    end subroutine printit_taskX_msg_i4_msg_i4_msg_i4

    subroutine printit_alltasks_1msg(task,message)

      use precision
      use mpi_vertex, only : myproc
      use configure
      implicit none
      integer(kind=ik) :: task
      character*(*)    :: message

      if (.not.(config%use_deactivate_stdout)) then
         if (task .ne. -1) then
            print '(a,i0,a,a)','Task ',myproc," ",message
         endif
      endif
    end subroutine printit_alltasks_1msg

end module print_stdout_mod

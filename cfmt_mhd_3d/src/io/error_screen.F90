!>
!> This module provides the function show_error_screen, which is
!> used to print error information
!>
!> \author A. Marek, MPA, March 2009
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
module error_screen
  use precision

  public :: show_error_screen

  interface show_error_screen

     module procedure show_val1_arr,show_val2_arr,show_val3_arr,     &
                      show_val4_arr,show_val5_arr,                   &
                      show_string_val0,show_string_index,            &
                      show_string_index2, show_string_val1_arr,      &
                      show_string_val2_arr,show_string_val3_arr,     &
                      show_string_val4_arr,show_string_val5_arr,     &
                      show_string_val1_float,show_string_val2_float, &
                      show_string_val3_float,show_string_val4_float, &
                      show_string_val5_float

  end interface


  contains

  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param index  index of array where error occured
  !> \param val    array with diagnostic values
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_val1_arr(index,val)
    use precision

    implicit none

    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val(:)
    integer(kind=ik) :: i

    do i=1,index
       print *,i,val(i)
    enddo

  end subroutine show_val1_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param index  index of array where error occured
  !> \param val1    array with diagnostic values
  !> \param val2    array with diagnostic values
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_val2_arr(index,val1,val2)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:)
    integer(kind=ik) :: i

    do i=1,index
       print *,i,val1(i),val2(i)
    enddo

  end subroutine show_val2_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param index  index of array where error occured
  !> \param val1    array with diagnostic values
  !> \param val2    array with diagnostic values
  !> \param val3    array with diagnostic values
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_val3_arr(index,val1,val2,val3)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:)

    integer(kind=ik) :: i

    do i=1,index
       print *,i,val1(i),val2(i),val3(i)
    enddo

  end subroutine show_val3_arr

  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param index  index of array where error occured
  !> \param val1    array with diagnostic values
  !> \param val2    array with diagnostic values
  !> \param val3    array with diagnostic values
  !> \param val4    array with diagnostic values
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_val4_arr(index,val1,val2,val3,val4)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:),val4(:)

    integer(kind=ik) :: i

    do i=1,index
       print *,i,val1(i),val2(i),val3(i),val4(i)
    enddo

  end subroutine show_val4_arr

  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param index  index of array where error occured
  !> \param val1    array with diagnostic values
  !> \param val2    array with diagnostic values
  !> \param val3    array with diagnostic values
  !> \param val4    array with diagnostic values
  !> \param val5    array with diagnostic values
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_val5_arr(index,val1,val2,val3,val4,val5)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:),val4(:),val5(:)

    integer(kind=ik) :: i

    do i=1,index
       print *,i,val1(i),val2(i),val3(i),val4(i),val5(i)
    enddo

  end subroutine show_val5_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val0(string1,string2)
    use precision

    implicit none

    character*(*),intent(in) :: string1,string2

    print *," "
    print *,trim(string1)," ",trim(string2)

  end subroutine show_string_val0


  subroutine show_string_index(string1,string2,index)
    use precision

    implicit none

    character*(*),intent(in) :: string1,string2
    integer(kind=ik), intent(in) :: index

    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index

  end subroutine show_string_index


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index1 first index to print
  !> \param index2 second index to print
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_index2(string1,string2,index1,index2)
    use precision

    implicit none

    character*(*),intent(in) :: string1,string2
    integer(kind=ik), intent(in) :: index1,index2

    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index1,index2

  end subroutine show_string_index2


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val     array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val1_arr(string1,string2,index,val)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val(:)
    
    character*(*),intent(in) :: string1,string2
    integer(kind=ik) :: i

    print *," "
    print *,trim(string1)," ",trim(string2)
    
    do i=1,index
       print *,i,val(i)
    enddo
      
  end subroutine show_string_val1_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val2_arr(string1,string2,index,val1,val2)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:)

    character*(*),intent(in) :: string1,string2
    integer(kind=ik) :: i

    print *," "
    print *,trim(string1)," ",trim(string2)

    do i=1,index
       print *,i,val1(i),val2(i)
    enddo

  end subroutine show_string_val2_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val3_arr(string1,string2,index,val1,val2,val3)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:)
    
    character*(*),intent(in) :: string1,string2
    integer(kind=ik) :: i

    print *," "
    print *,trim(string1)," ",trim(string2)
    
    do i=1,index
       print *,i,val1(i),val2(i),val3(i)
    enddo

  end subroutine show_string_val3_arr


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !> \param val4    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val4_arr(string1,string2,index,val1,val2,val3,val4)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:),val4(:)
    
    character*(*),intent(in) :: string1,string2
    integer(kind=ik) :: i

    print *," "
    print *,trim(string1)," ",trim(string2)
    
    do i=1,index
       print *,i,val1(i),val2(i),val3(i),val4(i)
    enddo

  end subroutine show_string_val4_arr
  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !> \param val4    array with diagnostic information
  !> \param val5    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val5_arr(string1,string2,index,val1,val2,val3,val4,val5)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1(:),val2(:),val3(:),val4(:),val5(:)
    
    character*(*),intent(in) :: string1,string2
    integer(kind=ik) :: i

    print *," "
    print *,trim(string1)," ",trim(string2)
    
    do i=1,index
       print *,i,val1(i),val2(i),val3(i),val4(i),val5(i)
    enddo

  end subroutine show_string_val5_arr

  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val     array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val1_float(string1,string2,index,val)
    use precision
    
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val

    character*(*),intent(in) :: string1,string2

    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index,val

      
  end subroutine show_string_val1_float


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val2_float(string1,string2,index,val1,val2)
    use precision

    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1,val2

    character*(*),intent(in) :: string1,string2

    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index,val1,val2

  end subroutine show_string_val2_float


  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val3_float(string1,string2,index,val1,val2,val3)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1,val2,val3
    
    character*(*),intent(in) :: string1,string2
    
    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index,val1,val2,val3

  end subroutine show_string_val3_float
  !> One instance of show_error_screen
  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !> \param val4    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val4_float(string1,string2,index,val1,val2,val3,val4)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1,val2,val3,val4
    
    character*(*),intent(in) :: string1,string2
    
    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index,val1,val2,val3,val4

  end subroutine show_string_val4_float

  !> \author A. Marek, MPA, March 2009
  !>
  !> \param string1 first string to print
  !> \param string2 second string to print
  !> \param index   index to print
  !> \param val1    array with diagnostic information
  !> \param val2    array with diagnostic information
  !> \param val3    array with diagnostic information
  !> \param val4    array with diagnostic information
  !> \param val5    array with diagnostic information
  !>
  !> \verbatim
  !>   SVN - Information  
  !>   $Revision$
  !>   $Date$
  !>   
  !> \endverbatim
  !>
  subroutine show_string_val5_float(string1,string2,index,val1,val2,val3,val4,val5)

    use precision
    implicit none
    integer(kind=ik),intent(in) :: index
    real(kind=rk),intent(in) :: val1,val2,val3,val4,val5
    
    character*(*),intent(in) :: string1,string2
    
    print *," "
    print *,trim(string1)," ",trim(string2)
    print *,index,val1,val2,val3,val4,val5

  end subroutine show_string_val5_float

end module error_screen

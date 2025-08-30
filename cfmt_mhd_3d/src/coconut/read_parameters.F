c     ==========================================================================

      MODULE gr_parameters

c     ==========================================================================
      
      IMPLICIT NONE

      CONTAINS



c     ==========================================================================

      SUBROUTINE read_parameters

c     ==========================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE perm_aux_cfc
      USE boundary_cfc

      integer (kind=ik) :: l


ccc   Runge-Kutta index
ccc   1: 1st order
ccc   2: 2nd order
ccc   3: 3rd order
      call get_integer_parameter ("rk_order=", rk_order)
      if (rk_order.lt.2.or.rk_order.gt.3) then
         print*
         print*, 'Nonfatal problem in:'
         print*, '   Parameter reader'
         print*
         print*, 'Diagnosis:'
         print*, '   Runge-Kutta order outside of range'
         print*
         print*, 'Action:'
         print*, '   Resetting Runge-Kutta order to 2; continuing'
         print*
         rk_order = 2
      endif

ccc   conversion_tolerance: tolerance for conversion to native fluid variables
ccc   determines maximal change in P per time step
      call get_real_parameter ("conversion_tolerance=", 
     &     conversion_tolerance)


ccc   tolerance for delta_u in metric solver
      call get_real_parameter ("metric_tolerance=",
     &     metric_tolerance)
 
      call get_integer_parameter ("newtonian_potential_form=",
     &     newtonian_potential_form)
      if (newtonian_potential_form.lt.1.or.
     &     newtonian_potential_form.gt.9) then
         print*
         print*, 'Fatal problem in:'
         print*, '   Parameter reader'
         print*
         print*, 'Diagnosis:'
         print*, '   Value for Newtonian potential form out of ' //
     &        'parameter range'
         print*
         print*, 'Action:'
         print*, '   Setting value for Newtonian potential form ' //
     &        'to 1; continuing'
         print*
         newtonian_potential_form = 1
      endif

      return

      END SUBROUTINE read_parameters

c     ===================================================================================



c     ===================================================================================

      SUBROUTINE get_string_parameter (name_string, value)

c     ===================================================================================

      USE precision
      USE abort

      IMPLICIT NONE

      character*(200) :: line_string
      character*(50) :: name
      character*(*) :: name_string, value
      integer (kind=ik) :: i, j, l, ll



      i = index (name_string, '=')
ccc   index of equal sign
      if(i.gt.0)then
ccc   name of parameter
         name = name_string (1 : i - 1)
      else
         name = name_string
      endif

ccc   store unoccupied unit handle in unit
      open (unit = 9, file = './parameters', status = 'unknown')
ccc   open paramter file
 10   format (a)
ccc   character string format

 15   continue
      read (9, 10, end = 19, err = 15) line_string
ccc   read line with format 10, goto 19 at end, goto 15 (redo) at error, store in line
      i = index (line_string, '=')
ccc   index of equal sign
      j = index (line_string, '#')
ccc   index of comment char

      if (i.eq.0.or.j.eq.1) goto 15
ccc   if the whole line is a comment or there is no
ccc   equal sign, then go on to the next line    

      if(j.gt.0.and.j.lt.i) goto 15
ccc   if there is an equal sign, but it is in a comment
ccc   then go on to the next line

      l = len (line_string)
ccc   if there is a comment in the line, make sure
ccc   we do not include it as part of the value

      if (j.gt.0) l = j - 1
ccc   strip out leading blanks

      ll = 0
 17   continue
      ll = ll + 1
      if (line_string (ll : ll).eq.' ') goto 17
      if (name.eq.line_string (ll : i - 1)) then
         value = line_string (i + 1 : l)
         goto 20
      endif
      goto 15

 19   continue
      print*
      print*, 'Fatal problem in:'
      print*, '   Parameter reader'
      print*
      print*, 'Diagnosis:'
      print*, '   Could not read parameter ', name
      print*
      print*, 'Action:'
      print*, '   Exiting'
      print*
      raise_abort("Could not read parameter '" // trim(name) // "'")

 20   continue

      close (9)

      return
      
      END SUBROUTINE get_string_parameter
c     ===================================================================================




c     ===================================================================================

      SUBROUTINE get_real_parameter (name, value)

c     ===================================================================================

      USE precision
      USE abort

      IMPLICIT NONE

      character*(*) :: name
      character*(50) :: name_tmp
      character*(200) :: value_string
      real (kind=rk) :: value
      integer (kind=ik) :: i

      call get_string_parameter (name, value_string)

ccc   make sure that there is a decimal point in the string
      if (index (value_string, '.').eq.0) then
      i = index (name, '=')
      if (i.eq.0) i = len (name)
         name_tmp = name (1 : i - 1)
         print*
         print*, 'Fatal problem in:'
         print*, '   Parameter reader'
         print*
         print*, 'Diagnosis:'
         print*, '   Could not find decimal point in real parameter ' //
     &        name_tmp
         print*
         print*, 'Action:'
         print*, '   Exiting'
         print*
         raise_abort("Could not find decimal point in real parameter")
      endif

ccc   use internal file to convert string to real
      read (value_string, 12) value
 12   format (e20.15)

      return

      END SUBROUTINE get_real_parameter

c     ===================================================================================



c     ===================================================================================

      SUBROUTINE get_integer_parameter (name, value)

c     ===================================================================================

      USE precision

      IMPLICIT NONE

      character*(*) :: name
      character*(50) :: value_string
      integer (kind=ik) :: value



      call get_string_parameter (name, value_string)

ccc   use internal file to convert string to integer
      read (value_string, 13) value
 13   format (i10)

      return

      END SUBROUTINE get_integer_parameter

c     ===================================================================================


c     ===================================================================================

      SUBROUTINE get_logical_parameter (name, value)

c     ===================================================================================

      IMPLICIT NONE

      character*(*) name
      character*(50) value_string
      logical value



      call get_string_parameter (name, value_string)
      read (value_string, 14) value
 14   format (l10)

      return

      END SUBROUTINE get_logical_parameter

c     ==========================================================================



      END MODULE gr_parameters

c     ==========================================================================

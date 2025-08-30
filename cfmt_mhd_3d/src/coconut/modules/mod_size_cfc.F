c     ==================================================================
      
      MODULE size_cfc

c     ==================================================================

      USE precision
      
      IMPLICIT NONE
      
      SAVE

      integer (kind=ik) :: m, n, o

      integer (kind=ik) :: nleg 


      integer (kind=ik) :: n_s, n_e, o_s, o_e, n_loc, o_loc

      integer (kind=ik) :: n_eq

      CONTAINS


c     ==================================================================
      
      SUBROUTINE init_size_cfc

c     ==================================================================

      USE configure

      IMPLICIT NONE

      m = config%qx

      n = config%qy

      o = config%qz

      if (config%nsdim .ge. 2) then
         nleg = 10
      else
         nleg = 0
      endif

      return

      END SUBROUTINE init_size_cfc

c     ==================================================================

      
      END MODULE size_cfc

c     ==================================================================

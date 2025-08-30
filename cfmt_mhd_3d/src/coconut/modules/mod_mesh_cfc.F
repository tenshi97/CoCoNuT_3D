c     ==================================================================

      MODULE mesh_cfc

c     ==================================================================

      USE precision
      USE size_cfc

      IMPLICIT NONE

      SAVE

      INTEGER (KIND=ik), ALLOCATABLE :: clevel (:), clevel_region (:)
      INTEGER (KIND=ik), ALLOCATABLE :: j_i (:), j_f (:) 
      INTEGER (KIND=ik) :: nregions, region_s, region_e
      
c#ifdef MPI_HYDRO
      LOGICAL, ALLOCATABLE :: need_intracell (:)
      INTEGER(KIND=ik), ALLOCATABLE :: left_root (:)
      INTEGER(KIND=ik), ALLOCATABLE :: right_root (:)

c     MPI communicators & groups:
      INTEGER(KIND=ik), ALLOCATABLE :: comm_super_cell (:)
      INTEGER(KIND=ik), ALLOCATABLE :: group_super_cell (:)
c#endif


      END MODULE mesh_cfc

c     ==================================================================

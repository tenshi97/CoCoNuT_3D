c     ==================================================================

      MODULE ancient_cfc

c     ==================================================================
  
      USE precision
      USE configure

      IMPLICIT NONE
  
      SAVE

      real(kind=rk), allocatable :: epsanc(:,:,:),
     &     preanc(:,:,:), wloanc(:,:,:) 

      real(kind=rk), allocatable :: d_cap_hat_anc(:,:,:),
     &     s_1_hat_anc(:,:,:), s_2_hat_anc(:,:,:), 
     &     s_3_hat_anc(:,:,:), tau_hat_anc(:,:,:), 
     &     d_cap_xnu_hat_anc(:,:,:,:)
      
      real(kind=rk), allocatable ::
     &     phi_anc(:,:,:),
     &     alpha_anc(:,:,:),
     &     beta_up_1_anc(:,:,:), 
     &     beta_up_2_anc(:,:,:), 
     &     beta_up_3_anc(:,:,:)


      CONTAINS
!>
!> \par This subroutine allocates the arrays from module savare_hy
!>  \author B. Mueller
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 1083 $
!>   $Date: 2010-05-11 14:48:13 +0200 (Tue, 11 May 2010) $
!>   
!> \endverbatim
!>  
c     ==================================================================

      SUBROUTINE allocate_ancient_cfc

c     ==================================================================

      USE precision
      USE abort
      
      USE mo_mpi
      
      USE configure

      IMPLICIT NONE
      
      integer(kind=ik) :: istat

      allocate(
     &     epsanc(config%qx,qy_s:qy_e,qz_s:qz_e),
     &     preanc(config%qx,qy_s:qy_e,qz_s:qz_e),
     &     wloanc(config%qx,qy_s:qy_e,qz_s:qz_e), stat=istat)
      
      allocate(
     &     d_cap_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e),
     &     s_1_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e),
     &     s_2_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e), 
     &     s_3_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e),
     &     tau_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e), 
     &     d_cap_xnu_hat_anc(config%qx,qy_s:qy_e,qz_s:qz_e,config%qn), 
     &     stat=istat)
      
      allocate(
     &     phi_anc(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),
     &     alpha_anc(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),
     &     beta_up_1_anc(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),
     &     beta_up_2_anc(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),
     &     beta_up_3_anc(0:config%qx+1,qy_s-1:qy_e+1,qz_s-1:qz_e+1),
     &     stat=istat)

      if (istat .ne. 0) then
         raise_abort("Allocation failed!.")
      end if

      
      END SUBROUTINE allocate_ancient_cfc

c     ==================================================================




c     ==================================================================

      SUBROUTINE deallocate_ancient_cfc
      
c     ==================================================================
      
      USE precision
      USE mo_mpi
      
      IMPLICIT NONE
      
      integer(kind=ik) :: istat

      deallocate(epsanc, preanc, wloanc)

      deallocate(d_cap_hat_anc, s_1_hat_anc, 
     &     s_2_hat_anc, s_3_hat_anc,
     &     tau_hat_anc, d_cap_xnu_hat_anc, stat=istat)

      deallocate(phi_anc, alpha_anc, beta_up_1_anc, beta_up_2_anc,
     &    beta_up_3_anc, stat=istat)

      return
      
      END SUBROUTINE deallocate_ancient_cfc

c     ==================================================================




      END MODULE ancient_cfc

c     ==================================================================

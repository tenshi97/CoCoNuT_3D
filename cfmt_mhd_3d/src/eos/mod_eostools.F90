!>
!> \verbatim
!> this module provides the some subroutines, which
!> are needed for the EoS implementation in the Vertex code
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
module eostools
  use precision
  use abort

  implicit none
  private
  public mapnuc, my_pack, my_pack_xi

#undef MAPNUC_OLD
  
contains
!> \verbatim
!> Purpose:
!>    - maps the typical nucleus (A,Z)_LS to two real nuclei (A,Z)_(1,2) 
!>        if possible and changes composition appropriately
!>
!>           Prescription: i)   out of all {(A,Z)} fulfilling
!>                                (1)      Z/(A-Z) <= [Z/(A-Z)]_LS
!>                              select (A,Z)_1 with
!>                                (2) dist(Z,Z_LS)=min
!>
!>                         ii)  out of all {(A,Z)} fulfilling
!>                                (1)      Z/(A-Z) > [Z/(A-Z)]_LS
!>                              select (A,Z)_2 with
!>                                (2) dist(Z,Z_LS)=min
!>                         iii) assume a gas OF SAME CHARGE AS BEFORE
!>                              of n,p,alfa,(A,Z)_(1,2)
!>                              with baryon-density n_b=rho/m_b
!>
!>                         iv)  n_n, n_p, n_alfa, n_e => n_A(1,2)
!>                    
!>            Advantage: binding energy changes very little
!>
!>    - leaves the composition  unchanged if no appropiate nucleus is 
!>          found from the collection introduced in module nucparam
!>
!> Note: since LS does not fulfil charge neutrality, the routine
!>       does not use Y_e to calculate "dne", but Z*X_h/A. This
!>       ensures that n_A(1,2) do not become negative. The fact that
!>       charge neutrality is not given should be fixed in LS;
!>       charge neutrality induces Y_e-X_p-0.5*X_a = Z*X_h/A.
!> \endverbatim
!>-----------------------------------------------------------------------
  subroutine mapnuc(xi,rho,xn,xp,xa,xh,z,a)
!DIR$ INLINEALWAYS mapnuc



    use precision
    use abort
    use configure
    use nucparam
!    use phycon ! forcheck
    
    implicit none
    
    real(kind=rk) ,intent(in) :: rho,xn,xp,xa,xh,z,a
    real(kind=rk) ,intent(inout) :: xi(:)

    real(kind=rk) :: dnn,dne,dnp,dnh,dna,dnb,adistmin,adist,rne_ls,rne
    real(kind=rk) :: adistmin2,deh
    integer(kind=ik) :: inu,inuc,nuc,inuc2,kn

    nuc=size(xi)
    
!    if ((a <= 8.0_rk) .or. (xh == 0._rk)) then
    if ((a <= 4.0_rk) .or. (xh == 0._rk)) then
       ! no heavy nuclei present (xh must be 0.0 if a=4)
       !         if (xh > 0.0) stop 'error mapnuc'
       xi(n_he4+1:nuc-1)=0.0
       xi(n_n)=xn+xh*(a-z)/a
#ifdef MAPNUC_OLD
       xi(n_n)=xn
#endif
       xi(n_p)=xp+xh*z/a
#ifdef MAPNUC_OLD
       xi(n_p)=xp
#endif
       xi(n_he4)=xa
!       if (abs(xi(nuc)-xi(2)-0.5*xi(3))>1e-7) write (*,*) abs(xi(nuc)-xi(2)-0.5*xi(3)),xi(1:3),xh,a,z
    else
         ! heavy nuclei present
       rne_ls=z/(a-z)
         
       adistmin=500._rk
       do inu=n_he4+1,nuc-1
          rne=pc_nuc(inu,1)/(pc_nuc(inu,2)-pc_nuc(inu,1))
          if (rne_ls >= rne) then
             adist=abs(a-pc_nuc(inu,2))
             if (adist < adistmin) then
                inuc=inu
                adistmin=adist
             endif
          endif
       enddo
       adistmin2=500._rk
       do inu=n_he4+1,nuc-1
          rne=pc_nuc(inu,1)/(pc_nuc(inu,2)-pc_nuc(inu,1))
          if (rne_ls < rne) then
             adist=abs(a-pc_nuc(inu,2))
             if (adist < adistmin2) then
                inuc2=inu
                adistmin2=adist
             endif
          endif
       enddo

#ifdef MAPNUC_OLD
       adistmin2=501._rk 
#endif
       if (adistmin < 500._rk) then
          if (adistmin2 < 500._rk) then
             xi(n_he4+1:nuc-1)=0.0_rk
             xi(n_n)=xn
             xi(n_p)=xp
             xi(n_he4)=xa
             dnb=rho/pc_mb
             dne=xi(nuc)*dnb
             dnn=xi(n_n)*rho/(1.*pc_mb)
             dnp=xi(n_p)*rho/(1.*pc_mb)
             dna=xi(n_he4)*rho/(4.*pc_mb)
             
             dnh=dnb-dnn-dnp-4.*dna
             !           deh=dne-dnp-2.*dna
             deh=z*xh/a*dnb
             
             xi(inuc2)=  (deh-dnh*pc_nuc(inuc,1)/pc_nuc(inuc,2))/        &
                  (pc_nuc(inuc2,1)/pc_nuc(inuc2,2)-                      &
                  pc_nuc(inuc,1)/pc_nuc(inuc,2)   )*pc_mb/rho
             xi(inuc) = -(deh-dnh*pc_nuc(inuc2,1)/pc_nuc(inuc2,2))/      &
                  (pc_nuc(inuc2,1)/pc_nuc(inuc2,2)-                      &
                  pc_nuc(inuc,1)/pc_nuc(inuc,2)   )*pc_mb/rho
             
             !cc TEST
             !            if ((xi(inuc2)<-1e-15).or.(xi(inuc)<-1e-15)) then
!               write (*,*) 'a',xi(inuc),xi(inuc2),xi(1),
!     &              xi(2),xi(3),xi(nuc),dnh*pc_mb/rho,deh*pc_mb/rho,
!     &              sum(xi(1:nuc-1)),pc_nuc(inuc,1),pc_nuc(inuc2,1),
!     &              pc_nuc(inuc,2),pc_nuc(inuc2,2),
!     &              xi(nuc)-(z/a*xh+xp+0.5*xa),a
!               write (*,*) 1.-sum(xi(1:3)),xh,z,a,rho
!               write (*,*) "negative mapnuc"
!            endif
!            if (sum(xi(1:nuc-1))>1.000001) then
!               write (*,*) 'b',xi(inuc),xi(inuc2),xi(1),
!     &              xi(2),xi(3),xi(nuc),dnh*pc_mb/rho,deh*pc_mb/rho,
!     &              sum(xi(1:nuc-1)),pc_nuc(inuc,1),pc_nuc(inuc2,1),
!     &              pc_nuc(inuc,2),pc_nuc(inuc2,2),
!     &              xi(nuc)-(z/a*xh+xp+0.5*xa)
!               write (*,*) "mapnuc too large"
!               stop
!            endif
!            if (sum(xi(1:nuc-1))<0.999999) then
!               write (*,*) 'c',xi(inuc),xi(inuc2),xi(1),
!     &              xi(2),xi(3),xi(nuc),dnh*pc_mb/rho,deh*pc_mb/rho,
!     &              sum(xi(1:nuc-1)),pc_nuc(inuc,1),pc_nuc(inuc2,1),
!     &              pc_nuc(inuc,2),pc_nuc(inuc2,2),
!     &              xi(nuc)-(z/a*xh+xp+0.5*xa)
!               write (*,*) "mapnuc large"
!               stop
!            endif
!
!            write (88,*) z*xh/a-pc_nuc(inuc,1)/pc_nuc(inuc,2)*xi(inuc)-
!     &                        pc_nuc(inuc2,1)/pc_nuc(inuc2,2)*xi(inuc2)

          else
             xi(n_he4+1:nuc-1)=0.0_rk
             xi(n_n)=xn
             xi(n_he4)=xa
             dnb=rho/pc_mb
             dne=xi(nuc)*dnb
             dnn=xi(n_n)*rho/(1.*pc_mb)
             dna=xi(n_p)*rho/(4.*pc_mb)
             
             dnh=(dnb-dnn-dne-2._rk*dna)/(pc_nuc(inuc,2)-pc_nuc(inuc,1))
             dnp=dne-pc_nuc(inuc,1)*dnh-2._rk*dna
             
             xi(inuc)=dnh*pc_nuc(inuc,2)*pc_mb/rho
             xi(n_p)   =dnp*pc_mb     /rho
          endif
          
       else
          if (config%restmass_version .gt. 1 ) then
             write(*,*) 'findnuc panic:',z,a,xi(nuc)
          endif
       endif
    endif
    
! -- attribute the error due to the dominant species, 
! in order to enforce SUM(X)=1.0
    kn = SUM(MAXLOC(xi(1:nuc-1)))  
    xi(kn) = 1.0_rk - ( SUM(xi(1:nuc-1)) - xi(kn) )

end subroutine mapnuc
!> \verbatim
!>
!>  Replace the Fortran intrinsic function pack, which turned
!> out to be buggy on some systems
!>
!> Author : A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  function my_pack(data_arr,log_arr) result(result_arr)

    use precision
    use abort
    implicit none


    real(kind=rk) :: data_arr(:), result_arr(size(data_arr))
    logical       :: log_arr(:)

    integer(kind=ik) :: i,j


    if(size(data_arr) .ne. size(log_arr)) then
       raise_abort("my_pack(): size error")
    endif

    result_arr(:) = 0.0_rk

    j=1

    do i=1,size(data_arr)
       if (log_arr(i)) then
          result_arr(j)=data_arr(i)
          j=j+1
       endif
    enddo

  end function my_pack
!> \verbatim
!>
!>  Replace the Fortran intrinsic function pack, which turned
!> out to be buggy on some systems
!>
!> Author : A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  function my_pack_xi(data_arr,log_arr) result(result_arr)

    use precision
    use abort
    implicit none


    real(kind=rk) :: data_arr(:,:), result_arr(size(data_arr,dim=1),size(data_arr,dim=2))
    logical       :: log_arr(:)

    integer(kind=ik) :: i,j,k


    if(size(data_arr, dim=1) .ne. size(log_arr)) then
       raise_abort("my_pack(): size error")
    endif

    result_arr(:,:) = 0.0_rk

    do k=1,size(data_arr,dim=2)
       j=1

       do i=1,size(data_arr,dim=1)
          if (log_arr(i)) then
             result_arr(j,k)=data_arr(i,k)
             j=j+1
          endif
          
       enddo
    enddo

  end function my_pack_xi

!> \verbatim
!>
!>  Replace the Fortran intrinsic function unpack, which turned
!> out to be buggy on some systems
!>
!> Author : A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  function my_unpack(data_arr,log_arr) result(result_arr)

    use precision
    use abort
    implicit none


    real(kind=rk) :: data_arr(:), result_arr(size(data_arr))
    logical       :: log_arr(:)

    integer(kind=ik) :: i,j


    if(size(data_arr) .ne. size(log_arr)) then
       raise_abort("my_pack(): size error")
    endif


    j=1

    do i=1,size(data_arr)
       if (log_arr(i)) then
          result_arr(j)=data_arr(i)
          j=j+1
       endif

    enddo


  end function my_unpack
  
end module eostools

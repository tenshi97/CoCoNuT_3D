#ifdef ONEMG_EOS
#define SWITCH_OFF_BOUNDARY_ERROR
#else
#undef SWITCH_OFF_BOUNDARY_ERROR
#endif

!>
!> \verbatim
!> this module provides the low-density NSE EoS
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
module low_density_nse
use precision
use abort

implicit none

private
public low_density_nse_eos

contains 
!> \verbatim
!>
!>  Compute low-density NSE -EoS
!>
!> Author : R. Buras
!> \endverbatim
!>
!>  \param dens density
!>  \param yel Ye
!>  \param temp temperature
!>  \param nnuc nuclear number densities
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine low_density_nse_eos(dens,yel,eden,temp,nnuc)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c TABLE for NSE composition in THJ regime, provided by THJ-KOK c
!c - Table needs den,ye,ene to calculate NSE-Composition, which c
!c   consists of n,p,alpha,and 11 symmetric and 3 asymmetric    c
!c   nuclei                                                     c
!c - the energy density is consistent with the eos_tj table     c
!c - check: Charge neutrality and mass normalization is taken care for c
!c - check: Warning! Composition can become as small as 1e-25 !!!      c
!c            Should be changed (evtl. in eos_tj)               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  use precision
  use abort
  use nucparam, only : n_n, n_p, n_he4, pc_nuc
  use phycon
  use eos_data_type
  use specfun, only : fastlog10
  use mo_mpi, only : myproc
  use print_stdout_mod
  
  use configure
  implicit none
  real(kind=rk), intent(in)               :: dens(:), yel(:), eden(:), temp(:)
  real(kind=rk), intent(inout)            :: nnuc(:,:) ! carefull nnuc as input
  ! and output is a number density, during interpolation, however, the mass
  ! fractions are stored in this array and later re-normalized into a number
  ! density

  real(kind=rk), dimension(lowden_table(1)%ntt) :: elarr, hearr, xharr

  real(kind=rk), dimension(lowden_table(1)%ntt,lowden_table(1)%nspec) :: &
                                                                      xiarr
  real(kind=rk), dimension(size(dens))    :: dens0, yel0, wr0, wr1, wy0,    &
                                             wy1, we0, we1, elden, logdens, &
                                             logeden
  integer(kind=ik), dimension(size(dens)) :: idens0, idens1, iyel0, &
                                             iyel1, ieden0, ieden1
  logical, dimension(size(dens)) :: innse

  integer(kind=ik)                :: ir, in, i, j, nn, NR, NT, NE, NQ


  if (config%low_den_nse_eos .gt. 2) then
     
     ir=size(dens)
     in=size(nnuc,dim=2)
      
     NR=lowden_table(1)%nro
     NT=lowden_table(1)%ntt
     NE=lowden_table(1)%nye
     NQ=lowden_table(1)%nspec


#ifdef CFC_TRANSPORT
     innse(:)=temp(:) .ge. config%tkok !.or. dens(:) .gt. 3.0e+7_rk
#else
     innse(:)=temp(:) .ge. config%tkok !no NSE below this temperature [MeV]
#endif


     logdens(:) = dens(:)
     call fastlog10(logdens, ir)

     logeden(:) = eden(:)
     call fastlog10(logeden, ir)

     do i=1,ir
        if (innse(i)) then
           ! check for kok-table boundaries
           if ((logdens(i) .lt. lowden_table(1)%romin) .or. &
               (logdens(i) .gt. lowden_table(1)%romax)) then
              write (*,*) "Zone: ",i, "rho: ",dens(i),lowden_table%romin, &
                       lowden_table%romax
#ifdef SWITCH_OFF_BOUNDARY_ERROR
              call printit_taskX(myproc,"Switching off NSE for this zone.")
              innse(i) = .false. 
#else
              raise_abort("low_density_nse_eos(): boundaries exceeded in density.")
#endif /* SWITCH_OFF_BOUNDARY_ERROR */
           endif

           if ((yel(i) .lt. lowden_table(1)%yemin) .or. &
                (yel(i) .gt.  lowden_table(1)%yemax)) then
              write (*,*) "Zone: ",i," ye: ",YEL(I),lowden_table(1)%yemin, &
                     lowden_table(1)%yemax
#ifdef SWITCH_OFF_BOUNDARY_ERROR
              call printit_taskX(myproc,"Switching off NSE for this zone.")
              innse(i) = .false. 
#else
              raise_abort("low_density_nse_eos(): boundaries exceeded in Ye.")
#endif /* SWITCH_OFF_BOUNDARY_ERROR */
           endif
           

           dens0(i) = max( 1.0_rk, min( real(NR,kind=rk), &
                1.0_rk+real((NR-1),kind=rk)*   &
                (logdens(i)-lowden_table(1)%romin) / &
                ((lowden_table(1)%romax)-            &
                (lowden_table(1)%romin)))) 

           yel0 (i) = max( 1.0_rk, min( real(NE,kind=rk), &
                1.0_rk+real((NE-1),kind=rk)*   &
                (       yel(i) -lowden_table(1)%yemin) / &
                (lowden_table(1)%yemax-lowden_table(1)%yemin))) 
           
           idens0(i) = int( dens0(i) )   ! calculate index of dens grid
           iyel0(i) = int(  yel0(i) )   ! calculate index of yel  grid
           
           idens1(i) = min(idens0(i) + 1,NR)
           iyel1(i) = min( iyel0(i) + 1,NE)

           wr1(i) = dens0(i) - real(idens0(i),kind=rk)
           wr0(i) = 1.0_rk - wr1(i)
           
           wy1(i) =  yel0(i) - real( iyel0(i),kind=rk)
           wy0(i) = 1.0_rk - wy1(i)
           
           elarr(:) = &
                wr0(i)* (  wy0(i)*real(lowden_table(1)%eltab(:,idens0(i),iyel0(i)),kind=rk)   &
                + wy1(i)*real(lowden_table(1)%eltab(:,idens0(i),iyel1(i)),kind=rk) ) &
                +wr1(i)* (  wy0(i)*real(lowden_table(1)%eltab(:,idens1(i),iyel0(i)),kind=rk)   &
                             + wy1(i)*real(lowden_table(1)%eltab(:,idens1(i),iyel1(i)),kind=rk) ) 

           if (config%low_den_nse_eos .eq. 4) then

              hearr(:)=      &
                   wr0(i)* (  wy0(i)*real(lowden_table(1)%hetab(:,idens0(i),iyel0(i)),kind=rk)  &
                   + wy1(i)*real(lowden_table(1)%hetab(:,idens0(i),iyel1(i)),kind=rk) )&
                   +wr1(i)* (  wy0(i)*real(lowden_table(1)%hetab(:,idens1(i),iyel0(i)),kind=rk)  &
                            + wy1(i)*real(lowden_table(1)%hetab(:,idens1(i),iyel1(i)),kind=rk) )
              
              xharr(:)=      &
                   wr0(i)* (  wy0(i)*real(lowden_table(1)%xhtab(:,idens0(i),iyel0(i)),kind=rk)  &
                   + wy1(i)*real(lowden_table(1)%xhtab(:,idens0(i),iyel1(i)),kind=rk) )&
                   + wr1(i)* (  wy0(i)*real(lowden_table(1)%xhtab(:,idens1(i),iyel0(i)),kind=rk)  &
                   + wy1(i)*real(lowden_table(1)%xhtab(:,idens1(i),iyel1(i)),kind=rk) )
           
           else

              xiarr(:,:)= &
                   wr0(i)* (  wy0(i)*real(lowden_table(1)%xitab(:,idens0(i),iyel0(i),:),kind=rk)   &
                   + wy1(i)*real(lowden_table(1)%xitab(:,idens0(i),iyel1(i),:),kind=rk) ) &
                   +wr1(i)* (  wy0(i)*real(lowden_table(1)%xitab(:,idens1(i),iyel0(i),:),kind=rk)   &
                   + wy1(i)*real(lowden_table(1)%xitab(:,idens1(i),iyel1(i),:),kind=rk) )
           endif ! low_den_nse_eos == 4

           elden(i) = logeden(i)

           ieden0(i)=1 ! this was missing, quick fix! check if works!
           
           do j=1,nt
              if (elarr(j).lt.elden(i)) ieden0(i)=j
           enddo
           ieden1(i)=ieden0(i)+1
           
           we1(i) = ( elden(i)         - elarr(ieden0(i)) ) / &
                ( elarr(ieden1(i)) - elarr(ieden0(i)) )
           we0(i) = 1.0_rk - we1(i)

           if ((we1(i) >= 0.0_rk) .and. (we1(i) <= 1.0_rk)) then
           
              nnuc(i,:) = 0.0_rk

              if (config%low_den_nse_eos .eq. 4) then
                 nnuc(i,n_he4) = &
                      we0(i)* hearr(ieden0(i)) + we1(i)* hearr(ieden1(i))
                 nnuc(i,lowden_table(1)%n_hv) = &
                      we0(i)* xharr(ieden0(i)) + we1(i)* xharr(ieden1(i))

                 ! neutrons and protons via charge neutrality
                 nnuc(i,n_n) = 1._rk - 0.5_rk       *nnuc(i,n_he4)        &
                      - (lowden_table(1)%aa-lowden_table(1)%zz)/ &
                      lowden_table(1)%aa*nnuc(i,lowden_table(1)%n_hv)         &
                      -  yel(i)

                 nnuc(i,n_p) = yel(i) - 0.5_rk  *nnuc(i,n_he4)            &
                      - lowden_table(1)%zz/                  &
                      lowden_table(1)%aa*nnuc(i,lowden_table(1)%n_hv)

                 ! composition >0
                 if (nnuc(i,n_n).lt.0._rk) then
                    nnuc(i,lowden_table(1)%n_hv)=(1._rk-yel(i)-0.5_rk*nnuc(i,n_he4))/        &
                         (lowden_table(1)%aa-lowden_table(1)%zz)*    &
                         lowden_table(1)%aa
                    nnuc(i,n_p)=1._rk-nnuc(i,lowden_table(1)%n_hv)-nnuc(i,n_he4)
                    nnuc(i,n_n)=0._rk
                    if (nnuc(i,lowden_table(1)%n_hv).lt.0._rk) then
                       nnuc(i,lowden_table(1)%n_hv)=0._rk
                       nnuc(i,n_he4)=2._rk*(1._rk-yel(i))
                       nnuc(i,n_p)=1._rk-nnuc(i,n_he4)
                    endif
                 endif
            
              else ! low_den_nse_eos == 4
                 nnuc(i,1:nq) = we0(i)* xiarr(ieden0(i),:) + we1(i)* xiarr(ieden1(i),:)
              endif ! low_den_nse_eos == 4


              ! tests
              if (minval(nnuc(i,:)) .lt. 0.0_rk)       raise_abort("low_density_nse_eos(): negnuc")
              if (abs(sum(nnuc(i,:))-1._rk) .gt. 1.e-13_rk) raise_abort("low_density_nse_eos(): sum error")


              if (config%low_den_nse_eos .eq. 4) then
                 if (ABS(nnuc(i,n_p)+0.5_rk*nnuc(i,n_he4)+lowden_table(1)%zz/lowden_table(1)%aa*nnuc(i,lowden_table(1)%n_hv) &
                      -yel(i)).gt.1.e-13_rk) then
                    raise_abort("low_density_nse_eos(): chargerr")
                 endif
              else

                 if (ABS(SUM(real(pc_nuc(:,1),kind=rk)/real(pc_nuc(:,2),kind=rk)*nnuc(i,:)) &
                      -yel(i)).gt.1.e-13_rk) then
                    raise_abort("low_density_nse_eos(): chargerr")
                 endif

              endif ! low_den_nse_eos
              
              ! renormalize
              do nn=1,in
                 nnuc(i,nn)=dens(i)*max(nnuc(i,nn),1.e-25_rk)/(pc_nuc(nn,2)*pc_mb)
              enddo
           else
              write(*,*) "In NSE with temp, out of NSE with energy for zone", i
           endif ! (we1(i) >= 0.0_rk) .and. (we1(i) <= 1.0_rk)

        endif ! innse(i)
     enddo

  endif ! low_den_nse_eos .gt. 2
     
end subroutine low_density_nse_eos

end module low_density_nse


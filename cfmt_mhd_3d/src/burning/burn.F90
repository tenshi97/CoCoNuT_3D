module burn_mod

  implicit none

contains

#undef SPLIT_BURN_LOOPS

!=======================================================================
  subroutine burn (jj, kk, selftime, childrentime)
!=======================================================================
! Autor          : Markus Rampp (MPA) 
! Modul          : $Id: sweeps.F,v 1.18 2005/09/28 11:33:25 rburas Exp $
! Version        : $Revision: 1.18 $
! Date           : $Date: 2005/09/28 11:33:25 $
!
! Purpose: compute change of composition due to nuclear burning
!           and recombination
!=======================================================================
    use precision
        
    use abort
    use bndinf_hy
    use intgrs_hy
    use hydro_hy
#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK)
    use gfloat_hy
#endif

    use vnew_hy
    use phycon
    use nucparam

#ifndef ONEMG_EOS
    use eos_sn2, only : lsrolo
#else
    use eos_sn2, only : lsrolo, lsrolo2
#endif

#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK)
    use burn_lim
    use cap
    use rateselec
    use nutrio_hy
#endif

#ifdef BURN_NETWORK
    use mod_get_rate
    use mod_tak
#endif      
    use cputim
    use configure
    use state, only : hydro

    implicit none
! LOCAL variables that are not in modules

    real(kind=rk), intent(out)          :: selftime(2), childrentime(2)
    real(kind=rk)                       :: selftime_start(2)
    integer(kind=ik), intent(in)        :: jj, kk
    integer(kind=ik)                    :: i, i4, j, k, kn
    real(kind=rk)                       :: t9, rhoc, xp_old, xn_old, rhonse

#ifdef BURN_SS
    real(kind=rk), dimension(size(rho)) :: vels
#endif
    real(kind=rk)                       :: dummy, dummy2, dummy3, dummy4, &
                                           dummy5, dummy6
#ifdef BURN_NETWORK
    real(kind=rk), dimension(size(rho)) :: sig, siga, sigb, &
                                           sigc, sigd, sige
    real(kind=rk)                       :: fac, fac2, fac3, z, z0, dz, eps
    logical, dimension(size(rho))       :: soliton, solitoni, lburn
    integer, dimension(size(rho))       :: iburn
    integer, dimension(size(rho),3)     :: i_burn_seq
    real(kind=rk), dimension(size(rho)) :: alph, alph1, alph2, alph3
    real(kind=rk), dimension(size(rho)) :: rho_burn, tem_burn, enu_burn
    real(kind=rk), dimension(size(rho),config%qn) :: xnu_burn
    real(kind=rk), dimension(size(rho),10) :: unu_tmp, tmp_tmp
    integer, dimension(size(rho))       :: soli_l, soli_r
    integer                             :: n_burn, burn_step
#endif
#ifdef TAK_RATES
    real(kind=rk), dimension(size(rho)) :: ecnes, ecmgs, ennes, enmgs
#endif

#ifdef TAK_SUBTIME
    real(kind=rk)                       :: dtrate
    real(kind=rk) zahs(size(rho),2),xnu(size(xnuc,dim=1), &
                  size(xnuc,dim=4))
    real(kind=rk), dimension(size(rho)) :: xhrepqs, ps, gamcs, ss, &
                                           cuqs, ceqs, cnqs, cpqs, &
                                           rhon, tem, eni
    logical                             :: err
#endif

    integer(kind=ik) :: flashing_index, last_intermediate_element
    real(kind=rk) :: escr(1:config%qx), offs

  selftime     = 0._rk
  childrentime = 0._rk
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

    if (config%restmass_version .gt. 0) escr(:)=energy(:,jj,kk)


    if (config%restmass_version .eq. 1) then
       do i=1,config%qx
!#if LOW_EOS_NSE_FLAG != 22
          offs = moffs + SUM(xnuc(i,jj,kk,1:n_he4)*mbar(1:n_he4)) &
                    + (1._rk-SUM(xnuc(i,jj,kk,1:n_he4)))*mbar(n_rep)
!#else
!     &               + (1.e0_rk-SUM(xnuc(i,jj,kk,1:n_he4)))*mbar(n_rep)
!     &             - (1.e0_rk-SUM(xnuc(i,jj,kk,n_d:n_he3)))*mbar(n_rep)
!#endif
          escr(i) = escr(i) + offs
       enddo

    endif


    if (config%restmass_version .eq. 2) then
       do i=1,config%qx
          offs = moffs + SUM(xnuc(i,jj,kk,1:config%qn-1)*mbar(1:config%qn-1))
          escr(i) = escr(i) + offs
       enddo

    endif


    if (config%restmass_version .eq. 3) then
       do i=1,config%qx
          offs = moffs + SUM(xnuc(i,jj,kk,1:config%qn)*mbar(1:config%qn))
          escr(i) = escr(i) + offs
       enddo

    endif

    do i = 1, nzn
       i4 = i + 4
#if defined(BURN_SS) || defined(NO_FLASHING_AT_SONIC_POINT)
       gamc (i4)  = gammac(i,jj,kk)
       p    (i4)  = press (i,jj,kk)
#endif
       rhonu(i4)   = densty(i,jj,kk)
       tmp  (i4)   = temp  (i,jj,kk)
#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK) || defined(NO_FLASHING_AT_SONIC_POINT)
       unu  (i4)   = velx  (i,jj,kk)
#endif
#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK)
       unu  (i4)   = velx  (i,jj,kk)
       enu  (i4)   = energy(i,jj,kk)
       xnnu(i4,config%qn) = xnuc(i,jj,kk,  config%qn)    !=Y_e Paco
#endif
       xnnu(i4,n_n)    = xnuc(i,jj,kk,n_n)
       xnnu(i4,n_p)    = xnuc(i,jj,kk,n_p)
       xnnu(i4,n_he4)  = xnuc(i,jj,kk,n_he4)
       xnnu(i4,n_si28) = xnuc(i,jj,kk,n_si28)
       xnnu(i4,n_ni56) = xnuc(i,jj,kk,n_ni56)
       xnnu(i4,n_c12)  = xnuc(i,jj,kk,n_c12)
       xnnu(i4,n_mg24) = xnuc(i,jj,kk,n_mg24)
       xnnu(i4,n_o16)  = xnuc(i,jj,kk,n_o16)
       xnnu(i4,n_ne20) = xnuc(i,jj,kk,n_ne20)
    enddo

#if defined(BURN_NETWORK) && !defined(BURN_NETW_NOA)

#define BURN_ASYM
#undef DISSOCIATION    
#define BURN_NE
#define BURN_MG

!******************************************************************
!     BURNING NETWORK (vectorized)
!******************************************************************

#ifdef ONEMG_EOS
    rhonse = lsrolo2
#else
    rhonse = lsrolo
#endif


    !     Soliton detection
    soliton(:)  = .false.
    solitoni(:) = .false.
    eps = 1.0e-5_rk

    soli_l(:) = 0_ik
    soli_r(:) = 0_ik
    
    unu_tmp(:,:)=spread(unu(:),dim=2,ncopies=10)
    kn = nzn4 - 5
    do j = -1,-5,-1
       do k = 1,5
          do i=10,kn
             if ( (1.0_rk-unu(i)/(unu_tmp(i+j,11+j)+eps)) .gt. 0.8_rk &
                   .and. (1.-unu(i)/(unu_tmp(i+k,k)+eps)) .gt. 0.8_rk &
                   .and. (unu(i)-unu_tmp(i+j,11+j))    .gt. 0.0_rk    &
                   .and. (unu(i)-unu_tmp(i+k,k))    .gt. 0.0_rk) then
                   solitoni(i) = .true.
                   soli_l(i)   =j
                   soli_r(i)   =k
                endif
             end do
          end do
       end do
       do i=10,kn
          if (solitoni(i)) then
             solitoni(i+soli_l(i)) = .true.
             solitoni(i+soli_r(i)) = .true.
          end if
       end do
       
!     This prevents from burning in zones that have higher 
!     temperatures due to an advection error.
       tmp_tmp(:,:)=spread(tmp(:),dim=2,ncopies=10)
       kn = nzn4 - 5
       do j = -1,-5,-1
          do k = 1,5
             do i = 10,kn
                if ( (tmp(i)/(tmp_tmp(i+j,11+j)+eps)-1.0_rk) .gt. 0.15_rk &
                   .and. (tmp(i)/(tmp_tmp(i+k,k)+eps)-1.0_rk) .gt. 0.15_rk    &
                   .and. (tmp(i)-tmp_tmp(i+j,11+j))    .gt. 0.0_rk       &
                   .and. (tmp(i)-tmp_tmp(i+k,k))    .gt. 0.0_rk          &
                   .and.  tmp(i)   .lt. thigh ) then
                   solitoni(i) = .true.
                   soli_l(i)   = j
                   soli_r(i)   = k
                endif
             enddo
          enddo
       enddo
       do i=10,kn
          if (solitoni(i)) then
             solitoni(i+soli_l(i)) = .true.
             solitoni(i+soli_r(i)) = .true.
          end if
       end do
       
       do i = 5,kn
          do k = i+1,i+5
             if (solitoni(i) .and. solitoni(k)) then
                do j = i,k
                   soliton(j) = .true.
                enddo
             endif
          enddo
       enddo
       

       lburn(:)=.false.

       !     gather arrays
       do i=5,nzn4
          lburn(i)=tmp(i)   .gt.  tlow           &
              .and. .not. soliton(i)             &
              .and. tmp(i)   .lt.  thigh         &
              .and. ishck(i-4,jj,kk) .lt.   1_ik &
              .and. rhonu(i)       .lt.   rhonse 
       end do



       n_burn=0
#ifdef NEC_COMPILER
!CDIR NODEP
#endif
       do i=5,nzn4
          if (lburn(i)) then
             n_burn           = n_burn+1
             rho_burn(n_burn) = rhonu(i)
             tem_burn(n_burn) = tmp(i)
#ifdef NEC_COMPILER
!CDIR EXPAND=25            
#endif
             xnu_burn(n_burn,1:config%qn) = xnnu(i,1:config%qn)
          end if
       end do

       if (n_burn .gt. 0_ik) then
!****************************************************************************
!     individual burning steps

! burn 12C and 12C to 4He and 20Ne 
!*     Y1+Y1-->Y2+Y3
!*     dY1/dt=-rate*Y1^2
!*     --> X1f=X1i/(1+X1i*rate*dt/A1)
!*     X2f=X2i+0.5*dX1*A2/A1
!*     X3f=X3i+0.5*dX1*A3/A1

          call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),1,0)
#ifdef BURN_ASYM
          call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         siga(1:n_burn),13,14)
          sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif 
         
          do i=1,n_burn
             if (xnu_burn(i,n_c12) .gt. 0.0_rk) then
                fac  = pc_nuc(n_c12,2)
                fac2 = 0.5_rk*pc_nuc(n_he4,2)/pc_nuc(n_c12,2)
                fac3 = 0.5_rk*pc_nuc(n_ne20,2)/pc_nuc(n_c12,2)
                
                dummy  = xnu_burn(i,n_he4)
                dummy2 = xnu_burn(i,n_ne20)
                dummy3 = xnu_burn(i,n_c12)            
                
                xnu_burn(i,n_c12) = xnu_burn(i,n_c12)/ &
                                  (1.0_rk+xnu_burn(i,n_c12)*(sig(i)/fac)*hydro%dt)        
               
                xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20)+ &
                                  fac3*(dummy3 - xnu_burn(i,n_c12))
                xnu_burn(i,n_he4) = xnu_burn(i,n_he4)+   &
                                  fac2*(dummy3 - xnu_burn(i,n_c12))
             end if
          end do
          


! burn 12C and 16O to 4He and 24Mg
!*     Y1+Y2-->Y3+Y4
!*     dY1/dt=-rate*Y1*Y2
!*     dY2/dt=-rate*Y1*Y2
!*     dY1/dt=dY2/dt=-dY3/dt=-dY4/dt
!*     z=Y1+Y2
!*     d2z/dt2=-rate*z*dz/dt
!*     -->zf=zi/(1+0.5*rate*zi*dt)
!*     X1f=X1i-0.5*dz*A1
!*     X2f=X2i-0.5*dz*A2
!*     X3f=X3i+0.5*dz*A3
!*     X4f=X4i+0.5*dz*A4

          call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),2,0)
#ifdef BURN_ASYM
          call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         siga(1:n_burn),15,16)
          sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif 

          do i = 1, n_burn
             if (tem_burn(i) .gt. 1.5e+9_rk .and.      &
                 xnu_burn(i,n_c12) .gt. 0.0_rk .and.   &
                 xnu_burn(i,n_o16) .gt. 0.0_rk) then
               
                dummy3 = xnu_burn(i,n_he4)
                dummy4 = xnu_burn(i,n_c12)
                dummy5 = xnu_burn(i,n_o16)
                dummy6 = xnu_burn(i,n_mg24)
                
                z0 = xnu_burn(i,n_c12)/pc_nuc(n_c12,2)+ &
                     xnu_burn(i,n_o16)/pc_nuc(n_o16,2)
!       z0 = min(2.*xnu_burn(i,n_c12)/pc_nuc(n_c12,2),2 .*xnu_burn(i,n_o16)/pc_nuc(n_o16,2))       
                z  = z0
                z  = z0/(1.0_rk+0.5_rk*sig(i)*z0*hydro%dt) 
                dz = z0 - z               
                xnu_burn(i,n_c12) =  &
                                  xnu_burn(i,n_c12) - 0.5_rk*dz*pc_nuc(n_c12,2)
                xnu_burn(i,n_o16) =  &
                                  xnu_burn(i,n_o16) - 0.5_rk*dz*pc_nuc(n_o16,2)
               
                xnu_burn(i,n_mg24) = &
                                   xnu_burn(i,n_mg24) + 0.5_rk*dz*pc_nuc(n_mg24,2)
                xnu_burn(i,n_he4) =  &
                                  xnu_burn(i,n_he4) + 0.5_rk*dz*pc_nuc(n_he4,2)

               
                if (xnu_burn(i,n_c12) .lt. 0.0_rk .or.  &
                    xnu_burn(i,n_o16) .lt. 0.0_rk) then
                   xnu_burn(i,n_he4)  = dummy3 
                   xnu_burn(i,n_c12)  = dummy4
                   xnu_burn(i,n_o16)  = dummy5  
                   xnu_burn(i,n_mg24) = dummy6 
                  
                   dz = min(2.0_rk*dummy4/pc_nuc(n_c12,2),2.* &
                            dummy5/pc_nuc(n_o16,2))        
                   xnu_burn(i,n_c12)  = dummy4 - 0.5_rk*dz*pc_nuc(n_c12,2)
                   xnu_burn(i,n_o16)  = dummy5 - 0.5_rk*dz*pc_nuc(n_o16,2)
                   
                   xnu_burn(i,n_mg24) = dummy6 + 0.5_rk*dz*pc_nuc(n_mg24,2)
                   xnu_burn(i,n_he4)  = dummy3 + 0.5_rk*dz*pc_nuc(n_he4,2)
                end if
                
             end if
          end do
         
#ifdef DISSOCIATION

! burn 20Ne to 4He and 16O
          call get_rate2(tem_burn(1:n_burn),rho_burn(1:n_burn),  &
                         sig(1:n_burn),21,22)           
         
          do i=1,n_burn
             if (xnu_burn(i,n_ne20) .gt. 0.0_rk) then
                fac2 = pc_nuc(n_he4,2)/pc_nuc(n_ne20,2)
                fac3 = pc_nuc(n_o16,2)/pc_nuc(n_ne20,2)
                
                dummy3 = xnu_burn(i,n_he4)
                dummy4 = xnu_burn(i,n_ne20)
                dummy5 = xnu_burn(i,n_o16)
                
                xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20)*exp(-sig(i)*hydro%dt)
                dummy = dummy4 - xnu_burn(i,n_ne20)
                
               
                xnu_burn(i,n_he4) = xnu_burn(i,n_he4) + dummy*fac2
                xnu_burn(i,n_o16) = xnu_burn(i,n_o16) + dummy*fac3
             endif
          end do
          
#endif /*DISSOCIATION*/


! burn 12C and 20Ne to 4He and 28Si
!*     Y1+Y2-->Y3+Y4
!*     dY1/dt=-rate*Y1*Y2
!*     dY2/dt=-rate*Y1*Y2
!*     dY1/dt=dY2/dt=-dY3/dt=-dY4/dt
!*     z=Y1+Y2
!*     d2z/dt2=-rate*z*dz/dt
!*     -->zf=zi/(1+0.5*rate*zi*dt)
!*     X1f=X1i-0.5*dz*A1
!*     X2f=X2i-0.5*dz*A2
!*     X3f=X3i+0.5*dz*A3
!*     X4f=X4i+0.5*dz*A4

                 
          call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),12,0)
#ifdef BURN_ASYM
          call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),17,18)
          sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif 

          do i=1,n_burn
             if (tem_burn(i) .gt. 1.5e+9_rk .and.    &
                 xnu_burn(i,n_c12)  .gt. 0.0_rk .and. &
                 xnu_burn(i,n_ne20) .gt. 0.0_rk) then
               
                dummy3 = xnu_burn(i,n_he4)
                dummy4 = xnu_burn(i,n_c12)
                dummy5 = xnu_burn(i,n_ne20)            
                dummy6 = xnu_burn(i,n_si28)
               
                z0 = xnu_burn(i,n_c12)/pc_nuc(n_c12,2)+   &
                     xnu_burn(i,n_ne20)/pc_nuc(n_ne20,2)
!            z0=min(2.*xnu_burn(i,n_c12)/pc_nuc(n_c12,2),
!     +           2 .*xnu_burn(i,n_ne20)/pc_nuc(n_ne20,2))      
                z  = z0/(1.0_rk+0.5_rk*sig(i)*z0*hydro%dt) 
                dz = z0 - z
                 
                xnu_burn(i,n_c12) =  &
                                 xnu_burn(i,n_c12) - 0.5_rk*dz*pc_nuc(n_c12,2)
                xnu_burn(i,n_ne20) = &
                                 xnu_burn(i,n_ne20) - 0.5_rk*dz*pc_nuc(n_ne20,2)
                 
                xnu_burn(i,n_si28) = &
                                 xnu_burn(i,n_si28) + 0.5_rk*dz*pc_nuc(n_si28,2)
                xnu_burn(i,n_he4) =  &
                                 xnu_burn(i,n_he4) + 0.5_rk*dz*pc_nuc(n_he4,2)
               
                if (xnu_burn(i,n_c12)  .lt. 0.0_rk .or. &
                    xnu_burn(i,n_ne20) .lt. 0.0_rk) then
                   xnu_burn(i,n_he4)  = dummy3 
                   xnu_burn(i,n_c12)  = dummy4
                   xnu_burn(i,n_ne20) = dummy5  
                   xnu_burn(i,n_si28) = dummy6 
                  
                   dz = min(2.0_rk*dummy4/pc_nuc(n_c12,2),  &
                            2.0_rk*dummy5/pc_nuc(n_ne20,2))       
                   xnu_burn(i,n_c12)  = dummy4 - 0.5_rk*dz*pc_nuc(n_c12,2)
                   xnu_burn(i,n_ne20) = dummy5 - 0.5_rk*dz*pc_nuc(n_ne20,2)
                   
                   xnu_burn(i,n_si28) = dummy6 + 0.5_rk*dz*pc_nuc(n_si28,2)
                   xnu_burn(i,n_he4)  = dummy3 + 0.5_rk*dz*pc_nuc(n_he4,2) 
                   
                endif
             endif
          end do
         
!     burn 16O to 28Si and 4He
!*     Y1+Y1-->Y2+Y3
!*     dY1/dt=-rate*Y1^2
!*     --> X1f=X1i/(1+X1i*rate*dt/A1)
!*     X2f=X2i+0.5*dX1*A2/A1
!*     X3f=X3i+0.5*dX1*A3/A1

          call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),3,0)
#ifdef BURN_ASYM
          call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         siga(1:n_burn),19,20)
          sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif 
         
          do i=1,n_burn
             if ( tem_burn(i) .gt. 1.5e+9_rk .and.   &
                  xnu_burn(i,n_o16) .gt. 0.0_rk) then
               
                fac  = pc_nuc(n_o16,2)
                fac2 = 0.5_rk*pc_nuc(n_he4,2)/pc_nuc(n_o16,2)
                fac3 = 0.5_rk*pc_nuc(n_si28,2)/pc_nuc(n_o16,2)
                
                dummy  = xnu_burn(i,n_he4)
                dummy2 = xnu_burn(i,n_si28)
                dummy3 = xnu_burn(i,n_o16)            
                
                xnu_burn(i,n_o16) = xnu_burn(i,n_o16)/ &
                                    (1.0_rk+xnu_burn(i,n_o16)*(sig(i)/fac)*hydro%dt)
               
                xnu_burn(i,n_si28) = xnu_burn(i,n_si28)+ &
                                    fac3*(dummy3 - xnu_burn(i,n_o16))
                xnu_burn(i,n_he4) = xnu_burn(i,n_he4)+ &
                                    fac2*(dummy3 - xnu_burn(i,n_o16))
             end if
          end do
          
#ifdef BURN_NE
! burn 20Ne to 28Si 
!*     Y1+Y1-->Y2
!*     dY1/dt=-rate*Y1^2
!*     --> X1f=X1i/(1+X1i*rate*dt/A1)
!*     X2f=X2i+dX1
      
          call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sig(1:n_burn),3,0)
         
#ifdef BURN_ASYM
          call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         siga(1:n_burn),19,20)
          sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif
         
          do i=1,n_burn
            if (tem_burn(i) .gt. 3.0e+9_rk .and. &
                xnu_burn(i,n_ne20) .gt. 0.0_rk) then
               
               fac  =pc_nuc(n_ne20,2)
               
               dummy2 = xnu_burn(i,n_si28)
               dummy3 = xnu_burn(i,n_ne20)
               
               xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20)/ &
                                    (1.0_rk+xnu_burn(i,n_ne20)*(sig(i)/fac)*hydro%dt)
               
               xnu_burn(i,n_si28) = xnu_burn(i,n_si28)+ &
                                    (dummy3 - xnu_burn(i,n_ne20))
            endif
         end do
#endif /*BURN_NE*/

#ifdef BURN_MG
! burn Mg24 to 28Si  
!*     Y1+Y1-->Y2
!*     dY1/dt=-rate*Y1^2
!*     --> X1f=X1i/(1+X1i*rate*dt/A1)
!*     X2f=X2i+dX1

         call get_rate3(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        sig(1:n_burn),3,0)
#ifdef BURN_ASYM
         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                       siga(1:n_burn),19,20)
         sig(1:n_burn) = sig(1:n_burn)+siga(1:n_burn)
#endif
         
         do i=1,n_burn
            if (tem_burn(i) .gt. 3.e+9_rk .and. &
                xnu_burn(i,n_mg24) .gt. 0.0_rk) then
               
               fac = pc_nuc(n_mg24,2)
               
               dummy2 = xnu_burn(i,n_si28)
               dummy3 = xnu_burn(i,n_mg24)
               
               xnu_burn(i,n_mg24) = xnu_burn(i,n_mg24)/ &
                                    (1.0_rk+xnu_burn(i,n_mg24)*(sig(i)/fac)*hydro%dt)
               
               xnu_burn(i,n_si28) = xnu_burn(i,n_si28)+ &
                                    (dummy3 - xnu_burn(i,n_mg24))
            endif
         end do
#endif /*BURN_MG*/

!******************************************************************
!     looks which of these reactions occur more rapidly and selects the
!     correct order!
!*     Y1+Y2-->Y3
!*     dY1/dt=-rate*Y1*Y2
!*     dY2/dt=-rate*Y1*Y2
!*     dY1/dt=dY2/dt=-dY3/dt
!*     z=Y1+Y2
!*     d2z/dt2=-rate*z*dz/dt
!*     -->zf=zi/(1+0.5*rate*zi*dt)
!*     X1f=X1i-0.5*dz*A1
!*     X2f=X2i-0.5*dz*A2
!*     X3f=X3i+0.5*dz*A3

! burn 16O and 4He to 20Ne
         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        siga(1:n_burn),8,9)
!     burn 20Ne and 4He to 24Mg
         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        sigb(1:n_burn),6,7)           
!     burn 24Mg and 4He to 28Si
         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                         sigc(1:n_burn),4,5)           
!******************************************************************
         
         do i=1,n_burn      
            if (siga(i) .gt. sigb(i) .and. sigb(i) .ge. sigc(i)) then
!     Fall 1:
               i_burn_seq(i,1) = 1
               i_burn_seq(i,2) = 2
               i_burn_seq(i,3) = 3
            else if (siga(i) .gt. sigc(i) .and. sigc(i) .gt. sigb(i)) then
               !     Fall 2:
               i_burn_seq(i,1) = 1
               i_burn_seq(i,2) = 3
               i_burn_seq(i,3) = 2
            else if (sigb(i) .ge. siga(i) .and. siga(i) .gt. sigc(i)) then
!     Fall 3:
               i_burn_seq(i,1) = 2
               i_burn_seq(i,2) = 1
               i_burn_seq(i,3) = 3
            else if (sigb(i) .ge. sigc(i) .and. sigc(i) .ge. siga(i)) then
!     Fall 4:
               i_burn_seq(i,1) = 2
               i_burn_seq(i,2) = 3
               i_burn_seq(i,3) = 1
            else if (sigc(i) .ge. siga(i) .and. siga(i) .gt. sigb(i)) then
!     Fall 5:
               i_burn_seq(i,1) = 3
               i_burn_seq(i,2) = 1
               i_burn_seq(i,3) = 2
            else if (sigc(i) .ge. sigb(i) .and. sigb(i) .ge. siga(i)) then
!     Fall 6:
               i_burn_seq(i,1) = 3
               i_burn_seq(i,2) = 2
               i_burn_seq(i,3) = 1
            end if
         end do


         do burn_step=1,3
!     Recalculate rates before each burning step (unnecessary if
!     the rates depend only on T and rho, are not changed in the
!     burning algorithm.
!         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn),
!     +        siga(1:n_burn),8,9)
!         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn),
!     +        sigb(1:n_burn),6,7)           
!         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn),
!     +        sigc(1:n_burn),4,5)           

            do i=1,n_burn
!     -------------------------------------------------------------------------------
               if (i_burn_seq(i,burn_step) .eq. 1) then
!     Reaction O16 + He4 -> Ne20:
                  if (xnu_burn(i,n_he4) .gt. 0.0_rk .and. &
                      xnu_burn(i,n_o16) .gt. 0.0_rk) then
                     
                     dummy3 = xnu_burn(i,n_he4)
                     dummy4 = xnu_burn(i,n_o16)
                     dummy5 = xnu_burn(i,n_ne20)
                     
                     z0 = xnu_burn(i,n_he4)/pc_nuc(n_he4,2)+ &
                          xnu_burn(i,n_o16)/pc_nuc(n_o16,2)
!     z=min(2.*xnu_burn(i,n_he4)/pc_nuc(n_he4,2),2*xnu_burn(i,n_o16)/pc_nuc(n_o16,2))
                     z = z0/(1.0_rk+0.5_rk*siga(i)*z0*hydro%dt) 
                     dz = z0-z
                     
                     xnu_burn(i,n_he4) = xnu_burn(i,n_he4)- &
                                         0.5_rk*dz*pc_nuc(n_he4,2)
                     xnu_burn(i,n_o16) = xnu_burn(i,n_o16)- &
                                         0.5_rk*dz*pc_nuc(n_o16,2)
                     
                     xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20)+ &
                                          0.5_rk*dz*pc_nuc(n_ne20,2)
                     
                     if (xnu_burn(i,n_he4) .lt. 0.0_rk .or. &
                         xnu_burn(i,n_o16) .lt. 0.0_rk) then
                        xnu_burn(i,n_he4) = dummy3
                        xnu_burn(i,n_o16)  = dummy4 
                        xnu_burn(i,n_ne20) = dummy5 
                        
                        dz = min(2.0_rk*dummy3/pc_nuc(n_he4,2), &
                                 2.0_rk*dummy4/pc_nuc(n_o16,2))
                        xnu_burn(i,n_he4) = dummy3- &
                                            0.5_rk*dz*pc_nuc(n_he4,2)
                        xnu_burn(i,n_o16) = dummy4- &
                                            0.5_rk*dz*pc_nuc(n_o16,2)
                        
                        xnu_burn(i,n_ne20) = dummy5+ &
                                             0.5_rk*dz*pc_nuc(n_ne20,2)        
                     endif
                  end if
#ifdef SPLIT_BURN_LOOPS
               end if
            end do
!     -------------------------------------------------------------------------------
         
            do i=1,n_burn
               if (i_burn_seq(i,burn_step) .eq. 2) then
#else            
               else if (i_burn_seq(i,burn_step) .eq. 2) then
#endif

!     Reaction Ne20 + He4 -> Mg24:
                  if (xnu_burn(i,n_he4)  .gt. 0.0_rk .and. &
                      xnu_burn(i,n_ne20) .gt. 0.0_rk) then
                     
                     dummy3 = xnu_burn(i,n_he4)
                     dummy4 = xnu_burn(i,n_ne20)
                     dummy5 = xnu_burn(i,n_mg24)            
                     
                     z0 = xnu_burn(i,n_he4)/pc_nuc(n_he4,2)+ &
                          xnu_burn(i,n_ne20)/pc_nuc(n_ne20,2)
!     z=min(2.*xnu_burn(i,n_he4)/pc_nuc(n_he4,2),2.*xnu_burn(i,n_ne20)/pc_nuc(n_ne20,2))
                     z = z0/(1.0_rk+0.5_rk*sigb(i)*z0*hydro%dt) 
                     dz = z0-z
                     
                     xnu_burn(i,n_he4) = xnu_burn(i,n_he4)- &
                                         0.5_rk*dz*pc_nuc(n_he4,2)
                     xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20)- &
                                          0.5_rk*dz*pc_nuc(n_ne20,2)
                     
                     xnu_burn(i,n_mg24) = xnu_burn(i,n_mg24)+ &
                                          0.5_rk*dz*pc_nuc(n_mg24,2)
                     
                     if (xnu_burn(i,n_he4)   .lt. 0.0_rk .or. &
                          xnu_burn(i,n_ne20) .lt. 0.0_rk) then
                         xnu_burn(i,n_he4) = dummy3
                         xnu_burn(i,n_ne20) = dummy4 
                         xnu_burn(i,n_mg24) = dummy5 
                        
                        dz = min(2.0_rk*dummy3/pc_nuc(n_he4,2), &
                                 2.0_rk*dummy4/pc_nuc(n_ne20,2))
                        xnu_burn(i,n_he4) = &
                                           dummy3-0.5_rk*dz*pc_nuc(n_he4,2)
                        xnu_burn(i,n_ne20) = &
                                           dummy4-0.5_rk*dz*pc_nuc(n_ne20,2)
                        xnu_burn(i,n_mg24) = &
                                           dummy5+0.5_rk*dz*pc_nuc(n_mg24,2)
                     end if
                  end if
#ifdef SPLIT_BURN_LOOPS
               end if
            end do
!     -------------------------------------------------------------------------------

            do i=1,n_burn
               if (i_burn_seq(i,burn_step) .eq. 3) then
#else
               else if (i_burn_seq(i,burn_step) .eq. 3) then
#endif

!     Reaction Mg24 + He4 -> Si28
                  if (xnu_burn(i,n_he4)  .gt. 0.0_rk .and. &
                      xnu_burn(i,n_mg24) .gt. 0.0_rk) then
                     
                     dummy3 = xnu_burn(i,n_he4)
                     dummy4 = xnu_burn(i,n_mg24)
                     dummy5 = xnu_burn(i,n_si28)            
                     
                     z0 = xnu_burn(i,n_he4)/pc_nuc(n_he4,2)+ &
                          xnu_burn(i,n_mg24)/pc_nuc(n_mg24,2)
!     z=min(2.*xnu_burn(i,n_he4)/pc_nuc(n_he4,2),2.*xnu_burn(i,n_mg24)/pc_nuc(n_mg24,2))
                     z = z0/(1.0_rk+0.5_rk*sigc(i)*z0*hydro%dt) 
                     dz = z0-z
                     
                     xnu_burn(i,n_he4)  = xnu_burn(i,n_he4)- &
                                          0.5_rk*dz*pc_nuc(n_he4,2)
                     xnu_burn(i,n_mg24) = xnu_burn(i,n_mg24)- &
                                          0.5_rk*dz*pc_nuc(n_mg24,2)
                     
                     xnu_burn(i,n_si28) = xnu_burn(i,n_si28)+ &
                                          0.5_rk*dz*pc_nuc(n_si28,2)
                     
                     if (xnu_burn(i,n_he4)  .lt. 0.0_rk .or. &
                         xnu_burn(i,n_mg24) .lt. 0.0_rk) then
                        xnu_burn(i,n_he4) = dummy3
                        xnu_burn(i,n_mg24) = dummy4 
                        xnu_burn(i,n_si28) = dummy5 
                        
                        dz = min(2.0_rk*dummy3/pc_nuc(n_he4,2),  &
                                 2.0_rk*dummy4/pc_nuc(n_mg24,2))
                        xnu_burn(i,n_he4) =                      &
                                          dummy3-0.5_rk*dz*pc_nuc(n_he4,2)
                        xnu_burn(i,n_mg24) = dummy4-0.5_rk*dz*   &
                                             pc_nuc(n_mg24,2)
                        
                        xnu_burn(i,n_si28) =                     &
                                              dummy5+0.5_rk*dz*pc_nuc(n_si28,2)
                     end if
                  end if          
               end if
            end do
         end do
!     -------------------------------------------------------------------------------

! burn 12C and 4He to 16O 

         call get_rate4(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        sig(1:n_burn),10,11)
         
         do i=1,n_burn
            if (xnu_burn(i,n_he4) .gt. 0.0_rk .and. &
                xnu_burn(i,n_c12) .gt. 0.0_rk) then
               
               dummy3 = xnu_burn(i,n_he4)
               dummy4 = xnu_burn(i,n_c12)
               dummy5 = xnu_burn(i,n_o16)
               
               z = xnu_burn(i,n_he4)/pc_nuc(n_he4,2)+ &
                   xnu_burn(i,n_c12)/pc_nuc(n_c12,2)
!     z=min(2.*xnu_burn(i,n_he4)/pc_nuc(n_he4,2),2.*xnu_burn(i,n_c12)/pc_nuc(n_c12,2))
               dummy = z 
               z = z/(1.0_rk+0.5_rk*sig(i)*z*hydro%dt) 
               z = dummy - z
               
               xnu_burn(i,n_he4) = xnu_burn(i,n_he4) - &
                                   0.5_rk*z*pc_nuc(n_he4,2)
               xnu_burn(i,n_c12) = xnu_burn(i,n_c12) -  &
                                   0.5_rk*z*pc_nuc(n_c12,2)
               
               xnu_burn(i,n_o16) = xnu_burn(i,n_o16) + &
                                   0.5_rk*z*pc_nuc(n_o16,2)
               
               if (xnu_burn(i,n_he4) .lt. 0.0_rk .or. &
                   xnu_burn(i,n_c12) .lt. 0.0_rk) then
                  xnu_burn(i,n_he4) = dummy3
                  xnu_burn(i,n_c12) = dummy4 
                  xnu_burn(i,n_o16) = dummy5 
               
                  z=min(2.0_rk*dummy3/pc_nuc(n_he4,2), &
                        2.0_rk*dummy4/pc_nuc(n_c12,2))
                  xnu_burn(i,n_he4) = dummy3 - 0.5_rk*z*pc_nuc(n_he4,2)
                  xnu_burn(i,n_c12) = dummy4 - 0.5_rk*z*pc_nuc(n_c12,2)
                  
                  xnu_burn(i,n_o16) = dummy5  + 0.5_rk*z*pc_nuc(n_o16,2)
               endif
            end if          
         end do
         
!******************************************************************
#ifdef DISSOCIATION
! burn 16O to 4He and 12C

         call get_rate2(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        sig(1:n_burn),23,24)           
         do i=1,n_burn
            if (xnu_burn(i,n_o16) .gt. 0.0_rk) then
               
               fac2              = pc_nuc(n_he4,2)/pc_nuc(n_o16,2)
               fac3              = pc_nuc(n_c12,2)/pc_nuc(n_o16,2)
               
               dummy3            = xnu_burn(i,n_he4)
               dummy4            = xnu_burn(i,n_o16)
               dummy5            = xnu_burn(i,n_c12)
               
               xnu_burn(i,n_o16) = xnu_burn(i,n_o16)*exp(-sig(i)*hydro%dt)
               dummy             = dummy4 - xnu_burn(i,n_o16)
               
               
               xnu_burn(i,n_he4) = xnu_burn(i,n_he4) + dummy*fac2
               xnu_burn(i,n_c12) = xnu_burn(i,n_c12) + dummy*fac3
            endif          
         end do


! burn 24Mg to 4He and 20Ne
                 
         call get_rate2(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                        sig(1:n_burn),25,26)           
         
         do i=1,n_burn
            if (xnu_burn(i,n_mg24) .gt. 0.0_rk) then
               fac2               = pc_nuc(n_he4,2)/pc_nuc(n_mg24,2)
               fac3               = pc_nuc(n_ne20,2)/pc_nuc(n_mg24,2)
               
               dummy3             = xnu_burn(i,n_he4)
               dummy4             = xnu_burn(i,n_mg24)
               dummy5             = xnu_burn(i,n_ne20)
               
               xnu_burn(i,n_mg24) = xnu_burn(i,n_mg24)*exp(-sig(i)*hydro%dt)
               dummy              = dummy4 - xnu_burn(i,n_mg24)
               
               
               xnu_burn(i,n_he4)  = xnu_burn(i,n_he4) + dummy*fac2
               xnu_burn(i,n_ne20) = xnu_burn(i,n_ne20) + dummy*fac3
            endif          
         end do


! burn 28Si to 4He and 24Mg

         call get_rate2(tem_burn(1:n_burn),rhonu(1:n_burn), &
                        sig(1:n_burn),27,28)           
         do i=1,n_burn
            if (xnu_burn(i,n_si28) .gt. 0.0) then
               fac2               = pc_nuc(n_he4,2)/pc_nuc(n_si28,2)
               fac3               = pc_nuc(n_mg24,2)/pc_nuc(n_si28,2)
               
               dummy3             = xnu_burn(i,n_he4)
               dummy4             = xnu_burn(i,n_si28)
               dummy5             = xnu_burn(i,n_mg24)
               
               xnu_burn(i,n_si28) = xnu_burn(i,n_si28)*exp(-sig(i)*hydro%dt)
               dummy              = dummy4 - xnu_burn(i,n_si28)
               

               xnu_burn(i,n_he4)  = xnu_burn(i,n_he4) + dummy*fac2
               xnu_burn(i,n_mg24) = xnu_burn(i,n_mg24) + dummy*fac3
            endif          
         end do
#endif
         
!******************************************************************

         if (config%use_flash_si) then
! burn 28Si to 56Ni

            call get_rate2(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                           sig(1:n_burn),27,28)
         
            do i=1,n_burn
               if (tem_burn(i).gt.  3.0e+9_rk.and. &
                   xnu_burn(i,n_si28).gt.0.0_rk) then
               
                  dummy4             = xnu_burn(i,n_si28)
                  dummy5             = xnu_burn(i,n_ni56)
                  
                  xnu_burn(i,n_si28) = xnu_burn(i,n_si28)*exp(-sig(i)*hydro%dt)
                  dummy              = dummy4 - xnu_burn(i,n_si28)
                  
                  xnu_burn(i,n_ni56) = xnu_burn(i,n_ni56) + dummy
               endif
            end do
         endif
         
         n_burn=0
         do i=5,nzn4
            if (lburn(i)) then
               n_burn        = n_burn+1
#ifdef NEC_COMPILER
!CDIR EXPAND=25            
#endif
               xnnu(i,1:config%qn) = xnu_burn(n_burn,1:config%qn)
            end if
         end do
         
#ifdef TRIPLEALPHA

         lburn=.false.
         
!     gather arrays
         do i=5,nzn4
            lburn(i)= rhonu(i) .lt.   rhonse  &
                 .and. .not. soliton(i)       &
                 .and. tmp(i)   .lt. thigh    &
                 .and. tmp(i)   .gt. 2.e+8_rk   
         end do
         
         n_burn=0
#ifdef NEC_COMPILER
!CDIR NODEP
#endif
         do i=5,nzn4
            if (lburn(i)) then
               n_burn                 = n_burn+1
               rho_burn(n_burn)       = rho(i)
               tem_burn(n_burn)       = tmp(i)
               xnu_burn(n_burn,n_he4) = xnnu(i,n_he4)
               xnu_burn(n_burn,n_c12) = xnnu(i,n_c12)
            end if
         end do

         if (n_burn .gt. 0_ik) then

! burn 3*4He to 12C  
!*     Y1+Y1+Y1-->Y2
!*     dY1/dt=-rate*Y1^3
!*     --> X1f=X1i/sqrt[(1+X1i^2*rate*dt/A1^2)]
!*     X2f=X2i+1/3*dX1*A2/A1=X2i+dX1

            call get_rate5(tem_burn(1:n_burn),rho_burn(1:n_burn), &
                           sig(1:n_burn),29,30)
            
            do i=1,n_burn
               if (xnu_burn(i,n_he4) .gt. 0.0_rk) then
                  
                  dummy  = xnu_burn(i,n_he4)
                  dummy3 = xnu_burn(i,n_c12)            
                  
                  xnu_burn(i,n_he4) = 1/(1.0_rk+xnu_burn(i,n_he4)* &
                                      xnu_burn(i,n_he4)*sig(i)*hydro%dt/ &
                                      (pc_nuc(n_he4,2)*pc_nuc(n_he4,2)))
                  xnu_burn(i,n_he4) = dummy*sqrt(xnu_burn(i,n_he4))
                  xnu_burn(i,n_c12) = xnu_burn(i,n_c12)+ &
                                      (dummy - xnu_burn(i,n_he4))  
               endif
            enddo
         end if
      
!     Write new compositions to xnnu-array:
         n_burn=0
         do i=5,nzn4
            if (lburn(i)) then
               n_burn         = n_burn+1
               xnnu(i,n_he4)  = xnu_burn(n_burn,n_he4)
               xnnu(i,n_c12 ) = xnu_burn(n_burn,n_c12 )
            end if
         end do
#endif

      end if ! (n_burn .gt. 0_ik)

#ifdef TAK_RATES 
#ifdef TAK_SUBTIME
      dtrate=hydro%dt*1.e-2_rk
      do j=1,100
#endif /* TAK_SUBTIME */

      do i=5, nzn4
         lburn(i)=rhonu(i) .lt.   rhonse                          &
             .and. .not. soliton(i)                               &
             .and. tmp(i)   .lt. thigh                            &
             .and. tmp(i)   .gt. 10.0_rk**ltaktemin               &
             .and. tmp(i)   .lt. 10.0_rk**ltaktemax               &
             .and. rhonu(i)*xnnu(i,config%qn) .lt. 10.0_rk**ltakroymax   &
             .and. rhonu(i)*xnnu(i,config%qn) .gt. 10.0_rk**ltakroymin
      end do

      n_burn=0
#ifdef NEC_COMPILER
!CDIR NODEP
#endif
      do i=5, nzn4
         if (lburn(i)) then
            n_burn                  = n_burn+1
            iburn(n_burn)           = i
            rho_burn(n_burn)        = rhonu(i)
            enu_burn(n_burn)        = enu(i)
            tem_burn(n_burn)        = tmp(i)
            xnu_burn(n_burn,n_mg24) = xnnu(i,n_mg24)
            xnu_burn(n_burn,n_ne20) = xnnu(i,n_ne20)
            xnu_burn(n_burn,config%qn)     = xnnu(i,config%qn)
         end if
      end do
      
      if (n_burn .gt. 0_ik) then
         
         call tak_vektor(rho_burn(1:n_burn), xnu_burn(1:n_burn,config%qn),  &
                tem_burn(1:n_burn), ecnes(1:n_burn), ecmgs(1:n_burn),&
                ennes(1:n_burn), enmgs(1:n_burn))

         do i=1,n_burn
            
!*     dYe=-SUMj(rate_ecapj*Xj/Aj)*dt

            xnu_burn(i,config%qn) = xnu_burn(i,config%qn) -                    &
                ( ecnes(i) * xnu_burn(i,n_ne20)/pc_nuc(n_ne20,2) &
                + ecmgs(i) * xnu_burn(i,n_mg24)/pc_nuc(n_mg24,2) ) * hydro%dt
            
!*     de=SUMj(rate_enuj)*dt [erg/g]

            enu_burn(i) = enu_burn(i) - (ennes(i)*xnu_burn(i,n_ne20) &
                                     + enmgs(i)*xnu_burn(i,n_mg24) ) &
                                     * pc_meverg * hydro%dt / pc_mb
         end do

#ifdef NEC_COMPILER
!CDIR NODEP
#endif
         do i=1,n_burn
            rhonu(iburn(i))       = rho_burn(i)
            enu(iburn(i))         = enu_burn(i)
            tmp(iburn(i))         = tem_burn(i)
            xnnu(iburn(i),n_mg24) = xnu_burn(i,n_mg24)
            xnnu(iburn(i),n_ne20) = xnu_burn(i,n_ne20)
            xnnu(iburn(i),config%qn)     = xnu_burn(i,config%qn)
         end do


#ifdef TAK_SUBTIME
         do i=5,ny+4
            i4=i - 4
            xnuc(i4,jj,kk,n_n)    = max(xnnu(i,n_n),0.0_rk)
            xnuc(i4,jj,kk,n_p)    = max(xnnu(i,n_p),0.0_rk)
            xnuc(i4,jj,kk,n_he4)  = max(xnnu(i,n_he4),0.0_rk)
            xnuc(i4,jj,kk,n_si28) = max(xnnu(i,n_si28),0.0_rk)
            xnuc(i4,jj,kk,n_ni56) = max(xnnu(i,n_ni56),0.0_rk)
            xnuc(i4,jj,kk,n_c12)  = max(xnnu(i,n_c12),0.0_rk)
            xnuc(i4,jj,kk,n_mg24) = max(xnnu(i,n_mg24),0.0_rk)
            xnuc(i4,jj,kk,n_o16)  = max(xnnu(i,n_o16),0.0_rk)
            xnuc(i4,jj,kk,n_ne20) = max(xnnu(i,n_ne20),0.0_rk)
            xnuc(i4,jj,kk,config%qn)     = max(xnnu(i,config%qn),0.0_rk)
         end do
#endif /* TAK_SUBTIME */

     
#ifdef TAK_SUBTIME
!     needs to be checked!
         eni(1:nzn)      = enu(5:nzn4)*rhonu(5:nzn4) 
         rhon(1:nzn)     = rhonu(5:nzn4)
         tem(1:nzn)      = tmp(5:nzn4)
         xnu(1:nzn,1:config%qn) = xnuc(5:nzn4,1,1,1:config%qn)
         
         err=.false.
         
         call eos(rhon, tem, xnu, xhrepqs,              &
                  zahs, eni, ps, gamcs, ss, cuqs, ceqs, &
                  cnqs, cpqs, mode=2, nsemode=0, ler=err)
 
         enddo
#endif /* TAK_SUBTIME */

      end if

#endif /* TAK_RATES */


!******************************************************************
#endif /* BURN_NETWORK && NOT BURN_NETW_NOA*/


#ifdef BURN_NETW_NOA
!     needs to be checked!

! burn 16O to 28Si
      do i = 5, nzn4
         if ( tmp(i)   .gt.  tlow &
              .and. tmp(i)   .lt.  thigh  &
              .and. ishck(i,1,1) .lt.   1 &
              .and. rhonu(i) .lt.   lsrolo   )  then

            if (xnnu (i,6) .eq. 0._rk .and. xnnu (i,8) .gt. 0._rk) then

            call get_rate(tmp(i),rhonu(i),sig,3,0,3)
           
            dummy2 = xnnu(i,4)
            dummy3 = xnnu (i,8)
 
            fac  =  pc_nuc(n_o16,2)
            fac2 = 1._rk

            xnnu (i,8) = xnnu (i,8)/(1._rk+xnnu(i,8)*(sig/fac*fac2)*hydro%dt)
            xnnu (i,4) = xnnu (i,4) + (dummy3 - xnnu (i,8))
            
            if ((xnnu(i,4)-dummy2).gt.0._rk) then
!               write(*,*)'O16 burn!','  zone:', i
!               write(*,*)' O16',xnnu (i,8),(xnnu(i,8)-dummy3)
!               write(*,*)' Si28',xnnu (i,4),(xnnu(i,4)-dummy2)
            endif
         endif
      endif
      enddo




! burn 12C and 12C to 24Mg 

            do i = 5, nzn4
               if ( tmp(i)   .gt.  tlow &
              .and. tmp(i)   .lt.  thigh  &
              .and. ishck(i,1,1) .lt.   1    &
              .and. rhonu(i) .lt.   lsrolo  )  then  

           if (xnnu (i,6).gt.0._rk .and. xnnu (i,8).eq.0._rk ) then

               call get_rate(tmp(i),rhonu(i),sig,1,0,3)

               dummy2 = xnnu(i,7)
               dummy3 = xnnu(i,6)            
 
               fac  = pc_nuc(n_c1212,2)
               fac2 = 1._rk

              xnnu (i,6) =xnnu (i,6)/(1._rk+xnnu(i,6)*(sig/fac*fac2)*hydro%dt)
               xnnu (i,7) = xnnu (i,7) + (dummy3 - xnnu (i,6))
            
               if ((xnnu(i,7)-dummy2).gt.0._rk) then
!                  write(*,*)'C12 burn!','  zone:', i
!                  write(*,*)' C12',xnnu (i,6),(xnnu(i,6)-dummy3)
!                  write(*,*)' Mg24',xnnu (i,7),(xnnu(i,7)-dummy2)
               endif
              endif

! burn 12C and 16O to 28Si

            if (xnnu (i,6).gt.0._rk .and. xnnu (i,8).gt.0._rk ) then

               call get_rate(tmp(i),rhonu(i),sig,2,0,3)
               fac  = 0.5_rk*(pc_nuc(n_c1212,2)+pc_nuc(n_o16,2))
               fac2 = 1._rk
               
!              Xc+o|i 
               dummy  = xnnu(i,6) + xnnu (i,8)
               dummy2 = dummy
               
!              Xc+o|f                
               dummy  = dummy/(1._rk+ dummy*(sig/fac*fac2)*hydro%dt)

               dummy4 = xnnu(i,6)
               dummy5 = xnnu(i,8)
               dummy6 = xnnu(i,4)


               xnnu(i,6) = xnnu(i,6)/dummy2*dummy
              
               xnnu(i,8) = xnnu(i,8)/dummy2*dummy

               
               xnnu (i,4) = xnnu (i,4) + (dummy2 - dummy)

            if ((xnnu(i,4)-dummy6).gt.0._rk) then
!               write(*,*)'C12+O16 burn!','  zone:', i
!               write(*,*)' C12',xnnu (i,6),(xnnu(i,6)-dummy4)
!               write(*,*)' O16',xnnu (i,8),(xnnu(i,8)-dummy5)
!               write(*,*)' Si28',xnnu(i,4),(xnnu(i,4)-dummy6)
            endif
            endif 
         endif
      enddo

! burn 20Ne to 24Mg

      do i = 5, nzn4
         if ( tmp(i)   .gt.  tlow &
             .and. tmp(i)   .lt.  thigh  &
             .and. ishck(i,1,1) .lt.   1    &
             .and. rhonu(i) .lt.   lsrolo  )  then  

            if (xnnu (i,9) .gt.0._rk) then

            call get_rate(tmp(i),rhonu(i),sig,6,7,4)           
            
            
            dummy5 = xnnu(i,7)
            dummy3 = xnnu (i,9)            

!              Xhe+ne|i 
               dummy  = xnnu (i,9)
               dummy2 = dummy
               
!              Xhe+ne|f                
               dummy  = dummy/(1._rk+ dummy*(sig/fac*fac2)*hydro%dt)


               xnnu(i,9) = xnnu(i,9) - (dummy2-dummy)
   
               xnnu (i,7) = xnnu (i,7) + (dummy2 - dummy)
                        
            if ((xnnu(i,7)-dummy2).gt.0._rk) then
!               write(*,*)'Ne20 burning!','  zone:', i
!               write(*,*)' Ne20',xnnu(i,9),(xnnu(i,9)-dummy3)
!               write(*,*)' Mg24',xnnu (i,7),(xnnu(i,7)-dummy5)
            endif
            endif
            endif
            enddo
             
! burn 24Mg to 28Si

      do i = 5, nzn4
         if ( tmp(i)   .gt.  tlow &
             .and. tmp(i)   .lt.  thigh &
             .and. ishck(i,1,1) .lt.   1   &
             .and. rhonu(i) .lt.   lsrolo  )  then  

            if (xnnu (i,7) .gt.0._rk) then

            call get_rate(tmp(i),rhonu(i),sig,4,5,4)           
            
           dummy5 = xnnu(i,7)
            dummy3 = xnnu (i,4)            

!              Xhe+mg|i 
               dummy  = xnnu (i,7)
               dummy2 = dummy
               
!              Xhe+mg|f                
               dummy  = dummy/(1._rk+ dummy*(sig/fac*fac2)*hydro%dt)


               xnnu(i,7) = xnnu(i,7) - (dummy2-dummy)
   
               xnnu (i,4) = xnnu (i,4) + (dummy2 - dummy)


            if ((xnnu(i,4)-dummy3).gt.0._rk) then
!               write(*,*)'Mg24 burning!','  zone:', i
!               write(*,*)' Mg24',xnnu(i,7),(xnnu(i,7)-dummy5)
!               write(*,*)' Si28',xnnu (i,4),(xnnu(i,4)-dummy3)
            endif
            endif
          
         endif
      enddo

! burn 28Si to 56Ni

      do i = 5, nzn4
         if (tmp(i).gt.4.5e9_rk.and.rhonu(i).lt.lsrolo) then
            dummy3 = xnnu(i,5)
            xnnu (i,5)=xnnu (i,5)+xnnu (i,4)
            xnnu (i,4)=0._rk
            if ((xnnu(i,5)-dummy3).gt.0._rk) then
!             write(*,*) 'Si28 burn!','  zone:',i
!             write(*,*)' Ni56',xnnu(i,5),(xnnu(i,5)-dummy2)
            endif
         endif
      enddo
#endif /*  BURN_NETW_NOA */


#if !defined(BURN_NETWORK) && !defined(BURN_NETW_NOA)

      if (config%use_flash_c) then
!******************************************************************
!     OLD BURNING
!******************************************************************
! burn 12C (and intermediate elemensts) to 24Mg
#define NUCLEAR_ELEMENT_TO_BURN_TO n_mg24
#define NUCLEAR_ELEMENT_TO_BURN_FROM n_c12
#define FIRST_IRON_ISOTOPE n_fe56
#define THRESHOLD_TEM 2.5e9_rk
#define NO_INTERMEDIATE
#undef  DEBUG_FLASH
#include "flashing.X90"
      endif

      if (config%use_flash_o) then
! burn 16O,20Ne,24Mg to 28Si
#define NUCLEAR_ELEMENT_TO_BURN_TO n_si28
#define NUCLEAR_ELEMENT_TO_BURN_FROM n_o16
#define FIRST_IRON_ISOTOPE n_fe56
#define THRESHOLD_TEM 3.5e9_rk
#undef  DEBUG_FLASH
#include "flashing.X90"
      endif

    
      if (config%use_flash_si) then  
! burn 28Si to 56Ni
#define NUCLEAR_ELEMENT_TO_BURN_TO n_ni56
#define NUCLEAR_ELEMENT_TO_BURN_FROM n_si28
#define FIRST_IRON_ISOTOPE n_fe56
#define THRESHOLD_TEM 4.5e9_rk
#define NO_INTERMEDIATE
#undef  DEBUG_FLASH
#include "flashing.X90"
      endif

#ifdef BURN_SS
!*******************************************************************
! supersonic burning: burn 16O,20Ne,24Mg,12C to 56Ni in the supersonic region
!*******************************************************************
      do i = 5, nzn4
       vels(i) = sqrt (gamc(i) * p(i) / rhonu(i))
         if (unu(i) .gt. 1.5_rk*vels(i) .and. rhonu(i) .lt. lsrolo) then
            dummy4 = xnnu(i,5)           
            xnnu (i,5)=xnnu (i,5)+xnnu (i,4)+xnnu (i,6) &
                      +xnnu (i,7)+xnnu (i,8)+xnnu (i,9)
            xnnu (i,4) =0._rk
            xnnu (i,6) =0._rk
            xnnu (i,7) = 0._rk
            xnnu (i,8) = 0._rk
            xnnu (i,9) = 0._rk
!            write(*,*) vels(i)
            if ((xnnu(i,10)-dummy4).gt.0._rk) then
!            write(*,*) i,xnnu(i,10),(xnnu(i,10)-dummy)
            endif
         endif
      enddo
!******************************************************************
#endif

#ifdef BURN_SS
!******************************************************************
! supersonic burning: burn 16O,20Ne,24Mg,12C to 56Ni in the supersonic region
!******************************************************************
      do i = 5, nzn4
         vels(i) = sqrt (gamc(i) * p(i) / rhonu(i))
         if (unu(i) .gt. 1.5_rk*vels(i) .and. rhonu(i) .lt. lsrolo) then
            dummy4 = xnnu(i,5)           
            xnnu (i,5)=xnnu (i,5)+xnnu (i,4)+xnnu (i,6) &
                      +xnnu (i,7)+xnnu (i,8)+xnnu (i,9)
            xnnu (i,4) = 0._rk
            xnnu (i,6) = 0._rk
            xnnu (i,7) = 0._rk
            xnnu (i,8) = 0._rk
            xnnu (i,9) = 0._rk
!            write(*,*) vels(i)
            if ((xnnu(i,10)-dummy4).gt.0._rk) then
!            write(*,*) i,xnnu(i,10),(xnnu(i,10)-dummy)
            endif
         endif
      enddo
!******************************************************************
#endif

      if (config%low_den_nse_eos .eq. 1) then
         if (config%restmass_version .gt. 0) then
            raise_abort("burn(): old_eos and no_rm do not work together!")
         endif

         ! recombine free nucleons and alphas to 56Ni
         do i = 5, nzn4
            t9=tmp(i)*1.e-9_rk
            rhoc=10._rk**(11.62_rk+1.5_rk*log10(t9)-39.17_rk/t9)
            if (rhonu(i) .lt. lsrolo .and. &
                 t9 .gt. 3._rk .and. rhonu(i) .gt. rhoc) then 
               xn_old         = xnnu(i,n_n)
               xp_old         = xnnu(i,n_p)
               
               xnnu(i,n_ni56)= xnnu(i,n_ni56)+xnnu(i,n_he4)+ &
                    2._rk*min(xn_old,xp_old)
               xnnu(i,n_he4) = 0._rk
               xnnu(i,n_n)   = max(xn_old-xp_old,0._rk)
               xnnu(i,n_p)   = max(xp_old-xn_old,0._rk)
               !            write(*,'(4e12.5)') xnnu(i,1),xnnu(i,2),xnnu(i,3),xnnu(i,5)
            endif
         enddo
      endif 

#endif /* NOT (BURN_NETWORK OR BURN_NETW_NOA) */


!-------  update composition

      do i = 1, nzn
         i4 = i + 4
#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK)
         energy(i,jj,kk)      = enu  (i4) ! Paco
         xnuc(i,jj,kk,config%qn)     = xnnu(i4,config%qn) ! Paco
#endif
         xnuc(i,jj,kk,n_n)    = xnnu(i4,n_n)
         xnuc(i,jj,kk,n_p)    = xnnu(i4,n_p)
         xnuc(i,jj,kk,n_he4)  = xnnu(i4,n_he4)
         xnuc(i,jj,kk,n_si28) = xnnu(i4,n_si28)
         xnuc(i,jj,kk,n_ni56) = xnnu(i4,n_ni56)
         xnuc(i,jj,kk,n_c12)  = xnnu(i4,n_c12)
         xnuc(i,jj,kk,n_mg24) = xnnu(i4,n_mg24)
         xnuc(i,jj,kk,n_o16)  = xnnu(i4,n_o16)
         xnuc(i,jj,kk,n_ne20) = xnnu(i4,n_ne20)
      end do

      if (config%restmass_version .eq. 1) then
         do i=1,config%qx
!#if LOW_EOS_NSE_FLAG != 22
            offs = moffs + SUM(xnuc(i,jj,kk,1:n_he4)*mbar(1:n_he4)) &
                      + (1._rk-SUM(xnuc(i,jj,kk,1:n_he4)))*mbar(n_rep)
!#else
!         offs = moffs + SUM(xnuc(i,jj,kk,1:n_he4)*mbar(1:n_he4))
!     &                - SUM(xnuc(i,jj,kk,n_d:n_he3)*mbar(n_d:n_he3))
!     &                + (1._rk-SUM(xnuc(i,jj,kk,1:n_he4)))*mbar(n_rep)
!     &                - (1._rk-SUM(xnuc(i,jj,kk,n_d:n_he3)))*mbar(n_rep)
!#endif
            escr(i) = escr(i) - offs
         enddo
         energy(:,jj,kk)=escr(:)
      endif


      if (config%restmass_version .eq. 2) then
         do i=1,config%qx
            offs = moffs + SUM(xnuc(i,jj,kk,1:config%qn-1)*mbar(1:config%qn-1))
            escr(i) = escr(i) - offs
         enddo
         energy(:,jj,kk)=escr(:)
      endif

      if (config%restmass_version .eq. 3) then
         do i=1,config%qx
            offs = moffs + SUM(xnuc(i,jj,kk,1:config%qn)*mbar(1:config%qn))
            escr(i) = escr(i) - offs
         enddo
         energy(:,jj,kk)=escr(:)
      endif


#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return
      end subroutine burn

      end module burn_mod

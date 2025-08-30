!>
!> This module provides the variables and functions which are needed
!> for calculating the approximative NSE treatment in the Boltzmann-gas
!> EoS.
!>
!>  Author: A. Marek, MPA, Jul. 2009
!>
module approximative_low_density_nse

  use precision

  implicit none
  private

  public approximative_nse

contains

#undef DEBG
function alpha_dissociation(t9) result(rho_critical)
  use precision
  use machcons, only : big10exp
  use specfun, only : fastlog10
  implicit none

  real(kind=rk), intent(in) :: t9(:)
  real(kind=rk), dimension(size(t9)) :: rho_critical, fit_argument

  fit_argument(:) = t9(:)

  call fastlog10(fit_argument, size(t9))

  fit_argument(:) = (11.62_rk+1.5_rk*fit_argument(:)-39.17_rk/t9(:))

  where(fit_argument(:) .gt. big10exp)

     fit_argument(:) = big10exp

  endwhere

  where(fit_argument(:) .lt. -big10exp)

     fit_argument(:) = -big10exp

  endwhere


  rho_critical(:) = 10._rk**(fit_argument(:))

end function alpha_dissociation

function nuclei_dissociation(t9) result(rho_critical)
  use precision
  use machcons, only : big10exp
  use specfun, only : fastlog10

  implicit none

  real(kind=rk), intent(in) :: t9(:)
  real(kind=rk), dimension(size(t9)) :: rho_critical, fit_argument

  fit_argument(:) = t9(:)

  call fastlog10(fit_argument, size(t9))

  fit_argument(:) = (10.60_rk+1.5_rk*fit_argument(:)-47.54_rk/t9(:))

  where(fit_argument(:) .gt. big10exp)

     fit_argument(:) = big10exp

  endwhere

  where(fit_argument(:) .lt. -big10exp)

     fit_argument(:) = -big10exp

  endwhere

  rho_critical(:) = 10._rk**fit_argument(:)

end function nuclei_dissociation


!> Driver for Calculating the nuclear composition by using an approximative
!> description of NSE as proposed in Rampp & Janka 2003
!>
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param temp temperatur [K]
!> \param dens density [g/ccm]
!> \param yel  Y_e
!> \param nnuc array of nuclear number densities
!> \param mnuc array of nuclear masses
!> \param etae degeneracy parameter electrons
!> \param etnuc degeneracy parameters nuclei
!> \param pres  pressure [erg/ccm]
!> \param entp  entropy
!> \param gamm  addiabatic index
!> \param mode  computational mode
!> \param ler   error flag
!>
  subroutine approximative_nse(temp,dens,yel,nnuc,mnuc,    &
                          etae,etnuc,eden,pres,entp,gamm,mode,ler)
    use precision
    use abort
    use phycon
    use nucparam
    implicit none

    integer(kind=ik), intent(in) :: mode
    real(kind=rk),    intent(in) :: dens(:),yel(:),mnuc(:)
    real(kind=rk),    intent(out):: etnuc(:,:),etae(:),pres(:),entp(:),        &
                             gamm(:)
    real(kind=rk),  intent(inout):: temp(:),eden(:),nnuc(:,:)
    logical, intent(out):: ler


    real(kind=rk), dimension(size(dens)) :: rhoc_alph,rhoc_nucl,temp_save,     &
                                            fac,facl,facr,t9,enbyrst,dtemp
    integer(kind=ik), dimension(size(dens)) :: info,nse_regime

    logical, dimension(size(dens)) :: not_done,no_dissoc,recomb

    real(kind=rk), dimension(size(dens),size(nnuc,dim=2)) ::   &
                   etnuc_wrk, nnuc_save, xi

    integer(kind=ik) k,kt,iwrk,nn,nuc,nseit
    real(kind=rk), parameter :: accur=1e-4_rk,taccur=1e-7_rk
    integer(kind=ik), parameter :: maxit=40

       kt=size(dens)
       nuc=size(nnuc,dim=2)

!         do nn=1,nuc
!            xi(:,nn)=(pc_nuc(nn,2)*pc_mb)*nnuc(:,nn)/dens(:)
!         enddo
!         do k=1,kt
!            if (abs(SUM(xi(k,1:nuc))-1.0).gt.1e-10) then
!               write(*,*) 'thj_nse 0:',k,kt,nuc,mode,
!     &              abs(SUM(xi(k,1:nuc))-1.0)
!               write(80,*) k,dens(k),xi(k,1:nuc)
!               stop
!            endif
!         enddo


       fac(:) = 0._rk


       select case(mode)
       case(1,4)
          ! do not try to use nse if temperature is
          ! smaller than 3.0e9_rk
          not_done(:)=temp(:)*pc_mevk .ge. 3.0e9_rk

          nse_regime(:)=0
          call eosnse(nse_regime,fac,not_done,temp,dens,yel,nnuc,mnuc,   &
               etae,etnuc,eden,pres,entp,gamm,mode,info,ler)
          ler=ANY(info(:).ne.0)
          if (ler) then
             write(*,*) 'BURN:'
             write(*,*) info
             write(*,*) not_done
          endif
          return
       case(2:3)
! test where we might enter the regime II at all by assuming a
!  fully recombined composition and exclude points accordingly
!  from composition iteration
          nnuc_save(:,:)=nnuc(:,:)
          temp_save(:)=temp(:)
          nse_regime(:)=1
          fac(:)=0.0_rk
          recomb(:)=temp(:)*pc_mevk.ge.3.0e9_rk
          call eosnse(nse_regime,fac,recomb(:),temp,dens,yel,nnuc,mnuc, &
               etae,etnuc,eden,pres,entp,gamm,mode,info,ler)
          nnuc(:,:)=nnuc_save(:,:)


          t9(:)=temp(:)*pc_mevk*1e-9_rk
          rhoc_alph(:)=alpha_dissociation(t9(:))

          rhoc_nucl(:)=nuclei_dissociation(t9(:))
          no_dissoc(:)=info(:).ne.0 .or. dens(:) > rhoc_alph(:)


          temp(:)=temp_save(:)

! call the eos for points which require no dissociation
!   (recombination is still possible !)
          nse_regime(:)=0
          call eosnse(nse_regime,fac,no_dissoc,temp,dens,yel,nnuc,      &
               mnuc,etae,etnuc,eden,pres,entp,gamm,mode,info,ler)

       end select

       enbyrst(:)=0.0_rk
       do nn=1,nuc
          enbyrst(:)=enbyrst(:)+nnuc(:,nn)*mnuc(nn)
       enddo
       enbyrst(:) = ( enbyrst(:) +      &
                    dens(:)/pc_mb*(-wc_mn+enullm+yel(:)*wc_me) ) *      &
                                                         pc_meverg


      t9(:)=temp(:)*pc_mevk*1e-9_rk
      rhoc_alph(:)=alpha_dissociation(t9(:))

#ifdef DEBG
      do nn=1,nuc
         xi(:,nn)=pc_nuc(nn,2)*nnuc(:,nn)/(dens(:)/pc_mb)
      enddo
      write(33,101) 'BURN 0:',dens(:).lt.rhoc_alph(:),not_done(:),      &
              no_dissoc(:),temp(:)*pc_mevk*1e-9,sum(xi(:,1:2),dim=2),   &
               xi(:,3), &
              sum(xi(:,4:nuc-1),dim=2),fac(:),info(:)
      write(33,101) ' '
#endif

!--------------------------------
! if possible (enough internal energy)
!    dissociate heavies into alphas

      nse_regime(:)=2
      nnuc_save(:,:)=nnuc(:,:)
      facl(:)=0.0_rk ; facr(:)=1.0_rk ; fac(:)=0.5_rk
!      fac(:)=max(min(1.0-rhoc_alph(:)/dens(:),1.0-1e-2),1e-2)


      not_done(:)=.not.no_dissoc(:)

#ifdef DEBG
      write(*,*) 'thj_nse BURN 2: ',kt,COUNT(not_done(:))
#endif

      nseit=0
      do while (COUNT(not_done).gt.0 .and. nseit.lt.maxit)
         dtemp(:)=temp(:)

         do nn=1,nuc
            where(not_done(:))
               nnuc(:,nn)=nnuc_save(:,nn)
            endwhere
         enddo


         do nn=1,nuc
            xi(:,nn)=pc_nuc(nn,2)*nnuc(:,nn)/(dens(:)/pc_mb)
         enddo
!         do k=1,kt
!            if (abs(SUM(xi(k,1:nuc))-1.0).gt.1e-10) then
!               write(*,*) 'thj_nse 1:',k,nseit,
!     &              abs(SUM(xi(k,1:nuc))-1.0)
!               stop
!            endif
!         enddo

         call eosnse(nse_regime,fac,not_done,temp,dens,yel,nnuc,mnuc,   &
                          etae,etnuc,eden,pres,entp,gamm,mode,info,ler)

         t9(:)=temp(:)*pc_mevk*1e-9_rk
         rhoc_alph(:)=alpha_dissociation(t9(:))



#ifdef DEBG
         do nn=1,nuc
               xi(:,nn)=pc_nuc(nn,2)*nnuc(:,nn)/(dens(:)/pc_mb)
         enddo
         write(33,101) 'BURN 2:',dens(:).lt.rhoc_alph(:),not_done(:),   &
              no_dissoc(:),     &
              temp(:)*pc_mevk*1e-9_rk,sum(xi(:,n_n:n_p),dim=2),xi(:,n_he4),        &
              sum(xi(:,n_he4+1:nuc-1),dim=2),fac(:),info(:)
#endif

         where(info(:) == -1 .and. not_done(:))
            facr(:)=fac(:)
!            temp(:)=temp(:)*0.9 ! only guess-value is changed
         endwhere

         where(dens(:).lt.rhoc_alph(:).and.not_done(:).and.info(:).ge.0)
            facl(:)=fac(:)
         endwhere
         where(dens(:).ge.rhoc_alph(:).and.not_done(:).and.info(:).ge.0)
            facr(:)=fac(:)
         endwhere

         where(facr(:)-facl(:) < accur .or. dtemp(:) < taccur)
            not_done(:)=.false.
         endwhere
         dtemp(:)=abs(temp(:)-dtemp(:))/temp(:)

         fac(:)=0.5_rk*(facl(:)+facr(:))


         nseit=nseit+1
      enddo
#ifdef DEBG
      write(*,*) 'thj_nse BURN 2: ',nseit-1
#endif

!--------------------------------
! if possible (temperature low enough)
!    or necessary (no temperature found in nse_regime=2)
!    recombine alphas (and nucleons) to 56Ni

      nse_regime(:)=1
      nnuc_save(:,:)=nnuc(:,:)
      facl(:)=0.0_rk ; facr(:)=1.0_rk ; fac(:)=0.5_rk

      t9(:)=temp(:)*pc_mevk*1e-9_rk
      rhoc_alph(:)=alpha_dissociation(t9(:))

      not_done(:)=(not_done(:) .or. dens(:) > rhoc_alph(:) .or. &
                  eden(:).le.enbyrst(:)) .and. recomb(:)

#ifdef DEBG
      write(*,*) 'thj_nse BURN 1: ',kt,COUNT(not_done(:))
#endif

      nseit=0
      do while (COUNT(not_done).gt.0 .and. nseit.lt.maxit)
         dtemp(:)=temp(:)
         do nn=1,nuc
            do k=1,kt
               if (not_done(k)) nnuc(k,nn)=nnuc_save(k,nn)
            enddo
         enddo

         call eosnse(nse_regime,fac,not_done,temp,dens,yel,nnuc,mnuc,   &
                          etae,etnuc,eden,pres,entp,gamm,mode,info,ler)


         t9(:)=temp(:)*pc_mevk*1e-9_rk
         rhoc_nucl(:)=nuclei_dissociation(t9(:))
         rhoc_alph(:)=alpha_dissociation(t9(:))

#ifdef DEBG
         do nn=1,nuc
               xi(:,nn)=pc_nuc(nn,2)*nnuc(:,nn)/(dens(:)/pc_mb)
         enddo
         write(33,101) 'BURN 1:',dens(:).lt.rhoc_alph(:),not_done(:),   &
              no_dissoc(:),     &
              temp(:)*pc_mevk*1e-9_rk ,sum(xi(:,n_n:n_p),dim=2),xi(:,n_he4),        &
              sum(xi(:,n_he4+1:nuc-1),dim=2),fac,info(:)
#endif

         where(info(:) == -1 .and. not_done(:))
            facl(:)=fac(:)
!            temp(:)=temp(:)*0.9
         endwhere

         where(dens(:).lt.rhoc_alph(:).and.not_done(:).and.info(:).ge.0)
            facl(:)=fac(:)
         endwhere
         where(dens(:).ge.rhoc_alph(:).and.not_done(:).and.info(:).ge.0)
            facr(:)=fac(:)
         endwhere

         dtemp(:)=abs(temp(:)-dtemp(:))/temp(:)
         where(facr(:)-facl(:) < accur .or. dtemp(:) < taccur)
            not_done(:)=.false.
         endwhere
         fac(:)=0.5_rk*(facl(:)+facr(:))


         nseit=nseit+1

      enddo
#ifdef DEBG
      write(*,*) 'thj_nse BURN 1: ',nseit-1
#endif


!--------------------------------
! if possible (enough internal energy)
!    dissociate alphas into nucleons

      nse_regime(:)=3
      nnuc_save(:,:)=nnuc(:,:)

      t9(:)=temp(:)*pc_mevk*1e-9_rk
      rhoc_nucl(:)=nuclei_dissociation(t9(:))

      not_done(:) = dens(:) < rhoc_nucl(:) .and. .not.no_dissoc
      where(not_done(:))
         facl(:)=0.0_rk ; facr(:)=1.0_rk
         fac(:)=max(min(1.0_rk-dens(:)/rhoc_nucl(:),1.0_rk-1e-2_rk),1e-2_rk)
      endwhere

#ifdef DEBG
      write(*,*) 'thj_nse BURN 3: ',kt,COUNT(not_done(:))
#endif

      nseit=0
      do while (COUNT(not_done).gt.0 .and. nseit.lt.maxit)
         dtemp(:)=temp(:)
         do nn=1,nuc
            do k=1,kt
               if (not_done(k)) nnuc(k,nn)=nnuc_save(k,nn)
            enddo
         enddo


         call eosnse(nse_regime,fac,not_done,temp,dens,yel,nnuc,mnuc,   &
                          etae,etnuc,eden,pres,entp,gamm,mode,info,ler)

         t9(:)=temp(:)*pc_mevk*1e-9_rk
         rhoc_nucl(:)=nuclei_dissociation(t9(:))

#ifdef DEBG
         do nn=1,nuc
               xi(:,nn)=pc_nuc(nn,2)*nnuc(:,nn)/(dens(:)/pc_mb)
         enddo
         write(33,101) 'BURN 3:',dens(:).lt.rhoc_nucl(:),not_done(:),   &
              no_dissoc(:),     &
              temp(:)*pc_mevk*1e-9_rk,sum(xi(:,n_n:n_p),dim=2),xi(:,n_he4),        &
              sum(xi(:,n_he4+1:nuc-1),dim=2),fac(:),info(:)
#endif

         where(info(:) == -1 .and. not_done(:))
            facr(:)=fac(:)
         endwhere

         where(dens(:).lt.rhoc_nucl(:).and.not_done(:).and.info(:).ge.0)
            facl(:)=fac(:)
         endwhere
         where(dens(:).ge.rhoc_nucl(:).and.not_done(:).and.info(:).ge.0)
            facr(:)=fac(:)
         endwhere

         dtemp(:)=abs(temp(:)-dtemp(:))/temp(:)
         where(facr(:)-facl(:) < accur .or. dtemp(:) < taccur)
            not_done(:)=.false.
         endwhere

         fac(:)=0.5_rk*(facl(:)+facr(:))

         nseit=nseit+1

      enddo
#ifdef DEBG
      write(*,*) 'thj_nse BURN 3: ',nseit-1
#endif

      ler=ANY(info(:).ne.0).or.COUNT(not_done(:)).gt.0
      if (ler) then
         write(*,*) 'BURN:'
         write(*,*) info
         write(*,*) not_done
         return
      endif
#ifdef DEBG
  101 format(a7,6l2,2(1pe11.4),8(1pe11.4),2I3)
#endif


end subroutine approximative_nse


!> Calculate the nuclear composition by using an approximatice
!> description of NSE as proposed in Rampp & Janka 2003
!>
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param nse_regime
!> \param fac
!> \param not_done array of zones that are not completed
!> \param temp temperatur [K]
!> \param dens density [g/ccm]
!> \param yel  Y_e
!> \param nnuc array of nuclear number densities
!> \param mnuc array of nuclear masses
!> \param etae degeneracy parameter electrons
!> \param etnuc degeneracy parameters nuclei
!> \param pres  pressure [erg/ccm]
!> \param entp  entropy
!> \param gamm  addiabatic index
!> \param mode  computational mode
!> \param info  array of info values
!> \param ler   error flag
!>
subroutine eosnse(nse_regime,fac,not_done,temp,dens,yel,nnuc,mnuc, &
                  etae,etnuc,eden,pres,entp,gamm,mode,info,ler)
  use precision
  use abort
  use phycon
  use nucparam
  use boltzman_gas_eos, only : boltzmann_gas

  implicit none

  ! arguments
  integer(kind=ik), intent(in) :: mode,nse_regime(:)
  real(kind=rk),    intent(in) :: dens(:),yel(:),mnuc(:),fac(:)
  real(kind=rk),    intent(out):: etnuc(:,:),etae(:),pres(:),entp(:), &
                                  gamm(:)
  real(kind=rk),  intent(inout):: temp(:),eden(:),nnuc(:,:)
  logical, intent(out):: ler
  logical, intent(in ):: not_done(:)
  integer(kind=ik), intent(out):: info(:)

  ! locals
  real(kind=rk), dimension(size(dens)) :: rh_wrk, fac_wrk
  integer(kind=ik), dimension(size(dens)) :: nse_wrk
  real(kind=rk), dimension(size(dens),size(nnuc,dim=2)) :: xi_wrk

  integer(kind=ik) :: iwrk,k,kt,nuc

  kt=size(dens)
  nuc=size(nnuc,dim=2)

  ler = .false.

  ! initialize info-array otherwise one sets implicitly the error flag not to
  ! false
  info(:) = 0

  if (mode == 2 .or. mode == 3) then

    ! gather
    iwrk=0
    do k=1,kt
       if (not_done(k)) then
          iwrk=iwrk+1
          rh_wrk(iwrk)=dens(k)
          xi_wrk(iwrk,:)=nnuc(k,:)*pc_nuc(:,2)*pc_mb/dens(k)
          fac_wrk(iwrk)=fac(k)
          nse_wrk(iwrk)=nse_regime(k)
       endif
    enddo

    if (iwrk > 0) then

      ! compute composition & temperature
      if (mode.eq.2.or.mode.eq.3) then
        call comp(nse_wrk(1:iwrk),xi_wrk(1:iwrk,1:nuc),fac_wrk(1:iwrk))
      endif

      ! "scatter" arrays
      iwrk=0
      do k=1,kt
        if (not_done(k)) then
           iwrk=iwrk+1
           nnuc(k,:) = rh_wrk(iwrk) * max(xi_wrk(iwrk,:),1e-25_rk)/(pc_nuc(:,2)*pc_mb)
        endif
      enddo
    endif
  endif

  call boltzmann_gas(temp(1:kt), dens(1:kt), yel(1:kt),         &
                     nnuc(1:kt,1:nuc), mnuc(1:nuc), etae(1:kt), &
                     etnuc(1:kt,1:nuc), eden(1:kt), pres(1:kt), &
                     entp(1:kt), gamm(1:kt), mode, info(1:kt))

end subroutine eosnse

!> Compute the compositional changes by nuclear dissociation and or
!> recombination
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param nse_regime
!> \param xi         mass fractions
!> \fac
!>
subroutine comp(nse_regime,xi,fac)

    use precision
    use nucparam

    implicit none

    integer(kind=ik), intent(in) :: nse_regime(:)
    real(kind=rk), intent(in) :: fac(:)
    real(kind=rk), intent(inout) :: xi(:,:)

    integer(kind=ik) :: i,kt,nuc
    real(kind=rk) :: xp_old,xn_old,xp_star,xn_star,xa_star

    if (config%low_den_nse_eos .eq. 2) then
       kt=size(nse_regime)
       nuc=size(xi,dim=2)

       do i = 1,kt
          if (nse_regime(i).eq.1) then
! recombine free nucleons and alphas to 56Ni
             xn_old=xi(i,n_n)
             xp_old=xi(i,n_p)

! recombine free nucleons to alphas
             xi(i,n_he4)=xi(i,n_he4)       &
                  +2.0_rk*min(xn_old,xp_old)
             xi(i,n_n) =max(xn_old-xp_old,0.0_rk)
             xi(i,n_p) =max(xp_old-xn_old,0.0_rk)

! partly recombine alphas to 56Ni
!  0 < fac < 1 gives the degree of recombination:
!   fac=0 is complete recombination
             xi(i,n_ni56)=xi(i,n_ni56)+(1.0_rk-fac(i))*xi(i,n_he4)
             xi(i,n_he4)=fac(i)*xi(i,n_he4)
          elseif (nse_regime(i).eq.2) then
             ! dissociate heavies into alphas (and free nucleons if necessary)
             !  0 < fac < 1 gives the degree of dissociation:
             !   fac=1 is complete dissociation
             xa_star= xi(i,n_he4)+fac(i)  &
                  *2.0_rk*SUM( min( pc_nuc(n_he4+1:nuc,2)-pc_nuc(n_he4+1:nuc,1),        &
                  pc_nuc(n_he4+1:nuc,1) )/pc_nuc(n_he4+1:nuc,2)      &
                  *xi(i,n_he4+1:nuc) )
             xn_star     = xi(i,n_n)     &
                  +fac(i)*SUM( max( 0.0_rk, &
                  pc_nuc(n_he4+1:nuc,2)-2._rk*pc_nuc(n_he4+1:nuc,1))/    &
                  pc_nuc(n_he4+1:nuc,2)*xi(i,n_he4+1:nuc) )
             xp_star     = xi(i,n_p)     &
                  +fac(i)*SUM( max( 0.0_rk, &
                  -pc_nuc(n_he4+1:nuc,2)+2._rk*pc_nuc(n_he4+1:nuc,1))/   &
                  pc_nuc(n_he4+1:nuc,2)*xi(i,n_he4+1:nuc) )

             xi(i,n_he4+1:nuc)=(1.0_rk-fac(i))*xi(i,n_he4+1:nuc)
! recombine (remaining) free nucleons to alphas
             xi(i,n_he4)=xa_star+2.0_rk*min(xn_star,xp_star)
             xi(i,n_n) =max(xn_star-xp_star,0.0_rk)
             xi(i,n_p) =max(xp_star-xn_star,0.0_rk)

          elseif (nse_regime(i).eq.3) then
! dissociate heavies and alphas into free nucleons
!  0 < fac < 1 gives the degree of dissociation:
!   fac=1 is complete dissociation
             xi(i,n_n)=xi(i,n_n) &
                  +fac(i)*SUM( (pc_nuc(n_he4:nuc,2)-pc_nuc(n_he4:nuc,1))/        &
                  pc_nuc(n_he4:nuc,2)*xi(i,n_he4:nuc) )
             xi(i,n_p)=xi(i,n_p) &
                  +fac(i)*SUM(  pc_nuc(n_he4:nuc,1)/pc_nuc(n_he4:nuc,2)  &
                  *xi(i,n_he4:nuc) )
             xi(i,n_he4:nuc)=(1.0_rk-fac(i))*xi(i,n_he4:nuc)
          endif
       enddo
    endif ! low_den_nse_eos == 2

end subroutine comp

end module approximative_low_density_nse

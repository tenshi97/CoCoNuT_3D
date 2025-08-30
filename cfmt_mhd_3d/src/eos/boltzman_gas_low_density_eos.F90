#undef IBM_OPT

#ifdef ONEMG_EOS

#undef COULOMB_CORRECTIONS
#define FIX_EOS_CONVERGENCE
#define SPIKE_DETEC

#else /* ONEMG_EOS */

#define COULOMB_CORRECTIONS
#undef FIX_EOS_CONVERGENCE
#undef SPIKE_DETEC

#endif /* ONEMG_EOS */

!>
!> This module provides the physical constants which are needed in
!> the calculation of the Boltzmann-gas EoS
!>
!>  Author: A. Marek, MPA, Jul. 2009
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
  module boltzmann_gas_eos_constants
    use precision
    use phycon

    implicit none
!- physical constants used by ThJ
  real(kind=rk), parameter :: pi    = pc_pi, &
                     c     = pc_cl, &           ! [cm/s]      &
                     hc    = 2.0_rk*pc_pi*wc_hc, & ! [MeVcm]     &
                     mtoe  = pc_meverg, &     ! [erg/MeV]   &
                     ktom  = pc_kmev, &       ! [MeV/K]     &
                     rme   = wc_me, &         ! [MeV]       &
                     mu    = pc_mb, &         ! [g] &
                     arad  = 8.5616e31_rk, &     ! [1/(MeV^3 cm^3)]    &
                     rbohr = 5.29177e-9_rk       ! [cm]

  end module boltzmann_gas_eos_constants
!>
!> This module provides the variables and procedures which are needed in
!> the calculation of the baryon contribution to the Boltzmann gas EoS
!>
!>  Author: A. Marek, MPA, Jul. 2009
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
  module baryon_contrib_boltzman_gas

    use precision

    implicit none
    private

    public baryon_contribution, nuclear_statistical_weight

  contains
!>
!> Calculate the cpntribution of the baryons the EoS
!>
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param temp temperatur [K]
!> \param dens density [g/cccm]
!> \param nnuc array of nuclear number densities
!> \param mnuc array of nuclear masses
!> \param gnuc array of statistical weights
!> \param etnuc degeneracy parameters nuclei
!> \param prby  pressure [erg/ccm]
!> \param enby  energy density
!> \param sby  entropy
!> \param dedtby  dEnergy_density/dtemperature
!> \param iter  iteration steps
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
    subroutine baryon_contribution (temp,dens,nnuc,mnuc,gnuc,       &
                         etnuc,enby,prby,sby,dedtby,iter)

!=======================================================================
!                            Baryonic EOS                              =
!      for a gas of nn species of nucleons and nuclei that are         =
!      treated as ideal, non-relativistic Boltzmann gases.             =
!      (Charge neutrality and baryon number conservation have to       =
!      be ensured by calling routine where the nuclear mixture         =
!      is defined.)                                                    =
!                                                                      =
!      INPUT:                                                          =
!            temp(1:np)       vector of temperatures [in MeV]          =
!            dens(1:np)       vector of densities (= mu*nb) [in g/cc]  =
!                             (baryon number density times atomic      =
!                             unit mass)                               =
!            nnuc(1:np,1:nn)  array of baryonic/nucleonic number       =
!                             densities [1/cc]                         =
!            mnuc(1:nn)       vector of baryonic/nucleonic             =
!                             rest-mass energies [in MeV]              =
!            gnuc(1:np,1:nn)  vector of baryonic/nucleonic             =
!                             statistical weights                      =
!            iter(1:nn)       mask (0=off) for points to be computed   =
!
!      OUTPUT:                                                         =
!            etnuc(1:np,1:nn) baryonic/nucleonic degeneracy            =
!                             parameters (= chem. pot./kT,             =
!                             without rest mass contribution)          =
!            enby(1:np)       baryonic energy density [in erg/cc]      =
!                             (without rest-mass energy)               =
!            prby(1:np)       baryonic pressure [in erg/cc]            =
!            sby(1:np)        baryonic entropy density [k_b/cc]        =
!            sby(1:np)        baryonic entropy density [k_b/cc]        =
!            dedtby(1:np)     dE/dT [erg/MeV/cc]
!
!
!            np               number of grid points                    =
!            nn               number of baryonic/nucleonic species     =
!
!----------------------------------------------------------------------=
!                                                                      =
!            Author:  H.-Thomas Janka, MPA;                            =
!                          Original-Version: June 99                   =
!                     (all rights reserved)                            =
!                                                                      =
!=======================================================================
      use precision

      use boltzmann_gas_eos_constants
      use specfun, only : fastlog, fastsqrt
      implicit none

      integer(kind=ik), intent(in) :: iter(:)
      real(kind=rk),    intent(in)  :: temp(:),dens(:),nnuc(:,:), &
                                       mnuc(:),gnuc(:,:)
      real(kind=rk),  intent(out)  :: etnuc(:,:), enby(:), prby(:), &
                                      sby(:), dedtby(:)

#ifdef IBM_OPT
      real(kind=rk) :: iter_float(size(iter))
#endif

!-----------------------------------------------------------------------
      integer(kind=ik) i,j,np,nn
      real(kind=rk) :: nsum(size(dens))
      real(kind=rk) :: arg(size(nnuc,dim=1), size(nnuc,dim=2))
      real(kind=rk), parameter :: eps = 1.0e-15_rk


      np=size(dens)
      nn=size(nnuc,dim=2)

#ifdef IBM_OPT
      iter_float(:)=-abs(real(iter(:),kind=rk))


      where(iter_float(:) .ne. 0._rk)
         nsum(:)=0.0_rk
         sby(:)=0.0_rk
      endwhere
#else /* IBM_OPT */
      do i = 1, np
         if (iter(i).ne.0) then
            nsum(i) = 0.0_rk
            sby(i)  = 0.0_rk
         endif
      enddo
#endif /* IBM_OPT */

      do  j = 1, nn

         do i = 1, np
            arg(i,j) = 2.0_rk*pi*mnuc(j)*temp(i)
         enddo

         call fastsqrt(arg(:,j), np)

         do i=1, np
            arg(i,j)     = (max(nnuc(i,j),eps)/gnuc(i,j))* &
                           (hc/arg(i,j))**3
         enddo

         call fastlog(arg(:,j), np)

         do  i = 1, np
#ifdef IBM_OPT

            nsum(i)=FSEL(iter_float(i),nsum(i),nsum(i)+nnuc(i,j) )

            etnuc(i,j)=FSEL(iter_float(i),etnuc(i,j), arg(i,j))

            sby(i)=FSEL(iter_float(i),sby(i), nnuc(i,j)*(2.5_rk-etnuc(i,j)))

#else /* IBM_OPT */

            if (iter(i).ne.0) then
               nsum(i) = nsum(i) + nnuc(i,j)
               etnuc(i,j) = arg(i,j)
               sby(i)  = sby(i) + nnuc(i,j)*(2.5_rk-etnuc(i,j))
            endif
#endif /* IBM_OPT */
         enddo
      enddo


#ifdef IBM_OPT
      where(iter_float(:) .ne. 0._rk)
         prby(:) = mtoe*nsum(:)*temp(:)
         enby(:) = 1.5_rk*prby(:)
         dedtby(:) = 1.5_rk*mtoe*nsum(:)
      endwhere
#else /* IBM_OPT */
      do i = 1, np
         if (iter(i).ne.0) then
            prby(i) = mtoe*nsum(i)*temp(i)
            enby(i) = 1.5_rk*prby(i)
            dedtby(i) = 1.5_rk*mtoe*nsum(i)
         endif
      enddo
#endif /* IBM_OPT */

    end subroutine baryon_contribution


!>
!> Calculate the nuclear statistical weights
!>
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param temk temperatur [K]
!> \param gnuc array of statistical weights
!> \param dgnucdt
!> \param np
!> \param nn
!> \param iter  iteration steps
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
    subroutine nuclear_statistical_weight(temk,gnuc,dgnucdt,np,nn,iter)
      use precision

      use nucparam
      use specfun, only : fastexp
      implicit none

      integer(kind=ik), intent(in) :: nn,np,iter(np)
      real(kind=rk), intent(in) :: temk(np)
      real(kind=rk), intent(out) :: gnuc(np,nn),dgnucdt(np,nn)

      real(kind=rk) :: t9(np), scrtch(np)

      real(kind=rk) fit1,fit2,fit3,fit4
      integer(kind=ik) i,j


      t9(:)=temk(:)*1e-9_rk
      do j=1,nn
         fit1=pc_nuc(j,4)
         fit2=pc_nuc(j,5)
         fit3=pc_nuc(j,6)
         fit4=pc_nuc(j,7)

         do i=1,np
            scrtch(i)=fit2/t9(i)+fit3+fit4*t9(i)
         enddo

         call fastexp(scrtch,np)

         do i=1,np
            if (iter(i) .ne. 0) then
               gnuc(i,j)    =  fit1*(1.0_rk+scrtch(i))
               dgnucdt(i,j) =  1e-9_rk*fit1*scrtch(i)*( -fit2/t9(i)**2+fit4 )
            endif
         enddo
      enddo


    end subroutine nuclear_statistical_weight

  end module baryon_contrib_boltzman_gas
!>
!> This module provides the variables and procedures which are needed in
!> the calculation of the Boltzmann-gas EoS
!>
!>  Author: A. Marek, MPA, Jul. 2009
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
  module boltzman_gas_eos

  use precision
  use phycon

  implicit none

  private
  public boltzmann_gas



  contains



!>
!> Compute the Boltzmann-gas EoS
!>
!>  Author: A. Marek, MPA, Jul. 2009, originally M. Rampp
!>
!> \param temp temperatur [K]
!> \param dens density [g/ccm]
!> \param yel  Ye
!> \param nnuc  array of nuclear number densities
!> \param mnuc  array of nuclear masses
!> \param etae degeneracy parameter electrons
!> \param etnuc degeneracy parameters nuclei
!> \param eden energy density
!> \param press pressure
!> \param entp entropy
!> \param gamm addiabatic index
!> \param mode calculation mode
!> \param info  info array
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
    subroutine boltzmann_gas(temp,dens,yel,nnuc,mnuc,etae,etnuc,eden, &
                            pres,entp,gamm,mode,info)

!=======================================================================
!      Combination of EOS for ideal baryon gas and ideal gas of        =
!      electrons+positrons+radiation, including Coulomb lattice        =
!      corrections for pressure, energy density, entropy, and          =
!      adiabatic index using a Wigner-Seitz approximation for a        =
!      body-centered cubic Coulomb lattice of nuclei embedded          =
!      by uniformly distributed ionization electrons.                  =
!                                                                      =
!      INPUT:                                                          =
!            temp(1:np)       vector of temperatures [in MeV]          =
!            dens(1:np)       vector of densities (= mu*nb) [in g/cc]  =
!                             (baryon number density times atomic      =
!                             unit mass)                               =
!            yel(1:np)        electron fraction (n(e-)-n(e+))/nb       =
!            nnuc(1:np,1:nn)  array of baryonic/nucleonic number       =
!                             densities [1/cc]                         =
!            mnuc(1:nn)       vector of baryonic/nucleonic             =
!                             rest-mass energies [in MeV]              =
!            mode             switch
!                             mode=1 temperature => energy
!                             mode=2 energy      => temperature
!                             mode=3 (equivalent mode=2)
!                             mode=4 entropy     => temperature
!      (Calling routine has to ensure charge neutrality and baryon     =
!       number conservation.)                                          =
!                                                                      =
!      OUTPUT:                                                         =
!            etae(1:np)       electron degeneracy parameter            =
!                             (= chem. pot./kT, incl. rest mass)       =
!            etnuc(1:np,1:nn) baryonic/nucleonic degeneracy            =
!                             parameters (= chem. pot./kT,             =
!                             without rest mass contribution)          =
!            eden(1:np)       energy density [in erg/cc]               =
!             (see below for the exact definitions concerning rest-mass)
!
!                             electron-positron pairs)                 =
!            pres(1:np)       pressure [in erg/cc]                     =
!            entp(1:np)       entropy density [k_b/cc]                 =
!            gamm(1:np)       adiabatic index = dln(P)/dln(rho)|_s     =
!
!            ler              error flag
!
!   Normalization of energy density:
!
!    *****************************************************************
!    *  eden = e_rel - n_b*(m_n*c^2) + n_e*(m_e*c^2) + n_b*(8.8 MeV) *
!    *****************************************************************
!
!             where: e_rel is the relativistic energy density
!                                            (including rest mass)
!                    n_b is the barion number-density
!                    m_n is the rest-mass of the neutron
!                    m_e is the rest-mass of the electron
!
!          the value 8.8 MeV (chosen as recommended by Swesty)
!              and used in producing the L&S-tables is approx. equal
!              the binding energy (per barion)
!              of Iron, which implies (assuming m_n=m_p)
!              eden ~> 0. for T->0
!
!            np               number of grid points                    =
!            nn               number of baryonic/nucleonic species     =
!----------------------------------------------------------------------=
!                                                                      =
!            Author:  H.-Thomas Janka, MPA;  Version: June 99          =
!                     (all rights reserved)                            =
!                                                                      =
!=======================================================================
      use precision
      use abort
      use error_screen
      use boltzmann_gas_eos_constants
      use machcons
#ifndef PROGRAM_remap
      use low_density_nse
#endif



      use lepton_radiation_low_den_eos, only :lepton_radiation_eos
      use eos_data_type, only : lept_radi_table

      use baryon_contrib_boltzman_gas, only :  baryon_contribution, nuclear_statistical_weight
      use nucparam, ONLY : massn,enullm !,nrep
      use specfun, only : fastlog10
      use configure

      implicit none

      integer(kind=ik), intent(in) :: mode
      real(kind=rk),    intent(in) :: dens(:),yel(:),mnuc(:)
      real(kind=rk),    intent(out):: etnuc(:,:), etae(:), pres(:), entp(:), &
                                      gamm(:)
      real(kind=rk),  intent(inout):: temp(:),eden(:),nnuc(:,:)
      integer(kind=ik), intent(out):: info(:)

!-----------------------------------------------------------------------
      real(kind=rk), dimension(size(dens)) :: enby, prby, sby, dedtby, &
                                              temk, roye, lete, lenlp, &
                                              lprlp, lslp, gamlp,      &
                                              dlpdlr, dlpdlt, dledlt,  &
                                              dedt, neli,dnsum, eclmb, &
                                              enbyrst, logtempl, logtempr
      real(kind=rk), dimension(size(nnuc,dim=1),size(nnuc,dim=2)) :: gnuc,dgnucdt
      integer(kind=ik) ::          i, np, nn
      real(kind=rk)    :: twt, fth, ee, zbar, pclmb, enlp, prlp, dledltl
!=======================================================================
      real(kind=rk), parameter :: eps = 1.0e-15_rk

! Parameters to match mass-offset of L&S EOS, see comment in
!   subroutine eos
      real(kind=rk), parameter :: masse=rme ! ,massu=mu*c**2/mtoe forcheck

! variables for newton search
      integer(kind=ik) iter(size(dens)),itsum,iit
      real(kind=rk), dimension(size(dens)) :: fun,dx,dlslp,dtemk,templ,tempr,damp
      real(kind=rk) :: dfun,dtem
      real(kind=rk) :: xacc,racc, scratch1, scratch2

      integer(kind=ik), parameter :: newtit=10, biseit=5
      logical lerr,ler

#ifdef SPIKE_DETEC
      real(kind=rk), dimension(size(dens),10) :: temk_tmp
      real(kind=rk) :: dt
      integer(kind=ik) :: j,k
#endif
#if defined(SPIKE_DETEC) || defined (CFC_TRANSPORT)
      real(kind=rk) :: medtem
      logical, dimension(size(dens)) :: ave_tem
#endif

      np=size(dens)
      nn=size(nnuc,dim=2)


      xacc=2.0_rk*emach
      racc=2.0_rk*emach

      lerr=.FALSE.
      info(1:np)=0

      twt = 2.0_rk/3.0_rk
      fth = 4.0_rk/3.0_rk
      ee  = mtoe*wc_hc**2/(rbohr*rme)

      iter(1:np)=1              ! flag grid points

#ifndef PROGRAM_remap
      if ((config%low_den_nse_eos .gt. 2) .and. ((mode.eq.3).or.(mode.eq.2))) then
         call low_density_nse_eos(dens(1:np),yel(1:np),eden(1:np),temp(1:np),nnuc(1:np,1:nn))
      endif
#endif

!c -- normalization of energy density
      enbyrst(:)=0.0_rk
      do i=1,nn
        enbyrst(:) = enbyrst(:)+nnuc(:,i)*mnuc(i)
      enddo
      enbyrst(:) = (enbyrst(:) + dens(:) / mu * (-massn + enullm + yel(:) * masse)) * mtoe
!      enbyrst(:) = 0.0
!      stop '++++++++++ CHECK tj +++++++++++'

      select case(mode)
      case(1)

         do i = 1, np
            temk(i) = temp(i)/ktom ! MeV ----> K
            roye(i) = dens(i)*yel(i) ! lepton table needs rho*Ye = ne*mu
            neli(i) = roye(i)/mu ! number density of ioniz. electrons
            dnsum(i) = SUM(nnuc(i,1:nn))
            if  (dnsum(i) .gt. 0.0_rk)  then
               zbar = neli(i)/dnsum(i)
            else
               zbar = 0.0_rk
            end if
#ifndef COULOMB_CORRECTIONS
            eclmb(i) = 0.0_rk
#else
            eclmb(i) = -1.44423_rk*ee*(zbar**twt)*neli(i)**fth
#endif /* COULOMB_CORRECTIONS */
         enddo
!     compute nuclear partition function (statistical weights)
         call NUCLEAR_STATISTICAL_WEIGHT(temk,gnuc,dgnucdt,np,nn,iter)

!     compute baryonic EOS:
         call baryon_contribution (temp,dens,nnuc,mnuc,gnuc,etnuc,enby,prby,sby,dedtby,iter)



         if (ANY(temk(1:np).lt.lept_radi_table(1)%ttmin .or. &
                 temk(1:np).gt.lept_radi_table(1)%ttmax)) then
            write(*,*) 'lepton_radiation_table> boundaries (temk) exceeded; mode= ',mode
            where(temk(1:np).lt.lept_radi_table(1)%ttmin.or. &
                  temk(1:np).gt.lept_radi_table(1)%ttmax)
               info(1:np)=-2
            elsewhere
               info(1:np)=0
            endwhere
            return
         endif

!     compute leptonic/radiation EOS:
         call  lepton_radiation_eos(temk,roye,lete,lenlp, &
                                                lprlp,lslp,gamlp,     &
                                                dlpdlr,dlpdlt,dledlt,dedt,iter)

         do i = 1, np
            enlp  = 10.0_rk**lenlp(i)
            prlp  = 10.0_rk**lprlp(i)
            pclmb = eclmb(i)/3.0_rk
            dledltl = dledlt(i)
            if (abs(dledltl) .lt. eps) dledltl = eps

            etae(i) = 10.0_rk**lete(i)
            eden(i) = enby(i) + enlp + eclmb(i) + enbyrst(i)
            pres(i) = prby(i) + prlp + pclmb
            entp(i) = sby(i)  + 10.0_rk**lslp(i)   &
                 + fth*eclmb(i)/(mtoe*temp(i))
            gamm(i) = (prlp*dlpdlr(i) + prby(i) + fth*pclmb +   &
                 (prlp*dlpdlt(i) + prby(i))**2/ &
                 (enlp*dledltl   + enby(i))   ) / pres(i)
         enddo


      case(2:3)                ! energy given
         do i = 1, np
            roye(i) = dens(i)*yel(i) ! lepton table needs rho*Ye = ne*mu
            neli(i) = roye(i)/mu ! number density of ioniz. electrons
            dnsum(i) = SUM(nnuc(i,1:nn))
            if  (dnsum(i) .gt. 0.0_rk)  then
               zbar = neli(i)/dnsum(i)
            else
               zbar = 0.0_rk
            end if
#ifndef COULOMB_CORRECTIONS
            eclmb(i) = 0.0_rk
#else
            eclmb(i) = -1.44423_rk*ee*(zbar**twt)*neli(i)**fth
#endif /* COULOMB_CORRECTIONS */
         enddo

         scratch1 = log10(ktom*lept_radi_table(1)%ttmin)
         scratch2 = log10(ktom*lept_radi_table(1)%ttmax)

         templ(:)= scratch1
         tempr(:)= scratch2

! -- begin BISECTION  iteration-loop
         do iit=1,biseit

            do i = 1, np
               if (iter(i).ne.0) then
                  temp(i)=10._rk**(0.5_rk*(templ(i)+tempr(i)))
                  temk(i) = temp(i)/ktom ! MeV ----> K
                  if (temk(i).lt.lept_radi_table(1)%ttmin .or.   &
                      temk(i).gt.lept_radi_table(1)%ttmax) info(i)=-1
               endif
            enddo
            if (ANY(info.eq.-1)) then
               do i = 1, np
                  if (info(i).eq.-1) then
                     write(*,*) i,temk(i),eden(i),yel(i),dens(i)
                     iter(i)=0
#ifdef FIX_EOS_CONVERGENCE
                     temk(i)=lept_radi_table(1)%ttmin*1.05_rk
#endif
                  endif
               enddo
#ifdef FIX_EOS_CONVERGENCE
               info(:)=0
#endif
               call show_error_screen("eostj","BISECTION iteration loop: boundaries (temk) exceeded; mode,iit",mode,iit)
!               ler=.true.
!               info=-1
!               return
            endif
            lerr=.FALSE.


!     compute nuclear partition function (statistical weights)
            call NUCLEAR_STATISTICAL_WEIGHT(temk,gnuc,dgnucdt,np,nn,iter)

!     compute baryonic EOS:
            call baryon_contribution (temp,dens,nnuc,mnuc,gnuc,etnuc,enby,prby,sby,dedtby,iter)

!     compute leptonic/radiation EOS:
            call lepton_radiation_eos(temk,roye,lete,lenlp, &
                                                  lprlp,lslp,gamlp,  &
                                                  dlpdlr,dlpdlt,dledlt,dedt,iter)


            logtempl(:) = temp(:)
            call fastlog10(logtempl, np)
            logtempr(:) = logtempl(:)

            do i=1,np
               if (iter(i).ne.0) then
                  fun(i)=( enby(i)+10.0_rk**lenlp(i)+eclmb(i)+enbyrst(i)   &
                         -eden(i) )
                  if (fun(i) .lt. 0.0_rk) then
                    templ(i)=logtempl(i)
                  else
                    tempr(i)=logtempr(i)
                  endif

               endif
            enddo

            itsum=SUM(iter(:))
         enddo
! -- end BISECTION iteration-loop


         iter(1:np)=1
! -- begin NEWTON iteration-loop
         do iit=1,newtit

            do i = 1, np
               if (iter(i).ne.0) then
                  temk(i) = temp(i)/ktom ! MeV ----> K
                  if (temk(i).lt.lept_radi_table(1)%ttmin .or.   &
                      temk(i).gt.lept_radi_table(1)%ttmax) info(i)=-1
               endif
            enddo
            if (ANY(info.eq.-1)) then
               do i = 1, np
                  if (info(i).eq.-1) then
                     write(0,*) i,temk(i),eden(i),yel(i),dens(i)
                     iter(i)=0
#ifdef FIX_EOS_CONVERGENCE
                     temk(i)=lept_radi_table(1)%ttmin*1.05_rk
#endif
                  endif
               enddo
#ifdef FIX_EOS_CONVERGENCE
               info(:)=0
#endif
               call show_error_screen("eostj","BISECTION iteration loop: boundaries (temk) exceeded; mode,iit",mode,iit)
!               ler=.true.
!               info=-1
!               return
            endif
            lerr=.FALSE.

!     compute nuclear partition function (statistical weights)
            call NUCLEAR_STATISTICAL_WEIGHT(temk,gnuc,dgnucdt,np,nn,iter)

!     compute baryonic EOS:
            call baryon_contribution (temp,dens,nnuc,mnuc,gnuc,etnuc,enby,prby,sby,dedtby,iter)

!     compute leptonic/radiation EOS:
            call lepton_radiation_eos(temk,roye,lete,lenlp, &
                                                  lprlp,lslp,gamlp,  &
                                                  dlpdlr,dlpdlt,dledlt,dedt,iter)

            do i=1,np
               if (iter(i).ne.0) then
                  fun(i)=( enby(i)+10.0_rk**lenlp(i)+eclmb(i)+enbyrst(i)   &
                         -eden(i) )
!                  dfun=dledlt(i)*10.0_rk**lenlp(i)/temp(i)+dedtby(i)
                  dfun=dedt(i)/ktom+dedtby(i)

!                  if (dx(i)*fun(i)/dfun .lt. 0.0_rk) then
!c                     damp(i)=0.99_rk
!                     damp(i)=1.0_rk
!                  else
!                     damp(i)=1.0_rk
!                  endif

                  if (dfun .gt. 0.0_rk) then
                     dx(i) = fun(i)/dfun
                  else
                     dx(i)=0.0_rk
                     lerr=.TRUE.
                  endif

                  if (dx(i) .gt. temp(i)) then
                     dx(i)=min(dx(i),temp(i)-lept_radi_table(1)%ttmin*ktom)
                     damp(i)=0.995_rk
                  else
                     damp(i)=1.0_rk
                  endif

                  temp(i)=temp(i)-dx(i)*damp(i)
                  if (abs(dx(i)/temp(i)).lt.xacc.or.    &
                       abs(fun(i)/eden(i)).lt.racc) iter(i)=0
!                  write(*,'(2I5,5e13.5)') i,iit,
!     &                 fun(i)/eden(i),dfun,dx(i)
               endif
            enddo
            if (lerr) then
               call show_error_screen("eostj","something wrong: mode= ",mode)
               info=-99
               return
            endif

            itsum=SUM(iter(:))
            if (itsum.eq.0) EXIT
         enddo
! -- end NEWTON iteration-loop


! try a final 'bisection' if Newton has not converged
         if (itsum.ne.0) then
            write(0,*) 'tj> final bisection used'
            do i = 1, np
               if (iter(i).ne.0) then
                  temp(i) = temp(i)+0.5_rk*dx(i)
                  temk(i) = temp(i)/ktom ! MeV ----> K
                  if (temk(i).lt.lept_radi_table(1)%ttmin .or.   &
                      temk(i).gt.lept_radi_table(1)%ttmax) info(i)=-1
               endif
            enddo
            if (ANY(info.eq.-1)) then
               do i = 1, np
                  if (info(i).eq.-1) then
                     write(0,*) i,temk(i),eden(i),yel(i),dens(i)
                     iter(i)=0
#ifdef FIX_EOS_CONVERGENCE
                     temk(i)=lept_radi_table(1)%ttmin*1.05_rk
#endif
                  endif
               enddo
#ifdef FIX_EOS_CONVERGENCE
               info(:)=0
#endif
               call show_error_screen("eostj","final BISECTION loop: boundaries (temk) exceeded; mode,iit",mode,iit)
!               ler=.true.
!               info=-1
!               return
            endif

!     compute nuclear partition function (statistical weights)
            call NUCLEAR_STATISTICAL_WEIGHT(temk,gnuc,dgnucdt,np,nn,iter)

!     compute baryonic EOS:
            call baryon_contribution (temp,dens,nnuc,mnuc,gnuc,etnuc,enby,prby,sby,dedtby,iter)

!     compute leptonic/radiation EOS:
            call lepton_radiation_eos(temk,roye,lete,lenlp, &
                                                  lprlp,lslp,gamlp,  &
                                                dlpdlr,dlpdlt,dledlt,dedt,iter)

            do i=1,np
               if (iter(i).ne.0) then
                  fun(i)=( enby(i)+10.0_rk**lenlp(i)+eclmb(i)+enbyrst(i)   &
                          -eden(i) )
                  if (abs(0.5_rk*dx(i)).lt.xacc.or.        &
                       abs(fun(i)/eden(i)).lt.racc) iter(i)=0
!                  write(*,*) i,fun(i)/eden(i),0.5*dx(i),fun(i),eden(i)
               endif
            enddo
            itsum=SUM(iter(:))
            if (itsum.ne.0) then
               do i=1,np
                  if (iter(i).ne.0) then
                     write(0,*) 'THJ:',i,0.5_rk*dx(i),  &
                       abs(fun(i)/eden(i)), fun(i),eden(i)

                  endif
               enddo
              call show_error_screen("eostj","ROOT not converged: mode , Number of points",mode,itsum)
               info=-5
               return
            endif
         endif

#ifdef SPIKE_DETEC
         ave_tem(:)=.true.
!         do while(any(ave_tem(:)))
! -- temperature spikes-----
            dT = 1.25d0
            ave_tem(:)=.false.
            temk_tmp(:,:)=spread(temk(:),dim=2,ncopies=10)
            do j=-2,-1
               do k=1,2
                  do i=6,np-5
                     if ( &
     &                    (temk(i) .gt. temk_tmp(i+j,11+j)*dT &
     &                    .and. temk(i) .gt. temk_tmp(i+k,k)*dT) &
     &                    .or. (temk_tmp(i+j,11+j) .gt. temk(i)*dT &
     &                    .and. temk_tmp(i+k,k) .gt. temk(i)*dT)) then
                        ave_tem(i)=.true.
!                     if (temk_tmp(i+j,11+j).gt.temk(i)*dT &
!     &                    .and. temk_tmp(i+k,k).gt.temk(i)*dT) then
!                        ave_tem(i)=.true.
                     endif
                  enddo
               end do
            enddo

            do i = 6,np-5
               if (ave_tem(i)) then

                  medtem = (temk(i)*0+temk_tmp(i-5,1)+temk_tmp(i-4,2)+ &
     &                temk_tmp(i-3,3)+temk_tmp(i-2,4)+temk_tmp(i-1,5)+ &
     &                temk_tmp(i+1,6)+temk_tmp(i+2,7)+temk_tmp(i+3,8)+ &
     &                temk_tmp(i+4,9)+temk_tmp(i+5,10))/10.0d0
!                  if (medtem.gt.tjttmin) then
                     temk(i) = medtem
                     temp(i) = medtem*ktom
!                  endif
!                  medtem = (pres(i)*0+pres(i-5)+pres(i-4)+ &
!     &                pres(i-3)+pres(i-2)+pres(i-1)+ &
!     &                pres(i+1)+pres(i+2)+pres(i+3)+ &
!     &                pres(i+4)+pres(i+5))/10.0d0
!                  pres(i) = medtem
               end if
            end do
!         enddo
! --------------------------
#endif

#ifdef CFC_TRANSPORT
!         ave_tem(:)=temk(:).lt.1.05*lept_radi_table(1)%ttmin
!         do while(any(ave_tem(:)))
            ave_tem(:)=.false.
            do i=1,np
               if (temk(i).lt.2*lept_radi_table(1)%ttmin) ave_tem(i)=.true.
            end do
            if (any(ave_tem(:))) then
               medtem=lept_radi_table(1)%ttmin
               do i=1,np
                  if (ave_tem(i)) then
                     temk(i)=max(medtem,temk(i))
                     temp(i)=temk(i)*ktom
                  endif
                  medtem=temk(i)
               end do
               do i=np,1,-1
                  if (ave_tem(i)) then
                     temk(i)=max(0.5*medtem+temk(i),temk(i))
                     temp(i)=temk(i)*ktom
                  endif
                  medtem=temk(i)
               end do
            end if
!         end do
#endif


! -- root found
         do i = 1, np
            enlp  = 10.0_rk**lenlp(i)
            prlp  = 10.0_rk**lprlp(i)
            pclmb = eclmb(i)/3.0_rk
            dledltl = dledlt(i)
!            if  (dledltl .eq. 0.0_rk)   dledltl = eps
            if (abs(dledltl) .lt. eps) dledltl = eps

            etae(i) = 10.0_rk**lete(i)
            pres(i) = prby(i) + prlp + pclmb
            entp(i) = sby(i)  + 10.0_rk**lslp(i)   &
                 + fth*eclmb(i)/(mtoe*temp(i))
            gamm(i) = (prlp*dlpdlr(i) + prby(i) + fth*pclmb +   &
                 (prlp*dlpdlt(i) + prby(i))**2/ &
                 (enlp*dledltl   + enby(i))   ) / pres(i)
         enddo


!--
      case(4)                ! entropy given
         do i = 1, np
            roye(i) = dens(i)*yel(i) ! lepton table needs rho*Ye = ne*mu
            neli(i) = roye(i)/mu ! number density of ioniz. electrons
            dnsum(i) = SUM(nnuc(i,1:nn))
            if  (dnsum(i) .gt. 0.0_rk)  then
               zbar = neli(i)/dnsum(i)
            else
               zbar = 0.0_rk
            end if
#ifndef COULOMB_CORRECTIONS
            eclmb(i) = 0.0_rk
#else
            eclmb(i) = -1.44423_rk*ee*(zbar**twt)*neli(i)**fth
#endif /* COULOMB_CORRECTIONS */
         enddo


! -- begin iteration-loop
         do iit=1,newtit

            do i = 1, np
               if (iter(i).ne.0) then
                  temk(i) = temp(i)/ktom ! MeV ----> K
                  if (temk(i).lt.lept_radi_table(1)%ttmin .or.   &
                      temk(i).gt.lept_radi_table(1)%ttmax) lerr=.TRUE.
               endif
            enddo
            if (lerr) then
               do i = 1, np
                  if (iter(i).ne.0) then
                     write(0,*) i,temk(i),entp(i),yel(i)
#ifdef ONEMG_EOS
                     temk(i)=lept_radi_table(1)%ttmin*1.05_rk
#endif
                  endif
               enddo
#ifdef ONEMG_EOS
               lerr=.false.
#else
               call show_error_screen("eostj","iteration loop: boundaries (temk) exceeded; mode,iit",mode,iit)
               ler=.true.
               return
#endif /* ONEMG_EOS */
            endif

!     compute nuclear partition function (statistical weights)
            call NUCLEAR_STATISTICAL_WEIGHT(temk,gnuc,dgnucdt,np,nn,iter)
!     compute baryonic EOS:
            call baryon_contribution (temp,dens,nnuc,mnuc,gnuc,etnuc,enby,prby,sby,dedtby,iter)

!     compute leptonic/radiation EOS:
            dtem=1.e-7_rk
            dtemk=temk(:)*(1._rk+dtem)
            call lepton_radiation_eos(dtemk,roye,lete, &
                                                  lenlp,lprlp,    &
                                                  dlslp,gamlp,dlpdlr,&
                                                  dlpdlt,dledlt,dedt,iter)
            call lepton_radiation_eos(temk,roye,lete, &
                                                  lenlp,lprlp,     &
                                                  lslp,gamlp,dlpdlr,&
                                                  dlpdlt,dledlt,dedt,iter)

            do i=1,np
               if (iter(i).ne.0) then
                  fun(i)=sby(i)+10.0_rk**lslp(i)   &
                       + fth*eclmb(i)/(mtoe*temp(i)) - entp(i)
                  fun(i)=fun(i)

                  dfun= SUM(1.5_rk*nnuc(i,:)/temp(i)       &
                            +dgnucdt(i,:)/ktom/gnuc(i,:))       &
                       +( 10._rk**dlslp(i)-10._rk**lslp(i) )/ &
                        ( (dtemk(i)-temk(i))*ktom )     &
                       - fth*eclmb(i)/(mtoe*temp(i)**2)
                  if (dfun .gt. 0.0_rk) then
                     dx(i) = fun(i)/dfun
                  else
                     dx(i)=0.0_rk
                     lerr=.TRUE.
                  endif
                  temp(i)=temp(i)-dx(i)
                  if (abs(dx(i)/temp(i)).lt.xacc.or.    &
                       abs(fun(i)/entp(i)).lt.racc) iter(i)=0
               endif
            enddo
            if (lerr) then
               call show_error_screen("eostj","something wrong; mode= ",mode)
               ler=.true.
               return
            endif
            itsum=SUM(iter(:))
            if (itsum.eq.0) EXIT
         enddo
! -- end iteration-loop

! -- root found
         do i = 1, np
            enlp  = 10.0_rk**lenlp(i)
            prlp  = 10.0_rk**lprlp(i)
            pclmb = eclmb(i)/3.0_rk
            dledltl = dledlt(i)
!            if  (dledltl .eq. 0.0_rk)   dledltl = eps
            if (abs(dledltl) .lt. eps) dledltl = eps
            etae(i) = 10.0_rk**lete(i)
            eden(i) = enby(i) + enlp + eclmb(i) + enbyrst(i)
            pres(i) = prby(i) + prlp + pclmb
            entp(i) = sby(i)  + 10.0_rk**lslp(i)   &
                 + fth*eclmb(i)/(mtoe*temp(i))
            gamm(i) = (prlp*dlpdlr(i) + prby(i) + fth*pclmb +   &
                 (prlp*dlpdlt(i) + prby(i))**2/ &
                 (enlp*dledltl   + enby(i))   ) / pres(i)
         enddo

      case default
        raise_abort("(): ")
      end select


    end subroutine boltzmann_gas



  end module boltzman_gas_eos

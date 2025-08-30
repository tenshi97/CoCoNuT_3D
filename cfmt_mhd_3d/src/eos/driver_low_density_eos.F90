!>
!> \verbatim
!> this module / subroutine is the driver for the low density EoS
!> that call eventually a) the baryonic EoS
!>                      b) lepton_radiation EoS
!>                      c) low-density NSE
!>                      d) approximative NSE
!>
!>  Author: A. Marek, MPA, Nov. 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module driver_low_density_eos

  use precision
  implicit none
    
  private

  public low_density_eos, in_low_density_eos
    
contains

!      subroutine init_low_density_eos(tbfile)
!        use precision
!        use abort
!        use lepton_radiation_low_den_eos, only : lepton_radiation_eos_init
!
!        implicit none
!        character(*) tbfile
!
!        call lepton_radiation_eos_init(tbfile)
!      end subroutine init_low_density_eos

!> \verbatim
!>
!>  Determine whether grid point is in low-density EoS    
!>
!> Author : L. Huedepohl
!>=======================================================================
!> \endverbatim
!>
!> \param rho density
!> \param tem temperature
!> \param Ye
!> \param is_in logical
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  elemental function in_low_density_eos(rho, tem, yee) result(is_in)
    use precision
    use eos_data_type, only : lept_radi_table

    implicit none

    real(kind=rk),  intent(in) :: rho, tem, yee
    logical           :: is_in
#ifndef ONEMG_EOS
    is_in = (                                     rho < lept_radi_table(1)%romax              ) .and. &
            (rho * yee < lept_radi_table(1)%royemax .and. rho * yee >= lept_radi_table(1)%royemin ) .and. &
            (tem       < lept_radi_table(1)%ttmax .and. tem       >= lept_radi_table(1)%ttmin )
#else
      ! Special case for ONeMg progenitors
      !
      !          T
      !          ^         |
      !          |         |  high. den. EoS
      !          |         |
      ! tem_cut >|         \-----\
      !          |               |
      !          | low. den. EoS |
      !          |               |
      !          0---------------|-----------> rho
      !                    ^     ^           
      !                    |     |_ romax_high_ld
      !                    |_ romax_low_ld
      !
    is_in = ((tem >= lept_radi_table(1)%tem_cut .and. rho < lept_radi_table(1)%romax_low_ld ) .or. &
             (tem <  lept_radi_table(1)%tem_cut .and. rho < lept_radi_table(1)%romax_high_ld)        ) .and. &
            (rho * yee < lept_radi_table(1)%royemax .and. rho * yee >= lept_radi_table(1)%royemin        ) .and. &
            (tem       < lept_radi_table(1)%ttmax .and. tem       >= lept_radi_table(1)%ttmin        )
#endif
  end function in_low_density_eos

!> \verbatim
!>
!>  Calculate low-density EoS
!>
!> Author : A. Marek, originally M. Rampp
!> \endverbatim
!>
!> \param rho density
!> \param tem temperature
!> \param xi composition
!> \param xhrep mass fraction of representative nucleus (not needed,
!>              just introduced to have the same interface as in
!>               high-density EoS)
!> \param za array of mass and charge number of representative heavy
!>           nucleus. (not needed, just introduced to have the same 
!>           interface as in high-density EoS)
!> \param ede energy density
!> \param pre pressure
!> \param ga  gamma
!> \param st  entropy
!> \param cu  neutrino chemical potential
!> \param ce  electron chemical potential
!> \param cn  neutron chemical potential
!> \param cp  proton chemical potential
!> \param mode calculation mode
!> \param nsemode nse calculation mode
!> \param error  error flag
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine low_density_eos(rho, tem, xi, xhrep, za, ede, pre, ga, st, cu, ce, cn, cp, mode, nsemode, error)
    use precision
    use abort


    use nucparam, only: pc_nuc
    use phycon
    use boltzman_gas_eos, only : boltzmann_gas

!#if LOW_EOS_NSE_FLAG == 2
    use approximative_low_density_nse, only : approximative_nse
!#endif
    use configure

    implicit none
    
    integer(kind=ik), intent(in)    :: mode, nsemode
    real(kind=rk),    intent(in)    :: rho(:)
    real(kind=rk),    intent(out)   :: pre(:),         &
                                       ga(:),          &
                                       za(:,:),        &
                                       xhrep(:),       &
                                       cu(:),          &
                                       ce(:),          &
                                       cn(:),          &
                                       cp(:)


    
    real(kind=rk),    intent(inout) :: tem(:), xi(:,:), ede(:), st(:)
    logical, intent(out)   :: error

    integer(kind=ik), dimension(size(rho))               :: info_wrk
    real(kind=rk), dimension(size(rho))                  :: tt_wrk
    real(kind=rk), dimension(size(xi,dim=2)-1)           :: mn_wrk
    real(kind=rk), dimension(size(rho),size(xi,dim=2)-1) :: cn_wrk,dn_wrk
    integer(kind=ik) :: k, nn, nuc

    nuc=size(xi,dim=2)
        
    error = .false.

    do nn=1,nuc-1
       mn_wrk(nn)=pc_nuc(nn,3)+pc_nuc(nn,2)*wc_mb
       ! write(*,'(I4,3F6.1)') nn,pc_nuc(nn,1),pc_nuc(nn,2), &
       !          (mn_tj(nn)-pc_nuc(nn,2)*wc_mn)/pc_nuc(nn,2)
    enddo
    
    tt_wrk = tem * pc_kmev
    if (mode .eq. 4) then
       st = st*rho/pc_mb
    endif
    
    dn_wrk(:,1:nuc-1) = spread(rho(:),dim=2,ncopies=nuc-1)*max(xi(:,1:nuc-1),1e-25_rk) &
                     / (spread(pc_nuc(1:nuc-1,2),dim=1,ncopies=size(rho))*pc_mb)

    if (config%low_den_nse_eos .eq. 2) then
       call approximative_nse(tt_wrk,rho,xi(:,nuc),                                 &
                       dn_wrk,mn_wrk,ce,cn_wrk,ede,pre,st,ga,mode,error)
    else
       call boltzmann_gas(tt_wrk,rho,xi(:,nuc),dn_wrk,mn_wrk,ce, &
                           cn_wrk,ede,pre,st,ga,mode,info_wrk)
       error=ANY(info_wrk/=0)
    endif
    
    if (error) return
        
    st  = st/rho*pc_mb

    ce = ce * tt_wrk
    cn = cn_wrk(:,1) * tt_wrk
    cp = cn_wrk(:,2) * tt_wrk
    cu = cp + ce - cn
    
    if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 4) &
         tem = tt_wrk * pc_mevk
    
    if (config%low_den_nse_eos .gt. 1) then
       do k = 1, size(rho)
          xi(k,1:nuc-1)=dn_wrk(k,1:nuc-1)*pc_nuc(1:nuc-1,2)/(rho(k)/pc_mb)
       enddo
    endif
    
    select case(nsemode)
    case (1) 
       xhrep(:) = 0.0_rk
       za(:,:) = 0.0_rk
    case(0)   
       ! do nothing
    case default
       raise_abort("tjeos(): nsemode not implemented")
        end select
      end subroutine low_density_eos


    end module driver_low_density_eos

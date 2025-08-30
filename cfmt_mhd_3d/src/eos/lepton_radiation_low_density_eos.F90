#undef IBM_OPT
!>
!> \verbatim
!> this module provides the variables and procedured
!> needed for computing the lepton-radiation part of the
!> low-density EoS
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
module lepton_radiation_low_den_eos

  use precision
  
  implicit none
  private
  public lepton_radiation_eos , lepton_radiation_eos_init
contains
!> \verbatim
!>
!>  Compute 
!>
!> Author : M. Rampp
!> \endverbatim
!>
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
    use precision
  
    implicit none
      

    REAL(kind=rk) :: d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
    INTEGER(kind=ik) ::i,j,l
    REAL(kind=rk) ::d1d2,cl(16),wt(16,16),x(16)
    SAVE wt
      DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*             &
           0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,        &
           1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,          &
           -6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,         &
           10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,          &
           -2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,         &
           2,-2,2*0,-1,1/

    d1d2=d1*d2

    do i=1,4
       x(i)=y(i)
       x(i+4)=y1(i)*d1
       x(i+8)=y2(i)*d2
       x(i+12)=y12(i)*d1d2
    enddo

    cl=MATMUL(wt,x)
  
    l=0
    do i=1,4
       do j=1,4
          l=l+1
          c(i,j)=cl(l)
       enddo
    enddo
    
    return
  END SUBROUTINE bcucof


!> \verbatim
!>
!>  Initialize the lepton-radiation EoS
!>
!> Author : M. Rampp
!> \endverbatim
!>
!> \param tbfile
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
 subroutine lepton_radiation_eos_init(tbfile)
   use precision
   use abort
   use eos_data_type, only : lept_radi_table
   use data_lep_rad_table, only : load_lepton_radiation_table,be ,bp
      
   implicit none
   integer(kind=ik) :: ierr,iro,ite,idl,idu,itl,itu
   character(*)    :: tbfile

   real(kind=rk) :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u     !helpers for intrp

   call load_lepton_radiation_table(tbfile)

!-- store matrices for bicubic interpolation


!CDIR IEXPAND (bcucof)
      do iro=2,lept_radi_table(1)%nro-2
         do ite=2,lept_radi_table(1)%ntt-2
            
            itl=ite
            itu=ite+1
            idl=iro
            idu=iro+1

            y(1)=lept_radi_table(1)%le(itl,idl)
            y(2)=lept_radi_table(1)%le(itu,idl)
            y(3)=lept_radi_table(1)%le(itu,idu)
            y(4)=lept_radi_table(1)%le(itl,idu)


            y1(1)=(lept_radi_table(1)%le(itl+1,idl)-  &
                   lept_radi_table(1)%le(itl-1,idl)) *&
                    0.5_rk*lept_radi_table(1)%dlttab

            y1(2)=(lept_radi_table(1)%le(itu+1,idl)-  &
                   lept_radi_table(1)%le(itu-1,idl)) *&
                   0.5_rk*lept_radi_table(1)%dlttab

            y1(3)=(lept_radi_table(1)%le(itu+1,idu)-  &
                   lept_radi_table(1)%le(itu-1,idu)) *&
                   0.5_rk*lept_radi_table(1)%dlttab

            y1(4)=(lept_radi_table(1)%le(itl+1,idu)-  &
                   lept_radi_table(1)%le(itl-1,idu))* &
                   0.5_rk*lept_radi_table(1)%dlttab

            y2(1)=(lept_radi_table(1)%le(itl,idl+1)-  &
                   lept_radi_table(1)%le(itl,idl-1))* &
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(2)=(lept_radi_table(1)%le(itu,idl+1)-  &
                   lept_radi_table(1)%le(itu,idl-1))* &
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(3)=(lept_radi_table(1)%le(itu,idu+1)-  &
                   lept_radi_table(1)%le(itu,idu-1)) *&
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(4)=(lept_radi_table(1)%le(itl,idu+1)-  &
                   lept_radi_table(1)%le(itl,idu-1))* &
                   0.5_rk*lept_radi_table(1)%dldtab


            y12(1)=(lept_radi_table(1)%le(itl+1,idl+1)-&
                    lept_radi_table(1)%le(itl+1,idl-1)-&
                    lept_radi_table(1)%le(itl-1,idl+1)+&
                    lept_radi_table(1)%le(itl-1,idl-1))*&
                    0.25_rk*lept_radi_table(1)%dldtab*  &
                          lept_radi_table(1)%dlttab

            y12(2)=(lept_radi_table(1)%le(itu+1,idl+1)- &
                    lept_radi_table(1)%le(itu+1,idl-1)- &
                    lept_radi_table(1)%le(itu-1,idl+1)+ &
                    lept_radi_table(1)%le(itu-1,idl-1))*&
                    0.25_rk*lept_radi_table(1)%dldtab*  &
                    lept_radi_table(1)%dlttab

            y12(3)=(lept_radi_table(1)%le(itu+1,idu+1) -&
                    lept_radi_table(1)%le(itu+1,idu-1)- &
                    lept_radi_table(1)%le(itu-1,idu+1)+ &
                    lept_radi_table(1)%le(itu-1,idu-1))*&
                    0.25_rk*lept_radi_table(1)%dldtab*  &
                          lept_radi_table(1)%dlttab

            y12(4)=(lept_radi_table(1)%le(itl+1,idu+1)- &
                    lept_radi_table(1)%le(itl+1,idu-1)- &
                    lept_radi_table(1)%le(itl-1,idu+1)+ &
                    lept_radi_table(1)%le(itl-1,idu-1))*&
                    0.25_rk*lept_radi_table(1)%dldtab*  &
                          lept_radi_table(1)%dlttab

            x1l=lept_radi_table(1)%ltt_ld(itl)
            x1u=lept_radi_table(1)%ltt_ld(itu)
            x2l=lept_radi_table(1)%lro_ld(idl)
            x2u=lept_radi_table(1)%lro_ld(idu)

            call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,be(1,1,ite,iro))


            y(1)=lept_radi_table(1)%lp(itl,idl)
            y(2)=lept_radi_table(1)%lp(itu,idl)
            y(3)=lept_radi_table(1)%lp(itu,idu)
            y(4)=lept_radi_table(1)%lp(itl,idu)

            y1(1)=(lept_radi_table(1)%lp(itl+1,idl)- &
                   lept_radi_table(1)%lp(itl-1,idl))*&
                   0.5_rk*lept_radi_table(1)%dlttab

            y1(2)=(lept_radi_table(1)%lp(itu+1,idl)- &
                   lept_radi_table(1)%lp(itu-1,idl))*&
                   0.5_rk*lept_radi_table(1)%dlttab

            y1(3)=(lept_radi_table(1)%lp(itu+1,idu)- &
                   lept_radi_table(1)%lp(itu-1,idu))*&
                   0.5_rk*lept_radi_table(1)%dlttab

            y1(4)=(lept_radi_table(1)%lp(itl+1,idu)- &
                   lept_radi_table(1)%lp(itl-1,idu))*&
                   0.5_rk*lept_radi_table(1)%dlttab

            y2(1)=(lept_radi_table(1)%lp(itl,idl+1)- &
                   lept_radi_table(1)%lp(itl,idl-1))*&
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(2)=(lept_radi_table(1)%lp(itu,idl+1)- &
                   lept_radi_table(1)%lp(itu,idl-1))*&
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(3)=(lept_radi_table(1)%lp(itu,idu+1)- &
                   lept_radi_table(1)%lp(itu,idu-1))*&
                   0.5_rk*lept_radi_table(1)%dldtab

            y2(4)=(lept_radi_table(1)%lp(itl,idu+1)- &
                   lept_radi_table(1)%lp(itl,idu-1))*&
                   0.5_rk*lept_radi_table(1)%dldtab

            y12(1)=(lept_radi_table(1)%lp(itl+1,idl+1)-&
                    lept_radi_table(1)%lp(itl+1,idl-1)-&
                    lept_radi_table(1)%lp(itl-1,idl+1)+&
                    lept_radi_table(1)%lp(itl-1,idl-1))* &
                    0.25_rk*lept_radi_table(1)%dldtab* &
                          lept_radi_table(1)%dlttab

           y12(2)=(lept_radi_table(1)%lp(itu+1,idl+1)-&
                   lept_radi_table(1)%lp(itu+1,idl-1)-&
                   lept_radi_table(1)%lp(itu-1,idl+1)+&
                   lept_radi_table(1)%lp(itu-1,idl-1))*&
                   0.25_rk*lept_radi_table(1)%dldtab* &
                          lept_radi_table(1)%dlttab

           y12(3)=(lept_radi_table(1)%lp(itu+1,idu+1)-&
                   lept_radi_table(1)%lp(itu+1,idu-1)-&
                   lept_radi_table(1)%lp(itu-1,idu+1)+&
                   lept_radi_table(1)%lp(itu-1,idu-1))*&
                   0.25_rk*lept_radi_table(1)%dldtab*&
                          lept_radi_table(1)%dlttab

           y12(4)=(lept_radi_table(1)%lp(itl+1,idu+1)-&
                   lept_radi_table(1)%lp(itl+1,idu-1)-&
                   lept_radi_table(1)%lp(itl-1,idu+1)+&
                   lept_radi_table(1)%lp(itl-1,idu-1))*&
                   0.25_rk*lept_radi_table(1)%dldtab*&
                          lept_radi_table(1)%dlttab

           x1l=lept_radi_table(1)%ltt_ld(itl)
           x1u=lept_radi_table(1)%ltt_ld(itu)
           x2l=lept_radi_table(1)%lro_ld(idl)
           x2u=lept_radi_table(1)%lro_ld(idu)
           call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,bp(1,1,ite,iro))

        enddo
     enddo

     lept_radi_table(1)%be => be
     lept_radi_table(1)%bp => bp

     
 end subroutine lepton_radiation_eos_init
  
!-----------------------------------------------------------------------
#define BICUBIC
!> \verbatim
!>
!>  Compute the lepton-radiation EoS
!>
!> Author : M. Rampp, IBM-version by A. Marek
!> \endverbatim
!>
!> \param temp temperature
!> \param dens density
!> \param leta
!> \param lenq
!> \param lprs pressure
!> \param lent entropy
!> \param gama adiabatic index
!> \param dlpdlr dPressure/dDensity
!> \param dlpdlt dPressure/dTemperature
!> \param dledt
!> \param dedt
!> \param iter iteration steps
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
   subroutine lepton_radiation_eos(temp,dens,leta,leng,lprs,lent,gama,     &
                         dlpdlr,dlpdlt,dledlt,dedt,iter)

!=======================================================================
!                    Lepton plus radiation EOS                         =
!      for leptons (electrons and positrons) treated as ideal          =
!      Fermi-gases and radiation tabulated in the temperature-         =
!      density range  (10^5 K .le. T .le. 10^12 K) and                 =
!      (10^-10 g/cc .le. rho*Ye .le. 10^12 g/cc) with 20 values        =
!      per decade, i.e. 141 T-values and 441 rho*Ye-values,            =
!      where rho*Ye/mu is the number density of ionization             =
!      electrons (for the atomic mass unit, mu = 1.66^-24 g            =
!      was used).                                                      =
!                                                                      =
!      INPUT:                                                          =
!            temp(1:np)       vector of temperatures [in K]            =
!            dens(1:np)       vector of density*Ye [in g/cc]           =
!                             (= rho*Ye = mu*nb*Ye = mu*ne)            =
!      TABULATED QUANTITIES (stored in module lreos_data)              =
!            lept_radi_table(1)%ntt               number of log10(T)-values (integer)      =
!            lept_radi_table(1)nro               number of log10(rho*Ye)-values           =
!                             (integer)                                =
!            lept_radi_table(1)%ltt_ld(1:nt)       log10(T[K]): T-values                    =
!            lept_radi_table(1)%lro_ld(1:nr)       log10(rho*Ye[g/cc]): rho*Ye-values       = 
!            lept_radi_table(1)%let(1:nt,1:nr)   log10(etae): electron degeneracy         =
!                             parameter (chem. pot./kT, incl.          =
!                             rest-mass energy)                        =
!            lept_radi_table(1)%le(1:nt,1:nr)    log10(E[erg/cc]): internal energy        =
!                             density of e- and e+ and radiation,      =
!                             plus rest-mass energy density of         =
!                             pairs                                    =
!            lept_radi_table(1)%lp(1:nt,1:nr)    log10(P[erg/cc]): pressure of e-         =
!                             and e+ and radiation                     =
!            lept_radi_table(1)%ls(1:nt,1:nr)    log10(s/k[1/cc]): entropy density/k      =
!                             for e- and e+ and radiation              =
!            glept_radi_table(1)%am(1:nt,1:nr)   dln(P)/dln(rho)|_s:  adiabatic           =
!                             index of leptons plus radiation          =
!      OUTPUT:                                                         =
!            leta(1:np)       interpolated logarithm of electron       =
!                             degeneracy parameter incl. electron      =
!                             rest mass contribution                   =
!            leng(1:np)       interpolated logarithm of internal       =
!                             energy density plus pair rest-mass       =
!                             energy [in erg/cc]                       =
!            lprs(1:np)       interpolated logarithm of pressure       =
!                             [in erg/cc]                              =
!            lent(1:np)       interpolated logarithm of entropy        =
!                             density/k [in 1/cc]                      =
!            gama(1:np)       interpolated adiabatic index             =
!            dlpdlr(1:np)     logarithmic derivative (*) of pressure   =
!                             for density                              =
!            dlpdlt(1:np)     logarithmic derivative (*) of pressure   =
!                             for temperature                          =
!            dledlt(1:np)     logarithmic derivative (*) of energy     =
!                             density for temperature                  =
!
!                             (*): if directive BICUBIC is set, the    =
!                                  derivatives are obtained by a       =
!                                  bicubic interpolation ala           =
!                                  Numerical Recipes (F77) Chap.~3.6.  =
!-
!            np               number of grid points                    =
!
!----------------------------------------------------------------------=
!                                                                      =
!            Author:  H.-Thomas Janka, MPA;  Version: June 99          =
!                     (all rights reserved)                            =
!                                                                      =
!            CRAY Version: (inlined) calls to bcuint_vec do not        =
!                          vectorize with f90                          =
!=======================================================================
     use precision
     use abort
     use eos_data_type, only : lept_radi_table
     use specfun, only : fastlog10

     implicit none

     integer(kind=ik) ,intent(in) :: iter(:)
     real(kind=rk) ,intent(in)  :: temp(:), dens(:)
     real(kind=rk) ,intent(out) :: leta(:), leng(:), lprs(:), lent(:), &
                                   gama(:), dlpdlr(:), dlpdlt(:),      &
                                   dledlt(:), dedt(:)

     real(kind=rk), dimension(size(iter)) :: ltemp, ldens

     integer(kind=ik) :: i, itl, idl, itu, idu, m, n, np
     real(kind=rk)    :: tsft, dsft, f1, f2, t, u



#undef IBM_OPT

#ifdef IBM_OPT
     real(kind=rk) :: iter_float(size(iter))
     real(kind=rk) :: itl_float,idl_float
#endif /* IBM_OPT */

!=======================================================================

     np=size(dens)

!      nt=lept_radi_table(1)%ntt
!      nr=lept_radi_table(1)%nro


#ifdef IBM_OPT
      iter_float(:)=-abs(real(iter(:),kind=rk))

#define POWER_FSEL(x,y) x=FSEL(iter_float(i),(x),(y))

#endif /* IBM_OPT */

#ifdef IBM_OPT
      do i=1,np
         POWER_FSEL(ltemp(i),log10(temp(i)))
         POWER_FSEL(ldens(i),log10(dens(i)))

         itl_float=real(itl,kind=rk)
         
         POWER_FSEL(itl_float,((ltemp(i)-lept_radi_table(1)%ltt_ld(1))* 
                             lept_radi_table(1)%dlttab)+1)
         itl=int(itl_float)

         idl_float=real(idl,kind=rk)

         POWER_FSEL(idl_float,((ldens(i)-lept_radi_table(1)%lro_ld(1))* 
                             lept_radi_table(1)%dldtab)+1)

         idl=int(idl_float)

         itu  = itl+1
         idu  = idl+1

         POWER_FSEL(tsft,(ltemp(i)-lept_radi_table(1)%ltt_ld(itl))*   
                         lept_radi_table(1)%dlttab )
         POWER_FSEL(dsft,(ldens(i)-lept_radi_table(1)%lro_ld(idl))*   
                         lept_radi_table(1)%dldtab)

         POWER_FSEL(f1,lept_radi_table(1)%let(itl,idl) +            
                 (lept_radi_table(1)%let(itu,idl) -            
                  lept_radi_table(1)%let(itl,idl))*tsft)

         POWER_FSEL(f2,lept_radi_table(1)%let(itl,idu) +            
                 (lept_radi_table(1)%let(itu,idu) -            
                  lept_radi_table(1)%let(itl,idu))*tsft)

         POWER_FSEL(leta(i),f1+ (f2-f1)*dsft)

         POWER_FSEL(f1,lept_radi_table(1)%le(itl,idl) +              
                 (lept_radi_table(1)%le(itu,idl) -              
                  lept_radi_table(1)%le(itl,idl))*tsft)

         POWER_FSEL(f2,lept_radi_table(1)%le(itl,idu) +         
                 (lept_radi_table(1)%le(itu,idu) -              
                  lept_radi_table(1)%le(itl,idu))*tsft)
         
         POWER_FSEL(leng(i),f1 + (f2-f1)*dsft)

         POWER_FSEL(f1,lept_radi_table(1)%ls(itl,idl) +         
                 (lept_radi_table(1)%ls(itu,idl) -             
                  lept_radi_table(1)%ls(itl,idl))*tsft)

         POWER_FSEL(f2,lept_radi_table(1)%lp(itl,idu) +        
                 (lept_radi_table(1)%lp(itu,idu) -             
                  lept_radi_table(1)%lp(itl,idu))*tsft)

         POWER_FSEL(lprs(i),f1 + (f2-f1)*dsft)

         POWER_FSEL(f1,lept_radi_table(1)%ls(itl,idl) +        
                 (lept_radi_table(1)%ls(itu,idl) -             
                  lept_radi_table(1)%ls(itl,idl))*tsft)

         POWER_FSEL(f2,lept_radi_table(1)%ls(itl,idu) +        
                 (lept_radi_table(1)%ls(itu,idu) -             
                  lept_radi_table(1)%ls(itl,idu))*tsft)


         POWER_FSEL(lent(i),f1 + (f2-f1)*dsft)

         POWER_FSEL(f1,lept_radi_table(1)%gam(itl,idl) +       
                 (lept_radi_table(1)%gam(itu,idl) -            
                  lept_radi_table(1)%gam(itl,idl))*tsft)

         POWER_FSEL(f2,lept_radi_table(1)%gam(itl,idu) +       
                 (lept_radi_table(1)%gam(itu,idu) -            
                  lept_radi_table(1)%gam(itl,idu))*tsft)

         POWER_FSEL(gama(i),f1 + (f2-f1)*dsft)

#ifdef BICUBIC
!            call bcuint_vec(tsft,dsft,dlttab,dldtab,be(1,1,itl,idl),
!     &           ansy,ansy1,ansy2)
!            dledlt(i)=ansy1
!            call bcuint_vec(tsft,dsft,dlttab,dldtab,bp(1,1,itl,idl),
!     &           ansy,ansy1,ansy2)
!            dlpdlt(i)=ansy1
!            dlpdlr(i)=ansy2

           t=tsft
           u=dsft
           m=itl
           n=idl

           POWER_FSEL(dledlt(i),(                              
            ( ( ( (3._rk*lept_radi_table(1)%be(4,4,m,n)*t+    
                   2._rk*lept_radi_table(1)%be(3,4,m,n))*t+   
                         lept_radi_table(1)%be(2,4,m,n) )  * u
               +  (3._rk*lept_radi_table(1)%be(4,3,m,n)*t+    
                   2._rk*lept_radi_table(1)%be(3,3,m,n))*t+   
                         lept_radi_table(1)%be(2,3,m,n) ) * u 
               +  (3._rk*lept_radi_table(1)%be(4,2,m,n)*t+    
                   2._rk*lept_radi_table(1)%be(3,2,m,n))*t+   
                         lept_radi_table(1)%be(2,2,m,n) )* u  
               +  (3._rk*lept_radi_table(1)%be(4,1,m,n)*t+    
                   2._rk*lept_radi_table(1)%be(3,1,m,n))*t+   
                         lept_radi_table(1)%be(2,1,m,n)       
                        ) * lept_radi_table(1)%dlttab)

           POWER_FSEL(dlpdlt(i),(                                      
            ( ( ( (3._rk*lept_radi_table(1)%bp(4,4,m,n)*t+    
                   2._rk*lept_radi_table(1)%bp(3,4,m,n))*t+   
                         lept_radi_table(1)%bp(2,4,m,n) )  * u
               +  (3._rk*lept_radi_table(1)%bp(4,3,m,n)*t+    
                   2._rk*lept_radi_table(1)%bp(3,3,m,n))*t+   
                         lept_radi_table(1)%bp(2,3,m,n) ) * u 
               +  (3._rk*lept_radi_table(1)%bp(4,2,m,n)*t+    
                   2._rk*lept_radi_table(1)%bp(3,2,m,n))*t+   
                         lept_radi_table(1)%bp(2,2,m,n) )* u  
               +  (3._rk*lept_radi_table(1)%bp(4,1,m,n)*t+    
                   2._rk*lept_radi_table(1)%bp(3,1,m,n))*t+   
                         lept_radi_table(1)%bp(2,1,m,n)       
                        ) * lept_radi_table(1)%dlttab)

           POWER_FSEL(dlpdlr(i),(                                      
            ( ( ( (3._rk*lept_radi_table(1)%bp(4,4,m,n)*u+    
                   2._rk*lept_radi_table(1)%bp(4,3,m,n))*u+   
                         lept_radi_table(1)%bp(4,2,m,n) )  * t
               +  (3._rk*lept_radi_table(1)%bp(3,4,m,n)*u+    
                   2._rk*lept_radi_table(1)%bp(3,3,m,n))*u+   
                         lept_radi_table(1)%bp(3,2,m,n) ) * t 
               +  (3._rk*lept_radi_table(1)%bp(2,4,m,n)*u+    
                   2._rk*lept_radi_table(1)%bp(2,3,m,n))*u+   
                         lept_radi_table(1)%bp(2,2,m,n) )* t  
               +  (3._rk*lept_radi_table(1)%bp(1,4,m,n)*u+    
                   2._rk*lept_radi_table(1)%bp(1,3,m,n))*u+   
                         lept_radi_table(1)%bp(1,2,m,n)       
                        ) * lept_radi_table(1)%dldtab)

           POWER_FSEL(f1,lept_radi_table(1)%le(itl,idl) +  
                  (lept_radi_table(1)%le(itl,idu) -  
                   lept_radi_table(1)%le(itl,idl))*dsft)

           POWER_FSEL(f2,lept_radi_table(1)%le(itu,idl) +  
                  (lept_radi_table(1)%le(itu,idu) -  
                   lept_radi_table(1)%le(itu,idl))*dsft)

           POWER_FSEL(dedt(i),(f2-f1)*lept_radi_table(1)%dlttab*10.0_rk**leng(i)/temp(i))

#else /* BICUBIV */

           POWER_FSEL(f1,lept_radi_table(1)%lp(itl,idl) + 
                  ( lept_radi_table(1)%lp(itl,idu) - 
                    lept_radi_table(1)%lp(itl,idl))*dsft)

           POWER_FSEL(f2,lept_radi_table(1)%lp(itu,idl) + 
                  ( lept_radi_table(1)%lp(itu,idu) - 
                    lept_radi_table(1)%lp(itu,idl))*dsft)

           POWER_FSEL(dlpldt(i),(f2-f1)*lept_radi_table(1)%dlttab)

           POWER_FSEL(f1,lept_radi_table(1)%le(itl,idl) + 
                  ( lept_radi_table(1)%le(itl,idu) - 
                    lept_radi_table(1)%le(itl,idl))*dsft)

           POWER_FSEL(f2,lept_radi_table(1)%le(itu,idl) + 
                  (lept_radi_table(1)%le(itu,idu) - 
                   lept_radi_table(1)%le(itu,idl))*dsft)

           POWER_FSEL(dledlt(i),(f2-f1)*lept_radi_table(1)%dlttab)

           POWER_FSEL(dedt(i),dledlt(i)*10.0_rk**leng(i)/temp(i))

           POWER_FSEL(f1,lept_radi_table(1)%lp(itl,idl) + 
                  (lept_radi_table(1)%lp(itu,idl) - 
                   lept_radi_table(1)%lp(itl,idl))*tsft)

           POWER_FSEL(f2,lept_radi_table(1)%lp(itl,idu) + 
                  (lept_radi_table(1)%lp(itu,idu) - 
                   lept_radi_table(1)%lp(itl,idu))*tsft)

           POWER_FSEL(dlpdlr(i),(f2-f1)*lept_radi_table(1)%dldtab)

#endif /* BICUBIC */


      enddo

#else /* IBM_OPT */

!CDIR IEXPAND (bcuint_vec)

      ltemp(:) = temp(:)
      ldens(:) = dens(:)

      call fastlog10(ltemp, np)
      call fastlog10(ldens, np)

     do i = 1, np

        if (iter(i) .ne.0) then
           itl  = int((ltemp(i)-lept_radi_table(1)%ltt_ld(1))* &
                             lept_radi_table(1)%dlttab)+1
           idl  = int((ldens(i)-lept_radi_table(1)%lro_ld(1))* &
                             lept_radi_table(1)%dldtab)+1
           itu  = itl+1
           idu  = idl+1

           tsft = (ltemp(i)-lept_radi_table(1)%ltt_ld(itl))*   &
                         lept_radi_table(1)%dlttab 
           dsft = (ldens(i)-lept_radi_table(1)%lro_ld(idl))*   &
                         lept_radi_table(1)%dldtab

           f1   = lept_radi_table(1)%let(itl,idl) +            &
                 (lept_radi_table(1)%let(itu,idl) -            &
                  lept_radi_table(1)%let(itl,idl))*tsft

           f2   = lept_radi_table(1)%let(itl,idu) +            &
                 (lept_radi_table(1)%let(itu,idu) -            &
                  lept_radi_table(1)%let(itl,idu))*tsft

           leta(i) = f1 + (f2-f1)*dsft

           f1   = lept_radi_table(1)%le(itl,idl) +              &
                 (lept_radi_table(1)%le(itu,idl) -              &
                  lept_radi_table(1)%le(itl,idl))*tsft

           f2   = lept_radi_table(1)%le(itl,idu) +              &
                 (lept_radi_table(1)%le(itu,idu) -              &
                  lept_radi_table(1)%le(itl,idu))*tsft

           leng(i) = f1 + (f2-f1)*dsft

           f1   = lept_radi_table(1)%lp(itl,idl) +             &
                 (lept_radi_table(1)%lp(itu,idl) -             &
                  lept_radi_table(1)%lp(itl,idl))*tsft

           f2   = lept_radi_table(1)%lp(itl,idu) +             &
                 (lept_radi_table(1)%lp(itu,idu) -             &
                  lept_radi_table(1)%lp(itl,idu))*tsft
           lprs(i) = f1 + (f2-f1)*dsft

           f1   = lept_radi_table(1)%ls(itl,idl) +             &
                 (lept_radi_table(1)%ls(itu,idl) -             &
                  lept_radi_table(1)%ls(itl,idl))*tsft

           f2   = lept_radi_table(1)%ls(itl,idu) +             & 
                 (lept_radi_table(1)%ls(itu,idu) -             &
                  lept_radi_table(1)%ls(itl,idu))*tsft
           lent(i) = f1 + (f2-f1)*dsft

           f1   = lept_radi_table(1)%gam(itl,idl) +            &
                 (lept_radi_table(1)%gam(itu,idl) -            &
                  lept_radi_table(1)%gam(itl,idl))*tsft

           f2   = lept_radi_table(1)%gam(itl,idu) +            &
                 (lept_radi_table(1)%gam(itu,idu) -            &
                  lept_radi_table(1)%gam(itl,idu))*tsft

           gama(i) = f1 + (f2-f1)*dsft

#ifdef BICUBIC
!            call bcuint_vec(tsft,dsft,dlttab,dldtab,be(1,1,itl,idl),
!     &           ansy,ansy1,ansy2)
!            dledlt(i)=ansy1
!            call bcuint_vec(tsft,dsft,dlttab,dldtab,bp(1,1,itl,idl),
!     &           ansy,ansy1,ansy2)
!            dlpdlt(i)=ansy1
!            dlpdlr(i)=ansy2

           t=tsft
           u=dsft
           m=itl
           n=idl

            dledlt(i)= (                                      &
            ( ( ( (3._rk*lept_radi_table(1)%be(4,4,m,n)*t+    &
                   2._rk*lept_radi_table(1)%be(3,4,m,n))*t+   &
                         lept_radi_table(1)%be(2,4,m,n) )  * u&
               +  (3._rk*lept_radi_table(1)%be(4,3,m,n)*t+    &
                   2._rk*lept_radi_table(1)%be(3,3,m,n))*t+   &
                         lept_radi_table(1)%be(2,3,m,n) ) * u &
               +  (3._rk*lept_radi_table(1)%be(4,2,m,n)*t+    &
                   2._rk*lept_radi_table(1)%be(3,2,m,n))*t+   &
                         lept_radi_table(1)%be(2,2,m,n) )* u  &
               +  (3._rk*lept_radi_table(1)%be(4,1,m,n)*t+    &
                   2._rk*lept_radi_table(1)%be(3,1,m,n))*t+   &
                         lept_radi_table(1)%be(2,1,m,n)       &
                        ) * lept_radi_table(1)%dlttab

            dlpdlt(i)= (                                      &
            ( ( ( (3._rk*lept_radi_table(1)%bp(4,4,m,n)*t+    &
                   2._rk*lept_radi_table(1)%bp(3,4,m,n))*t+   &
                         lept_radi_table(1)%bp(2,4,m,n) )  * u&
               +  (3._rk*lept_radi_table(1)%bp(4,3,m,n)*t+    &
                   2._rk*lept_radi_table(1)%bp(3,3,m,n))*t+   &
                         lept_radi_table(1)%bp(2,3,m,n) ) * u &
               +  (3._rk*lept_radi_table(1)%bp(4,2,m,n)*t+    &
                   2._rk*lept_radi_table(1)%bp(3,2,m,n))*t+   &
                         lept_radi_table(1)%bp(2,2,m,n) )* u  &
               +  (3._rk*lept_radi_table(1)%bp(4,1,m,n)*t+    &
                   2._rk*lept_radi_table(1)%bp(3,1,m,n))*t+   &
                         lept_radi_table(1)%bp(2,1,m,n)       &
                        ) * lept_radi_table(1)%dlttab

            dlpdlr(i)= (                                      &
            ( ( ( (3._rk*lept_radi_table(1)%bp(4,4,m,n)*u+    &
                   2._rk*lept_radi_table(1)%bp(4,3,m,n))*u+   &
                         lept_radi_table(1)%bp(4,2,m,n) )  * t&
               +  (3._rk*lept_radi_table(1)%bp(3,4,m,n)*u+    &
                   2._rk*lept_radi_table(1)%bp(3,3,m,n))*u+   &
                         lept_radi_table(1)%bp(3,2,m,n) ) * t &
               +  (3._rk*lept_radi_table(1)%bp(2,4,m,n)*u+    &
                   2._rk*lept_radi_table(1)%bp(2,3,m,n))*u+   &
                         lept_radi_table(1)%bp(2,2,m,n) )* t  &
               +  (3._rk*lept_radi_table(1)%bp(1,4,m,n)*u+    &
                   2._rk*lept_radi_table(1)%bp(1,3,m,n))*u+   &
                         lept_radi_table(1)%bp(1,2,m,n)       &
                        ) * lept_radi_table(1)%dldtab

            f1   = lept_radi_table(1)%le(itl,idl) +  &
                  (lept_radi_table(1)%le(itl,idu) -  &
                   lept_radi_table(1)%le(itl,idl))*dsft

            f2   = lept_radi_table(1)%le(itu,idl) +  &
                  (lept_radi_table(1)%le(itu,idu) -  &
                   lept_radi_table(1)%le(itu,idl))*dsft

            dedt(i) = (f2-f1)*lept_radi_table(1)%dlttab*10.0_rk**leng(i)/temp(i)

#else /* BICUBIC */
            f1   =  lept_radi_table(1)%lp(itl,idl) + &
                  ( lept_radi_table(1)%lp(itl,idu) - &
                    lept_radi_table(1)%lp(itl,idl))*dsft

            f2   =  lept_radi_table(1)%lp(itu,idl) + &
                  ( lept_radi_table(1)%lp(itu,idu) - &
                    lept_radi_table(1)%lp(itu,idl))*dsft

            dlpdlt(i) = (f2-f1)*lept_radi_table(1)%dlttab

            f1   =  lept_radi_table(1)%le(itl,idl) + &
                  ( lept_radi_table(1)%le(itl,idu) - &
                    lept_radi_table(1)%le(itl,idl))*dsft

            f2   = lept_radi_table(1)%le(itu,idl) + &
                  (lept_radi_table(1)%le(itu,idu) - &
                   lept_radi_table(1)%le(itu,idl))*dsft

            dledlt(i) = (f2-f1)*lept_radi_table(1)%dlttab
            dedt(i) = dledlt(i)*10.0_rk**leng(i)/temp(i)

            f1   = lept_radi_table(1)%lp(itl,idl) + &
                  (lept_radi_table(1)%lp(itu,idl) - &
                   lept_radi_table(1)%lp(itl,idl))*tsft

            f2   = lept_radi_table(1)%lp(itl,idu) + &
                  (lept_radi_table(1)%lp(itu,idu) - &
                   lept_radi_table(1)%lp(itl,idu))*tsft

            dlpdlr(i) = (f2-f1)*lept_radi_table(1)%dldtab
#endif /* BICUBIC */

         endif
      enddo

#endif /* IBM_OPT */
    end subroutine lepton_radiation_eos

  end module lepton_radiation_low_den_eos

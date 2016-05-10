! Commented out by roj April 2006
! #include <misc.h>
! #include <params.h>

module pkg_cldoptics

!---------------------------------------------------------------------------------
! Purpose:
!
! Compute cloud optical properties: liquid and ice partical size; emissivity
!
! Author: Byron Boville  Sept 06, 2002, assembled from existing subroutines
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8

! comment out by roj April 2006 --- variables added - directly ...

!  use ppgrid,        only: pcols, pver, pverp

       	integer pcols      ! number of columns (max)
   	integer pver       ! number of vertical levels
  	integer pverp      ! pver + 1
 
 	parameter (pcols  = 1)
!	parameter (pver  = 26)		! Roj, Oct 2008 - modified subroutines to pass pver as INPUT!!  Also removed two surboutine I wasn't using
!	parameter (pverp  = pver +1)

!
! end of roj's new stuff
	
  private
  public :: cldefr, cldems, reitab, reltab

contains

!===============================================================================
  subroutine cldefr(lchnk   ,ncol1, ncol    , &
       landfrac,t       ,rel     ,rei     ,ps      ,pmid    , landm, icefrac, snowh, pver)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud water and ice particle size 
! 
! Method: 
! use empirical formulas to construct effective radii
! 
! Author: J.T. Kiehl, B. A. Boville, P. Rasch
! 
!-----------------------------------------------------------------------

    implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in) :: ncol1, ncol,pver                  ! number of atmospheric columns

    real(r8), intent(in) :: landfrac(pcols)      ! Land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! Ice fraction
    real(r8), intent(in) :: t(pcols,pver)        ! Temperature
    real(r8), intent(in) :: ps(pcols)            ! Surface pressure
    real(r8), intent(in) :: pmid(pcols,pver)     ! Midpoint pressures
    real(r8), intent(in) :: landm(pcols)
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
!
! Output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
    real(r8), intent(out) :: rei(pcols,pver)      ! Ice effective drop size (microns)
!

!++pjr
! following Kiehl
         call reltab(ncol1, ncol, t, landfrac, landm, icefrac, rel, snowh,pver)

! following Kristjansson and Mitchell
         call reitab(ncol1, ncol, t, rei, pver)
!--pjr
!
!
    return
  end subroutine cldefr

!===============================================================================
  subroutine cldems(lchnk   ,ncol1, ncol    ,clwp    ,fice    ,rei     ,emis,pver    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud emissivity using cloud liquid water path (g/m**2)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J.T. Kiehl
! 
!-----------------------------------------------------------------------

    implicit none
!------------------------------Parameters-------------------------------
!
    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361)
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: lchnk, pver                   ! chunk identifier
    integer, intent(in) :: ncol1, ncol             ! number of atmospheric columns

    real(r8), intent(in) :: clwp(pcols,pver)       ! cloud liquid water path (g/m**2)
    real(r8), intent(in) :: rei(pcols,pver)        ! ice effective drop size (microns)
    real(r8), intent(in) :: fice(pcols,pver)       ! fractional ice content within cloud
!
! Output arguments
!
    real(r8), intent(out) :: emis(pcols,pver)       ! cloud emissivity (fraction)
!
!---------------------------Local workspace-----------------------------
!
    integer i,k                 ! longitude, level indices
    real(r8) kabs                   ! longwave absorption coeff (m**2/g)
    real(r8) kabsi                  ! ice absorption coefficient
!
!-----------------------------------------------------------------------
!
    do k=1,pver
       do i=ncol1,ncol
          kabsi = 0.005 + 1./rei(i,k)
          kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
       end do
    end do
!
    return
  end subroutine cldems


!===============================================================================
  subroutine reltab(ncol1, ncol, t, landfrac, landm, icefrac, rel, snowh, pver)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud water size
! 
! Method: 
! analytic formula following the formulation originally developed by J. T. Kiehl
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
 	

! added constant ... Roj April 2006	
	real(r8), parameter :: tmelt = 273.16
!   use physconst,          only: tmelt
!   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol1, ncol, pver
    real(r8), intent(in) :: landfrac(pcols)      ! Land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! Ice fraction
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    real(r8), intent(in) :: landm(pcols)         ! Land fraction ramping to zero over ocean
    real(r8), intent(in) :: t(pcols,pver)        ! Temperature

!
! Output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
!
!---------------------------Local workspace-----------------------------
!
    integer i,k               ! Lon, lev indices
    real(r8) rliqland         ! liquid drop size if over land
    real(r8) rliqocean        ! liquid drop size if over ocean
    real(r8) rliqice          ! liquid drop size if over sea ice
!
!-----------------------------------------------------------------------
!
#ifdef CRM
    rliqocean = 13.0_r8
    rliqice   = 13.0_r8
    rliqland  = 7.0_r8
#else
    rliqocean = 14.0_r8
    rliqice   = 14.0_r8
    rliqland  = 8.0_r8
#endif
    do k=1,pver
       do i=ncol1,ncol
! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * min(1.0_r8,max(0.0_r8,(tmelt-t(i,k))*0.05))
          ! Modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,snowh(i)*10.))
          ! Ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,1.0-landm(i)))
          ! Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * min(1.0_r8,max(0.0_r8,icefrac(i)))
! end jrm
       end do
    end do
  end subroutine reltab

!===============================================================================
  subroutine reitab(ncol1, ncol, t, re,pver)
    !

    integer, intent(in) :: ncol1,ncol,pver
    real(r8), intent(out) :: re(pcols,pver)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8) retab(95)
    real(r8) corr
    integer i
    integer k
    integer index
    !
    !       Tabulated values of re(T) in the temperature interval
    !       180 K -- 274 K; hexagonal columns assumed:
    !
    data retab / 						&
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,	&
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,	&
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 	&
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,	&
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,	&
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, 	&
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, 	&
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,	&
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,	&
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,	&
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,	&
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,	&
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, 	&
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,	&
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,	&
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639/	
    !
    save retab
    !
    do k=1,pver
       do i=ncol1,ncol
          index = int(t(i,k)-179.)
          index = min(max(index,1),94)
          corr = t(i,k) - int(t(i,k))
          re(i,k) = retab(index)*(1.-corr)		&
               +retab(index+1)*corr
          !           re(i,k) = amax1(amin1(re(i,k),30.),10.)
       end do
    end do
    !
    return
  end subroutine reitab

end module pkg_cldoptics

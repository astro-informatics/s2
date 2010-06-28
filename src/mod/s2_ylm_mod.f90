!------------------------------------------------------------------------------
! s2_ylm_mod
!
!! Provides functionality to compute spherical harmonic functions.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 February 2008
!
! Revisions:
!   February 2008 - Jason McEwen
!------------------------------------------------------------------------------

module s2_ylm_mod

  use s2_types_mod
  use s2_error_mod
  use s2_dl_mod
  use s2_sky_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------
  
  public :: &
       s2_ylm_sky, &
       s2_ylm_eval_leg
  

  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  ! None.
  
  
  !---------------------------------------
  ! Global variables
  !---------------------------------------
  
  !! Evalute ylms using associated Legendre functions.
  integer, public, parameter :: S2_YLM_EVAL_METHOD_LEG = 1
  
  !! Evaluate ylms using Wigner functions.
  integer, public, parameter :: S2_YLM_EVAL_METHOD_WIG = 2

  !! Evaluate real part of ylms.
  integer, public, parameter :: S2_YLM_REALITY_REAL = 1

  !! Evaluate imaginary part of ylms.
  integer, public, parameter :: S2_YLM_REALITY_IMAG = 2

  !! Evaluate absolute value of ylms.
  integer, public, parameter :: S2_YLM_REALITY_ABS = 3


  !---------------------------------------
  ! Data types
  !---------------------------------------
  
  ! None.
  
  
  !----------------------------------------------------------------------------
  
  contains
  

    !--------------------------------------------------------------------------
    ! s2_ylm_sky
    !
    !! Evaluate spherical harmonic function Ylm over full sky.
    !!
    !! Variables:
    !!   - el: Harmonic scale l parameter.
    !!   - m: Harmonic azimuthal scale m parameter.
    !!   - nside: Healpix resolution parameter to construct sky at.
    !!   - reality: Reality type to specify whether to compute the real part,
    !!     imaginary part or absolute value of Ylm (since s2_sky_mod only 
    !!     handles real functions on the sphere).
    !!   - pix_scheme_in: Pixelisation scheme for map.  
    !!   - method_in: Method used to compute Ylm (S2_YLM_EVAL_METHOD_LEG for
    !!     assocaite Legendre functions; S2_YLM_EVAL_METHOD_WIG for Wigner dlmn
    !!     functions.
    !!   - ylm: Spherical harmonic function evaluated on the full sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 February 2008
    !
    ! Revisions:
    !   February 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_ylm_sky(el, m, nside, reality, pix_scheme_in, method_in) &
         result(ylm)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      integer, intent(in) :: el, m
      integer, intent(in) :: nside
      integer, intent(in) :: reality
      integer, intent(in), optional :: pix_scheme_in
      integer, intent(in), optional :: method_in
      type(s2_sky) :: ylm

      real(s2_sp), allocatable :: map(:)
      integer :: ipix, npix, fail
      integer :: pix_scheme, method
      real(s2_dp) :: theta, phi

      ! Initialise parameters.
      fail = 0
      npix = 12*nside**2
      if(present(pix_scheme_in)) then
         pix_scheme = pix_scheme_in
      else
         pix_scheme = S2_SKY_RING
      end if
      if(present(method_in)) then
         method = method_in
      else
         method = S2_YLM_EVAL_METHOD_LEG
      end if

      ! Allocate memory.
      allocate(map(0:npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_ylm_sky')
      end if

      ! Compute Ylm values over sky.
      do ipix = 0,npix-1

         ! Compute theta and phi.
         if(pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(nside, ipix, theta, phi)
         else if(pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_ylm_sky')
         end if
         
         ! Evaluate Ylm.
         if(method == S2_YLM_EVAL_METHOD_LEG) then

            select case (reality)

              case (S2_YLM_REALITY_REAL)
                 map(ipix) = real(s2_ylm_eval_leg(el, m, theta, phi))

              case (S2_YLM_REALITY_IMAG)
                 map(ipix) = aimag(s2_ylm_eval_leg(el, m, theta, phi))

              case (S2_YLM_REALITY_ABS)
                 map(ipix) = abs(s2_ylm_eval_leg(el, m, theta, phi))

            end select

         else if(method == S2_YLM_EVAL_METHOD_WIG) then

            select case (reality)

              case (S2_YLM_REALITY_REAL)
                 map(ipix) = real(s2_ylm_eval_wig(el, m, theta, phi))

              case (S2_YLM_REALITY_IMAG)
                 map(ipix) = aimag(s2_ylm_eval_wig(el, m, theta, phi))

              case (S2_YLM_REALITY_ABS)
                 map(ipix) = abs(s2_ylm_eval_wig(el, m, theta, phi))

            end select

         else
            call s2_error(S2_ERROR_YLM_ARG_INVALID, 's2_ylm_sky', &
                 comment_add='Invalid method')
         end if

      end do

      ! Construct sky object from map values.
      ylm = s2_sky_init(map(0:npix-1), nside, pix_scheme)

      ! Free memory.
      deallocate(map)

    end function s2_ylm_sky


    !--------------------------------------------------------------------------
    ! s2_ylm_eval_wig
    !
    !! Evaluate spherical harmonic function Ylm at specified value of theta
    !! and phi using Wigner dlmn function evaluations.
    !!
    !! Notes:
    !!   - This function is considerably slower thatn s2_ylm_eval_leg but is
    !!     included for a consistency check.
    !!
    !! Variables:
    !!   - el: Harmonic scale l parameter.
    !!   - m: Harmonic azimuthal scale m parameter.
    !!   - theta: Spherical coordiante theta value to evalute Ylm at.
    !!   - phi: Spherical coordiante phi value to evalute Ylm at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 February 2008
    !
    ! Revisions:
    !   February 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_ylm_eval_wig(el, m, theta, phi) result(ylm)

      integer, intent(in) :: el, m
      real(s2_dp), intent(in) :: theta, phi
      complex(s2_dpc) :: ylm
      
			real(s2_dp), pointer :: dl(:,:)
      integer :: fail = 0
      complex(s2_dpc) :: I

      I = cmplx(0d0, 1d0)

      allocate(dl(-el:el,-el:el), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_ylm_eval_wig')
      end if
      call s2_dl_beta_operator(dl, theta, el)

      ylm = sqrt( (2d0*el+1d0)/(4d0*PI )) * dl(m,0) &
           * exp(I*m*phi)

      deallocate(dl)

    end function s2_ylm_eval_wig


    !--------------------------------------------------------------------------
    ! s2_ylm_eval_leg
    !
    !! Evaluate spherical harmonic function Ylm at specified value of theta
    !! and phi using associated Legendre function evaluations.
    !!
    !! Variables:
    !!   - el: Harmonic scale l parameter.
    !!   - m: Harmonic azimuthal scale m parameter.
    !!   - theta: Spherical coordiante theta value to evalute Ylm at.
    !!   - phi: Spherical coordiante phi value to evalute Ylm at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 February 2008
    !
    ! Revisions:
    !   February 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_ylm_eval_leg(el, m, theta, phi) result(ylm)

      integer, intent(in) :: el, m
      real(s2_dp), intent(in) :: theta, phi
      complex(s2_dpc) :: ylm

      integer :: mabs
      complex(s2_dpc) :: I

      I = cmplx(0d0, 1d0)

      if(m<0) then
         mabs = -m
      else
         mabs = m
      end if

      ylm = sqrt( (2d0*el+1d0)/(4d0*PI) ) &
           * sqrt( exp(logfact(el-mabs)-logfact(el+mabs)) ) &
           * s2_ylm_plm(el,mabs,cos(theta)) &
           * exp(I*mabs*phi)  

      if(m<0) ylm = (-1)**m * conjg(ylm)      

    end function s2_ylm_eval_leg


    !--------------------------------------------------------------------------
    ! s2_ylm_plm
    !
    !! Computes the associated Legendre function for x, l and m.
    !! Adapted from numerical recipes.
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Computes the associated Legendre polynomial P_m^l(x).
    !!   - Includes (-1)^m factor inside 
    !!
    !! Variables:
    !!   - el: Legendre function l parameter.
    !!   - m: Legendre function m parameter.
    !!   - x: Point to evaluate specified Legendre funtion at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_ylm_plm(el,m,x) result(plm)
            
      implicit none
      
      integer, intent(in) :: el, m
      real(s2_dp), intent(in) :: x
      real(s2_dp) :: plm

      integer :: i, ll
      real(s2_dp) ::  fact, pll, pmm, pmmp1, somx2
      
      if(m<0 .or. m>el .or. abs(x)>1d0) then
         call s2_error(S2_ERROR_YLM_ARG_INVALID, 's2_ylm_plm')
      end if

      ! Compute Pmm.
      pmm=1.                     
      if(m.gt.0) then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
         enddo
      endif

      if(el.eq.m) then
         plm=pmm
      else
         ! Compute Pm,m+1.
         pmmp1=x*(2*m+1)*pmm      
         if(el.eq.m+1) then
            plm=pmmp1
         else                    
            ! Compute Pm,l for l > m+1.
            do ll=m+2,el
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            enddo
            plm=pll
         endif
      endif
      return

    end function s2_ylm_plm
    

    !--------------------------------------------------------------------------
    ! logfact
    !
    !! Computes the natural logarithm of an (integer) factorial.
    !!
    !! Variables:
    !!  - n: Integer to compute factorial of [input].
		!!  - logfactn: Natural logarithm of factorial value computed [output].
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2007
    !
    ! Revisions:
    !   October 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function logfact(n) result(logfactn)

      integer, intent(in) :: n
      real(s2_dp) :: logfactn

      real(s2_dp) :: y, temp, sum, c(6), loggamma, x
      integer :: nn

      if (n < 0) then

         call s2_error(S2_ERROR_YLM_ARG_INVALID, 'logfact', &
              comment_add='Factorial argument negative')

      else

         ! The engine of this function actually calculates the gamma function,
         ! for which the real argument is x = n + 1.

         x = real(n, s2_dp) + 1.0

         ! Table of fitting constants.

         c(1) = 76.18009172947146
         c(2) = - 86.50532032941677
         c(3) = 24.01409824083091
         c(4) = - 1.231739572450155
         c(5) = 0.1208650973866179e-2
         c(6) = - 0.5395239384953e-5

         ! Add up fit.

         temp = x + 5.5 - (x + 0.5) * log(x + 5.5);
         sum = 1.000000000190015
         y = x

         do nn = 1, 6
            y = y + 1.0;
            sum = sum + c(nn) / y;
         end do

         loggamma = - temp + log(2.5066282746310005 * sum / x);

      end if

      ! Finally make explicit the conversion back to log of the factorial.

      logfactn = loggamma

    end function logfact


end module s2_ylm_mod

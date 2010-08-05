!------------------------------------------------------------------------------
! s2_sky_mod -- S2 library sky class
!
!! Provides functionality to support and manipulate a function defined on 
!! the sky(/sphere).  Support is provided to representing the sky in both 
!! real (map) space and in harmonic (alm) space, and also to convert between 
!! the two representations.  One may also scale, dilate, rotate, downsample 
!! and add functions defined on the sky.  The sky pixelisation used is
!! currently based on HEALPix.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module s2_sky_mod

  use s2_types_mod, only: s2_sp, s2_spc, s2_dp, s2_dpc, pi
  use s2_error_mod
  use s2_vect_mod
  use s2_pl_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_sky_init, &
    s2_sky_free, s2_sky_remove_map, &
    s2_sky_compute_alm, s2_sky_compute_alm_iter, &
    s2_sky_compute_map, s2_sky_irregular_invsht, &
    s2_sky_map_convert, &
    s2_sky_der, s2_sky_der_discrete_phi, s2_sky_der_discrete_phi_fovop, &
    s2_sky_der_discrete_theta, s2_sky_der_discrete_theta_fovop, &
    s2_sky_der_discrete_grad, &
    s2_sky_conv, s2_sky_conv_space, s2_sky_convpt_space, &
    s2_sky_conv_space_fovop, s2_sky_convpt_space_weights, &
    s2_sky_offset, s2_sky_scale,  &
    s2_sky_add, s2_sky_add_alm, s2_sky_product, s2_sky_thres, s2_sky_thres_abs, &
    s2_sky_error_twonorm, s2_sky_rms, & ! s2_sky_error_onenorm, s2_sky_error_pnorm, &
    s2_sky_dilate, &
    s2_sky_rotate, s2_sky_rotate_alm, &
    s2_sky_power_map, s2_sky_power_alm, &
    s2_sky_azimuthal_bl, &
    s2_sky_admiss, &
    s2_sky_admiss_dil, s2_sky_fov, &
    s2_sky_extract_ab, s2_sky_extract_ab_fsht, &
    s2_sky_downsample, s2_sky_upsample, &
    s2_sky_draw_dot, &
    s2_sky_write_map_file, s2_sky_write_matmap_file, &
    s2_sky_write_alm_file, s2_sky_write_matalm_file, &
    s2_sky_write_file, s2_sky_io_fits_write, &
    s2_sky_set_lmax, &
    s2_sky_set_nside, &
    s2_sky_get_init, &
    s2_sky_get_nside, &
    s2_sky_get_npix, &
    s2_sky_get_lmax, &
    s2_sky_get_mmax, &
    s2_sky_get_map, &
    s2_sky_get_map_pix, &
    s2_sky_get_alm, &
    s2_sky_get_cl, &
    s2_sky_get_cm, &
    s2_sky_get_pix_scheme, &
    s2_sky_get_map_status, &
    s2_sky_get_alm_status, &
    s2_sky_get_n_param, &
    s2_sky_get_param


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  interface s2_sky_init

     module procedure &
       s2_sky_init_map, &
       s2_sky_init_alm, &
       s2_sky_init_ab, &
       s2_sky_init_fun, &
       s2_sky_init_fun_alm, &
       s2_sky_init_file, &
       s2_sky_init_copy
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Healpix pixelisation ring type.
  integer, public, parameter :: S2_SKY_RING = 1

  !! Healpix pixelisation nest type.
  integer, public, parameter :: S2_SKY_NEST = 2
  
  !! Specifier to indicate function passed to construct sky
  !! is defined on the sphere.
  integer, public, parameter :: S2_SKY_FUN_TYPE_SPHERE = 1

  !! Specifier to indicate function passed to construct sky
  !! is defined on the plane (and will be numerically projected onto the 
  !! sphere).
  integer, public, parameter :: S2_SKY_FUN_TYPE_PLANE = 2

  !! Default function type if fun_type flag is not specified.
  integer, public, parameter :: S2_SKY_DEFAULT_FTYPE = S2_SKY_FUN_TYPE_SPHERE

  !! Map fits file type.
  integer, public, parameter :: S2_SKY_FILE_TYPE_MAP = 1

  !! Full s2_sky file type.
  integer, public, parameter :: S2_SKY_FILE_TYPE_SKY = 2

  !! Alm fits file type.
  integer, public, parameter :: S2_SKY_FILE_TYPE_ALM = 3

  !! Circle method type for setting field-of-view.
  integer, public, parameter :: S2_SKY_FOV_METHOD_CIRCLE = 1

  !! Square method type for setting field-of-view.
  integer, public, parameter :: S2_SKY_FOV_METHOD_SQUARE = 2

  !! Theta derivative type.
  integer, public, parameter :: S2_SKY_DER_TYPE_THETA = 1

  !! Phi derivative type.
  integer, public, parameter :: S2_SKY_DER_TYPE_PHI = 2

  !! Gradient derivative type.
  integer, public, parameter :: S2_SKY_DER_TYPE_GRAD = 3


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_sky
    private
    logical :: init = .false.
    integer :: nside = 0
    integer :: npix = 0
    integer :: lmax = 0
    integer :: mmax = 0
    integer :: pix_scheme = S2_SKY_RING
    real(s2_sp), allocatable :: map(:)
    complex(s2_spc), allocatable :: alm(:,:)
    logical :: map_status = .false.
    logical :: alm_status = .false.
    real(s2_sp), allocatable :: param(:)
    integer :: n_param = 0
  end type s2_sky


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_sky_init_empty
    !
    !! Initialised an empty sky at a specified resolution.  Neither the map 
    !! or alm sky representations are allocated.  They are allocated when
    !! required.
    !!
    !! Variables:
    !!  - [nside]: Healpix nside.
    !!  - [pix_scheme]: Pixelisation scheme for map.
    !!  - [lmax]: Alm lmax.
    !!  - [mmax]: Alm lmin.
    !!  - sky: The returned sky initialised with resolutions.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_empty(nside, pix_scheme, lmax, mmax) &
      result(sky)

      use pix_tools, only: nside2npix
      
      integer, intent(in), optional :: nside, pix_scheme, lmax, mmax
      type(s2_sky) :: sky

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_empty')
        return
      end if

      if(present(nside)) then
        sky%nside = nside
        sky%npix = nside2npix(nside)
      end if
      if(present(pix_scheme)) sky%pix_scheme = pix_scheme
      if(present(lmax)) then
        sky%lmax = lmax
        if(present(mmax)) then
          sky%mmax = mmax
        else
          sky%mmax = sky%lmax
        end if
      end if
      
      if(present(nside) .and. present(lmax)) then
         call s2_sky_valid_sizes(sky)
      end if

      sky%init = .true.

    end function s2_sky_init_empty


    !--------------------------------------------------------------------------
    ! s2_sky_init_map
    !
    !! Initalise a sky from a specified map.
    !!    
    !! Variables:
    !!  - map: The map array used to initialised the sky.
    !!  - nside: Healpix nside.
    !!  - pix_scheme: Pixelisation scheme for map.  The passed map is assumed 
    !!    to be in this format.
    !!  - [lmax]: Alm lmax.
    !!  - [mmax]: Alm mmax.
    !!  - sky: The returned sky initialised with the specified map.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_map(map, nside, pix_scheme, lmax, mmax) result(sky)
    
      real(s2_sp), intent(in) :: map(:)
      integer, intent(in) :: nside, pix_scheme
      integer, intent(in), optional :: lmax, mmax
      type(s2_sky) :: sky

      integer :: fail

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_map')
        return
      end if

      ! Initialise empty sky.
      sky = s2_sky_init_empty(nside, pix_scheme, lmax, mmax)

      ! Check size of passed map is consistent with other sizes passed.
      if(size(map) /= sky%npix) then 
        call s2_error(S2_ERROR_INIT_FAIL, 's2_sky_init_map', &
          comment_add='Passed map invalid size')
      end if

      ! Allocate space.
      allocate(sky%map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_map')
      end if

      ! Save map.
      sky%map = map

      ! Set status and initialised flags.
      sky%map_status = .true.
      sky%init = .true.   ! Should already be set from s2_sky_init_empty
                          ! but set regardless.

    end function s2_sky_init_map


    !--------------------------------------------------------------------------
    ! s2_sky_init_alm
    !    
    !! Initalise a sky from a specified alm.
    !!
    !! Variables:
    !!  - alm: The alm array used to initialised the sky.
    !!  - lmax: Alm lmax.
    !!  - mmax: Alm mmax.
    !!  - [nside]: Healpix nside.
    !!  - [pix_scheme]: Pixelisation scheme for map.  No map is initialised 
    !!    but if created this will be the default pixelisation scheme.
    !!  - sky: The returned sky initialised with the specified alm.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_alm(alm, lmax, mmax, nside, pix_scheme) result(sky)
    
      complex(s2_spc), intent(in) :: alm(:,:)
      integer, intent(in) :: lmax, mmax
      integer, intent(in), optional :: nside, pix_scheme
      type(s2_sky) :: sky

      integer :: fail

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_alm')
        return
      end if

      ! Initialise empty sky.
      sky = s2_sky_init_empty(nside, pix_scheme, lmax, mmax)

      ! Check size of passed alm is consistent with other sizes passed.
      if(size(alm,1) /= sky%lmax+1 .or.size(alm,2) /= sky%mmax+1 ) then 
        call s2_error(S2_ERROR_INIT_FAIL, 's2_sky_init_alm', &
          comment_add='Passed alm invalid size')
      end if

      ! Allocate space.
      allocate(sky%alm(0:sky%lmax, 0:sky%mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_alm')
      end if
      
      ! Save map.
      sky%alm = alm

      ! Set status and initialised flags.
      sky%alm_status = .true.
      sky%init = .true.   ! Should already be set from s2_sky_init_empty
                          ! but set regardless.

    end function s2_sky_init_alm


    !--------------------------------------------------------------------------
    ! s2_sky_init_ab
    !
    !! Initialise a sky from the ecp (equi-sampled) alpha-beta array defined on
    !! the sphere.
    !!
    !! Notes:
    !!   - Linear interpolation performed if interp=.true.
    !!
    !! Variables:
    !!   - x_ab: Ecp (equi-sampled) alpha-beta array defining a function on 
    !!     the sphere to construct sky from.
    !!   - interp: Logical specifying whether interpolation is to be performed.
    !!   - nside: Healpix nside.
    !!   - pix_scheme: Pixelisation scheme for map. 
    !!   - [lmax]: Alm lmax.
    !!   - [mmax]: Alm mmax.
    !!   - [sdw]: Logical to specify SDW pixelisation.
    !!   - sky: The returned sky initialised from the ecp (equi-sampled) 
    !!    alpha-beta array defined over the sphere.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_ab(x_ab, interp, nside, pix_scheme, lmax, mmax, sdw) &
      result(sky)

      use pix_tools, only: pix2ang_ring, pix2ang_nest, nside2npix

      real(s2_sp), intent(in) :: x_ab(:,:)
      logical, intent(in) :: interp
      integer, intent(in) :: nside, pix_scheme
      integer, intent(in), optional :: lmax, mmax
			logical, intent(in), optional :: sdw
      type(s2_sky) :: sky

      real(s2_dp) :: alpha, beta
      integer :: npix, ipix, fail
      real(s2_sp), allocatable :: map(:)

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_ab')
        return
      end if

      ! Set sizes and allocate memory for temporary map constructed.
      npix = nside2npix(nside)
      allocate(map(0:npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_ab')
      end if

      ! Construct map.      
      do ipix = 0,npix-1

          ! Get theta (beta) and phi (alpha) angles corresponding to pixel.
          if(pix_scheme == S2_SKY_RING) then
             call pix2ang_ring(nside, ipix, beta, alpha)
          else if(pix_scheme == S2_SKY_NEST) then
             call pix2ang_nest(nside, ipix, beta, alpha)
          else
             call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_init_ab')
          end if

          map(ipix) = s2_sky_interp_ab(x_ab, alpha, beta, interp, sdw)

      end do

      ! Initialise sky from map.
      sky = s2_sky_init_map(map, nside, pix_scheme, lmax, mmax)

      ! Free temporary map made.  (Copy is now stored in sky.)
      deallocate(map)

    end function s2_sky_init_ab


    !--------------------------------------------------------------------------
    ! s2_sky_init_fun
    !
    !! Initialise a sky from a template function.
    !!
    !! Variables:
    !!  - nside: Healpix nside.
    !!  - pix_scheme: Pixelisation scheme for map.  
    !!  - [lmax]: Alm lmax.
    !!  - [mmax]: Alm mmax.
    !!  - [param]: Parameter array specifying analytic parameters for the
    !!    function to be evaluated.
    !!  - [fun_type_in]: Integer specifing the function type (whether function 
    !!     defined directly on sphere or whether defined on plane and to be 
    !!     numericall projected onto the sphere).
    !!  - sky: The initialised sky with the template function evaluated over 
    !!    the map.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_fun(fun, nside, pix_scheme, lmax, mmax, &
      param, fun_type_in) result(sky)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      integer, intent(in) :: nside, pix_scheme
      integer, intent(in), optional :: lmax, mmax
      real(s2_sp), intent(in), optional :: param(:)
      integer, intent(in), optional :: fun_type_in
      type(s2_sky) :: sky
      interface 
        function fun(theta, phi, param) result(val)
          use s2_types_mod
          real(s2_sp), intent(in) :: theta, phi
          real(s2_sp), intent(in), optional :: param(:)
          real(s2_sp) :: val
        end function fun
      end interface

      integer :: fail, ipix
      real(s2_dp) :: theta, phi
      real(s2_sp) :: r, norm_preserving_factor
      integer :: fun_type

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_fun')
        return
      end if

      ! Initialise empty sky.
      sky = s2_sky_init_empty(nside, pix_scheme, lmax, mmax)

      ! Allocate space.
      allocate(sky%map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_fun')
      end if
      
      ! Set function type.
      if(present(fun_type_in)) then
         fun_type = fun_type_in
      else
         fun_type = S2_SKY_DEFAULT_FTYPE
      end if

      ! Compute sky.
      if(fun_type == S2_SKY_FUN_TYPE_SPHERE) then

         ! Evaluate fun on sphere for each pixel.
         do ipix = 0,sky%npix-1

            if(sky%pix_scheme == S2_SKY_RING) then
               call pix2ang_ring(nside, ipix, theta, phi)
            else if(sky%pix_scheme == S2_SKY_NEST) then
               call pix2ang_nest(nside, ipix, theta, phi)
            else
               call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_init_fun')
            end if
            
            sky%map(ipix) = fun(real(theta,s2_sp), real(phi,s2_sp), param)
            
         end do
      
      else if(fun_type == S2_SKY_FUN_TYPE_PLANE) then

         ! Evaluate fun on plane and project onto sphere for each pixel.
         do ipix = 0,sky%npix-1
            
            ! Convert pixel index to angular coordinates theta and phi.
            if(sky%pix_scheme == S2_SKY_RING) then
               call pix2ang_ring(nside, ipix, theta, phi)
            else if(sky%pix_scheme == S2_SKY_NEST) then
               call pix2ang_nest(nside, ipix, theta, phi)
            else
               call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_init_fun')
            end if
            
            r = 2 * tan(real(theta, s2_sp) / 2.0d0)
            norm_preserving_factor = 2 / (1 + cos(real(theta, s2_sp)))
            
            sky%map(ipix) = norm_preserving_factor &
              * fun(r, real(phi, s2_sp), param)
        
         end do
         
      else

         call s2_error(S2_ERROR_SKY_FTYPE_INVALID, 's2_sky_init_fun')

      end if
    
      ! Save parameters if present.
      if(present(param)) then
         sky%n_param = size(param)
         allocate(sky%param(1:sky%n_param), stat=fail)
         sky%param = param
         if(fail /= 0) then 
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_fun')
         end if
      end if

      ! Set status and initialised flags.     
      sky%map_status = .true.
      sky%init = .true.
    
    end function s2_sky_init_fun


    !--------------------------------------------------------------------------
    ! s2_sky_init_fun_alm
    !
    !! Initialise a sky from a template function defined in harmonic space.
    !!
    !! Variables:
    !!  - lmax: Alm lmax.
    !!  - mmax: Alm mmax.
    !!  - azisum: Logical specifying whether the template function is 
    !!    azimuthally symmetric, in which case non-zero harmonic 
    !!    cofficients are computed for m=0 only.
    !!  - [nside]: Healpix nside.
    !!  - [pix_scheme]: Pixelisation scheme for map.  
    !!  - [param]: Parameter array specifying analytic parameters for the
    !!    function to be evaluated.
    !!  - sky: The initialised sky with the template function evaluated over 
    !!    the map.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_fun_alm(fun, lmax, mmax, azisym, nside, pix_scheme, &
      param) result(sky)

      integer, intent(in) :: lmax, mmax
      logical, intent(in) :: azisym
      integer, intent(in), optional :: nside, pix_scheme
      real(s2_sp), intent(in), optional :: param(:)
      type(s2_sky) :: sky
      interface 
        function fun(el, m, param) result(val)
          use s2_types_mod
          integer, intent(in) :: el, m
          real(s2_sp), intent(in), optional :: param(:)
          complex(s2_spc) :: val
        end function fun
      end interface

      integer :: fail = 0
      integer :: el, m

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_fun_alm')
        return
      end if

      ! Initialise empty sky.
      sky = s2_sky_init_empty(nside, pix_scheme, lmax, mmax)

      ! Allocate space.
      allocate(sky%alm(0:sky%lmax,0:sky%mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_fun_alm')
      end if
      sky%alm(0:sky%lmax,0:sky%mmax) = cmplx(0e0, 0e0)

      ! Compute alms.
      do el = 0,lmax
         if (azisym) then
            m = 0
            sky%alm(el,m) = fun(el, m, param)
         else
            do m = 0,min(el,mmax)
               sky%alm(el,m) = fun(el, m, param)
            end do
         end if
      end do

      ! Save parameters if present.
      if(present(param)) then
         sky%n_param = size(param)
         allocate(sky%param(1:sky%n_param), stat=fail)
         sky%param = param
         if(fail /= 0) then 
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_fun_alm')
         end if
      end if

      ! Set status and initialised flags.     
      sky%alm_status = .true.
      sky%init = .true.

    end function s2_sky_init_fun_alm


    !--------------------------------------------------------------------------
    ! s2_sky_init_file
    ! 
    !! Wrapper to initialise a sky from either a map fits file or a full
    !! s2_sky fits file.
    !!
    !! Variables:
    !!  - filename: The fits filename to read the sky from.
    !!  - file_type: Fit type specifier to specify whether to read from a 
    !!    fits map file or from a fits full s2_sky file.
    !!  - [extension]: The fits extension the map to be read is stored in, 
    !!    If not specified the default is 1. 
    !!    (If file_type == S2_SKY_FILE_TYPE_MAP.)
    !!  - sky: The initialised sky with the map, nside, etc. defined from the
    !!    fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !   November 2005 - Modified by Jason McEwen to incorporate alm fits file.
    !--------------------------------------------------------------------------

    function s2_sky_init_file(filename, file_type, extension) result(sky)

      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_type
      integer, intent(in), optional :: extension
      type(s2_sky) :: sky
    
      select case(file_type)

         case(S2_SKY_FILE_TYPE_SKY)
            sky = s2_sky_init_file_sky(filename)

         case(S2_SKY_FILE_TYPE_MAP)
            sky = s2_sky_init_file_map(filename, extension)

         case(S2_SKY_FILE_TYPE_ALM)
            sky = s2_sky_init_file_alm(filename)

         case default
            call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky_init_file', &
              comment_add='Invalid file type specifier')

      end select

    end function s2_sky_init_file


    !--------------------------------------------------------------------------
    ! s2_sky_init_file_map
    ! 
    !! Initialise a sky from a map file.
    !!
    !! Variables:
    !!  - filename: The fits filename to read the map from.
    !!  - [extension_in]: The fits extension the map to be read is stored in, 
    !!    If not specified the default is 1.
    !!  - sky: The initialised sky with the map, nside, etc. defined from the
    !!    map contained in the fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_file_map(filename, extension_in) result(sky)
    
      use fitstools, only: getsize_fits, input_map

      character(len=*), intent(in) :: filename
      integer, intent(in), optional :: extension_in
      type(s2_sky) :: sky
    
      integer :: nmaps, polarisation, extension
      integer :: fail
      real(s2_sp), allocatable :: map_temp(:,:)

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_map_file')
        return
      end if

      ! Get sizes from fits file.
      sky%npix = getsize_fits(filename, nmaps, sky%pix_scheme, &
        nside=sky%nside, polarisation=polarisation)

      ! Pixel scheme unknown.
      if(sky%pix_scheme /= S2_SKY_RING &
          .and. sky%pix_scheme /= S2_SKY_NEST) then
        call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_init_file_map')
      end if

      ! Functionality to support a polarised map not yet incorporated.
      if(polarisation /= 0) then
        call s2_error(S2_ERROR_SKY_POL_DEF, 's2_sky_init_map_file')
      end if

      ! Check and set extension status.
      if(present(extension_in)) then
        extension = extension_in
        if(extension > nmaps) then
          call s2_error(S2_ERROR_SKY_EXT_INVALID, 's2_sky_init_map_file')
        end if
      else
         ! If extension not passed then set to 1.
         extension = 1
      end if
 
      ! Allocate space.
      allocate(map_temp(0:sky%npix, nmaps), stat=fail)
      allocate(sky%map(0:sky%npix), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_map_file')
      end if

      ! Load all maps from fits file.
      call input_map(filename, map_temp, sky%npix, nmaps)

      ! Copy map data to sky (only extension specified).
      sky%map = map_temp(:,extension)

      ! Free local dynamic storage space.
      deallocate(map_temp)

      ! Set status and initialised flags.      
      sky%map_status = .true.
      sky%init = .true.

    end function s2_sky_init_file_map


    !--------------------------------------------------------------------------
    ! s2_sky_init_file_alm
    ! 
    !! Initialise a sky from a fits alm file.
    !!
    !! Variables:
    !!  - filename: The fits filename to read the alms from.
    !!  - sky: The initialised sky with the alms, lmax and mmax defined from 
    !!    the alms contained in the fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2005
    !
    ! Revisions:
    !   November 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_file_alm(filename) result(sky)

      use fitstools, only: number_of_alms, fits2alms

      character(len=*), intent(in) :: filename
      type(s2_sky) :: sky

      integer :: nalm, ncol, next, fail
      integer :: lmax, mmax
      integer :: l, m, i
      real(s2_sp), allocatable :: alm_input(:,:,:)
      integer, parameter :: HEADER_LEN = 180
      character(len=80) :: header(HEADER_LEN,1)

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_file_alm')
        return
      end if

      ! Set/get sizes.
      nalm = number_of_alms(filename)
      ncol = 3
      next = 1

      ! Allocate space for input alm array.
      allocate(alm_input(1:nalm, 1:(ncol+1), 1:next), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_file_alm')
      end if

      ! Get alms from file.
      call fits2alms(filename, nalm, alm_input, ncol, &
        header, HEADER_LEN, next)

      ! Set lmax and mmax from maximum values contained in alm_input
      ! l and m columns.
      lmax = maxval(alm_input(:,1,1))
      mmax = maxval(alm_input(:,2,1))

      ! Alternatively could use last l and m.
      ! lmax = alm_input(nalm,1,1)
      ! mmax = alm_input(nalm,2,1)

      ! Alternatively could use formula to get l from lm but would
      ! have to get lm_max from file (differs to nalm)

      ! Allocate space for alms.
      allocate(sky%alm(0:lmax, 0:mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_file_alm')
      end if
      sky%alm = 0e0

      ! Copy input alm array values to sky object alm array.
      i = 1
      do l = 0,lmax
         do m = 0,min(l,mmax)
            sky%alm(l,m) = cmplx(alm_input(i,3,1), alm_input(i,4,1))
            i = i + 1
         end do
      end do
      i = i-1

      ! Initialise other sky attributes.
      sky%lmax = lmax
      sky%mmax = mmax
      sky%alm_status = .true.
      sky%init = .true. 

      ! Free temporary memory.
      deallocate(alm_input)

    end function s2_sky_init_file_alm


    !--------------------------------------------------------------------------
    ! s2_sky_init_file_sky
    !
    !! Wrapper to initialise a sky data structure from a full s2_sky file.
    !! The sky structure is read and initialised by the routine 
    !! s2_sky_io_fits_read.
    !!
    !! Notes:
    !!   - Cannot be part of s2_sky_init interface since has same signature
    !!     as s2_sky_init_map_file.
    !!
    !! Variables:
    !!   - filename: Name of s2_sky fits file containing the sky data to 
    !!     be read.
    !!   - sky: Returned sky structure initialised with the data contained in
    !!     the input s2_sky fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_file_sky(filename) result(sky)

      character(len=*), intent(in) :: filename
      type(s2_sky) :: sky

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_read')
        return
      end if

      ! Read s2_sky file.
      call s2_sky_io_fits_read(filename, sky)

      ! All status flags set in s2_sky_io_fits_read routine.
      
    end function s2_sky_init_file_sky

  
    !--------------------------------------------------------------------------
    ! s2_sky_init_copy
    !
    !! Initialise a sky as a copy of another sky.
    !!
    !! Variables:
    !!   - orig: The original sky to copy.
    !!   - copy: The initialised map copied from the original.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_init_copy(orig) result(copy)

      type(s2_sky), intent(in) :: orig
      type(s2_sky) :: copy

      integer :: fail

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_init_copy')
        return
      end if

      ! Copy object attributes.

      copy%nside = orig%nside
      copy%npix = orig%npix
      copy%lmax = orig%lmax
      copy%mmax = orig%mmax
      copy%pix_scheme = orig%pix_scheme
      copy%map_status = orig%map_status
      copy%alm_status = orig%alm_status
      copy%n_param = orig%n_param

      if(copy%map_status) then
         allocate(copy%map(0:copy%npix-1), stat=fail)
         if(fail /= 0) then
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_copy')
         end if
         copy%map(0:copy%npix-1) = orig%map(0:copy%npix-1)
      end if

      if(copy%alm_status) then
         allocate(copy%alm(0:copy%lmax, 0:copy%mmax), stat=fail)
         if(fail /= 0) then
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_copy')
         end if
         copy%alm(0:copy%lmax, 0:copy%mmax) = orig%alm(0:copy%lmax, 0:copy%mmax)
      end if

      if(copy%n_param /= 0) then
         allocate(copy%param(1:copy%n_param), stat=fail)
         if(fail /= 0) then
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_copy')
         end if
         copy%param(1:copy%n_param) = orig%param(1:copy%n_param)
      end if

      copy%init = .true.

    end function s2_sky_init_copy


    !--------------------------------------------------------------------------
    ! s2_sky_free
    !
    !! Free all data associated with an initialised sky and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - sky: The sky to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_free(sky)

      type(s2_sky), intent(inout) :: sky

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_free')
      end if 

      ! Free space.
      if(allocated(sky%map)) deallocate(sky%map)
      if(allocated(sky%alm)) deallocate(sky%alm)
      if(allocated(sky%param)) deallocate(sky%param)
 
      ! Reset attributes.
      sky%nside = 0
      sky%npix = 0
      sky%lmax = 0
      sky%mmax = 0
      sky%map_status = .false.
      sky%alm_status = .false.
      sky%n_param = 0

      sky%init = .false.

    end subroutine s2_sky_free


 		!--------------------------------------------------------------------------
    ! s2_sky_remove_map
    !
    !! Free all data associated with a map and reset all map related 
    !! attributes.
    !!
    !! Variables:
    !!   - sky: The sky containing the map to be removed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2007
    !
    ! Revisions:
    !   November 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

 		subroutine s2_sky_remove_map(sky)

      type(s2_sky), intent(inout) :: sky

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_remove_map')
      end if 

      ! Free space.
      if(allocated(sky%map)) deallocate(sky%map)
 
      ! Reset attributes.
      sky%nside = 0
      sky%npix = 0
      sky%map_status = .false.

    end subroutine s2_sky_remove_map


    !--------------------------------------------------------------------------
    ! s2_sky_compute_alm
    !
    !! Compute the alms for a sky from the map representation.  Errors occur
    !! if the initialised sky does not have a map defined.
    !!
    !! Variables:
    !!   - sky: The sky to compute the alms for.
    !!   - [lmax]: If specified overwrites the lmax of the sky.
    !!   - [mmax]: If specified overwrites the mmax of the sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_compute_alm(sky, lmax, mmax, message)

      use alm_tools, only: map2alm

      type(s2_sky), intent(inout) :: sky
      integer, intent(in), optional :: lmax, mmax
      logical, intent(in), optional :: message

      integer :: fail
      logical :: message_use
      complex(s2_spc), allocatable :: alm_temp(:,:,:)
      real(s2_dp), allocatable :: w8ring(:,:)
      !real(s2_dp) :: cos_theta_cut = -1.0d0
      real(s2_dp) :: zbounds(1:2)

      message_use = .true.

      ! Set message status.
      if(present(message)) message_use = message

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_compute_alm')
      end if 

      ! Redefine sky lmax and mmax if present and different to before.
      if(present(lmax)) then
         call s2_sky_set_lmax(sky, lmax, mmax)
      end if

      ! Check if alm already defined.
      if(sky%alm_status) then
         ! If alm already calculated at current lmax and mmax then done.
         ! (if lmax and mmax changed and map was previously defined it will 
         ! have been removed by s2_sky_set_lmax).
         call s2_error(S2_ERROR_SKY_ALM_DEF, 's2_sky_compute_alm')
         return
      end if
      
      ! Check that lmax and mmax are defined (i.e. not 0) and sizes valid.
      if(sky%lmax == 0 .or. sky%mmax == 0) then
         call s2_error(S2_ERROR_SKY_SIZE_NOT_DEF, 's2_sky_compute_alm')
      end if
      call s2_sky_valid_sizes(sky)

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_compute_alm')
      end if

      ! Check map in ring pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Calculate alm if passed all previous tests...
      
      ! Allocate space (one polarisation only).
      allocate(alm_temp(1, 0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(sky%alm(0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(w8ring(1:2*sky%nside, 1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_compute_alm')
      end if

      ! Compute alm.
      if(message_use) then
          write(*,'(a,a,i4,a,i4)') 's2_sky_compute_alm> ', &
               'Computing alm from map at nside =', &
               sky%nside, ', lmax =', sky%lmax
      end if

      w8ring = 1.0d0
      zbounds(1) = -1.0d0
      zbounds(2) = 1.0d0

      call map2alm(sky%nside, sky%lmax, sky%mmax, sky%map, &
        alm_temp, zbounds, w8ring)

      ! Copy alm data to sky (one polarisation only).
      sky%alm(0:sky%lmax, 0:sky%mmax) = alm_temp(1,0:sky%lmax,0:sky%mmax)

      ! Free local dynamic storage space.
      deallocate(alm_temp)

      ! Set alm_status.
      sky%alm_status = .true.

    end subroutine s2_sky_compute_alm


    !--------------------------------------------------------------------------
    ! s2_sky_compute_alm_iter
    !
    !! Compute the alms for a sky from the map representation using iteration.
    !! Errors occur if the initialised sky does not have a map defined.
    !!
    !! Variables:
    !!   - sky: The sky to compute the alms for.
    !!   - iter_order: Number of iterations to perform (iter_order=0 gives 
    !!     same results as s2_sky_compute_alm
    !!   - [lmax]: If specified overwrites the lmax of the sky.
    !!   - [mmax]: If specified overwrites the mmax of the sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2007
    !
    ! Revisions:
    !   November 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_compute_alm_iter(sky, iter_order, lmax, mmax, message)

      use alm_tools, only: map2alm, alm2map

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: iter_order
      integer, intent(in), optional :: lmax, mmax
      logical, intent(in), optional :: message

      integer :: fail, iter
      logical :: message_use
      complex(s2_spc), allocatable :: alm_temp(:,:,:)
      complex(s2_spc), allocatable :: alm_temp2(:,:,:)
      real(s2_sp), allocatable :: map_temp(:)
      real(s2_dp), allocatable :: w8ring(:,:)
      real(s2_dp) :: zbounds(1:2)
      !real(s2_dp) :: cos_theta_cut = -1.0d0

      message_use = .true.

      ! Set message status.
      if(present(message)) message_use = message

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_compute_alm_iter')
      end if 

      ! Redefine sky lmax and mmax if present and different to before.
      if(present(lmax)) then
         call s2_sky_set_lmax(sky, lmax, mmax)
      end if

      ! Check if alm already defined.
      if(sky%alm_status) then
         ! If alm already calculated at current lmax and mmax then done.
         ! (if lmax and mmax changed and map was previously defined it will 
         ! have been removed by s2_sky_set_lmax).
         call s2_error(S2_ERROR_SKY_ALM_DEF, 's2_sky_compute_alm_iter')
         return
      end if
      
      ! Check that lmax and mmax are defined (i.e. not 0) and sizes valid.
      if(sky%lmax == 0 .or. sky%mmax == 0) then
         call s2_error(S2_ERROR_SKY_SIZE_NOT_DEF, 's2_sky_compute_alm_iter')
      end if
      call s2_sky_valid_sizes(sky)

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_compute_alm_iter')
      end if

      ! Check map in ring pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Calculate alm if passed all previous tests...
      
      ! Allocate space (one polarisation only).
      allocate(alm_temp(1, 0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(alm_temp2(1, 0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(map_temp(0:sky%npix-1), stat=fail)
      allocate(sky%alm(0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(w8ring(1:2*sky%nside, 1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_compute_alm_iter')
      end if

      ! Compute alm.

      if(message_use) then
          write(*,'(a,a,i4,a,i4)') 's2_sky_compute_alm_iter> ', &
               'Computing alm from map at nside =', &
               sky%nside, ', lmax =', sky%lmax
      end if

      w8ring = 1.0d0
      zbounds(1) = -1d0 
      zbounds(2) = 1d0

      alm_temp(1, 0:sky%lmax, 0:sky%mmax) = 0e0
      alm_temp2(1, 0:sky%lmax, 0:sky%mmax) = 0e0
      map_temp(0:sky%npix-1) = sky%map(0:sky%npix-1)

      do iter = 0, iter_order

				!write(*,*) 'iter=', iter

        call map2alm(sky%nside, sky%lmax, sky%mmax, map_temp, &
          alm_temp, zbounds, w8ring)

       alm_temp2(1, 0:sky%lmax, 0:sky%mmax) = alm_temp2(1, 0:sky%lmax, 0:sky%mmax) &
				 + alm_temp(1, 0:sky%lmax, 0:sky%mmax)

        if (iter == iter_order) exit

        ! alm_2 -> map : reconstructed map
        call alm2map(sky%nside, sky%lmax, sky%mmax, alm_temp2, map_temp)

        ! map = map_2 - map : residual map
        map_temp(0:sky%npix-1) = sky%map(0:sky%npix-1) - map_temp(0:sky%npix-1)

     end do

      ! Copy alm data to sky (one polarisation only).
      sky%alm = alm_temp2(1,:,:)

      ! Free local dynamic storage space.
      deallocate(alm_temp, alm_temp2, map_temp)
      deallocate(w8ring)

      ! Set alm_status.
      sky%alm_status = .true.

    end subroutine s2_sky_compute_alm_iter


    !--------------------------------------------------------------------------
    ! s2_sky_compute_map
    ! 
    !! Compute the map for a sky from the alm representation.  Errors occur
    !! if the initialised sky does not have an alm defined.
    !!
    !! Variables:
    !!   - sky: The sky to compute the map for.
    !!   - [nside]: If specified overwrites the nside of the sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_compute_map(sky, nside, message)

      use alm_tools, only: alm2map

      type(s2_sky), intent(inout) :: sky
      integer, intent(in), optional :: nside
      logical, intent(in), optional :: message

      integer :: fail  
      logical :: message_use
      complex(s2_spc), allocatable :: alm_temp(:,:,:)

      ! Set message status.
      message_use = .true.
      if(present(message)) message_use = message

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_compute_map')
      end if 

      ! Redefine sky nside if present.
      if(present(nside)) then
         call s2_sky_set_nside(sky, nside)
      end if

      ! Check if map already defined.
      if(sky%map_status) then
         ! If map already calculated at current nside then done
         ! (if nside changed and map was previously defined it will have 
         ! been removed by s2_sky_set_nside).
         call s2_error(S2_ERROR_SKY_MAP_DEF, 's2_sky_compute_map')
         return
      end if

      ! Check that nside defined (i.e. not 0) and sizes valid.
      if(sky%nside == 0) then
         call s2_error(S2_ERROR_SKY_SIZE_NOT_DEF, 's2_sky_compute_map')
      end if
      call s2_sky_valid_sizes(sky)

      ! Check alm defined
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_compute_map')
      end if

      ! Calculate map if passed all previous tests...

      ! Allocate space (one polarisation only).
      allocate(alm_temp(1, 0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(sky%map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_compute_map')
      end if

      ! Compute map.
      alm_temp(1,:,:) = sky%alm    ! Copy sky to alm data (one pol only).

      if(message_use) then
         write(*,'(a,a,i4,a,i4)') 's2_sky_compute_map> ', &
              'Computing map from alms at nside =', &
              sky%nside, ', lmax =', sky%lmax
      end if
         
      call alm2map(sky%nside, sky%lmax, sky%mmax, alm_temp, sky%map)
      
      ! Set pix_scheme.
      sky%pix_scheme = S2_SKY_RING

      ! Free local dynamic storage space.
      deallocate(alm_temp)

      ! Set map_status.
      sky%map_status = .true.

    end subroutine s2_sky_compute_map


    !--------------------------------------------------------------------------
    ! s2_sky_irregular_invsht
    ! 
    !! Compute real space values of a sky, on an irregular grid, from its
    !! spherical harmonic coefficients.  Errors occur
    !! if the initialised sky does not have an alm already defined.
    !!
    !! Notes:
    !!   - Implementation is based on separation of vairables, however true 
    !!     separation of variables is not possible since irregular (theta,phi)
    !!     grid may not be separable (for example, every point inthe grid may 
    !!     have different (theta,phi)).
    !!
    !! Variables:
    !!   - sky: The sky to compute the map for.
    !!   - f(0:N-1): The output function computed, defined over the irregular 
    !!     grid.
    !!   - N: The number of points on the irregular grid.
    !!   - thetas(0:N-1): The theta positions of the irregular grid.
    !!   - phis(0:N-1): The phi positions of the irregular grid.
    !!   - [azisym]: Logical specifying whether the function is
    !!     azimuthally symmetric, in which case summations over m are
    !!     truncated to m=0 only.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_irregular_invsht(sky, f, N, thetas, phis, azisym)

      use s2_dl_mod, only: s2_dl_beta_operator

      type(s2_sky), intent(in) :: sky
      integer, intent(in) :: N
      real(s2_dp), intent(in) :: thetas(0:N-1)
      real(s2_dp), intent(in) :: phis(0:N-1)
      real(s2_dp), intent(out) :: f(0:N-1)
      logical, intent(in), optional :: azisym

      complex(s2_dpc), allocatable :: fmtheta(:,:)
      real(s2_dp), pointer :: dl(:,:)
      integer :: lmax, mmax, itheta, el, m
      integer :: fail = 0
      complex(s2_dpc) :: I
      
      I = cmplx(0d0, 1d0)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_irregular_invsht')
      end if

      ! Check alms defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_irregular_invsht')
         return
      end if

      ! Allocate temporary space.
      lmax = sky%lmax
      mmax = sky%mmax
      if(present(azisym)) then 
         if(azisym) then
            mmax = 0
         end if
      end if
      allocate(fmtheta(0:mmax,0:N-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_irregular_invsht')
      end if
      fmtheta(0:mmax,0:N-1) = 0e0

      ! Compute f_m(theta).
      do itheta = 0,N-1
         do el = 0,lmax
            ! Compute dlmns.
            allocate(dl(-el:el,-el:el), stat=fail)
            if(fail /= 0) then
               call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_irregular_invsht')
            end if
            call s2_dl_beta_operator(dl, thetas(itheta), el)

            do m = 0,min(el,mmax)
               fmtheta(m, itheta) = fmtheta(m, itheta)  + &
                    sqrt((2.0*el+1.0)/(4.0*PI)) * dl(m,0) * sky%alm(el,m)
            end do

            ! Free dlmns.
            deallocate(dl)
         end do

      end do

      ! Compute f(theta,phi).
      f(0:N-1) = 0e0
      do itheta = 0,N-1
         f(itheta) =  fmtheta(0,itheta)
         do m = 1,mmax
            f(itheta) = f(itheta) &
                 + 2 * real(exp(I*m*phis(itheta)) * fmtheta(m,itheta),s2_dp)
         end do
      end do

      ! Free memory.
      deallocate(fmtheta)

    end subroutine s2_sky_irregular_invsht


    !--------------------------------------------------------------------------
    ! s2_sky_map_convert
    !
    !! Convert the map pixelisation scheme to that specified.  If the map is
    !! already in the required pixelisation scheme then do nothing.
    !!
    !! Variables:
    !!   - sky: Sky to convert pixelisation scheme.
    !!   - pix_scheme: Required pixelisation scheme.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_map_convert(sky, pix_scheme)

      use pix_tools, only: convert_ring2nest, convert_nest2ring

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: pix_scheme

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_map_convert')
      end if 
  
      ! Check not already in required pixelisation scheme.
      if(pix_scheme == sky%pix_scheme) then
         ! Actually don't print warning message
         ! write(*,'(a)') &
         ! 's2_sky_map_convert> Map already in required pixelisation scheme'
         return
      end if

      sky%pix_scheme = pix_scheme
      if(sky%pix_scheme == S2_SKY_RING) then
         call convert_nest2ring(sky%nside, sky%map(0:sky%npix-1))
      else if(sky%pix_scheme == S2_SKY_NEST) then
         call convert_ring2nest(sky%nside, sky%map(0:sky%npix-1))
      else
         call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_map_convert')
      end if

    end subroutine s2_sky_map_convert


    !--------------------------------------------------------------------------
    ! s2_sky_der
    !
    !! Compute continuous derivative on sphere.
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!
    !! Variables:
    !!   - sky: Sky to compute derivate of.
    !!   - der_type: Type of derivative to return 
    !!     (S2_SKY_DER_TYPE_THETA: der = dT/dtheta; 
    !!     S2_SKY_DER_TYPE_PHI: der = 1/sin(theta) * dT/dphi;
    !!     S2_SKY_DER_TYPE_GRAD: der = sqrt( (dT/dtheta)**2 + (1/sin(theta) * dT/dphi)**2) )
    !!   - lmax: Alm lmax (if harmonic coefficients not already computed).
    !!   - mmax: Alm mmax (if harmonic coefficients not already computed).
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_der(sky, der_type, lmax, mmax) result(der)

      use alm_tools, only: alm2map_der

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: der_type
      integer, intent(in), optional :: lmax, mmax     
      type(s2_sky) :: der

      complex(s2_spc), allocatable :: alm_temp(:,:,:)
      real(s2_sp), allocatable :: der1(:,:)
      real(s2_sp), allocatable :: map(:)
      integer :: fail = 0

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der')
      end if 

      ! Compute alm if not already computed.
      if(.not. sky%alm_status) call s2_sky_compute_alm(sky, lmax, mmax)

      ! Allocate space.
      allocate(alm_temp(1, 0:sky%lmax, 0:sky%mmax), stat=fail)
      allocate(der1(0:sky%npix-1, 1:2), stat=fail)
      allocate(map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der')
      end if

      ! Compute first derivatives.
      alm_temp(1,0:sky%lmax, 0:sky%mmax) = sky%alm(0:sky%lmax, 0:sky%mmax)
      call alm2map_der(sky%nside, sky%lmax, sky%mmax, alm_temp, map, der1)
      
      ! Discard map; overwrite with desired derivative.
      select case(der_type)
          
          case(S2_SKY_DER_TYPE_THETA)
             map(0:sky%npix-1) = der1(0:sky%npix-1, 1)

          case(S2_SKY_DER_TYPE_PHI)
             map(0:sky%npix-1) = der1(0:sky%npix-1, 2)

          case(S2_SKY_DER_TYPE_GRAD)
             map(0:sky%npix-1) = &
                  sqrt(der1(0:sky%npix-1, 1)**2 + der1(0:sky%npix-1, 2)**2)

          case default
             call s2_error(S2_ERROR_SKY_DER_TYPE_INVALID, 's2_sky_der')

      end select
     
      ! Initialise output sky with new map.
      der = s2_sky_init_map(map, sky%nside, S2_SKY_RING, sky%lmax, sky%mmax)

      ! Free memory.
      deallocate(alm_temp)
      deallocate(der1)
      deallocate(map)

    end function s2_sky_der


    !--------------------------------------------------------------------------
    ! s2_sky_der_discrete_phi
    !
    !! Compute discrete derivative on sphere with respect to phi.
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!
    !! Variables:
    !!   - sky: Sky to compute derivate of.
    !!   - divbysin: Logical specifying whether to divide by sin(theta) 
    !!     ( if true, return 1/sin(theta) dT/dphi; else return dT/dphi).
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_der_discrete_phi(sky, divbysin) result (der)

      use pix_tools, only: pix2ang_ring, ring2nest, nest2ring, next_in_line_nest

      type(s2_sky), intent(inout) :: sky
      logical, intent(in) :: divbysin
      type(s2_sky) :: der

      real(s2_sp), allocatable :: map(:)
      integer :: ipix, ipix_adj, ipix_nest, ipix_adj_nest, fail=0
      real(s2_dp) :: theta, phi, theta_adj, phi_adj, dot, dphi
      real(s2_dp), parameter :: TOL = 1d-8

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der_discrete_phi')
      end if 

      ! Check map computed.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_der_discrete_phi')
      end if
      
      ! Ensure in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Allocate space for temporary map.
      allocate(map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_phi')
      end if

      ! Compute derivative.
      do ipix = 0,sky%npix-1
         
         ! Find adjacent pixel on ring.
         call ring2nest(sky%nside, ipix, ipix_nest)
         call next_in_line_nest(sky%nside, ipix_nest, ipix_adj_nest)
         call nest2ring(sky%nside, ipix_adj_nest, ipix_adj)

         ! Compute finite difference.
         map(ipix) = sky%map(ipix_adj) - sky%map(ipix)

         ! Divide by dphi.
         call pix2ang_ring(sky%nside, ipix, theta, phi)
         call pix2ang_ring(sky%nside, ipix_adj, theta_adj, phi_adj)
         dot = cos(phi)*cos(phi_adj) + sin(phi)*sin(phi_adj)
         if (dot > 1.0) dot = 1.0
         if (dot < -1.0) dot = -1.0 
         dphi = acos(dot)
         if(abs(dphi) < TOL) then
            call s2_error(S2_ERROR_ARTH, 's2_sky_der_discrete_phi')
         end if
         map(ipix) = map(ipix) / dphi
        
         ! Divide by sin(theta).
         if(divbysin) then
            map(ipix) = map(ipix) / sin(theta)
         end if

      end do

      ! Initialise output sky with new map.
      der = s2_sky_init_map(map, sky%nside, S2_SKY_RING)

      ! Free memory.
      deallocate(map)
      
    end function s2_sky_der_discrete_phi


    !--------------------------------------------------------------------------
    ! s2_sky_der_discrete_phi_fovop
    !
    !! Compute a sparse matrix representation of the discrete derivative 
    !! operator on sphere with respect to phi.
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!
    !! Variables:
    !!   - sky: Sky to extract polar cap.
    !!   - divbysin: Logical specifying whether to divide by sin(theta) 
    !!     ( if true, return 1/sin(theta) dT/dphi; else return dT/dphi).
    !!   - theta_fov: Size of field of view (fov), note that fov extends
    !!     from 0 to theta_fov/2 .
    !!   - nop: Number of non-zero enties in the derivative operator.
    !!   - op(0:nop-1,0:2): Discrete phi derivative operator (specifies 
    !!     indices and values of entries in operator matrix).
    !!   - nsphere: Number of pixels on the sphere within the fov.
    !!   - xmap(0:nsphere-1): Vector of pixels values on the sphere within 
    !!     the fov.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_der_discrete_phi_fovop(sky, divbysin, &
         theta_fov, nop, op, nsphere, xmap)

      use pix_tools, only: pix2ang_ring, ring2nest, nest2ring, next_in_line_nest

      type(s2_sky), intent(inout) :: sky 
      logical, intent(in) :: divbysin
      real(s2_dp), intent(in) :: theta_fov
      integer, intent(out) :: nop
      real(s2_dp), allocatable, intent(out) :: op(:,:)
      integer, intent(out) :: nsphere
      real(s2_sp), allocatable, intent(out) :: xmap(:)

      integer :: ipix, iop
      integer :: ipix_adj, ipix_nest, ipix_adj_nest
      integer :: fail = 0
      real(s2_dp) :: theta, phi, scale
      real(s2_dp) :: theta_adj, phi_adj, dot, dphi
      real(s2_dp), allocatable :: opx(:,:)
      real(s2_dp), parameter :: TOL = 1d-8

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der_discrete_phi_fovop')
      end if 

      ! Check map computed.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_der_discrete_phi_fovop')
      end if

      ! Ensure map in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Count number of pixels within FOV.
      nsphere = 0
      do ipix = 0,sky%npix-1 
         call pix2ang_ring(sky%nside, ipix, theta, phi)
         if(theta > theta_fov/2.0) then
            nsphere = ipix
            exit            
         end if
      end do

      ! Get xmap vector.
      allocate(xmap(0:nsphere-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_phi_fovop')
      end if
      xmap(0:nsphere-1) = sky%map(0:nsphere-1)

      ! Allocate maximum required space for sparse representation of operator.
      allocate(opx(0:nsphere*nsphere-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_phi_fovop')
      end if  

      ! Compute sparse matrix representation of convolution operator.
      iop = 0
      do ipix = 0,nsphere-1

         ! Get theta and phi angles corresponding to pixel.
         call pix2ang_ring(sky%nside, ipix, theta, phi)
        
         if(theta <= theta_fov/2.0) then

            ! Find adjacent pixel on ring.
            call ring2nest(sky%nside, ipix, ipix_nest)
            call next_in_line_nest(sky%nside, ipix_nest, ipix_adj_nest)
            call nest2ring(sky%nside, ipix_adj_nest, ipix_adj)
            
            ! Compute dphi scaling.
            call pix2ang_ring(sky%nside, ipix, theta, phi)
            call pix2ang_ring(sky%nside, ipix_adj, theta_adj, phi_adj)
            dot = cos(phi)*cos(phi_adj) + sin(phi)*sin(phi_adj)
            if (dot > 1.0) dot = 1.0
            if (dot < -1.0) dot = -1.0 
            dphi = acos(dot)
            if(abs(dphi) < TOL) then
               call s2_error(S2_ERROR_ARTH, 's2_sky_der_discrete_phi')
            end if
            scale = 1d0 / dphi

            ! Compute optional sin(theta) scaling.
            if(divbysin) then
               scale = scale / sin(theta)
            end if

            ! Set entries of sparse matrix...

            opx(iop, 0) = ipix
            opx(iop, 1) = ipix_adj
            opx(iop, 2) = scale
            iop = iop + 1

            opx(iop, 0) = ipix
            opx(iop, 1) = ipix
            opx(iop, 2) = - scale
            iop = iop + 1

         end if

      end do
      nop = iop

      ! Copy required operator points.
      allocate(op(0:nop-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_phi_fovop')
      end if   
      op(0:nop-1, 0:2) = opx(0:nop-1, 0:2)

      ! Free memory.
      deallocate(opx)

    end subroutine s2_sky_der_discrete_phi_fovop


    !--------------------------------------------------------------------------
    ! s2_sky_der_discrete_theta
    !
    !! Compute discrete derivative on sphere with respect to theta (using 
    !! convolutional or nearest interpolation).
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!
    !! Variables:
    !!   - sky: Sky to compute derivate of.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !!   - nearest: Logical to specify nearest interpolation rather than
    !!     convolutional.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_der_discrete_theta(sky, support_theta, kernel, param, &
         inclusive, nearest) result (der)

      use pix_tools, only: pix2ang_ring, ring_num, ring2z, ang2pix_ring

      type(s2_sky), intent(inout) :: sky
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      logical, intent(in), optional :: nearest
      type(s2_sky) :: der
      interface 
         function kernel(theta, param) result(val)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta
           real(s2_dp), intent(in), optional :: param(:)
           real(s2_dp) :: val
         end function kernel
      end interface

      integer :: fail = 0
      integer :: ipix, iring
      real(s2_dp) :: theta, theta_adj, phi
      real(s2_dp) :: z, zadj
      real(s2_dp) :: val_adj, val, dtheta
      real(s2_sp), allocatable :: map(:)
      integer :: ipix_adj
      logical :: nearest_use
      real(s2_dp), parameter :: TOL = 1d-8

      if(present(nearest)) then
         nearest_use = nearest
      else
         nearest_use = .false.
      end if

      ! Check object initialised.
      if(.not. sky%init) then
         call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der_discrete_theta')
      end if

      ! Check map computed.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_der_discrete_theta')
      end if

      ! Ensure in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Allocate space for temporary map.
      allocate(map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_theta')
      end if

      ! Compute derivative.
      do ipix = 0,sky%npix-1

         ! Compute ring number.
         call pix2ang_ring(sky%nside, ipix, theta, phi)               
         z = cos(theta)
         iring = ring_num(sky%nside, z)

         ! If final ring set derivative to zero.
         ! (Note iring ranges from [1, 4*Nside-1] inclusive.)
         if (iring == 4*sky%nside-1) then

            map(ipix) = 0.0

         else

            ! Compute theta for next ring.
            zadj = ring2z(sky%nside, iring+1)
            ! Eliminate numerical noise that could otherwise cause instability.
            if (zadj > 1.0) zadj = 1.0
            if (zadj < -1.0) zadj = -1.0
            theta_adj = acos(zadj)
            dtheta = theta_adj - theta
            if(abs(dtheta) < TOL) then
               call s2_error(S2_ERROR_ARTH, 's2_sky_der_discrete_theta')
            end if

            if (nearest_use) then

               ! Find nearest pixel to use in computing finite difference.
               call ang2pix_ring(sky%nside, theta_adj, phi, ipix_adj)
               map(ipix) = sky%map(ipix_adj) - sky%map(ipix)

            else

               ! Compute map value at (theta_adj, phi) by convolutional interpolation.
               val_adj = s2_sky_convpt_space(sky, support_theta, kernel, &
                    theta_adj, phi, param, inclusive)
               
               ! Compute map value at (theta, phi) by convolutional interpolation.
               val = s2_sky_convpt_space(sky, support_theta, kernel, &
                    theta, phi, param, inclusive)

               ! Compute finite difference.
               map(ipix) = (val_adj - val) / dtheta

            end if

         end if

      end do

      ! Initialise output sky with new map.
      der = s2_sky_init_map(map, sky%nside, S2_SKY_RING)

      ! Free memory.
      deallocate(map)

    end function s2_sky_der_discrete_theta


    !--------------------------------------------------------------------------
    ! s2_sky_der_discrete_theta_fovop
    !
    !! Compute a sparse matrix representation of the discrete derivative
    !! operator on the sphere with respect to theta (using 
    !! convolutional or nearest interpolation).
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!
    !! Variables:
    !!   - sky: Sky to compute derivate of.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !!   - nearest: Logical to specify nearest interpolation rather than
    !!     convolutional.
    !!   - theta_fov: Size of field of view (fov), note that fov extends
    !!     from 0 to theta_fov/2 .
    !!   - nop: Number of non-zero enties in the derivative operator.
    !!   - op(0:nop-1,0:2): Discrete theta derivative operator (specifies 
    !!     indices and values of entries in operator matrix).
    !!   - nsphere: Number of pixels on the sphere within the fov.
    !!   - xmap(0:nsphere-1): Vector of pixels values on the sphere within 
    !!     the fov.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_der_discrete_theta_fovop(sky, &
         support_theta, kernel, param, inclusive, nearest, &
         theta_fov, nop, op, nsphere, xmap)

      use pix_tools, only: pix2ang_ring, ring_num, ring2z, ang2pix_ring

      type(s2_sky), intent(inout) :: sky 
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      logical, intent(in), optional :: nearest
      real(s2_dp), intent(in) :: theta_fov
      integer, intent(out) :: nop
      real(s2_dp), allocatable, intent(out) :: op(:,:)
      integer, intent(out) :: nsphere
      real(s2_sp), allocatable, intent(out) :: xmap(:)
      interface 
         function kernel(theta, param) result(val)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta
           real(s2_dp), intent(in), optional :: param(:)
           real(s2_dp) :: val
         end function kernel
      end interface

      integer :: ipix, ipix_adj, iring, iop, iw, iwa
      integer :: fail = 0
      real(s2_dp) :: theta, phi, theta_adj, dtheta
      real(s2_dp) :: icurr, wcurr
      real(s2_dp) :: z, zadj
      real(s2_dp), allocatable :: opx(:,:)
      integer :: nweights, nweights_adj
      integer, allocatable :: indices(:), indices_adj(:)
      real(s2_dp), allocatable :: weights(:), weights_adj(:)
      logical :: nearest_use
      real(s2_dp), parameter :: TOL = 1d-8

      if(present(nearest)) then
         nearest_use = nearest
      else
         nearest_use = .false.
      end if

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der_discrete_theta_fovop')
      end if 

      ! Check map computed.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_der_discrete_theta_fovop')
      end if

      ! Ensure map in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Count number of pixels within FOV.
      nsphere = 0
      do ipix = 0,sky%npix-1 
         call pix2ang_ring(sky%nside, ipix, theta, phi)
         if(theta > theta_fov/2.0) then
            nsphere = ipix
            exit            
         end if
      end do

      ! Get xmap vector.
      allocate(xmap(0:nsphere-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_theta_fovop')
      end if
      xmap(0:nsphere-1) = sky%map(0:nsphere-1)

      ! Allocate maximum required space for sparse representation of operator.
      allocate(opx(0:nsphere*nsphere-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_theta_fovop')
      end if  

      ! Compute sparse matrix representation of convolution operator.
      iop = 0
      do ipix = 0,nsphere-1

         ! Get theta and phi angles corresponding to pixel.
         call pix2ang_ring(sky%nside, ipix, theta, phi)
        
         if(theta <= theta_fov/2.0) then

            ! Compute ring number.
            z = cos(theta)
            iring = ring_num(sky%nside, z)

            ! If final ring set derivative to zero, i.e. all zero 
            ! row of matrix.
            ! (Note iring ranges from [1, 4*Nside-1] inclusive.)
            if (iring /= 4*sky%nside-1) then

               ! Compute theta for next ring.
               zadj = ring2z(sky%nside, iring+1)
               ! Eliminate numerical noise that could otherwise cause instability.
               if (zadj > 1.0) zadj = 1.0
               if (zadj < -1.0) zadj = -1.0
               theta_adj = acos(zadj)
               dtheta = theta_adj - theta
               if(abs(dtheta) < TOL) then
                  call s2_error(S2_ERROR_ARTH, 's2_sky_der_discrete_theta_fovop')
               end if

               if (nearest_use) then

                  ! Find nearest pixel to use in computing finite difference.
                  call ang2pix_ring(sky%nside, theta_adj, phi, ipix_adj)

                  if(theta_adj <= theta_fov/2.0) then
                     opx(iop, 0) = ipix
                     opx(iop, 1) = ipix_adj
                     opx(iop, 2) = 1d0 / dtheta
                     iop = iop + 1

                     opx(iop, 0) = ipix
                     opx(iop, 1) = ipix
                     opx(iop, 2) = -1d0 / dtheta
                     iop = iop + 1
                  end if
               else

                  ! Compute convolution weights.
                  call s2_sky_convpt_space_weights(nweights_adj, indices_adj, &
                       weights_adj, &
                       sky%nside, S2_SKY_RING, support_theta, kernel, &
                       theta_adj, phi, param, inclusive)
                  call s2_sky_convpt_space_weights(nweights, indices, weights, &
                       sky%nside, S2_SKY_RING, support_theta, kernel, &
                       theta, phi, param, inclusive)

                  ! Add entries due to adjacent point and given 
                  ! point if it coincides.
                  do iwa = 0,nweights_adj-1
                     if(indices_adj(iwa) < nsphere) then   
                        ! Ensure only include pixels within fov.             
                        icurr = indices_adj(iwa)
                        wcurr = weights_adj(iwa)
                        do iw = 0,nweights-1
                           if(indices(iw) == icurr) then
                              wcurr = wcurr - weights(iw)
                              indices(iw) = -1
                           end if
                        end do
                        opx(iop, 0) = ipix
                        opx(iop, 1) = icurr
                        opx(iop, 2) = wcurr / dtheta
                        iop = iop + 1
                     end if
                  end do

                  ! Add entries due to remaining given points.
                  do iw = 0,nweights-1
                     if(indices(iw) < nsphere) then   
                        if(indices(iw) /= -1) then
                           ! Ensure only include pixels within fov.             
                           opx(iop, 0) = ipix
                           opx(iop, 1) = indices(iw)
                           opx(iop, 2) = -weights(iw) / dtheta
                           iop = iop + 1
                        end if
                     end if
                  end do

                  ! Free temporary memory.
                  deallocate(indices_adj, weights_adj)
                  deallocate(indices, weights)

               end if
            end if

         end if

      end do
      nop = iop

      ! Copy required operator points.
      allocate(op(0:nop-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_theta_fovop')
      end if   
      op(0:nop-1, 0:2) = opx(0:nop-1, 0:2)

      ! Free memory.
      deallocate(opx)

    end subroutine s2_sky_der_discrete_theta_fovop


    !--------------------------------------------------------------------------
    ! s2_sky_der_discrete_grad
    !
    !! Compute discrete gradient on sphere with respect to theta (using 
    !! convolutional or nearest interpolation for dtheta).
    !!
    !! Notes:
    !!   - Memory allocated herein for der sphere; must be freed by calling
    !!     routine.
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!
    !! Variables:
    !!   - sky: Sky to compute derivate of.
    !!   - divbysin: Logical specifying whether to divide by sin(theta) 
    !!     ( if true, return 1/sin(theta) dT/dphi; else return dT/dphi).
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !!   - nearest: Logical to specify nearest interpolation rather than
    !!     convolutional.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_sky_der_discrete_grad(sky, divbysin, support_theta, kernel, &
         param, inclusive, nearest) result (der)

      type(s2_sky), intent(inout) :: sky
      logical, intent(in) :: divbysin
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      logical, intent(in), optional :: nearest
      type(s2_sky) :: der
      interface 
         function kernel(theta, param) result(val)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta
           real(s2_dp), intent(in), optional :: param(:)
           real(s2_dp) :: val
         end function kernel
      end interface

      real(s2_sp), allocatable :: map(:)
      integer :: fail = 0
      integer :: ipix
      type(s2_sky) :: dphi
      type(s2_sky) :: dtheta

      ! Check object initialised.
      if(.not. sky%init) then
         call s2_error(S2_ERROR_NOT_INIT, 's2_sky_der_discrete_grad')
      end if

      ! Check map computed.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_der_discrete_grad')
      end if

      ! Ensure in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Compute partials.
      dphi = s2_sky_der_discrete_phi(sky, divbysin)
      dtheta = s2_sky_der_discrete_theta(sky, support_theta, kernel, param, &
           inclusive, nearest)

      ! Allocate space for temporary map.
      allocate(map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_der_discrete_grad')
      end if

      ! Compute gradient.
      map(0:sky%npix-1) = &
           sqrt(dphi%map(0:sky%npix-1)**2 + dtheta%map(0:sky%npix-1)**2)

      ! Initialise output sky with new map.
      der = s2_sky_init_map(map, sky%nside, S2_SKY_RING)

      ! Free memory.
      deallocate(map)
      call s2_sky_free(dphi)
      call s2_sky_free(dtheta)

    end function s2_sky_der_discrete_grad


    !--------------------------------------------------------------------------
    ! s2_sky_conv
    !
    !! Convolve a sky with a beam.
    !!
    !! Notes:
    !!   - Sky alms must already be defined since if not lmax and mmax may not
    !!     be defined.
    !!
    !! Variables:
    !!   - sky: Sky to convolve with beam.  Convolved sky on ouput.
    !!   - beam: Beam to convolve with sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_conv(sky, beam)

      type(s2_sky), intent(inout) :: sky
      type(s2_pl), intent(in) :: beam

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_conv')
      end if 

      ! Check alms computed.
      ! Must already be computed since lmax and mmax may otherwise be unknown.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_conv')
      end if

      ! Convolve sky with beam (overwrite alms).
      call s2_pl_conv(beam, sky%alm)

      ! If map oringally defined then must clear and recompute.
      if(sky%map_status) then
         sky%map_status = .false.
         sky%map = 0.0e0
         deallocate(sky%map)
         call s2_sky_compute_map(sky)
      end if
         
    end subroutine s2_sky_conv


    !--------------------------------------------------------------------------
    ! s2_sky_conv_space
    !
    !! Compute the convolution in real space of a sky with a kernel.
    !!
    !! Notes:
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!
    !! Variables:
    !!   - sky: Sky to compute convolution.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_conv_space(sky, support_theta, kernel, param, inclusive)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      type(s2_sky), intent(inout) :: sky
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      interface 
        function kernel(theta, param) result(val)
          use s2_types_mod
          real(s2_dp), intent(in) :: theta
          real(s2_dp), intent(in), optional :: param(:)
          real(s2_dp) :: val
        end function kernel
      end interface

      real(s2_sp), allocatable :: map(:)
      integer :: ipix, fail=0
      real(s2_dp) :: theta0
      real(s2_dp) :: phi0

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_conv_space')
      end if 

      ! Check map computed.
      ! Must already be computed since lmax and mmax may otherwise be unknown.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_conv_space')
      end if

      ! Allocate space for temporary map.
      allocate(map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_conv_space')
      end if

      ! Perform convolution.
      do ipix = 0,sky%npix-1

         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta0, phi0)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta0, phi0)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_conv_space')
         end if

         map(ipix) = s2_sky_convpt_space(sky, support_theta, kernel, &
              theta0, phi0, param, inclusive)

      end do

      ! Copy temporary map.
      sky%map(0:sky%npix-1) = map(0:sky%npix-1)

      ! Free temporary memory.
      deallocate(map)
    
    end subroutine s2_sky_conv_space
       

    !--------------------------------------------------------------------------
    ! s2_sky_convpt_space
    !
    !! Compute the convolution in real space of a sky with a kernel at
    !! a given position.
    !!
    !! Notes:
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!
    !! Variables:
    !!   - sky: Sky to compute convolution.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - theta0: Theta position of point to compute convolution at.
    !!   - phi0: Phi position of point to compute convolution at.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_convpt_space(sky, support_theta, kernel, &
         theta0, phi0, param, inclusive) result(val)

      use pix_tools, only: query_disc, pix2ang_ring, pix2ang_nest

      type(s2_sky), intent(in) :: sky
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in) :: theta0
      real(s2_dp), intent(in) :: phi0
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      real(s2_dp) :: val
      interface 
        function kernel(theta, param) result(val)
          use s2_types_mod
          real(s2_dp), intent(in) :: theta
          real(s2_dp), intent(in), optional :: param(:)
          real(s2_dp) :: val
        end function kernel
      end interface

      type(s2_vect) :: vec0, vec
      real(s2_sp) :: x0_sp(1:3)
      real(s2_dp) :: x0_dp(1:3)
      real(s2_dp) :: theta, phi, dot, ang, normalisation
      integer :: nest, ndisc, ipix, idisc, fail=0
      integer, allocatable :: disc_ipix(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_convpt_space')
      end if 

      ! Check map computed.
      ! Must already be computed since lmax and mmax may otherwise be unknown.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_convpt_space')
      end if

      ! Compute cartesial coordinates of (theta,phi).
      vec0 = s2_vect_init(1.0, real(theta0,s2_sp), real(phi0,s2_sp))
      call s2_vect_convert(vec0, S2_VECT_TYPE_CART)
      x0_sp = s2_vect_get_x(vec0)
      x0_dp(1:3) = x0_sp(1:3)

      ! Query disc.
      if(sky%pix_scheme == S2_SKY_RING) then
         nest = 0
      else if(sky%pix_scheme == S2_SKY_NEST) then
         nest = 1
      else
         call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_convpt_space')
      end if
      allocate(disc_ipix(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_convpt_space')
      end if      
      call query_disc(sky%nside, x0_dp, support_theta/2.0, disc_ipix, ndisc, nest, inclusive)

      ! Compute convolution.
      val = 0.0
      normalisation = 0.0
      do idisc = 0, ndisc-1

         ! Get index of pixel in disk.
         ipix = disc_ipix(idisc)

         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta, phi)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_convpt_space')
         end if

         ! Compute angle between this point and x0.
         vec = s2_vect_init(1.0, real(theta,s2_sp), real(phi,s2_sp))
         dot = s2_vect_dot(vec0, vec)         
         if (dot > 1d0) dot = 1d0  ! Remove numerical noise that could cause acos to fail.
         ang = acos(dot)
         call s2_vect_free(vec)

         ! Compute contribution to convolution.
         val = val + kernel(ang, param) * sky%map(ipix)

         normalisation = normalisation + kernel(ang, param)
         
      end do

      ! Normalise.
      val = val / normalisation

      ! Free memory.
      deallocate(disc_ipix)
      call s2_vect_free(vec0)

    end function s2_sky_convpt_space


    !--------------------------------------------------------------------------
    ! s2_sky_conv_space_fovop
    !
    !! Compute a sparse matrix representation of the spatial
    !! convolution operator for a polar cap within a specified
    !! field-of-view.
    !!
    !! Notes:
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!   - Memory allocated herein; must be freed by calling routine.
    !!
    !! Variables:
    !!   - sky: Sky to extract polar cap.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - theta_fov: Size of field of view (fov), note that fov extends
    !!     from 0 to theta_fov/2 .
    !!   - nop: Number of non-zero enties in the convolution operator.
    !!   - op(0:nop-1,0:2): Convolution operator (specifies indices and values 
    !!     of entries in operator matrix).
    !!   - nsphere: Number of pixels on the sphere within the fov.
    !!   - xmap(0:nsphere-1): Vector of pixels values on the sphere within 
    !!     the fov.
    !!   - kernel: Convolution kernel.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_conv_space_fovop(sky, support_theta, theta_fov, nop, op, &
         nsphere, xmap, kernel, param, inclusive)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      type(s2_sky), intent(inout) :: sky  ! Inout since may need to change pix_scheme.
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in) :: theta_fov
      integer, intent(out) :: nop
      real(s2_dp), allocatable, intent(out) :: op(:,:)
      integer, intent(out) :: nsphere
      real(s2_sp), allocatable, intent(out) :: xmap(:)
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      interface 
        function kernel(theta, param) result(val)
          use s2_types_mod
          real(s2_dp), intent(in) :: theta
          real(s2_dp), intent(in), optional :: param(:)
          real(s2_dp) :: val
        end function kernel
      end interface

      integer :: ipix, iop
      integer :: fail = 0
      real(s2_dp) :: theta, phi
      real(s2_dp), allocatable :: opx(:,:)
      integer :: nweights, iweight
      integer, allocatable :: indices(:)
      real(s2_dp), allocatable :: weights(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_conv_space_fovop')
      end if 

      ! Check map computed.
      ! Must already be computed since lmax and mmax may otherwise be unknown.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_conv_space_fovop')
      end if

      ! Ensure map in RING pixelisation scheme.
      call s2_sky_map_convert(sky, S2_SKY_RING)

      ! Count number of pixels within FOV.
      nsphere = 0
      do ipix = 0,sky%npix-1 
         call pix2ang_ring(sky%nside, ipix, theta, phi)
         if(theta > theta_fov/2.0) then
            nsphere = ipix
            exit            
         end if
      end do

      ! Get xmap vector.
      allocate(xmap(0:nsphere-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_conv_space_fovop')
      end if
      xmap(0:nsphere-1) = sky%map(0:nsphere-1)

      ! Allocate maximum required space for sparse representation of operator.
      allocate(opx(0:nsphere*nsphere-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_conv_space_fovop')
      end if  

      ! Compute sparse matrix representation of convolution operator.
      iop = 0
      do ipix = 0,nsphere-1

         ! Get theta and phi angles corresponding to pixel.
         call pix2ang_ring(sky%nside, ipix, theta, phi)
        
         if(theta <= theta_fov/2.0) then

            call s2_sky_convpt_space_weights(nweights, indices, weights, &
                 sky%nside, S2_SKY_RING, support_theta, kernel, &
                 theta, phi, param)

            do iweight = 0,nweights-1
               if(indices(iweight) < nsphere) then   
                  ! Ensure only include pixels within fov.             
                  opx(iop, 0) = ipix
                  opx(iop, 1) = indices(iweight)
                  opx(iop, 2) = weights(iweight)
                  iop = iop + 1
               end if
            end do

            deallocate(indices, weights)

         end if

      end do
      nop = iop

      ! Copy required operator points.
      allocate(op(0:nop-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_conv_space_fovop')
      end if   
      op(0:nop-1, 0:2) = opx(0:nop-1, 0:2)

      ! Free memory.
      deallocate(opx)

    end subroutine s2_sky_conv_space_fovop


    !--------------------------------------------------------------------------
    ! s2_sky_convpt_space_weights
    !
    !! Compute the weights and indices to perform convolution in real
    !! space at a given position.  i.e. compute one row of the matrix
    !! convolution operator.
    !!
    !! Notes:
    !!   - Support of the kernel is the full support, i.e. kernel has support
    !!     with theta in [-support_theta/2, support_theta/2].
    !!   - Memory allocated herein; must be freed by calling routine.
    !!
    !! Variables:
    !!   - nweights: Number of weights 
    !!   - indices(0:nweights-1): Indices of points on sphere required in 
    !!     computation of convolution.
    !!   - weights(0:nweights-1): Weight for each pixel.
    !!   - nside: Resolution of Healpix sphere considered in computing pixel 
    !!     positions.
    !!   - pix_schmeme: Pixelisation scheme of pixel indices.
    !!   - support_theta: Full support of the convolution kernel (see note).
    !!   - kernel: Convolution kernel.
    !!   - theta0: Theta position of point to compute convolution at.
    !!   - phi0: Phi position of point to compute convolution at.
    !!   - param: Optional parameters to pass to the kernel function.
    !!   - inclusive: If set to 1, all the pixels overlapping (even
    !!     partially) with the disc are listed, otherwise only those
    !!     whose center lies within the disc are listed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_convpt_space_weights(nweights, indices, weights, &
         nside, pix_scheme, support_theta, kernel, &
         theta0, phi0, param, inclusive)

      use pix_tools, only: query_disc, pix2ang_ring, pix2ang_nest, nside2npix

      integer, intent(out) :: nweights
      integer, intent(out), allocatable :: indices(:)
      real(s2_dp), intent(out), allocatable :: weights(:)
      integer, intent(in) :: nside
      integer, intent(in) :: pix_scheme
      real(s2_dp), intent(in) :: support_theta
      real(s2_dp), intent(in) :: theta0
      real(s2_dp), intent(in) :: phi0
      real(s2_dp), intent(in), optional :: param(:)
      integer, intent(in), optional :: inclusive
      real(s2_dp) :: val
      interface 
        function kernel(theta, param) result(val)
          use s2_types_mod
          real(s2_dp), intent(in) :: theta
          real(s2_dp), intent(in), optional :: param(:)
          real(s2_dp) :: val
        end function kernel
      end interface

      type(s2_vect) :: vec0, vec
      real(s2_sp) :: x0_sp(1:3)
      real(s2_dp) :: x0_dp(1:3)
      real(s2_dp) :: theta, phi, dot, ang, normalisation
      integer :: nest, ndisc, ipix, idisc, fail=0, npix
      integer, allocatable :: disc_ipix(:)

      ! Compute cartesial coordinates of (theta,phi).
      vec0 = s2_vect_init(1.0, real(theta0,s2_sp), real(phi0,s2_sp))
      call s2_vect_convert(vec0, S2_VECT_TYPE_CART)
      x0_sp = s2_vect_get_x(vec0)
      x0_dp(1:3) = x0_sp(1:3)

      ! Query disc.
      if(pix_scheme == S2_SKY_RING) then
         nest = 0
      else if(pix_scheme == S2_SKY_NEST) then
         nest = 1
      else
         call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_convpt_space_weights')
      end if
      npix = nside2npix(nside) 
      allocate(disc_ipix(0:npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_convpt_space_weights')
      end if      
      call query_disc(nside, x0_dp, support_theta/2.0, disc_ipix, ndisc, nest, inclusive)

      ! Allocate space.
      nweights = ndisc
      allocate(indices(0:nweights-1), stat=fail)
      allocate(weights(0:nweights-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_convpt_space_weights')
      end if
      
      ! Compute convolution weights.
      indices(0:nweights-1) = disc_ipix(0:nweights-1)
      normalisation = 0.0
      do idisc = 0, ndisc-1

         ! Get index of pixel in disk.
         ipix = disc_ipix(idisc)

         ! Get theta and phi angles corresponding to pixel.
         if(pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(nside, ipix, theta, phi)
         else if(pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_convpt_space_weights')
         end if

         ! Compute angle between this point and x0.
         vec = s2_vect_init(1.0, real(theta,s2_sp), real(phi,s2_sp))
         dot = s2_vect_dot(vec0, vec)         
         if (dot > 1d0) dot = 1d0  ! Remove numerical noise that could cause acos to fail.
         ang = acos(dot)
         call s2_vect_free(vec)

         ! Compute contribution to convolution.
         weights(idisc) = kernel(ang, param)
         normalisation = normalisation + kernel(ang, param)
         
      end do

      ! Normalise.
      weights(0:nweights-1) = weights(0:nweights-1) / normalisation

      ! Free memory.
      deallocate(disc_ipix)
      call s2_vect_free(vec0)

    end subroutine s2_sky_convpt_space_weights


    !--------------------------------------------------------------------------
    ! s2_sky_offset
    !
    !! Add an offset to the sky.
    !!
    !! Variables:
    !!   - sky: Sky to add offset to.
    !!   - offset: Offset to add.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2006
    !
    ! Revisions:
    !   June 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_offset(sky, offset)

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: offset

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_offset')
      end if 

      ! Check map present.
      if(.not. sky%map_status) call s2_sky_compute_map(sky)

      ! Add offset.
      sky%map = sky%map + offset
      
      ! Recompute alms if were present previously.
      if(sky%alm_status) then
         sky%alm_status = .false.
         deallocate(sky%alm)
         call s2_sky_compute_alm(sky)
      end if
      
    end subroutine s2_sky_offset


    !--------------------------------------------------------------------------
    ! s2_sky_scale
    !
    !! Scale both the sky map and alm.
    !!
    !! Variables:
    !!   - sky: Sky to scale.
    !!   - scale: Scale to apply
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_scale(sky, scale)

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: scale

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_scale')
      end if 

      if(sky%map_status) sky%map = sky%map * scale
      if(sky%alm_status) sky%alm = sky%alm * scale

    end subroutine s2_sky_scale


    !--------------------------------------------------------------------------
    ! s2_sky_add
    !
    !! Add the maps of two skys.  
    !!
    !! Variables:
    !!  - a: First sky to be added.
    !!  - b: Second sky to be added
    !!  - [subtract]: Logcal specifying whether to subtract, rather than 
    !!    add maps.
    !!  - c: The output *initialised* sky whos map is the sum of the maps of 
    !!    the other two passed skies.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_add(a, b, subtract) result(c)

      type(s2_sky), intent(inout) :: a, b  
      logical, intent(in), optional :: subtract
        ! Need to be inout incase have to compute map.
      type(s2_sky) :: c

      logical :: subtract_use
      integer :: nside, npix, pix_scheme, lmax, mmax, fail
      real(s2_sp), allocatable :: map_temp(:)

      subtract_use = .false.
      if(present(subtract)) subtract_use = subtract

      ! Check objects initialised.
      if((.not. a%init) .or.(.not. b%init)) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_add')
      end if 

      ! Check skies to be added have same nside.
      if(a%nside /= b%nside) then
         call s2_error(S2_ERROR_SKY_NON_CONFORM, 's2_sky_add', &
           comment_add='Skies to be added have different nside')
      end if
      
      ! Check skies to be added are in the same pixelisation scheme.
      if(a%pix_scheme /= b%pix_scheme) then
         call s2_error(S2_ERROR_SKY_PIX_DIFF, 's2_sky_add', &
           comment_add='Converting sky b to pixelisation scheme of sky a')
         ! Convert pix scheme of b to a.
         call s2_sky_map_convert(b, a%pix_scheme)
      end if

      ! Compute maps if not already computed.
      if(.not. a%map_status) call s2_sky_compute_map(a)
      if(.not. b%map_status) call s2_sky_compute_map(b)

      ! Define parameters for resultant sky.
      nside = a%nside
      npix = a%npix
      pix_scheme = a%pix_scheme
      lmax = min(a%lmax, b%lmax)  ! If alm not defined then will be zero.
      mmax = min(a%mmax, b%mmax)  ! Ditto.

      ! Generate added map
      allocate(map_temp(0:npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_add')
      end if

      if(subtract_use) then
         map_temp = a%map - b%map
      else
         map_temp = a%map + b%map
      end if

      ! Initialise output sky with new map.
      c = s2_sky_init_map(map_temp, nside, pix_scheme, lmax, mmax)

      ! Free memory used for temporary storage.
      deallocate(map_temp)

    end function s2_sky_add


    !--------------------------------------------------------------------------
    ! s2_sky_add_alm
    !
    !! Add the alms of two skys.  
    !!
    !! Variables:
    !!  - a: First sky to be added.
    !!  - b: Second sky to be added
    !!  - [subtract]: Logcal specifying whether to subtract, rather than 
    !!    add alms.
    !!  - c: The output *initialised* sky whos alms are the sum of the alms of 
    !!    the other two passed skies.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_add_alm(a, b, subtract) result(c)

      type(s2_sky), intent(inout) :: a, b  
      logical, intent(in), optional :: subtract
        ! Need to be inout incase have to compute alms.
      type(s2_sky) :: c

      logical :: subtract_use
      integer :: nside, npix, pix_scheme, lmax, mmax, fail
      complex(s2_spc), allocatable :: alm_temp(:,:)

      subtract_use = .false.
      if(present(subtract)) subtract_use = subtract

      ! Check objects initialised.
      if((.not. a%init) .or.(.not. b%init)) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_add_alm')
      end if 

      ! Check skies to be added have same lmax and mmax.
      if(a%lmax /= b%lmax .or. a%mmax /= b%mmax) then
         call s2_error(S2_ERROR_SKY_NON_CONFORM, 's2_sky_add_alm', &
           comment_add='Skies to be added have different lmax or mmax')
      end if
     
      ! Compute alms if not already computed.
      if(.not. a%alm_status) call s2_sky_compute_alm(a)
      if(.not. b%alm_status) call s2_sky_compute_alm(b)

      ! Define parameters for resultant sky.
      nside = a%nside
      npix = a%npix
      pix_scheme = a%pix_scheme
      lmax = a%lmax
      mmax = a%mmax

      ! Generate added map
      allocate(alm_temp(0:lmax,0:mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_add_alm')
      end if

      if(subtract_use) then
         alm_temp = a%alm(0:lmax,0:mmax) - b%alm(0:lmax,0:mmax)
      else
         alm_temp = a%alm(0:lmax,0:mmax) + b%alm(0:lmax,0:mmax)
      end if

      ! Initialise output sky with new alms.
      c = s2_sky_init_alm(alm_temp, lmax, mmax)

      ! Free memory used for temporary storage.
      deallocate(alm_temp)

    end function s2_sky_add_alm


    !--------------------------------------------------------------------------
    ! s2_sky_product
    !
    !! Multiply two sky maps together.
    !!
    !! Variables:
    !!  - a: First sky.
    !!  - b: Second sky.
    !!  - [divide]: Logcal specifying whether to dive, rather than 
    !!    multiply maps (if so result is a/b).
    !!  - c: The output *initialised* sky whos map is the product of the maps
    !!    of the other two passed skies.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2005
    !
    ! Revisions:
    !   August 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_product(a, b, divide) result(c)

      type(s2_sky), intent(inout) :: a, b  
      logical, intent(in), optional :: divide
      type(s2_sky) :: c

      logical :: divide_use
        ! Don't set divide_use value here, if do it becomes a saved variable. 
      integer :: nside, npix, pix_scheme, lmax, mmax, fail
      real(s2_sp), allocatable :: map_temp(:)

      divide_use = .false.
      if(present(divide)) divide_use = divide

      ! Check objects initialised.
      if((.not. a%init) .or.(.not. b%init)) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_product')
      end if 

      ! Check skies to be added have same nside.
      if(a%nside /= b%nside) then
         call s2_error(S2_ERROR_SKY_NON_CONFORM, 's2_sky_product', &
           comment_add='Skies to be multiplied have different nside')
      end if
      
      ! Check skies to be added are in the same pixelisation scheme.
      if(a%pix_scheme /= b%pix_scheme) then
!**         call s2_error(S2_ERROR_SKY_PIX_DIFF, 's2_sky_product', &
!           comment_add='Converting sky b to pixelisation scheme of sky a')
         ! Convert pix scheme of b to a.
         call s2_sky_map_convert(b, a%pix_scheme)
      end if

      ! Compute maps if not already computed.
      if(.not. a%map_status) call s2_sky_compute_map(a)
      if(.not. b%map_status) call s2_sky_compute_map(b)

      ! Define parameters for resultant sky.
      nside = a%nside
      npix = a%npix
      pix_scheme = a%pix_scheme
      lmax = min(a%lmax, b%lmax)  ! If alm not defined then will be zero.
      mmax = min(a%mmax, b%mmax)  ! Ditto.

      ! Allocate space for product map.
      allocate(map_temp(0:npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_product')
      end if

      ! Compute product map.
      if(divide_use) then
         map_temp(0:npix-1) = a%map(0:npix-1) / b%map(0:npix-1)
      else
         map_temp(0:npix-1) = a%map(0:npix-1) * b%map(0:npix-1)
      end if

      ! Initialise output sky with new map.
      c = s2_sky_init_map(map_temp, nside, pix_scheme, lmax, mmax)

      ! Free memory used for temporary storage.
      deallocate(map_temp)

    end function s2_sky_product
    

    !--------------------------------------------------------------------------
    ! s2_sky_thres
    !
    !! Threshold all values outside range to specified limits.
    !!
    !! Variables:
    !!  - sky: Sky to threshold.
    !!  - thres_lower: Lower threshold value.
    !!  - thres_upper: Upper threshold value.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2007
    !
    ! Revisions:
    !   April 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_thres(sky, thres_lower, thres_upper)

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: thres_lower, thres_upper

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_thres')
      end if 

      ! Compute map if not already computed.
      if(.not. sky%map_status) call s2_sky_compute_map(sky)

      ! Perform thresholding.
      where(sky%map(0:sky%npix-1) <= thres_lower)
         sky%map(0:sky%npix-1) = thres_lower
      end where

      where(sky%map(0:sky%npix-1) >= thres_upper)
         sky%map(0:sky%npix-1) = thres_upper
      end where
      
      ! Recompute alms if were present previously.
      if(sky%alm_status) then
         sky%alm_status = .false.
         deallocate(sky%alm)
         call s2_sky_compute_alm(sky)
      end if
      
    end subroutine s2_sky_thres


    !--------------------------------------------------------------------------
    ! s2_sky_thres_abs
    !
    !! Threshold all values below (in absolute value) specified threshold.
    !!
    !! Variables:
    !!  - sky: Sky to threshold.
    !!  - thres: Threshold value.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2007
    !
    ! Revisions:
    !   April 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_thres_abs(sky, thres)

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: thres

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_thres')
      end if 

      ! Compute map if not already computed.
      if(.not. sky%map_status) call s2_sky_compute_map(sky)

      ! Perform thresholding.
      where(abs(sky%map(0:sky%npix-1)) <= thres)
         sky%map(0:sky%npix-1) = 0e0
      end where

      ! Recompute alms if were present previously.
      if(sky%alm_status) then
         sky%alm_status = .false.
         deallocate(sky%alm)
         call s2_sky_compute_alm(sky)
      end if

    end subroutine s2_sky_thres_abs


    !--------------------------------------------------------------------------
    ! s2_sky_error_twonorm
    !
    !! Compute two-norm error (i.e. mean square error) or two sky maps.
    !!
    !! Variables:
    !!  - a: First sky.
    !!  - b: Second sky.
    !!  - error: Two-norm error computed by 
    !!    er = sqrt( int( (a(w)-b(w))^2 )dw ).
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2007
    !
    ! Revisions:
    !   April 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_error_twonorm(a, b) result(er)

      type(s2_sky), intent(inout) :: a, b  
      real(s2_sp) :: er

      er = 0e0

      ! Check objects initialised.
      if((.not. a%init) .or.(.not. b%init)) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_error_twonorm')
      end if 

      ! Check skies to be added have same nside.
      if(a%nside /= b%nside) then
         call s2_error(S2_ERROR_SKY_NON_CONFORM, 's2_sky_error_twonorm', &
           comment_add='Skies to be multiplied have different nside')
      end if
      
      ! Check skies to be added are in the same pixelisation scheme.
      if(a%pix_scheme /= b%pix_scheme) then
!**         call s2_error(S2_ERROR_SKY_PIX_DIFF, 's2_sky_error_twonorm', &
!           comment_add='Converting sky b to pixelisation scheme of sky a')
         ! Convert pix scheme of b to a.
         call s2_sky_map_convert(b, a%pix_scheme)
      end if

      ! Compute maps if not already computed.
      if(.not. a%map_status) call s2_sky_compute_map(a)
      if(.not. b%map_status) call s2_sky_compute_map(b)

      ! Compute two norm error.
      er = sum((a%map(0:a%npix-1) - b%map(0:b%npix-1))**2)
      er = sqrt(er * 4e0 * pi / real(a%npix,s2_sp))

    end function s2_sky_error_twonorm
   

    !--------------------------------------------------------------------------
    ! s2_sky_rms
    !
    !! Compute root-mean squared value of sky (same as error_twonorm with 
    !! second sky set to zero).
    !!
    !! Variables:
    !!  - a: Sky to compute RMS value of.
    !!  - rms: RMS value of sky computed by 
    !!    er = sqrt( int( (a(w))^2 )dw ).
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2007
    !
    ! Revisions:
    !   April 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_rms(a) result(rms)

      type(s2_sky), intent(inout) :: a
      real(s2_sp) :: rms

      rms = 0e0

      ! Check object initialised.
      if(.not. a%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_rms')
      end if 

      ! Compute map if not already computed.
      if(.not. a%map_status) call s2_sky_compute_map(a)

      ! Compute two norm error.
      rms = sum(a%map(0:a%npix-1)**2)
      rms = sqrt(rms * 4e0 * pi / real(a%npix,s2_sp))

    end function s2_sky_rms
   

    !--------------------------------------------------------------------------
    ! s2_sky_dilate
    !
    !! Dilate the map of a sky.  The dilation may be performed in two 
    !! orthogonal directions.  To perform the usual isotropic dilation set
    !! a=b.
    !!
    !! Variables:
    !!   - sky: Sky to dilate. 
    !!   - a: First dilation scale.
    !!   - b: Second dilation scale.
    !!   - norm_preserve: Logical specifying whether a norm preserving 
    !!     `cocycle' factor is to be applied.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_dilate(sky, a, b, norm_preserve)

      use pix_tools, only: pix2ang_ring, pix2ang_nest, &
        ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: a, b
      logical, intent(in) :: norm_preserve

      integer :: ipix, ipix_dil, fail
      real(s2_dp) :: theta, phi, theta_dil, phi_dil
      real(s2_sp) :: D, R, A2, cocycle
      real(s2_sp), allocatable :: map_dil(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_dilate')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_dilate')
         return
      end if

      ! Allocate temporary storage space for dilated map.
      allocate(map_dil(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_dilate')
      end if

      do ipix = 0,sky%npix-1

         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta, phi)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_dilate')
         end if

         ! Compute dilated theta and phi values.
         D = 1/a**2 * (cos(phi))**2 + 1/b**2 * (sin(phi))**2
         theta_dil = 2 * atan( tan(theta / 2.0e0) * sqrt(D))

         ! Must account for correct quadrant that phi lies in.
         ! Otherwise inverse will always give IV and I quadrants.
         phi = mod(phi, 2.0d0*pi)            ! Restrict phi to range [0,2*pi).
         if(phi >= 0.0e0 .and. phi < pi/2.0e0) then
            phi_dil = atan( a/b * tan(phi) )
         elseif(phi >= pi/2.0e0 .and. phi < pi) then
            phi_dil = atan( a/b * tan(phi) ) + pi
         elseif(phi >= pi .and. phi < 3.0e0*pi/2.0e0) then
            phi_dil = atan( a/b * tan(phi) ) + pi
         elseif(phi >= 3.0e0*pi/2.0e0) then
            phi_dil = atan( a/b * tan(phi) ) + 2.0e0*pi
         end if

         ! Compute index corresponding to dilated angles.
         if(sky%pix_scheme == S2_SKY_RING) then
            call ang2pix_ring(sky%nside, theta_dil, phi_dil, ipix_dil)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call ang2pix_nest(sky%nside, theta_dil, phi_dil, ipix_dil)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_dilate')
         end if
         
         ! Compute norm preserving cocycle term.
         if(norm_preserve) then
            R = ( b**2 - (a**2)*(tan(phi)**2) ) &
              / ( b**2 + (a**2)*(tan(phi)**2) )
            A2 = a**2 / 2.0e0 * (R + 1) &
              + b**2 / 2.0e0 * (1 - R)
            cocycle = 2.0e0 * a / (b * cos(phi)**2) * (R + 1) &
              * A2 / (A2+1 + (A2-1)*cos(theta))**2
            ! Note A2 in code here = A^2 from theory.
         else
            cocycle = 1.0e0
         end if

         ! Compute dilated map.
         map_dil(ipix) = sqrt(cocycle) * sky%map(ipix_dil)

      end do

      ! Overwrite old map with dilated version.
      sky%map = map_dil

      ! Free temporary storage space used.
      deallocate(map_dil)

      ! Ensure alm not present.
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.
      
    end subroutine s2_sky_dilate


    !--------------------------------------------------------------------------
    ! s2_sky_rotate
    !
    !! Rotate the map of a sky.  The rotation is perform in real space.
    !! The rotation is defined by the rotation routine in the s2_vect_mod
    !! class.
    !!
    !! Variables
    !!  - sky: The sky rotated.
    !!  - alpha: Alpha Euler angle of the rotation.
    !!  - beta: Beta Euler angle of the rotation.
    !!  - gamma: Gamma Euler angle of the rotation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_rotate(sky, alpha, beta, gamma)

      use pix_tools, only: pix2ang_ring, pix2ang_nest, &
        ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(inout) :: sky
      real(s2_sp) :: alpha, beta, gamma
      integer :: ipix, ipix_rot, fail
      real(s2_dp) :: theta, phi, theta_rot, phi_rot
      type(s2_vect) :: vect
      real(s2_sp), allocatable :: map_rot(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_rotate')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_rotate')
         return
      end if

     ! Allocate temporary storage space for rotated map.
      allocate(map_rot(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_rotate')
      end if

      do ipix = 0,sky%npix-1
         
         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta, phi)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_rotate')
         end if

         ! Rotate coordinates.
         ! (Very inefficient since compute rotation matric every time).
         vect = s2_vect_init(1.0e0, real(theta,s2_sp), real(phi,s2_sp))
!         call s2_vect_rotate(vect, alpha, beta, gamma)
         call s2_vect_rotate(vect, -gamma, -beta, -alpha)
         theta_rot = s2_vect_get_theta(vect)
         phi_rot = s2_vect_get_phi(vect)
         call s2_vect_free(vect)

         ! Compute index corresponding to rotated coordinates.
         if(sky%pix_scheme == S2_SKY_RING) then
            call ang2pix_ring(sky%nside, theta_rot, phi_rot, ipix_rot)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call ang2pix_nest(sky%nside, theta_rot, phi_rot, ipix_rot)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_rotate')
         end if

         ! Compute rotated map.
         map_rot(ipix) = sky%map(ipix_rot)

      end do

      ! Overwrite old map with rotated version.
      sky%map(0:sky%npix-1) = map_rot(0:sky%npix-1)

      ! Free temporary storage space used.
      deallocate(map_rot)

      ! Ensure alm not present.
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.

    end subroutine s2_sky_rotate


    !--------------------------------------------------------------------------
    ! s2_sky_rotate_alm
    !
    !! Rotate the alms of a sky.  The rotation is perform in harmomnic space.
    !! The rotation is defined by the rotation routine in the s2_vect_mod
    !! class but is performing using dlmn Wigner coefficients.
    !!
    !! Variables
    !!  - sky: The sky rotated.
    !!  - alpha: Alpha Euler angle of the rotation.
    !!  - beta: Beta Euler angle of the rotation.
    !!  - gamma: Gamma Euler angle of the rotation.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_rotate_alm(sky, alpha, beta, gamma, azisym_in)

      use s2_dl_mod

      type(s2_sky), intent(inout) :: sky
      real(s2_dp) :: alpha, beta, gamma
      logical, optional :: azisym_in

      complex(s2_spc), allocatable :: alm_rot(:,:)
      integer :: fail, el, m, n
      real(s2_dp) :: nsign
      real(s2_dp), pointer :: dl(:,:)
      complex(s2_dpc) :: I
      logical :: azisym = .false.

      I = cmplx(0d0, 1d0)
      
      ! Parse inputs.
      if(present(azisym_in)) then
         azisym = azisym_in
      end if

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_rotate_alm')
      end if

      ! Check alm defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_rotate_alm')
         return
      end if

      ! Allocate temporary storage space for rotated alms.
      ! (Note that the mmax is lmax for the rotated function.)
      allocate(alm_rot(0:sky%lmax, 0:sky%lmax), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_rotate_alm')
      end if

      ! Rotate harmonic coefficients.
      do el = 0,sky%lmax

         ! Compute dlms.
         allocate(dl(-el:el,-el:el), stat=fail)
         if(fail /= 0) then
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_rotate_alm')
         end if
         call s2_dl_beta_operator(dl, beta, el)

         ! Perform rotation.
         do m = 0,el
            alm_rot(el,m) = dl(m,0) * exp(-I*m*alpha) * sky%alm(el,0)
            if (.not. azisym) then
               nsign = 1d0
               do n = 1,min(el,sky%mmax)
                  nsign = -nsign
                  !nsign = (-1)**(n)
                  alm_rot(el, m) = alm_rot(el,m) &
                       + dl(m,n) * exp(-I*m*alpha) * exp(-I*n*gamma) * sky%alm(el,n) &
                       + dl(m,-n) * exp(-I*m*alpha) * exp(+I*n*gamma) * nsign *conjg(sky%alm(el,n))
               end do
            end if
         end do

         ! Clear dlmns.
         deallocate(dl)

      end do

      ! Overwrite old map with rotated version.
      sky%alm(0:sky%lmax,0:sky%lmax) = alm_rot(0:sky%lmax,0:sky%lmax)
      sky%mmax = sky%lmax

      ! Free temporary storage space used.
      deallocate(alm_rot)

      ! Ensure map not present.
      if(allocated(sky%map)) deallocate(sky%map)
      sky%map_status = .false.

    end subroutine s2_sky_rotate_alm


    !--------------------------------------------------------------------------
    ! s2_sky_power_map
    !
    !! Compute power of function on sphere in map space.
    !!
    !! Variables:
    !!   - sky: Sky to compute power of.
    !!   - power: Full power computed.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2006
    !
    ! Revisions:
    !   May 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_power_map(sky) result(power)

      type(s2_sky), intent(in) :: sky
      real(s2_sp) :: power

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_power_map')
      end if

      ! Check map defined.
      ! If map not computed stop since lmax and mmax may not be defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_power_map')
      end if

      power = sum(sky%map(0:sky%npix-1)**2) * 4*pi/real(sky%npix,s2_sp)

    end function s2_sky_power_map

    
    !--------------------------------------------------------------------------
    ! s2_sky_power_alm
    !
    !! Compute power of function on sphere in alm space.
    !!
    !! Variables:
    !!   - sky: Sky to compute power of.
    !!   - power: Full power computed.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2006
    !
    ! Revisions:
    !   May 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_power_alm(sky) result(power)

      type(s2_sky), intent(in) :: sky
      real(s2_sp) :: power

      type(s2_pl) :: cl
      integer :: l, m

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_power_alm')
      end if

      ! Check alm defined.
      ! If alm not computed stop since lmax and mmax may not be defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_power_alm')
      end if

!      ! Compute power directly.
!      power = 0e0
!      do l = 0,sky%lmax
!
!         power = power + abs(sky%alm(l,0))**2
!
!         ! Note m starting from 1 in summation.
!         do m = 1,min(l, sky%mmax)
!            power = power + 2e0*abs(sky%alm(l,m))**2.0e0
!         end do
!
!      end do

      ! Compute power from cl.      
      cl = s2_sky_get_cl(sky)
      power = s2_pl_power(cl)

      ! Free memory.
      call s2_pl_free(cl)

    end function s2_sky_power_alm


    !--------------------------------------------------------------------------
    ! s2_sky_azimuthal_bl
    !
    !! Find azimuthal band limit of sky.  Finds lowest m' value such that 
    !! cutoff_prop*100 percent of the cm power is contained in the alms 
    !! with m index below m'.
    !!
    !! Notes:
    !!   - Only successful if alms for the sky are already computed.  Cannot 
    !      automatically compute alms if not already computed since lmax and 
    !!     mmax may not be defined.
    !!
    !! Variables
    !!  - sky: Sky to find azimuthal band limit for.
    !!  - [cutoff_prop_in]: Cutoff proportion to use when finding azimuthal 
    !!    band limit.  If not specified default value is used.
    !!  - mmax_min: The azimuthal band limit found.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_azimuthal_bl(sky, cutoff_prop_in) result(mmax_min)

      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(in), optional :: cutoff_prop_in
      integer :: mmax_min

      type(s2_pl) :: cl
      real(s2_sp), allocatable :: cl_vals(:), cm_vals(:), cm_cummulative(:)
      real(s2_sp) :: cutoff, cutoff_prop
      integer :: m, fail
      real(s2_sp), parameter :: CUTOFF_PROP_DEFAULT = 0.95e0

      if(present(cutoff_prop_in)) then
         ! Check cutoff proportion passed is in valid range.
         if(cutoff_prop_in <= 0.0e0 .or. cutoff_prop_in > 1.0e0) then
            call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_azimuthal_bl', &
              comment_add='Invalid cutoff proportion')
         end if
         cutoff_prop = cutoff_prop_in
      else
         cutoff_prop = CUTOFF_PROP_DEFAULT
      end if

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_azimuthal_bl')
      end if

      ! Check alm defined.
      ! If alm not computed stop since lmax and mmax may not be defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_azimuthal_bl')
      end if

      ! Allocate space for cl and cm spectra.
      allocate(cm_vals(0:sky%mmax), stat=fail)
      allocate(cm_cummulative(0:sky%mmax), stat=fail)
      allocate(cl_vals(0:sky%lmax), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_azimuthal_bl')
      end if

      ! Compute cls and cms.
      cl = s2_sky_get_cl(sky)
      call s2_pl_get_spec(cl, cl_vals)     
      call s2_sky_get_cm(sky, cm_vals(0:sky%mmax))

      ! Compute cummulative cm array,
      cm_cummulative = cm_vals
      cm_cummulative(0) = (2*sky%lmax+1) * cm_cummulative(0)
      do m = 1,sky%mmax
!         cm_cummulative(m) = sum(cm_vals(0:m))    !must start m index from 0 if use this technique
         cm_cummulative(m) = (2*(sky%lmax-m)+2) * cm_cummulative(m) &
           + cm_cummulative(m-1)
      end do

      ! Set cutoff as proportion of total sum of cms.
      cutoff = cm_cummulative(sky%mmax) * cutoff_prop
      
      ! Step through m's and find first index in cumulative that has required 
      ! portion below.
      do m = 0,sky%mmax
         mmax_min = m
         if(cm_cummulative(m) >= cutoff) exit   ! Stop loop, found mmax_min.
      end do

      ! Free memory.
      call s2_pl_free(cl)
      deallocate(cl_vals, cm_vals, cm_cummulative)

    end function s2_sky_azimuthal_bl


    !--------------------------------------------------------------------------
    ! s2_sky_admiss
    !    
    !! Calculte the normalised numerical admissibility of a wavelet defined 
    !! on the sphere.
    !!
    !! Notes:
    !!   - To developers: Have a consistent area element in every term of
    !!     discretised integral thus it cancels and may be removed.
    !!
    !! Variables:
    !!   - sky: Sky to calculate wavelet admissibility of.
    !!   - admiss: Normalised numerical admissibility value calculated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_admiss(sky) result(admiss)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      type(s2_sky), intent(in) :: sky
      real(s2_sp) :: admiss

      integer :: ipix, fail
      real(s2_sp) :: c, c_abswav
      real(s2_dp) :: phi_unused
      real(s2_dp), allocatable :: theta(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_admiss')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_admiss')
      end if

      allocate(theta(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_admiss')
      end if

      ! Evaluate theta for each sample over the sphere.
      do ipix = 0,sky%npix-1

         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta(ipix), phi_unused)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta(ipix), phi_unused)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_rotate')
         end if

      end do
      
      c = sum(sky%map/(1.0e0+cos(real(theta,s2_sp))))
      c_abswav = sum(abs(sky%map)/(1.0e0+cos(real(theta,s2_sp))))
      
      admiss = c / c_abswav

      deallocate(theta)

    end function s2_sky_admiss
    

    !--------------------------------------------------------------------------
    ! s2_sky_admiss_dil
    !
    !! Compute admissibility for function defined on sphere for a range of
    !! dilations.
    !!
    !! Variables:
    !!   - sky: Original (undilated) sky to calculated dilations of.
    !!   - dilation1: Array of first dilation components.
    !!   - dilation2: Array of second dilation components.
    !!   - admiss: Array of admissibility values calculated corresponding to
    !!     the dilated version of the passed function defined on the sphere.
    !!   - [norm_preserve_in]: Logical to specify whether dilations are to 
    !!     preserve the 2-norm.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_admiss_dil(sky, dilation1, dilation2, admiss, &
      norm_preserve_in)

      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(in) :: dilation1(:)
      real(s2_sp), intent(in) :: dilation2(:)
      real(s2_sp), intent(out) :: admiss(:)
      logical, intent(in), optional :: norm_preserve_in

      logical, parameter :: NORM_PRESERVE_DEFAULT = .true.
      integer :: idil, ndil
      type(s2_sky) :: temp_sky
      logical :: norm_preserve

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_admiss_dil')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_admiss_dil')
      end if

      ndil = size(dilation1)

      ! Check dilation arrays consistent.
      if(ndil /= size(dilation2)) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_admiss_dil', &
           comment_add='Dilation components inconsistent size')
      end if

      ! Check admiss size consistent.
      if(ndil /= size(admiss)) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_admiss_dil', &
           comment_add='Admissibility array inconsistent size')
      end if

      ! Set norm_preserve.
      if(present(norm_preserve_in)) then
         norm_preserve = norm_preserve_in
      else
         norm_preserve = NORM_PRESERVE_DEFAULT
      end if

      ! Calculate admissibility for each dilation.
      do idil = 1,ndil

         temp_sky = s2_sky_init(sky)
         call s2_sky_dilate(temp_sky, dilation1(idil), dilation2(idil), &
           norm_preserve)    
         admiss(idil) = s2_sky_admiss(temp_sky)
         call s2_sky_free(temp_sky)

      end do

    end subroutine s2_sky_admiss_dil


    !--------------------------------------------------------------------------
    ! s2_sky_fov
    !
    !! Zero portion of sky outside of field-of-view (fov).  Two methods exist: 
    !! (i) S2_SKY_FOV_METHOD_CIRCLE: restrict fov based on theta_fiv only;
    !! (ii) S2_SKY_FOV_METHOD_SQUARE: restrict fov to ensure only planar grid 
    !! within theta_fov remains.
    !!
    !! Variables:
    !!   - sky: Sky to restrict fov.
    !!   - theta_fov: Fov to restruct to (note theta_max = theta_fov/2).
    !!   - method: Method to use when determining fov.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_fov(sky, theta_fov, method)

      use pix_tools, only: pix2ang_ring, pix2ang_nest
      use s2_vect_mod

      type(s2_sky), intent(inout) :: sky
      real(s2_dp), intent(in) :: theta_fov
      integer, intent(in) :: method

      integer :: ipix
      real(s2_dp) :: theta, phi
      real(s2_sp) :: x(3)
      real(s2_sp) :: L
      type(s2_vect) :: vec

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_fov')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_fov')
      end if
     
      ! Zero maps values outside fov.
      L = sqrt(2.0) * sin(theta_fov/2.0)
      do ipix = 0,sky%npix-1

         ! Get theta and phi angles corresponding to pixel.
         if(sky%pix_scheme == S2_SKY_RING) then
            call pix2ang_ring(sky%nside, ipix, theta, phi)
         else if(sky%pix_scheme == S2_SKY_NEST) then
            call pix2ang_nest(sky%nside, ipix, theta, phi)
         else
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_fov')
         end if

         ! Zero outside fov.
         select case(method)

             case(S2_SKY_FOV_METHOD_CIRCLE)
                if (theta > theta_fov/2.0) sky%map(ipix) = 0e0

             case(S2_SKY_FOV_METHOD_SQUARE)
                vec = s2_vect_init(1.0, real(theta,s2_sp), real(phi,s2_sp))
                call s2_vect_convert(vec, S2_VECT_TYPE_CART)
                x = s2_vect_get_x(vec)
                call s2_vect_free(vec)
                if( (abs(x(1)) > L/2.0) .or. (abs(x(2)) > L/2.0) ) then
                   sky%map(ipix) = 0e0
                end if

             case default
                call s2_error(S2_ERROR_SKY_FOV_METHOD_INVALID, 's2_sky_fov', &
                     comment_add='Invalid fov method type specifier')
            
         end select

      end do

      ! Ensure alm not defined.
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.

    end subroutine s2_sky_fov


    !--------------------------------------------------------------------------
    ! s2_sky_extract_ab_fsht
    !
    !! Extra an ecp (equispaced) sampled theta-phi array over the sphere
    !! for the grid used for the Fast Spherical Harmonic Transform.
    !!
    !! Notes:
    !!   - No interpolation performed.
    !!   - The size of xtp must already be specified and allocated.
    !!
    !! Variables:
    !!   - sky: Sky to extract theta-phi array representation of.
    !!   - xtp: The extracted discrete theta-phi data array corresponding
    !!     to the given sky.
    !  
    !! @author J. D. McEwen
    !! @version 0.1 - April 2008
    !
    ! Revisions:
    !   April 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_extract_ab_fsht(sky, xtp, L)

      use pix_tools, only: ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(out) :: xtp(0:L,0:2*L)
      integer, intent(in) :: L

      integer :: itheta, iphi, ipix
      real(s2_dp) :: theta, phi

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_extract_ab_fsht')
      end if

      do itheta = 0,L

         theta = 2*pi*(2.*itheta+1)/real(2*(2*L+1),s2_dp)
         theta = mod(theta, PI)
         if(theta>=3.1415926) theta = 3.1415926

         do iphi = 0,2*L

            phi = 2*pi*(2.*iphi+1)/real(2*(2*L+1),s2_dp)
            phi = mod(phi, 2*PI)

            ! Compute index corresponding to theta (beta) and phi (alpha)
            ! angles.
            if(sky%pix_scheme == S2_SKY_RING) then
               call ang2pix_ring(sky%nside, theta, phi, ipix)
            else if(sky%pix_scheme == S2_SKY_NEST) then
               call ang2pix_nest(sky%nside, theta, phi, ipix)
            else
               call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_extract_ab_fsht')
            end if
            
            xtp(itheta, iphi) = sky%map(ipix)
write(*,*) 'xtp(', itheta+1, ',', iphi+1, ') = ', xtp(itheta, iphi), ';'

         end do
         
      end do
      
    end subroutine s2_sky_extract_ab_fsht


    !--------------------------------------------------------------------------
    ! s2_sky_extract_ab_s2dw
    !
    !! Extra an ecp (equispaced) sampled theta-phi array over the sphere
    !! for the grid used for S2DW.
    !!
    !! Notes:
    !!   - No interpolation performed.
    !!   - The size of xtp must already be specified and allocated.
    !!
    !! Variables:
    !!   - sky: Sky to extract theta-phi array representation of.
    !!   - xtp: The extracted discrete theta-phi data array corresponding
    !!     to the given sky.
    !  
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_extract_ab_s2dw(sky, xtp, B)

      use pix_tools, only: ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(in) :: sky
      real(s2_dp), intent(out) :: xtp(0:2*B-1,0:2*B-2)
      integer, intent(in) :: B

      integer :: itheta, iphi, ipix
      real(s2_dp) :: theta, phi

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_extract_ab_fsht')
      end if

      do itheta = 0,2*B-1

         theta = pi*(2*itheta+1)/real(4*B,s2_dp)
         theta = mod(theta, PI)

         do iphi = 0,2*B-2

            phi = 2*pi*iphi/real(2*B-1,s2_dp)
            phi = mod(phi, 2*PI)

            ! Compute index corresponding to theta (beta) and phi (alpha)
            ! angles.
            if(sky%pix_scheme == S2_SKY_RING) then
               call ang2pix_ring(sky%nside, theta, phi, ipix)
            else if(sky%pix_scheme == S2_SKY_NEST) then
               call ang2pix_nest(sky%nside, theta, phi, ipix)
            else
               call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_extract_ab_fsht')
            end if
            
            xtp(itheta, iphi) = sky%map(ipix)

         end do
         
      end do
      
    end subroutine s2_sky_extract_ab_s2dw


    !--------------------------------------------------------------------------
    ! s2_sky_extract_ab
    !
    !! Extra an ecp (equispaced) sampled alpha-beta array over the sphere
    !! corresponding to the given sky.
    !!
    !! Notes:
    !!   - No interpolation performed.
    !!   - The size of x_ab must already be specified and allocated.
    !!
    !! Variables:
    !!   - sky: Sky to extract alpha-beta array representation of.
    !!   - x_ab: The extracted discrete alpha-beta data array corresponding
    !!     to the given sky.
    !  
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   March 2004 - Written by Jason McEwen for cswt-0.1
    !   November 2004 - Incorporated in s2_sky_mod
    !--------------------------------------------------------------------------

    subroutine s2_sky_extract_ab(sky, x_ab)

      use pix_tools, only: ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(out) :: x_ab(:,:)

      integer :: n_alpha, n_beta, i_alpha, i_beta, ipix
      real(s2_dp) :: alpha, beta

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_extract_ab')
      end if

      n_alpha = size(x_ab, 1)
      n_beta = size(x_ab, 2)

      do i_alpha = 1,n_alpha
         alpha = (i_alpha-1) * 2.0d0*pi/real(n_alpha, s2_dp)

         do i_beta = 1,n_beta
            beta = (i_beta-1) * pi / real(n_beta, s2_dp)
            
            ! Compute index corresponding to theta (beta) and phi (alpha)
            ! angles.
            if(sky%pix_scheme == S2_SKY_RING) then
               call ang2pix_ring(sky%nside, beta, alpha, ipix)
            else if(sky%pix_scheme == S2_SKY_NEST) then
               call ang2pix_nest(sky%nside, beta, alpha, ipix)
            else
               call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_extract_ab')
            end if
            
            x_ab(i_alpha, i_beta) = sky%map(ipix)
            
         end do
         
      end do
      
    end subroutine s2_sky_extract_ab


    !--------------------------------------------------------------------------
    ! s2_sky_inter_ab
    !
    !! Interpolate a pixelised ab array for any continuous value 
    !! in the appropriate angle ranges.
    !!
    !! Notes:
    !!   - Linear interpolation used.
    !!   - Angles outside the valid range for alpha [0,2*pi] and
    !!     beta [0,pi] are mapped onto the mod of 
    !!     the angle for the relevant range.
    !!
    !! Variables:
    !!   - x_ab: The discrete ab data array to be
    !!     interpolated. 
    !!   - alpha: The alpha continuous angle to interpolate at.
    !!   - beta: The beta continuous angle to interpolate at.
    !!   - interp: Logical to specify whether interpolation is to be 
    !!     performed.  If not, then the nearest (lower) sample is simply taken.
    !!   - [sdw]: Logical to specify SDW pixelisation.
    !  
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   March 2004 - Written by Jason McEwen for cswt-0.1
    !   November 2004 - Incorporated in s2_sky_mod
    !--------------------------------------------------------------------------

    function s2_sky_interp_ab(x_ab, alpha, beta, interp, sdw) &
      result(interp_val)
      
      real(s2_sp), intent(in) :: x_ab(:,:)
      real(s2_dp), intent(inout) :: alpha, beta
      logical, intent(in) :: interp
			logical, intent(in), optional :: sdw
      real(s2_sp) :: interp_val
    
      integer :: n_a, n_b
      integer :: i_a, i_b
      real(s2_sp) :: i_a_interp, i_b_interp
      
      n_a = size(x_ab,1)
      n_b = size(x_ab,2)
      
      ! Ensure angles in valid range.
      alpha = mod(alpha, 2.0d0*pi)
      beta = mod(beta, real(pi,s2_dp))
      
      ! Find bounding indices.
      i_a_interp = alpha * n_a / (2.0e0 * pi)
      i_a = floor(i_a_interp)
      
			if(present(sdw) .and. sdw) then
				i_b_interp = max(beta * real(n_b,s2_sp) / pi - 0d5, 0d0)
			else
      	i_b_interp = beta * n_b / pi
			end if
      i_b = floor(i_b_interp)
      

      if(interp) then
      
        ! Calculate interpolated value.  (See log book p26-28 for derivation of
        ! formula.)
        ! Note extra plus 1 in x_abg indicies since index 1:n.

        ! Ensure index at most second to last sample so interpolation 
        ! valid.  If i_x_interp between n_x-1 and n_x then interpolation 
        ! rule works for extrapolation.    
        if(i_a == n_a-1 .and. n_a > 1) i_a = n_a-2
        if(i_b == n_b-1 .and. n_b > 1) i_b = n_b-2
      
        interp_val = &
          !
          & x_ab(i_a+1+1, i_b+1+1) &
          & * (i_a_interp - i_a) * (i_b_interp - i_b) &
          !
          & - x_ab(i_a+1, i_b+1+1) &
          & * (i_a_interp - (i_a+1)) * (i_b_interp - i_b) &
          !
          & - x_ab(i_a+1+1, i_b+1) &
          & * (i_a_interp - i_a) * (i_b_interp - (i_b+1)) &
          !
          & + x_ab(i_a+1, i_b+1) &
          & * (i_a_interp - (i_a+1)) * (i_b_interp - (i_b+1))
    
      else
      
        ! Use closest sample instead of interpolating.
        interp_val = x_ab(i_a+1, i_b+1)

      end if
    
    end function s2_sky_interp_ab
    

    !--------------------------------------------------------------------------
    ! s2_sky_downsample
    !      
    !! Downsample a sky map nside to the specified (lower!) nside.
    !!
    !! Variables:
    !!   - sky: Sky to be downsampled. 
    !!   - nside_down: The new nside of the downsampled sky.
    !!   - mask: Logical specifying whether sky to downsample is a mask 
    !!     (in which case pixels are multiplied to downsample rather than
    !!     averaged).
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_downsample(sky, nside_down, mask)

      use pix_tools, only: nside2npix

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: nside_down
      logical, optional, intent(in) :: mask

      integer :: npix_down

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_downsample')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_downsample')
         return
      end if

      ! Check nside_down valid size.
      if(nside_down > sky%nside) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_downsample', &
           comment_add='nside_down larger than original size')
      end if
      npix_down = nside2npix(nside_down) ! Also checks nside_down power of 2.

      ! Downsample by two until reach required size.
      do
         if(sky%nside <= nside_down) exit   ! Finished downsampling.
         call s2_sky_down_nsideby2(sky, mask)
      end do

      ! Ensure alm not present.
      ! (Should already be removed each time downsample nside by 2 but
      ! remove here as precaution.)
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.

    end subroutine s2_sky_downsample


    !--------------------------------------------------------------------------
    ! s2_sky_down_nsideby2
    !
    !! Downsample a sky map nside by a factor of two.
    !!
    !! Variables:
    !!   - sky: Sky to be downsampled. 
    !!   - mask_in: Logical specifying whether sky to downsample is a mask 
    !!     (in which case pixels are multiplied to downsample rather than
    !!     averaged).
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_down_nsideby2(sky, mask_in)

      use pix_tools, only: nside2npix

      type(s2_sky), intent(inout) :: sky
      logical, optional, intent(in) :: mask_in

      integer :: nside_down, npix_down, ipix, ipix_down, fail, pix_scheme_orig
      real(s2_sp), allocatable :: map_down(:)
      logical :: mask 

      mask = .false.
      if(present(mask_in)) mask = mask_in

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_down_nsideby2')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_down_nsideby2')
         return
      end if

      ! Define new sizes.
      nside_down = sky%nside / 2
      npix_down = nside2npix(nside_down)

     ! Allocate temporary storage space for dilated map.
      allocate(map_down(0:npix_down-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_down_nsideby2')
      end if
      
      ! Initialise new map with zeros.

      if(mask) then
         map_down = 1.0e0
      else
         map_down = 0.0e0
      end if

      ! Ensure in nested pixelisation scheme.
      pix_scheme_orig = sky%pix_scheme
      call s2_sky_map_convert(sky, S2_SKY_NEST)

      ! Create downsampled map.
      do ipix = 0,sky%npix-1
        ipix_down = ipix/4    ! Integer division.
        if(mask) then
           map_down(ipix_down) = map_down(ipix_down) * sky%map(ipix)
        else
           map_down(ipix_down) = map_down(ipix_down) + sky%map(ipix)
        end if
      end do
      if(.not. mask) map_down = map_down / 4.0e0

      ! Update sky.
      sky%nside = nside_down
      sky%npix = npix_down
      deallocate(sky%map)
      allocate(sky%map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_down_nsideby2')
      end if
      sky%map = map_down

      ! Convert back to original pixelisation scheme.
      call s2_sky_map_convert(sky, pix_scheme_orig)

      ! Free temporary storage space used.
      deallocate(map_down)

      ! Ensure alm not present.
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.

    end subroutine s2_sky_down_nsideby2
   

    !--------------------------------------------------------------------------
    ! s2_sky_upsample
    !      
    !! Up-sample a sky map nside to the specified (greater!) nside.
    !! Map is up-sampled by simply setting all sub pixels in a nested region
    !! to the value of the parent pixel).
    !!
    !! Variables:
    !!   - sky: Sky to be upsampled.
    !!   - nside_down: The new nside of the downsampled sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2005
    !
    ! Revisions:
    !   August 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_upsample(sky, nside_up)

      use pix_tools, only: nside2npix

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: nside_up

      integer :: npix_up

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_upsample')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_upsample')
         return
      end if

      ! Check nside_down valid size.
      if(nside_up < sky%nside) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_upsample', &
           comment_add='nside_up smaller than original size')
      end if
      npix_up = nside2npix(nside_up) ! Also checks nside_up power of 2.

      ! Upsample by two until reach required size.
      do
         if(sky%nside >= nside_up) exit   ! Finished up-sampling.
         call s2_sky_up_nsideby2(sky)
      end do

      ! Ensure alm not present.
      ! (Should already be removed each time upsample nside by 2 but
      ! remove here as precaution.)
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.

    end subroutine s2_sky_upsample


    !--------------------------------------------------------------------------
    ! s2_sky_up_nsideby2
    !
    !! Up-sample a sky map nside by a factor of two (simply set all sub
    !! pixels in a nested region to the value of the parent pixel).
    !!
    !! Variables:
    !!   - sky: Sky to be up-sampled.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2005
    !
    ! Revisions:
    !   August 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_up_nsideby2(sky)

      use pix_tools, only: nside2npix

      type(s2_sky), intent(inout) :: sky

      integer :: nside_up, npix_up, ipix, ipix_up, fail, pix_scheme_orig
      real(s2_sp), allocatable :: map_up(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_up_nsideby2')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_up_nsideby2')
         return
      end if

      ! Define new sizes.
      nside_up = sky%nside * 2
      npix_up = nside2npix(nside_up)

     ! Allocate temporary storage space for dilated map.
      allocate(map_up(0:npix_up-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_up_nsideby2')
      end if
      
      ! Initialise new map with zeros.
      map_up = 0.0e0
      
      ! Ensure in nested pixelisation scheme.
      pix_scheme_orig = sky%pix_scheme
      call s2_sky_map_convert(sky, S2_SKY_NEST)

      ! Create up-sampled map.

      do ipix = 0,sky%npix-1
        ipix_up = 4*ipix
        map_up(ipix_up:ipix_up+3) = sky%map(ipix)
      end do

      ! Update sky.
      sky%nside = nside_up
      sky%npix = npix_up
      deallocate(sky%map)
      allocate(sky%map(0:sky%npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_up_nsideby2')
      end if
      sky%map(0:sky%npix-1) = map_up(0:sky%npix-1)

      ! Convert back to original pixelisation scheme.
      call s2_sky_map_convert(sky, pix_scheme_orig)

      ! Free temporary storage space used.
      deallocate(map_up)

      ! Ensure alm not present.
      if(allocated(sky%alm)) deallocate(sky%alm)
      sky%alm_status = .false.
      
    end subroutine s2_sky_up_nsideby2


    !--------------------------------------------------------------------------
    ! s2_sky_valid_sizes
    !
    !! Check the nside, lmax and mmax resolutions are valid.  Invalid if any
    !! are negative.  Also warning given if lmax > 3*nside of mmax < lmax
    !!
    !! Variables:
    !!   - sky: Sky to check sizes of.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_valid_sizes(sky)
      
      type(s2_sky), intent(in) :: sky
      
      if(sky%nside < 0 .or. sky%lmax < 0 .or. sky%mmax < 0) then
        call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_valid_sizes')
      end if

      if(sky%lmax > 3*sky%nside .or. sky%mmax < sky%lmax) then
         call s2_error(S2_ERROR_SKY_SIZE_WARNING, 's2_sky_valid_sizes')
      end if

    end subroutine s2_sky_valid_sizes


    !--------------------------------------------------------------------------
    ! s2_sky_draw_dot
    !
    !! Draw dots on sky map at specified positions.
    !!
    !! Variables:
    !!   - sky: Sky to draw dots on map of.
    !!   - alpha: Alpha array of dot positions.
    !!   - beta: Beta array of dot positions.  (Must have same size as alpha 
    !!     array.)
    !!   - [large]: Logical to indicate whether to draw large dots.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_draw_dot(sky, alpha, beta, large)

      use pix_tools, only: ang2pix_nest, neighbours_nest

      type(s2_sky), intent(inout) :: sky
      real(s2_sp), intent(in) :: alpha(:)
      real(s2_sp), intent(in) :: beta(:)
      logical, intent(in), optional :: large

      real(s2_sp), parameter :: DOT_VAL = -1.6375e30
      integer :: orig_pix_scheme, i, j, ipix
      integer :: neigh(8), n_neigh      
      integer :: neigh2(8), n_neigh2      

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_draw_dot')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_draw_dot')
         return
      end if

      ! Check size of alpha and beta arrays are consistent.
      if(size(alpha) /= size(beta)) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_draw_dot', &
           comment_add='Size of alpha and beta arrays is inconsistent')
      end if

      ! Get original pixelisation scheme.
      orig_pix_scheme = sky%pix_scheme

      ! Make sure in nested pixelisation scheme for now.
      call s2_sky_map_convert(sky, S2_SKY_NEST)

      ! Draw on each dot.
      do i = 1,size(alpha)

         call ang2pix_nest(sky%nside, real(beta(i),s2_dp), &
           real(alpha(i),s2_dp), ipix)

         call neighbours_nest(sky%nside, ipix, neigh, n_neigh)

         sky%map(neigh(1:n_neigh)) = DOT_VAL
         sky%map(ipix) = DOT_VAL

         ! If require large dot find second level neighbours.
         if(present(large) .and. large) then
            do j = 1,n_neigh
               call neighbours_nest(sky%nside, neigh(j), neigh2, n_neigh2)
               sky%map(neigh2(1:n_neigh2)) = DOT_VAL
            end do
         end if

      end do

      ! Return to original pixelisation scheme.
      call s2_sky_map_convert(sky, orig_pix_scheme)

    end subroutine s2_sky_draw_dot


    !--------------------------------------------------------------------------
    ! IO routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_sky_init_file
    ! 
    !! Wrapper to write a s2 sky various types of files.
    !!
    !! Variables:
    !!   - sky: Sky containing the map and/or alms to write to a file.
    !!   - filename: Name of the output file.
    !!   - file_type: Fit type specifier to specify whether to write a 
    !!     fits map file, fits alm file or a fits full s2_sky file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_write_file(sky, filename, file_type, comment)

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_type
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_write_file')
      end if

      ! Write file.
      select case(file_type)

         case(S2_SKY_FILE_TYPE_SKY)
            call s2_sky_io_fits_write(filename, sky, comment)

         case(S2_SKY_FILE_TYPE_MAP)
            call s2_sky_write_map_file(sky, filename, comment)

         case(S2_SKY_FILE_TYPE_ALM)
            call s2_sky_write_alm_file(sky, filename, comment)

         case default
            call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky_write_file', &
              comment_add='Invalid file type specifier')

      end select

    end subroutine s2_sky_write_file


    !--------------------------------------------------------------------------
    ! s2_sky_write_map_file
    !
    !! Write a sky map to a fits file.
    !!
    !! Variables:
    !!   - sky: Sky containing the map to write to a fits file.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_write_map_file(sky, filename, comment)

      use head_fits, only: add_card
      use fitstools, only: write_bintab

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer, parameter :: HEADER_LEN = 180
      character(len=80) :: header(HEADER_LEN)
      real(s2_sp), allocatable :: map_temp(:,:)
      integer :: i, fail

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_write_map')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_write_map')
      end if

      ! Initialise empty header.
      do i = 1,HEADER_LEN
        header(i) = ""
      end do

      ! Write header.

      call add_card(header, "COMMENT", &
        "-----------------------------------------------")
      call add_card(header, "COMMENT", &
        "     Sky Map Pixelisation Specific Keywords    ")
      call add_card(header, "COMMENT", &
        "-----------------------------------------------")

      call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")

      if(sky%pix_scheme == S2_SKY_RING) then
        call add_card(header, "ORDERING", "RING",&
          "Pixel ordering scheme, either RING or NESTED")
      else if(sky%pix_scheme == S2_SKY_NEST) then
        call add_card(header, "ORDERING", "NESTED",&
          "Pixel ordering scheme, either RING or NESTED")
      else
        call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_sky_write_map_file')
      end if

      call add_card(header, "NSIDE" ,sky%nside, &
        "Resolution parameter for HEALPIX")
      call add_card(header, "FIRSTPIX", 0, "First pixel # (0 based)")
      call add_card(header, "LASTPIX", sky%npix-1, "Last pixel # (0 based)")

      call add_card(header) ! Blank line

      call add_card(header, "COMMENT", &
        "-----------------------------------------------")
      call add_card(header, "COMMENT", &
        "      S2 Specific Keywords/Comments           ")
      call add_card(header, "COMMENT", &
        "-----------------------------------------------")  

      if(present(comment)) call add_card(header, "COMMENT", trim(comment)) 

      ! If wish to add other sky attributes to file in future do so
      ! here, e.g. lmax.

      call add_card(header) ! blank line

      call add_card(header,"COMMENT","*************************************")

      ! Write map.
      allocate(map_temp(0:sky%npix-1,1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_write_map_file')
      end if

      ! Copy map to temporary storage array of correct dimension for
      ! writing via healpix routines.
      map_temp(:,1) = sky%map(0:sky%npix-1)
      call write_bintab(map_temp, sky%npix, 1, header, HEADER_LEN, filename)
      deallocate(map_temp)
      
    end subroutine s2_sky_write_map_file


    !--------------------------------------------------------------------------
    ! s2_sky_write_alm_file
    !
    !! Write sky alms to a fits file.
    !!
    !! Variables:
    !!   - sky: Sky containing the alm to write to a fits file.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2005
    !
    ! Revisions:
    !   November 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_write_alm_file(sky, filename, comment)

      use head_fits, only: add_card
      use fitstools, only: alms2fits

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: nalm, ncol, next, fail
      integer :: l, m, i
      real(s2_sp), allocatable :: alm_output(:,:,:)
      integer, parameter :: HEADER_LEN = 180
      character(len=80) :: header(HEADER_LEN,1)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_write_alm_file')
      end if

      ! Check alms defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_write_alm_file')
      end if

      ! Define sizes.
      nalm = (sky%lmax+1)*(sky%mmax+2)/2
      ncol = 3
      next = 1

      ! Allocate space for output alm array.
      allocate(alm_output(1:nalm, 1:(ncol+1), 1:next), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_write_alm_file')
      end if

      ! Format alm array for output.
      i = 1
      do l = 0,sky%lmax
         do m = 0,min(l,sky%mmax)
            alm_output(i,1,1) = l
            alm_output(i,2,1) = m
            alm_output(i,3,1) = real(sky%alm(l,m), s2_sp)
            alm_output(i,4,1) = real(aimag(sky%alm(l,m)), s2_sp)
            i = i + 1
         end do
      end do
      i = i-1

      ! Initialise empty header.
      do i = 1,HEADER_LEN
        header(i,1) = ""
      end do

      ! Write header.
      call add_card(header(:,1)) ! Blank line
      call add_card(header(:,1), "COMMENT", &
        "-----------------------------------------------")
      call add_card(header(:,1), "COMMENT", &
        "      S2 Specific Keywords/Comments           ")
      call add_card(header(:,1), "COMMENT", &
        "-----------------------------------------------")  

      call add_card(header(:,1), "LMAX" ,sky%lmax, &
        "Maximum harmonic l of alms")
      call add_card(header(:,1), "MMAX" ,sky%mmax, &
        "Maximum harmonic m of alms")

      if(present(comment)) then
         call add_card(header(:,1), "COMMENT", trim(comment)) 
      end if

      ! If wish to add other sky attributes to file in future do so
      ! here, e.g. lmax.

      call add_card(header(:,1)) ! blank line
      call add_card(header(:,1),"COMMENT", &
        "***********************************************")

      ! Write alms to output file.
      call alms2fits(filename, nalm, alm_output, ncol, &
        header, HEADER_LEN, next)

      ! Free memory.
      deallocate(alm_output)

    end subroutine s2_sky_write_alm_file


    !--------------------------------------------------------------------------
    ! s2_sky_write_matmap_file
    !
    !! Write a sky map to a matlab map file.  Note that the matlab map is
    !! defined on an equiangular grid.
    !!
    !! Variables:
    !!   - sky: Sky containing the map to write to a file.
    !!   - filename: Name of the output matlab map file.
    !!   - B: Band limit corresponding to equi-angular grid.
    !!   - [comment]: Optional additional comment to be added to the file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_write_matmap_file(sky, filename, B, comment)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in) :: filename
      integer, intent(in) :: B
      character(len=*), intent(in), optional :: comment

      real(s2_dp) :: theta, phi
      integer :: fileid, ipix, itheta, iphi
      integer :: fail = 0
      real(s2_dp), allocatable :: xtp(:,:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_write_matmap_file')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_write_matmap_file')
      end if

      ! Extract equi-angular sampled sphere to write to file.
      allocate(xtp(0:2*B-1,0:2*B-2), stat=fail)
      if(fail /= 0) then 
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
              's2_sky_write_matmap_file')
      end if
      call s2_sky_extract_ab_s2dw(sky, xtp, B)

      ! Open file.
      fileid = 11
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')

      ! Write to file.
      if(present(comment)) write(fileid,'(a,a)') '% ', trim(comment)
      do itheta = 0,2*B-1
         theta = pi*(2*itheta+1)/real(4*B,s2_dp)
         theta = mod(theta, PI)
         do iphi = 0,2*B-2
            phi = 2*pi*iphi/real(2*B-1,s2_dp)
            phi = mod(phi, 2*PI)
            write(fileid,'(3e28.20)') theta, phi, xtp(itheta, iphi)
         end do
      end do

      ! Close file.
      close(fileid)

      ! Free memory.
      deallocate(xtp)

    end subroutine s2_sky_write_matmap_file


    !--------------------------------------------------------------------------
    ! s2_sky_write_matalm_file
    !
    !! Write alms to a matlab alm file.  
    !!
    !! Variables:
    !!   - sky: Sky containing the alms to write to a file.
    !!   - filename: Name of the output matlab alm file.
    !!   - [comment]: Optional additional comment to be added to the file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_write_matalm_file(sky, filename, comment)

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: fileid, el, m

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_write_matalm_file')
      end if

      ! Check alms defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_write_matalm_file')
      end if

      ! Open file.
      fileid = 12
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')

      ! Write to file.
      if(present(comment)) write(fileid,'(a,a)') '% ', trim(comment)
      do m = 0,sky%lmax
         do el = 0,sky%lmax
            if (m <= sky%mmax) then
               write(fileid,'(2e28.20)') &
                    real(sky%alm(el,m), s2_dp), real(aimag(sky%alm(el,m)), s2_dp)
            else
               write(fileid,'(2e28.20)') real(0.0, s2_dp), real(0.0, s2_dp)
            end if
         end do
      end do

      ! Close file.
      close(fileid)

    end subroutine s2_sky_write_matalm_file


    !--------------------------------------------------------------------------
    ! s2_sky_io_fits_write
    !
    !! Write a s2_sky object to a fits file.  All s2_sky data is writtend to
    !! the file (i.e. both the map and alms are written if present, as well as
    !! any sky function parameters and all other s2_sky variables.
    !!
    !! Variables:
    !!   - filename: Name of the output fits file to write the sky structure
    !!     data to.
    !!   - sky: The sky structure containing the data to be written to the
    !!     output fits file.
    !!   - [comment]: Optional comment string to be added to the output fits 
    !!     file header if present.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_io_fits_write(filename, sky, comment)

      character(len=*), intent(in) :: filename
      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      integer :: group,dim1
      logical :: simple, extend, file_exists
      integer :: naxis
      integer :: naxes(1), naxes_alm(2)
      integer :: tfields, nrows, varidat
      character(len=32) :: ttype(1), tform(1), tunit(1), extname
      integer :: frow, felem, colnum

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_io_fits_write')
      end if

      ! Check at least one of map or alms defiend.
      ! This error should never be able to occur, since for sky to be 
      ! initialised either a map or alm must be defined but check regardless.
      if(.not. sky%map_status .and. .not. sky%alm_status) then
         call s2_error(S2_ERROR_NOT_INIT, 's2_sky_io_fits_write', &
           comment_add='Neither map or alm defined')
      end if

      ! Define FITS parameters.

      bitpix=-32 ! Real single precision.
      status=0   ! Initialse error status to zero.

      ! Check if file already exists.
      call s2_sky_io_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call s2_error(S2_ERROR_SKY_FILE_EXISTS, &
              's2_sky_io_fits_write')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.
      call ftpcom(unit, &
        '  Sph_sky data file created by sph-0.1',status)
      call ftpcom(unit, &
        '  May contain both a map and alms and other s2_sky data',status)
      call ftpcom(unit, &
        '  Primary extension empty',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if
      call ftpkyj(unit,'NSIDE', sky%nside, &
        'HEALPix nside',status)
      call ftpkyj(unit,'NPIX',sky%npix, &
        'HEALPix npix',status)
      call ftpkyj(unit,'LMAX',sky%lmax, &
        'max spherical harmonic l considered',status)
      call ftpkyj(unit,'MMAX',sky%mmax, &
        'max spherical harmonic m considered',status)
      call ftpkyj(unit,'NPARAM',sky%n_param, &
        'number of parameters present',status)
      call ftpkyj(unit,'PIXSCHM',sky%pix_scheme, &
        'HEALPix pixelisation scheme of map',status)
      if(sky%pix_scheme == S2_SKY_RING) then
         call ftpkys(unit,'PIXSTR','RING', &
           'String describing HEALPix pixelisation scheme',status)
      elseif(sky%pix_scheme == S2_SKY_NEST) then
         call ftpkys(unit,'PIXSTR','NEST', &
           'String describing HEALPix pixelisation scheme',status)
      end if
      call ftpkyl(unit,'MAPSTAT',sky%map_status, &
        'Set if map is present',status)
      call ftpkyl(unit,'ALMSTAT',sky%alm_status, &
        'Set if alms are present',status)

      ! If map is present then write to binary table in new extension.
      if(sky%map_status) then

         ! Insert binary table extension for map.
         extname='MAP'
         ttype(1)='MAPVALS'
         tform(1)='1E'
         tunit(1)=''
         tfields=1
         nrows=sky%npix
         varidat=0
         call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,&
           varidat,status)

         ! Write map values to binary table.
         frow=1
         felem=1
         colnum=1
         call ftpcle(unit,colnum,frow,felem,nrows,&
           sky%map(0:sky%npix-1), status)

      end if

      ! If alms are present then write to 2D images in new extensions.
      ! Write real and complex parts into separate images.
      if(sky%alm_status) then

         ! Define parameters for alms images.
         naxis=2
         naxes_alm(1)=sky%lmax + 1
         naxes_alm(2)=sky%mmax + 1

         ! Insert a new image extension for real part of alms.
        call ftiimg(unit,bitpix,naxis,naxes_alm,status)

        ! Write additional header keywords.
        call ftpkyj(unit,'LMAX',sky%lmax, &
          'max spherical harmonic l considered',status)
        call ftpkyj(unit,'MMAX',sky%mmax, &
          'max spherical harmonic m considered',status)
        call ftpkys(unit,'EXTNAME','ALM_RL','entension name',status)
      
        ! Write real part of alms as a 2D image.
        group = 1
        dim1 = sky%lmax + 1
        call ftp2de(unit, group, dim1, sky%lmax+1, sky%mmax+1, &
          real(sky%alm(:,:),s2_sp), status)

        ! Insert a new image extension for imaginary part of alms.
        call ftiimg(unit,bitpix,naxis,naxes_alm,status)

        ! Write additional header keywords.
        call ftpkyj(unit,'LMAX',sky%lmax, &
          'max spherical harmonic l considered',status)
        call ftpkyj(unit,'MMAX',sky%mmax, &
          'max spherical harmonic m considered',status)
        call ftpkys(unit,'EXTNAME','ALM_IM','entension name',status)
      
        ! Write imaginary part of alms as a 2D image.
        group = 1
        dim1 = sky%lmax + 1
        call ftp2de(unit, group, dim1, sky%lmax+1, sky%mmax+1, &
          aimag(sky%alm(:,:)), status)

      end if

      ! If present write parameters to binary table in new extension.
      if(sky%n_param > 0) then

         ! Insert binary table extension for map.
         extname='PARAM'
         ttype(1)='PARAMVALS'
         tform(1)='1E'
         tunit(1)=''
         tfields=1
         nrows=sky%n_param
         varidat=0
         call ftibin(unit,nrows,tfields,ttype,tform,tunit,&
           extname,varidat,status)

         ! Write map values to binary table.
         frow=1
         felem=1
         colnum=1
         call ftpcle(unit,colnum,frow,felem,nrows,&
           sky%param(1:sky%n_param), status)

      end if

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2_sky_io_fits_error_check(status, .true.)

    end subroutine s2_sky_io_fits_write

  
    !--------------------------------------------------------------------------
    ! s2_sky_io_fits_read
    !
    !! Reads a full s2_sky fits file and allocates a new sky structure with
    !! the data read.
    !!
    !! Variables:
    !!   - filename: Name of s2_sky fits file containing the sky data to 
    !!     be read.
    !!   - sky: Returned sky structure initialised with the data contained in
    !!     the input s2_sky fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_io_fits_read(filename, sky)

      character(len=*), intent(in) :: filename
      type(s2_sky), intent(out) :: sky

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: ihdu, hdutype, naxis, n_check
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: colnum, frow, felem, nelem
      integer :: group, dim1
      real(s2_sp) :: nullval

      integer :: nside, npix, lmax, mmax, n_param, pix_scheme, fail
      logical :: map_status, alm_status
      real(s2_sp), allocatable :: alm_rl(:,:), alm_im(:,:)

      ! Check object not already initialised.
      if(sky%init) then
        call s2_error(S2_ERROR_INIT, 's2_sky_io_fits_read')
        return
      end if

      ! Initialse error status to zero.
      status=0   

      ! Check if file already exists.
      call s2_sky_io_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call s2_error(S2_ERROR_SKY_FILE_INVALID, &
           's2_sky_io_fits_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)


      ! --------------------------------------
      ! Read primary header
      ! --------------------------------------

      call ftgkyj(unit, 'NSIDE', nside, comment, status)
      call ftgkyj(unit, 'NPIX', npix, comment, status)
      call ftgkyj(unit, 'LMAX', lmax, comment, status)
      call ftgkyj(unit, 'MMAX', mmax, comment, status)
      call ftgkyj(unit, 'NPARAM', n_param, comment, status)
      call ftgkyj(unit, 'PIXSCHM', pix_scheme, comment, status)
      call ftgkyl(unit, 'MAPSTAT', map_status, comment, status)
      call ftgkyl(unit, 'ALMSTAT', alm_status, comment, status)

      ! Initialise sky data structure.
      ! Map, alm and param saved as read in.
      sky = s2_sky_init_empty(nside, pix_scheme, lmax, mmax)
      sky%map_status = map_status
      sky%alm_status = alm_status
      sky%n_param = n_param
      sky%init = .true.      ! Actually already set by s2_sky_init_empty.

      ! Check correct number of HDUs in input file.
      hdunum = 1 ! Due to primary header.
      if(map_status) hdunum = hdunum + 1
      if(alm_status) hdunum = hdunum + 2  ! One for real and one for imaginary.
      if(n_param > 0) hdunum = hdunum + 1
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
         call s2_error(S2_ERROR_SKY_FILE_INVALID, &
           's2_sky_io_fits_read', &
           comment_add='Invalid number of headers')
      end if

      ! Set hdu index ready to read from next extension.
      ihdu = 2 
      

      ! --------------------------------------
      ! Read map if present
      ! --------------------------------------

      if(map_status) then

         ! Allocate space for map.
         allocate(sky%map(0:npix-1), stat=fail)
         if(fail /= 0) then 
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
                 's2_sky_io_fits_read')
         end if
         ! Initialise with zeros.
         sky%map = 0.0e0 

         ! Move to next ihdu extension (i.e. map extension).
         call ftmahd(unit, ihdu, hdutype, status)

         ! Check correct hdutype (i.e. binary table).
         if(hdutype /= 2) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Map not stored in binary table')
         end if

         ! Read header NAXIS2 and check same as npix.
         call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
         if(naxis/=npix) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent number of map pixels')
         end if

         ! Read map values from binary table.
         frow=1
         felem=1
         nelem=npix
         nullval = -999  ! Arbitrary since will stop and return error 
                         ! if null values detected.
         colnum=1      
         call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           sky%map(0:npix-1),anynull,status)
         if(anynull) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Null map values contained in file')
         end if

         ! Increment hdu index ready for next block.
         ihdu = ihdu + 1

      end if
      
      
      ! --------------------------------------
      ! Read alms if present
      ! --------------------------------------
     
      if(alm_status) then

         ! Allocate space for alms.
         allocate(sky%alm(0:lmax,0:mmax), stat=fail)
         allocate(alm_rl(0:lmax,0:mmax), stat=fail)
         allocate(alm_im(0:lmax,0:mmax), stat=fail)
         if(fail /= 0) then 
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
                 's2_sky_io_fits_read')
         end if
         ! Initialise with zeros.
         sky%alm = cmplx(0.0e0, 0.0e0)
         alm_rl = 0.0e0
         alm_im = 0.0e0

         ! Move to next ihdu extension (i.e. alm_rl extension).
         call ftmahd(unit, ihdu, hdutype, status)

         ! Check correct hdutype (i.e. image).
         if(hdutype /= 0) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Real part of alms not stored in image')
         end if

         ! Read header and check correct sizes.
         call ftgkyj(unit, 'NAXIS1', n_check, comment, status)
         if(n_check/=lmax+1) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent lmax size alm_rl extension')
         end if
         call ftgkyj(unit, 'NAXIS2', n_check, comment, status)
         if(n_check/=mmax+1) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent mmax size alm_rl extension')
         end if

         ! Read coefficients as 2D data cube.
         group = 1
         nullval = -999
         dim1 = lmax+1
         call ftg2de(unit, group, nullval, dim1, lmax+1, mmax+1, &
           alm_rl(0:lmax,0:mmax), anynull, status)
         if(anynull) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Null alm_rl values contained in file')
         end if

         ! Move to next ihdu extension (i.e. alm_im extension).
         ihdu = ihdu + 1
         call ftmahd(unit, ihdu, hdutype, status)

         ! Check correct hdutype (i.e. image).
         if(hdutype /= 0) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Real part of alms not stored in image')
         end if

         ! Read header and check correct sizes.
         call ftgkyj(unit, 'NAXIS1', n_check, comment, status)
         if(n_check/=lmax+1) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent lmax size alm_im extension')
         end if
         call ftgkyj(unit, 'NAXIS2', n_check, comment, status)
         if(n_check/=mmax+1) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent mmax size alm_im extension')
         end if

         ! Read coefficients as 2D data cube.
         group = 1
         nullval = -999
         dim1 = lmax+1
         call ftg2de(unit, group, nullval, dim1, lmax+1, mmax+1, &
           alm_im(0:lmax,0:mmax), anynull, status)
         if(anynull) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Null alm_im values contained in file')
         end if

         ! Create complex alm data type form real and imaginary 
         ! components.
         sky%alm = cmplx(alm_rl, alm_im)
         deallocate(alm_rl, alm_im)

         ! Increment hdu index ready for next block.
         ihdu = ihdu + 1

      end if


      ! --------------------------------------
      ! Read param if present
      ! --------------------------------------
     
      if(n_param > 0) then

         ! Allocate space for map.
         allocate(sky%param(1:n_param), stat=fail)
         if(fail /= 0) then 
            call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
                 's2_sky_io_fits_read')
         end if
         ! Initialise with zeros.
        sky% param = 0.0e0 

         ! Move to next ihdu extension (i.e. param extension).
         call ftmahd(unit, ihdu, hdutype, status)

         ! Check correct hdutype (i.e. binary table).
         if(hdutype /= 2) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Param data not stored in binary table')
         end if

         ! Read header NAXIS2 and check same as n_param.
         call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
         if(naxis/=n_param) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Inconsistent size for  param array')
         end if

         ! Read param values from binary table.
         frow=1
         felem=1
         nelem=n_param
         nullval = -999  ! Arbitrary since will stop and return error 
                         ! if null values detected.
         colnum=1      
         call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           sky%param(1:n_param),anynull,status)
         if(anynull) then
            call s2_error(S2_ERROR_SKY_FILE_INVALID, &
              's2_sky_io_fits_read', &
              comment_add='Null param values contained in file')
         end if

         ! Last possible extension so don't need to increment hdu index.

      end if


      ! --------------------------------------
      ! Tidy up
      ! --------------------------------------

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2_sky_io_fits_error_check(status, .true.)

    end subroutine s2_sky_io_fits_read


    !--------------------------------------------------------------------------
    ! s2_sky_io_fits_error_check
    !
    !! Check if a fits error has occured and print error message.  Halt
    !! program execution if halt flag is set.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - status: Fits integer status code.
    !!   - halt: Logical to indicate whether to halt program execution if an 
    !!     error is detected.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_io_fits_error_check(status, halt)

      integer, intent(inout) :: status
      logical, intent(in) :: halt

      character(len=30) :: errtext
      character(len=80) :: errmessage

      !  Check if status is OK (no error); if so, simply return.
      if (status .le. 0) return

      ! The FTGERR subroutine returns a descriptive 30-character text 
      ! string that corresponds to the integer error status number.  
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      ! The FTGMSG subroutine retrieves the oldest message from
      ! the stack and shifts any remaining messages on the stack down one
      ! position.  FTGMSG is called repeatedly until a blank message is
      ! returned, which indicates that the stack is empty.  Each error message
      ! may be up to 80 characters in length. 
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write(*,*) trim(errmessage)
          call ftgmsg(errmessage)
      end do

      if(halt) stop

    end subroutine s2_sky_io_fits_error_check


    !--------------------------------------------------------------------------
    ! s2_sky_io_fits_exists
    !
    !! Check if a fits file exists.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to check existence of.
    !!   - status: Fits integer status code.
    !!   - exists: Logical indicating whether the fits file already exists.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_io_fits_exists(filename, status, exists)

      character(len=*), intent(in) :: filename
      integer, intent(inout) :: status
      logical, intent(out) :: exists

      integer :: unit, blocksize
      logical :: halt

      ! Simply return if status is already greater than zero.
      if (status .gt. 0) return

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      call ftopen(unit, filename, 1, blocksize, status)

      ! Check status of opening file.
      if(status == 0) then

        ! File was opened.  Close it and set exists flag accordingly.
        call ftclos(unit, status)
        exists = .true.

      else if (status == 104) then
        
        ! File does not exist.  Reset status and set exists flag accordingly.
         status = 0
         exists = .false.

      else

        ! Some other error occured while opening file.
        halt = .false.
        call s2_sky_io_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine s2_sky_io_fits_exists


    !--------------------------------------------------------------------------
    ! s2_sky_io_fits_del
    !
    !! Delete a fits file.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to detele.
    !!   - status: Fits integer status code.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_io_fits_del(filename, status)

      character(len=*), intent(in) :: filename
      integer, intent(inout) ::  status

      integer :: unit, blocksize

      ! Simply return if status is greater than zero.
      if (status .gt. 0)return

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Try to open the file, to see if it exists.
      call ftopen(unit,filename,1,blocksize,status)

      if(status .eq. 0) then
         ! File was opened;  so now delete it.
         call ftdelt(unit,status)
      else if(status .eq. 103) then
         ! File doesn't exist, so just reset status to zero and clear errors.
          status=0
          call ftcmsg
      else
         ! There was some other error opening the file; delete the file anyway.
         status=0
         call ftcmsg
         call ftdelt(unit,status)
      end if

      ! Free the unit number for later reuse.
      call ftfiou(unit, status)

    end subroutine s2_sky_io_fits_del


    !--------------------------------------------------------------------------
    ! Set routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_sky_set_lmax
    !
    !! Set a sky lmax.  If alm is calculated it is removed.
    !!
    !! Variables:
    !!   - sky: Sky to set lmax of.
    !!   - lmax: New lmax to set.
    !!   - [mmax_in]: New mmax to set.  If not specified defaults to new lmax.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_set_lmax(sky, lmax, mmax_in)

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: lmax
      integer, intent(in), optional :: mmax_in
 
      integer :: mmax

      if(present(mmax_in)) then
         mmax = mmax_in
      else
         mmax = lmax
      end if
      
      if(lmax == sky%lmax .and. mmax == sky%mmax) then
         ! If same as current lmax and mmax values then do nothing.
         return
      end if

      ! Otherwise set new lmax and mmax.
      sky%lmax = lmax
      sky%mmax = mmax

      ! Set alm_status to not calculated for this new resolution.
      sky%alm_status = .false.
      if(allocated(sky%alm)) deallocate(sky%alm)

      call s2_sky_valid_sizes(sky)
      
    end subroutine s2_sky_set_lmax


    !--------------------------------------------------------------------------
    ! s2_sky_set_nside
    !
    !! Set a sky nside.  If map is calculated it is removed.
    !!
    !! Variables:
    !!   - sky: Sky to set nside of.
    !!   - nside: New nside to set.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_set_nside(sky, nside)

      use pix_tools, only: nside2npix

      type(s2_sky), intent(inout) :: sky
      integer, intent(in) :: nside

      if(nside == sky%nside) then
         ! If nside if the same as previous value then do nothing.
         return
      end if

      ! Otherwise set new nside.
      sky%nside = nside

      ! Calculate npix for the new nside.
      sky%npix = nside2npix(nside)

      ! Set map_status to not calculated for this new resolution.
      sky%map_status = .false.
      if(allocated(sky%map)) deallocate(sky%map)

      call s2_sky_valid_sizes(sky)

    end subroutine s2_sky_set_nside


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_sky_get_init
    !
    !! Get init variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_init(sky) result(init)
      
      type(s2_sky), intent(in) :: sky
      logical :: init

      init = sky%init

    end function s2_sky_get_init


    !--------------------------------------------------------------------------
    ! s2_sky_get_nside
    !
    !! Get nside variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - nside: Object nside variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_nside(sky) result(nside)
      
      type(s2_sky), intent(in) :: sky
      integer :: nside

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_nside')
      end if

      nside = sky%nside

    end function s2_sky_get_nside


    !--------------------------------------------------------------------------
    ! s2_sky_get_npix
    !
    !! Get nside variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - npix: Object npix variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_npix(sky) result(npix)
      
      type(s2_sky), intent(in) :: sky
      integer :: npix

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_npix')
      end if

      npix = sky%npix

    end function s2_sky_get_npix


    !--------------------------------------------------------------------------
    ! s2_sky_get_lmax
    !
    !! Get lmax variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - lmax: Object lmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_lmax(sky) result(lmax)
      
      type(s2_sky), intent(in) :: sky
      integer :: lmax

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_lmax')
      end if

      lmax = sky%lmax

    end function s2_sky_get_lmax


    !--------------------------------------------------------------------------
    ! s2_sky_get_mmax
    !
    !! Get mmax variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - mmax: Object mmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_mmax(sky) result(mmax)
      
      type(s2_sky), intent(in) :: sky
      integer :: mmax

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_mmax')
      end if

      mmax = sky%mmax

    end function s2_sky_get_mmax


    !--------------------------------------------------------------------------
    ! s2_sky_get_pix_scheme
    !
    !! Get pix_scheme variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - pix_scheme: Object pix_scheme variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_pix_scheme(sky) result(pix_scheme)
      
      type(s2_sky), intent(in) :: sky
      integer :: pix_scheme

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_pix_scheme')
      end if

      pix_scheme = sky%pix_scheme

    end function s2_sky_get_pix_scheme


    !--------------------------------------------------------------------------
    ! s2_sky_get_map
    !
    !! Get map variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - map(:): Object map variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_get_map(sky, map)
      
      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(out) :: map(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_map')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_get_map')
      end if
      
      ! Check sizes consistent.
      if(size(map) /= sky%npix) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_get_map', comment_add='Inconsistent size for maps')
      end if

      map = sky%map

    end subroutine s2_sky_get_map


    !--------------------------------------------------------------------------
    ! s2_sky_get_map_pix
    !
    !! Get map pixel from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - map_pix: Object map pixel variable returned.
    !!   - ipix: The index of the pixel to get (in range [0:sky%npix-1]).
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_map_pix(sky, ipix) result(map_pix)
      
      type(s2_sky), intent(in) :: sky
      integer, intent(in) :: ipix
      real(s2_sp) :: map_pix

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_map_pix')
      end if

      ! Check map defined.
      if(.not. sky%map_status) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_sky_get_map_pix')
      end if

      ! Check index valid.
      if(ipix < 0 .or. ipix >= sky%npix) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_get_map_pix', comment_add='Pixel index out of range')
      end if

      map_pix = sky%map(ipix)

    end function s2_sky_get_map_pix


    !--------------------------------------------------------------------------
    ! s2_sky_get_alm
    !
    !! Get alm variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - alm(:,:): Object alm variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_get_alm(sky, alm)
      
      type(s2_sky), intent(in) :: sky
      complex(s2_spc), intent(out) :: alm(:,:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_alm')
      end if

      ! Check alm defined.
      ! If alm not computed stop since lmax and mmax may not be defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_get_alm')
      end if

      ! Check sizes consistent.
      if(size(alm,1) /= sky%lmax+1 .or. size(alm,2) /= sky%mmax+1) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_get_alm', comment_add='Inconsistent size for alms')
      end if

      alm = sky%alm

    end subroutine s2_sky_get_alm


    !--------------------------------------------------------------------------
    ! s2_sky_get_cl
    !
    !! Compute cl for the given sky.  
    !!
    !! Notes:
    !!   - Only successful if alms for the sky are already computed.  Cannot 
    !      automatically compute alms if not already computed since lmax and 
    !!     mmax may not be defined.
    !!   - Initialses a new pl object which must be freed by calling routine.
    !!   - Map assumed real (thus alms only stored from m>=0) hence cl computed
    !!     by: cl = 1/(2*l+1) * ( |a_{l,0}|^2 + 2 * \sum_{m=1}^l |a_{l,m}|^2 )
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - cl: Cls computed for the sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_cl(sky) result(cl)

      type(s2_sky), intent(in) :: sky
      type(s2_pl) :: cl

      real(s2_sp), allocatable :: cl_vals(:)
      integer :: l, m, fail
      
      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_cl')
      end if

      ! Check alm defined.
      ! If alm not computed stop since lmax and mmax may not be defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_get_cl')
      end if

      ! Allocate space for cl values.
      allocate(cl_vals(0:sky%lmax), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_get_cl')
      end if
      
      ! Compute cl values.

      cl_vals = 0

      do l = 0,sky%lmax

         ! Note m starting from 1 in summation.
         do m = 1,min(l, sky%mmax)
            cl_vals(l) = cl_vals(l) + abs(sky%alm(l,m))**2.0e0
         end do

         cl_vals(l) = (2.0e0 * cl_vals(l) + abs(sky%alm(l,0))**2.0e0) &
              / real(2.0e0*min(l,sky%mmax)+1, s2_sp)

      end do

      ! Initialise a pl object with cl values.
      cl = s2_pl_init(cl_vals)

      ! Free temporary memory used.
      deallocate(cl_vals)

    end function s2_sky_get_cl


    !--------------------------------------------------------------------------
    ! s2_sky_get_cm
    !
    !! Compute cm for the given sky.  
    !!
    !! Notes:
    !!   - Only successful if alms for the sky are already computed.  Cannot 
    !      automatically compute alms if not already computed since lmax and 
    !!     mmax may not be defined.
    !!   - Map assumed real (thus alms only stored from m>=0) hence cm computed
    !!     by: cm = 1/(2*(lmax-m)+2-dm0) 
    !!            * \sum_{l=m}^lmax [ 2*(1-dm0)*|a_{l,m}|^2 + dm0*|a_{l,0}|^2 ]
    !!     where dm0 is the Kronecker delta.
    !!   - Predominately used for determined azimuthal band limit of a sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - cm: Cms computed for the sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_get_cm(sky, cm)

      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(out) :: cm(0:)

      integer :: l,m
      real(s2_sp) :: delta_m0 = 1.0e0

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_cm')
      end if

      ! Check alm defined.
      ! If alm not computed stop since lmax and mmax may not be defined.
      if(.not. sky%alm_status) then
         call s2_error(S2_ERROR_SKY_ALM_NOT_DEF, 's2_sky_get_cm')
      end if

      ! Check cm of correct size.
      if(size(cm) /= sky%mmax+1) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_sky_get_cm', &
           comment_add='Size of cm array invalid')
      end if

      ! Compute cm values.

      cm = 0.0e0

      delta_m0 = 1.0e0

      do m = 0,sky%mmax

         do l = m,sky%lmax
            cm(m) = cm(m) &
               + 2.0e0 * (1.0e0-delta_m0) * abs(sky%alm(l,m))**2.0e0 &
               + delta_m0 * abs(sky%alm(l,0))**2.0e0
         end do

         cm(m) = cm(m) / real(2.0e0*(sky%lmax - m) + 2 - delta_m0, s2_sp)

         delta_m0 = 0.0e0

      end do

    end subroutine s2_sky_get_cm


    !--------------------------------------------------------------------------
    ! s2_sky_get_map_status
    !
    !! Get map_status variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - map_status: Object map_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_map_status(sky) result(map_status)
      
      type(s2_sky), intent(in) :: sky
      logical :: map_status

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_map_status')
      end if

      map_status = sky%map_status

    end function s2_sky_get_map_status


    !--------------------------------------------------------------------------
    ! s2_sky_get_alm_status
    !
    !! Get map_status variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get variable of.
    !!   - alm_status: Object alm_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_alm_status(sky) result(alm_status)
      
      type(s2_sky), intent(in) :: sky
      logical :: alm_status

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_alm_status')
      end if

      alm_status = sky%alm_status

    end function s2_sky_get_alm_status


    !--------------------------------------------------------------------------
    ! s2_sky_get_n_param
    !
    !! Get n_param variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get variable of.
    !!   - n_param: Object n_param variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_sky_get_n_param(sky) result(n_param)
      
      type(s2_sky), intent(in) :: sky
      integer :: n_param

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_n_param')
      end if

      n_param = sky%n_param

    end function s2_sky_get_n_param


    !--------------------------------------------------------------------------
    ! s2_sky_get_param
    !
    !! Get param variable from the passed sky.
    !!
    !! Variables:
    !!   - sky: Sky object to get the variable of.
    !!   - param(:): Object param variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_sky_get_param(sky, param)
      
      type(s2_sky), intent(in) :: sky
      real(s2_sp), intent(out) :: param(:)

      ! Check object initialised.
      if(.not. sky%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_sky_get_param')
      end if

      if(size(param) /= sky%n_param .or. sky%n_param <= 0) then
         call s2_error(S2_ERROR_SKY_SIZE_INVALID, &
           's2_sky_get_param', comment_add='Inconsistent size for param')
      end if

      param = sky%param

    end subroutine s2_sky_get_param


end module s2_sky_mod

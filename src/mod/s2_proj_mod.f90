!------------------------------------------------------------------------------
! s2_sky_mod -- S2 library sky class
!
!! Provides functionality to project a sky object onto a planar image.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module s2_proj_mod

  use s2_types_mod, only: s2_sp, s2_spc, s2_dp, s2_dpc, pi
  use s2_error_mod
  use s2_sky_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_proj_init, &
    s2_proj_free, &
    s2_proj_operator_nearest_neighbour, &
    s2_proj_operator_kernel, &
    s2_proj_write_image_file


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  ! No interfaces.


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! No projection method.
  integer, public, parameter :: S2_PROJ_METHOD_NONE = 0

  !! Project using nearest neighbour method.
  integer, public, parameter :: S2_PROJ_METHOD_NEAREST_NEIGHBOUR = 1

  !! Project using harmonic interpolation method.
  integer, public, parameter :: S2_PROJ_METHOD_HARMONIC_INTERP = 2

  !! Project using convolution kernel interpolation method.
  integer, public, parameter :: S2_PROJ_METHOD_KERNEL = 3

  !! Project upper hemisphere to plane.
  integer, public, parameter :: S2_PROJ_FIELD_HEMISPHERE_UPPER = 0

  !! Project lower hemisphere to plane.
  integer, public, parameter :: S2_PROJ_FIELD_HEMISPHERE_LOWER = 1


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_proj
    private
    logical :: init = .false.
    type(s2_sky) :: parent
    real(s2_sp) :: theta_fov = 0e0
    integer :: N = 0
    real(s2_sp) :: dx = 0e0
    real(s2_sp) :: umax = 0e0
    real(s2_sp) :: image_size = 0e0
    integer :: method = S2_PROJ_METHOD_NONE
    integer :: field = S2_PROJ_FIELD_HEMISPHERE_UPPER
    real(s2_sp), allocatable :: image(:,:)
  end type s2_proj


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_proj_init
    !
    !! Initialie a proj object from a sky.
    !!
    !! Variables:
    !!   - sky: Sky to project onto planar image.
    !!   - theta_fov: Theta field of view.
    !!   - method: Method to use to project sky onto planar image.
    !!   - field: Hemisphere to project onto plane.
    !!   - [nside]: Healpix nside (required by some projection methods if not 
    !!     provided in sky object).
    !!   - [lmax]: Alm lmax (required by some projection methods if not 
    !!     provided in sky object).
    !!   - proj: The returned sky initialised with resolutions and projected 
    !!     planar image.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_proj_init(sky, theta_fov, method, field, nside, lmax) &
         result(proj)

      type(s2_sky), intent(inout) :: sky
        ! inout since may need to compute map or alms
      real(s2_sp), intent(in) :: theta_fov
      integer, intent(in) :: method
      integer, intent(in), optional :: field
      integer, intent(in), optional :: nside
      integer, intent(in), optional :: lmax
      type(s2_proj) :: proj
      
      integer :: N
      integer :: lmax_use
      real(s2_sp) :: dx
      real(s2_sp) :: umax
      real(s2_sp) :: image_size
      integer :: fail = 0

      ! Check object not already initialised.
      if(proj%init) then
        call s2_error(S2_ERROR_INIT, 's2_proj_init')
        return
      end if

      if (present(lmax)) then
         lmax_use = lmax
      else
         lmax_use = s2_sky_get_lmax(sky)
      end if

      ! Compute projection parameters.
      call s2_proj_compute_params(theta_fov, lmax_use, N, dx, umax, image_size) 

      ! Allocate space for image.
      allocate(proj%image(0:N-1, 0:N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_init')
      end if
      proj%image(0:N-1, 0:N-1) = 0e0

      ! Set attributes.
      proj%parent = s2_sky_init(sky)
      proj%theta_fov = theta_fov
      proj%N = N
      proj%dx = dx
      proj%umax = umax
      proj%image_size = image_size
      proj%method = method
      if(present(field)) proj%field = field

      ! Create image.
      select case(method)

         case(S2_PROJ_METHOD_NEAREST_NEIGHBOUR)
            call s2_proj_projection_nearest_neighbour(sky, proj, nside)

         case(S2_PROJ_METHOD_HARMONIC_INTERP)
            call s2_proj_projection_harmonic_interp(sky, proj, lmax_use)

         case(S2_PROJ_METHOD_KERNEL)
            call s2_proj_projection_kernel(sky, proj, nside)

         case default
            call s2_error(S2_ERROR_PROJ_METHOD_INVALID, 's2_proj_init', &
              comment_add='Invalid method type specifier')

      end select

      ! Set as initialised.
      proj%init = .true.

    end function s2_proj_init


    !--------------------------------------------------------------------------
    ! s2_proj_free
    !
    !! Free all data associated with an initialised proj object and reset all
    !! other attributes.
    !!
    !! Variables:
    !!   - proj: The proj object to be freed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_free(proj)

      type(s2_proj), intent(inout) :: proj

      ! Check object initialised.
      if(.not. proj%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_proj_free')
      end if 

      ! Free space.
      if(allocated(proj%image)) deallocate(proj%image)
      
      ! Reset attributes.
      proj%theta_fov = 0e0
      proj%N = 0
      proj%dx = 0e0
      proj%umax = 0e0
      proj%image_size = 0e0
      proj%method = S2_PROJ_METHOD_NONE
      proj%field = S2_PROJ_FIELD_HEMISPHERE_UPPER

      proj%init = .false.

    end subroutine s2_proj_free


    !--------------------------------------------------------------------------
    ! s2_proj_compute_params
    !
    !! Compute size parameters of planar image from field of view and harmonic
    !! band-limit.
    !!
    !! Variables:
    !!   - theta_fov: Input theta field of view.
    !!   - lmax: Input harmonic band-limit.
    !!   - N: Computed number of pixels in each dimension of square image.
    !!   - dx: Computed spacing between pixels in image.
    !!   - umax: Computed band-limit in plane.
    !!   - image_size: Compute size of image.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_compute_params(theta_fov, lmax, N, dx, umax, image_size) 
      
      real(s2_sp), intent(in) :: theta_fov
      integer, intent(in) :: lmax
      integer, intent(out) :: N
      real(s2_sp), intent(out) :: dx
      real(s2_sp), intent(out) :: umax
      real(s2_sp), intent(out) :: image_size

      real(s2_sp) :: image_size_ideal

      if(lmax == 0) then
        call s2_error(S2_ERROR_SKY_SIZE_NOT_DEF, 's2_proj_compute_params', &
             comment_add='lmax is zero')
      end if

      umax = real(lmax,s2_sp) / (2.0 * pi * cos(theta_fov/2.0))
      dx = 1.0 / (2.0 * umax)
      image_size_ideal = 2.0 * sin(theta_fov/2.0)
      N = ceiling(image_size_ideal / dx)
      ! If natural N is not integer, increase N to ensure achieve at least 
      ! required fov.

      image_size = N * dx

    end subroutine s2_proj_compute_params


    !--------------------------------------------------------------------------
    ! s2_proj_projection_nearest_neighbour
    !
    !! Project sky to planar image using nearest neighbour method.
    !!
    !! Notes:
    !!   - Space to store image must be allocated before calling this routine.
    !!
    !! Variables:
    !!   - sky: Sky to project to planar image.
    !!   - proj: Projected sky.
    !!   - [nside]: Optional Healpix nside to consider in map of sky when 
    !!     projecting image (if not provided taken from sky).
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_projection_nearest_neighbour(sky, proj, nside)

      use s2_vect_mod
      use pix_tools, only: ang2pix_ring, ang2pix_nest

      type(s2_sky), intent(inout) :: sky 
        ! Inout since may need to compute map.
      type(s2_proj), intent(inout) :: proj
      integer, intent(in), optional :: nside

      integer :: N, nside_use, pix_scheme, ipix, i, j, k
      integer :: fail = 0
      real(s2_sp), allocatable :: grid(:)
      real(s2_sp) :: x(1:3)
      real(s2_dp) :: theta, phi
      type(s2_vect) :: vec

      ! Check object not already initialised.
      if(proj%init) then
        call s2_error(S2_ERROR_INIT, 's2_proj_projection_nearest_neighbour')
        return
      end if

      ! Only projection of upper hemisphere supported at present.
      if(proj%field /= S2_PROJ_FIELD_HEMISPHERE_UPPER) then
         call s2_error(S2_ERROR_PROJ_FIELD_INVALID, 's2_proj_projection_nearest_neighbour', &
              comment_add='Only upper hemisphere supported at present.')
      end if

      ! Get local parameters.
      N = proj%N
      if(present(nside)) then
         nside_use = nside
      else
         nside_use = s2_sky_get_nside(sky)
      end if
      pix_scheme = s2_sky_get_pix_scheme(sky)

      ! Define planar grid.
      allocate(grid(0:N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_projection_nearest_neighbour')
      end if
      grid = (/  ((k+0.5)*proj%dx - proj%image_size/2.0, k = 0,N-1) /) 

      ! Ensure sky map is defined and compute otherwise.
      call s2_sky_compute_map(sky, nside_use)

      ! Project on planar grid.
      do i = 0,N-1
         x(1) = grid(i)
         do j =0,N-1             
            x(2) = grid(j)
            if(x(1)**2 + x(2)**2 > 1.0) then
               x(3) = 0 
                 ! Account for points outside disk with unit radius 
                 ! (will ensure map to theta=pi/2, which is outside fov).
            else
               x(3) = sqrt(1.0 - x(1)**2 - x(2)**2)
            end if

            ! Compute spherical corrdinates of planar grid point.
            vec = s2_vect_init(x)
            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
            theta = s2_vect_get_theta(vec)
            phi = s2_vect_get_phi(vec)
            call s2_vect_free(vec)

           if(theta <= proj%theta_fov/2.0) then

              ! Find pixel index corresponding to position on sphere.
              if(pix_scheme == S2_SKY_RING) then
                 call ang2pix_ring(nside_use, theta, phi, ipix)
              else if(pix_scheme == S2_SKY_NEST) then
                 call ang2pix_nest(nside_use, theta, phi, ipix)
              else
                 call s2_error(S2_ERROR_SKY_PIX_INVALID, &
                      's2_proj_projection_nearest_neighbour')
              end if

              ! Set image value.
              proj%image(i,j) = s2_sky_get_map_pix(sky, ipix)

           else
              proj%image(i,j) = 0.0
           end if

         end do
      end do

      ! Free memory.
      deallocate(grid)
      
    end subroutine s2_proj_projection_nearest_neighbour


    !--------------------------------------------------------------------------
    ! s2_proj_operator_nearest_neighbour
    !
    !! Compute sparse representaion of linear operation to perform nearest 
    !! neighbour projection.  The map vector containing the pixels on the 
    !! sphere within the fov is also returned.
    !!
    !! Notes:
    !!   - Space for op and xmap allocate herein; must be freed by calling 
    !!     routine.
    !!
    !! Variables:
    !!   - proj: Projected sky.
    !!   - [nside]: Optional Healpix nside to consider in map of sky when 
    !!     projecting image (if not provided taken from sky).
    !!   - nop: Number of non-zero enties in the projection operator.
    !!   - op(0:nop-1,0:1): Projection operator (specifies indices of
    !!     unit entries in operator matrix).
    !!   - nsphere: Number of pixels on the sphere within the fov.
    !!   - xmap(0:nsphere-1): Vector of pixels values on the sphere within 
    !!     the fov.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_operator_nearest_neighbour(proj, nside, nop, op, nsphere, xmap)

      use s2_vect_mod
      use pix_tools, only: ang2pix_ring, ang2pix_nest, pix2ang_ring

        ! Inout since may need to compute map.
      type(s2_proj), intent(inout) :: proj
      integer, intent(in), optional :: nside
      integer, intent(out) :: nop
      integer, allocatable, intent(out) :: op(:,:)
      integer, intent(out) :: nsphere
      real(s2_sp), allocatable, intent(out) :: xmap(:)

      integer :: nside_use, pix_scheme, npix
      integer :: fail = 0
      real(s2_dp) :: theta, phi
      integer :: ipix, i, j, k, iang, iop
      real(s2_sp) :: x(1:3)
      type(s2_vect) :: vec
      real(s2_sp), allocatable :: map(:)
      integer, allocatable :: opx(:,:)
      real(s2_sp), allocatable :: grid(:)

      ! Check object not already initialised.
      if(.not. proj%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_proj_operator_nearest_neighbour')
        return
      end if

      ! Only projection of upper hemisphere supported at present.
      if(proj%field /= S2_PROJ_FIELD_HEMISPHERE_UPPER) then
         call s2_error(S2_ERROR_PROJ_FIELD_INVALID, 's2_proj_operator_nearest_neighbour', &
              comment_add='Only upper hemisphere supported at present.')
      end if

      ! Get local parameters.
      if(present(nside)) then
         nside_use = nside
      else
         nside_use = s2_sky_get_nside(proj%parent)
      end if

      ! Ensure sky map is defined and compute otherwise.
      call s2_sky_compute_map(proj%parent, nside_use)

      ! Ensure map in RING pixelisation scheme.
      call s2_sky_map_convert(proj%parent, S2_SKY_RING)
      pix_scheme = s2_sky_get_pix_scheme(proj%parent)
      npix = s2_sky_get_npix(proj%parent)

      ! Define planar grid.
      allocate(grid(0:proj%N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      grid = (/  ((k+0.5)*proj%dx - proj%image_size/2.0, k = 0,proj%N-1) /) 

      ! Get xmap vector.
      allocate(map(0:npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      call s2_sky_get_map(proj%parent, map)

      nsphere = 0
      do ipix = 0,npix-1 
         call pix2ang_ring(nside_use, ipix, theta, phi)
         if(theta > proj%theta_fov/2.0) then
            nsphere = ipix
            exit            
         end if
      end do

      allocate(xmap(0:nsphere-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      xmap(0:nsphere-1) = map(0:nsphere-1)

      ! Allocate maximum required space for sparse representation of operator.
      allocate(opx(0:proj%N*proj%N-1, 0:1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if      

      ! Construct operator.
      iang = 0
      iop = 0
      do i = 0,proj%N-1
         x(1) = grid(i)
         do j =0,proj%N-1             
            x(2) = grid(j)
            if(x(1)**2 + x(2)**2 > 1.0) then
               x(3) = 0 
                 ! Account for points outside disk with unit radius 
                 ! (will ensure map to theta=pi/2, which is outside fov).
            else
               x(3) = sqrt(1.0 - x(1)**2 - x(2)**2)
            end if

            ! Compute spherical corrdinates of planar grid point.
            vec = s2_vect_init(x)
            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
            theta = s2_vect_get_theta(vec)
            phi = s2_vect_get_phi(vec)
            call s2_vect_free(vec)

           if(theta <= proj%theta_fov/2.0) then

              ! Find pixel index corresponding to position on sphere.
              if(pix_scheme == S2_SKY_RING) then
                 call ang2pix_ring(nside_use, theta, phi, ipix)
              else if(pix_scheme == S2_SKY_NEST) then
                 call ang2pix_nest(nside_use, theta, phi, ipix)
              else
                 call s2_error(S2_ERROR_SKY_PIX_INVALID, &
                      's2_proj_operator_nearest_neighbour')
              end if
              if(ipix < nsphere) then
                 ! Ensure only include pixels within fov.
                 ! (Otherwise could find theta within fov, but nearest pixel outside fov).
                 opx(iop, 0) = iang
                 opx(iop, 1) = ipix
                 iop = iop + 1
              end if
           end if

           iang = iang + 1
         end do
      end do
      nop = iop

      ! Copy required operator points.
      allocate(op(0:nop-1, 0:1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if   
      op(0:nop-1, 0:1) = opx(0:nop-1, 0:1)

      ! Free memory.
      deallocate(grid)
      deallocate(map)
      deallocate(opx)

    end subroutine s2_proj_operator_nearest_neighbour


    !--------------------------------------------------------------------------
    ! s2_proj_projection_kernel
    !
    !! Project sky to planar image using convolution kernel projection method.
    !!
    !! Notes:
    !!   - Space to store image must be allocated before calling this routine.
    !!
    !! Variables:
    !!   - sky: Sky to project to planar image.
    !!   - proj: Projected sky.
    !!   - [nside]: Optional Healpix nside to consider in map of sky when 
    !!     projecting image (if not provided taken from sky).
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_projection_kernel(sky, proj, nside)

      use s2_vect_mod
      use pix_tools, only: ang2pix_ring, ang2pix_nest
      
      interface 
         function kernel(theta, param) result(val)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta
           real(s2_dp), intent(in), optional :: param(:)
           real(s2_dp) :: val
         end function kernel
      end interface

      type(s2_sky), intent(inout) :: sky 
        ! Inout since may need to compute map.
      type(s2_proj), intent(inout) :: proj
      integer, intent(in), optional :: nside

      integer :: N, nside_use, pix_scheme, i, j, k
      integer :: fail = 0
      real(s2_sp), allocatable :: grid(:)
      real(s2_sp) :: x(1:3)
      real(s2_dp) :: theta, phi
      type(s2_vect) :: vec
      real(s2_dp) :: sigma(1), support_theta

      ! Check object not already initialised.
      if(proj%init) then
        call s2_error(S2_ERROR_INIT, 's2_proj_projection_nearest_neighbour')
        return
      end if

      ! Only projection of upper hemisphere supported at present.
      if(proj%field /= S2_PROJ_FIELD_HEMISPHERE_UPPER) then
         call s2_error(S2_ERROR_PROJ_FIELD_INVALID, 's2_proj_projection_nearest_neighbour', &
              comment_add='Only upper hemisphere supported at present.')
      end if

      ! Get local parameters.
      N = proj%N
      if(present(nside)) then
         nside_use = nside
      else
         nside_use = s2_sky_get_nside(sky)
      end if
      pix_scheme = s2_sky_get_pix_scheme(sky)

      ! Define planar grid.
      allocate(grid(0:N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_projection_nearest_neighbour')
      end if
      grid = (/  ((k+0.5)*proj%dx - proj%image_size/2.0, k = 0,N-1) /) 

      ! Ensure sky map is defined and compute otherwise.
      call s2_sky_compute_map(sky, nside_use)

      ! Project on planar grid.
      do i = 0,N-1
         x(1) = grid(i)
         do j =0,N-1             
            x(2) = grid(j)
            if(x(1)**2 + x(2)**2 > 1.0) then
               x(3) = 0 
                 ! Account for points outside disk with unit radius 
                 ! (will ensure map to theta=pi/2, which is outside fov).
            else
               x(3) = sqrt(1.0 - x(1)**2 - x(2)**2)
            end if

            ! Compute spherical corrdinates of planar grid point.
            vec = s2_vect_init(x)
            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
            theta = s2_vect_get_theta(vec)
            phi = s2_vect_get_phi(vec)
            call s2_vect_free(vec)

           if(theta <= proj%theta_fov/2.0) then

              sigma(1) = 0.02
              support_theta = 4 * sigma(1)
              proj%image(i,j) = s2_sky_convpt_space(sky, support_theta, kernel, &
                   theta, phi, sigma)

           else
              proj%image(i,j) = 0.0
           end if

         end do
      end do

      ! Free memory.
      deallocate(grid)
      
    end subroutine s2_proj_projection_kernel


    !--------------------------------------------------------------------------
    ! s2_proj_operator_kernel
    !
    !! Compute sparse representaion of linear operation to perform convolution  
    !! kernel projection.  The map vector containing the pixels on the 
    !! sphere within the fov is also returned.
    !!
    !! Notes:
    !!   - Space for op and xmap allocate herein; must be freed by calling 
    !!     routine.
    !!
    !! Variables:
    !!   - proj: Projected sky.
    !!   - [nside]: Optional Healpix nside to consider in map of sky when 
    !!     projecting image (if not provided taken from sky).
    !!   - nop: Number of non-zero enties in the projection operator.
    !!   - op(0:nop-1,0:2): Projection operator (specifies indices and values 
    !!     of entries in operator matrix).
    !!   - nsphere: Number of pixels on the sphere within the fov.
    !!   - xmap(0:nsphere-1): Vector of pixels values on the sphere within 
    !!     the fov.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   July 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_operator_kernel(proj, nside, nop, op, nsphere, xmap)

      use s2_vect_mod
      use pix_tools, only: ang2pix_ring, ang2pix_nest, pix2ang_ring

      interface 
         function kernel(theta, param) result(val)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta
           real(s2_dp), intent(in), optional :: param(:)
           real(s2_dp) :: val
         end function kernel
      end interface

        ! Inout since may need to compute map.
      type(s2_proj), intent(inout) :: proj
      integer, intent(in), optional :: nside
      integer, intent(out) :: nop
      real(s2_dp), allocatable, intent(out) :: op(:,:)
      integer, intent(out) :: nsphere
      real(s2_sp), allocatable, intent(out) :: xmap(:)

      integer :: nside_use, pix_scheme, npix
      integer :: fail = 0
      real(s2_dp) :: theta, phi
      integer :: ipix, i, j, k, iang, iop
      real(s2_sp) :: x(1:3)
      type(s2_vect) :: vec


      real(s2_dp) :: sigma(1), support_theta
      integer :: nweights, iweight
      integer, allocatable :: indices(:)
      real(s2_dp), allocatable :: weights(:)

      real(s2_sp), allocatable :: map(:)
      real(s2_dp), allocatable :: opx(:,:)
      real(s2_sp), allocatable :: grid(:)

      ! Check object not already initialised.
      if(.not. proj%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_proj_operator_nearest_neighbour')
        return
      end if

      ! Only projection of upper hemisphere supported at present.
      if(proj%field /= S2_PROJ_FIELD_HEMISPHERE_UPPER) then
         call s2_error(S2_ERROR_PROJ_FIELD_INVALID, 's2_proj_operator_nearest_neighbour', &
              comment_add='Only upper hemisphere supported at present.')
      end if

      ! Get local parameters.
      if(present(nside)) then
         nside_use = nside
      else
         nside_use = s2_sky_get_nside(proj%parent)
      end if

      ! Ensure sky map is defined and compute otherwise.
      call s2_sky_compute_map(proj%parent, nside_use)

      ! Ensure map in RING pixelisation scheme.
      call s2_sky_map_convert(proj%parent, S2_SKY_RING)
      pix_scheme = s2_sky_get_pix_scheme(proj%parent)
      npix = s2_sky_get_npix(proj%parent)

      ! Define planar grid.
      allocate(grid(0:proj%N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      grid = (/  ((k+0.5)*proj%dx - proj%image_size/2.0, k = 0,proj%N-1) /) 

      ! Get xmap vector.
      allocate(map(0:npix-1), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      call s2_sky_get_map(proj%parent, map)

      nsphere = 0
      do ipix = 0,npix-1 
         call pix2ang_ring(nside_use, ipix, theta, phi)
         if(theta > proj%theta_fov/2.0) then
            nsphere = ipix
            exit            
         end if
      end do

      allocate(xmap(0:nsphere-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if
      xmap(0:nsphere-1) = map(0:nsphere-1)

      ! Allocate maximum required space for sparse representation of operator.
      allocate(opx(0:proj%N*proj%N*proj%N*proj%N-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if      

      ! Construct operator.
      iang = 0
      iop = 0
      do i = 0,proj%N-1
         x(1) = grid(i)
         do j =0,proj%N-1             
            x(2) = grid(j)
            if(x(1)**2 + x(2)**2 > 1.0) then
               x(3) = 0 
                 ! Account for points outside disk with unit radius 
                 ! (will ensure map to theta=pi/2, which is outside fov).
            else
               x(3) = sqrt(1.0 - x(1)**2 - x(2)**2)
            end if

            ! Compute spherical corrdinates of planar grid point.
            vec = s2_vect_init(x)
            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
            theta = s2_vect_get_theta(vec)
            phi = s2_vect_get_phi(vec)
            call s2_vect_free(vec)

           if(theta <= proj%theta_fov/2.0) then

              sigma(1) = 0.02
              support_theta = 4 * sigma(1)
              call s2_sky_convpt_space_weights(nweights, indices, weights, &
                   nside_use, S2_SKY_RING, support_theta, kernel, &
                   theta, phi, sigma)

              do iweight = 0,nweights-1
                 if(indices(iweight) < nsphere) then   
                    ! Ensure only include pixels within fov.             
                    opx(iop, 0) = iang
                    opx(iop, 1) = indices(iweight)
                    opx(iop, 2) = weights(iweight)
                    iop = iop + 1
                 end if
              end do

              deallocate(indices, weights)

           end if

           iang = iang + 1
         end do
      end do
      nop = iop

      ! Copy required operator points.
      allocate(op(0:nop-1, 0:2), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_proj_operator_nearest_neighbour')
      end if   
      op(0:nop-1, 0:2) = opx(0:nop-1, 0:2)

      ! Free memory.
      deallocate(grid)
      deallocate(map)
      deallocate(opx)

    end subroutine s2_proj_operator_kernel


    !--------------------------------------------------------------------------
    ! s2_proj_projection_harmonic_interp
    !
    !! Project sky to planar image using harmonic interpolation method.
    !!
    !! Notes:
    !!   - Space to store image must be allocated before calling this routine.
    !!   - The associate Legendre functions computed using Numerical Recipes
    !!     routines appear to be stable to lmax~128 (otherwise return NaN), 
    !!     while using Wigner dlmn matrices is relatively slow, but can run 
    !!     at very high band-limits.  The current impleentation is based on 
    !!     Wigner dlmn matrices.
    !!
    !! Variables:
    !!   - sky: Sky to project to planar image.
    !!   - proj: Projected sky.
    !!   - [lmax]: Optional harmonic band-limit lmax to consider in harmonic
    !!     representation of sky when projecting image (if not provided taken 
    !!     from sky).
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_projection_harmonic_interp(sky, proj, lmax)

      use s2_vect_mod
      use s2_ylm_mod

      type(s2_sky), intent(inout) :: sky 
        ! Inout since may need to compute alms.
      type(s2_proj), intent(inout) :: proj
      integer, intent(in), optional :: lmax

      integer :: N, lmax_use, i, j, k, el, m
      integer :: fail = 0
      real(s2_sp), allocatable :: grid(:)
!!$      complex(s2_spc), allocatable :: alm(:,:)
      real(s2_sp) :: x(1:3)
      real(s2_dp) :: theta, phi
      type(s2_vect) :: vec


      integer :: iang
      real(s2_dp), allocatable :: thetas(:), phis(:), f(:)

      ! Check object not already initialised.
      if(proj%init) then
        call s2_error(S2_ERROR_INIT, 's2_proj_projection_harmonic_interp')
        return
      end if

      ! Only projection of upper hemisphere supported at present.
      if(proj%field /= S2_PROJ_FIELD_HEMISPHERE_UPPER) then
         call s2_error(S2_ERROR_PROJ_FIELD_INVALID, 's2_proj_projection_nearest_neighbour', &
              comment_add='Only upper hemisphere supported at present.')
      end if

      ! Get local parameters.
      N = proj%N
      if(present(lmax)) then
         lmax_use = lmax
      else
         lmax_use = s2_sky_get_lmax(sky)
      end if

      ! Define planar grid.
      allocate(grid(0:N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
             's2_proj_projection_harmonic_interp')
      end if
      grid = (/  ((k+0.5)*proj%dx - proj%image_size/2.0, k = 0,N-1) /) 

      ! Ensure sky alms are defined and compute otherwise.
      call s2_sky_compute_alm(sky, lmax_use, lmax_use)

!!$      ! Get alms.
!!$      allocate(alm(0:lmax_use, 0:lmax_use), stat=fail)
!!$      if(fail /= 0) then
!!$        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
!!$             's2_proj_projection_harmonic_interp')
!!$      end if
!!$      call s2_sky_get_alm(sky, alm)

      ! Compute (theta,phi) grid points.
      allocate(thetas(0:N*N-1), stat=fail)
      allocate(phis(0:N*N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
             's2_proj_projection_harmonic_interp')
      end if
      iang = 0
      do i = 0,N-1
         x(1) = grid(i)
         do j = 0,N-1             
            x(2) = grid(j)
            if(x(1)**2 + x(2)**2 > 1.0) then
               x(3) = 0 
                 ! Account for points outside disk with unit radius 
                 ! (will ensure map to theta=pi/2, which is outside fov).
            else
               x(3) = sqrt(1.0 - x(1)**2 - x(2)**2)
            end if

            ! Compute spherical corrdinates of planar grid point.
            vec = s2_vect_init(x)
            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
            thetas(iang) = s2_vect_get_theta(vec)
            phis(iang) = s2_vect_get_phi(vec)
            iang = iang + 1
            call s2_vect_free(vec)

         end do
      end do

      ! Compute irregular inverse spherical harmonic transform.
      allocate(f(0:N*N-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
             's2_proj_projection_harmonic_interp')
      end if
      call s2_sky_irregular_invsht(sky, f, N*N, thetas, phis)

      ! Extact image.
      iang = 0
      do i = 0,N-1
         do j = 0,N-1 
            if (thetas(iang) <= proj%theta_fov/2.0) then
               proj%image(i,j) = f(iang)
            else
               proj%image(i,j) = 0.0
            end if
            iang = iang + 1
         end do
      end do

!!$      ! Project on planar grid.
!!$      do i = 0,N-1
!!$         x(1) = grid(i)
!!$         do j = 0,N-1             
!!$            x(2) = grid(j)
!!$            x(3) = sqrt(1.0 - x(1)**2 + x(2)**2)
!!$
!!$            ! Compute spherical corrdinates of planar grid point.
!!$            vec = s2_vect_init(x)
!!$            call s2_vect_convert(vec, S2_VECT_TYPE_S2)
!!$            theta = s2_vect_get_theta(vec)
!!$            phi = s2_vect_get_phi(vec)
!!$            call s2_vect_free(vec)
!!$
!!$            ! Compute projected image
!!$            proj%image(i,j) = 0e0
!!$            do el = 0,lmax_use
!!$!write(*,*) 'el=', el
!!$!write(*,*) proj%image(i,j)
!!$               proj%image(i,j) = proj%image(i,j) + &
!!$                    real(alm(el,0)) * s2_ylm_eval_leg(el, 0, theta, phi)
!!$!write(*,*) proj%image(i,j)
!!$
!!$               do m = 1,min(el,lmax_use)
!!$!write(*,*) 'm=',m
!!$                  proj%image(i,j) = proj%image(i,j) + &
!!$                       2.0 * real(alm(el,m) * s2_ylm_eval_leg(el, m, theta, phi))
!!$!write(*,*) proj%image(i,j)
!!$               end do
!!$!write(*,*)
!!$!write(*,*) 
!!$            end do
!!$
!!$         end do
!!$      end do

      ! Free memory.
      deallocate(grid)
!!$      deallocate(alm)
      deallocate(f)
      deallocate(thetas, phis)
      
    end subroutine s2_proj_projection_harmonic_interp


    !--------------------------------------------------------------------------
    ! s2_proj_write_image_file
    !
    !! Write a projected image to a file. 
    !!
    !! Variables:
    !!   - proj: Projected sky.
    !!   - filename: Name of the output file.
    !!   - [comment]: Optional additional comment to be added to the file
    !!     header.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control.
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_proj_write_image_file(proj, filename, comment)

      use pix_tools, only: pix2ang_ring, pix2ang_nest

      type(s2_proj), intent(in) :: proj
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: fileid, i, j
      integer :: fail = 0

      ! Check object initialised.
      if(.not. proj%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_proj_write_image_file')
      end if

      ! Open file.
      fileid = 22
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')

      ! Write to file.
      if(present(comment)) write(fileid,'(a,a)') '% ', trim(comment)
      do i = 0,proj%N-1
         do j = 0,proj%N-1                       
            write(fileid,'(e28.20)') proj%image(i,j)
         end do
      end do

      ! Close file.
      close(fileid)

    end subroutine s2_proj_write_image_file


end module s2_proj_mod



function kernel(theta, param) result(val)

  use s2_types_mod
  real(s2_dp), intent(in) :: theta
  real(s2_dp), intent(in), optional :: param(:)
  real(s2_dp) :: val

  real(s2_dp) :: sigma

  sigma = param(1)
  val = exp(-theta**2 / (2.0*sigma**2))

end function kernel

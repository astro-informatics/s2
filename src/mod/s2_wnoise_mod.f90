!------------------------------------------------------------------------------
! s2_wnoise_mod -- S2 library wnoise class
!
!! Provides functionality to support white noise realisations on the sky.
!! Noise may either be constructed from a uniform standard devaition constant 
!! over the sky or from a standard deviation map that varies over the sky 
!! depending on the number of observations at a particular position on the sky.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_wnoise_mod

  use s2_types_mod
  use s2_error_mod
  use s2_distn_mod
  use s2_sky_mod
  use s2_pl_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_wnoise_init, &
    s2_wnoise_free, &
    s2_wnoise_gen_sky, &
    s2_wnoise_conv, &
    s2_wnoise_compute_alm, &
    s2_wnoise_map_convert, &
    s2_wnoise_downsample, &
    s2_wnoise_write_sky_file, &
    s2_wnoise_write_nobs_file, &
    s2_wnoise_write_std_file, &
    s2_wnoise_get_init, &
    s2_wnoise_get_nside, &
    s2_wnoise_get_type, &
    s2_wnoise_get_sky, &
    s2_wnoise_get_nobs, &
    s2_wnoise_get_sigma0, &
    s2_wnoise_get_std_const, &
    s2_wnoise_get_std_sky, &
    s2_wnoise_get_beam_status


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  interface s2_wnoise_init
    module procedure &
      s2_wnoise_init_const, &
      s2_wnoise_init_sky_file, &
      s2_wnoise_init_sky, &
      s2_wnoise_init_copy
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! To specify constant noise std over the sky.
  integer, public, parameter :: S2_WNOISE_TYPE_STD_CONST = 1
  !! To specify variable noise std over the sky.
  integer, public, parameter :: S2_WNOISE_TYPE_STD_SKY = 2


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_wnoise
     private
     logical :: init = .false.
     integer :: seed
     integer :: nside = 0
     integer :: type = S2_WNOISE_TYPE_STD_CONST
     type(s2_sky) :: sky
     type(s2_sky) :: nobs
     type(s2_sky) :: std_sky
     real(s2_sp) :: sigma0 = 0.0e0
     real(s2_sp) :: std_const = 0.0e0
     logical :: beam_status = .false.
  end type s2_wnoise


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_wnoise_init_const
    !
    !! Initialise wnoise with a constant std over the sky.
    !!
    !! Variables:
    !!   - std: Standard deviation of noise (constant over sky).
    !!   - nside: Healpix resolution of sky.
    !!   - seed: Seed for generating samples from uniform distribution.
    !!   - wnoise: Wnoise object initialised.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_init_const(std, nside, seed) result(wnoise)

      real(s2_sp), intent(in) :: std
      integer, intent(in) :: nside, seed
      type(s2_wnoise) :: wnoise
      
      ! Check object not already initialised.
      if(wnoise%init) then
        call s2_error(S2_ERROR_INIT, 's2_wnoise_init_const')
        return
      end if

      wnoise%type = S2_WNOISE_TYPE_STD_CONST
      wnoise%nside = nside
      wnoise%std_const = std
      wnoise%seed = seed
      wnoise%init = .true.

      call s2_wnoise_gen_sky(wnoise)

    end function s2_wnoise_init_const


    !--------------------------------------------------------------------------
    ! s2_wnoise_init_sky_file
    !
    !! Initialise a wnoise with a std that varies over the sky.  The std map  
    !! is calculated from the `number of observations' field contained in the
    !! fits file read.
    !!
    !! Notes:
    !!   - The wnoise healpix sky resolution will be set to the resolution of 
    !!     the nobs field contained in the fits file.  If a lower resolution
    !!     is later required, then the map resolution map be downsampled using 
    !!     the routine s2_wnoise_downsample.
    !!
    !! Variables:
    !!   - filename_nobs: Name of fits file containing nobs field.
    !!   - extension: Fits file extension specifying extension of nobs field in
    !!     fits file.
    !!   - sigma0: Noise parameter used to calculate std_sky from nobs field.
    !!   - seed: Seed for generating samples from uniform distribution.
    !!   - wnoise: Wnoise object initialised.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_init_sky_file(filename_nobs, extension, sigma0, seed) &
      result(wnoise)
      
      character(len=*), intent(in) :: filename_nobs
      integer, intent(in) :: extension
      real(s2_sp), intent(in) :: sigma0
      integer, intent(in) :: seed
      type(s2_wnoise) :: wnoise
      
      ! Check object not already initialised.
      if(wnoise%init) then
        call s2_error(S2_ERROR_INIT, 's2_wnoise_init_sky_file')
        return
      end if

      wnoise%type = S2_WNOISE_TYPE_STD_SKY
      wnoise%nobs = s2_sky_init(filename_nobs, S2_SKY_FILE_TYPE_MAP, &
        extension)
      wnoise%nside = s2_sky_get_nside(wnoise%nobs)
      wnoise%sigma0 = sigma0
      wnoise%seed = seed
      wnoise%init = .true.

      call s2_wnoise_compute_std_sky(wnoise)
      call s2_wnoise_gen_sky(wnoise)

    end function s2_wnoise_init_sky_file


    !--------------------------------------------------------------------------
    ! s2_wnoise_init_sky
    !
    !! Initialise a wnoise with a std that varies over the sky.  The std map  
    !! is calculated from the `number of observations' field contained in the
    !! passed s2_sky object.
    !!
    !! Notes:
    !!   - The wnoise healpix sky resolution will be set to the resolution of 
    !!     the nobs field contained in the fits file.  If a lower resolution
    !!     is later required, then the map resolution map be downsampled using 
    !!     the routine s2_wnoise_downsample.
    !!
    !! Variables:
    !!   - nobs: Sky object containing nobs field.
    !!   - sigma0: Noise parameter used to calculate std_sky from nobs field.
    !!   - seed: Seed for generating samples from uniform distribution.
    !!   - wnoise: Wnoise object initialised.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_init_sky(nobs, sigma0, seed) result(wnoise)

      type(s2_sky), intent(in) :: nobs
      real(s2_sp), intent(in) :: sigma0
      integer, intent(in) :: seed
      type(s2_wnoise) :: wnoise

      ! Check object not already initialised.
      if(wnoise%init) then
        call s2_error(S2_ERROR_INIT, 's2_wnoise_init_sky')
        return
      end if

      wnoise%type = S2_WNOISE_TYPE_STD_SKY
      wnoise%nobs = s2_sky_init(nobs)
      wnoise%nside = s2_sky_get_nside(wnoise%nobs)
      wnoise%sigma0 = sigma0
      wnoise%seed = seed
      wnoise%init = .true.

      call s2_wnoise_compute_std_sky(wnoise)
      call s2_wnoise_gen_sky(wnoise)

    end function s2_wnoise_init_sky


    !--------------------------------------------------------------------------
    ! s2_wnoise_init_copy
    !    
    !! Initialise a wnoise object as a copy of an original wnoise object. 
    !!
    !! Variables:
    !!   - orig: Original wnoise object to copy.
    !!   - copy: Initialised wnoise object as a copy of orig.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_init_copy(orig) result(copy)
 
      type(s2_wnoise), intent(in) :: orig
      type(s2_wnoise) :: copy

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_wnoise_init_copy')
        return
      end if

      ! Copy object atrributes.

      copy%seed = orig%seed
      copy%nside = orig%nside
      copy%type = orig%type
      copy%sky = s2_sky_init(orig%sky)

      ! Only copy nobs and std_sky skies if they are initialised 
      ! (i.e. if tyoe is S2_WNOISE_TYPE_STD_SKY).
      if(s2_sky_get_init(orig%nobs)) then
         copy%nobs = s2_sky_init(orig%nobs)
      end if
      if(s2_sky_get_init(orig%std_sky)) then
         copy%std_sky = s2_sky_init(orig%std_sky)
      end if
      
      copy%sigma0 = orig%sigma0
      copy%std_const = orig%std_const
      copy%beam_status = orig%beam_status

      copy%init = .true.

    end function s2_wnoise_init_copy
    

    !--------------------------------------------------------------------------
    ! s2_wnoise_free 
    !    
    !! Free all data associated with an initialised wnoise and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_free(wnoise)

      type(s2_wnoise), intent(inout) :: wnoise
      
      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_free')
      end if 
      
      if(s2_sky_get_init(wnoise%sky)) call s2_sky_free(wnoise%sky)
      if(s2_sky_get_init(wnoise%nobs)) call s2_sky_free(wnoise%nobs)
      if(s2_sky_get_init(wnoise%std_sky)) call s2_sky_free(wnoise%std_sky)

      wnoise%nside = 0
      wnoise%type = S2_WNOISE_TYPE_STD_CONST
      wnoise%sigma0 = 0.0e0
      wnoise%std_const = 0.0e0
      wnoise%beam_status = .false.

      wnoise%init = .false.

    end subroutine s2_wnoise_free


    !--------------------------------------------------------------------------
    ! s2_wnoise_compute_std_sky
    !    
    !! Compute the std sky map from nobs field: 
    !!   std_map(ipix) = sigma0/sqrt(nobs(ipix))
    !!    
    !! Notes:
    !!   - Object must be initialised before this routine is call.  If to be 
    !!     called from init function then must set init status first.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to calculate std_sky of.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_compute_std_sky(wnoise)

      type(s2_wnoise), intent(inout) :: wnoise

      real(s2_sp), allocatable :: temp_map(:)
      integer :: fail

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_free')
      end if 

      ! Get nobs map.
      allocate(temp_map(0:s2_sky_get_npix(wnoise%nobs)-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_wnoise_compute_std_sky')
      end if
      call s2_sky_get_map(wnoise%nobs, temp_map)

      ! Overwrite nobs_map with standard deviation that varies over sky
      ! and save as s2_sky type called std_sky.
      temp_map  = wnoise%sigma0 / sqrt(temp_map)

      ! Free wnoise%std_sky if already constructed.
      ! (When downsample must recompute std_sky so may already be constructed.)
      if(s2_sky_get_init(wnoise%std_sky)) call s2_sky_free(wnoise%std_sky)

      wnoise%std_sky = s2_sky_init(temp_map, s2_sky_get_nside(wnoise%nobs), &
        s2_sky_get_pix_scheme(wnoise%nobs))

      ! Free temporary memory used.
      deallocate(temp_map)

    end subroutine s2_wnoise_compute_std_sky


    !--------------------------------------------------------------------------
    ! s2_wnoise_gen_sky
    !    
    !! Generate a realisation of the wnoise sky satisfying the noise properties
    !! specified by the attributes of the wnoise object.
    !!    
    !! Notes:
    !!   - Object must be initialised before this routine is call.  If to be 
    !!     called from init function then must set init status first.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to generate sky realisation for.
    !!   - seed_in: Seed used to generate samples.  If not present then the 
    !!     wnoise seed attribute is used.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_gen_sky(wnoise, seed_in)
      
      type(s2_wnoise), intent(inout) :: wnoise
      integer, intent(in), optional :: seed_in
 
      real(s2_sp), allocatable :: sky_map(:)
      real(s2_sp) :: std
      integer :: fail, ipix, npix, pix_scheme, seed

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_gen_sky')
      end if 

      ! Set seed to use.
      if(present(seed_in)) then
         seed = seed_in
      else
         seed = wnoise%seed
      end if

      ! If sky already defined then free to create new sky realisation.
      if(s2_sky_get_init(wnoise%sky)) call s2_sky_free(wnoise%sky)

      ! Set constant noise dispersion if noise of this type.  Also set 
      ! pix_scheme depending on noise type.
      if(wnoise%type == S2_WNOISE_TYPE_STD_CONST) then
         std = wnoise%std_const
         pix_scheme = S2_SKY_RING  ! Say (could choose either).
      else if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         pix_scheme = s2_sky_get_pix_scheme(wnoise%nobs)
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, 's2_wnoise_gen_sky')
         return
      end if

      ! Compute npix for healpix nside resolution.
      npix = 12 * wnoise%nside**2

      ! Generate white noise over sky.
      allocate(sky_map(0:npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_wnoise_gen_sky')
      end if
      sky_map = 0.0e0
      do ipix = 0,npix-1
         ! Set variable noise dispersion for each pixel if noise of this type.
         if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
            std = s2_sky_get_map_pix(wnoise%std_sky, ipix)
         end if
         sky_map(ipix) = s2_distn_sample_gauss(seed, 0.0e0, std) 
      end do
      
      ! Save sky map as s2_sky type.
      wnoise%sky = s2_sky_init(sky_map, wnoise%nside, pix_scheme)

      ! Set beam status as false since new sky just calculated.
      wnoise%beam_status = .false.

      ! Free temporary storage.
      deallocate(sky_map)

    end subroutine s2_wnoise_gen_sky


    !--------------------------------------------------------------------------
    ! s2_wnoise_conv
    !
    !! Apply a beam by convolving it with the wnoise sky.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object containing sky to apply beam to.
    !!   - beam: Beam to apply.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_wnoise_conv(wnoise, beam)

      type(s2_wnoise), intent(inout) :: wnoise
      type(s2_pl), intent(in) :: beam

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_conv')
      end if 

      ! Convolve beam with wnoise sky.
      call s2_sky_conv(wnoise%sky, beam)

      ! Set beam status.
      wnoise%beam_status = .true.

    end subroutine s2_wnoise_conv
 

    !--------------------------------------------------------------------------
    ! s2_wnoise_compute_alm
    !
    !! Compute the alms of the wnoise sky at the lmax and mmax specified.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to compute alms of sky.
    !!   - lmax: Healpix spherical harmonic lmax.
    !!   - mmax: Healpix spherical harmonic mmax.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_wnoise_compute_alm(wnoise, lmax, mmax)

      type(s2_wnoise), intent(inout) :: wnoise
      integer, intent(in) :: lmax, mmax

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_compute_alm')
      end if 

      ! Compute wnoise sky alms.
      call s2_sky_compute_alm(wnoise%sky, lmax, mmax)

    end subroutine s2_wnoise_compute_alm


    !--------------------------------------------------------------------------
    ! s2_wnoise_map_convert    
    !
    !! Convert a wnoise sky map to the specified pixelisation scheme.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to convert sky of.
    !!   - pix_scheme: Pixelisation scheme to conver to sky map to.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_wnoise_map_convert(wnoise, pix_scheme)

      type(s2_wnoise), intent(inout) :: wnoise
      integer, intent(in) :: pix_scheme
      
      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_map_convert')
      end if 
      
      call s2_sky_map_convert(wnoise%sky, pix_scheme)
      call s2_sky_map_convert(wnoise%nobs, pix_scheme)
      call s2_sky_map_convert(wnoise%std_sky, pix_scheme)

    end subroutine s2_wnoise_map_convert


    !--------------------------------------------------------------------------
    ! s2_wnoise_downsample
    !    
    !! Downsample wnoise sky map to the specified nside_down.
    !!
    !! Notes:
    !!   - Actually generate a new wnoise realisation at the lower resolution 
    !!     (based on the new std map).
    !!
    !! Variables:
    !!   - wnoise: Wnoise containing sky to downsample.
    !!   - nside_down: New Healpix nside to downsample to.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_downsample(wnoise, nside_down)

      type(s2_wnoise), intent(inout) :: wnoise
      integer, intent(in) :: nside_down

      wnoise%nside = nside_down

      ! For full sky case downsample nobs field and compute new
      if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         call s2_sky_downsample(wnoise%nobs, nside_down)        
         call s2_wnoise_compute_std_sky(wnoise)
      end if

      call s2_wnoise_gen_sky(wnoise)

    end subroutine s2_wnoise_downsample


    !--------------------------------------------------------------------------
    ! s2_wnoise_write_sky_file
    !    
    !! Write the wnoise sky to an output fits file.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to write sky of.
    !!   - filename: Output fits filename.
    !!   - [comment]: Optional comment to be appended to the fits file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_write_sky_file(wnoise, filename, comment)

      type(s2_wnoise), intent(in) :: wnoise
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_write_sky_file')
      end if 

      call s2_sky_write_map_file(wnoise%sky, filename, comment)

    end subroutine s2_wnoise_write_sky_file


    !--------------------------------------------------------------------------
    ! s2_wnoise_write_nobs_file
    !    
    !! Write the wnoise nobs file, providing the wnoise type is
    !! S2_WNOISE_TYPE_STD_SKY (i.e. the noise std varies over the full sky).
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to write nobs field of.
    !!   - filename: Output fits filename.
    !!   - [comment]: Optional comment to be appended to the fits file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_write_nobs_file(wnoise, filename, comment)

      type(s2_wnoise), intent(in) :: wnoise
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_write_nobs_file')
      end if 

      if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         call s2_sky_write_map_file(wnoise%nobs, filename, comment)
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, &
           's2_wnoise_write_nobs_file', &
           comment_add='Std constant over full sky, nobs not used',&
           halt_in=.false.)
      end if

    end subroutine s2_wnoise_write_nobs_file


    !--------------------------------------------------------------------------
    ! s2_wnoise_write_std_file
    !    
    !! Write the wnoise std file, providing the wnoise type is
    !! S2_WNOISE_TYPE_STD_SKY (i.e. the noise std varies over the full sky).
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to write std_sky of.
    !!   - filename: Output fits filename.
    !!   - [comment]: Optional comment to be appended to the fits file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_wnoise_write_std_file(wnoise, filename, comment)

      type(s2_wnoise), intent(in) :: wnoise
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_write_std_file')
      end if 

      if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         call s2_sky_write_map_file(wnoise%std_sky, filename, comment)
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, &
           's2_wnoise_write_std_sky_file', &
           comment_add='Std constant over full sky', &
           halt_in=.false.)
      end if

    end subroutine s2_wnoise_write_std_file


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------
   
    !--------------------------------------------------------------------------
    ! s2_wnoise_get_init
    !
    !! Get init variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_init(wnoise) result(init)
      
      type(s2_wnoise), intent(in) :: wnoise
      logical :: init

      init = wnoise%init

    end function s2_wnoise_get_init


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_nside
    !
    !! Get nside variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - nside: Object nside variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_nside(wnoise) result(nside)
      
      type(s2_wnoise), intent(in) :: wnoise
      integer :: nside

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_nside')
      end if

      nside = wnoise%nside

    end function s2_wnoise_get_nside


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_type
    !
    !! Get type variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - type: Object type variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_type(wnoise) result(type)
      
      type(s2_wnoise), intent(in) :: wnoise
      integer :: type

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_type')
      end if

      type = wnoise%type

    end function s2_wnoise_get_type


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_sky
    !
    !! Get sky variable from the passed wnoise object.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the wnoise sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - sky: Object sky variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_sky(wnoise) result(sky)
      
      type(s2_wnoise), intent(in) :: wnoise
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_sky')
      end if

      sky = s2_sky_init(wnoise%sky)

    end function s2_wnoise_get_sky
 

    !--------------------------------------------------------------------------
    ! s2_wnoise_get_nobs
    !
    !! Get nobs variable from the passed wnoise object.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the wnoise nobs sky.
    !!   - The returned nobs sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - nobs: Object nobs variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_nobs(wnoise) result(nobs)
      
      type(s2_wnoise), intent(in) :: wnoise
      type(s2_sky) :: nobs

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_nobs')
      end if

      if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         nobs = s2_sky_init(wnoise%nobs)
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, &
           's2_wnoise_get_nobs', &
           comment_add='Std constant over full sky')
      end if

    end function s2_wnoise_get_nobs


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_sigma0
    !
    !! Get sigma0 variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - sigma0: Object sigma0 variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_sigma0(wnoise) result(sigma0)
      
      type(s2_wnoise), intent(in) :: wnoise
      real(s2_sp) :: sigma0

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_sigma0')
      end if

      sigma0 = wnoise%sigma0

    end function s2_wnoise_get_sigma0


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_std_const
    !
    !! Get std_const variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - std_const: Object std_const variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_std_const(wnoise) result(std_const)
      
      type(s2_wnoise), intent(in) :: wnoise
      real(s2_sp) :: std_const

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_std_const')
      end if

      if(wnoise%type == S2_WNOISE_TYPE_STD_CONST) then
         std_const = wnoise%std_const
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, &
           's2_wnoise_get_std_const', &
           comment_add='Std not constant over full sky')
      end if

    end function s2_wnoise_get_std_const


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_std_sky
    !
    !! Get std_sky variable from the passed wnoise object.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the wnoise std sky.
    !!   - The returned std sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - std_sky: Object std_sky variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_std_sky(wnoise) result(std_sky)
      
      type(s2_wnoise), intent(in) :: wnoise
      type(s2_sky) :: std_sky

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_std_sky')
      end if

      if(wnoise%type == S2_WNOISE_TYPE_STD_SKY) then
         std_sky = s2_sky_init(wnoise%std_sky)
      else
         call s2_error(S2_ERROR_WNOISE_TYPE_INVALID, &
           's2_wnoise_get_std_sky', &
           comment_add='Std constant over full sky')
      end if

    end function s2_wnoise_get_std_sky


    !--------------------------------------------------------------------------
    ! s2_wnoise_get_beam_status
    !
    !! Get beam_status variable from the passed wnoise object.
    !!
    !! Variables:
    !!   - wnoise: Wnoise object to get the variable of.
    !!   - beam_status: Object beam_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_wnoise_get_beam_status(wnoise) result(beam_status)
      
      type(s2_wnoise), intent(in) :: wnoise
      logical :: beam_status

      ! Check object initialised.
      if(.not. wnoise%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_wnoise_get_beam_status')
      end if

      beam_status = wnoise%beam_status

    end function s2_wnoise_get_beam_status


end module s2_wnoise_mod

  


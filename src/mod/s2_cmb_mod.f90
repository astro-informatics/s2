!------------------------------------------------------------------------------
! s2_cmb_mod -- S2 library cmb class
!
!! Provides functionality to create a Gaussian simuated CMB map.  The map 
!! is realised from Gaussain alms that satisfy the specified CMB power 
!! spectrum.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_cmb_mod

  use s2_types_mod
  use s2_error_mod
  use s2_distn_mod
  use s2_pl_mod
  use s2_sky_mod
  use s2_wnoise_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_cmb_init, &
    s2_cmb_free, &
    s2_cmb_gen_sky, &
    s2_cmb_map_convert, &
    s2_cmb_add_noise, &
    s2_cmb_write_sky, &
    s2_cmb_get_init, &
    s2_cmb_get_nside, &
    s2_cmb_get_npix, &
    s2_cmb_get_lmax, &
    s2_cmb_get_sky, &
    s2_cmb_get_clt, &
    s2_cmb_get_noise_added, &
    s2_cmb_get_beam_applied


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface s2_cmb_init
     module procedure &
       s2_cmb_init_pl, &
       s2_cmb_init_array, &
       s2_cmb_init_file, &
       s2_cmb_init_copy
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_cmb
     private
     logical :: init = .false.
     integer :: seed
     integer :: nside = 0
     type(s2_sky) :: sky
     type(s2_pl) :: clt
     logical :: noise_added = .false.
     logical :: beam_applied = .false.
  end type s2_cmb


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_cmb_init_pl
    !
    !! Initialise a Gaussian cmb form a clt spectrum pl object.
    !! If beam is present then the realised cmb alms are convolved with the 
    !! beam before the cmb map is generated.
    !!
    !! Variables:
    !!   - clt: Clt spectrum that fully defines the statistical 
    !!     characteristics of a cmb realisation as a pl object.
    !!   - nside: Sky Healpix resolution.
    !!   - seed: Seed for generating samples from a Gaussian distribution.
    !!   - [beam]: If present then the amb alms are convolved with this beam
    !!     before a map representation is constructed.
    !!   - [compute_map]: Logical to specify whether to compute the cmb map from
    !!     the alms (default is true).
    !!   - cmb: Cmb initialised.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !    August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_cmb_init_pl(clt, nside, seed, beam, compute_map) result(cmb)

      type(s2_pl), intent(in) :: clt
      integer, intent(in) :: nside, seed
      type(s2_pl), intent(in), optional :: beam
      logical, intent(in), optional :: compute_map
      type(s2_cmb) :: cmb

      ! Check object not already initialised.
      if(cmb%init) then
        call s2_error(S2_ERROR_INIT, 's2_cmb_init_pl')
        return
      end if

      cmb%nside = nside
      cmb%seed = seed
      cmb%clt = s2_pl_init(clt)
      cmb%init = .true.
      if(present(beam)) cmb%beam_applied = .true.
      call s2_cmb_gen_sky(cmb, beam=beam, compute_map=compute_map)

    end function s2_cmb_init_pl



    !--------------------------------------------------------------------------
    ! s2_cmb_init_array
    !
    !! Initialise a Gaussian cmb form a clt spectrum array.
    !! If beam is present then the realised cmb alms are convolved with the 
    !! beam before the cmb map is generated.
    !!
    !! Variables:
    !!   - clt_spec: Clt spectrum that fully defines the statistical 
    !!     characteristics of a cmb realisation.
    !!   - nside: Sky Healpix resolution.
    !!   - seed: Seed for generating samples from a Gaussian distribution.
    !!   - [beam]: If present then the amb alms are convolved with this beam
    !!     before a map representation is constructed.
    !!   - [compute_map]: Logical to specify whether to compute the cmb map from
    !!     the alms (default is true).
    !!   - cmb: Cmb initialised.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !    August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_cmb_init_array(clt_spec, nside, seed, beam, compute_map) &
      result(cmb)

      real(s2_sp), intent(in) :: clt_spec(:)
      integer, intent(in) :: nside, seed
      type(s2_pl), intent(in), optional :: beam
      logical, intent(in), optional :: compute_map
      type(s2_cmb) :: cmb

      ! Check object not already initialised.
      if(cmb%init) then
        call s2_error(S2_ERROR_INIT, 's2_cmb_init_array')
        return
      end if

      cmb%nside = nside
      cmb%seed = seed
      cmb%clt = s2_pl_init(clt_spec)
      cmb%init = .true.
      if(present(beam)) cmb%beam_applied = .true.
      call s2_cmb_gen_sky(cmb, beam=beam, compute_map=compute_map)

    end function s2_cmb_init_array


    !--------------------------------------------------------------------------
    ! s2_cmb_init_file
    !
    !! Initialise a Gaussian cmb from a file containing a clt spectrum.
    !! If beam is present then the realised cmb alms are convolved with the 
    !! beam before the cmb map is generated.
    !!
    !! Variables:
    !!   - filename_clt: Name of the file containing the clt spectrum.
    !!   - nside: Sky Healpix resolution.
    !!   - lmin: Minimum l value read from file.  (The first value read will 
    !!     be saved a this l.  Clt for all lower values of l is set to zero. 
    !!     Typically lmin=2 so monopole and dipole not included.)
    !!   - lmax: Maximum l value of clt.
    !!   - ncomment: Number of comment lines in file to ignore before start 
    !!     reading clt values.
    !!   - seed: Seed for generating samples from a Gaussian distribution.
    !!   - [beam]: If present then the amb alms are convolved with this beam
    !!     before a map representation is constructed.
    !!   - [compute_map]: Logical to specify whether to compute the cmb map from
    !!     the alms (default is true).
    !!   - cmb: Cmb initialised.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_init_file(filename_clt, nside, lmin, lmax, ncomment, &
      seed, beam, scale_cl, compute_map) result(cmb)

      character(len=*), intent(in) :: filename_clt
      integer, intent(in) :: nside, lmin, lmax, ncomment, seed
      type(s2_pl), intent(in), optional :: beam
      logical, intent(in), optional :: scale_cl
      logical, intent(in), optional :: compute_map
      type(s2_cmb) :: cmb

      logical :: scale_cl_use = .true.

!real(s2_sp) :: power

      if(present(scale_cl)) scale_cl_use = scale_cl
      
      ! Check object not already initialised.
      if(cmb%init) then
        call s2_error(S2_ERROR_INIT, 's2_cmb_init_file')
        return
      end if
      
      ! Initialise object variables and generate alms and map 
      ! (i.e. generate sky).
      cmb%nside = nside
      cmb%seed = seed
      cmb%clt = s2_pl_init(filename_clt, lmin, lmax, ncomment, scale_cl_use)

!power = s2_pl_power(cmb%clt)
!write(*,*) 'cmb power cl: ', power / 1e6

      cmb%init = .true.
      if(present(beam)) cmb%beam_applied = .true.
      call s2_cmb_gen_sky(cmb, beam=beam, compute_map=compute_map)

!power = s2_sky_power_map(cmb%sky)
!write(*,*) 'cmb power sky: ', power

!power = s2_sky_power_alm(cmb%sky)
!write(*,*) 'cmb power alm: ', power

    end function s2_cmb_init_file


    !--------------------------------------------------------------------------
    ! s2_cmb_init_copy
    !
    !! Initialise a cmb from a copy of an original cmb.
    !!
    !! Variables:
    !!   - orig: Original cmb to be copied.
    !!   - copy: Initialised cmb with attributes copied from orig.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_init_copy(orig) result(copy)

      type(s2_cmb), intent(in) :: orig
      type(s2_cmb) :: copy

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_cmb_init_copy')
        return
      end if

      ! Copy object atrributes.

      copy%sky = s2_sky_init(orig%sky)
      copy%clt = s2_pl_init(orig%clt)

      copy%seed = orig%seed
      copy%nside = orig%nside
      copy%noise_added = orig%noise_added
      copy%beam_applied = orig%beam_applied

      copy%init = .true.

    end function s2_cmb_init_copy


    !--------------------------------------------------------------------------
    ! s2_cmb_free
    !
    !! Free all data associated with an initialised cmb and reset all other 
    !! attributes.
    !
    !! Variables:
    !!   - cmb: Cmb to be freed.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_free(cmb)

      type(s2_cmb), intent(inout) :: cmb

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_free')
      end if 

      if(s2_sky_get_init(cmb%sky)) call s2_sky_free(cmb%sky)
      if(s2_pl_get_init(cmb%clt)) call s2_pl_free(cmb%clt)
      
      cmb%nside = 0
      cmb%noise_added = .false.
      cmb%beam_applied = .false.
      cmb%init = .false.

    end subroutine s2_cmb_free


    !--------------------------------------------------------------------------
    ! s2_cmb_gen_sky
    !
    !! Generate a sky realisation of the cmb from the clt spectrum.  
    !! First Gaussian alm values that satisfy the clt spectrum are generated, 
    !! these are convolved with a beam (if present), before the map realisation
    !! is computed from an inverse spherical harmonic transform.
    !!
    !! Notes:
    !!   - Object must be initialised before this routine is call.  If to be 
    !!     called from init function then must set init status first.
    !!   - Map is converted here to have units of mK. -- Not any more!
    !!
    !! Variables:
    !!   - cmb: Cmb to generate sky map of.
    !!   - [seed]: Seed for sampling Gaussian distribution.  If not present 
    !!     the cmb seed attribute will be used.
    !!   - [beam]: Beam function.  No beam is convolved if not present.
    !!   - [compute_map]: Logical to specify whether to compute the cmb map from
    !!     the alms (default is true).
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    ! 
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_gen_sky(cmb, seed, beam, compute_map)

      type(s2_cmb), intent(inout) :: cmb
      integer, intent(in), optional :: seed
      type(s2_pl), intent(in), optional :: beam
      logical, intent(in), optional :: compute_map

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_gen_sky')
      end if 

      call s2_cmb_create_alm(cmb, seed)
      if(present(beam)) call s2_cmb_apply_beam(cmb, beam)

      if(present(compute_map)) then
        if(compute_map) call s2_cmb_compute_map(cmb)
      else
        call s2_cmb_compute_map(cmb)
      end if

      ! Convert map to from uK to mK.
!write(*,*) 'WARNING: NO LONGER CONVERTING TO mK IN S2_CMB'
!      call s2_sky_scale(cmb%sky, 1.0e-3)

    end subroutine s2_cmb_gen_sky


    !--------------------------------------------------------------------------
    ! s2_cmb_apply_beam
    !
    !! Convolve beam with cmb alms.  The cmb alms are overwritten with the 
    !! convolved alms.
    !!
    !! Variables:
    !!   - cmb: Cmb to apply beam to.
    !!   - beam: Beam to apply.
    !   
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_apply_beam(cmb, beam)

      type(s2_cmb), intent(inout) :: cmb
      type(s2_pl), intent(in) :: beam

      integer, parameter :: DEFAUL_SKY_PIX_SCHEME = S2_SKY_RING
      integer :: lmax, mmax, pix_scheme, fail
      complex(s2_spc), allocatable :: alm(:,:)

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_apply_beam')
      end if 

      ! Check sky initialised.
      if(.not. s2_sky_get_init(cmb%sky)) then
         call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_apply_beam', &
           comment_add='Sky not initialised')
         return
      end if

      ! Set local size parameters before wipe.
      lmax = s2_sky_get_lmax(cmb%sky)
      mmax = s2_sky_get_mmax(cmb%sky)
      pix_scheme = s2_sky_get_pix_scheme(cmb%sky)
      
      ! Allocate temp alm array to calculate.
      allocate(alm(0:lmax,0:mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_cmb_apply_beam')
      end if

      ! Get the alm and convolve with beam.
      call s2_sky_get_alm(cmb%sky, alm)
      call s2_pl_conv(beam, alm)               ! Alm overwritten on output.

      ! Set convolved alms as new sky.
      call s2_sky_free(cmb%sky)
      cmb%sky = s2_sky_init(alm, lmax, mmax, cmb%nside, &
              pix_scheme)

      ! Free temporary alm memory used.
      deallocate(alm)

    end subroutine s2_cmb_apply_beam


    !--------------------------------------------------------------------------
    ! s2_cmb_create_alm
    !
    !! Create Gaussian cmb alms that satify the clt spectrum.
    !!
    !! Variables:
    !!   - cmb: Cmb to create alms for.
    !!   - [seed_in]: Seed to generate samples from a Gaussian distribution.  If
    !!     not present then the cmb seed attribute is used.
    !   
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_create_alm(cmb, seed_in)

      type(s2_cmb), intent(inout) :: cmb
      integer, intent(in), optional :: seed_in

      complex(s2_spc), allocatable :: alm(:,:)
      real(s2_sp) :: hsqrt2, std
      integer :: lmax, mmax, seed, fail
      integer :: l, m
      integer, parameter :: DEFAUL_SKY_PIX_SCHEME = S2_SKY_RING

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_create_alm')
      end if 

      ! Set seed to use.
      if(present(seed_in)) then
         seed = seed_in
      else
         seed = cmb%seed
      end if
      
      ! Ensure seed negative to initialise random deviate generators.
      if (seed > 0) seed = -seed

      ! Remove any old sky.
      if(s2_sky_get_init(cmb%sky)) call s2_sky_free(cmb%sky)

      ! Set local lmax and mmax variables.
      lmax = s2_pl_get_lmax(cmb%clt)
      mmax = lmax

      ! Allocate temp alm array to calculate.
      allocate(alm(0:lmax,0:mmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky_init_map')
      end if
      alm(0:lmax,0:mmax) = cmplx(0e0, 0e0)

      ! Create Gaussian alm based on clt spectrum.
      hsqrt2 = 1.0e0 / sqrt(2.0e0)
      do l = 0,lmax
         
         if(s2_pl_get_spec_l(cmb%clt, l) < 0) then
           std = 0e0
         else
           std  = sqrt(s2_pl_get_spec_l(cmb%clt, l))
         end if

         ! m = 0 case
         alm(l,0) = cmplx( &
           s2_distn_sample_gauss(seed, 0.0e0, std), &
           0.0e0 )
       
         ! m > 0 case
         do m = 1,l
            alm(l,m) = cmplx( &
              s2_distn_sample_gauss(seed, 0.0e0, std * hsqrt2), &
              s2_distn_sample_gauss(seed, 0.0e0, std * hsqrt2) )
         end do
      
      end do

      ! Initialise sky from alm.
      cmb%sky = s2_sky_init(alm(0:lmax,0:mmax), lmax, mmax, cmb%nside, &
        DEFAUL_SKY_PIX_SCHEME)

      ! Free temporary storage space used.
      deallocate(alm)

    end subroutine s2_cmb_create_alm


    !--------------------------------------------------------------------------
    ! s2_cmb_compute_map
    !
    !! Compute the cmb sky map from alms.
    !!
    !! Variables:
    !!   - cmb: Cmb to compute map of.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_compute_map(cmb)

      type(s2_cmb), intent(inout) :: cmb
      
      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_compute_map')
      end if 
      
      call s2_sky_compute_map(cmb%sky)

    end subroutine s2_cmb_compute_map


    !--------------------------------------------------------------------------
    ! s2_cmb_map_convert
    !
    !! Convert sky map pixelisation scheme between nested and ring.
    !!
    !! Variables:
    !!   - cmb: Cmb to convert pixelisation scheme of.
    !!   - pix_scheme: Pixelisation scheme to convert to.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_map_convert(cmb, pix_scheme)

      type(s2_cmb), intent(inout) :: cmb
      integer, intent(in) :: pix_scheme
      
      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_map_convert')
      end if 
      
      call s2_sky_map_convert(cmb%sky, pix_scheme)

    end subroutine s2_cmb_map_convert


    !--------------------------------------------------------------------------
    ! s2_cmb_add_noise
    !
    !! Add noise to the cmb sky map.  The cmb noise_added attribute flag is 
    !! also set.
    !!
    !! Notes:
    !!   - When noise is added to the cmb sky, the old sky is overwritten with
    !!     the new one (the sum of the old one and the noise sky).  Since the
    !!     old alms correspond only to the old simulated cmb map, and not the
    !!     map plus noise, they are wiped from the resultant sky saved.
    !!
    !! Variables:
    !!   - cmb: Cmb to add noise to.
    !!   - wnoise: White noise object containing noise realisation to add to 
    !!     cmb sky map.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_cmb_add_noise(cmb, wnoise)

      type(s2_cmb), intent(inout) :: cmb
      type(s2_wnoise), intent(in) :: wnoise  
        ! Inout since must be able to compute sky map if doesn't already exist.

      type(s2_sky) :: wnoise_sky

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_add_noise')
      end if

      wnoise_sky = s2_wnoise_get_sky(wnoise)
      cmb%sky = s2_sky_add(cmb%sky, wnoise_sky)

      cmb%noise_added = .true.

      call s2_sky_free(wnoise_sky)

    end subroutine s2_cmb_add_noise


    !--------------------------------------------------------------------------
    ! s2_cmb_write_sky
    !
    !! Write the cmb sky to an output fits file.
    !!
    !! Variables:
    !!   - cmb: Cmb to be written to a fits file.
    !!   - filename: Name of the outout fits file.
    !!   - [comment]: Optional comment to add to header of output file.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine s2_cmb_write_sky(cmb, filename, comment)

      type(s2_cmb), intent(in) :: cmb
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_write_sky')
      end if 

      call s2_sky_write_map_file(cmb%sky, filename, comment)

    end subroutine s2_cmb_write_sky


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------
       
    !--------------------------------------------------------------------------
    ! s2_cmb_get_init
    !
    !! Get init variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_init(cmb) result(init)

      type(s2_cmb), intent(in) :: cmb
      logical :: init

      init = cmb%init

    end function s2_cmb_get_init


    !--------------------------------------------------------------------------
    ! s2_cmb_get_nside
    !
    !! Get nside variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - nside: Object nside variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_nside(cmb) result(nside)

      type(s2_cmb), intent(in) :: cmb
      integer :: nside

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_nside')
      end if 

      nside = cmb%nside

    end function s2_cmb_get_nside


    !--------------------------------------------------------------------------
    ! s2_cmb_get_npix
    !
    !! Get npix variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - npix: Object npix variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_npix(cmb) result(npix)

      type(s2_cmb), intent(in) :: cmb
      integer :: npix

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_npix')
      end if 

      npix = s2_sky_get_npix(cmb%sky)

    end function s2_cmb_get_npix


    !--------------------------------------------------------------------------
    ! s2_cmb_get_lmax
    !
    !! Get lmax variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - lmax: Object lmax variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_lmax(cmb) result(lmax)

      type(s2_cmb), intent(in) :: cmb
      integer :: lmax

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_lmax')
      end if 

      lmax = s2_pl_get_lmax(cmb%clt)

    end function s2_cmb_get_lmax


    !--------------------------------------------------------------------------
    ! s2_cmb_get_sky
    !
    !! Get sky variable from the passed cmb.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the cmb sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - sky: Object sky variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_sky(cmb) result(sky)

      type(s2_cmb), intent(in) :: cmb
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_sky')
      end if 

      ! Make a copy for the returned sky.
      sky = s2_sky_init(cmb%sky)

    end function s2_cmb_get_sky


    !--------------------------------------------------------------------------
    ! s2_cmb_get_clt
    !
    !! Get clt variable from the passed cmb.
    !!
    !! Notes:
    !!   - Initialises a new clt as a copy of the cmb clt.
    !!   - The returned clt is subsequently independed of the clt stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - clt: Object clt variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_clt(cmb) result(clt)

      type(s2_cmb), intent(in) :: cmb
      type(s2_pl) :: clt

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_clt')
      end if 

      clt = s2_pl_init(cmb%clt)

    end function s2_cmb_get_clt


    !--------------------------------------------------------------------------
    ! s2_cmb_get_noise_added
    !
    !! Get noise_added variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - noise_added: Object noise_added variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_noise_added(cmb) result(noise_added)

      type(s2_cmb), intent(in) :: cmb
      logical :: noise_added

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_noise_added')
      end if 

      noise_added = cmb%noise_added

    end function s2_cmb_get_noise_added


    !--------------------------------------------------------------------------
    ! s2_cmb_get_beam_applied
    !
    !! Get beam_applied variable from the passed cmb.
    !!
    !! Variables:
    !!   - cmb: Cmb object to get the variable of.
    !!   - beam_applied: Object beam_applied variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function s2_cmb_get_beam_applied(cmb) result(beam_applied)

      type(s2_cmb), intent(in) :: cmb
      logical :: beam_applied

      ! Check object initialised.
      if(.not. cmb%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_cmb_get_beam_applied')
      end if 

      beam_applied = cmb%beam_applied

    end function s2_cmb_get_beam_applied


end module s2_cmb_mod

!------------------------------------------------------------------------------
! s2_gcmbcoad
!
!! Generate co-added Gaussian CMB and noise realisations.  Realisations are
!! simulated for each WMAP channel by convolving with the approriate beam
!! function and adding appropriate noise for the specific receiver.  The
!! signals from each receiver are combined to give a coadded map, using the
!! processing pipeline specified by Komatsu et al. (2003).
!!
!! Usage: s2_gcmbcoad
!!   - [-help]: Display usage information.
!!   - [-isim_start isim_start]: Index of first simulation to compute.
!!   - [-isim_end isim_end (optional)]: Index fo last simulation to compute.
!!     (If not specified only one simualtion is made for index isim_start.)
!!   - [-dir_cl dir_cl]: Directory containing CMB power specftrum file.
!!   - [-dir_nobs dir_nobs (set to dir_cl if not specified)]: Directory
!!     containing WMAP map files with NOBS field in extension 2.
!!   - [-dir_beam dir_beam (set to dir_cl if not specified)]: Directory
!!     containting beam transfer function files.
!!   - [-mask filename_mask (optional)]: Mask file to apply.  If not
!!     sepcified no mask is applied.
!!   - [-lmax lmax]: Maximum harmonic l to consider.
!!   - [-out_cmb_prefix filename_out_cmb_prefix]: Prefix of output
!!     coadded cmb map file names (simulation index and .fits automatically
!!     appended to end of name).
!!   - [-out_noise_prefix filename_out_noise_prefix]: Prefix of output
!!     coadded noise map file names (simulation index and .fits automatically
!!     appended to end of name).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - October 2005
!
! Revisions:
!   October 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_gcmbcoad

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2_cmb_mod
  use s2_wnoise_mod

  implicit none

  interface
    function s2_gcmbcoad_weightedsum(skies, weights) result(sky_coadded)
      use s2_sky_mod
      use s2_error_mod
      implicit none
      type(s2_sky), intent(inout) :: skies(1:)
      type(s2_sky), intent(inout) :: weights(1:)
      type(s2_sky) :: sky_coadded
    end function s2_gcmbcoad_weightedsum
  end interface
  
  integer, parameter :: NUM_BANDS = 8
  integer, parameter :: NOBS_EXT = 2
  integer, parameter :: LMAX_LIMIT = 1000
  integer, parameter :: NSIDE_FINAL = 256
  character(len=*), parameter :: msgpfx = 's2_gcmbcoad> '

#ifdef WMAP5
  real(s2_sp), parameter :: SIGMA0(NUM_BANDS) = (/ &
       2.254, 2.141, 3.314, &
       2.953, 5.899, 6.565, &
       6.926, 6.761 /)

  character(len=S2_STRING_LEN), parameter :: filename_nobs(NUM_BANDS) = (/ &
    'wmap_imap_r9_5yr_Q1_v3.fits', &
    'wmap_imap_r9_5yr_Q2_v3.fits', &
    'wmap_imap_r9_5yr_V1_v3.fits', &
    'wmap_imap_r9_5yr_V2_v3.fits', &
    'wmap_imap_r9_5yr_W1_v3.fits', &
    'wmap_imap_r9_5yr_W2_v3.fits', &
    'wmap_imap_r9_5yr_W3_v3.fits', &
    'wmap_imap_r9_5yr_W4_v3.fits'   /)

  character(len=S2_STRING_LEN), parameter :: filename_beam(NUM_BANDS) = (/ &
    'wmap_Q1_ampl_bl_5yr_v3.txt', &
    'wmap_Q2_ampl_bl_5yr_v3.txt', &
    'wmap_V1_ampl_bl_5yr_v3.txt', &
    'wmap_V2_ampl_bl_5yr_v3.txt', &
    'wmap_W1_ampl_bl_5yr_v3.txt', &
    'wmap_W2_ampl_bl_5yr_v3.txt', &
    'wmap_W3_ampl_bl_5yr_v3.txt', &
    'wmap_W4_ampl_bl_5yr_v3.txt'    /)

  character(S2_STRING_LEN) :: filename_mask
  character(S2_STRING_LEN) :: filename_cl = 'wmap_lcdm_sz_lens_wmap5_cl_v3.dat'
  character(S2_STRING_LEN) :: filename_out
  character(S2_STRING_LEN) :: filename_out_cmb_prefix = 'gcmb5_coadded_'
  character(S2_STRING_LEN) :: filename_out_noise_prefix = 'wnoise5_coadded_'
  character(S2_STRING_LEN) :: comment_out = 'seed not specified'
  character(len=S2_STRING_LEN) :: filename_current
  character(len=S2_STRING_LEN) :: dir_cl = 'data_input'
  character(len=S2_STRING_LEN) :: dir_nobs = 'data_input'
  character(len=S2_STRING_LEN) :: dir_beam = 'data_input'

!!$#elifdef WMAP3
!!$  real(s2_sp), parameter :: SIGMA0(NUM_BANDS) = (/ &
!!$    2.2449, 2.1347, 3.3040, &
!!$    2.9458, 5.8833, 6.5324, &
!!$    6.8849, 6.7441 /)
!!$
!!$  character(len=S2_STRING_LEN), parameter :: filename_nobs(NUM_BANDS) = (/ &
!!$    'wmap_imap_r9_3yr_Q1_v2.fits', &
!!$    'wmap_imap_r9_3yr_Q2_v2.fits', &
!!$    'wmap_imap_r9_3yr_V1_v2.fits', &
!!$    'wmap_imap_r9_3yr_V2_v2.fits', &
!!$    'wmap_imap_r9_3yr_W1_v2.fits', &
!!$    'wmap_imap_r9_3yr_W2_v2.fits', &
!!$    'wmap_imap_r9_3yr_W3_v2.fits', &
!!$    'wmap_imap_r9_3yr_W4_v2.fits'   /)
!!$
!!$  character(len=S2_STRING_LEN), parameter :: filename_beam(NUM_BANDS) = (/ &
!!$    'wmap_q1_ampl_bl_3yr_v2.txt', &
!!$    'wmap_q2_ampl_bl_3yr_v2.txt', &
!!$    'wmap_v1_ampl_bl_3yr_v2.txt', &
!!$    'wmap_v2_ampl_bl_3yr_v2.txt', &
!!$    'wmap_w1_ampl_bl_3yr_v2.txt', &
!!$    'wmap_w2_ampl_bl_3yr_v2.txt', &
!!$    'wmap_w3_ampl_bl_3yr_v2.txt', &
!!$    'wmap_w4_ampl_bl_3yr_v2.txt'    /)
!!$
!!$  character(S2_STRING_LEN) :: filename_mask
!!$  character(S2_STRING_LEN) :: filename_cl = 'wmap_comb_tt_powspec_3yr_v2.txt'
!!$  character(S2_STRING_LEN) :: filename_out
!!$  character(S2_STRING_LEN) :: filename_out_cmb_prefix = 'gcmb3_coadded_'
!!$  character(S2_STRING_LEN) :: filename_out_noise_prefix = 'wnoise3_coadded_'
!!$  character(S2_STRING_LEN) :: comment_out = 'seed not specified'
!!$  character(len=S2_STRING_LEN) :: filename_current
!!$  character(len=S2_STRING_LEN) :: dir_cl = 'data_input'
!!$  character(len=S2_STRING_LEN) :: dir_nobs = 'data_input'
!!$  character(len=S2_STRING_LEN) :: dir_beam = 'data_input'

#else

  ! See http://cmbdata.gsfc.nasa.gov/product/map/IMaps_cleaned.cfm for sigma0
  ! values and definition.
  real(s2_sp), parameter :: SIGMA0(NUM_BANDS) = (/ &
    2.26677, 2.15567, 3.28789, &
    2.93683, 5.85196, 6.53276, &
    6.88032, 6.72537 /)

  character(len=S2_STRING_LEN), parameter :: filename_nobs(NUM_BANDS) = (/ &
    'wmap_q1_cleanimap_yr1_v1.fits', &
    'wmap_q2_cleanimap_yr1_v1.fits', &
    'wmap_v1_cleanimap_yr1_v1.fits', &
    'wmap_v2_cleanimap_yr1_v1.fits', &
    'wmap_w1_cleanimap_yr1_v1.fits', &
    'wmap_w2_cleanimap_yr1_v1.fits', &
    'wmap_w3_cleanimap_yr1_v1.fits', &
    'wmap_w4_cleanimap_yr1_v1.fits'   /)

 character(len=S2_STRING_LEN), parameter :: filename_beam(NUM_BANDS) = (/ &
    'map_q1_ampl_bl_yr1_v1.txt', &
    'map_q2_ampl_bl_yr1_v1.txt', &
    'map_v1_ampl_bl_yr1_v1.txt', &
    'map_v2_ampl_bl_yr1_v1.txt', &
    'map_w1_ampl_bl_yr1_v1.txt', &
    'map_w2_ampl_bl_yr1_v1.txt', &
    'map_w3_ampl_bl_yr1_v1.txt', &
    'map_w4_ampl_bl_yr1_v1.txt'    /)

  character(S2_STRING_LEN) :: filename_mask
  character(S2_STRING_LEN) :: filename_cl = 'wmap_lcdm_pl_model_yr1_v1.txt'
  character(S2_STRING_LEN) :: filename_out
  character(S2_STRING_LEN) :: filename_out_cmb_prefix = 'gcmb_coadded_'
  character(S2_STRING_LEN) :: filename_out_noise_prefix = 'wnoise_coadded_'
  character(S2_STRING_LEN) :: comment_out = 'seed not specified'
  character(len=S2_STRING_LEN) :: filename_current
  character(len=S2_STRING_LEN) :: dir_cl = 'data_input'
  character(len=S2_STRING_LEN) :: dir_nobs = 'data_input'
  character(len=S2_STRING_LEN) :: dir_beam = 'data_input'

#endif

  integer :: isim_start = 1, isim_end = 1
  integer :: lmax, lmin_cl, lmin_beam
  integer :: iband, isim, fail
  integer :: ncomment_cl, ncomment_beam
  integer :: seed_cmb = 1, seed_noise = 1
  logical :: scale_in_cl = .true.      
  logical :: scale_in_beam = .false.
  logical :: line_nos_in_cl = .true.
  logical :: line_nos_in_beam = .true.
  logical :: apply_mask = .false.
  
  type(s2_pl) :: cl
  type(s2_pl) :: beam(NUM_BANDS)
  type(s2_sky) :: weights(NUM_BANDS)
  type(s2_sky) :: weights_sum, weights_temp
  type(s2_sky) :: mask
  type(s2_sky) :: nobs(NUM_BANDS)
  type(s2_cmb) :: cmb
  type(s2_sky) :: cmb_sky
  type(s2_sky) :: cmb_sky_beam(NUM_BANDS)
  type(s2_sky) :: sky_temp
  type(s2_sky) :: cmb_coadded, cmb_coadded_unscaled
  type(s2_wnoise) :: wnoise
  type(s2_sky) :: wnoise_sky(NUM_BANDS)
  type(s2_sky) :: wnoise_coadded, wnoise_coadded_unscaled
  real(s2_sp), allocatable :: nobs_map(:), weights_map(:)

#ifdef DEBUG
  integer :: l
#endif

#ifdef MAKE_COADDED_DATA_MAP
  type(s2_sky) :: data_map(NUM_BANDS)
  type(s2_sky) :: data_map_coadded_unscaled
  type(s2_sky) :: data_map_coadded
#endif  

#ifdef MAKE_COADDED_NOISE_DISP_MAP
  real(s2_sp), allocatable :: ones_map(:)
  type(s2_sky) :: sky_ones
  type(s2_sky) :: noise_disp
#endif


  !----------------------------------------------------------------
  ! Initialise variables.
  !----------------------------------------------------------------

  ! Set default parameter values.
  ! (Not necessary to reset some variables but all done here for
  ! completeness and easy of reference.)
  isim_start = 1
  isim_end = 1
  lmax = 1000
  lmin_cl = 2
  lmin_beam = 0
  ncomment_cl = 29
#ifdef WMAP3
  ncomment_cl = 45
#endif
#ifdef WMAP5
  ncomment_cl = 0
#endif

  ncomment_beam = 8
  scale_in_cl = .true.
  scale_in_beam = .false.
  line_nos_in_cl = .true.
  line_nos_in_beam = .true.
  seed_cmb = 1
  seed_noise = 1
  apply_mask = .false.
  
  ! Parse options.
  call parse_options()

  ! Check lmax not outside limit.
  if(lmax > LMAX_LIMIT) then
    call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_gcmbcoad', &
      comment_add='lmax outside range specified in input files')
  end if

  ! Check isim numbers valid.
  if(isim_start > isim_end) then
    call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_gcmbcoad', &
      comment_add='Simulation number range invalid')
  end if


  !----------------------------------------------------------------
  ! Read in power spectra of cmb and beams
  ! (and mask if required).
  !----------------------------------------------------------------

  ! Read in cl spectrum.
  write(filename_current,'(a,a,a)') trim(dir_cl), '/', trim(filename_cl)
  cl = s2_pl_init(filename_current, lmin_cl, lmax, &
    ncomment_cl, scale_in_cl, line_nos_in_cl)

#ifdef DEBUG
    write(*,*) 'Cl spectrum file: ', trim(filename_current)
    do l=0,lmax
      write(*,*) 'l=', l, ': ', s2_pl_get_spec_l(cl, l)
    end do
    write(*,*)
#endif

  ! Read in power spectrum for each beam.
  do iband = 1, NUM_BANDS
    write(filename_current,'(a,a,a)') trim(dir_beam), '/', &
      trim(filename_beam(iband))
    beam(iband) = s2_pl_init(filename_current, lmin_beam, lmax, &
      ncomment_beam, scale_in_beam, line_nos_in_beam)

#ifdef DEBUG
    write(*,*) 'Band file: ', trim(filename_current)
    write(*,*) 'Band no: ', iband
    do l=0,lmax
      write(*,*) 'l=', l, ': ', s2_pl_get_spec_l(beam(iband), l)
    end do
    write(*,*)
#endif

  end do

  ! Read in mask.
  if(apply_mask) then
    mask = s2_sky_init(filename_mask, S2_SKY_FILE_TYPE_MAP)
  end if


  !----------------------------------------------------------------
  ! Read in Nobs map for each receiver and compute noise dispersion
  ! and weight maps.
  !----------------------------------------------------------------
  
  do iband = 1, NUM_BANDS

    ! Read in Nobs field.
    write(filename_current,'(a,a,a)') trim(dir_nobs), '/', &
      trim(filename_nobs(iband))

    nobs(iband) = s2_sky_init(filename_current, S2_SKY_FILE_TYPE_MAP, NOBS_EXT)

    ! Convert nobs field to ring pixelisation.
    ! (Not stictly necessary to do here, since routines will do it when
    ! necessary, but better to just perform the conversion once here.)
    call s2_sky_map_convert(nobs(iband), S2_SKY_RING)

    ! Get map values.
    allocate(nobs_map(0:s2_sky_get_npix(nobs(iband))-1), stat=fail)
    allocate(weights_map(0:s2_sky_get_npix(nobs(iband))-1), stat=fail)
    if(fail /= 0) then
      call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_gcmbcoad')
    end if
    call s2_sky_get_map(nobs(iband), nobs_map)

    ! Compute weights for each band.
    weights_map = nobs_map / SIGMA0(iband)**2
    weights(iband) = s2_sky_init(weights_map, &
      s2_sky_get_nside(nobs(iband)), s2_sky_get_pix_scheme(nobs(iband)))

    ! Compute weight sum.
    if(iband == 1) then
      weights_sum = s2_sky_init(weights(iband))
    else
      weights_temp = s2_sky_init(weights_sum)
      call s2_sky_free(weights_sum)
      weights_sum = s2_sky_add(weights(iband), weights_temp)
      call s2_sky_free(weights_temp)
    end if
    deallocate(nobs_map, weights_map)

  end do


  !----------------------------------------------------------------
  ! Make coadded data map if required.
  !----------------------------------------------------------------

#ifdef MAKE_COADDED_DATA_MAP

  ! Read in data maps.
  do iband = 1, NUM_BANDS
    write(filename_current,'(a,a,a)') trim(dir_nobs), '/', &
      trim(filename_nobs(iband))
write(*,*) trim(filename_current)
    data_map(iband) = s2_sky_init(filename_current, S2_SKY_FILE_TYPE_MAP, 1)
    call s2_sky_map_convert(data_map(iband), S2_SKY_RING)
  end do

  ! Combine to produce coadded map
  data_map_coadded_unscaled = s2_gcmbcoad_weightedsum(data_map, weights)
  data_map_coadded = s2_sky_product(data_map_coadded_unscaled, weights_sum, divide=.true.)

  ! Downsample map.
  call s2_sky_downsample(data_map_coadded, NSIDE_FINAL)

  ! Write coadded data map.
  call s2_sky_write_map_file(data_map_coadded, 'wmap5_coadded.fits')

  ! Free memory used.
  do iband = 1, NUM_BANDS
    call s2_sky_free(data_map(iband))
  end do
  call s2_sky_free(data_map_coadded_unscaled)
  call s2_sky_free(data_map_coadded)

#endif


  !----------------------------------------------------------------
  ! Make coadded noise dispersion map if required.
  !----------------------------------------------------------------

#ifdef MAKE_COADDED_NOISE_DISP_MAP

  allocate(ones_map(0:s2_sky_get_npix(weights_sum)-1), stat=fail)
  if(fail /= 0) then
    call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_gcmbcoad')
  end if
  ones_map(0:s2_sky_get_npix(weights_sum)-1) = 1e0

  sky_ones = s2_sky_init(ones_map, s2_sky_get_nside(weights_sum), &
    s2_sky_get_pix_scheme(weights_sum))

  noise_disp = s2_sky_product(sky_ones, weights_sum, divide=.true.)

  call s2_sky_write_map_file(noise_disp, 'wnoise_coadded_disp.fits')

  call s2_sky_free(sky_ones)
  call s2_sky_free(noise_disp)
  deallocate(ones_map)

#endif 


  !----------------------------------------------------------------
  ! Create Gaussian CMB realisation and then compute co-added map.
  !----------------------------------------------------------------
  
  do isim = isim_start, isim_end

    write(*,'(a,a,i5)') msgpfx, 'Running simulation number ', isim

    ! Use system clock to set seed if not already specified.
!     if(.not. seed_set) then
!       call system_clock(seed_cmb)
!     end if

    ! Set seed to simulation number.
    seed_cmb = isim

    ! Construct simulated cmb map.
    ! Don't need to compute map yet since have to convolve with various beams.
    cmb = s2_cmb_init(cl, s2_sky_get_nside(nobs(1)), seed_cmb, &
      compute_map=.false.)

    ! Get cmb sky map.
    cmb_sky = s2_cmb_get_sky(cmb)

    ! Convert cmb from uK to mK to produce final output map in mK
    ! (so same units as cmb band maps).
    ! See http://cmbdata.gsfc.nasa.gov/
!    call s2_sky_scale(cmb_sky, 1.0e-3)

    ! Perform weighted sum over bands.
    do iband = 1, NUM_BANDS

      ! Convolve with beam for current band.
      cmb_sky_beam(iband) = s2_sky_init(cmb_sky)
      call s2_sky_conv(cmb_sky_beam(iband), beam(iband))

      ! Generate noise map for current band.
!       call system_clock(seed_noise)
      seed_noise = isim * iband  ! Set deterministic noise seed.
      wnoise = s2_wnoise_init(nobs(iband), SIGMA0(iband), seed_noise)
      wnoise_sky(iband) = s2_wnoise_get_sky(wnoise)
      call s2_wnoise_free(wnoise)

      ! Compute map from alms after convolved cmb with beam.
      call s2_sky_compute_map(cmb_sky_beam(iband))

    end do

    ! Combine bands to produce coadded map.
    cmb_coadded_unscaled = s2_gcmbcoad_weightedsum(cmb_sky_beam, weights)
    wnoise_coadded_unscaled = s2_gcmbcoad_weightedsum(wnoise_sky, weights)
    cmb_coadded = s2_sky_product(cmb_coadded_unscaled, weights_sum, divide=.true.)
    wnoise_coadded = s2_sky_product(wnoise_coadded_unscaled, weights_sum, divide=.true.)

    ! Apply mask.
    if(apply_mask) then
      sky_temp = s2_sky_init(cmb_coadded)
      cmb_coadded = s2_sky_product(sky_temp, mask)
      call s2_sky_free(sky_temp)
      sky_temp = s2_sky_init(wnoise_coadded)
      wnoise_coadded = s2_sky_product(sky_temp, mask)
      call s2_sky_free(sky_temp)
    end if

    ! Downsample to nside=256.
    call s2_sky_downsample(cmb_coadded, NSIDE_FINAL)
    call s2_sky_downsample(wnoise_coadded, NSIDE_FINAL)

    ! Save simulated map.
    write(filename_out,'(a,i6.6,a)') trim(filename_out_cmb_prefix), isim, '.fits'
    write(comment_out,'(a,i6.6)') 'iseed_cmb=', seed_cmb
    call s2_sky_write_map_file(cmb_coadded, filename_out, comment_out)
    write(filename_out,'(a,i6.6,a)') trim(filename_out_noise_prefix), isim, '.fits'
    call s2_sky_write_map_file(wnoise_coadded, filename_out)

    ! Free temporary memory used in loop.
    call s2_cmb_free(cmb)
    call s2_sky_free(cmb_sky)
    call s2_sky_free(cmb_coadded)
    call s2_sky_free(cmb_coadded_unscaled)
    call s2_sky_free(wnoise_coadded)
    call s2_sky_free(wnoise_coadded_unscaled)
    do iband = 1, NUM_BANDS
      call s2_sky_free(cmb_sky_beam(iband))
      call s2_sky_free(wnoise_sky(iband))
    end do
  
  end do

  
  !----------------------------------------------------------------
  ! Free memory.
  !----------------------------------------------------------------

  do iband = 1, NUM_BANDS
    call s2_sky_free(nobs(iband))
    call s2_pl_free(beam(iband))   
  end do
  call s2_pl_free(cl)
  call s2_sky_free(weights_sum)
  if(apply_mask) call s2_sky_free(mask)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - October 2005
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: n, i
      character(len=S2_STRING_LEN) :: opt
      character(len=S2_STRING_LEN) :: arg
      
      n = nArguments()
     
      do i = 1,n,2
        
        call getArgument(i,opt)
     
        if (i == n .and. trim(opt) /= '-help') then
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a,a)') 'Usage: s2_gcmbcoad ', &
              '[-isim_start isim_start]'
            write(*,'(a,a)') '                   ', &
              '[-isim_end isim_end (optional)]'
            write(*,'(a,a)') '                   ', &
              '[-dir_cl dir_cl]'
            write(*,'(a,a)') '                   ', &
              '[-dir_nobs dir_nobs (set to dir_cl if not specified)]'
            write(*,'(a,a)') '                   ', &            
              '[-dir_beam dir_beam (set to dir_cl if not specified)]'
            write(*,'(a,a)') '                   ', &
              '[-mask filename_mask (optional)]'
            write(*,'(a,a)') '                   ', &
              '[-lmax lmax]'
            write(*,'(a,a)') '                   ', &
              '[-out_cmb_prefix filename_out_cmb_prefix]'
            write(*,'(a,a)') '                   ', &
              '[-out_noise_prefix filename_out_noise_prefix]'              
            stop
          
          case ('-isim_start')
            read(arg,*) isim_start
            isim_end = isim_start    ! If just making one don't need
                                     ! to specify isim_end.

          case ('-isim_end')
            read(arg,*) isim_end

          case ('-dir_cl')
              dir_cl = trim(arg)
              dir_nobs = dir_cl
              dir_beam = dir_cl

          case ('-dir_nobs')
              dir_nobs = trim(arg)

          case ('-dir_beam')
              dir_beam = trim(arg)
          
          case ('-mask')
            filename_mask = trim(arg)
            apply_mask = .true.

          case ('-lmax')
            read(arg,*) lmax

          case ('-out_cmb_prefix')
            filename_out_cmb_prefix = trim(arg)

          case ('-out_noise_prefix')
            filename_out_noise_prefix = trim(arg)
            
          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options

end program s2_gcmbcoad


!--------------------------------------------------------------------------
! s2_gcmbcoad_weightedsum
!
!! Compute a weighted sum of sky maps.
!!
!! Notes:
!!   - Intent of skies and weights variables is inout since may need to
!!     change pixelisation scheme to be consistent in s2_sky_product.
!!
!! Variables:
!!   - skies: Array of skies to sum.
!!   - weights: Weight maps to apply to each corresponding sky in the
!!     summation.
!!   - weighted_sky_sum: Resultant weighted sum map.
!
!! @author J. D. McEwen
!! @version 0.1 October 2005
!
! Revisions:
!   October 2005 - Written by Jason McEwen
!--------------------------------------------------------------------------

function s2_gcmbcoad_weightedsum(skies, weights) result(weighted_sky_sum)

  use s2_sky_mod
  use s2_error_mod

  implicit none

  type(s2_sky), intent(inout) :: skies(1:)
  type(s2_sky), intent(inout) :: weights(1:)
  type(s2_sky) :: weighted_sky_sum

  integer :: iband, nbands
  type(s2_sky) :: weighted_sky_iband, weighted_sky_sum_temp
  
  ! Check sizes consistent.
  if(size(weights) /= size(skies)) then
    call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_gcmbcoad_weightedsum', &
      comment_add='Inconsistend number of sky and weight maps')
  end if 

  nbands = size(weights)

  do iband = 1,nbands

    weighted_sky_iband = s2_sky_product(weights(iband), skies(iband))

    ! Compute sum of weighted maps.
    if(iband == 1) then
      weighted_sky_sum = s2_sky_init(weighted_sky_iband)
    else
      weighted_sky_sum_temp = s2_sky_init(weighted_sky_sum)
      call s2_sky_free(weighted_sky_sum)
      weighted_sky_sum = s2_sky_add(weighted_sky_iband, weighted_sky_sum_temp)
      call s2_sky_free(weighted_sky_sum_temp)
    end if

    ! Free temporary memory used in loop.
    call s2_sky_free(weighted_sky_iband)

  end do
  
end function s2_gcmbcoad_weightedsum

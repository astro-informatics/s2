!------------------------------------------------------------------------------
! s2_gcmb
!
!! Generate a Gaussian CMB realisation from an input power spectrum.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_cl]: Name of cl input ascii file.
!!   - [-scale_cl scale_cl (logical)]: Logical to specify whether to scale
!!     the input cl values by 2*pi/(l(l+1)).
!!   - [-out filename_out]: Name of output cmb map file.
!!   - [-file_type_out file_type_out_str]: String specifying output file type.
!!   - [-nside nside]: Healpix nside resolution to generate map for.
!!   - [-lmax lmax]: Maximum hamonic l to consider in power spectrum when
!!     generating map.
!!   - [-lmin lmin]: Minimum l of power spectrum contained in input file
!!     (others are set to zero).  Often this is 2.
!!   - [-ncomment ncomment]: Number of comment lines in input file to ignore.
!!   - [-seed seed]: Integer seed for random number generator used to
!!     computed Gaussian realised CMB. If not specified then use system clock.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - July 2005
!
! Revisions:
!   July 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_gcmb

  use s2_types_mod
  use s2_cmb_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2_error_mod
  use s2_vect_mod

  implicit none

  type(s2_cmb) :: cmb
  integer :: nside, lmin, lmax, ncomment, seed, fail = 0
  character(S2_STRING_LEN) :: filename_cl
  character(S2_STRING_LEN) :: filename_out = 'out.fits'
  logical :: scale_cl = .true.
  logical :: seed_set = .false.
  real(s2_sp) :: noise_var = 0.0e0, beam_fwhm = 0.0e0
  real(s2_sp), allocatable :: npl(:)
  logical :: beam_present = .false.
  logical :: noise_present = .false.
  type(s2_pl) :: beam_pl, noise_pl, background_pl, background_pl_tmp
  type(s2_sky) :: cmb_sky

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type_out = S2_SKY_FILE_TYPE_MAP
  character(len=S2_STRING_LEN) :: file_type_out_str = MAP_FILE
  logical :: compute_map = .false.

  ! Set default parameter values.
  nside = 512
  lmax = 1024
  lmin = 2
  ncomment = 29   ! Value required to read WMAP cl files.
  seed = 1
  scale_cl = .true.
  
  ! Parse options.
  call parse_options()

  ! Set out sky file type.
  select case (trim(file_type_out_str))
    case (MAP_FILE)
       file_type_out = S2_SKY_FILE_TYPE_MAP
    case (ALM_FILE)
       file_type_out = S2_SKY_FILE_TYPE_ALM
    case (SKY_FILE)
       file_type_out = S2_SKY_FILE_TYPE_SKY
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_axiconv', &
         comment_add='Invalid file type option')
  end select

  ! Use system clock to set seed if not already specified.
  if(.not. seed_set) then
    call system_clock(seed)
  end if

  ! Read CMB power spectrum from file.
  background_pl = s2_pl_init(filename_cl, lmin, &
       lmax, ncomment, scale_cl)

  ! Apply beam. 
  if(beam_present) then

     ! Convert beam_fwhm to radians.
     beam_fwhm = s2_vect_arcmin_to_rad(beam_fwhm)

     ! Construct beam.
     beam_pl = s2_pl_init_guassian(beam_fwhm, lmax)

     ! Convolve background spectrum with beam.
     call s2_pl_conv(background_pl, beam_pl)

  end if

  ! Add noise.
  if(noise_present) then

     ! Generate noise pl vector.
     allocate(npl(0:lmax), stat=fail)
     if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_gcmb')
     end if
     npl(0:1) = 0e0
     npl(2:lmax) = noise_var             
     noise_pl = s2_pl_init(npl)
     deallocate(npl)

     ! Add noise spectrum.
     background_pl_tmp = s2_pl_add(background_pl, noise_pl)
     call s2_pl_free(background_pl)
     background_pl = s2_pl_init(background_pl_tmp)
     call s2_pl_free(background_pl_tmp)

  end if

  ! Simulate CMB.
  if (file_type_out /= S2_SKY_FILE_TYPE_ALM) compute_map = .true.
  cmb = s2_cmb_init(background_pl, nside, seed, compute_map=compute_map) 

  ! Write output file.
  cmb_sky = s2_cmb_get_sky(cmb)
  call s2_sky_write_file(cmb_sky, filename_out, file_type_out)
  call s2_sky_free(cmb_sky)

  ! Free memory.
  call s2_cmb_free(cmb)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - November 2004
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
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2_gcmb [-inp filename_cl]'
            write(*,'(a)') '               [-scale_cl scale_cl (logical)]'
            write(*,'(a)') '               [-out filename_out]'
            write(*,'(a,a)') '               ', &
                 '[-file_type_out file_type_out_str (sky; map; alm)]'
            write(*,'(a)') '               [-nside nside]'
            write(*,'(a)') '               [-lmax lmax]'
            write(*,'(a)') '               [-lmin lmin]'
            write(*,'(a)') '               [-beam_fwhm beam_fwhm]'
            write(*,'(a)') '               [-noise_var noise_var]'
            write(*,'(a)') '               [-ncomment ncomment]'
            write(*,'(a)') '               [-seed seed]'
            stop
          
          case ('-inp')
            filename_cl = trim(arg)

          case ('-scale_cl')
            read(arg,*) scale_cl
            
          case ('-out')
            filename_out = trim(arg)

          case ('-file_type_out')
            file_type_out_str = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax

          case ('-lmin')
            read(arg,*) lmin

          case ('-beam_fwhm')
            read(arg,*) beam_fwhm
            beam_present = .true.

          case ('-noise_var')
            read(arg,*) noise_var
            noise_present = .true.

          case ('-ncomment')
            read(arg,*) ncomment

          case ('-seed')
            read(arg,*) seed
            seed_set = .true.

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_gcmb

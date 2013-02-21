!------------------------------------------------------------------------------
! s2_almbeam
!
!! Read a HEALPix alm fits file and apply (optionally) a beam and
!! (optionally) a pixel window function.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_alm]: Name of HEAPix alm fits file to read.
!!   - [-out filename_out]: Name of output HEALPix alm file to write.
!!   - [-beam filename_beam (optional)]: Name of input beam fits file.
!!   - [-fwhm fwhm (arcmin) (optional)]: Full-width-half-maximum (in arcmin) 
!!     of Gaussian beam.
!!   - [-nside_pixel_window (optional)]: HEALPix nside of map if to be computed.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - February 2013
!
! Revisions:
!   February 2013 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_almbeam

  use s2_types_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2_vect_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_alm, filename_beam, filename_out
  type(s2_sky) :: sky
  type(s2_pl) :: beam, pixel_window
  logical :: gaussian_beam = .true.
  logical :: apply_beam = .false.
  logical :: apply_pixel_window = .false.
  integer :: nside_pixel_window
  real(s2_sp) :: fwhm, fwhm_arcmin = 13.2e0

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with alms read in from alm fits file.
  sky = s2_sky_init(filename_alm, S2_SKY_FILE_TYPE_ALM)

  ! Construct and apply beam.
  if (apply_beam) then

     ! Construct beam.
     if (gaussian_beam) then

        fwhm = s2_vect_arcmin_to_rad(fwhm_arcmin)
        beam = s2_pl_init_guassian(fwhm, s2_sky_get_lmax(sky))

     else

        beam = s2_pl_init(filename_beam)

     end if

     ! Apply beam.
     call s2_sky_conv(sky, beam)

     ! Free beam.
     call s2_pl_free(beam)

  end if

  ! Construct and apply pixel window function.
  if (apply_pixel_window) then

     ! Construct pixel window function.
     pixel_window =  &
          s2_pl_init_pixel_window(nside_pixel_window, &
                                  s2_sky_get_lmax(sky))

     ! Apply pixel window function.
     call s2_sky_conv(sky, pixel_window)

     ! Free pixel window function.
     call s2_pl_free(pixel_window)

  end if

  ! Write beamed sky to file.
  call s2_sky_write_alm_file(sky, filename_out)

  ! Free memory.
  call s2_sky_free(sky)



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
            write(*,'(a)') 'Usage: s2_almbeam [-inp filename_alm]'
            write(*,'(a)') '                  [-out filename_out]'
            write(*,'(a)') '                  [-beam filename_beam (optional)]'
            write(*,'(a)') '                  [-fwhm fwhm (in arcmin) (optional)]'
            write(*,'(a)') '                  [-nside_pixel_window nside (optional)]'
            stop

          case ('-inp')
            filename_alm = trim(arg)

          case ('-out')
            filename_out = trim(arg)
          
          case ('-beam')
            filename_beam = trim(arg)
            gaussian_beam = .false.
            apply_beam = .true.

          case ('-fwhm')
            read(arg,*) fwhm_arcmin
            gaussian_beam = .true.
            apply_beam = .true.

          case ('-nside_pixel_window')
            read(arg,*) nside_pixel_window
            apply_pixel_window = .true.

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program s2_almbeam

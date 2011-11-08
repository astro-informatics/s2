!------------------------------------------------------------------------------
! s2_alm2sky
!
!! Read a HEALPix alm fits file and write the read map to a full s2_sky fits
!! file.  Compute map if requested.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-alm filename_alm]: Name of map HEAPix alm fits file to read.
!!   - [-sky filename_sky]: Name of full s2_sky fits file to write.
!!   - [-compute_map compute_map]: Logical specifying whether to compute map.
!!   - [-nside]: HEALPix nside of map if to be computed.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2005
!
! Revisions:
!   November 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_alm2sky

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_sky, filename_alm
  type(s2_sky) :: sky
  logical :: compute_map = .false.
  integer :: nside = 64

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with alms read in from alm fits file.
  sky = s2_sky_init(filename_alm, S2_SKY_FILE_TYPE_ALM)

  ! Compute map if required.
  if(compute_map) then
     call s2_sky_compute_map(sky, nside)
  end if

  ! Write sky to s2_sky file.
  call s2_sky_io_fits_write(filename_sky, sky)

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
            write(*,'(a)') 'Usage: s2_alm2sky [-alm filename_alm]'
            write(*,'(a)') '                  [-sky filename_sky]'
            write(*,'(a)') '                  [-compute_map compute_map]'
            write(*,'(a)') '                  [-nside nside]'
            stop

          case ('-alm')
            filename_alm = trim(arg)
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-compute_map')
            read(arg,*) compute_map

          case ('-nside')
            read(arg,*) nside
            compute_map = .true.

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program s2_alm2sky

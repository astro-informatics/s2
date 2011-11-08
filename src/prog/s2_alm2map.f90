!------------------------------------------------------------------------------
! s2_alm2map
!
!! Read a HEALPix alm fits file, compute map and then write the map to a 
!! HEALPix map fits file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-alm filename_alm]: Name of map HEAPix alm fits file to read.
!!   - [-map filename_map]: Name of HEALPix map fits file to write.
!!   - [-nside]: HEALPix nside of map to be computed.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2005
!
! Revisions:
!   November 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_alm2map

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_map, filename_alm
  type(s2_sky) :: sky
  integer :: nside = 64

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with alms read in from alm fits file.
  sky = s2_sky_init(filename_alm, S2_SKY_FILE_TYPE_ALM)

  ! Compute map.
  call s2_sky_compute_map(sky, nside)

  ! Write sky to Healpix map fits file.
  call s2_sky_write_map_file(sky, filename_map)

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
            write(*,'(a)') 'Usage: s2_alm2map [-alm filename_alm]'
            write(*,'(a)') '                  [-map filename_map]'
            write(*,'(a)') '                  [-nside nside]'
            stop

          case ('-alm')
            filename_alm = trim(arg)
          
          case ('-map')
            filename_map = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case default
            print '("Unknown option ",a," ignored")', trim(opt)
        end select
      end do

    end subroutine parse_options


end program s2_alm2map

!------------------------------------------------------------------------------
! s2_skyup
!
!! Up-sample a sky.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of input file containing sky to up-sample.
!!   - [-out filename_out]: Name of output file containing up-sampled sky.
!!   - [-nside nside_new]: New nside of up-sampled map.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [-ext ext (optional)]: File extension for map files.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - August 2005
!
! Revisions:
!   August 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skyup

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  integer :: nside_new
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyup', &
         comment_add='Invalid file type option')

  end select

  ! Initialse sky with map read in from fits file.
  sky = s2_sky_init(filename_in, file_type, ext)

  ! Up-sample the map.
  call s2_sky_upsample(sky, nside_new)
  
  ! Save up-sampled sky map.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(sky, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, sky)
       
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyup', &
            comment_add='Invalid file type specifier')
       
  end select

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
            write(*,'(a)') 'Usage: s2_skyup [-inp filename_in]'
            write(*,'(a)') '                 [-out filename_out]'
            write(*,'(a)') '                 [-nside nside_new]'
            write(*,'(a)') '                 [-file_type file_type_str]'
            write(*,'(a)') '                 [-ext ext (optional)]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-nside')
            read(arg,*) nside_new

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_skyup





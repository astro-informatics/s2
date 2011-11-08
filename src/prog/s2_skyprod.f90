!------------------------------------------------------------------------------
! s2_skyprod
!
!! Take the product of two maps of the sky.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-sky1 filename1]: Name of file containing first sky.
!!   - [-sky2 filename2]: Name of file containing second sky.
!!   - [-out filename_out]: Name of output sky file.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [-ext ext (optional)]: File extension for map files.
!!   - [div divide (logical)]: Logical specifying whether to divide rather 
!!     than multiply maps (if so result is a/b).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - August 2005
!
! Revisions:
!   August 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skyprod

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: filename1, filename2, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  logical :: divide = .false.
  type(s2_sky) :: sky1, sky2, sky_new

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyprod', &
         comment_add='Invalid file type option')

  end select

  ! Initialse sky with map read in from fits file.
  sky1 = s2_sky_init(filename1, file_type, ext)
  sky2 = s2_sky_init(filename2, file_type, ext)

  ! Create product map.
  sky_new = s2_sky_product(sky1, sky2, divide)

  ! Save new sky map.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(sky_new, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, sky_new)
       
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyprod', &
            comment_add='Invalid file type specifier')
       
  end select

  ! Free memory.
  call s2_sky_free(sky1)
  call s2_sky_free(sky2)
  call s2_sky_free(sky_new)


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
            write(*,'(a)') 'Usage: s2_skyprod [-sky1 filename1]'
            write(*,'(a)') '                   [-sky2 filename2]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-file_type file_type_str]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            write(*,'(a)') '                   [-div divide (logical)]'
            stop
          
          case ('-sky1')
            filename1 = trim(arg)

          case ('-sky2')
            filename2 = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case ('-div')
            read(arg,*) divide

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_skyprod





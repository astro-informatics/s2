!------------------------------------------------------------------------------
! s2_skythres
!
!! Threshold sky.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-sky filename]: Name of file containing sky to threshold.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [-ext ext (optional)]: File extension for map files.
!!   - [-lower thres_lower]: Lower threshold value.
!!   - [-upper thres_upper]: Upper threshold value.
!!   - [-abs_thres abs_thres (logical)]: Logical to specify whether to 
!!     threhold based on absolute values (in which case lower threshold is 
!!     used only).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2007
!
! Revisions:
!   April 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_skythres

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: filename, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  type(s2_sky) :: sky
  real(s2_sp) :: thres_lower=0e0, thres_upper=5e0
  logical :: abs_thres = .false.

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skythres', &
         comment_add='Invalid file type option')

  end select

  ! Initialse sky with map read in from fits file.
  sky = s2_sky_init(filename, file_type, ext)

  ! Perform thresholding.
  if(abs_thres) then
     call s2_sky_thres_abs(sky, thres_lower)
  else
     call s2_sky_thres(sky, thres_lower, thres_upper)
  end if

  ! Write output.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(sky, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, sky)

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skythres', &
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
            write(*,'(a)') 'Usage: s2_skythres [-sky filename]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-file_type file_type_str]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            write(*,'(a)') '                   [-lower thres_lower]'
            write(*,'(a)') '                   [-upper thres_upper]'
            write(*,'(a)') '                   [-abs_thres abs_thres (logical)]'
            stop
          
          case ('-sky')
            filename = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case ('-out')
            filename_out = trim(arg)

          case ('-lower')
            read(arg,*) thres_lower

          case ('-upper')
            read(arg,*) thres_upper

          case ('-abs_thres')
            read(arg,*) abs_thres

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_skythres





!------------------------------------------------------------------------------
! s2_almrot
!
!! Rotate a sky in harmonic space.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of file containing input sky alms.
!!   - [-out filename_out]: Name of output file for rotated sky alms.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [alpha alpha]: Alpha Euler angle of rotation (in degrees).
!!   - [beta beta]: Beta Euler angle of rotation (in degrees).
!!   - [gamma gamma]: Gamma Euler angle of rotation (in degrees).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version Under svn version control.
!
! Revisions:
!   June 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_almrot

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_ALM
  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = ALM_FILE
  real(s2_dp) :: alpha=0d0, beta=0d0, gamma=0d0
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (ALM_FILE)
       file_type = S2_SKY_FILE_TYPE_ALM

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_almrot', &
         comment_add='Invalid file type option')

  end select

  ! Initialse sky with map read in from fits file.
  sky = s2_sky_init(filename_in, file_type)

  ! Rotate sky.
  call s2_sky_rotate_alm(sky, alpha*pi/180d0, beta*pi/180d0, gamma*pi/180d0)

  ! Save new sky map.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_ALM)
       call s2_sky_write_alm_file(sky, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, sky)
       
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_almrot', &
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
            write(*,'(a)') 'Usage: s2_almrot [-inp filename_in]'
            write(*,'(a)') '                 [-out filename_out]'
            write(*,'(a)') '                 [-file_type file_type_str (alm; sky)]'
            write(*,'(a)') '                 [-alpha alpha (degrees)]'
            write(*,'(a)') '                 [-beta beta (degrees)]'
            write(*,'(a)') '                 [-gamma gamma (degrees)]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-alpha')
            read(arg,*) alpha

          case ('-beta')
            read(arg,*) beta

          case ('-gamma')
            read(arg,*) gamma

          case default
            print '("Unknown option ",a," ignored")', trim(opt) 

        end select
      end do

    end subroutine parse_options


end program s2_almrot





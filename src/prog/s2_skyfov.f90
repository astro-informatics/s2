!------------------------------------------------------------------------------
! s2_sky2proj
!
!!



 
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of file containing input sky .
!!   - [-file_type file_type_str]: String specifying input file type.
!!   - [-out filename_out]: Name of output file for projected image.
!!   - [method method]: Projection method to use.
!!   - [theta_fov theta_fov]: Field of view (in degrees).
!!   - [nside nside]: Healpix nside to use when projecting image (required by 
!!     some projection methods if not provided in sky object).
!!   - [lmax lmax]: Harmonic band-limit to use when projecting image
!!     (required by some projection methods if not provided in sky object).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skyfov

  use s2_types_mod
  use s2_sky_mod
  use s2_proj_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_ALM

  character(len=*), parameter ::  METHOD_CIRCLE = 'circle'
  character(len=*), parameter ::  METHOD_SQUARE = 'square'
  character(len=S2_STRING_LEN) :: method_str = METHOD_CIRCLE
  integer :: method = S2_SKY_FOV_METHOD_CIRCLE

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE

  real(s2_dp) :: theta_fov = 90.0
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()
  theta_fov = theta_fov * pi / 180e0

  ! Set sky file type.
  select case (trim(file_type_str))
    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (ALM_FILE)
       file_type = S2_SKY_FILE_TYPE_ALM

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyfov', &
         comment_add='Invalid file type option')
  end select

  ! Set method type.
  select case (trim(method_str))
    case (METHOD_CIRCLE)
       method = S2_SKY_FOV_METHOD_CIRCLE

    case (METHOD_SQUARE)
       method = S2_SKY_FOV_METHOD_SQUARE

    case default
       call s2_error(S2_ERROR_SKY_FOV_METHOD_INVALID, 's2_skyfov', &
         comment_add='Invalid fov method type option')
  end select

  ! Initialse sky from file.
  sky = s2_sky_init(filename_in, file_type)

  ! Restrict sky map to fov.
  call s2_sky_fov(sky, theta_fov, method)

  ! Write resultant sky to file.
  call s2_sky_write_file(sky, trim(filename_out), file_type)

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
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2_sky2proj [-inp filename_in]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-file_type file_type_str (sky; map; alm)]'
            write(*,'(a)') '                   [-method method (circle; square)]'
            write(*,'(a)') '                   [-theta_fov theta_fov (in degrees)]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-method')
            method_str = trim(arg)

          case ('-theta_fov')
            read(arg,*) theta_fov

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_skyfov





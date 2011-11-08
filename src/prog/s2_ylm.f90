!------------------------------------------------------------------------------
! s2_ylm
!
!! Compute spherical harmonic function on full sky.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-el el]: Harmonic scale l.
!!   - [-m m]: Harmonic azimuthal scale m.
!!   - [-nside nside]:  Healpix resolution parameter.
!!   - [-out filename_out]: Name of output file containing Ylm sky.
!!   - [-method method (leg or wig; optional)]: Method used to compute Ylm.
!!   - [-reality reality (abs, real or imag; optional)]: Reality type to 
!!     specify whether to compute the real part, imaginary part or absolute
!!     value of Ylm.
!!   - [-file_type file_type_str (optional)]: String specifying file type.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - February 2008
!
! Revisions:
!   February 2008 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_ylm

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_ylm_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  character(len=*), parameter ::  METHOD_LEG = 'leg'
  character(len=*), parameter ::  METHOD_WIG = 'wig'
  character(len=S2_STRING_LEN) :: method_str = METHOD_LEG
  integer :: method = S2_YLM_EVAL_METHOD_LEG
  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=*), parameter ::  REALITY_REAL = 'real'
  character(len=*), parameter ::  REALITY_IMAG = 'imag'
  character(len=*), parameter ::  REALITY_ABS = 'abs'
  character(len=S2_STRING_LEN) :: reality_str = REALITY_ABS
  integer :: reality = S2_YLM_REALITY_ABS
  integer :: nside, el, m
  type(s2_sky) :: ylm

  ! Set default parameters.
  nside = 128
  el = 1
  m = 1
  filename_out = 'ylm.fits'

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_ylm', &
         comment_add='Invalid file type option')

  end select

  ! Set method used to compute Ylm.
  select case (trim(method_str))

    case (METHOD_LEG)
       method = S2_YLM_EVAL_METHOD_LEG

    case (METHOD_WIG)
       method = S2_YLM_EVAL_METHOD_WIG

    case default
       call s2_error(S2_ERROR_YLM_ARG_INVALID, 's2_ylm', &
         comment_add='Invalid method')

  end select

  ! Set reality parameter.
  select case (trim(reality_str))

    case (REALITY_REAL)
       reality = S2_YLM_REALITY_REAL

    case (REALITY_IMAG)
       reality = S2_YLM_REALITY_IMAG

    case (REALITY_ABS)
       reality = S2_YLM_REALITY_ABS

    case default
       call s2_error(S2_ERROR_YLM_ARG_INVALID, 's2_ylm', &
         comment_add='Invalid reality parameter')

  end select
  
  ! Compute spherical harmonic function.
  ylm = s2_ylm_sky(el, m, nside, reality, S2_SKY_RING, method)

  ! Write file.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(ylm, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, ylm)
       
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_ylm', &
            comment_add='Invalid file type specifier')
       
  end select

  ! Free memory.
  call s2_sky_free(ylm)


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
            write(*,'(a)') 'Usage: s2_ylm [-el el]'
            write(*,'(a)') '              [-m m]'
            write(*,'(a)') '              [-nside nside]'
            write(*,'(a)') '              [-out filename_out]'
            write(*,'(a)') '              [-method method (leg or wig; optional)]'
            write(*,'(a)') '              [-reality reality (abs, real or imag; optional)]'
            write(*,'(a)') '              [-file_type file_type_str (optional)]'
            stop
          
          case ('-el')
            read(arg,*) el

          case ('-m')
            read(arg,*) m

          case ('-nside')
            read(arg,*) nside

          case ('-out')
            filename_out = trim(arg)

          case ('-method')
            method_str = trim(arg)

          case ('-reality')
            reality_str = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_ylm





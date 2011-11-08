!------------------------------------------------------------------------------
! s2_sky2fsht
!
!! Read a HEALPix pixelised sky file (map or sky format) and extra an 
!! ecp (equispaced) sampled theta-phi array over the sphere for the grid 
!! used for the Fast Spherical Harmonic Transform.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-sky filename_sky]: Name of full file to read.
!!   - [-L L]: Harmonic band-limit L for FSHT pixelisation (related to number 
!!     of samples on the sphere for the FSHT ECP pixelisation).
!!   - [-file_type file_type_str]: String to specify file type (either map 
!!     or sky.
!!   - [-ext ext (optional)]: Extension number of Healpix fits file.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2008
!
! Revisions:
!   April 2008 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_sky2fsht

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1, fail = 0
  character(len=S2_STRING_LEN) :: filename_sky
  type(s2_sky) :: sky
  real(s2_sp), allocatable :: xtp(:,:)
  integer :: L = 50

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky2fsht', &
         comment_add='Invalid file type option')

  end select

  ! Allocate memory.
  allocate(xtp(0:L,0:2*L), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_sky2fsht')
  end if

  ! Read full sky object.
  sky = s2_sky_init(filename_sky, file_type)

  ! Extract theta-phi grid for FSHT.
  call s2_sky_extract_ab_fsht(sky, xtp(0:L,0:2*L), L)


  ! Currently xtp writting to standard out when extracted but could save here
  ! if required.


  ! Free memory.
  call s2_sky_free(sky)
  deallocate(xtp)


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
            write(*,'(a)') 'Usage: s2_sky2fsht [-sky filename_sky]'
            write(*,'(a)') '                   [-L L]'
            write(*,'(a)') '                   [-file_type file_type_str]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            stop
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-L')
            read(arg,*) L

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_sky2fsht

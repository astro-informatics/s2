!------------------------------------------------------------------------------
! s2_about
!
!! Display information about the S2 package.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2011
!
! Revisions:
!   November 2011 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_about

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  character(len=*), parameter ::  RING_STR = 'ring'
  character(len=*), parameter ::  NEST_STR = 'nest'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  character(len=S2_STRING_LEN) :: order_str
  integer :: order
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()

  ! Display info.
  write(*,'(a)') "=========================================================="
  write(*,'(a)') "S2 package to analyse functions on the sphere"
  write(*,'(a)') "By Jason McEwen"

  write(*,'(a)') "See www.jasonmcewen.org for more information."
  write(*,'(a)') "See LICENSE.txt for license details."

  write(*,'(a,a)') "Version: ", S2_VERSION
  write(*,'(a,a)') "Build: ", S2_BUILD
  write(*,'(a)') "=========================================================="


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
            write(*,'(a)') 'Usage: s2_about'
            stop

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_about





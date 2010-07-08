!------------------------------------------------------------------------------
! s2_skyconvsp
!
!! Compute the convolution in real space of a sky with a kernel.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of file containing input sky.
!!   - [-out filename_out]: Name of output file for convolved sky.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [-sigma sigma]: Sigma of the Gaussian kernel (only a Gaussian kernel 
!!     is supported at present).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skyconvsp

  use s2_types_mod
  use s2_sky_mod
  use s2_proj_mod
  use s2_error_mod

  implicit none

  interface 
     function kernel(theta, param) result(val)
       use s2_types_mod
       real(s2_dp), intent(in) :: theta
       real(s2_dp), intent(in), optional :: param(:)
       real(s2_dp) :: val
     end function kernel
  end interface

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE

  real(s2_dp) :: sigma(1), support_theta
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))
    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (ALM_FILE)
       file_type = S2_SKY_FILE_TYPE_ALM
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyconvsp', &
         comment_add='Invalid file type option')

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyconvsp', &
         comment_add='Invalid file type option')
  end select

  ! Initialse sky with map or alms read from file.
  sky = s2_sky_init(filename_in, file_type)

  ! Perform convolution.
!  sigma(1) = 0.03
  support_theta = 2.0 * sigma(1)
  call s2_sky_conv_space(sky, support_theta, kernel, sigma)

  ! Save output file.
  call s2_sky_write_file(sky, filename_out, file_type)

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
            write(*,'(a)') 'Usage: s2_skyconvsp [-inp filename_in]'
            write(*,'(a)') '                    [-out filename_out]'
            write(*,'(a)') '                    [-file_type file_type_str (sky; map)]'
            write(*,'(a)') '                    [-sigma sigma]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-sigma')
            read(arg,*) sigma(1)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_skyconvsp


!---------------------------------------------------------------------
! kernel
!
!! Gaussian convolution kernel function
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   July 2010 - Written by Jason McEwen 
!---------------------------------------------------------------------

function kernel(theta, param) result(val)

  use s2_types_mod
  real(s2_dp), intent(in) :: theta
  real(s2_dp), intent(in), optional :: param(:)
  real(s2_dp) :: val

  real(s2_dp) :: sigma

  sigma = param(1)
  val = exp(-theta**2 / (2.0*sigma**2))

end function kernel

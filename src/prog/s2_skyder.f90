!------------------------------------------------------------------------------
! s2_skyder
!
!! Compute continuous or discrete derivatives on the sphere.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of file containing input sky.
!!   - [-out filename_out]: Name of output file for convolved sky.
!!   - [-file_type file_type_str]: String specifying file types.
!!   - [-der_type der_type_str (theta; phi; grad)]:  Derivative to compute.
!!   - [-discrete discrete (true; false)]: Whether to compute a discrete or
!!     continuous derivative.
!!   - [-apply_sin apply_sin (true; false)]: Whether to apply 1/sin(theta) 
!!     when computing dT/dphi.
!!   - [-lmax lmax]: Harmonic band limit considered for continuous 
!!     derivatives.
!!   - [-mmax mmax]: Azimuthal harmonic band limit considered for continuous 
!!     derivatives.
!!   - [-sigma sigma]: Sigma of the Gaussian kernel (only a Gaussian kernel 
!!     is supported at present).
!!   - [-op_file filename_op]: Filename of operator to save (optional).
!!   - [-theta_fov theta_fov (in degrees)]: Field-of-view operator defined 
!!     on.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   August 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skyder

  use s2_types_mod
  use s2_sky_mod
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

  character(len=*), parameter ::  DER_THETA = 'theta'
  character(len=*), parameter ::  DER_PHI = 'phi'
  character(len=*), parameter ::  DER_GRAD = 'grad'
  character(len=S2_STRING_LEN) :: der_type_str = DER_GRAD
  integer :: der_type = S2_SKY_DER_TYPE_GRAD

  type(s2_sky) :: sky
  type(s2_sky) :: der
  integer :: lmax = 128, mmax = 128
  logical :: discrete = .false.
  logical :: apply_sin = .false.
  real(s2_dp) :: sigma(1), support_theta

  logical :: save_op = .false.
  real(s2_dp) :: theta_fov
  character(len=S2_STRING_LEN) :: filename_op
  integer :: nsphere, nop, j, fileid
  real(s2_dp), allocatable :: op(:,:)
  real(s2_sp), allocatable :: xmap(:)

  ! Parse input parameters.
  sigma(1) = 0.02
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
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skyder', &
         comment_add='Invalid file type option')
  end select

  ! Set der type.
  select case (trim(der_type_str))
    case (DER_THETA)
       der_type = S2_SKY_DER_TYPE_THETA

    case (DER_PHI)
       der_type = S2_SKY_DER_TYPE_PHI

    case (DER_GRAD)
       der_type = S2_SKY_DER_TYPE_GRAD

    case default
       call s2_error(S2_ERROR_SKY_DER_TYPE_INVALID, 's2_skyder')
  end select

  ! Initialse sky with map or alms read from file.
  sky = s2_sky_init(filename_in, file_type)

  if(discrete) then

     select case(der_type)
          
          case(S2_SKY_DER_TYPE_THETA)
             support_theta = 4.0 * sigma(1)
             der = s2_sky_der_discrete_theta(sky, support_theta, &
                  kernel, sigma, nearest=.false.)

          case(S2_SKY_DER_TYPE_PHI)
             der = s2_sky_der_discrete_phi(sky, apply_sin)

          case(S2_SKY_DER_TYPE_GRAD)
             support_theta = 4.0 * sigma(1)
             der = s2_sky_der_discrete_grad(sky, apply_sin, &
                  support_theta, kernel, sigma, nearest=.false.)

          case default
             call s2_error(S2_ERROR_SKY_DER_TYPE_INVALID, 's2_skyder')

      end select

  else
     
     ! Compute alms.
     if(.not. s2_sky_get_alm_status(sky)) then
        call s2_sky_compute_alm(sky, lmax, mmax)
     end if

     ! Compute derivative.
     der = s2_sky_der(sky, der_type)

  end if

  ! Save output file.
  call s2_sky_write_file(der, filename_out, S2_SKY_FILE_TYPE_MAP)

  ! Write the derivative operator to file.
  if (save_op) then

     ! Compute operator.
     select case(der_type)
          
         case(S2_SKY_DER_TYPE_THETA)            
            support_theta = 4.0 * sigma(1)
            call s2_sky_der_discrete_theta_fovop(sky, &
                 support_theta, kernel, sigma, nearest=.false., &
                 theta_fov=theta_fov, nop=nop, op=op, nsphere=nsphere, xmap=xmap)

          case(S2_SKY_DER_TYPE_PHI)
             call s2_sky_der_discrete_phi_fovop(sky, apply_sin, theta_fov, nop, op, &
                  nsphere, xmap)

          case default
             call s2_error(S2_ERROR_SKY_DER_TYPE_INVALID, 's2_skyder', &
                  comment_add='Gradient operator not supported at present.')

      end select

     ! Write to file.
     fileid = 51
     open(unit=fileid, file=trim(filename_op), status='new', action='write', &
          form='formatted')
     do j = 0,nop-1
        write(fileid,'(2i20,e20.10)') nint(op(j,0)), &
             nint(op(j,1)), op(j,2)
     end do
     close(fileid)

     ! Free memory.
     deallocate(xmap)
     deallocate(op)

  end if

  ! Free memory.
  call s2_sky_free(sky)
  call s2_sky_free(der)

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
            write(*,'(a)') 'Usage: s2_skyder [-inp filename_in]'
            write(*,'(a)') '                 [-file_type file_type_str (sky; map; alm)]'
            write(*,'(a)') '                 [-out filename_out]'
            write(*,'(a)') '                 [-der_type der_type_str (theta; phi; grad)]'
            write(*,'(a)') '                 [-discrete discrete (true; false)]'
            write(*,'(a)') '                 [-apply_sin apply_sin (true; false)]'
            write(*,'(a)') '                 [-lmax lmax]'
            write(*,'(a)') '                 [-mmax mmax]'
            write(*,'(a)') '                 [-sigma sigma]'
            write(*,'(a)') '                 [-op_file filename_op]'
            write(*,'(a)') '                 [-theta_fov theta_fov (in degrees)]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-der_type')
            der_type_str = trim(arg)

          case ('-discrete')
            read(arg,*) discrete

          case ('-apply_sin')
            read(arg,*) apply_sin

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case ('-sigma')
            read(arg,*) sigma(1)

          case ('-op_file')
            filename_op = trim(arg)
            save_op = .true.

          case ('-theta_fov')
            read(arg,*) theta_fov

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_skyder


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

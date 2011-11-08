!------------------------------------------------------------------------------
! s2_sky2proj
!
!! Project a sky onto a planar image.
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
!!   - [-op_file filename_op (optional)]: File to save matrix representation 
!!     of projection operator.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_sky2proj

  use s2_types_mod
  use s2_sky_mod
  use s2_proj_mod
  use s2_error_mod

  implicit none

  interface 
     function conv_kernel(theta, param) result(val)
       use s2_types_mod
       real(s2_dp), intent(in) :: theta
       real(s2_dp), intent(in), optional :: param(:)
       real(s2_dp) :: val
     end function conv_kernel
  end interface

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_ALM

  character(len=*), parameter ::  METHOD_NEAREST = 'nearest'
  character(len=*), parameter ::  METHOD_HARMONIC = 'harmonic'
  character(len=*), parameter ::  METHOD_KERNEL = 'kernel'
  integer :: method = S2_PROJ_METHOD_NEAREST_NEIGHBOUR

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = ALM_FILE
  character(len=S2_STRING_LEN) :: method_str = METHOD_NEAREST

  real(s2_sp) :: theta_fov
  integer :: nside = 0, lmax = 0
  type(s2_sky) :: sky
  type(s2_proj) :: proj
  real(s2_dp) :: sigma_conv = 0.02, sigma(1), support_theta

  logical :: save_op = .false.
  logical :: save_convop = .false.
  logical :: save_xmap = .false.
  logical :: save_ang = .false.
  character(len=S2_STRING_LEN) :: filename_op, filename_convop, filename_xmap
  character(len=S2_STRING_LEN) :: filename_ang
  integer :: nsphere, nop, j, fileid
  integer, allocatable :: op(:,:)
  real(s2_dp), allocatable :: op_dp(:,:)
  real(s2_sp), allocatable :: xmap(:)
  real(s2_dp), allocatable :: thetas(:), phis(:)

  ! Parse input parameters.
  call parse_options()
  theta_fov = theta_fov * pi / 180e0
  sigma(1) = sigma_conv

  ! Set sky file type.
  select case (trim(file_type_str))
    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (ALM_FILE)
       file_type = S2_SKY_FILE_TYPE_ALM

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky2proj', &
         comment_add='Invalid file type option')
  end select

  ! Set method type.
  select case (trim(method_str))
    case (METHOD_NEAREST)
       method = S2_PROJ_METHOD_NEAREST_NEIGHBOUR

    case (METHOD_HARMONIC)
       method = S2_PROJ_METHOD_HARMONIC_INTERP

    case (METHOD_KERNEL)
       method = S2_PROJ_METHOD_KERNEL

    case default
       call s2_error(S2_ERROR_PROJ_METHOD_INVALID, 's2_sky2proj', &
         comment_add='Invalid projection type option')
  end select

  ! Initialse sky with map or alms read from file.
  sky = s2_sky_init(filename_in, file_type)

  ! Project sky onto image plane.
  if (nside == 0 .and. lmax == 0) then
     proj = s2_proj_init(sky, theta_fov, method)
  elseif (lmax == 0) then
     proj = s2_proj_init(sky, theta_fov, method, nside=nside)
  elseif (nside == 0) then
     proj = s2_proj_init(sky, theta_fov, method, lmax=lmax)
  else
     proj = s2_proj_init(sky, theta_fov, method, nside=nside, lmax=lmax)
  end if

  ! Write image to file.
  call s2_proj_write_image_file(proj, trim(filename_out))

  ! Write the projection operator to file.
  if (save_op .or. save_xmap .or. save_ang) then

      select case(method)

         case(S2_PROJ_METHOD_NEAREST_NEIGHBOUR)

            if (nside == 0) then
               call s2_proj_operator_nearest_neighbour(proj, nop=nop, op=op, &
                    nsphere=nsphere, xmap=xmap)
            else
               call s2_proj_operator_nearest_neighbour(proj, nside, nop, op, nsphere, xmap)
            end if

            if(save_op) then

               fileid = 43
               open(unit=fileid, file=trim(filename_op), status='new', action='write', &
                    form='formatted')

               do j = 0,nop-1
                  write(fileid,'(2i20)') op(j,0), op(j,1)
               end do

               close(fileid)

            end if

            if(save_xmap) then

               fileid = 44
               open(unit=fileid, file=trim(filename_xmap), status='new', action='write', &
                    form='formatted')

               do j = 0,nsphere-1
                  write(fileid,'(e20.10)') xmap(j)
               end do

               close(fileid)

            end if

            if(save_ang) then

               fileid = 45
               open(unit=fileid, file=trim(filename_ang), status='new', action='write', &
                    form='formatted')

               if(allocated(xmap)) deallocate(xmap)
               call s2_proj_get_xmap(proj, nside, nsphere, xmap, thetas, phis)
               do j = 0,nsphere-1
                  write(fileid,'(e20.10,e20.10)') thetas(j), phis(j)
               end do

               close(fileid)
               deallocate(thetas, phis)

            end if

            ! Free memory.
            deallocate(xmap, op)

         case(S2_PROJ_METHOD_KERNEL)

            if (nside == 0) then
               call s2_proj_operator_kernel(proj, nop=nop, op=op_dp, &
                    nsphere=nsphere, xmap=xmap)
            else
               call s2_proj_operator_kernel(proj, nside, nop, op_dp, nsphere, xmap)
            end if

            if(save_op) then

               fileid = 43
               open(unit=fileid, file=trim(filename_op), status='new', action='write', &
                    form='formatted')

               do j = 0,nop-1
                  write(fileid,'(2i20,e20.10)') nint(op_dp(j,0)), &
                       nint(op_dp(j,1)), op_dp(j,2)
               end do

               close(fileid)

            end if

            if(save_xmap) then

               fileid = 44
               open(unit=fileid, file=trim(filename_xmap), status='new', action='write', &
                    form='formatted')

               do j = 0,nsphere-1
                  write(fileid,'(e20.10)') xmap(j)
               end do

               close(fileid)

            end if

            if(save_ang) then

               fileid = 45
               open(unit=fileid, file=trim(filename_ang), status='new', action='write', &
                    form='formatted')

               if(allocated(xmap)) deallocate(xmap)
               call s2_proj_get_xmap(proj, nside, nsphere, xmap, thetas, phis)
               do j = 0,nsphere-1
                  write(fileid,'(e20.10,e20.10)') thetas(j), phis(j)
               end do

               close(fileid)
               deallocate(thetas, phis)

            end if

            ! Free memory.
            deallocate(xmap, op_dp)

         case default
            call s2_error(S2_ERROR_PROJ_METHOD_INVALID, 's2_sky2proj', &
                 comment_add='Output projection operator only defined for nearst neighbour and kernel projection.')

      end select

  end if

  ! Write the convolution operator to file.
  if (save_convop) then

     support_theta = 4 * sigma(1)
     call s2_sky_conv_space_fovop(sky, support_theta, real(theta_fov,s2_dp), nop, op_dp, &
         nsphere, xmap, conv_kernel, sigma)

     fileid = 43
     open(unit=fileid, file=trim(filename_convop), status='new', action='write', &
          form='formatted')
     do j = 0,nop-1
        write(fileid,'(2i20,e20.10)') nint(op_dp(j,0)), &
             nint(op_dp(j,1)), op_dp(j,2)
     end do
     close(fileid)

     ! Free memory.
     deallocate(xmap, op_dp)

  end if

  ! Free memory.
  call s2_sky_free(sky)
  call s2_proj_free(proj)


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
            write(*,'(a)') 'Usage: s2_sky2proj [-inp filename_in]'
            write(*,'(a)') '                   [-file_type file_type_str (sky; map; alm)]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-method method (nearest; harmonic)]'
            write(*,'(a)') '                   [-theta_fov theta_fov (in degrees)]'
            write(*,'(a)') '                   [-nside nside (optional)]'
            write(*,'(a)') '                   [-lmax lmax (optional)]'
            write(*,'(a)') '                   [-op_file filename_op (optional)]'
            write(*,'(a)') '                   [-convop_file filenameconv_op (optional)]'
            write(*,'(a)') '                   [-sigma_conv sigma_conv (optional)]'
            write(*,'(a)') '                   [-xmap_file filename_xmap (optional)]'
            write(*,'(a)') '                   [-ang_file filename_ang (optional)]'
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

          case ('-nside')
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax

          case ('-op_file')
            filename_op = trim(arg)
            save_op = .true.

          case ('-convop_file')
            filename_convop = trim(arg)
            save_convop = .true.

          case ('-sigma_conv')
            read(arg,*) sigma_conv

          case ('-xmap_file')
            filename_xmap = trim(arg)
            save_xmap = .true.

          case ('-ang_file')
            filename_ang = trim(arg)
            save_ang = .true.

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_sky2proj




function conv_kernel(theta, param) result(val)

  use s2_types_mod
  real(s2_dp), intent(in) :: theta
  real(s2_dp), intent(in), optional :: param(:)
  real(s2_dp) :: val

  real(s2_dp) :: sigma

  sigma = param(1)
  val = exp(-theta**2 / (2.0*sigma**2))

end function conv_kernel


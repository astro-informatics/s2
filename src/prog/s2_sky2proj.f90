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

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_ALM

  character(len=*), parameter ::  METHOD_NEAREST = 'nearest'
  character(len=*), parameter ::  METHOD_HARMONIC = 'harmonic'
  integer :: method = S2_PROJ_METHOD_NEAREST_NEIGHBOUR

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = ALM_FILE
  character(len=S2_STRING_LEN) :: method_str = METHOD_NEAREST

  real(s2_sp) :: theta_fov
  integer :: nside = 0, lmax = 0
  type(s2_sky) :: sky
  type(s2_proj) :: proj

  logical :: save_op = .false.
  character(len=S2_STRING_LEN) :: filename_op
  integer :: nsphere, nop, j, fileid
  integer, allocatable :: op(:,:)
  real(s2_sp), allocatable :: xmap(:)

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
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky2proj', &
         comment_add='Invalid file type option')
  end select

  ! Set method type.
  select case (trim(method_str))
    case (METHOD_NEAREST)
       method = S2_PROJ_METHOD_NEAREST_NEIGHBOUR

    case (METHOD_HARMONIC)
       method = S2_PROJ_METHOD_HARMONIC_INTERP

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
  if (save_op) then

     ! Output projection operator file only defined for nearest 
     ! neighbour projection method at present.
     if( method /= S2_PROJ_METHOD_NEAREST_NEIGHBOUR) then
        call s2_error(S2_ERROR_PROJ_METHOD_INVALID, 's2_sky2proj', &
             comment_add='Output projection operator only defined for nearst neighbour projection.')
     end if

     if (nside == 0) then
        call s2_proj_operator_nearest_neighbour(proj, nop=nop, op=op, nsphere=nsphere, xmap=xmap)
     else
        call s2_proj_operator_nearest_neighbour(proj, nside, nop, op, nsphere, xmap)
     end if

     ! Open file.
     fileid = 43
     open(unit=fileid, file=trim(filename_op), status='new', action='write', &
          form='formatted')

     write(fileid,'(a)') 'op = [...'
     do j = 0,nop-1
        write(fileid,'(i20,a,i20,a)') op(j,0), ', ', op(j,1), '; ...'
     end do
     write(fileid,'(a)') ']'

     write(fileid,'(a)') ''

     write(fileid,'(a)') 'xmap = [...'
     do j = 0,nsphere-1
        write(fileid,'(e20.10,a)') xmap(j), '; ...'
     end do
     write(fileid,'(a)') ']'

     ! Close file.
     close(fileid)

     ! Free memory.
     deallocate(xmap, op)

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
          write(*,*) 'option ', opt, ' has no argument'
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

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_sky2proj





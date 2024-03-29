!------------------------------------------------------------------------------
! s2_map2matmap
!
!! Read a HEALPix map fits file and extract an ecp (equispaced)
!! sampled theta-phi array over the sphere for the grid used for FSSHT.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_inp]: Name of map HEAPix fits file to read.
!!   - [-out filename_out]: Name of matlab map file to write.
!!   - [-ext extension]: Optional extension of HEALpix map fits file to read 
!!     map from.
!!   - [-B B]: Band-limit of sphere in written matlab map file.
!!   - [-grid_type grid_type]: Type of ECP grid to covert to (s2dw; mw; mwss).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_map2matmap

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_inp, filename_out
  type(s2_sky) :: sky
  integer :: ext = 1
  integer :: B = 128
  character(len=*), parameter :: grid_type_str_s2dw = 's2dw'
  character(len=*), parameter :: grid_type_str_mw = 'mw'
  character(len=*), parameter :: grid_type_str_mwss = 'mwss'
  character(len=S2_STRING_LEN) :: grid_type_str = grid_type_str_mw

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with map read in from map fits file.
  sky = s2_sky_init(filename_inp, S2_SKY_FILE_TYPE_MAP, ext)

  ! Write healpix map to matlab file.
  select case (trim(grid_type_str))

    case (grid_type_str_s2dw)
       call s2_sky_write_matmap_file(sky, filename_out, B, &
            grid_type=S2_SKY_ABGRID_TYPE_S2DW)

    case (grid_type_str_mw)
       call s2_sky_write_matmap_file(sky, filename_out, B, &
            grid_type=S2_SKY_ABGRID_TYPE_MW)

    case (grid_type_str_mwss)
       call s2_sky_write_matmap_file(sky, filename_out, B, &
            grid_type=S2_SKY_ABGRID_TYPE_MWSS)

    case default
       call s2_error(S2_ERROR_SKY_ABGRID_TYPE_INVALID, &
            's2_map2matmap')

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
            write(*,'(a)') 'Usage: s2_map2matmap [-inp filename_inp]'
            write(*,'(a)') '                     [-out filename_out]'
            write(*,'(a)') '                     [-ext ext (optional)]'
            write(*,'(a)') '                     [-B B]'
            write(*,'(a)') '                     [-grid_type grid_type (s2dw; mw; mwss)]'
            stop
          
          case ('-inp')
            filename_inp = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-ext')
            read(arg,*) ext

         case ('-B')
            read(arg,*) B

          case ('-grid_type')
            grid_type_str = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_map2matmap

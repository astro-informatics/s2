!------------------------------------------------------------------------------
! s2_skymask
!
!! Apply a mask (containing only ones and zeros) to a sky.  The output map is
!! the product of the mask and sky maps.  If the display status is set then
!! masked pixels of the output map are overwritten with a magic number that
!! appears grey when plotted.  Output maps producted with the display option
!! set should *only* be used for display purposes, and *not* for any
!! subsequent analysis.
!!
!! Usage: s2_skymask
!!   - [-help]: Display usage information.
!!   - [-sky filename_sky]: Name of input file containing sky to mask.
!!   - [-mask filename_mask]: Name of input file containined mask to apply.
!!   - [-display display]:  Logical specify whether to produce an output
!!     masked map for displya purposes only, in which case the masked
!!     values are set to a magic number so that they appear grey when plotted.
!!   - [-out filename_out]: Name of output file to write containing masked sky.
!!   - [-file_type file_type_str (map or sky)]: String specifying file types.
!!   - [-ext ext (optional)]: File extension for map files.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - August 2005
!
! Revisions:
!   August 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_skymask

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  real(s2_sp), parameter :: FITS_DISPLAY_GREY_MAGIC_NUMBER = -1.6375e30
  real(s2_sp), parameter :: ZERO_TOL = 1e-5

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: filename_sky, filename_mask, filename_out
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  type(s2_sky) :: sky, mask, prod
  real(s2_sp), allocatable :: mask_map(:), prod_map(:)
  logical :: display = .true.
  integer :: npix, nside, pix_scheme, lmax, mmax, fail

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skymask', &
         comment_add='Invalid file type option')

  end select

  ! Initialse sky with map read in from fits file.
  sky = s2_sky_init(filename_sky, file_type, ext)
  mask = s2_sky_init(filename_mask, file_type, ext)

  ! Take product of sky and mask.
  ! Checks skies are of same size and same pixel scheme
  ! (converts pixel scheme if required).
  prod = s2_sky_product(sky, mask)

  ! If display set, then overwrite masked pixels of product map with
  ! magic number that is displayed as grey.
  if(display) then

    ! Get product map and original mask map.
    npix = s2_sky_get_npix(sky)
    nside = s2_sky_get_nside(sky)
    pix_scheme = s2_sky_get_pix_scheme(sky)
    lmax = s2_sky_get_lmax(sky)
    mmax = s2_sky_get_mmax(sky)
    allocate(mask_map(0:npix-1), stat=fail)
    allocate(prod_map(0:npix-1), stat=fail)
    if(fail /= 0) then
      call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_skymask')
    end if
    call s2_sky_get_map(mask, mask_map)
    call s2_sky_get_map(prod, prod_map)

    ! Overwrite masked pixels with magic number.
    where(mask_map < ZERO_TOL)
      prod_map = FITS_DISPLAY_GREY_MAGIC_NUMBER
    end where

    ! Clear old prod sky and replace with new version.
    call s2_sky_free(prod)
    prod = s2_sky_init(prod_map, nside, pix_scheme, lmax, mmax)

    ! Free temporary memory used in this block.
    deallocate(mask_map, prod_map)
  
  end if

  ! Save output prod map.
  select case(file_type)

    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(prod, filename_out)

    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, prod)
       
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_skymask', &
            comment_add='Invalid file type specifier')
       
  end select

  ! Free memory.
  call s2_sky_free(sky)
  call s2_sky_free(mask)
  call s2_sky_free(prod)


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
            write(*,'(a)') 'Usage: s2_skymask [-sky filename_sky]'
            write(*,'(a)') '                   [-mask filename_mask]'
            write(*,'(a)') '                   [-display display]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-file_type file_type_str (map or sky)]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            stop
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-mask')
            filename_mask = trim(arg)

          case ('-display')
            read(arg,*) display

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type_str = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_skymask





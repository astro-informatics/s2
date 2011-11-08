!------------------------------------------------------------------------------
! s2_sky2map
!
!! Read a full s2_sky file and write the map to a HEALPix fits file.
!!
!! Notes:
!!   - If the sky is stored only as alms then the map is computed at the nside 
!!     specified in the sky file.  If no nside is specified in the sky file 
!!     then one may be specified as an input argument to this program or else 
!!     an error will occur.
!! 
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-sky filename_sky]: Name of full s2_sky fits file.
!!   - [-map filename_map]: Name of map HEAPix fits file to write.
!!   - [-nside nside]: Optional nside of output map fits file (see note above).
!!   - [-beta beta]: Optional beta rotation to rotate north pole to center
!!      of map for optimal viewing.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_sky2map

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_sky, filename_map
  type(s2_sky) :: sky
  integer :: nside
  logical :: new_nside = .false., beta_center = .false. 
  logical :: norm_pres = .false.
  logical :: dil1_status = .false., dil2_status = .false.
  real(s2_sp) :: dil1, dil2

  ! Parse input parameters.
  call parse_options()

  ! Read full sky object.
  sky = s2_sky_init(filename_sky, S2_SKY_FILE_TYPE_SKY)

  ! Compute map if not present. 
  ! If map is present then these calls won't do anything.
  ! If alms are not defined then no problem since routines will return 
  ! before checking if alms are present if map is defined.
  if(new_nside) then
     ! Map will be calculated at new nside.
     ! If nside is same as current sky%nside and map is already computed
     ! then nothing will happen (routine will simply return).
     call s2_sky_compute_map(sky, nside)
  else
     if(.not. s2_sky_get_map_status(sky)) then  ! Not necessary, but
                                                 ! eliminates warning message.
        call s2_sky_compute_map(sky)
     end if
  end if

  ! Dilate sky if required.
  if(dil1_status .or. dil2_status) then

     ! If only one dilation parameter set then set as symmetric dilation.
     if(.not. dil1_status) dil1 = dil2
     if(.not. dil2_status) dil2 = dil1

     ! Dilate sky.
     call s2_sky_dilate(sky, dil1, dil2, norm_pres)

  end if


  ! Rotate map by beta if required.
  if(beta_center) then
     call s2_sky_rotate(sky, 0.0e0, pi/2.0e0, 0.0e0)
  end if


  ! Write healpix map to fits file.
  call s2_sky_write_map_file(sky, filename_map)

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
            write(*,'(a)') 'Usage: s2_sky2map [-sky filename_sky]'
            write(*,'(a)') '                   [-map filename_map]'
            write(*,'(a)') '                   [-nside nside (optional)]'
            write(*,'(a,a)') '                   ', &
              '[-beta_center beta_center (optional)]'
            write(*,'(a)') '                   [-dil1 dil1]'
            write(*,'(a)') '                   [-dil2 dil2]'
            write(*,'(a)') '                   [-npres norm_pres]'
         stop
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-map')
            filename_map = trim(arg)

          case ('-nside')
            read(arg,*) nside
            new_nside = .true.

          case ('-beta_center')
            read(arg,*) beta_center

          case ('-dil1')
            read(arg,*) dil1
            dil1_status = .true.

          case ('-dil2')
            read(arg,*) dil2
            dil2_status = .true.

          case ('-npres')
            read(arg,*) norm_pres

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_sky2map

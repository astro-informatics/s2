!------------------------------------------------------------------------------
! s2_sky2cl
!
!! Read a map from a fits file, compute the alms of the map and then the cl
!! spectrum.  Write the computed cl spectrum values to the standard output.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of sky fits file containing map to compute
!!     cls of.
!!   - [-out filename_out]: Name of ouput pl file.
!!   - [-beam filename_beam]: Name of beam file to convolve with only if
!!     present.
!!   - [-type file_type (''map'' or ''sky'')]: String to specify type of 
!!     input file.
!!   - [-lmax lmax]: lmax of alm representation of map.
!!   - [-mmax mmax]: mmax of alm representation of map.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - March 2005
!
! Revisions:
!   March 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_sky2cl

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_pl_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  character(len=S2_STRING_LEN) :: filename_in, filename_out, filename_beam
  character(len=S2_STRING_LEN) :: file_type = SKY_FILE
  type(s2_sky) :: sky
  type(s2_pl) :: cl, beam
  integer :: lmax=0, mmax=0
  logical :: new_lmax = .false., new_mmax = .false.
  logical :: dil1_status = .false., dil2_status = .false.
  logical :: norm_pres = .false.
  logical :: beam_status = .false.
  real(s2_sp) :: dil1, dil2

  ! Read command line options.
  call parse_options()

  ! Read specified sky file.
  select case (trim(file_type))

    case (MAP_FILE)
       sky = s2_sky_init(filename_in, S2_SKY_FILE_TYPE_MAP) 

    case (SKY_FILE)
       sky = s2_sky_init(filename_in, S2_SKY_FILE_TYPE_SKY) 

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_sky2cl', &
         comment_add='Invalid file type option')

  end select

  ! Dilate sky if required.
  if(dil1_status .or. dil2_status) then

     ! If only one dilation parameter set then set as symmetric dilation.
     if(.not. dil1_status) dil1 = dil2
     if(.not. dil2_status) dil2 = dil1

     ! If dilate will definitely need to recompute alms.  
     ! Set lmax and mmax and status accordingly.
     if(.not. new_lmax) then
        lmax = s2_sky_get_lmax(sky)
        mmax = s2_sky_get_mmax(sky)
        new_lmax = .true.
        new_mmax = .true.
     end if

     ! Dilate sky.
     call s2_sky_dilate(sky, dil1, dil2, norm_pres)

  end if

  ! Compute alms if required.
  ! If alms are already present and no new lmax or mmax given, then
  ! s2_sky_compute_alm will return without computing any new alms.
  if(new_lmax .and. new_mmax) then
     call s2_sky_compute_alm(sky, lmax, mmax)
  elseif(new_lmax) then
     mmax = lmax
     call s2_sky_compute_alm(sky, lmax, mmax)
  else
     call s2_sky_compute_alm(sky)
  end if

  ! Get cls.
  cl = s2_sky_get_cl(sky)

  ! Convolve with beam if beam is present.
  if(beam_status) then
     beam = s2_pl_init(filename_beam)
     call s2_pl_conv(cl, beam)
  end if

  ! Write cls to output file.
  call s2_pl_io_fits_write(filename_out, cl)

  ! Free memory.
  call s2_pl_free(cl)
  if(beam_status) call s2_pl_free(beam)
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
            write(*,'(a)') 'Usage: s2_sky2cl [-inp filename_in]'
            write(*,'(a)') '                  [-out filename_out]'
            write(*,'(a)') '                  [-beam filename_beam]'
            write(*,'(a,a)') '                  ', &
                 '[-type file_type (''map'' or ''sky'')]'
            write(*,'(a)') '                  [-lmax lmax]'
            write(*,'(a)') '                  [-mmax mmax]' 
            write(*,'(a)') '                  [-dil1 dil1]'
            write(*,'(a)') '                  [-dil2 dil2]'
            write(*,'(a)') '                  [-npres norm_pres]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-beam')
            filename_beam = trim(arg)
            beam_status = .true.

          case ('-type')
            file_type = trim(arg)

          case ('-lmax')
            read(arg,*) lmax
            new_lmax = .true.

          case ('-mmax')
            read(arg,*) mmax
            new_mmax = .true.

          case ('-dil1')
            read(arg,*) dil1
            dil1_status = .true.

          case ('-dil2')
            read(arg,*) dil2
            dil2_status = .true.

          case ('-npres')
            read(arg,*) norm_pres

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_sky2cl

!------------------------------------------------------------------------------
! s2_map2sky
!
!! Read a HEALPix map fits file and write the read map to a full s2_sky fits
!! file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-map filename_map]: Name of map HEAPix fits file to read.
!!   - [-sky filename_sky]: Name of full s2_sky fits file to write.
!!   - [-ext extension]: Optional extension of HEALpix map fits file to read 
!!     map from.
!!   - [-alm compute_alm]: Logical specifying whether to compute alms.
!!   - [-lmax lmax]: Maximum harmonic l to consider when computing alms.
!!   - [-mmax mmax]: Maximum harmonic m to consider when computing alms.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_map2sky

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_sky, filename_map
  type(s2_sky) :: sky
  integer :: ext = 1
  logical :: compute_alm = .false.
  integer :: lmax = 128, mmax = 128

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with map read in from map fits file.
  sky = s2_sky_init(filename_map, S2_SKY_FILE_TYPE_MAP, ext)

  ! Compute alms if required.
  if(compute_alm) then
     call s2_sky_compute_alm(sky, lmax, mmax)
  end if

  ! Write sky to full s2_sky file.
  call s2_sky_io_fits_write(filename_sky, sky)

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
            write(*,'(a)') 'Usage: s2_map2sky [-map filename_map]'
            write(*,'(a)') '                   [-sky filename_sky]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            write(*,'(a)') '                   [-alm compute_alm]'
            write(*,'(a)') '                   [-lmax lmax]'
            write(*,'(a)') '                   [-mmax mmax]'
            stop
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-map')
            filename_map = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case ('-alm')
            read(arg,*) compute_alm

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_map2sky

!------------------------------------------------------------------------------
! s2_map2alm
!
!! Read a HEALPix map fits file and write the read map to a HEALPix fits alm
!! file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-map filename_map]: Name of map HEAPix fits file to read.
!!   - [-alm filename_alm]: Name of alm HEALPix fits file to write.
!!   - [-ext extension]: Optional extension of HEALpix map fits file to read 
!!     map from.
!!   - [-lmax lmax]: Maximum harmonic l to consider when computing alms.
!!   - [-mmax mmax]: Maximum harmonic m to consider when computing alms.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2005
!
! Revisions:
!   November 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_map2alm

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_alm, filename_map
  type(s2_sky) :: sky
  integer :: ext = 1
  integer :: lmax = 128, mmax = 128, iter_order = 0

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with map read in from map fits file.
  sky = s2_sky_init(filename_map, S2_SKY_FILE_TYPE_MAP, ext)

  ! Compute alms.
  call s2_sky_compute_alm_iter(sky, iter_order, lmax, mmax)

  ! Write sky to alm file.
  call s2_sky_write_alm_file(sky, filename_alm)

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
            write(*,'(a)') 'Usage: s2_map2alm [-map filename_map]'
            write(*,'(a)') '                  [-alm filename_alm]'
            write(*,'(a)') '                  [-ext ext (optional)]'
            write(*,'(a)') '                  [-lmax lmax]'
            write(*,'(a)') '                  [-mmax mmax]'
            write(*,'(a)') '                  [-iter iter_order]'
            stop
          
          case ('-map')
            filename_map = trim(arg)

          case ('-alm')
            filename_alm = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case ('-iter')
            read(arg,*) iter_order

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_map2alm

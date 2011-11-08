!------------------------------------------------------------------------------
! s2_sky2alm
!
!! Read a HEALPix map fits file and write the read map to a full s2_sky fits
!! file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-sky filename_sky]: Name of full s2_sky fits file to read.
!!   - [-alm filename_alm]: Name of map HEAPix fits file to write.
!!   - [-alm_compute compute_alm]: Logical specifying whether to compute alms.
!!   - [-lmax lmax]: Maximum harmonic l to consider when computing alms.
!!   - [-mmax mmax]: Maximum harmonic m to consider when computing alms.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2005
!
! Revisions:
!   November 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_sky2alm

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_sky, filename_alm
  type(s2_sky) :: sky
  logical :: compute_alm = .false.
  integer :: lmax = 128, mmax = 128

  ! Parse input parameters.
  call parse_options()

  ! Read full sky object.
  sky = s2_sky_init(filename_sky, S2_SKY_FILE_TYPE_SKY)

  ! Compute alms if required.
  ! Alm may already exist, in which case don't need to compute.
  if(compute_alm) then
     call s2_sky_compute_alm(sky, lmax, mmax)
  end if

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
            write(*,'(a)') 'Usage: s2_sky2alm [-sky filename_sky]'
            write(*,'(a)') '                  [-alm filename_alm]'
            write(*,'(a)') '                  [-compute_alm compute_alm]'
            write(*,'(a)') '                  [-lmax lmax]'
            write(*,'(a)') '                  [-mmax mmax]'
            stop
          
          case ('-sky')
            filename_sky = trim(arg)

          case ('-alm')
            filename_alm = trim(arg)

          case ('-compute_alm')
            read(arg,*) compute_alm

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_sky2alm

!------------------------------------------------------------------------------
! s2_gcmb
!
!! Generate a Gaussian CMB realisation from an input power spectrum.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_cl]: Name of cl input ascii file.
!!   - [-scale_cl scale_cl (logical)]: Logical to specify whether to scale
!!     the input cl values by 2*pi/(l(l+1)).
!!   - [-out filename_out]: Name of output cmb map file.
!!   - [-nside nside]: Healpix nside resolution to generate map for.
!!   - [-lmax lmax]: Maximum hamonic l to consider in power spectrum when
!!     generating map.
!!   - [-lmin lmin]: Minimum l of power spectrum contained in input file
!!     (others are set to zero).  Often this is 2.
!!   - [-ncomment ncomment]: Number of comment lines in input file to ignore.
!!   - [-seed seed]: Integer seed for random number generator used to
!!     computed Gaussian realised CMB. If not specified then use system clock.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - July 2005
!
! Revisions:
!   July 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_gcmb

  use s2_types_mod
  use s2_cmb_mod

  implicit none

  type(s2_cmb) :: cmb
  integer :: nside, lmin, lmax, ncomment, seed
  character(S2_STRING_LEN) :: filename_cl
  character(S2_STRING_LEN) :: filename_out = 'out.fits'
  logical :: scale_cl = .true.
  logical :: seed_set = .false.

  ! Set default parameter values.
  nside = 512
  lmax = 1024
  lmin = 2
  ncomment = 29   ! Value required to read WMAP cl files.
  seed = 1
  scale_cl = .true.
  
  ! Parse options.
  call parse_options()

  ! Use system clock to set seed if not already specified.
  if(.not. seed_set) then
    call system_clock(seed)
  end if

  ! Construct simulated cmb map.
  cmb = s2_cmb_init(filename_cl, nside, lmin, lmax, ncomment, &
      seed, scale_cl=scale_cl)

  ! Write computed map to output file.
  call s2_cmb_write_sky(cmb, filename_out)

  ! Free memory.
  call s2_cmb_free(cmb)


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
            write(*,'(a)') 'Usage: s2_gcmb [-inp filename_cl]'
            write(*,'(a)') '                [-scale_cl scale_cl (logical)]'
            write(*,'(a)') '                [-out filename_out]'
            write(*,'(a)') '                [-nside nside]'
            write(*,'(a)') '                [-lmax lmax]'
            write(*,'(a)') '                [-lmin lmin]'
            write(*,'(a)') '                [-ncomment ncomment]'
            write(*,'(a)') '                [-seed seed]'
            stop
          
          case ('-inp')
            filename_cl = trim(arg)

          case ('-scale_cl')
            read(arg,*) scale_cl
            
          case ('-out')
            filename_out = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax

          case ('-lmin')
            read(arg,*) lmin

          case ('-ncomment')
            read(arg,*) ncomment

          case ('-seed')
            read(arg,*) seed
            seed_set = .true.

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_gcmb
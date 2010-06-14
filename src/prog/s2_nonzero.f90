!------------------------------------------------------------------------------
! s2_nonzero
!
!! Count the number of non-zero values in a map.  Two techniques are used.
!! One assumes the map is binary and sums the values.  The other counts the
!! number of values above a threshold.  The results are printed to the screen.
!!
!! Usage: cswt_mask_nonzero
!!   - [-help]: Display usage information.
!!   - [-inp filename_mask]: Name of sky fits file containing map to count.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - December 2004
!
! Revisions:
!   December 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_nonzero

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename
  type(s2_sky) :: sky
  real(s2_sp), allocatable :: map(:)
  integer :: npix, i, num_nonzero_count = 0
  real(s2_sp) :: num_nonzero_sum
  real(s2_sp), parameter :: ZERO_TOL = 0.01

  ! Read filename.
  call parse_options()

  ! Read specified sky file.
  sky = s2_sky_init(filename, S2_SKY_FILE_TYPE_MAP)

  ! Allcoate space for map.
  npix = s2_sky_get_npix(sky)
  allocate(map(0:npix-1))

  ! Get map.
  call s2_sky_get_map(sky, map)
  
  ! Count number of ones in map.

  ! Summing approach (assuming pixels exactly 1).
  num_nonzero_sum = sum(map)

  ! Counting approach (count every pixel greater than ZERO_TOL).
  do i = 0,npix-1
     if(abs(map(i)) > ZERO_TOL) num_nonzero_count = num_nonzero_count + 1
  end do

  ! Print results.
  write(*,'(a,a)') 'Filename: ', trim(filename)
  write(*,'(a,i46)') 'Total number of pixels: ', npix
  write(*,'(a,i25)') 'Number of non-zero values (count technique): ', &
    num_nonzero_count
  write(*,'(a,f21.4)') 'Proportion of non-zero values (count technique): ', &
    num_nonzero_count / real(npix,s2_sp)
  write(*,'(a,a,f11.1)') 'Number of one values, assuming binary map ', &
    '(sum technique): ', num_nonzero_sum
  write(*,'(a,a,f7.4)') 'Proportion of one values, assuming binary map ', &
    '(sum technique): ', num_nonzero_sum / real(npix,s2_sp)

  ! Free memory.
  call s2_sky_free(sky)
  deallocate(map)


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
            write(*,'(a)') 'Usage: s2_nonzero [-inp filename]'
            stop
          
          case ('-inp')
            filename = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_nonzero

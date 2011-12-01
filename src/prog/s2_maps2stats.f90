!------------------------------------------------------------------------------
! s2_maps2stats
!



!! Compute the axisymmetric convolution of a sky with a kernel in
!! harmonic space.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of file containing input sky.
!!   - [-inp filename_in]: Name of file containing convolution kernel.
!!   - [-out filename_out]: Name of output file for convolved sky.
!!   - [-file_type_in file_type_in_str]: String specifying input file type.
!!   - [-file_type_ker file_type_ker_str]: String specifying kernel file type.
!!   - [-file_type_out file_type_out_str]: String specifying output file type.
!!   - [-lmax lmax]: Maximum harmonic l to consider.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   December 2011 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_maps2stats

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_mean, filename_std

  type(s2_sky) :: sky
  type(s2_sky) :: mean, mean_tmp
  type(s2_sky) :: mom2, mom2_tmp
  type(s2_sky) :: std

  integer :: nside, nside_check
  integer :: fileid = 31, iostat
  character(len=S2_STRING_LEN) :: line

  integer ::  n = 0

!!$  interface 
!!$     function square(x) result(val)
!!$       use s2_types_mod
!!$       real(s2_sp), intent(in) :: x
!!$       real(s2_sp) :: val
!!$     end function square
!!$  end interface
!!$  interface 
!!$     function square_root(x) result(val)
!!$       use s2_types_mod
!!$       real(s2_sp), intent(in) :: x
!!$       real(s2_sp) :: val
!!$     end function square_root
!!$  end interface


  
  ! Parse input parameters.
  call parse_options()

  ! Open file containing list of maps.
  open(fileid, file=filename_in, form='formatted', status='old')

  ! Read maps and compute moments.
  do
     read(fileid,'(a)',iostat=iostat) line
     if (iostat < 0) exit

     ! Print status update.
     write(*,'(a,a,a,a)') 's2_maps2stats> ', &
          'Processing map ', &
          trim(line), '...'

     ! Read map.
     sky = s2_sky_init(trim(line), S2_SKY_FILE_TYPE_MAP)
     n = n + 1

     if (n == 1) then

        ! Initialise moments if first map.
        mean = s2_sky_init(sky)
        mom2 = s2_sky_init(sky)
        call s2_sky_fun(mom2, square)

     else

        ! Update mean.
        mean_tmp = s2_sky_add(mean, sky)
        call s2_sky_free(mean)
        mean = s2_sky_init(mean_tmp)
        call s2_sky_free(mean_tmp)

        ! Update second moment.
        call s2_sky_fun(sky, square)
        mom2_tmp = s2_sky_add(mom2, sky)
        call s2_sky_free(mom2)
        mom2 = s2_sky_init(mom2_tmp)
        call s2_sky_free(mom2_tmp)

     end if

     ! Free current map.
     call s2_sky_free(sky)

  end do

  ! Close file.
  close(fileid)

  ! Compute mean and save map.
  call  s2_sky_scale(mean, 1.0/real(n,s2_sp))
  call s2_sky_write_file(mean, filename_mean, S2_SKY_FILE_TYPE_MAP)

  ! Compute standard deviation and save map.
  call  s2_sky_scale(mom2, 1.0/real(n-1,s2_sp))
  call s2_sky_fun(mean, square)
  std = s2_sky_add(mom2, mean, subtract=.true.)
  call s2_sky_fun(std, square_root)
  call s2_sky_write_file(std, filename_std, S2_SKY_FILE_TYPE_MAP)

  ! Free memory.
  call s2_sky_free(mean)
  call s2_sky_free(mom2)
  call s2_sky_free(std)


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
            write(*,'(a)') 'Usage: s2_maps2stats [-inp filename_in]'
            write(*,'(a)') '                     [-out_mean filename_mean]'
            write(*,'(a)') '                     [-out_std filename_std]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out_mean')
            filename_mean = trim(arg)

          case ('-out_std')
            filename_std = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


    !---------------------------------------------------------------------
    ! square
    !
    !! Function to square a map pixel value.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - December 2012
    !
    ! Revisions:
    !   December 2012 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    function square(x) result (val)

      use s2_types_mod
      
      implicit none

      real(s2_sp), intent(in) :: x
      real(s2_sp) :: val

      val = x * x

    end function square


    !---------------------------------------------------------------------
    ! square_root
    !
    !! Function to square root a map pixel value.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - December 2012
    !
    ! Revisions:
    !   December 2012 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    function square_root(x) result (val)

      use s2_types_mod
      
      implicit none

      real(s2_sp), intent(in) :: x
      real(s2_sp) :: val

      val = sqrt(x)

    end function square_root


end program s2_maps2stats


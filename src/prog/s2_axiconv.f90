!------------------------------------------------------------------------------
! s2_axiconv
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
!   November 2011 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_axiconv

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type_in = S2_SKY_FILE_TYPE_ALM
  integer :: file_type_ker = S2_SKY_FILE_TYPE_ALM
  integer :: file_type_out = S2_SKY_FILE_TYPE_ALM

  character(len=S2_STRING_LEN) :: filename_in, filename_ker, filename_out
  character(len=S2_STRING_LEN) :: file_type_in_str = ALM_FILE
  character(len=S2_STRING_LEN) :: file_type_ker_str = ALM_FILE
  character(len=S2_STRING_LEN) :: file_type_out_str = ALM_FILE

  type(s2_sky) :: sky_in
  type(s2_sky) :: kernel
  type(s2_sky) :: sky_out
  integer :: lmax, lmax_sky, lmax_ker
  integer :: l, m, fail = 0

  complex(s2_spc), allocatable :: sky_in_alm(:,:)
  complex(s2_spc), allocatable :: kernel_alm(:,:)
  complex(s2_spc), allocatable :: sky_out_alm(:,:)

  ! Parse input parameters.
  call parse_options()

  ! Set input sky file type.
  select case (trim(file_type_in_str))
    case (MAP_FILE)
       file_type_in = S2_SKY_FILE_TYPE_MAP
    case (ALM_FILE)
       file_type_in = S2_SKY_FILE_TYPE_ALM
    case (SKY_FILE)
       file_type_in = S2_SKY_FILE_TYPE_SKY
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_axiconv', &
         comment_add='Invalid file type option')
  end select

  ! Set kernel sky file type.
  select case (trim(file_type_ker_str))
    case (MAP_FILE)
       file_type_ker = S2_SKY_FILE_TYPE_MAP
    case (ALM_FILE)
       file_type_ker = S2_SKY_FILE_TYPE_ALM
    case (SKY_FILE)
       file_type_ker = S2_SKY_FILE_TYPE_SKY
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_axiconv', &
         comment_add='Invalid file type option')
  end select

  ! Set out sky file type.
  select case (trim(file_type_out_str))
    case (MAP_FILE)
       file_type_out = S2_SKY_FILE_TYPE_MAP
    case (ALM_FILE)
       file_type_out = S2_SKY_FILE_TYPE_ALM
    case (SKY_FILE)
       file_type_out = S2_SKY_FILE_TYPE_SKY
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_axiconv', &
         comment_add='Invalid file type option')
  end select

  ! Read sky and kernel from files.
  sky_in = s2_sky_init(filename_in, file_type_in)
  kernel = s2_sky_init(filename_ker, file_type_ker)

  ! Check lmax consistent if defined.
  lmax_sky = s2_sky_get_lmax(sky_in)
  lmax_ker = s2_sky_get_lmax(kernel)
  if (lmax_sky /= 0 .and. lmax_sky /= lmax) then
     call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_axiconv', &
          comment_add='Inconsistent lmax for sky')
  end if
  if (lmax_ker /= 0 .and. lmax_ker /= lmax) then
     call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_axiconv', &
          comment_add='Inconsistent lmax for kernel')
  end if

  ! Compute harmonic coefficients of sky and kernel if not already
  ! computed.
  if (.not. s2_sky_get_alm_status(sky_in)) then
     call s2_sky_compute_alm(sky_in, lmax, lmax)
  end if
  if (.not. s2_sky_get_alm_status(kernel)) then
     call s2_sky_compute_alm(kernel, lmax, lmax)
  end if

  ! Get copies of alms.
  allocate(sky_in_alm(0:lmax, 0:lmax), stat=fail)
  allocate(kernel_alm(0:lmax, 0:lmax), stat=fail)
  allocate(sky_out_alm(0:lmax, 0:lmax), stat=fail)
  sky_in_alm(0:lmax, 0:lmax) = 0.0e0
  kernel_alm(0:lmax, 0:lmax) = 0.0e0
  sky_out_alm(0:lmax, 0:lmax) = 0.0e0
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
          's2_axiconv')
  end if
  call s2_sky_get_alm(sky_in, sky_in_alm)
  call s2_sky_get_alm(kernel, kernel_alm)

  ! Perform convolution.
  do l = 0,lmax 
     do m = 0,lmax
        sky_out_alm(l,m) = sqrt(4e0*pi/real(2*l+1,s2_sp)) &
             * kernel_alm(l,0) * sky_in_alm(l,m)
     end do
  end do

  ! Initialise convolved sky.
  sky_out = s2_sky_init(sky_out_alm, lmax, lmax, &
       s2_sky_get_nside(kernel), s2_sky_get_pix_scheme(sky_in))

  ! Compute map if required for output.
  if (file_type_out == S2_SKY_FILE_TYPE_MAP) then
     call s2_sky_compute_map(sky_out)
  end if
  
  ! Save output file.
  call s2_sky_write_file(sky_out, filename_out, file_type_out)

  ! Free memory.
  call s2_sky_free(sky_in)
  call s2_sky_free(kernel)
  call s2_sky_free(sky_out)
  deallocate(sky_in_alm)
  deallocate(kernel_alm)
  deallocate(sky_out_alm)


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
            write(*,'(a)') 'Usage: s2_axiconv [-inp filename_in]'
            write(*,'(a)') '                  [-ker filename_ker]'
            write(*,'(a)') '                  [-out filename_out]'
            write(*,'(a,a)') '                  ', &
                 '[-file_type_in file_type_in_str (sky; map; alm)]'
            write(*,'(a,a)') '                  ', &
                 '[-file_type_ker file_type_ker_str (sky; map; alm)]'
            write(*,'(a,a)') '                  ', &
                 '[-file_type_out file_type_out_str (sky; map; alm)]'
            write(*,'(a)') '                  [-lmax lmax]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-ker')
            filename_ker = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type_in')
            file_type_in_str = trim(arg)

          case ('-file_type_ker')
            file_type_ker_str = trim(arg)

          case ('-file_type_out')
            file_type_out_str = trim(arg)

          case ('-lmax')
            read(arg,*) lmax

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_axiconv


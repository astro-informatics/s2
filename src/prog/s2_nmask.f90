!------------------------------------------------------------------------------
! s2_nmask
!
!! Generate noisy masks and computed harmonic space noise covariance matrix
!! (not yet outputted) from a binary mask.
!!
!! Usage: s2_nmask
!!   - [-help]: Display usage information.
!!   - [-mask filename_mask]: Name of input file containined mask construct 
!!     noisy mask from.
!!   - [-lmax lmax]: Harmonic lmax for computing covariance matrix.
!!   - [-nsigma nsigma]: Noise standard deviation.
!!   - [-out filename_out]: Name of output noisy mask.
!!   - [-file_type file_type_str (map or sky)]: String specifying file types.
!!   - [-ext ext (optional)]: File extension for map files.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - February 2008
!
! Revisions:
!   February 2008 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_nmask

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_distn_mod
  use s2_ylm_mod

  implicit none

  real(s2_sp), parameter :: ZERO_TOL = 1e-5

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type = S2_SKY_FILE_TYPE_MAP, ext = 1
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  character(len=S2_STRING_LEN) :: filename_mask, filename_out

  type(s2_sky) :: mask, nmask
  real(s2_sp), allocatable :: mask_map(:), nmask_map(:)
  type(s2_sky) :: ylm1_real, ylm1_imag, ylm2_real, ylm2_imag
  real(s2_sp), allocatable :: ylm1_real_map(:), ylm1_imag_map(:)
  real(s2_sp), allocatable :: ylm2_real_map(:), ylm2_imag_map(:)
  complex(s2_spc), allocatable :: cov(:,:,:,:)

  integer :: ipix, npix, nside, pix_scheme, lmax, fail, seed
  integer :: el1, m1, el2, m2
  real(s2_sp) :: nsigma = 1e0
  complex(s2_dpc), parameter :: I = cmplx(0d0,1d0)

  ! Parse input parameters.
  call parse_options()

  ! Set sky file type.
  select case (trim(file_type_str))

    case (MAP_FILE)
       file_type = S2_SKY_FILE_TYPE_MAP

    case (SKY_FILE)
       file_type = S2_SKY_FILE_TYPE_SKY

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_nmask', &
         comment_add='Invalid file type option')

  end select


  !----------------------------------------------------------------------------
  ! Construct noisy mask
  !----------------------------------------------------------------------------

  mask = s2_sky_init(filename_mask, file_type, ext)

  npix = s2_sky_get_npix(mask)
  nside = s2_sky_get_nside(mask)
  pix_scheme = s2_sky_get_pix_scheme(mask)
  fail = 0
  allocate(mask_map(0:npix-1), stat=fail)
  allocate(nmask_map(0:npix-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_nmask')
  end if

  call s2_sky_get_map(mask, mask_map)

  do ipix = 0,npix-1
     if(abs(mask_map(ipix)) > ZERO_TOL) then
        nmask_map(ipix) = 0e0
     else
        seed = 1
        nmask_map(ipix) = s2_distn_sample_gauss(seed, 0e0, nsigma)
     end if
  end do

  nmask = s2_sky_init(nmask_map, nside, pix_scheme)

  ! Save noisy mask map.
  select case(file_type)
    case(S2_SKY_FILE_TYPE_MAP)
       call s2_sky_write_map_file(nmask, filename_out)
    case(S2_SKY_FILE_TYPE_SKY)
       call s2_sky_io_fits_write(filename_out, nmask)
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_nmask', &
            comment_add='Invalid file type specifier')
  end select


  !----------------------------------------------------------------------------
  ! Compute harmonic space noise covariance matrix
  !----------------------------------------------------------------------------

  fail = 0
  allocate(ylm1_real_map(0:npix-1), stat=fail)
  allocate(ylm1_imag_map(0:npix-1), stat=fail)
  allocate(ylm2_real_map(0:npix-1), stat=fail)
  allocate(ylm2_imag_map(0:npix-1), stat=fail)
  allocate(cov(0:lmax,-lmax:lmax,0:lmax,-lmax:lmax), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_nmask')
  end if
  cov(0:lmax,-lmax:lmax,0:lmax,-lmax:lmax) = cmplx(0e0, 0e0)
!s1=0
!s2=0
  do el1 = 0,lmax
 write(*,*) 'el1=', el1
     do m1 = -el1,el1
!s1=s1+1
 write(*,*) 'm1=', m1
        ! Compute Ylm1 over full sky.
        ylm1_real = s2_ylm_sky(el1, m1, nside, S2_YLM_REALITY_REAL, &
             pix_scheme, S2_YLM_EVAL_METHOD_LEG)
        call s2_sky_get_map(ylm1_real, ylm1_real_map)
        call s2_sky_free(ylm1_real)
        ylm1_imag = s2_ylm_sky(el1, m1, nside, S2_YLM_REALITY_IMAG, &
             pix_scheme, S2_YLM_EVAL_METHOD_LEG)
        call s2_sky_get_map(ylm1_imag, ylm1_imag_map)
        call s2_sky_free(ylm1_imag)

        do el2 = 0,lmax
           do m2 = -el2,el2
!s2=s2+1
              ! Compute Ylm2 over full sky.
              ylm2_real = s2_ylm_sky(el2, m2, nside, &
                   S2_YLM_REALITY_REAL, &
                   pix_scheme, S2_YLM_EVAL_METHOD_LEG)
              call s2_sky_get_map(ylm2_real, ylm2_real_map)
              call s2_sky_free(ylm2_real)
              ylm2_imag = s2_ylm_sky(el2, m2, nside, &
                   S2_YLM_REALITY_IMAG, &
                   pix_scheme, S2_YLM_EVAL_METHOD_LEG)
              call s2_sky_get_map(ylm2_imag, ylm2_imag_map)
              call s2_sky_free(ylm2_imag)

              ! Sum over pixels to compute cov.
              do ipix=0,npix-1

                 ! Only take contributions from zero region of original
                 ! mask (i.e. non-zero region of noise mask).
                 if(abs(mask_map(ipix)) < ZERO_TOL) then
!cov(s1,s2)=...
                    cov(el1,m1,el2,m2) = cov(el1,m1,el2,m2) &
                     + ( &
                          ( ylm1_real_map(ipix) * ylm2_real_map(ipix) &
                            + ylm1_imag_map(ipix) * ylm2_imag_map(ipix) ) &
                       +I*( ylm1_imag_map(ipix) * ylm2_real_map(ipix) &
                            - ylm1_real_map(ipix) * ylm2_imag_map(ipix) ) &
                       ) &
                       * (nsigma**2) &
                       * 4*pi/real(npix,s2_sp)

                  end if

              end do

           end do
        end do

     end do
  end do

  ! Free memory.
  deallocate(mask_map, nmask_map)
  deallocate(ylm1_real_map, ylm1_imag_map, ylm2_real_map, ylm2_imag_map)
  deallocate(cov)
  call s2_sky_free(mask)
  call s2_sky_free(nmask)


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
            write(*,'(a)') 'Usage: s2_nmask [-mask filename_mask]'
            write(*,'(a)') '                [-lmax lmax]'
            write(*,'(a)') '                [-nsigma nsigma]'
            write(*,'(a)') '                [-out filename_out]'
            write(*,'(a)') '                [-file_type file_type_str (map or sky)]'
            write(*,'(a)') '                [-ext ext (optional)]'
            stop

          case ('-mask')
            filename_mask = trim(arg)

          case ('-lmax')
            read(arg,*) lmax

          case ('-nsigma')
            read(arg,*) nsigma

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


end program s2_nmask





!------------------------------------------------------------------------------
! s2_grf2
!
!! Simulate two Gaussian random fields (GRFs) on the sky that satisfy
!! individual and cross power spectra.
!!
!! Notes:
!!   - Not tested at all.
!!
!! Usage: s2_grf2
!!   - [-help]: Display usage information.
!!   - [-inp_al filename_al]: Name of a input ascii file containing power
!!     spectrum.
!!   - [-inp_bl filename_bl]: Name of b input ascii file containing power
!!     spectrum.
!!   - [-inp_xl filename_xl]: Name of x input ascii file containing power
!!     spectrum.
!!   - [-scale_al scale_al]: Logical to specify whether to scale
!!     the input cl values of al by 2*pi/(l(l+1)).
!!   - [-scale_bl scale_bl]: Logical to specify whether to scale
!!     the input cl values of bl by 2*pi/(l(l+1)).
!!   - [-scale_xl scale_xl]: Logical to specify whether to scale
!!     the input cl values of xl by 2*pi/(l(l+1)).
!!   - [-ncomment_al ncomment_al]: Number of comment lines in input a power
!!     spectrum file to ignore.
!!   - [-ncomment_bl ncomment_bl]: Number of comment lines in input b power
!!     spectrum file to ignore.
!!   - [-ncomment_xl ncomment_xl]: Number of comment lines in input x power
!!     spectrum file to ignore.
!!   - [-nside nside]: Healpix nside resolution to generate maps at.
!!   - [-lmax lmax]: Maximum hamonic l to consider in power spectrum when
!!     generating map.
!!   - [-lmin lmin]: Minimum l of power spectrum contained in input file
!!     (others are set to zero).  Often this is 2.
!!   - [-seed seed]: Integer seed for random number generator used to
!!     computed GRFs.  If not specified then use system clock.
!!   - [-outa filename_out_a]: Name of sky a output file.
!!   - [-outb filename_out_b]: Name of sky b output file.
!!   - [-file file_type_str (map or sky)]: Type of output files to write.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - August 2005
!
! Revisions:
!   August 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program s2_grf2

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_distn_mod
  use s2_pl_mod

  implicit none

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  character(len=S2_STRING_LEN) :: file_type_str = MAP_FILE
  character(len=S2_STRING_LEN) :: filename_out_a = 'a.tmp'
  character(len=S2_STRING_LEN) :: filename_out_b = 'b.tmp'
  character(len=S2_STRING_LEN) :: filename_al, filename_bl, filename_xl

  logical :: scale_al = .false., scale_bl = .false., scale_xl = .false.
  integer :: ncomment_al = 0, ncomment_bl = 0, ncomment_xl = 0
  logical :: line_nos

  integer :: seed = 0
  logical :: seed_set = .false.
  
  real(s2_sp) :: MEAN_ZERO = 0e0
  real(s2_sp) :: STD_UNIT = 1e0

  integer :: l, m, fail
  integer :: nside, lmin, lmax
  logical :: nside_set = .false.
  real(s2_dp) :: a, b, c, xr, yr, xi, yi
  type(s2_pl) :: a_cl, b_cl, x_cl
  real(s2_dp), allocatable :: al(:), bl(:), xl(:)
  real(s2_sp), allocatable :: al_sp(:), bl_sp(:), xl_sp(:)
  complex(s2_spc), allocatable :: alm(:,:), blm(:,:)
  type(s2_sky) :: a_sky, b_sky

  ! Init variables.
  lmin = 2
  lmax = 194
  ncomment_al = 0   ! 29  ! Value required to read WMAP cl files.
  ncomment_bl = 0
  ncomment_xl = 0
  seed = 1
  line_nos = .false.

  ! Parse input parameters.
  call parse_options()

  ! Use system clock to set seed if not already specified.
  if(.not. seed_set) then
    call system_clock(seed)
  end if

  ! Allocate space.
  allocate(al(0:lmax), stat=fail)
  allocate(bl(0:lmax), stat=fail)
  allocate(xl(0:lmax), stat=fail)
  allocate(al_sp(0:lmax), stat=fail)
  allocate(bl_sp(0:lmax), stat=fail)
  allocate(xl_sp(0:lmax), stat=fail)
  allocate(alm(0:lmax,0:lmax), stat=fail)
  allocate(blm(0:lmax,0:lmax), stat=fail)
  if(fail /= 0) then
    call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_grf2')
  end if

  ! Initialise variables allocated to zero.
  alm = cmplx(0d0,0d0)
  blm = cmplx(0d0,0d0)

  ! Read input.
  a_cl = s2_pl_init(filename_al, lmin, lmax, ncomment_al, scale_al, line_nos)
  b_cl = s2_pl_init(filename_bl, lmin, lmax, ncomment_bl, scale_bl, line_nos)
  x_cl = s2_pl_init(filename_xl, lmin, lmax, ncomment_xl, scale_xl, line_nos)

  ! Get power spectrum values from pl objects and cast as doubles.
  call s2_pl_get_spec(a_cl, al_sp)
  call s2_pl_get_spec(b_cl, bl_sp)
  call s2_pl_get_spec(x_cl, xl_sp)
  al = real(al_sp, s2_dp)
  bl = real(bl_sp, s2_dp)
  xl = real(xl_sp, s2_dp)

#ifdef DEBUG
  write(*,'(a16,a12,a12)') 'al', 'bl', 'xl'
  do l = 0, lmax
    write(*,'(i4,e12.4,e12.4,e12.4)') l, al(l), bl(l), xl(l)
  end do
#endif

  ! To avoid divide by zero for l=0 and 1 only consider from l=2 now.
  lmin = 2

  ! Compute harmonic coefficients of simulated fields.
  do l = lmin,lmax

    ! Check not taking sqrt of negative values.
    if(al(l) * bl(l) < xl(l)**2d0) then
      call s2_error(S2_ERROR_ARTH, 's2_grf2', &
        comment_add='Attempt to take sqrt of negative value')
    end if
    
    c = ( sqrt(al(l) * bl(l)) + sqrt(al(l) * bl(l) - xl(l)**2d0) ) / xl(l)
    b = sqrt( bl(l) / (2d0*(1d0 + c**2d0)) )
    a = sqrt( al(l) / (2d0*(1d0 + c**2d0)) )

    do m = 0,l

      xr = real(s2_distn_sample_gauss(seed, MEAN_ZERO, STD_UNIT), s2_dp)
      yr = real(s2_distn_sample_gauss(seed, MEAN_ZERO, STD_UNIT), s2_dp)
      xi = real(s2_distn_sample_gauss(seed, MEAN_ZERO, STD_UNIT), s2_dp)
      yi = real(s2_distn_sample_gauss(seed, MEAN_ZERO, STD_UNIT), s2_dp)

      alm(l,m) = real(cmplx(  a*xr + a*c*yr,   a*xi + a*c*yi), s2_spc)
      blm(l,m) = real(cmplx(b*c*xr +   b*yr, b*c*xi +   b*yi), s2_spc)

    end do

  end do

  ! Initialise sky structures from alms.
  a_sky = s2_sky_init(alm, lmax, lmax)
  b_sky = s2_sky_init(blm, lmax, lmax)

  ! Compute maps if nside set.
  if(nside_set) then
    call s2_sky_compute_map(a_sky, nside)
    call s2_sky_compute_map(b_sky, nside)
  end if 
   
  ! Save skies in file of appropriate format.
  select case (trim(file_type_str))

    case (MAP_FILE)
      if(nside_set) then
        call s2_sky_write_map_file(a_sky, filename_out_a)
        call s2_sky_write_map_file(b_sky, filename_out_b)
      else
        call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_grf2', &
          comment_add='Cannot save map file since map not computed')
      end if

    case (SKY_FILE)     
      call s2_sky_io_fits_write(filename_out_a, a_sky)
      call s2_sky_io_fits_write(filename_out_b, b_sky)

    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_grf2', &
         comment_add='Invalid file type option')

  end select

  ! Free memory.
  call s2_sky_free(a_sky)
  call s2_sky_free(b_sky)
  call s2_pl_free(a_cl)
  call s2_pl_free(b_cl)
  call s2_pl_free(x_cl)
  deallocate(alm, blm)
  deallocate(al, bl, xl)
  deallocate(al_sp, bl_sp, xl_sp)


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
            write(*,'(a)') 'Usage: s2_grf2 [-inp_al filename_al]'
            write(*,'(a)') '                [-inp_bl filename_bl]'
            write(*,'(a)') '                [-inp_xl filename_xl]'
            write(*,'(a)') '                [-scale_al scale_al]'
            write(*,'(a)') '                [-scale_bl scale_bl]'
            write(*,'(a)') '                [-scale_xl scale_xl]'
            write(*,'(a)') '                [-ncomment_al ncomment_al]'
            write(*,'(a)') '                [-ncomment_bl ncomment_bl]'
            write(*,'(a)') '                [-ncomment_xl ncomment_xl]'
            write(*,'(a)') '                [-line_nos line_nos]'
            write(*,'(a)') '                [-nside nside]'
            write(*,'(a)') '                [-lmax lmax]'
            write(*,'(a)') '                [-lmin lmin]'
            write(*,'(a)') '                [-seed seed]'
            write(*,'(a)') '                [-outa filename_out_a]'
            write(*,'(a)') '                [-outb filename_out_b]'
            write(*,'(a)') '                [-file file_type_str (map or sky)]'
            stop

          case ('-inp_al')
            filename_al = trim(arg)

          case ('-inp_bl')
            filename_bl = trim(arg)

          case ('-inp_xl')
            filename_xl = trim(arg)

          case ('-scale_al')
            read(arg,*) scale_al

          case ('-scale_bl')
            read(arg,*) scale_bl

          case ('-scale_xl')
            read(arg,*) scale_xl

          case ('-ncomment_al')
            read(arg,*) ncomment_al

          case ('-ncomment_bl')
            read(arg,*) ncomment_bl
          
          case ('-ncomment_xl')
            read(arg,*) ncomment_xl

          case ('-line_nos')
            read(arg,*) line_nos
            
          case ('-nside')
            nside_set = .true.
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax

          case ('-lmin')
            read(arg,*) lmin

          case ('-seed')
            read(arg,*) seed
            seed_set = .true.

          case ('-outa')
            filename_out_a = trim(arg)

          case ('-outb')
            filename_out_b = trim(arg)

          case ('-file')
            read(arg,*) file_type_str
                     
          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_grf2

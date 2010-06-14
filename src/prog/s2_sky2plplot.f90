!------------------------------------------------------------------------------
! s2_sky2plplot
!
!! Read a full fits s2_sky file and compute and plot the cl spectrum (dilate 
!! sky also if required, see notes below).
!!
!! Notes:
!!   - When constructing optimal filters, template sky is saved for dilation 
!!     of one, when actual template used may be at different scale.  Thus
!!     functionality incorporated here to dilate sky before computing and 
!!     plotting cl spectrum.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_sky]: Name of full s2_sky fits file containing sky to
!!      compute and plot pls of.
!!   - [-out filename_out]: Name of output postscript to print.
!!   - [-scale scale]: Logical to specify whether to scale y axis values to 
!!     plot dls rather than cls.
!!   - [-title title]: Optional title to print on plot.
!!   - [-lmax lmax]: Spherical harmonic lmax if computing alms.
!!   - [-mmax mmax]: Spherical harmonic mmax if computing alms.
!!   - [-dil1 dil1]: First dilation value if dilating sky before computing 
!!     alms.
!!   - [-dil2 dil2]: Second dilation value if dilating sky before computing 
!!     alms.
!!   - [-npres norm_pres]: Logical to specify whether norm preserving dilation 
!!     is to be performed, if indeed performing a dilation.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_sky2plplot

  use s2_types_mod
  use s2_pl_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_sky, filename_out
  logical :: scale = .true.
  character(len=S2_STRING_LEN) :: title = ''
  type(s2_sky) :: sky
  type(s2_pl) :: cl
  integer :: lmax=0, mmax=0
  logical :: new_lmax = .false., new_mmax = .false.
  logical :: dil1_status = .false., dil2_status = .false.
  logical :: norm_pres = .false.
  real(s2_sp) :: dil1, dil2

  ! Read command line options.
  call parse_options()

  ! Read full sky object.
  sky = s2_sky_init(filename_sky, S2_SKY_FILE_TYPE_SKY)

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

  ! Get cls from sky.
  cl = s2_sky_get_cl(sky)

  ! Produce plot.
  call s2_pl_plot(cl, trim(filename_out), scale, trim(title))

  ! Free memory.
  call s2_pl_free(cl)
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
            write(*,'(a)') 'Usage: s2_sky2plplot [-inp filename_sky]'
            write(*,'(a)') '                      [-out filename_out]'
            write(*,'(a)') '                      [-scale scale]'
            write(*,'(a)') '                      [-title title]'
            write(*,'(a)') '                      [-lmax lmax]'
            write(*,'(a)') '                      [-mmax mmax]'
            write(*,'(a)') '                      [-dil1 dil1]'
            write(*,'(a)') '                      [-dil2 dil2]'
            write(*,'(a)') '                      [-npres norm_pres]'
            stop
          
          case ('-inp')
            filename_sky = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-scale')
            read(arg,*) scale

          case ('-title')
            title = trim(arg)

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


end program s2_sky2plplot

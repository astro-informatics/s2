!------------------------------------------------------------------------------
! s2_pl_mod -- S2 library pl class
!
!! Provides functionality to support and manipulate a p(l) (l spectrum) 
!! function.
!    
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_pl_mod

  use s2_types_mod
  use s2_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_pl_init, &
    s2_pl_init_guassian, &
    s2_pl_free, &
    s2_pl_add, &
    s2_pl_conv, &
    s2_pl_power, &
    s2_pl_plot, &
    s2_pl_io_ascii_write, &
    s2_pl_io_fits_write, &
    s2_pl_get_init, &
    s2_pl_get_lmax, &
    s2_pl_get_spec, &
    s2_pl_get_spec_l


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  interface s2_pl_init
     module procedure &
       s2_pl_init_array, &
       s2_pl_init_array_const, &
       s2_pl_init_file_ascii, &
       s2_pl_init_file_fits, &
       s2_pl_init_copy, &
       s2_pl_init_copy_truncate
  end interface
  
  interface s2_pl_conv
     module procedure &
       s2_pl_conv_alm, &
       s2_pl_conv_pl
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_pl
     private
     logical :: init = .false.
     integer :: lmax = 0
     real(s2_sp), allocatable :: pl(:)
  end type s2_pl



  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_pl_init_array
    !
    !! Initiliase a pl from a array of data values.
    !!
    !! Variables:
    !!   - data: Data array to initialise pl spectrum.
    !!   - pl: Initialised pl.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_array(data) result(pl)

      real(s2_sp), intent(in) :: data(:)
      type(s2_pl) :: pl

      integer :: fail

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_array')
        return
      end if

      ! Set size.
      pl%lmax = size(data) - 1

      ! Allocate space.
      allocate(pl%pl(0:pl%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_array')
      end if

      ! Save spectrum.
      pl%pl = data

      ! Set initialised.
      pl%init = .true.

    end function s2_pl_init_array


    !--------------------------------------------------------------------------
    ! s2_pl_init_array_const
    !
    !! Initiliase a pl with a constant value over all l.
    !!
    !! Variables:
    !!   - const: Constant values used to  array to initialise pl spectrum.
    !!   - lmax: Lmax of the constructed pl spectrum.
    !!   - pl: Initialised pl.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 February 2005
    !
    ! Revisions:
    !   February 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_array_const(const, lmax) result(pl)

      real(s2_sp), intent(in) :: const
      integer, intent(in) :: lmax
      type(s2_pl) :: pl

      integer :: l, fail

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_array')
        return
      end if

      ! Set size.
      pl%lmax = lmax

      ! Allocate space.
      allocate(pl%pl(0:pl%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_array')
      end if

      ! Save spectrum.
      do l = 0,lmax
         pl%pl(l) = const
      end do

      ! Set initialised.
      pl%init = .true.

    end function s2_pl_init_array_const


    !--------------------------------------------------------------------------
    ! s2_pl_init_file_ascii
    !
    !! Initialise a pl from a ascii file.
    !!
    !! Variables:
    !!   - filename: Name of ascii file to read spectrum from.
    !!   - lmin: Minimum l to read (first value read starts at this index ).  
    !!     (Pl values below lmin are set to zero.)
    !!   - lmax: Pl lmax.
    !!   - ncomment: Number of initial comment lines in ascii file to ignore.
    !!   - [scale_in]: Logical to specify whether to scale the real in pl 
    !!     values.  If true the values are scaled by 2*pi/(l*(l+1)) 
    !!     (necessary to scale input cmb cl spectrums to get actual pl 
    !!     values).  If not present defaults to false.
    !!   - [line_nos_in]: Logical to specify whether each line in the file
    !!     contains a preceeding l count.  Defaults to true.
    !!   - pl: Initialised pl.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_file_ascii(filename, lmin, lmax, ncomment, &
      scale_in, line_nos_in) result(pl)
      
      character(len=*), intent(in) :: filename
      integer, intent(in) :: lmin, lmax, ncomment
      logical, intent(in), optional :: scale_in
      logical, intent(in), optional :: line_nos_in
      type(s2_pl) :: pl

      integer :: l, fileid = 10, fail
      real(s2_dp) :: ldum
      character(len=80) :: line
      logical, parameter :: DEFAULT_SCALE_STATUS = .false.
      logical :: scale
      logical :: line_nos = .true.

! real(s2_sp) :: plval_dp

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_file_ascii')
        return
      end if

      if(present(scale_in)) then
         scale = scale_in
      else
         scale = DEFAULT_SCALE_STATUS
      end if

      if(present(line_nos_in)) line_nos = line_nos_in

      ! Set size.
      pl%lmax = lmax

      ! Open file.
      open(fileid, file=filename, form='formatted', status='old')

      ! Ignore comment lines.
      do l = 1,ncomment
         read(fileid,*) line
      end do
      
      ! Allocate space.
      allocate(pl%pl(0:pl%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_file_ascii')
      end if

      ! Initialise to zero so pl values below lmin are zero.
      pl%pl = 0.0e0

      ! Read data.
      do l = lmin,lmax

         if(line_nos) then
            read(fileid,*) ldum, pl%pl(l)
         else
            read(fileid,*) pl%pl(l)
         end if

! read(fileid,*) ldum, line      
! write(*,*) trim(line)
! read(line,*) plval_dp
! pl%pl(l) = real(plval_dp, s2_sp)
! write(*,*) ldum, pl%pl(l)
      end do

      ! Close file.
      close(fileid)

      ! Scale if required.
      if(scale) then
         ! Ignore monopole to avoid divide by zero problem. (This will 
         ! usually be set to zero anyway.)
         do l = 2,lmax  
            pl%pl(l) = pl%pl(l) * 2.0e0*PI / real((l*(l+1)))
         end do
      end if

      ! Set initialised.
      pl%init = .true.

    end function s2_pl_init_file_ascii
   

    !--------------------------------------------------------------------------
    ! s2_pl_init_file_fits
    !
    !! Wrapper to initialise a pl data structure from a s2_pl fits file.
    !! The pl structure is read and initialised by the routine 
    !! s2_pl_io_fits_read.
    !!
    !! Variables:
    !!   - filename: Name of s2_pl fits file containing the pl data to 
    !!     be read.
    !!   - pl: Returned pl structure initialised with the data contained in
    !!     the input s2_pl fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_file_fits(filename) result(pl)

      character(len=*), intent(in) :: filename
      type(s2_pl) :: pl

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_file_fits')
        return
      end if

      ! Read s2_pl file.
      call s2_pl_io_fits_read(filename, pl)

      ! All status flags set in s2_pl_io_fits_read routine.
      
    end function s2_pl_init_file_fits


    !--------------------------------------------------------------------------
    ! s2_pl_init_copy
    !
    !! Initialse a pl from a copy of another pl.
    !!
    !! Variables:
    !!   - orig: Original pl to copy.
    !!   - copy: Copied pl.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_copy(orig) result(copy)

      type(s2_pl), intent(in) :: orig
      type(s2_pl) :: copy

      integer :: fail

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_copy')
        return
      end if

      ! Copy object attributes.

      copy%lmax = orig%lmax

      allocate(copy%pl(0:copy%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_copy')
      end if

      copy%pl = orig%pl
      copy%init = .true.
      
    end function s2_pl_init_copy
  

    !--------------------------------------------------------------------------
    ! s2_pl_init_copy_truncate
    !
    !! Initialse a pl from a copy of another pl but truncate to a lower lmax.
    !!
    !! Variables:
    !!   - orig: Original pl to copy.
    !!   - lmax: Lmax of copy pl to truncate to.
    !!   - copy: Copied pl.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 Januray 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_copy_truncate(orig, lmax) result(copy)

      type(s2_pl), intent(in) :: orig
      integer, intent(in) :: lmax
      type(s2_pl) :: copy

      integer :: fail

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_init_copy_truncate')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_copy_truncate')
        return
      end if

      ! Check new lmax less than lmax of orig.
      if(lmax < 0 .or. lmax > orig%lmax) then
        call s2_error(S2_ERROR_PL_SIZE_INVALID, &
          's2_pl_init_copy_truncate', &
          comment_add='New lmax out of range')
     end if

     ! Construct copy.

      copy%lmax = lmax

      allocate(copy%pl(0:copy%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_copy_truncate')
      end if

      copy%pl = orig%pl(0:copy%lmax)
      copy%init = .true.

    end function s2_pl_init_copy_truncate


    !--------------------------------------------------------------------------
    ! s2_pl_init_gaussian
    !
    !! Initiliase a pl with a Gaussian beam with specified full width half max.
    !!
    !! Notes:
    !!   - Gaussian beam defined as exp(-l(l+1)*sigma2), where 
    !!     sigma2 = fwhm**2.0e0 / (8.0e0 * log(2.0e0)).
    !!
    !! Variables:
    !!   - fwhm: Full width half maximum of Gaussian beam (in radians).
    !!   - lmax: Lmax of the constructed pl spectrum.
    !!   - pl: Initialised pl with Gaussian beam.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_init_guassian(fwhm, lmax) result(pl)

      real(s2_sp), intent(in) :: fwhm
      integer, intent(in) :: lmax
      type(s2_pl) :: pl
      
      real(s2_sp) :: sigma2
      integer :: l, fail

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_init_gaussian')
        return
      end if

      ! Set size.
      pl%lmax = lmax

      ! Allocate space.
      allocate(pl%pl(0:pl%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_gaussian')
      end if

      ! Compute variance from full width half max.
      sigma2 = fwhm**2.0e0 / (8.0e0 * log(2.0e0))

      ! Compute spectrum values.
      do l = 0,lmax
         pl%pl(l) = exp(-l*(l+1)*sigma2 / 2.0)
      end do

      ! Set initialised.
      pl%init = .true.

    end function s2_pl_init_guassian


    !--------------------------------------------------------------------------
    ! s2_pl_free
    !
    !! Free all data associated with an initialised pl and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - pl: The pl to free.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_free(pl)

      type(s2_pl), intent(inout) :: pl

      ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_free')
      end if

      if(allocated(pl%pl)) deallocate(pl%pl)
      pl%lmax = 0
      pl%init = .false.

    end subroutine s2_pl_free


    !--------------------------------------------------------------------------
    ! s2_pl_add
    !
    !! All two pl spectrums.
    !!
    !! Notes:
    !!   - The returned pl_sum object is initialised herein and must be freed 
    !!     by the calling routine.
    !!
    !! Variables:
    !!   - pl1: First pl spectrum to be added.
    !!   - pl2: Second pl spectru to be added.
    !!   - pl_sum: Summed spectrum.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_add(pl1, pl2) result(pl_sum)

      type(s2_pl), intent(in) :: pl1, pl2
      type(s2_pl) :: pl_sum

      integer :: fail

      ! Check objects initialised.
      if(.not. pl1%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_add')
      end if
       if(.not. pl2%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_add')
      end if

      ! Check sizes consistent.
      if(pl1%lmax /= pl2%lmax) then
        call s2_error(S2_ERROR_PL_SIZE_INVALID, 's2_pl_add', &
          comment_add='Pl spectra to be added have inconsistent size')
      end if

      ! Initialse pl_sum object and allocate space.
      pl_sum%lmax = pl1%lmax
      allocate(pl_sum%pl(0:pl_sum%lmax), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_init_copy')
      end if

      ! Sum spectra.
      pl_sum%pl = pl1%pl + pl2%pl

      ! Initialise status.
      pl_sum%init = .true.

    end function s2_pl_add


    !--------------------------------------------------------------------------
    ! s2_pl_conv_alm
    !
    !! Convolve an alm with the pl.  The original alm is overwritten with 
    !! the convolved alm on ouput.  If length of alm is shorter than pl then
    !! neglect larger pl values.  If length of alm is greater than pl then
    !! effectivly pad the pl with zeros, i.e. set high l alm values to zero.
    !!
    !! Variables:
    !!   - pl: Pl to convolve.
    !!   - alm: Alm to convolve.  Overwritten with new alm on output.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_conv_alm(pl, alm)

      type(s2_pl), intent(in) :: pl
      complex(s2_spc), intent(inout) :: alm(0:,0:)

      integer :: lmax_alm, l, m

      ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_conv')
      end if

      lmax_alm = size(alm,1) - 1

      do l = 0,min(pl%lmax, lmax_alm)
         do m = 0,l
            alm(l,m) = alm(l,m) * pl%pl(l)
         end do
      end do

      ! If alm extends to higher frequency than pl%lmax then effectively 
      ! pad the pl spectrum with zeros and convolve (i.e. set high frequency 
      ! alm terms to zero).
      if(lmax_alm > pl%lmax) then

         ! Display warning.
         call s2_error(S2_ERROR_PL_LMAX_LOW, 's2_pl_conv', &
           comment_add='Padding pl spectrum with zeros')

         do l = pl%lmax+1,lmax_alm
            do m = 0,l
               alm(l,m) = cmplx(0.0e0,0.0e0)
            end do
         end do
      end if

    end subroutine s2_pl_conv_alm


    !--------------------------------------------------------------------------
    ! s2_pl_conv_pl
    !
    !! Convolve two power spectra.  Pl1 is overwritten with the convolved 
    !! spectrum on output.  (Spectra must have same lmax.)
    !!
    !! Variables:
    !!   - pl1: First pl to convolve.  Overwritten with new spectrum on output.
    !!   - pl2: Second pl to convolve.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_conv_pl(pl1, pl2)

      type(s2_pl), intent(inout) :: pl1
      type(s2_pl), intent(in) :: pl2

      integer :: l

      ! Check objects initialised.
      if(.not. pl1%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_conv')
      end if
      if(.not. pl2%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_conv')
      end if

      ! Check pls have same lmax.
      if(pl1%lmax /= pl2%lmax) then
         call s2_error(S2_ERROR_PL_SIZE_INVALID, 's2_pl_conv', &
           comment_add='Spectra to convolve have inconsistent lmax')
         stop
      end if

      ! Perform convolution.
      do l = 0,pl1%lmax
         pl1%pl(l) = pl1%pl(l) * pl2%pl(l)
      end do

    end subroutine s2_pl_conv_pl


    !--------------------------------------------------------------------------
    ! s2_pl_power
    !
    !! Compute full power from a pl spectrum
    !!
    !! Variables:
    !!   - pl: Pl structure containing power spectrum to compute total power 
    !!     of.
    !!   - power: Full power computed.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 May 2006
    !
    ! Revisions:
    !   May 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_power(pl) result(power)

      type(s2_pl), intent(in) :: pl
      real(s2_sp) :: power

      integer :: l

      ! Check objects initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_power')
      end if

      ! Compute power.
      power = pl%pl(0)
      do l = 1,pl%lmax
         power = power + (2*l+1) * pl%pl(l)
      end do

    end function s2_pl_power


    !--------------------------------------------------------------------------
    ! s2_pl_plot
    !
    !! Plot the pl angular power spectrum as a postscript.
    !!
    !! Notes:
    !!   - l=0 term (dc term) is not plotted.  In almost all cases this will 
    !!     be zero anyway.
    !!
    !! Variables:
    !!   - pl: Pl structure containing power spectrum to plot.
    !!   - filename: Name of output postscript file (including .ps extension 
    !!     but excluding local /cps device flag).
    !!   - [scale_to_dl_in]: Logical to specify whether to scale the pl cl 
    !!     spectrum to a dl spectrum, where d(l) = c(l) * l * (l+1) / (2*pi).
    !!     (Default is true.)
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 February 2005
    !
    ! Revisions:
    !   February 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_plot(pl, filename, scale_to_dl_in, title)

      type(s2_pl), intent(in) :: pl
      character(len=*), intent(in) :: filename
      logical, intent(in), optional :: scale_to_dl_in
      character(len=*), intent(in), optional :: title

      character(len=S2_STRING_LEN) :: filename_local
      real(s2_sp), allocatable :: l_axis(:)
      real(s2_sp), allocatable :: pl_plot(:)
      integer :: l, fail
      logical :: scale_to_dl = .true.

      ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_plot')
      end if

      ! Create l index array to plot against.
      allocate(l_axis(1:pl%lmax), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_plot')
      end if
      do l = 1,pl%lmax
         l_axis(l) = real(l,s2_sp)
      end do
      
      ! Create copy of spectrum in case need to scale.
      allocate(pl_plot(1:pl%lmax), stat=fail)
      if(fail /= 0) then
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_pl_plot')
      end if
      pl_plot = pl%pl(1:pl%lmax)

      ! Scale to dl if requested.
      if(present(scale_to_dl_in)) scale_to_dl = scale_to_dl_in
      if(scale_to_dl) then
         do l = 1,pl%lmax
            pl_plot(l) = pl_plot(l) * l * (l+1) / (2.0e0*pi)
         end do
      end if

#ifndef NO_PGPLOT
      ! Plot spectrum.
      write(filename_local, '(a,a)') trim(filename), '/cps'
      call pgopen(trim(filename_local))
      call pgenv(0.0e0, real(pl%lmax,s2_sp), 0.0e0, 1.2e0*maxval(pl_plot), &
        0, 0)
      call pgline(pl%lmax, l_axis, pl_plot)

      ! Add labels.
      call pgmtxt('B', 2.0e0, 0.5e0, 0.5e0, '\fsl')
      if(scale_to_dl) then
         call pgmtxt('L', 2.0e0, 0.5e0, 0.5e0, &
           '\fsl\fn(\fsl\fn+1)C\d\fsl\fn\u/2\gp')
      else
         call pgmtxt('L', 2.0e0, 0.5e0, 0.5e0, 'C\u\fsl')
      end if
      if(present(title)) then 
         call pgmtxt('T', 2.0e0, 0.5e0, 0.5e0, trim(title))
      end if

      ! Close plot.
      call pgclos
#else
      call s2_error(S2_ERROR_PL_PLOT_FAIL, 's2_pl_plot', &
         comment_add='Not compiled with pgplot')
#endif

      ! Free memory.
      deallocate(l_axis)
      deallocate(pl_plot)
      
    end subroutine s2_pl_plot


    !--------------------------------------------------------------------------
    ! File IO routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_pl_io_ascii_write
    !
    !! Write a s2_pl object to an ascii file.
    !!
    !! Variables:
    !!   - filename: Name of the output fits file to write the pl structure
    !!     data to.
    !!   - pl: The pl structure containing the data to be written to the
    !!     output fits file.
    !!   - [dl_scale]: Logical specifying whether to scale to a dl spectrum.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_ascii_write(filename, pl, dl_scale)

      character(len=*), intent(in) :: filename
      type(s2_pl), intent(in) :: pl
      logical, intent(in), optional :: dl_scale

      integer :: fileid = 21, l
      logical :: scale = .false.

      if(present(dl_scale)) scale = dl_scale
      
      open(unit=fileid, file=filename, status='replace', action='write')

!write(fileid,*) 'cl=['

      do l = 0, pl%lmax

         if(scale) then
            write(fileid,'(i4,e20.10)') l, pl%pl(l) * l * (l+1) / (2.0e0*pi)
         else
            write(fileid,'(i4,e20.10)') l, pl%pl(l)
!            write(fileid,'(e20.10,a)') pl%pl(l), ';'
         end if

      end do

!write(fileid,*) '];'

      close(fileid)

    end subroutine s2_pl_io_ascii_write


    !--------------------------------------------------------------------------
    ! s2_pl_io_fits_write
    !
    !! Write a s2_pl object to a fits file.
    !!
    !! Variables:
    !!   - filename: Name of the output fits file to write the pl structure
    !!     data to.
    !!   - pl: The pl structure containing the data to be written to the
    !!     output fits file.
    !!   - [comment]: Optional comment string to be added to the output fits 
    !!     file header if present.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_fits_write(filename, pl, comment)

      character(len=*), intent(in) :: filename
      type(s2_pl), intent(in) :: pl
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      logical :: simple, extend, file_exists
      integer :: naxis
      integer :: naxes(1)
      integer :: tfields, nrows, varidat
      character(len=32) :: ttype(1), tform(1), tunit(1), extname
      integer :: frow, felem, colnum

      ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_io_fits_write')
      end if

      ! Define FITS parameters.

      bitpix=-32 ! Real single precision.
      status=0   ! Initialse error status to zero.

      ! Check if file already exists.
      call s2_pl_io_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call s2_error(S2_ERROR_PL_FILE_EXISTS, &
              's2_pl_io_fits_write')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.
      call ftpcom(unit, &
        '  Pl harmonic power spectrum data file created by sph-0.1',status)
      call ftpcom(unit, &
        '  Primary extension empty',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if
      call ftpkyj(unit,'LMAX', pl%lmax, &
        'max spherical harmonic l considered',status)

      ! Insert binary table extension for spectrum values.
      extname='CL'
      ttype(1)='CLVALS'
      tform(1)='1E'
      tunit(1)=''
      tfields=1
      nrows=pl%lmax+1
      varidat=0
      call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

      ! Write map values to binary table.
      frow=1
      felem=1
      colnum=1
      call ftpcle(unit,colnum,frow,felem,nrows,pl%pl(0:pl%lmax), status)

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2_pl_io_fits_error_check(status, .true.)

    end subroutine s2_pl_io_fits_write


    !--------------------------------------------------------------------------
    ! s2_pl_io_fits_read
    !
    !! Read a fits pl file and allocate a new pl structure with the data read.
    !!
    !! Variables:
    !!   - filename: Name of s2_pl file containing the pl spectrum data to
    !!     be read.
    !!   - pl: Returned pl structure initialised with the data contained in
    !!     the input s2_pl fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_fits_read(filename, pl)

      character(len=*), intent(in) :: filename
      type(s2_pl), intent(out) :: pl

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: ihdu, hdutype, naxis
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: colnum, frow, felem, nelem
      real(s2_sp) :: nullval

      integer :: fail

      ! Check object not already initialised.
      if(pl%init) then
        call s2_error(S2_ERROR_INIT, 's2_pl_io_fits_read')
        return
      end if

      ! Initialse error status to zero.
      status=0   

      ! Check if file already exists.
      call s2_pl_io_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call s2_error(S2_ERROR_PL_FILE_INVALID, &
           's2_pl_io_fits_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)

      ! Read primary header lmax.
      call ftgkyj(unit, 'LMAX', pl%lmax, comment, status)

      ! Check correct number of HDUs in input file.
      hdunum = 2  ! One for primary header and one for inary data table.
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
         call s2_error(S2_ERROR_PL_FILE_INVALID, &
           's2_pl_io_fits_read', &
           comment_add='Invalid number of extensions')
      end if

      ! Allocate space for values.
      allocate(pl%pl(0:pl%lmax), stat=fail)
      if(fail /= 0) then 
         call s2_error(S2_ERROR_MEM_ALLOC_FAIL, &
           's2_pl_io_fits_read')
      end if
      ! Initialise with zeros.
      pl%pl = 0.0e0
      
      ! Set hdu index ready to read from next extension.
      ihdu = 2 

      ! Move to next ihdu extension (i.e. data values extension).
      call ftmahd(unit, ihdu, hdutype, status)

      ! Check correct hdutype (i.e. binary table).
      if(hdutype /= 2) then
         call s2_error(S2_ERROR_PL_FILE_INVALID, &
           's2_pl_io_fits_read', &
           comment_add='Values not stored in binary table')
      end if

      ! Read header NAXIS2 and check same as pl%lmax + 1.
      call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
      if(naxis/=pl%lmax + 1) then
         call s2_error(S2_ERROR_PL_FILE_INVALID, &
           's2_pl_io_fits_read', &
           comment_add='Inconsistent number of data samples')
      end if

      ! Read map values from binary table.
      frow=1
      felem=1
      nelem=pl%lmax+1
      nullval = -999  ! Arbitrary since will stop and return error 
                      ! if null values detected.
      colnum=1      
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           pl%pl(0:pl%lmax),anynull,status)
      if(anynull) then
         call s2_error(S2_ERROR_PL_FILE_INVALID, &
           's2_pl_io_fits_read', &
           comment_add='Null data values contained in file')
      end if

      ! Set pl initialised status.
      pl%init = .true.

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2_pl_io_fits_error_check(status, .true.)

    end subroutine s2_pl_io_fits_read


    !--------------------------------------------------------------------------
    ! s2_pl_io_fits_error_check
    !
    !! Check if a fits error has occured and print error message.  Halt
    !! program execution if halt flag is set.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - status: Fits integer status code.
    !!   - halt: Logical to indicate whether to halt program execution if an 
    !!     error is detected.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_fits_error_check(status, halt)

      integer, intent(inout) :: status
      logical, intent(in) :: halt

      character(len=30) :: errtext
      character(len=80) :: errmessage

      !  Check if status is OK (no error); if so, simply return.
      if (status .le. 0) return

      ! The FTGERR subroutine returns a descriptive 30-character text 
      ! string that corresponds to the integer error status number.  
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      ! The FTGMSG subroutine retrieves the oldest message from
      ! the stack and shifts any remaining messages on the stack down one
      ! position.  FTGMSG is called repeatedly until a blank message is
      ! returned, which indicates that the stack is empty.  Each error message
      ! may be up to 80 characters in length. 
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write(*,*) trim(errmessage)
          call ftgmsg(errmessage)
      end do

      if(halt) stop

    end subroutine s2_pl_io_fits_error_check


    !--------------------------------------------------------------------------
    ! s2_pl_io_fits_exists
    !
    !! Check if a fits file exists.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to check existence of.
    !!   - status: Fits integer status code.
    !!   - exists: Logical indicating whether the fits file already exists.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_fits_exists(filename, status, exists)

      character(len=*), intent(in) :: filename
      integer, intent(inout) :: status
      logical, intent(out) :: exists

      integer :: unit, blocksize
      logical :: halt

      ! Simply return if status is already greater than zero.
      if (status .gt. 0) return

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      call ftopen(unit, filename, 1, blocksize, status)

      ! Check status of opening file.
      if(status == 0) then

        ! File was opened.  Close it and set exists flag accordingly.
        call ftclos(unit, status)
        exists = .true.

      else if (status == 104) then
        
        ! File does not exist.  Reset status and set exists flag accordingly.
         status = 0
         exists = .false.

      else

        ! Some other error occured while opening file.
        halt = .false.
        call s2_pl_io_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine s2_pl_io_fits_exists


    !--------------------------------------------------------------------------
    ! s2_pl_io_fits_del
    !
    !! Delete a fits file.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to detele.
    !!   - status: Fits integer status code.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_io_fits_del(filename, status)

      character(len=*), intent(in) :: filename
      integer, intent(inout) ::  status

      integer :: unit, blocksize

      ! Simply return if status is greater than zero.
      if (status .gt. 0)return

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Try to open the file, to see if it exists.
      call ftopen(unit,filename,1,blocksize,status)

      if(status .eq. 0) then
         ! File was opened;  so now delete it.
         call ftdelt(unit,status)
      else if(status .eq. 103) then
         ! File doesn't exist, so just reset status to zero and clear errors.
          status=0
          call ftcmsg
      else
         ! There was some other error opening the file; delete the file anyway.
         status=0
         call ftcmsg
         call ftdelt(unit,status)
      end if

      ! Free the unit number for later reuse.
      call ftfiou(unit, status)

    end subroutine s2_pl_io_fits_del


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_pl_get_init
    !
    !! Get init variable from the passed pl.
    !!
    !! Variables:
    !!   - pl: Pl object to get the variable of.
    !!   - init: Object init variable returned.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_get_init(pl) result(init)

      type(s2_pl), intent(in) :: pl
      logical :: init

      init = pl%init

    end function s2_pl_get_init


    !--------------------------------------------------------------------------
    ! s2_pl_get_lmax
    !
    !! Get lmax variable from the passed pl.
    !!
    !! Variables:
    !!   - pl: Pl object to get the variable of.
    !!   - lmax: Object lmax variable returned.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Winit variable from the passed pl.ritten by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_get_lmax(pl) result(lmax)

      type(s2_pl), intent(in) :: pl
      integer :: lmax

      ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_get_lmax')
      end if

      lmax = pl%lmax

    end function s2_pl_get_lmax
    

    !--------------------------------------------------------------------------
    ! s2_pl_get_spec
    !
    !! Get pl spectrum variable from the passed pl.
    !!
    !! Variables:
    !!   - pl: Pl object to get the variable of.
    !!   - pl_spec: Object pl spec variable returned.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_pl_get_spec(pl, pl_spec)

      type(s2_pl), intent(in) :: pl
      real(s2_sp), intent(out) :: pl_spec(:)

     ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_get_pl')
      end if

      if(size(pl_spec) /= pl%lmax + 1) then
         call s2_error(S2_ERROR_PL_SIZE_INVALID, &
           's2_pl_get_spec', comment_add='Inconsistent size for spectrums')
      end if

      pl_spec = pl%pl

    end subroutine s2_pl_get_spec


    !--------------------------------------------------------------------------
    ! s2_pl_get_spec_l
    !
    !! Get spectrum value for a given l from the pl.
    !!
    !! Variables:
    !!   - pl: Pl object to get the variable of.
    !!   - l: l value to get value of.
    !!   - pl_val: Spectrum value for specified l.
    !    
    !! @author J. D. McEwen 
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_pl_get_spec_l(pl, l) result(pl_val)

      type(s2_pl), intent(in) :: pl
      integer, intent(in) :: l
      real(s2_sp) :: pl_val

     ! Check object initialised.
      if(.not. pl%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_pl_get_pl')
      end if

      if(l > pl%lmax) then
         call s2_error(S2_ERROR_PL_SIZE_INVALID, &
           's2_pl_get_spec_l', comment_add='Requested index outside bound')
      end if

      pl_val = pl%pl(l)

    end function s2_pl_get_spec_l


end module s2_pl_mod

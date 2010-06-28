!------------------------------------------------------------------------------
! s2_error_mod  -- S2 library error class
! 
!! Functionality to handle errors that may occur in the s2 library.  Public
!! s2 error codes are defined, with corresponding private error comments and 
!! default halt execution status.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_error_mod

  use s2_types_mod, only: S2_STRING_LEN

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: s2_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: S2_ERROR_NUM = 34

  integer, public, parameter :: &
    S2_ERROR_NONE = 0, &
    S2_ERROR_INIT = 1, &
    S2_ERROR_NOT_INIT = 2, &
    S2_ERROR_INIT_FAIL = 3, &
    S2_ERROR_MEM_ALLOC_FAIL = 4, &
    S2_ERROR_ARTH = 5, &
    S2_ERROR_SKY_SIZE_WARNING = 6, &
    S2_ERROR_SKY_SIZE_INVALID = 7, &
    S2_ERROR_SKY_SIZE_NOT_DEF = 8, &
    S2_ERROR_SKY_POL_DEF = 9, &
    S2_ERROR_SKY_MAP_DEF = 10, &
    S2_ERROR_SKY_MAP_NOT_DEF = 11, &
    S2_ERROR_SKY_ALM_DEF = 12, &
    S2_ERROR_SKY_ALM_NOT_DEF = 13, &
    S2_ERROR_SKY_EXT_INVALID = 14, &
    S2_ERROR_SKY_PIX_INVALID = 15, &
    S2_ERROR_SKY_NON_CONFORM = 16, &
    S2_ERROR_SKY_PIX_DIFF = 17, &
    S2_ERROR_SKY_FTYPE_INVALID = 18, &
    S2_ERROR_SKY_FILE_EXISTS = 19, &
    S2_ERROR_SKY_FILE_INVALID = 20, &
    S2_ERROR_SKY_FOV_METHOD_INVALID = 21, &
    S2_ERROR_DISNT_BND_INVALID = 22, &
    S2_ERROR_PL_LMAX_LOW = 23, &
    S2_ERROR_PL_SIZE_INVALID = 24, &
    S2_ERROR_PL_FILE_EXISTS = 25, &
    S2_ERROR_PL_FILE_INVALID = 26, &
    S2_ERROR_PL_PLOT_FAIL = 27, &
    S2_ERROR_WNOISE_TYPE_INVALID = 28, &
    S2_ERROR_VECT_TYPE_INVALID = 29, &
    S2_ERROR_VECT_CART_DIM_INVALID = 30, &
    S2_ERROR_YLM_ARG_INVALID = 31, &
    S2_ERROR_PROJ_METHOD_INVALID = 32, &
    S2_ERROR_PROJ_FIELD_INVALID = 33



  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.
  !! Comment associated with each error type.
  character(len=S2_STRING_LEN), parameter :: &
    error_comment(S2_ERROR_NUM) = &
      (/ & 
      'No error                                                                 ', &
      'Attempt to initialise object that has already been initialised           ', &
      'Object not initialised                                                   ', &
      'Object initialisation failed                                             ', &
      'Memory allocation failed                                                 ', &
      'Arithmetic exception                                                     ', &
      'Warning: Sky sizes not in reccommended range (0<lmax<=3*nside, mmax=lmax)', &
      'Invalid sky sizes                                                        ', &
      'Sky size not defined                                                     ', &
      'Warning: Functionality to handle polarisation not yet incorporated       ', &
      'Warning: Map already defined                                             ', &
      'Map not defined when required                                            ', &
      'Warning: Alm already defined                                             ', &
      'Alm not defined when required                                            ', &
      'Fits file extension greater than nmaps                                   ', &
      'Invalid pixelisation scheme                                              ', &
      'Skies do not conform                                                     ', &
      'Warning: Skies in different pixelisation schemes                         ', &
      'Invalid function type                                                    ', &
      'Full sky (map and alm) fits file already exists                          ', &
      'Full sky (map and alm) fits file invalid                                 ', &
      'Field-of-view method type invalid                                        ', &
      'Warning: Invalid bounds when sampling uniform distribution               ', &
      'Warning: Alm greater lmax than pl spectrum                               ', &
      'Invalid pl spectrum sizes                                                ', &
      'Pl file aleady exists                                                    ', &
      'Pl file invalid                                                          ', &
      'Pl plot failed                                                           ', &
      'Invalid noise type                                                       ', &
      'Invalid vector coordinate type                                           ', &
      'Invalid number of dimensions for cartesian vector                        ', &
      'Invalid agruments                                                        ', &
      'Invalid projection method                                                ', &
      'Invalid projection field                                                 ' &
      /) 
  
  !! Default program halt status of each error type.
  logical, parameter :: &
    halt_default(S2_ERROR_NUM) = &
      (/ &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .false., &
      .true.,  &
      .true.,  &
      .false., &
      .false., &
      .true.,  &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .false., &
      .true.,  &   
      .true.,  &   
      .true.,  &   
      .true.,  &   
      .false., &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.  /)
  
  
  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_error
    !
    !! Display error message corresponding to error_code and halt program 
    !! execution if required.
    !!
    !! Variables:
    !!   - error_code: Integer error code.
    !!   - [procedure]: Procedure name where s2_error called from.  Displayed 
    !!     when error message printed to screen.
    !!   - [comment_add]: If present, additional comment to append to default 
    !!     error comment.
    !!   - [comment_out]: If present the error comment is copied to comment_out
    !!     on output.
    !!   - [halt_in]: If present overrides default halt value.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_error(error_code, procedure, comment_add, &
      comment_out, halt_in)

      integer, intent(in) :: error_code
      character(len=*), intent(in), optional :: procedure, comment_add
      character(len=*), intent(inout), optional :: comment_out
      logical, intent(in), optional :: halt_in

      logical :: halt
      character(len=*), parameter :: comment_prefix = 'S2_ERROR: '

      !---------------------------------------
      ! Display error message
      !---------------------------------------

      if(present(procedure)) then

        if(present(comment_add)) then
          write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            '''', &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            ''''
        end if
 
     else

        if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
            trim(error_comment(error_code+1)), &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
        end if

      end if

      ! Copy error comment if comment_out present.
      if(present(comment_out)) comment_out = error_comment(error_code+1)

      !---------------------------------------
      ! Halt program execution if required
      !---------------------------------------
      
      if( present(halt_in) ) then
        halt = halt_in
      else
        halt = halt_default(error_code+1)
      end if

      if( halt ) then
        write(*,'(a,a,a,a,a)') comment_prefix, &
          '  Halting program execution ', &
          'due to error ''', trim(error_comment(error_code+1)), ''''
        stop
      end if

    end subroutine s2_error


end module s2_error_mod

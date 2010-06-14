!------------------------------------------------------------------------------
! s2_types_mod -- S2 library types class
!
!! Definition of intrinsic types and constants used in the s2 library.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_types_mod

  implicit none

  private


  ! --------------------------------------
  ! Intrinsic type definitions
  ! --------------------------------------

  !! Type definition for single precision real.
  integer, public, parameter :: s2_sp  = SELECTED_REAL_KIND(5,30)
  !! Type definition for double precision real.
  integer, public, parameter :: s2_dp  = SELECTED_REAL_KIND(12,200)

  !! Type definition for single precisison complex.
  integer, public, parameter :: s2_spc = KIND((1.0_s2_sp, 1.0_s2_sp))
  !! Type definition for double precision complex.
  integer, public, parameter :: s2_dpc = KIND((1.0_s2_dp, 1.0_s2_dp))


  ! --------------------------------------
  ! Constants
  ! --------------------------------------

  !! String buffer length.
  integer, public, parameter :: S2_STRING_LEN = 256

  !! PI definition.
  real(s2_sp), public, parameter :: PI = 3.141592653589793238462643383279502884197


end module s2_types_mod

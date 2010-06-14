!------------------------------------------------------------------------------
! s2_vect_mod -- S2 library vect class
!
!! Provide functionality to support and rotate vectors in R^3, represented in
!! either cartesian or spherical coordinates.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2_vect_mod

  use s2_types_mod
  use s2_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------
  
  public :: &
    s2_vect_init, &
    s2_vect_free, &
    s2_vect_convert, &
    s2_vect_rad_to_arcmin, &
    s2_vect_arcmin_to_rad, &
    s2_vect_dot, &
    s2_vect_rotate, &
    s2_vect_set_unit, &
    s2_vect_get_init, &
    s2_vect_get_x, &
    s2_vect_get_r, &
    s2_vect_get_theta, &
    s2_vect_get_phi, &
    s2_vect_get_type


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface s2_vect_init
     module procedure &
       s2_vect_init_cart, &
       s2_vect_init_sph, &
       s2_vect_init_copy
  end interface
  

  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Cartesian coordinate status type.
  integer, public, parameter :: S2_VECT_TYPE_CART = 0
  !! Spherical coordinate status type.
  integer, public, parameter :: S2_VECT_TYPE_S2  = 1
  !! Dimension of 3d cartesian vector.
  integer, parameter :: S2_VECT_CART_DIM = 3


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_vect
    private
    logical :: init = .false.
    real(s2_sp) :: x(S2_VECT_CART_DIM) = 0.0e0
    real(s2_sp) :: r = 0.0e0
    real(s2_sp) :: theta = 0.0e0
    real(s2_sp) :: phi = 0.e0
    integer :: type = S2_VECT_TYPE_CART
  end type s2_vect


  !----------------------------------------------------------------------------

  contains 


    !--------------------------------------------------------------------------
    ! s2_vect_init_cart
    !
    !! Initialise a vector in cartesian coordinates.
    !!
    !! Variables:
    !!   - x: array of 3 cartesian values (x,y,z).
    !!   - vect: Initialised vector.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_init_cart(x) result(vect)

      real(s2_sp), intent(in) :: x(:)
      type(s2_vect) :: vect

      ! Check object not already initialised.
      if(vect%init) then
        call s2_error(S2_ERROR_INIT, 's2_vect_init_cart')
        return
      end if

      if(size(x) /= S2_VECT_CART_DIM) then
         call s2_error(S2_ERROR_VECT_CART_DIM_INVALID, 's2_vect_init_cart')
      end if

      vect%type = S2_VECT_TYPE_CART
      vect%x = x
      vect%init = .true.

    end function s2_vect_init_cart


    !--------------------------------------------------------------------------
    ! s2_vect_init_sph
    !
    !! Initialise a vector in spherical coordinates.
    !!
    !! Variables:
    !!   - r: R spherical coordinate.
    !!   - theta: Theta spherical coordinate.
    !!   - phi: Phi spherical coordinate.
    !!   - vect: Initialised vector.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_init_sph(r, theta, phi) result(vect)

      real(s2_sp), intent(in) :: r, theta, phi
      type(s2_vect) :: vect

      ! Check object not already initialised.
      if(vect%init) then
        call s2_error(S2_ERROR_INIT, 's2_vect_init_sph')
        return
      end if

      vect%type = S2_VECT_TYPE_S2
      vect%r = r
      vect%theta = theta
      vect%phi = phi
      vect%init = .true.

    end function s2_vect_init_sph


    !--------------------------------------------------------------------------
    ! s2_vect_init_copy
    !    
    !! Initialise a vector by making a copy of another vector.
    !!
    !! Variables:
    !!   - orig: Original vector to copy.
    !!   - copy: The copy vector initialised.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_init_copy(orig) result(copy)

      type(s2_vect), intent(in) :: orig
      type(s2_vect) :: copy

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2_error(S2_ERROR_INIT, 's2_vect_init_copy')
        return
      end if
      
      ! Check orig object initialised.
      if(.not. orig%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_init_copy')
        return
      end if

      ! Copy object attributes.
      copy%x = orig%x
      copy%r = orig%r
      copy%theta = orig%theta
      copy%phi = orig%phi
      copy%type = orig%type
      
      copy%init = .true.
      
    end function s2_vect_init_copy


    !--------------------------------------------------------------------------
    ! s2_vect_free
    !
    !! Reset vector values.  (No memory to be freed.)
    !!
    !! Variables:
    !!   - vect: Vector to free.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_free(vect)

      type(s2_vect), intent(inout) :: vect

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_free')
        return
      end if

      vect%init = .false.
      vect%x = 0.0e0
      vect%r = 0.0e0
      vect%theta = 0.0e0
      vect%phi = 0.0e0
      vect%type = S2_VECT_TYPE_CART

    end subroutine s2_vect_free


    !--------------------------------------------------------------------------
    ! s2_vect_convert
    ! 
    !! Convert a vector to the coordinates system specified by type.  If the
    !! vector is already in the soordinate system then do nothing.
    !!
    !! Variables:
    !!   - vect: Vector to convert coordinate type.
    !!   - type: Coordinate type to convert to.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_convert(vect, type)
      
      type(s2_vect), intent(inout) :: vect
      integer, intent(in) :: type

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_convert')
        return
      end if
      
      ! If vector already in current coordinates then do nothing.
      if(vect%type == type) return

      ! Convert vector.
      if(type == S2_VECT_TYPE_CART) then
         call s2_vect_compute_cart(vect)
         vect%type = type
      else if(type == S2_VECT_TYPE_S2) then
         call s2_vect_compute_sph(vect)
         vect%type = type
      else
         call s2_error(S2_ERROR_VECT_TYPE_INVALID)
      end if

    end subroutine s2_vect_convert


    !--------------------------------------------------------------------------
    ! s2_vect_compute_sph
    !
    !! Compute spherical coordinates of a vector from cartesian coordinates.
    !!
    !! Variables:
    !!   - vect: Vector to compute spherical coordinates of.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_compute_sph(vect)

      type(s2_vect), intent(inout) :: vect

      real(s2_sp), parameter :: tol = 1.0d-5
    
      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_compute_sph')
        return
      end if

      vect%theta = 0.0e0
      vect%phi = 0.0e0
      vect%r = sqrt(vect%x(1)**2 + vect%x(2)**2 + vect%x(3)**2)

      ! If r = 0, then theta and phi arbitrary.
      ! Leave set to 0.
      if(abs(vect%r - tol) > 0) then

        vect%theta = atan2(sqrt(vect%x(1)**2+vect%x(2)**2), vect%x(3))

        ! If theta = 0, the phi arbitrary.
        ! Leave set to 0.
        if(abs(vect%theta - tol) > 0) then
          vect%phi = atan2(vect%x(2), vect%x(1))
          vect%phi = mod(vect%phi, 2*PI)
        end if

      end if
   
    end subroutine s2_vect_compute_sph


    !--------------------------------------------------------------------------
    ! s2_vect_compute_cart
    !
    !! Compute cartesian coordinates of a vector from spherical coordinates.
    !!
    !! Variables:
    !!   - vect: Vector to compute cartesian coordinates of.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_compute_cart(vect)

      type(s2_vect), intent(inout) :: vect

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_compute_cart')
        return
      end if

      vect%x(1) = vect%r * sin(vect%theta) * cos(vect%phi)
      vect%x(2) = vect%r * sin(vect%theta) * sin(vect%phi)
      vect%x(3) = vect%r * cos(vect%theta)

    end subroutine s2_vect_compute_cart


    !--------------------------------------------------------------------------
    ! s2_vect_rad_to_arcmin
    !
    !! Convert an angle from radians to arcmins.
    !!
    !! Variables:
    !!   - theta_rad: Angle in radians to be converted.
    !!   - theta_arcmin: Equivalent angle in arcmins.
    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_rad_to_arcmin(theta_rad) result(theta_arcmin)

      real(s2_sp), intent(in) :: theta_rad
      real(s2_sp) :: theta_arcmin

      theta_arcmin = (180.0e0 * 60.0e0 / pi) * theta_rad

    end function s2_vect_rad_to_arcmin


    !--------------------------------------------------------------------------
    ! s2_vect_arcmin_to_rad
    !
    !! Convert an angle from arcmins to radians.
    !!
    !! Variables:
    !!   - theta_arcmin: Angle in arcmins to be converted.
    !!   - theta_rad: Equivalent angle in radians.

    !
    !! @author J. D. McEwen
    !! @version 0.1 November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_arcmin_to_rad(theta_arcmin) result(theta_rad)

      real(s2_sp), intent(in) :: theta_arcmin
      real(s2_sp) :: theta_rad

      theta_rad = (pi / (60.0e0 * 180.0e0)) * theta_arcmin

    end function s2_vect_arcmin_to_rad


    !--------------------------------------------------------------------------
    ! s2_vect_dot
    !
    !! Compute the dot product of two vectors.
    !!
    !! Variables:
    !!   - a: First vector of dot product.
    !!   - b: Second vector of dot product.
    !!   - dot: Dot product of a and b computed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2006
    !
    ! Revisions:
    !   April 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_dot(a, b) result(dot)

      type(s2_vect), intent(inout) :: a   ! inout so can compute cart if nec,
      type(s2_vect), intent(inout) :: b
      real(s2_sp) :: dot

      real(s2_sp) :: x_a(S2_VECT_CART_DIM)
      real(s2_sp) :: x_b(S2_VECT_CART_DIM)

      ! Check both vectors initialised.
      if(.not. a%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_dot')
        return
      end if
      if(.not. b%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_dot')
        return
      end if

      ! Get cartesian coordinates (without changing representating
      ! of a and b).
      x_a = s2_vect_get_x(a)
      x_b = s2_vect_get_x(b)

      dot =  x_a(1)*x_b(1) + x_a(2)*x_b(2) + x_a(3)*x_b(3)

    end function s2_vect_dot


    !--------------------------------------------------------------------------
    ! s2_vect_rotate
    !
    !! Rotate a vector by the specified Euler angles.
    !!
    !! Variables:
    !!   - vect: The vector to be rotated on input and the rotated vector on 
    !!     output.
    !!   - alpha: Alpha Euler angle.
    !!   - beta: Beta Euler angle.
    !!   - gamma: Gamma Euler angle.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotate(vect, alpha, beta, gamma)

      type(s2_vect), intent(inout) :: vect
      real(s2_sp), intent(in) :: alpha, beta, gamma

      real(s2_sp) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      integer :: type_orig

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_rotate')
        return
      end if

      type_orig = vect%type

      ! Ensure vector in cartesian form.
      call s2_vect_convert(vect, S2_VECT_TYPE_CART)
      
      ! Generate rotation matrix.
      call s2_vect_rotmat_point(R, alpha, beta, gamma)
!     call s2_vect_rotmat_axis(R, alpha, beta, gamma)

      ! Rotations appear to differ.  Different definitions.
      !write(*,*) 'R: ', R(1,:)     
      !write(*,*) 'R: ', R(2,:)     
      !write(*,*) 'R: ', R(3,:)     

      ! This routine will be marginally faster.
!      call s2_vect_rotmat_direct(R, alpha, beta, gamma)
      
      !write(*,*) 'R: ', R(1,:)     
      !write(*,*) 'R: ', R(2,:)     
      !write(*,*) 'R: ', R(3,:)     
      !write(*,*)

      ! Rotate vector.
      vect%x = matmul(R, vect%x)

      ! Return vector to original coordinate type.
      call s2_vect_convert(vect, type_orig)

    end subroutine s2_vect_rotate
   


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_direct
    !
    !! Generate a rotation matrix to rotate a vector (in cartesian coordinates)
    !! by the speficied Euler angles.  The rotation matrix is coded directly 
    !! rather than multiplying rotation matrices that rotate about various 
    !! axis.
    !
    ! Notes:
    !   - ** Definition differs to that of rotmat.  Investigate further. **
    !!
    !! Variables:
    !!   - R: Rotation matrix generated.
    !!   - alpha: Alpha Euler angle.
    !!   - beta: Beta Euler angle.
    !!   - gamma: Gamma Euler angle.
    !
    !! @author D. J. Mortlock
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Sourced from general library written by 
    !                    Daniel Mortlock
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_direct(R, alpha, beta, gamma)

      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: alpha, beta, gamma

      real(s2_sp) :: sinalpha, sinbeta, singamma, cosalpha, cosbeta, cosgamma

      sinalpha = sin(alpha)
      sinbeta = sin(beta)
      singamma = sin(gamma)
      cosalpha = cos(alpha)
      cosbeta = cos(beta)
      cosgamma = cos(gamma)

      R(1, 1) = cosalpha * cosbeta * cosgamma - sinalpha * singamma
      R(1, 2) = - cosalpha * cosbeta * singamma - sinalpha * cosgamma
      R(1, 3) = cosalpha * sinbeta

      R(2, 1) = sinalpha * cosbeta * cosgamma + cosalpha * singamma
      R(2, 2) = - sinalpha * cosbeta * singamma + cosalpha * cosgamma
      R(2, 3) = sinalpha * sinbeta

      R(3, 1) = - sinbeta * cosgamma
      R(3, 2) = sinbeta * singamma
      R(3, 3) = cosbeta

    end subroutine s2_vect_rotmat_direct



    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_point
    !
    !! Generate a rotation matrix to rotate a vector (in cartesian coordinates)
    !! by the specified Euler angles.  The point is rotated in the right hand 
    !! sense (axis rotated in left hand sense).    
    !!
    !! Notes:
    !!   - The rotation matrix is generated from:
    !!       R = Rz(gamma) * Ry(beta) * Rz(alpha)
    !!
    !! Variables:
    !!   - R: Rotation matrix generated.
    !!   - alpha: Alpha Euler angle.
    !!   - beta: Beta Euler angle.
    !!   - gamma: Gamma Euler angle.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_point(R, alpha, beta, gamma)
   
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: alpha, beta, gamma

      real(s2_sp), dimension(S2_VECT_CART_DIM, S2_VECT_CART_DIM) :: &
        Rz_gamma, Ry_beta, Rz_alpha

      call s2_vect_rotmat_z_point(Rz_gamma, gamma)
      call s2_vect_rotmat_y_point(Ry_beta, beta)
      call s2_vect_rotmat_z_point(Rz_alpha, alpha)

      R = matmul(matmul(Rz_alpha, Ry_beta), Rz_gamma)

    end subroutine s2_vect_rotmat_point


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_x_point
    !
    !! Generate a rotation matrix to rotate about the x axis by the specified 
    !! angle.  The point is rotated in the right hand sense (axis 
    !! rotated in left hand sense). 
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about x axis.
    !!   - ang: Angle to rotate about x axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_x_point(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0
      R(1,1) = 1.0
      R(2,2) = cos_ang
      R(2,3) = -sin_ang 
      R(3,2) = sin_ang
      R(3,3) = cos_ang

    end subroutine s2_vect_rotmat_x_point


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_y_point
    !
    !! Generate a rotation matrix to rotate about the y axis by the specified 
    !! angle.   The point is rotated in the right hand sense (axis 
    !! rotated in left hand sense). 
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about y axis.
    !!   - ang: Angle to rotate about y axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_y_point(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0   
      R(1,1) = cos_ang
      R(1,3) = sin_ang
      R(2,2) = 1.0
      R(3,1) = -sin_ang
      R(3,3) = cos_ang
      
    end subroutine s2_vect_rotmat_y_point


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_z_point
    !
    !! Generate a rotation matrix to rotate about the z axis by the specified 
    !! angle.  The point is rotated in the right hand sense (axis 
    !! rotated in left hand sense). 
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about z axis.
    !!   - ang: Angle to rotate about z axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_z_point(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0
      R(1,1) = cos_ang
      R(1,2) = -sin_ang
      R(2,1) = sin_ang
      R(2,2) = cos_ang
      R(3,3) = 1.0e0

    end subroutine s2_vect_rotmat_z_point


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_axis
    !
    !! Generate a rotation matrix to rotate a vector (in cartesian coordinates)
    !! by the specified Euler angles.  The axis is rotated in the right hand 
    !! sense (point rotated in left hand sense).    
    !!
    !! Notes:
    !!   - The rotation matrix is generated from:
    !!       R = Rz(gamma) * Ry(beta) * Rz(alpha)
    !!
    !! Variables:
    !!   - R: Rotation matrix generated.
    !!   - alpha: Alpha Euler angle.
    !!   - beta: Beta Euler angle.
    !!   - gamma: Gamma Euler angle.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_axis(R, alpha, beta, gamma)
   
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: alpha, beta, gamma

      real(s2_sp), dimension(S2_VECT_CART_DIM, S2_VECT_CART_DIM) :: &
        Rz_gamma, Ry_beta, Rz_alpha

      call s2_vect_rotmat_z_axis(Rz_gamma, gamma)
      call s2_vect_rotmat_y_axis(Ry_beta, beta)
      call s2_vect_rotmat_z_axis(Rz_alpha, alpha)

      R = matmul(matmul(Rz_gamma, Ry_beta), Rz_alpha)

    end subroutine s2_vect_rotmat_axis


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_x_axis
    !
    !! Generate a rotation matrix to rotate about the x axis by the specified 
    !! angle.  The axis is rotated in the right hand sense (point
    !! rotated in left hand sense). 
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about x axis.
    !!   - ang: Angle to rotate about x axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_x_axis(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0
      R(1,1) = 1.0
      R(2,2) = cos_ang
      R(2,3) = sin_ang 
      R(3,2) = -sin_ang
      R(3,3) = cos_ang

    end subroutine s2_vect_rotmat_x_axis


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_y_axis
    !
    !! Generate a rotation matrix to rotate about the y axis by the specified 
    !! angle.  The axis is rotated in the right hand sense (point
    !! rotated in left hand sense).  
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about y axis.
    !!   - ang: Angle to rotate about y axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_y_axis(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0   
      R(1,1) = cos_ang
      R(1,3) = -sin_ang
      R(2,2) = 1.0
      R(3,1) = sin_ang
      R(3,3) = cos_ang
      
    end subroutine s2_vect_rotmat_y_axis


    !--------------------------------------------------------------------------
    ! s2_vect_rotmat_z_axis
    !
    !! Generate a rotation matrix to rotate about the z axis by the specified 
    !! angle.   The axis is rotated in the right hand sense (point
    !! rotated in left hand sense).  
    !!
    !! Variables:
    !!   - R: Rotation matrix generated to rotate about z axis.
    !!   - ang: Angle to rotate about z axis by.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_rotmat_z_axis(R, ang)
      
      real(s2_sp), intent(out) :: R(S2_VECT_CART_DIM, S2_VECT_CART_DIM)
      real(s2_sp), intent(in) :: ang

      real(s2_sp) :: sin_ang, cos_ang

      sin_ang = sin(ang)
      cos_ang = cos(ang)

      R = 0.0e0
      R(1,1) = cos_ang
      R(1,2) = sin_ang
      R(2,1) = -sin_ang
      R(2,2) = cos_ang
      R(3,3) = 1.0e0

    end subroutine s2_vect_rotmat_z_axis


    !--------------------------------------------------------------------------
    ! s2_vect_set_unit
    !
    !! Scale the vector magnitude so it has unit length (orientation
    !! maintained).
    !!
    !! Variables:
    !!   - vect: Vector to scale to unit vector.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_vect_set_unit(vect)
      
      type(s2_vect), intent(inout) :: vect
      
      integer :: type_orig

      type_orig = vect%type

      ! Ensure vector in spherical coordinate form.
      call s2_vect_convert(vect, S2_VECT_TYPE_S2)
      
      ! Convert to unit vector.
      vect%r = 1.0e0

      ! Return vector to original coordinate type.
      call s2_vect_convert(vect, type_orig)

    end subroutine s2_vect_set_unit


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2_vect_get_init
    !
    !! Get init variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_init(vect) result(init)
      
      type(s2_vect), intent(in) :: vect
      logical :: init

      init = vect%init

    end function s2_vect_get_init


    !--------------------------------------------------------------------------
    ! s2_vect_get_x
    !
    !! Get x variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - x: Object x coordinate variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_x(vect) result(x)
      
      type(s2_vect), intent(inout) :: vect
      real(s2_sp) :: x(S2_VECT_CART_DIM)

      integer :: type_orig

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_get_x')
      end if

      ! Ensure in cartesian form.
      type_orig = vect%type
      call s2_vect_convert(vect, S2_VECT_TYPE_CART)

      x = vect%x

      ! Convert back to original type.
      call s2_vect_convert(vect, type_orig)

    end function s2_vect_get_x


    !--------------------------------------------------------------------------
    ! s2_vect_get_r
    !
    !! Get r variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - r: Object r coordinate variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_r(vect) result(r)
      
      type(s2_vect), intent(inout) :: vect
      real(s2_sp) :: r

      integer :: type_orig

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_get_r')
      end if

      ! Ensure in spherical coordinate form.
      type_orig = vect%type
      call s2_vect_convert(vect, S2_VECT_TYPE_S2)

      r = vect%r

      ! Convert back to original type.
      call s2_vect_convert(vect, type_orig)

    end function s2_vect_get_r


    !--------------------------------------------------------------------------
    ! s2_vect_get_theta
    !
    !! Get theta variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - theta: Object theta coordinate variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_theta(vect) result(theta)
      
      type(s2_vect), intent(inout) :: vect
      real(s2_sp) :: theta

      integer :: type_orig

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_get_theta')
      end if

      ! Ensure in spherical coordinate form.
      type_orig = vect%type
      call s2_vect_convert(vect, S2_VECT_TYPE_S2)

      theta = vect%theta

      ! Convert back to original type.
      call s2_vect_convert(vect, type_orig)

    end function s2_vect_get_theta


    !--------------------------------------------------------------------------
    ! s2_vect_get_phi
    !
    !! Get phi variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - phi: Object phi coordinate variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_phi(vect) result(phi)
      
      type(s2_vect), intent(inout) :: vect
      real(s2_sp) :: phi
 
      integer :: type_orig

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_get_phi')
      end if

      ! Ensure in spherical coordinate form.
      type_orig = vect%type
      call s2_vect_convert(vect, S2_VECT_TYPE_S2)

      phi = vect%phi

      ! Convert back to original type.
      call s2_vect_convert(vect, type_orig)

    end function s2_vect_get_phi


    !--------------------------------------------------------------------------
    ! s2_vect_get_type
    !
    !! Get type variable from the passed vect.
    !!
    !! Variables:
    !!   - vect: Vect object to get the variable of.
    !!   - type: Object type coordinate variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_vect_get_type(vect) result(type)
      
      type(s2_vect), intent(in) :: vect
      integer :: type

      ! Check object initialised.
      if(.not. vect%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_vect_get_type')
      end if

      type = vect%type

    end function s2_vect_get_type


end module s2_vect_mod

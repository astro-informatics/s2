!------------------------------------------------------------------------------
! s2_graph_mod -- S2 library graph class
!
!! Provides functionality to create a graph representation of a (potentially
!! masked) sky.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk) and A. Stradis
!
! Revisions:
!   December 2011 - Written by Jason McEwen and Athamos Stradis
!------------------------------------------------------------------------------

module s2_graph_mod

  use s2_types_mod, only: s2_sp, s2_spc, s2_dp, s2_dpc, pi
  use s2_error_mod
  use s2_vect_mod
  use s2_sky_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2_graph_init, &
    s2_graph_free


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  ! No interfaces.


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! No global variables.
  integer, public, parameter :: CRAXY_NUMBER = 3

  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_graph
    private
    logical :: init = .false.
    integer :: nside = 0
    integer :: pix_scheme = S2_SKY_NEST
    integer :: nvertices = 0
    integer, allocatable :: ivertices(:)
    type(s2_vect), allocatable :: vertices(:)
    real(s2_dp), allocatable :: adj(:,:)
    real(s2_dp), allocatable :: deg(:)
    type(s2_sky) :: mask
    logical :: mask_status = .false.
  end type s2_graph


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_graph_init
    !
    !! Initialie a graph object.
    !!
    !! ...
    !--------------------------------------------------------------------------

    function s2_graph_init(nside, mask) result(graph)

      use pix_tools, only: nside2npix

      integer, intent(in) :: nside
      type(s2_sky), intent(inout), optional :: mask
      type(s2_graph) :: graph

      integer :: npix, ipix, nneighbours, ineighbour
      integer :: neighbours(0:7)      
      real(s2_dp) :: theta, phi, theta_n, phi_n, dot, ang
      type(s2_vect) :: vec, vec_n
      integer :: fail = 0

      ! Check object not already initialised.
      if(graph%init) then
        call s2_error(S2_ERROR_INIT, 's2_graph_init')
        return
      end if

      ! Check mask consistent with passed nside and pix_scheme.
      if(present(mask)) graph%mask_status = .true.
      if(graph%mask_status) then

         ! Check mask initialised.         
         if(.not. s2_sky_get_init(mask)) then
            call s2_error(S2_ERROR_NOT_INIT, 's2_graph_init')
         end if
         
         ! Check mask has map defined.
         if(.not. s2_sky_get_init(mask)) then
            call s2_error(S2_ERROR_NOT_INIT, 's2_graph_init')
         end if
         if(.not. s2_sky_get_map_status(mask)) then
            call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_graph_init')
         end if
         
         ! Check resolution.
         if(nside /= s2_sky_get_nside(mask)) then
            call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_graph_init', &
                 comment_add='Mask has invalid nside.')
         end if

         ! Check mask in nested pixelisation scheme.
         call s2_sky_map_convert(mask, S2_SKY_NEST)
     
      end if
      npix = nside2npix(nside)

      ! Allocate space for graph object.
      allocate(graph%ivertices(0:npix-1), stat=fail)
      allocate(graph%vertices(0:npix-1), stat=fail)
      allocate(graph%adj(0:npix-1,0:npix-1), stat=fail)
      allocate(graph%deg(0:npix-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_graph_init')
      end if      

      ! Initialise graph attributes.
      graph%nside = nside
      graph%nvertices = npix

      ! Compute adjacency matrix.
      do ipix = 0, npix-1

         ! Compute position of current pixel and save.
         call pix2ang_nest(nside, ipix, theta, phi)
         vec = s2_vect_init(1.0, real(theta,s2_sp), real(phi,s2_sp))
         graph%ivertices(ipix) = ipix
         graph%vertices(ipix) = s2_vect_init(vec)

         ! Find neighbours.
         call neighbour_nest(nside, ipix, neighbours(0:7), nneighbours)

         ! Set adjacency matrix entries for all neighbours.
         do ineighbour = 0, nneighbours-1

            ! Compute position of neighbour.
            call pix2ang_nest(nside, neighbours(ineighbour), theta_n, phi_n)
            vec_n = s2_vect_init(1.0, real(theta_n,s2_sp), real(phi_n,s2_sp))

            ! Compute anglar separation between pixels.
            dot = s2_vect_dot(vec, vec_n)         
            if (dot > 1d0) dot = 1d0  ! Remove numerical noise that could cause acos to fail.
            if (dot < -1d0) dot = -1d0 
            ang = acos(dot)
            call s2_vect_free(vec_n)

            ! Set adjacency matrix value.
            graph%adj(ipix, neighbours(ineighbour)) = dot
         end do

         ! Compute diagonal degree matrix components.
         graph%deg(ipix) = sum(graph%adj(ipix,0:npix-1))
         
         call s2_vect_free(vec)

      end do

      ! Set as initialised.
      graph%init = .true.

    end function s2_graph_init


    !--------------------------------------------------------------------------
    ! s2_graph_free
    !
    !! Free all data associated with an initialised graph object and reset all
    !! other attributes.
    !!
    !! ...
    !--------------------------------------------------------------------------

    subroutine s2_graph_free(graph)

      type(s2_graph), intent(inout) :: graph

      ! Check object initialised.
      if(.not. graph%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_graph_free')
      end if 

      ! Free space.
      if(allocated(graph%ivertices)) deallocate(graph%ivertices)
      if(allocated(graph%vertices)) deallocate(graph%vertices)
      if(allocated(graph%adj)) deallocate(graph%adj)
      if(allocated(graph%deg)) deallocate(graph%deg)
      if(graph%mask_status) call s2_sky_free(graph%mask)
      
      ! Reset attributes.
      graph%nside = 0
      graph%pix_scheme = S2_SKY_NEST
      graph%nvertices = 0
      graph%mask_status = .false.

      graph%init = .false.

    end subroutine s2_graph_free


end module s2_graph_mod

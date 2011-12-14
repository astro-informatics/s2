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


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2_graph
    private
    logical :: init = .false.
    integer :: nside = 0
    integer :: pix_scheme = S2_SKY_RING
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

    function s2_graph_init(nside, pix_scheme, mask) result(graph)

      use pix_tools, only: nside2npix

      integer, intent(in) :: nside
      integer, intent(in) :: pix_scheme
      type(s2_sky), intent(in), optional :: mask
      type(s2_graph) :: graph

      integer :: npix, ipix

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

         ! Check pixelisation scheme.
         if(pix_scheme /= s2_sky_get_pix_scheme(mask)) then
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_graph_init', &
                 comment_add='Mask has invalid pix_scheme.')
         end if    
     
      end if
      npix = nside2npix(nside)




      ! Construct graph adjacency and degree matrices...
      
      ! Use the healpix function neighbours_nest to find adjacent pixel positions.







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
      graph%pix_scheme = S2_SKY_RING
      graph%nvertices = 0
      graph%mask_status = .false.

      graph%init = .false.

    end subroutine s2_graph_free


end module s2_graph_mod

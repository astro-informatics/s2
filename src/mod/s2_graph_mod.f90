!------------------------------------------------------------------------------
! s2_graph_mod -- S2 library graph class
!
!! Provides functionality to create a graph representation of a (potentially
!! masked) sky.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   September 2010 - Written by Jason McEwen 
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
    real(s2_dp) :: support_radius = 0d0
    integer :: nvertices = 0
    integer, allocatable :: ivertices(:)
    type(s2_vect), allocatable :: vertices(:)
    real(s2_dp), allocatable :: adj_matrix(:,:)
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
    !! Variables:
    !!   - 
    !!   -
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   September 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_graph_init(nside, pix_scheme, support_radius, mask) result(graph)

      use pix_tools, only: nside2npix

      integer, intent(in) :: nside
      integer, intent(in) :: pix_scheme
      real(s2_dp), intent(in) :: support_radius
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
         ! Check resolution.
         if(nside /= s2_sky_get_nsidex(mask)) then
            call s2_error(S2_ERROR_SKY_SIZE_INVALID, 's2_graph_init', &
                 comment_add='Mask has invalid nside.')
            return
         end if
         ! Check pixelisation scheme.
         if(pix_scheme /= s2_sky_get_pix_scheme(mask)) then
            call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2_graph_init', &
                 comment_add='Mask has invalid pix_scheme.')
            return
         end if         
      end if
      npix = nside2npix(nside)

      ! Count vertices.
      if(graph%mask_status) then         
         ! Get mask map values.

         ! Threshold all values above TOL to one, below to ZERO.


         ! Sum resultant values.



         ! Round to nearest integer.
         
         graph%nvertices = 




      else
         graph%nvertices = npix
      end if


      ! Allocate space to store vertices.


      ! Allocate temporary space to store sparse representation of adjacency matrix.

      

      ! Extract vertices and adjacency matrix.
      do ipix = 0,npix-1

!** TODO: don't keep every vertex if mask present.
         keep = .true.





         do idisc = 0,ndisc - 1

         end do
      end do



      ! Free temporary memory.



      ! Set as initialised.
      graph%init = .true.

    end function s2_graph_init


    !--------------------------------------------------------------------------
    ! s2_graph_free
    !
    !! Free all data associated with an initialised graph object and reset all
    !! other attributes.
    !!
    !! Variables:
    !!   - graph: The graph object to be freed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   September 2010 - Written by Jason McEwen
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
      if(allocated(graph%adj_matrix)) deallocate(graph%adj_matrix)
      if(graph%mask_status) call s2_sky_free(graph%mask)
      
      ! Reset attributes.
      graph%nside = 0
      graph%pix_scheme = S2_SKY_RING
      graph%support_radius = 0d0
      graph%nvertices = 0
      graph%mask_status = .false.

      graph%init = .false.

    end subroutine s2_graph_free






    ! s2_graph_io_
    ! s2_graph_vfun ??

end module s2_graph_mod

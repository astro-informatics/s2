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

  use s2_types_mod
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
    s2_graph_free, &
    s2_graph_write_file, &
    s2_graph_vals, &
    s2_graph_get_nvertices


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
    !! Variable:
    !!   - nside: Healpix resolution of corresponding sphere.
    !!   - mask: Optimal mask defining pixels to exclude from graph.
    !!   - graph: Graph object constructed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   January 2012 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_graph_init(nside, mask) result(graph)

      use pix_tools, only: nside2npix, pix2ang_nest, neighbours_nest

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
         call neighbours_nest(nside, ipix, neighbours(0:7), nneighbours)

         ! Set adjacency matrix entries for all neighbours.
         do ineighbour = 0, nneighbours-1

            ! Compute position of neighbour.
            call pix2ang_nest(nside, neighbours(ineighbour), theta_n, phi_n)
            vec_n = s2_vect_init(1.0, real(theta_n,s2_sp), real(phi_n,s2_sp))

            ! Compute anglar separation between pixels.
            dot = s2_vect_dot(vec, vec_n)         
            if (dot > 1d0) dot = 1d0  ! Remove numerical noise.
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
    !! Variables:
    !!   - graph: The graph to be freed.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   January 2012 - Written by Jason McEwen
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


    !--------------------------------------------------------------------------
    ! s2_graph_write_file
    !
    !! Write graph representation to data files.
    !!
    !! Variables:
    !!   - graph: The graph to be written to files.
    !!   - filename_prefix: Prefix of the output filenames.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   January 2012 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2_graph_write_file(graph, filename_prefix)

      type(s2_graph), intent(inout) :: graph
      character(len=*), intent(in) :: filename_prefix

      character(len=S2_STRING_LEN) :: filename
      integer :: fileid, v, vv

      ! Check object initialised.
      if(.not. graph%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_graph_write_file')
      end if 

      ! Write vertices
      fileid = 11
      write(filename, '(a,a)') trim(filename_prefix), '_vertices.dat'
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')
      do v = 0,graph%nvertices-1             
         write(fileid,'(3e20.10)') s2_vect_get_x(graph%vertices(v))
      end do
      close(fileid)

      ! Write adjacency matrix.
      write(filename, '(a,a)') trim(filename_prefix), '_adj.dat'
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')
      do v = 0,graph%nvertices-1             
         do vv = 0,graph%nvertices-1             
            write(fileid,'(e20.10)') graph%adj(vv,v)
         end do
      end do
      close(fileid)

      ! Write degree matrix (diagonal component).
      write(filename, '(a,a)') trim(filename_prefix), '_deg.dat'
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')
      do v = 0,graph%nvertices-1             
         write(fileid,'(e20.10)') graph%deg(v)
      end do
      close(fileid)

    end subroutine s2_graph_write_file


    !--------------------------------------------------------------------------
    ! s2_graph_vals
    !
    !! Extract the map values of a sky on the graph vertices.
    !!
    !! Notes:
    !!   - gvals is allocated herein and must be freed by the calling routine.
    !!
    !! Variables:
    !!   - graph: Graph representation of sky.
    !!   - sky: Sky containing map with defined values.
    !!   - gvals: The graph values, defined on each vertex of the
    !!     graph, that are extracted from the map.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   January 2012 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_graph_vals(graph, sky) result(gvals)

      type(s2_graph), intent(in) :: graph
      type(s2_sky), intent(inout) :: sky
      real(s2_sp), allocatable :: gvals(:)

      integer :: v
      integer :: fail = 0

      ! Check objects initialised.
      if((.not. graph%init) .or. (.not. s2_sky_get_init(sky))) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_graph_vals')
      end if 

      ! Check sky map defined.
      if (.not. s2_sky_get_map_status(sky)) then
         call s2_error(S2_ERROR_SKY_MAP_NOT_DEF, 's2_graph_vals')
      end if

      ! Check sky size consistent with graph.
      if (graph%nside /= s2_sky_get_nside(sky)) then
         call s2_error(S2_ERROR_SKY_NON_CONFORM, 's2_graph_vals', &
              comment_add='Sky map and graph have inconsistent nside')
      end if
      
      ! Check sky pixel scheme consistent with graph and convert if
      ! not.
      if (graph%pix_scheme /= s2_sky_get_pix_scheme(sky)) then
         call s2_sky_map_convert(sky, graph%pix_scheme)
      end if

      ! Allocate memory for gvals.
      allocate(gvals(0:graph%nvertices-1), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_graph_vals')
      end if      

      ! Extract map values.
      do v = 0,graph%nvertices-1    
         gvals(v) = s2_sky_get_map_pix(sky, graph%ivertices(v))
      end do

    end function s2_graph_vals


    !--------------------------------------------------------------------------
    ! s2_graph_get_nvertices
    !
    !! Get nvertices from the passed graph.
    !!
    !! Variables:
    !!   - graph: Graph object to get attribute of.
    !!   - nvertices: Number of vertices in the graph.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   January 2012 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2_graph_get_nvertices(graph) result(nvertices)

      type(s2_graph), intent(in) :: graph
      integer :: nvertices

      ! Check object initialised.
      if(.not. graph%init) then
        call s2_error(S2_ERROR_NOT_INIT, 's2_graph_get_nvertices')
      end if 

      nvertices = graph%nvertices

    end function s2_graph_get_nvertices


end module s2_graph_mod

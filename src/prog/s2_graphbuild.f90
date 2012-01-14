!------------------------------------------------------------------------------
! s2_graph
!
!! ...
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk) and A. Stradis
!
! Revisions:
!   December 2011 - Written by Jason McEwen and Athamos Stradis
!------------------------------------------------------------------------------

program s2_graphbuild

  use s2_types_mod
  use s2_sky_mod
  use s2_graph_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_mask, filename_map
  character(len=S2_STRING_LEN) :: filename_out_prefix, filename
  type(s2_sky) :: mask
  type(s2_sky) :: sky
  type(s2_graph) :: graph
  logical :: mask_present = .false.
  logical :: map_present = .false.
  integer :: nside = 16
  integer :: v
  real(s2_sp), allocatable :: gvals(:)
  integer :: fileid = 20
  
  ! Parse input parameters.
  filename_out_prefix = 'graph'
  call parse_options()

  ! Construct graph.
  if (mask_present) then
     mask = s2_sky_init(filename_mask, S2_SKY_FILE_TYPE_MAP)
     graph = s2_graph_init(nside, mask)
  else
     graph = s2_graph_init(nside)
  end if

  ! Save graph.
  call s2_graph_write_file(graph, filename_out_prefix)

  ! Extract values from map.
  if (map_present) then

     ! Extract values.
     sky = s2_sky_init(filename_map, S2_SKY_FILE_TYPE_MAP)
     call s2_graph_vals(graph, sky, gvals)

     ! Save values.
     write(filename, '(a,a)') trim(filename_out_prefix), '_vals.dat'
     open(unit=fileid, file=trim(filename), status='new', action='write', &
          form='formatted')
     do v = 0, s2_graph_get_nvertices(graph)-1             
        write(fileid,'(e20.10)') gvals(v)
     end do
     close(fileid)

     ! Free values.
     deallocate(gvals)

  end if

  ! Free memory.
  if (mask_present) call s2_sky_free(mask)
  if (map_present) call s2_sky_free(sky)
  call s2_graph_free(graph)


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
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2_graphbuild [-nside nside]'
            write(*,'(a)') '                     [-map filename_map (optional)]'
            write(*,'(a)') '                     [-mask filename_mask (optional)]'
            write(*,'(a)') '                     [-out filename_out_prefix]'
            stop

          case ('-nside')
            read(arg,*) nside

          case ('-map')
            filename_map = trim(arg)
            map_present = .true.

          case ('-mask')
            filename_mask = trim(arg)
            mask_present = .true.

          case ('-out')
            filename_out_prefix = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_graphbuild

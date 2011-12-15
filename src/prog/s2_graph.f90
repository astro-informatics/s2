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

program s2_graph

  use s2_types_mod
  use s2_sky_mod
  use s2_graph_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_mask, filename_out
  type(s2_sky) :: mask
  logical :: mask_present = .false.
  integer :: nside = 32
  integer :: pix_scheme = S2_SKY_NEST

  ! Parse input parameters.
  call parse_options()

  ! Load mask.
  if (mask_present) mask = s2_sky_init(filename_mask, S2_SKY_FILE_TYPE_MAP)




  write(*,*) 'Build graph representation of healpix map...'

  

  ! Free memory.
  if (mask_present) call s2_sky_free(mask)


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
            write(*,'(a)') 'Usage: s2_graph [-nside nside]'
            write(*,'(a)') '                [-pix_scheme pix_scheme]'
            write(*,'(a)') '                [-mask filename_mask (optional)]'
            write(*,'(a)') '                [-out filename_out]'
            stop

          case ('-nside')
            read(arg,*) nside

          case ('-pix_scheme')
            read(arg,*) pix_scheme
          
          case ('-mask')
            filename_mask = trim(arg)
            mask_present = .true.

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_graph
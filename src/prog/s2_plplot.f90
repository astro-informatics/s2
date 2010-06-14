!------------------------------------------------------------------------------
! s2_plplot
!
!! Read a fit pl file and plot the values to and output postscript file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of pl fits file containing pl values.
!!   - [-out filename_out]: Name of output postscript to print.
!!   - [-scale scale]: Logical to specify whether to scale y axis values to 
!!     plot dls rather than cls.
!!   - [-title title]: Optional title to print on plot.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_plplot

  use s2_types_mod
  use s2_pl_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  logical :: scale = .true.
  character(len=S2_STRING_LEN) :: title = ''
  type(s2_pl) :: pl
 
  ! Read command line options.
  call parse_options()

  ! Read pl fits file.
  pl = s2_pl_init(filename_in)

  ! Produce plot.
  call s2_pl_plot(pl, trim(filename_out), scale, trim(title))

  ! Free memory.
  call s2_pl_free(pl)


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
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2_plplot [-inp filename_in]'
            write(*,'(a)') '                  [-out filename_out]'
            write(*,'(a)') '                  [-scale scale]'
            write(*,'(a)') '                  [-title title]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-scale')
            read(arg,*) scale

          case ('-title')
            title = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_plplot

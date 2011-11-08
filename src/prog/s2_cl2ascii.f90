!------------------------------------------------------------------------------
! s2_cl2ascii
!
!! Read a fits cl file and write out an ascii file (spectrum may be scaled to 
!! dls depending on dl_status flag).
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Input fits cl file name.
!!   - [-out filename_out]: Output acsii cl/dl file name.
!!   - [-dl_status dl_status]: Logical specifying whether to scale cl values 
!!     to dl values when writting to output ascii file.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - June 2005
!
! Revisions:
!   June 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_cl2ascii

  use s2_types_mod
  use s2_pl_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  type(s2_pl) :: cl
  logical :: dl_status = .false.

  filename_out = 'out.txt'

  ! Read command line options.
  call parse_options()

  ! Read in cl file.
  cl = s2_pl_init(filename_in)

  ! Write cl values to ascii file.
  call s2_pl_io_ascii_write(filename_out, cl, dl_status)

  ! Free memory.
  call s2_pl_free(cl)


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
            write(*,'(a)') 'Usage: s2_cl2ascii [-inp filename_in]'
            write(*,'(a)') '                    [-out filename_out]'
            write(*,'(a)') '                    [-dl_status dl_status]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-dl_status')
            read(arg,*) dl_status

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2_cl2ascii

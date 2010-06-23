!------------------------------------------------------------------------------
! s2_alm2matalm
!
!! Read a HEALPix alm fits file and write the read alms to a matlab alm file
!! that can be consumed by s2ea.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_inp]: Name of HEAPix alm fits file to read.
!!   - [-out filename_out]: Name of matlab alm file to write.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   June 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_alm2matalm

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_inp, filename_out
  type(s2_sky) :: sky

  ! Parse input parameters.
  call parse_options()

  ! Initialse sky with alms read in from alm fits file.
  sky = s2_sky_init(filename_inp, S2_SKY_FILE_TYPE_ALM)

  ! Write matlab alm file.
  call s2_sky_write_matalm_file(sky, filename_out)

  ! Free memory.
  call s2_sky_free(sky)


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
            write(*,'(a)') 'Usage: s2_alm2matalm [-inp filename_inp]'
            write(*,'(a)') '                     [-out filename_out]'
            stop

          case ('-inp')
            filename_inp = trim(arg)
          
          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2_alm2matalm

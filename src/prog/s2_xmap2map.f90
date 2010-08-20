!------------------------------------------------------------------------------
! s2_xmap2map
!
!! Convert an xmap file written by s2ea_io_writemap to a Healpix map file.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_inp]: Name of map HEAPix fits file to read.
!!   - [-out filename_out]: Name of matlab map file to write.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   August 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_xmap2map

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_inp, filename_out
  type(s2_sky) :: sky
  real(s2_sp), allocatable :: map(:)
  integer :: npix, nside, ipix
  integer :: fileid = 16, fail = 0

  character(len=1), parameter :: COMMENT_CHAR = '#'
  character(len=S2_STRING_LEN) :: line, line2

  ! Parse input parameters.
  call parse_options()

  ! Open file
  open(fileid, file=filename_inp, &
       form='formatted', status='old')

  ! Ignore leading comment lines.
  line = COMMENT_CHAR
  do while(line(1:1) == COMMENT_CHAR)
     read(fileid,'(a)') line
  end do

  ! Read number of files.
  read(line, *) line2, nside
  npix = 12*nside**2

  ! Allocate space for map.  
  allocate(map(0:npix-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'spsim_analysis')
  end if

  ! Read values.
  do ipix = 0,npix-1
     read(fileid,*) map(ipix)
  end do

  ! Close file
  close(fileid)

  ! Initialise sky.
  sky = s2_sky_init(map, nside, S2_SKY_RING)

  ! Write sky to file.
  call s2_sky_write_map_file(sky, trim(filename_out))

  ! Free memory.
  deallocate(map)
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
            write(*,'(a)') 'Usage: s2_xmap2map [-inp filename_inp]'
            write(*,'(a)') '                   [-out filename_out]'
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


end program s2_xmap2map

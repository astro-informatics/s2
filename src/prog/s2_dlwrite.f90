!------------------------------------------------------------------------------
! s2_dlwrite
!
!! Write out Wigner d function values for theta=pi/2.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-lmax lmax]: Lmax to write values out for.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2007
!
! Revisions:
!   April 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2_dlwrite

  use s2_types_mod
  use s2_error_mod
  use s2_dl_mod, only: s2_dl_beta_operator

  implicit none

  real(s2_dp) :: PION2 = 1.570796326794896619231321691639751442099d0
  real(s2_dp), pointer :: dl(:, :) => null() 
  integer :: lmax = 12
  integer :: l, m, mm, fail=0

  ! Parse input parameters.
  call parse_options()

  ! Compute and write out dl values.
  do l = 0,lmax   

     ! Allocate space for dl.
     allocate(dl(-lmax:lmax, -lmax:lmax), stat=fail)
     if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2_dlwrite')
     end if
  
     ! Compute dl values.
     call s2_dl_beta_operator(dl, PION2, l)

     ! Write out
!     do m = -lmax,lmax
!        do mm = -lmax,lmax
     do m = -l,l
        do mm = -l,l
           write(*,*) 'd(', l+1, ',', m+lmax+1, ',', mm+lmax+1, ')=', dl(m,mm), ';'
        end do
     end do

     ! Deallocate dl.
     ! (Allocate and free dl each loop to ensure fresh values.)
     deallocate(dl)

  end do


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
            write(*,'(a)') 'Usage: s2_dlwrite [-lmax lmax]'
            stop

         case ('-lmax')
            read(arg,*) lmax

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program s2_dlwrite





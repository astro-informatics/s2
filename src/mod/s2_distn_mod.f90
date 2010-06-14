!------------------------------------------------------------------------------
! s2_distn_mod
!
!! Functionality to sample from uniform and Gaussian distributions.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Jason McEwen
!------------------------------------------------------------------------------

module s2_distn_mod

  use s2_types_mod, only: s2_sp, s2_dp
  use s2_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       s2_distn_sample_gauss, &
       s2_distn_sample_uniform, &
       gasdev2_dp


  !---------------------------------------
  ! Interfaces
  !---------------------------------------
  
  ! None.


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Data types
  !---------------------------------------

  ! None.


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2_distn_sample_gauss
    !
    !! Generate sample from Gaussian distribution with mean 'mean' and 
    !! standard deviation 'std' (based on gasdev2 which generates sample from
    !! Gaussian distribution with mean 0 and standard deviation 1).
    !!
    !! Variables:
    !!   - seed: Integer seed for generator.
    !!   - [mean]: Mean of Gaussian distribution to sample.
    !!   - [std]: Standard deviation of Gaussian distribution to sample.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_distn_sample_gauss(seed, mean, std) result(sample)

      integer, intent(in) :: seed
      real(s2_sp), intent(in), optional :: mean, std
      real(s2_sp) :: sample

      sample = gasdev2(seed)

      if(present(std)) sample = sample * std
      if(present(mean)) sample = sample + mean

    end function s2_distn_sample_gauss


    !--------------------------------------------------------------------------
    ! s2_distn_sample_uniform
    !
    !! Generate sample from uniform distribution in range [lower, upper) 
    !! (based on ran2 which generates sample from uniform distribution in 
    !! range [0,1)).
    !!
    !! Variables:
    !!   - seed: Integer seed for generator.
    !!   - [lower]: Lower bound of uniform distribution to sample.
    !!   - [upper]: Upper bound of uniform distribution to sample.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function s2_distn_sample_uniform(seed, lower, upper) result(sample)

      integer, intent(in) :: seed
      real(s2_sp), intent(in), optional :: lower, upper
      real(s2_sp) :: sample

      real(s2_sp) :: TOL = 1e-5

      sample = ran2(seed)

      ! Convert to uniform deviate sample in range [lower,upper) if variables
      ! present.
      if(present(lower) .and. present(upper)) then
         sample = sample * (upper-lower) + lower
         ! If lower and upper the same then return value for sample.
         if(abs(lower - upper) < TOL) sample = lower
      else if(present(lower) .or. present(upper)) then
         ! Display warning.  
         call s2_error(S2_ERROR_DISNT_BND_INVALID, &
           's2_distn_sample_uniform', comment_add = &
           'Only one of range boundaries specified, defaulting to [0,1)')
      end if

    end function s2_distn_sample_uniform


    !--------------------------------------------------------------------------
    ! gasdev2
    !
    !! Generate sample from Gaussian distribution of mean 0 and standard 
    !! deviation 1 given seed.
    !!
    !! Notes:
    !!   - gaussian deviate inpsired from gasdev (Num Rec 1992, chap 7.3),
    !!     the only difference is the use of RAN2 instead of RAN1
    !!
    !! Variables:
    !!   - idum: Seed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    ! 
    ! Revisions:
    !   March 2004 - Copied from Mike Hobson
    !--------------------------------------------------------------------------

    function gasdev2(idum)

      INTEGER idum
      REAL   gasdev2
      INTEGER iset
      REAL   fac,gset,rsq,v1,v2  !,ran2 -- removed since now internal function
      ! external ran2 -- removed since now internal function
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1        v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev2=v2*fac
         iset=1
      else
         gasdev2=gset
         iset=0
      endif
      return
      
    end function gasdev2


    !--------------------------------------------------------------------------
    ! ran2
    !
    !! Generate uniform deviate in range [0,1) given seed.
    !!
    !! Notes: 
    !!   - uniform deviate (Num rec 1992, chap 7.1), original routine
    !!     said to be 'perfect' ...
    !!
    !! Variables:
    !!   - idum: Seed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    ! 
    ! Revisions:
    !   March 2004 - Copied from Mike Hobson
    !--------------------------------------------------------------------------

    function ran2(idum)

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL   ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11       continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
         
    end function ran2


    !--------------------------------------------------------------------------
    ! gasdev2_dp
    !
    !! Generate sample from Gaussian distribution of mean 0 and standard 
    !! deviation 1 given seed.  (Using double precision.)
    !!
    !! Notes:
    !!   - gaussian deviate inpsired from gasdev (Num Rec 1992, chap 7.3),
    !!     the only difference is the use of RAN2 instead of RAN1
    !!
    !! Variables:
    !!   - idum: Seed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    ! 
    ! Revisions:
    !   March 2004 - Copied from Mike Hobson
    !--------------------------------------------------------------------------

    function gasdev2_dp(idum)

      INTEGER :: idum
      REAL(s2_dp) :: gasdev2_dp
      INTEGER :: iset
      REAL(s2_dp) :: fac,gset,rsq,v1,v2  !,ran2 -- removed since now internal function
      ! external ran2 -- removed since now internal function
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
2        v1=2.*ran2_dp(idum)-1.
         v2=2.*ran2_dp(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 2
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev2_dp=v2*fac
         iset=1
      else
         gasdev2_dp=gset
         iset=0
      endif
      return
      
    end function gasdev2_dp


    !--------------------------------------------------------------------------
    ! ran2_dp
    !
    !! Generate uniform deviate in range [0,1) given seed.
    !! (Using double precision.)
    !!
    !! Notes: 
    !!   - uniform deviate (Num rec 1992, chap 7.1), original routine
    !!     said to be 'perfect' ...
    !!
    !! Variables:
    !!   - idum: Seed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    ! 
    ! Revisions:
    !   March 2004 - Copied from Mike Hobson
    !--------------------------------------------------------------------------

    function ran2_dp(idum)

      INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL(s2_dp) :: ran2_dp,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 12 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
12       continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2_dp=min(AM*iy,RNMX)
      return
         
    end function ran2_dp


end module s2_distn_mod

!====================================================================
!                          RODRIGO BOMBARDI
!
!	Subroutines that calculate the WIlcoxon Signed Rank test
!
!====================================================================
MODULE Wilcoxon_Signed_Rank

implicit none

contains

!========================================================================
! Giving two arrays of forecast error the subroutine calculates the exact
! test of equality skill using the Wilcoxon Signed-Rank Test at 5%
! significance level.
!  array1  = An n-element single-precision floating-point vector.
!  array2  = An n-element single-precision floating-point vector.
!  rsig    = Wilcoxon Sgned_Rank significance test flag.
!            rsig = 1.0 -> reject null hypotheis
!            rsig = 0.0 -> Do not reject null hypotheis
!  missval = Value for missing value in case the test cannot be
!            calculated.
!========================================================================

  SUBROUTINE wilcoxon(array1,array2,rsig,zscore,missval)

    IMPLICIT NONE
    integer :: tt, increment, mtot, nmiss, newdim, nd
    real    :: sigmaw, sig, sig1, sig2, mt
    real, intent(out) :: rsig, zscore
    real, intent(in)  :: missval, array1(:),array2(:)
    real, dimension(:), allocatable :: tmp1,tmp2,original,&
                                    absolute,ranked,signed
  
    mtot = SIZE(array1)
!------------- checking for pairs with difference equal zero ------------
    nmiss=0
    do tt=1,mtot
       if(array1(tt)-array2(tt) .eq. 0.0) nmiss=nmiss+1
    enddo
    newdim=mtot-nmiss
    if(newdim .ge. 6)then   ! statistical significance requires at
                             ! least 5 degrees of greedom
                             ! N-1 = 5 => N=6 
       allocate(tmp1(newdim),tmp2(newdim),original(newdim))
       allocate(absolute(newdim),ranked(newdim),signed(newdim))

       if(nmiss .ne. 0)then
       !---- removing pairs with difference equal zero
          nd=1
          do tt=1,mtot
             if(array1(tt)-array2(tt) .ne. 0.0)then
                tmp1(nd)=array1(tt)
                tmp2(nd)=array2(tt)
                nd=nd+1
             endif
          enddo
       else
          tmp1(:)=array1(:)
          tmp2(:)=array2(:)
       endif
        
!-------------- calculating differences and ranks ----------------------
       original(:)=tmp1(:)-tmp2(:)
       !---- the sign function
       absolute(:)=abs(tmp1(:)-tmp2(:))
       call cranks(absolute(:),ranked(:))
       !---- calculating signed ranked scores
       where( original(:) .lt. absolute(:))
         signed=(-1.)*ranked
       elsewhere
         signed=ranked
       end where 

!-------------- Wilcoxon Signed-Rank Test ------------------------------
       sigmaw=sqrt(newdim*(newdim+1.)*(2.*newdim+1.)/24.)
       mt=newdim*(newdim+1.)/4.
       sig1=0.
       sig2=0.
       do tt=1,newdim
          if(signed(tt) .gt. 0)sig1=sig1+signed(tt)
          if(signed(tt) .lt. 0)sig2=sig2+abs(signed(tt))
       enddo
       sig=max(sig1,sig2)
       zscore=(sig-mt)/sigmaw
       if(newdim .ge. 10)then !---- z score
          if(abs(zscore) .gt. 1.96) rsig=1.
          if(abs(zscore) .le. 1.96) rsig=0.
       else
          if(newdim .eq. 6)then 
             if(sig .gt. 21.) rsig=1.
             if(sig .le. 21.) rsig=0.
          endif
          if(newdim .eq. 7)then
             if(sig .gt. 24.) rsig=1.
             if(sig .le. 24.) rsig=0.
          endif
          if(newdim .eq. 8)then
             if(sig .gt. 30.) rsig=1.
             if(sig .le. 30.) rsig=0.
          endif
          if(newdim .eq. 9)then
             if(sig .gt. 35.) rsig=1.
             if(sig .le. 35.) rsig=0.
          endif
       endif

       deallocate(tmp1,tmp2,original)
       deallocate(absolute,ranked,signed)

    else
       rsig=missval
    endif

  END SUBROUTINE wilcoxon
  
!========================================================================
! Giving an array the subroutine sorts that array and returns the
! ranks of that array. Ties between values receive the average of the ranks
!========================================================================
 
  SUBROUTINE cranks(array,ranks)
    IMPLICIT NONE
    INTEGER :: mtot, jt, beg, ned
    REAL, INTENT(in)  :: array(:)
    REAL, INTENT(out) :: ranks(:)
    REAL, DIMENSION(:), ALLOCATABLE :: tmp(:)
  
    mtot = SIZE(array)
    ALLOCATE(tmp(mtot))
    tmp(mtot)=mtot
    jt=1
    do while(jt .lt. mtot)
       if(array(jt+1) .ne. array(jt))then   !Not a tie
          tmp(jt)=jt
          jt=jt+1
       else                               ! a tie
          beg=jt
       do while(array(jt+1).eq.array(jt) .and. jt.lt.mtot)
          jt=jt+1
          ned=jt
       end do
       tmp(beg:ned)=0.5*(beg+ned)
       if(jt+1 .eq. mtot .and. array(jt+1).eq. &
          array(jt)) tmp(beg:ned+1)=0.5*(beg+ned+1)
          jt=jt+1
       endif            
    enddo
    ranks(:)=tmp(:)
!  return
  END SUBROUTINE cranks

END MODULE Wilcoxon_Signed_Rank

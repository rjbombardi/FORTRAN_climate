!====================================================================
!                          RODRIGO BOMBARDI
!
!	Subroutines that perform basic statistical calculations
! such as  sort, rank, median, correlation, rank correlation as 
! as well as statistical significance tests
!
!====================================================================
MODULE stats

implicit none

contains

!========================================================================
! Giving an array data the subroutine sorts the array data from lowest
! to highest using the Shell method
!========================================================================
  SUBROUTINE Shell_Sort(array)
 
    IMPLICIT NONE
    INTEGER :: i, j, increment
    REAL :: tmp
    REAL, INTENT(in out) :: array(:)
 
    increment = SIZE(array) / 2
    DO WHILE (increment > 0)
       DO i = increment+1, SIZE(array)
          j = i
          tmp = array(i)
          DO WHILE (j >= increment+1 .AND. array(j-increment) > tmp)
             array(j) = array(j-increment)
             j = j - increment
          END DO
          array(j) = tmp
       END DO
       IF (increment == 2) THEN
          increment = 1
       ELSE
          increment = increment * 5 / 11
       END IF      
    END DO
  END SUBROUTINE Shell_Sort

!========================================================================
! Giving an array data the subroutine sorts that array and returns the
! ranks of that array for the calculation of rank correlation. If
! Ties between vallues 
!========================================================================
 
  SUBROUTINE crank(array)
    IMPLICIT NONE
    INTEGER :: mtot, jt, beg, ned
    REAL, INTENT(in out) :: array(:)
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
    array(:)=tmp(:)
!  return
  END SUBROUTINE crank

! --------------------------------------------------------------------
! REAL FUNCTION  Median() :
!    This function receives an array X of N entries, copies its value
! to a local array tmp(), sorts tmp() and computes the median.
!    The returned value is of REAL type.
! --------------------------------------------------------------------

  REAL FUNCTION  Median(array)
      IMPLICIT  NONE
      INTEGER                         :: i, mtot
      REAL, DIMENSION(:), ALLOCATABLE :: tmp
      REAL, INTENT(in)                :: array(:)

      mtot = SIZE(array)
      ALLOCATE(tmp(mtot))
      DO i = 1, mtot                       ! make a copy
         tmp(i) = array(i)
      END DO
      CALL  Shell_Sort(tmp)                ! sort the copy
      IF (MOD(mtot,2) == 0) THEN           ! compute the median
         Median = 0.5*(tmp(mtot/2) + tmp(mtot/2+1))
      ELSE
         Median = tmp(mtot/2+1)
      END IF
  END FUNCTION  Median


! --------------------------------------------------------------------
! REAL FUNCTION  MAD() :
!    This function receives an array X of N entries and alculates 
! computes the Median ABSOLUTE DEVIATION. It is a measure of 
! variance for skeewed data.
!    The returned value is of REAL type.
! --------------------------------------------------------------------

  REAL FUNCTION  MAD(array)
      IMPLICIT  NONE
      INTEGER                         :: i, mtot
      REAL, DIMENSION(:), ALLOCATABLE :: tmp,tmp2
      REAL, INTENT(in)                :: array(:)
      REAL                            :: med

      mtot = SIZE(array)
      ALLOCATE(tmp(mtot),tmp2(mtot))
      DO i = 1, mtot                       ! make a copy
         tmp(i) = array(i)
      END DO
      med=median(tmp(:))
      tmp2(:)=abs(tmp(:)-med)
      MAD=median(tmp2(:))

  END FUNCTION  MAD



! --------------------------------------------------------------------
! REAL FUNCTION  Correlate() :
! This function receives two arrays X and Y of same N entries and
! calculates the correlation between them. The returned value is of 
! REAL type.
! --------------------------------------------------------------------
   REAL FUNCTION correlate(array1,array2)
     IMPLICIT  NONE
     INTEGER                         :: tt, mtot
     REAL                            :: nomr,var1,var2
     REAL, DIMENSION(:), ALLOCATABLE :: diff1,diff2
     REAL, INTENT(in)                :: array1(:),array2(:)

!------------------------------------------------------------------------
!                Calculating Correlation Terms
!------------------------------------------------------------------------
     mtot = SIZE(array1)
     ALLOCATE(diff1(mtot))
     ALLOCATE(diff2(mtot))
!--------------------------- x - xb -------------------------------------
     do tt=1,mtot
        diff1(tt)=array1(tt)-sum(array1(:))/mtot
     end do
!--------------------------- y - yb -------------------------------------
     do tt=1,mtot
        diff2(tt)=array2(tt)-sum(array2(:))/mtot
     end do
!------------------------ correlation -----------------------------------
     nomr=sum(diff1(:)*diff2(:))
     var1=sqrt(sum((diff1(:))**2))
     var2=sqrt(sum((diff2(:))**2))
     correlate=nomr/(var1*var2)

  END FUNCTION  correlate

!========================================================================
! Given two arrays of same size this function calculates the rank
! correlation between the two arrays
!========================================================================
 
!  SUBROUTINE r_correlate(array1,array2,rs)
  REAL FUNCTION r_correlate(array1,array2)
    IMPLICIT NONE 
    INTEGER                         :: mtot, jt
    REAL, INTENT(in)                :: array1(:),array2(:)
    REAL, DIMENSION(:), ALLOCATABLE :: sort1,sort2, crank1, crank2, &
                                     tmp1, tmp2
!    REAL                            :: diff, rs
    REAL                            :: diff

    mtot = SIZE(array1)
    ALLOCATE(sort1(mtot),sort2(mtot),crank1(mtot),crank2(mtot))
    ALLOCATE(tmp1(mtot),tmp2(mtot))
!----------------- Sorting Inout Data -------------------------
    sort1(:)=array1(:)
    sort2(:)=array2(:)
    CALL Shell_Sort(sort1)
    CALL Shell_Sort(sort2)
!----------------- Calculating Ranks -------------------------
    crank1(:)=sort1(:)
    crank2(:)=sort2(:)
    CALL crank(crank1)
    CALL crank(crank2)
!----------------- Rearranging Ranks -------------------------
    tmp1(:)=array1(:)
    tmp2(:)=array2(:)
    do jt=1,mtot
       where(array1(:) .eq. sort1(jt))
          tmp1=crank1(jt)
       end where 
       where(array2(:) .eq. sort2(jt))
          tmp2=crank2(jt)
       end where 
    end do
!----------------- Rank correlation coefficient -------------------------
    diff=sum((tmp1(:)-tmp2(:))**2) !Sum the squared difference of ranks.
!    rs=1.-6.*diff/(mtot*(mtot**2-1.))
    r_correlate=1.-6.*diff/(mtot*(mtot**2-1.))
 
!  END SUBROUTINE r_correlate
  END FUNCTION r_correlate

!========================================================================
! Given the value of the rank correlation and the number of data points
! (mtot) this subroutine verifies whether the rank correlation (rs) is
! statistically significant at the 5% level
! rsig = 1 => reject null hypothesis
! rsig = 0 => Do NOT reject null hypothesis
!========================================================================

  SUBROUTINE rank_sig(rs,mtot,rsig)

    real,   intent(in)  :: rs
    real,   intent(out) :: rsig
    integer,intent(in)  :: mtot
    real                :: tcut95, tcut98, tmp2,tmp3
    integer             :: it, tmp1, degfree
    
    degfree=mtot-2
    if(degfree .lt. 5 .or. &
       degfree .gt. 100) STOP "# deg freedom must be >= 5 or <= 100"
    if(degfree .le. 100 .and. degfree .ge. 5)then
       open(10, file='/homes/bombardi/help_files/rank_sig_1_100.txt', &
                 status='old')
       do it=1,71
          READ(10,*) tmp1,tmp2,tmp3
          if(degfree .le. 50 .and. tmp1 .eq. degfree) tcut95=tmp2
          if(degfree .gt. 50 .and. &
            (tmp1 .eq. degfree .or. tmp1+1 .eq. degfree)) tcut95=tmp2
       enddo
       close(10) 
    end if
    if(abs(rs) .ge. tcut95) rsig=1.0
    if(abs(rs) .lt. tcut95) rsig=0.0

  END SUBROUTINE rank_sig

!========================================================================
! Given the value of the correlation and the number of data points
! (mtot) this subroutine verifies whether rank correlation (corr) is
! statistically significant at the 5% level
! rsig = 1 => reject null hypothesis
! rsig = 0 => Do NOT reject null hypothesis
!========================================================================

  SUBROUTINE corr_sig(corr,mtot,rsig)

    real,   intent(in)  :: corr
    real,   intent(out) :: rsig
    integer,intent(in)  :: mtot
    real                :: tcut95, tmp2, tstu
    integer             :: it, tmp1, degfree

    degfree=mtot-2
    if(degfree .lt. 1000)then
       open(10, file='/homes/bombardi/help_files/t_student_975.txt', &
                 status='old')
       do it=1,35
          READ(10,*) tmp1,tmp2
          if(degfree .le. 30 .and. tmp1 .eq. degfree) tcut95=tmp2
          if(degfree .gt. 30 .and. degfree .lt. 40 .and. &
             tmp1 .eq. 30) tcut95=tmp2
          if(degfree .ge. 40 .and. degfree .lt. 60 .and. &
             tmp1 .eq. 40) tcut95=tmp2
          if(degfree .ge. 60 .and. degfree .lt. 80 .and. &
             tmp1 .eq. 60) tcut95=tmp2
          if(degfree .ge. 80 .and. degfree .lt. 100 .and. &
             tmp1 .eq. 80) tcut95=tmp2
          if(degfree .ge. 100 .and. degfree .lt. 1000  .and. &
             tmp1 .eq. 100) tcut95=tmp2
       enddo
       close(10) 
    end if
    if(degfree .ge. 1000) tcut95=1.96

    tstu=(corr**2)/sqrt((1-corr**2)/(degfree))
    if(abs(tstu) .ge. tcut95) rsig=1.0
    if(abs(tstu) .lt. tcut95) rsig=0.0

  END SUBROUTINE corr_sig

!========================================================================
! Given two arrays with "potentially" some missing data, this subroutine
! calculates rearranges the arrays and calculates the correlation or
! the rank correlation and its statistical significance 
!========================================================================

  SUBROUTINE arrange(array1,array2,rs,rsig,rtype,missval)
 
    INTEGER           :: mtot,nmiss,newdim,nd,tt
    REAL, intent(out) :: rs, rsig
    REAL, INTENT(in)  :: array1(:), array2(:), missval
    REAL,dimension(:), allocatable :: tmp1,tmp2
    CHARACTER*4         :: rtype

      mtot = SIZE(array1)
!-------------------- checking for missing data -------------------------
      nmiss=0
      do tt=1,mtot
         if(array1(tt) .eq. missval .or. &
            array2(tt) .eq. missval) nmiss=nmiss+1
      enddo
      newdim=mtot-nmiss
      if(newdim .ge. 7)then    ! statistical significance requires at
                               ! least 5 degress of freedom (n-2)=5
                               ! => n=7
         allocate(tmp1(newdim),tmp2(newdim))
         if(nmiss .ne. 0)then
            nd=1
            do tt=1,mtot
               if(array1(tt) .ne. missval .and. &
                  array2(tt) .ne. missval)then
                  tmp1(nd)=array1(tt)
                  tmp2(nd)=array2(tt)
                  nd=nd+1
               endif
            enddo
         else
            tmp1(:)=array1(:)
            tmp2(:)=array2(:)
         endif
         if(rtype .eq. "rank") then
            rs=r_correlate(tmp1(:),tmp2(:))
            call rank_sig(rs,newdim,rsig)
         endif
         if(rtype .eq. "corr") then
            rs=correlate(tmp1(:),tmp2(:))
            call corr_sig(rs,newdim,rsig)
         endif
      else
         rs=missval
         rsig=missval
      endif

  END SUBROUTINE arrange


! --------------------------------------------------------------------
! REAL FUNCTION  rmse() :
! This function receives two arrays X and Y of same N entries and
! calculates the Root Mean Square Error between them. The returned
! value is of REAL type.
! --------------------------------------------------------------------
   REAL FUNCTION rmse(sim,obs,missval)
     IMPLICIT  NONE
     INTEGER                         :: mtot,nmiss,newdim,nd,tt
     REAL                            :: nomr
     REAL, DIMENSION(:), ALLOCATABLE :: diff, tmp1, tmp2
     REAL, INTENT(in)                :: sim(:),obs(:),missval

!------------------------------------------------------------------------
!                Calculating Terms
!------------------------------------------------------------------------
     mtot = SIZE(sim)
!-------------------- checking for missing data -------------------------
     nmiss=0
     do tt=1,mtot
        if(sim(tt) .eq. missval .or. &
           obs(tt) .eq. missval) nmiss=nmiss+1
     enddo
     newdim=mtot-nmiss

!--------------------- verifying missing data ---------------------------
     if(newdim .ge. 3)then    ! rmse requires at least 3 records

        allocate(tmp1(newdim),tmp2(newdim),diff(newdim))
        if(nmiss .ne. 0)then
           nd=1
           do tt=1,mtot
              if(sim(tt) .ne. missval .and. &
                 obs(tt) .ne. missval)then
                 tmp1(nd)=sim(tt)
                 tmp2(nd)=obs(tt)
                 nd=nd+1
              endif
           enddo
        else
           tmp1(:)=sim(:)
           tmp2(:)=obs(:)
        end if
!------------------------ (sim - obs)^2 ---------------------------------
        do tt=1,newdim
           diff(tt)=(tmp1(tt)-tmp2(tt))**2
        end do
!------------------------------ RMSE ------------------------------------
        nomr=sum(diff(:))/newdim
        rmse=sqrt(nomr)
     else
        rmse=missval
     endif
     
  END FUNCTION  rmse


! --------------------------------------------------------------------
! REAL FUNCTION  bias() :
! This function receives two arrays X and Y of same N entries and
! calculates the Bias between them. The returned value is of REAL type.
! --------------------------------------------------------------------
   REAL FUNCTION bias(sim,obs,missval)
     IMPLICIT  NONE
     INTEGER                         :: mtot,nmiss,newdim,nd,tt
     REAL, DIMENSION(:), ALLOCATABLE :: diff, tmp1, tmp2
     REAL, INTENT(in)                :: sim(:),obs(:),missval

!------------------------------------------------------------------------
!                Calculating Terms
!------------------------------------------------------------------------
     mtot = SIZE(sim)
!-------------------- checking for missing data -------------------------
     nmiss=0
     do tt=1,mtot
        if(sim(tt) .eq. missval .or. &
           obs(tt) .eq. missval) nmiss=nmiss+1
     enddo
     newdim=mtot-nmiss

!--------------------- verifying missing data ---------------------------
     if(newdim .ge. 3)then    ! rmse requires at least 3 records

        allocate(tmp1(newdim),tmp2(newdim),diff(newdim))
        if(nmiss .ne. 0)then
           nd=1
           do tt=1,mtot
              if(sim(tt) .ne. missval .and. &
                 obs(tt) .ne. missval)then
                 tmp1(nd)=sim(tt)
                 tmp2(nd)=obs(tt)
                 nd=nd+1
              endif
           enddo
        else
           tmp1(:)=sim(:)
           tmp2(:)=obs(:)
        end if
!------------------------ (sim - obs)^2 ---------------------------------
        do tt=1,newdim
           diff(tt)=tmp1(tt)-tmp2(tt)
        end do
!------------------------------ RMSE ------------------------------------
        bias=sum(diff(:))/newdim
     else
        bias=missval
     endif
     
  END FUNCTION  bias


! --------------------------------------------------------------------
! REAL FUNCTION  mean() :
! This function receives an arrays X of N entries and
! calculates the Mean.
! The returned value is of REAL type.
! --------------------------------------------------------------------
   REAL FUNCTION mean(sim,missval)
     IMPLICIT  NONE
     INTEGER                         :: mtot,nmiss,newdim,nd,tt
     REAL, DIMENSION(:), ALLOCATABLE :: tmp1
     REAL, INTENT(in)                :: sim(:),missval
!------------------------------------------------------------------------
!                Calculating Terms
!------------------------------------------------------------------------
     mtot = SIZE(sim)
!-------------------- checking for missing data -------------------------
     nmiss=0
     do tt=1,mtot
        if(sim(tt) .eq. missval) nmiss=nmiss+1
     enddo
     newdim=mtot-nmiss

!--------------------- verifying missing data ---------------------------
     if(newdim .ge. 3)then    ! rmse requires at least 3 records

        allocate(tmp1(newdim))
        if(nmiss .ne. 0)then
           nd=1
           do tt=1,mtot
              if(sim(tt) .ne. missval)then
                 tmp1(nd)=sim(tt)
                 nd=nd+1
              endif
           enddo
        else
           tmp1(:)=sim(:)
        end if
!----------------------------- average -----------------------------------
        mean = sum(tmp1(:))/newdim
     else
        mean=missval
     endif

     
  END FUNCTION mean

! --------------------------------------------------------------------
! REAL FUNCTION  stdev() :
! This function receives an arrays X of N entries and
! calculates the Standard Mean Deviation.
! The returned value is of REAL type.
! --------------------------------------------------------------------
   REAL FUNCTION stdev(sim,missval)
     IMPLICIT  NONE
     INTEGER                         :: mtot,nmiss,newdim,nd,tt
     REAL, DIMENSION(:), ALLOCATABLE :: tmp1
     REAL, INTENT(in)                :: sim(:),missval
     REAL                            :: avg
!------------------------------------------------------------------------
!                Calculating Terms
!------------------------------------------------------------------------
     mtot = SIZE(sim)
!-------------------- checking for missing data -------------------------
     nmiss=0
     do tt=1,mtot
        if(sim(tt) .eq. missval) nmiss=nmiss+1
     enddo
     newdim=mtot-nmiss

!--------------------- verifying missing data ---------------------------
     if(newdim .ge. 3)then    ! rmse requires at least 3 records

        allocate(tmp1(newdim))
        if(nmiss .ne. 0)then
           nd=1
           do tt=1,mtot
              if(sim(tt) .ne. missval)then
                 tmp1(nd)=sim(tt)
                 nd=nd+1
              endif
           enddo
        else
           tmp1(:)=sim(:)
        end if
!----------------------------- average -----------------------------------
        avg = sum(tmp1(:))/newdim
!------------------------------ STDEV ------------------------------------
        stdev=sqrt(sum((tmp1(:)-avg)**2)/newdim)
     else
        stdev=missval
     endif
     
  END FUNCTION stdev


END MODULE stats

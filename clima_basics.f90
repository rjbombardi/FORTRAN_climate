!====================================================================
!                          RODRIGO BOMBARDI
!
!	Subroutines that perform basic calculations in climate
! sciences. For example: generating arrays of time (hours, days,
! months, and years); calculating trends; calculating the mean
! annual cycle and the smoothed mean annual cycle; and calculating
! anomalies by removing the annual cycle and/or trends
!
!====================================================================
module clima_basics

implicit none

contains
!====================================================================
!	Subroutine that generates arrays of months and years for
! monthly datasets
!====================================================================
!       month0          == first month of the dataset
!       year0           == first year of the dataset
!	mtot		== total number of points in time
!	month           == vector with monthly data
!	year            == vector with yearly data
!       month_days      == vector with number of days in each month
!====================================================================
     subroutine dates_month(mtot,month0,year0,month,year,month_days)

     integer, intent(in)      :: mtot, month0, year0
     integer, dimension(mtot) :: month, year, month_days
     integer                  :: tt
     integer, dimension(12)   :: mm
!---------------------------------------------------------------------
! Generating monthly and yearly arrays
!---------------------------------------------------------------------
     month(1) = month0
     year(1)  = year0
     do tt=2,mtot
        month(tt) = month(tt-1) + 1
        year(tt)  = year(tt-1)
        if(month(tt).gt.12)then
           month(tt) = 1
           year(tt) = year(tt-1) + 1
        end if
     end do
     mm(:)=[31,28,31,30,31,30,31,31,30,31,30,31]
     do tt=1,12
        where(month.eq.tt)
              month_days=mm(tt)
        end where
     end do
     return
!===================================================================
!   End of Subroutine    
!===================================================================
     end subroutine dates_month

!====================================================================
!	Subroutine that generates arrays of pentads and years
! for pentad datasets
!====================================================================
!       month0          == first month of the dataset
!       year0           == first year of the dataset
!	mtot		== total number of points in time
!	month           == vector with monthly data
!	year            == vector with yearly data
!       month_days      == vector with number of days in each month
!====================================================================
     subroutine dates_pentad(mtot,pentad0,year0,pentad,year)

     integer, intent(in)      :: mtot, pentad0, year0
     integer, intent(out)     :: pentad(mtot), year(mtot)
     integer                  :: tt
!---------------------------------------------------------------------
! Generating pentads and year arrays
!---------------------------------------------------------------------
     pentad(1) = pentad0
     year(1)  = year0
     do tt=2,mtot
        pentad(tt) = pentad(tt-1) + 1
        year(tt)  = year(tt-1)
        if(pentad(tt).gt.73)then
           pentad(tt) = 1
           year(tt) = year(tt-1) + 1
        end if
     end do
!===================================================================
!   End of Subroutine    
!===================================================================
     end subroutine dates_pentad



!==================================================================
!	Subroutines that generates arrays of hours, days, months,
! and years for hourly/daily datasets considering leap years
!==================================================================
     subroutine dates_day_hour(dh, mtot, day0, month0, year0, &
                               hour, day, month, year)
!============================ ATENTION ============================
!       The dataset MUST start at hour 00. Otherwise, the dataset
! MUST have daily resolution (dh = 24)
!       dh              == hourly interval. Examples:
!                          dh=1  --> One record each hour
!                          dh=3  --> One record every 3 hours
!                          dh=24 --> One record a day (MAXIMUM)
!       day0            == first day of the dataset
!       month0          == first month of the dataset
!       year0           == first year of the dataset
!	mtot		== total number of points in time
!	hour            == vector with hourly data
!	day             == vector with daily data
!	month           == vector with monthly data
!	year            == vector with yearly data
!==================================================================
!---------------------- Declaring Variables -----------------------
      integer, intent(in)      :: dh, mtot, day0, month0, year0
!      integer, dimension(mtot) :: hour, day, month, year
      integer, intent(out) :: hour(mtot), day(mtot), month(mtot), &
                              year(mtot)
      integer :: tot, tt, fin, beg, dt, mt, yt
      integer, dimension(:), allocatable :: hh
      integer, dimension(12)             :: mm
!==================================================================
! Generating hourly, daily, monthly, and yearly arrays
!==================================================================
!---- allocating memory
        tot=int(24./dh)
        allocate (hh(tot))
!---------------------------------------------------------------------
! Generating hourly arrays
!---------------------------------------------------------------------
        if(tot.eq.1) hh(1)=0
        if(tot.gt.1)then
           hh(1)=0
           do tt=2,tot
              hh(tt)=hh(tt-1)+dh
           end do
        end if

        if(mod(mtot,tot).eq.0) fin=0
        if(mod(mtot,tot).ne.0) fin=tot
        do tt=1,mtot-fin,tot 
           hour(tt:tt+tot-1)=hh(:)
        end do
!---- Correcting for incomplete datasets
        beg=mod(mtot,tot)
        if(beg.gt.0) hour(mtot-beg+1:mtot)=hh(1:beg)
!---------------------------------------------------------------------
! Generating daily, monthly, and yearly arrays
!---------------------------------------------------------------------
        dt=day0
        mt=month0
        yt=year0
!---- Checking for loop years and defyning monthly days
        if((mod(yt,4).eq.0 .and. mod(yt,100).ne.0) .or. &
           (mod(yt,4).eq.0 .and. mod(yt,400).eq.0))then
            mm(:)=[31,29,31,30,31,30,31,31,30,31,30,31]
        else
            mm(:)=[31,28,31,30,31,30,31,31,30,31,30,31]
        end if 

        do tt=1,mtot-fin,tot 
           day(tt:tt+tot-1)   = dt
           month(tt:tt+tot-1) = mt
           year(tt:tt+tot-1)  = yt
           dt=dt+1
           if(dt.gt.mm(mt))then !---- Checking for month changes
              dt=1
              mt=mt+1
              if(mt.gt.12)then     !---- Checking for year changes
                 mt=1
                 yt=yt+1
!---- Checking for loop years and defyning monthly days again
!---- since the year changed
                 if((mod(yt,4).eq.0 .and. mod(yt,100).ne.0) .or. &
                    (mod(yt,4).eq.0 .and. mod(yt,400).eq.0))then
                     mm(:)=[31,29,31,30,31,30,31,31,30,31,30,31]
                 else
                     mm(:)=[31,28,31,30,31,30,31,31,30,31,30,31]
                 end if
              end if
           end if  
        end do
!===================================================================
!   End of Subroutine    
!===================================================================
    end subroutine dates_day_hour

!==================================================================
!	Subroutines that generates arrays of hours, days, months,
! and years for hourly/daily datasets for datasets without Feb 29
!==================================================================
     subroutine dates_noleap(dh, mtot, day0, month0, year0, &
                               hour, day, month, year)
!============================ ATENTION ============================
!       The dataset MUST start at hour 00. Otherwise, the dataset
! MUST have daily resolution (dh = 24)
!       dh              == hourly interval. Examples:
!                          dh=1  --> One record each hour
!                          dh=3  --> One record every 3 hours
!                          dh=24 --> One record a day (MAXIMUM)
!       day0            == first day of the dataset
!       month0          == first month of the dataset
!       year0           == first year of the dataset
!	mtot		== total number of points in time
!	hour            == vector with hourly data
!	day             == vector with daily data
!	month           == vector with monthly data
!	year            == vector with yearly data
!==================================================================
!---------------------- Declaring Variables -----------------------
      integer, intent(in)      :: dh, mtot, day0, month0, year0
!      integer, dimension(mtot) :: hour, day, month, year
      integer, intent(out) :: hour(mtot), day(mtot), month(mtot), &
                              year(mtot)
      integer :: tot, tt, fin, beg, dt, mt, yt
      integer, dimension(:), allocatable :: hh
      integer, dimension(12)             :: mm
!==================================================================
! Generating hourly, daily, monthly, and yearly arrays
!==================================================================

!---- allocating memory
        tot=int(24./dh)
        allocate (hh(tot))
!---------------------------------------------------------------------
! Generating hourly arrays
!---------------------------------------------------------------------
        if(tot.eq.1) hh(1)=0
        if(tot.gt.1)then
           hh(1)=0
           do tt=2,tot
              hh(tt)=hh(tt-1)+dh
           end do
        end if

        if(mod(mtot,tot).eq.0) fin=0
        if(mod(mtot,tot).ne.0) fin=tot
        do tt=1,mtot-fin,tot 
           hour(tt:tt+tot-1)=hh(:)
        end do
!---- Correcting for incomplete datasets
        beg=mod(mtot,tot)
        if(beg.gt.0) hour(mtot-beg+1:mtot)=hh(1:beg)
!---------------------------------------------------------------------
! Generating daily, monthly, and yearly arrays
!---------------------------------------------------------------------
        dt=day0
        mt=month0
        yt=year0
        mm(:)=[31,28,31,30,31,30,31,31,30,31,30,31]
        do tt=1,mtot-fin,tot 
           day(tt:tt+tot-1)   = dt
           month(tt:tt+tot-1) = mt
           year(tt:tt+tot-1)  = yt
           dt=dt+1
           if(dt.gt.mm(mt))then !---- Checking for month changes
              dt=1
              mt=mt+1
              if(mt.gt.12)then     !---- Checking for year changes
                 mt=1
                 yt=yt+1
              end if
           end if  
        end do
!===================================================================
!   End of Subroutine    
!===================================================================
    end subroutine dates_noleap

!==================================================================
!	Subroutines that calculates daily data anomalies. Starts
! by calculating the mean annual cycle and smoothing the mean annual
! cycle.
!==================================================================
    subroutine anomalies_3d(bin, bin_anomal, missval, &  ! real
                           lon0, dlon, lat0, dlat,&      ! real
               nlon,nlat,ntmp,tot,nyrs,npass)      ! integer

!============================ SETTINGS ============================

    ! VARIABLE: REAL
    ! bin        = original dataset
    ! bin_anomal = anomalies
    ! missval    = missing value
    ! lon0       = initial longitude
    ! dlon       = longitude interval
    ! lat0       = initial latitude
    ! dlat       = latitude interval

    ! VARIABLE: INTEGER
    ! nlon     = number of points of longitude
    ! nlat     = number of points of latitude
    ! ntime    = number of points in time
    ! nvars    = number of variables
    ! tot      = number of points in time for one year
    ! nyrs     = number of years
    ! npass    = number of times the filter is applied
    
!====================== Defining Variables ========================
        
        integer, intent(in) :: nlon,nlat,tot,nyrs,ntmp,npass
        real, intent(in)    :: lon0,dlon,lat0,dlat,missval
        integer :: beg,ned,yt,np,tt
        integer, parameter :: nvars=2
        real    :: fac
        real, dimension(nlon,nlat,tot)  :: annualmean,annualcycle,tmp
        real, dimension(nlon,nlat,ntmp) :: bin,bin_anomal

!---------------------------------------------------------------------
! Calculating the mean annual cycle
!---------------------------------------------------------------------

        annualmean(:,:,:)=0.0
        annualcycle(:,:,:)=0.0
        beg=1
!---- multiplication is more efficient than division
        fac=1./nyrs
!---- caluclating the mean for each day/pentad/month
        do yt=1,nyrs
           ned=beg+tot-1
           tmp(:,:,1:tot)=bin(:,:,beg:ned)
           where (tmp.eq.missval)
                  tmp = 0.0
           elsewhere
                  tmp = tmp
           end where
           annualmean(:,:,1:tot) = annualmean(:,:,1:tot) + &
                           tmp(:,:,1:tot)*fac ! average
           beg=ned+1
        end do
        write(*,*) "Annual cycle calculated."
!---------------------------------------------------------------------
! Smoothing the mean annual cycle
!---------------------------------------------------------------------
        tmp(:,:,:)=annualmean(:,:,:)
        np=0
        do while(np.lt.npass)
!---- the beginning of the series
           annualcycle(:,:,1) = 0.25*tmp(:,:,tot) + &
                                0.50*tmp(:,:,1)   + &
                                0.25*tmp(:,:,2)
!---- the middle of the series
           do tt=2,tot-1
              annualcycle(:,:,tt) = 0.25*tmp(:,:,tt-1) + &
                                    0.50*tmp(:,:,tt)   + &
                                    0.25*tmp(:,:,tt+1)
           end do
!---- the end of the series
           annualcycle(:,:,tot) = 0.25*tmp(:,:,tot-1) + &
                                  0.50*tmp(:,:,tot)   + &
                                  0.25*tmp(:,:,1)
           tmp(:,:,:) = annualcycle(:,:,:)
           np=np+1 
        enddo
        write(*,*) "Annual cycle smoothed."
!---------------------------------------------------------------------
! Removing smoothed annual cycle from original data
!---------------------------------------------------------------------
        bin_anomal(:,:,:)=bin(:,:,:)
        beg=1
        do tt=1,ntmp,tot
           ned=beg+tot-1
           write(*,*) tt,beg,ned, ntmp
           if(ned.le.ntmp)then !---- complete years
              bin_anomal(:,:,beg:ned) = bin_anomal(:,:,beg:ned) - &
                         annualcycle(:,:,1:tot)
           else !---- in case the last year is incomplete
              bin_anomal(:,:,beg:ntmp) = bin_anomal(:,:,beg:ntmp) - &
                         annualcycle(:,:,1:tot-(ned-ntmp))
           endif
           beg=ned+1
        enddo
        write(*,*) "Annual cycle removed."

    return
!===================================================================
!   End of subroutine    
!===================================================================
    end subroutine anomalies_3d



!==================================================================
!	Subroutines that calculates only the mean annual cycle
!==================================================================
    subroutine annual_cycle(bin, annualmean, annualcycle, & ! real
                                                missval, &  ! real
                            nlon,nlat,ntmp,tot,nyrs,npass)  ! integer

!====================== Defining Variables ========================
        
        integer, intent(in) :: nlon,nlat,ntmp,tot,nyrs,npass
        real, intent(in)    :: missval
        integer :: beg,ned,np,tt
        real    :: fac
        real, dimension(nlon,nlat,tot)  :: annualcycle,annualmean,tmp
        real, dimension(nlon,nlat,ntmp) :: bin
!---------------------------------------------------------------------
! Calculating the mean annual cycle
!---------------------------------------------------------------------
        annualmean(:,:,:)=0.0
        annualcycle(:,:,:)=0.0
        beg=1
!---- multiplication is more efficient than division
        fac=1./nyrs
!---- caluclating the mean for each day/pentad/month
        do tt=1,nyrs
           ned=beg+tot-1
           tmp(:,:,1:tot)=bin(:,:,beg:ned)
           where (tmp.eq.missval)
                  tmp = 0.0
           elsewhere
                  tmp = tmp
           end where
           annualmean(:,:,1:tot) = annualmean(:,:,1:tot) + &
                           tmp(:,:,1:tot)*fac ! average
           beg=ned+1
        end do
        write(*,*) "Annual cycle calculated."

!---------------------------------------------------------------------
! Smoothing the mean annual cycle
!---------------------------------------------------------------------
        tmp(:,:,:)=annualmean(:,:,:)
        np=0
        do while(np.lt.npass)
!---- the beginning of the series
           annualcycle(:,:,1) = 0.50*tmp(:,:,1)   + &
                                0.50*tmp(:,:,2)
!---- the middle of the series
           do tt=2,tot-1
              annualcycle(:,:,tt) = 0.25*tmp(:,:,tt-1) + &
                                    0.50*tmp(:,:,tt)   + &
                                    0.25*tmp(:,:,tt+1)
           end do
!---- the end of the series
           annualcycle(:,:,tot) = 0.50*tmp(:,:,tot-1) + &
                                  0.50*tmp(:,:,tot)
           tmp(:,:,:) = annualcycle(:,:,:)
           np=np+1 
        enddo
        write(*,*) "Annual cycle smoothed."

        return
!===================================================================
!   End of subroutine    
!===================================================================
    end subroutine annual_cycle

!=========================================================================
!                        RODRIGO J BOMNARDI
! 
! Given two paired array data (x,y) this subroutine performs a linear
! regression analysis.
! The output from the program is the slope and y-intercept of the least-
! squares best fit straight line through the data points.
!
!=========================================================================

  subroutine linreg(tmpx,tmpy,alpha,beta,r2,missval)
  implicit none

    ! r2    = squared correlation coefficient
    ! beta  = y-intercept of least-squares best fit line
    ! ntmp  = number of data points
    ! sumx  = sum of x
    ! sumx2 = sum of x**2
    ! sumxy = sum of x * y
    ! sumy  = sum of y
    ! sumy2 = sum of y**2
    ! tmpx  = input x data
    ! tmpy  = input y data

    integer             :: tt,ntmp,nmiss,newdim
    real                :: sumx, sumx2, sumxy, sumy, sumy2
    real, intent(out)   :: r2 , beta, alpha
    real, intent(in)    :: tmpx(:), tmpy(:), missval

!-------------------- checking for missing data -------------------------
    ntmp = size(tmpx)
    nmiss=0
    do tt=1,ntmp
       if(tmpx(tt) .eq. missval .or. &
          tmpy(tt) .eq. missval) nmiss=nmiss+1
    enddo
    newdim=ntmp-nmiss
!-------------------- minimum mean square errors ------------------------
    sumxy = 0.0
    sumx  = 0.0
    sumy  = 0.0
    sumx2 = 0.0
    sumy2 = 0.0
    if(newdim .ge. 3)then    ! requires at least 3 records
       do tt=1,ntmp
          if(tmpx(tt) .ne. missval .and. tmpy(tt) .ne. missval)then
             sumx  = sumx  + tmpx(tt)             ! compute sum of x
             sumx2 = sumx2 + tmpx(tt)*tmpx(tt)    ! compute sum of x**2
             sumxy = sumxy + tmpx(tt)*tmpy(tt)    ! compute sum of x * y
             sumy  = sumy  + tmpy(tt)             ! compute sum of y
             sumy2 = sumy2 + tmpy(tt)*tmpy(tt)    ! compute sum of y**2
          endif
       end do
!----------------------------- compute slope ----------------------------
       alpha = (newdim*sumxy - sumx*sumy)/(newdim * sumx2 - sumx**2)
!------------------------- compute y-intercept --------------------------
       beta = (sumy * sumx2 - sumx * sumxy) / (newdim * sumx2 - sumx**2) 
!------------------- squared correlation coefficient --------------------
       r2 = (sumxy-sumx*sumy/newdim)/ &
             sqrt((sumx2-sumx**2/newdim)*(sumy2-sumy**2/newdim))
    else
       alpha=missval
       beta=missval
       r2=missval
    endif

  end subroutine linreg

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

!=========================================================================
!                        RODRIGO J BOMNARDI
! 
! Given two paired array data (x,y) and a level of significance this
! subroutine performs a bootstrap significance test for the linear
! coefficient of the regression analysis.
! The output from the program is the value of the coefficient at the
! level of significance provided.
!
!=========================================================================

  subroutine linreg_sig(tmpx,tmpy,t_cut,missval,r_sig)
  implicit none

    ! ntmp  = number of data points
    ! tmpx  = input x data
    ! tmpy  = input y data
    ! t_cut = level of significance
    ! r_sig = value of the linear coefficient at the level of 
    !         significance provided

    integer, parameter :: mtot=1000000, ntot=1000
    integer            :: s(1),i,ntmp,tt,jt,it,nt,ties
    integer,dimension(:), allocatable :: tmp(:)
    real                  :: r2 , beta, alpha
    real, dimension(mtot) :: rand
    real, dimension(ntot) :: test
    real, intent(out)     :: r_sig
    real, intent(in)      :: tmpx(:), tmpy(:), t_cut, missval

    ntmp = size(tmpx)
    allocate(tmp(ntmp))
    i=1
    call system_clock(count=s(1)) 
    call random_seed(size=i) 
    call random_seed(put=s(1:1)) 
    call random_number(rand)

    tmp(1)=rand(1)
    ties=0
    tt=1
    jt=1
    nt=1
    do while (tt .le. mtot .and. nt .lt. ntot)
       if(ties .eq. 0) then
          jt=jt+1
          tt=tt+1
          tmp(jt)=1+int((ntmp)*rand(tt))
       else
          tt=tt+1
          tmp(jt)=1+int((ntmp)*rand(tt))
       endif   
       ties=0
       do it=1,jt-1
          if(tmp(jt) .eq. tmp(it)) ties=1
       enddo
       if(jt .eq. ntmp)then
          call linreg(tmpx(tmp(:)),tmpy(:),alpha,beta,r2,missval)
          test(nt)=alpha
          jt=1
          nt=nt+1
       endif
    end do 

    call Shell_Sort(test(:))
    r_sig=test(int(t_cut*ntot/100.))

  end subroutine linreg_sig


!=========================================================================
!                        RODRIGO J BOMNARDI
! 
! Given an array data this subroutine detrends the data using Minimum Squares.
!
!=========================================================================

  subroutine detrend(tmp,missval)
  implicit none

    ! ntmp  = number of data points
    ! tmp   = input data

    integer             :: tt,ii,ntmp,nmiss,newdim
    real                :: r2 , beta, alpha
    real, allocatable, dimension(:) :: tmpx,tmpy
    real, intent(in out)            :: tmp(:)
    real, intent(in)               :: missval

!-------------------- checking for missing data -------------------------
    ntmp = size(tmp)
    nmiss=0
    do tt=1,ntmp
       if(tmp(tt) .eq. missval) nmiss=nmiss+1
    enddo
    newdim=ntmp-nmiss
!-------------------------- organizing data -----------------------------
    allocate(tmpx(newdim),tmpy(newdim))
    do tt=1,newdim
       tmpx(tt)=tt
    enddo
!---------------------- removing missing data ---------------------------
    if(newdim .ge. 3)then    ! requires at least 3 records
       if(newdim .ne. ntmp)then
          ii=1
          do tt=1,ntmp
             if(tmp(tt) .ne. missval) then
                tmpy(ii)=tmp(tt)
                ii=ii+1
             endif
          enddo
       else
          tmpy=tmp
       endif
!---------------------- calculating regression --------------------------
       call linreg(tmpx,tmpy,alpha,beta,r2,missval)
!--------------------- re-inserting missing data ------------------------
       if(newdim .ne. ntmp)then
          ii=1
          do tt=1,ntmp
             if(tmp(tt) .ne. missval) then
                tmp(tt)=tmpy(ii)-tt*alpha ! removing trend
                ii=ii+1
             endif
          enddo
       else
          do tt=1,ntmp
             tmp(tt)= tmpy(tt)-tt*alpha   ! removing trend
          enddo
       endif
       deallocate(tmpx,tmpy)

    endif


  end subroutine detrend

!=========================================================================
! 
! Given a day and a month the function calculates the Julian day
!
!=========================================================================

  INTEGER FUNCTION julian(day,month,year)

    IMPLICIT NONE 
    INTEGER, INTENT(in)    :: day,month,year
    INTEGER, DIMENSION(12) :: mm

    if((mod(year,4).eq.0 .and. mod(year,100).ne.0) .or. &
       (mod(year,4).eq.0 .and. mod(year,400).eq.0))then
        mm(:)=[31,29,31,30,31,30,31,31,30,31,30,31]
    else
        mm(:)=[31,28,31,30,31,30,31,31,30,31,30,31]
    end if 

    if(month .eq. 1) julian=day
    if(month .gt. 1) julian=sum(mm(1:month-1))+day
 
  END FUNCTION julian


!===================================================================
!   End of Module    
!===================================================================
end module clima_basics
!===================================================================
!    END    END    END    END    END    END    END    END    END     
!===================================================================

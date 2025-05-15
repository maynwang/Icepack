module gemdrv_sbc_rpntls

!copyright (C) 2001  MSC-RPN COMM  %%%MC2%%%

   !!======================================================================
   !!         ***  from CONCEPTS MODULE  sbc_rpntls  ***
   !!=====================================================================
   !! History ICEPACK:   0.0  !  11-2018  (M. Plante, with subroutines from above) 
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC

   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision   
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   
   REAL (kind=8), PUBLIC :: tforc_1, &  ! Julian day (decimal) atm field 1
  &                 tforc_2     ! Julian day, atm field 2
   CHARACTER (LEN=16), PUBLIC :: current_atmf, & ! RPN formatted date (atm field no 2)
  &                              Mod_runstrt_S   ! RPN formatted date (simulation starting time)

   PUBLIC incdatsd, prsdate, pdfjdate2, inter_field3
          

CONTAINS

      SUBROUTINE incdatsd(newdate, olddate, dt)
      IMPLICIT NONE
      character(LEN=16) :: newdate,olddate
      real*8 :: dt

      real*8 ::jolddate,jnewdate

      integer :: newyy,newmo,newdd,newhh,newmm,newss
      integer :: oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign      
!----------------------------------------------------------------------------      
      call prsdate(oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign,olddate,1)

      call pdfjdate(jolddate,oldyy,oldmo,olddd,oldhh,oldmm,oldss)
      jnewdate=jolddate+dt
      print *, 'jnewdate: ', jnewdate
      print *, 'joldate: ', jolddate
      print *, 'dt: ', dt
      print *, 'prsdate vars: ', oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign  
      print *, 'jolddate: ', jolddate
      call pdfcdate(newyy,newmo,newdd,newhh,newmm,newss,jnewdate)
      print *, 'pdfcdate outputs: ',newyy,newmo,newdd,newhh,newmm,newss,jnewdate
      write(newdate,12) newyy,newmo,newdd,newhh,newmm,newss
      print *, 'newdate: ', newdate
 12   format(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
      return

      END SUBROUTINE incdatsd


      SUBROUTINE prsdate(yy,mo,dd,hh,mm,ss,sign,date,mode)
      IMPLICIT NONE

      integer :: yy,mo,dd,hh,mm,ss,sign,mode
      character(LEN=*) :: date
      character(LEN=16) :: tmpdate
      character(LEN=4) :: cyy
      character(LEN=2) :: cmo,cdd,chh,cmm,css
!----------------------------------------------------------------

      if (mode == 1) then

      if (date(1:1) == '-') then
         sign = -1
         tmpdate=date(2:16)
      else
         if (date(1:1) == ' ') then
            sign = 1
            tmpdate=date(2:16)
         else
            sign = 1
            tmpdate=date(1:15)
         endif
      endif

      cyy=tmpdate(1:4)
      cmo=tmpdate(5:6)
      cdd=tmpdate(7:8)
      chh=tmpdate(10:11)
      cmm=tmpdate(12:13)
      css=tmpdate(14:15)
PRINT *, "Reading value from input file...date = ", date
      read(cyy,'(I4)') yy
      read(cmo,'(I2)') mo
      read(cdd,'(I2)') dd
      read(chh,'(I2)') hh
      read(cmm,'(I2)') mm
      read(css,'(I2)') ss

      elseif (mode == 2) then

      write(date,10) yy,mo,dd,hh,mm,ss
 10   format(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)

      else

      write(date,12) yy,mo,dd,hh
 12   format(i4.4,i2.2,i2.2,i2.2,'_000')

      endif
      
      RETURN

      END SUBROUTINE prsdate

      SUBROUTINE pdfjdate2 (jdate,yyyy,mo,dd,hh,mm,ss)
      implicit none
      real*8 jdate
      integer yyyy,mo,dd,hh,mm,ss
!
!  calculate julian calendar day
!  see cacm letter to editor by fliegel and flandern 1968
!  page 657
!
      integer jd,jyy,jmo,jdd
      real*8 one,sec_in_day
      parameter (one=1.0d0, sec_in_day=one/86400.0d0)

      jd(jyy,jmo,jdd)=jdd-32075+1461*(jyy+4800+(jmo-14)/12)/4     &
     &     +  367*(jmo-2-(jmo-14)/12*12)/12 - 3     &
     &     *((jyy+4900+(jmo-14)/12)/100)/4


      jdate = jd(yyyy,mo,dd) - 2433646 ! good from 1951 onwards
      if (jdate.lt.0) then
         print*, 'Negative Julian day in pdfjdate2 --- ABORT ---'
         stop
      endif

      jdate = jdate + dble(hh*3600+mm*60+ss)*sec_in_day

      return

      END SUBROUTINE pdfjdate2

      SUBROUTINE pdfjdate(jdate,yyyy,mo,dd,hh,mm,ss)
      IMPLICIT NONE

      real*8 :: jdate
      integer ::  yyyy,mo,dd,hh,mm,ss
!
!  calculate julian calendar day
!  see cacm letter to editor by fliegel and flandern 1968
!  page 657
!
      integer jd,jyy,jmo,jdd

! Old fortran: a  statement function
      jd(jyy,jmo,jdd)=jdd-32075+1461*(jyy+4800+(jmo-14)/12)/4    &
     &     +  367*(jmo-2-(jmo-14)/12*12)/12 - 3         &
     &     *((jyy+4900+(jmo-14)/12)/100)/4
!------------------------------------------------------------------------

      jdate = jd(yyyy,mo,dd)
      jdate = jdate + (hh*3600+mm*60+ss)/86400.0
      return

      END SUBROUTINE pdfjdate

      SUBROUTINE pdfcdate(yyyy,mo,dd,hh,mm,ss,jdate)
      IMPLICIT NONE
      real*8 jdate
      integer yyyy,mo,dd,hh,mm,ss,seconds
      real*8 :: f,rj
!--------------------------------------------------------------------

      rj = int(jdate)
      f = jdate - rj
      seconds = nint(f * 86400.0)
      
      ss = mod(seconds, 60)
      mm = mod(seconds - ss,3600)/60
      
      
      hh = (seconds-60*mm-ss) / 3600
      if (hh.eq.24) then
         hh = 0
         seconds = seconds - 86400
         rj = rj+1.0
      endif
      mm = (seconds - hh * 3600 - ss) / 60
      
      ! Comment because datec is an internal function only accessible to ECCC machine
      ! call datec(int(rj),yyyy,mo,dd)
      ! Instead, hard code the date conversion:

      !      DeltaDays = rj - 2460310 ! This is the delta days since jan 1 2024

      yyyy = 2024

      if (rj - 2460310 .lt. 32) then

          mo = 01

          dd = rj - 2460310

      elseif (rj - 2460310 .lt. 61) then

          mo = 02

          dd = rj - 2460310 - 31

      elseif (rj - 2460310 .lt. 92) then

          mo = 03

          dd = rj - 2460310 - 31-29

      elseif (rj - 2460310 .lt. 122) then

          mo = 04

          dd = rj - 2460310 - 31-29-30

      endif


      return

      END SUBROUTINE pdfcdate


      SUBROUTINE inter_field3 (b,fdta1,fdta2,fnow,ni,nj)
      implicit none
      real(wp), intent(in)                      :: b
      integer,  intent(in)                      :: ni,nj
      real(wp), INTENT(in),    DIMENSION(ni,nj) :: fdta1, fdta2
      real(wp), INTENT(inout), DIMENSION(ni,nj) :: fnow

      real(wp)                                  :: fac1

      fac1 = 1. - b
      fnow(:,:) = fac1*fdta1(:,:) + b*fdta2(:,:)
      return

      END SUBROUTINE inter_field3


end module gemdrv_sbc_rpntls

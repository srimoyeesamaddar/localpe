CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE                                                                 JVE
CJVE    J. Vincent Eccles                                            JVE
CJVE    1487 Lynnwood                                                JVE
CJVE    Logan, UT 84341                                              JVE
CJVE    435 753-3819                                                 JVE
CJVE                                                                 JVE
CJVE    All rights reserved.                                         JVE
CJVE    No license given.                                            JVE
CJVE    No warranty implied.  This software is worthless anyway.     JVE
CJVE                                                                 JVE
CJVE    This software was developed on my own time for my own        JVE
CJVE    purposes.  No funds for development of this software         JVE
CJVE    where obtained from any company nor any government.          JVE
CJVE                                                                 JVE
CJVE    The software is useless.  You are permitted to distribute    JVE
CJVE    these subroutines to any one.  If you alter the software     JVE
CJVE    please note the modification in the header and send me a     JVE
CJVE    copy to update these tools.                                  JVE
CJVE                                                                 JVE
CJVE    Please report bugs (with command file examples of bugs) or   JVE
CJVE    suggested advances to Vince at address above.                JVE 
CJVE                                                                 JVE 
CJVE    AUTHORIZED DISTRIBUTION LIST:                                JVE
CJVE                                                                 JVE
CJVE#BEGIN DISTRIBUTION LIST *************************************** JVE
CJVE                                                                 JVE
CJVE    Anybody                                                      JVE
CJVE    Anywhere                                                     JVE
CJVE                                                                 JVE 
CJVE#END DISTRIBUTION LIST ***************************************** JVE
CJVE                                                                 JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE Last modified: 09/24/04  on                                     JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  These subroutines and functions are for easy date manipulation.
C  They were written by Vince Eccles for free distribution.  The
C  come with no guarantee to work properly.  I believe all dates
C  assume a 4 digit year.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INTEGER FUNCTION IGETDOC(IY,IDOY)
C
C  IGETDOC returns a day-of-century value given a year and day-of-year.
C
      IMPLICIT NONE
      INTEGER IDOY,IY,MY
C
      MY=MOD(IY,100)
      IGETDOC=MY*365+MY/4-MY/400+IDOY
      RETURN
      END
C
      INTEGER FUNCTION IGETDOY(IY,IM,ID)
C
C  IGETDOY returns a day-of-year value given a day-of-month, month, year.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,I
      INTEGER IMM(12)
      DATA IMM/31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER IN_FEB
      EXTERNAL IN_FEB
C
      IMM(2)=IN_FEB(IY)
      IDOY=ID
      IF (IM.GT.12) THEN
        WRITE(*,*)'Bad month in IGETDOY:',IM
        STOP
      ENDIF
      DO I=1,IM-1
        IDOY=IDOY+IMM(I)
      ENDDO
      IGETDOY=IDOY
C
      RETURN
      END
C
      INTEGER FUNCTION IGETUNIXDAY(IY,IM,ID)
C
C  IGETUNIXDAY returns a unix day value given a day-of-month, 
C  month, year.  January 1, 1970 is date that the unix  
C  calander was adopted.  
C  A continuous unix day count has a day 1 
C  on January 1, 1970 AD
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,I
      INTEGER IGETDOY
      EXTERNAL IGETDOY
C
      IDOY=IGETDOY(IY,IM,ID)
      I=IY-1
      IGETUNIXDAY=I*365+I/4-I/100+I/400+IDOY-719162
      RETURN
      END
C
      SUBROUTINE GETUNIXDATE(IY,IM,ID,IUNIXDAY)
C
C  GETGDATE returns a month, day-of-month, and year given 
C  a gregorian day.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IUNIXDAY,IUDAY,IDOY
      INTEGER IGETUNIXDAY
      EXTERNAL IGETUNIXDAY
C
      ID=1
      IM=1
      IY=IUNIXDAY/365+1970
      IY=IY-(IY/4/365)+(IY/100/365)-(IY/400/365)
      IUDAY=IGETUNIXDAY(IY,IM,ID)
      DOWHILE(IUDAY.LE.IUNIXDAY)
        IY=IY+1
        IUDAY=IGETUNIXDAY(IY,IM,ID)
      ENDDO
      DOWHILE(IUDAY.GT.IUNIXDAY)
        IY=IY-1
        IUDAY=IGETUNIXDAY(IY,IM,ID)
      ENDDO
    1 IUDAY=IGETUNIXDAY(IY,IM,ID)
      IDOY=IUNIXDAY-IUDAY+1
      IF (IDOY.LE.0) THEN
        IY=IY-1
        GOTO 1
      ENDIF
      CALL GETDATE(IY,IM,ID,IDOY)
      IUDAY=IGETUNIXDAY(IY,IM,ID)
C
      RETURN
      END
C
      INTEGER FUNCTION IGETMJD(IY,IM,ID)
C
C  IGETMJD returns a Modified Julian Day (MJD) value given a 
C  day-of-month, month, year.  The Julian Date (JD) is an integer
C  counter of days beginning at noon 1 January 4713 BC, which
C  is Julian Day 0.  MJD begins at midnight and MJD drops the first 
C  two (most significant) numbers.
C  
C             MJD = JD - 2400000.5
C
C        January 30,  2000 is 51543 +1 MJD.
C        January  1, -4713 is 0 +1 MJD.
C        January  1,  2004 is 53005 +1 MJD.
C        January  1,  2009 is 54831 +1 MJD.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY
      INTEGER IGETDOY,IGETMJD1
      EXTERNAL IGETDOY,IGETMJD1
C
      IDOY=IGETDOY(IY,IM,ID)
      IGETMJD=IGETMJD1(IY,IDOY)
C
      RETURN
      END
C
      INTEGER FUNCTION IGETMJD1(IY,IDOY)
C
C  IGETMJD returns a Modified Julian Day (MJD) value given a 
C  day-of-month, month, year.  The Julian Date (JD) is an integer
C  counter of days beginning at noon 1 January 4713 BC, which
C  is Julian Day 1.  MJD begins at midnight and MJD drops the first 
C  two (most significant) numbers.
C  
C             MJD = JD - 2400000.5
C
C        January 30,  2000 is 51543 +1 MJD.
C        January  1, -4713 is 0 +1 MJD.
C        January  1,  2004 is 53005 +1 MJD.
C        January  1,  2009 is 54831 +1 MJD.
C
      IMPLICIT NONE
      INTEGER IY,IDOY,I
      INTEGER I1,I2
C
      I=IY-1
      I1=I*365
      I2=I/4-I/100+I/400  
      IGETMJD1=I1+I2+IDOY-678576
C
      RETURN
      END
C
      SUBROUTINE GETMJDDATE(IY,IM,ID,MJD)
C
C  GETMJDDATE returns a month, day-of-month, and year given 
C  a modified Julian Day (MJD).
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,MJD,MJDDAY,IDOY
      INTEGER IGETMJD
      EXTERNAL IGETMJD
C
      ID=1
      IM=1
      IY=MJD/365+1970
      IY=IY-(IY/4/365)+(IY/100/365)-(IY/400/365)
      MJDDAY=IGETMJD(IY,IM,ID)
      DOWHILE(MJDDAY.LE.MJD)
        IY=IY+1
        MJDDAY=IGETMJD(IY,IM,ID)
      ENDDO
      DOWHILE(MJDDAY.GT.MJD)
        IY=IY-1
        MJDDAY=IGETMJD(IY,IM,ID)
      ENDDO
    1 MJDDAY=IGETMJD(IY,IM,ID)
      IDOY=MJD-MJDDAY+1
      IF (IDOY.LE.0) THEN
        IY=IY-1
        GOTO 1
      ENDIF
      CALL GETDATE(IY,IM,ID,IDOY)
      MJDDAY=IGETMJD(IY,IM,ID)
C
      RETURN
      END
C
      INTEGER FUNCTION IGETGREG(IY,IM,ID)
C
C  IGETGREG returns a gregorian day value given a day-of-month, 
C  month, year.  October 15, 1582 is date that the gregorian
C  calander was adopted.  The years goes as
C
C        2 B.C. 1 B.C. 1 A.D. 2 A.D.    Gregorian Years
C       -1      0      1      2         Astronomical Years
C  
C  A continuous gregorian day count can have a day 1 anywhere
C  so January 1, 1 A.D.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,I
      INTEGER IGETDOY
      EXTERNAL IGETDOY
C
      IDOY=IGETDOY(IY,IM,ID)
      I=IY-1
      IGETGREG=I*365+I/4-I/100+I/400+IDOY
      RETURN
      END
C
      SUBROUTINE GETGDATE(IY,IM,ID,IGREG)
C
C  GETGDATE returns a month, day-of-month, and year given 
C  a gregorian day.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IGREG,IG,IDOY
      INTEGER IGETGREG
      EXTERNAL IGETGREG
C
      ID=1
      IM=1
      IY=IGREG/365+1
      IY=IY-(IY/4/365)+(IY/100/365)-(IY/400/365)
      IG=IGETGREG(IY,IM,ID)
      DOWHILE(IG.LE.IGREG)
        IY=IY+1
        IG=IGETGREG(IY,IM,ID)
      ENDDO
      DOWHILE(IG.GT.IGREG)
        IY=IY-1
        IG=IGETGREG(IY,IM,ID)
      ENDDO
    1 IG=IGETGREG(IY,IM,ID)
      IDOY=IGREG-IG+1
      IF (IDOY.LE.0) THEN
        IY=IY-1
        GOTO 1
      ENDIF
      CALL GETDATE(IY,IM,ID,IDOY)
      IG=IGETGREG(IY,IM,ID)
C
      RETURN
      END
C
      SUBROUTINE UNIXSEC2DATE(IY,IM,ID,IH,IMI,IS,IUNIXSECONDS)
C         This subroutine accepts unix seconds and supplies the DATE
C         (seconds since 1970 --> DATE)
C         Mike Howsden & Vince Eccles     19JUL2001
      INTEGER*8 IUNIXSECONDS,I8
      INTEGER IUTSEC,IY,IM,ID,IH,IMI,IS,IUNIXDAY
C
      I8=(IUNIXSECONDS/86400)+1
      IF(I8.LT.0) I8=I8-1
      IUTSEC=(IUNIXSECONDS-(I8-1)*86400)
      IUNIXDAY=I8
C                   Get date from iunixseconds.
C                   IGREG=gregorian date of 1,1,1970
      CALL GETUNIXDATE(IY,IM,ID,IUNIXDAY)
      IH=IUTSEC/3600
      IMI=(IUTSEC-IH*3600)/60
      IS=(IUTSEC-IH*3600-IMI*60)
C
      RETURN
      END
C
      INTEGER*8 FUNCTION IGETUNIXSECONDS(IY,IM,ID,IH,IMI,IS)
C        This function returns the UNIXSECONDS given the DATE
C        (DATE --> seconds since 1970)
C         Mike Howsden & Vince Eccles     19JUL2001
      INTEGER IGETUNIXDAY
      EXTERNAL IGETUNIXDAY
      INTEGER*8 IUNIXDAY
C                               get the gregorian day
      IUNIXDAY=IGETUNIXDAY(IY,IM,ID)
C                               get seconds since 1970
      IGETUNIXSECONDS=(IUNIXDAY-1)*86400+IH*3600+IMI*60+IS
C
      RETURN
      END
C
      SUBROUTINE GETDATE(IY,IM,ID,IDOY)
C
C  GETDATE returns a month and day-of-month given a year and day-of-year.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,IDO
      INTEGER IMM(12)
      DATA IMM/31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER IN_FEB
      EXTERNAL IN_FEB
C
      IMM(2)=IN_FEB(IY)
      IDO=IDOY
      DO IM=1,11
        IF (IDO.GT.IMM(IM)) THEN
          IDO=IDO-IMM(IM)
        ELSE
          GOTO 1
        END IF
      ENDDO
    1 ID=IDO
C
      IF (ID.EQ.32.OR.IDOY.GT.366) THEN
        WRITE(*,*)'BAD DATE REQUESTED IN GETDATE.'
        WRITE(*,*)'  IDOY=',IDOY
        WRITE(*,*)'    ID=',ID
        WRITE(*,*)'    IM=',IM
        WRITE(*,*)'    IY=',IY
        STOP
      ENDIF
C
      RETURN
      END
C
      INTEGER FUNCTION IN_FEB(IY)
C
C  IN_FEB returns the number of days in February for a given year.
C
      IMPLICIT NONE
      INTEGER IY
C
      IN_FEB=28
      IF( IY/4*4 .EQ. IY ) THEN
        IF( IY/100*100 .EQ. IY ) THEN
          IF( IY/400*400 .EQ. IY ) IN_FEB=29
        ELSE
          IN_FEB=29
        ENDIF
      ENDIF
C
      RETURN
      END
C
      INTEGER FUNCTION IN_YEAR(IY)
C
C  IN_YEAR returns the numbers of days in the year given a year.
C
      IMPLICIT NONE
      INTEGER IY
C
      IN_YEAR=365
      IF( IY/4*4 .EQ. IY ) THEN
        IF( IY/100*100 .EQ. IY ) THEN
          IF( IY/400*400 .EQ. IY ) IN_YEAR=366
        ELSE
          IN_YEAR=366
        ENDIF
      ENDIF
C
      RETURN
      END
C
      INTEGER FUNCTION IN_MONTH(IY,IM)
C
      IMPLICIT NONE
      INTEGER IY,IM
      INTEGER IMM(12)
      DATA IMM/31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER IN_FEB
      EXTERNAL IN_FEB
C
C  IN_MONTH returns the number of days in the month given the month and
C  year.
C
      IF (IM.EQ.2) THEN
        IN_MONTH=IN_FEB(IY)
      ELSE
        IN_MONTH=IMM(IM)
      ENDIF
C
      RETURN
      END
C
      INTEGER FUNCTION IGETMON(CMON)
C  
C  IGETMON returns the integer month given a the first 3 letters of the
C  months name (in caps).
C
      IMPLICIT NONE
      INTEGER I
      CHARACTER*3 CMON,CCMON(12)
      DATA CCMON /'JAN','FEB','MAR','APR','MAY','JUN',
     &'JUL','AUG','SEP','OCT','NOV','DEC'/
C
      DO I=1,12
        IF(CMON.EQ.CCMON(I)) THEN
          IGETMON=I
          RETURN
        ENDIF
      ENDDO
      WRITE(0,*)'opps in IGETMON:',CMON
C
      STOP
      END
C
      SUBROUTINE CGETMON(IMON,CMON)
C
C  IGETMON returns the integer month given a the first 3 letters of the
C  months name (in caps).
C
      IMPLICIT NONE
      INTEGER IMON
      CHARACTER*3 CMON,CCMON(12)
      DATA CCMON /'JAN','FEB','MAR','APR','MAY','JUN',
     &'JUL','AUG','SEP','OCT','NOV','DEC'/
C
      IF (IMON.LT.1.OR.IMON.GT.12) THEN
        WRITE(0,*)'bad month in CGETMON',IMON
        STOP
      ENDIF
      CMON=CCMON(IMON)
C
      RETURN
      END
C
      INTEGER FUNCTION IGETDOW(IY,IM,ID)
C
C  IGETDOW returns a day-of-week number given day-of-month, month, year.
C
C   igetdow s,m,t,w,t,f,s
C   igetdow 1,2,3,4,5,6,7
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,NODS0
      INTEGER IGETDOY
      EXTERNAL IGETDOY
C
      IDOY=IGETDOY(IY,IM,ID)
      IF (IM.GT.2) THEN
        NODS0=(IY-1)*365+(IY-1)/4+IDOY-1
      ELSE
        NODS0=(IY-1)*365+(IY-1)/4+IDOY-1
      END IF
      IGETDOW=MOD(NODS0,7)+1
C
      RETURN
      END
C
      SUBROUTINE UPDATE(IY,IM,ID)
C
C  UPDATE returns a day,month,year value one day advanced from the input
C  value.  Takes into account end-of-months, end-of-years, leap years.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID,IDOY,NDOY
      INTEGER IGETDOY,IN_YEAR
      EXTERNAL IGETDOY,IN_YEAR
C
      IDOY=IGETDOY(IY,IM,ID)+1
      NDOY=IN_YEAR(IY)
      IF (IDOY.GT.NDOY) THEN
        IDOY=1
        IY=IY+1
      ENDIF
      CALL GETDATE(IY,IM,ID,IDOY)
C
      RETURN
      END
C
      REAL*8 FUNCTION GET_LUNAR_AGE (IY,IM,ID,IUTSEC,LAT,LON)
C
C  PHASEMOON returns the moon for a particular date.
C
C        0 is new moon.
C        6 is 1st quarter.
C       12 is full moon.
C       18 is last quarter
C       24 is new moon.
C
C  The position is also given but is not accurate yet.
C
      IMPLICIT NONE
      INTEGER IY,IM,ID
      INTEGER LAT,LON,IUTSEC
C
      INTEGER IG,LONS
      INTEGER IGETGREG
C
      REAL*8 G0,GG,GP,TWOPI,UT
      REAL*8 ANGLE
C
      TWOPI=6.283185307D0
      G0=   725755.226388889D0    ! 19 Jan 1988 5:26 is new moon
      GP=   29.5324596774195D0    ! is period of new moon
C     G0=G0-GP*245923D0           ! zero point is long ago just to get positive nos.
C
      IG=IGETGREG(IY,IM,ID)
      UT=DBLE(IUTSEC)/3600.D0
      GG=DBLE(IG)+UT/24.0D0
C
      ANGLE=360.D0*(GG-G0)/GP
      LAT=0
      LONS=NINT(15.D0*(12.D0-UT))+360
      LON=MOD(LONS-NINT(ANGLE)+360*9000,360)
C
      GET_LUNAR_AGE=MOD((GG-G0)/GP*24.D0,24.D0)
C
      RETURN
      END
C
C ----------------------------------------------------
C
      SUBROUTINE CHANGETIME(IY,ID,IU)
C
C  given a negative day or time or too large a day or
C  time, this subroutine walks in the direction that
C  the date and time need to walk to make them proper
C  dates and times.
C
      INTEGER IY,ID,IU
      INTEGER ND,IN_YEAR
      EXTERNAL IN_YEAR
C
    1 ND=IN_YEAR(IY)
    2 IF (ID.LE.0) THEN
        IY=IY-1
        ND=IN_YEAR(IY)
        ID=ND-ID
        GO TO 2
      ELSE IF (ID.GT.ND) THEN
        IY=IY+1
        ID=ID-ND
        GOTO 1
      ENDIF
C
      IF (IU.LT.0) THEN
        IU=IU+288
        ID=ID-1
        GOTO 2
      ELSE IF (IU.GE.288) THEN
        IU=IU-288
        ID=ID+1
        GOTO 2
      ENDIF
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GETBESTTIME(IY,ID,IS)
    1 IF (IS.LT.86400) THEN
        IF (IS.GE.0) RETURN
        IS=IS+86400
        CALL DNDAY(IY,ID)
        GOTO 1
      ELSE
        IS=IS-86400
        CALL UPDAY(IY,ID)
        GOTO 1
      ENDIF
      RETURN
      END
      SUBROUTINE DNDAY(IY,ID)
      ID=ID-1
      IF (ID.LE.0) THEN
        IY=IY-1
        ND=365
        IF( IY/4*4 .EQ. IY ) THEN
          IF( IY/100*100 .EQ. IY ) THEN
            IF( IY/400*400 .EQ. IY ) ND=366
          ELSE
            ND=366
          ENDIF
        ENDIF
        ID=ND
      ENDIF
      RETURN
      END
      SUBROUTINE UPDAY(IY,ID)
      ID=ID+1
      ND=365
      IF( IY/4*4 .EQ. IY ) THEN
        IF( IY/100*100 .EQ. IY ) THEN
          IF( IY/400*400 .EQ. IY ) ND=366
        ELSE
          ND=366
        ENDIF
      ENDIF
      IF (ID.GT.ND) THEN
        IY=IY+1
        ID=1
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE                                                                 JVE
CJVE    J. Vincent Eccles                                            JVE
CJVE    1487 Lynnwood                                                JVE
CJVE    Logan, UT 84341                                              JVE
CJVE    435 753-3819                                                 JVE
CJVE                                                                 JVE
CJVE    All rights reserved.                                         JVE
CJVE    No license given.                                            JVE
CJVE    No warranty implied.  This software is worthless anyway.     JVE
CJVE                                                                 JVE
CJVE    This software was developed on my own time for my own        JVE
CJVE    purposes.  No funds for development of this software         JVE
CJVE    where obtained from any company nor any government.          JVE
CJVE                                                                 JVE
CJVE    The software is useless.  You are permitted to distribute    JVE
CJVE    these subroutines to any one.  If you alter the software     JVE
CJVE    please note the modification in the header and send me a     JVE
CJVE    copy to update these tools.                                  JVE
CJVE                                                                 JVE
CJVE    Please report bugs (with command file examples of bugs) or   JVE
CJVE    suggested advances to Vince at address above.                JVE 
CJVE                                                                 JVE 
CJVE    AUTHORIZED DISTRIBUTION LIST:                                JVE
CJVE                                                                 JVE
CJVE#BEGIN DISTRIBUTION LIST *************************************** JVE
CJVE                                                                 JVE
CJVE    Anybody                                                      JVE
CJVE    Anywhere                                                     JVE
CJVE                                                                 JVE 
CJVE#END DISTRIBUTION LIST ***************************************** JVE
CJVE                                                                 JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE
CJVE Last modified: 09/24/04  on                                     JVE
CJVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE_JVE

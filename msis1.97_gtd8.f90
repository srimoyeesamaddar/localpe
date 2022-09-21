! ===============================================================================
! MSIS® (NRL-SOF-014-1) SOFTWARE

! MSIS® is a registered trademark of the Government of the United States of 
! America, as represented by the Secretary of the Navy. Unauthorized use of 
! the trademark is prohibited. 

! The MSIS® Software (hereinafter Software) is property of the United States 
! Government, as represented by the Secretary of the Navy. The Government of 
! the United States of America, as represented by the Secretary of the Navy, 
! herein grants a non-exclusive, non-transferable license to the Software for 
! academic, non-commercial, purposes only. A user of the Software shall not: 
! (i) use the Software for any non-academic, commercial purposes, (ii) make 
! any modification or improvement to the Software, (iii) disseminate the 
! Software or any supporting data to any other person or entity who will use 
! the Software for any non-academic, commercial purposes, or (iv) copy the 
! Software or any documentation related thereto except for (a) distribution 
! among the user’s personal computer systems, archival, or emergency repair 
! purposes, or (b) distribution for non-commercial, academic purposes, without 
! first obtaining the written consent of IP Counsel for the Naval Research 
! Laboratory. 

! As the owner of MSIS®, the United States, the United States Department of 
! Defense, and their employees: (1) Disclaim any warranties, express, or 
! implied, including but not limited to any implied warranties of 
! merchantability, fitness for a particular purpose, title or non-infringement, 
! (2) Do not assume any legal liability or responsibility for the accuracy, 
! completeness, or usefulness of the software, (3) Do not represent that use of 
! the software would not infringe privately owned rights, (4) Do not warrant 
! that the software will function uninterrupted, that is error-free or that any 
! errors will be corrected.

! BY USING THIS SOFTWARE YOU ARE AGREEING TO THE ABOVE TERMS AND CONDITIONS.  
! ===============================================================================

!!! 
!!! NRLMSIS 1.97:
!!! Beta version of Neutral Atmosphere Empirical Model from the surface to 
!!! lower exosphere
!!!
!!! GTD8: Legacy wrapper with input and output arguments used in NRLMSISE-00
!!!
!!! Doug Drob (douglas.drob@nrl.navy.mil)
!!! John Emmert (john.emmert@nrl.navy.mil)
!!!
!!! ***Do not redistribute without contacting authors***
!!! =========================================================================
!     INPUT VARIABLES:
!        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
!              (Year ignored in current model)
!        SEC - UT(SEC)
!        ALT - ALTITUDE(KM)
!        GLAT - GEODETIC LATITUDE(DEG)
!        GLONG - GEODETIC LONGITUDE(DEG)
!        STL - LOCAL SOLAR TIME(Ignored in 1.97; calculated from GLAT and GLONG)
!        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
!        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
!        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!           - ARRAY CONTAINING:
!             (1) DAILY AP
!             (2) 3 HR AP INDEX FOR CURRENT TIME
!             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!                    TO CURRENT TIME
!             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
!                    TO CURRENT TIME
!        MASS - MASS NUMBER (Not used in 1.97; temperature and all densities
!                            are always calculated)
!
!     NOTES ON INPUT VARIABLES: 
!        F107 and F107A values used to generate the model correspond
!        to the 10.7 cm radio flux at the actual distance of the Earth
!        from the Sun rather than the radio flux at 1 AU. 
!
!     OUTPUT VARIABLES:
!        D(1) - HE NUMBER DENSITY(CM-3)
!        D(2) - O NUMBER DENSITY(CM-3)
!        D(3) - N2 NUMBER DENSITY(CM-3)
!        D(4) - O2 NUMBER DENSITY(CM-3)
!        D(5) - AR NUMBER DENSITY(CM-3)
!        D(6) - TOTAL MASS DENSITY (including anomalous oxygen)(GM/CM3)
!        D(7) - H NUMBER DENSITY(CM-3)
!        D(8) - N NUMBER DENSITY(CM-3)
!        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
!        T(1) - EXOSPHERIC TEMPERATURE(Not available in 1.97)
!        T(2) - TEMPERATURE AT ALTITUDE(K)
!
!     SWITCHES
!        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC8(SW),
!        WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
!        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
!        FOR THE FOLLOWING VARIATIONS
!               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
!               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
!               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
!               7 - DIURNAL               8 - SEMIDIURNAL
!               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
!              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
!              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
!              15-25 - NOT USED IN 1.97
!        To get current values of SW: CALL TRETRV8(SW)

!!! =========================================================================

subroutine gtd8(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

    use pymsis, only:  initmsis, msiscalc
    implicit none

    ! MSIS Legacy subroutine arguments
    integer,intent(in)  :: iyd
    real(4),intent(in)  :: sec
    real(4),intent(in)  :: alt
    real(4),intent(in)  :: glat
    real(4),intent(in)  :: glong
    real(4),intent(in)  :: stl
    real(4),intent(in)  :: f107a
    real(4),intent(in)  :: f107
    real(4),intent(in)  :: ap(7)
    integer,intent(in)  :: mass
    real(4),intent(out) :: d(9),t(2)

    ! MSIS 1.97 subroutine arguments
    real(8)             :: xday
    real(8)             :: xutsec
    real(8)             :: xalt
    real(8)             :: xlat
    real(8)             :: xlon
    real(8)             :: xsfluxavg,xsflux
    real(8)             :: xap(1:7)
    real(8)             :: xtn
    real(8)             :: xdn(1:10)

    !Initialization flag
    logical             :: initflag = .false.

    if (.not. initflag) then
      call initmsis()
      initflag = .true.
    endif

    ! Cast the input arguments to the legacy format
    xday = dble(mod(iyd,1000))
    xutsec = dble(sec)
    xalt = dble(alt)
    xlat = dble(glat)
    xlon = dble(glong)
    xsfluxavg = dble(f107a)
    xsflux = dble(f107)
    xap = dble(ap)

    ! Call the new subroutine
    call msiscalc(xday,xutsec,xalt,xlat,xlon,xsfluxavg,xsflux,xap,xtn,xdn)

    ! Cast the output arguments to the legacy format (mks to cgs)
    t(2) = xtn
    d(1) = dble(xdn(5))*1.0e-6  ! [He]
    d(2) = dble(xdn(4))*1.0e-6  ! [O]
    d(3) = dble(xdn(2))*1.0e-6  ! [N2]
    d(4) = dble(xdn(3))*1.0e-6  ! [O2]
    d(5) = dble(xdn(7))*1.0e-6  ! [O2]
    d(6) = dble(xdn(1))*1.0e-3  ! [Rho]
    d(7) = dble(xdn(6))*1.0e-6  ! [H]
    d(8) = dble(xdn(8))*1.0e-6  ! [N]
    d(9) = dble(xdn(9))*1.0e-6  ! [Anomalous O]

    return

end subroutine gtd8

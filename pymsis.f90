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
!!! pymsis.f90: Subroutines to initialize the model, set switches, and 
!!!             evaluate the model using a new interface
!!!
!!! Doug Drob (douglas.drob@nrl.navy.mil)
!!! John Emmert (john.emmert@nrl.navy.mil)
!!!
!!! ***Do not redistribute without contacting authors***
!!! =========================================================================

module pymsis

  use constants
  implicit none
  
  real(4)    :: swleg(1:25), swc(1:25), sav(1:25)  !Legacy switch arrays

contains

!!! ======================================================
!!! Entry point for setting switches
!!! ======================================================

  subroutine setsw(switchin)

    use constants
    use msis,only: sw

    implicit none
    
    logical,intent(in) :: switchin(0:maxnbf-1)

    sw = switchin
    return
    
  end subroutine setsw
  
!!! ======================================================
!!! Legacy switches and mapping to new switches
!!! ======================================================
  
  subroutine tselec8(sv)
  
    use msis,only: sw
    implicit none
    real(4), intent(in)  :: sv(1:25)
    integer   :: i
    
    do i = 1, 25
        sav(i) = sv(i)
        swleg(i) = amod(sv(i),2.)
        if(abs(sv(i)).eq.1.or.abs(sv(i)).eq.2.) then
          swc(i)=1.
        else
          swc(i)=0.
        endif
    enddo
    
    !Map legacy switches to new switches
    !Main effects
    sw(0)                           = .true.                !Global term must be on
    sw(csfx:csfx+nsfx-1)            = (swleg(1) .eq. 1.0)   !Solar flux
    sw(1:6)                         = (swleg(2) .eq. 1.0)   !Time independent
    sw(306:307)                     = (swleg(2) .eq. 1.0)   !Time independent (extra, solar-flux modulated terms)
    sw((/7,8,11,12,15,16,19,20/))   = (swleg(3) .eq. 1.0)   !Symmetric annual
    sw((/21,22,25,26,29,30,33,34/)) = (swleg(4) .eq. 1.0)   !Symmetric semiannual
    sw((/9,10,13,14,17,18/))        = (swleg(5) .eq. 1.0)   !Asymmetric annual
    sw((/23,24,27,28,31,32/))       = (swleg(6) .eq. 1.0)   !Asymmetric semiannual
    sw(35:94)                       = (swleg(7) .eq. 1.0)   !Diurnal
    sw(300:305)                     = (swleg(7) .eq. 1.0)   !Solar zenith angle
    sw(95:144)                      = (swleg(8) .eq. 1.0)   !Semidiurnal
    sw(145:184)                     = (swleg(14) .eq. 1.0)  !Terdiurnal
    sw(cmag:cmag+1)                 = .false.               !Geomagnetic activity mode master switch
    if((swleg(9) .gt. 0) .or. (swleg(13) .eq. 1)) sw(cmag:cmag+1) = (/.true.,.true./)  !Daily mode master switch
    if(swleg(9) .lt. 0)                           sw(cmag:cmag+1) = (/.false.,.true./) !Storm-time mode master switch
    sw(cmag+2:cmag+12)              = (swleg(9) .eq. 1.0)   !Daily geomagnetic activity terms
    sw(cmag+28:cmag+40)             = (swleg(9) .eq. -1.0)  !Storm-time geomagnetic activity terms
    sw(cspw:csfx-1)                 = ((swleg(11) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Longitudinal
    sw(cut:cut+nut-1)               = ((swleg(12) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !UT/Lon
    sw(cmag+13:cmag+25)             = ((swleg(13) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Mixed UT/Lon/Geomag (Daily mode terms)
    sw(cmag+41:cmag+53)             = ((swleg(13) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Mixed UT/Lon/Geomag (Storm-time mode terms)
    !Cross terms
    sw(csfxmod:csfxmod+nsfxmod-1)   = (swc(1) .eq. 1.0)     !Solar activity modulation
    if (swc(1) .eq. 0) then                                 !Solar activity modulation
        sw(303:305) = .false.                                   !Solar zenith angle
        sw(306:307) = .false.                                   !Time independent
        sw(447) = .false.                                       !UT/Lon
        sw(454) = .false.                                       !UT/Lon
    endif
    if (swc(3) .eq. 0) then                                 !Symmetric annual
        sw(201:204) = .false.                                   !SPW1 (2,1)
        sw(209:212) = .false.                                   !SPW1 (4,1)
        sw(217:220) = .false.                                   !SPW1 (6,1)
        sw(255:258) = .false.                                   !SPW2 (2,2)
        sw(263:266) = .false.                                   !SPW2 (4,2)
        sw(271:274) = .false.                                   !SPW2 (6,2)
    endif
    if (swc(4) .eq. 0) then                                 !Symmetric semiannual
        sw(225:228) = .false.                                   !SPW1 (2,1)
        sw(233:236) = .false.                                   !SPW1 (4,1)
        sw(241:244) = .false.                                   !SPW1 (6,1)
        sw(275:278) = .false.                                   !SPW2 (2,2)
        sw(283:286) = .false.                                   !SPW2 (4,2)
        sw(291:294) = .false.                                   !SPW2 (6,2)
    endif
    if (swc(5) .eq. 0) then                                 !Asymmetric annual
        sw(47:50) = .false.                                     !Diurnal (1,1)
        sw(51:54) = .false.                                     !Diurnal (2,1) !In MSISE-00, swc(5) is applied to all annual modulated tides
        sw(55:58) = .false.                                     !Diurnal (3,1)
        sw(59:62) = .false.                                     !Diurnal (4,1)
        sw(63:66) = .false.                                     !Diurnal (5,1)
        sw(67:70) = .false.                                     !Diurnal (6,1)
        sw(105:108) = .false.                                   !Semidiurnal (2,2)
        sw(109:112) = .false.                                   !Semidiurnal (3,2)
        sw(113:116) = .false.                                   !Semidiurnal (4,2)
        sw(117:120) = .false.                                   !Semidiurnal (5,2)
        sw(121:124) = .false.                                   !Semidiurnal (6,2)
        sw(153:156) = .false.                                   !Terdiurnal (3,3)
        sw(157:160) = .false.                                   !Terdiurnal (4,3)
        sw(161:164) = .false.                                   !Terdiurnal (5,3)
        sw(165:168) = .false.                                   !Terdiurnal (6,3)
        sw(197:200) = .false.                                   !SPW1 (1,1)
        sw(205:208) = .false.                                   !SPW1 (3,1)
        sw(213:216) = .false.                                   !SPW1 (5,1)
        sw(259:262) = .false.                                   !SPW2 (3,2)
        sw(267:270) = .false.                                   !SPW2 (5,2)
        sw(394:397) = .false.                                   !Geomag (Daily mode terms)
        sw(407:410) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
        sw(422:425) = .false.                                   !Geomag (Storm-time mode terms)
        sw(435:438) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
        sw(446)     = .false.                                   !UT/Lon
    endif
    if (swc(6) .eq. 0) then                                 !Asymmetric semiannual
        sw(221:224) = .false.                                   !SPW1 (1,1)
        sw(229:232) = .false.                                   !SPW1 (3,1)
        sw(237:240) = .false.                                   !SPW1 (5,1)
        sw(279:282) = .false.                                   !SPW2 (3,2)
        sw(287:290) = .false.                                   !SPW2 (5,2)
    endif
    if (swc(7) .eq. 0) then                                 !Diurnal
        sw(398:401) = .false.                                   !Geomag (Daily mode terms)
        sw(426:429) = .false.                                   !Geomag (Storm-time mode terms)
    endif
    if (swc(11) .eq. 0) then                                !Longitude
        sw(402:410) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
        sw(430:438) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
        sw(452:453) = .false.                                   !UT/Lon
    endif
    if (swc(12) .eq. 0) then                                !UT/Lon
        sw(411:414) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
        sw(439:440) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
    endif
    
  end subroutine tselec8

  subroutine tretrv8(svv)

    implicit none
    real(4), intent(out)  :: svv(1:25)
    integer   :: i

    do i = 1, 25
        svv(i) = sav(i)
    enddo
  
  end subroutine tretrv8

  ! =========================================================================
  ! The main MSIS subroutine entry point
  ! =========================================================================

  subroutine msiscalc(day,utsec,zalt,lat,lon,sfluxavg,sflux,ap,tn,dn)

    use constants
    use tfn
    use dfn
    use msis,only: globe, eta
    implicit none

    real(8),intent(in)  :: day
    real(8),intent(in)  :: utsec
    real(8),intent(in)  :: zalt
    real(8),intent(in)  :: lat
    real(8),intent(in)  :: lon
    real(8),intent(in)  :: sfluxavg,sflux,ap(1:7)
    real(8),intent(out) :: tn,dn(1:10)
    real(8),external    :: hgt2gph
    real(8),external    :: dilog
  
    real(8)              :: z
    real(8)              :: Vz,Wz,lnPz,lndtotz,delz
    integer              :: ix
    integer              :: i,j,kmax
    integer              :: massid

    real(8),save         :: lastday = 1.0e-32
    real(8),save         :: lastutsec = 1.0e-32
    real(8),save         :: lastlat = 1.0e-32
    real(8),save         :: lastlon = 1.0e-32
    real(8),save         :: lastz = 1.0e-32
    real(8),save         :: lastsflux = 1.0e-32
    real(8),save         :: lastsfluxavg = 1.0e-32
    real(8),save         :: lastapd = 1.0e-32
    real(8),save         :: gf(0:maxnbf-1)
    integer,save         :: iz
    real(8),save         :: Sz(-5:0,2:6)
    type(tnparm),save    :: tpro
    type(dnparm),save    :: dpro(1:10)

    z = hgt2gph(lat,zalt*1000.0d0)/1000.0d0

    ! If only altitude changes then update the local spline
    ! weights
    if (z .lt. zetaB) then
        if (z .ne. lastz) then
            if (z .lt. zetaF) then
                kmax = 5
            else
                kmax = 6
            endif
            call bspline(z,nodes,nd+2,kmax,eta,Sz,iz)
            lastz = z
        endif
    endif

     ! If location, time or solar flux conditions change
     ! the (re)compute the profile parameters
     if ((day .ne. lastday) .or. (utsec .ne. lastutsec) .or. &
         (lat .ne. lastlat) .or. (lon .ne. lastlon) .or.     &
         (sflux .ne. lastsflux) .or. (sfluxavg .ne. lastsfluxavg) .or. &
         (ap(1) .ne. lastapd)) then
		call globe(day,utsec,lat,lon,sfluxavg,sflux,ap,gf)
        call tfnparm(gf,tpro)
        do massid = 2,9
            call dfnparm(massid,gf,tpro,dpro(massid))
        enddo
        lastday = day
        lastutsec = utsec
        lastlat = lat
        lastlon = lon
        lastsflux = sflux
        lastsfluxavg = sfluxavg
        lastapd = ap(1)
    endif

    ! Compute the temperature at altitude
    tn = tfnx(z,iz,Sz(-3:0,4),tpro)

     ! Total number density and temperature integration terms
    delz = z - zetaB
    if (z .lt. zetaF) then
        i = max(iz-4,0)
        if (iz .lt. 4) then
           j = -iz
        else
           j = -4
        endif
        Vz = dot_product(tpro%beta(i:iz),Sz(j:0,5)) + tpro%cVS
        Wz = 0d0
        lnPz = lnP0 - fh*(Vz - tpro%Vzeta0)
        lndtotz = lnPz - log(boltz*tn)
    else
        if (z .lt. zetaB) then
            Vz = dot_product(tpro%beta(iz-4:iz),Sz(-4:0,5)) + tpro%cVS
            Wz = dot_product(tpro%gamma(iz-5:iz),Sz(-5:0,6)) + tpro%cVS*delz + tpro%cWS
        else
            Vz = (delz + log(tn/tpro%tex)/tpro%sigma)/tpro%tex + tpro%cVB
            Wz = (0.5d0*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
                 + tpro%cVB*delz + tpro%cWB
        endif
    endif
        
    ! Compute the species number densities at altitude
    dn(2:10) = 1.0d0
    do massid = 2,9
        dn(massid) = dfnx(z,tn,lndtotz,Vz,Wz,tpro,dpro(massid))
    enddo
    dn(1) = dot_product(dn(2:9),amu(2:9))*amu2kg

    return
	
  end subroutine msiscalc

  
  ! =========================================================================
  ! Load in an model parameter set estimate
  ! =========================================================================

  subroutine initmsis

    use constants
    use msis,only: initparmspace,loadparmset,haveparmspace

    implicit none

    character(128)      	  :: betain = 'msis1.97.bin'   !Change path to parameter file here, as necessary

    ! -------------------------------------------------------------
    ! Initialize parameter arrays
    ! -------------------------------------------------------------
    if (.not. haveparmspace) call initparmspace()  ! allocate parameter arrays

    ! -------------------------------------------------------------
    ! Load a priori parameter set
    ! -------------------------------------------------------------
    call loadparmset(betain)

    return

  end subroutine initmsis

end module pymsis

!!!
!!! Authors:  Douglas P. Drob and John Emmert, NRL 7632
!!!
!!! =================================================
!!! MSIS horizontal expansion functions
!!! =================================================

module msis

  use constants

  implicit none
  
  logical                :: haveparmspace = .false.
  logical                :: zsfx(0:mbf) = .false.    !Flags zonal mean terms to be modulated by F1 (MSISE-00 legacy multiplier)
  logical                :: tsfx(0:mbf) = .false.    !Flags tide terms to be modulated by F2 (MSISE-00 legacy multiplier)
  logical                :: psfx(0:mbf) = .false.    !Flags SPW terms to be modulated by F3 (MSISE-00 legacy multiplier)
  logical                :: smod(0:nl) = .false.    !Flags which temperature levels get solar flux modulation; loadparmset turns flags on based on parameter values
  logical                :: sw(0:maxnbf-1) = .true.  !Switch array for globe subroutine.

  real(8)                :: plg(0:maxn,0:maxn)
  real(8)                :: cdoy(2), sdoy(2)
  real(8)                :: clst(3), slst(3)
  real(8)                :: clon(2), slon(2)
  real(8)                :: lastlat = -999.9
  real(8)                :: lastdoy = -999.9
  real(8)                :: lastlst = -999.9
  real(8)                :: lastlon = -999.9

  real(8)                :: sfluxavgref = 150.0 ! Reference F10.7 value (=150 in NRLMSISE-00)

  ! The basis definition and model parameters
  type basissubset
    sequence
    character(8)              :: name
    integer                   :: bl,nl
    real(8),allocatable       :: beta(:,:)
    logical,allocatable       :: active(:,:)
    integer,allocatable       :: fitb(:,:)
  end type basissubset

  type (basissubset)     :: TN
  type (basissubset)     :: PR
  type (basissubset)     :: N2
  type (basissubset)     :: O2
  type (basissubset)     :: O1
  type (basissubset)     :: HE
  type (basissubset)     :: H1
  type (basissubset)     :: AR
  type (basissubset)     :: N1
  type (basissubset)     :: OA   !Anomalous O
  type (basissubset)     :: NO
  
  ! Cubic B-spline quadrature locations and weights
  real(8)                :: qzm(0:nl),qzp(0:nl)
  real(8)                :: mwght(0:p-1,0:nl),pwght(0:p-1,0:nl)
  
  ! Reciprocal node difference arrays (constant values needed for B-spline calculations)
  real(8)                :: eta(0:30,2:6) = 0d0
  real(8)                :: etaO1(0:30,2:6) = 0d0
  real(8)                :: etaNO(0:30,2:6) = 0d0
  
contains

! ===============================================
! Initialize (allocate) the model parameter space
! ===============================================

    subroutine initparmspace()

        implicit none
        integer  :: n,m,j,k

        ! --------------------------------------
        ! Model formulation parameter subsets
        ! --------------------------------------
        call initsubset(TN,0,nl,maxnbf,'TN')
        call initsubset(PR,0,nl,maxnbf,'PR')

        call initsubset(N2,0,nls,maxnbf,'N2')
        call initsubset(O2,0,nls,maxnbf,'O2')
        call initsubset(HE,0,nls,maxnbf,'HE')
        call initsubset(H1,0,nls,maxnbf,'H1')
        call initsubset(AR,0,nls,maxnbf,'AR')
        call initsubset(N1,0,nls,maxnbf,'N1')
        call initsubset(OA,0,nls,maxnbf,'OA')
        ! Atomic O: 3 for mix to diffusive, 3 Chapman parameters, 8 independent splines
        call initsubset(O1,0,nls+nsplO1,maxnbf,'O1')
        ! Nitric Oxide: 3 for mix to diffusive, 3 Chapman parameters, 8 independent splines
        call initsubset(NO,0,nls+nsplNO,maxnbf,'NO')
        
        ! --------------------------------------
        ! Set solar flux modulation flags
        ! --------------------------------------
        zsfx(:) = .false.
        tsfx(:) = .false.
        psfx(:) = .false.
        ! F1, solar flux modulation of the zonal mean asymmetric annual terms
        zsfx(9:10) = .true.    !Pl(1,0) annual terms
        zsfx(13:14) = .true.   !Pl(3,0) annual terms
        zsfx(17:18) = .true.   !Pl(5,0) annual terms
        ! F2, solar sflux modulation of the tides
        tsfx(ctide:cspw-1) = .true.
        ! F3, solar sflux modulation of stationary planetary wave 1
        psfx(cspw:cspw+59) = .true. 

        ! --------------------------------------------
        ! Calculate reciprocal node difference arrays
        ! --------------------------------------------
        do k = 2, 6
            do j = 0, nl
                eta(j,k) = 1d0 / (nodes(j+k-1) - nodes(j))
            enddo
        enddo
        do k = 2, 4
            do j = 0, ndO1-k+1
                etaO1(j,k) = 1d0 / (nodesO1(j+k-1) - nodesO1(j))
            enddo
            do j = 0, ndNO-k+1
                etaNO(j,k) = 1d0 / (nodesNO(j+k-1) - nodesNO(j))
            enddo
        enddo

        haveparmspace = .true.

        return

        ! -----------------------------------------
        ! subset subroutine 
        ! -----------------------------------------
        contains

            subroutine initsubset(subset,bl,nl,maxnbf,name)
                implicit none
                type (basissubset),intent(inout) :: subset
                integer,intent(in)               :: bl
                integer,intent(in)               :: nl
                integer,intent(in)               :: maxnbf
                character(2)                     :: name
                integer                          :: iz
                subset%name = name
                subset%bl = bl
                subset%nl = nl
                allocate(subset%beta(0:maxnbf-1,bl:nl), &
                         subset%active(0:maxnbf-1,bl:nl), &
                         subset%fitb(0:maxnbf-1,bl:nl))
                subset%beta = 0.0d0
                subset%beta(cmag+1,:) = 1.0d0  !This nonlinear geomagnetic activity parameter must be nonzero
                subset%active = .false.
                subset%fitb = 0
                return
            end subroutine initsubset

    end subroutine initparmspace

!!! =====================================================
!!! Read in a parameter file
!!! =====================================================

  subroutine loadparmset(name)

    implicit none

    character(128),intent(in)  :: name
    integer                    :: i,k
    logical                    :: havefile

    inquire(file=trim(name),exist=havefile)
    if (havefile) then
       open(unit=1,file=trim(name),status='old',access='stream')
    else
       print *,"MSIS parameter set ",trim(name)," not found. Stopping."
       stop
    endif

    smod(:) = .false.
    
    do i = 0,nl
        read(1) (TN%beta(k,i),k=0,maxnbf-1)
        !Set solar flux modulation flags; if on for a given vertical parameter, then sfluxmod is called by tfnparm
        if ((Tn%beta(csfxmod+0,i) .ne. 0) .or. &
           (Tn%beta(csfxmod+1,i) .ne. 0) .or. &
           (Tn%beta(csfxmod+2,i) .ne. 0)) smod(i) = .true.
    enddo

    read(1) (PR%beta(k,0),k=0,maxnbf-1)

    do i = N2%bl,N2%nl
        read(1) (N2%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = O2%bl,O2%nl
        read(1) (O2%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = O1%bl,O1%nl
        read(1) (O1%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = HE%bl,HE%nl
        read(1) (HE%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = H1%bl,H1%nl
        read(1) (H1%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = AR%bl,AR%nl
        read(1) (AR%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = N1%bl,N1%nl
        read(1) (N1%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = OA%bl,OA%nl
        read(1) (OA%beta(k,i),k=0,maxnbf-1)
    enddo
    
    do i = NO%bl,NO%nl
        read(1) (NO%beta(k,i),k=0,maxnbf-1)
    enddo

    close(1)

    call pressparm()

    return

end subroutine loadparmset

!!! ==========================================================
!!! Calculate pressure parameters from temperature parameters
!!! ==========================================================

  subroutine pressparm()

    implicit none

    integer                    :: j,b,iz
    real(8)                    :: lnz

    !Integrate pressure on nodes up to the last fully mixed level
    do j = 0,mbf
        lnz = 0.0
        do b = 0,3
            lnz = lnz + TN%beta(j,b)*gwht(b)*fh
        enddo
        PR%beta(j,1) = -lnz
        do iz = 1,izfmx
            lnz = 0.0
            do b = 0,3
                lnz = lnz + TN%beta(j,iz+b)*gwht(b)*fh
            enddo
            PR%beta(j,iz+1) = PR%beta(j,iz) - lnz
        enddo
    enddo

    return

end subroutine pressparm

!!! ========================================================
!!! Calculate horizontal and time-dependent basis functions
!!! (Same purpose as NRLMSISE-00 "GLOBE7" subroutine)
!!! ========================================================

  subroutine globe(doy,utsec,lat,lon,sfluxavg,sflux,ap,bf)

    implicit none
    real(8),intent(in)   :: doy
    real(8),intent(in)   :: utsec
    real(8),intent(in)   :: lat
    real(8),intent(in)   :: lon
    real(8),intent(in)   :: sfluxavg,sflux,ap(1:7)
    real(8),intent(out)  :: bf(0:maxnbf-1)

    real(8)              :: lst
    real(8)              :: slat, clat, clat2, clat4, slat2
    real(8)              :: cosdoy,sindoy
    real(8)              :: coslon,sinlon
    real(8)              :: pl
    real(8)              :: coslst,sinlst
    real(8)              :: dfa,df
    real(8)              :: theta
    real(8)              :: sza
    integer              :: n,m,l,s,c

    ! ---------------------------------
    ! Legendre polynomials
    ! ---------------------------------
    if (lat .ne. lastlat) then

        clat = sin(lat*deg2rad)  ! clat <=> sin, Legendre polyomial defined in colat
        slat = cos(lat*deg2rad)  ! slat <=> cos, Legendre polyomial defined in colat
        clat2 = clat*clat
        clat4 = clat2*clat2
        slat2 = slat*slat

        plg(0, 0) = 1d0
        plg(1, 0) = clat
        plg(2, 0) = 0.5d0 * (3d0 * clat2 - 1d0)
        plg(3, 0) = 0.5d0 * (5d0 * clat * clat2 - 3d0 * clat)
        plg(4, 0) = (35d0 * clat4 - 30d0 * clat2 + 3d0)/8d0
        plg(5, 0) = (63d0 * clat2 * clat2 * clat - 70d0 * clat2 * clat + 15d0 * clat)/8d0
        plg(6, 0) = (11d0 * clat * plg(5, 0) - 5d0 * plg(4, 0))/6d0

        plg(1, 1) = slat
        plg(2, 1) = 3d0 * clat * slat
        plg(3, 1) = 1.5d0 * (5d0 * clat2 - 1d0) * slat
        plg(4, 1) = 2.5d0 * (7d0 * clat2 * clat - 3d0 * clat) * slat
        plg(5, 1) = 1.875d0 * (21d0 * clat4 - 14d0 * clat2 + 1d0) * slat
        plg(6, 1) = (11d0 * clat * plg(5, 1) - 6d0 * plg(4, 1))/5d0

        plg(2, 2) = 3d0 * slat2
        plg(3, 2) = 15d0 * slat2 * clat
        plg(4, 2) = 7.5d0 * (7d0 * clat2 - 1d0) * slat2
        plg(5, 2) = 3d0 * clat * plg(4, 2) - 2d0 * plg(3, 2)
        plg(6, 2) = (11d0 * clat * plg(5, 2) - 7d0 * plg(4, 2))/4d0

        plg(3, 3) = 15d0 * slat2 * slat
        plg(4, 3) = 105d0 * slat2 * slat * clat
        plg(5, 3) = (9d0 * clat * plg(4, 3) - 7d0 * plg(3, 3))/2d0
        plg(6, 3) = (11d0 * clat * plg(5, 3) - 8d0 * plg(4, 3))/3d0

        lastlat = lat

    endif

    ! ------------------------------------------------------------
    ! Fourier harmonics of local time, day of year, and longitude
    ! ------------------------------------------------------------

    lst = mod(utsec/3600.0 + lon/15.0 + 24.0,24.0)
    if (lst .ne. lastlst) then
       clst(1) = cos(lst2rad*lst)
       slst(1) = sin(lst2rad*lst)
       clst(2) = cos(lst2rad*lst*2.0)
       slst(2) = sin(lst2rad*lst*2.0)
       clst(3) = cos(lst2rad*lst*3.0)
       slst(3) = sin(lst2rad*lst*3.0)
       lastlst = lst
    endif

    if (doy .ne. lastdoy) then
       cdoy(1) = cos(doy2rad*doy)
       sdoy(1) = sin(doy2rad*doy)
       cdoy(2) = cos(doy2rad*doy*2.0)
       sdoy(2) = sin(doy2rad*doy*2.0)
       lastdoy = doy
    endif

    if (lon .ne. lastlon) then
        clon(1) = cos(deg2rad*lon)
        slon(1) = sin(deg2rad*lon)
        clon(2) = cos(deg2rad*lon*2.0)
        slon(2) = sin(deg2rad*lon*2.0)
        lastlon = lon
    endif

    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/
    !                  Coupled Linear Terms
    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    ! ======== Zonal Average (modulated) ==========
    ! time independent, annual, and semiannual
    bf(:) = 0.0
    c = ctimeind
    do n = 0,amaxn
       bf(c) = plg(n,0)
       c = c + 1
    enddo
    do s = 1,amaxs
       cosdoy = cdoy(s)
       sindoy = sdoy(s)
       do n = 0,amaxn
          pl = plg(n,0)
          bf(c) = pl*cosdoy
          bf(c+1) = pl*sindoy
          c = c + 2
       enddo
    enddo

    ! ============= Local Time Variation / Tides (modulated) ==============
    !  Diurnal, semidiurnal, terdiurnal with Annual and semiannual modulation
    if (c .ne. ctide) stop 'problem with basis definitions'
    do l = 1,tmaxl
       coslst = clst(l)
       sinlst = slst(l)
       do n = l,tmaxn
          pl = plg(n,l)
          bf(c) = pl*coslst
          bf(c+1) = pl*sinlst
          c = c + 2
       enddo
       do s = 1,tmaxs
          cosdoy = cdoy(s)
          sindoy = sdoy(s)
          do n = l,tmaxn
             pl = plg(n,l)
             bf(c) = pl*coslst*cosdoy
             bf(c+1) = pl*sinlst*cosdoy
             bf(c+2) = pl*coslst*sindoy
             bf(c+3) = pl*sinlst*sindoy
             c = c + 4
          enddo
       enddo
    enddo

    ! ============= Planetary waves ===================
    if (c .ne. cspw) stop 'problem with basis definitions'
    do m = 1,pmaxm
       coslon = clon(m)
       sinlon = slon(m)
       do n = m,pmaxn
          pl = plg(n,m)
          bf(c) = pl*coslon
          bf(c+1) = pl*sinlon
          c = c + 2
       enddo
       do s = 1,pmaxs
          cosdoy = cdoy(s)
          sindoy = sdoy(s)
          do n = m,pmaxn
             pl = plg(n,m)
             bf(c) = pl*coslon*cosdoy
             bf(c+1) = pl*sinlon*cosdoy
             bf(c+2) = pl*coslon*sindoy
             bf(c+3) = pl*sinlon*sindoy
             c = c + 4
          enddo
       enddo
    enddo
    
    ! ============= Linear solar flux terms ===================
    if (c .ne. csfx) stop 'problem with basis definitions'
    dfa = sfluxavg - sfluxavgref
    df = sflux - sfluxavg
    bf(c) = dfa            ! sfxon(0)
    bf(c+1) = dfa*dfa      ! sfxon(1)
    bf(c+2) = df           ! sfxon(2)
    bf(c+3) = df*df        ! sfxon(3)
    bf(c+4) = df*dfa       ! sfxon(4)
    c = c + nsfx

    ! ============= Additional terms ===================
    if (c .ne. cextra) stop 'problem with basis definitions'
    sza = solzen(doy,lst,lat,lon)
    bf(c)   = -0.5*tanh((sza-90.0)/4.)  !Solar zenith angle logistic function (transition width 2 deg)
    bf(c+1) = -0.5*tanh((sza-90.0)/20.)  !Solar zenith angle logistic function (transition width 10 deg)
    bf(c+2) = cos(sza*deg2rad)                       !cosine of solar zenith angle
    bf(c+3) = dfa * bf(c)       !Solar flux modulation of logistic sza term
    bf(c+4) = dfa * bf(c+1)     !Solar flux modulation of logistic sza term
    bf(c+5) = dfa * bf(c+2)     !Solar flux modulation of cos(sza) term
    bf(c+6) = dfa * plg(2,0)    !Solar flux modulation of P(2,0) term
    bf(c+7) = dfa * plg(4,0)    !Solar flux modulation of P(4,0) term
    
    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/
    !          Nonlinear terms
    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    c = cnonlin

    ! ============= Solar flux modulation terms ===================
    if (c .ne. csfxmod) stop 'problem with basis definitions'
    bf(c) = dfa  
    bf(c+1) = dfa*dfa
    bf(c+2) = df 
    bf(c+3) = df*df
    bf(c+4) = df*dfa
    c = c + nsfxmod

    ! ======== Terms needed for legacy geomagnetic activity dependence =============
    if (c .ne. cmag) stop 'problem with basis set'
    bf(c:c+6) = ap - 4.0
    bf(c+8) =   doy2rad*doy
    bf(c+9) =   lst2rad*lst
    bf(c+10) =  deg2rad*lon
    bf(c+11) =  lst2rad*utsec/3600.0
    bf(c+12) =  abs(lat)
    c = c + 13
    do m = 0,1
        do n = 0,amaxn
           bf(c) = plg(n,m)
           c = c + 1
        enddo
    enddo

    ! ======== Terms needed for legacy UT dependence =============
    c = cut
    bf(c) =   lst2rad*utsec/3600.0
    bf(c+1) = doy2rad*doy
    bf(c+2) = dfa
    bf(c+3) = deg2rad*lon
    bf(c+4) = plg(1,0)
    bf(c+5) = plg(3,0)
    bf(c+6) = plg(5,0)
    bf(c+7) = plg(3,2)
    bf(c+8) = plg(5,2)

    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/
    !         Apply Switches
    ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    where(.not. sw(0:mbf)) bf(0:mbf) = 0.0
    
    return

  end subroutine globe

! ================================================================
! Nonlinear solar flux modulation cross terms
! ================================================================

real(8) function sfluxmod(iz,gf,parmset,dffact)

    implicit none

    integer,intent(in)            :: iz
    real(8),intent(in)            :: gf(0:maxnbf-1)
    type(basissubset),intent(in)  :: parmset
    real(8),intent(in)            :: dffact  !Turns on or adjusts the delta-F terms added to F1 and F2 (eqns. A22b and A22c in Hedin (1987)).

    real(8)                       :: f1,f2,f3,sum
    integer                       :: j

    if (sw(csfxmod)) then      !Apply switch again to ensure modulation is turned off
      !f1 = parmset%beta(csfxmod,iz) * gf(csfxmod)
      !f1 = parmset%beta(csfxmod,iz) * gf(csfxmod) + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) )/ parmset%beta(0,iz)  !MSISE-00 version; temperature solar flux terms must be normalized
      f1 = parmset%beta(csfxmod,iz) * gf(csfxmod) + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f1 = 0.0d0
    endif

    if (sw(csfxmod+1)) then      !Apply switch again to ensure modulation is turned off
      !f2 = parmset%beta(csfxmod+1,iz) * gf(csfxmod)
      !f2 = parmset%beta(csfxmod+1,iz) * gf(csfxmod) + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) )/ parmset%beta(0,iz)  !MSISE-00 version; temperature solar flux terms must be normalized
      f2 = parmset%beta(csfxmod+1,iz) * gf(csfxmod) + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f2 = 0.0d0
    endif

    if (sw(csfxmod+2)) then      !Apply switch again to ensure modulation is turned off
      f3 = parmset%beta(csfxmod+2,iz) * gf(csfxmod)
    else
      f3 = 0.0d0
    endif

    sum = 0.0
    do j = 0, mbf

!        if (.not. parmset%active(j,iz)) cycle

        !Solar flux modulation of zonal mean intra-annual terms
        if (zsfx(j)) then
            sum = sum + parmset%beta(j,iz)*gf(j)*f1
            cycle
        endif

        !Solar flux modulation of tides
        if (tsfx(j)) then
            sum = sum + parmset%beta(j,iz)*gf(j)*f2
            cycle
        endif
        
        !Solar flux modulation of longitude dependence
        if (psfx(j)) then
            sum = sum + parmset%beta(j,iz)*gf(j)*f3
            cycle
        endif

    enddo

    sfluxmod = sum

    return

end function sfluxmod

! ================================================================
! MSIS legacy ap dependence (daily ap mode and ap history mode)
! including mixed ap/UT/Longitude terms.
! ================================================================
! From the orginal NRLMSISE-00
!
!     t(9)=apdf*(p(33)+p(46)*plg(3,1)+p(35)*plg(5,1)+
!      (p(101)*plg(2,1)+p(102)*plg(4,1)+p(103)*plg(6,1))*cd14*swc(5)+
!      (p(122)*plg(2,2)+p(123)*plg(4,2)+p(124)*plg(6,2))*swc(7)*
!      cos(hr*(tloc-p(125))))
!    t(13)= apdf*swc(11)*(1.+p(121)*plg(2,1))*
!     ((p(61)*plg(3,2)+p(62)*plg(5,2)+p(63)*plg(7,2))*
!          cos(dgtr*(long-p(64))))
!      +apdf*swc(11)*swc(5)*
!      (p(116)*plg(2,2)+p(117)*plg(4,2)+p(118)*plg(6,2))*
!          cd14*cos(dgtr*(long-p(119)))
!      + apdf*swc(12)*
!      (p(84)*plg(2,1)+p(85)*plg(4,1)+p(86)*plg(6,1))*
!          cos(sr*(sec-p(76)))
!
! where the parameters have been reindexed as
!
!    bf(cmag:cmag+6) = ap - 4.0       !ap array
!    bf(cmag+8) = doy2rad*doy
!    bf(cmag+9) = lst2rad*lst
!    bf(cmag+10) = deg2rad*lon
!    bf(cmag+11) = lst2rad*utsec/3600
!    bf(cmag+12) = abs(lat)
!    bf(cmag+13) = plg(0,0)
!    bf(cmag+14) = plg(1,0)
!    ...
!    bf(cmag+26) = plg(6,1)
!
! Master switch control is as follows
!   sw(cmag) .nor. sw(cmag+1)   Do nothing: Return zero
!   sw(cmag) .and. sw(cmag+1)   Daily Ap mode
!   sw(cmag) .neqv. sw(cmag+1)  3-hour ap history mode

real(8) function geomag(p0,bf,plg)

    implicit none

    real(8),intent(in)  :: p0(0:nmag-1)
    real(8),intent(in)  :: bf(0:12)
    real(8),intent(in)  :: plg(0:6,0:1)

    real(8)             :: p(0:nmag-1)    !Copy of parameters used to apply switches
    logical             :: sw1(0:nmag-1)  !Copy of switches
    real(8)             :: delA, gbeta, ex, sumex, G(1:6)
    integer(4)          :: i
    
    !Return zero if both master switches are off    
    if (.not. (sw(cmag) .or. sw(cmag+1))) then
        geomag = 0
        return
    endif

    !Copy parameters
    p = p0
    sw1 = sw(cmag:cmag+nmag-1)

    !Daily Ap mode
    if (sw1(0) .eqv. sw1(1)) then
        where(.not. sw1(2:25)) p(2:25) = 0.0  !Apply switches
        p(8) = p0(8)       !Need doy phase term
        delA = G0fn(bf(0),p(0),p(1))
        geomag = ( p(2)*plg(0,0) +  p(3)*plg(2,0) +  p(4)*plg(4,0)                    &  ! time independent
          + (p(5)*plg(1,0) + p(6)*plg(3,0) + p(7)*plg(5,0)) * cos(bf(8) - p(8))       &  ! doy modulation
          + (p(9)*plg(1,1) + p(10)*plg(3,1) + p(11)*plg(5,1)) * cos(bf(9) - p(12))    &  ! local time modulation
          + (1.0d0 + p(13)*plg(1,0)) *                                                &
            (p(14)*plg(2,1) + p(15)*plg(4,1) + p(16)*plg(6,1)) * cos(bf(10) - p(17))   &  ! longitude effect
          + (p(18)*plg(1,1) + p(19)*plg(3,1) + p(20)*plg(5,1)) * cos(bf(10) - p(21)) * &
            cos(bf(8) - p(8))                                                         &  ! longitude with doy modulaiton
          + (p(22)*plg(1,0) + p(23)*plg(3,0) + p(24)*plg(5,0)) * cos(bf(11) - p(25)) ) &  ! universal time
          *delA
    !3-hour ap history mode
    else
        where(.not. sw1(30:)) p(30:) = 0.0    !Apply switches
        p(36) = p0(36)       !Need doy phase term
        gbeta = p(28)/(1 + p(29)*(45d0 - bf(12)))
        ex = exp(-10800d0*gbeta)
        sumex = 1 + (1 - ex**19d0) * ex**(0.5d0) / (1 - ex)
        do i = 1, 6
            G(i) = G0fn(bf(i),p(26),p(27))
        enddo
        delA = ( G(1)                                                                &
                     + ( G(2)*ex + G(3)*ex*ex + G(4)*ex**3d0                         &
                          +(G(5)*ex**4d0 + G(6)*ex**12d0)*(1-ex**8d0)/(1-ex) ) ) / sumex
        geomag = ( p(30)*plg(0,0) +  p(31)*plg(2,0) +  p(32)*plg(4,0)                    &  ! time independent
          + (p(33)*plg(1,0) + p(34)*plg(3,0) + p(35)*plg(5,0)) * cos(bf(8) - p(36))       &  ! doy modulation
          + (p(37)*plg(1,1) + p(38)*plg(3,1) + p(39)*plg(5,1)) * cos(bf(9) - p(40))    &  ! local time modulation
          + (1.0d0 + p(41)*plg(1,0)) *                                                &
            (p(42)*plg(2,1) + p(43)*plg(4,1) + p(44)*plg(6,1)) * cos(bf(10) - p(45))   &  ! longitude effect
          + (p(46)*plg(1,1) + p(47)*plg(3,1) + p(48)*plg(5,1)) * cos(bf(10) - p(49)) * &
            cos(bf(8) - p(36))                                                         &  ! longitude with doy modulaiton
          + (p(50)*plg(1,0) + p(51)*plg(3,0) + p(52)*plg(5,0)) * cos(bf(11) - p(53)) ) &  ! universal time
          *delA
    endif

    return

    contains

    real(8) function G0fn(a,k00r,k00s)
        real(8),intent(in)  :: a, k00r, k00s
        G0fn = a + (k00r - 1.d0) * (a + (exp(-a*k00s) - 1.d0)/k00s)
        return
    end function G0fn

  end function geomag


! ================================================================
! MSIS legacy UT dependence
! ================================================================
! From the orginal NRLMSISE-00
!    T(12)=(1.+P(96)*PLG(2,1))*(1.+P(82)*DFA*SWC(1))*
!   $(1.+P(120)*PLG(2,1)*SWC(5)*CD14)*
!   $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
!   $     COS(SR*(SEC-P(72))))
!    T(12)=T(12)+SWC(11)*
!   $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
!   $     COS(SR*(SEC-P(80))+2.*DGTR*LONG)*(1.+P(138)*DFA*SWC(1))
!
! The basis functions have been reindexed:
    !bf(cut) =   lst2rad*utsec/3600.0
    !bf(cut+1) = doy2rad*doy
    !bf(cut+2) = dfa
    !bf(cut+3) = deg2rad*lon
    !bf(cut+4) = plg(1,0)
    !bf(cut+5) = plg(3,0)
    !bf(cut+6) = plg(5,0)
    !bf(cut+7) = plg(3,2)
    !bf(cut+8) = plg(5,2)
! as have the model parameters

real(8) function utdep(p0,bf)

    implicit none

    real(8),intent(in)  :: p0(0:nut-1)
    real(8),intent(in)  :: bf(0:8)

    real(8)             :: p(0:nut-1)    !Copy of parameters used to apply switches
    logical             :: sw1(0:nut-1)  !Copy of switches

    !Copy parameters
    p = p0
    sw1 = sw(cut:cut+nut-1)
    where(.not. sw1(3:nut-1)) p(3:nut-1) = 0.0  !Apply switches

    utdep = cos(bf(0)-p(0)) *                          &
            (1 + p(3)*bf(4)*cos(bf(1)-p(1)))  *        &
            (1 + p(4)*bf(2)) * (1 + p(5)*bf(4)) *      &
            (p(6)*bf(4) + p(7)*bf(5) + p(8)*bf(6)) +   &
            cos(bf(0)-p(2)+2*bf(3)) * (p(9)*bf(7) + p(10)*bf(8)) * (1 + p(11)*bf(2))

    return

  end function utdep


  ! --------------------------------------------------
  ! Solar zenith angle approximation
  ! --------------------------------------------------

  real(8) function solzen(ddd,lst,lat,lon)

      implicit none
      real(8),intent(in)   :: ddd
      real(8),intent(in)   :: lst
      real(8),intent(in)   :: lat
      real(8),intent(in)   :: lon

      real(8)              :: wlon,dec
      real(8)              :: teqnx,tf,teqt
      real(8)              :: rlat,phi,cosx

      real(8),parameter    :: pi = 3.141592653589793
      real(8),parameter    :: deg2rad = 0.017453292519943295
      real(8),parameter    :: humr = pi/12.0
      real(8),parameter    :: dumr = pi/182.5
      real(8),parameter    :: p(5) = &
        (/0.017203534,0.034407068,0.051610602,0.068814136,0.103221204/)

      wlon = 360.0 - lon
      teqnx = ddd + (lst + wlon / 15.0) / 24.0 + 0.9369
      teqnx = ddd + 0.9369

    ! declination

      dec = 23.256 * sin(p(1) * (teqnx - 82.242)) + 0.381 * sin(p(2)*(teqnx - 44.855))  &
           + 0.167 * sin(p(3) * (teqnx - 23.355)) - 0.013 * sin(p(4)*(teqnx + 11.97)) &
           + 0.011 * sin(p(5) * (teqnx - 10.410)) + 0.339137
      dec = dec * deg2rad

      !equation of time

      tf = teqnx - 0.5
      teqt = -7.38 * sin(p(1) * (tf -  4.0)) - 9.87 * sin(p(2) * (tf +  9.0)) &
            + 0.27 * sin(p(3) * (tf - 53.0)) -  0.2 * cos(p(4) * (tf - 17.0))

      phi = humr * ( lst - 12.0) + teqt * deg2rad / 4.0
      rlat = lat * deg2rad

      cosx = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(phi)
      if (abs(cosx) .gt. 1.0) cosx = sign(1.0d0,cosx)

      solzen = acos(cosx) / deg2rad
      return

    end function solzen

end module msis

!!! =====================================================
!!! Returns array of nonzero b-spline values, for all
!!! orders up to specified order (maximum 6)
!!! =====================================================

subroutine bspline(x,nodes,nd,kmax,eta,S,i)
    implicit none

    real(8),intent(in)        :: x           !Location at which splines are to be evaluated
    real(8),intent(in)        :: nodes(0:30) !Spline node locations
    integer,intent(in)        :: nd          !Number of spline nodes minus one (0:nd)
    integer,intent(in)        :: kmax        !Maximum order (up to 6 allowed) of evaluated splines
	real(8),intent(in)  	  :: eta(0:30,2:6)  !Array of precomputed weights for recursion (reciprocals of node differences)
	real(8),intent(out)       :: S(-5:0,2:6) !Array of b-spline values (spline index relative to i (-5:0), spline order (2:6))
    integer,intent(out)       :: i           !Index of last nonzero bspline

    integer				      :: j,k,l
    integer       		      :: low,high
    real(8)                   :: w(-4:0) !Weights for recursion relation
 
	S(:,:) = 0.0d0
	
	!Find index of last (rightmost) nonzero spline
	if (x .ge. nodes(nd)) then
		i = nd
		return
	endif
	if (x .le. nodes(0)) then
		i = -1
		return
	endif
    low = 0
    high = nd
    i = (low + high)/2
    do while (x .lt. nodes(i) .or. x .ge. nodes(i + 1))
       if (x .lt. nodes(i)) then
          high = i
       else
          low = i
       endif
       i = (low + high)/2
    end do
 	
	!Initialize with linear splines
	S(0,2) = (x - nodes(i)) * eta(i,2)
	if (i .gt. 0) S(-1,2) = 1 - S(0,2)
    if (i .ge. nd-1) S(0,2) = 0d0   !Reset out-of-bounds spline to zero
	
	do k = 3, kmax		!Loop through orders (starting with linear (k=2) splines
		!Calculate recursion weights
		w(:) = 0.0d0
		do l = 0, -(k-2), -1
			j = i + l
			if (j .lt. 0) cycle    !Skip out-of-bounds splines
			w(l) = (x - nodes(j)) * eta(j,k)
		enddo
		!Calculate splines
		if (i .lt. (nd-k+1)) S(0,k) = w(0)*S(0,k-1)	    !First spline (from right)
		do l = -1, -(k-2), -1   !Loop through middle splines, right to left
			if ( ((i+l) .lt. 0) .or. ((i+l) .ge. (nd-k+1)) ) cycle
			S(l,k) = w(l)*S(l,k-1) + (1-w(l+1))*S(l+1,k-1)
		enddo
		if ((i-k+1) .ge. 0) S(-k+1,k) = (1-w(-k+2))*S(-k+2,k-1) !Last spline
	enddo
	
	return
	
  end subroutine bspline

!!! =====================================================
!!! Function to calculate dilogarithm in the domain [0,1), 
!!! using the algorithm described by Ginsberg, E. S., and 
!!! D. Zaborowski (1975), The Dilogarithm function of a real
!!! argument, Commun. ACM, 18, 200–202. Retains terms up
!!! to order 3 in the expansion, which results in relative
!!! errors less than 1E-5. JTE, 3/7/17, 12/12/18
!!! =====================================================

  real(8) function dilog(x0)
    implicit none
    real(8),intent(in)        :: x0
    real(8)                   :: x, xx, x4, lnx
    real(8),parameter    :: Pi = 3.1415926535897932384626433832795d0
    real(8),parameter    :: Pi2_6 = Pi*Pi / 6d0
    x = x0
    if (x .gt. 0.5d0) then
      lnx = log(x)
      x = 1d0 - x          !Reflect argument into [0,0.5] range
      xx = x*x
      x4 = 4d0*x
      dilog = pi2_6 - lnx*log(x) &
              - (4d0*xx*(23d0/16d0 + x/36d0 + xx/576d0 + xx*x/3600d0) &
                 + x4 + 3d0*(1d0 - xx)*lnx) / (1d0 + x4 + xx)
    else
      xx = x*x
      x4 = 4d0*x
      dilog = (4d0*xx*(23d0/16d0 + x/36d0 + xx/576d0 + xx*x/3600d0) &
               + x4 + 3d0*(1d0 - xx)*log(1d0 - x)) / (1d0 + x4 + xx)
    endif

    return
  end function dilog

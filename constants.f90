!!!
!!! Authors:  Douglas P. Drob and John Emmert, NRL 7632
!!!
!!! =================================================
!!! MSIS constants
!!! =================================================

module constants

    implicit none

    ! -----------------------------------------------------------------------------
    ! Trigonometric constants
    ! -----------------------------------------------------------------------------
    real(8),parameter    :: Pi = 3.1415926535897932384626433832795d0
    real(8),parameter    :: deg2rad = Pi / 180d0
    real(8),parameter    :: doy2rad = 2d0*Pi / 365d0
    real(8),parameter    :: lst2rad = Pi / 12d0
    real(8),parameter    :: tanh1 = 0.761594155955765d0  ! tanh(1.0)

    ! -----------------------------------------------------------------------------
    ! Thermodynamic constants (mks units, CODATA 2014 and Mohr et al., [2012])
    ! -----------------------------------------------------------------------------
                                            !      mbar  Mass      N2         O2            O          He     H    Ar                N       Anom O                   NO
    real(8),parameter      :: amu(0:10) =    (/28.96546,  0.0, 28.0134,   31.9988, 31.9988/2.0,     4.0026, 1.0,   39.948, 28.0134/2.0, 31.9988/2.0, (28.0134+31.9988)/2/)
    real(8),parameter      :: volmix(0:10) = (/     1.0,  0.0, 0.780848, 0.209390,         1.0, 0.00000524, 1.0, 0.009332,         1.0,         1.0,                 1.0/)
    real(8),parameter      :: boltz = 1.38064852d-23   ! mks units
    real(8),parameter      :: amu2kg = 1.660539040d-27
    real(8),parameter      :: mbar = 28.96546d0
    real(8),parameter      :: rgas = 8.3144598
    real(8),parameter      :: g0 = 9.7976432222d0
    real(8),parameter      :: lnP0 = 11.5080482  ! Natural log of global average surface pressure (hardwired for now)
    real(8),parameter      :: g0divR = g0/rgas 
    real(8),parameter      :: fh = mbar*g0/rgas
    real(8),parameter      :: fho1 = 31.9988/2.0*g0/rgas

    !! -----------------------------------------------------------------------------
    !! Thermodynamic constants (mks units, NRLMSISE-00 Values)
    !! -----------------------------------------------------------------------------
    !                                        !      mbar  Mass      N2         O2            O          He     H    Ar                N       Anom O                   NO
    !real(8),parameter      :: amu(0:10) =    (/28.9502,  0.0,      28.0,     32.0,        16.0,           4.0, 1.0,       40.0,        14.0,        16.0,                30.0/)
    !real(8),parameter      :: volmix(0:10) = (/     1.0,  0.0, 0.781105, 0.209547,         1.0, 0.00000524177, 1.0, 0.00934318,         1.0,         1.0,                 1.0/)
    !real(8),parameter      :: boltz = 1.3806d-23
    !real(8),parameter      :: amu2kg = 1.66d-27
    !real(8),parameter      :: mbar = 28.9502d0
    !real(8),parameter      :: rgas = 8.314d0
    !real(8),parameter      :: g0 = 9.7976432222d0
    !real(8),parameter      :: lnp0 = 11.5390  ! Natural log of global average surface pressure
    !real(8),parameter      :: fh = mbar*g0/rgas
    !real(8),parameter      :: fho1 = 32.0/2.0*g0/rgas

    ! --------------------------------------------
    ! Vertical profile parameters
    ! --------------------------------------------

    integer,parameter      :: nd = 27    !Number of temperature profile nodes
    integer,parameter      :: p = 4      !Spline order
    integer,parameter      :: nl = nd - p  !Last temperature profile level index
    integer,parameter      :: nls = 9    ! Last parameter index for each species (excluding O, NO splines)

    integer,parameter      :: nspec = 11  !Number of species including temperature
    integer,parameter      :: ngrp = 1024 !Maximum number of group IDs.
    
    real(8),parameter      :: nodes(0:nd+2) = &  !Nodes for temperature profile splines
        (/ -15., -10.,  -5.,   0.,   5., 10., 15., 20.,  25.,  30.,  35.,  40., 45., 50., &
            55.,  60.,  65.,  70.,  75., 80., 85., 92.5, 102.5, 112.5, 122.5, 132.5, 142.5, &
           152.5, 162.5, 172.5/)

    real(8),parameter      :: bwalt = 122.5d0     ! Reference geopotential height for Bates Profile
    real(8),parameter      :: zetaF = 70.0d0      ! Fully mixed below this, uses constant mixing ratios
    real(8),parameter      :: zetaB = bwalt       ! Bates Profile above this altitude
    real(8),parameter      :: zetaA = 85.0d0      ! Default reference height for active minor species
    
    integer,parameter      :: izfmx = 13          ! fully mixed below this spline index
    integer,parameter      :: izax = 17           ! Spline index at zetaA
    integer,parameter      :: izfx = 14           ! Spline index at zeta F
    integer,parameter      :: itex = nl           ! Index of Bates exospheric temperature
    integer,parameter      :: itgb0 = nl - 1      ! Index of Bates temperature gradient at lower boundary
    integer,parameter      :: itb0 = nl - 2       ! Index of Bates temperature at lower boundary
    
    real(8),parameter      :: nodesO1(0:13) = &   !Nodes for O1 splines (Domain 50-85 km)
        (/ 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 92.5, 102.5, 112.5/)
    integer,parameter      :: ndO1 = 13
    integer,parameter      :: nsplO1 = ndO1-5     !Number of unconstrained spline parameters for O1 (there are 2 additional C1-constrained splines)
    real(8),parameter      :: zetarefO1 = zetaA      !Joining height for O1 splines, and reference height for O1 density
    
    real(8),parameter      :: nodesNO(0:13) = &   !Nodes for NO splines (Domain 70-122.5 km)
        (/ 47.5, 55., 62.5, 70., 77.5, 85., 92.5, 100., 107.5, 115., 122.5, 130., 137.5, 145./)
    integer,parameter      :: ndNO = 13
    integer,parameter      :: nsplNO = ndNO-5     !Number of unconstrained spline parameters for NO (there are 2 additional C1-constrained splines)
    real(8),parameter      :: zetarefNO = zetaB      !Joining height for NO splines, and reference height for NO density
    
    ! Matrix for matching b-spline boundary conditions
    real(8),parameter      :: c2tn(3,3) = &
        reshape((/1.0d0, -10.d0,  33.333333333333336d0, &
                  1.0d0,   0.d0, -16.666666666666668d0, &
                  1.0d0,  10.d0,  33.333333333333336d0/), &
                (/3,3/))

    !C1 Continuity for O; Last 2 splines are constrained
    real(8),parameter      :: c1o1(2,2) = &
        reshape((/ 1.75d0,               -2.916666573405061d0, &
                  -1.624999900076852d0,  21.458332647194382d0 /), &
                (/2,2/))
    real(8),parameter      :: c1o1adj(2) = (/0.257142857142857d0, -0.102857142686844d0/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
    
    !C1 Continuity for NO; Last 2 splines are constrained
    real(8),parameter      :: c1NO(2,2) = &
        reshape((/ 1.5d0,  -3.75d0, &
                   0.0d0,  15.0d0 /), &
                (/2,2/))
    real(8),parameter      :: c1NOadj(2) = (/0.166666666666667d0, -0.066666666666667d0/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
    
    ! Anomalous Oxygen parameters (legacy profile from NRLMSISE-00
    real(8),parameter      :: zetarefOA = zetaB   !Reference height for anomalous oxygen density
    real(8),parameter      :: TOA = 4000.         !Temperature of anomalous oxygen density (K)
    real(8),parameter      :: HOA = boltz * TOA / (16.0*amu2kg*g0) / 1000.  !Hydrostatic scale height of anomalous oxygen density (km)
    

    ! --------------------------------------------
    ! Horizontal basis function parameters
    ! --------------------------------------------

    integer,parameter      :: maxnbf = 512   ! Number of basis functions to be allocated

    integer,parameter      :: maxn = 6       ! Maximum latitude (Legendre) spectral degree
    integer,parameter      :: maxl = 3       ! Maximum local time (tidal) spectral order
    integer,parameter      :: maxm = 2       ! Maximum longitude (stationary planetary wave) order
    integer,parameter      :: maxs = 2       ! Maximimum day of year (intra-annual) Fourier order

    integer,parameter      :: amaxn = 6      ! Maximum Legendre degree used in time independent and intra-annual zonal mean terms
    integer,parameter      :: amaxs = 2      ! Maximum intra-annual order used in zonal mean terms

    integer,parameter      :: tmaxl = 3      ! Maximum tidal order used
    integer,parameter      :: tmaxn = 6      ! Maximum Legendre degree coupled with tides
    integer,parameter      :: tmaxs = 2      ! Maximum intra-annual order coupled with tides

    integer,parameter      :: pmaxm = 2      ! Maximum stationary planetary wave order used
    integer,parameter      :: pmaxn = 6      ! Maximum Legendre degree coupled with SPW
    integer,parameter      :: pmaxs = 2      ! Maximum intra-annual order coupled with SPW

    integer,parameter      :: nsfx = 5       ! Number of linear solar flux terms
    integer,parameter      :: nsfxmod = 5    ! Number of nonlinear modulating solar flux terms (legacy NRLMSISE-00 terms)
    
    integer,parameter      :: nmag = 54      ! Number of terms in NRLMSISE-00 legacy geomagnetic parameterization
    integer,parameter      :: nut = 12       ! Number of terms in NRLMSISE-00 legacy UT parameterization

    integer,parameter      :: ctimeind = 0   ! Starting index of time-independent terms
    integer,parameter      :: cintann = ctimeind + (amaxn+1)   ! Starting index of zonal mean intra-annual terms
    integer,parameter      :: ctide = cintann + ((amaxn+1)*2*amaxs)   ! Starting index of zonal mean intra-annual terms
    integer,parameter      :: cspw = ctide + (4*tmaxs+2)*(tmaxl*(tmaxn+1)-(tmaxl*(tmaxl+1))/2) ! Starting index of SPW terms
    integer,parameter      :: csfx = cspw + (4*pmaxs+2)*(pmaxm*(pmaxn+1)-(pmaxm*(pmaxm+1))/2)   ! Starting index of linear solar flux terms
    integer,parameter      :: cextra = csfx + nsfx   ! Starting index of time-independent terms
    integer,parameter      :: mbf = 383           ! Last index of linear terms
    integer,parameter      :: cnonlin = mbf + 1   ! Starting index of nonlinear terms
    integer,parameter      :: csfxmod = cnonlin   ! Starting index of modulating solar flux terms
    integer,parameter      :: cmag = csfxmod + nsfxmod    ! Starting index of daily geomagnetic terms
    integer,parameter      :: cut = cmag + nmag    ! Starting index of UT terms
    
    ! --------------------------------------------------------------------------------
    ! Weights for pressure coefficient calculation (from temperature coefficients)
    ! --------------------------------------------------------------------------------
    
    real(8),parameter       :: gwht(0:3) = &
           (/ 5.0d0/24.0d0, 55.0d0/24.0d0, &
             55.0d0/24.0d0, 5.0d0/24.0d0/)

    ! --------------------------------------------------------------------------------
    ! Constants need for analytical integration by parts of piecewise mass profile
    ! --------------------------------------------------------------------------------

    real(8),parameter       :: wbeta(0:nl) = (nodes(4:nd)-nodes(0:nl)) / 4d0 !Weights for 1st spline integration
    real(8),parameter       :: wgamma(0:nl) = (nodes(5:nd+1)-nodes(0:nl)) / 5d0 !Weights for 2nd spline integration
    
    !Non-zero bspline values at zetaB (5th and 6th order)
    real(8),parameter       :: S5zetaB(0:3) = (/0.041666666666667d0, 0.458333333333333d0, 0.458333333333333d0, &
                                                0.041666666666667d0/)
    real(8),parameter       :: S6zetaB(0:4) = (/0.008771929824561d0, 0.216228070175439d0, 0.550000000000000d0, &
                                                0.216666666666667d0, 0.008333333333333d0/)

    !Weights for calculating temperature gradient at zetaA
    real(8),parameter       :: wghtAxdz(0:2) = (/-0.102857142857d0, 0.0495238095238d0, 0.053333333333d0/)

    !Non-zero bspline values at zetaA (4th, 5th and 6th order)
    real(8),parameter       :: S4zetaA(0:2) = (/0.257142857142857d0, 0.653968253968254d0, 0.088888888888889d0/)
    real(8),parameter       :: S5zetaA(0:3) = (/0.085714285714286d0, 0.587590187590188d0, 0.313020313020313d0, &
                                                0.013675213675214d0/)
    real(8),parameter       :: S6zetaA(0:4) = (/0.023376623376623d0, 0.378732378732379d0, 0.500743700743701d0, &
                                                0.095538448479625d0, 0.001608848667672d0/)

    !Non-zero bspline values at zetaF (4th and 5th order)
    real(8),parameter       :: S4zetaF(0:2) = (/0.166666666666667d0, 0.666666666666667d0, 0.166666666666667d0/)
    real(8),parameter       :: S5zetaF(0:3) = (/0.041666666666667d0, 0.458333333333333d0, 0.458333333333333d0, &
                                                0.041666666666667d0/)

    !Non-zero bspline values at zeta=0 (5th order)
    real(8),parameter       :: S5zeta0(0:2) = (/0.458333333333333d0, 0.458333333333333d0, 0.041666666666667d0/)

    end module constants

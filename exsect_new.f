C Subroutine EXSECT
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Adapted from Banks & Nagy 2-stream code by Stan Solomon, 1988
C Added high-energy relativistic cross section correction, SCS, 1999
C Updated comments, SCS, 2002
C Included in GLOW v. 0.97, SCS, 2005
C
C Calculates electron impact cross sections
C
C Definitions:
C SIGS   elastic cross sections for each species, energy; cm2
C PE     elastic backscatter probabilities for each species, energy
C PI     inelastic  "
C SIGA   energy loss cross section for each species, loss, energy; cm2
C SEC    secondary production xsect for species, Esec, Epri; cm2
C SIGEX  excitation xsect for each state, species, energy; cm2
C     O states:   1D,   1S, 3s5S, 3s3S, 3p5P, 3p3P, 3d3D, 3s'3D
C     O2 states:   a,    b, AA'c,    B,  9.9, Ryds,  vib
C     N2 states: ABW,   B',    C, aa'w,  1Pu,   b', Ryds,  vib
C SIGIX  ionization xsect for each state, species, energy; cm2
C     O states:   4S,  2Do,  2Po
C     O2 states:   X,    a,    A,    b,    B,   c,  37eV
C     N2 states:   X,    A,    B,    D,    C, 40eV
C IIMAX  number of bins for secondary production for each primary energy
C WW     energy threshold for each excited state, species; eV
C WW, AO, OMEG, ANU, BB: revised excitation cross section parameters,
C        from Green & Stolarski (1972) formula (W, A, omega, nu, gamma)
C AUTO   autoionization coefs (= 0 as autoion. included in ion xsects)
C THI    energy threshold for each ionized state, species; eV
C AK, AJ, TS, TA, TB, GAMS, GAMB:  Jackman et al (1977) ioniz. params
C ENER   energy grid; eV
C DEL    energy grid spacing; eV
C NNN    number of excited states for each species
C NINN   number of ionized states for each species
C NUM    number of points on elastic data trid for each species
C EC     data energy grid of elastic xsects and backscatter ratios
C        for each species; eV
C CC     elastic xsects on data grid for each species, cm2
C CE     elastic backscat. probs on data grid for each species; cm2
C CI     inelastic "
C
C Array dimensions:
C NBINS  number of energy levels
C NMAJ   number of major species
C NEI    number of slots for excited and ionized states
C
C
      SUBROUTINE EXSECT (ENER, DEL)
C
      INCLUDE 'glow.h'
      PARAMETER (NMAJ=3)
      PARAMETER (NEI=10)
C  Added NNEW. 11/20/11, jy
      PARAMETER (NNEW=12)
C  
      COMMON /CXSECT/ SIGS(NMAJ,NBINS), PE(NMAJ,NBINS), PI(NMAJ,NBINS),
     >                SIGA(NMAJ,NBINS,NBINS), SEC(NMAJ,NBINS,NBINS),
     >                SIGEX(NEI,NMAJ,NBINS), SIGIX(NEI,NMAJ,NBINS),
     >                IIMAXX(NBINS)
C
      COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ),
     >                ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ),
     >                THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ),
     >                TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ),
     >                GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)
c&&
      COMMON /CXNEW/ WWN(NNEW),COA(NNEW),COB(NNEW),COC(NNEW),
     >               COD(NNEW),COE(NNEW),COF(NNEW),COG(NNEW),
     >               COH(NNEW),SIGNE(NNEW,NBINS)


C
      DIMENSION ENER(NBINS), DEL(NBINS), SIGI(NBINS), T12(NBINS),
     >          RATIO(NBINS), NNN(NMAJ), NINN(NMAJ), NUM(NMAJ),
     >          EC(31,NMAJ), CC(31,NMAJ), CE(31,NMAJ), CI(31,NMAJ)
C
      DATA QQN/6.51E-14/, NNN/8,7,8/, NINN/3,7,6/, NUM/31,28,28/
C
      DATA WW  /1.96, 4.17, 9.29, 9.53,10.76,10.97,12.07,12.54, 0.,0.,
     >          0.98, 1.64, 4.50, 8.44, 9.90,13.50, 0.25, 0.00, 0.,0.,
     >          6.17, 8.16,11.03, 8.40,12.85,14.00,13.75, 1.85, 0.,0./
      DATA AO /.0100,.0042,.1793,.3565,.0817,.0245,.0293,.1221, 0.,0.,
     >         .0797,.0211,.0215,.3400,.0657,1.110,3.480, 0.00, 0.,0.,
     >         2.770,.1140,.1790,.0999,.8760,.6010,1.890,1.350, 0.,0./
      DATA OMEG/1.00, 1.00, 3.00, 0.75, 3.00, 0.85, 0.75, 0.75, 0.,0.,
     >          2.00, 2.00, 1.15, 0.75, 0.75, 0.75, 7.00, 0.00, 0.,0.,
     >          3.00, 3.00, 3.00, 1.00, 0.75, 0.75, 0.75, 8.00, 0.,0./
      DATA ANU /2.00, 1.04, 2.53, 0.54, 2.43, 2.87, 0.93, 0.72, 0.,0.,
     >          6.18, 4.14, 1.00, 1.05, 1.60, 3.00,10.87, 0.00, 0.,0.,
     >          4.53, 4.78, 4.32, 4.05, 1.47, 1.27, 3.00, 1.58, 0.,0./
      DATA BB / 1.00, 0.50, 1.02, 0.01, 4.19, 4.88, 0.66, 0.17, 0.,0.,
     >          0.53, 0.51, 0.98, 0.99, 1.86, 1.00, 1.00, 0.00, 0.,0.,
     >          1.42, 3.54,12.70, 5.20, 0.86, 0.45, 1.00, 1.00, 0.,0./
      DATA AUTO/30*0.0/
      DATA THI/13.60,16.90,18.50, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >         12.10,16.10,16.90,18.20,20.00,23.00,37.00, 0.,0.,0.,
     >         15.58,16.73,18.75,22.00,23.60,40.00, 0.00, 0.,0.,0./
      DATA AK / 1.13, 1.25, 0.67, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >          0.47, 1.13, 1.13, 1.01, 0.65, 0.95, 0.59, 0.,0.,0.,
     >          2.42, 1.06, 0.55, 0.37, 0.37, 0.53, 0.00, 0.,0.,0./
      DATA AJ / 1.81, 1.79, 1.78, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >          3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 0.,0.,0.,
     >          1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 0.00, 0.,0.,0./
      DATA TS / 6.41, 6.41, 6.41, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >          1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 0.,0.,0.,
     >          4.71, 4.71, 4.71, 4.71, 4.71, 4.71, 0.00, 0.,0.,0./
      DATA TA /3450.,3450.,3450.,   0.,   0.,   0.,   0., 0.,0.,0.,
     >         1000.,1000.,1000.,1000.,1000.,1000.,1000., 0.,0.,0.,
     >         1000.,1000.,1000.,1000.,1000.,1000.,   0., 0.,0.,0./
      DATA TB/162.00,162.0,162.0, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >         24.20,32.20,33.80,36.40,40.60,46.00,74.00, 0.,0.,0.,
     >         31.16,33.46,37.50,44.00,47.20,80.00, 0.00, 0.,0.,0./
      DATA GAMS/13.00,13.0,13.00, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >          18.50,18.5,18.50,18.50,18.50,18.50,18.50, 0.,0.,0.,
     >          13.80,13.8,13.80,13.80,13.80,13.80, 0.00, 0.,0.,0./
      DATA GAMB/-.815,-.815,-.815,0.00, 0.00, 0.00, 0.00, 0.,0.,0.,
     >          12.10,16.10,16.90,18.2,20.30,23.00,37.00, 0.,0.,0.,
     >          15.58,16.73,18.75,22.0,23.60,40.00, 0.00, 0.,0.,0./
      DATA EC /      1.00,     2.00,     4.00,     6.00,     8.00,
     >              10.00,    12.00,    14.00,    16.00,    18.00,
     >              20.00,    30.00,    40.00,    50.00,    60.00,
     >              70.00,    80.00,    90.00,   100.00,   150.00,
     >             200.00,   300.00,   500.00,  1000.00,  2000.00,
     >            3000.00,  5000.00, 10000.00, 20000.00, 40000.00,
     >           50000.00,
     >               1.00,     2.00,     3.00,     5.00,     7.00,
     >              10.00,    15.00,    20.00,    30.00,    40.00,
     >              50.00,    70.00,   100.00,   150.00,   200.00,
     >             300.00,   400.00,   500.00,   600.00,   700.00,
     >            1000.00,  2000.00,  3000.00,  5000.00, 10000.00,
     >           20000.00, 40000.00, 50000.00,     0.00,     0.00,
     >               0.00,
     >               1.00,     2.00,     2.50,     3.00,     4.00,
     >               5.00,     6.00,     8.00,    10.00,    15.00,
     >              20.00,    30.00,    40.00,    50.00,    70.00,
     >             100.00,   200.00,   300.00,   500.00,   700.00,
     >            1000.00,  2000.00,  3000.00,  5000.00, 10000.00,
     >           20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0/
      DATA CC /  5.00E-16, 6.00E-16, 7.50E-16, 7.60E-16, 7.70E-16,
     >           7.80E-16, 7.50E-16, 7.20E-16, 6.90E-16, 6.70E-16,
     >           6.50E-16, 5.60E-16, 4.60E-16, 4.00E-16, 3.50E-16,
     >           3.20E-16, 2.90E-16, 2.70E-16, 2.50E-16, 1.90E-16,
     >           1.50E-16, 1.20E-16, 8.00E-17, 5.00E-17, 3.02E-17,
     >           1.99E-17, 1.20E-17, 6.08E-18, 3.06E-18, 1.55E-18,
     >           1.24E-18,
     >           5.50E-16, 6.90E-16, 7.50E-16, 8.50E-16, 9.60E-16,
     >           1.00E-15, 1.00E-15, 9.00E-16, 8.30E-16, 7.70E-16,
     >           6.90E-16, 5.70E-16, 4.40E-16, 3.30E-16, 2.70E-16,
     >           2.10E-16, 1.80E-16, 1.60E-16, 1.40E-16, 1.30E-16,
     >           1.10E-16, 7.00E-17, 5.00E-17, 3.00E-17, 1.53E-17,
     >           7.72E-18, 3.90E-18, 3.13E-18, 0.00E+00, 0.00E+00,
     >           0.00E+00,
     >           9.00E-16, 2.27E-15, 2.52E-15, 1.93E-15, 1.32E-15,
     >           1.15E-15, 1.16E-15, 1.17E-15, 1.18E-15, 1.14E-15,
     >           1.13E-15, 9.50E-16, 8.60E-16, 7.30E-16, 5.90E-16,
     >           4.70E-16, 3.30E-16, 2.50E-16, 1.60E-16, 1.30E-16,
     >           1.10E-16, 6.35E-17, 4.18E-17, 2.54E-17, 1.28E-17,
     >           6.44E-18, 3.27E-18, 2.62E-18, 0.00E+00, 0.00E+00, 0.0/
      DATA CE /   0.50000,  0.49500,  0.46800,  0.43600,  0.42000,
     >            0.40500,  0.37000,  0.36000,  0.34000,  0.33000,
     >            0.32000,  0.27000,  0.24000,  0.22000,  0.20000,
     >            0.18000,  0.17000,  0.16000,  0.15000,  0.13000,
     >            0.11500,  0.09000,  0.06800,  0.04600,  0.02400,
     >            0.01660,  0.01000,  0.00510,  0.00255,  0.00125,
     >            0.00100,
     >            0.50000,  0.50000,  0.49000,  0.44500,  0.42700,
     >            0.40500,  0.36800,  0.34300,  0.31600,  0.28900,
     >            0.25800,  0.22000,  0.18400,  0.16400,  0.13300,
     >            0.11000,  0.10000,  0.09200,  0.08500,  0.08000,
     >            0.06800,  0.03700,  0.02600,  0.01600,  0.00800,
     >            0.00400,  0.00200,  0.00160,  0.00000,  0.00000,
     >            0.00000,
     >            0.50000,  0.50000,  0.50000,  0.49000,  0.46800,
     >            0.44500,  0.43600,  0.42000,  0.40500,  0.36800,
     >            0.34300,  0.31600,  0.28900,  0.25800,  0.22000,
     >            0.18400,  0.14000,  0.11000,  0.08400,  0.07400,
     >            0.06300,  0.03400,  0.02400,  0.01500,  0.00740,
     >            0.00370,  0.00180,  0.00140,  0.00000,  0.00000, 0.0/
      DATA CI /   0.60000,  0.60000,  0.60000,  0.60000,  0.60000,
     >            0.60000,  0.55000,  0.46000,  0.40000,  0.36000,
     >            0.32000,  0.22000,  0.15000,  0.10000,  0.08200,
     >            0.07000,  0.06100,  0.05400,  0.05000,  0.04400,
     >            0.03800,  0.02800,  0.02000,  0.01050,  0.00600,
     >            0.00400,  0.00250,  0.00130,  0.00060,  0.00030,
     >            0.00025,
     >            0.50000,  0.50000,  0.50000,  0.50000,  0.48000,
     >            0.44000,  0.36000,  0.28000,  0.20000,  0.14000,
     >            0.10000,  0.07000,  0.05000,  0.04600,  0.04300,
     >            0.03700,  0.03200,  0.02800,  0.02400,  0.02100,
     >            0.01600,  0.00900,  0.00620,  0.00400,  0.00200,
     >            0.00100,  0.00050,  0.00040,  0.00000,  0.00000,
     >            0.00000,
     >            0.50000,  0.50000,  0.50000,  0.50000,  0.50000,
     >            0.50000,  0.50000,  0.50000,  0.50000,  0.50000,
     >            0.44000,  0.30000,  0.20000,  0.13000,  0.09000,
     >            0.06000,  0.05000,  0.04200,  0.03200,  0.02500,
     >            0.02000,  0.01100,  0.00800,  0.00500,  0.00250,
     >            0.00120,  0.00060,  0.00050,  0.00000,  0.00000, 0.0/
C ADDED 11/20/11, JY      
      DATA ERYD /13.61/
      
      DATA WWN  /6.14, 7.30, 7.36, 8.16, 
     >    11.03, 11.28,11.52,11.75,
     >    11.97, 8.2,  8.4,  8.9/
    
      DATA COA /9.57530e-16, 148.927e-16,  2.85419e-16, 1.39180e-16,
     >5.24082e-16, 43.9056e-16,  4.72802e-16, 0.182832e-16,
     >15.3774e-16, 0.245721e-16, 74.1954e-16, 26.2727e-16/
      DATA COB /2.41725, 3.37549,  3.05068, 3.80445,
     >2.07034, 3.78085,  2.58944, 1.77218, 
     >3.79464, 2.41297,  5.70981, 6.13567/ 
      DATA COC /3.48174, 2.62388,  6.78624, 6.83396, 
     >2.99729, 2.93550,  2.62034, 3.86083, 
     >1.64863, 8.19281,  4.66609, 5.09806/
      DATA COD /6.42020, 1.48676,  1.31338,  1.35693,  
     >1.31644, 1.31263,  1.62154,  3.50790,   
     >1.39476, 3.44784,  0.938729,  0.806361/ 
      DATA COE /1.05341e-16, 0.0832431e-16, 2.38554e-16, 1.39227e-16,  
     >0.0149362e-16,  0.00963622e-16, 2.77677e-16,  9.43598e-16, 
     >0.00373335e-16, 0.0376054e-16,  0.356180e-16, 0./
      DATA COF /2.06570,  1.89009,   20.1055,   25.9739,  
     >0.164347, 0.0386813, 10.3283,   14.6259,   
     >4.05230,  13.6749,   2.49686,   1.0/
      DATA COG /6.74498,  21.5238,   11.6571,   11.5631,  
     >37.0967,  36.8368,   8.19713,   8.34705,  
     >12.5569,  13.5317,   12.8333,   1.0/
      DATA COH /1.02078,  2.48385,   4.35359,   1.35315, 
     >9.51190,  8.54345,   1.62767,   1.81088, 
     >7.56624,  1.37048,   0.851090,  1.0/

C
C
C Interpolate elastic cross sections and backscatter ratios:
C
      DO 90 IJ=1,NMAJ
        DO  80 IV=1,NBINS
          EX=ENER(IV)
          DO 50 II=1,NUM(IJ)
            IF (EC(II,IJ) .GT. EX) GOTO 60
   50     CONTINUE
          SIGS(IJ,IV)=CC(NUM(IJ),IJ)*(EC(NUM(IJ),IJ)/EX)**0.8
          IF(IJ.EQ.1) SIGS(IJ,IV)=CC(NUM(IJ),IJ)*(EC(NUM(IJ),IJ)/EX)**2
          PE(IJ,IV) = CE(NUM(IJ),IJ)* (EC(NUM(IJ),IJ)/EX)
          PI(IJ,IV) = CI(NUM(IJ),IJ)* (EC(NUM(IJ),IJ)/EX)
          GOTO 80
   60     I=II-1
          IF (I .LE. 0) THEN
            SIGS(IJ,IV)=CC(II,IJ)
            PE(IJ,IV)=CE(II,IJ)
            PI(IJ,IV)=CI(II,IJ)
          ELSE
            FAC = ALOG (EX/EC(I,IJ)) / ALOG (EC(II,IJ)/EC(I,IJ))
            SIGS(IJ,IV) = EXP (ALOG (CC(I,IJ))
     >                         + ALOG (CC(II,IJ)/CC(I,IJ)) * FAC)
            PE(IJ,IV) = EXP (ALOG (CE(I,IJ))
     >                       + ALOG (CE(II,IJ)/CE(I,IJ)) * FAC)
            PI(IJ,IV) = EXP (ALOG (CI(I,IJ))
     >                       + ALOG (CI(II,IJ)/CI(I,IJ)) * FAC)
          ENDIF
   80   CONTINUE
   90 CONTINUE
C
C
C Calculate electron impact excitation and ionization cross sections:
C
      DO 140 I=1,NMAJ
        DO 140 K=1,NEI
          DO 140 J=1,NBINS
          IF (ENER(J).GT.WW(K,I) .AND. WW(K,I).GT.0.001) THEN
            WE = WW(K,I) / ENER(J)
            SIGEX(K,I,J) = QQN * AO(K,I)
     >                       * (WE**OMEG(K,I) / WW(K,I)**2)
     >                       * (1.0 - WE**BB(K,I)) ** ANU(K,I)
            IF (SIGEX(K,I,J) .LT. 1.E-30) SIGEX(K,I,J) = 0.0
          ELSE
            SIGEX(K,I,J) = 0.0
          ENDIF
          IF (ENER(J).GT.THI(K,I) .AND. THI(K,I).GT.0.001) THEN
            AE = AK(K,I)/ENER(J) * ALOG(ENER(J)/AJ(K,I))
            GAMMA = GAMS(K,I) * ENER(J) / (ENER(J)+GAMB(K,I))
            T0 = TS(K,I) - (TA(K,I)/(ENER(J)+TB(K,I)))
            SIGIX(K,I,J) = 1.E-16 * AE * GAMMA
     >                    * ( ATAN(((ENER(J)-THI(K,I))/2.-T0)/GAMMA)
     >                       +ATAN(T0/GAMMA) )
            IF (SIGIX(K,I,J) .LT. 1.E-30) SIGIX(K,I,J) = 0.0
          ELSE
            SIGIX(K,I,J) = 0.0
          ENDIF
  140 CONTINUE
C
C
C Obtain high-energy correction factors:
C
      CALL HEXC(ENER,SIGIX,RATIO)
      DO 141 J=1,NBINS
        DO 141 I=1,NMAJ
          DO 141 K=1,NEI
            SIGIX(K,I,J)=SIGIX(K,I,J)/RATIO(J)
 141     CONTINUE
CADDED 11/20/11, JY
C   New cross-sections based on Johnson 2005 and Malone 2009
	    DO 142 J=1,NBINS
	    IF (ENER(J).GT.WWN(I) .AND. WWN(I).GT.0.001) THEN
              WEN=ENER(J)-WWN(I)
              SIGNE(I,J)=COA(I)*((WEN/ERYD)**(COB(I)))*
     >                 (1.0 + WEN/COC(I))**(-1.*(COB(I) + COD(I)))
     >                 +COE(I)*((WEN/ERYD)**(COF(I)))*
     >                 (1.0 + WEN/COG(I))**(-1.*(COF(I) + COH(I)))
          IF (SIGNE(I,J) .LT. 1.E-30) SIGNE(I,J) = 0.0
          ELSE
            SIGNE(I,J) = 0.0
          ENDIF
 142   CONTINUE

C  Redefine SIGEX, Zero the hi ener NANS in W and B' 11/20/11 JY
      DO 143 I=1,NBINS
      IF (ENER(I).GT.425.) SIGNE(4,I)=0.
      IF (ENER(I).GT.1150.) SIGNE(3,I)=0.

            SIGEX(1,3,I)=SIGNE(1,I)+SIGNE(2,I)+SIGNE(3,I)
	    SIGEX(2,3,I)=SIGNE(4,I)
            SIGEX(3,3,I)=SIGNE(5,I)+SIGNE(6,I)+SIGNE(7,I)
     >		+SIGNE(8,I)+SIGNE(9,I)
	    SIGEX(4,3,I)=SIGNE(10,I)+SIGNE(11,I)+SIGNE(12,I)
 143     CONTINUE

 


C Zero energy loss xsect and secondary production xsect arrays:
C
      DO 145 I1=1,NMAJ
      DO 145 I2=1,NBINS
      DO 145 I3=1,NBINS
        SIGA(I1,I2,I3)=0.0
        SEC(I1,I2,I3)=0.0
  145 CONTINUE
C
C
C Loop over energy:
C
      DO 500 JY=1,NBINS
C
      KUK=0
      KUK1=0
      ETJ=ENER(JY)
      DETJ=DEL(JY)
C
C
C Loop over species:
C
      DO 400 I=1,NMAJ
C
C
C Loop over excited states:
C
      DO 200 J=1,NNN(I)
C
C
C Calculate energy loss from JY to J-K for each species.
C The cross secton is divided proportionally between bin INV and bin
C INV-1, the two bins closest to J-K:
C
      SIGG = SIGEX(J,I,JY)
      ETA = ETJ - WW(J,I)
      IF (ETA .GT. 0.) THEN
        IE = INV (ETA,JY,ENER)
        IEE = IE - 1
        IF (IEE .LT. 1) IEE=IE
        K = JY - IE
        KK = JY - IEE
        IF (KK .GE. KUK) KUK = KK
        IF (IE .EQ. JY) THEN
          IF (JY .EQ. 1) THEN
            SIGA(I,1,JY) = SIGA(I,1,JY) + SIGG
          ELSE
            SIGA(I,1,JY) = SIGA(I,1,JY) + SIGG * (DETJ/DEL(JY-1))
     >                    * WW(J,I) / (ENER(JY)-ENER(JY-1))
          ENDIF
        ELSE
          IF (IE .EQ. 1) THEN
C           SIGA(I,K,JY) = SIGA(I,K,JY) + SIGG * 2.*ETA*DETJ/DEL(1)**2
            SIGA(I,K,JY) = SIGA(I,K,JY) + SIGG * DETJ/DEL(1)
          ELSE
            FF = (ENER(IE)-ETA) / (ENER(IE)-ENER(IEE))
            FF = 1.0 - ABS(FF)
            SIGA(I,K,JY) = SIGA(I,K,JY) + SIGG * FF * DETJ/DEL(IE)
            SIGA(I,KK,JY)= SIGA(I,KK,JY) + SIGG * (1.0-FF)*DETJ/DEL(IEE)
cc            if (i .eq. 3) then
c               if (jy .eq. 67) then
c                  write(*,*) 'e',ener(jy),k,kk,ie,iee,sigg,ff,
c     >                       detj/del(ie),detj/del(iee)
c               endif
c            endif
          ENDIF
          WAG=WW(J,I)-THI(1,I)
          IF (WAG .GT. 0. .AND. AUTO(J,I) .GT. 0.) THEN
            IBZ = INV (WAG,JY,ENER)
            SEC(I,IBZ,JY) = SEC(I,IBZ,JY)+SIGG*(DETJ/DEL(IBZ))*AUTO(J,I)
            IF(IBZ.GE.KUK1)KUK1=IBZ
          ENDIF
        ENDIF
      ENDIF
  200 CONTINUE
C
C
C Loop over ion states:
C
      DO 300 ML=1,NINN(I)
C
      DO 210 II=1,NBINS
      SIGI(II) = 0.0
      T12(II) = 0.0
  210 CONTINUE
C
C
C Calculate cross-section for production of secondaries into each
C bin from 1 to ITMAX and store in SIGI(II).  Apply relativistic correction.
C Also store the average energy of the secondaries in T12(II):
C
      WAG = THI(ML,I)
      TMAX = (ETJ-WAG) / 2.
      if (tmax .gt. 1.e6) tmax=1.e6
      IF (TMAX .LE. 0.) GOTO 300
      ITMAX = INV (TMAX,JY,ENER)
      IF (ITMAX .GE. KUK1) KUK1 = ITMAX + 1
      TMT = ENER(1) + DEL(1) / 2.0
      IF (TMAX .LT. TMT) TMT = TMAX
      SIGI(1) = SIGION(I,ML,ETJ,0.0,TMT,T12(1)) / RATIO(JY)
      TMT = ENER(1) + DEL(1) / 2.
      IF (TMAX .GT. TMT) THEN
        IF (TMAX .LE. ENER(2)) ITMAX = 2
        DO 220 II=2,ITMAX
        E1 = ENER(II) - DEL(II) / 2.
        E2 = E1 + DEL(II)
        IF (E2 .GT. TMAX) E2 = TMAX
        IF(E1 .LE. E2)SIGI(II)=SIGION(I,ML,ETJ,E1,E2,T12(II))/RATIO(JY)
  220   CONTINUE
      ENDIF
C
C
C Add the secondary production cross-section to SEC; calculate
C the ionization energy loss cross-section and add to SIGA:
C
      DO 250 II=1,ITMAX
      SEC(I,II,JY) = SEC(I,II,JY) + SIGI(II) * DETJ / DEL(II)
      WTH1 = T12(II) + WAG
      ETA = ETJ - WTH1
      IF (ETA .GT. 0.) THEN
        IE = INV (ETA,JY,ENER)
        IEE = IE - 1
        K = JY - IE
        KK = JY - IEE
        IF (IEE .LT. 1) IEE = IE
        IF (IE .EQ. JY) THEN
          SIGA(I,1,JY) = SIGA(I,1,JY) + SIGI(II) * (DETJ/DEL(JY-1))
     >                  * WTH1 / (ENER(JY)-ENER(JY-1))
        ELSE
          IF (KK .GE. KUK) KUK = KK
          IF (IE .EQ. 1) THEN
C           SIGA(I,K,JY) = SIGA(I,K,JY)+2.*ETA*SIGI(II)*DETJ/DEL(1)**2
            SIGA(I,K,JY) = SIGA(I,K,JY)+SIGI(II)*DETJ/DEL(1)
c            if (i .eq. 2) then
c               if (jy .eq. 67) then
c                  write(*,*) ener(jy),k,ie,ii,sigi(ii),detj/del(1)
c               endif
c            endif
          ELSE
            FF = (ENER(IE)-ETA) / (ENER(IE)-ENER(IEE))
            FF = 1.0 - ABS(FF)
            SIGA(I,K,JY) = SIGA(I,K,JY) + FF*SIGI(II)*DETJ/DEL(IE)
            SIGA(I,KK,JY)=SIGA(I,KK,JY)+(1.0-FF)*SIGI(II)*DETJ/DEL(IEE)
c            if (i .eq. 3) then
c               if (jy .eq. 67) then
c                  write(*,*) ener(jy),k,kk,ie,iee,ii,sigi(ii),ff
c               endif
c            endif
          ENDIF
        ENDIF
      ENDIF
C
  250 CONTINUE
C
  300 CONTINUE
C
  400 CONTINUE
C
      IIMAXX(JY) = KUK1
C
  500 CONTINUE
C
C
      goto 999
      write(*,*) 'Writing the e-impact cross sections...'


      open(unit=6,file='sigexO.dat')
      write(6,455)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 455  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(6,692)
 692  FORMAT ('     ENER     O(1D)     O(1S)     O(3s5S)   O(3s3S)'
     >        '   O(3p5P)   O(3p3P)   O(3d3D)   O(3sp3D)')
      
      do 920 i=1,1
c        do 920 k=1,NEI
          do 920 j=1,NBINS
          write (6,544) ENER(j),(SIGEX(k,i,j),k=1,nei)
 544      format (1x, f11.5, 10e10.2)
 920  continue
      close(6)

      write(*,*) '6 is closed'

      open(unit=7,file='sigexO2.dat')
      write(7,456)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 456  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(7,693)
 693  FORMAT ('     ENER     O2(a)     O2(b)     O2(AApc)  O2(B)'
     >       '     O2(9.9)   O2(Ryds)  O2(vib)   ')
      do 921 i=2,2
c        do 921 k=1,NEI
          do 921 j=1,NBINS
          write (7,545) ENER(j),(SIGEX(k,i,j),k=1,nei)
 545      format (1x, f11.5, 10e10.2)
 921  continue
      close(7)

      write(*,*) '7 is closed'

      open(unit=8,file='sigexN2.dat')
      write(8,457)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 457  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(8,694)
 694  FORMAT ('     ENER     N2(ABW)   N2(Bp)    N2(C)     N2(aaw)   ',
     >        'N2(1Pu)   N2(bp)    N2(Ryds)  N2(vib) ')     
      do 922 i=3,3
c        do 922 k=1,NEI
          do 922 j=1,NBINS
          write (8,546) ENER(j),(SIGEX(k,i,j),k=1,nei)
 546      format (1x, f11.5, 10e10.2)
 922   continue
      close(8)

      write(*,*) '8 is closed'

      open(unit=9,file='sigixO+.dat')
      write(9,555)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 555  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(9,792)
 792  FORMAT ('     ENER     O+(4S)    O+(2Do)   O+(2Po)   ')
      do 820 i=1,1
c        do 820 k=1,NEI
          do 820 j=1,NBINS
          write (9,644) ENER(j),(SIGIX(k,i,j),k=1,nei)
 644      format (1x, f11.5, 10e10.2)
 820  continue
      close(9)

      write(*,*) '9 is closed'
 
      open(unit=10,file='sigixO2+.dat')
      write(10,556)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 556  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(10,793)
 793  FORMAT ('     ENER     O2+(X)    O2+(a)    O2+(A)    O2+(b)'
     >       '    O2+(B)    O2+(c)    O2+(37eV)')   
      do 821 i=2,2
c        do 821 k=1,NEI
          do 821 j=1,NBINS
          write (10,645) ENER(j),(SIGIX(k,i,j),k=1,nei)
 645      format (1x, f11.5, 10e10.2)
 821  continue
      close(10)

      write(*,*) '10 is closed'

      open(unit=21,file='sigixN2+.dat')
      write(21,557)  IDATE, UT, GLAT, GLON, F107, F107A, AP
 557  FORMAT ('     Date=' ,i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
      WRITE(21,794)
 794  FORMAT ('     ENER     N2+(X)    N2+(A)    N2+(B)'
     >       '    N2+(D)    N2+(C)    N2+(40eV)')
      do 822 i=3,3
c       do 822 k=1,NEI
          do 822 j=1,NBINS
          write (21,646) ENER(j),(SIGIX(k,i,j),k=1,nei)
 646      format (1x, f11.5, 10e10.2)
  822  continue

      close(unit=21)
      
      Write(*,*) 'Writing siga...'

      open(unit=12,file='siga.dat')
      write(12,*) nbins
      write (12,647) (ENER(j),j=1,nbins)
      write (12,647) (del(j),j=1,nbins)
       do 823 i=1,nmaj
        do 823 k=1,Nbins
           write (12,647) (SIGa(i,k,j),j=1,nbins)
 823    continue
       close(unit=12)

      open(unit=13,file='sec.dat')
      write(13,*) nbins
      write (13,647) (ENER(j),j=1,nbins)
      write (13,647) (del(j),j=1,nbins)
       do 824 i=1,nmaj
        do 824 k=1,Nbins
          write (13,647) (sec(i,k,j),j=1,nbins)
  824      continue
       close(13)


       write(*,*) 'Wrote all the files...'

 647      format (1x, 10e10.2)

 999  continue
      RETURN
      END
C
C
C
C
C Function SIGION calculates ionization cross section for species I,
C state ML, primary energy E, secondary energy from E1 to E2 
C
      FUNCTION SIGION(I,ML,E,E1,E2,T12)
C
      PARAMETER (NMAJ=3)
      PARAMETER (NEI=10)
C
      COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ),
     >                ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ),
     >                THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ),
     >                TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ),
     >                GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)
C
      DOUBLE PRECISION ABB, ABC, ABD
      DATA QQ/1.E-16/
C
C
      IF (E .LE. THI(ML,I)) GOTO 30
C
      AK1=AK(ML,I)
      AJ1=AJ(ML,I)
      TS1=TS(ML,I)
      TA1=TA(ML,I)
      TB1=TB(ML,I)
      GAMS1=GAMS(ML,I)
      GAMB1=GAMB(ML,I)
      S=QQ*AK1*ALOG(E/AJ1)
      A=S/E
      TZ=TS1-TA1/(E+TB1)
      GG=(GAMS1*E)/(E+GAMB1)
      TTL=(E-THI(ML,I))/2.0
      TTL1=TTL-0.01
      IF(E1.GE.TTL1)GO TO 30
      IF(E2.GE.TTL)E2=TTL
      ABB=(E2-TZ)/GG
      ABC=(E1-TZ)/GG
      AL2=GG*GG*(ABB*ABB+1.0)
      AL1=GG*GG*(ABC*ABC+1.0)
      ABD=DATAN(ABB)-DATAN(ABC)
      T12=TZ+0.5*GG*(ALOG(AL2)-ALOG(AL1))/ABD
      SIGION=A*GG*ABD
      RETURN
C
   30 SIGION=0.0
      RETURN
C
      END
C
C
C
C
C Function INV finds the bin number closest to energy ETA on grid ENER.
C Bin INV or INV-1 will contain ETA.
C
      FUNCTION INV (ETA, JY, ENER)
C
      INCLUDE 'glow.h'
C
      DIMENSION ENER(NBINS)
C
      IF (ETA .LT. 0.) THEN
        INV = -1
      ELSE
        DO 30 IV=1,JY
          IF (ETA .LE. ENER(IV)) GOTO 40
   30   CONTINUE
        IV = JY
   40   INV = IV
      ENDIF
C
      RETURN
      END
C
C
C
C
C Subroutine HEXC
C
C High Energy Cross Section Correction
C Calculates ratio of low energy (non-relativistic) to high energy
C (relativistic) ionization cross sections, based on N2.
C Extends to 1 GeV.
C
C Originally coded by Ann Windnagel, 11/98
C Re-written by Stan Solomon, 2/99
C Re-designed with table lookup, SCS, 4/99
C Updated comments, SCS, 4/02
C References:
C   Porter et al., J. Chem. Phys., 65, 154, 1976.
C   Rieke and Prepejchal, Phys. Rev. A, 6, 1507, 1990.
C   Saksena et al., Int. Jour. of Mass Spec. & Ion Proc., 171, L1, 1997.


      SUBROUTINE HEXC(ENER,SIGIX,RATIO)

      INCLUDE 'glow.h'
      PARAMETER (NMAJ=3)
      PARAMETER (NEI=10)

      DIMENSION ENER(NBINS), SIGIX(NEI,NMAJ,NBINS), RATIO(NBINS)
      DIMENSION TOTX(NBINS), TOTNEW(NBINS), EGR(13), SGR(13)
      DATA EGR/1.E4,      2.E4,      5.E4,      1.E5,      2.E5,
     >         3.E5,      5.E5,      1.E6,      2.E6,      5.E6,
     >         1.E7,      1.E8,      1.E9/
      DATA SGR/1.20E-17,  7.03E-18,  3.37E-18,  1.96E-18,  1.26E-18,
     >         1.05E-18,  9.50E-19,  9.00E-19,  9.00E-19,  9.40E-19,
     >         1.00E-18,  1.26E-18,  1.59E-18/


C Calculate total low-energy cross section for N2:

      DO 20 K = 1,NBINS
        TOTX(K) = 0.
        DO 20 I = 1,NEI
          TOTX(K) = TOTX(K) + SIGIX(I,3,K)
   20 CONTINUE


C Calculate high-energy cross section for N2, using tabulated values:

      DO 70 K=1,NBINS
        IF (ENER(K) .GE. EGR(1)) THEN
        DO 60 KG=1,12
          IF (ENER(K) .GE. EGR(KG) .AND. ENER(K) .LT. EGR(KG+1))
     >    TOTNEW(K)=TERPOO(ENER(K),EGR(KG),EGR(KG+1),SGR(KG),SGR(KG+1))
   60   CONTINUE
        ELSE
          TOTNEW(K)=TOTX(K)
        ENDIF
   70 CONTINUE


C Calculate ratio (=1 < 10 keV):

      DO 90 K = 1,NBINS
        IF (ENER(K) .GE. EGR(1)) THEN
          RATIO(K) = TOTX(K)/TOTNEW(K) 
C         IF (RATIO(K) .GT. 1.) RATIO(K) = 1.
        ELSE
          RATIO(K) = 1.
        ENDIF
   90 CONTINUE

      RETURN

      END



      FUNCTION TERPOO(X,X1,X2,Y1,Y2)
      TERPOO = EXP ( ALOG(Y1) + ALOG(X/X1)*ALOG(Y2/Y1)/ALOG(X2/X1) )
      RETURN
      END

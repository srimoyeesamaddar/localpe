
;C Calculates electron impact cross sections
;C
;C Definitions:
;C SIGS   elastic cross sections for each species, energy; cm2
;C PE     elastic backscatter probabilities for each species, energy
;C PI     inelastic  "
;C SIGA   energy loss cross section for each species, loss, energy; cm2
;C SEC    secondary production xsect for species, Esec, Epri; cm2
;C SIGEX  excitation xsect for each state, species, energy; cm2
;C SIGIX  ionization xsect for each state, species, energy; cm2
;C IIMAX  number of bins for secondary production for each primary energy
;C WW     energy threshold for each excited state, species; eV
;C WW, AO, OMEG, ANU, BB: revised excitation cross section parameters,
;C        from Green & Stolarski (1972) formula (W, A, omega, nu, gamma)
;C AUTO   autoionization coefs (= 0 as autoion. included in ion xsects)
;C THI    energy threshold for each ionized state, species; eV
;C AK, AJ, TS, TA, TB, GAMS, GAMB:  Jackman et al (1977) ioniz. params
;C ENER   energy grid; eV
;C DEL    energy grid spacing; eV
;C NNN    number of excited states for each species
;C NINN   number of ionized states for each species
;C NUM    number of points on elastic data trid for each species
;C EC     data energy grid of elastic xsects and backscatter ratios
;C        for each species; eV
;C CC     elastic xsects on data grid for each species, cm2
;C CE     elastic backscat. probs on data grid for each species; cm2
;C CI     inelastic "
;C
;C Array dimensions:
;C NBINS  number of energy levels
;C NMAJ   number of major species
;C NEI    number of slots for excited and ionized states

    
 ; NMAJ   number of major species
 ; NEI    number of slots for excited and ionized states
 ;ww=dblarr(nei,nmaj)
    nbins=190
    NMAJ=3
    NEI=10
    JMAX=120
    LMAX=123

    SIGS=dblarr(NMAJ,NBINS)
    PE=dblarr(NMAJ,NBINS)
    PI=dblarr(NMAJ,NBINS)
    SIGA=dblarr(NMAJ,NBINS,NBINS)
    SEC=dblarr(NMAJ,NBINS,NBINS)
    SIGEX=dblarr(NEI,NMAJ,NBINS)
    SIGIX=dblarr(NEI,NMAJ,NBINS)
    IIMAXX=dblarr(NBINS) 

    WW=dblarr(NEI,NMAJ)
    AO=dblarr(NEI,NMAJ)
    OMEG=dblarr(NEI,NMAJ)
    ANU=dblarr(NEI,NMAJ)
    BB=dblarr(NEI,NMAJ)
    AUTO=dblarr(NEI,NMAJ)
    THI=dblarr(NEI,NMAJ)
    AK=dblarr(NEI,NMAJ)
    AJ=dblarr(NEI,NMAJ)
    TS=dblarr(NEI,NMAJ)
    TA=dblarr(NEI,NMAJ)
    TB=dblarr(NEI,NMAJ)
    GAMS=dblarr(NEI,NMAJ)
    GAMB=dblarr(NEI,NMAJ)

    ENER=dblarr(NBINS)
    DEL=dblarr(NBINS)
    SIGI=dblarr(NBINS)
    T12=dblarr(NBINS)
    ;RATIO=dblarr(NBINS)
    NNN=dblarr(NMAJ)
    NINN=dblarr(NMAJ)
    NUM=dblarr(NMAJ)
    EC=dblarr(31,NMAJ)
    CC=dblarr(31,NMAJ)
    CE=dblarr(31,NMAJ)
    CI=dblarr(31,NMAJ)


    QQN=6.51E-14
    NNN=[8,7,8]
    NINN=[3,7,6]
    NUM=[31,28,28]
    ratio=dblarr(nbins)
    RATIO(*)=1.0

;C     O states:   1D,   1S, 3s5S, 3s3S, 3p5P, 3p3P, 3d3D, 3s'3D
;C     O2 states:   a,    b, AA'c,    B,  9.9, Ryds,  vib
;C     N2 states: A,BW,   B',    C, aa'w,  1Pu,   b', Ryds,  vib

    WW=[[1.96, 4.17, 9.29, 9.53,10.76,10.97,12.07,12.54, 0.,0.],[0.98, 1.64, 4.50, 8.44, 9.90,13.50, 0.25, 0.00, 0.,0.],$
    [6.17, 8.16,11.03, 8.40,12.85,14.00,13.75, 1.85, 0.,0.]]

    AO=[[.0100,.0042,.1793,.3565,.0817,.0245,.0293,.1221, 0.,0.],$
              [.0797,.0211,.0215,.3400,.0657,1.110,3.480, 0.00, 0.,0.],$
              [2.770,.1140,.1790,.0999,.8760,.6010,1.890,1.350, 0.,0.]]
   
    OMEG=[[1.00, 1.00, 3.00, 0.75, 3.00, 0.85, 0.75, 0.75, 0.,0.],$
               [2.00, 2.00, 1.15, 0.75, 0.75, 0.75, 7.00, 0.00, 0.,0.],$
               [3.00, 3.00, 3.00, 1.00, 0.75, 0.75, 0.75, 8.00, 0.,0.]]

    ANU=[[2.00, 1.04, 2.53, 0.54, 2.43, 2.87, 0.93, 0.72, 0.,0.],$
               [6.18, 4.14, 1.00, 1.05, 1.60, 3.00,10.87, 0.00, 0.,0.],$
               [4.53, 4.78, 4.32, 4.05, 1.47, 1.27, 3.00, 1.58, 0.,0.]]
    BB=[[1.00, 0.50, 1.02, 0.01, 4.19, 4.88, 0.66, 0.17, 0.,0.],$
               [0.53, 0.51, 0.98, 0.99, 1.86, 1.00, 1.00, 0.00, 0.,0.],$
               [1.42, 3.54,12.70, 5.20, 0.86, 0.45, 1.00, 1.00, 0.,0.]]
    AUTO= fltarr(nei,nmaj)
    THI=[[13.60,16.90,18.50, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
              [12.10,16.10,16.90,18.20,20.00,23.00,37.00, 0.,0.,0.],$
              [15.58,16.73,18.75,22.00,23.60,40.00, 0.00, 0.,0.,0.]]
     AK=[[ 1.13, 1.25, 0.67, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
               [0.47, 1.13, 1.13, 1.01, 0.65, 0.95, 0.59, 0.,0.,0.],$
               [2.42, 1.06, 0.55, 0.37, 0.37, 0.53, 0.00, 0.,0.,0.]]
     AJ=[[ 1.81, 1.79, 1.78, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
               [3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 0.,0.,0.],$
               [1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 0.00, 0.,0.,0.]]
    TS=[[ 6.41, 6.41, 6.41, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
              [ 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 0.,0.,0.],$
               [4.71, 4.71, 4.71, 4.71, 4.71, 4.71, 0.00, 0.,0.,0.]]
    TA=[[3450.,3450.,3450.,   0.,   0.,   0.,   0., 0.,0.,0.],$
              [1000.,1000.,1000.,1000.,1000.,1000.,1000., 0.,0.,0.],$
              [1000.,1000.,1000.,1000.,1000.,1000.,   0., 0.,0.,0.]]
     TB= [[162.00,162.0,162.0, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
              [24.20,32.20,33.80,36.40,40.60,46.00,74.00, 0.,0.,0.],$
             [31.16,33.46,37.50,44.00,47.20,80.00, 0.00, 0.,0.,0.]]
    GAMS= [[13.00,13.0,13.00, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
               [18.50,18.5,18.50,18.50,18.50,18.50,18.50, 0.,0.,0.],$
               [13.80,13.8,13.80,13.80,13.80,13.80, 0.00, 0.,0.,0.]]
    GAMB= [[-.815,-.815,-.815,0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
              [ 12.10,16.10,16.90,18.2,20.30,23.00,37.00, 0.,0.,0.],$
              [ 15.58,16.73,18.75,22.0,23.60,40.00, 0.00, 0.,0.,0.]]
     EC =    [[  1.00,     2.00,     4.00,     6.00,     8.00,$
                   10.00,    12.00,    14.00,    16.00,    18.00,$
                   20.00,    30.00,    40.00,    50.00,    60.00,$
                   70.00,    80.00,    90.00,   100.00,   150.00,$
                  200.00,   300.00,   500.00,  1000.00,  2000.00,$
                 3000.00,  5000.00, 10000.00, 20000.00, 40000.00,$
                50000.00],$
                  [  1.00,     2.00,     3.00,     5.00,     7.00,$
                   10.00,    15.00,    20.00,    30.00,    40.00,$
                   50.00,    70.00,   100.00,   150.00,   200.00,$
                  300.00,   400.00,   500.00,   600.00,   700.00,$
                 1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
                20000.00, 40000.00, 50000.00,     0.00,     0.00,$
                    0.00],$
                 [   1.00,     2.00,     2.50,     3.00,     4.00,$
                    5.00,     6.00,     8.00,    10.00,    15.00,$
                   20.00,    30.00,    40.00,    50.00,    70.00,$
                  100.00,   200.00,   300.00,   500.00,   700.00,$
                 1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
                20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0]]
    CC=   [[ 5.00E-16, 6.00E-16, 7.50E-16, 7.60E-16, 7.70E-16,$
                7.80E-16, 7.50E-16, 7.20E-16, 6.90E-16, 6.70E-16,$
                6.50E-16, 5.60E-16, 4.60E-16, 4.00E-16, 3.50E-16,$
                3.20E-16, 2.90E-16, 2.70E-16, 2.50E-16, 1.90E-16,$
                1.50E-16, 1.20E-16, 8.00E-17, 5.00E-17, 3.02E-17,$
                1.99E-17, 1.20E-17, 6.08E-18, 3.06E-18, 1.55E-18,$
                1.24E-18],$
                [5.50E-16, 6.90E-16, 7.50E-16, 8.50E-16, 9.60E-16,$
                1.00E-15, 1.00E-15, 9.00E-16, 8.30E-16, 7.70E-16,$
                6.90E-16, 5.70E-16, 4.40E-16, 3.30E-16, 2.70E-16,$
                2.10E-16, 1.80E-16, 1.60E-16, 1.40E-16, 1.30E-16,$
                1.10E-16, 7.00E-17, 5.00E-17, 3.00E-17, 1.53E-17,$
                7.72E-18, 3.90E-18, 3.13E-18, 0.00E+00, 0.00E+00,$
                0.00E+00],$
               [ 9.00E-16, 2.27E-15, 2.52E-15, 1.93E-15, 1.32E-15,$
                1.15E-15, 1.16E-15, 1.17E-15, 1.18E-15, 1.14E-15,$
                1.13E-15, 9.50E-16, 8.60E-16, 7.30E-16, 5.90E-16,$
                4.70E-16, 3.30E-16, 2.50E-16, 1.60E-16, 1.30E-16,$
                1.10E-16, 6.35E-17, 4.18E-17, 2.54E-17, 1.28E-17,$
                6.44E-18, 3.27E-18, 2.62E-18, 0.00E+00, 0.00E+00, 0.0]]
    CE =  [[ 0.50000,  0.49500,  0.46800,  0.43600,  0.42000,$
                 0.40500,  0.37000,  0.36000,  0.34000,  0.33000,$
                 0.32000,  0.27000,  0.24000,  0.22000,  0.20000,$
                 0.18000,  0.17000,  0.16000,  0.15000,  0.13000,$
                 0.11500,  0.09000,  0.06800,  0.04600,  0.02400,$
                 0.01660,  0.01000,  0.00510,  0.00255,  0.00125,$
                 0.00100],$
                 [0.50000,  0.50000,  0.49000,  0.44500,  0.42700,$
                 0.40500,  0.36800,  0.34300,  0.31600,  0.28900,$
                 0.25800,  0.22000,  0.18400,  0.16400,  0.13300,$
                 0.11000,  0.10000,  0.09200,  0.08500,  0.08000,$
                 0.06800,  0.03700,  0.02600,  0.01600,  0.00800,$
                 0.00400,  0.00200,  0.00160,  0.00000,  0.00000,$
                 0.00000],$
                [ 0.50000,  0.50000,  0.50000,  0.49000,  0.46800,$
                 0.44500,  0.43600,  0.42000,  0.40500,  0.36800,$
                 0.34300,  0.31600,  0.28900,  0.25800,  0.22000,$
                 0.18400,  0.14000,  0.11000,  0.08400,  0.07400,$
                 0.06300,  0.03400,  0.02400,  0.01500,  0.00740,$
                 0.00370,  0.00180,  0.00140,  0.00000,  0.00000, 0.0]]
     CI=   [ [0.60000,  0.60000,  0.60000,  0.60000,  0.60000,$
                 0.60000,  0.55000,  0.46000,  0.40000,  0.36000,$
                 0.32000,  0.22000,  0.15000,  0.10000,  0.08200,$
                 0.07000,  0.06100,  0.05400,  0.05000,  0.04400,$
                 0.03800,  0.02800,  0.02000,  0.01050,  0.00600,$
                 0.00400,  0.00250,  0.00130,  0.00060,  0.00030,$
                 0.00025],$
                 [0.50000,  0.50000,  0.50000,  0.50000,  0.48000,$
                 0.44000,  0.36000,  0.28000,  0.20000,  0.14000,$
                 0.10000,  0.07000,  0.05000,  0.04600,  0.04300,$
                 0.03700,  0.03200,  0.02800,  0.02400,  0.02100,$
                 0.01600,  0.00900,  0.00620,  0.00400,  0.00200,$
                 0.00100,  0.00050,  0.00040,  0.00000,  0.00000,$
                 0.00000],$
                 [0.50000,  0.50000,  0.50000,  0.50000,  0.50000,$
                 0.50000,  0.50000,  0.50000,  0.50000,  0.50000,$
                 0.44000,  0.30000,  0.20000,  0.13000,  0.09000,$
                 0.06000,  0.05000,  0.04200,  0.03200,  0.02500,$
                 0.02000,  0.01100,  0.00800,  0.00500,  0.00250,$
                 0.00120,  0.00060,  0.00050,  0.00000,  0.00000, 0.0]]

  
; set up energy array

     FOR N=0,NBINS-1 DO BEGIN
        IF (N LE 21) THEN BEGIN
          ENER(N) = 0.5 * FLOAT(N)
        ENDIF ELSE BEGIN
          ENER(N) = EXP (0.05 * FLOAT(N+26))
        ENDELSE
     ENDFOR
      DEL(0) = 0.5
      FOR N=1,NBINS-1 DO BEGIN
        DEL(N) = ENER(N)-ENER(N-1)
    ENDFOR
      FOR N=0,NBINS-1 DO BEGIN
        ENER(N) = ENER(N) - DEL(N)/2.0
     ENDFOR


 ;C Interpolate elastic cross sections and backscatter ratios:
      for IJ=0,NMAJ-1 do begin
        for IV=0,NBINS-1 do begin
          EX=ENER(IV)
          for II=0,NUM(IJ)-1 do begin
            IF (EC(II,IJ) GT EX) then GOTO,err 
          endfor
          SIGS(IJ,IV)=CC(NUM(IJ),IJ)*((EC(NUM(IJ),IJ)/EX)^0.8)
          IF(IJ EQ 0) then SIGS(IJ,IV)=CC(NUM(IJ),IJ)*((EC(NUM(IJ),IJ)/EX)^2)
          PE(IJ,IV) = CE(NUM(IJ),IJ)* (EC(NUM(IJ),IJ)/EX)
          PI(IJ,IV) = CI(NUM(IJ),IJ)* (EC(NUM(IJ),IJ)/EX)
          GOTO,EIGHTY
   err:     I=II-1
          IF (I LE 0) THEN begin
            SIGS(IJ,IV)=CC(II,IJ)
            PE(IJ,IV)=CE(II,IJ)
            PI(IJ,IV)=CI(II,IJ)
            endif ELSE begin
            FAC = ALOG(EX/EC(I,IJ)) / ALOG(EC(II,IJ)/EC(I,IJ))
            SIGS(IJ,IV) = EXP(ALOG (CC(I,IJ))+ ALOG (CC(II,IJ)/CC(I,IJ)) * FAC)
            PE(IJ,IV) = EXP(ALOG (CE(I,IJ))+ ALOG (CE(II,IJ)/CE(I,IJ)) * FAC)
            PI(IJ,IV) = EXP(ALOG (CI(I,IJ))+ ALOG (CI(II,IJ)/CI(I,IJ)) * FAC)
         Endelse
   EIGHTY: 
         endfor
         endfor



; Calculate electron impact excitation and ionization cross sections:

      for I=0,NMAJ-1 do begin
        for K=0,NEI-1 do begin
          for J=0,NBINS-1 do begin
            IF ((ENER(J) GT WW(K,I)) AND (WW(K,I) GT 0.001)) THEN begin
              WE = WW(K,I) / ENER(J)
              SIGEX(K,I,J) = QQN * AO(K,I)*(WE^(OMEG(K,I)) / (WW(K,I)^2)) * [(1.0 - (WE^BB(K,I)))^ ANU(K,I)]
              IF (SIGEX(K,I,J) LT 1E-30) THEN SIGEX(K,I,J) = 0.0 
            endif ELSE begin
              SIGEX(K,I,J) = 0.0
            Endelse
            IF [(ENER(J) GT THI(K,I)) AND (THI(K,I) GT 0.001)] THEN begin
              AE = AK(K,I)/ENER(J) * ALOG(ENER(J)/AJ(K,I))
              GAMMA = GAMS(K,I) * ENER(J) / (ENER(J)+GAMB(K,I))
              T0 = TS(K,I) - (TA(K,I)/(ENER(J)+TB(K,I)))
              SIGIX(K,I,J) = 1E-16 * AE * GAMMA* ( ATAN(((ENER(J)-THI(K,I))/2.-T0)/GAMMA)+ATAN(T0/GAMMA) )
              IF (SIGIX(K,I,J) LT 1E-30) then SIGIX(K,I,J) = 0.0
            ENDIF ELSE BEGIN
              SIGIX(K,I,J) = 0.0
            ENDELSE
         ENDFOR
       endfor
     ENDFOR

;***************8n2bw reduced by .5      
;;     N2 states: ABW,   B',    C, aa'w,  1Pu,   b', Ryds,  vib
;C,W,B',B,A.

restore,'john_mal_xs.sav'
;restore,'broad_xs.sav'
;sigex(0,2,*)=xsect(4,*)+xsect(1,*)+xsect(3,*)
;sigex(1,2,*)=xsect(2,*)
;sigex(2,2,*)=xsect_c(0,*)+xsect_c(1,*)+xsect_c(2,*)+xsect_c(3,*)+xsect_c(4,*);C

;sigex(3,2,*)=sigex(3,2,*)/2.

;sigex(6,2,*)=sigex(6,2,*)/2.





sec=fltarr(nmaj,nbins,nbins) ; (species,lower energy, upper energy) cross section for cascade forming PE with lower energy by 
; PE with higher energy and collisions (ionization) with species, cm^2
siga=sec ; (species, lower energy, higher energy) cross section for energy loss by PE of higher energy to PE of lower energy
                                ; through collisions (ionizatin & excitation) with species, cm^2

for inmaj=0,nmaj-1 do begin
for iprim=nbins-1,0,-1 do begin

for iionstate=0,ninn(inmaj)-1 do begin

newenergylo= ener(iprim)-del(iprim)/2.0 - thi(iionstate,inmaj)
newenergyhi= ener(iprim)+del(iprim)/2.0 - thi(iionstate,inmaj)

wbin,ener,del,newenergylo,newenergyhi,w
sec(inmaj,*,iprim) =sec(inmaj,*,iprim) + sigix(iionstate,inmaj,iprim)*w
siga(inmaj,*,iprim)=siga(inmaj,*,iprim)+ sigix(iionstate,inmaj,iprim)*w

endfor

for iexstate=0,nnn(inmaj)-1 do begin

newenergylo= ener(iprim)-del(iprim)/2.0 - ww(iexstate,inmaj)
newenergyhi= ener(iprim)+del(iprim)/2.0 - ww(iexstate,inmaj)

wbin,ener,del,newenergylo,newenergyhi,w
siga(inmaj,*,iprim)=siga(inmaj,*,iprim) + sigex(iexstate,inmaj,iprim)*w

endfor
endfor
endfor

  ;    RETURN

file='exsect_cs.sav'
save,file=file,/var
;stop
END


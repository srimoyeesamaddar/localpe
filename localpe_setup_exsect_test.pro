;pro localpe_setup_exsect,demo=demo
;
;  ; Calculates electron impact cross sections
;  ; cross section fit parameters based on work of Green and Tawada 1972 and Jackman et al 1977
;  ; fits are redone by Solomon and taken from GLOW package
;
;  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
;  ; see this file for definitions
;
;
;  ; define array dimensions
;  ; nbins=190  ; number of electron energy bins
;  ;nmaj=3                 ; number of major species in atmosphere (o, o2, and n2)
;  nei=10                  ; number of product states for excitation, ionization
;
;  ;nei=14  ;Making 14 states to add N2**
;
;  ; cross section arrays
;  sigex=fltarr(nei,nmaj,nbins)
;  sigix=fltarr(nei,nmaj,nbins)
;  iimaxx=fltarr(nbins)
;
;  ener=fltarr(nbins)
;  del=fltarr(nbins)
;  sigi=fltarr(nbins)
;  t12=fltarr(nbins)
;
;  qqn=6.51e-14
;  nnn=[8,7,8]
;  ;nnn=[8,7,14]  ;Making 14 states for N2 **
;  ninn=[3,7,6]
;  num=[31,28,28]
;  ratio=fltarr(nbins)+1.0
;
;
;
;  ;C     O states:   1D,   1S, 3s5S, 3s3S, 3p5P, 3p3P, 3d3D, 3s'3D
;  ;C     O2 states:   a,    b, AA'c,    B,  9.9, Ryds,  vib
;  ;C     N2 states: A,BW,   B',    C, aa'w,  1Pu,   b', Ryds,  vib
;  eistates=$
;    [['1D',   '1S',  '3s5S',  '3s3S', '3p5P', '3p3P', '3d3D', '3s3D',' ',' '],$
;    ['a',     'b',    'AA',     'c',      'B',         '9.9',    'Ryds',   'vib',    ' ',' '],$
;    ['A',    'BW', "B'",     'C',     "aa'w",   '1Pu',   "b'",       'Ryds',  'vib',' ']]
;
;  ;NOTE:Srimoyee- Adding 4 more extra zeros in each parameters to account for 4 extra states of N2**
;  ;;
;  ;    ww=[[1.96, 4.17, 9.29, 9.53,10.76,10.97,12.07,12.54, 0.,0.],[0.98, 1.64, 4.50, 8.44, 9.90,13.50, 0.25, 0.00, 0.,0.],$
;  ;       [6.17, 8.16,11.03, 8.40,12.85,14.00,13.75, 1.85, 0.,0.]] ;N2 'ww' for original model
;
;
;  ww=[[1.96, 4.17, 9.29, 9.53,10.76,10.97,12.07,12.54, 0.,0., 0.,0., 0.,0.],[0.98, 1.64, 4.50, 8.44, 9.90,13.50, 0.25, 0.00, 0.,0., 0.,0., 0.,0.],$
;    [6.169,7.353,7.362,8.165,8.399,8.549,11.032,11.875,12.255,12.500,12.935,12.854,8.890,12.910]]
;
;
;
;
;  ao=[[.0100,.0042,.1793,.3565,.0817,.0245,.0293,.1221, 0.,0.],$
;    [.0797,.0211,.0215,.3400,.0657,1.110,3.480, 0.00, 0.,0.],$
;    [2.770,.1140,.1790,.0999,.8760,.6010,1.890,1.350, 0.,0.]]
;
;  omeg=[[1.00, 1.00, 3.00, 0.75, 3.00, 0.85, 0.75, 0.75, 0.,0.],$
;    [2.00, 2.00, 1.15, 0.75, 0.75, 0.75, 7.00, 0.00, 0.,0.],$
;    [3.00, 3.00, 3.00, 1.00, 0.75, 0.75, 0.75, 8.00, 0.,0.]]
;
;  anu=[[2.00, 1.04, 2.53, 0.54, 2.43, 2.87, 0.93, 0.72, 0.,0.],$
;    [6.18, 4.14, 1.00, 1.05, 1.60, 3.00,10.87, 0.00, 0.,0.],$
;    [4.53, 4.78, 4.32, 4.05, 1.47, 1.27, 3.00, 1.58, 0.,0.]]
;  bb=[[1.00, 0.50, 1.02, 0.01, 4.19, 4.88, 0.66, 0.17, 0.,0.],$
;    [0.53, 0.51, 0.98, 0.99, 1.86, 1.00, 1.00, 0.00, 0.,0.],$
;    [1.42, 3.54,12.70, 5.20, 0.86, 0.45, 1.00, 1.00, 0.,0.]]
;
;  auto= fltarr(nei,nmaj)
;  thi=[[13.60,16.90,18.50, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [12.10,16.10,16.90,18.20,20.00,23.00,37.00, 0.,0.,0.],$
;    [15.58,16.73,18.75,22.00,23.60,40.00, 0.00, 0.,0.,0.]]
;
;  ak=[[ 1.13, 1.25, 0.67, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [0.47, 1.13, 1.13, 1.01, 0.65, 0.95, 0.59, 0.,0.,0.],$
;    [2.42, 1.06, 0.55, 0.37, 0.37, 0.53, 0.00, 0.,0.,0.]]
;
;  aj=[[ 1.81, 1.79, 1.78, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 0.,0.,0.],$
;    [1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 0.00, 0.,0.,0.]]
;
;  ts=[[ 6.41, 6.41, 6.41, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [ 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 0.,0.,0.],$
;    [4.71, 4.71, 4.71, 4.71, 4.71, 4.71, 0.00, 0.,0.,0.]]
;
;
;  ta=[[3450.,3450.,3450.,   0.,   0.,   0.,   0., 0.,0.,0.],$
;    [1000.,1000.,1000.,1000.,1000.,1000.,1000., 0.,0.,0.],$
;    [1000.,1000.,1000.,1000.,1000.,1000.,   0., 0.,0.,0.]]
;
;  tb= [[162.00,162.0,162.0, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [24.20,32.20,33.80,36.40,40.60,46.00,74.00, 0.,0.,0.],$
;    [31.16,33.46,37.50,44.00,47.20,80.00, 0.00, 0.,0.,0.]]
;
;  gams= [[13.00,13.0,13.00, 0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [18.50,18.5,18.50,18.50,18.50,18.50,18.50, 0.,0.,0.],$
;    [13.80,13.8,13.80,13.80,13.80,13.80, 0.00, 0.,0.,0.]]
;
;  gamb= [[-.815,-.815,-.815,0.00, 0.00, 0.00, 0.00, 0.,0.,0.],$
;    [ 12.10,16.10,16.90,18.2,20.30,23.00,37.00, 0.,0.,0.],$
;    [ 15.58,16.73,18.75,22.0,23.60,40.00, 0.00, 0.,0.,0.]]
;
;  ;Not changing these cuz don't know what they are **
;  ec =    [[  1.00,     2.00,     4.00,     6.00,     8.00,$
;    10.00,    12.00,    14.00,    16.00,    18.00,$
;    20.00,    30.00,    40.00,    50.00,    60.00,$
;    70.00,    80.00,    90.00,   100.00,   150.00,$
;    200.00,   300.00,   500.00,  1000.00,  2000.00,$
;    3000.00,  5000.00, 10000.00, 20000.00, 40000.00,$
;    50000.00],$
;    [  1.00,     2.00,     3.00,     5.00,     7.00,$
;    10.00,    15.00,    20.00,    30.00,    40.00,$
;    50.00,    70.00,   100.00,   150.00,   200.00,$
;    300.00,   400.00,   500.00,   600.00,   700.00,$
;    1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
;    20000.00, 40000.00, 50000.00,     0.00,     0.00,$
;    0.00],$
;    [   1.00,     2.00,     2.50,     3.00,     4.00,$
;    5.00,     6.00,     8.00,    10.00,    15.00,$
;    20.00,    30.00,    40.00,    50.00,    70.00,$
;    100.00,   200.00,   300.00,   500.00,   700.00,$
;    1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
;    20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0]]
;  cc=   [[ 5.00e-16, 6.00e-16, 7.50e-16, 7.60e-16, 7.70e-16,$
;    7.80e-16, 7.50e-16, 7.20e-16, 6.90e-16, 6.70e-16,$
;    6.50e-16, 5.60e-16, 4.60e-16, 4.00e-16, 3.50e-16,$
;    3.20e-16, 2.90e-16, 2.70e-16, 2.50e-16, 1.90e-16,$
;    1.50e-16, 1.20e-16, 8.00e-17, 5.00e-17, 3.02e-17,$
;    1.99e-17, 1.20e-17, 6.08e-18, 3.06e-18, 1.55e-18,$
;    1.24e-18],$
;    [5.50e-16, 6.90e-16, 7.50e-16, 8.50e-16, 9.60e-16,$
;    1.00e-15, 1.00e-15, 9.00e-16, 8.30e-16, 7.70e-16,$
;    6.90e-16, 5.70e-16, 4.40e-16, 3.30e-16, 2.70e-16,$
;    2.10e-16, 1.80e-16, 1.60e-16, 1.40e-16, 1.30e-16,$
;    1.10e-16, 7.00e-17, 5.00e-17, 3.00e-17, 1.53e-17,$
;    7.72e-18, 3.90e-18, 3.13e-18, 0.00e+00, 0.00e+00,$
;    0.00e+00],$
;    [ 9.00e-16, 2.27e-15, 2.52e-15, 1.93e-15, 1.32e-15,$
;    1.15e-15, 1.16e-15, 1.17e-15, 1.18e-15, 1.14e-15,$
;    1.13e-15, 9.50e-16, 8.60e-16, 7.30e-16, 5.90e-16,$
;    4.70e-16, 3.30e-16, 2.50e-16, 1.60e-16, 1.30e-16,$
;    1.10e-16, 6.35e-17, 4.18e-17, 2.54e-17, 1.28e-17,$
;    6.44e-18, 3.27e-18, 2.62e-18, 0.00e+00, 0.00e+00, 0.0]]
;  ce =  [[ 0.50000,  0.49500,  0.46800,  0.43600,  0.42000,$
;    0.40500,  0.37000,  0.36000,  0.34000,  0.33000,$
;    0.32000,  0.27000,  0.24000,  0.22000,  0.20000,$
;    0.18000,  0.17000,  0.16000,  0.15000,  0.13000,$
;    0.11500,  0.09000,  0.06800,  0.04600,  0.02400,$
;    0.01660,  0.01000,  0.00510,  0.00255,  0.00125,$
;    0.00100],$
;    [0.50000,  0.50000,  0.49000,  0.44500,  0.42700,$
;    0.40500,  0.36800,  0.34300,  0.31600,  0.28900,$
;    0.25800,  0.22000,  0.18400,  0.16400,  0.13300,$
;    0.11000,  0.10000,  0.09200,  0.08500,  0.08000,$
;    0.06800,  0.03700,  0.02600,  0.01600,  0.00800,$
;    0.00400,  0.00200,  0.00160,  0.00000,  0.00000,$
;    0.00000],$
;    [ 0.50000,  0.50000,  0.50000,  0.49000,  0.46800,$
;    0.44500,  0.43600,  0.42000,  0.40500,  0.36800,$
;    0.34300,  0.31600,  0.28900,  0.25800,  0.22000,$
;    0.18400,  0.14000,  0.11000,  0.08400,  0.07400,$
;    0.06300,  0.03400,  0.02400,  0.01500,  0.00740,$
;    0.00370,  0.00180,  0.00140,  0.00000,  0.00000, 0.0]]
;  ci=   [ [0.60000,  0.60000,  0.60000,  0.60000,  0.60000,$
;    0.60000,  0.55000,  0.46000,  0.40000,  0.36000,$
;    0.32000,  0.22000,  0.15000,  0.10000,  0.08200,$
;    0.07000,  0.06100,  0.05400,  0.05000,  0.04400,$
;    0.03800,  0.02800,  0.02000,  0.01050,  0.00600,$
;    0.00400,  0.00250,  0.00130,  0.00060,  0.00030,$
;    0.00025],$
;    [0.50000,  0.50000,  0.50000,  0.50000,  0.48000,$
;    0.44000,  0.36000,  0.28000,  0.20000,  0.14000,$
;    0.10000,  0.07000,  0.05000,  0.04600,  0.04300,$
;    0.03700,  0.03200,  0.02800,  0.02400,  0.02100,$
;    0.01600,  0.00900,  0.00620,  0.00400,  0.00200,$
;    0.00100,  0.00050,  0.00040,  0.00000,  0.00000,$
;    0.00000],$
;    [0.50000,  0.50000,  0.50000,  0.50000,  0.50000,$
;    0.50000,  0.50000,  0.50000,  0.50000,  0.50000,$
;    0.44000,  0.30000,  0.20000,  0.13000,  0.09000,$
;    0.06000,  0.05000,  0.04200,  0.03200,  0.02500,$
;    0.02000,  0.01100,  0.00800,  0.00500,  0.00250,$
;    0.00120,  0.00060,  0.00050,  0.00000,  0.00000, 0.0]]
;
;
;  ; ADDED 11/20/11, JY
;  ERYD=13.61
;
;  WWN  =[6.14, 7.30, 7.36, 8.16, $
;    11.03, 11.28,11.52,11.75,$
;    11.97, 8.2,  8.4,  8.9]
;
;  COA =[9.57530e-16, 148.927e-16,  2.85419e-16, 1.39180e-16,$
;    5.24082e-16, 43.9056e-16,  4.72802e-16, 0.182832e-16,$
;    15.3774e-16, 0.245721e-16, 74.1954e-16, 26.2727e-16]
;
;  COB=[2.41725, 3.37549,  3.05068, 3.80445,$
;    2.07034, 3.78085,  2.58944, 1.77218, $
;    3.79464, 2.41297,  5.70981, 6.13567]
;
;  COC=[3.48174, 2.62388,  6.78624, 6.83396, $
;    2.99729, 2.93550,  2.62034, 3.86083, $
;    1.64863, 8.19281,  4.66609, 5.09806]
;
;  COD=[6.42020, 1.48676,  1.31338,  1.35693,  $
;    1.31644, 1.31263,  1.62154,  3.50790,   $
;    1.39476, 3.44784,  0.938729,  0.806361]
;
;  COE=[1.05341e-16, 0.0832431e-16, 2.38554e-16, 1.39227e-16,  $
;    0.0149362e-16,  0.00963622e-16, 2.77677e-16,  9.43598e-16, $
;    0.00373335e-16, 0.0376054e-16,  0.356180e-16, 0.]
;
;  COF=[2.06570,  1.89009,   20.1055,   25.9739,  $
;    0.164347, 0.0386813, 10.3283,   14.6259,   $
;    4.05230,  13.6749,   2.49686,   1.0]
;
;  COG=[6.74498,  21.5238,   11.6571,   11.5631,  $
;    37.0967,  36.8368,   8.19713,   8.34705,  $
;    12.5569,  13.5317,   12.8333,   1.0]
;
;  COH=[1.02078,  2.48385,   4.35359,   1.35315, $
;    9.51190,  8.54345,   1.62767,   1.81088, $
;    7.56624,  1.37048,   0.851090,  1.0]
;
;  ;                 Nothing changed till here****
;
;
;  ; set up energy array
;  ; there are two options, the normal bins are variable and size and 0.05*energy in width
;  ;  up to 10s of keV energies. This means the bin size is much larger than the ionization
;  ;  energy. The program should be Ok with that, but for debugging and other purposes,
;  ;  there is code blocked off by a goto statement that will create a higher resolution
;  ;  energy array where the width never gets larger than 4eV
;  ;
;
;  for n=0,nbins-1 do begin
;    if (n le 20) then begin
;      ener(n) = 0.5 * float(n+1)
;    endif else begin
;      ener(n) = exp(0.05 * float(n+27))
;    endelse
;  endfor
;
;  del(0) = 0.5
;  del(1:*)=ener(1:*)-ener(0:*)
;  ener=ener-del/2.0
;
;  ; Calculate electron impact excitation and ionization cross sections:
;
;  for i=0,nmaj-1 do begin
;    for k=0,9 do begin     ;changing nei to 9 here instead of 14
;      for j=0,nbins-1 do begin
;        if ((ener(j) gt ww(k,i)) and (ww(k,i) gt 0.001)) then begin
;          we = ww(k,i) / ener(j)
;          sigex(k,i,j) = qqn * ao(k,i)*(we^(omeg(k,i)) / (ww(k,i)^2)) * [(1.0 - (we^bb(k,i)))^ anu(k,i)] > 1e-30
;        endif else begin
;          sigex(k,i,j) = 0.0
;        endelse
;        if [(ener(j) gt thi(k,i)) and (thi(k,i) gt 0.001)] then begin
;          ae = ak(k,i)/ener(j) * alog(ener(j)/aj(k,i))
;          gamma = gams(k,i) * ener(j) / (ener(j)+gamb(k,i))
;          t0 = ts(k,i) - (ta(k,i)/(ener(j)+tb(k,i)))
;          sigix(k,i,j) = 1e-16 * ae * gamma* ( atan(((ener(j)-thi(k,i))/2.-t0)/gamma)+atan(t0/gamma) ) > 1e-30
;        endif else begin
;          sigix(k,i,j) = 0.0
;        endelse
;      endfor
;    endfor
;  endfor
;  ;______________________________________________Adding my fitting for N2 electron impact excitation here_______________________________________________________________________________
;
;
;
;
;  ; N2 Excitation States
;  ;                    predissociation
;  ;0  A 3Σu+   6.169
;  ;1  B 3Пg    7.353
;  ;2  W 3∆u    7.362
;  ;3  B’ 3Σu-  8.165
;  ;4  a’ 1Σu-  8.399
;  ;5  a 1Пg    8.549     12%
;  ;6  C 3Пu    11.032
;  ;7  E 3Σg+   11.875
;  ;8  a’’1Σg+  12.255
;  ;9  b 1Пu    12.500    95%
;  ;10  c’4 1Σu+ 12.935   10%
;  ;11  b’ 1Σu+  12.854   84%
;  ;12  w 1∆u    8.890
;  ;13  c 1Пu    12.910   99%
;
;
;
;  ;Eth=[6.169,7.353,7.362,8.165,8.399,8.549,11.032,11.875,12.255,12.500,12.935,12.854,8.890,12.910]
;  ;
;  ;
;  ;p=[ [9.57530,2.41725,3.48174,6.42020,1.05341,2.06570,6.74498,1.02078],$
;  ;  [148.9270,3.375490,2.623880,1.486760,0.083243,1.890090,21.52380,2.483850],$
;  ;  [2.85419,3.05068,6.78624,1.31338,2.38554,20.10550,11.65710,4.35359],$
;  ;  [99.6172,7.9575,0.7916,1.0916,0.5619,2.9250,7.7983,1.2620],$
;  ;  [99.6141,8.1527,0.6147,0.8563,7.5378,4.9560,4.9845,0.8058],$
;  ;  [99.5923,8.4896,0.8923,1.7735,0.5584,1.4904,11.3971,0.8712],$
;  ;  [499.8625,4.9213,3.0122,2.0885,0.1251,-1.4981,10.2382,-0.4407],$
;  ;  [98.4972,19.0442,8.6910,1.4275,0.0032,-0.0748,8.1990,19.9242],$
;  ;  [98.8627,14.6712,8.2492,0.5516,0.0824,0.7443,6.1166,3.1550],$
;  ;  [3.1219,3.0088,4.0583,0.3174,0.1606,1.1773,21.0129,0.7159],$
;  ;  [0.0728,1.2415,33.5635,0.4428,4.9739,3.3763,0.0134,1.9084],$
;  ;  [0.1237,1.3607,22.3639,0.4218,4.9800,4.2974,0.1170,1.8511],$
;  ;  [20.3014,5.9437,5.1718,0.8282,0,0,0,0],$
;  ;  [0.0817,1.1999,27.6399,0.5989,0,0,0,0] ]
;  ;
;  ;
;  ;ER=13.6056980659
;  ;
;  ;
;  ;; Calculate electron impact excitation cross sections:
;  ;
;  ;for k=0,nei-1 do begin
;  ;  for j=0,nbins-1 do begin
;  ;    if ((ener(j) gt Eth(k)) and (Eth(k) gt 0.001)) then begin
;  ;
;  ;
;  ;
;  ;      XE=ener(j)-Eth(k)
;  ;
;  ;      if k lt 12 then $
;  ;        sigex(k,2,j)= ( p(0,k) * (XE/ER)^p(1,k) / (  1.0+   (XE/p(2,k))^(p(1,k)+p(3,k))     )     )  +    $
;  ;        ( p(4,k) * (XE/ER)^p(5,k) / (  1.0+   (XE/p(6,k))^(p(5,k)+p(7,k))     )     )    > 1e-30   else $
;  ;
;  ;        sigex(k,2,j)= ( p(0,k) * (XE/ER)^p(1,k) / (  1.0+   (XE/p(2,k))^(p(1,k)+p(3,k))     )     )      > 1e-30
;  ;
;  ;
;  ;    endif else begin
;  ;      sigex(k,2,j) = 0.0
;  ;    endelse
;  ;
;  ;  endfor
;  ;endfor
;  ;sigex=sigex*1.e-16
;  ;
;  ;
;  ;
;  ;
;  ;
;  ;______________________________________________________________________________________________________________________________________________________________________________________
;
;
;
;  ; for when we want to add in Justin's changes
;  ; ****note that I is not defined in this loop!
;  ;ADDED 11/20/11, JY
;  ;   New cross-sections based on Johnson 2005 and Malone 2009
;  ;for J=0,NBINS-1 do begin
;  ;      IF (ENER(J).GT.WWN(I) .AND. WWN(I).GT.0.001) THEN
;  ;              WEN=ENER(J)-WWN(I)
;  ;              SIGNE(I,J)=COA(I)*((WEN/ERYD)**(COB(I)))*
;  ;     >                 (1.0 + WEN/COC(I))**(-1.*(COB(I) + COD(I)))
;  ;     >                 +COE(I)*((WEN/ERYD)**(COF(I)))*
;  ;     >                 (1.0 + WEN/COG(I))**(-1.*(COF(I) + COH(I)))
;  ;          IF (SIGNE(I,J) .LT. 1.E-30) SIGNE(I,J) = 0.0
;  ;          ELSE
;  ;            SIGNE(I,J) = 0.0
;  ;          ENDIF
;  ;endfor
;
;  ;C  Redefine SIGEX, Zero the hi ener NANS in W and B' 11/20/11 JY
;  ;      DO 143 I=1,NBINS
;  ;      IF (ENER(I).GT.425.) SIGNE(4,I)=0.
;  ;      IF (ENER(I).GT.1150.) SIGNE(3,I)=0.
;
;  ;            SIGEX(1,3,I)=SIGNE(1,I)+SIGNE(2,I)+SIGNE(3,I)
;  ;      SIGEX(2,3,I)=SIGNE(4,I)
;  ;            SIGEX(3,3,I)=SIGNE(5,I)+SIGNE(6,I)+SIGNE(7,I)
;  ;     >    +SIGNE(8,I)+SIGNE(9,I)
;  ;      SIGEX(4,3,I)=SIGNE(10,I)+SIGNE(11,I)+SIGNE(12,I)
;  ; 143     CONTINUE
;
;
;
;  sec=fltarr(nmaj,nbins,nbins) ; (species,lower energy, upper energy) cross section for cascade forming PE with lower energy by
;  ; PE with higher energy and collisions (ionization) with species, cm^2
;  siga=sec ; (species, lost energy, higher energy) cross section for energy loss by PE of higher energy to PE of lower energy
;  ; through collisions (ionizatin & excitation) with species, cm^2, note it is different in format / use from
;  ; sec in that the second dimension is energy lost by electron rather than energy of electron after collision
;  sigloss=fltarr(nmaj,nbins) ; = total loss cross ection
;  ; we consider production and loss due to cascade in energy through ionization & excitation and production of secondaries
;  ; through ionization
;
;  
;  for iprim=0,nbins-1 do begin    ;Loop over energy
;  
;        for inmaj=0,nmaj-1 do begin   ;Loop over species
;          
;          
;          for iionstate=0,ninn(inmaj)-1 do begin
;
;            wag=thi[iionstate,inmaj] ; ionization threshold for state
;            ;               newenergylo= ener(iprim)-del(iprim)/2.0 - thi(iionstate,inmaj)
;            ;               newenergyhi= ener(iprim)+del(iprim)/2.0 - thi(iionstate,inmaj)
;            ;               newenergydel=newenergyhi-newenergylo
;            ;               newenergylo=newenergylo>0.
;            maxaveenergy = ener(iprim) - wag
;            halfenergy = maxaveenergy / 2.0
;
;            itmax=(where(ener[0:iprim] gt halfenergy,/null))[0]
;
;            if halfenergy gt 0.0 then begin
;
;              ; first, calculate the cross section for production of secondaries
;              ; into the bins corresponding to those bins which can contain
;              ; secondaries (bottom 1/2 of possible energies after ionization)
;
;
;              for indx=0,itmax-1 do begin
;                e1=((ener[indx]-del[indx]/2.0)<halfenergy)>0.0
;                e2=(e1 + del[indx])<halfenergy
;                sigsecprod,inmaj,iionstate,ener[iprim],e1,e2,meansecenergy,sigsp
;                ; production of secondaries
;                sec(inmaj,indx,iprim) = sec(inmaj,indx,iprim) + sigsp * del[iprim]/del[indx]
;
;                t12=meansecenergy; ener[indx] ; assumed average energy of the secondaries that are produced
;                wth1 = t12 + wag ; mean energy lossed by primaries
;                eta = ener[iprim] - wth1 ; new energy of primary
;
;                if eta gt 0. then begin
;
;                  ;                      find the bin number closest to energy ETA on grid ENER.
;                  ie=(where(ener[0:iprim] gt eta,/null))[0]
;                  iee=(ie-1)>0  ; loss of energy can span two energy bins
;
;                  k=iprim-ie-1  ; this is the index corresponding to the amount of energy lost by the primary
;                  kk=(iprim-iee-1)>0 ; above - 1
;
;                  if (ie eq iprim) then begin
;
;                    if (iprim eq 0) then $
;                      siga(inmaj,0,iprim)  =  siga(inmaj,0  ,iprim) + sigsp $
;                    else $
;                      siga(inmaj,0,iprim) = siga(inmaj,0,iprim) + sigsp * (del[iprim]/del[iprim-1]))* wth1 / (ener[iprim]-ener[iprim-1])
;
;                  endif else begin
;
;                    if (ie eq 0) then begin
;                      siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp* 2*eta* del[iprim]/(del[0]^2)
;                      siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp* del[iprim]/del[0]
;                    endif else begin
;                      frac =1.0- abs((ener[ie]-eta) / (ener[ie]-ener[iee]))
;                      siga(inmaj,k,iprim) =  siga(inmaj,k,iprim)  + sigsp * (del[iprim]/del[ie])*  frac
;                      siga(inmaj,kk,iprim) = siga(inmaj,kk,iprim) + sigsp * (del[iprim]/del[iee]) * (1.0- frac)
;                    endelse
;
;                  endelse
;
;
;                endif
;
;
;              endfor
;
;            endif
;
;
;
;          endfor   ;ionisation states
;          
;          
;    
;;          _________________________________________________________________________________________________________
;    
;           for iexstate=0,nnn(inmaj)-1 do begin   
;
;
;                sigsp=sigex[iexstate,inmaj,iprim]
;                wag=ww[iexstate,inmaj] ; excitation threshold for state
;;                newenergylo= ener(iprim)-del(iprim)/2.0 - wag
;;                newenergyhi= ener(iprim)+del(iprim)/2.0 - wag
;;                newenergydel=newenergyhi-newenergylo
;;                newenergylo=newenergylo>0.
;;                meannewenergy=(newenergylo+newenergyhi)/2.0
;
;
;                eta = ener[iprim] - wag ; new energy of primary
;                if eta gt 0. then begin
;
;;                      find the bin number closest to energy ETA on grid ENER.
;                        ie=(where(ener[0:iprim] gt eta,/null))[0]
;                        iee=(ie-1)>0  ; loss of energy can span two energy bins
;  
;                        k=iprim-ie-1  ; this is the index corresponding to the amount of energy lost by the primary
;                        kk=(iprim-iee-1)>0 ; above - 1
;                        
;                        if (ie eq iprim) then begin
;                          
;                          if (iprim eq 0) then $
;                             siga(inmaj,0,iprim)  =  siga(inmaj,0  ,iprim) + sigsp $
;                          else $
;                            siga(inmaj,0,iprim) = siga(inmaj,0,iprim) + sigsp * (del[iprim]/del[iprim-1]))* wag / (ener[iprim]-ener[iprim-1])
;                       
;                        endif else begin
;                          
;                          if (ie eq 0) then begin
;                              siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp* 2*eta* del[iprim]/(del[0]^2) 
;                              siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp* del[iprim]/del[0]
;                          endif else begin
;                            frac =1.0- abs((ener[ie]-eta) / (ener[ie]-ener[iee]))
;                           siga(inmaj,k,iprim) =  siga(inmaj,k,iprim)  + sigsp * (del[iprim]/del[ie])*  frac
;                           siga(inmaj,kk,iprim) = siga(inmaj,kk,iprim) + sigsp * (del[iprim]/del[iee]) * (1.0- frac) 
;                          endelse
;                           
;                        endelse
;                        
;                        
;                endif
;
;           endfor  ;excitation states
;        
;;        _______________________________________________________________________________________________________________
;
;
;; calculate total loss cross section (tsa in GLOW)
;          for k=0,(iprim-1)<(nbins-1) do $
;               sigloss[inmaj,iprim] = sigloss[inmaj,iprim]+siga[inmaj,k,iprim]*del[iprim-k-1]/del[iprim]
;
;        
;        endfor     ;species
;  
;  
;  
;  endfor           ;primary energy bins
;  
;  
;
;
;
;  return
;END      
;           
;           
;            
;    
;
;pro sigsecprod,imaj,istate,eprim,e1,e2,meansecenergy,sigsecprod,full_integral=full_integral
;
;  ; pro SIGSECPROD calculate the  cross section for secondary electron production by electron
;  ; impact ionization for species imaj, state index
;  ; istate, primary energy eprim, secondary energy from e1 to e2 (e1 < e2).
;  ; the equation for this calculation is equation 4 from Green and Sawada JASTP, 1972
;
;  ; keyword full_integral integrates the entire cross section from 0 to tm where tm is the max energy of a secondary
;  ; see eqn 5 of Green and Sawada
;
;  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
;  ; see this file for definitions
;
;  a=(ak[istate,imaj]/eprim)*alog(eprim/aj[istate,imaj])*1e-16 ; G&S Eqn 9
;  gamma=gams[istate,imaj]*eprim/(eprim+gamb[istate,imaj]) ; G&S Eqn 8
;  t0=ts[istate,imaj]-(ta[istate,imaj]/(eprim+tb[istate,imaj])) ; G&S Eqn 7
;  tm=(0.5 * (eprim-thi[istate,imaj]))>0. ; G&S Eqn 6
;
;  ABB=(E2-T0)/gamma
;  ABC=(E1-T0)/gamma
;  AL2=gamma*gamma*(ABB*ABB+1.0)
;  AL1=gamma*gamma*(ABC*ABC+1.0)
;  ABD=ATAN(ABB)-ATAN(ABC)
;  T12=T0+0.5*gamma*(ALOG(AL2)-ALOG(AL1))/ABD
;  meansecenergy=t12
;
;  if keyword_set(full_integral) then begin
;
;    sigsecprod= a * gamma * ( atan((tm - t0)/gamma) + atan(t0/gamma) )
;
;  endif else begin
;
;    sigsecprod= (a * gamma * ( atan((e2 - t0)/gamma))) - (a * gamma * ( atan((e1 - t0)/gamma)))
;
;  endelse
;
;  return
;end
;
;; get rid of once this all works
;; now calculate the cross section for cascade in energy due to ionization
;; these are the possibilities in the upper half of possible energy after
;; ionization, but note that the cross section is calculated the same way
;; as the ionization cross section is symmetric (see Green and Sawada, 1972)
;;wbin,ener,del,halfenergy,maxaveenergy,w
;;ind=where(w gt 0.,n_ind)
;;for i=0,n_ind-1 do begin
;;     indx=ind[i]
;;     e1=maxaveenergy - (ener[indx]-del[indx]/2.0)
;;     if e1 gt 0. then begin
;;     e2=e1 + del[indx]
;;    sigion,inmaj,iionstate,ener[iprim],e1,e2,sigi
;;     sec(inmaj,indx,iprim)= sec(inmaj,indx,iprim) + sigi*w[indx] ;* del[iprim]/del[indx] * mask[indx]
;;    siga(inmaj,indx,iprim)=siga(inmaj,indx,iprim) + sigi*w[indx] ;* del[iprim]/del[indx] * mask[indx]
;;     endif
;;     endfor
;;if iprim mod 100 eq 0 then print,iprim,nbins,inmaj
;;if iprim eq 3000 then stop
;;sige = sigex[iexstate,inmaj,iprim]
;
;

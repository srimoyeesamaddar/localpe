pro ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi

; This procedure uses tabulated values of elastic cross section parameters for o, o2, and n2
; and places them onto the chosen energy grid

; this procedure is converted directly from fortran to IDL from the GLOW model
; of Stan Solomon, specifically etrans.f; see GLOW documentation for further details

EC=[[1.00,     2.00,     4.00,     6.00,     8.00,$
                    10.00,    12.00,    14.00,    16.00,    18.00,$
                    20.00,    30.00,    40.00,    50.00,    60.00,$
                    70.00,    80.00,    90.00,   100.00,   150.00,$
                   200.00,   300.00,   500.00,  1000.00,  2000.00,$
                  3000.00,  5000.00, 10000.00, 20000.00, 40000.00,$
                 50000.00],$
                     [1.00,     2.00,     3.00,     5.00,     7.00,$
                    10.00,    15.00,    20.00,    30.00,    40.00,$
                    50.00,    70.00,   100.00,   150.00,   200.00,$
                   300.00,   400.00,   500.00,   600.00,   700.00,$
                  1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
                 20000.00, 40000.00, 50000.00,     0.00,     0.00,$
                     0.00],$
                     [1.00,     2.00,     2.50,     3.00,     4.00,$
                     5.00,     6.00,     8.00,    10.00,    15.00,$
                    20.00,    30.00,    40.00,    50.00,    70.00,$
                   100.00,   200.00,   300.00,   500.00,   700.00,$
                  1000.00,  2000.00,  3000.00,  5000.00, 10000.00,$
                 20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0]]
                 
CC=[[ 5.00E-16, 6.00E-16, 7.50E-16, 7.60E-16, 7.70E-16,$
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
                 [9.00E-16, 2.27E-15, 2.52E-15, 1.93E-15, 1.32E-15,$
                 1.15E-15, 1.16E-15, 1.17E-15, 1.18E-15, 1.14E-15,$
                 1.13E-15, 9.50E-16, 8.60E-16, 7.30E-16, 5.90E-16,$
                 4.70E-16, 3.30E-16, 2.50E-16, 1.60E-16, 1.30E-16,$
                 1.10E-16, 6.35E-17, 4.18E-17, 2.54E-17, 1.28E-17,$
                 6.44E-18, 3.27E-18, 2.62E-18, 0.00E+00, 0.00E+00, 0.0]]
                 
CE=[[0.50000,  0.49500,  0.46800,  0.43600,  0.42000,$
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
                  [0.50000,  0.50000,  0.50000,  0.49000,  0.46800,$
                  0.44500,  0.43600,  0.42000,  0.40500,  0.36800,$
                  0.34300,  0.31600,  0.28900,  0.25800,  0.22000,$
                  0.18400,  0.14000,  0.11000,  0.08400,  0.07400,$
                  0.06300,  0.03400,  0.02400,  0.01500,  0.00740,$
                  0.00370,  0.00180,  0.00140,  0.00000,  0.00000, 0.0]]
     
CI=[[ 0.60000,  0.60000,  0.60000,  0.60000,  0.60000,$
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
 
; Interpolate elastic cross sections and backscatter ratios:
 
 sigs=fltarr(nmaj,nbins)
 pe=fltarr(nmaj,nbins)
 pi=fltarr(nmaj,nbins)
 a10ener=alog10(ener)
 
 num=lonarr(nmaj)
 
 for i=0,nmaj-1 do begin
  num[i]=total(ci[*,i] gt 0.)
  sigs[i,*]=interpol(alog10(cc[0:num[i]-1,i]),alog10(ec[0:num[i]-1,i]),a10ener,/spline)
   pe[i,*]=interpol(alog10(ce[0:num[i]-1,i]),alog10(ec[0:num[i]-1,i]),a10ener,/spline)
    pi[i,*]=interpol(alog10(ci[0:num[i]-1,i]),alog10(ec[0:num[i]-1,i]),a10ener,/spline)
    sigs[i,*]=10.^sigs[i,*]
    pe[i,*]=10.^pe[i,*]
    pi[i,*]=10.^pi[i,*]
 endfor
 
 return
 end
 
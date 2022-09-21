
;###########################################################################################################
;###########################################################################################################

;Srimoyee- NEW COMPILATION CODE STARTS HERE


pro localpe_setup_exsect_new,demo=demo

  ; Calculates electron impact cross sections
  ; New compilation of cross section 

  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
  ; see this file for definitions


  ; define array dimensions
  ; nbins=190  ; number of electron energy bins
  ;nmaj=3                 ; number of major species in atmosphere (o, o2, and n2)
  nei=20                  ; maximum number of product states for excitation, ionization

  

  ; cross section arrays
  sigex=fltarr(nei,nmaj,nbins)
  sigix=fltarr(nei,nmaj,nbins)
  iimaxx=fltarr(nbins)

  ener=fltarr(nbins)
  del=fltarr(nbins)
  sigi=fltarr(nbins)
  t12=fltarr(nbins)

  nnn=[14,8,20]  ;Excitation states 
  ninn=[1,2,2]   ;Ionization states
  
;  Excitation thresholds for different states

  ww=[[9.5,1.96, 12.08, 12.53, 14.11, 4.18, 15.65, 10.98, 10.73, 9.14, 12.75, 13, 16.08, 14, 0.,0., 0.,0., 0.,0.],$
      [0.2,4.5,6.12,10.,10.3,16.,0.977,1.627,0.0,0.0, 0.,0., 0.,0.,0., 0.,0., 0.,0.,0.],$
      [6.169,7.353,7.362,8.165,8.399,8.549,11.032,11.875,12.255,12.500,12.935,12.854,16.4,17.4,11,8.890,12.910,23.7,40,1.85]]

;    Ionisation thresholds for different states

  thi=[[13.60,0.0],$
       [12.10,20.0],$
       [15.581,22.0]]
    

;_____________________________________________________________________________________________________________________________
  ; set up energy array
  ; there are two options, the normal bins are variable and size and 0.05*energy in width
  ;  up to 10s of keV energies. This means the bin size is much larger than the ionization
  ;  energy. The program should be Ok with that, but for debugging and other purposes,
  ;  there is code blocked off by a goto statement that will create a higher resolution
  ;  energy array where the width never gets larger than 4eV
  ;

  for n=0,nbins-1 do begin
    if (n le 20) then begin
      ener(n) = 0.5 * float(n+1)
    endif else begin
      ener(n) = exp(0.05 * float(n+27))
    endelse
  endfor

  del(0) = 0.5
  del(1:*)=ener(1:*)-ener(0:*)
  ener=ener-del/2.0
  
 ;_________________________________________________Calculate electron impact excitation and ionization cross sections____________________________________________________________________________
 


  ;#1 N2 Excitation States
  ;                    predissociation
  ;0  A 3Σu+   6.169
  ;1  B 3Пg    7.353
  ;2  W 3∆u    7.362
  ;3  B’ 3Σu-  8.165
  ;4  a’ 1Σu-  8.399
  ;5  a 1Пg    8.549     12%
  ;6  C 3Пu    11.032    50%
  ;7  E 3Σg+   11.875
  ;8  a’’1Σg+  12.255
  ;9  b 1Пu    12.500    95%
  ;10  c’4 1Σu+ 12.935   10%
  ;11  b’ 1Σu+  12.854   84%
  ;15  w 1∆u    8.890
  ;16  c 1Пu    12.910   99%
  ;
  ;Other loss channels added to account for dissociation due to electron impact- Strickland 1996
  ;12  15.8 eV peak 16.4 eV  100%
  ;13  17.3 eV peak  17.4eV 100%
  ;14  triplet manifold  11eV 100%
  ;
  ;17  VUV       23.7 eV 100%
  ;18  Ryd atoms 40 eV 100%
  ;
  ;
  ;19 vib : NOT changing this state
  
  
  ;#2  O2 Excitation States
  

  ;0 vib          0.2 eV
  ;1 AA'c         4.5 eV     prediss
  ;2 SR           6.12       prediss
  ;3 LB           10         prediss
  ;4 SB           10.3       prediss
  ;5 Ryds         16         prediss
  ;6 a1del_g      0.977
  ;7 b1Sigma_g    1.627


  ;#3  O Excitation States
  
; 0   3s^3S^o      9.5 eV 
; 1   2p^4 ^1D     1.96
; 2   3d^3D^o      12.08
; 3   3s'^3D^o     12.53
; 4   3s''^3P^o    14.11 
; 5   2p^4 ^1S     4.18 
; 6   2p^5P^o      15.65
; 7  3p^5P        10.73
; 8  4d^3D^o      12.75
; 9  5d^3D^o      13
; 10  Rydbergs     14
; 11   3p^3P        10.98
; 12   3s^5S^o      9.14
; 13 4d'^3P^o     16.08

p_O=[[0.1186, 1.0238,13.8980,1.1150,0.0274,0.8784,52.6275,0.5161],$
     [0.0469,14.0079,6.5948,0.9170,495.2045,3.7869,2.4143,0.7713],$
     [0.1185,0.8416,3.7233,0.4250,0.0179,1.7335,20.1638,1.0148],$
     [0.0478,0.7445,37.1416,0.6746,0.1312,0.8873,4.2165,0.8589],$
     [0.0277,1.4506,32.4527,0.8741,0.1872,1.0028,2.6481,0.3287],$
     [0.0032,0.6338,20.0968,2.4291,0.1812,1.0276,3.8687,0.9775],$
     [0.0037,5.1219,18.5104,1.6939,0.1107,0.6327,35.0869,0.8162],$
     [0.3638,1.1726,6.6056,1.0698,0.0,0.0,1.0,0.0],$
     [0.1538,4.8487,7.2615,1.6631,0.0324,0.2460,11.9178,2.5099],$
     [98.4005,7.1143,4.8595,2.2804,0.0,0.0,1.0,0.0],$
     [0.0180,0.8349,34.5364,0.9812,9.8181,3.8669,1.6838,3.3599],$
     [0.0045,1.6012,27.9098,0.9244,0.0108,0.5718,12.6714,0.6722],$
     [0.0006,7.0204,23.9377,0.5416,0.0,0.0,1.0,0.0],$
     [1.1437,1.0298,7.3643,0.4891,9.4479,11.5893,0.2458,10.4649]]     

  
  p_N2=[ [9.57530,2.41725,3.48174,6.42020,1.05341,2.06570,6.74498,1.02078],$
         [148.9270,3.375490,2.623880,1.486760,0.083243,1.890090,21.52380,2.483850],$
         [2.85419,3.05068,6.78624,1.31338,2.38554,20.10550,11.65710,4.35359],$
         [99.6172,7.9575,0.7916,1.0916,0.5619,2.9250,7.7983,1.2620],$
         [99.6141,8.1527,0.6147,0.8563,7.5378,4.9560,4.9845,0.8058],$
         [99.5923,8.4896,0.8923,1.7735,0.5584,1.4904,11.3971,0.8712],$
         [499.8625,4.9213,3.0122,2.0885,0.1251,-1.4981,10.2382,-0.4407],$
         [98.4972,19.0442,8.6910,1.4275,0.0032,-0.0748,8.1990,19.9242],$
         [98.8627,14.6712,8.2492,0.5516,0.0824,0.7443,6.1166,3.1550],$
         [3.1219,3.0088,4.0583,0.3174,0.1606,1.1773,21.0129,0.7159],$
         [0.0728,1.2415,33.5635,0.4428,4.9739,3.3763,0.0134,1.9084],$
         [0.1237,1.3607,22.3639,0.4218,4.9800,4.2974,0.1170,1.8511],$
         [0.3908,1.1866,15.6934,0.5019,0.0326,1.4785,2.4783,13.1506],$
         [0.1848,1.0904,14.5947,0.4722,1.3188,1.6956,0.7298,8.9221],$
         [99.0223,6.5109,4.9607,1.4943,2.9187,2.0930,3.7883,3.6823],$
         [20.3014,5.9437,5.1718,0.8282,0,0,1,0],$
         [0.0817,1.1999,27.6399,0.5989,0,0,1,0] ,$
         [0.0527,1.2101,58.9775,0.7044,0,0,1,0],$
         [0.0237,0.5292,100.2524,1.0136,0,0,1,0]]
  
  
  p_O2=[[12.9375,7.5566, 10.1148,13.5803,37.1013,5.8420,6.2180,0.8618],$
    [25.4093,7.5201,5.7868,50.1096,0.9881,1.9706,6.5458,1.2318],$
    [99.5800,7.8081,0.4631,0.7547,0.8460,0.5494,26.5331,0.7434],$
    [99.5343,9.6625,4.7581,0.6902,0.0624,0.2952,31.1376,0.8421],$
    [99.5481,9.0000,3.2491,0.5718,0.0144,0.5148,41.2214,0.7297],$
    [3.0477,1.5724,12.6191,0.6797,9.9951,1.1743,0.0018,4.1789],$
    [1.7862,2.3891,5.4963,1.7106,0.0,0.0,1.0,0.0],$
    [5.7594,3.8908,4.2622,1.7539,0.0,0.0,1.0,0.0]]
    
    p=fltarr(8,nei,nmaj)
    p(*,0:13,0)=p_O
    p(*,0:7,1)=p_O2
    p(*,0:18,2)=p_N2
  
  ER=13.6056980659
  
  
  ; Calculate electron impact excitation cross sections:
   for i=0,nmaj-1 do begin
       for k=0,nnn[i]-1 do begin
          for j=0,nbins-1 do begin
                 if ((ener(j) gt ww(k,i)) and (ww(k,i) gt 0.001)) then begin
  
                      XE=ener(j)-ww(k,i)
                      sigex(k,i,j)= ( ( p(0,k,i) * (XE/ER)^p(1,k,i) / (  1.0+   (XE/p(2,k,i))^(p(1,k,i)+p(3,k,i))     )     )  +    $
                                1    ( p(4,k,i) * (XE/ER)^p(5,k,i) / (  1.0+   (XE/p(6,k,i))^(p(5,k)+p(7,k,i))     )     )  )    *1.e-16 > 1.e-30    
                 
                 endif else begin
                      sigex(k,i,j) = 0.0
               endelse
  
          endfor
       endfor
   endfor
  

;___________________________________________________________________________________________________________________________________________ 
;  Keeping the vibrational excitation state of N2 same as in GLOW model
  qqn=6.51e-14
  
  ao=[[.0100,.0042,.1793,.3565,.0817,.0245,.0293,.1221, 0.,0.],$
      [.0797,.0211,.0215,.3400,.0657,1.110,3.480, 0.00, 0.,0.],$
      [2.770,.1140,.1790,.0999,.8760,.6010,1.890,1.350, 0.,0.]]

 omeg=[[1.00, 1.00, 3.00, 0.75, 3.00, 0.85, 0.75, 0.75, 0.,0.],$
       [2.00, 2.00, 1.15, 0.75, 0.75, 0.75, 7.00, 0.00, 0.,0.],$
       [3.00, 3.00, 3.00, 1.00, 0.75, 0.75, 0.75, 8.00, 0.,0.]]

anu= [[2.00, 1.04, 2.53, 0.54, 2.43, 2.87, 0.93, 0.72, 0.,0.],$
      [6.18, 4.14, 1.00, 1.05, 1.60, 3.00,10.87, 0.00, 0.,0.],$
      [4.53, 4.78, 4.32, 4.05, 1.47, 1.27, 3.00, 1.58, 0.,0.]]
      
bb=[[1.00, 0.50, 1.02, 0.01, 4.19, 4.88, 0.66, 0.17, 0.,0.],$
    [0.53, 0.51, 0.98, 0.99, 1.86, 1.00, 1.00, 0.00, 0.,0.],$
    [1.42, 3.54,12.70, 5.20, 0.86, 0.45, 1.00, 1.00, 0.,0.]]

sigex(19,2,*)=0.0
for j=0,nbins-1 do begin
        if ((ener(j) gt ww(19,2)) and (ww(19,2) gt 0.001)) then begin
          we = ww(19,2) / ener(j)
          sigex(19,2,j) = qqn * ao(7,2)*(we^(omeg(7,2)) / (ww(19,2)^2)) * [(1.0 - (we^bb(7,2)))^ anu(7,2)] > 1.e-30
        endif else begin
          sigex(19,2,j) = 0.0
        endelse
endfor

;___________________________________________________________________________________________________________________________________________

;N2, O and O2 ionisation cross-section parameters:
;N2/O2 ionisation:
;;Single dissociation -Non dissociative
;;N+N_plus -Dissociative
;O ionisation

               
 ak=[[ 4.8787, 0.0],$
     [1.9634, 0.3287],$
     [6.4298, 1.4726]]

 aj=[[11.1297, 0.0],$
     [8.9237, 0.5764],$
     [14.2167, 22.3963]]

 ts=[[ 17.4974,0.],$
    [-48.6810, 13.5476],$
    [-49.7604, -44.8376]]

 ta=[[3450.8,0.],$
    [-7.6249,8373.8],$
    [830.0339,-6192.0]]

 tb= [[121.6107,0.],$
     [-35.1641, 156.8474],$
     [2182.8,398.1121]]

 gams= [[8.1247,0.],$
       [65.5480,47.8684],$
       [41.3117,32.8228]]

gamb= [[ -4.1235,0.],$
       [ -40.3355,-63.1845],$
       [ -23.3569,-42.3218]]
       
       
; Calculate electron impact ionization cross sections:

for i=0,nmaj-1 do begin
       for k=0,ninn[i]-1  do begin     
         for j=0,nbins-1 do begin
           
           if [(ener(j) gt thi(k,i)) and (thi(k,i) gt 0.001)] then begin
             ae = ak(k,i)/ener(j) * alog(ener(j)/aj(k,i))
             gamma = gams(k,i) * ener(j) / (ener(j)+gamb(k,i))
             t0 = ts(k,i) - (ta(k,i)/(ener(j)+tb(k,i)))
             sigix(k,i,j) = 1.e-16 * ae * gamma* ( atan(((ener(j)-thi(k,i))/2.-t0)/gamma)+atan(t0/gamma) ) > 1.e-30
           endif else begin
             sigix(k,i,j) = 0.0
           endelse
         endfor
       endfor
     endfor

;  ______________________________________________________________________________________________________________________



  sec=fltarr(nmaj,nbins,nbins) ; (species,lower energy, upper energy) cross section for cascade forming PE with lower energy by
  ; PE with higher energy and collisions (ionization) with species, cm^2
  siga=sec ; (species, lost energy, higher energy) cross section for energy loss by PE of higher energy to PE of lower energy
  ; through collisions (ionizatin & excitation) with species, cm^2, note it is different in format / use from
  ; sec in that the second dimension is energy lost by electron rather than energy of electron after collision
  sigloss=fltarr(nmaj,nbins) ; = total loss cross section
  ; we consider production and loss due to cascade in energy through ionization & excitation and production of secondaries
  ; through ionization

  ; start the key loops
  for iprim=nbins-1,0,-1 do begin
    kuk=0 ; track max value (index of) of secondary energy so we don't have to loop over all energies later
    for inmaj=0,nmaj-1 do begin


      mask=fltarr(nbins)+1.0 ; mask is used to make sure we get the production/loss within a large bin correct
      ;mask[iprim]=0. ; there is loss from a bin due to the creation of secondaries, but for no other

      for iionstate=0,ninn(inmaj)-1 do begin

        newenergylo= ener(iprim)-del(iprim)/2.0 - thi(iionstate,inmaj)
        newenergyhi= ener(iprim)+del(iprim)/2.0 - thi(iionstate,inmaj)
        newenergydel=newenergyhi-newenergylo
        newenergylo=newenergylo>0.
        maxaveenergy = ener(iprim) - thi(iionstate,inmaj)
        halfenergy = maxaveenergy / 2.0

        if maxaveenergy gt 1.0 then begin

          ; first, calculate the cross section for production of secondaries
          ; into the bins corresponding to those bins which can contain
          ; secondaries (bottom 1/2 of possible energies after ionization)
          wbin,ener,del,0.0,halfenergy,w
          ind=where(w gt 0.,n_ind)
          abin,ener,halfenergy,mw,f
          n_ind=mw
          for i=0,n_ind-1 do begin
            indx=i;ind[i]
            e1=((ener[indx]-del[indx]/2.0)<halfenergy)>0.0
            e2=(e1 + del[indx])<halfenergy
            sigsecprod,inmaj,iionstate,ener[iprim],e1,e2,meansecenergy,sigsp
            ; production of secondaries
            sec(inmaj,indx,iprim) = sec(inmaj,indx,iprim) + sigsp * del[iprim]/del[indx]

            ; loss of primary energy through ionization / creation of secondaries
            wag=thi[iionstate,inmaj] ; ionization threshold for state
            t12=meansecenergy; ener[indx] ; assumed average energy of the secondaries that are produced

            wth1 = t12 + wag ; mean energy lossed by primaries
            eta = ener[iprim] - wth1 ; new energy of primary
            abin,ener,eta,ie,frac
            iee=(ie-1)>0  ; loss of energy can span two energy bins
            k=iprim-ie-1  ; this is the index corresponding to the amount of energy lost by the primary
            kk=(iprim-iee-1)>0 ; above - 1

            ; if the primary does not change energy bins, then put the energy loss in bin 0
            ; this is neccessary because siga is essentailly tracking energy loss by bin - not absolute energy
            if ie ge iprim or (k lt 0) or (kk lt 0) then begin
              k=0
              kk=0
              ie=(iprim-1)>0
              iee=ie
              frac=wth1/(ener[iprim]-ener[ie])
              ;     if ~finite(frac) then frac=1.
              siga(inmaj,0,iprim)  =  siga(inmaj,0  ,iprim) + sigsp*del[iprim]/del[ie]*frac
              ;      if ~finite(frac) then print,'Fraction not finite',frac
              ;if iprim eq 66 then print,sigsp*del[iprim]/del[ie]*frac
            endif else begin
              siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp*del[iprim]/del[ie  ]*frac
              siga(inmaj,kk,iprim) = siga(inmaj,kk,iprim) + sigsp*del[iprim]/del[iee]*(1.0-frac)
              ;      if ~finite(frac) then print,'Fraction infinite',frac
            endelse
            ;stop
            ;if ener[indx] ge 15. then stop
            kuk=kuk>k  ; keep track of the largest value that k ever gets to
          endfor

        endif

      endfor ; ion states

      ; cascade due to excitation
      ; use a similar proportional division approach to that above for ionization cascade

      for iexstate=0,nnn(inmaj)-1 do begin

        sigsp=sigex[iexstate,inmaj,iprim]
        wag=ww[iexstate,inmaj] ; excitation threshold for state
        newenergylo= ener(iprim)-del(iprim)/2.0 - wag
        newenergyhi= ener(iprim)+del(iprim)/2.0 - wag
        newenergydel=newenergyhi-newenergylo
        newenergylo=newenergylo>0.
        meannewenergy=(newenergylo+newenergyhi)/2.0

        if newenergyhi ge 0. then begin

          eta = ener[iprim] - wag ; new energy of primary
          abin,ener,eta,ie,frac
          iee=(ie-1)>0  ; loss of energy can span two energy bins
          k=iprim-ie-1  ; this is the index corresponding to the amount of energy lost by the primary
          kk=(iprim-iee-1)>0 ; above - 1

          ; if the primary does not change energy bins, then put the energy loss in bin 0
          ; this is neccessary because siga is essentailly tracking energy loss by bin - not absolute energy
          ;     print, ener[iprim],ener[ie];
          if ie ge iprim or (k lt 0) or (kk lt 0) then begin
            k=0
            kk=0
            ie=(iprim-1)>0
            iee=ie

            frac=wag/(ener[iprim]-ener[ie])
            ;        if ~finite(frac) then frac=1.
            siga(inmaj,0,iprim)  =  siga(inmaj,0  ,iprim) + sigsp*del[iprim]/del[ie]*frac
            ;        if ~finite(siga(inmaj,0,iprim)) then print,'0 not nfinite',frac

          endif else begin
            siga(inmaj,k,iprim)  =  siga(inmaj,k  ,iprim) + sigsp*del[iprim]/del[ie]*frac
            siga(inmaj,kk,iprim) = siga(inmaj,kk,iprim) + sigsp*del[iprim]/del[iee]*(1.0-frac)
            ;        if ~finite(siga(inmaj,k,iprim)) then print,'k infinite',frac
            ;        if ~finite(siga(inmaj,kk,iprim)) then print,'kk infinite',frac
          endelse

          kuk=kuk>k  ; keep track of the largest value that k ever gets to

        endif

      endfor ; excitation state

      ; calculate total loss cross section (tsa in GLOW)
      for k=0,(iprim-1)<(nbins-1) do sigloss[inmaj,iprim] = sigloss[inmaj,iprim]+siga[inmaj,k,iprim]*del[iprim-k-1]/del[iprim]

    endfor ; species

  endfor ; primary energy

  first_exsect=1
  ;stop
  print,'Have completed setup_exsect...'

  save,siga,sigloss,sec,sigix,sigex,ener,filename='/home/srimoyee/Desktop/nrl_files/sav_files/sigasiglossNC.sav'

  return


END

pro sigsecprod,imaj,istate,eprim,e1,e2,meansecenergy,sigsecprod,full_integral=full_integral

  ; pro SIGSECPROD calculate the  cross section for secondary electron production by electron
  ; impact ionization for species imaj, state index
  ; istate, primary energy eprim, secondary energy from e1 to e2 (e1 < e2).
  ; the equation for this calculation is equation 4 from Green and Sawada JASTP, 1972

  ; keyword full_integral integrates the entire cross section from 0 to tm where tm is the max energy of a secondary
  ; see eqn 5 of Green and Sawada

  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
  ; see this file for definitions

  a=(ak[istate,imaj]/eprim)*alog(eprim/aj[istate,imaj])*1e-16 ; G&S Eqn 9
  gamma=gams[istate,imaj]*eprim/(eprim+gamb[istate,imaj]) ; G&S Eqn 8
  t0=ts[istate,imaj]-(ta[istate,imaj]/(eprim+tb[istate,imaj])) ; G&S Eqn 7
  tm=(0.5 * (eprim-thi[istate,imaj]))>0. ; G&S Eqn 6

  ABB=(E2-T0)/gamma
  ABC=(E1-T0)/gamma
  AL2=gamma*gamma*(ABB*ABB+1.0)
  AL1=gamma*gamma*(ABC*ABC+1.0)
  ABD=ATAN(ABB)-ATAN(ABC)
  T12=T0+0.5*gamma*(ALOG(AL2)-ALOG(AL1))/ABD
  meansecenergy=t12




  if keyword_set(full_integral) then begin

    sigsecprod= a * gamma * ( atan((tm - t0)/gamma) + atan(t0/gamma) )
    print,'set'

  endif else begin

    sigsecprod= (a * gamma * ( atan((e2 - t0)/gamma))) - (a * gamma * ( atan((e1 - t0)/gamma)))
    ;                  print,'not set'
  endelse

  


  return
end




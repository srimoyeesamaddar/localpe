; N2 plots for presentation
; Comparing the NRL and EUVAC binning

;NOTE:Plot code for sample_localpe_v6.pro
;
;____________________________________________________________________________________________________________________________________________________________________
; N2 Excitation States
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

;States added to account for dissociation due to electron impact- Strickland 1996
;12  15.8 eV peak 16.4 eV  100%
;13  17.3 eV peak  17.4eV 100%
;14  triplet manifold  11eV 100%

;17  VUV       23.7 eV 100%
;18  Ryd atoms 40 eV 100%


;N2 ionisation
;0 Single dissociation -Non dissociative
;1 N+N_plus -Dissociative

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
;____________________________________________________________________________________________________________________________________________________________________

pro dissstates_NEW,Sp_exctrate,nei_N2,nei_O2,Sp_ionzrate,nii_N2,nii_O2,Sp_photoki,$
  N2_dissext,N2_dissionz,O2_dissext,O2_dissionz
  ; Counts only the states needed for dissociation excitation and ionisation

  if n_elements(where(~finite(Sp_exctrate),/null)) gt 0 then stop
  if n_elements(where(~finite(Sp_ionzrate),/null)) gt 0 then stop
  if n_elements(where(~finite(Sp_photoki),/null)) gt 0 then stop


  jmax=(size(Sp_exctrate))[3]
  nbins=(size(Sp_exctrate))[4]
  lmax=(size(Sp_exctrate))[5]

  ;N2 states
  N2_dissext=fltarr(nei_N2,jmax,nbins,lmax)
  N2_dissionz=fltarr(nii_N2,jmax,nbins,lmax)
  

  N2_dissext[0,*,*,*]=  reform(Sp_exctrate[5,2,*,*,*])
  N2_dissext[1,*,*,*]=  reform(Sp_exctrate[6,2,*,*,*])
  N2_dissext[2,*,*,*]=  reform(Sp_exctrate[9,2,*,*,*])
  N2_dissext[3,*,*,*]=  reform(Sp_exctrate[10,2,*,*,*])
  N2_dissext[4,*,*,*]=  reform(Sp_exctrate[11,2,*,*,*])
  N2_dissext[5,*,*,*]=  reform(Sp_exctrate[12,2,*,*,*])
  N2_dissext[6,*,*,*]=  reform(Sp_exctrate[13,2,*,*,*])
  N2_dissext[7,*,*,*]=  reform(Sp_exctrate[14,2,*,*,*])
  N2_dissext[8,*,*,*]=  reform(Sp_exctrate[16,2,*,*,*])
  N2_dissext[9,*,*,*]=  reform(Sp_exctrate[17,2,*,*,*])
  N2_dissext[10,*,*,*]= reform(Sp_exctrate[18,2,*,*,*])

  N2_dissionz[0,*,*,*]= reform(Sp_ionzrate[1,2,*,*,*])



  ;O2 states

  O2_dissext=fltarr(nei_O2,jmax,nbins,lmax)
  O2_dissionz=fltarr(nii_O2,jmax,nbins,lmax)
 

  O2_dissext[0,*,*,*]=  reform(Sp_exctrate[1,1,*,*,*])
  O2_dissext[1,*,*,*]=  reform(Sp_exctrate[2,1,*,*,*])
  O2_dissext[2,*,*,*]=  reform(Sp_exctrate[3,1,*,*,*])
  O2_dissext[3,*,*,*]=  reform(Sp_exctrate[4,1,*,*,*])
  O2_dissext[4,*,*,*]=  reform(Sp_exctrate[5,1,*,*,*])



  O2_dissionz[0,*,*,*]= reform(Sp_ionzrate[1,1,*,*,*])

  



end

;____________________________________________________________________________________________________________________________________________________________________



;pro dissstates_OLD,Sp_exctrate,nei_N2,nei_O2,Sp_ionzrate,nii_N2,nii_O2,Sp_photoki,pii_N2,pii_O2,$
;  N2_dissext,N2_dissionz,N2_phtndissionz,O2_dissext,O2_dissionz,O2_phtndissionz
;  ; Counts only the states needed for dissociation excitation and ionisation
;  if n_elements(where(~finite(Sp_exctrate),/null)) gt 0 then stop
;  if n_elements(where(~finite(Sp_ionzrate),/null)) gt 0 then stop
;  if n_elements(where(~finite(Sp_photoki),/null)) gt 0 then stop
;
;  jmax=(size(Sp_exctrate))[3]
;  nbins=(size(Sp_exctrate))[4]
;  lmax=(size(Sp_exctrate))[5]
;
;  ;N2 states
;  N2_dissext=fltarr(nei_N2,jmax,nbins,lmax)
;  N2_dissionz=fltarr(nii_N2,jmax,nbins,lmax)
;  N2_phtndissionz=fltarr(pii_N2,jmax,lmax)
;
;  N2_dissext[0,*,*,*]=  reform(Sp_exctrate[4,2,*,*,*])
;  N2_dissext[1,*,*,*]=  reform(Sp_exctrate[5,2,*,*,*])
;  N2_dissext[2,*,*,*]=  reform(Sp_exctrate[6,2,*,*,*])
;
;  N2_dissionz[0,*,*,*]= reform(Sp_ionzrate[3,2,*,*,*])
;  N2_dissionz[1,*,*,*]= reform(Sp_ionzrate[4,2,*,*,*])
;  N2_dissionz[2,*,*,*]= reform(Sp_ionzrate[5,2,*,*,*])
;
;  N2_phtndissionz[0,*,*]=  reform(Sp_photoki[2,3,*,*])
;  N2_phtndissionz[1,*,*]=  reform(Sp_photoki[2,4,*,*])
;  N2_phtndissionz[2,*,*]=  reform(Sp_photoki[2,5,*,*])
;  N2_phtndissionz[3,*,*]=  reform(Sp_photoki[2,6,*,*])
;  N2_phtndissionz[4,*,*]=  reform(Sp_photoki[2,7,*,*])
;
;
;  ;O2 states
;
;  O2_dissext=fltarr(nei_O2,jmax,nbins,lmax)
;  O2_dissionz=fltarr(nii_O2,jmax,nbins,lmax)
;  O2_phtndissionz=fltarr(pii_O2,jmax,lmax)
;
;  O2_dissext[0,*,*,*]=  reform(Sp_exctrate[2,1,*,*,*])
;  O2_dissext[1,*,*,*]=  reform(Sp_exctrate[3,1,*,*,*])
;  O2_dissext[2,*,*,*]=  reform(Sp_exctrate[4,1,*,*,*])
;  O2_dissext[3,*,*,*]=  reform(Sp_exctrate[5,1,*,*,*])
;
;
;  O2_dissionz[0,*,*,*]= reform(Sp_ionzrate[4,1,*,*,*])
;  O2_dissionz[1,*,*,*]= reform(Sp_ionzrate[5,1,*,*,*])
;  O2_dissionz[2,*,*,*]= reform(Sp_ionzrate[6,1,*,*,*])
;
;  O2_phtndissionz[0,*,*]=  reform(Sp_photoki[1,3,*,*])
;  O2_phtndissionz[1,*,*]=  reform(Sp_photoki[1,4,*,*])
;  O2_phtndissionz[2,*,*]=  reform(Sp_photoki[1,5,*,*])
;  O2_phtndissionz[3,*,*]=  reform(Sp_photoki[1,6,*,*])
;  O2_phtndissionz[4,*,*]=  reform(Sp_photoki[1,7,*,*])
;
;
;
;end
;
;

;____________________________________________________________________________________________________________________________________________________________________


pro totionzrate,Sp_ionzrate,Sp_photoi,eiionz_1,eiionz_2;,photoi_1
  nei=(size(Sp_ionzrate))[1]
  nmaj=(size(Sp_ionzrate))[2]
  jmax=(size(Sp_ionzrate))[3]
  nbins=(size(Sp_ionzrate))[4]
  lmax=(size(Sp_ionzrate))[5]




  eiionz_1 = fltarr(nmaj,jmax)
  eiionz_2 = fltarr(nmaj,jmax,lmax)

  ;  photoi_1=reform(total(Sp_photoi,2))


  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude
        for sp=0,nmaj-1 do begin   ;species
          for s=0,nei-1 do begin  ;states

            eiionz_1[sp,j]=eiionz_1[sp,j]+Sp_ionzrate[s,sp,j,e,l]
            eiionz_2[sp,j,l]=eiionz_2[sp,j,l]+Sp_ionzrate[s,sp,j,e,l]



          endfor

        endfor
      endfor
    endfor

  endfor




end



;____________________________________________________________________________________________________________________________________________________________________


pro disscontrb,Sp_exctrate,Sp_ionzrate,$
  Sp_exctrate_1,Sp_exctrate_2,Sp_exctrate_3,Sp_exctrate_5,$
  Sp_ionzrate_1,Sp_ionzrate_2,Sp_ionzrate_3,Sp_ionzrate_5,$
  Sp_prcnttotdisionz,Sp_prcnttotdisexct,$
  Sp_contrdexct,Sp_contrdi,$
  St_contrdexct

  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  nii=(size(Sp_ionzrate))[1]

  pii=(size(Sp_phtndissionz))[1]


  Sp_exctrate_1 = fltarr(nei,jmax,nbins)    ;DE
  Sp_exctrate_2 = fltarr(jmax,nbins)
  Sp_exctrate_3 = fltarr(jmax)
  Sp_exctrate_4 = fltarr(nei,jmax)
  Sp_exctrate_5 = fltarr(jmax,lmax)

  Sp_ionzrate_1 = fltarr(nii,jmax,nbins)   ;DI
  Sp_ionzrate_2 = fltarr(jmax,nbins)
  Sp_ionzrate_3 = fltarr(jmax)
  Sp_ionzrate_4 = fltarr(nii,jmax)
  Sp_ionzrate_5 = fltarr(jmax,lmax)



  

  ;DE

  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states



          Sp_exctrate_1[s,j,e]=Sp_exctrate_1[s,j,e]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_2[j,e]=Sp_exctrate_2[j,e]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_3[j]=Sp_exctrate_3[j]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_4[s,j]=Sp_exctrate_4[s,j]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_5[j,l]=Sp_exctrate_5[j,l]+Sp_exctrate[s,j,e,l]




        endfor
      endfor
    endfor

  endfor

  ;DI
  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nii-1 do begin  ;states



          Sp_ionzrate_1[s,j,e]=Sp_ionzrate_1[s,j,e]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_2[j,e]=Sp_ionzrate_2[j,e]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_3[j]=Sp_ionzrate_3[j]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_4[s,j]=Sp_ionzrate_4[s,j]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_5[j,l]=Sp_ionzrate_5[j,l]+Sp_ionzrate[s,j,e,l]




        endfor
      endfor
    endfor

  endfor



  

  ; Percent contribution to dissociation-New- without photon DI

  Sp_prcnttotdisionz=100*Sp_ionzrate_3/(Sp_exctrate_3+Sp_ionzrate_3)
  Sp_prcnttotdisexct=100*Sp_exctrate_3/(Sp_exctrate_3+Sp_ionzrate_3)


  Sp_prcnttotdisionz[where(~finite(Sp_prcnttotdisionz),/null)]=0.
  Sp_prcnttotdisexct[where(~finite(Sp_prcnttotdisexct),/null)]=0.




  ;      Contribution to dissociation of each state for each energy and altitude

  Sp_contrdexct=fltarr(nei,jmax,nbins)
  Sp_contrdi=fltarr(nii,jmax,nbins)
  ;  toteidiss=reform(total(Sp_ionzrate_1,1))+reform(total(Sp_exctrate_1,1))
  ;
  ;;  for s=0, nei-1 do begin
  ;for e=0, nbins-1 do begin
  ;;     Sp_contrdexct[s,*,*]= 100.*reform( Sp_exctrate_1[s,*,*])/toteidiss
  ;    Sp_contrdexct[*,*,e]= 100.*reform( Sp_exctrate_1[*,*,e])/reform(total(Sp_exctrate_1,3))
  ;    Sp_contrdi[*,*,e]= 100.*reform( Sp_ionzrate_1[*,*,e])/reform(total(Sp_ionzrate_1,3))
  ;  endfor
  ;
  ;
  ;;  for s=0, nii-1 do begin
  ;;    Sp_contrdi[s,*,*]= 100.*reform(Sp_ionzrate_1[s,*,*])/toteidiss
  ;;
  ;;  endfor


  ;>>>>>>>>>>.......>>>>>>>>>>>>>>>>>>>>>>>>>>>..
  ; 3 rd testing to calculate percentage contribution

  toteidiss=reform(total(Sp_ionzrate_1,1))+reform(total(Sp_exctrate_1,1))
  toteidiss=reform(total(toteidiss,2))

  for s=0, nei-1 do begin
    for e=0, nbins-1 do begin
      Sp_contrdexct[s,*,e]= 100.*reform( Sp_exctrate_1[s,*,e])/toteidiss

    endfor
  endfor

  for s=0, nii-1 do begin
    for e=0, nbins-1 do begin

      Sp_contrdi[s,*,e]= 100.*reform( Sp_ionzrate_1[s,*,e])/toteidiss
    endfor
  endfor



  ;>>>>>>>>>>.......>>>>>>>>>>>>>>>>>>>>>>>>>>>..

  Sp_contrdexct[where(~finite(Sp_contrdexct),/null)]=0.
  Sp_contrdi[where(~finite(Sp_contrdi),/null)]=0.

  ;_________________________________________________________________________
  ;  Contribution of each state to total dissociative excitation


  totexctdiss=reform(total(Sp_exctrate_4,1))
  St_contrdexct=fltarr(nei,jmax)

  for s=0, nei-1 do begin
    for j=0, jmax-1 do begin
      St_contrdexct[s,j]= 100.*reform(Sp_exctrate_4[s,j])/totexctdiss[j]

    endfor
  endfor


  St_contrdexct[where(~finite(St_contrdexct),/null)]=0.


end


;____________________________________________________________________________________________________________________________________________________________________


pro average_ener,Sp_exctrate,ener,avg_ener_3,avg_ener_2,avg_ener_1,Sp_rate_1
  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]



  Sp_tot_3 = fltarr(nei,jmax,lmax)
  Sp_tot_2 = fltarr(nei,jmax)
  Sp_tot_1 = fltarr(jmax)

  Sp_rate_3 = fltarr(nei,jmax,lmax)
  Sp_rate_2 = fltarr(nei,jmax)
  Sp_rate_1 = fltarr(jmax)

  avg_ener_3 = fltarr(nei,jmax,lmax)
  avg_ener_2 = fltarr(nei,jmax)
  avg_ener_1 = fltarr(jmax)




  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states

          Sp_tot_3[s,j,l]=Sp_tot_3[s,j,l]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_tot_2[s,j]=Sp_tot_2[s,j]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_tot_1[j]=Sp_tot_1[j]+Sp_exctrate[s,j,e,l]*ener[e]


          Sp_rate_3[s,j,l]=Sp_rate_3[s,j,l]+Sp_exctrate[s,j,e,l]
          Sp_rate_2[s,j]=  Sp_rate_2[s,j]+Sp_exctrate[s,j,e,l]
          Sp_rate_1[j]=    Sp_rate_1[j]+Sp_exctrate[s,j,e,l]




        endfor
      endfor
    endfor

  endfor

  avg_ener_3=Sp_tot_3/Sp_rate_3
  avg_ener_2=Sp_tot_2/Sp_rate_2
  avg_ener_1=Sp_tot_1/Sp_rate_1




end


;____________________________________________________________________________________________________________________________________________________________________


pro max_ener,Sp_exctrate,Sp_disionzrate,ener,maxener,sdener,avg_ener_3,$
  cumprcntdiss,prcntdiss,totatmpercentdiss

  ;Calculates energy of max dissociation and 95% dissociation

  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  nii= (size(Sp_disionzrate))[1]



  ; DE


  Sp_erate_1 = fltarr(nei,jmax,nbins)
  Sp_erate_2 = fltarr(jmax,nbins)
  Sp_erate_3 = fltarr(jmax)
  Sp_erate_4 = fltarr(nei,jmax)
  Sp_erate_5 = fltarr(nei,jmax,lmax)

  Sp_etot_1 = fltarr(nei,jmax,nbins)
  Sp_etot_2 = fltarr(jmax,nbins)
  Sp_etot_3 = fltarr(jmax)
  Sp_etot_4 = fltarr(nei,jmax)
  Sp_etot_5 = fltarr(nei,jmax,lmax)




  ;DI


  Sp_irate_2 = fltarr(jmax,nbins)
  Sp_irate_3 = fltarr(jmax)
  Sp_itot_2 = fltarr(jmax,nbins)
  Sp_itot_3 = fltarr(jmax)



  maxener = fltarr(jmax)
  sdener = fltarr(jmax)


  cumprcntdiss=fltarr(jmax,nbins) ;how much % dissociation we have gotten to at a particular ener bin
  prcntdiss=fltarr(jmax,nbins)
  totatmpercentdiss =fltarr(jmax,nbins)


  ;  DE
  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states



          Sp_erate_1[s,j,e]=Sp_erate_1[s,j,e]+Sp_exctrate[s,j,e,l]
          Sp_erate_2[j,e]=Sp_erate_2[j,e]+Sp_exctrate[s,j,e,l]
          Sp_erate_3[j]=Sp_erate_3[j]+Sp_exctrate[s,j,e,l]
          Sp_erate_4[s,j]=Sp_erate_4[s,j]+Sp_exctrate[s,j,e,l]

          Sp_etot_1[s,j,e]=Sp_etot_1[s,j,e]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_etot_2[j,e]=Sp_etot_2[j,e]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_etot_3[j]=Sp_etot_3[j]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_etot_4[s,j]=Sp_etot_4[s,j]+Sp_exctrate[s,j,e,l]*ener[e]


        endfor
      endfor
    endfor

  endfor


  ;  DI
  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nii-1 do begin  ;states




          Sp_irate_2[j,e]=Sp_irate_2[j,e]+Sp_disionzrate[s,j,e,l]
          Sp_irate_3[j]=Sp_irate_3[j]+Sp_disionzrate[s,j,e,l]


          Sp_itot_2[j,e]=Sp_itot_2[j,e]+Sp_disionzrate[s,j,e,l]*ener[e]
          Sp_itot_3[j]=Sp_itot_3[j]+Sp_disionzrate[s,j,e,l]*ener[e]

        endfor
      endfor
    endfor

  endfor

  ;  Total dissociation rate
  Sp_rate_2=Sp_irate_2+Sp_erate_2
  Sp_tot_2=Sp_itot_2+Sp_etot_2

  Sp_rate_3=Sp_irate_3+Sp_erate_3
  Sp_tot_3=Sp_itot_3+Sp_etot_3

  ;  Average Energy
  avg_ener_3=Sp_tot_3/Sp_rate_3
  avg_ener_3[where(~finite(avg_ener_3),/null)]=0.


  totatmpercentdiss=Sp_rate_2/total(Sp_rate_2)

  for j=0,jmax-1 do begin  ;altitude

    mi1=max(reform(Sp_rate_2[j,*]),ind1)
    maxener[j]=ener[ind1]
    ;    if j eq 75 then stop
    ;print, j,ener[ind1]
    Sp_tot=fltarr(nbins)
    Sp_tot[0]=Sp_tot_2[j,0]

    rateratio=reform(Sp_rate_2[j,*])/total(reform(Sp_rate_2[j,*]));N2_tot/reform(total(N2_tot))
    prcntdiss[j,*]=rateratio

    cumprcntdiss[j,0]=rateratio[0]

    for e=1,nbins-1 do begin

      Sp_tot[e]=total(Sp_tot_2[j,0:e],2)
      cumprcntdiss[j,e]= total(rateratio[0:e])
    endfor

    jj=where(  ( (total(Sp_tot_2[j,*])*0.68)-Sp_tot )  ge 0,/null)
    ;         print,ener[jj[-1]]
    ;         stop
    sdener[j]=ener[jj[-1]]

  endfor




end

;____________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________

const=1.e-16

sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"

fname=[sav_loc+'SQX9_NC_17JUL2022.sav',$    ;2 X9 spectra with SQ binning   
       sav_loc+'X9_NC_17JUL2022.sav']       
  



; file has:
; Sp_exctrate1_transp(nstate,nmaj,jmax,nbins,lmax),
; Sp_ionzrate1_transp(nstate,nmaj,jmax,nbins,lmax),
; Sp_exctrate1_local(nstate,nmaj,jmax,nbins,lmax),
; Sp_ionzrate1_local(nstate,nmaj,jmax,nbins,lmax),
; Sp_photoi(nmaj,jmax,lmax),
; Sp_photoki(nmaj,pstate,jmax,lmax),
; tflux_all(jmax,nbins,lmax),
; pespec_all(jmax,nbins,lmax),
; zz,ener,del,wv1,wv2,sza,lat,lon,tau,sigix,sigex,zmaj


;____________________________________________________________________________________________________________________________________________________________________



plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"
;____________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________
;
;Quiet New
restore,fname[0]
tflux_all_QN=tflux_all
pespec_all_QN=pespec_all
Sp_photoi_QNT=Sp_photoi
N2_NDI_3_QN=reform( total(reform(total(reform(Sp_ionzrate1_transp[0,2,*,*,*]),2)),2)  )
O2_NDI_3_QN=reform( total(reform(total(reform(Sp_ionzrate1_transp[0,1,*,*,*]),2)),2)  )
N2_NDI_6_QN= reform(total(reform(Sp_ionzrate1_transp[0,2,*,*,*]),2))
O2_NDI_6_QN=reform(total(reform(Sp_ionzrate1_transp[0,1,*,*,*]),2))


zmaj_QN=zmaj
;tflux_QN=reform(total(tflux_all,3))
pespec_QN=reform(total(pespec_all,3))
tflux_2_QN= reform(total(tflux_all_QN,3))



dissstates_NEW,Sp_exctrate1_transp,11,5,Sp_ionzrate1_transp,1,1,Sp_photoki,$
  N2_dissext_QNT,N2_dissionz_QNT,O2_dissext_QNT,O2_dissionz_QNT




;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_1_QNT,eiionz_2_QNT;,photoi_QNT

;N2
disscontrb,N2_dissext_QNT,N2_dissionz_QNT,N2_exctrate_1_QNT,N2_exctrate_2_QNT,N2_exctrate_3_QNT,N2_exctrate_5_QNT,$
  N2_ionzrate_1_QNT,N2_ionzrate_2_QNT,N2_ionzrate_3_QNT,N2_ionzrate_5_QNT,$
  N2_prcnttotdisionz_QNT,N2_prcnttotdisexct_QNT,$
  N2_contrdexct_QNT,N2_contrdi_QNT,N2_stcontrdexct_QNT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT
;
max_ener,N2_dissext_QNT,N2_dissionz_QNT,ener,N2_maxener_QNT,N2_sdener_QNT,N2_avg_ener_3_QNT,N2_cumprcntdiss_QNT,N2_prcntdiss_QNT,N2_totatmpercentdiss_QNT


;O2
disscontrb,O2_dissext_QNT,O2_dissionz_QNT,O2_exctrate_1_QNT,O2_exctrate_2_QNT,O2_exctrate_3_QNT,O2_exctrate_5_QNT,$
  O2_ionzrate_1_QNT,O2_ionzrate_2_QNT,O2_ionzrate_3_QNT,O2_ionzrate_5_QNT,$
  O2_prcnttotdisionz_QNT,O2_prcnttotdisexct_QNT,$
  O2_contrdexct_QNT,O2_contrdi_QNT,O2_stcontrdexct_QNT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,O2_dissext_QNT,O2_dissionz_QNT,ener,O2_maxener_QNT,O2_sdener_QNT,O2_avg_ener_3_QNT,O2_cumprcntdiss_QNT,O2_prcntdiss_QNT,O2_totatmpercentdiss_QNT
;
;O

O_dissext_QNT=reform(Sp_exctrate1_transp[0:13,0,*,*,*])
O_exctrate_QNT=reform(total(reform(total(reform(total(O_dissext_QNT,1)),2)),2))

O_dissionz_QNT=(reform(Sp_ionzrate1_transp[*,0,*,*,*]))[0,*,*,*]
;O_phtndissionz_QNT=reform( Sp_photoki(0,0:5,*,*))

disscontrb,O_dissext_QNT,O_dissionz_QNT,$
  O_exctrate_1_QNT,O_exctrate_2_QNT,O_exctrate_3_QNT,O_exctrate_5_QNT,$
  O_ionzrate_1_QNT,O_ionzrate_2_QNT,O_ionzrate_3_QNT,O_ionzrate_5_QNT,$
  O_prcnttotdisionz_QNT,O_prcnttotdisexct_QNT,$
  O_contrdexct_QNT,O_contrdi_QNT,O_stcontrdexct_QNT


max_ener,O_dissext_QNT,O_dissionz_QNT,ener,O_maxener_QNT,O_sdener_QNT,O_avg_ener_3_QNT,O_cumprcntdiss_QNT,O_prcntdiss_QNT,O_totatmpercentdiss_QNT

;_____________________________________________________________________________________

;X9 New
restore,fname[1]
tflux_all_X9N=tflux_all
pespec_all_X9N=pespec_all
Sp_photoi_X9NT=Sp_photoi
N2_NDI_3_X9N=reform( total(reform(total(reform(Sp_ionzrate1_transp[0,2,*,*,*]),2)),2)  )
O2_NDI_3_X9N=reform( total(reform(total(reform(Sp_ionzrate1_transp[0,1,*,*,*]),2)),2)  )
N2_NDI_6_X9N= reform(total(reform(Sp_ionzrate1_transp[0,2,*,*,*]),2))
O2_NDI_6_X9N=reform(total(reform(Sp_ionzrate1_transp[0,1,*,*,*]),2))

zmaj_X9N=zmaj
sigex_NEW=sigex
sigix_NEW=sigix


;tflux_X9N=reform(total(tflux_all,3))
pespec_X9N=reform(total(pespec_all,3))
tflux_2_X9N= reform(total(tflux_all_X9N,3))

;Dissociatiove excitation and ionisation cross-sections-NEW
eii=[where(ener le 10,/null)] ; for this state dissociation does not occur below 10 eV
sigex_NEW(5,2,eii)=0.

;Zero out predissociation of O2 below dissociation limit- AA'c 5.583 eV and SR 7.083 eV**

eii=[where(ener le 5.583,/null)]
sigex_NEW(1,1,eii)=0.

eii=[where(ener le 7.083,/null)]
sigex_NEW(2,1,eii)=0.



totdisscrss_N2= 0.12*sigex_NEW[5,2,*]+0.5*sigex_NEW[6,2,*]+0.95*sigex_NEW[9,2,*]+0.10*sigex_NEW[10,2,*]+0.84*sigex_NEW[11,2,*]+sigex_NEW[12,2,*]+$
  sigex_NEW[13,2,*]+sigex_NEW[14,2,*]+ 0.99*sigex_NEW[16,2,*]+sigex_NEW[17,2,*]+sigex_NEW[18,2,*]+sigix_NEW[1,2,*]
totdissexccrss_N2= 0.12*sigex_NEW[5,2,*]+0.5*sigex_NEW[6,2,*]+0.95*sigex_NEW[9,2,*]+0.10*sigex_NEW[10,2,*]+0.84*sigex_NEW[11,2,*]+sigex_NEW[12,2,*]+$
  sigex_NEW[13,2,*]+sigex_NEW[14,2,*]+ 0.99*sigex_NEW[16,2,*]+sigex_NEW[17,2,*]+sigex_NEW[18,2,*]
totdissionzcrss_N2= sigix_NEW[1,2,*]

;totdisscrss_O2= sigex_NEW[1,1,*]+sigex_NEW[2,1,*]+sigex_NEW[3,1,*]+sigex_NEW[4,1,*]+sigex_NEW[5,1,*]+sigix_NEW[1,1,*]
;totdissexccrss_O2= sigex_NEW[1,1,*]+sigex_NEW[2,1,*]+sigex_NEW[3,1,*]+sigex_NEW[4,1,*]+sigex_NEW[5,1,*]
;totdissionzcrss_O2= sigix_NEW[1,1,*]
;
;
dissstates_NEW,Sp_exctrate1_transp,11,5,Sp_ionzrate1_transp,1,1,Sp_photoki,$
  N2_dissext_X9NT,N2_dissionz_X9NT,O2_dissext_X9NT,O2_dissionz_X9NT
;
;
;
;
;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_1_X9NT,eiionz_2_X9NT;,photoi_QNT

;;N2
;
;;;Adding the b1piu and c1piu into 1 to compare with GLOW -COMMENT IT LATER **
;;N2_dissext_X9NT(2,*,*,*)=reform(N2_dissext_X9NT(2,*,*,*))+reform(N2_dissext_X9NT(8,*,*,*))
;;N2_dissext_X9NT=[N2_dissext_X9NT(0:7,*,*,*),N2_dissext_X9NT(9:10,*,*,*)]


disscontrb,N2_dissext_X9NT,N2_dissionz_X9NT,N2_exctrate_1_X9NT,N2_exctrate_2_X9NT,N2_exctrate_3_X9NT,N2_exctrate_5_X9NT,$
  N2_ionzrate_1_X9NT,N2_ionzrate_2_X9NT,N2_ionzrate_3_X9NT,N2_ionzrate_5_X9NT,$
  N2_prcnttotdisionz_X9NT,N2_prcnttotdisexct_X9NT,$
  N2_contrdexct_X9NT,N2_contrdi_X9NT,N2_stcontrdexct_X9NT


;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,N2_dissext_X9NT,N2_dissionz_X9NT,ener,N2_maxener_X9NT,N2_sdener_X9NT,N2_avg_ener_3_X9NT,N2_cumprcntdiss_X9NT,N2_prcntdiss_X9NT,N2_totatmpercentdiss_X9NT


;O2
disscontrb,O2_dissext_X9NT,O2_dissionz_X9NT,O2_exctrate_1_X9NT,O2_exctrate_2_X9NT,O2_exctrate_3_X9NT,O2_exctrate_5_X9NT,$
  O2_ionzrate_1_X9NT,O2_ionzrate_2_X9NT,O2_ionzrate_3_X9NT,O2_ionzrate_5_X9NT,$
  O2_prcnttotdisionz_X9NT,O2_prcnttotdisexct_X9NT,$
  O2_contrdexct_X9NT,O2_contrdi_X9NT,O2_stcontrdexct_X9NT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,O2_dissext_X9NT,O2_dissionz_X9NT,ener,O2_maxener_X9NT,O2_sdener_X9NT,O2_avg_ener_3_X9NT,O2_cumprcntdiss_X9NT,O2_prcntdiss_X9NT,O2_totatmpercentdiss_X9NT


;O

O_dissext_X9NT=reform(Sp_exctrate1_transp[0:13,0,*,*,*])
O_exctrate_X9NT=reform(total(reform(total(reform(total(O_dissext_X9NT,1)),2)),2))

O_dissionz_X9NT=(reform(Sp_ionzrate1_transp[*,0,*,*,*]))[0,*,*,*]
;O_phtndissionz_X9NT=reform( Sp_photoki(0,0:5,*,*))
disscontrb,O_dissext_X9NT,O_dissionz_X9NT,$
  O_exctrate_1_X9NT,O_exctrate_2_X9NT,O_exctrate_3_X9NT,O_exctrate_5_X9NT,$
  O_ionzrate_1_X9NT,O_ionzrate_2_X9NT,O_ionzrate_3_X9NT,O_ionzrate_5_X9NT,$
  O_prcnttotdisionz_X9NT,O_prcnttotdisexct_X9NT,$
  O_contrdexct_X9NT,O_contrdi_X9NT,O_stcontrdexct_X9NT


max_ener,O_dissext_X9NT,O_dissionz_X9NT,ener,O_maxener_X9NT,O_sdener_X9NT,O_avg_ener_3_X9NT,O_cumprcntdiss_X9NT,O_prcntdiss_X9NT,O_totatmpercentdiss_X9NT




;____________________________________________________________________________________________________________________________________________________________________






;Universal Plot variables
xg=0 ;gridstyles
yg=0
xtl=1.  ;ticklen
ytl=1.
xsg=1
ysg=1
xstl=1.
ystl=1.
fntsz=14
lgdsz=12
fntstyl=1 ;Bold 0-Normal
transp=100

colors=['red','dodger blue','spring green','fuchsia','yellow','deep pink','firebrick','medium blue','orange','purple','dark khaki','blue','gold']
lstyl=[0,2,3,4]
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]

postleft=[0.1,0.1,0.5,0.90]
postright=[0.55,0.1,0.95,0.90]

thck=5.
leg_loc=[0.47, 0.85]
leg_loc_left=[0.37, 0.87]
maj=['O ','O!L2!N ','N!L2!N ']










;_________________________________________________________________________________________________________________________________

;Solar flux plot

restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5
ssflux_M5=peakm5
ssflux_X9=meanpeakx9
wv1_NRL=wave_gcm1
wv2_NRL=wave_gcm2

;Quiet spectrum with NRLFLARE resolution
restore,'/home/srimoyee/Desktop/nrl_files/sav_files/SQ_newspectra.sav'
;file has: wave_gcm1,wave_gcm2,meanpeakx9
ssflux_SQ=meanpeakx9
wv1_SQ=wave_gcm1
wv2_SQ=wave_gcm2


;1 solar flux
xr=[0.05,105.]
yr=[1.e2,1.e12]
xt='Wavelength (nm)';$\AA$
yt='Flux (Photons cm!U-2!N s!U-1!N)'
titl="Input Solar Spectra"
w2a=window(window_title=titl,dimension=[450,450])

s1=plot(wv1_NRL*0.1, ssflux_X9, $
  xtitle=xt,ytitle=yt,title=titl,name='X9 Flare',$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,margin=[0.25,0.20,0.15,0.15],$
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,histogram=1,/current) ;layout=[3,1,1],



s3=plot(wv1_SQ*0.1, ssflux_SQ,name="SQ'05",$
  histogram=1,color=colors[1],thick=3,linestyle=lstyl[0],/overplot)

leg10 = LEGEND(TARGET=[s1,s3], POSITION=[0.2,0.9],font_size=lgdsz,/normal,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.15,font_style=fntstyl)

h=6.62607004e-34 ;Js
c=3.e8 ;m s-1

;NRL energy unit
NRLX9_eneunit =  (h*c*ssflux_x9*1.e14)/(0.5*(wv1_NRL+wv2_NRL)) ;W m-2
SQ_eneunit =  (h*c*ssflux_SQ*1.e14)/(0.5*(wv1_SQ+wv2_SQ)) ;W m-2

yt='Solar Flux (W m!U-2!N  )'
xt='Wavelength (nm)';($\AA$)'
xr=[0.05,107.5]
yr=[1.e-8,1.e-1]


w2b = WINDOW(DIMENSIONS=[450,450])

n1=plot(wv1_nrl*0.1,nrlx9_eneunit,name='X9 Flare',$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.25,0.20,0.15,0.15],$
  ;  ytickvalues=[1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.e0],ytickname=['1.e-5','1.e-4','1.e-3','1.e-2','1.e-1','1'],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,histogram=1,$;
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=1,/current)

n3=plot(wv1_SQ*0.1,sq_eneunit,name="SQ'05",$
  thick=3,color=colors[1],linestyle=lstyl[0],histogram=1,/overplot)

leg10 = LEGEND(TARGET=[n1,n3], POSITION=[0.2,0.9],font_size=lgdsz,/normal,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.15,font_style=fntstyl)


;_____________________________________________________________________________________________________
;----------------------------FLARE--------------------------------------------------------------------

;Total rates - Photon Ionisation, Ionisation, DI,DE, Dissoc.
;titl='Total Dissociation rate of N!L2!Nyt=' Altitude (km)'

yt=' Altitude (km)'


yr1=[50,195.]
xr1=[1.e-4,1.e5]

plotmargin=[0.2,0.2,0.2,0.2]

titl=['0.05-0.4 nm','0.4-0.8 nm','0.8-1.8 nm']

;#1
w7a = WINDOW(DIMENSIONS=[450,450])
xt='Rate (cm!U-3!N s!U-1!N)'
;x9_photoi_tot=Sp_photoi_X9NT[2,*,0]+Sp_photoi_X9NT[2,*,1]+Sp_photoi_X9NT[2,*,2]+Sp_photoi_X9NT[2,*,3]+Sp_photoi_X9NT[2,*,4]+Sp_photoi_X9NT[2,*,5]
;x9_photoi_tot=Sp_photoi_X9NT[2,*,6]+Sp_photoi_X9NT[2,*,7]+Sp_photoi_X9NT[2,*,8]
x9_photoi_tot=Sp_photoi_X9NT[2,*,9]+Sp_photoi_X9NT[2,*,10]+Sp_photoi_X9NT[2,*,11]
wvlind=2
d1a=plot(x9_photoi_tot,zz,name='Pi(X9 Flare)',$
        xtitle=xt,ytitle=yt,title='Photoionization Rate !C in'+titl[wvlind],$
        xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
        ytickvalues=[50,100,150,200],$
        thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
        margin=plotmargin,/current)

d1b=plot(Sp_photoi_QNT[2,*,wvlind],zz,name="Pi(SQ'05)",$
        thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
        
leg7a = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


w7b = WINDOW(DIMENSIONS=[450,450])
;x9_eiionz_tot=  eiionz_2_X9NT[2,*,0]+eiionz_2_X9NT[2,*,1]+eiionz_2_X9NT[2,*,2]+eiionz_2_X9NT[2,*,3]+eiionz_2_X9NT[2,*,4]+eiionz_2_X9NT[2,*,5]
;x9_eiionz_tot=  eiionz_2_X9NT[2,*,6]+eiionz_2_X9NT[2,*,7]+eiionz_2_X9NT[2,*,8]
x9_eiionz_tot=  eiionz_2_X9NT[2,*,9]+eiionz_2_X9NT[2,*,10]+eiionz_2_X9NT[2,*,11]
d1a=plot(x9_eiionz_tot,zz,name='Pe(X9 Flare)',$
        xtitle=xt,ytitle=yt,title='Photoelectron Ionization !CRate in '+titl[wvlind],$
        xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
        ytickvalues=[50,100,150,200],$
        thick=thck,color=colors[3],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
        margin=plotmargin,/current)

d1b=plot(eiionz_2_QNT[2,*,wvlind],zz,name="Pe(SQ'05)",$
  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)

leg7b = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

w7c = WINDOW(DIMENSIONS=[450,450])
;x9_ndi_tot=N2_NDI_6_X9N[*,0]+N2_NDI_6_X9N[*,1]+N2_NDI_6_X9N[*,2]+N2_NDI_6_X9N[*,3]+N2_NDI_6_X9N[*,4]+N2_NDI_6_X9N[*,5]
;x9_ndi_tot=N2_NDI_6_X9N[*,6]+N2_NDI_6_X9N[*,7]+N2_NDI_6_X9N[*,8]
x9_ndi_tot=N2_NDI_6_X9N[*,9]+N2_NDI_6_X9N[*,10]+N2_NDI_6_X9N[*,11]

d1a=plot(x9_ndi_tot,zz,name='N!L2!N!U+!N(X9 Flare)',$
         xtitle=xt,ytitle=yt,title='Photoelectron Non-Dissociative !CIonization Rate in'+titl[wvlind],$
         xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
         ytickvalues=[50,100,150,200],$
         thick=thck,color=colors[8],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
         margin=plotmargin,/current)

d1b=plot(N2_NDI_6_QN[*,wvlind],zz,name="N!L2!N!U+!N(SQ'05)",$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)

leg7c = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)




w7d = WINDOW(DIMENSIONS=[450,450])

;x9_ionzrate_5_tot=N2_ionzrate_5_X9NT[*,0]+N2_ionzrate_5_X9NT[*,1]+N2_ionzrate_5_X9NT[*,2]+N2_ionzrate_5_X9NT[*,3]+N2_ionzrate_5_X9NT[*,4]+N2_ionzrate_5_X9NT[*,5]
;x9_ionzrate_5_tot=N2_ionzrate_5_X9NT[*,6]+N2_ionzrate_5_X9NT[*,7]+N2_ionzrate_5_X9NT[*,8]
x9_ionzrate_5_tot=N2_ionzrate_5_X9NT[*,9]+N2_ionzrate_5_X9NT[*,10]+N2_ionzrate_5_X9NT[*,11]


  d1a=plot(x9_ionzrate_5_tot,zz,name='D.I(X9 Flare)',$
    xtitle=xt,ytitle=yt,title='Photoelectron Dissociative!C Ionization Rate in '+titl[wvlind],$
    xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
    ytickvalues=[50,100,150,200],$
    thick=thck,color=colors[2],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
    margin=plotmargin,/current)


d1c=plot(N2_ionzrate_5_QNT[*,wvlind],zz,name="D.I(SQ'05)",$
  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)

leg7d = LEGEND(TARGET=[d1a,d1c], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

w7e = WINDOW(DIMENSIONS=[450,450])
;x9_exctrate_5_tot=N2_exctrate_5_X9NT[*,0]+N2_exctrate_5_X9NT[*,1]+N2_exctrate_5_X9NT[*,2]+N2_exctrate_5_X9NT[*,3]+N2_exctrate_5_X9NT[*,4]+N2_exctrate_5_X9NT[*,5]
;x9_exctrate_5_tot=N2_exctrate_5_X9NT[*,6]+N2_exctrate_5_X9NT[*,7]+N2_exctrate_5_X9NT[*,8]
x9_exctrate_5_tot=N2_exctrate_5_X9NT[*,9]+N2_exctrate_5_X9NT[*,10]+N2_exctrate_5_X9NT[*,11]

d1a=plot(x9_exctrate_5_tot,zz,name='D.E(X9 Flare)',$
    xtitle=xt,ytitle=yt,title='Photoelectron Dissociative!C Excitation Rate in '+titl[wvlind],$
    xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
    ytickvalues=[50,100,150,200],$
    thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
    margin=plotmargin,/current)


d1d=plot(N2_exctrate_5_QNT[*,wvlind],zz,name="D.E(SQ'05)",$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

leg7e = LEGEND(TARGET=[d1a,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


w7f = WINDOW(DIMENSIONS=[450,450])
x9_dissoc_tot=x9_ionzrate_5_tot+x9_exctrate_5_tot
 
d1a=plot(x9_dissoc_tot,zz,name='Dissoc(X9 Flare)',$
    xtitle=xt,ytitle=yt,title='Photoelectron Dissociation!C Rate in '+titl[wvlind],$
    xstyle=1,ystyle=1,yrange=yr1,xrange=xr1,xlog=1,ylog=0,$;
    ytickvalues=[50,100,150,200],$
    thick=thck,color=colors[7],linestyle=lstyl[1],font_size=fntsz,font_style=1,$
    margin=plotmargin,/current)

d1e=plot(N2_ionzrate_5_QNT[*,wvlind]+N2_exctrate_5_QNT[*,wvlind],zz,name="Dissoc(SQ'05)",$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)

leg7f = LEGEND(TARGET=[d1a,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


;________________________________________________________________________________________________________________________________________________

;Total rates - Photon Ionisation, Ionisation, DI,DE, Dissoc.
;fntsz=14
;lgdsz=12
yt=' Altitude (km)'
thck=5

yr=[55,195]

w3a = WINDOW(DIMENSIONS=[450,450])

;N2

titl='Comparisons of Photoelectron !CImpact Rates of N!L2!N'
xt="Rate (cm!U-3!N s!U-1!N)!C (X9 Flare-Dashed,SQ'05-Bold)";

d1a=plot((reform(total(Sp_photoi_X9NT,3)))[2,*],zz,name='Pi',$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.2,0.2,0.2,0.2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=[1.e-1,1.e6],xlog=1,ylog=0,$;xrange=[1.e-4,1.e5]
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(eiionz_1_X9NT[2,*],zz,name='Pe',$
  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

d1c=plot(N2_NDI_3_X9N,zz,name='N!L2!N!U+!N',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1d=plot(N2_ionzrate_3_X9NT,zz,name= 'DI',$
  thick=thck,color=colors[2],linestyle=lstyl[1],/overplot)

d1e=plot(N2_exctrate_3_X9NT,zz,name='DE',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1f=plot(N2_ionzrate_3_X9NT+N2_exctrate_3_X9NT,zz,name='Dissoc',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

;Quiet

q1a=plot((reform(total(Sp_photoi_QNT,3)))[2,*],zz,name='Pi',$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
q1b=plot(eiionz_1_QNT[2,*],zz,name='Pe',$
  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)

q1c=plot(N2_NDI_3_QN,zz,name='N!L2!N!U+!N',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)


q1d=plot(N2_ionzrate_3_QNT,zz,name='EI DI',$
  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)

q1e=plot(N2_exctrate_3_QNT,zz,name='EI DE',$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
q1f=plot(N2_ionzrate_3_QNT+N2_exctrate_3_QNT,zz,name='EI Dissoc',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)


leg3a = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

n2_r1=eiionz_1_X9NT[2,*]/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r2=N2_NDI_3_X9N/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r3=N2_ionzrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r4=N2_exctrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r5=(N2_ionzrate_3_X9NT+N2_exctrate_3_X9NT)/(reform(total(Sp_photoi_X9NT,3)))[2,*]

n2_r6=eiionz_1_QNT[2,*]/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r7=N2_NDI_3_QN/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r8=N2_ionzrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r9=N2_exctrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r10=(N2_ionzrate_3_QNT+N2_exctrate_3_QNT)/(reform(total(Sp_photoi_QNT,3)))[2,*]


N2_r1[where(~finite(N2_r1),/null)]=0.
N2_r2[where(~finite(N2_r2),/null)]=0.
N2_r3[where(~finite(N2_r3),/null)]=0.
N2_r4[where(~finite(N2_r4),/null)]=0.
N2_r5[where(~finite(N2_r5),/null)]=0.

N2_r6[where(~finite(N2_r6),/null)]=0.
N2_r7[where(~finite(N2_r7),/null)]=0.
N2_r8[where(~finite(N2_r8),/null)]=0.
N2_r9[where(~finite(N2_r9),/null)]=0.
N2_r10[where(~finite(N2_r10),/null)]=0.

w3b = WINDOW(DIMENSIONS=[450,450])
titl='Comparisons of Photoelectron!C Impact Ratios of N!L2!N'
xt="Ratio !C (X9 Flare-Dashed,SQ'05-Bold)"
e1a=plot(n2_r1,zz,name='Pe$\slash$Pi Total',$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.2,0.2,0.2,0.2],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.01,1000],xlog=1,ylog=0,xtickvalues=[0.01,0.1,1,10,100,1000],xtickname=['0.01','0.1','1','10','10!U2!N','10!U3!N'],$;
  thick=thck,color=colors[3],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

e1b=plot(N2_r2,zz,name='Pe$\slash$Pi N!L2!N !U+!N ',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)


e1c=plot(N2_r3,zz,name='Pe$\slash$Pi DI ',$
  thick=thck,color=colors[2],linestyle=lstyl[1],/overplot)

e1d=plot(N2_r4,zz,name='Pe$\slash$Pi DE ',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

e1e=plot(N2_r5,zz,name='Pe$\slash$Pi Dissoc ',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)


f1a=plot(n2_r6,zz,name=' Pe$\slash$Pi Total',$
  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)

f1b=plot(N2_r7,zz,name='Pe$\slash$Pi N!L2!N!U+!N',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)

f1c=plot(N2_r8,zz,name='Pe$\slash$Pi DI ',$
  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)

f1d=plot(N2_r9,zz,name='Pe$\slash$Pi DE ',$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

f1e=plot(N2_r10,zz,name='Pe$\slash$Pi Dissoc ',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
;
;

leg3b = LEGEND(TARGET=[e1a,e1b,e1c,e1d,e1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)




end

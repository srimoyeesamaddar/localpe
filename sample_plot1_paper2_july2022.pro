; N2 plots for Paper2
; Comparing the flare and quiet EUVAC spectra

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

pro dissstates_NEW,Sp_exctrate,nei_N2,nei_O2,Sp_ionzrate,nii_N2,nii_O2,Sp_photoki,pii_N2,pii_O2,$
  N2_dissext,N2_dissionz,N2_phtndissionz,O2_dissext,O2_dissionz,O2_phtndissionz
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
  N2_phtndissionz=fltarr(pii_N2,jmax,lmax)

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

  N2_phtndissionz[0,*,*]=  reform(Sp_photoki[2,3,*,*])
  N2_phtndissionz[1,*,*]=  reform(Sp_photoki[2,4,*,*])
  N2_phtndissionz[2,*,*]=  reform(Sp_photoki[2,5,*,*])
  N2_phtndissionz[3,*,*]=  reform(Sp_photoki[2,6,*,*])
  N2_phtndissionz[4,*,*]=  reform(Sp_photoki[2,7,*,*])


  ;O2 states

  O2_dissext=fltarr(nei_O2,jmax,nbins,lmax)
  O2_dissionz=fltarr(nii_O2,jmax,nbins,lmax)
  O2_phtndissionz=fltarr(pii_O2,jmax,lmax)

  O2_dissext[0,*,*,*]=  reform(Sp_exctrate[1,1,*,*,*])
  O2_dissext[1,*,*,*]=  reform(Sp_exctrate[2,1,*,*,*])
  O2_dissext[2,*,*,*]=  reform(Sp_exctrate[3,1,*,*,*])
  O2_dissext[3,*,*,*]=  reform(Sp_exctrate[4,1,*,*,*])
  O2_dissext[4,*,*,*]=  reform(Sp_exctrate[5,1,*,*,*])



  O2_dissionz[0,*,*,*]= reform(Sp_ionzrate[1,1,*,*,*])

  O2_phtndissionz[0,*,*]=  reform(Sp_photoki[1,3,*,*])
  O2_phtndissionz[1,*,*]=  reform(Sp_photoki[1,4,*,*])
  O2_phtndissionz[2,*,*]=  reform(Sp_photoki[1,5,*,*])
  O2_phtndissionz[3,*,*]=  reform(Sp_photoki[1,6,*,*])
  O2_phtndissionz[4,*,*]=  reform(Sp_photoki[1,7,*,*])



end

;____________________________________________________________________________________________________________________________________________________________________



pro dissstates_OLD,Sp_exctrate,nei_N2,nei_O2,Sp_ionzrate,nii_N2,nii_O2,Sp_photoki,pii_N2,pii_O2,$
  N2_dissext,N2_dissionz,N2_phtndissionz,O2_dissext,O2_dissionz,O2_phtndissionz
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
  N2_phtndissionz=fltarr(pii_N2,jmax,lmax)

  N2_dissext[0,*,*,*]=  reform(Sp_exctrate[4,2,*,*,*])
  N2_dissext[1,*,*,*]=  reform(Sp_exctrate[5,2,*,*,*])
  N2_dissext[2,*,*,*]=  reform(Sp_exctrate[6,2,*,*,*])

  N2_dissionz[0,*,*,*]= reform(Sp_ionzrate[3,2,*,*,*])
  N2_dissionz[1,*,*,*]= reform(Sp_ionzrate[4,2,*,*,*])
  N2_dissionz[2,*,*,*]= reform(Sp_ionzrate[5,2,*,*,*])

  N2_phtndissionz[0,*,*]=  reform(Sp_photoki[2,3,*,*])
  N2_phtndissionz[1,*,*]=  reform(Sp_photoki[2,4,*,*])
  N2_phtndissionz[2,*,*]=  reform(Sp_photoki[2,5,*,*])
  N2_phtndissionz[3,*,*]=  reform(Sp_photoki[2,6,*,*])
  N2_phtndissionz[4,*,*]=  reform(Sp_photoki[2,7,*,*])


  ;O2 states

  O2_dissext=fltarr(nei_O2,jmax,nbins,lmax)
  O2_dissionz=fltarr(nii_O2,jmax,nbins,lmax)
  O2_phtndissionz=fltarr(pii_O2,jmax,lmax)

  O2_dissext[0,*,*,*]=  reform(Sp_exctrate[2,1,*,*,*])
  O2_dissext[1,*,*,*]=  reform(Sp_exctrate[3,1,*,*,*])
  O2_dissext[2,*,*,*]=  reform(Sp_exctrate[4,1,*,*,*])
  O2_dissext[3,*,*,*]=  reform(Sp_exctrate[5,1,*,*,*])


  O2_dissionz[0,*,*,*]= reform(Sp_ionzrate[4,1,*,*,*])
  O2_dissionz[1,*,*,*]= reform(Sp_ionzrate[5,1,*,*,*])
  O2_dissionz[2,*,*,*]= reform(Sp_ionzrate[6,1,*,*,*])

  O2_phtndissionz[0,*,*]=  reform(Sp_photoki[1,3,*,*])
  O2_phtndissionz[1,*,*]=  reform(Sp_photoki[1,4,*,*])
  O2_phtndissionz[2,*,*]=  reform(Sp_photoki[1,5,*,*])
  O2_phtndissionz[3,*,*]=  reform(Sp_photoki[1,6,*,*])
  O2_phtndissionz[4,*,*]=  reform(Sp_photoki[1,7,*,*])



end



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


pro disscontrb,Sp_exctrate,Sp_ionzrate,Sp_phtndissionz,$
  Sp_exctrate_1,Sp_exctrate_2,Sp_exctrate_3,Sp_exctrate_5,$
  Sp_ionzrate_1,Sp_ionzrate_2,Sp_ionzrate_3,Sp_ionzrate_5,$
  Sp_phtndirate_3,Sp_phtndirate_4,$
  Sp_prcnttotdisionz,Sp_prcnttotdisexct,Sp_prcnttotpdi,$
  Sp_contrdexct,Sp_contrdi

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



  Sp_phtndirate_3=fltarr(pii,jmax,lmax)      ;Photon DI
  Sp_phtndirate_4=fltarr(pii,jmax)


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



  ;Photon DI

  for l=0,lmax-1 do begin  ;wavelength bins


    for j=0,jmax-1 do begin  ;altitude

      for s=0,pii-1 do begin  ;states




        Sp_phtndirate_3[j]=Sp_phtndirate_3[j]+Sp_phtndissionz[s,j,l]
        Sp_phtndirate_4[s,j]=Sp_phtndirate_4[s,j]+Sp_phtndissionz[s,j,l]



      endfor
    endfor


  endfor


  ; Percent contribution to dissociation-New

  Sp_prcnttotdisionz=100*Sp_ionzrate_3/(Sp_exctrate_3+Sp_ionzrate_3+Sp_phtndirate_3)
  Sp_prcnttotdisexct=100*Sp_exctrate_3/(Sp_exctrate_3+Sp_ionzrate_3+Sp_phtndirate_3)
  Sp_prcnttotpdi=100*Sp_phtndirate_3/(Sp_exctrate_3+Sp_ionzrate_3+Sp_phtndirate_3)

  Sp_prcnttotdisionz[where(~finite(Sp_prcnttotdisionz),/null)]=0.
  Sp_prcnttotdisexct[where(~finite(Sp_prcnttotdisexct),/null)]=0.
  Sp_prcnttotpdi[where(~finite(Sp_prcnttotpdi),/null)]=0.

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
  avg_ener_2=Sp_tot_2/Sp_rate_2021
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

    Sp_tot=fltarr(nbins)
    Sp_tot[0]=Sp_tot_2[j,0]

    rateratio=reform(Sp_rate_2[j,*])/total(reform(Sp_rate_2[j,*]));N2_tot/reform(total(N2_tot))
    prcntdiss[j,*]=rateratio

    cumprcntdiss[j,0]=rateratio[0]

    for e=1,nbins-1 do begin

      Sp_tot[e]=total(Sp_tot_2[j,0:e],2)
      cumprcntdiss[j,e]= total(rateratio[0:e])
    endfor

    jj=where(  ( (total(Sp_tot_2[j,*])*0.95)-Sp_tot )  ge 0,/null)
    ;         print,ener[jj[-1]]
    ;         stop
    sdener[j]=ener[jj[-1]]

  endfor




end

;____________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________

const=1.e-16

sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"

fname=[sav_loc+'QUIET_NC_17JUL2022.sav',$    ;0 Quiet, New
  sav_loc+'X9_NC_17JUL2022.sav',$       ;1 X9, New
  sav_loc+'M5_NC_19NOV2021.sav',$       ;2 M5, New
  sav_loc+'QUIET_OC_19NOV2021.sav',$    ;3 Quiet, Old
  sav_loc+'X9_OC_19NOV2021.sav',$       ;4 X9, Old
  sav_loc+'M5_OC_19NOV2021.sav']        ;5 M5, Old




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
zmaj_QN=zmaj
;tflux_QN=reform(total(tflux_all,3))
pespec_QN=reform(total(pespec_all,3))
tflux_2_QN= reform(total(tflux_all_QN,3))


xr=[9., 1.e4]
dissstates_NEW,Sp_exctrate1_transp,11,5,Sp_ionzrate1_transp,1,1,Sp_photoki,5,5,$
  N2_dissext_QNT,N2_dissionz_QNT,N2_phtndissionz_QNT,O2_dissext_QNT,O2_dissionz_QNT,O2_phtndissionz_QNT




;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_1_QNT,eiionz_2_QNT;,photoi_QNT

;N2
disscontrb,N2_dissext_QNT,N2_dissionz_QNT,N2_phtndissionz_QNT,N2_exctrate_1_QNT,N2_exctrate_2_QNT,N2_exctrate_3_QNT,N2_exctrate_5_QNT,$
  N2_ionzrate_1_QNT,N2_ionzrate_2_QNT,N2_ionzrate_3_QNT,N2_ionzrate_5_QNT,N2_phtndirate_3_QNT,N2_phtndirate_4_QNT,$
  N2_prcnttotdisionz_QNT,N2_prcnttotdisexct_QNT,N2_prcnttotpdi_QNT,$
  N2_contrdexct_QNT,N2_contrdi_QNT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,N2_dissext_QNT,N2_dissionz_QNT,ener,N2_maxener_QNT,N2_sdener_QNT,N2_avg_ener_3_QNT,N2_cumprcntdiss_QNT,N2_prcntdiss_QNT,N2_totatmpercentdiss_QNT


;O2
disscontrb,O2_dissext_QNT,O2_dissionz_QNT,O2_phtndissionz_QNT,O2_exctrate_1_QNT,O2_exctrate_2_QNT,O2_exctrate_3_QNT,O2_exctrate_5_QNT,$
  O2_ionzrate_1_QNT,O2_ionzrate_2_QNT,O2_ionzrate_3_QNT,O2_ionzrate_5_QNT,O2_phtndirate_3_QNT,O2_phtndirate_4_QNT,$
  O2_prcnttotdisionz_QNT,O2_prcnttotdisexct_QNT,O2_prcnttotpdi_QNT,$
  O2_contrdexct_QNT,O2_contrdi_QNT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,O2_dissext_QNT,O2_dissionz_QNT,ener,O2_maxener_QNT,O2_sdener_QNT,O2_avg_ener_3_QNT,O2_cumprcntdiss_QNT,O2_prcntdiss_QNT,O2_totatmpercentdiss_QNT

;O

O_dissext_QNT=reform(Sp_exctrate1_transp[0:13,0,*,*,*])
O_exctrate_QNT=reform(total(reform(total(reform(total(O_dissext_QNT,1)),2)),2))

O_dissionz_QNT=(reform(Sp_ionzrate1_transp[*,0,*,*,*]))[0,*,*,*]
O_phtndissionz_QNT=reform( Sp_photoki(0,0:5,*,*))

disscontrb,O_dissext_QNT,O_dissionz_QNT,O_phtndissionz_QNT,$
  O_exctrate_1_QNT,O_exctrate_2_QNT,O_exctrate_3_QNT,O_exctrate_5_QNT,$
  O_ionzrate_1_QNT,O_ionzrate_2_QNT,O_ionzrate_3_QNT,O_ionzrate_5_QNT,$
  O_phtndirate_3_QNT,O_phtndirate_4_QNT,$
  O_prcnttotdisionz_QNT,O_prcnttotdisexct_QNT,O_prcnttotpdi_QNT,$
  O_contrdexct_QNT,O_contrdi_QNT


max_ener,O_dissext_QNT,O_dissionz_QNT,ener,O_maxener_QNT,O_sdener_QNT,O_avg_ener_3_QNT,O_cumprcntdiss_QNT,O_prcntdiss_QNT,O_totatmpercentdiss_QNT

;_____________________________________________________________________________________

;X9 New
restore,fname[1]
tflux_all_X9N=tflux_all
pespec_all_X9N=pespec_all
Sp_photoi_X9NT=Sp_photoi
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

totdisscrss_O2= sigex_NEW[1,1,*]+sigex_NEW[2,1,*]+sigex_NEW[3,1,*]+sigex_NEW[4,1,*]+sigex_NEW[5,1,*]+sigix_NEW[1,1,*]
totdissexccrss_O2= sigex_NEW[1,1,*]+sigex_NEW[2,1,*]+sigex_NEW[3,1,*]+sigex_NEW[4,1,*]+sigex_NEW[5,1,*]
totdissionzcrss_O2= sigix_NEW[1,1,*]


dissstates_NEW,Sp_exctrate1_transp,11,5,Sp_ionzrate1_transp,1,1,Sp_photoki,5,5,$
  N2_dissext_X9NT,N2_dissionz_X9NT,N2_phtndissionz_X9NT,O2_dissext_X9NT,O2_dissionz_X9NT,O2_phtndissionz_X9NT




;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_1_X9NT,eiionz_2_X9NT;,photoi_QNT

;N2

;;Adding the b1piu and c1piu into 1 to compare with GLOW -COMMENT IT LATER **
;N2_dissext_X9NT(2,*,*,*)=reform(N2_dissext_X9NT(2,*,*,*))+reform(N2_dissext_X9NT(8,*,*,*))
;N2_dissext_X9NT=[N2_dissext_X9NT(0:7,*,*,*),N2_dissext_X9NT(9:10,*,*,*)]


disscontrb,N2_dissext_X9NT,N2_dissionz_X9NT,N2_phtndissionz_X9NT,N2_exctrate_1_X9NT,N2_exctrate_2_X9NT,N2_exctrate_3_X9NT,N2_exctrate_5_X9NT,$
  N2_ionzrate_1_X9NT,N2_ionzrate_2_X9NT,N2_ionzrate_3_X9NT,N2_ionzrate_5_X9NT,N2_phtndirate_3_X9NT,N2_phtndirate_4_X9NT,$
  N2_prcnttotdisionz_X9NT,N2_prcnttotdisexct_X9NT,N2_prcnttotpdi_X9NT,$
  N2_contrdexct_X9NT,N2_contrdi_X9NT


;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,N2_dissext_X9NT,N2_dissionz_X9NT,ener,N2_maxener_X9NT,N2_sdener_X9NT,N2_avg_ener_3_X9NT,N2_cumprcntdiss_X9NT,N2_prcntdiss_X9NT,N2_totatmpercentdiss_X9NT


;O2
disscontrb,O2_dissext_X9NT,O2_dissionz_X9NT,O2_phtndissionz_X9NT,O2_exctrate_1_X9NT,O2_exctrate_2_X9NT,O2_exctrate_3_X9NT,O2_exctrate_5_X9NT,$
  O2_ionzrate_1_X9NT,O2_ionzrate_2_X9NT,O2_ionzrate_3_X9NT,O2_ionzrate_5_X9NT,O2_phtndirate_3_X9NT,O2_phtndirate_4_X9NT,$
  O2_prcnttotdisionz_X9NT,O2_prcnttotdisexct_X9NT,O2_prcnttotpdi_X9NT,$
  O2_contrdexct_X9NT,O2_contrdi_X9NT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,O2_dissext_X9NT,O2_dissionz_X9NT,ener,O2_maxener_X9NT,O2_sdener_X9NT,O2_avg_ener_3_X9NT,O2_cumprcntdiss_X9NT,O2_prcntdiss_X9NT,O2_totatmpercentdiss_X9NT


;O

O_dissext_X9NT=reform(Sp_exctrate1_transp[0:13,0,*,*,*])
O_exctrate_X9NT=reform(total(reform(total(reform(total(O_dissext_X9NT,1)),2)),2))

O_dissionz_X9NT=(reform(Sp_ionzrate1_transp[*,0,*,*,*]))[0,*,*,*]
O_phtndissionz_X9NT=reform( Sp_photoki(0,0:5,*,*))
disscontrb,O_dissext_X9NT,O_dissionz_X9NT,O_phtndissionz_X9NT,$
  O_exctrate_1_X9NT,O_exctrate_2_X9NT,O_exctrate_3_X9NT,O_exctrate_5_X9NT,$
  O_ionzrate_1_X9NT,O_ionzrate_2_X9NT,O_ionzrate_3_X9NT,O_ionzrate_5_X9NT,$
  O_phtndirate_3_X9NT,O_phtndirate_4_X9NT,$
  O_prcnttotdisionz_X9NT,O_prcnttotdisexct_X9NT,O_prcnttotpdi_X9NT,$
  O_contrdexct_X9NT,O_contrdi_X9NT


max_ener,O_dissext_X9NT,O_dissionz_X9NT,ener,O_maxener_X9NT,O_sdener_X9NT,O_avg_ener_3_X9NT,O_cumprcntdiss_X9NT,O_prcntdiss_X9NT,O_totatmpercentdiss_X9NT




;_____________________________________________________________________________________________________________

;X9 Old cross-section
restore,fname[4]
tflux_all_X9O=tflux_all
pespec_all_X9O=pespec_all
Sp_photoi_X9OT=Sp_photoi
sigex_OLD=sigex
sigix_OLD=sigix
zmaj_X9O=zmaj

tflux_X9O=reform(total(tflux_all,3))
pespec_X9O=reform(total(pespec_all,3))
tflux_2_X9O= reform(total(tflux_all_X9N,3))

;Zero out predissociation cross-section of O2 below dissociation limit- AA'c 5.583 eV **

eii=[where(ener le 5.583,/null)]
sigex_OLD(2,1,eii)=0.

;Zero out predissociation rates of O2 below dissociation limit- AA'c 5.583 eV *

eii=where(ener le 5.583,/null)
Sp_exctrate1_transp(2,1,*,eii,*)=0.



;Dissociatiove excitation and ionisation cross-sections-OLD

totdisscrss_N2_OLD= sigex_OLD[4,2,*]+sigex_OLD[5,2,*]+sigex_OLD[6,2,*]+$
  sigix_OLD[3,2,*]+sigix_OLD[4,2,*]+sigix_OLD[5,2,*]
totdissexccrss_N2_OLD= sigex_OLD[4,2,*]+sigex_OLD[5,2,*]+sigex_OLD[6,2,*]
totdissionzcrss_N2_OLD= sigix_OLD[3,2,*]+sigix_OLD[4,2,*]+sigix_OLD[5,2,*]



totdisscrss_O2_OLD= sigex_OLD[2,1,*]+sigex_OLD[3,1,*]+sigex_OLD[4,1,*]+sigex_OLD[5,1,*]+$
  sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+sigix_OLD[6,1,*]
totdissexccrss_O2_OLD= sigex_OLD[2,1,*]+sigex_OLD[3,1,*]+sigex_OLD[4,1,*]+sigex_OLD[5,1,*]
totdissionzcrss_O2_OLD=sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+sigix_OLD[6,1,*]


dissstates_OLD,Sp_exctrate1_transp,3,4,Sp_ionzrate1_transp,3,3,Sp_photoki,5,5,$
  N2_dissext_X9OT,N2_dissionz_X9OT,N2_phtndissionz_X9OT,O2_dissext_X9OT,O2_dissionz_X9OT,O2_phtndissionz_X9OT


;Total ionisation rates

totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_1_X9OT,eiionz_2_X9OT;,photoi_QNT

;N2
disscontrb,N2_dissext_X9OT,N2_dissionz_X9OT,N2_phtndissionz_X9OT,N2_exctrate_1_X9OT,N2_exctrate_2_X9OT,N2_exctrate_3_X9OT,N2_exctrate_5_X9OT,$
  N2_ionzrate_1_X9OT,N2_ionzrate_2_X9OT,N2_ionzrate_3_X9OT,N2_ionzrate_5_X9OT,N2_phtndirate_3_X9OT,N2_phtndirate_4_X9OT,$
  N2_prcnttotdisionz_X9OT,N2_prcnttotdisexct_X9OT,N2_prcnttotpdi_X9OT,$
  N2_contrdexct_X9OT,N2_contrdi_X9OT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,N2_dissext_X9OT,N2_dissionz_X9OT,ener,N2_maxener_X9OT,N2_sdener_X9OT,N2_avg_ener_3_X9OT,N2_cumprcntdiss_X9OT,N2_prcntdiss_X9OT,N2_totatmpercentdiss_X9OT


;O2
disscontrb,O2_dissext_X9OT,O2_dissionz_X9OT,O2_phtndissionz_X9OT,O2_exctrate_1_X9OT,O2_exctrate_2_X9OT,O2_exctrate_3_X9OT,O2_exctrate_5_X9OT,$
  O2_ionzrate_1_X9OT,O2_ionzrate_2_X9OT,O2_ionzrate_3_X9OT,O2_ionzrate_5_X9OT,O2_phtndirate_3_X9OT,O2_phtndirate_4_X9OT,$
  O2_prcnttotdisionz_X9OT,O2_prcnttotdisexct_X9OT,O2_prcnttotpdi_X9OT,$
  O2_contrdexct_X9OT,O2_contrdi_X9OT



;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT

max_ener,O2_dissext_X9OT,O2_dissionz_X9OT,ener,O2_maxener_X9OT,O2_sdener_X9OT,O2_avg_ener_3_X9OT,O2_cumprcntdiss_X9OT,O2_prcntdiss_X9OT,O2_totatmpercentdiss_X9OT



;O

O_dissext_X9OT=reform(Sp_exctrate1_transp[0:7,0,*,*,*])
O_dissionz_X9OT=reform(Sp_ionzrate1_transp[0:2,0,*,*,*])
O_phtndissionz_X9OT=reform( Sp_photoki(0,0:5,*,*))
disscontrb,O_dissext_X9OT,O_dissionz_X9OT,O_phtndissionz_X9OT,$
  O_exctrate_1_X9OT,O_exctrate_2_X9OT,O_exctrate_3_X9OT,O_exctrate_5_X9OT,$
  O_ionzrate_1_X9OT,O_ionzrate_2_X9OT,O_ionzrate_3_X9OT,O_ionzrate_5_X9OT,$
  O_phtndirate_3_X9OT,O_phtndirate_4_X9OT,$
  O_prcnttotdisionz_X9OT,O_prcnttotdisexct_X9OT,O_prcnttotpdi_X9OT,$
  O_contrdexct_X9OT,O_contrdi_X9OT


max_ener,O_dissext_X9OT,O_dissionz_X9OT,ener,O_maxener_X9OT,O_sdener_X9OT,O_avg_ener_3_X9OT,O_cumprcntdiss_X9OT,O_prcntdiss_X9OT,O_totatmpercentdiss_X9OT


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









;_____________________________________________________________________________________________________
;_____________________________________________________________________________________________________
;Percentage contribution to dissociation as function of altitude

xt='Percentage contribution'
yt='Altitude (km)'

xr=[0.,100.]
yr=[55.,195]
;fntsz=14
;lgdsz=12
thck=5
;
leg_nam1=['EIE- X9 Flare','EII-X9 Flare','PI-X9 Flare',$
  'EIE-Quiet','EII-Quiet','PI-Quiet',$
  'EIE- M5 Flare','EII-M5 Flare','PI-M5 Flare']


;X9 flare

titl=[];'Percentage contribution to !C dissociation of  N!L2!N ';(Transport)
w20a = WINDOW(DIMENSIONS=[450,450])


f1=plot(N2_prcnttotdisexct_X9NT,zz,name='EI Excitation',$;leg_nam1[0]
  xtitle=xt,ytitle=yt,title=titl,margin=[0.25,0.2,0.25,0.2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

f2=plot(N2_prcnttotdisionz_X9NT,zz,name='EI Ionization',$;leg_nam1[1]
  thick=thck,color='dark blue',linestyle=lstyl[1],/overplot)
  ;leg20b = LEGEND(TARGET=[f4,f5], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)



leg20a = LEGEND(TARGET=[f1,f2], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.15,font_style=fntstyl)



;____________________________________________________________________________________________________________________________________________________________________




















;____________________________________________________________________________________________________________________________________________________________________
;Dissociation cross-sections plots

ener_ZM=[10.0,10.5, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
crss_ZM=[1.27, 1.87, 2.49, 4.92, 11.0, 23.5, 33.6, 44.1, 62.3, 77.4, 94.6, 121.0, 139.0, 150.0, 159.0, 175.0, 186.0, 194.0, 197.0, 198.0, 196.0]*1.e-18

yt='Cross-section (10!U-16!N cm!U2!N)'
xt='Energy (eV)'
xr=[9., 3.e4]
yr=[0.0,3.]

leg_nam=['Present',"SQ'05"]


;N2
w1a = WINDOW(DIMENSIONS=[450,450])
;
;fntsz=10
;lgdsz=10
transp=100
titl=[];'Dissociative Excitation of N!L2!N '
c1a= plot(ener,totdissexccrss_N2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.2,0.2,0.2,0.2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

c1b=plot(ener,totdissexccrss_N2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg1a = LEGEND(TARGET=[c1a,c1b], POSITION=[0.27,0.9],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

w1b = WINDOW(DIMENSIONS=[450,450])
titl=[];'Dissociative Ionization of N!L2!N '
yr=[0.,1.0]
c2a= plot(ener,totdissionzcrss_N2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.2,0.2,0.2,0.2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
 thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

c2b=plot(ener,totdissionzcrss_N2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg2a = LEGEND(TARGET=[c2a,c2b],  POSITION=[0.6,0.9],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

w1c = WINDOW(DIMENSIONS=[450,450])
titl=[];'Total Dissociation of N!L2!N '
yr=[0.0,3.]
c3a= plot(ener,totdisscrss_N2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.2,0.2,0.2,0.2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

c3b=plot(ener,totdisscrss_N2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

c3c=plot(ener_ZM,crss_ZM/const,name='ZM 1978',$
  thick=thck,color='orange',linestyle=lstyl[0],/overplot)


leg3a = LEGEND(TARGET=[c3a,c3b,c3c], POSITION=[0.95,0.9],/normal, font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

;O2
xr=[1., 3.e4]
yr=[0.0,3.]

titl='Dissociative Excitation of O!L2!N '
w1d = WINDOW(DIMENSIONS=[450,450])
d1a= plot(ener,totdissexccrss_O2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  layout=[3,2,4],thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(ener,totdissexccrss_O2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg4a = LEGEND(TARGET=[d1a,d1b], POSITION=[0.27,0.5],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.05,font_style=fntstyl); POSITION=leg_loc,/normal

titl='Dissociative Ionization of O!L2!N '
yr=[0.0,1.5]
d2a= plot(ener,totdissionzcrss_O2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  layout=[3,2,5],thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

d2b=plot(ener,totdissionzcrss_O2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg5a = LEGEND(TARGET=[d2a,d2b], POSITION=[0.6,0.5],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.05,font_style=fntstyl)

titl='Total Dissociation of O!L2!N '
yr=[0.0,3.]
d3a= plot(ener,totdisscrss_O2/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  layout=[3,2,6],thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

d3b=plot(ener,totdisscrss_O2_OLD/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg6a = LEGEND(TARGET=[d3a,d3b] ,POSITION=[0.95,0.5],/normal, font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.05,font_style=fntstyl)


;O
w1d = WINDOW(DIMENSIONS=[800,400])
xr=[1., 3.e4]
yr=[0.0,1.0]

titl='Electron Impact Excitaion of O'

c1= plot(ener,total(reform(sigex_NEW[*,0,*]),1)/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  layout=[2,1,1],thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

c2=plot(ener,total(reform(sigex_OLD[*,0,*]),1)/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg1d = LEGEND(TARGET=[c1,c2], POSITION=[0.9,0.85],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.05,font_style=fntstyl)

titl='Electron Impact Ionization of O'
yr=[0.0,1.5]
d1= plot(ener,total(reform(sigix_NEW[*,0,*]),1)/const,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  layout=[2,1,2], thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,font_style=fntstyl,/current)

d2=plot(ener,total(reform(sigix_OLD[*,0,*]),1)/const,name=leg_nam[1],$
  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)

leg2d = LEGEND(TARGET=[d1,d2], POSITION=[0.9,0.35],/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.05,font_style=fntstyl)

;____________________________________________________________________________________________________________________________________________________________________













































;____________________________________________________________________________________________________________________________________________________________________
;
;Total rates - Photon Ionisation, Ionisation, DI,DE, Dissoc.
;fntsz=12
;lgdsz=12
yt=' Altitude (km)'
thck=5

yr=[55,195]



w6a = WINDOW(DIMENSIONS=[450,450])
mrgn=[0.2,0.2,0.2,0.2]

;N2

titl=[];'N!L2!N
xt='Rate (cm!U-3!N s!U-1!N)!C (Flare-Dashed,Quiet-Bold)';

d1a=plot((reform(total(Sp_photoi_X9NT,3)))[2,*],zz,name='Photoionization (Pi)',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$ ;,layout=[2,3,1]
  xstyle=1,ystyle=1,yrange=yr,xrange=[1.e-1,1.e5],xlog=1,ylog=0,$;xrange=[1.e-4,1.e5]
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(eiionz_1_X9NT[2,*],zz,name='EI Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)


d1e=plot(N2_ionzrate_3_X9NT+N2_exctrate_3_X9NT,zz,name='EI Dissoc',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

;Quiet

q1a=plot((reform(total(Sp_photoi_QNT,3)))[2,*],zz,name='Photoionization',$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
q1b=plot(eiionz_1_QNT[2,*],zz,name='Pe',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)

q1e=plot(N2_ionzrate_3_QNT+N2_exctrate_3_QNT,zz,name='EI Dissoc',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)


leg6a = LEGEND(TARGET=[d1a,d1b,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

n2_r1=eiionz_1_X9NT[2,*]/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r2=N2_ionzrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r3=N2_exctrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[2,*]
n2_r4=(N2_ionzrate_3_X9NT+N2_exctrate_3_X9NT)/(reform(total(Sp_photoi_X9NT,3)))[2,*]

n2_r5=eiionz_1_QNT[2,*]/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r6=N2_ionzrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r7=N2_exctrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[2,*]
n2_r8=(N2_ionzrate_3_QNT+N2_exctrate_3_QNT)/(reform(total(Sp_photoi_QNT,3)))[2,*]


N2_r1[where(~finite(N2_r1),/null)]=0.
N2_r2[where(~finite(N2_r2),/null)]=0.
N2_r3[where(~finite(N2_r3),/null)]=0.
N2_r4[where(~finite(N2_r4),/null)]=0.

N2_r5[where(~finite(N2_r5),/null)]=0.
N2_r6[where(~finite(N2_r6),/null)]=0.
N2_r7[where(~finite(N2_r7),/null)]=0.
N2_r8[where(~finite(N2_r8),/null)]=0.

w6aa = WINDOW(DIMENSIONS=[450,450])
xt='Ratio !C (Flare-Dashed,Quiet-Bold)'
e1a=plot(n2_r1,zz,name='Pe$\slash$Pi',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$;layout=[2,3,2],
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.01,1000],xlog=1,ylog=0,xtickvalues=[0.01,0.1,1,10,100,1000],xtickname=['0.01','0.1','1','10','10!U2!N','10!U3!N'],$;
  thick=thck,color=colors[8],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)


e1d=plot(N2_r4,zz,name='EI Dissoc$\slash$Pi ',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)


f1a=plot(n2_r5,zz,name=' Pe$\slash$Pi ',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
f1d=plot(N2_r8,zz,name='EI Dissoc$\slash$Pi ',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
leg6c = LEGEND(TARGET=[e1a,e1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

;____________________________________________________________________________________________________
;
;O2
titl=[];'O!L2!N'

w6b = WINDOW(DIMENSIONS=[450,450])
xt='Rate (cm!U-3!N s!U-1!N) !C (Flare-Dashed,Quiet-Bold) '

d1a=plot((reform(total(Sp_photoi_X9NT,3)))[1,*],zz,name='Photoionization (Pi)',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$;layout=[2,3,3],
  xstyle=1,ystyle=1,yrange=yr,xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(eiionz_1_X9NT[1,*],zz,name='EI Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)


d1e=plot(O2_ionzrate_3_X9NT+O2_exctrate_3_X9NT,zz,name='EI Dissoc',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

;Quiet

q1a=plot((reform(total(Sp_photoi_QNT,3)))[1,*],zz,name='Photoionization',$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

q1b=plot(eiionz_1_QNT[1,*],zz,name='Pe',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)


q1e=plot(O2_ionzrate_3_QNT+O2_exctrate_3_QNT,zz,name='EI Dissoc',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)


leg6a = LEGEND(TARGET=[d1a,d1b,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

o2_r1=eiionz_1_X9NT[1,*]/(reform(total(Sp_photoi_X9NT,3)))[1,*]
o2_r2=O2_ionzrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[1,*]
o2_r3=O2_exctrate_3_X9NT/(reform(total(Sp_photoi_X9NT,3)))[1,*]
o2_r4=(O2_ionzrate_3_X9NT+O2_exctrate_3_X9NT)/(reform(total(Sp_photoi_X9NT,3)))[1,*]

o2_r5=eiionz_1_QNT[1,*]/(reform(total(Sp_photoi_QNT,3)))[1,*]
o2_r6=O2_ionzrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[1,*]
o2_r7=O2_exctrate_3_QNT/(reform(total(Sp_photoi_QNT,3)))[1,*]
o2_r8=(O2_ionzrate_3_QNT+O2_exctrate_3_QNT)/(reform(total(Sp_photoi_QNT,3)))[1,*]

O2_r1[where(~finite(O2_r1),/null)]=0.
O2_r2[where(~finite(O2_r2),/null)]=0.
O2_r3[where(~finite(O2_r3),/null)]=0.
O2_r4[where(~finite(O2_r4),/null)]=0.

O2_r5[where(~finite(O2_r5),/null)]=0.
O2_r6[where(~finite(O2_r6),/null)]=0.
O2_r7[where(~finite(O2_r7),/null)]=0.
O2_r8[where(~finite(O2_r8),/null)]=0.

w6bb = WINDOW(DIMENSIONS=[450,450])
xt='Ratio !C (Flare-Dashed,Quiet-Bold)'
e1a=plot(o2_r1,zz,name='Pe$\slash$Pi',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$;layout=[2,3,4]
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.01,1000],xlog=1,ylog=0,xtickvalues=[0.01,0.1,1,10,100,1000],xtickname=['0.01','0.1','1','10','10!U2!N','10!U3!N'],$;
  thick=thck,color=colors[8],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)



e1d=plot(O2_r4,zz,name='EI Dissoc$\slash$Pi ',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)


f1a=plot(o2_r5,zz,name=' Pe$\slash$Pi ',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)




f1d=plot(O2_r8,zz,name='EI Dissoc$\slash$Pi ',$
  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)



leg6c = LEGEND(TARGET=[e1a,e1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


;____________________________________________________________________________________________________
;
;O
titl=[];'O'
yt=' Altitude (km)'
xt='Rate (cm!U-3!N s!U-1!N) !C (Flare-Dashed,Quiet-Bold) '

w6c = WINDOW(DIMENSIONS=[450,450])

d1a=plot((reform(total(Sp_photoi_X9NT,3)))[0,*],zz,name='Photoionization (Pi)',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$ ;layout=[2,3,5],
  xstyle=1,ystyle=1,yrange=yr,xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(eiionz_1_X9NT[0,*],zz,name='EI Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1c=plot(O_exctrate_X9NT,zz,name='EI Excitation',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

;Quiet

q1a=plot((reform(total(Sp_photoi_QNT,3)))[0,*],zz,name='Photoionization',$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

q1b=plot(eiionz_1_QNT[0,*],zz,name='Pe',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)

q1c=plot(O_exctrate_QNT,zz,name='EI Excitation',$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)


leg8a = LEGEND(TARGET=[d1a,d1b,d1c], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

o_r1=eiionz_1_X9NT[0,*]/(reform(total(Sp_photoi_X9NT,3)))[0,*]
o_r2=O_exctrate_X9NT/(reform(total(Sp_photoi_X9NT,3)))[0,*]

o_r3=eiionz_1_QNT[0,*]/(reform(total(Sp_photoi_QNT,3)))[0,*]
o_r4=O_exctrate_QNT/(reform(total(Sp_photoi_QNT,3)))[0,*]

O_r1[where(~finite(O_r1),/null)]=0.
O_r2[where(~finite(O_r2),/null)]=0.
O_r3[where(~finite(O_r3),/null)]=0.
O_r4[where(~finite(O_r4),/null)]=0.
w6cc = WINDOW(DIMENSIONS=[450,450])

xt='Ratio !C (Flare-Dashed,Quiet-Bold)'
e1a=plot(o_r1,zz,name='Pe $\slash$ Pi',$
  xtitle=xt,ytitle=yt,title=titl,margin=mrgn,$ ;layout=[2,3,6],
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.1,10000],xlog=1,ylog=0,xtickvalues=[0.1,1,10,100,1000,1.e4],xtickname=['0.1','1','10','10!U2!N','10!U3!N','10!U4!N'],$;
  thick=thck,color=colors[8],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

e1b=plot(O_r2,zz,name='EI Excitation$\slash$Pi ',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)


f1a=plot(o_r3,zz,name='Pe$\slash$Pi ',$
  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)

f1b=plot(O_r4,zz,name='EI Exct and Pi (Quiet) ',$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)


leg6c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


;fntsz=14


;__________________________________________________________________________________________________________________________________________________________________























































;____________________________________________________________________________________________________________________________________________________________________
;Average energy , max energy and  95% dissociation energy


xt=' Energy (eV)'
yt='Altitude (km)'


;xr=[1.,1.e5]
xr=[1.,2.5e4]
yr=[55,195.]


xt=' Energy (eV)'
yt='Altitude (km)'


;xr=[1.,1.e5]
xr=[10.,2.5e4]
yr=[55,195.]

w1 =  WINDOW(DIMENSIONS=[450,450])
;fntsz=14
;lgdsz=12
;1c. N2

titl=[];'Energy of 95% dissociation of N!L2!N '


d1a=plot(N2_sdener_X9NT,zz,name='X9 Flare',$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,margin=[0.2,0.2,0.2,0.2],$;layout=[3,2,3],
  thick=5,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

d1b=plot(N2_sdener_QNT,zz,name='Quiet',$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

leg11a = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)

;____________________________________________________________________________________________________________________________________________________________________































;____________________________________________________________________________________________________________________________________________________________________

;___________________________________________60 km__________________________________________________________
;__________________________________________________________________________________________________________

;fntsz=13
;lgdsz=11

plot_margin = [0.2, 0.2, 0.2, 0.2]

xr=[1.,1.e4]
yr1=[1.e-10,1.e3] ;Dissociation Rate

yr2=[0,5.5]
yr3=[0,7.0]

yr4=[1.e-10,1.e12]
yr5=[1.e-19,1.e-15]

leg_nam2=['X9Flare', 'Quiet']
leg_nam3=['PD(X9 Flare)','PD(Quiet)','CPD(X9 Flare)','CPD(Quiet)']
leg_nam4=['PE Flux(X9 Flare)','PE Flux(Quiet)',$
  'N!L2!N Crosssection','O!L2!N Crosssection']




w10a = WINDOW(DIMENSIONS=[450,450])
htind= 80;80;20

titl=[];'N!L2!N@ '+strtrim(string(fix(zz[htind])),2)+'km'


;thck=5

xt='Energy (eV)'



p1 = plot(ener, tflux_2_X9N[htind,*],name=leg_nam2[0], $
  xtitle=xt,ytitle="Photoelectron Flux",title=titl,$
  axis_style=1,yrange=[1.e-5,1.e15],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,/current)

p1b=plot(ener,tflux_2_QN[htind,*],name=leg_nam2[1],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

p2 = plot(ener,  N2_cumprcntdiss_X9NT[htind,*]*100. ,name='CPD '+leg_nam2[0],$
  axis_style=0,yrange=[0.,100.],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
  thick=thck,color='blue',linestyle=lstyl[1],font_size=fntsz,/current)

p2b=plot(ener,N2_cumprcntdiss_QNT[htind,*]*100.,name='CPD '+leg_nam2[1],$
    thick=2,color='blue',linestyle=lstyl[0],/overplot)




a2 = axis('y',target=p2, $
  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
  tickfont_size=fntsz,$
  tickfont_style=fntstyl,$
  textpos=1, $                           ; text faces outward
  tickdir=1, $                           ; ticks face inward
  title="Cumulative % Dissociation (CPD)")

l11a=legend(target=[p1,p1b,p2,p2b],font_size=lgdsz,transparency=transp,vertical_spacing=0.005,font_style=fntstyl)


























;_________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________

;Solar flux plot

restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5
ssflux_M5=peakm5
ssflux_X9=meanpeakx9
wv1_NRL=wave_gcm1
wv2_NRL=wave_gcm2

;Quiet spectrum with NRLFLARE resolution
restore,'/home/srimoyee/Desktop/nrl_files/sav_files/NRL_QuietSpectrum_July2022.sav'
;File has :nrl_wave1, nrl_wave2, nrl_current
ssflux_NRLQuiet=nrl_current


;1 solar flux
xr=[0.05,175.]
yr=[1.e-1,1.e12]
xt='Wavelength (nm)';$\AA$
yt='Flux (Photons cm!U-2!N s!U-1!N)'
titl=[];"Input Solar Spectra"
w2a=window(window_title=titl,dimension=[500,500])

s1=plot(wv1_NRL*0.1, ssflux_X9, $
  xtitle=xt,ytitle=yt,title=titl,name='X9 Flare',$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,margin=[0.2,0.2,0.2,0.2],$
  thick=3,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,histogram=1,/current) ;layout=[3,1,1],

s3=plot(wv1_NRL*0.1, ssflux_NRLQuiet,name='Quiet',$
  histogram=1,color=colors[1],thick=3,linestyle=lstyl[0],/overplot)

leg10 = LEGEND(TARGET=[s1,s3], POSITION=[0.2,0.9],font_size=lgdsz,/normal,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.15,font_style=fntstyl)

h=6.62607004e-34 ;Js
c=3.e8 ;m s-1

;NRL energy unit
NRLX9_eneunit =  (h*c*ssflux_x9*1.e14)/(0.5*(wv1_NRL+wv2_NRL)) ;W m-2
NRLM5_eneunit =  (h*c*ssflux_m5*1.e14)/(0.5*(wv1_NRL+wv2_NRL)) ;W m-2
NRLQuiet_eneunit =  (h*c*ssflux_NRLQuiet*1.e14)/(0.5*(wv1_NRL+wv2_NRL)) ;W m-2

yt='Solar Flux (W m!U-2!N  )'
xt='Wavelength (nm)';($\AA$)'
xr=[0.05,107.5]
yr=[1.e-10,1.e-2]
titl='NRL Solar Spectrum'

w2b = WINDOW(DIMENSIONS=[450,450])

n1=plot(wv1_nrl*0.1,nrlx9_eneunit,name='X9 Flare',$
  xtitle=xt,ytitle=yt,title=titl,margin=[0.18,0.15,0.15,0.15],$
  ;  ytickvalues=[1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.e0],ytickname=['1.e-5','1.e-4','1.e-3','1.e-2','1.e-1','1'],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,histogram=1,$;
  thick=3,color=colors[1],linestyle=lstyl[1],font_size=fntsz,font_style=1,/current)
n3=plot(wv1_nrl*0.1,nrlQuiet_eneunit,name='Quiet',$
    thick=3,color=colors[1],linestyle=lstyl[0],histogram=1,/overplot)

leg10 = LEGEND(TARGET=[n1,n3], POSITION=[0.2,0.9],font_size=lgdsz,/normal,transparency=transp,vertical_spacing=0.005,SAMPLE_WIDTH=0.15,font_style=fntstyl)

;____________________________________________________________________________________________________________________________________________________________________























end

; Testing the dissociation rates of N2 due to electron impact
;NOTE:Plot code for sample_localpe_v6.pro
;
;____________________________________________________________________________________________________________________________________________________________________
; O2 Excitation States
; 
; Strickland 1996 :
;                         predissociation
;0 vib           0.3 eV
;1 Ryds          16 eV    100%
;2 AA'c          4.5 eV   100%
;3 a1∆_u         1 eV
;4 b1Σ_g+        1.6 eV   
;5 Long Band     10 eV    100%
;6 1^3П_g        7.6 eV   100%
;7 B^3Σ_u        8.3 eV   100%
;8 Second Band   10.3eV   100%
;9 8.9           8.9 eV   100%

;GLOW O2 states:   
;0 a 
;1 b
;2 AA'c                    100%
;3 B,                      100%
;4 9.9                    100%
;5 Ryds                   100%
;6 vib


;____________________________________________________________________________________________________________________________________________________________________

pro newstates_NEW,Sp_exctrate,nei,Sp_ionzrate,nii,Sp_dissext,Sp_dissionz

  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  Sp_dissext=fltarr(nei,jmax,nbins,lmax)
  Sp_dissionz=fltarr(nii,jmax,nbins,lmax)

  Sp_dissext[0,*,*,*]= Sp_exctrate[1,*,*,*]
  Sp_dissext[1,*,*,*]= Sp_exctrate[2,*,*,*]
  Sp_dissext[2,*,*,*]= Sp_exctrate[5,*,*,*]
  Sp_dissext[3,*,*,*]= Sp_exctrate[6,*,*,*]
  Sp_dissext[4,*,*,*]= Sp_exctrate[7,*,*,*]
  Sp_dissext[5,*,*,*]= Sp_exctrate[8,*,*,*]
  Sp_dissext[6,*,*,*]= Sp_exctrate[9,*,*,*]
 

  Sp_dissionz[0,*,*,*]= Sp_ionzrate[4,*,*,*]
  Sp_dissionz[1,*,*,*]= Sp_ionzrate[5,*,*,*]
  Sp_dissionz[2,*,*,*]= Sp_ionzrate[6,*,*,*]

end


pro newstates_OLD,Sp_exctrate,nei,Sp_ionzrate,nii,Sp_dissext,Sp_dissionz

  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  Sp_dissext=fltarr(nei,jmax,nbins,lmax)
  Sp_dissionz=fltarr(nii,jmax,nbins,lmax)

  Sp_dissext[0,*,*,*]= Sp_exctrate[2,*,*,*]
  Sp_dissext[1,*,*,*]= Sp_exctrate[3,*,*,*]
  Sp_dissext[2,*,*,*]= Sp_exctrate[4,*,*,*]
  Sp_dissext[3,*,*,*]= Sp_exctrate[5,*,*,*]
  
  Sp_dissionz[0,*,*,*]= Sp_ionzrate[4,*,*,*]
  Sp_dissionz[1,*,*,*]= Sp_ionzrate[5,*,*,*]
  Sp_dissionz[2,*,*,*]= Sp_ionzrate[6,*,*,*]

end





pro disscontrb,Sp_exctrate,Sp_ionzrate,Sp_exctrate_1,Sp_exctrate_2,$
               Sp_exctrate_3,Sp_ionzrate_1,Sp_ionzrate_2,Sp_ionzrate_3,$
               Sp_prcnttotdisionz,Sp_prcnttotdisexct
  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  nii=(size(Sp_ionzrate))[1]



  Sp_exctrate_1 = fltarr(nei,jmax,nbins)
  Sp_exctrate_2 = fltarr(jmax,nbins)
  Sp_exctrate_3 = fltarr(jmax)
  Sp_exctrate_4 = fltarr(nei,jmax)

  Sp_ionzrate_1 = fltarr(nii,jmax,nbins)
  Sp_ionzrate_2 = fltarr(jmax,nbins)
  Sp_ionzrate_3 = fltarr(jmax)
  Sp_ionzrate_4 = fltarr(nii,jmax)




  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states



          Sp_exctrate_1[s,j,e]=Sp_exctrate_1[s,j,e]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_2[j,e]=Sp_exctrate_2[j,e]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_3[j]=Sp_exctrate_3[j]+Sp_exctrate[s,j,e,l]
          Sp_exctrate_4[s,j]=Sp_exctrate_4[s,j]+Sp_exctrate[s,j,e,l]





        endfor
      endfor
    endfor

  endfor


  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nii-1 do begin  ;states



          Sp_ionzrate_1[s,j,e]=Sp_ionzrate_1[s,j,e]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_2[j,e]=Sp_ionzrate_2[j,e]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_3[j]=Sp_ionzrate_3[j]+Sp_ionzrate[s,j,e,l]
          Sp_ionzrate_4[s,j]=Sp_ionzrate_4[s,j]+Sp_ionzrate[s,j,e,l]





        endfor
      endfor
    endfor

  endfor





  ; Percent contribution to dissociation-New

  Sp_prcnttotdisionz=100*Sp_ionzrate_3/(Sp_exctrate_3+Sp_ionzrate_3)
  Sp_prcnttotdisexct=100*Sp_exctrate_3/(Sp_exctrate_3+Sp_ionzrate_3)

  Sp_prcnttotdisionz[where(~finite(Sp_prcnttotdisionz),/null)]=0.
  Sp_prcnttotdisexct[where(~finite(Sp_prcnttotdisexct),/null)]=0.
  
         
end




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




pro max_ener,Sp_exctrate,ener,maxener,maxener_states,sdener,sdener_states,$
             Sp_rate_1,Sp_rate_2,cumprcntdiss,prcntdiss,totatmpercentdiss

  ;Calculates energy of max dissociation and 95% dissociation

  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]




  Sp_rate_1 = fltarr(nei,jmax,nbins)
  Sp_rate_2 = fltarr(jmax,nbins)

  Sp_tot_1 = fltarr(nei,jmax,nbins)
  Sp_tot_2 = fltarr(jmax,nbins)


  maxener = fltarr(jmax)
  maxener_states = fltarr(nei,jmax)

  sdener = fltarr(jmax)
  sdener_states = fltarr(nei,jmax)

  cumprcntdiss=fltarr(jmax,nbins) ;how much % dissociation we have gotten to at a particular ener bin
  prcntdiss=fltarr(jmax,nbins)
  totatmpercentdiss =fltarr(jmax,nbins)
  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states



          Sp_rate_1[s,j,e]=Sp_rate_1[s,j,e]+Sp_exctrate[s,j,e,l]
          Sp_rate_2[j,e]=Sp_rate_2[j,e]+Sp_exctrate[s,j,e,l]

          Sp_tot_1[s,j,e]=Sp_tot_1[s,j,e]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_tot_2[j,e]=Sp_tot_2[j,e]+Sp_exctrate[s,j,e,l]*ener[e]






        endfor
      endfor
    endfor

  endfor

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



  for j=0,jmax-1 do begin  ;altitude

    for s=0,nei-1 do begin  ;states
      mi2=max(reform(Sp_rate_1[s,j,*]),ind2)
      maxener_states[s,j]=ener[ind2]

      Sp_tot=fltarr(nbins)
      Sp_tot[0]=Sp_tot_1[s,j,0]

      for e=1,nbins-1 do Sp_tot[e]=total(Sp_tot_1[s,j,0:e],3)

      jj=where(  ( (total(Sp_tot_1[s,j,*])*0.95)-Sp_tot )  ge 0,/null)


      sdener_states[s,j]=ener[jj[-1]]



    endfor
  endfor




end

;____________________________________________________________________________________________________________________________________________________________________


sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"

fname=[sav_loc+'N2_FLAREOLD_Aug2021.sav',$        ;0 flare,old
       sav_loc+'N2_FLARENEW_Aug2021.sav',$        ;1 flare,new

       sav_loc+'N2_QUIETOLD_Aug2021.sav',$        ;2 quiet,old
       sav_loc+'N2_QUIETNEW_Aug2021.sav']         ;3 quiet,new


; file has:
;  Sp_exctrate1_transp,Sp_ionzrate1_transp,Sp_exctrate1_local,Sp_ionzrate1_local,
;  tflux_all,pespec_all,zz,ener,del,wv1,wv2,sza,lat,lon,tau,
;  sigix,sigex
;N2 ionisation
;Single dissociation -Non dissociative
;N+N_plus -Dissociative


plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"

;____________________________________________________________________________________________________________________________________________________________________

restore,fname[0]
newstates_OLD,Sp_exctrate1_transp,4,Sp_ionzrate1_transp,3,Sp_dissext_FOT,Sp_dissionz_FOT
newstates_OLD,Sp_exctrate1_local,4,Sp_ionzrate1_local,3,Sp_dissext_FOL,Sp_dissionz_FOL
tflux_all_FO=tflux_all
pespec_all_FO=pespec_all
tflux_FO=reform(total(tflux_all,3))
pespec_FO=reform(total(pespec_all,3))
sigex_OLD=sigex
sigix_OLD=sigix
totdisscrss_OLD= sigex_OLD[2,1,*]+sigex_OLD[3,1,*]+sigex_OLD[4,1,*]+sigex_OLD[5,1,*]+$
                 sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+ sigix_OLD[6,1,*]



;disscontrb,Sp_dissext_FOT,Sp_dissionz_FOT,Sp_exctrate_1_FOT,Sp_exctrate_2_FOT,Sp_exctrate_3_FOT,$
;          Sp_ionzrate_1_FOT,Sp_ionzrate_2_FOT,Sp_ionzrate_3_FOT,Sp_prcnttotdisionz_FOT,Sp_prcnttotdisexct_FOT
;
;
;disscontrb,Sp_dissext_FOL,Sp_dissionz_FOL,Sp_exctrate_1_FOL,Sp_exctrate_2_FOL,Sp_exctrate_3_FOL,$
;          Sp_ionzrate_1_FOL,Sp_ionzrate_2_FOL,Sp_ionzrate_3_FOL,Sp_prcnttotdisionz_FOL,Sp_prcnttotdisexct_FOL

;average_ener,Sp_dissext_FOT,ener,avg_ener_3_FOT,avg_ener_2_FOT,avg_ener_1_FOT,Sp_totrate_FOT
;average_ener,Sp_dissext_FOL,ener,avg_ener_3_FOL,avg_ener_2_FOL,avg_ener_1_FOL,Sp_totrate_FOL

max_ener,Sp_dissext_FOT,ener,maxener_FOT,maxener_states_FOT,sdener_FOT,sdener_states_FOT,Sp_rate_1_FOT,Sp_rate_2_FOT,cumprcntdiss_FOT,prcntdiss_FOT,totatmpercentdiss_FOT
max_ener,Sp_dissext_FOL,ener,maxener_FOL,maxener_states_FOL,sdener_FOL,sdener_states_FOL,Sp_rate_1_FOL,Sp_rate_2_FOL,cumprcntdiss_FOL,prcntdiss_FOL,totatmpercentdiss_FOL

restore,fname[1]
newstates_NEW,Sp_exctrate1_transp,7,Sp_ionzrate1_transp,3,Sp_dissext_FNT,Sp_dissionz_FNT
newstates_NEW,Sp_exctrate1_local,7,Sp_ionzrate1_local,3,Sp_dissext_FNL,Sp_dissionz_FNL
tflux_all_FN=tflux_all
pespec_all_FN=pespec_all
tflux_FN=reform(total(tflux_all,3))
pespec_FN=reform(total(pespec_all,3))
sigex_NEW=sigex
sigix_NEW=sigix

totdisscrss_NEW= sigex_NEW[1,1,*]+sigex_NEW[2,1,*]+sigex_NEW[5,1,*]+sigex_NEW[6,1,*]+sigex_NEW[7,1,*]+$
                 sigex_NEW[8,1,*]+ sigex_NEW[9,1,*]+sigix_NEW[4,1,*]+sigix_NEW[5,1,*]+ sigix_NEW[6,1,*]

;disscontrb,Sp_dissext_FNT,Sp_dissionz_FNT,Sp_exctrate_1_FNT,Sp_exctrate_2_FNT,Sp_exctrate_3_FNT,$
;  Sp_ionzrate_1_FNT,Sp_ionzrate_2_FNT,Sp_ionzrate_3_FNT,Sp_prcnttotdisionz_FNT,Sp_prcnttotdisexct_FNT
;
;
;disscontrb,Sp_dissext_FNL,Sp_dissionz_FNL,Sp_exctrate_1_FNL,Sp_exctrate_2_FNL,Sp_exctrate_3_FNL,$
;  Sp_ionzrate_1_FNL,Sp_ionzrate_2_FNL,Sp_ionzrate_3_FNL,Sp_prcnttotdisionz_FNL,Sp_prcnttotdisexct_FNL


;average_ener,Sp_dissext_FNT,ener,avg_ener_3_FNT,avg_ener_2_FNT,avg_ener_1_FNT,Sp_totrate_FNT
;average_ener,Sp_dissext_FNL,ener,avg_ener_3_FNL,avg_ener_2_FNL,avg_ener_1_FNL,Sp_totrate_FNL

max_ener,Sp_dissext_FNT,ener,maxener_FNT,maxener_states_FNT,sdener_FNT,sdener_states_FNT,Sp_rate_1_FNT,Sp_rate_2_FNT,cumprcntdiss_FNT,prcntdiss_FNT,totatmpercentdiss_FNT
max_ener,Sp_dissext_FNL,ener,maxener_FNL,maxener_states_FNL,sdener_FNL,sdener_states_FNL,Sp_rate_1_FNL,Sp_rate_2_FNL,cumprcntdiss_FNL,prcntdiss_FNL,totatmpercentdiss_FNL


restore,fname[2]
newstates_OLD,Sp_exctrate1_transp,4,Sp_ionzrate1_transp,3,Sp_dissext_QOT,Sp_dissionz_QOT
newstates_OLD,Sp_exctrate1_local,4,Sp_ionzrate1_local,3,Sp_dissext_QOL,Sp_dissionz_QOL
tflux_all_QO=tflux_all
pespec_all_QO=pespec_all
tflux_QO=reform(total(tflux_all,3))
pespec_QO=reform(total(pespec_all,3))

;disscontrb,Sp_dissext_QOT,Sp_dissionz_QOT,Sp_exctrate_1_QOT,Sp_exctrate_2_QOT,Sp_exctrate_3_QOT,$
;  Sp_ionzrate_1_QOT,Sp_ionzrate_2_QOT,Sp_ionzrate_3_QOT,Sp_prcnttotdisionz_QOT,Sp_prcnttotdisexct_QOT
;
;disscontrb,Sp_dissext_QOL,Sp_dissionz_QOL,Sp_exctrate_1_QOL,Sp_exctrate_2_QOL,Sp_exctrate_3_QOL,$
;          Sp_ionzrate_1_QOL,Sp_ionzrate_2_QOL,Sp_ionzrate_3_QOL,Sp_prcnttotdisionz_QOL,Sp_prcnttotdisexct_QOL

;average_ener,Sp_dissext_QOT,ener,avg_ener_3_QOT,avg_ener_2_QOT,avg_ener_1_QOT,Sp_totrate_QOT
;average_ener,Sp_dissext_QOL,ener,avg_ener_3_QOL,avg_ener_2_QOL,avg_ener_1_QOL,Sp_totrate_QOL

max_ener,Sp_dissext_QOT,ener,maxener_QOT,maxener_states_QOT,sdener_QOT,sdener_states_QOT,Sp_rate_1_QOT,Sp_rate_2_QOT,cumprcntdiss_QOT,prcntdiss_QOT,totatmpercentdiss_QOT
max_ener,Sp_dissext_QOL,ener,maxener_QOL,maxener_states_QOL,sdener_QOL,sdener_states_QOL,Sp_rate_1_QOL,Sp_rate_2_QOL,cumprcntdiss_QOL,prcntdiss_QOL,totatmpercentdiss_QOL



restore,fname[3]
newstates_NEW,Sp_exctrate1_transp,7,Sp_ionzrate1_transp,3,Sp_dissext_QNT,Sp_dissionz_QNT
newstates_NEW,Sp_exctrate1_local,7,Sp_ionzrate1_local,3,Sp_dissext_QNL,Sp_dissionz_QNL
tflux_all_QN=tflux_all
pespec_all_QN=pespec_all
tflux_QN=reform(total(tflux_all,3))
pespec_QN=reform(total(pespec_all,3))


;disscontrb,Sp_dissext_QNT,Sp_dissionz_QNT,Sp_exctrate_1_QNT,Sp_exctrate_2_QNT,Sp_exctrate_3_QNT,$
;  Sp_ionzrate_1_QNT,Sp_ionzrate_2_QNT,Sp_ionzrate_3_QNT,Sp_prcnttotdisionz_QNT,Sp_prcnttotdisexct_QNT
;
; disscontrb,Sp_dissext_QNL,Sp_dissionz_QNL,Sp_exctrate_1_QNL,Sp_exctrate_2_QNL,Sp_exctrate_3_QNL,$
;    Sp_ionzrate_1_QNL,Sp_ionzrate_2_QNL,Sp_ionzrate_3_QNL,Sp_prcnttotdisionz_QNL,Sp_prcnttotdisexct_QNL

;average_ener,Sp_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,Sp_totrate_QNT
;average_ener,Sp_dissext_QNL,ener,avg_ener_3_QNL,avg_ener_2_QNL,avg_ener_1_QNL,Sp_totrate_QNL

max_ener,Sp_dissext_QNT,ener,maxener_QNT,maxener_states_QNT,sdener_QNT,sdener_states_QNT,Sp_rate_1_QNT,Sp_rate_2_QNT,cumprcntdiss_QNT,prcntdiss_QNT,totatmpercentdiss_QNT
max_ener,Sp_dissext_QNL,ener,maxener_QNL,maxener_states_QNL,sdener_QNL,sdener_states_QNL,Sp_rate_1_QNL,Sp_rate_2_QNL,cumprcntdiss_QNL,prcntdiss_QNL,totatmpercentdiss_QNL


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
fntsz=18
lgdsz=15
transp=50

colors=['red','dodger blue','spring green','purple','gray','deep pink','firebrick','medium blue','orange','fuchsia','dark khaki','dark slate gray','gold']
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

;;xt='Percentage contribution'
;;yt='Altitude (km)'
;;
;;xr=[0.,100.]
;; yr=[min(zz),max(zz)]
;;
;;;
;;leg_nam1=['Excitation-Flare&Strickland','Ionization-Flare&Strickland','Excitation-Flare&Glow','Ionization-Flare&Glow',$
;;           'Excitation-Quiet&Strickland','Ionization-Quiet&Strickland','Excitation-Quiet&Glow','Ionization-Quiet&Glow']
;;titl='Percentage contribution to dissociation of  O!L2!N (Transport)'
;;w1 = WINDOW(DIMENSIONS=[1000,800])
;;
;;
;;f1=plot(Sp_prcnttotdisexct_FNT,zz,name=leg_nam1[0],$
;;  xtitle=xt,ytitle=yt,title=titl,$
;;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;;  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)
;;
;;f2=plot(Sp_prcnttotdisionz_FNT,zz,name=leg_nam1[1],$
;;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;;
;;f3=plot(Sp_prcnttotdisexct_FOT,zz,name=leg_nam1[2],$
;;    thick=thck,color=colors[0],linestyle=lstyl[1],/overplot)
;;
;;
;;f4=plot(Sp_prcnttotdisionz_FOT,zz,name=leg_nam1[3],$
;;    thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;;
;;
;;q1=plot(Sp_prcnttotdisexct_QNT,zz,name=leg_nam1[4],$
;;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;;
;;q2=plot(Sp_prcnttotdisionz_QNT,zz,name=leg_nam1[5],$
;;    thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;;
;;q3=plot(Sp_prcnttotdisexct_QOT,zz,name=leg_nam1[6],$
;;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;;
;;q4=plot(Sp_prcnttotdisionz_QOT,zz,name=leg_nam1[7],$
;;    thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;;
;;
;;leg1 = LEGEND(TARGET=[f1,f2,f3,f4,q1,q2,q3,q4], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;;
;
;
;
;
;titl='Percentage contribution to dissociation of  O!L2!N (Local)'
;w2 = WINDOW(DIMENSIONS=[1000,800])
;
;
;f1=plot(Sp_prcnttotdisexct_FNL,zz,name=leg_nam1[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)
;
;f2=plot(Sp_prcnttotdisionz_FNL,zz,name=leg_nam1[1],$
;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;
;f3=plot(Sp_prcnttotdisexct_FOL,zz,name=leg_nam1[2],$
;  thick=thck,color=colors[0],linestyle=lstyl[1],/overplot)
;
;
;f4=plot(Sp_prcnttotdisionz_FOL,zz,name=leg_nam1[3],$
;  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;
;
;q1=plot(Sp_prcnttotdisexct_QNL,zz,name=leg_nam1[4],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;q2=plot(Sp_prcnttotdisionz_QNL,zz,name=leg_nam1[5],$
;  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;
;q3=plot(Sp_prcnttotdisexct_QOL,zz,name=leg_nam1[6],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;q4=plot(Sp_prcnttotdisionz_QOL,zz,name=leg_nam1[7],$
;  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
;
;leg2 = LEGEND(TARGET=[f1,f2,f3,f4,q1,q2,q3,q4], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;



;____________________________________________________________________________________________________________________________________________________________________

;;yt='Dissociative Excitation cross-section (cm!U2!N)'
;;xt='Energy (eV)'
;;xr=[1., 5.e4]
;;yr=[];[0.,100.01]
;;titl='Dissociative Excitation cross-section'
;;
;;
;;w3 = WINDOW(DIMENSIONS=[1000,800])
;;
;;
;;c1= plot(ener,totdisscrss_NEW- (sigix_NEW[4,1,*]+sigix_NEW[5,1,*]+ sigix_NEW[6,1,*]),name='Strickland 1996',$
;;  xtitle=xt,ytitle=yt,title=titl,$
;;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;;  thick=thck,color=colors[1],linestyle=lstyl[0],font_size=fntsz,/current)
;;  
;;c2=plot(ener,totdisscrss_OLD- (sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+ sigix_OLD[6,1,*]),name='GLOW',$
;;    thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;;
;;leg3 = LEGEND(TARGET=[c1,c2], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;;
;;;Calculating percentage contribution to dissociation-Strickland
;;
;;
;;;0 vib           0.3 eV
;;;1 Ryds          16 eV    100%
;;;2 AA'c          4.5 eV   100%
;;;3 a1∆_u         1 eV
;;;4 b1Σ_g+        1.6 eV
;;;5 Long Band     10 eV    100%
;;;6 1^3П_g        7.6 eV   100%
;;;7 B^3Σ_u        8.3 eV   100%
;;;8 Second Band   10.3eV   100%
;;;9 8.9           8.9 eV   100%
;;
;;
;;
;;percentdiss_1=100.*sigex_NEW[1,1,*]/totdisscrss_NEW
;;percentdiss_2=100.*sigex_NEW[2,1,*]/totdisscrss_NEW
;;percentdiss_3=100.*sigex_NEW[5,1,*]/totdisscrss_NEW
;;percentdiss_4=100.*sigex_NEW[6,1,*]/totdisscrss_NEW
;;
;;percentdiss_5=100.*sigex_NEW[7,1,*]/totdisscrss_NEW
;;percentdiss_6=100.*sigex_NEW[8,1,*]/totdisscrss_NEW
;;percentdiss_7=100.*sigex_NEW[9,1,*]/totdisscrss_NEW
;;
;;percentdiss_ionz=100.*(sigix_NEW[4,2,*]+sigix_NEW[5,2,*]+sigix_NEW[6,2,*])/totdisscrss_NEW
;;
;;
;;
;;
;;percentdiss_1[where(~finite(percentdiss_1),/null    )] =0.
;;percentdiss_2[where(~finite(percentdiss_2),/null    )] =0.
;;percentdiss_3[where(~finite(percentdiss_3),/null    )] =0.
;;percentdiss_4[where(~finite(percentdiss_4),/null    )] =0.
;;percentdiss_5[where(~finite(percentdiss_5),/null    )] =0.
;;percentdiss_6[where(~finite(percentdiss_6),/null    )] =0.
;;percentdiss_7[where(~finite(percentdiss_7),/null    )] =0.
;;percentdiss_ionz[where(~finite(percentdiss_ionz),/null    )] =0.
;;
;;
;;leg_nam=['Rydbergs', "A+A'+c",'Long Band','1!U3!N $\Pi$ !Lg!N',"B'!U3!N $\Sigma$ !Lu!N !U-!N",$
;;         "Second Band",'8.9 eV']
;;
;;yt='Percentage contribution'
;;xt='Energy (eV)'
;;xr=[1., 2.5e4]
;;yr=[0.,100.01]
;;titl='Percentage contribution of states to dissociation cross-section (Strickland 1996)'
;;
;;
;;
;;
;;w4 = WINDOW(DIMENSIONS=[1000,800])
;;
;;c1= plot(ener,percentdiss_1,name=leg_nam[0],$
;;  xtitle=xt,ytitle=yt,title=titl,$
;;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;;  thick=thck,color=colors[4],linestyle=lstyl[0],font_size=fntsz,/current)
;;
;;c2=plot(ener,percentdiss_2,name=leg_nam[1],$
;;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;;
;;c3= plot(ener,percentdiss_3,name=leg_nam[2],$
;;  thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)
;;
;;c4=plot(ener,percentdiss_4,name=leg_nam[3],$
;;  thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)
;;
;;c5=plot(ener,percentdiss_5,name=leg_nam[4],$
;;  thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)
;;
;;c6=plot(ener,percentdiss_6,name=leg_nam[5],$
;;  thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)
;;
;;c7=plot(ener,percentdiss_7,name=leg_nam[6],$
;;  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;;
;;c8=plot(ener,percentdiss_ionz,name='Dissociative Ionisation',$
;;  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;;
;;
;;leg4 = LEGEND(TARGET=[c1,c2,c3,c4,c5,c6,c7,c8], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;;
;;
;
;
;
;
;;Calculating percentage contribution to dissociation-GLOW
;
;
;;GLOW O2 states:
;;0 a
;;1 b
;;2 AA'c                    100%
;;3 B,                      100%
;;4 9.9                    100%
;;5 Ryds                   100%
;;6 vib
;
;
;
;percentdiss_1=100.*sigex_OLD[2,1,*]/totdisscrss_OLD
;percentdiss_2=100.*sigex_OLD[3,1,*]/totdisscrss_OLD
;percentdiss_3=100.*sigex_OLD[4,1,*]/totdisscrss_OLD
;percentdiss_4=100.*sigex_OLD[5,1,*]/totdisscrss_OLD
;
;percentdiss_ionz=100.*(sigix_OLD[4,2,*]+sigix_OLD[5,2,*]+sigix_OLD[6,2,*])/totdisscrss_OLD
;
;
;
;
;percentdiss_1[where(~finite(percentdiss_1),/null    )] =0.
;percentdiss_2[where(~finite(percentdiss_2),/null    )] =0.
;percentdiss_3[where(~finite(percentdiss_3),/null    )] =0.
;percentdiss_4[where(~finite(percentdiss_4),/null    )] =0.
;percentdiss_ionz[where(~finite(percentdiss_ionz),/null    )] =0.
;
;
;
;
;
;leg_nam=[ "A+A'+c","B'!U3!N $\Sigma$ !Lu!N !U-!N",'9.9 eV','Rydbergs']
;
;yt='Percentage contribution'
;xt='Energy (eV)'
;xr=[1., 2.5e4]
;yr=[0.,100.01]
;titl='Percentage contribution of states to dissociation cross-section (GLOW)'
;
;
;
;
;w5 = WINDOW(DIMENSIONS=[1000,800])
;
;c1= plot(ener,percentdiss_1,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;c2=plot(ener,percentdiss_2,name=leg_nam[1],$
;  thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)
;
;c3= plot(ener,percentdiss_3,name=leg_nam[2],$
;  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;
;c4=plot(ener,percentdiss_4,name=leg_nam[3],$
;  thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;c5=plot(ener,percentdiss_ionz,name='Dissociative Ionisation',$
;  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg5 = LEGEND(TARGET=[c1,c2,c3,c4,c5], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;







;_____________________________________________________________________________________________________


;titl='Total Dissociation rate of O!L2!N
;yt=' Altitude (km)'
;xt='Photoelectron Dissociation rate (cm!U-3!N s!U-1!N)'
;
;w6 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(Sp_totrate_FOT,zz,name='GLOW,Flare',$
;  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
;  xstyle=1,ystyle=1,yrange=[],xrange=[],xlog=1,ylog=0,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(Sp_totrate_FNT,zz,name='Strickland,Flare',$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(Sp_totrate_QOT,zz,name='GLOW,Quiet',$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(Sp_totrate_QNT,zz,name='Strickland,Flare',$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;titl='Ratio of change in  dissociation rate of O!L2!N
;xt='Ratio (GLOW/Strickalnd)'
;e1a=plot(Sp_totrate_FOT/Sp_totrate_FNT,zz,name='Flare',$
;  xtitle=xt,ytitle=[],title=titl,position=postright,$
;  xstyle=1,ystyle=1,yrange=[],xrange=[1.,2.],xlog=0,ylog=0,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;e1b=plot(Sp_totrate_QOT/Sps_totrate_QNT,zz,name='Quiet',$
;  thick=2,color=colors[-2],linestyle=lstyl[0],/overplot)
;
;
;leg6 = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal)
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;leg_nam2=['Flare,Glow','Flare,Strickland ',$
;  'Quiet,Glow','Quiet,Strickland']
;
;
;xt='Average Energy (eV)'
;yt='Altitude (km)'
;
;
;xr=[1.,10000.]
;yr=[min(zz),max(zz)]
;titl='Average energy of dissociation of O!L2!N (Transport) '
;
;
;w7 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(avg_ener_1_FOT,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(avg_ener_1_FNT,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(avg_ener_1_QOT,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(avg_ener_1_QNT,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg7 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;

;titl='Average energy of dissociation of O!L2!N (Local) '
;
;
;w8 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(avg_ener_1_FOL,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(avg_ener_1_FNL,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(avg_ener_1_QOL,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(avg_ener_1_QNL,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg8 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;____________________________________________________________________________________________________________________________________________________________________

;leg_nam2=['Flare,Glow','Flare,Strickland ',$
;  'Quiet,Glow','Quiet,Strickland']
;
;
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;xr=[1.,1.e4]
;;xr=[10.,10000.]
;yr=[min(zz),max(zz)]
;titl='Energy of maximum dissociation of O!L2!N (Transport) '
;
;
;w9 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(maxener_FOT,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(maxener_FNT,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(maxener_QOT,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(maxener_QNT,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg9 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;titl='Energy of maximum dissociation of O!L2!N (Local) '
;
;
;w10 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(maxener_FOL,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(maxener_FNL,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(maxener_QOL,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(maxener_QNL,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg10 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;leg_nam2=['Flare,Glow','Flare,Strickland ',$
;  'Quiet,Glow','Quiet,Strickland']
;
;
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;xr=[1.,1.e5]
;;xr=[10.,10000.]
;yr=[min(zz),max(zz)]
;titl='Energy of 95% dissociation of O!L2!N (Transport) '
;
;
;w11 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(sdener_FOT,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(sdener_FNT,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(sdener_QOT,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(sdener_QNT,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg11 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;titl='Energy of 95% dissociation of O!L2!N (Local) '
;
;
;w12 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(sdener_FOL,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(sdener_FNL,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(sdener_QOL,zz,name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(sdener_QNL,zz,name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;leg12 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)



;_____________________________________________________________________________________________________

leg_nam2=['Flare,Glow','Flare,Strickland ',$
  'Quiet,Glow','Quiet,Strickland']

xt=' Energy (eV)'
yt='Photoelectron Dissociation rate (cm!U-3!N s!U-1!N)'

xr=[1.,1.e4]
yr=[]

htind=150

titl='Dissociation rate of O!L2!N at altitiude '+string(zz[htind])+'km'



w10a = WINDOW(DIMENSIONS=[1000,1000])

d1a=plot(ener,Sp_rate_2_FOT[htind,*],name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

d1b=plot(ener,Sp_rate_2_FNT[htind,*],name=leg_nam2[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1c=plot(ener,Sp_rate_2_QOT[htind,*],name=leg_nam2[2],$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

d1d=plot(ener,Sp_rate_2_QNT[htind,*],name=leg_nam2[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)




;w10a = WINDOW(DIMENSIONS=[1000,1000])
;
;d1a=plot(ener,Sp_rate_2_FOL[htind,*],name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
;  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(ener,Sp_rate_2_FNL[htind,*],name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(ener,Sp_rate_2_QOL[htind,*],name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,Sp_rate_2_QNL[htind,*],name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;



;_____________________________________________________________________________________________________


xt='Energy (eV)'
xr=[0.1,1.e4]
leg_nam3=['Photoelectron Flux(GLOW,Flare)','Photoelectron Flux(Strick,Flare)','Photoelectron Flux(GLOW,Quiet)','Photoelectron Flux(Strick,Quiet)',$
  'Cross-section(GLOW)','Cross-section(Strick)']


titl='Dissociation of O!L2!N at altitiude '+string(zz[htind])+'km'
w11 = WINDOW(DIMENSIONS=[1000,1000])


plot_margin = [0.15, 0.25, 0.15, 0.15]

p1 = plot(ener, reform(tflux_FO[htind,*]),name=leg_nam3[0], $
  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N eV!U-1!N )",title=titl,$
  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)


p1b=plot(ener,reform(tflux_FN[htind,*]),name=leg_nam3[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

p1c=plot(ener,reform(tflux_QO[htind,*]),name=leg_nam3[2],$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

p1d=plot(ener,reform(tflux_QN[htind,*]),name=leg_nam3[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)


p2 = plot(ener, totdisscrss_OLD- (sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+ sigix_OLD[6,1,*]),name=leg_nam3[4] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)

p2b=plot(ener,totdisscrss_NEW- (sigix_NEW[4,1,*]+sigix_NEW[5,1,*]+ sigix_NEW[6,1,*]),name=leg_nam3[5],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)



a2 = axis('y',target=p2, $
  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
  tickfont_size=fntsz,$
  textpos=1, $                           ; text faces outward
  tickdir=1, $                           ; ticks face inward
  title='Cross-section (cm!U2!N)')

l11a=legend(target=[p1,p1b,p1c],font_size=lgdsz,transparency=transp)
l11b=legend(target=[p1d,p2,p2b],font_size=lgdsz,transparency=transp)









;p1 = plot(ener, reform(pespec_FO[htind,*]),name=leg_nam3[0], $
;  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N eV!U-1!N )",title=titl,$
;  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;p1b=plot(ener,reform(pespec_FN[htind,*]),name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;p1c=plot(ener,reform(pespec_QO[htind,*]),name=leg_nam3[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;p1d=plot(ener,reform(pespec_QN[htind,*]),name=leg_nam3[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;p2 = plot(ener, totdisscrss_OLD- (sigix_OLD[4,1,*]+sigix_OLD[5,1,*]+ sigix_OLD[6,1,*]),name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,totdisscrss_NEW- (sigix_NEW[4,1,*]+sigix_NEW[5,1,*]+ sigix_NEW[6,1,*]),name=leg_nam3[5],$
;  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;a2 = axis('y',target=p2, $
;  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;  tickfont_size=fntsz,$
;  textpos=1, $                           ; text faces outward
;  tickdir=1, $                           ; ticks face inward
;  title='Cross-section (cm!U2!N)')
;
;l11a=legend(target=[p1,p1b,p1c],font_size=lgdsz,transparency=transp)
;l11b=legend(target=[p1d,p2,p2b],font_size=lgdsz,transparency=transp)



;_____________________________________________________________________________________________________

thck=4

xt='Energy (eV)'
xr=[1.,1.e4]
leg_nam3=['Percent Diss(GLOW,Flare)','Percent Diss(Strick,Flare)','Percent Diss(GLOW,Quiet)','Percent Diss(Strick,Quiet)',$
  'Cumulative Percent Diss(GLOW,Flare)','Cumulative Percent Diss(Strick,Flare)','Cumulative Percent Diss(GLOW,Quiet)','Cumulative Percent Diss(Strick,Quiet)']


sym=['+','*']
titl='Dissociation of O!L2!N at altitiude '+string(zz[htind])+'km'

w11 = WINDOW(DIMENSIONS=[1000,1000])

plot_margin = [0.15, 0.25, 0.15, 0.15]

p1 = plot(ener, prcntdiss_FOT[htind,*]*100.,name=leg_nam3[0], $
  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

p1b=plot(ener,prcntdiss_FNT[htind,*]*100.,name=leg_nam3[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

p2 = plot(ener, cumprcntdiss_FOT[htind,*]*100.,name=leg_nam3[4] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],linestyle=lstyl[1],font_size=fntsz,/current)

p2b=plot(ener,cumprcntdiss_FNT[htind,*]*100.,name=leg_nam3[5],$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

ar1 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
t2 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Cumulative', /normal,font_size=fntsz,/overplot)

a2 = axis('y',target=p2, $
  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
  tickfont_size=fntsz,$
  textpos=1, $                           ; text faces outward
  tickdir=1, $                           ; ticks face inward
  title='Percentage Dissociation (Cumulative)')

l11a=legend(target=[p1,p1b,p2,p2b],font_size=lgdsz,transparency=transp)




;p1 = plot(ener, prcntdiss_FOL[htind,*]*100.,name=leg_nam3[0], $
;  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
;  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;p1b=plot(ener,prcntdiss_FNL[htind,*]*100.,name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;p2 = plot(ener, cumprcntdiss_FOL[htind,*]*100.,name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[6],linestyle=lstyl[1],font_size=fntsz,/current)
;
;p2b=plot(ener,cumprcntdiss_FNL[htind,*]*100.,name=leg_nam3[5],$
;  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)
;
;ar1 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
;t2 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Cumulative', /normal,font_size=fntsz,/overplot)
;
;a2 = axis('y',target=p2, $
;  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;  tickfont_size=fntsz,$
;  textpos=1, $                           ; text faces outward
;  tickdir=1, $                           ; ticks face inward
;  title='Percentage Dissociation (Cumulative)')
;
;l11a=legend(target=[p1,p1b,p2,p2b],font_size=lgdsz,transparency=transp)



w12 = WINDOW(DIMENSIONS=[1000,1000])

plot_margin = [0.15, 0.25, 0.15, 0.15]

p1 = plot(ener, prcntdiss_QOT[htind,*]*100.,name=leg_nam3[2], $
  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=2,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)

p1b=plot(ener,prcntdiss_QNT[htind,*]*100.,name=leg_nam3[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

p2 = plot(ener, cumprcntdiss_QOT[htind,*]*100.,name=leg_nam3[6] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=2,color=colors[6],linestyle=lstyl[0],font_size=fntsz,/current)

p2b=plot(ener,cumprcntdiss_QNT[htind,*]*100.,name=leg_nam3[7],$
    thick=2,color=colors[7],linestyle=lstyl[0],/overplot)

ar1 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
t2 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Cumulative', /normal,font_size=fntsz,/overplot)

a2 = axis('y',target=p2, $
  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
  tickfont_size=fntsz,$
  textpos=1, $                           ; text faces outward
  tickdir=1, $                           ; ticks face inward
  title='Percentage Dissociation (Cumulative)')

l11a=legend(target=[p1,p1b,p2,p2b],font_size=lgdsz,transparency=transp)




;p1 = plot(ener, prcntdiss_QOL[htind,*]*100.,name=leg_nam3[2], $
;  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
;  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=2,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p1b=plot(ener,prcntdiss_QNL[htind,*]*100.,name=leg_nam3[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;p2 = plot(ener, cumprcntdiss_QOL[htind,*]*100.,name=leg_nam3[6] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=2,color=colors[6],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,cumprcntdiss_QNL[htind,*]*100.,name=leg_nam3[7],$
;  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
;
;ar1 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
;t2 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Cumulative', /normal,font_size=fntsz,/overplot)
;
;a2 = axis('y',target=p2, $
;  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;  tickfont_size=fntsz,$
;  textpos=1, $                           ; text faces outward
;  tickdir=1, $                           ; ticks face inward
;  title='Percentage Dissociation (Cumulative)')
;
;l11a=legend(target=[p1,p1b,p2,p2b],font_size=lgdsz,transparency=transp)


;_________________________________________________________________________________________________________________________________

;hind=[17,20,30,45,60,150]
;
;leg_nam4=[string(fix(zz[hind[0]]))+'km',$
;  string(fix(zz[hind[1]]))+'km',$
;  string(fix(zz[hind[2]]))+'km',$
;  string(fix(zz[hind[3]]))+'km',$
;  string(fix(zz[hind[4]]))+'km',$
;  string(fix(zz[hind[5]]))+'km']
;
;xt=' Energy (eV)'
;yt='Percentage dissociation'
;
;xr=[1.,1.e4]
;;yr=[0.001,100.]
;
;
;totatmpercentdiss=totatmpercentdiss_FOT
;titl='Percentage Dissociative Excitationof N!L2!N !C(Flare/GLOW Cross-section/Transport)'
;
;
;
;w12 = WINDOW(DIMENSIONS=[1000,1000])
;
;d1a=plot(ener,totatmpercentdiss[hind[0],*],name=leg_nam4[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
;  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[2],linestyle=lstyl[0],font_size=fntsz,/current)
;
;d1b=plot(ener,totatmpercentdiss[hind[1],*],name=leg_nam4[1],$
;  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;
;d1c=plot(ener,totatmpercentdiss[hind[2],*],name=leg_nam4[2],$
;  thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,totatmpercentdiss[hind[3],*],name=leg_nam4[3],$
;  thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)
;
;d1e=plot(ener,totatmpercentdiss[hind[4],*],name=leg_nam4[4],$
;  thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)
;
;d1f=plot(ener,totatmpercentdiss[hind[5],*],name=leg_nam4[5],$
;  thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)
;
;leg12 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;








end

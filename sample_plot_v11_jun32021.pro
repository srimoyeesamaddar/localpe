; Testing the dissociation rates due to electron impact
;NOTE:Plot code for sample_localpe_v6.pro
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
;____________________________________________________________________________________________________________________________________________________________________

pro newstates_NEW,N2_exctrate,nei,N2_ionzrate,nii,N2_dissext,N2_dissionz

  jmax=(size(N2_exctrate))[2]
  nbins=(size(N2_exctrate))[3]
  lmax=(size(N2_exctrate))[4]

  N2_dissext=fltarr(nei,jmax,nbins,lmax)
  N2_dissionz=fltarr(nii,jmax,nbins,lmax)

  N2_dissext[0,*,*,*]= N2_exctrate[5,*,*,*]
  N2_dissext[1,*,*,*]=N2_exctrate[6,*,*,*]
  N2_dissext[2,*,*,*]=N2_exctrate[9,*,*,*]
  N2_dissext[3,*,*,*]=N2_exctrate[10,*,*,*]
  N2_dissext[4,*,*,*]=N2_exctrate[11,*,*,*]
  N2_dissext[5,*,*,*]=N2_exctrate[12,*,*,*]
  N2_dissext[6,*,*,*]=N2_exctrate[13,*,*,*]
  N2_dissext[7,*,*,*]=N2_exctrate[14,*,*,*]
  N2_dissext[8,*,*,*]=N2_exctrate[16,*,*,*]
  N2_dissext[9,*,*,*]=N2_exctrate[17,*,*,*]
  N2_dissext[10,*,*,*]=N2_exctrate[18,*,*,*]


  N2_dissionz[0,*,*,*]= N2_ionzrate[3,*,*,*]
  N2_dissionz[1,*,*,*]=N2_ionzrate[4,*,*,*]
  N2_dissionz[2,*,*,*]=N2_ionzrate[5,*,*,*]

end


pro newstates_OLD,N2_exctrate,nei,N2_ionzrate,nii,N2_dissext,N2_dissionz

  jmax=(size(N2_exctrate))[2]
  nbins=(size(N2_exctrate))[3]
  lmax=(size(N2_exctrate))[4]

  N2_dissext=fltarr(nei,jmax,nbins,lmax)
  N2_dissionz=fltarr(nii,jmax,nbins,lmax)

  N2_dissext[0,*,*,*]= N2_exctrate[4,*,*,*]
  N2_dissext[1,*,*,*]=N2_exctrate[5,*,*,*]
  N2_dissext[2,*,*,*]=N2_exctrate[6,*,*,*]

  N2_dissionz[0,*,*,*]= N2_ionzrate[3,*,*,*]
  N2_dissionz[1,*,*,*]=N2_ionzrate[4,*,*,*]
  N2_dissionz[2,*,*,*]=N2_ionzrate[5,*,*,*]

end





pro disscontrb,N2_exctrate,N2_ionzrate,N2_exctrate_1,N2_exctrate_2,N2_exctrate_3,N2_ionzrate_1,N2_ionzrate_2,N2_ionzrate_3,N2_prcnttotdisionz,N2_prcnttotdisexct
  nei=(size(N2_exctrate))[1]
  jmax=(size(N2_exctrate))[2]
  nbins=(size(N2_exctrate))[3]
  lmax=(size(N2_exctrate))[4]

  nii=(size(N2_ionzrate))[1]



  N2_exctrate_1 = fltarr(nei,jmax,nbins)
  N2_exctrate_2 = fltarr(jmax,nbins)
  N2_exctrate_3 = fltarr(jmax)

  N2_ionzrate_1 = fltarr(nii,jmax,nbins)
  N2_ionzrate_2 = fltarr(jmax,nbins)
  N2_ionzrate_3 = fltarr(jmax)





  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states



          N2_exctrate_1[s,j,e]=N2_exctrate_1[s,j,e]+N2_exctrate[s,j,e,l]
          N2_exctrate_2[j,e]=N2_exctrate_2[j,e]+N2_exctrate[s,j,e,l]
          N2_exctrate_3[j]=N2_exctrate_3[j]+N2_exctrate[s,j,e,l]






        endfor
      endfor
    endfor

  endfor


  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nii-1 do begin  ;states



          N2_ionzrate_1[s,j,e]=N2_ionzrate_1[s,j,e]+N2_ionzrate[s,j,e,l]
          N2_ionzrate_2[j,e]=N2_ionzrate_2[j,e]+N2_ionzrate[s,j,e,l]
          N2_ionzrate_3[j]=N2_ionzrate_3[j]+N2_ionzrate[s,j,e,l]






        endfor
      endfor
    endfor

  endfor





  ; Percent contribution to dissociation-New

  N2_prcnttotdisionz=100*N2_ionzrate_3/(N2_exctrate_3+N2_ionzrate_3)
  N2_prcnttotdisexct=100*N2_exctrate_3/(N2_exctrate_3+N2_ionzrate_3)

  N2_prcnttotdisionz[where(~finite(N2_prcnttotdisionz),/null)]=0.
  N2_prcnttotdisexct[where(~finite(N2_prcnttotdisexct),/null)]=0.

end




pro average_ener,N2_exctrate,ener,avg_ener_3,avg_ener_2,avg_ener_1
  nei=(size(N2_exctrate))[1]
  jmax=(size(N2_exctrate))[2]
  nbins=(size(N2_exctrate))[3]
  lmax=(size(N2_exctrate))[4]



  N2_tot_3 = fltarr(nei,jmax,lmax)
  N2_tot_2 = fltarr(nei,jmax)
  N2_tot_1 = fltarr(jmax)

  N2_rate_3 = fltarr(nei,jmax,lmax)
  N2_rate_2 = fltarr(nei,jmax)
  N2_rate_1 = fltarr(jmax)

  avg_ener_3 = fltarr(nei,jmax,lmax)
  avg_ener_2 = fltarr(nei,jmax)
  avg_ener_1 = fltarr(jmax)




  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states

          N2_tot_3[s,j,l]=N2_tot_3[s,j,l]+N2_exctrate[s,j,e,l]*ener[e]
          N2_tot_2[s,j]=N2_tot_2[s,j]+N2_exctrate[s,j,e,l]*ener[e]
          N2_tot_1[j]=N2_tot_1[j]+N2_exctrate[s,j,e,l]*ener[e]


          N2_rate_3[s,j,l]=N2_rate_3[s,j,l]+N2_exctrate[s,j,e,l]
          N2_rate_2[s,j]=N2_rate_2[s,j]+N2_exctrate[s,j,e,l]
          N2_rate_1[j]=N2_rate_1[j]+N2_exctrate[s,j,e,l]




        endfor
      endfor
    endfor

  endfor

  avg_ener_3=N2_tot_3/N2_rate_3
  avg_ener_2=N2_tot_2/N2_rate_2
  avg_ener_1=N2_tot_1/N2_rate_1




end




pro max_ener,N2_exctrate,ener,maxener,maxener_states,sdener,sdener_states,N2_rate_1,N2_rate_2,cumprcntdiss,prcntdiss,totatmpercentdiss

;Calculates energy of max dissociation and 95% dissociation
  
  nei=(size(N2_exctrate))[1]
  jmax=(size(N2_exctrate))[2]
  nbins=(size(N2_exctrate))[3]
  lmax=(size(N2_exctrate))[4]




  N2_rate_1 = fltarr(nei,jmax,nbins)
  N2_rate_2 = fltarr(jmax,nbins)

  N2_tot_1 = fltarr(nei,jmax,nbins)
  N2_tot_2 = fltarr(jmax,nbins)


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



          N2_rate_1[s,j,e]=N2_rate_1[s,j,e]+N2_exctrate[s,j,e,l]
          N2_rate_2[j,e]=N2_rate_2[j,e]+N2_exctrate[s,j,e,l]

          N2_tot_1[s,j,e]=N2_tot_1[s,j,e]+N2_exctrate[s,j,e,l]*ener[e]
          N2_tot_2[j,e]=N2_tot_2[j,e]+N2_exctrate[s,j,e,l]*ener[e]






        endfor
      endfor
    endfor

  endfor

totatmpercentdiss=N2_rate_2/total(N2_rate_2)

  for j=0,jmax-1 do begin  ;altitude

    mi1=max(reform(N2_rate_2[j,*]),ind1)
    maxener[j]=ener[ind1]

    N2_tot=fltarr(nbins)
    N2_tot[0]=N2_tot_2[j,0]

    rateratio=reform(N2_rate_2[j,*])/total(reform(N2_rate_2[j,*]));N2_tot/reform(total(N2_tot))
    prcntdiss[j,*]=rateratio
    
    cumprcntdiss[j,0]=rateratio[0]
    
    for e=1,nbins-1 do begin
      
      N2_tot[e]=total(N2_tot_2[j,0:e],2)
      cumprcntdiss[j,e]= total(rateratio[0:e])
    endfor
    
    jj=where(  ( (total(N2_tot_2[j,*])*0.95)-N2_tot )  ge 0,/null)
    ;         print,ener[jj[-1]]
    ;         stop
    sdener[j]=ener[jj[-1]]

  endfor



  for j=0,jmax-1 do begin  ;altitude

    for s=0,nei-1 do begin  ;states
      mi2=max(reform(N2_rate_1[s,j,*]),ind2)
      maxener_states[s,j]=ener[ind2]

      N2_tot=fltarr(nbins)
      N2_tot[0]=N2_tot_1[s,j,0]

      for e=1,nbins-1 do N2_tot[e]=total(N2_tot_1[s,j,0:e],3)

      jj=where(  ( (total(N2_tot_1[s,j,*])*0.95)-N2_tot )  ge 0,/null)


      sdener_states[s,j]=ener[jj[-1]]



    endfor
  endfor




end

;____________________________________________________________________________________________________________________________________________________________________


sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"

fname=[sav_loc+'N2exctrate1_Jul2021_FLAREOLD.sav',$        ;0 flare,old 
       sav_loc+'N2exctrate1_Jul2021_FLARENEW.sav',$             ;1 flare,new

       sav_loc+'N2exctrate1_Jul2021_QUIETOLD.sav',$             ;2 quiet,old
       sav_loc+'N2exctrate1_Jul2021_QUIETNEW.sav']             ;3 quiet,new

  
; file has:  
;  N2_exctrate1_transp,N2_ionzrate1_transp,N2_exctrate1_local,N2_ionzrate1_local,
;  tflux_all,pespec_all,zz,ener,del,wv1,wv2,sza,lat,lon,tau,
;  sigix,sigex
 

plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"

;____________________________________________________________________________________________________________________________________________________________________

restore,fname[0]
newstates_OLD,N2_exctrate1_transp,3,N2_ionzrate1_transp,3,N2_dissext_FOT,N2_dissionz_FOT
newstates_OLD,N2_exctrate1_local,3,N2_ionzrate1_local,3,N2_dissext_FOL,N2_dissionz_FOL
tflux_all_FO=tflux_all
pespec_all_FO=pespec_all
tflux_FO=reform(total(tflux_all,3))
pespec_FO=reform(total(pespec_all,3))
sigex_OLD=sigex
totdisscrss_OLD= sigex_OLD[4,2,*]+sigex_OLD[5,2,*]+sigex_OLD[6,2,*]



;disscontrb,N2_dissext_FOT,N2_dissionz_FOT,N2_exctrate_1_FOT,N2_exctrate_2_FOT,N2_exctrate_3_FOT,$
;          N2_ionzrate_1_FOT,N2_ionzrate_2_FOT,N2_ionzrate_3_FOT,N2_prcnttotdisionz_FOT,N2_prcnttotdisexct_FOT
;
;
;disscontrb,N2_dissext_FOL,N2_dissionz_FOL,N2_exctrate_1_FOL,N2_exctrate_2_FOL,N2_exctrate_3_FOL,$
;          N2_ionzrate_1_FOL,N2_ionzrate_2_FOL,N2_ionzrate_3_FOL,N2_prcnttotdisionz_FOL,N2_prcnttotdisexct_FOL
          
;average_ener,N2_dissext_FOT,ener,avg_ener_3_FOT,avg_ener_2_FOT,avg_ener_1_FOT
;average_ener,N2_dissext_FOL,ener,avg_ener_3_FOL,avg_ener_2_FOL,avg_ener_1_FOL

max_ener,N2_dissext_FOT,ener,maxener_FOT,maxener_states_FOT,sdener_FOT,sdener_states_FOT,N2_rate_1_FOT,N2_rate_2_FOT,cumprcntdiss_FOT,prcntdiss_FOT,totatmpercentdiss_FOT
max_ener,N2_dissext_FOL,ener,maxener_FOL,maxener_states_FOL,sdener_FOL,sdener_states_FOL,N2_rate_1_FOL,N2_rate_2_FOL,cumprcntdiss_FOL,prcntdiss_FOL,totatmpercentdiss_FOL

restore,fname[1]
newstates_NEW,N2_exctrate1_transp,11,N2_ionzrate1_transp,3,N2_dissext_FNT,N2_dissionz_FNT
newstates_NEW,N2_exctrate1_local,11,N2_ionzrate1_local,3,N2_dissext_FNL,N2_dissionz_FNL
tflux_all_FN=tflux_all
pespec_all_FN=pespec_all
tflux_FN=reform(total(tflux_all,3))
pespec_FN=reform(total(pespec_all,3))

totdisscrss_NEW= 0.12*sigex[5,2,*]+0.5*sigex[6,2,*]+0.95*sigex[9,2,*]+0.10*sigex[10,2,*]+0.84*sigex[11,2,*]+sigex[12,2,*]+sigex[13,2,*]+sigex[14,2,*]+$
  0.99*sigex[16,2,*]+sigex[17,2,*]+sigex[18,2,*]

;disscontrb,N2_dissext_FNT,N2_dissionz_FNT,N2_exctrate_1_FNT,N2_exctrate_2_FNT,N2_exctrate_3_FNT,$
;  N2_ionzrate_1_FNT,N2_ionzrate_2_FNT,N2_ionzrate_3_FNT,N2_prcnttotdisionz_FNT,N2_prcnttotdisexct_FNT
;
;
;disscontrb,N2_dissext_FNL,N2_dissionz_FNL,N2_exctrate_1_FNL,N2_exctrate_2_FNL,N2_exctrate_3_FNL,$
;  N2_ionzrate_1_FNL,N2_ionzrate_2_FNL,N2_ionzrate_3_FNL,N2_prcnttotdisionz_FNL,N2_prcnttotdisexct_FNL


;average_ener,N2_dissext_FNT,ener,avg_ener_3_FNT,avg_ener_2_FNT,avg_ener_1_FNT
;average_ener,N2_dissext_FNL,ener,avg_ener_3_FNL,avg_ener_2_FNL,avg_ener_1_FNL

max_ener,N2_dissext_FNT,ener,maxener_FNT,maxener_states_FNT,sdener_FNT,sdener_states_FNT,N2_rate_1_FNT,N2_rate_2_FNT,cumprcntdiss_FNT,prcntdiss_FNT,totatmpercentdiss_FNT
max_ener,N2_dissext_FNL,ener,maxener_FNL,maxener_states_FNL,sdener_FNL,sdener_states_FNL,N2_rate_1_FNL,N2_rate_2_FNL,cumprcntdiss_FNL,prcntdiss_FNL,totatmpercentdiss_FNL


restore,fname[2]
newstates_OLD,N2_exctrate1_transp,3,N2_ionzrate1_transp,3,N2_dissext_QOT,N2_dissionz_QOT
newstates_OLD,N2_exctrate1_local,3,N2_ionzrate1_local,3,N2_dissext_QOL,N2_dissionz_QOL
tflux_all_QO=tflux_all
pespec_all_QO=pespec_all
tflux_QO=reform(total(tflux_all,3))
pespec_QO=reform(total(pespec_all,3))

;disscontrb,N2_dissext_QOT,N2_dissionz_QOT,N2_exctrate_1_QOT,N2_exctrate_2_QOT,N2_exctrate_3_QOT,$
;  N2_ionzrate_1_QOT,N2_ionzrate_2_QOT,N2_ionzrate_3_QOT,N2_prcnttotdisionz_QOT,N2_prcnttotdisexct_QOT
;
;disscontrb,N2_dissext_QOL,N2_dissionz_QOL,N2_exctrate_1_QOL,N2_exctrate_2_QOL,N2_exctrate_3_QOL,$
;          N2_ionzrate_1_QOL,N2_ionzrate_2_QOL,N2_ionzrate_3_QOL,N2_prcnttotdisionz_QOL,N2_prcnttotdisexct_QOL
          
;average_ener,N2_dissext_QOT,ener,avg_ener_3_QOT,avg_ener_2_QOT,avg_ener_1_QOT
;average_ener,N2_dissext_QOL,ener,avg_ener_3_QOL,avg_ener_2_QOL,avg_ener_1_QOL

max_ener,N2_dissext_QOT,ener,maxener_QOT,maxener_states_QOT,sdener_QOT,sdener_states_QOT,N2_rate_1_QOT,N2_rate_2_QOT,cumprcntdiss_QOT,prcntdiss_QOT,totatmpercentdiss_QOT
max_ener,N2_dissext_QOL,ener,maxener_QOL,maxener_states_QOL,sdener_QOL,sdener_states_QOL,N2_rate_1_QOL,N2_rate_2_QOL,cumprcntdiss_QOL,prcntdiss_QOL,totatmpercentdiss_QOL



restore,fname[3]
newstates_NEW,N2_exctrate1_transp,11,N2_ionzrate1_transp,3,N2_dissext_QNT,N2_dissionz_QNT
newstates_NEW,N2_exctrate1_local,11,N2_ionzrate1_local,3,N2_dissext_QNL,N2_dissionz_QNL
tflux_all_QN=tflux_all
pespec_all_QN=pespec_all
tflux_QN=reform(total(tflux_all,3))
pespec_QN=reform(total(pespec_all,3))


;disscontrb,N2_dissext_QNT,N2_dissionz_QNT,N2_exctrate_1_QNT,N2_exctrate_2_QNT,N2_exctrate_3_QNT,$
;  N2_ionzrate_1_QNT,N2_ionzrate_2_QNT,N2_ionzrate_3_QNT,N2_prcnttotdisionz_QNT,N2_prcnttotdisexct_QNT
;
; disscontrb,N2_dissext_QNL,N2_dissionz_QNL,N2_exctrate_1_QNL,N2_exctrate_2_QNL,N2_exctrate_3_QNL,$
;    N2_ionzrate_1_QNL,N2_ionzrate_2_QNL,N2_ionzrate_3_QNL,N2_prcnttotdisionz_QNL,N2_prcnttotdisexct_QNL

;average_ener,N2_dissext_QNT,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT
;average_ener,N2_dissext_QNL,ener,avg_ener_3_QNL,avg_ener_2_QNL,avg_ener_1_QNL

max_ener,N2_dissext_QNT,ener,maxener_QNT,maxener_states_QNT,sdener_QNT,sdener_states_QNT,N2_rate_1_QNT,N2_rate_2_QNT,cumprcntdiss_QNT,prcntdiss_QNT,totatmpercentdiss_QNT
max_ener,N2_dissext_QNL,ener,maxener_QNL,maxener_states_QNL,sdener_QNL,sdener_states_QNL,N2_rate_1_QNL,N2_rate_2_QNL,cumprcntdiss_QNL,prcntdiss_QNL,totatmpercentdiss_QNL


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
fntsz=22
lgdsz=18
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



leg_nam=['a!U1!N $\Pi$ !Lg!N', 'C!U3!N $\Pi$ !Lu!N','b!U1!N $\Pi$ !Lu!N',"c'!L4!N !U1!N $\Sigma$ !Lu!N !U+!N",$
           "b'!U1!N $\Sigma$ !Lu!N !U+!N",'15.8 eV','17.3 eV','Triplet manifold','c!U1!N $\Pi$ !Lu!N','VUV',$
           'Ryd atoms','Dissociative Ionisation','Total Dissociative excitation']

  ;____________________________________________________________________________________________________________________________________________________________________


  ;N2 Electron Impact Dissociation percentage contribution plot
;
;
; 
;xt='Percentage contribution'
;yt='Altitude (km)'
;
;;#1
;;xr=[0.,100.]
; yr=[min(zz),max(zz)]
;titl='Percentage contribution to dissociation of different N!L2!N states'
;
;
;w1 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(N2_percentdiss[5,*],zz,name=leg_nam[0],$
;         xtitle=xt,ytitle=yt,title=titl,$
;         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
;         xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;         thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;d1b=plot(N2_percentdiss[6,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)
;
;d1c=plot(N2_percentdiss[9,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)
;
;d1d=plot(N2_percentdiss[10,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;
;d1e=plot(N2_percentdiss[11,*],zz,name=leg_nam[4],$
;    thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d1f=plot(N2_percentdiss[12,*],zz,name=leg_nam[5],$
;    thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)
;
;d1g=plot(N2_percentdiss[13,*],zz,name=leg_nam[6],$
;    thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)
;
;d1h=plot(N2_percentdiss[14,*],zz,name=leg_nam[7],$
;    thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)
;    
;d1i=plot(N2_percentdiss[16,*],zz,name=leg_nam[8],$
;    thick=thck,color=colors[8],linestyle=lstyl[0],/overplot)
;
;d1j=plot(N2_percentdiss[17,*],zz,name=leg_nam[9],$
;    thick=thck,color=colors[9],linestyle=lstyl[0],/overplot)
;
;d1k=plot(N2_percentdiss[18,*],zz,name=leg_nam[10],$
;    thick=thck,color=colors[10],linestyle=lstyl[0],/overplot)
;
;d1l=plot(N2_percentdiss[20,*],zz,name=leg_nam[11],$
;    thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;    
;
;d1m=plot(N2_percentdisexct,zz,name=leg_nam[12],$
;      thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;leg1 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f,d1g,d1h,d1i,d1j,d1k,d1l,d1m], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)

;_____________________________________________________________________________________________________

;xt='Percentage contribution'
;yt='Altitude (km)'
;
;xr=[0.,100.]
; yr=[min(zz),max(zz)]
;;titl='Percentage contribution to dissociation of different N!L2!N states'
;;
;leg_nam1=['Excitation-Flare&New','Ionization-Flare&New','Excitation-Flare&Glow','Ionization-Flare&Glow',$
;           'Excitation-Quiet&New','Ionization-Quiet&New','Excitation-Quiet&Glow','Ionization-Quiet&Glow']
;titl='Percentage contribution to dissociation of  N!L2!N (Transport)'
;w2 = WINDOW(DIMENSIONS=[1000,800])
;
;
;f1=plot(N2_prcnttotdisexct_FNT,zz,name=leg_nam1[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)
;
;f2=plot(N2_prcnttotdisionz_FNT,zz,name=leg_nam1[1],$
;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;  
;f3=plot(N2_prcnttotdisexct_FOT,zz,name=leg_nam1[2],$
;    thick=thck,color=colors[0],linestyle=lstyl[1],/overplot)
;
;
;f4=plot(N2_prcnttotdisionz_FOT,zz,name=leg_nam1[3],$
;    thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;    
;    
;q1=plot(N2_prcnttotdisexct_QNT,zz,name=leg_nam1[4],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;    
;q2=plot(N2_prcnttotdisionz_QNT,zz,name=leg_nam1[5],$
;    thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;
;q3=plot(N2_prcnttotdisexct_QOT,zz,name=leg_nam1[6],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;q4=plot(N2_prcnttotdisionz_QOT,zz,name=leg_nam1[7],$
;    thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
;
;leg2 = LEGEND(TARGET=[f1,f2,f3,f4,q1,q2,q3,q4], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;
;
;
;titl='Percentage contribution to dissociation of  N!L2!N (Local)'
;w3 = WINDOW(DIMENSIONS=[1000,800])
;
;
;f1=plot(N2_prcnttotdisexct_FNL,zz,name=leg_nam1[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)
;
;f2=plot(N2_prcnttotdisionz_FNL,zz,name=leg_nam1[1],$
;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;
;f3=plot(N2_prcnttotdisexct_FOL,zz,name=leg_nam1[2],$
;  thick=thck,color=colors[0],linestyle=lstyl[1],/overplot)
;
;
;f4=plot(N2_prcnttotdisionz_FOL,zz,name=leg_nam1[3],$
;  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;
;
;q1=plot(N2_prcnttotdisexct_QNL,zz,name=leg_nam1[4],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;q2=plot(N2_prcnttotdisionz_QNL,zz,name=leg_nam1[5],$
;  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;
;q3=plot(N2_prcnttotdisexct_QOL,zz,name=leg_nam1[6],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;q4=plot(N2_prcnttotdisionz_QOL,zz,name=leg_nam1[7],$
;  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
;
;leg3 = LEGEND(TARGET=[f1,f2,f3,f4,q1,q2,q3,q4], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)





;____________________________________________________________________________________________________________________________________________________________________

restore,sav_loc+'N2cross_secJun2021.sav'
;  file has:sigix,sigex,ener

eii=[where(ener le 10,/null)]
sigex(5,2,eii)=0.


;Calculating percentage contribution to dissociation

totdisscrss= 0.12*sigex[5,2,*]+0.5*sigex[6,2,*]+0.95*sigex[9,2,*]+0.10*sigex[10,2,*]+0.84*sigex[11,2,*]+sigex[12,2,*]+sigex[13,2,*]+sigex[14,2,*]+$
  0.99*sigex[16,2,*]+sigex[17,2,*]+sigex[18,2,*] +sigix[3,2,*]+sigix[4,2,*]+sigix[5,2,*];tot_ionz


percentdiss_a1Pig=100.*(0.12*sigex[5,2,*])/totdisscrss
percentdiss_C3Piu=100.*(0.5*sigex[6,2,*])/totdisscrss
percentdiss_b1Piu=100.*(0.95*sigex[9,2,*])/totdisscrss
percentdiss_cprm1Sigmau=100.*(0.10*sigex[10,2,*])/totdisscrss

percentdiss_bprm1Sigmau=100.*(0.84*sigex[11,2,*])/totdisscrss
percentdiss_15=100.*sigex[12,2,*]/totdisscrss
percentdiss_17=100.*sigex[13,2,*]/totdisscrss
percentdiss_trip=100.*sigex[14,2,*]/totdisscrss

percentdiss_c1Piu=100.*(0.99*sigex[16,2,*])/totdisscrss
percentdiss_vuv=100.*sigex[17,2,*]/totdisscrss
percentdiss_Ryd=100.*sigex[18,2,*]/totdisscrss

percentdiss_ionz=100.*(sigix[3,2,*]+sigix[4,2,*]+sigix[5,2,*])/totdisscrss




percentdiss_a1Pig[where(~finite(percentdiss_a1Pig),/null    )] =0.
percentdiss_C3Piu[where(~finite(percentdiss_C3Piu),/null    )] =0.
percentdiss_b1Piu[where(~finite(percentdiss_b1Piu),/null    )] =0.
percentdiss_cprm1Sigmau[where(~finite(percentdiss_cprm1Sigmau),/null    )] =0.
percentdiss_bprm1Sigmau[where(~finite(percentdiss_bprm1Sigmau),/null    )] =0.
percentdiss_15[where(~finite(percentdiss_15),/null    )] =0.
percentdiss_17[where(~finite(percentdiss_17),/null    )] =0.
percentdiss_trip[where(~finite(percentdiss_trip),/null    )] =0.
percentdiss_c1Piu[where(~finite(percentdiss_c1Piu),/null    )] =0.
percentdiss_vuv[where(~finite(percentdiss_vuv),/null    )] =0.
percentdiss_Ryd[where(~finite(percentdiss_Ryd),/null    )] =0.
percentdiss_ionz[where(~finite(percentdiss_ionz),/null    )] =0.

;____________________________________________________________________________________________________________________________________________________________________

;yt='Percentage contribution'
;xt='Energy (eV)'
;xr=[9., 2.5e4]
;yr=[0.,100.01]
;titl='Percentage contribution of states to dissociation cross-section'
;
;
;
;
;w4 = WINDOW(DIMENSIONS=[1000,800])
;
;c1= plot(ener,percentdiss_a1Pig,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;c2=plot(ener,percentdiss_C3Piu,name=leg_nam[1],$
;  thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)
;
;c3= plot(ener,percentdiss_b1Piu,name=leg_nam[2],$
;  thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)
;
;c4=plot(ener,percentdiss_cprm1Sigmau,name=leg_nam[3],$
;  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;
;c5=plot(ener,percentdiss_bprm1Sigmau,name=leg_nam[4],$
;  thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;c6=plot(ener,percentdiss_15,name=leg_nam[5],$
;  thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)
;
;c7=plot(ener,percentdiss_17,name=leg_nam[6],$
;  thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)
;
;c8=plot(ener,percentdiss_trip,name=leg_nam[7],$
;  thick=thck,color=colors[7],linestyle=lstysigex_OLD=sigex

;
;c9=plot(ener,percentdiss_c1Piu,name=leg_nam[8],$
;  thick=thck,color=colors[8],linestyle=lstyl[0],/overplot)
;
;c10=plot(ener,percentdiss_vuv,name=leg_nam[9],$
;  thick=thck,color=colors[9],linestyle=lstyl[0],/overplot)
;
;c11=plot(ener,percentdiss_Ryd,name=leg_nam[10],$
;  thick=thck,color=colors[10],linestyle=lstyl[0],/overplot)
;
;c12=plot(ener,percentdiss_ionz,name=leg_nam[11],$
;  thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)
;
;
;leg4 = LEGEND(TARGET=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)

;____________________________________________________________________________________________________________________________________________________________________

;leg_nam2=['Flare PE MODEL-with Glow Cross-section','Flare PE Model with New Cross-section',$
;  'Quiet PE MODEL-with Glow Cross-section','Quiet PE Model with New Cross-section']
;
;
;xt='Average Energy (eV)'
;yt='Altitude (km)'
;
;
;xr=[10.,10000.]
;yr=[min(zz),max(zz)]
;titl='Average energy of dissociation of N!L2!N (Transport) '
;
;
;w5 = WINDOW(DIMENSIONS=[1000,800])
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
;leg5 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;titl='Average energy of dissociation of N!L2!N (Local) '
;
;
;w6 = WINDOW(DIMENSIONS=[1000,800])
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
;leg6 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)

;____________________________________________________________________________________________________________________________________________________________________

;leg_nam2=['Flare PE MODEL-with Glow Cross-section','Flare PE Model with New Cross-section',$
;  'Quiet PE MODEL-with Glow Cross-section','Quiet PE Model with New Cross-section']
;
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;xr=[0.1,1.e4]
;;xr=[10.,10000.]
;yr=[min(zz),max(zz)]
;titl='Energy of maximum dissociation of N!L2!N (Transport) '
;
;
;w7 = WINDOW(DIMENSIONS=[1000,800])
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
;leg7 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;titl='Energy of maximum dissociation of N!L2!N (Local) '
;
;
;w8 = WINDOW(DIMENSIONS=[1000,800])
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
;leg8 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;


;____________________________________________________________________________________________________________________________________________________________________

;leg_nam2=['Flare PE MODEL-with Glow Cross-section','Flare PE Model with New Cross-section',$
;  'Quiet PE MODEL-with Glow Cross-section','Quiet PE Model with New Cross-section']
;
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;xr=[10.,1.e5]
;;xr=[10.,10000.]
;yr=[min(zz),max(zz)]
;titl='Energy of 95% dissociation of N!L2!N (Transport) '
;
;
;w9 = WINDOW(DIMENSIONS=[1000,800])
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
;leg9 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;
;
;titl='Energy of 95% dissociation of N!L2!N (Local) '
;
;
;w10 = WINDOW(DIMENSIONS=[1000,800])
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
;leg10 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;


;_____________________________________________________________________________________________________

leg_nam2=['Flare PE MODEL-with Glow Cross-section','Flare PE Model with New Cross-section',$
    'Quiet PE MODEL-with Glow Cross-section','Quiet PE Model with New Cross-section']
  
xt=' Energy (eV)'
yt='Photoelectron Dissociation rate (cm!U-3!N s!U-1!N)'

xr=[1.,1.e4]
yr=[]

htind=150

titl='Dissociation rate of N!L2!N at altitiude '+string(zz[htind])+'km'



;w10a = WINDOW(DIMENSIONS=[1000,1000])
;
;d1a=plot(ener,N2_rate_2_FOT[htind,*],name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
;  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(ener,N2_rate_2_FNT[htind,*],name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(ener,N2_rate_2_QOT[htind,*],name=leg_nam2[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,N2_rate_2_QNT[htind,*],name=leg_nam2[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)




w10a = WINDOW(DIMENSIONS=[1000,1000])

d1a=plot(ener,N2_rate_2_FOL[htind,*],name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

d1b=plot(ener,N2_rate_2_FNL[htind,*],name=leg_nam2[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1c=plot(ener,N2_rate_2_QOL[htind,*],name=leg_nam2[2],$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

d1d=plot(ener,N2_rate_2_QNL[htind,*],name=leg_nam2[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)




;_____________________________________________________________________________________________________


xt='Energy (eV)'
xr=[0.1,1.e4]
leg_nam3=['Photoelectron Flux(GLOW,Flare)','Photoelectron Flux(NEW,Flare)','Photoelectron Flux(GLOW,Quiet)','Photoelectron Flux(NEW,Quiet)',$
  'Cross-section(GLOW)','Cross-section(NEW)']

  
titl='Dissociation of N!L2!N at altitiude '+string(zz[htind])+'km'
w11 = WINDOW(DIMENSIONS=[1000,1000])


plot_margin = [0.15, 0.25, 0.15, 0.15]

;p1 = plot(ener, reform(tflux_FO[htind,*]),name=leg_nam3[0], $
;  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N eV!U-1!N )",title=titl,$
;  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;p1b=plot(ener,reform(tflux_FN[htind,*]),name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;p1c=plot(ener,reform(tflux_QO[htind,*]),name=leg_nam3[2],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;p1d=plot(ener,reform(tflux_QN[htind,*]),name=leg_nam3[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;  
;
;p2 = plot(ener, totdisscrss_OLD,name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,totdisscrss_NEW,name=leg_nam3[5],$
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
;








p1 = plot(ener, reform(pespec_FO[htind,*]),name=leg_nam3[0], $
  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N eV!U-1!N )",title=titl,$
  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)


p1b=plot(ener,reform(pespec_FN[htind,*]),name=leg_nam3[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

p1c=plot(ener,reform(pespec_QO[htind,*]),name=leg_nam3[2],$
  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)

p1d=plot(ener,reform(pespec_QN[htind,*]),name=leg_nam3[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)


p2 = plot(ener, totdisscrss_OLD,name=leg_nam3[4] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)

p2b=plot(ener,totdisscrss_NEW,name=leg_nam3[5],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)



a2 = axis('y',target=p2, $
  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
  tickfont_size=fntsz,$
  textpos=1, $                           ; text faces outward
  tickdir=1, $                           ; ticks face inward
  title='Cross-section (cm!U2!N)')

l11a=legend(target=[p1,p1b,p1c],font_size=lgdsz,transparency=transp)
l11b=legend(target=[p1d,p2,p2b],font_size=lgdsz,transparency=transp)



;_____________________________________________________________________________________________________

thck=4

xt='Energy (eV)'
xr=[10.,1.e4]
leg_nam3=['Percent Diss(GLOW,Flare)','Percent Diss(NEW,Flare)','Percent Diss(GLOW,Quiet)','Percent Diss(NEW,Quiet)',$
  'Cumulative Percent Diss(GLOW,Flare)','Cumulative Percent Diss(NEW,Flare)','Cumulative Percent Diss(GLOW,Quiet)','Cumulative Percent Diss(NEW,Quiet)']


sym=['+','*']
titl='Dissociation of N!L2!N at altitiude '+string(zz[htind])+'km'

w11 = WINDOW(DIMENSIONS=[1000,1000])

;plot_margin = [0.15, 0.25, 0.15, 0.15]
;
;p1 = plot(ener, prcntdiss_FOT[htind,*]*100.,name=leg_nam3[0], $
;  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
;  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;p1b=plot(ener,prcntdiss_FNT[htind,*]*100.,name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;p2 = plot(ener, cumprcntdiss_FOT[htind,*]*100.,name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[6],linestyle=lstyl[1],font_size=fntsz,/current)
;
;p2b=plot(ener,cumprcntdiss_FNT[htind,*]*100.,name=leg_nam3[5],$
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




p1 = plot(ener, prcntdiss_FOL[htind,*]*100.,name=leg_nam3[0], $
  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

p1b=plot(ener,prcntdiss_FNL[htind,*]*100.,name=leg_nam3[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

p2 = plot(ener, cumprcntdiss_FOL[htind,*]*100.,name=leg_nam3[4] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],linestyle=lstyl[1],font_size=fntsz,/current)

p2b=plot(ener,cumprcntdiss_FNL[htind,*]*100.,name=leg_nam3[5],$
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



w12 = WINDOW(DIMENSIONS=[1000,1000])

plot_margin = [0.15, 0.25, 0.15, 0.15]

;p1 = plot(ener, prcntdiss_QOT[htind,*]*100.,name=leg_nam3[2], $
;  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
;  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=2,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p1b=plot(ener,prcntdiss_QNT[htind,*]*100.,name=leg_nam3[3],$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;p2 = plot(ener, cumprcntdiss_QOT[htind,*]*100.,name=leg_nam3[6] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=2,color=colors[6],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,cumprcntdiss_QNT[htind,*]*100.,name=leg_nam3[7],$
;    thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
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
;



p1 = plot(ener, prcntdiss_QOL[htind,*]*100.,name=leg_nam3[2], $
  xtitle=xt,ytitle="Percentage dissociation",title=titl,$
  axis_style=1,yrange=[],ylog=0,xrange=xr,xlog=1,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=2,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)

p1b=plot(ener,prcntdiss_QNL[htind,*]*100.,name=leg_nam3[3],$
  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)

p2 = plot(ener, cumprcntdiss_QOL[htind,*]*100.,name=leg_nam3[6] ,$
  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=0,margin=plot_margin,$
  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=2,color=colors[6],linestyle=lstyl[0],font_size=fntsz,/current)

p2b=plot(ener,cumprcntdiss_QNL[htind,*]*100.,name=leg_nam3[7],$
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


;_________________________________________________________________________________________________________________________________

hind=[17,20,30,45,60,150]

leg_nam4=[string(fix(zz[hind[0]]))+'km',$
  string(fix(zz[hind[1]]))+'km',$
  string(fix(zz[hind[2]]))+'km',$
  string(fix(zz[hind[3]]))+'km',$
  string(fix(zz[hind[4]]))+'km',$
  string(fix(zz[hind[5]]))+'km']

xt=' Energy (eV)'
yt='Percentage dissociation'

xr=[1.,1.e4]
;yr=[0.001,100.]


totatmpercentdiss=totatmpercentdiss_FOT
titl='Percentage Dissociative Excitationof N!L2!N !C(Flare/GLOW Cross-section/Transport)'



w12 = WINDOW(DIMENSIONS=[1000,1000])

d1a=plot(ener,totatmpercentdiss[hind[0],*],name=leg_nam4[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],font_size=fntsz,/current)

d1b=plot(ener,totatmpercentdiss[hind[1],*],name=leg_nam4[1],$
  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)

d1c=plot(ener,totatmpercentdiss[hind[2],*],name=leg_nam4[2],$
  thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)

d1d=plot(ener,totatmpercentdiss[hind[3],*],name=leg_nam4[3],$
  thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)

d1e=plot(ener,totatmpercentdiss[hind[4],*],name=leg_nam4[4],$
  thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)

d1f=plot(ener,totatmpercentdiss[hind[5],*],name=leg_nam4[5],$
  thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)

leg12 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)









end

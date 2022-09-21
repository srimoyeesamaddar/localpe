; N2 plots for presentation 
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


;____________________________________________________________________________________________________________________________________________________________________

pro dissstates_NEW,Sp_exctrate,nei,Sp_ionzrate,nii,Sp_dissext,Sp_dissionz
  ; Counts only the states needed for dissociation excitation and ionisation

  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  Sp_dissext=fltarr(nei,jmax,nbins,lmax)
  Sp_dissionz=fltarr(nii,jmax,nbins,lmax)

  Sp_dissext[0,*,*,*]=  Sp_exctrate[5,*,*,*]
  Sp_dissext[1,*,*,*]=  Sp_exctrate[6,*,*,*]
  Sp_dissext[2,*,*,*]=  Sp_exctrate[9,*,*,*]
  Sp_dissext[3,*,*,*]=  Sp_exctrate[10,*,*,*]
  Sp_dissext[4,*,*,*]=  Sp_exctrate[11,*,*,*]
  Sp_dissext[5,*,*,*]=  Sp_exctrate[12,*,*,*]
  Sp_dissext[6,*,*,*]=  Sp_exctrate[13,*,*,*]
  Sp_dissext[7,*,*,*]=  Sp_exctrate[14,*,*,*]
  Sp_dissext[8,*,*,*]=  Sp_exctrate[16,*,*,*]
  Sp_dissext[9,*,*,*]=  Sp_exctrate[17,*,*,*]
  Sp_dissext[10,*,*,*]= Sp_exctrate[18,*,*,*]



  Sp_dissionz[0,*,*,*]= Sp_ionzrate[1,*,*,*]


end

;____________________________________________________________________________________________________________________________________________________________________


pro dissstates_OLD,Sp_exctrate,nei,Sp_ionzrate,nii,Sp_dissext,Sp_dissionz

  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  Sp_dissext=fltarr(nei,jmax,nbins,lmax)
  Sp_dissionz=fltarr(nii,jmax,nbins,lmax)

  Sp_dissext[0,*,*,*]= Sp_exctrate[4,*,*,*]
  Sp_dissext[1,*,*,*]= Sp_exctrate[5,*,*,*]
  Sp_dissext[2,*,*,*]= Sp_exctrate[6,*,*,*]


  Sp_dissionz[0,*,*,*]= Sp_ionzrate[3,*,*,*]
  Sp_dissionz[1,*,*,*]= Sp_ionzrate[4,*,*,*]
  Sp_dissionz[2,*,*,*]= Sp_ionzrate[5,*,*,*]

end

;____________________________________________________________________________________________________________________________________________________________________


pro totionzrate,Sp_ionzrate,Sp_photoi,eiionz_1,eiionz_2,photoi_1
  nei=(size(Sp_ionzrate))[1]
  jmax=(size(Sp_ionzrate))[2]
  nbins=(size(Sp_ionzrate))[3]
  lmax=(size(Sp_ionzrate))[4]



  
  eiionz_1 = fltarr(jmax)
  eiionz_2 = fltarr(jmax,lmax)
  
  photoi_1=reform(total(Sp_photoi,2))


  for l=0,lmax-1 do begin  ;wavelength bins

    for e=0,nbins-1 do begin
      for j=0,jmax-1 do begin  ;altitude

        for s=0,nei-1 do begin  ;states

          eiionz_1[j]=eiionz_1[j]+Sp_ionzrate[s,j,e,l]
          eiionz_2[j,l]=eiionz_2[j,l]+Sp_ionzrate[s,j,e,l]
          
          



        endfor
      endfor
    endfor

  endfor

  


end



;____________________________________________________________________________________________________________________________________________________________________



pro disscontrb,Sp_exctrate,Sp_ionzrate,Sp_exctrate_1,Sp_exctrate_2,$
  Sp_exctrate_3,Sp_exctrate_5,Sp_ionzrate_1,Sp_ionzrate_2,Sp_ionzrate_3,Sp_ionzrate_5,$
  Sp_prcnttotdisionz,Sp_prcnttotdisexct
  
  nei=(size(Sp_exctrate))[1]
  jmax=(size(Sp_exctrate))[2]
  nbins=(size(Sp_exctrate))[3]
  lmax=(size(Sp_exctrate))[4]

  nii=(size(Sp_ionzrate))[1]



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





  ; Percent contribution to dissociation-New

  Sp_prcnttotdisionz=100*Sp_ionzrate_3/(Sp_exctrate_3+Sp_ionzrate_3)
  Sp_prcnttotdisexct=100*Sp_exctrate_3/(Sp_exctrate_3+Sp_ionzrate_3)

  Sp_prcnttotdisionz[where(~finite(Sp_prcnttotdisionz),/null)]=0.
  Sp_prcnttotdisexct[where(~finite(Sp_prcnttotdisexct),/null)]=0.


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


pro max_ener,Sp_exctrate,Sp_disionzrate,ener,maxener,sdener,$
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

  Sp_etot_1 = fltarr(nei,jmax,nbins)
  Sp_etot_2 = fltarr(jmax,nbins)




;DI


  Sp_irate_2 = fltarr(jmax,nbins)
  Sp_itot_2 = fltarr(jmax,nbins)

  
  
  
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

          Sp_etot_1[s,j,e]=Sp_etot_1[s,j,e]+Sp_exctrate[s,j,e,l]*ener[e]
          Sp_etot_2[j,e]=Sp_etot_2[j,e]+Sp_exctrate[s,j,e,l]*ener[e]

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

          
          Sp_itot_2[j,e]=Sp_itot_2[j,e]+Sp_disionzrate[s,j,e,l]*ener[e]

        endfor
      endfor
    endfor

  endfor
  
;  Total dissociation rate
  Sp_rate_2=Sp_irate_2+Sp_erate_2
  Sp_tot_2=Sp_itot_2+Sp_etot_2
  
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


sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"

fname=[sav_loc+'N2_FLAREOLD_Sept2021.sav',$        ;0 flare,old
  sav_loc+'N2_FLARENEW_Sept2021.sav']       ;1 flare,new

  

; file has:
;  Sp_exctrate1_transp,Sp_ionzrate1_transp,Sp_exctrate1_local,Sp_ionzrate1_local,Sp_photoi
;  tflux_all,pespec_all,zz,ener,del,wv1,wv2,sza,lat,lon,tau,
;  sigix,sigex


;____________________________________________________________________________________________________________________________________________________________________

;Cosby 1993 dissociative excitation cross-section data
cosby_ener=[18.5, 23.5, 28.5, 38.5, 48.5, 73.5, 88.5, 148.5];
cosby_crss=[17.4, 66.5, 80.7, 86.5, 101.7, 95.5, 84.0, 90.2 ]*1.e-18;



plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"
;____________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________

restore,fname[0]
dissstates_OLD,Sp_exctrate1_transp,3,Sp_ionzrate1_transp,3,Sp_dissext_FOT,Sp_dissionz_FOT
;stop
tflux_all_FO=tflux_all
pespec_all_FO=pespec_all
tflux_FO=reform(total(tflux_all,3))
pespec_FO=reform(total(pespec_all,3))
sigex_OLD=sigex
sigix_OLD=sigix
Sp_photoi_FOT=Sp_photoi

;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_FOT,eiionz_2_FOT,photoi_FOT

;Dissociative excitation and ionisation cross-sections-OLD
totdisscrss_OLD= sigex_OLD[4,2,*]+sigex_OLD[5,2,*]+sigex_OLD[6,2,*]
totdissionzcrss_OLD= sigix_OLD[3,2,*]+sigix_OLD[4,2,*]+sigix_OLD[5,2,*]

disscontrb,Sp_dissext_FOT,Sp_dissionz_FOT,Sp_exctrate_1_FOT,Sp_exctrate_2_FOT,Sp_exctrate_3_FOT,Sp_exctrate_5_FOT,$
  Sp_ionzrate_1_FOT,Sp_ionzrate_2_FOT,Sp_ionzrate_3_FOT,Sp_ionzrate_5_FOT,Sp_prcnttotdisionz_FOT,Sp_prcnttotdisexct_FOT



average_ener,Sp_dissext_FOT,ener,avg_ener_3_FOT,avg_ener_2_FOT,avg_ener_1_FOT,Sp_totrate_FOT

max_ener,Sp_dissext_FOT,Sp_dissionz_FOT,ener,maxener_FOT,sdener_FOT,cumprcntdiss_FOT,prcntdiss_FOT,totatmpercentdiss_FOT



;_____________________________________________________________________________________
restore,fname[1]
dissstates_NEW,Sp_exctrate1_transp,11,Sp_ionzrate1_transp,1,Sp_dissext_FNT,Sp_dissionz_FNT

tflux_all_FN=tflux_all
pespec_all_FN=pespec_all
tflux_FN=reform(total(tflux_all,3))
pespec_FN=reform(total(pespec_all,3))
sigex_NEW=sigex
sigix_NEW=sigix
Sp_photoi_FNT=Sp_photoi


eii=[where(ener le 10,/null)] ; for this state dissociation does not occur below 10 eV
sigex_NEW(5,2,eii)=0.

;Dissociatiove excitation and ionisation cross-sections-NEW

totdisscrss_NEW= 0.12*sigex_NEW[5,2,*]+0.5*sigex_NEW[6,2,*]+0.95*sigex_NEW[9,2,*]+0.10*sigex_NEW[10,2,*]+0.84*sigex_NEW[11,2,*]+sigex_NEW[12,2,*]+$
  sigex_NEW[13,2,*]+sigex_NEW[14,2,*]+ 0.99*sigex_NEW[16,2,*]+sigex_NEW[17,2,*]+sigex_NEW[18,2,*]

totdissionzcrss_NEW= sigix_NEW[1,2,*]


;Total ionisation rates
totionzrate,Sp_ionzrate1_transp,Sp_photoi,eiionz_FNT,eiionz_2_FNT,photoi_FNT

disscontrb,Sp_dissext_FNT,Sp_dissionz_FNT,Sp_exctrate_1_FNT,Sp_exctrate_2_FNT,Sp_exctrate_3_FNT,Sp_exctrate_5_FNT,$
  Sp_ionzrate_1_FNT,Sp_ionzrate_2_FNT,Sp_ionzrate_3_FNT,Sp_ionzrate_5_FNT,Sp_prcnttotdisionz_FNT,Sp_prcnttotdisexct_FNT




average_ener,Sp_dissext_FNT,ener,avg_ener_3_FNT,avg_ener_2_FNT,avg_ener_1_FNT,Sp_totrate_FNT

max_ener,Sp_dissext_FNT,Sp_dissionz_FNT,ener,maxener_FNT,sdener_FNT,cumprcntdiss_FNT,prcntdiss_FNT,totatmpercentdiss_FNT



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

colors=['red','dodger blue','spring green','fuchsia','gray','deep pink','firebrick','medium blue','orange','purple','dark khaki','dark slate gray','gold']
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
;Comparisons of dissociation cross-sections of N2
yt='Cross-section (10!U16!N cm!U2!N)'
xt='Energy (eV)'

xr=[5.,5.e4]
yr=[0., 2.]

const=1.e-16
leg_nam1=['GLOW ','New Compilation','Cosby 1993']
titl='Total dissociation cross-section due to !C electron impact excitation of N!L2!N'
w1 = WINDOW(DIMENSIONS=[1000,800])


f1=plot(ener, totdisscrss_OLD/const,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

f2=plot(ener, totdisscrss_NEW/const,name=leg_nam1[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

f3=plot(cosby_ener,cosby_crss/const,name=leg_nam1[2],$
  thick=thck,color=colors[3],linestyle=lstyl[1],symbol='*',sym_thick=2,/overplot)


leg1 = LEGEND(TARGET=[f1,f2,f3], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


leg_nam2=['GLOW ','Tian and Vidal 1998']
titl='Total dissociation cross-section due to !C electron impact ionisation of N!L2!N'
yr=[0., 0.7]


w2 = WINDOW(DIMENSIONS=[1000,800])


f1=plot(ener, totdissionzcrss_OLD/const,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

f2=plot(ener, totdissionzcrss_NEW/const,name=leg_nam2[1],$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)



leg2 = LEGEND(TARGET=[f1,f2], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)




;_____________________________________________________________________________________________________






xt='Percentage contribution'
yt='Altitude (km)'

xr=[0.,100.]
yr=[50.,max(zz)]

;
;leg_nam1=['Excitation-Flare&New','Ionization-Flare&New','Excitation-Flare&Glow','Ionization-Flare&Glow',$
;  'Excitation-Quiet&New','Ionization-Quiet&SNew','Excitation-Quiet&Glow','Ionization-Quiet&Glow']


leg_nam1=['Dissociative Excitation','Dissociative Ionization']

titl='Percentage contribution to dissociation of  N!L2!N ';(Transport)
w20 = WINDOW(DIMENSIONS=[800,800])


f1=plot(Sp_prcnttotdisexct_FNT,zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)

f2=plot(Sp_prcnttotdisionz_FNT,zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

;f3=plot(Sp_prcnttotdisexct_FOT,zz,name=leg_nam1[2],$
;  thick=thck,color=colors[0],linestyle=lstyl[1],/overplot)
;
;
;f4=plot(Sp_prcnttotdisionz_FOT,zz,name=leg_nam1[3],$
;  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)



;ar1 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
;t1 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Excitation', /normal,font_size=fntsz,/overplot)


;ar2 = ARROW([0.7,0.6], [0.5,0.7], COLOR=colors[2], /normal, /overplot)
;t2 = TEXT([0.4, 0.22, 0.1], [0.8, 0.58, 0.2],'Ionization', /normal,font_size=fntsz,/overplot)


leg20a = LEGEND(TARGET=[f1,f2], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)



;____________________________________________________________________________________________________________________________________________________________________
;Plots of cross-sections


leg_nam=['a!U1!N $\Pi$ !Lg!N', 'C!U3!N $\Pi$ !Lu!N','b!U1!N $\Pi$ !Lu!N',"c'!L4!N !U1!N $\Sigma$ !Lu!N !U+!N",$
  "b'!U1!N $\Sigma$ !Lu!N !U+!N",'15.8 eV','17.3 eV','Triplet manifold','c!U1!N $\Pi$ !Lu!N','VUV',$
  'Ryd atoms','Dissociative Ionisation','Total Dissociative excitation']

  yt='Cross-section (10!U-16!N cm!U2!N)'
  xt='Energy (eV)'
  xr=[9., 1.e4]
  yr=[0.,100.01]



  d1a=plot(Sp_photoi_FNT[*,wvlind[1]],zz,name='Photoionisation (Pi)',$
    xtitle=xt,ytitle=yt,title=titl[wvlind[1]],$
    xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
    ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,$
    layout=[2,3,3],/current)



;____________________________________________________________________________________________________________________________________________________________________



;Calculating percentage contribution to dissociation

totdisscrss= 0.12*sigex_NEW[5,2,*]+0.5*sigex_NEW[6,2,*]+0.95*sigex_NEW[9,2,*]+0.10*sigex_NEW[10,2,*]+0.84*sigex_NEW[11,2,*]+$
  sigex_NEW[12,2,*]+sigex_NEW[13,2,*]+sigex_NEW[14,2,*]+0.99*sigex_NEW[16,2,*]+sigex_NEW[17,2,*]+sigex_NEW[18,2,*]+$
  sigix_NEW[1,2,*];tot_ionz


percentdiss_a1Pig=100.*(0.12*sigex_NEW[5,2,*])/totdisscrss
percentdiss_C3Piu=100.*(0.5*sigex_NEW[6,2,*])/totdisscrss
percentdiss_b1Piu=100.*(0.95*sigex_NEW[9,2,*])/totdisscrss
percentdiss_cprm1Sigmau=100.*(0.10*sigex_NEW[10,2,*])/totdisscrss

percentdiss_bprm1Sigmau=100.*(0.84*sigex_NEW[11,2,*])/totdisscrss
percentdiss_15=100.*sigex_NEW[12,2,*]/totdisscrss
percentdiss_17=100.*sigex_NEW[13,2,*]/totdisscrss
percentdiss_trip=100.*sigex_NEW[14,2,*]/totdisscrss

percentdiss_c1Piu=100.*(0.99*sigex_NEW[16,2,*])/totdisscrss
percentdiss_vuv=100.*sigex_NEW[17,2,*]/totdisscrss
percentdiss_Ryd=100.*sigex_NEW[18,2,*]/totdisscrss

percentdiss_ionz=100.*(sigix_NEW[1,2,*])/totdisscrss




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

;______________________________________________________________________________________________

yt='Percentage contribution'
xt='Energy (eV)'
xr=[9., 1.e4]
yr=[0.,100.01]
titl='Percentage contribution of states to dissociation cross-section'


leg_nam=['a!U1!N $\Pi$ !Lg!N', 'C!U3!N $\Pi$ !Lu!N','b!U1!N $\Pi$ !Lu!N',"c'!L4!N !U1!N $\Sigma$ !Lu!N !U+!N",$
  "b'!U1!N $\Sigma$ !Lu!N !U+!N",'15.8 eV','17.3 eV','Triplet manifold','c!U1!N $\Pi$ !Lu!N','VUV',$
  'Ryd atoms','Dissociative Ionisation','Total Dissociative excitation']


w4 = WINDOW(DIMENSIONS=[1000,800])

c1= plot(ener,percentdiss_a1Pig,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)

c2=plot(ener,percentdiss_C3Piu,name=leg_nam[1],$
  thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)

c3= plot(ener,percentdiss_b1Piu,name=leg_nam[2],$
  thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)

c4=plot(ener,percentdiss_cprm1Sigmau,name=leg_nam[3],$
  thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)

c5=plot(ener,percentdiss_bprm1Sigmau,name=leg_nam[4],$
  thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)

c6=plot(ener,percentdiss_15,name=leg_nam[5],$
  thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)

c7=plot(ener,percentdiss_17,name=leg_nam[6],$
  thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)

c8=plot(ener,percentdiss_trip,name=leg_nam[7],$
  thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)


c9=plot(ener,percentdiss_c1Piu,name=leg_nam[8],$
  thick=thck,color=colors[8],linestyle=lstyl[0],/overplot)

c10=plot(ener,percentdiss_vuv,name=leg_nam[9],$
  thick=thck,color=colors[9],linestyle=lstyl[0],/overplot)

c11=plot(ener,percentdiss_Ryd,name=leg_nam[10],$
  thick=thck,color=colors[10],linestyle=lstyl[0],/overplot)

c12=plot(ener,percentdiss_ionz,name=leg_nam[11],$
  thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)


leg4a = LEGEND(TARGET=[c1,c2,c3,c4,c5], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
leg4b = LEGEND(TARGET=[c6,c7,c8,c9,c10,c11,c12], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)

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

;Total rates - Photon Ionisation, Ionisation, DI,DE, Dissoc.
;titl='Total Dissociation rate of N!L2!N
yt=' Altitude (km)'
xt='Rate (cm!U-3!N s!U-1!N)'

w6 = WINDOW(DIMENSIONS=[1000,100])

d1a=plot(photoi_FNT,zz,name='Photoionisation (Pi)',$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.01,1.e5],xlog=1,ylog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)

d1b=plot(eiionz_FNT,zz,name='Electron Impact Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1c=plot(Sp_ionzrate_3_FNT,zz,name='EI Dissociative Ionisation',$
    thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

d1d=plot(Sp_exctrate_3_FNT,zz,name='EI Dissociative Excitation',$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1e=plot(Sp_ionzrate_3_FNT+Sp_exctrate_3_FNT,zz,name='EI Dissociation',$
    thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


xt='Ratio'
e1a=plot(eiionz_FNT/photoi_FNT,zz,name='Ratio of Pe and Pi',$
  xtitle=xt,ytitle=[],title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.1,1000],xlog=1,ylog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,/current)

e1b=plot((Sp_ionzrate_3_FNT+Sp_exctrate_3_FNT)/photoi_FNT,zz,name='Ratio of EI Dissoc and Pi ',$
    thick=thck,color=colors[12],linestyle=lstyl[1],/overplot)

leg10c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)
;____________________________________________________________________________________________________________________________________________________________________



;_____________________________________________________________________________________________________

;Total rates - Photon Ionisation, Ionisation, DI,DE, Dissoc.
;titl='Total Dissociation rate of N!L2!N
yt=' Altitude (km)'
xt='Rate (cm!U-3!N s!U-1!N)'
fntsz= 12
;wvlind=[0,1,2]
;wvlind=[6,7,8]
wvlind=[9,10,11]
titl=['0.5-1 angstrom','1-1.5 angstrom','1.5-2 angstrom',$
         '2-2.5 angstrom','2.5-3 angstrom','3-4 angstrom',$
         '4-5 angstrom','5-6 angstrom','6-8 angstrom',$
         '8-10 angstrom','10-14 angstrom','14-18 angstrom']

;#1
w7 = WINDOW(DIMENSIONS=[1000,1000])

d1a=plot(Sp_photoi_FNT[*,wvlind[0]],zz,name='Photoionisation (Pi)',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[0]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,$
    layout=[2,3,1],/current)

d1b=plot(eiionz_2_FNT[*,wvlind[0]],zz,name='Electron Impact Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1c=plot(Sp_ionzrate_5_FNT[*,wvlind[0]],zz,name='EI Dissociative Ionisation',$
  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

d1d=plot(Sp_exctrate_5_FNT[*,wvlind[0]],zz,name='EI Dissociative Excitation',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1e=plot(Sp_ionzrate_5_FNT[*,wvlind[0]]+Sp_exctrate_5_FNT[*,wvlind[0]],zz,name='EI Dissociation',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


xt='Ratio'
e1a=plot(reform(eiionz_2_FNT[*,wvlind[0]])/reform(Sp_photoi_FNT[*,wvlind[0]]),zz,name='Ratio of Pe and Pi',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[0]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.1,1000],xlog=1,ylog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,$
  layout=[2,3,2],/current)

e1b=plot(reform(Sp_ionzrate_5_FNT[*,wvlind[0]]+Sp_exctrate_5_FNT[*,wvlind[0]])/reform(Sp_photoi_FNT[*,wvlind[0]]),zz,name='Ratio of EI Dissoc and Pi ',$
  thick=thck,color=colors[12],linestyle=lstyl[1],/overplot)

;leg10c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


;#2

d1a=plot(Sp_photoi_FNT[*,wvlind[1]],zz,name='Photoionisation (Pi)',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[1]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,$
  layout=[2,3,3],/current)

d1b=plot(eiionz_2_FNT[*,wvlind[1]],zz,name='Electron Impact Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1c=plot(Sp_ionzrate_5_FNT[*,wvlind[1]],zz,name='EI Dissociative Ionisation',$
  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

d1d=plot(Sp_exctrate_5_FNT[*,wvlind[1]],zz,name='EI Dissociative Excitation',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1e=plot(Sp_ionzrate_5_FNT[*,wvlind[1]]+Sp_exctrate_5_FNT[*,wvlind[1]],zz,name='EI Dissociation',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


xt='Ratio'
e1a=plot(reform(eiionz_2_FNT[*,wvlind[1]])/reform(Sp_photoi_FNT[*,wvlind[1]]),zz,name='Ratio of Pe and Pi',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[1]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.1,1000],xlog=1,ylog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,$
  layout=[2,3,4],/current)
e1b=plot(reform(Sp_ionzrate_5_FNT[*,wvlind[1]]+Sp_exctrate_5_FNT[*,wvlind[1]])/reform(Sp_photoi_FNT[*,wvlind[1]]),zz,name='Ratio of EI Dissoc and Pi ',$
  thick=thck,color=colors[12],linestyle=lstyl[1],/overplot)

;leg10c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)




;#3

d1a=plot(Sp_photoi_FNT[*,wvlind[2]],zz,name='Photoionisation (Pi)',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[2]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[1.e-4,1.e5],xlog=1,ylog=0,$;
  ;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,$
  layout=[2,3,5],/current)

d1b=plot(eiionz_2_FNT[*,wvlind[2]],zz,name='Electron Impact Ionisation (Pe)',$
  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)

d1c=plot(Sp_ionzrate_5_FNT[*,wvlind[2]],zz,name='EI Dissociative Ionisation',$
  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)

d1d=plot(Sp_exctrate_5_FNT[*,wvlind[2]],zz,name='EI Dissociative Excitation',$
  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

d1e=plot(Sp_ionzrate_5_FNT[*,wvlind[2]]+Sp_exctrate_5_FNT[*,wvlind[2]],zz,name='EI Dissociation',$
  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)

leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)


xt='Ratio'
e1a=plot(reform(eiionz_2_FNT[*,wvlind[2]])/reform(Sp_photoi_FNT[*,wvlind[2]]),zz,name='Ratio of Pe and Pi',$
  xtitle=xt,ytitle=yt,title=titl[wvlind[2]],$
  xstyle=1,ystyle=1,yrange=[50.,200.],xrange=[0.1,1000],xlog=1,ylog=0,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,$
  layout=[2,3,6],/current)
e1b=plot(reform(Sp_ionzrate_5_FNT[*,wvlind[2]]+Sp_exctrate_5_FNT[*,wvlind[2]])/reform(Sp_photoi_FNT[*,wvlind[2]]),zz,name='Ratio of EI Dissoc and Pi ',$
  thick=thck,color=colors[12],linestyle=lstyl[1],/overplot)

leg10c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal,font_size=lgdsz,transparency=transp)

;____________________________________________________________________________________________________________________________________________________________________
;Energy of 95% dissociation


xt=' Energy (eV)'
yt='Altitude (km)'

fntsz=15
;xr=[1.,1.e5]
xr=[10.,1.e5]
yr=[min(zz),max(zz)]
titl='Energy of 95% dissociation of N!L2!N '


w11 = WINDOW(DIMENSIONS=[700,700])

d1a=plot(sdener_FNT,zz,$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],linestyle=lstyl[1],font_size=fntsz,/current)




;_____________________________________________________________________________________________________


end

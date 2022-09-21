;pro sample_plot_v6

  ;Plot code for sample_localpe_v6.pro

  sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"
  fname=[sav_loc+'m5_time0.sav',$
    sav_loc+'m5_time316.sav',$
    sav_loc+'x9_time0.sav',$
    sav_loc+'x9_time720.sav']
    
plt_loc="/home/srimoyee/Desktop/nrl_files/NRLFeline_plots/"
  ;____________________________________________________________________________________________________________________________________________________________________
;
;  restore,fname[0]
;
;  O_photoi=fltarr(n_elements(fname),n_elements(zz))
;  O2_photoi=fltarr(n_elements(fname),n_elements(zz))
;  N2_photoi=fltarr(n_elements(fname),n_elements(zz))
;
;  O_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))
;  O2_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))
;  N2_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))
;
;  O_pepi=fltarr(n_elements(fname),n_elements(zz))
;  O2_pepi=fltarr(n_elements(fname),n_elements(zz))
;  N2_pepi=fltarr(n_elements(fname),n_elements(zz))
;
;  spectrum=fltarr(n_elements(fname),n_elements(wv1))
;
;  for i=0, n_elements(fname)-1 do begin
;    O_photoi[i,*]=photoi[0,*]
;    O2_photoi[i,*]=photoi[1,*]
;    N2_photoi[i,*]=photoi[2,*]
;
;    O_eiionz_transp[i,*]=eiionz_transp[0,*]
;    O2_eiionz_transp[i,*]=eiionz_transp[1,*]
;    N2_eiionz_transp[i,*]=eiionz_transp[2,*]
;
;    O_pepi[i,*]=O_eiionz_transp[i,*]/O_photoi[i,*]
;    O2_pepi[i,*]=O2_eiionz_transp[i,*]/O2_photoi[i,*]
;    N2_pepi[i,*]=N2_eiionz_transp[i,*]/N2_photoi[i,*]
;
;    O_pepi[i,where(~finite(reform(O_pepi[i,*])),/null)]=0.
;    O2_pepi[i,where(~finite(reform(O2_pepi[i,*])),/null)]=0.
;    N2_pepi[i,where(~finite(reform(N2_pepi[i,*])),/null)]=0.
;
;    spectrum[i,*]=ssflux
;
;    if i+1 lt n_elements(fname) then restore, fname[i+1]
;
;
;  endfor
;
;  ;____________________________________________________________________________________________________________________________________________________________________

 ;Universal Plot variables
  xg=0 ;gridstyles
  yg=0
  xtl=1.  ;ticklen
  ytl=1.
  xsg=1
  ysg=1
  xstl=1.
  ystl=1.

  colors=['red','dark green','dodger blue','purple']
  lstyl=[0,2,3,4]
  postop=[0.1,0.5,0.45,0.95]
  posbot=[0.1,0.1,0.45,0.45]
  thck=5.
  leg_loc=[0.90, 0.85]
  leg_loc_left=[0.37, 0.87]
  maj=['O ','O!L2!N ','N!L2!N ']
  leg_nam=['M5 time=0','M5 time=316', 'X9 time=0', 'X9 time=720']

  head_tilt='NRL spectrum '

;
;  ;____________________________________________________________________________________________________________________________________________________________________
;
;
;  ;O Electron Impact Ionisation plot
;
;  xr=[1.e-5,1.e5]
;  yr=[min(zz),max(zz)]
;  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[0]+'Electron impact ionization'
;
;
;  w1 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;    
;  O_ei1=plot(O_eiionz_transp[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O_ei2=plot(O_eiionz_transp[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O_ei3=plot(O_eiionz_transp[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O_ei4=plot(O_eiionz_transp[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg1 = LEGEND(TARGET=[O_ei1,O_ei2,O_ei3,O_ei4], POSITION=leg_loc,/normal)
;
;  w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;  ;_________________________________________________________________________
;
;  ;O Photoionisation plot
;
;  xr=[1.e-5,1.e5]
;  yr=[min(zz),max(zz)]
;  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[0]+'Photoionization'
;
;
;  w2 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;  O_pi1=plot(O_photoi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O_pi2=plot(O_photoi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O_pi3=plot(O_photoi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O_pi4=plot(O_photoi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg2 = LEGEND(TARGET=[O_pi1,O_pi2,O_pi3,O_pi4], POSITION=leg_loc,/normal)
;
;  w2.Save, plt_loc+"O_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;  ;_________________________________________________________________________
;
;  ;O Pe/Pi plot
;
;
;  xr=[0.1,1000.]
;  yr=[min(zz),max(zz)]
;  xt='Pe/Pi'
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[0]+'Pe/Pi'
;
;
;  w3 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;  O_pp1=plot(O_pepi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O_pp2=plot(O_pepi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O_pp3=plot(O_pepi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O_pp4=plot(O_pepi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg3 = LEGEND(TARGET=[O_pp1,O_pp2,O_pp3,O_pp4], POSITION=leg_loc,/normal)
;
;  w3.Save, plt_loc+"O_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;  ;____________________________________________________________________________________________________________________________________________________________________
;
;
;
;  ;O2 Electron Impact Ionisation plot
;
;  xr=[1.e-5,1.e5]
;  yr=[min(zz),max(zz)]
;  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[1]+'Electron impact ionization'
;
;
;  w4 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;  O2_ei1=plot(O2_eiionz_transp[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O2_ei2=plot(O2_eiionz_transp[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O2_ei3=plot(O2_eiionz_transp[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O2_ei4=plot(O2_eiionz_transp[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg4 = LEGEND(TARGET=[O2_ei1,O2_ei2,O2_ei3,O2_ei4], POSITION=leg_loc,/normal)
;
;  w4.Save, plt_loc+"O2_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;  ;_________________________________________________________________________
;
;  ;O2 Photoionisation plot
;
;  xr=[1.e-5,1.e5]
;  yr=[min(zz),max(zz)]
;  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[1]+'Photoionization'
;
;
;  w5 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;  O2_pi1=plot(O2_photoi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O2_pi2=plot(O2_photoi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O2_pi3=plot(O2_photoi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O2_pi4=plot(O2_photoi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg2 = LEGEND(TARGET=[O2_pi1,O2_pi2,O2_pi3,O2_pi4], POSITION=leg_loc,/normal)
;
;  w5.Save, plt_loc+"O2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;  ;_________________________________________________________________________
;
;  ;O2 Pe/Pi plot
;
;
;  xr=[0.1,1000.]
;  yr=[min(zz),max(zz)]
;  xt='Pe/Pi'
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[1]+'Pe/Pi'
;
;
;  w6 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;  O2_pp1=plot(O2_pepi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  O2_pp2=plot(O2_pepi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  O2_pp3=plot(O2_pepi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;  O2_pp4=plot(O2_pepi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;  leg3 = LEGEND(TARGET=[O2_pp1,O2_pp2,O2_pp3,O2_pp4], POSITION=leg_loc,/normal)
;
;  w6.Save, plt_loc+"O2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;____________________________________________________________________________________________________________________________________________________________________
;
; 
;;N2 Electron Impact Ionisation plot
;
;  xr=[1.e-5,1.e5]
;  yr=[min(zz),max(zz)]
;  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
;  yt='Altitude (km)'
;  titl=head_tilt+ maj[2]+'Electron impact ionization'
;
;
;w7 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;N2_ei1=plot(N2_eiionz_transp[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;N2_ei2=plot(N2_eiionz_transp[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;N2_ei3=plot(N2_eiionz_transp[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;N2_ei4=plot(N2_eiionz_transp[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;leg7 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4], POSITION=leg_loc,/normal)
;
;w7.Save, plt_loc+"N2_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;;_________________________________________________________________________
;
;;N2 Photoionisation plot
;
;xr=[1.e-5,1.e5]
;yr=[min(zz),max(zz)]
;xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
;yt='Altitude (km)'
;titl=head_tilt+ maj[2]+'Photoionization'
;
;
;w8 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;N2_pi1=plot(N2_photoi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;N2_pi2=plot(N2_photoi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;N2_pi3=plot(N2_photoi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;N2_pi4=plot(N2_photoi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;leg8 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4], POSITION=leg_loc,/normal)
;
;w8.Save, plt_loc+"N2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;_________________________________________________________________________
;
;;N2 Pe/Pi plot
;
;
;xr=[0.1,1000.]
;yr=[min(zz),max(zz)]
;xt='Pe/Pi'
;yt='Altitude (km)'
;titl=head_tilt+ maj[2]+'Pe/Pi'
;
;
;w9 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
;
;N2_pp1=plot(N2_pepi[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;N2_pp2=plot(N2_pepi[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;N2_pp3=plot(N2_pepi[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)
;
;N2_pp4=plot(N2_pepi[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[3],/overplot)
;
;leg9 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4], POSITION=leg_loc,/normal)
;
;w9.Save, plt_loc+"N2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;
;
;
;;Solar flux
;
;xr=[0.5,1.75e3]
;;yr=[1.e4,1.e12]
;yr=[min(spectrum),max(spectrum)]
;xt='Wavelength ($\AA$)'
;yt='Solar Flux (Photons cm!U-2!N s!U-1!N))'
;titl=head_tilt+ 'Solar Flux'
;
;
;w10=window(window_title='Solar Flux',dimension=[800,800])
;
;s1=plot(wv1, spectrum[0,*], name=leg_nam[0],$
;         xtitle=xt,ytitle=yt,title=titl,$
;         xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,$
;         xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;         thick=thck,histogram=1,color=colors[0],linestyle=lstyl[0],/current) 
;
;s2=plot(wv1, spectrum[1,*], name=leg_nam[1],$
;         thick=thck,histogram=1,color=colors[1],linestyle=lstyl[1],/overplot)
;
;s3=plot(wv1, spectrum[2,*], name=leg_nam[2],$
;         thick=thck,histogram=1,color=colors[2],linestyle=lstyl[2],/overplot)
;         
;s4=plot(wv1, spectrum[3,*], name=leg_nam[3],$
;         thick=thck,histogram=1,color=colors[3],linestyle=lstyl[3],/overplot)
;
;a1 = ARROW([1.9,1.9], [max(spectrum),min(spectrum)], COLOR='blue', /DATA, /CURRENT)
;t1 = TEXT(0.27,0.67,'Fe line', /normal,font_size=14, font_color='blue')
;
;leg10 = LEGEND(TARGET=[s1,s2,s3,s4], POSITION=leg_loc_left,/normal)
;
;w10.Save, plt_loc+"solar flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  
;____________________________________________________________________________________________________________________________________________________________________

;__________________________________________________________________________________________________________________________________________________

;plotting the cross-sections used for NRL spectra and SQ'05

SQ_crss=sav_loc+'new_crosssec_SQ05.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
restore, SQ_crss
SQ_wv=spec_wave
SQ_sigi_o=sigi_o
SQ_sigi_o2=sigi_o2
SQ_sigi_n2=sigi_n2
SQ_sigab_o=sigab_o
SQ_sigab_o2=sigab_o2
SQ_sigab_n2=sigab_n2

HF_crss=sav_loc+'new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
restore, HF_crss

HF_wv=spec_wave
HF_sigi_o=sigi_o
HF_sigi_o2=sigi_o2
HF_sigi_n2=sigi_n2
HF_sigab_o=sigab_o
HF_sigab_o2=sigab_o2
HF_sigab_n2=sigab_n2

;Cross-section plots


xr=[min(HF_wv),max(HF_wv)]
yt='Cross-section (cm!U2!N)'
xt='Wavelength ($\AA$)'
xylog=[1,1]
ititl='Ionisation Cross-section'
atitl='Absorption Cross-section'
nam=['SQ 2005','Henk+Fennely']

; O ionisation
w8 = WINDOW(WINDOW_TITLE=maj[0]+ititl,DIMENSIONS=[800,800])

csiO_1=plot(SQ_wv,SQ_sigi_o,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+ititl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-16],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigi_o),max(HF_sigi_o)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csiO_2=plot(HF_wv,HF_sigi_o,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg8 = LEGEND(TARGET=[csiO_1,csiO_2], POSITION=leg_loc_left,/normal)

w8.Save, plt_loc+"O_ioncrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;  O2 ionisation
w9 = WINDOW(WINDOW_TITLE=maj[1]+ititl,DIMENSIONS=[800,800])

csiO2_1=plot(SQ_wv,SQ_sigi_o2,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[1]+ititl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-16],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigi_o2),max(HF_sigi_o2)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csiO2_2=plot(HF_wv,HF_sigi_o2,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg9 = LEGEND(TARGET=[csiO2_1,csiO2_2], POSITION=leg_loc_left,/normal)

w9.Save, plt_loc+"O2_ioncrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;  N2 ionisation
w10 = WINDOW(WINDOW_TITLE=maj[2]+ititl,DIMENSIONS=[800,800])

csiN2_1=plot(SQ_wv,SQ_sigi_n2,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[2]+ititl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-16],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigi_n2),max(HF_sigi_n2)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csiN2_2=plot(HF_wv,HF_sigi_o2,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg10 = LEGEND(TARGET=[csiN2_1,csiN2_2], POSITION=leg_loc_left,/normal)

w10.Save, plt_loc+"N2_ioncrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;  O absorption
w11 = WINDOW(WINDOW_TITLE=maj[0]+atitl,DIMENSIONS=[800,800])

csaO_1=plot(SQ_wv,SQ_sigab_o,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+atitl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-16],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigab_o),max(HF_sigab_o)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csaO_2=plot(HF_wv,HF_sigab_o,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg11 = LEGEND(TARGET=[csaO_1,csaO_2], POSITION=leg_loc_left,/normal)

w11.Save, plt_loc+"O_abcrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;  O2 absorption
w12 = WINDOW(WINDOW_TITLE=maj[1]+atitl,DIMENSIONS=[800,800])

csaO2_1=plot(SQ_wv,SQ_sigab_o2,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[1]+atitl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-16],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigab_o2),max(HF_sigab_o2)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csaO2_2=plot(HF_wv,HF_sigab_o2,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg12 = LEGEND(TARGET=[csaO2_1,csaO2_2], POSITION=leg_loc_left,/normal)

w12.Save, plt_loc+"O2_abcrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;  N2 absorption
w13 = WINDOW(WINDOW_TITLE=maj[2]+atitl,DIMENSIONS=[800,800])

csaN2_1=plot(SQ_wv,SQ_sigab_n2,name=nam[0],$
  xtitle=xt,ytitle=yt,title=maj[2]+atitl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=[1.e-22,1.e-15],xlog=xylog[0],ylog=xylog[1],$;,xrange=xr,yrange=[min(HF_sigab_n2),max(HF_sigab_n2)],
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[3],/current)

csaN2_2=plot(HF_wv,HF_sigab_n2,name=nam[1],$
  thick=thck,color=colors[2],linestyle=lstyl[2],/overplot)

leg13 = LEGEND(TARGET=[csaN2_1,csaN2_2], POSITION=leg_loc_left,/normal)

w13.Save, plt_loc+"N2_abcrss_comp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;__________________________________________________________________________________________________________________________________________________
close,/all

end

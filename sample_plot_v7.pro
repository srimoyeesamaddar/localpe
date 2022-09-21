; This version tests the new NRL spectrum provided in
; newspectra.sav file and the code is sample_localpe_v6.pro


  ;Plot code for sample_localpe_v6.pro

  sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"
  fname=[sav_loc+'peakm5.sav',$
    sav_loc+'peakx9.sav']
    
plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"
 ;____________________________________________________________________________________________________________________________________________________________________

  restore,fname[0]

  O_photoi=fltarr(n_elements(fname),n_elements(zz))
  O2_photoi=fltarr(n_elements(fname),n_elements(zz))
  N2_photoi=fltarr(n_elements(fname),n_elements(zz))

  O_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))
  O2_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))
  N2_eiionz_transp=fltarr(n_elements(fname),n_elements(zz))

  O_pepi=fltarr(n_elements(fname),n_elements(zz))
  O2_pepi=fltarr(n_elements(fname),n_elements(zz))
  N2_pepi=fltarr(n_elements(fname),n_elements(zz))

  spectrum=fltarr(n_elements(fname),n_elements(wv1))

  for i=0, n_elements(fname)-1 do begin
    O_photoi[i,*]=photoi[0,*]
    O2_photoi[i,*]=photoi[1,*]
    N2_photoi[i,*]=photoi[2,*]

    O_eiionz_transp[i,*]=eiionz_transp[0,*]
    O2_eiionz_transp[i,*]=eiionz_transp[1,*]
    N2_eiionz_transp[i,*]=eiionz_transp[2,*]

    O_pepi[i,*]=O_eiionz_transp[i,*]/O_photoi[i,*]
    O2_pepi[i,*]=O2_eiionz_transp[i,*]/O2_photoi[i,*]
    N2_pepi[i,*]=N2_eiionz_transp[i,*]/N2_photoi[i,*]

    O_pepi[i,where(~finite(reform(O_pepi[i,*])),/null)]=0.
    O2_pepi[i,where(~finite(reform(O2_pepi[i,*])),/null)]=0.
    N2_pepi[i,where(~finite(reform(N2_pepi[i,*])),/null)]=0.

    spectrum[i,*]=ssflux

    if i+1 lt n_elements(fname) then restore, fname[i+1]


  endfor

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

  colors=['red','dodger blue','dark green','purple']
  lstyl=[0,2,3,4]
  postop=[0.1,0.5,0.45,0.95]
  posbot=[0.1,0.1,0.45,0.45]
  thck=5.
  leg_loc=[0.90, 0.85]
  leg_loc_left=[0.37, 0.87]
  maj=['O ','O!L2!N ','N!L2!N ']
  leg_nam=['M5 flare', 'X9 flare']

  head_tilt='NRL spectrum '


  ;____________________________________________________________________________________________________________________________________________________________________


  ;O Electron Impact Ionisation plot

  xr=[1.e-5,1.e5]
  yr=[min(zz),max(zz)]
  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
  yt='Altitude (km)'
  titl=head_tilt+ maj[0]+'Electron impact ionization'


  w1 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])
    
  O_ei1=plot(O_eiionz_transp[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O_ei2=plot(O_eiionz_transp[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

  
  leg1 = LEGEND(TARGET=[O_ei1,O_ei2], POSITION=leg_loc,/normal)

  w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  ;_________________________________________________________________________

  ;O Photoionisation plot

  xr=[1.e-5,1.e5]
  yr=[min(zz),max(zz)]
  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
  yt='Altitude (km)'
  titl=head_tilt+ maj[0]+'Photoionization'


  w2 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

  O_pi1=plot(O_photoi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O_pi2=plot(O_photoi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

 leg2 = LEGEND(TARGET=[O_pi1,O_pi2], POSITION=leg_loc,/normal)

  w2.Save, plt_loc+"O_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


  ;_________________________________________________________________________

  ;O Pe/Pi plot


  xr=[0.1,1000.]
  yr=[min(zz),max(zz)]
  xt='Pe/Pi'
  yt='Altitude (km)'
  titl=head_tilt+ maj[0]+'Pe/Pi'


  w3 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

  O_pp1=plot(O_pepi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O_pp2=plot(O_pepi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

  leg3 = LEGEND(TARGET=[O_pp1,O_pp2], POSITION=leg_loc,/normal)

  w3.Save, plt_loc+"O_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


  ;____________________________________________________________________________________________________________________________________________________________________



  ;O2 Electron Impact Ionisation plot

  xr=[1.e-5,1.e5]
  yr=[min(zz),max(zz)]
  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
  yt='Altitude (km)'
  titl=head_tilt+ maj[1]+'Electron impact ionization'


  w4 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

  O2_ei1=plot(O2_eiionz_transp[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O2_ei2=plot(O2_eiionz_transp[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

  leg4 = LEGEND(TARGET=[O2_ei1,O2_ei2], POSITION=leg_loc,/normal)

  w4.Save, plt_loc+"O2_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  ;_________________________________________________________________________

  ;O2 Photoionisation plot

  xr=[1.e-5,1.e5]
  yr=[min(zz),max(zz)]
  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
  yt='Altitude (km)'
  titl=head_tilt+ maj[1]+'Photoionization'


  w5 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

  O2_pi1=plot(O2_photoi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O2_pi2=plot(O2_photoi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

  
  leg2 = LEGEND(TARGET=[O2_pi1,O2_pi2], POSITION=leg_loc,/normal)

  w5.Save, plt_loc+"O2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


  ;_________________________________________________________________________

  ;O2 Pe/Pi plot


  xr=[0.1,1000.]
  yr=[min(zz),max(zz)]
  xt='Pe/Pi'
  yt='Altitude (km)'
  titl=head_tilt+ maj[1]+'Pe/Pi'


  w6 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

  O2_pp1=plot(O2_pepi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

  O2_pp2=plot(O2_pepi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)

  leg3 = LEGEND(TARGET=[O2_pp1,O2_pp2], POSITION=leg_loc,/normal)

  w6.Save, plt_loc+"O2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;____________________________________________________________________________________________________________________________________________________________________

 
;N2 Electron Impact Ionisation plot

  xr=[1.e-5,1.e5]
  yr=[min(zz),max(zz)]
  xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '
  yt='Altitude (km)'
  titl=head_tilt+ maj[2]+'Electron impact ionization'


w7 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

N2_ei1=plot(N2_eiionz_transp[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

N2_ei2=plot(N2_eiionz_transp[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)


leg7 = LEGEND(TARGET=[N2_ei1,N2_ei2], POSITION=leg_loc,/normal)

w7.Save, plt_loc+"N2_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;N2 Photoionisation plot

xr=[1.e-5,1.e5]
yr=[min(zz),max(zz)]
xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
yt='Altitude (km)'
titl=head_tilt+ maj[2]+'Photoionization'


w8 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

N2_pi1=plot(N2_photoi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

N2_pi2=plot(N2_photoi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)


leg8 = LEGEND(TARGET=[N2_pi1,N2_pi2], POSITION=leg_loc,/normal)

w8.Save, plt_loc+"N2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;_________________________________________________________________________

;N2 Pe/Pi plot


xr=[0.1,1000.]
yr=[min(zz),max(zz)]
xt='Pe/Pi'
yt='Altitude (km)'
titl=head_tilt+ maj[2]+'Pe/Pi'


w9 = WINDOW(WINDOW_TITLE=head_tilt,DIMENSIONS=[800,800])

N2_pp1=plot(N2_pepi[0,*],zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[0],linestyle=lstyl[0],/current)

N2_pp2=plot(N2_pepi[1,*],zz,name=leg_nam[1],$
    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)


leg9 = LEGEND(TARGET=[N2_pp1,N2_pp2], POSITION=leg_loc,/normal)

w9.Save, plt_loc+"N2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;____________________________________________________________________________________________________________________________________________________________________




;Solar flux

xr=[0.5,1.75e3]
;yr=[1.e4,1.e12]
yr=[min(spectrum),max(spectrum)]
xt='Wavelength ($\AA$)'
yt='Solar Flux (Photons cm!U-2!N s!U-1!N))'
titl=head_tilt+ 'Solar Flux'


w10=window(window_title='Solar Flux',dimension=[800,800])

s1=plot(wv1, spectrum[0,*], name=leg_nam[0],$
         xtitle=xt,ytitle=yt,title=titl,$
         xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,$
         xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
         thick=thck,histogram=1,color=colors[0],linestyle=lstyl[0],/current) 

s2=plot(wv1, spectrum[1,*], name=leg_nam[1],$
         thick=thck,histogram=1,color=colors[1],linestyle=lstyl[1],/overplot)


a1 = ARROW([1.9,1.9], [max(spectrum),min(spectrum)], COLOR=colors[3], /DATA, /CURRENT)
t1 = TEXT(0.27,0.67,'Fe line', /normal,font_size=14, font_color=colors[3])

leg10 = LEGEND(TARGET=[s1,s2], POSITION=leg_loc_left,/normal)

w10.Save, plt_loc+"solar flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  
;____________________________________________________________________________________________________________________________________________________________________


end

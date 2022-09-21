;Theoritical calculation of photoionisation and verification of the model photoionisation
;June 29, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'
restore,filename=floc+'SQ_peakm5.sav'
;restore,filename=floc+'SQ_peakx9.sav'
SQ_photoi_wv= photoi_wv
SQ_wv1= wv1
SQ_wv2= wv2
SQ_ssflux= ssflux

restore,filename=floc+'peakm5.sav'
;restore, filename=floc+'peakx9.sav'
; file has:photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;;Making SQ solar flux from NRL flux
;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5
;
;wavelength_low=[0.5,4.,8.,18.,32.,70.,155.,224.,290.,320.,540.,650.,650.,798.,798.,798.,913.,$
;  913.,913.,975.,987.,1027.]
;
;wavelength_high=[4.,8.,18.,32.,70.,155.,224.,290.,320.,540.,650.,798.,798.,913.,913.,913.,$
;  975.,975.,975.,987.,1027.,1050.]
;
;peakm5_new= fltarr(n_elements(wavelength_low))
;meanpeakx9_new= fltarr(n_elements(wavelength_low))
;
;peakm5_new[0] = total(peakm5[0:5])
;peakm5_new[1] = total(peakm5[6:8])
;peakm5_new[2] = total(peakm5[9:11])
;peakm5_new[3:-1]=peakm5[12:30]
;
;meanpeakx9_new[0] = total(meanpeakx9[0:5])
;meanpeakx9_new[1] = total(meanpeakx9[6:8])
;meanpeakx9_new[2] = total(meanpeakx9[9:11])
;meanpeakx9_new[3:-1]=meanpeakx9[12:30]
;
;
;meanpeakx9=meanpeakx9_new
;peakm5=peakm5_new
;
;wave_gcm1=wavelength_low 
;wave_gcm2=wavelength_high
;
;save,wave_gcm1,wave_gcm2,meanpeakx9,peakm5, filename= '/home/srimoyee/Desktop/nrl_files/sav_files/SQ_newspectra.sav'

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________



plt_loc="/home/srimoyee/Desktop/nrl_files/jul12_plots/m5_"

;Universal Plot variables
xg=0 ;gridstyles
yg=0
xtl=1.  ;ticklen
ytl=1.
xsg=1
ysg=1
xstl=1.
ystl=1.

colors=['red','dark green','dodger blue','purple',"dark goldenrod","midnight blue","orchid","dim gray"]
lstyl=[0,2,3,4]
fnt_sz=15
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]
thck=3.
ln_thck=5
sym_thck=5
leg_loc=[0.92, 0.88]
leg_loc2=[0.95, 0.60]
leg_loc3=[0.95, 0.26]
leg_loc_left=[0.4, 0.88]
maj=['O ','O!L2!N ','N!L2!N ']
leg_nam=['NRL','SQ']

sym=["*","+","tu"]

;___________________________________________________________________________________________________________________________________________________________________________________


leg_1=['.5-1 $\AA$','1-1.5 $\AA$','1.5-2 $\AA$','2-2.5 $\AA$','2.5-3 $\AA$','3-4 $\AA$s',$
           '4-5 $\AA$','5-6 $\AA$','6-8 $\AA$',$
           '8-10 $\AA$','10-14 $\AA$','14-18 $\AA$']
leg_2=['0.5-4 $\AA$','4-8 $\AA$','8-18 $\AA$']

; solar flux

xr=[0.,1027.]
yr=[min(ssflux),max(ssflux)]
xt='Wavelength ($\AA$)'
yt='Flux (Photons cm!U-2!N s!U-1!N))'
titl='M5 Solar Flux'


w1=window(window_title=titl,dimension=[800,800])

s1=plot(wv1, ssflux, $
  xtitle=xt,ytitle=yt,title=titl,name=leg_nam[0],$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=1,color=colors[0],font_size=fnt_sz,/current)

s2=plot(SQ_wv1, SQ_ssflux,name=leg_nam[1],$
        histogram=1,color=colors[2],thick=thck,/overplot)
  
leg1 = LEGEND(TARGET=[s1,s2], POSITION=leg_loc,font_size=fnt_sz,/normal)
      
       
w1.Save, plt_loc+"flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot



  yr=[min(zz),max(zz)]
  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
  yt='Altitude (km)'
  titl='Photoionization'
 
  
  

;  O Photoionisation plot

  xr=[1.e-9,1.e3]  

  w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

  O_pi1=plot(photoi_wv[*,0,0],zz,name=leg_1[0],$
    xtitle=xt,ytitle=yt,title=maj[0]+titl+leg_2[0],$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,/current)

  
  O_pi2=plot(photoi_wv[*,0,1],zz,name=leg_1[1],$
    thick=thck,color=colors[1],/overplot)

  O_pi3=plot(photoi_wv[*,0,2],zz,name=leg_1[2],$
    thick=thck,color=colors[7],/overplot)

  O_pi4=plot(photoi_wv[*,0,3],zz,name=leg_1[3],$
    thick=thck,color=colors[3],/overplot)

  O_pi5=plot(photoi_wv[*,0,4],zz,name=leg_1[4],$
    thick=thck,color=colors[4],/overplot)

  O_pi6=plot(photoi_wv[*,0,5],zz,name=leg_1[5],$
    thick=thck,color=colors[5],/overplot)

  O_pi7=plot(total(photoi_wv[*,0,0:5],3),zz,name=leg_2[0],$
    thick=thck,color=colors[0],/overplot)
    
  O_pi8=plot(SQ_photoi_wv[*,0,0],zz,name="SQ "+leg_2[0],$
    thick=thck,color=colors[2],/overplot)

  leg2 = LEGEND(TARGET=[O_pi1,O_pi2,O_pi3,O_pi4,O_pi5,O_pi6,O_pi7,O_pi8], POSITION=leg_loc,font_size=fnt_sz,/normal)


  w2.Save, plt_loc+"O_photoi_0.5-4.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

  ;O2 Photoionisation plot



  w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

  O2_pi1=plot(photoi_wv[*,1,0],zz,name=leg_1[0],$
    xtitle=xt,ytitle=yt,title=maj[1]+titl+leg_2[0],$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,/current)

  
  O2_pi2=plot(photoi_wv[*,1,1],zz,name=leg_1[1],$
    thick=thck,color=colors[1],/overplot)

  O2_pi3=plot(photoi_wv[*,1,2],zz,name=leg_1[2],$
    thick=thck,color=colors[7],/overplot)

  O2_pi4=plot(photoi_wv[*,1,3],zz,name=leg_1[3],$
    thick=thck,color=colors[3],/overplot)

  O2_pi5=plot(photoi_wv[*,1,4],zz,name=leg_1[4],$
    thick=thck,color=colors[4],/overplot)

  O2_pi6=plot(photoi_wv[*,1,5],zz,name=leg_1[5],$
    thick=thck,color=colors[5],/overplot)

  O2_pi7=plot(total(photoi_wv[*,1,0:5],3),zz,name=leg_2[0],$
    thick=thck,color=colors[0],/overplot)

  O2_pi8=plot(SQ_photoi_wv[*,1,0],zz,name="SQ "+leg_2[0],$
    thick=thck,color=colors[2],/overplot)
    
  leg3 = LEGEND(TARGET=[O2_pi1,O2_pi2,O2_pi3,O2_pi4,O2_pi5,O2_pi6,O2_pi7,O2_pi8], POSITION=leg_loc,font_size=fnt_sz,/normal)


  w3.Save, plt_loc+"O2_photoi_0.5-4.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;N2 Photoionisation plot



  w4 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

  N2_pi1=plot(photoi_wv[*,2,0],zz,name=leg_1[0],$
    xtitle=xt,ytitle=yt,title=maj[2]+titl+leg_2[0],$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,/current)

  
  N2_pi2=plot(photoi_wv[*,2,1],zz,name=leg_1[1],$
    thick=thck,color=colors[1],/overplot)

  N2_pi3=plot(photoi_wv[*,2,2],zz,name=leg_1[2],$
    thick=thck,color=colors[7],/overplot)

  N2_pi4=plot(photoi_wv[*,2,3],zz,name=leg_1[3],$
    thick=thck,color=colors[3],/overplot)

  N2_pi5=plot(photoi_wv[*,2,4],zz,name=leg_1[4],$
    thick=thck,color=colors[4],/overplot)

  N2_pi6=plot(photoi_wv[*,2,5],zz,name=leg_1[5],$
    thick=thck,color=colors[5],/overplot)

  N2_pi7=plot(total(photoi_wv[*,2,0:5],3),zz,name=leg_2[0],$
    thick=thck,color=colors[0],/overplot)

  N2_pi8=plot(SQ_photoi_wv[*,2,0],zz,name="SQ "+leg_2[0],$
    thick=thck,color=colors[2],/overplot)

  
  leg4 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6,N2_pi7,N2_pi8], POSITION=leg_loc,font_size=fnt_sz,/normal)


  w4.Save, plt_loc+"N2_photoi_0.5-4.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  
  
;___________________________________________________________________________________________________________________________________________________________________________________
  
  


;  O Photoionisation plot


w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O_pi1=plot(photoi_wv[*,0,6],zz,name=leg_1[6],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl+leg_2[1],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


O_pi2=plot(photoi_wv[*,0,7],zz,name=leg_1[7],$
  thick=thck,color=colors[1],/overplot)

O_pi3=plot(photoi_wv[*,0,8],zz,name=leg_1[8],$
  thick=thck,color=colors[7],/overplot)

O_pi4=plot(total(photoi_wv[*,0,6:8],3),zz,name=leg_2[1],$
  thick=thck,color=colors[0],/overplot)

O_pi5=plot(SQ_photoi_wv[*,0,1],zz,name="SQ "+leg_2[1],$
    thick=thck,color=colors[2],/overplot)


leg5 = LEGEND(TARGET=[O_pi1,O_pi2,O_pi3,O_pi4,O_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w5.Save, plt_loc+"O_photoi_4-8.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Photoionisation plot



w6 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O2_pi1=plot(photoi_wv[*,1,6],zz,name=leg_1[6],$
  xtitle=xt,ytitle=yt,title=maj[1]+titl+leg_2[1],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


O2_pi2=plot(photoi_wv[*,1,7],zz,name=leg_1[7],$
  thick=thck,color=colors[1],/overplot)

O2_pi3=plot(photoi_wv[*,1,8],zz,name=leg_1[8],$
  thick=thck,color=colors[7],/overplot)

O2_pi4=plot(total(photoi_wv[*,1,6:8],3),zz,name=leg_2[1],$
  thick=thck,color=colors[0],/overplot)
  
O2_pi5=plot(SQ_photoi_wv[*,1,1],zz,name="SQ "+leg_2[1],$
    thick=thck,color=colors[2],/overplot)

leg6 = LEGEND(TARGET=[O2_pi1,O2_pi2,O2_pi3,O2_pi4,O2_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w6.Save, plt_loc+"O2_photoi_4-8.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;N2 Photoionisation plot


w7 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

N2_pi1=plot(photoi_wv[*,2,6],zz,name=leg_1[6],$
  xtitle=xt,ytitle=yt,title=maj[2]+titl+leg_2[1],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


N2_pi2=plot(photoi_wv[*,2,7],zz,name=leg_1[7],$
  thick=thck,color=colors[1],/overplot)

N2_pi3=plot(photoi_wv[*,2,8],zz,name=leg_1[8],$
  thick=thck,color=colors[7],/overplot)

N2_pi4=plot(total(photoi_wv[*,2,6:8],3),zz,name=leg_2[1],$
  thick=thck,color=colors[0],/overplot)
  
N2_pi5=plot(SQ_photoi_wv[*,2,1],zz,name="SQ "+leg_2[1],$
    thick=thck,color=colors[2],/overplot)

leg7 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w7.Save, plt_loc+"N2_photoi_4-8.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;___________________________________________________________________________________________________________________________________________________________________________________

;  O Photoionisation plot



w8 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O_pi1=plot(photoi_wv[*,0,9],zz,name=leg_1[6],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl+leg_2[2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


O_pi2=plot(photoi_wv[*,0,10],zz,name=leg_1[7],$
  thick=thck,color=colors[1],/overplot)

O_pi3=plot(photoi_wv[*,0,11],zz,name=leg_1[8],$
  thick=thck,color=colors[7],/overplot)

O_pi4=plot(total(photoi_wv[*,0,9:11],3),zz,name=leg_2[2],$
  thick=thck,color=colors[0],/overplot)
  
 O_pi5=plot(SQ_photoi_wv[*,0,2],zz,name="SQ "+leg_2[2],$
    thick=thck,color=colors[2],/overplot)

leg8 = LEGEND(TARGET=[O_pi1,O_pi2,O_pi3,O_pi4,O_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w8.Save, plt_loc+"O_photoi_8-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Photoionisation plot



w9 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O2_pi1=plot(photoi_wv[*,1,9],zz,name=leg_1[9],$
  xtitle=xt,ytitle=yt,title=maj[1]+titl+leg_2[2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


O2_pi2=plot(photoi_wv[*,1,10],zz,name=leg_1[10],$
  thick=thck,color=colors[1],/overplot)

O2_pi3=plot(photoi_wv[*,1,11],zz,name=leg_1[11],$
  thick=thck,color=colors[7],/overplot)

O2_pi4=plot(total(photoi_wv[*,1,9:11],3),zz,name=leg_2[2],$
  thick=thck,color=colors[0],/overplot)

O2_pi5=plot(SQ_photoi_wv[*,1,2],zz,name="SQ "+leg_2[2],$
    thick=thck,color=colors[2],/overplot)

leg9 = LEGEND(TARGET=[O2_pi1,O2_pi2,O2_pi3,O2_pi4,O2_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w9.Save, plt_loc+"O2_photoi_8-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;N2 Photoionisation plot



w10 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

N2_pi1=plot(photoi_wv[*,2,9],zz,name=leg_1[9],$
  xtitle=xt,ytitle=yt,title=maj[2]+titl+leg_2[2],$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,/current)


N2_pi2=plot(photoi_wv[*,2,10],zz,name=leg_1[10],$
  thick=thck,color=colors[1],/overplot)

N2_pi3=plot(photoi_wv[*,2,11],zz,name=leg_1[11],$
  thick=thck,color=colors[7],/overplot)

N2_pi4=plot(total(photoi_wv[*,2,9:11],3),zz,name=leg_2[2],$
  thick=thck,color=colors[0],/overplot)

N2_pi5=plot(SQ_photoi_wv[*,2,2],zz,name="SQ "+leg_2[2],$
    thick=thck,color=colors[2],/overplot)

leg10 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


w10.Save, plt_loc+"N2_photoi_8-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


end
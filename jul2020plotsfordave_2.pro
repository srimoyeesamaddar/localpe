;Theoritical calculation of photoionisation and verification of the model photoionisation
;June 29, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'

restore,filename=floc+'SQ_peakm5.sav'
; file has:photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2
SQ_m5_photoi_wv= photoi_wv
SQ_m5_wv1= wv1
SQ_m5_wv2= wv2
SQ_m5_ssflux= ssflux

restore,filename=floc+'SQ_peakx9.sav'
SQ_x9_photoi_wv= photoi_wv
SQ_x9_wv1= wv1
SQ_x9_wv2= wv2
SQ_x9_ssflux= ssflux


restore,filename=floc+'peakm5.sav'
m5_photoi_wv= photoi_wv
m5_wv1= wv1
m5_wv2= wv2
m5_ssflux= ssflux

restore, filename=floc+'peakx9.sav'
x9_photoi_wv= photoi_wv
x9_wv1= wv1
x9_wv2= wv2
x9_ssflux= ssflux




;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________



plt_loc="/home/srimoyee/Desktop/nrl_files/jul12_plots/x9_"

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
posleft=[0.05,0.1,0.45,0.9]
posright=[0.55,0.1,0.95,0.9]
thck=3.
ln_thck=5
sym_thck=5
leg_loc=[0.92, 0.88]

leg_loc_left=[0.05, 0.9]
leg_loc_right=[0.9, 0.9]

maj=['O ','O!L2!N ','N!L2!N ']
leg_nam=['NRL (0.5-18 $\AA$)','SQ (0.5-18 $\AA$)']

sym=["*","+","tu"]

;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot

tit_fl=['M5 flare -','X9 flare -']

yr=[min(zz),max(zz)]
xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
yt='Altitude (km)'
titl='Photoionization'




;  O Photoionisation plot

xr=[1.e-9,1.e4]

w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O_pi1=plot(total(m5_photoi_wv[*,0,0:11],3),zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,position=posleft,/current)



O_pi2=plot(total(SQ_m5_photoi_wv[*,0,0:2],3),zz,name=leg_nam[1],$
  thick=thck,color=colors[0],/overplot)
  
leg1 = LEGEND(TARGET=[O_pi1,O_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
  
  
O_pi3=plot(total(x9_photoi_wv[*,0,0:11],3),zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,position=posright,/current)


O_pi4=plot(total(SQ_x9_photoi_wv[*,0,0:2],3),zz,name=leg_nam[1],$
    thick=thck,color=colors[0],/overplot)

leg2 = LEGEND(TARGET=[O_pi3,O_pi4], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Photoionisation plot



w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O2_pi1=plot(total(m5_photoi_wv[*,1,0:11],3),zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,position=posleft,/current)



O2_pi2=plot(total(SQ_m5_photoi_wv[*,1,0:2],3),zz,name=leg_nam[1],$
  thick=thck,color=colors[0],/overplot)
  
leg3 = LEGEND(TARGET=[O2_pi1,O2_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
  
  
O2_pi3=plot(total(x9_photoi_wv[*,1,0:11],3),zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,position=posright,/current)


O2_pi4=plot(total(SQ_x9_photoi_wv[*,1,0:2],3),zz,name=leg_nam[1],$
    thick=thck,color=colors[0],/overplot)

leg4 = LEGEND(TARGET=[O2_pi3,O2_pi4], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w2.Save, plt_loc+"O2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;N2 Photoionisation plot


w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

N2_pi1=plot(total(m5_photoi_wv[*,2,0:11],3),zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[6],font_size=fnt_sz,position=posleft,/current)

N2_pi2=plot(total(SQ_m5_photoi_wv[*,2,0:2],3),zz,name=leg_nam[1],$
  thick=thck,color=colors[0],/overplot)
  
leg5 = LEGEND(TARGET=[N2_pi1,N2_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
  
  
N2_pi3=plot(total(x9_photoi_wv[*,2,0:11],3),zz,name=leg_nam[0],$
    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
    thick=thck,color=colors[6],font_size=fnt_sz,position=posright,/current)


N2_pi4=plot(total(SQ_x9_photoi_wv[*,2,0:2],3),zz,name=leg_nam[1],$
    thick=thck,color=colors[0],/overplot)

leg6 = LEGEND(TARGET=[N2_pi3,N2_pi4], POSITION=leg_loc_right,font_size=fnt_sz,/normal)

;w3.Save, plt_loc+"N2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;___________________________________________________________________________________________________________________________________________________________________________________








end
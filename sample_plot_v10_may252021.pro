; Testing the dissociation rates due to electron impact


;Plot code for sample_localpe_v6.pro

sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"
fname=[sav_loc+'N2predissrate_May2021.sav',$
  sav_loc+'N2predissrateOLD_May2021.sav']

; file has:
;eionz_wv_transp ,photi_wv ,N2_wv_prediss,zz,wv1,wv2,sza,lat,lon,tau

plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"
;____________________________________________________________________________________________________________________________________________________________________

restore,fname[0]
N2_wv_prediss_NEW=N2_wv_prediss
;N2_wv_prediss_NEW=N2_wv_prediss
;N2_wv_prediss_NEW=N2_wv_prediss


restore,fname[1]
N2_wv_prediss_OLD=N2_wv_prediss

;____________________________________________________________________________________________________________________________________________________________________
; percentchange=100*(N2_wv_prediss_OLD-N2_wv_prediss_NEW)/N2_wv_prediss_OLD

percentchange=dblarr(161,19)
for i=0,18 do begin
;  percentchange[*,i]=100*(N2_wv_prediss_OLD[*,i]-N2_wv_prediss_NEW[*,i])/N2_wv_prediss_OLD[*,i]
  percentchange[*,i]=(N2_wv_prediss_OLD[*,i])/N2_wv_prediss_NEW[*,i]
endfor
percentchange[where(~finite(percentchange),/null)]=0.

;totpercentchange=100*(reform(total(N2_wv_prediss_OLD[*,0:18],2))-reform(total(N2_wv_prediss_NEW[*,0:18],2)))/reform(total(N2_wv_prediss_OLD[*,0:18],2))
totpercentchange=(reform(total(N2_wv_prediss_OLD[*,0:18],2)))/reform(total(N2_wv_prediss_NEW[*,0:18],2))
totpercentchange[where(~finite(totpercentchange),/null)]=0.

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

postleft=[0.1,0.1,0.5,0.90]
postright=[0.55,0.1,0.95,0.90]

thck=5.
leg_loc=[0.47, 0.85]
leg_loc_left=[0.37, 0.87]
leg_transp=50
maj=['O ','O!L2!N ','N!L2!N ']

leg_nam=['With GLOW cross-section', 'With New compilation']

head_tilt=['0.5-1 $\AA$', '1-1.5 $\AA$','1.5-2 $\AA$','2-2.5 $\AA$','2.5-3 $\AA$','3- 4 $\AA$',$
           '4- 5 $\AA$','5-6 $\AA$','6-8 $\AA$','8-10 $\AA$','10-14 $\AA$','14-18 $\AA$']
;____________________________________________________________________________________________________________________________________________________________________


;N2 Electron Impact Dissociation plot


yr=[50.,max(zz)]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
yt='Altitude (km)'

;#1
xr=[1.e-5,1.e2]
titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[0]


w1 = WINDOW(DIMENSIONS=[1000,800])

d1a=plot(N2_wv_prediss_OLD[*,0],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=[min(zz),150],xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d1b=plot(N2_wv_prediss_NEW[*,0],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg1 = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[0]
d1c=plot(percentchange[*,0],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=[min(zz),150],xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)



;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#2
xr=[1.e-5,1.e5]
titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[1]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w2 = WINDOW(DIMENSIONS=[1000,800])

d2a=plot(N2_wv_prediss_OLD[*,1],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d2b=plot(N2_wv_prediss_NEW[*,1],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg2 = LEGEND(TARGET=[d2a,d2b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[1]
d2c=plot(percentchange[*,1],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________


;#3

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[2]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w3 = WINDOW(DIMENSIONS=[1000,800])

d3a=plot(N2_wv_prediss_OLD[*,2],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d3b=plot(N2_wv_prediss_NEW[*,2],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg3 = LEGEND(TARGET=[d3a,d3b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[2]
d3c=plot(percentchange[*,2],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#4

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[3]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w4 = WINDOW(DIMENSIONS=[1000,800])

d4a=plot(N2_wv_prediss_OLD[*,3],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d4b=plot(N2_wv_prediss_NEW[*,3],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg4 = LEGEND(TARGET=[d4a,d4b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[3]
d4c=plot(percentchange[*,3],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#5

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[4]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w5 = WINDOW(DIMENSIONS=[1000,800])

d5a=plot(N2_wv_prediss_OLD[*,4],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d5b=plot(N2_wv_prediss_NEW[*,4],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg5 = LEGEND(TARGET=[d5a,d5b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[4]
d5c=plot(percentchange[*,4],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#6

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[5]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w6 = WINDOW(DIMENSIONS=[1000,800])

d6a=plot(N2_wv_prediss_OLD[*,5],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d6b=plot(N2_wv_prediss_NEW[*,5],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg6 = LEGEND(TARGET=[d6a,d6b], POSITION=leg_loc,/normal,transparency=leg_transp)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[5]
d6c=plot(percentchange[*,5],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)

;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#7

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[6]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w7 = WINDOW(DIMENSIONS=[1000,800])

d7a=plot(N2_wv_prediss_OLD[*,6],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d7b=plot(N2_wv_prediss_NEW[*,6],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg7 = LEGEND(TARGET=[d7a,d7b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[6]
d7c=plot(percentchange[*,6],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#8

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[7]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w8 = WINDOW(DIMENSIONS=[1000,800])

d8a=plot(N2_wv_prediss_OLD[*,7],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d8b=plot(N2_wv_prediss_NEW[*,7],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg8 = LEGEND(TARGET=[d8a,d8b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[7]
d8c=plot(percentchange[*,7],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________


;#9

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[8]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w9 = WINDOW(DIMENSIONS=[1000,800])

d9a=plot(N2_wv_prediss_OLD[*,8],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d9b=plot(N2_wv_prediss_NEW[*,8],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg9 = LEGEND(TARGET=[d9a,d9b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[8]
d9c=plot(percentchange[*,8],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#10

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[9]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w10 = WINDOW(DIMENSIONS=[1000,800])

d10a=plot(N2_wv_prediss_OLD[*,9],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d10b=plot(N2_wv_prediss_NEW[*,9],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg10 = LEGEND(TARGET=[d10a,d10b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[9]
d10c=plot(percentchange[*,9],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#11

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[10]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w11 = WINDOW(DIMENSIONS=[1000,800])

d11a=plot(N2_wv_prediss_OLD[*,10],zz,name=leg_nam[0],position=postleft,$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d11b=plot(N2_wv_prediss_NEW[*,10],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg11 = LEGEND(TARGET=[d11a,d11b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[10]
d11c=plot(percentchange[*,10],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#12

titl='Electron impact dissociation of molecular N!L2!N !C'+ head_tilt[11]
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w12 = WINDOW(DIMENSIONS=[1000,800])

d12a=plot(N2_wv_prediss_OLD[*,11],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

d12b=plot(N2_wv_prediss_NEW[*,11],zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


leg12 = LEGEND(TARGET=[d12a,d12b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'+ head_tilt[11]
d12c=plot(percentchange[*,11],zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)


;w1.Save, plt_loc+"O_eiionz_transp.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_________________________________________________________________________

;#Final

xr=[1.e0,1.e5]
titl='Electron impact dissociation of molecular N!L2!N !C'
xt='Electron impact dissociation (Pe)!C !C(cm!U-3!N s!U-1!N) '
w12 = WINDOW(DIMENSIONS=[1000,800])

df1a=plot(total(N2_wv_prediss_OLD,2),zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],linestyle=lstyl[0],/current)

df1b=plot(total(N2_wv_prediss_NEW,2),zz,name=leg_nam[1],$
  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)


legdf = LEGEND(TARGET=[df1a,df1b], POSITION=leg_loc,/normal)

xt='Ratio change '
titl='Ratio(GLOW/NEW) of rate'
df1c=plot(totpercentchange,zz,$
  xtitle=xt,title=titl,position=postright,$
  xstyle=1,ystyle=1,yrange=yr,xrange=[0,10],xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],linestyle=lstyl[0],/current)




end

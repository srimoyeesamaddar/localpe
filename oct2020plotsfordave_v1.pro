;Theoritical calculation of photoionisation and verification of the model photoionisation-FINAL RESULTS with ols data
;Oct 27, 2020
;checking calculatios for each bin and only low solar zenith angle


;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'



;______________________________________________________________________________________________________


restore,filename=floc+'peakx9_oct2020.sav'

; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau


nrl_wv1= wv1
nrl_wv2= wv2

x9_tau=tau   ;[altitude, wavelength]

x9_szalo_photoi= photi_wv
x9_szalo_eiionz_loc= eionz_wv_loc
x9_szalo_eiionz_transp= eionz_wv_transp

;______________________________________________________________________________________________________


restore,filename=floc+'peakm5_oct2020.sav'

; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau


m5_tau=tau   ;[altitude, wavelength]

m5_szalo_photoi= photi_wv
m5_szalo_eiionz_loc= eionz_wv_loc
m5_szalo_eiionz_transp= eionz_wv_transp



close,/all
;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________


plt_loc="/home/srimoyee/Desktop/nrl_files/jul12_plots/x9_"

;Universal Plot variables
;xg=0 ;gridstyles
;yg=0
;xtl=1.  ;ticklen
;ytl=1.
;xsg=-1
;ysg=-1
;xstl=0.5
;ystl=0.5

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
thck=3.
ln_thck=5
sym_thck=2
sym_incr=10
leg_loc=[0.1, 0.9]

leg_loc_left=[0.04, 0.9]
leg_loc_right=[0.98, 0.9]
leg_transp=50
maj=['O ','O!L2!N ','N!L2!N ']

sym=["*","+","tu"]


;___________________________________________________________________________________________________________________________________________________________________________________

postpleft=[0.05,0.55,0.33,0.9]
postpmid=[0.39,0.55,0.68,0.9]
postpright=[0.73,0.55,0.98,0.9]

posbtleft=[0.05,0.15,0.33,0.45]
posbtmid=[0.39,0.15,0.68,0.45]
posbtright=[0.73,0.15,0.98,0.45]


leg_nam1=['0.5-1.0 A','1.0-1.5 A','1.5-2.0 A','2.0-2.5 A','2.5-3.0 A','3.0-4.0 A']
leg_nam2=['4-5 A','5-6 A','6-8 A']
leg_nam3=['8-10 A','10-14 A','14-18 A']
leg_nam='0.5-18 A'

tit_fl=['M5 flare -','X9 flare -']


;_______________________________________________________________________________________________________________________

;___________________________________________________N2 plot (1st 6bins)__________________________________________________


;____________________________________________________M5 Flare____________________________________________________________






;___________________________________________________Photoionisation plot__________________________________________________

sp=2
xr=[1.e-5,1.e2]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)
yt='Altitude (km)'

titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(m5_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(m5_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(m5_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[1.,1000.]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,0]/m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,1]/m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,2]/m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(m5_szalo_eiionz_transp[sp,*,3]/m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(m5_szalo_eiionz_transp[sp,*,4]/m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(m5_szalo_eiionz_transp[sp,*,5]/m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(x9_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(x9_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(x9_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[1.,1000.]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,0]/x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,1]/x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,2]/x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(x9_szalo_eiionz_transp[sp,*,3]/x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(x9_szalo_eiionz_transp[sp,*,4]/x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(x9_szalo_eiionz_transp[sp,*,5]/x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________











;___________________________________________________N2 plot (2nd 3bins)____________________________________________________________________________________


;____________________________________________________M5 Flare______________________________________________________________________________________________






;___________________________________________________Photoionisation plot___________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'



w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot____________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________________

xr=[1.,1000.]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,6]/m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,7]/m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,8]/m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare________________________________________________________________________________________________

;___________________________________________________Photoionisation plot_________________________________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionisation '

N2_pi1=plot(x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'

N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,6]/x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,7]/x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,8]/x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________












;___________________________________________________N2 plot (3rd 3bins)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,9]/m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,10]/m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,11]/m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e5]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'


N2_pi1=plot(x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,9]/x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,10]/x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,11]/x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________









;___________________________________________________N2 plot (Total)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e5]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w4 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)


leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[0.1,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3)/total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)


;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e5]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'


N2_pi1=plot(total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)


;leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'

N2_ei1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3)/total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)




;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________












;_______________________________________________________________________________________________________________________

;___________________________________________________O2 plot (1st 6bins)__________________________________________________


;____________________________________________________M5 Flare____________________________________________________________






;___________________________________________________Photoionisation plot__________________________________________________

sp=1
xr=[1.e-4,1.e3]


;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)
yt='Altitude (km)'

titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(m5_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(m5_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(m5_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,0]/m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,1]/m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,2]/m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(m5_szalo_eiionz_transp[sp,*,3]/m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(m5_szalo_eiionz_transp[sp,*,4]/m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(m5_szalo_eiionz_transp[sp,*,5]/m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________


xr=[1.e-4,1.e3]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'

N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(x9_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(x9_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(x9_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,0]/x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,1]/x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,2]/x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(x9_szalo_eiionz_transp[sp,*,3]/x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(x9_szalo_eiionz_transp[sp,*,4]/x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(x9_szalo_eiionz_transp[sp,*,5]/x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________











;___________________________________________________O2 plot (2nd 3bins)____________________________________________________________________________________


;____________________________________________________M5 Flare______________________________________________________________________________________________






;___________________________________________________Photoionisation plot___________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot____________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,6]/m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,7]/m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,8]/m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare________________________________________________________________________________________________

;___________________________________________________Photoionisation plot_________________________________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelctron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,6]/x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,7]/x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,8]/x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________












;___________________________________________________O2 plot (3rd 3bins)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,9]/m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,10]/m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,11]/m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,9]/x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,10]/x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,11]/x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________









;___________________________________________________O2 plot (Total)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e5]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w4 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)


leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[0.1,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3)/total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)


;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e5]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'


N2_pi1=plot(total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)


;leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl= 'Photoelectron Ionisation'


N2_ei1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3)/total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)




;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________













;_______________________________________________________________________________________________________________________

;___________________________________________________O plot (1st 6bins)__________________________________________________


;____________________________________________________M5 Flare____________________________________________________________






;___________________________________________________Photoionisation plot__________________________________________________

sp=0
xr=[1.e-9,1.e1]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)


titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(m5_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(m5_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(m5_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,0]/m5_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,1]/m5_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,2]/m5_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(m5_szalo_eiionz_transp[sp,*,3]/m5_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(m5_szalo_eiionz_transp[sp,*,4]/m5_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(m5_szalo_eiionz_transp[sp,*,5]/m5_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________

xr=[1.e-9,1.e1]

;xr=[1.e-4,1.e1]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photionisation'


N2_pi1=plot(x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(x9_szalo_eiionz_transp[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(x9_szalo_eiionz_transp[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(x9_szalo_eiionz_transp[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,0]/x9_szalo_photoi[sp,*,0],zz,name=leg_nam1[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,1]/x9_szalo_photoi[sp,*,1],zz,name=leg_nam1[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,2]/x9_szalo_photoi[sp,*,2],zz,name=leg_nam1[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(x9_szalo_eiionz_transp[sp,*,3]/x9_szalo_photoi[sp,*,3],zz,name=leg_nam1[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(x9_szalo_eiionz_transp[sp,*,4]/x9_szalo_photoi[sp,*,4],zz,name=leg_nam1[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(x9_szalo_eiionz_transp[sp,*,5]/x9_szalo_photoi[sp,*,5],zz,name=leg_nam1[5],$
  thick=thck,color=colors[7],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________











;___________________________________________________O plot (2nd 3bins)____________________________________________________________________________________


;____________________________________________________M5 Flare______________________________________________________________________________________________






;___________________________________________________Photoionisation plot___________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot____________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,6]/m5_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,7]/m5_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,8]/m5_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare________________________________________________________________________________________________

;___________________________________________________Photoionisation plot_________________________________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'


N2_pi1=plot(x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,6]/x9_szalo_photoi[sp,*,6],zz,name=leg_nam2[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,7]/x9_szalo_photoi[sp,*,7],zz,name=leg_nam2[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,8]/x9_szalo_photoi[sp,*,8],zz,name=leg_nam2[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________












;___________________________________________________O plot (3rd 3bins)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)

N2_pi2=plot(m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

N2_ei2=plot(m5_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(m5_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)





;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_szalo_eiionz_transp[sp,*,9]/m5_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)

N2_pp2=plot(m5_szalo_eiionz_transp[sp,*,10]/m5_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(m5_szalo_eiionz_transp[sp,*,11]/m5_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photionisation'


N2_pi1=plot(x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)

N2_pi2=plot(x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)




;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_szalo_eiionz_transp[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_szalo_eiionz_transp[sp,*,9]/x9_szalo_photoi[sp,*,9],zz,name=leg_nam3[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[sp,*,10]/x9_szalo_photoi[sp,*,10],zz,name=leg_nam3[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[sp,*,11]/x9_szalo_photoi[sp,*,11],zz,name=leg_nam3[2],$
  thick=thck,color=colors[4],/overplot)


;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________









;___________________________________________________O plot (Total)_____________________________________________________________________________________________


;____________________________________________________M5 Flare________________________________________________________________________________________________________






;___________________________________________________Photoionisation plot_____________________________________________________________________________________________


xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)

titl='Photoionization'


w4 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)


leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,orientation=1,/normal)


;___________________________________________________Photoelectron Ionisation plot_______________________________________________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)

;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot___________________________________________________________________________________________________________

xr=[1.,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(total(m5_szalo_eiionz_transp[sp,*,0:11],3)/total(m5_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)


;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare_______________________________________________________________________________________________________

;___________________________________________________Photoionisation plot________________________________________________________________________________________________


xr=[1.e-4,1.e4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)


;leg1 = LEGEND(TARGET=[N2_pi1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot_________________________________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'

N2_ei1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot_____________________________________________________________________________________________________________

xr=[1.,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(total(x9_szalo_eiionz_transp[sp,*,0:11],3)/total(x9_szalo_photoi[sp,*,0:11],3),zz,name=leg_nam,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)




;leg3 = LEGEND(TARGET=[N2_pp1], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_______________________________________________________________________________________________________________________________________________________________________
;________________________________________________________________________________________________________________________________________________________________________











end






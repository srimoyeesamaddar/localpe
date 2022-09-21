


floc='/home/srimoyee/Desktop/nrl_files/sav_files/'

file_m5= floc+'NRL_m5.sav'
restore, file_m5
;file has:
; photoi, eiionz_transp, zz, ssflux, wv_lo,wv_hi,sza_deg

m5_photoi= photoi
m5_eiionz_transp= eiionz_transp

file_x9= floc+'NRL_x9.sav'
restore, file_x9
x9_photoi= photoi
x9_eiionz_transp= eiionz_transp


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



tit_fl=['M5 flare -','X9 flare -']


;_______________________________________________________________________________________________________________________

;___________________________________________________N2 plot ____________________________________________________________


;____________________________________________________M5 Flare____________________________________________________________






;___________________________________________________Photoionisation plot__________________________________________________

sp=2
xr=[1.e-4,1.e5]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)
yt='Altitude (km)'

titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=postpleft,/current)






;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=postpmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[0.01,1000.]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_eiionz_transp[sp,*]/m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=postpright,/current)




;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________


xr=[1.e-4,1.e5]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posbtleft,/current)



;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[0.01,1000.]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_eiionz_transp[sp,*]/x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posbtright,/current)



;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________







;_______________________________________________________________________________________________________________________

;___________________________________________________O2 plot ____________________________________________________________


;____________________________________________________M5 Flare____________________________________________________________






;___________________________________________________Photoionisation plot__________________________________________________

sp=1
xr=[1.e-4,1.e5]


;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)
yt='Altitude (km)'

titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=postpleft,/current)


;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=postpmid,/current)



;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[0.01,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_eiionz_transp[sp,*]/m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=postpright,/current)



;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________


xr=[1.e-4,1.e5]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionisation'

N2_pi1=plot(x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=posbtleft,/current)



;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'

N2_ei1=plot(x9_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[0.01,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_eiionz_transp[sp,*]/x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],font_size=fnt_sz,position=posbtright,/current)



;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________






;_______________________________________________________________________________________________________________________

;___________________________________________________O plot _____________________________________________________________


;____________________________________________________M5 Flare____________________________________________________________





;___________________________________________________Photoionisation plot__________________________________________________

sp=0
xr=[1.e-4,1.e4]

;xr=[]
;yr=[]

xt=[];'Photoionisation (Pi) (cm!U-3!N s!U-1!N)


titl='Photoionization'


w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())


N2_pi1=plot(m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpleft,/current)


;___________________________________________________Photoelectron Ionisation plot___________________________________________________


xt=[];'Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(m5_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot______________________________________________________________________

xr=[0.01,1000]
xt=[];'Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(m5_eiionz_transp[sp,*]/m5_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=postpright,/current)




;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;______________________________________________________________________________________________________________________________________________


;________________________________________________________X9 flare______________________________________________________________________________

;___________________________________________________Photoionisation plot_______________________________________________________________________

xr=[1.e-4,1.e4]

;xr=[1.e-4,1.e1]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photionisation'


N2_pi1=plot(x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtleft,/current)



;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________Photoelectron Ionisation plot__________________________________________________________________


xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionisation'


N2_ei1=plot(x9_eiionz_transp[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtmid,/current)




;leg2 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________Pe/Pi plot________________________________________________________________________________________

xr=[0.01,1000]
xt='Pe/Pi '
titl='Pe/Pi'


N2_pp1=plot(x9_eiionz_transp[sp,*]/x9_photoi[sp,*],zz,$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[sp]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posbtright,/current)



;leg3 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________________________________________________________________________________________
;_________________________________________________________________________________________________________________________________________________________














end












;Theoritical calculation of photoionisation and verification of the model photoionisation
;Aug 24, 2020
;checking calculatios for each bin and only low solar zenith angle


;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'


;______________________________________________________________________________________________________

;
;restore,filename=floc+'SQpeakm5_wv_szalo.sav'
;restore,filename=floc+'SQpeakx9_wv_szalo.sav'
restore,filename=floc+'SQpeakx9_wv_szalo_test.sav'

;restore,filename=floc+'SQpeakx9_wv_szalo_35eV.sav'

; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau

SQx9_szalo_photoi= photi_wv
SQx9_szalo_eiionz_loc= eionz_wv_loc
SQx9_szalo_eiionz_transp= eionz_wv_transp
x9_lat= fix(lat)
x9_lon= fix(lon)
x9_sza=  fix(round(sza*!radeg))

SQ_wv1= wv1
SQ_wv2= wv2
SQ_tau = tau



;
;restore,filename=floc+'SQpeakx9_wv_szalo_v6.sav'  ;with new cross-section bins and states
;
;SQx9_szalo_photoi= photi_wv
;SQx9_szalo_eiionz_loc= eionz_wv_loc
;SQx9_szalo_eiionz_transp= eionz_wv_transp
;
;SQ_wv1= wv1
;SQ_wv2= wv2

;______________________________________________________________________________________________________


;restore,filename=floc+'peakx9_wv_szalo_v8.sav'
;restore,filename=floc+'peakm5_wv_szalo.sav'
restore,filename=floc+'peakx9_ssflux1.sav';'peakx9_wv_szalo_2states.sav';peakx9_wv_szalo_2states.sav

; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau


nrl_wv1= wv1
nrl_wv2= wv2
nrl_tau=tau   ;[altitude, wavelength]

x9_szalo_photoi= photi_wv
x9_szalo_eiionz_loc= eionz_wv_loc
x9_szalo_eiionz_transp= eionz_wv_transp




close,/all
;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________



;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/SQ_newspectra.sav'
;
;;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5
;
;ssflux=meanpeakx9
;
;jmax=n_elements(zz)
;dz=(zz(1:*)-zz(0:-2))
;dz=[dz[0],dz]
;; Energy conservation tests for photoionisation and photoelectrons
;SQ_wv1=nrl_wv1
;SQ_wv2=nrl_wv2
;
;;the height integrated photoionisation energy should add up to the solar flux,
;;atleast for the lower wavlength bins
;
;photoi_ht = fltarr(n_elements(zz),n_elements(SQ_wv1))  ;height integrated photoionisation energy
;
;;photoi_wv:Array[ species,alt, wavelength]
;;photoi_stot=total(x9_szalo_photoi,1)  ;total photoization in all species
;photoi_stot=reform(x9_szalo_photoi[2,*,*])  ;N2 photoization rate
;
;for l=0, n_elements(SQ_wv1)-1 do begin
;  photoi_ht(*,l)=photoi_stot[*,l]*dz*1.e5 ;height integrated ionisation rate
;endfor
;photoi_tot = total(photoi_ht,1)
;
;;h=6.62607004e-34
;;c=3.e8
;
;E=12397./(0.5*(SQ_wv1+SQ_wv2))   ;eV cm-2 s-1
;ssflux_E=ssflux*E
;
;ph_E= 25 ;25eV to ionizie N2
;E=ph_E
;photoi_E=photoi_tot*E;eV cm-2 s-1    *6.242e18
;
;;_________________________________________________________________________________________________
;;
;;the height integrated photoelectron ionisation energy
;
;
;eiionz_ht = fltarr(n_elements(zz),n_elements(SQ_wv1))  ;height integrated photoionisation energy
;
;;eiionz_transp:Array[species,alt,  wavelength]
;;eiionz_stot=total(x9_szalo_eiionz_transp,1)  ;total photoization in all species
;eiionz_stot=reform(x9_szalo_eiionz_transp[2,*,*])  ;N2 photoization rate
;
;
;for l=0, n_elements(SQ_wv1)-1 do begin
;  eiionz_ht(*,l)=eiionz_stot[*,l]*dz*1.e5 ;height integrated ionisation rate
;endfor
;eiionz_tot = total(eiionz_ht,1)  ;function of wavlength
;
;E=(12397./(0.5*(SQ_wv1+SQ_wv2)))-ph_E;-500.
;
;eiionz_E=eiionz_tot*E;eV cm-2 s-1
;
;;for i=3, 21 do $
;;  print,SQ_wv1[i],SQ_wv2[i],ssflux_E[i],photoi_E[i],eiionz_E[i]
;;print,ssflux_E[0:11]
;;print,photoi_E[0:11]
;;print,eiionz_E[0:9]
;;
;;
;;s=[4.64679e+08,  4.10889e+09,  2.11353e+10,  2.06691e+10 , 2.94863e+10,  9.34729e+10 , 1.09602e+11,  1.51535e+11,  3.83339e+11,  4.21854e+11,  1.97834e+12,  2.35365e+12,$
;;   1.79600e+13 , 1.27308e+13,  2.02927e+13,  2.57466e+13,  4.05461e+13,  3.23691e+13,  2.48097e+13]
;;
;;p=[432926.,  6.95858e+06,  5.02719e+07,  6.35841e+07,  1.11389e+08,  4.51587e+08,  6.85419e+08,  1.16429e+09 , 3.76635e+09,  5.33064e+09,  3.31953e+10,  5.21237e+10,$
;;   6.35841e+07,  1.11389e+08,  4.51587e+08,  6.85419e+08,  1.16429e+09,  3.76635e+09,  5.33064e+09]
;;
;;e=[ 1.45963e+11 , 8.94228e+11,  3.28687e+12,  2.52929e+12,  2.98737e+12,  7.71300e+12,  7.19096e+12,  8.31525e+12,  1.73364e+13,  1.54176e+13,$
;;    2.52929e+12,  2.98737e+12,  7.71300e+12,  7.19096e+12,  8.31525e+12,  1.73364e+13,  1.54176e+13]
;;    
;;    
;;    wv=[0.5,1.0,1.5,2.0,2.5,3.0,4,5,6,8,10,14,18,32,70,155,224,290,320]
;
 
  

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
posleft=[0.05,0.15,0.33,0.9]
posmid=[0.39,0.15,0.68,0.9]
posright=[0.73,0.15,0.98,0.9]
thck=3.
ln_thck=5
sym_thck=2
sym_incr=10
leg_loc=[0.74, 0.90]

leg_loc_left=[0.04, 0.9]
leg_loc_right=[0.98, 0.9]
leg_transp=50
maj=['O ','O!L2!N ','N!L2!N ']

sym=["*","+","tu"]

;w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())
;
;N2_pi1=plot(wv,s,name='Solar flux',$
;  xtitle='Wavelength ($\AA$)',ytitle='eV cm!U-2!N s!U-1!N ',$
;  xstyle=1,ystyle=1,yrange=[1.e4,1.e16],xrange=[0.1,1.e3],xlog=1,ylog=1,histogram=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],font_size=fnt_sz,/current)
;
;N2_pi2=plot(wv,p,name='N!L2!N Photoionisation Energy',histogram=1,$
;  thick=thck,color=colors[3],/overplot)
;;
;N2_pi3=plot(wv,e,name='N!L2!N Photoelectronionisation  Energy',histogram=1,$
;  thick=thck,color=colors[4],/overplot,symbol=sym[0])
;  
;  N2_pi4=plot(wv,e+p,name='N!L2!N Total P-Ionz+E-Ionz  Energy',histogram=1,$
;    thick=thck,color=colors[0],/overplot)
;
;
;leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot

tit_fl=['M5 flare -','X9 flare -']
;tit_fl=['X9 flare -','M5 flare -']






leg_nam=['0.5-1.0','1.0-1.5','1.5-2.0','2.0-2.5','2.5-3.0','3.0-4.0','NRL 0.5-4','SQ 0.5-4']

;leg_nam=['8-10','10-14','14-18','NRL 8-18','SQ 8-18']
;leg_nam=['4-5','5-6','6-8','NRL 4-8','SQ 4-8']


;_________________________________________________________________________



;  N2
;
xr=[1.e-10,1.e-3]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
yt='Altitude (km)'
titl='Photoionization'

;X9 flare

w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

N2_pi1=plot(x9_szalo_photoi[2,*,0],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posleft,/current)

N2_pi2=plot(x9_szalo_photoi[2,*,1],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[2,*,2],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)

N2_pi4=plot(x9_szalo_photoi[2,*,3],zz,name=leg_nam[3],$
  thick=thck,color=colors[5],/overplot)

N2_pi5=plot(x9_szalo_photoi[2,*,4],zz,name=leg_nam[4],$
  thick=thck,color=colors[6],/overplot)

N2_pi6=plot(x9_szalo_photoi[2,*,5],zz,name=leg_nam[5],$
  thick=thck,color=colors[7],/overplot)

;
;N2_pi7=plot(total(x9_szalo_photoi[2,*,0:5],3),zz,name=leg_nam[6],$
;    thick=thck,color=colors[2],/overplot)
;
;N2_pi8=plot(SQx9_szalo_photoi[2,*,0],SQ_tau[*,0],name=leg_nam[7],$
;  thick=thck,color=colors[0],/overplot)

leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3,N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)




xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(x9_szalo_eiionz_transp[2,*,0],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[2,*,1],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[2,*,2],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)

N2_ei4=plot(x9_szalo_eiionz_transp[2,*,3],zz,name=leg_nam[3],$
  thick=thck,color=colors[5],/overplot)

N2_ei5=plot(x9_szalo_eiionz_transp[2,*,4],zz,name=leg_nam[4],$
  thick=thck,color=colors[6],/overplot)

N2_ei6=plot(x9_szalo_eiionz_transp[2,*,5],zz,name=leg_nam[5],$
  thick=thck,color=colors[7],/overplot)

;
;N2_ei7=plot(total(x9_szalo_eiionz_transp[2,*,0:5],3),zz,name=leg_nam[6],$
;    thick=thck,color=colors[2],/overplot)
;
;N2_ei8=plot(SQx9_szalo_eiionz_transp[2,*,0],SQ_tau[*,0],name=leg_nam[7],$
;  thick=thck,color=colors[0],/overplot)


leg5 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei4,N2_ei5,N2_ei6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)
;leg5 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3,N2_ei7,N2_ei8], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________

xr=[10.,1000]
xt='Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(x9_szalo_eiionz_transp[2,*,0]/x9_szalo_photoi[2,*,0],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=1,$;yrange=yr,,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[2,*,1]/x9_szalo_photoi[2,*,1],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[2,*,2]/x9_szalo_photoi[2,*,2],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)

N2_pp4=plot(x9_szalo_eiionz_transp[2,*,3]/x9_szalo_photoi[2,*,3],zz,name=leg_nam[3],$
  thick=thck,color=colors[5],/overplot)

N2_pp5=plot(x9_szalo_eiionz_transp[2,*,4]/x9_szalo_photoi[2,*,4],zz,name=leg_nam[4],$
  thick=thck,color=colors[6],/overplot)

N2_pp6=plot(x9_szalo_eiionz_transp[2,*,5]/x9_szalo_photoi[2,*,5],zz,name=leg_nam[5],$
  thick=thck,color=colors[7],/overplot)

;N2_pp7=plot(total(x9_szalo_eiionz_transp[2,*,0:5],3)/total(x9_szalo_photoi[2,*,0:5],3),zz,name=leg_nam[6],$
;  thick=thck,color=colors[2],/overplot)

;N2_pp8=plot(SQx9_szalo_eiionz_transp[2,*,0]/SQx9_szalo_photoi[2,*,0],SQ_tau[*,0],name=leg_nam[7],$
;  thick=thck,color=colors[0],/overplot)


leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp4,N2_pp5,N2_pp6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)
;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3,N2_pp7,N2_pp8], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________________________________________________________________________________________________________________________________________

;;SQ  N2
;
;xr=[10.,10000.]
;yr=[50.,max(zz)]
;xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)'
;titl='Photoionization'
;
;;X9 flare
;
;w6 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())
;
;N2_pi1=plot(reform(total(SQx9_szalo_photoi[2,*,*],3)),zz,$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],font_size=fnt_sz,/current)
;
;
;___________________________________________________________________________________________________________________________________________________________________________________

;;SQ Pe/Pi ratio plot for 4 bins
;;4 32-70 A
;;5 70-155A
;;6 155-224 A
;;8 290-320 A
;
;posupleft= [0.05,0.6,0.45,0.9]
;posupright= [0.5,0.6,0.9,0.9]
;posdoleft= [0.05,0.15,0.45,0.5]
;posdoright= [0.5,0.15,0.9,0.5]
;
;leg_nam=['32-70 A','70-155A','155-224 A','290-320 A']
;
;yr=[4,0]
;xt='Pe/Pi '
;yt= 'Optical Depth'
;
;
;w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())
;
;xr=[1.,100]
;N2_pp1=plot(SQx9_szalo_eiionz_transp[2,*,4]/SQx9_szalo_photoi[2,*,4],tau[*,4],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+leg_nam[0],$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],font_size=fnt_sz,position=posupleft,/current)
;
;
;N2_pp2=plot(SQx9_szalo_eiionz_transp[2,*,5]/SQx9_szalo_photoi[2,*,5],tau[*,5],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+leg_nam[1],$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],font_size=fnt_sz,position=posupright,/current)
;
; xr=[0.1,10.0]
; N2_pp3=plot(SQx9_szalo_eiionz_transp[2,*,6]/SQx9_szalo_photoi[2,*,6],tau[*,6],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+leg_nam[2],$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[1],font_size=fnt_sz,position=posdoleft,/current)
;
;
;xr=[0.01,1.0]
;N2_pp4=plot(SQx9_szalo_eiionz_transp[2,*,8]/SQx9_szalo_photoi[2,*,8],tau[*,8],$
;    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+leg_nam[3],$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[1],font_size=fnt_sz,position=posdoright,/current)


;___________________________________________________________________________________________________________________________________________________________________________________

;  leg_nam=['SQ cross-sections','H+F cross-sections']
;  xr=[1e-3,1e5]
;yr=[min(zz),max(zz)]
;  xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
;  titl='Photoionization'
;
;
;
;  w7 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())
;
;N2_pi1=plot(SQx9_szalo_photoi[2,*,0],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[1],font_size=fnt_sz,position=posleft,/current)
;
;N2_pi2=plot(SQx9_szalo_photoi_v6[2,*,0],zz,name=leg_nam[1],$
;    thick=thck,color=colors[3],/overplot)
;
;leg7 = LEGEND(TARGET=[N2_pi1,N2_pi2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)
;
;
;  xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
;  titl='Photoelectron Ionization'
;
;  N2_ei1=plot(SQx9_szalo_eiionz_transp[2,*,0],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[1],font_size=fnt_sz,position=posmid,/current)
;
;  N2_ei2=plot(SQx9_szalo_eiionz_transp_v6[2,*,0],zz,name=leg_nam[1],$
;    thick=thck,color=colors[3],/overplot)
;
;
;  leg8 = LEGEND(TARGET=[N2_ei1,N2_ei2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)
;
;
;  xr=[10.,1000]
;  xt='Pe/Pi '
;  titl='Pe/Pi'
;
;  N2_pp1=plot(SQx9_szalo_eiionz_transp[2,*,0]/SQx9_szalo_photoi[2,*,0],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[1],font_size=fnt_sz,position=posright,/current)
;
;  N2_pp2=plot(SQx9_szalo_eiionz_transp_v6[2,*,0]/SQx9_szalo_photoi_v6[2,*,0],zz,name=leg_nam[1],$
;    thick=thck,color=colors[3],/overplot)
;
;
;  leg9 = LEGEND(TARGET=[N2_pp1,N2_pp2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)


;___________________________________________________________________________________________________________________________________________________________________________________



;;Pi=12397./[(.5+1.)/2.,(1.+1.5)/2.,(1.5+2.)/2.,(2.+2.5)/2.,(2.5+3.)/2.,(3.+4.)/2.,(0.5+4.)/2.]
;Pi=[12397./((0.5+4)/2)]
;;Pi=12397./[(4.+5.)/2.,(5.+6.)/2.,(6.+8.)/2.,(4.+8.)/2.];;Pi=12397./[(.5+1.)/2.,(1.+1.5)/2.,(1.5+2.)/2.,(2.+2.5)/2.,(2.5+3.)/2.,(3.+4.)/2.,(0.5+4.)/2.]
;Pi=[12397./((0.5+4)/2)]
;;Pi=12397./[(4.+5.)/2.,(5.+6.)/2.,(6.+8.)/2.,(4.+8.)/2.]
;
;;Pi=12397./[(8.+10.)/2.,(10.+14.)/2.,(14.+18.)/2.,(8.+18.)/2.]
;
;auger_wvln=30.; N2
;auger_energy=12397./auger_wvln
;Ionz_en=[15.60,29.0]
;;Ionz_en= [15.60,16.70, 18.80, 25.3,29.0,33.40,36.80,37.8, 43.6,413.233]
;
;PePi= fltarr(n_elements(Ionz_en),n_elements(Pi))
;
;for i=0, n_elements(Pi)-1 do begin
;
;  for j= 0, n_elements(Ionz_en)-1 do begin
;    if (j eq 0) then begin
;;      print,":)"
;       Pepi(j,i)=(Pi(i)-auger_energy-Ionz_en(j))/Ionz_en(j)
;    endif else $
;                 Pepi(j,i)=(Pi(i)-Ionz_en(j))/Ionz_en(j)
;  endfor
;endfor
;end
;

;
;;Pi=12397./[(8.+10.)/2.,(10.+14.)/2.,(14.+18.)/2.,(8.+18.)/2.]
;
;auger_wvln=30.; N2
;auger_energy=12397./auger_wvln
;Ionz_en=[15.60,29.0]
;;Ionz_en= [15.60,16.70, 18.80, 25.3,29.0,33.40,36.80,37.8, 43.6,413.233]
;
;PePi= fltarr(n_elements(Ionz_en),n_elements(Pi))
;
;for i=0, n_elements(Pi)-1 do begin
;
;  for j= 0, n_elements(Ionz_en)-1 do begin
;    if (j eq 0) then begin
;;      print,":)"
;       Pepi(j,i)=(Pi(i)-auger_energy-Ionz_en(j))/Ionz_en(j)
;    endif else $
;                 Pepi(j,i)=(Pi(i)-Ionz_en(j))/Ionz_en(j)
;  endfor
;endfor
end
;




;Seems like there are enormous differences between the SQ binning and ours in the photoelectron ionization rate.
;Like a factor of 10 at some altitudes, like between 70-100 km. And yet, the differences at these same altitudes in
;the N2 photoionization rates are small. Why is that? Are we sure that we are getting the same pe/pi as SQ05?
; Ultimately, we will have to plot pe/pi on an optical depth vertical axis like Figure 3 of their paper.
;  They then picked off the value at tau = 1 as their value for the table. We need to show that we can reproduce that
;  figure and indeed do it for each bin of SQ to show that we can reproduce their table A2-A4.
;Theoritical calculation of photoionisation and verification of the model photoionisation
;Aug 24, 2020
;checking calculatios for each bin and only low solar zenith angle


;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'


;______________________________________________________________________________________________________


;restore,filename=floc+'SQpeakm5_wv_szalo.sav'
;restore,filename=floc+'SQpeakx9_wv_szalo_v2.sav'
;restore,filename=floc+'SQpeakx9_wv_szalo_35eV.sav'
restore,filename=floc+'SQpeakx9_wv_szalo_test.sav'


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
;
;restore,filename=floc+'SQpeakx9_wv_szalo_v6.sav'  ;with new cross-section bins and states
;
;SQx9_szalo_photoi_v6= photi_wv
;SQx9_szalo_eiionz_loc_v6= eionz_wv_loc
;SQx9_szalo_eiionz_transp_v6= eionz_wv_transp

;SQ_wv1= wv1
;SQ_wv2= wv2
;______________________________________________________________________________________________________


;restore,filename=floc+'peakx9_wv_szalo_v8.sav'
;restore,filename=floc+'peakm5_wv_szalo.sav'
restore,filename=floc+'peakx9_ssflux1.sav';'peakx9_wv_szalo_2states.sav'


; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon, tau


nrl_wv1= wv1
nrl_wv2= wv2
nrl_tau =tau

x9_szalo_photoi= photi_wv
x9_szalo_eiionz_loc= eionz_wv_loc
x9_szalo_eiionz_transp= eionz_wv_transp


close,/all

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;
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
;
;
;;the height integrated photoionisation energy should add up to the solar flux,
;;atleast for the lower wavlength bins
;
;photoi_ht = fltarr(n_elements(zz),n_elements(nrl_wv1))  ;height integrated photoionisation energy
;
;;photoi_wv:Array[alt, species, wavelength]
;photoi_stot=total(x9_szalo_photoi,1)  ;total photoization in all species
;
;  for l=0, n_elements(nrl_wv1)-1 do begin
;    photoi_ht(*,l)=photoi_stot[*,l]*dz*1.e5 ;height integrated ionisation rate
;  endfor
;  photoi_tot = total(photoi_ht,1)
;h=6.62607004e-34
;c=3.e8
;E=h*c/(0.5*1.e-10*(nrl_wv1+nrl_wv2))   ;J cm-2 s-1
;photoi_tot=photoi_tot;*E;*6.242e18
;ssflux_tot=ssflux;*E
;
;print,total(photoi_tot),total(ssflux_tot)
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

;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot

tit_fl=['M5 flare -','X9 flare -']
;tit_fl=['X9 flare -','M5 flare -']



leg_nam=['4-5','5-6','6-8','NRL 4-8','SQ 4-8']

;_________________________________________________________________________



  ;  N2

xr=[1.e-14,1.e-4]
yr=[min(zz),max(zz)]

;xr=[]
;yr=[]


yt='Altitude (km)'

xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)
titl='Photoionization'

;X9 flare

w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

N2_pi1=plot(x9_szalo_photoi[2,*,6],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posleft,/current)

N2_pi2=plot(x9_szalo_photoi[2,*,7],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)
;
N2_pi3=plot(x9_szalo_photoi[2,*,8],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)

;
;
;N2_pi4=plot(total(x9_szalo_photoi[2,*,6:8],3),zz,name=leg_nam[3],$
;  thick=thck,color=colors[2],/overplot)
;
;N2_pi5=plot(SQx9_szalo_photoi[2,*,1],SQ_tau[*,1],name=leg_nam[4],$
;  thick=thck,color=colors[0],/overplot)

leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)





xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) '
titl='Photoelectron Ionization'

N2_ei1=plot(x9_szalo_eiionz_transp[2,*,6],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posmid,/current)

N2_ei2=plot(x9_szalo_eiionz_transp[2,*,7],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)

N2_ei3=plot(x9_szalo_eiionz_transp[2,*,8],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)

;
;
;N2_ei4=plot(total(x9_szalo_eiionz_transp[2,*,6:8],3),zz,name=leg_nam[3],$
;  thick=thck,color=colors[2],/overplot)

;N2_ei5=plot(SQx9_szalo_eiionz_transp[2,*,1],SQ_tau[*,1],name=leg_nam[4],$
;  thick=thck,color=colors[0],/overplot)


leg5 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;_________________________________________________________________________

xr=[10.,1000.]
xt='Pe/Pi '
titl='Pe/Pi'

N2_pp1=plot(x9_szalo_eiionz_transp[2,*,6]/x9_szalo_photoi[2,*,6],zz,name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],font_size=fnt_sz,position=posright,/current)

N2_pp2=plot(x9_szalo_eiionz_transp[2,*,7]/x9_szalo_photoi[2,*,7],zz,name=leg_nam[1],$
  thick=thck,color=colors[3],/overplot)

N2_pp3=plot(x9_szalo_eiionz_transp[2,*,8]/x9_szalo_photoi[2,*,8],zz,name=leg_nam[2],$
  thick=thck,color=colors[4],/overplot)


;
;N2_pp4=plot(total(x9_szalo_eiionz_transp[2,*,6:8],3)/total(x9_szalo_photoi[2,*,6:8],3),zz,name=leg_nam[3],$
;  thick=thck,color=colors[2],/overplot)

;N2_pp5=plot(SQx9_szalo_eiionz_transp[2,*,1]/SQx9_szalo_photoi[2,*,1],SQ_tau[*,1],name=leg_nam[4],$
;  thick=thck,color=colors[0],/overplot)


leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;___________________________________________________________________________________________________________________________________________________________________________________


;
end
;




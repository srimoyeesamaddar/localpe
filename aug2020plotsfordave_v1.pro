;Theoritical calculation of photoionisation and verification of the model photoionisation
;Aug 17, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'

m5_lat=intarr(3)
m5_lon=intarr(3)
m5_sza=intarr(3)

x9_lat=intarr(3)
x9_lon=intarr(3)
x9_sza=intarr(3)

restore,filename=floc+'SQpeakm5_wv_szalo.sav'
; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon

SQm5_szalo_photoi= photi_wv
SQm5_szalo_eiionz_loc= eionz_wv_loc
SQm5_szalo_eiionz_transp= eionz_wv_transp


SQ_wv1= wv1
SQ_wv2= wv2

m5_lat[0]=fix(lat)
m5_lon[0]= fix(lon)
m5_sza[0]=fix( round( sza*!radeg))


restore,filename=floc+'SQpeakm5_wv_szamed.sav'

SQm5_szamed_photoi= photi_wv
SQm5_szamed_eiionz_loc= eionz_wv_loc
SQm5_szamed_eiionz_transp= eionz_wv_transp


m5_lat[1]= fix(lat)
m5_lon[1]= fix(lon)
m5_sza[1]= fix(round( sza*!radeg))


restore,filename=floc+'SQpeakm5_wv_szahi.sav'


SQm5_szahi_photoi= photi_wv
SQm5_szahi_eiionz_loc= eionz_wv_loc
SQm5_szahi_eiionz_transp= eionz_wv_transp


m5_lat[2]= fix(lat)
m5_lon[2]= fix(lon)
m5_sza[2]= fix(round( sza*!radeg))

;______________________________________________________________________________________________________

restore,filename=floc+'SQpeakx9_wv_szalo_v2.sav'
; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon

SQx9_szalo_photoi= photi_wv
SQx9_szalo_eiionz_loc= eionz_wv_loc
SQx9_szalo_eiionz_transp= eionz_wv_transp

x9_lat[0]= fix(lat)
x9_lon[0]= fix(lon)
x9_sza[0]=  fix(round(sza*!radeg))


restore,filename=floc+'SQpeakx9_wv_szamed.sav'

SQx9_szamed_photoi= photi_wv
SQx9_szamed_eiionz_loc= eionz_wv_loc
SQx9_szamed_eiionz_transp= eionz_wv_transp


x9_lat[1]= fix(lat)
x9_lon[1]= fix(lon)
x9_sza[1]=fix( round( sza*!radeg))


restore,filename=floc+'SQpeakx9_wv_szahi.sav'


SQx9_szahi_photoi= photi_wv
SQx9_szahi_eiionz_loc= eionz_wv_loc
SQx9_szahi_eiionz_transp= eionz_wv_transp


x9_lat[2]= fix(lat)
x9_lon[2]= fix(lon)
x9_sza[2]= fix (round(sza*!radeg))
;______________________________________________________________________________________________________


restore,filename=floc+'peakm5_wv_szalo.sav'
; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon

m5_szalo_photoi= photi_wv
m5_szalo_eiionz_loc= eionz_wv_loc
m5_szalo_eiionz_transp= eionz_wv_transp


nrl_wv1= wv1
nrl_wv2= wv2


restore,filename=floc+'peakm5_wv_szamed.sav'

m5_szamed_photoi= photi_wv
m5_szamed_eiionz_loc= eionz_wv_loc
m5_szamed_eiionz_transp= eionz_wv_transp


restore,filename=floc+'peakm5_wv_szahi.sav'

m5_szahi_photoi= photi_wv
m5_szahi_eiionz_loc= eionz_wv_loc
m5_szahi_eiionz_transp= eionz_wv_transp

;______________________________________________________________________________________________________


restore,filename=floc+'peakx9_wv_szalo_v2.sav'
; file has:eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon

x9_szalo_photoi= photi_wv
x9_szalo_eiionz_loc= eionz_wv_loc
x9_szalo_eiionz_transp= eionz_wv_transp

restore,filename=floc+'peakx9_wv_szamed.sav'

x9_szamed_photoi= photi_wv
x9_szamed_eiionz_loc= eionz_wv_loc
x9_szamed_eiionz_transp= eionz_wv_transp


restore,filename=floc+'peakx9_wv_szahi.sav'

x9_szahi_photoi= photi_wv
x9_szahi_eiionz_loc= eionz_wv_loc
x9_szahi_eiionz_transp= eionz_wv_transp




SQm5_szalo_pepi= total(SQm5_szalo_eiionz_transp,3)/total(SQm5_szalo_photoi,3)
SQm5_szalo_pepi[where(~finite(SQm5_szalo_pepi),/null)]=0.

SQx9_szalo_pepi= total(SQx9_szalo_eiionz_transp,3)/total(SQx9_szalo_photoi,3)
;SQx9_szalo_pepi[where(~finite(SQx9_szalo_pepi),/null)]=0.


close,/all
;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;Making new arrays combining NRL and SQ bins for final results of photoelectron model


  n_wvl=n_elements(nrl_wv1[0:11])+n_elements(SQ_wv1[3:21])
  wv_lo= fltarr(n_wvl)
  wv_hi= fltarr(n_wvl)

   wv_lo=[nrl_wv1[0:11],SQ_wv1[3:21]]
   wv_hi=[nrl_wv2[0:11],SQ_wv2[3:21]]

; M5 flare

;      Low sza
   
   m5_photoi_szalo=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_loc_szalo=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_transp_szalo=fltarr(3,n_elements(zz),n_wvl)
;   m5_pepi_szalo= fltarr(3,n_elements(zz),n_wvl)
       
   m5_photoi_szalo[*,*,0:11] =m5_szalo_photoi[*,*,0:11]
   m5_photoi_szalo[*,*,12:-1]=SQm5_szalo_photoi[*,*,3:21]
   
   m5_eiionz_loc_szalo[*,*,0:11] = m5_szalo_eiionz_loc[*,*,0:11]
   m5_eiionz_loc_szalo[*,*,12:-1] = SQm5_szalo_eiionz_loc[*,*,3:21]
   
   m5_eiionz_transp_szalo[*,*,0:11] =m5_szalo_eiionz_transp[*,*,0:11]     
   m5_eiionz_transp_szalo[*,*,12:-1] = SQm5_szalo_eiionz_transp[*,*,3:21]
   
   m5_pepi_szalo=total(m5_eiionz_transp_szalo,3)/total(m5_photoi_szalo,3)
   m5_pepi_szalo[where(~finite(m5_pepi_szalo),/null)]=0.
   


   ;      Medium sza
   
   m5_photoi_szamed=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_loc_szamed=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_transp_szamed=fltarr(3,n_elements(zz),n_wvl)
;   m5_pepi_szamed= fltarr(3,n_elements(zz),n_wvl)

   m5_photoi_szamed[*,*,0:11] =m5_szamed_photoi[*,*,0:11]
   m5_photoi_szamed[*,*,12:-1]=SQm5_szamed_photoi[*,*,3:21]

   m5_eiionz_loc_szamed[*,*,0:11] = m5_szamed_eiionz_loc[*,*,0:11]
   m5_eiionz_loc_szamed[*,*,12:-1] = SQm5_szamed_eiionz_loc[*,*,3:21]

   m5_eiionz_transp_szamed[*,*,0:11] =m5_szamed_eiionz_transp[*,*,0:11]
   m5_eiionz_transp_szamed[*,*,12:-1] = SQm5_szamed_eiionz_transp[*,*,3:21]
   
   m5_pepi_szamed= total(m5_eiionz_transp_szamed,3)/total( m5_photoi_szamed,3)
   m5_pepi_szamed[where(~finite(m5_pepi_szamed),/null)]=0.
   
       
   ;      High sza
   
   
   m5_photoi_szahi=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_loc_szahi=fltarr(3,n_elements(zz),n_wvl)
   m5_eiionz_transp_szahi=fltarr(3,n_elements(zz),n_wvl)
;   m5_pepi_szahi= fltarr(3,n_elements(zz),n_wvl)
   
   m5_photoi_szahi[*,*,0:11] =m5_szahi_photoi[*,*,0:11]
   m5_photoi_szahi[*,*,12:-1]=SQm5_szahi_photoi[*,*,3:21]

   m5_eiionz_loc_szahi[*,*,0:11] = m5_szahi_eiionz_loc[*,*,0:11]
   m5_eiionz_loc_szahi[*,*,12:-1] = SQm5_szahi_eiionz_loc[*,*,3:21]

   m5_eiionz_transp_szahi[*,*,0:11] =m5_szahi_eiionz_transp[*,*,0:11]
   m5_eiionz_transp_szahi[*,*,12:-1] = SQm5_szahi_eiionz_transp[*,*,3:21]

   m5_pepi_szahi= total(m5_eiionz_transp_szahi,3)/ total(m5_photoi_szahi,3)
   m5_pepi_szahi[where(~finite(m5_pepi_szahi),/null)]=0.

   ; X9 flare

   ;      Low sza

   x9_photoi_szalo=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_loc_szalo=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_transp_szalo=fltarr(3,n_elements(zz),n_wvl)
;   x9_pepi_szalo= fltarr(3,n_elements(zz),n_wvl)

   x9_photoi_szalo[*,*,0:11] =x9_szalo_photoi[*,*,0:11]
   x9_photoi_szalo[*,*,12:-1]=SQx9_szalo_photoi[*,*,3:21]

   x9_eiionz_loc_szalo[*,*,0:11] = x9_szalo_eiionz_loc[*,*,0:11]
   x9_eiionz_loc_szalo[*,*,12:-1] = SQx9_szalo_eiionz_loc[*,*,3:21]

   x9_eiionz_transp_szalo[*,*,0:11] =x9_szalo_eiionz_transp[*,*,0:11]
   x9_eiionz_transp_szalo[*,*,12:-1] = SQx9_szalo_eiionz_transp[*,*,3:21]

   x9_pepi_szalo=total( x9_eiionz_transp_szalo,3)/total( x9_photoi_szalo,3)
;   x9_pepi_szalo[where(~finite(x9_pepi_szalo),/null)]=0.

   ;      Medium sza

   x9_photoi_szamed=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_loc_szamed=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_transp_szamed=fltarr(3,n_elements(zz),n_wvl)
;   x9_pepi_szamed=fltarr(3,n_elements(zz),n_wvl)

   x9_photoi_szamed[*,*,0:11] =x9_szamed_photoi[*,*,0:11]
   x9_photoi_szamed[*,*,12:-1]=SQx9_szamed_photoi[*,*,3:21]

   x9_eiionz_loc_szamed[*,*,0:11] = x9_szamed_eiionz_loc[*,*,0:11]
   x9_eiionz_loc_szamed[*,*,12:-1] = SQx9_szamed_eiionz_loc[*,*,3:21]

   x9_eiionz_transp_szamed[*,*,0:11] =x9_szamed_eiionz_transp[*,*,0:11]
   x9_eiionz_transp_szamed[*,*,12:-1] = SQx9_szamed_eiionz_transp[*,*,3:21]

   x9_pepi_szamed=total(x9_eiionz_transp_szamed,3)/total(  x9_photoi_szamed,3)
   x9_pepi_szamed[where(~finite(x9_pepi_szamed),/null)]=0.


   ;      High sza


   x9_photoi_szahi=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_loc_szahi=fltarr(3,n_elements(zz),n_wvl)
   x9_eiionz_transp_szahi=fltarr(3,n_elements(zz),n_wvl)
;   x9_pepi_szahi= fltarr(3,n_elements(zz),n_wvl)
   
   x9_photoi_szahi[*,*,0:11] =x9_szahi_photoi[*,*,0:11]
   x9_photoi_szahi[*,*,12:-1]=SQx9_szahi_photoi[*,*,3:21]

   x9_eiionz_loc_szahi[*,*,0:11] = x9_szahi_eiionz_loc[*,*,0:11]
   x9_eiionz_loc_szahi[*,*,12:-1] = SQx9_szahi_eiionz_loc[*,*,3:21]

   x9_eiionz_transp_szahi[*,*,0:11] =x9_szahi_eiionz_transp[*,*,0:11]
   x9_eiionz_transp_szahi[*,*,12:-1] = SQx9_szahi_eiionz_transp[*,*,3:21]

   x9_pepi_szahi=total(x9_eiionz_transp_szahi,3)/ total( x9_photoi_szahi,3)
   x9_pepi_szahi[where(~finite(x9_pepi_szahi),/null)]=0.
   


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
posleft=[0.05,0.15,0.45,0.9]
posright=[0.55,0.15,0.95,0.9]
thck=3.
ln_thck=5
sym_thck=2
sym_incr=10
leg_loc=[0.92, 0.88]

leg_loc_left=[0.04, 0.9]
leg_loc_right=[0.98, 0.9]

maj=['O ','O!L2!N ','N!L2!N ']

  sym=["*","+","tu"]

;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot

tit_fl=['M5 flare -','X9 flare -']

leg_nam=['lat= '+strtrim(string(m5_lat[0]),1)+' lon= '+strtrim(string(m5_lon[0]),1)+' sza= '+strtrim(string(m5_sza[0]),1),$
         'lat= '+strtrim(string(m5_lat[1]),1)+' lon= '+strtrim(string(m5_lon[1]),1)+' sza= '+strtrim(string(m5_sza[1]),1),$
         'lat= '+strtrim(string(m5_lat[2]),1)+' lon= '+strtrim(string(m5_lon[2]),1)+' sza= '+strtrim(string(m5_sza[2]),1),$
         'lat= '+strtrim(string(x9_lat[0]),1)+' lon= '+strtrim(string(x9_lon[0]),1)+' sza= '+strtrim(string(x9_sza[0]),1),$
         'lat= '+strtrim(string(x9_lat[1]),1)+' lon= '+strtrim(string(x9_lon[1]),1)+' sza= '+strtrim(string(x9_sza[1]),1),$
         'lat= '+strtrim(string(x9_lat[2]),1)+' lon= '+strtrim(string(x9_lon[2]),1)+' sza= '+strtrim(string(x9_sza[2]),1)]

yr=[min(zz),max(zz)]
xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N)'
yt='Altitude (km)'
titl='Photoionization'




;;  O Photoionisation plot
;
;xr=[1.e-9,1.e4]
;
;w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;O_pi1=plot(total(m5_photoi_szalo[0,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;
;O_pi2=plot(total(m5_photoi_szamed[0,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;  
;  
;O_pi3=plot(total(m5_photoi_szahi[0,*,*],3),zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],/overplot)
;
;
;leg1 = LEGEND(TARGET=[O_pi1,O_pi2,O_pi3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;O_pi4=plot(total(x9_photoi_szalo[0,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;O_pi5=plot(total(x9_photoi_szamed[0,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O_pi6=plot(total(x9_photoi_szahi[0,*,*],3),zz,name=leg_nam[5],$
;    thick=thck,color=colors[2],/overplot)
;
;
;leg2 = LEGEND(TARGET=[O_pi4,O_pi5,O_pi6], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;_________________________________________________________________________
;
;;O2 Photoionisation plot
;
;
;
;w2 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;O2_pi1=plot(total(m5_photoi_szalo[1,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;
;O2_pi2=plot(total(m5_photoi_szamed[1,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;  
;  
;O2_pi3=plot(total(m5_photoi_szahi[1,*,*],3),zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],/overplot)
;
;leg3 = LEGEND(TARGET=[O2_pi1,O2_pi2,O2_pi3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;O2_pi4=plot(total(x9_photoi_szalo[1,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;O2_pi5=plot(total(x9_photoi_szamed[1,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O2_pi6=plot(total(x9_photoi_szahi[1,*,*],3),zz,name=leg_nam[5],$
;    thick=thck,color=colors[2],/overplot)
;
;
;leg4 = LEGEND(TARGET=[O2_pi4,O2_pi5,O2_pi6], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w2.Save, plt_loc+"O2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;_________________________________________________________________________
;
;;N2 Photoionisation plot
;
;
;w3 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;N2_pi1=plot(total(m5_photoi_szalo[2,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;
;N2_pi2=plot(total(m5_photoi_szamed[2,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;  
;  
;N2_pi3=plot(total(m5_photoi_szahi[2,*,*],3),zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],/overplot)
;
;leg5 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;N2_pi4=plot(total(x9_photoi_szalo[2,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;N2_pi5=plot(total(x9_photoi_szamed[2,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;N2_pi6=plot(total(x9_photoi_szahi[2,*,*],3),zz,name=leg_nam[5],$
;    thick=thck,color=colors[2],/overplot)2.4964066e+08
;
;leg6 = LEGEND(TARGET=[N2_pi4,N2_pi5,N2_pi6], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;;w3.Save, plt_loc+"N2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;___________________________________________________________________________________________________________________________________________________________________________________
;
;
;
;;Photoelectron Ionisation plot
;
;
;
;;yr=[min(zz),max(zz)]
;xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) !C (solid line- local, starred line- transport)'
;;yt='Altitude (km)'
;titl='Photoelectron Ionization'
;
;
;
;
;;  O Photoelectron ionisation plot
;
;;xr=[1.e-9,1.e4]
;
;w4 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;O_ei1=plot(total(m5_eiionz_loc_szalo[0,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;O_ei2=plot(total(m5_eiionz_loc_szamed[0,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;O_ei3=plot(total(m5_eiionz_loc_szahi[0,*,*],3),zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;
;O_ei4=plot(total(m5_eiionz_transp_szalo[0,*,*],3),zz,name=leg_nam[1],$
;    thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O_ei5=plot(total(m5_eiionz_transp_szamed[0,*,*],3),zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O_ei6=plot(total(m5_eiionz_transp_szahi[0,*,*],3),zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;leg1 = LEGEND(TARGET=[O_ei1,O_ei2,O_ei3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;
;
;O_ei7=plot(total(x9_eiionz_loc_szalo[0,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;O_ei8=plot(total(x9_eiionz_loc_szamed[0,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O_ei9=plot(total(x9_eiionz_loc_szahi[0,*,*],3),zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;
;
;O_ei10=plot(total(x9_eiionz_transp_szalo[0,*,*],3),zz,name=leg_nam[3],$
;    thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O_ei11=plot(total(x9_eiionz_transp_szamed[0,*,*],3),zz,name=leg_nam[4],$
;    thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O_ei12=plot(total(x9_eiionz_transp_szahi[0,*,*],3),zz,name=leg_nam[5],$
;    thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;
;leg2 = LEGEND(TARGET=[O_ei7,O_ei8,O_ei9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;_________________________________________________________________________
;
;;O2 Photoelectron ionisation plot
;
;w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;;xr=[1.e-9,1.e4]
;
;O2_ei1=plot(total(m5_eiionz_loc_szalo[1,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;O2_ei2=plot(total(m5_eiionz_loc_szamed[1,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;O2_ei3=plot(total(m5_eiionz_loc_szahi[1,*,*],3),zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;
;O2_ei4=plot(total(m5_eiionz_transp_szalo[1,*,*],3),zz,name=leg_nam[0],$
;    thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O2_ei5=plot(total(m5_eiionz_transp_szamed[1,*,*],3),zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O2_ei6=plot(total(m5_eiionz_transp_szahi[1,*,*],3),zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;leg3 = LEGEND(TARGET=[O2_ei1,O2_ei2,O2_ei3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;O2_ei7=plot(total(x9_eiionz_loc_szalo[1,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;O2_ei8=plot(total(x9_eiionz_loc_szamed[1,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O2_ei9=plot(total(x9_eiionz_loc_szahi[1,*,*],3),zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;
;
;O2_ei10=plot(total(x9_eiionz_transp_szalo[1,*,*],3),zz,name=leg_nam[3],$
;    thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O2_ei11=plot(total(x9_eiionz_transp_szamed[1,*,*],3),zz,name=leg_nam[4],$
;    thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;O2_ei12=plot(total(x9_eiionz_transp_szahi[1,*,*],3),zz,name=leg_nam[5],$
;    thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;
;leg4 = LEGEND(TARGET=[O2_ei7,O2_ei8,O2_ei9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;_________________________________________________________________________________________________
;
;; N2 photoelectron ionisation
;
;w6= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;
;N2_ei1=plot(total(m5_eiionz_loc_szalo[2,*,*],3),zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;N2_ei2=plot(total(m5_eiionz_loc_szamed[2,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;N2_ei3=plot(total(m5_eiionz_loc_szahi[2,*,*],3),zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;
;N2_ei4=plot(total(m5_eiionz_transp_szalo[2,*,*],3),zz,name=leg_nam[0],$
;  thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;N2_ei5=plot(total(m5_eiionz_transp_szamed[2,*,*],3),zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;N2_ei6=plot(total(m5_eiionz_transp_szahi[2,*,*],3),zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;leg5 = LEGEND(TARGET=[N2_ei1,N2_ei2,N2_ei3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;N2_ei7=plot(total(x9_eiionz_loc_szalo[2,*,*],3),zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;N2_ei8=plot(total(x9_eiionz_loc_szamed[2,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;N2_ei9=plot(total(x9_eiionz_loc_szahi[2,*,*],3),zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;N2_ei10=plot(total(x9_eiionz_transp_szalo[2,*,*],3),zz,name=leg_nam[3],$
;  thick=thck,color=colors[0],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;
;
;N2_ei11=plot(total(x9_eiionz_transp_szamed[2,*,*],3),zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;N2_ei12=plot(total(x9_eiionz_transp_szahi[2,*,*],3),zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],symbol= sym[0],sym_thick=sym_thck,sym_increment=sym_incr,/overplot)
;
;
;leg6 = LEGEND(TARGET=[N2_ei7,N2_ei8,N2_ei9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;___________________________________________________________________________________________________________________________________________________________________________________

;Pe/Pi ratio plot

;
;
;yr=[min(zz),max(zz)]
;xt='Pe/Pi '
;yt='Altitude (km)'
;titl='Pe/Pi'
;
;
;
;
;;  O Pe/Pi plot
;
;xr=[1.e-3,1.e3]
;
;w7 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;O_pp1=plot(m5_pepi_szalo[0,*],zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;O_pp2=plot( m5_pepi_szamed[0,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;O_pp3=plot(m5_pepi_szahi[0,*],zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;
;
;leg1 = LEGEND(TARGET=[O_pp1,O_pp2,O_pp3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;
;
;O_pp7=plot(x9_pepi_szalo[0,*],zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;O_pp8=plot(x9_pepi_szamed[0,*],zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O_pp9=plot(x9_pepi_szahi[0,*],zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;
;
;
;leg2 = LEGEND(TARGET=[O_pp7,O_pp8,O_pp9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;_________________________________________________________________________
;
;;O2 Pe/Pi plot
;
;w8 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;;xr=[1.e-9,1.e4]
;
;O2_pp1=plot(m5_pepi_szalo[1,*],zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;O2_pp2=plot(m5_pepi_szamed[1,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;O2_pp3=plot(m5_pepi_szahi[1,*],zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;
;leg3 = LEGEND(TARGET=[O2_pp1,O2_pp2,O2_pp3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;O2_pp7=plot(x9_pepi_szalo[1,*],zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;O2_pp8=plot(x9_pepi_szamed[1,*],zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;O2_pp9=plot(x9_pepi_szahi[1,*],zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;
;
;
;leg4 = LEGEND(TARGET=[O2_pp7,O2_pp8,O2_pp9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;_________________________________________________________________________________________________
;
;; N2 Pe/Pi plot
;
;w9 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;
;N2_pp1=plot(m5_pepi_szalo[2,*],zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)
;
;
;N2_pp2=plot(m5_pepi_szamed[2,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;N2_pp3=plot(m5_pepi_szahi[2,*],zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;leg5 = LEGEND(TARGET=[N2_pp1,N2_pp2,N2_pp3], POSITION=leg_loc_left,font_size=fnt_sz,/normal)
;
;
;N2_pp7=plot(x9_pepi_szalo[2,*],zz,name=leg_nam[3],$
;  xtitle=xt,ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)
;
;
;N2_pp8=plot(x9_pepi_szamed[2,*],zz,name=leg_nam[4],$
;  thick=thck,color=colors[1],/overplot)
;
;N2_pp9=plot(x9_pepi_szahi[2,*],zz,name=leg_nam[5],$
;  thick=thck,color=colors[2],/overplot)
;
;leg6 = LEGEND(TARGET=[N2_pp7,N2_pp8,N2_pp9], POSITION=leg_loc_right,font_size=fnt_sz,/normal)



;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;;Comparing new NRL bins results to SQ results
;
leg_name =['NRL','SQ05']

;Photoionisation plot

yt='Altitude (km)'
xt='Photoionisation (Pi) (cm!U-3!N s!U-1!N) !C '
titl='Photoionization'
;yr=[min(zz),max(zz)]
xr=[10.,10000.]


;  O Photoionisation plot


w10 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O_pi1=plot(total(m5_photoi_szalo[0,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O_pi2=plot(total(SQm5_szalo_photoi[0,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg7 = LEGEND(TARGET=[O_pi1,O_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


O_pi4=plot(total(x9_photoi_szalo[0,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O_pi5=plot(total(SQx9_szalo_photoi[0,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg8 = LEGEND(TARGET=[O_pi4,O_pi5], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Photoionisation plot



w11 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])



O2_pi1=plot(total(m5_photoi_szalo[1,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O2_pi2=plot(total(SQm5_szalo_photoi[1,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg9 = LEGEND(TARGET=[O2_pi1,O2_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


O2_pi4=plot(total(x9_photoi_szalo[1,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O2_pi5=plot(total(SQx9_szalo_photoi[1,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg10 = LEGEND(TARGET=[O2_pi4,O2_pi5], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w2.Save, plt_loc+"O2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;N2 Photoionisation plot


w11 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])


N2_pi1=plot(total(m5_photoi_szalo[2,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

N2_pi2=plot(total(SQm5_szalo_photoi[2,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg11 = LEGEND(TARGET=[N2_pi1,N2_pi2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


N2_pi4=plot(total(x9_photoi_szalo[2,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

N2_pi5=plot(total(SQx9_szalo_photoi[2,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg12 = LEGEND(TARGET=[N2_pi4,N2_pi5], POSITION=leg_loc_right,font_size=fnt_sz,/normal)

;w3.Save, plt_loc+"N2_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;___________________________________________________________________________________________________________________________________________________________________________________



;Photoelectron Ionisation plot



;yr=[min(zz),max(zz)]
xr=[10.,1.e5]
;yt='Altitude (km)'
xt='Photoelectron Ionisation (Pe) (cm!U-3!N s!U-1!N) !C '
titl='Photoelectron Ionization'


;  O Photoelectron ionisation plot

w12 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O_ei1=plot(total(m5_eiionz_transp_szalo[0,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O_ei2=plot(total(SQm5_szalo_eiionz_transp[0,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg13 = LEGEND(TARGET=[O_ei1,O_ei2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)



O_ei7=plot(total(x9_eiionz_transp_szalo[0,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O_ei8=plot(total(SQx9_szalo_eiionz_transp[0,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg14 = LEGEND(TARGET=[O_ei7,O_ei8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Photoelectron ionisation plot

w13 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])

O2_ei1=plot(total(m5_eiionz_transp_szalo[1,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O2_ei2=plot(total(SQm5_szalo_eiionz_transp[1,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg15 = LEGEND(TARGET=[O2_ei1,O2_ei2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


O2_ei7=plot(total(x9_eiionz_transp_szalo[1,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O2_ei8=plot(total(SQx9_szalo_eiionz_transp[1,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg16 = LEGEND(TARGET=[O2_ei7,O2_ei8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;
;;_________________________________________________________________________________________________
;
; N2 photoelectron ionisation

w14= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])



N2_ei1=plot(total(m5_eiionz_transp_szalo[2,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

N2_ei2=plot(total(SQm5_szalo_eiionz_transp[2,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg17 = LEGEND(TARGET=[N2_ei1,N2_ei2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


N2_ei7=plot(total(x9_eiionz_transp_szalo[2,*,*],3),zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

N2_ei8=plot(total(SQx9_szalo_eiionz_transp[2,*,*],3),zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg18 = LEGEND(TARGET=[N2_ei7,N2_ei8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)
;
;
;;___________________________________________________________________________________________________________________________________________________________________________________

;Pe/Pi ratio plot


xt='Pe/Pi !C  '
;yt='Altitude (km)'
titl='Pe/Pi'
xr=[0.01,1000.]



;  O Pe/Pi plot


w15 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])


O_pp1=plot( m5_pepi_szalo[0,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O_pp2=plot(SQm5_szalo_pepi[0,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg19 = LEGEND(TARGET=[O_pp1,O_pp2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)



O_pp7=plot(x9_pepi_szalo[0,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[0]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O_pp8=plot(SQx9_szalo_pepi[0,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg20 = LEGEND(TARGET=[O_pp7,O_pp8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________

;O2 Pe/Pi plot

w16 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])


O2_pp1=plot(m5_pepi_szalo[1,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

O2_pp2=plot(SQm5_szalo_pepi[1,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg21 = LEGEND(TARGET=[O2_pp1,O2_pp2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


O2_pp7=plot(x9_pepi_szalo[1,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[1]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

O2_pp8=plot(SQx9_szalo_pepi[1,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg22 = LEGEND(TARGET=[O2_pp7,O2_pp8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)


;w1.Save, plt_loc+"O_photoi_0.5-18.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;_________________________________________________________________________________________________

; N2 Pe/Pi plot

w17 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])


N2_pp1=plot(m5_pepi_szalo[2,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[0],ytitle=yt,title=tit_fl[0]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posleft,/current)

N2_pp2=plot(SQm5_szalo_pepi[2,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)

leg23 = LEGEND(TARGET=[N2_pp1,N2_pp2], POSITION=leg_loc_left,font_size=fnt_sz,/normal)


N2_pp7=plot(x9_pepi_szalo[2,*],zz,name=leg_name[0],$
  xtitle=xt+ leg_nam[3],ytitle=yt,title=tit_fl[1]+maj[2]+titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],font_size=fnt_sz,position=posright,/current)

N2_pp8=plot(SQx9_szalo_pepi[2,*],zz,name=leg_name[1],$
  thick=thck,color=colors[1],/overplot)
  
leg24 = LEGEND(TARGET=[N2_pp7,N2_pp8], POSITION=leg_loc_right,font_size=fnt_sz,/normal)



end





;Seems like there are enormous differences between the SQ binning and ours in the photoelectron ionization rate. 
;Like a factor of 10 at some altitudes, like between 70-100 km. And yet, the differences at these same altitudes in 
;the N2 photoionization rates are small. Why is that? Are we sure that we are getting the same pe/pi as SQ05?
; Ultimately, we will have to plot pe/pi on an optical depth vertical axis like Figure 3 of their paper.
;  They then picked off the value at tau = 1 as their value for the table. We need to show that we can reproduce that 
;  figure and indeed do it for each bin of SQ to show that we can reproduce their table A2-A4.
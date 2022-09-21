
@localpe_setup_exsect
@ace_etransport
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software

; common localpe,nbins,nmaj,nei,nst,jmax,lmax,sigs,pe,pi,siga,sec,sigex,sigix,tpot,prob,$
;        ener,del,peflux,primary,sigloss,iimax,$
;        wv1,wv2,ssflux,$
;       sigabs,sigionx,auger_energy,auger_wvln,$
;        zmaj,zz,z,zt,$
;        lat,sza,idate,utsec,f107,f107a,ap,$
;        first_neutral,first_ssflux,first_pxsect,first_exsect, $
;        zcol,tau,flux,photoki,photoi,aprod,aloss,$
;        ww,a0,omeg,anu,bb,auto,thi,ak,aj,gams,gamb,ts,tb, ta, eistates
        

idate=1999266
lat=0.

f107=137.8
f107a=157.077
ap=28

localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi
endif
localpe_ssflux

dum=min(abs(ener-20.),ia20)
dum=min(abs(ener-100.),ia100)
dum=min(abs(ener-1000.),ia1000)

n_hours=14
glow20=fltarr(n_hours)
glow100=fltarr(n_hours)
glow1000=fltarr(n_hours)
ace20=fltarr(n_hours)
ace100=fltarr(n_hours)
ace1000=fltarr(n_hours)
asza=fltarr(n_hours)

for i=0,n_hours-1 do begin

hour=float(i)+6.
utsec=3600.*hour
localpe_neutralatm
asza[i]=sza
 localpe_photoionz
 localpe_approxeden,zz,f107,eden,etemp
 eden=eden/8.0
localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo
toaflux=fltarr(nbins)
;ace_etransport,toaflux,etemp,eden,upflux,downflux
ace_etransport
tflux=(upflux+downflux)/0.5

z1=110.
z2=250.
dum1l=min(abs(local_zz-z1),ilz1)
dum1l=min(abs(local_zz-z2),ilz2)

cs=2
ct=3

window,i
plot,local_ener,local_pespec[ilz1,*],xr=[1,5e4],/yl,yr=[1e-10,1e12],/xl,psym=10,/nodata,$
charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg

oplot,local_ener,local_pespec[ilz1,*],psym=10,color=!ct.navy,thick=ct*2
oplot,local_ener,local_pespec[ilz2,*],psym=10,color=!ct.navy,thick=ct
oplot,local_ener,tflux[ilz1,*],psym=10,color=!ct.blue,thick=ct*2
oplot,local_ener,tflux[ilz2,*],psym=10,color=!ct.blue,thick=ct

if i eq 0 then begin

; get the original glow run
atm_folder='D:\Users\bailey\Documents\localpe\msisgfiles\'
str_lat=strtrim(fix(abs(lat)),2)
sign_lat='N'
if lat lt 0. then sign_lat='S'
str_idate=strtrim(long(idate),2)
peflux_file=atm_folder+str_idate+'_'+str_lat+sign_lat+'_EUVAC_pespecmt_z.dat'

openr,lun,peflux_file,/get_lun
line=0
s='a string'
glow_nbins=0
for ihour=0,5 do begin
;read the # of energy bins (nbins) and # of alt bins (jmax)
readf,lun,glow_nbins,glow_jmax
glow_ener=(glow_del=(glow_uflx=(glow_phitop=(fltarr(glow_nbins)))))
glow_zz=fltarr(glow_jmax)

;read the electron energies (in eV)
readf,lun,glow_ener
;read the width of each energy bin
readf,lun,glow_del
;read the altitudes (in cm)
readf,lun,glow_zz
a=fltarr(glow_nbins)

dum=min(abs(glow_ener-20.),ig20)
dum=min(abs(glow_ener-100.),ig100)
dum=min(abs(glow_ener-1000.),ig1000)

;read the flux spectrum for each altitude
;each row of zflx is the spectra for a different altitude
glow_zflx=fltarr(glow_jmax,glow_nbins)
for ii=0,glow_jmax-1 do begin
    readf,lun,a
    glow_zflx(ii,*)=a
endfor
;read the flux at the highest altitude
readf,lun,glow_uflx
;read the flux escaping from the top
readf,lun,glow_phitop

endfor

endif 

readf,lun,glow_nbins,glow_jmax
readf,lun,glow_ener
readf,lun,glow_del
readf,lun,glow_zz
glow_zflx=fltarr(glow_jmax,glow_nbins)
for ii=0,glow_jmax-1 do begin
    readf,lun,a
    glow_zflx(ii,*)=a
end
;read the flux at the highest altitude
readf,lun,glow_uflx
;read the flux escaping from the top
readf,lun,glow_phitop

;change zz to km
glow_zz=glow_zz*1e-5
dum1g=min(abs(glow_zz-z1),igz1)
dum2g=min(abs(glow_zz-z2),igz2)

;neglect any negative fluxes
glow_zflx=glow_zflx>1e-15

oplot,glow_ener,glow_zflx[igz1,*],psym=10,color=!snoe.ct.red,thick=ct*2
oplot,glow_ener,glow_zflx[igz2,*],psym=10,color=!snoe.ct.red,thick=ct

glow20[i]=glow_zflx[igz1,ig20]
glow100[i]=glow_zflx[igz1,ig100]
glow1000[i]=glow_zflx[igz1,ig1000]
ace20[i]=tflux[ilz1,ia20]
ace100[i]=tflux[ilz1,ia100]
ace1000[i]=tflux[ilz1,ia1000]

endfor ; hour

window,13
hour=findgen(12)+6.
plot,hour,glow20,/nodata,xr=[0,24],/yl,yr=[1e0,1e9],xtitle='Hour',ytitle='Electron Flux'
oplot,hour,glow20,color=!ct.red
oplot,hour,glow100,color=!ct.red
oplot,hour,glow1000,color=!ct.red
oplot,hour,ace20,color=!ct.navy
oplot,hour,ace100,color=!ct.navy
oplot,hour,ace1000,color=!ct.navy


free_lun,lun

end

; leftovers
;window,1
;plot,local_ener,local_pespec[ilz1,*],xr=[10,30],yr=[0,3e8],psym=10,/nodata,$
;charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg
;oplot,local_ener,local_pespec[ilz2,*],psym=10,color=!ct.navy,thick=ct
;oplot,local_ener,tflux[ilz2,*],psym=10,color=!ct.blue,thick=ct
;oplot,glow_ener,glow_zflx[igz2,*],psym=10,color=!snoe.ct.red,thick=ct

;window,2
;plot,local_ener,local_pespec[ilz1,*],xr=[0,30],yr=[0,3e7],psym=10,/nodata,$
;charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg
;oplot,local_ener,local_pespec[ilz1,*],psym=10,color=!ct.navy,thick=ct
;oplot,local_ener,tflux[ilz1,*],psym=2,color=!ct.blue,thick=ct
;oplot,glow_ener,glow_zflx[igz1,*],psym=10,color=!snoe.ct.red,thick=ct

;window,3
;plot,local_ener,local_pespec[ilz1,*],xr=[300,600],/yl,yr=[1e2,1e5],psym=10,/nodata,$
;charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg

;oplot,local_ener,local_pespec[ilz1,*],psym=10,color=!ct.navy,thick=ct*2
;oplot,local_ener,local_pespec[ilz2,*],psym=10,color=!ct.navy,thick=ct
;oplot,local_ener,tflux[ilz1,*],psym=10,color=!ct.blue,thick=ct*2
;oplot,local_ener,tflux[ilz2,*],psym=10,color=!ct.blue,thick=ct
;oplot,glow_ener,glow_zflx[igz1,*],psym=10,color=!snoe.ct.red,thick=ct*2
;oplot,glow_ener,glow_zflx[igz2,*],psym=10,color=!snoe.ct.red,thick=ct
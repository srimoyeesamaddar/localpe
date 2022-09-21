
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
        

idate=1999072
lat=0.

f107=142.7
f107a=128.551
ap=6

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



hour=12.
utsec=3600.*hour
localpe_neutralatm
;asza[i]=sza
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

window,0
plot,local_ener,local_pespec[ilz1,*],xr=[1,5e4],/yl,yr=[1e-10,1e12],/xl,psym=10,/nodata;,$
;charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg

oplot,local_ener,local_pespec[ilz1,*],psym=10,color=!ct.navy,thick=ct*2
oplot,local_ener,local_pespec[ilz2,*],psym=10,color=!ct.navy,thick=ct
oplot,local_ener,tflux[ilz1,*],psym=10,color=!ct.blue,thick=ct*2
oplot,local_ener,tflux[ilz2,*],psym=10,color=!ct.blue,thick=ct

end


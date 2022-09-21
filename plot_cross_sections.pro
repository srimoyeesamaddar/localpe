 common localpe,nbins,nmaj,nei,nst,jmax,lmax,sigs,pe,pi,siga,sec,sigex,sigix,tpot,prob,$
        ener,del,peflux,primary,$
        wv1,wv2,ssflux,$
       sigabs,sigionx,auger_energy,auger_wvln,$
        zmaj,zz,z,zt,$
        lat,sza,idate,utsec,f107,f107a,ap,$
        first_neutral,first_ssflux,first_pxsect,first_exsect
 
 window,3
 ct=3
 cs=2
 plot,ener,sigix[0,0,*],xtitle='Electron Energy (eV)',ytitle='Cross Section (cm!u2!n)',$
 charsize=cs,charthick=ct,xthick=ct,ythick=ct,/nodata,back=!snoe.ct.fg,color=!snoe.ct.bg,$
 /xlog,/ylog,yr=[1e-20,1e-15],xr=[1,1e5]
 
 for imaj=0,nmaj-1 do begin
 for istate=0,nei-1 do begin
 
 if max(sigix[istate,imaj,*]) gt 0. then begin
 oplot,ener,sigix[istate,imaj,*],color=!snoe.ct.navy,thick=ct
 endif
 
 endfor
 endfor
 
 end
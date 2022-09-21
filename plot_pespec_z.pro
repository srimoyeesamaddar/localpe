

;Filename: read_pespec_z.pro
;Author: Melissa Rice
;Date: Summer 2003

; Renamed plot_pespec_z.pro, SMB July 8, 2008

;Purpose: This program reads the output file 'pespec_z.dat' from
;msisg_peflux and plots photoelectron spectra at four
;different altitudes

;open the data file
openr,lun,'pespec_z.dat',/get_lun
line=0
s='a string'
nbins=0
;read the # of energy bins (nbins) and # of alt bins (jmax)
readf,lun,nbins,jmax
ener=(del=(uflx=(phitop=(fltarr(nbins)))))
zz=fltarr(jmax)
;read the electron energies (in eV)
readf,lun,ener
;read the width of each energy bin
readf,lun,del
;read the altitudes (in cm)
readf,lun,zz
a=fltarr(nbins)

;read the flux spectrum for each altitude
;each row of zflx is the spectra for a different altitude
zflx=fltarr(jmax,nbins)
for i=0,jmax-1 do begin
    readf,lun,a
    zflx(i,*)=a
end
;read the flux at the highest altitude
readf,lun,uflx
;read the flux escaping from the top
readf,lun,phitop
;change zz to km
zz=zz*1e-5
stop
;neglect any negative fluxes
zflx=zflx>1e-10
free_lun,lun

;create the plot
;set_plot,'ps'
;device,filename='spec.ps'
;device,color=256
title='Photoelectron Spectra for Max Solar Activity (2000-7-12)'
ytitle='Energy (eV)'
xtitle='Flux (electrons/(cm^2 s eV))'
plot,ener,zflx(56,*),/xlog,/ylog,xrange=[3,1e3],yrange=[100,1e11],/xst,/yst,title=title,xtitle=xtitle,ytitle=ytitle,thick=6,color=-252
oplot,ener,zflx(18,*),color=-30,thick=6
oplot,ener,zflx(40,*),color=!snoe.ct.green,thick=6
oplot,ener,zflx(76,*),color=!snoe.ct.red,thick=6
;c=0.
;plot,zz,zflx(*,c)
;xyouts,ener(c),0.,'*'
;xyouts,200,10^9.5,'Energy:',color=-30,size=1.5
xyouts,200,10^9.5,'106 km',color=-30,size=1.5
xyouts,200,10^8.75,'150 km',color=!snoe.ct.green,size=1.5
xyouts,200,10^8.25,color=-252,'203 km',size=1.5
xyouts,200,10^7.25,'303 km',color=!snoe.ct.red,size=1.5
;device,/close

;set_plot,'x'

end

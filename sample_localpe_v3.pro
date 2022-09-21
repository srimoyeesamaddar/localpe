;Running the ACEPE code for each flux bin in the solar spectrum
; for calcualting pe/pi  factors in each solar flux bins


;Main photoelectron calling function

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software

;............Model constants....................
nbins=180
nmaj=3
;lmax=139
;lmax=156
;lmax=123
lmax=22 ;Solomon and Qian 2005
lon=0.0
;altitude bins
z0=40.
zz=findgen(161)*2.+z0
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]


;..............Inputs..................


;idate=1999266
idate=20031104 ; change this
lat=0.

;f107=137.8
;f107a=157.077

f107=100.
f107a=100.
ap=15
;binsize=178


tic

localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi

endif

;..............1 run with all bins for photoionisation....................
;.........................................................................

;maskk=0
;mask=0
photoi_run=fltarr(nmaj,jmax,lmax) 
photoi_run_loc=fltarr(nmaj,jmax,lmax) 
primary_run_loc=fltarr(lmax,jmax,nbins)
primary_run=fltarr(lmax,jmax,nbins)
;localpe_ssflux,maskk=maskk,mask
;
;dum=min(abs(ener-20.),ia20)
;dum=min(abs(ener-100.),ia100)
;dum=min(abs(ener-1000.),ia1000)
;
;n_hours=14
;glow20=fltarr(n_hours)
;glow100=fltarr(n_hours)
;glow1000=fltarr(n_hours)
;ace20=fltarr(n_hours)
;ace100=fltarr(n_hours)
;ace1000=fltarr(n_hours)
;asza=fltarr(n_hours)
;
;;for i=0,n_hours-1 do begin
;;  for i=0, 3 do begin
;i=6.
;hour=float(i)+6.
;utsec=3600.*hour
;localpe_neutralatm
;asza[i]=sza
;localpe_photoionz
;
;
;localpe_approxeden,zz,f107,eden,etemp
;;eden=eden/8.0
;localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo
;
;ace_etransport
;
;;photoi_run=fltarr(jmax,nmaj,lmax) ; alt,species,bins
;
;;photoi_run=photoi_wv


;............................................................................

;............................................................................




;set up a mask for solar flux
maskk=1
n_mask=lmax  ; total number of runs for masking(mask=1) one bin in each run
mask=fltarr(lmax,n_mask)
eiionz_run=fltarr(nmaj,jmax,lmax)  ;from transport calculation
eiionzt_run= fltarr(nmaj,jmax,lmax)



for run=0, n_mask-1 do mask[run,run]=1.


for run=0,n_mask-1 do begin  ;n_mask-
  localpe_ssflux,maskk=maskk,reform(mask[*,run])

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

  ;for i=0,n_hours-1 do begin
  ;  for i=0, 3 do begin
  i=6.
  hour=float(i)+6.
  utsec=3600.*hour
  localpe_neutralatm
  asza[i]=sza
  localpe_photoionz


  localpe_approxeden,zz,f107,eden,etemp
  ;eden=eden/8.0
  localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo
 photoi_run_loc[*,*,run]=photoi
 primary_run_loc(run,*,*)=primary

localpe_photoionz
  ace_etransport
   

;   eiionz_run[*,*,run]=total(eiionz_wv,3); local calculation
;   eiionzt_run[*,*,run]=total(eiionz_transp_wv,3) ; transport calculation  jmax,nmaj,lmax
;  
   photoi_run[*,*,run]=photoi 
   primary_run(run,*,*)=primary
   eiionz_run[*,*,run]=eiionz; local calculation
   eiionzt_run[*,*,run]=eiionz_transp ; transport calculation  jmax,nmaj,lmax
  print,'Done', run

endfor
pepi=eiionz_run/photoi_run  ;alt,species,bins
pepi_t=eiionzt_run/photoi_run

;pepi_bins=reform(total(eiionzt_run,1))/reform(total(photoi_run,1)) ;[species (nmaj), bins(lmax)]
;eiionz_bins=reform(total(eiionzt_run,1))
;photoi_bins=reform(total(photoi_run,1))
ssflux_run=ssflux

toc  ;total time elapsed


loc_sav='C:\Users\srimo\Desktop\nrl_files\sav_files\'
fname='Solomon&Qian2005' ;change this each time




save, photoi_run,eiionz_run,eiionzt_run,pepi,pepi_t,zz,ssflux_run,wv1,wv2,tau,$
  filename=loc_sav+fname+'.sav'



;....................................................
;.....................Plots..........................
loc_plt='C:\Users\srimo\Desktop\nrl_files\plots2\'

alt_ind= 10


;cgps_open, filename=loc_plt+fname+'.ps'
;gspath='C:\Program Files\gs\gs9.27'

;cgdisplay,1000,1000, /free
;cgplot,eiionzt_run(alt_ind,0,*)/ photoi_run(alt_ind,0,*), 0.05*(wv1+wv2) ,$
;       /xs,/ys,/xl,yr=[80.,250.],$ ;,xr=[0.001,1.e3]
;       color='red', charthick=2, thick=5.,$
;       xtitle='Photoelectron/Photoionization',$
;       ytitle='Wavelength (nm)',$
;       title='Pe/Pi transport for 1 altitude'
;cgplot,eiionzt_run(alt_ind,0,*)/ photoi_run(alt_ind,0,*), 0.05*(wv1+wv2),/overplot,$
;       color='charcoal',thick=5.
;cgplot, eiionzt_run(alt_ind,0,*)/ photoi_run(alt_ind,0,*), 0.05*(wv1+wv2),/overplot,$
;        color='dodger blue',thick=5.
;cglegend,titles=['O','O$\down2$','N$\down2$'],$
;         psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
;         alignment=1,/box,Location=[0.9, 0.9]
;         
;cgdisplay,1000,1000, /free
;cgplot,pepi_bins(0,*), 0.05*(wv1+wv2) ,$
;         /xs,/ys,/xl,yr=[80.,250.],$;,xr=[0.001,1.e3]
;         color='red', charthick=2, thick=5.,$
;         xtitle='Photoelectron/Photoionization',$
;         ytitle='Wavelength (nm)',$
;         title='Pe/Pi transport for all altitudes'
;cgplot,pepi_bins(1,*), 0.05*(wv1+wv2),/overplot,$
;         color='charcoal',thick=5.
;cgplot, pepi_bins(2,*), 0.05*(wv1+wv2),/overplot,$
;         color='dodger blue',thick=5.
;cglegend,titles=['O','O$\down2$','N$\down2$'],$
;         psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
;         alignment=1,/box,Location=[0.9, 0.9]
;         
;;cgps2pdf,loc_plt+fname+'.ps',$
;;         loc_plt+fname+'.pdf',$
;;                gs_path=gspath
;;       cgps_close,/pdf
;



end

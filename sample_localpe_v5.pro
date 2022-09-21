;  running ACEPE for all the solar flux in the file:
; 'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\x28_custombins_myver.sav'
  




;Main photoelectron calling function

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software

;............Model constants....................
nbins=180
nmaj=3
;lmax=139
lmax=156
;lmax=123
;lmax=22 ;Solomon and Qian 2005
lon=0.0
;altitude bins
z0=40.
zz=findgen(161)*2.+z0
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]


;..............Inputs..................


;idate=1999266
idate=20031104
lat=0.

;f107=137.8
;f107a=157.077

f107=70.
f107a=70.
ap=15
;binsize=178
pepi=fltarr(3,11,jmax)
ssflux_fi=fltarr(11,lmax)

;.....................................................
floc='C:\Users\srimo\Desktop\nrl_files\sav_files\'
;ploc=
;fname='Solomon&Qian2005_allrun'
fname='X28_334'
;.....................................................

tic

localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi

endif

;DEFAULT solar flask mask inputs
maskk=0  ;deafult values ;make maskk=1 when using one flux bin at a time
mask=0
fi=[10,20,50,90,100,110,120,150,200,220,250]
;Solar flux inputs
for fii=0, 10 do begin
localpe_ssflux,maskk=maskk,mask,fi[fii]

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


;photoelectron flux transport calculation
ace_etransport

pepi[0,fii,*]=eiionz_transp(0,*)/photoi(0,*)
pepi[1,fii,*]=eiionz_transp(1,*)/photoi(1,*)
pepi[2,fii,*]=eiionz_transp(2,*)/photoi(2,*)
ssflux_fi(fii,*)=ssflux

toc
print,fi[fii],' Done '
endfor
save, pepi,zz,ssflux_fi,wv1,wv2, filename=floc+fname+'.sav'
end
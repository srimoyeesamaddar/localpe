;Theoritical calculation of photoionisation and verification of the model cross-section data
;June 29, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'
restore,filename=floc+'peakm5.sav'
;restore, filename=floc+'peakx9.sav'
; file has:photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2
photoi_mod= photoi
photoi_wv_mod=photoi_wv
eiionz_transp_mod=eiionz_transp
zz_mod= zz
ssflux_mod=ssflux
flux_mod=flux
wv1_mod=wv1
wv2_mod=wv2

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________
; Base Cross-section files

restore, '/home/srimoyee/Desktop/nrl_files/sav_files/ion&abs_cross_files.sav' ;from Henke & Fennely
; file has:
;wave_o,wave_o2,wave_n2
;totionz_o,totionz_o2,totionz_n2
;totabs_o,totabs_o2,totabs_n2

wave_o_base = wave_o  
wave_o2_base= wave_o2
wave_n2_base =wave_n2 
totionz_o_base =totionz_o
totionz_o2_base= totionz_o2
totionz_n2_base =totionz_n2
totabs_o_base =totabs_o
totabs_o2_base= totabs_o2
totabs_n2_base =totabs_n2


restore, '/home/srimoyee/Desktop/nrl_files/sav_files/baseline_states.sav'
;  file has:
;    O_crss_states,O2_crss_states,N2_crss_states
;    first column is wavelength

O_crss_states_base =O_crss_states
O2_crss_states_base = O2_crss_states
N2_crss_states_base = N2_crss_states

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model Cross-section files 

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2

restore, probstate_crss
restore, ionabs_crss


O_prob_state_mod = O_prob_state
O2_prob_state_mod= O2_prob_state
N2_prob_state_mod= N2_prob_state

spec_wave_mod =spec_wave
sigi_o_mod= sigi_o
sigi_o2_mod = sigi_o2
sigi_n2_mod= sigi_n2
sigab_o_mod =sigab_o
sigab_o2_mod = sigab_o2
sigab_n2_mod = sigab_n2
;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;  SQ cross-section files

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_SQ05.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec_SQ05.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
;
restore, probstate_crss
restore, ionabs_crss


O_prob_state_SQ = O_prob_state
O2_prob_state_SQ= O2_prob_state
N2_prob_state_SQ= N2_prob_state

spec_wave_SQ =spec_wave
sigi_o_SQ = sigi_o
sigi_o2_SQ = sigi_o2
sigi_n2_SQ= sigi_n2
sigab_o_SQ =sigab_o
sigab_o2_SQ = sigab_o2
sigab_n2_SQ = sigab_n2
;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________


;  Hinteregger Solar Flux

file_hint= '/home/srimoyee/Desktop/bailey_bins/sc21refw.dat'
openr,lun_hint, file_hint, /get_lun
n_hint=1661   ;wavelength bins

hintwvln=dblarr(n_hint)  ;in angstrom
hintflux=dblarr(n_hint)  ;in  photons/cm2/s
hintlable=strarr(n_hint)  ;flux line
hintclass=intarr(n_hint)
r_correlation= fltarr(n_hint)

format='$(f8.2,f8.1,a11,i1,f10.6)'
l=''
for i=0,n_hint-1 do begin
  readf,lun_hint,format,a,b,l,c,d
  hintwvln(i)=a
  hintflux(i)=b*1.e9 /1000.  ; to be in units of photons/cm2/s
  hintlable(i)=l
  hintclass(i)=c
  r_correlation(i)=d
endfor

free_lun,lun_hint
hintener=12397./hintwvln

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;Inputs
z0=40.
zz=findgen(17)*10.+z0
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]

nmaj=3

lat=20.
lon=0.
f107= 70.
f107a =70.
ap=15

idate=2017249   ;x9 flare
i=6  ;i th hour
hour=float(i)+6.
utsec=3600.*hour
doy=idate-(long(idate/1000)*1000)
utdoy=doy+utsec/(3600.*24.) ; ut in days since Jan 1 of year (input to zenith.pro)
sza=zenith(utdoy,lat,lon);105*!dtor;
localhour=(((utsec/240. + lon)/15. + 24. ) mod 24.)

;__________________________________________________________________________________________________________________________________________________________________________________

;#1 Using the new msis beta to calculate the values

;##1 make an file with the inputs to the MSIS model

openw, lun_msisinp,'/home/srimoyee/Desktop/nrl_files/sav_files/msis_inputs.dat',/get_lun
printf,lun_msisinp,'iyd', 'sec', 'number of altitudes','alt','glat','glong',$
  'stl','f107a','f107','ap','mass'
printf,lun_msisinp,idate
printf,lun_msisinp,utsec
printf,lun_msisinp,n_elements(zz)
printf,lun_msisinp,zz
printf,lun_msisinp,lat
printf,lun_msisinp,lon
printf,lun_msisinp,localhour
printf,lun_msisinp,f107a
printf,lun_msisinp,f107
printf,lun_msisinp,ap
printf,lun_msisinp,1.

close, lun_msisinp
;__________________________________________________________________________________________

;##2 make the msis executable file
spawn,$
  "gfortran -ffree-line-length-200 -o msis_idl_interface.exe hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 msis_idl_interface.f90"
;__________________________________________________________________________________________

;##3 call the fortran wrapper using spawn command

msis_filename=['/home/srimoyee/Desktop/localpe/msis_idl_interface.exe']  ; fortran executable file/wrapper. KEEP The MSIS files in the same directory as localpe

spawn,msis_filename,result,/noshell


;__________________________________________________________________________________________

;##4 open the output file generated by fortran MSIS to get the varibles

msis_outdata=dblarr(4,n_elements(zz))  ;4 columns: O, O2 and N2 densitiy and temperature

msis_outfile='/home/srimoyee/Desktop/nrl_files/sav_files/msis_outputs.dat'
openr, lun_msis, msis_outfile, /get_lun
skip_lun, lun_msis, 1, /lines ; 2 header lines
readf, lun_msis , msis_outdata
free_lun, lun_msis

;##5 save the file data into the local variables

zt=msis_outdata(3,*)
zo=msis_outdata(0,*)
zo2=msis_outdata(1,*)
zn2=msis_outdata(2,*)


close,/all

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;Calculating scale height H

z=zz*1.e5
am=[16., 32., 28.] ; amu for O, O2 anf N2
kb= 1.38064852e-16   ; Boltzmann's constant ergs per k

re=6378.e5 ; radius of the Earth in cm
g=980.  ;acceleration due to gravity in cm^2/sec
gr=g*(re/(re+z))^2 ; account for small change in g with altitude

amave=(.2*32 + .79*28 + .01*40.) ; average mass of atmosphere, used for lower altitudes
gperamu=1.66054*1.e-24  ;coverting amu to grams


;    Assumption: at low altitudes, every constituent follows
;                a single scale height, above that, they follow
;                their own, assume change occurs at 100km


ilow=z le 100.*1e5
ihi=1-ilow

H=fltarr(nmaj,n_elements(zz))  ;in km
for i=0,nmaj-1 do begin
  amass=amave*ilow + am[i]*ihi     ;mass for every height

  H[i,*]=(kb*zt/(amass*gperamu*gr))/1.e5


endfor

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________
;Trying to see the flux and photoionisation deviation in 1 bin
i_o=where((wave_o_base ge 320.) and (wave_o_base lt 540.),n_o)

st_o=where((reform(o_crss_states_base[0,*]) ge 320.) and (reform(o_crss_states_base[0,*]) lt 540),n_st)
st_wave= reform(o_crss_states_base[0,st_o])

lmax=n_o
nst =10
lmax_st= n_st

; SIGABS  photoabsorption cross sections,0- O, 1-O2, 2-N2; cm2
; SIGIONx  photoionization cross sections,0- O,1- O2,2- N2; cm2

sigionx=fltarr(nmaj,lmax)
sigabs=fltarr(nmaj,lmax)

sigabs(0,*)=totabs_o_base[i_o]
sigabs(1,*)=totabs_o2_base[i_o]
sigabs(2,*)=totabs_n2_base[i_o]

sigionx(0,*)=totionz_o_base[i_o]
sigionx(1,*)=totionz_o2_base[i_o]
sigionx(2,*)=totionz_n2_base[i_o]


;PROB    branching ratios for each state, species, and wavelength bin:
nst =10


prob=fltarr(nst,nmaj,lmax_st) ;states species and energies, 0-O, 1-O2, 2- N2


;0-O
n_cols=(size(O_prob_state))[1]  ;number of states

for l=0,lmax_st-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,0,l)=O_crss_states_base(k,st_o[l])
endfor

;1-O2
n_cols=(size(O2_prob_state))[1]  ;number of states

for l=0,lmax_st-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,1,l)=O2_crss_states_base(k,st_o[l])
endfor


;2-N2
n_cols=(size(N2_prob_state))[1]  ;number of states

for l=0,lmax_st-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,2,l)=N2_crss_states_base(k,st_o[l])
endfor


;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;Calculating theoritical flux

flux_th = fltarr(jmax,lmax)
tau_th=fltarr(jmax,lmax)
pi_o=fltarr(jmax,lmax)
pi_o2=fltarr(jmax,lmax)
pi_n2=fltarr(jmax,lmax)
pi=fltarr(jmax,lmax,nst,nmaj)


for i=0 , lmax-1 do begin

  tau_th(*,i)= (H(0,*)*sigabs(0,i)*zo + H(1,*)*sigabs(1,i)*zo2 + H(2,*)*sigabs(2,i)*zn2)*1.e5; optical depth
  flux_th(*,i) = ssflux(18)*exp(-1.*tau_th[*,i])  ;flux 320-540 A
  
  min_o=min(abs(st_wave-wave_o_base[i_o[i]]),imin)
  print, st_wave[imin],wave_o_base[i_o[i]]
  
  for j=0, jmax-1 do begin

    for l=0,nst-1 do begin
      pi(j,i,l,0)=flux_th(j,i)*prob(l,0,imin)*zo(j)
      pi(j,i,l,1)=flux_th(j,i)*prob(l,1,imin)*zo2(j)
      pi(j,i,l,2)=flux_th(j,i)*prob(l,2,imin)*zn2(j)


    endfor

    pi_o(j,i)= total(pi(j,i,*,0))
    pi_o2(j,i)=total(pi(j,i,*,1))
    pi_n2(j,i)=total(pi(j,i,*,2))

  endfor




endfor


;***************************************************************************************************************************

;sigionx[0,12:30]=sigi_o[3:-1]
;sigionx[1,12:30]=sigi_o2[3:-1]
;sigionx[2,12:30]=sigi_n2[3:-1]


;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

plt_loc="/home/srimoyee/Desktop/nrl_files/jun30_plots/"

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
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]
thck=3.
ln_thck=5
sym_thck=5
leg_loc=[0.95, 0.70]
leg_loc2=[0.95, 0.60]
leg_loc3=[0.95, 0.26]
leg_loc_left=[0.37, 0.87]
maj=['O ','O!L2!N ','N!L2!N ']
leg_nam=['HF','SQ']

sym=["*","+","tu"]

;_________________________________________________________________________________________________________________

; Plotting the original base cross-section files

xr=[0.,1050.]
yr=[min(sigab_o/1.e-18),60.]
xt='Wavelength ($\AA$)'
yt='Cross-section (cm!U2!N )'
titl='Cross-section'
xr=[17.,660.]

;#O

yr=[0.,20.]
arrow_yr=[yr[1],yr[0]]

w5=window(window_title=titl,dimension=[800,800])

ab_o=plot(wave_o_base, totabs_o_base/1.e-18,$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],/current)

ab_o1=scatterplot(spec_wave_SQ, sigab_o_SQ/1.e-18,$
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

ab_o2=scatterplot(spec_wave_SQ, sigab_o_SQ/1.e-18,$
  sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)

o_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

;w5.Save, plt_loc+"Oabs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;#O2

yr=[0.,30.]
arrow_yr=[yr[1],yr[0]]


w6=window(window_title=titl,dimension=[800,800])

ab_o2=plot(wave_o2_base, totabs_o2_base/1.e-18,$
  xtitle=xt,ytitle=yt,title=maj[1]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],/current)

ab_o21=scatterplot(spec_wave_SQ, sigab_o2_SQ/1.e-18,$
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

ab_o22=scatterplot(spec_wave_mod, sigab_o2_mod/1.e-18,$
  sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)

o2_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)


;w6.Save, plt_loc+"O2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;#N2

yr=[0.,30.]
arrow_yr=[yr[1],yr[0]]

w7=window(window_title=titl,dimension=[800,800])

ab_n2=plot(wave_n2_base, totabs_n2_base/1.e-18, $
  xtitle=xt,ytitle=yt,title=maj[2]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[4],/current)

ab_n21=scatterplot(spec_wave_SQ, sigab_n2_SQ/1.e-18, $
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

ab_n22=scatterplot(spec_wave_mod, sigab_n2_mod/1.e-18, $
  sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)

n2_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)



;w7.Save, plt_loc+"N2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;_____________________________________________________________________________________________________

;  Hinteregger Solar Flux



ii_hint1=where((hintwvln ge 645.) and (hintwvln lt 1027.),n_hint1)
ii_hint2=where((hintwvln ge 798.) and (hintwvln lt 913.),n_hint2)
ii_hint3=where((hintwvln ge 913.) and (hintwvln le 975.),n_hint3)

y_lab1=make_array(n_hint1,value=1.e10)
y_lab2=make_array(n_hint2,value=1.e10)
y_lab3=make_array(n_hint3,value=1.e10)



xr=[17.,660.]
yr=[0.,5.e11]
arrow_yr=[1.e10,1.]

xt='Wavelength ($\AA$)'
yt='Flux (Photons cm!U-2!N s!U-1!N))'
titl='Hinteregger Solar Flux'




w8=window(window_title=titl,dimension=[800,800])

h1=plot(hintwvln, hintflux, $
  xtitle=xt,ytitle=yt,title=titl,$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)

ss1=scatterplot((wv1_mod+wv2_mod)/2., ssflux_mod,$
  sym_color=colors[4],sym_thick=thck,/overplot)

;h_1 = ARROW([650,650], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
;h_2 = ARROW([798,798], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)
;h_3 = ARROW([913,913], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
;h_4 = ARROW([975,975], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)



;t1 = TEXT(hintwvln[ii_hint1],y_lab1,hintlable[ii_hint1], /data,font_size=10,orientation=90,font_color=colors[3])
;
;
;;w8.Save, plt_loc+"hint_flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________



end
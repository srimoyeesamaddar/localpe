;SQ binning reproduction
;Jul 10, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________


;Model input files
restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5

wv1=wave_gcm1
wv2=wave_gcm2
;ssflux=meanpeakx9
ssflux=peakm5

;__________________________________________________________________________________________________________________________________________________________________________________

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




;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________


;SQ'05
wavelength_low=[0.5,4.,8.,18.,32.,70.,155.,224.,290.,320.,540.,650.,650.,798.,798.,798.,913.,$
  913.,913.,975.,987.,1027.]

wavelength_high=[4.,8.,18.,32.,70.,155.,224.,290.,320.,540.,650.,798.,798.,913.,913.,913.,$
  975.,975.,975.,987.,1027.,1050.]

spec_wave=0.5*(wavelength_low+wavelength_high)

; Absorption cross-sections
sigab_o=[0.0023,0.0170,0.1125,0.1050,0.3247,1.319,3.7832,6.0239,7.7205,10.7175,13.1253,8.5159,4.7889,$
  3.0031,4.1048,3.7947,0.0,0.0,0.0,0.0,0.0,0.0]*1e-18  ;in cm^2
sigab_o2=[0.0045,0.034,0.2251,0.2101,0.6460,2.6319,7.6283,13.2125,16.8233,20.3066,27.0314,23.5669,24.9102,$
  10.498,10.9075,13.3122,13.395,14.4042,32.5038,18.7145,1.6320,1.1500]*1e-18  ;in cm^2
sigab_n2=[0.0025,0.0201,0.1409,1.1370,0.3459,1.5273,5.0859,9.9375,11.7383,19.6514,23.0931,23.0346,54.5252,$
  2.1434,13.1062,71.6931,2.1775,14.4390,115.257,2.5465,0.0,0.0]*1e-18  ;in cm^2

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;the original base cross-section files H+F
restore, '/home/srimoyee/Desktop/nrl_files/sav_files/ion&abs_cross_files.sav' ;from Henke & Fennely
; file has:
;wave_o,wave_o2,wave_n2
;totionz_o,totionz_o2,totionz_n2
;totabs_o,totabs_o2,totabs_n2

;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

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

; Interpolating the absorption cross-sections to hinterreger bins 
intp_abs_o=exp(interpol(alog(totabs_o),alog(wave_o),alog(hintwvln)))
intp_abs_o2=exp(interpol(alog(totabs_o2),alog(wave_o2),alog(hintwvln)))
intp_abs_n2=exp(interpol(alog(totabs_n2),alog(wave_n2),alog(hintwvln)))

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

colors=['red','dark green','dodger blue','purple',"dark goldenrod"]
lstyl=[0,2,3,4]
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]
thck=3.
ln_thck=3
sym_thck=5
leg_loc=[0.95, 0.70]
leg_loc2=[0.95, 0.60]
leg_loc3=[0.95, 0.26]
leg_loc_left=[0.37, 0.87]
maj=['O ','O!L2!N ','N!L2!N ']
leg_nam=['HF','SQ','Interpolated']

sym=["*","+","tu"]










;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

; Plotting the original base cross-section files


xt='Wavelength ($\AA$)'
yt='Cross-section (cm!U2!N )'
titl='Cross-section'

xr=[600.,980.]
yr=[1.e-19,1.e-15]


arrow_yr=[1.e-15,1.e-19]

;#O

w1=window(window_title=titl,dimension=[800,800])

ab_o=plot(wave_o, totabs_o, name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$ ;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],/current)
  
;intp_o=plot(hintwvln, intp_abs_o, name=leg_nam[2],$
;            color=colors[0],thick=thck,/overplot)

ab_o1=scatterplot(spec_wave, sigab_o, name=leg_nam[1],$
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)


o_1 = ARROW([650,650], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)
o_2 = ARROW([798,798], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)
o_3 = ARROW([913,913], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)
o_4 = ARROW([975,975], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)

;leg1 = LEGEND(TARGET=[ab_o,intp_o], POSITION=leg_loc_left,/normal)

;w1.Save, plt_loc+"Oabs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;#O2

yr=[1.e-18,1.e-16]
arrow_yr=[1.e-16,1.e-18]


w2=window(window_title=titl,dimension=[800,800])

ab_o2=plot(wave_o2, totabs_o2, name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[1]+titl,$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],/current)
;
;intp_o2=plot(hintwvln, intp_abs_o2, name=leg_nam[2],$
;    color=colors[0],thick=thck,/overplot)

ab_o21=scatterplot(spec_wave, sigab_o2, name=leg_nam[1],$
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

o2_1 = ARROW([650,650], arrow_yr, COLOR=colors[0],line_thick=ln_thck, /DATA, /CURRENT)
o2_2 = ARROW([798,798], arrow_yr, COLOR=colors[0],line_thick=ln_thck, /DATA, /CURRENT)
o2_3 = ARROW([913,913], arrow_yr, COLOR=colors[0],line_thick=ln_thck, /DATA, /CURRENT)
o2_4 = ARROW([975,975], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)

;leg2 = LEGEND(TARGET=[ab_o2,intp_o2], POSITION=leg_loc_left,/normal)


;w2.Save, plt_loc+"O2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;#N2

yr=[1.e-19,1.e-15]
arrow_yr=[1.e-15,1.e-19]

w3=window(window_title=titl,dimension=[800,800])

ab_n2=plot(wave_n2, totabs_n2, name=leg_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[2]+titl,$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[4],/current)
  
;intp_n2=plot(hintwvln, intp_abs_n2, name=leg_nam[2],$
;    color=colors[0],thick=thck,/overplot)


ab_n21=scatterplot(spec_wave, sigab_n2, name=leg_nam[1],$
  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)


n2_1 = ARROW([650,650], arrow_yr, COLOR=colors[0],line_thick=ln_thck, /DATA, /CURRENT)
n2_2 = ARROW([798,798], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)
n2_3 = ARROW([913,913], arrow_yr, COLOR=colors[0],line_thick=ln_thck, /DATA, /CURRENT)
n2_4 = ARROW([975,975], arrow_yr, COLOR=colors[0],line_thick=ln_thck,/DATA, /CURRENT)

;leg3 = LEGEND(TARGET=[ab_n2,intp_n2], POSITION=leg_loc_left,/normal)


;w3.Save, plt_loc+"N2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;  Hinteregger Solar Flux



ii_hint1=where((hintwvln ge 600.) and (hintwvln lt 980.),n_hint1)
ii_hint2=where((hintwvln ge 798.) and (hintwvln lt 913.),n_hint2)
ii_hint3=where((hintwvln ge 913.) and (hintwvln le 975.),n_hint3)

y_lab1=make_array(n_hint1,value=1.e10)
y_lab2=make_array(n_hint2,value=1.e10)
y_lab3=make_array(n_hint3,value=1.e10)



yr=[1.e4,1.e12]

;yr=[1.e4,1.e12]
;      yr=[min(flux[i,*]),max(flux[i,*])]
xt='Wavelength ($\AA$)'
yt='Flux (Photons cm!U-2!N s!U-1!N))'
arrow_yr=[1.e12,1.e4]

titl='Hinteregger Solar Flux'


w4=window(window_title=titl,dimension=[800,800])

h1=plot(hintwvln, hintflux, $
  xtitle=xt,ytitle=yt,title=titl,$
  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)

ss1=scatterplot((wv1+wv2)/2., ssflux,$
  sym_color=colors[4],sym_thick=thck,/overplot)

h_1 = ARROW([650,650], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
h_2 = ARROW([798,798], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)
h_3 = ARROW([913,913], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
h_4 = ARROW([975,975], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)



t1 = TEXT(hintwvln[ii_hint1],y_lab1,hintlable[ii_hint1], /data,font_size=10,orientation=90,font_color=colors[3])


;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;  Hinteregger * absorption cross section Solar Flux
 
xt='Wavelength ($\AA$)'
yt='Flux * Absorption Cross-section(Photons  s!U-1!N))'
titl='Hinteregger Solar Flux * Absorption Cross-section'

yr=[1.e-12,1.e-7]
arrow_yr=[1.e-7,1.e-12]

w5=window(window_title=titl,dimension=[800,800])

hab1=plot(hintwvln, intp_abs_o*hintflux, $
  xtitle=xt,ytitle=yt,title=titl+ "-"+ maj[0],$
  xrange=xr, yrange=yr,ylog=1,ystyle=1,xlog=0,xstyle=1,$  
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],layout=[1,3,1],/current)


h_1 = ARROW([650,650], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
h_2 = ARROW([798,798], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)
h_3 = ARROW([913,913], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
h_4 = ARROW([975,975], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)



;w6=window(window_title=titl,dimension=[800,800])

hab1=plot(hintwvln, intp_abs_o2*hintflux, $
  xtitle=xt,ytitle=yt,title=titl+ "-"+ maj[1],$
  xrange=xr, yrange=yr,ylog=1,ystyle=1,xlog=0,xstyle=1,$  
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=0,color=colors[1],linestyle=lstyl[0],layout=[1,3,2],/current)





;w7=window(window_title=titl,dimension=[800,800])

hab1=plot(hintwvln, intp_abs_n2*hintflux, $
  xtitle=xt,ytitle=yt,title=titl+ "-"+ maj[2],$
  xrange=xr, yrange=yr,ylog=1,ystyle=1,xlog=0,xstyle=1,$  
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,histogram=0,color=colors[2],linestyle=lstyl[0],layout=[1,3,3],/current)









end
;trying out new bins for solar flux


sav_loc='/home/srimoyee/Desktop/nrl_files/sav_files/'
filename='ion&abs_cross_files.sav'
;file has: wave_o,wave_o2,wave_n2,totionz_o,totionz_o2,totionz_n2,totabs_o,totabs_o2,totabs_n2
restore, sav_loc+filename

vert=7.;vertical altitude resolution in km
H= 7. ;scale height in km

n=n_elements(where(wave_o le 600.))  ; doing this for wavelengths less than equal to 500 angstrom
;n=n_elements(wave_o)

;_________________________________________________________________________________________________________________________________________________
;#1st method
;ln(sigma2(lambda2))-ln(sigma1(lambda1))= vert/H  ;calculating in ln scale
;
;logabs_o=alog(totabs_o[0])+findgen(n)*(vert/H)
;logabs_o2=alog(totabs_o2[0])+findgen(n)*(vert/H)
;logabs_n2=alog(totabs_n2[0])+findgen(n)*(vert/H)
;
;abs_o=exp(logabs_o)  
;abs_o2=exp(logabs_o2)
;abs_n2=exp(logabs_n2)
;
;
;o_ind=where(wave_o le 600.)
;o2_ind=where(wave_o2 le 600.)
;n2_ind=where(wave_n2 le 600.)
;
;tababs_o = fltarr(n)
;tababs_o2 = fltarr(n)
;tababs_n2 = fltarr(n)
;
;wv_o = fltarr(n)
;wv_o2 = fltarr(n)
;wv_n2 = fltarr(n)
;
;for i=1, n-1 do begin
;  r_o=min(abs(totabs_o[o_ind]-abs_o[i]),oi)
;  r_o2=min(abs(totabs_o2[o2_ind]-abs_o2[i]),o2i)
;  r_n2=min(abs(totabs_n2[n2_ind]-abs_n2[i]),n2i)
;  
;  tababs_o[i]=totabs_o[o_ind[oi]]]
;  tababs_o2[i]=totabs_o2[o2_ind[o2i]]]
;  tababs_n2[i]=totabs_n2[n2_ind[n2i]]]
;  
;  wv_o[i]=wave_o[o_ind[oi]]
;  wv_o2[i]=wave_o2[o2_ind[o2i]]
;  wv_n2[i]=wave_n2[n2_ind[n2i]]
;  
;endfor
;tababs_o[0]=totabs_o[0]   ;setting the first bin for the calulated and table data to be equal
;tababs_o2[0]=totabs_o2[0]
;tababs_n2[0]=totabs_n2[0]
;
;wv_o[0]=wave_o[0]
;wv_o2[0]=wave_o2[0]
;wv_n2[0]=wave_n2[0]


;_________________________________________________________________________________________________________________________________________________

; 2nd method 
;sigma2(lambda2)/sigma1(lambda1)= exp(vert/H)  ;calculating in ln scale

abs_o=(exp(vert/H))^findgen(n)*totabs_o[0]
abs_o2=(exp(vert/H))^findgen(n)*totabs_o2[0]
abs_n2=(exp(vert/H))^findgen(n)*totabs_n2[0]


o_ind=where(wave_o le 600.)
o2_ind=where(wave_o2 le 600.)
n2_ind=where(wave_n2 le 600.)

tababs_o = fltarr(n)
tababs_o2 = fltarr(n)
tababs_n2 = fltarr(n)

wv_o = fltarr(n)
wv_o2 = fltarr(n)
wv_n2 = fltarr(n)

for i=1, n-1 do begin
  r_o=min(abs(totabs_o[o_ind]-abs_o[i]),oi)
  r_o2=min(abs(totabs_o2[o2_ind]-abs_o2[i]),o2i)
  r_n2=min(abs(totabs_n2[n2_ind]-abs_n2[i]),n2i)
  
  tababs_o[i]=totabs_o[o_ind[oi]]
  tababs_o2[i]=totabs_o2[o2_ind[o2i]]
  tababs_n2[i]=totabs_n2[n2_ind[n2i]]

  wv_o[i]=wave_o[o_ind[oi]]
  wv_o2[i]=wave_o2[o2_ind[o2i]]
  wv_n2[i]=wave_n2[n2_ind[n2i]]

endfor

tababs_o[0]=totabs_o[0]   ;setting the first bin for the calulated and table data to be equal
tababs_o2[0]=totabs_o2[0]
tababs_n2[0]=totabs_n2[0]

wv_o[0]=wave_o[0]
wv_o2[0]=wave_o2[0]
wv_n2[0]=wave_n2[0]
;_________________________________________________________________________________________________________________________________________________


; OLD/IDL  msis model for n0 value 
idate=2016205
lat=0.
lon=0.
f107a=80.
f107=80.
ap=15
i=6
hour=float(i)+6.
utsec=3600.*hour
doy=idate-(long(idate/1000)*1000)
utdoy=doy+utsec/(3600.*24.) ; ut in days since Jan 1 of year (input to zenith.pro)
localhour=(((utsec/240. + lon)/15. + 24. ) mod 24.)
zz=60. ;in km
;#1 OLD/IDL  msis method
;
; function MSIS2K,Zin,YYYYDDD,UTSEC,LAT,LON,SLTHR,F10AVG,F10,APINDX,$
;                mass=massno,meters=meters,switches=switches,drag=drg,$
;                t_inf_in=t_inf_in  ,s_bates_in=s_bates_in  ,t_lb_in=t_lb_in,$
;                t_inf_out=t_inf_out,s_bates_out=s_bates_out,t_lb_out=t_lb_out

msis=msis2k(zz,idate, utsec,lat,lon,localhour,f107a,f107,ap)
zo=double(msis.o[0])
zo2=double(msis.o2[0])
zn2=double(msis.n2[0])

n0=2.7e19
;________________________________________________________________________________________________________________________________________________
;Height of maximum absorption: z=z0+H*ln(H*n0*sigma*sec(chi)) 
chi=0.*!dtor ;zenith angle, FOR NOW NOT USING Zenith function
z0=0  ; FOR NOW; assuming density is constant below zz
;zmax=z0+H*alog(H*1.e5*(zo*abs_o+zo2*abs_o2+zn2*abs_n2)/cos(chi))
 
zmax=z0+H*alog(H*1.e5*n0*abs_n2/cos(chi))
 
;________________________________________________________________________________________________________________________________________________
  
head_var=['calculated cross-section (cm2)','H+F cross-section','H+F wavelength','H+F energy','Altitude(Tau=1)']
tab_headO= 'O data for Vertical resolution of ' + strtrim(string(vert))+ ' km'
tab_headO2= 'O2 data for Vertical resolution of ' + strtrim(string(vert))+ ' km'
tab_headN2= 'N2 data for Vertical resolution of ' + strtrim(string(vert))+ ' km'

;output data file
outfile1=sav_loc+'vert7_newbins_gen_o_v2.csv'
outfile2=sav_loc+'vert7_newbins_gen_o2_v2.csv'
outfile3=sav_loc+'vert7_newbins_gen_n2_v2.csv'



WRITE_CSV, outfile1, abs_o,tababs_o,wv_o,12397./wv_o, zmax,$
  HEADER=head_var,Table_header=tab_headO

WRITE_CSV, outfile2, abs_o2,tababs_o2,wv_o2,12397./wv_o2, zmax, $
  HEADER=head_var,Table_header=tab_headO2

WRITE_CSV, outfile3,abs_n2,tababs_n2,wv_n2,12397./wv_n2, zmax, $
    HEADER=head_var,Table_header=tab_headN2


;________________________________________________________________________________________________________
;_____________________________________Plots______________________________________________________________

plt_loc='/home/srimoyee/Desktop/nrl_files/nrl_newbins_gen/'

;Universal Plot variables
xg=0 ;gridstyles
yg=0
xtl=1.  ;ticklen
ytl=1.
xsg=1
ysg=1
xstl=1.
ystl=1.

colors=['red','dark green','dodger blue','purple']
lstyl=[0,2,3,4]
sym=["star","diamond"]
sym_fil=0
sym_thck=3.
sym_col=['dodger blue','yellow','purple','red']
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]
thck=5.
leg_loc=[0.90, 0.85]
leg_loc_left=[0.37, 0.87]
maj=['O ','O!L2!N ','N!L2!N ']


head_titl='Henke +Fennely '
;_____________________________________________________________________________________________________

;xt='Wavelength ($\AA$)'
;yt='Cross-sections (cm!U2!N))'
;
;;O absorption cross-section
;
;;xr=[min(wave_o),max(wave_o)]
;xr=[min(wave_o),600.]
;;yr=[min(totabs_o),max(totabs_o)]
;yr=[0.,2.e-17]
;titl=maj[0]+ 'Absorption Cross-sections'
;
;
;w1=window(window_title=head_titl,dimension=[800,800])
;
;abc1=plot(wave_o, totabs_o, name=maj[0],$
;         xtitle=xt,ytitle=yt,title=titl,$
;         xrange=xr, yrange=yr,ylog=0,ystyle=1,xlog=0,xstyle=1,$;
;         xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;         thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)
;
;w1.Save, plt_loc+"o_abscross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;O2 absorption cross-section
;
;;xr=[min(wave_o2),max(wave_o2)]
;xr=[min(wave_o2),600.]
;;yr=[min(totabs_o2),max(totabs_o2)]
;yr=[0.,3.e-17]
;titl=maj[1]+ 'Absorption Cross-sections'
;
;
;w2=window(window_title=head_titl,dimension=[800,800])
;
;abc2=plot(wave_o2, totabs_o2, name=maj[1],$
;         xtitle=xt,ytitle=yt,title=titl,$
;         xrange=xr, yrange=yr,ylog=0,ystyle=1,xlog=0,xstyle=1,$
;         xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;         thick=thck,histogram=0,color=colors[1],linestyle=lstyl[0],/current)
;
;w2.Save, plt_loc+"o2_abscross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;
;;N2 absorption cross-section
;
;;xr=[min(wave_n2),max(wave_n2)]
;xr=[min(wave_n2),600.]
;;yr=[min(totabs_n2),max(totabs_n2)]
;yr=[0.,3.e-17]
;titl=maj[2]+ 'Absorption Cross-sections'
;
;
;w3=window(window_title=head_titl,dimension=[800,800])
;
;abc3=plot(wave_n2, totabs_n2, name=maj[2],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[2],linestyle=lstyl[0],/current)
;
;w3.Save, plt_loc+"n2_abscross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;_________________________________________________________________________________________________________________________________________________
;_______________________________________________Verification Plots________________________________________________________________________________
;
;xt='Wavelength ($\AA$)'
;yt='Cross-sections (cm!U2!N))'
;
;;O absorption cross-section
;
;
;xr=[min(wave_o),600.]
;yr=[0.,2.e-17]
;titl=maj[0]+ 'Absorption Cross-sections'
;
;
;w4=window(window_title=head_titl,dimension=[800,800])
;
;abc11=plot(wave_o, totabs_o, name='H+F',$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr, yrange=yr,ylog=0,ystyle=1,xlog=0,xstyle=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)
;  
;abc12=scatterplot(wv_o[0:31],tababs_o[0:31],name='Calculated',$
;    sym_thick=sym_thck,sym_color=sym_col[0],symbol=sym[0],/overplot)
;    
;abc13=scatterplot(wv_o[0:31],abs_o[0:31],name='From H+F table',$
;    sym_thick=sym_thck,sym_color=sym_col[1],symbol=sym[1],/overplot)
;
;leg4 = LEGEND(TARGET=[abc11,abc12,abc13], POSITION=leg_loc_left,/normal)
;w4.Save, plt_loc+"o_abscross_ver3.5.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;O2 absorption cross-section
;
;
;xr=[min(wave_o2),600.]
;yr=[0.,3.e-17]
;titl=maj[1]+ 'Absorption Cross-sections'
;
;
;w5=window(window_title=head_titl,dimension=[800,800])
;
;abc21=plot(wave_o2, totabs_o2, name='H+F',$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr, yrange=yr,ylog=0,ystyle=1,xlog=0,xstyle=1,$
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[1],linestyle=lstyl[0],/current)
;  
;abc22=scatterplot(wv_o2[0:29],tababs_o2[0:29],name='Calculated',$
;                  sym_thick=sym_thck,sym_color=sym_col[0],symbol=sym[0],/overplot)
;
;abc23=scatterplot(wv_o2[0:29],abs_o2[0:29],name='From H+F table',$
;                  sym_thick=sym_thck,sym_color=sym_col[1],symbol=sym[1],/overplot)
;
;leg5 = LEGEND(TARGET=[abc21,abc22,abc23], POSITION=leg_loc_left,/normal)
;w5.Save, plt_loc+"o2_abscross_ver3.5.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;
;;N2 absorption cross-section
;
;
;xr=[min(wave_n2),600.]
;yr=[0.,3.e-17]
;titl=maj[2]+ 'Absorption Cross-sections'
;
;
;w6=window(window_title=head_titl,dimension=[800,800])
;
;abc31=plot(wave_n2, totabs_n2, name='H+F',$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[2],linestyle=lstyl[0],/current)
;
;abc32=scatterplot(wv_n2[0:30],tababs_n2[0:30],name='Calculated',$
;    sym_thick=sym_thck,sym_color=sym_col[2],symbol=sym[0],/overplot)
;
;abc33=scatterplot(wv_n2[0:30],abs_n2[0:30],name='From H+F table',$
;    sym_thick=sym_thck,sym_color=sym_col[3],symbol=sym[1],/overplot)
;
;
;leg6 = LEGEND(TARGET=[abc31,abc32,abc33], POSITION=leg_loc_left,/normal)
;
;w6.Save, plt_loc+"n2_abscross_ver3.5.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT





;_________________________________________________________________________________________________________________________________________________

end

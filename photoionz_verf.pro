;Theoritical calculation of photoionisation and verification of the model photoionisation
;June 29, 2020

;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________

;Model data file
floc='/home/srimoyee/Desktop/nrl_files/sav_files/'
restore,filename=floc+'peakm5.sav'
;restore, filename=floc+'peakx9.sav'
; file has:photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2
flux_m5= flux

wv2[25]=975.
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



;Cross-section files and data

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2


;____________________________________________________________________________________________
;
restore, probstate_crss
restore, ionabs_crss
lmax= n_elements(spec_wave)



; SIGABS  photoabsorption cross sections,0- O, 1-O2, 2-N2; cm2
; SIGIONx  photoionization cross sections,0- O,1- O2,2- N2; cm2

sigionx=fltarr(nmaj,lmax)
sigabs=fltarr(nmaj,lmax)

sigabs(0,*)=sigab_o
sigabs(1,*)=sigab_o2
sigabs(2,*)=sigab_n2

sigionx(0,*)=sigi_o
sigionx(1,*)=sigi_o2
sigionx(2,*)=sigi_n2
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
; Absorption cross-sections
;sigab_o=[0.0023,0.0170,0.1125,0.1050,0.3247,1.319,3.7832,6.0239,7.7205,10.7175,13.1253,8.5159,4.7889,$
;  3.0031,4.1048,3.7947,0.0,0.0,0.0,0.0,0.0,0.0]*1e-18  ;in cm^2
;sigab_o2=[0.0045,0.034,0.2251,0.2101,0.6460,2.6319,7.6283,13.2125,16.8233,20.3066,27.0314,23.5669,24.9102,$
;  10.498,10.9075,13.3122,13.395,14.4042,32.5038,18.7145,1.6320,1.1500]*1e-18  ;in cm^2
;sigab_n2=[0.0025,0.0201,0.1409,1.1370,0.3459,1.5273,5.0859,9.9375,11.7383,19.6514,23.0931,23.0346,54.5252,$
;  2.1434,13.1062,71.6931,2.1775,14.4390,115.257,2.5465,0.0,0.0]*1e-18  ;in cm^2
;  
;  sigabs[0,12:30]=sigab_o[3:-1]
;  sigabs[1,12:30]=sigab_o2[3:-1]
;  sigabs[2,12:30]=sigab_n2[3:-1]
  

;;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

;PROB    branching ratios for each state, species, and wavelength bin:
nst =10

prob=fltarr(nst,nmaj,lmax) ;states species and energies, 0-O, 1-O2, 2- N2


;0-O
n_cols=(size(O_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,0,l)=O_prob_state(k,l)
endfor

;1-O2
n_cols=(size(O2_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,1,l)=O2_prob_state(k,l)
endfor


;2-N2
n_cols=(size(N2_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
  for k=1,n_cols-1 do $
    prob(k-1,2,l)=N2_prob_state(k,l)
endfor
;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

;Calculating theoritical flux

flux_th = fltarr(jmax,lmax)

tau_th=fltarr(jmax,lmax)

pi=fltarr(jmax,lmax,nst,nmaj)
tau_th=fltarr(jmax,lmax)
pi_o=fltarr(jmax,lmax)
pi_o2=fltarr(jmax,lmax)
pi_n2=fltarr(jmax,lmax)

for i=0 , lmax-1 do begin
  
tau_th(*,i)= (H(0,*)*sigabs(0,i)*zo + H(1,*)*sigabs(1,i)*zo2 + H(2,*)*sigabs(2,i)*zn2)*1.e5; optical depth
flux_th(*,i) = ssflux(i)*exp(-1.*tau_th[*,i])  
 
 for j=0, jmax-1 do begin
  
    for l=0,nst-1 do begin
          pi(j,i,l,0)=flux_th(j,i)*prob(l,0,i)*sigionx(0,i)*zo(j)
          pi(j,i,l,1)=flux_th(j,i)*prob(l,1,i)*sigionx(1,i)*zo2(j)
          pi(j,i,l,2)=flux_th(j,i)*prob(l,2,i)*sigionx(2,i)*zn2(j)
          
          
     endfor
     
    pi_o(j,i)= total(pi(j,i,*,0))
    pi_o2(j,i)=total(pi(j,i,*,1))
    pi_n2(j,i)=total(pi(j,i,*,2))
    
  endfor
 

endfor


;***************************************************************************************************************************
;  ; Input cross-section files
;
; probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_SQ05.sav'
;;  file has : O_prob_state,O2_prob_state,N2_prob_state
; ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec_SQ05.sav'
;;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
;
;restore, probstate_crss
;restore, ionabs_crss
;
;sigionx[0,12:30]=sigi_o[3:-1]
;sigionx[1,12:30]=sigi_o2[3:-1]
;sigionx[2,12:30]=sigi_n2[3:-1]
;
;
;;Calculating theoritical flux
;
;wl=fix(12)
;
;for i=3 , 21 do begin
;  
;  for j=0, jmax-1 do begin
;   for l=1,3 do begin
;      pi_o(j,wl)=pi_o(j,wl)+(flux(j,wl)*o_prob_state(l,i)*sigionx(0,wl)*zo(j))
;     
;
;    endfor
;  endfor
;  wl++
;  print, wl
;endfor
;
;wl=fix(12)
;for i=3 , 21 do begin
;  for j=0, jmax-1 do begin
;    for l=1,2 do begin
;      pi_o2(j,wl)=pi_o2(j,wl)+(flux(j,wl)*o2_prob_state(l,i)*sigionx(1,wl)*zo2(j))
;      pi_n2(j,wl)=pi_n2(j,wl)+(flux(j,wl)*n2_prob_state(l,i)*sigionx(2,wl)*zn2(j))
;
;
;    endfor
;  endfor
;  wl++
;endfor

 ;___________________________________________________________________________________________________________________________________________________________________________________
 ;___________________________________________________________________________________________________________________________________________________________________________________

; Checking flux with model

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


; flux

; xr=[min(wv1),max(wv1)]
; yr=[0.,5.e11]
;xt='Wavelength ($\AA$)'
;yt='Flux (Photons cm!U-2!N s!U-1!N))'
;restore, filename=floc+"photoionz_verf_flux.sav"
;
;wv2[25]=975.
;for i=0, jmax-1 do begin
;  
;     titl='Flux (M5)'+strtrim(string(fix(zz[i])),1)+"km"
;
;
;      ;yr=[1.e4,1.e12]
;;      yr=[min(flux[i,*]),max(flux[i,*])]
;      
;  
;     w1=window(window_title=titl,dimension=[800,800])
;
;     s1=plot(wv1, flux_th[i,*], name=leg_nam[0],$
;             xtitle=xt,ytitle=yt,title=titl,$
;             yrange=yr, ylog=1,ystyle=1,xlog=1,xstyle=1,$ ;xrange=xr,
;             xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;             thick=thck,histogram=1,color=colors[0],linestyle=lstyl[0],/current)
;
;     s2=plot(wv1, flux_th_HF[i,*], name=leg_nam[1],$
;             thick=thck,histogram=1,color=colors[1],linestyle=lstyl[1],/overplot)
;
; 
;     leg1 = LEGEND(TARGET=[s1,s2], POSITION=leg_loc_left,/normal)
;
;;     w1.Save, plt_loc+"m5_flux_"+strtrim(string(fix(zz[i])),1)+"km.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;  
;endfor

;flux_th_HF=flux_th
;save, wv1, wv2,flux_th_HF, zz, filename=floc+"photoionz_verf_flux.sav"
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

;Cross-section plots

leg_nam=['H+F','SQ05']


xr=[0.,1050.]
yr=[min(sigab_o/1.e-18),60.]
xt='Wavelength ($\AA$)'
yt='Cross-section (cm!U2!N )'
titl='Cross-section'

;#O

;w2=window(window_title=titl,dimension=[800,800])
;
;ab1=scatterplot((wv1+wv2)/2., sigabs[0,*]/1.e-18, name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
;  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/current)
;  
; ab2=scatterplot(spec_wave, sigab_o/1.e-18, name=leg_nam[1],$
;          sym_color=colors[1],symbol=sym[0],sym_thick=sym_thck,/overplot)
;    
;    
;leg2 = LEGEND(TARGET=[ab1,ab2], POSITION=leg_loc_left,/normal)
;
;w2.Save, plt_loc+"Oabs_cross_aw.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;#O2
;
;yr=[min(sigab_o2/1.e-18),max(sigab_o2/1.e-18)]
;w3=window(window_title=titl,dimension=[800,800])
;
;ab3=scatterplot((wv1+wv2)/2., sigabs[1,*]/1.e-18, name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=maj[1]+titl,$
;  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/current)
;  
; ab4=scatterplot(spec_wave, sigab_o2/1.e-18, name=leg_nam[1],$
;    sym_color=colors[1],symbol=sym[0],sym_thick=sym_thck,/overplot)
;
;
;  leg3= LEGEND(TARGET=[ab3,ab4], POSITION=leg_loc_left,/normal)
;
;w3.Save, plt_loc+"O2abs_cross_aw.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;;#N2
;
;yr=[min(sigab_n2/1.e-18),max(sigab_n2/1.e-18)]
;w4=window(window_title=titl,dimension=[800,800])
;
;ab5=scatterplot((wv1+wv2)/2., sigabs[2,*]/1.e-18, name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=maj[2]+titl,$
;  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/current)
;
;
;ab6=scatterplot(spec_wave, sigab_o/1.e-18, name=leg_nam[1],$
;    sym_color=colors[1],symbol=sym[0],sym_thick=sym_thck,/overplot)
;
;
;  leg4= LEGEND(TARGET=[ab5,ab6], POSITION=leg_loc_left,/normal)
;  
;w4.Save, plt_loc+"N2abs_cross_aw.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT



;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________
 
; Plotting the original base cross-section files


                             
             
        
          
            
           
      


restore, '/home/srimoyee/Desktop/nrl_files/sav_files/ion&abs_cross_files.sav' ;from Henke & Fennely
; file has:
;wave_o,wave_o2,wave_n2
;totionz_o,totionz_o2,totionz_n2
;totabs_o,totabs_o2,totabs_n2

;xr=[17.,660.]

xr=[640.,990.]
;yr=[min(totabs_o/1.e-18),max(totabs_o/1.e-18)]
yr=[0.,50.]
;arrow_yr=[max(totabs_o/1.e-18),min(totabs_o/1.e-18)]
arrow_yr=[yr[1],yr[0]]
;#O

w5=window(window_title=titl,dimension=[800,800])

ab_o=plot(wave_o, totabs_o/1.e-18, name=maj[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[1],/current)

ab_o1=scatterplot(spec_wave, sigab_o/1.e-18, name=leg_nam[1],$
    sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)
    
ab_o2=scatterplot((wv1+wv2)/2., sigabs(0,*)/1.e-18, name=leg_nam[1],$
    sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)

;o_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

o_1 = ARROW([650,650 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_2 = ARROW([798,798], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_3 = ARROW([913,913], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o_4 = ARROW([975,975], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

;w5.Save, plt_loc+"Oabs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;#O2
;
yr=[0.,50.]
arrow_yr=[yr[1],yr[0]]


w6=window(window_title=titl,dimension=[800,800])

ab_o2=plot(wave_o2, totabs_o2/1.e-18, name=maj[1],$
  xtitle=xt,ytitle=yt,title=maj[1]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[2],/current)

ab_o21=scatterplot(spec_wave, sigab_o2/1.e-18, name=leg_nam[1],$
    sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)
    
ab_o22=scatterplot((wv1+wv2)/2., sigabs(1,*)/1.e-18, name=leg_nam[1],$
    sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)
;
;o2_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;o2_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

o2_1 = ARROW([650,650 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_2 = ARROW([798,798], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_3 = ARROW([913,913], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
o2_4 = ARROW([975,975], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

;w6.Save, plt_loc+"O2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;#N2

yr=[0.,120.]
;yr=[min(totabs_n2/1.e-18),max(totabs_n2/1.e-18)]
arrow_yr=[yr[1],yr[0]]

w7=window(window_title=titl,dimension=[800,800])

ab_n2=plot(wave_n2, totabs_n2/1.e-18, name=maj[2],$
  xtitle=xt,ytitle=yt,title=maj[2]+titl,$
  xrange=xr,yrange=yr, ylog=0,ystyle=1,xlog=0,xstyle=1,$ ;xrange=xr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[4],/current)
  
ab_n21=scatterplot(spec_wave, sigab_n2/1.e-18, name=leg_nam[1],$
    sym_color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)
    
ab_n22=scatterplot((wv1+wv2)/2., sigabs(2,*)/1.e-18, name=leg_nam[1],$
    sym_color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)

;n2_1 = ARROW([18,18 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_2 = ARROW([32,32], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_3 = ARROW([70,70], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_4 = ARROW([155,155], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_5 = ARROW([224,224], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_6 = ARROW([290,290 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_7 = ARROW([320,320], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_8 = ARROW([540,540], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
;n2_8 = ARROW([650,650], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

n2_1 = ARROW([650,650 ], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_2 = ARROW([798,798], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_3 = ARROW([913,913], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)
n2_4 = ARROW([975,975], arrow_yr, COLOR=colors[5],line_thick=ln_thck,/DATA, /CURRENT)

;w7.Save, plt_loc+"N2abs_basecross.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
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


;  NRLEUV test model............................

file_euv= '/home/srimoyee/Desktop/bailey_bins/nrleuv_bastille_flare.dat'
openr,lun_euv, file_euv, /get_lun
skip_lun, lun_euv, 14,/lines
n_euv=2398  ;wavelength bins
euvwvln=dblarr(n_euv)  ;in angstrom
euvflux=dblarr(n_euv)  ;in  photons/cm2/s

format='$(f8.2,f8.1,a11,i1,f10.6)'

for i=0,n_euv-1 do begin
  readf,lun_euv,a,b,c
  euvwvln(i)=a
  euvflux(i)=c ; to be in units of photons/cm2/s
endfor
free_lun,lun_euv
euvener=12397./euvwvln

ii_hint1=where((hintwvln ge 645.) and (hintwvln lt 1027.),n_hint1)
ii_hint2=where((hintwvln ge 798.) and (hintwvln lt 913.),n_hint2)
ii_hint3=where((hintwvln ge 913.) and (hintwvln le 975.),n_hint3)

y_lab1=make_array(n_hint1,value=1.e10)
y_lab2=make_array(n_hint2,value=1.e10)
y_lab3=make_array(n_hint3,value=1.e10)

; flux

xr=[640.,990.]
yr=[0.,5.e11]

;yr=[1.e4,1.e12]
;      yr=[min(flux[i,*]),max(flux[i,*])]
xt='Wavelength ($\AA$)'
yt='Flux (Photons cm!U-2!N s!U-1!N))'
arrow_yr=[1.e10,1.]

titl='Hinteregger Solar Flux'

 
w8=window(window_title=titl,dimension=[800,800])

h1=plot(hintwvln, hintflux, $
             xtitle=xt,ytitle=yt,title=titl,$
             xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$ 
             xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
             thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)

ss1=scatterplot((wv1+wv2)/2., ssflux,$
         sym_color=colors[4],sym_thick=thck,/overplot)

;h_1 = ARROW([650,650], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
;h_2 = ARROW([798,798], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)
;h_3 = ARROW([913,913], arrow_yr, COLOR=colors[3],line_thick=ln_thck, /DATA, /CURRENT)
;h_4 = ARROW([975,975], arrow_yr, COLOR=colors[3],line_thick=ln_thck,/DATA, /CURRENT)



;t1 = TEXT(hintwvln[ii_hint1],y_lab1,hintlable[ii_hint1], /data,font_size=10,orientation=90,font_color=colors[3])
;
;w9=window(window_title=titl,dimension=[800,800])
;xr=[798.,913.]
;h2=plot(hintwvln, hintflux, $
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)
;
;
;;ss2=plot(wv1, ssflux,$
;;    color=colors[1],histogram=1,thick=thck,/overplot)
;
;t2 = TEXT(hintwvln[ii_hint2],y_lab2,hintlable[ii_hint2], /data,font_size=10,orientation=90,font_color=colors[3])
;
;
;w10=window(window_title=titl,dimension=[800,800])
;xr=[913.,975.]
;h3=plot(hintwvln, hintflux, $
;  xtitle=xt,ytitle=yt,title=titl,$
;  xrange=xr,yrange=yr, ylog=1,ystyle=1,xlog=0,xstyle=1,$
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,histogram=0,color=colors[0],linestyle=lstyl[0],/current)
;
;;ss3=plot(wv1, ssflux,$
;;    color=colors[1],histogram=1,thick=thck,/overplot)
;
;t3 = TEXT(hintwvln[ii_hint3],y_lab3,hintlable[ii_hint3], /data,font_size=10,orientation=90,font_color=colors[3])
;
;;w8.Save, plt_loc+"hint_flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


;___________________________________________________________________________________________________________________________________________________________________________________
;___________________________________________________________________________________________________________________________________________________________________________________

  ;O Photoionisation plot
  
;  pi_o_hf=pi_o
;  pi_o2_hf=pi_o2
;  pi_n2_hf=pi_n2
;save, pi_o_hf,pi_o2_hf,pi_n2_hf,filename="m5_photoi_hf.sav"
;restore,"m5_photoi_hf.sav"
;  
;  yr=[min(zz),max(zz)]
;  xr=[1.e-9,10000.]
;  xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
;  yt='Altitude (km)'
;  titl='Photoionization'
; 
;for i=0,29 do begin
;  
;  tit="("+strtrim(string(wv1[i]),1)+"-"+strtrim(string(wv2[i]),1)+"$\AA$)"
;  
;;  xr=[min(photoi_wv[*,0,i]),max(photoi_wv[*,0,i])]
;  
;  w1 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=[800,800])
;
;  O_pi1=plot(pi_o_hf[*,i],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=maj[0]+titl+tit,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],layout=[1,3,1],/current)
;
;  O_pi2=plot(pi_o[*,i],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  
;  leg1 = LEGEND(TARGET=[O_pi1,O_pi2], POSITION=leg_loc,/normal)
;
;
;
;;_________________________________________________________________________
;
;  ;O2 Photoionisation plot
;
;;  xr=[min(photoi_wv[*,1,i]),max(photoi_wv[*,1,i])]          
; 
;  
;
;  O2_pi1=plot(pi_o2_hf[*,i],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=maj[1]+titl+tit,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],layout=[1,3,2],/current)
;
;  O2_pi2=plot(pi_o2[*,i],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;  
;
;;  leg2 = LEGEND(TARGET=[O2_pi1,O2_pi2], POSITION=leg_loc2,/normal)
;
;
;
;;_________________________________________________________________________
;
;;N2 Photoionisation plot
;
;;xr=[min(photoi_wv[*,2,i]),max(photoi_wv[*,2,i])]         
;
;N2_pi1=plot(pi_n2_hf[*,i],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=maj[2]+titl+tit,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],layout=[1,3,3],/current)
;
;N2_pi2=plot(pi_n2[*,i],zz,name=leg_nam[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;
;
;;leg3 = LEGEND(TARGET=[N2_pi1,N2_pi2], POSITION=leg_loc3,/normal)
;
;;w1.Save, plt_loc+"photoi_" +strtrim(string(i),1)+".jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;endfor

end


;____________________________________________________________________________________________________________________________________________________________________

;_______________________________________________Main photoelectron calling function__________________________________________________________________________________

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software

;________________________________________________Model constants and Inputs_____________________________________________________________________________________________________

nbins=176

nmaj=3
;lat=5.    ;x9 flare
;lon=0.0


lat=20.  ;m5 flare
lon=105.0



;ALTITUDE BINS
z0=40.
zz=findgen(161)*1.+z0
;zz=findgen(17)*10.+z0
;zz=findgen(400)*0.5+z0
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]

;idate=2002089 ;SQ'05
idate =2016205  ;m5 flare
;idate=2017249   ;x9 flare


floc='/home/srimoyee/Desktop/nrl_files/sav_files/'
;____________________________________________________________________________________________________________________________________________________________________
;________________________________________________F10.7 and ap values _____________________________________________________________________________________________________


restore,floc+'f107datafile.sav'
;this file was created from f107_values.pro from
;data available from
;https://lasp.colorado.edu/lisird/data/penticton_radio_flux/
; file has: month, day, year, hour, minute, second,yyyyddd,f107d,f107_81
f107_ii=where(yyyyddd eq idate,count)
if count gt 0 then begin
  f107= f107d[f107_ii[0]]
  f107a=f107_81[f107_ii[0]]
  print, f107,f107a
endif else begin
  f107=70.
  f107a=70.
  print,'F107 values not found for the given date, taking default value of (f107, f107-81day) ',f107, f107a
endelse

;OR USE YOUR OWN VALUES
;f107=70.
;f107a=70.
;____________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________

ap=15


;____________________________________________________________________________________________________________________________________________________________________

;_________________________________________________INPUT SOLAR FLUX___________________________________________________________________________________________________

;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
restore, '/home/srimoyee/Desktop/nrl_files/sav_files/SQ_newspectra.sav'

;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5

wv1=wave_gcm1
wv2=wave_gcm2
;ssflux=meanpeakx9

ssflux=peakm5
;if (hi eq 0) then ssflux[hi]= 2.e9
;if ((hi ge 0) and (hi le 5)) then ssflux[hi]= 2.e9
mask=fltarr(n_elements(wv1))
mask[hi]=1.    ;hi
ssflux=ssflux*mask

lmax=n_elements(wv1)

;____________________________________________________________________________________________________________________________________________________________________

;_________________________________________________PHOTOIONISATION CROSS-SECTION SETUP _______________________________________________________________________________

;CALL THE NEW CROSS SECTIONS FILE GENERATION HERE
;
;
;
;  Input cross-section files
;
probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_SQ05.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec_SQ05.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2




;____________________________________________________________________________________________

;probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states.sav'
;;  file has : O_prob_state,O2_prob_state,N2_prob_state
;ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
;____________________________________________________________________________________________


localpe_setup_pxsect, probstate_crss, ionabs_crss,lmax, sigionx,sigabs,prob,auger_wvln,auger_energy,tpot, nst
; INPUT: probstate_crss, ionabs_crss,lmax
;OUTPUT: sigionx,sigabs,prob,auger_wvln,auger_energy,tpot, nst


;____________________________________________________________________________________________________________________________________________________________________

;_________________________________________________ELECTRON CROSS-SECTION SETUP ___________________________________________________________________________________________________


if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi

endif


dum=min(abs(ener-20.),ia20)
dum=min(abs(ener-100.),ia100)
dum=min(abs(ener-1000.),ia1000)





;____________________________________________________________________________________________________________________________________________________________________

;____________________________________________________________________________________________________________________________________________________


;output variables - for calculating parameters in each bin

;pepi=fltarr(nmaj,lmax)
;tau1h=fltarr(lmax)
;ii_check=intarr(12)




eionz_wv_loc = fltarr(nmaj,jmax,lmax)
eionz_wv_transp = fltarr(nmaj,jmax,lmax)
photi_wv = fltarr(nmaj,jmax,lmax)

tic
;_______________________________________________________________________________________________________________________________________________________________________


;____________________________________________________________________________________________________________________________________________________________________

;_________________________________________________NEUTRAL ATMOSPHERE AND ELECTRON DENSITY AND TEMPERATURE___________________________________________________________________________________________________

;f you need to loop over time in a day
;    n_hours=14
;    glow20=fltarr(n_hours)
;    glow100=fltarr(n_hours)
;    glow1000=fltarr(n_hours)
;    ace20=fltarr(n_hours)
;    ace100=fltarr(n_hours)
;    ace1000=fltarr(n_hours)
;    asza=fltarr(n_hours)




;for i=0,n_hours-1 do begin
;     i=6  ;i th hour
hour= 5  ; m5 flare UT time
;     hour =12  ; x9 flare UT time
;    hour=float(i)+6.
utsec=3600.*hour
;________________________________________________________________________________


;doy=ace_date(idate,/to_doy,/from_yyyyddd)
doy=idate-(long(idate/1000)*1000)
utdoy=doy+utsec/(3600.*24.) ; ut in days since Jan 1 of year (input to zenith.pro)
sza=zenith(utdoy,lat,lon);105*!dtor;
localhour=(((utsec/240. + lon)/15. + 24. ) mod 24.)
;___________________________________________________________________________________________________________________________________________________


localpe_neutralatm,idate, utsec, zz,lat,lon,localhour, f107a,f107,ap,sza,zmaj,zcol
;INPUT: idate, utsec, zz,lat,lon,localhour, f107a,f107,ap,sza
;OUTPUT: zmaj,zcol

localpe_approxeden,zz,f107,f107a,idate,lat,lon,utsec,eden,etemp
;INPUT: zz,f107,f107a,idate,lat,lon,utsec
;OUTPUT: eden,etemp


;______________________________________________________________________________________________________________________________________________________________________

;_____________________________________________________________PHOTOIONISATION__________________________________________________________________________________________

localpe_photoionz,nbins,zcol,sigabs,sigionx,prob,tpot,auger_energy,wv1,wv2 ,tau, flux, photoi, photoi_wv, primary
;INPUT: nbins,zcol,sigabs,sigionx,prob,tpot,auger_energy,wv1,wv2
;OUTPUT: tau, flux, photoi, photoi_wv, primary





;__________________________________________________________________________________________________________________________






 
          
for hi=0,0 do begin  ;solar flux loop

  localpe_ssflux,hi

  

  localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo

  ace_etransport

  ;     Printing out parameters for optical depth tau=1
  ii=where(reform(tau(*,hi)) ge 1.,/null)
  min_tau=min(tau(ii,hi),min_i)

  ;min_dif=min(abs(reform(tau(*,hi))-1.),min_i)
  ;      print, ii[min_i]
  ;      ii_check[hi]=ii[min_i]
  ;print, min_i
  print,'Optical Depth: ' ,min_tau
  print,'Height of unit optical depth:',hi,'.',zz[ii[min_i]]
  print,'O:',eiionz[0,ii[min_i]]/photoi[0,ii[min_i]],eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]
  print,'O2:',eiionz[1,ii[min_i]]/photoi[1,ii[min_i]],eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]
  print,'N2:',eiionz[2,ii[min_i]]/photoi[2,ii[min_i]],eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]


  ;      pepi[0,hi]=eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]
  ;      pepi[1,hi]=eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]
  ;      pepi[2,hi]=eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]
  ;      tau1h[hi]=zz[ii[min_i]]


  ;      min_difO=min(abs(reform(tau_wv(0,*,i))-1.),Oi)
  ;      min_difO2=min(abs(reform(tau_wv(1,*,i))-1.),O2i)
  ;      min_difN2=min(abs(reform(tau_wv(2,*,i))-1.),N2i)
  ;
  ;      print,'O:',eiionz[0,Oi]/photoi[0,Oi]
  ;      print,'O2:',eiionz[1,O2i]/photoi[1,O2i]
  ;      print,'N2:',eiionz[2,N2i]/photoi[2,N2i]
  ;      print,""
  ;         print,hi
  ;
  eionz_wv_loc[0,*,hi] =eiionz[0,*]
  eionz_wv_loc[1,*,hi] =eiionz[1,*]
  eionz_wv_loc[2,*,hi]=eiionz[2,*]

  eionz_wv_transp[0,*,hi] =eiionz_transp[0,*]
  eionz_wv_transp[1,*,hi] =eiionz_transp[1,*]
  eionz_wv_transp[2,*,hi]=eiionz_transp[2,*]

  photi_wv[0,*,hi] =photoi[0,*]
  photi_wv[1,*,hi] =photoi[1,*]
  photi_wv[2,*,hi]=photoi[2,*]





  ;
endfor  ;solar flux loop end
;
;pepi[where(~finite(pepi),/null)]=0.



toc
;____________________________________________________________________________________________________________________________________________________________________

;; Energy conservation tests for photoionisation and photoelectrons
;
;
;;the height integrated photoionisation energy should add up to the solar flux,
;;atleast for the lower wavlength bins
;
;photoi_ht=dblarr(lmax)  ;height integrated photoionisation energy
;
;;photoi_wv:Array[alt, species, wavelength]
;photoi_stot=total(photoi_wv,2)  ;total photoization in all species
;
;  for l=0, lmax-1 do begin
;    photoi_ht(l)=total(reform(photoi_stot(*,l)*dz*1.e5)) ;height integrated ionisation rate
;  endfor
;h=6.62607004e-34
;c=3.e8
;E=h*c/(0.5*1.e-10*(wv1+wv2))   ;J cm-2 s-1
;photoi_ht=photoi_ht*E;*6.242e18
;ssflux_ht=ssflux*E
;
;print,total(photoi_ht),total(ssflux_ht)






;;_____________________________________________________________________________________________________________________________
;
;energy conservation test  for transport code

eiionz_ht= fltarr(nbins)
eiionz_stot=total(eiionz_transp_wv,2)  ;species total   ; eiionz_wv:Array[alt, species, wavelength]

for i=0, nbins-1 do $
  eiionz_ht[i]=total(reform(eiionz_stot(*,i)*dz*1.e5)) ;height integrated electron ionisation rate

eiionz_ht=eiionz_ht;*1.60218e-19*ener  ;J cm-2 s-1
;
;  eiionz_tot=fltarr(lmax)  ;putting the electron ionisation energy in solar flux bins
;  prim_tot=fltarr(jmax,lmax)
;
;  for l=0,lmax-1 do begin
;    ii=where(((12397./ener ge wv1[l])and (12397./ener lt wv2[l])),cnt,/null)
;    if cnt gt 0 then begin $
;         eiionz_tot[l]=total(eiionz_ht(ii))
;         prim_tot[*,l]=total(primary(*,ii))
;    endif
;  endfor

;  for l=0, lmax-1 do begin
;    ii=where(((ener ge wv1[l])and (ener lt wv2[l])),cnt,/null)
;    prim_tot[*,l]=total(primary(*,ii))
;
;
;  endfor
tflux_ht=fltarr(nbins)
tflux_ht= total(tflux,1)
;  for i=0, nbins-1 do $
;      tflux_ht(i)=total(reform(tflux(*,i)*del(i)))
;

pespec_ht=fltarr(nbins)
pespec_ht= total(pespec,1)

prim_ht=fltarr(nbins)
for i=0, nbins-1 do $
  prim_ht(i)=total(reform(primary(*,i)*dz*1.e5))


cgdisplay,2000,2000, /free

cgplot,ener,eiionz_ht,color='red',thick=1,/xs,/ys,/yl,/xl,yr=[1.e2,1.e10],psym=-46,$
  ytitle='photons cm-2 s-1 ',xtitle='Energy',xr=[0.01,1.e5]
cgplot, ener,tflux_ht*del,color='dodger blue',/overplot,thick=3
cgplot, ener,prim_ht,color='green',/overplot,thick=3
cglegend, titles=['Height Integrated Photoelectron ','Tflux','Primary'],$
  psyms=[15], symcolors=['red','dodger blue','green'],$
  alignment=1,/box,Location=[0.93, 0.9]


;____________________________________________________________________________________________________________________________________________________________________

; output- for running each solar flux bin at a time

;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_m5time0.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_m5time316.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_x9time0.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_x9time720.sav'




;save, photoi, eiionz_transp,pepi, zz, ssflux, wv1,wv2,zcol,eden,etemp, filename=floc+'SQ05_zcol_neut3.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_org.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_iri.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_neut1.sav'
;save,pepi,tau1h,wv1,wv2,sza,filename=floc+'SQ05_sza6_4.sav'

;save,eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'SQpeakm5_wv_szalo_test.sav'
;save,eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'peakm5_wv_szalo_2states.sav'

;output- for running the entire solar flux in one go
;
;pepi=fltarr(nmaj,lmax)
;pepi=eiionz_transp/photoi
;pepi[where(~finite(pepi),/null)]=0.

;save, photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2, zcol,filename=floc+'peakm5.sav'
;save, photoi, photoi_wv,eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'peakx9.sav'
;
;save,eionz_wv, zz, ssflux, flux,wv1,wv2,filename=floc+'peakm5_eiionz_v2.sav'
;save, eionz_wv, zz, ssflux, flux,wv1,wv2, filename=floc+'peakx9_eiionz_v2.sav'

;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'m5_time0.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'m5_time316.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'x9_time0.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'x9_pt5-4.sav'




;___________________________________________________________________________________________________
;z=zz*1.e5
;am=[16., 32., 28.] ; amu for O, O2 anf N2
;kb= 1.38064852e-16   ; Boltzmann's constant ergs per k
;
;re=6378.e5 ; radius of the Earth in cm
;g=980.  ;acceleration due to gravity in cm^2/sec
;gr=g*(re/(re+z))^2 ; account for small change in g with altitude
;
;amave=(.2*32 + .79*28 + .01*40.) ; average mass of atmosphere, used for lower altitudes
;gperamu=1.66054*1.e-24  ;coverting amu to grams
;
;
;;    Assumption: at low altitudes, every constituent follows
;;                a single scale height, above that, they follow
;;                their own, assume change occurs at 100km
;
;
;ilow=z le 100.*1e5
;ihi=1-ilow
;
;H=fltarr(3,n_elements(zz))
;for i=0,2 do begin
;   amass=amave*ilow + am[i]*ihi     ;mass for every height
;
;   H[i,*]=(kb*zt/(amass*gperamu*gr))/1.e5
;
;
;endfor
;___________________________________________________________________________________________________

;Theoritical calculation of scale height H

;z_0=0.
;H=7.
;n0=1.90120e19
;n0_1=1.90120e19;2.7e19
;chi=0
;abs_n2=[2.48E-05,1.47E-04,4.57E-04,1.02E-03,1.90E-03,4.00E-03,8.55E-03,$
;  1.55E-02,3.12E-02,6.32E-02,0.14,0.3]*1e-18
;
;;zmax=[43.,55.5,63.5,69.1,73.3,78.7,84.,88.1,93,98,103.5,109.]  ;km
;;
;;;nz=n0*exp((zmax-z_0)/H)
;;
;;
;;n0_2=[1.63198e+19,3.01390e+19,5.72532e+19,1.13414e+20,1.85562e+20,2.44977e+20,6.30649e+20,1.83762e+21,1.57760e+21,$
;;    9.80887e+20,1.19087e+21,5.36025e+20]
;zmax_z1=fltarr(12)
;;zmax_z2=fltarr(12)
;
;
;;
;zmax_z=z_0+H*alog(H*1.e5*n0*abs_n2/cos(chi))
;;
;for i=0,11 do begin
;  zmax_z1[i]=z_0+H[2,ii_check[i]]*alog(H[2,ii_check[i]]*1.e5*n0_1*abs_n2[i])
;  print, H[2,ii_check[i]]
;;  zmax_z2[i]=z_0+H[2,ii_check[i]]*alog(H[2,ii_check[i]]*1.e5*n0_2[i]*abs_n2[i])
;
;endfor
;;zmax=[45.,58,65.5,71,75,79,83.5,87,90.5,94.5,99.,103.]
;;n_0=fltarr(12)
;;for i=0,11 do $
;;   n_0[i]=(1/(H[2,ii_check[i]]*1.e5*abs_n2[i]))*exp((zmax[i]-z_0)/H[2,ii_check[i]])
;
;
;




;____________________________________________________________________________________________________





;____________________________________________________________________________________________________________________________________________________________________

end
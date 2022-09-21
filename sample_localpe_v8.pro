; This version tests for different soalr zenith angles

;____________________________________________________________________________________________________________________________________________________________________

;_______________________________________________Main photoelectron calling function__________________________________________________________________________________

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software

;________________________________________________Model constants and Inputs_____________________________________________________________________________________________________
nbins=178
nmaj=3
 
;
;lat_arr=0.+findgen(17)*5.; m5 flare
;lon_arr=20.+findgen(16)*10.

lat_arr=0.+findgen(16)*5.;x9 flare
lon_arr=0.+findgen(18)*5.

;altitude bins
z0=40.
zz=findgen(181)*1.+z0

jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]

;idate=2002089 ;SQ'05
;idate =2016205  ;m5 flare
idate=2017249   ;x9 flare


floc='/home/srimoyee/Desktop/nrl_files/sav_files/'

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

ap=15
lmax= 22
;_______________________________________________________________________________________________________________________________________________________________________

;output variables - for calculating parameters in each bin

;pepi=fltarr(nmaj,46)
;tau1h=fltarr(46)
;ii_check=intarr(12)
;eionz_wv= fltarr(jmax,nmaj,46)

ij = 0
n_sza= n_elements(lat_arr)*n_elements(lon_arr)

latlon_arr = strarr(n_sza)
photoi_wv_sza=fltarr(jmax,nmaj,lmax,n_sza)
photoi_sza=fltarr(nmaj,jmax,n_sza)
eiionz_sza=fltarr(nmaj,jmax,n_sza)
eiionz_transp_sza=fltarr(nmaj,jmax,n_sza)
sza_arr=fltarr(n_sza)

tic
;_______________________________________________________________________________________________________________________________________________________________________


localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi

endif


dum=min(abs(ener-20.),ia20)
dum=min(abs(ener-100.),ia100)
dum=min(abs(ener-1000.),ia1000)
;____________________________________________________________If you need to loop over time in a day________________________________________________________________________________
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
;hour= 5  ; m5 flare UT time
hour =12  ; x9 flare UT time
;    hour=float(i)+6.
utsec=3600.*hour

localpe_ssflux

 for j=0 ,n_elements(lat_arr)-1 do begin

 for i=0,n_elements(lon_arr)-1 do begin    
  lat = lat_arr[j]
  lon= lon_arr[i]
  
  

  localpe_neutralatm
  ;    asza[i]=sza

  localpe_approxeden;,zz,f107,f107a,eden,etemp,idate,lat,lon,utsec

  ;______________________________________________________________________  Solar flux__________________________________________________________________________________________

  ;______________________________________________________________________________________________________________________________________________________________________-

  ;

  ;for hi=0, 11 do begin  ;solar flux loop

  ;       localpe_ssflux;,hi

  localpe_photoionz

  localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo

  ace_etransport

  ;     Printing out parameters for optical depth tau=1
  ;      ii=where(reform(tau(*,hi)) ge 1.,/null)
  ;      min_tau=min(tau(ii,hi),min_i)
  ;
  ;      ;min_dif=min(abs(reform(tau(*,hi))-1.),min_i)
  ;      print, ii[min_i]
  ;;      ii_check[hi]=ii[min_i]
  ;      ;print, min_i
  ;      print,'Optical Depth: ' ,min_tau
  ;      print,'Height of unit optical depth:',hi,'.',zz[ii[min_i]]
  ;      print,'O:',eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]
  ;      print,'O2:',eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]
  ;      print,'N2:',eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]
  ;
  ;
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
  ;      eionz_wv[*,0,hi] =eiionz_transp[0,*]
  ;      eionz_wv[*,1,hi] =eiionz_transp[1,*]
  ;      eionz_wv[*,2,hi]=eiionz_transp[2,*]
  ;
  ; endfor  ;solar flux loop end
  ;
  ;pepi[where(~finite(pepi),/null)]=0.
  ;_________________________________________________________________________
  ;longitude loop


  photoi_wv_sza[*,*,*,ij]= photoi_wv
  photoi_sza[*,*,ij]= photoi[*,*]
  eiionz_sza[*,*,ij]=eiionz
  eiionz_transp_sza[*,*,ij]= eiionz_transp
  sza_arr[ij]= sza
  latlon_arr[ij]=strtrim(string(fix(lat)))+'-'+strtrim(string(fix(lon)))
  
  
  
  print,ij,'Done'
  ij++

  endfor   ;longitude loop

endfor    ;latitude loop

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
;
;;_____________________________________________________________________________________________________________________________
;;
;  ;energy conservation test  for transport code
;
;  eiionz_ht= fltarr(nbins)
;  eiionz_stot=total(eiionz_transp_wv,2)  ;species total   ; eiionz_wv:Array[alt, species, wavelength]
;
;  for i=0, nbins-1 do $
;    eiionz_ht[i]=total(reform(eiionz_stot(*,i)*dz*1.e5)) ;height integrated electron ionisation rate
;
;  eiionz_ht=eiionz_ht*1.60218e-19*ener  ;J cm-2 s-1
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
;  tflux_ht=fltarr(nbins)
;  for i=0, nbins-1 do $
;      tflux_ht(i)=total(reform(tflux(*,i)*del(i)))
;

;    cgdisplay,1000,1000, /free
;
;    cgplot,12397./ener,eiionz_ht,color='red',thick=1,/xs,/ys,psym=-46,$
;      ytitle='photons cm-2 s-1 ',xtitle='Energy',xr=[0.,2000.],yr=[0.,1e9]
;    cgplot, 12397./ener,tflux_ht,color='dodger blue',/overplot,thick=3
;    cglegend, titles=['Height Integrated Photoelectron ','Tflux'],$
;      psyms=[15], symcolors=['red','dodger blue'],$
;      alignment=1,/box,Location=[0.93, 0.9]


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



;output- for running the entire solar flux in one go
;
;pepi=fltarr(nmaj,lmax)
;pepi=eiionz_transp/photoi
;pepi[where(~finite(pepi),/null)]=0.



;save, photoi_wv_sza,photoi_sza,eiionz_sza,eiionz_transp_sza,sza_arr, lotlon_arr,zz, wv1,wv2, ssflux,$
;  filename = floc+"SQ_peak_m5_sza.sav"

save, photoi_wv_sza,photoi_sza,eiionz_sza,eiionz_transp_sza,sza_arr,lotlon_arr, zz, wv1,wv2, ssflux,$
    filename = floc+"SQ_peak_x9_sza.sav"

;_________________________________________________________________________________________________________________________________________________________________________________________


end
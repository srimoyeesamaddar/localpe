pro localpe_ssflux;,hi
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
; see this file for definitions

;___________________________________________________________________________________________________________________________________________________________________         
;#1 NRL spectrum 

;filename='/home/srimoyee/Desktop/nrl_files/sample_spectra.sav'
;filename='/home/srimoyee/Desktop/nrl_files/sav_files/irradiance_23jul2016_myver.sav'
;filename='/home/srimoyee/Desktop/nrl_files/sav_files/irradiance_04nov2003_myver.sav'

;restore, filename
;wv1=wavelength_low
;wv2=wavelength_high
;ssflux=spectrum2
;___________________________________________________________________________________________________________________________________________________________________

;#2 NRl Spectrum
;
;filename='/home/srimoyee/Desktop/nrl_files/sav_files/x28_custombins_myver.sav'
;restore, filename
;wv1=wavelength_low
;wv2=wavelength_high
;ssflux=spectrum[90,*] ;for x_28custombins_myver.sav file only
;___________________________________________________________________________________________________________________________________________________________________
;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;restore, '/home/srimoyee/Desktop/nrl_files/sav_files/SQ_newspectra.sav'
;
;;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5

;;meanpeakx9[0:6]=meanpeakx9[0:6]/1.e5
;;meanpeakx9[7:8]=meanpeakx9[7:8]/1.e4
;;meanpeakx9[9:12]=meanpeakx9[9:12]/1.e2
;;meanpeakx9[13:14]=meanpeakx9[13:14]/1.e1



;meanpeakx9[0:5]=meanpeakx9[0:5]/1.e5
;meanpeakx9[6:8]=meanpeakx9[6:8]/1.e4
;meanpeakx9[9:11]=meanpeakx9[9:11]/1.e3
;meanpeakx9[12]=meanpeakx9[12]/1.e2
;meanpeakx9[13:-1]=meanpeakx9[13:-1]/10.
;
;wv1=wave_gcm1
;wv2=wave_gcm2
;ssflux=meanpeakx9

;ssflux=peakm5
;if (hi eq 0) then ssflux[hi]= 2.e9
;if ((hi ge 0) and (hi le 5)) then ssflux[hi]= 2.e9
;


;Quiet spectrum with NRLFLARE resolution

restore,'/home/srimoyee/Desktop/nrl_files/sav_files/NRL_QuietSpectrum_July2022.sav'

;File has :nrl_wave1, nrl_wave2, nrl_current
wv1=nrl_wave1
wv2=nrl_wave2
ssflux=nrl_current



;mask=fltarr(n_elements(wv1))
;mask[hi]=1.    
;ssflux=ssflux*mask

;ssflux=mask  ;checking with solar flux values all set to 1
;___________________________________________________________________________________________________________________________________________________________________


; #4 EUVAC- from Solomon and Qian 2005
;filename='/home/srimoyee/Desktop/nrl_files/sav_files/Solomon_Qian2005.sav'
;;file has:A_fac, fref, wavelength_low, wavelength_high
;restore, filename
;
;f107=70
;f107a=70
;P107 = (F107+F107A)/2.
;n_wvl=n_elements(a_fac)
;ssflux= fltarr(n_wvl)
;for L=0,n_wvl-1 do begin 
;          SSFLUX(L) = fref(L) * (1. + A_fac(L)*(P107-80.))
;          IF (SSFLUX(L) LT 0.8*fref(L)) then SSFLUX(L) = 0.8*fref(L)
;endfor
;
;wv1=wavelength_low
;wv2=wavelength_high

;mask=fltarr(22)
;mask[hi]=1.
;ssflux=ssflux*mask
;......................................................................
;......................................................................
;;Adding EUVAC spectra to higher wavlengths in NRL spectra for testing
;;NRL does not have values in bins higer than 450 A ,.i.e., from index 81  
;
;;EUVAC model solar spectrum
;n_wvl=123
;waves=fltarr(n_wvl)
;wavel=fltarr(n_wvl)
;rflux=fltarr(n_wvl)  ;reference flux
;a=fltarr(n_wvl)
;
;euv_data=fltarr(4,n_wvl)
;openr,lun_euvac,'ssflux_euvac.dat',/get_lun
;skip_lun, lun_euvac, 1, /lines ; 1 header line
;readf,lun_euvac,euv_data
;waves=reform(euv_data[0,*])
;wavel=reform(euv_data[1,*])
;rflux=reform(euv_data[2,*])
;a=reform(euv_data[3,*])
;free_lun,lun_euvac
;
;P107 = (F107+F107A)/2.
;euv_sflux=fltarr(n_wvl) ; EUV solar spectrum
;
;for L=0,n_wvl-1 do begin
;          EUV_SFLUX(L) = RFLUX(L) * (1. + A(L)*(P107-80.))
;          IF (EUV_SFLUX(L) LT 0.8*RFLUX(L)) then EUV_SFLUX(L) = 0.8*RFLUX(L)
;endfor
;
;;Adding EUVAC to NRL spectrum for higher wavelengths
;ssflux[81:-1]=euv_sflux[48:-1]  ;FINAL SPECTRUM



; Making some bins..................for this keep maskk=0 
;ssflux_copy=ssflux
;ssflux=ssflux*0.0

;#1 soft x-rays less than 300A
;ssflux[0:65]=ssflux_copy[0:65]

;#2  He II 304A
;ssflux[66]=ssflux_copy[66] 
;
;;#3 greater than 304A
;ssflux[67:-1]=ssflux_copy[67:-1]  
;
;;#4 Lyman beta
;ssflux[138]=ssflux_copy[138]
;
;;#5 Lyman alpha
;ssflux[144]=ssflux_copy[144]

;#6 Auger N2
;ssflux[40]=ssflux_copy[40]

;___________________________________________________________________________________________________________________________________________________________________
;; SQ05 with and without Fe line
;file='/home/srimoyee/Desktop/nrl_files/sav_files/x9_stanbands_myver.sav'
;;file='/home/srimoyee/Desktop/nrl_files/sav_files/m5_stanbands_myver.sav'
;; file has: wave_lo, wave_hi,spectrum
;
;restore, file
;;for SQ' 05 + Fe line
;wv1=wave_lo
;wv2=wave_hi
;ssflux=spectrum[*,720] 
;mask=fltarr(38)
;mask[hi]=1.
;ssflux=ssflux*mask

;___________________________________________________________________________________________________________________________________________________________________

; FISM Spectra

;filenam= '/home/srimoyee/Desktop/nrl_files/sav_files/fism_1606_253sep17.txt'
;
;openr,lun_fism, filenam, /get_lun
;n_wv=32   ;wavelength bins
;wave1=fltarr(n_wv)  ;in nm ?
;wave2=fltarr(n_wv)  ;in nm ?
;fismflux=dblarr(n_wv)  ;in  photons/cm2/s
;skip_lun, lun_fism, 2, /LINES
;count=0
;
;while not eof(lun_fism) do begin
;  readf,lun_fism,a,b,c
;  wave1(count)=a*10.0
;  wave2(count)=b*10.0
;  fismflux(count)=c
;  count++
;endwhile
;free_lun,lun_fism
;
;wv1=wave1
;wv2=wave2
;ssflux=fismflux
;
;
;;mask=fltarr(n_elements(wv1))
;;mask[hi]=1.
;;ssflux=ssflux*mask
;
;___________________________________________________________________________________________________________________________________________________________________
lmax=n_elements(wv1)




return
end
  
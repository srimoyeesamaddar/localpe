

fname="/home/srimoyee/Desktop/nrl_files/sav_files/X9_NC_17JUL2022.sav"      ; X9, New
restore,fname

; file has:
; Sp_exctrate1_transp(nstate,nmaj,jmax,nbins,lmax),
; Sp_ionzrate1_transp(nstate,nmaj,jmax,nbins,lmax),
; Sp_exctrate1_local(nstate,nmaj,jmax,nbins,lmax),
; Sp_ionzrate1_local(nstate,nmaj,jmax,nbins,lmax),
; Sp_photoi(nmaj,jmax,lmax),
; Sp_photoki(nmaj,pstate,jmax,lmax),
; tflux_all(jmax,nbins,lmax),
; pespec_all(jmax,nbins,lmax),
; zz,ener,del,wv1,wv2,sza,lat,lon,tau,sigix,sigex,zmaj

out_fname="/home/srimoyee/Desktop/nrl_files/sav_files/N2_vibcross_JUL2022.csv" 
N2vib_crss=reform(sigex[19,2,*])
  
  WRITE_CSV, out_fname, ener,N2vib_crss
  
  end
;idate =2016205  ;m5 flare
idate=2017249   ;x9 flare



;hour =5 ; m5 flare UT time
hour=12 ; x9 flare UT time

doy=idate-(long(idate/1000)*1000)
utsec=3600.*hour
utdoy=doy+utsec/(3600.*24.) ; ut in days since Jan 1 of year (input to zenith.pro)



;
;lat=0.+findgen(17)*5.; m5 flare
;lon=20.+findgen(16)*5.

lat=0.+findgen(16)*5.;x9 flare
lon=0.+findgen(18)*5.

ij=0

sza=fltarr(n_elements(lon),n_elements(lat))
for j=0, n_elements(lat)-1 do begin  ;
   for i=0,n_elements(lon)-1 do begin

    sza[i,j]=zenith(utdoy,lat[j],lon[i])*!radeg;
    

    print,ij,lat[j],lon[i],sza[i,j]
    ij++
  endfor
  
endfor
;min_sza=min(sza,ii)
;print, ii, min_sza
end
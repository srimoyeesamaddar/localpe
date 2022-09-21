pro compare_x,ener,del,siga,sec,sigloss

openr,lun,'siga.dat',/get_lun
readf,lun,glow_nbins
dum=fltarr(glow_nbins)
glow_ener=fltarr(glow_nbins)
glow_del=fltarr(glow_nbins)
glow_siga=fltarr(3,glow_nbins,glow_nbins)
readf,lun,glow_ener
readf,lun,glow_del
for i=0,2 do begin
for j=0,glow_nbins-1 do begin
  readf,lun,dum
  glow_siga[i,j,*]=dum
endfor
endfor
free_lun,lun
openr,lun,'sec.dat',/get_lun
readf,lun,glow_nbins
dum=fltarr(glow_nbins)
glow_ener=fltarr(glow_nbins)
glow_del=fltarr(glow_nbins)
glow_sec=fltarr(3,glow_nbins,glow_nbins)
readf,lun,glow_ener
readf,lun,glow_del
for i=0,2 do begin
for j=0,glow_nbins-1 do begin
  readf,lun,dum
  glow_sec[i,j,*]=dum
endfor
endfor
free_lun,lun
openr,lun,'sigloss.dat',/get_lun
readf,lun,glow_nbins
dum=fltarr(glow_nbins)
glow_ener=fltarr(glow_nbins)
glow_del=fltarr(glow_nbins)
glow_sigloss=fltarr(3,glow_nbins)
readf,lun,glow_ener
readf,lun,glow_del
for i=0,2 do begin
  readf,lun,dum
  glow_sigloss[i,*]=dum
endfor
free_lun,lun

dum=min(abs(ener-20.),i20)
dum=min(abs(glow_ener-20),glow_i20)

dum=min(abs(ener-100.),i100)
dum=min(abs(glow_ener-100),glow_i100)

dum=min(abs(ener-1000.),i1000)
dum=min(abs(glow_ener-1000),glow_i1000)

dum=min(abs(ener-10000.),i10000)
dum=min(abs(glow_ener-10000),glow_i10000)

j=0

window,10
plot,ener,siga[j,*,i100],/xlog,xr=[1e-1,1e2],/ylog,/nodata,yr=[1e-25,1e-15]
oplot,ener,siga[j,*,i20],color=!ct.navy,psym=-4
oplot,glow_ener,glow_siga[j,*,glow_i20],color=!ct.red,psym=-4
oplot,ener,siga[j,*,i100],color=!ct.navy,psym=-1
oplot,glow_ener,glow_siga[j,*,glow_i100],color=!ct.red,psym=-1
oplot,ener,siga[j,*,i1000],color=!ct.navy,psym=-2
oplot,glow_ener,glow_siga[j,*,glow_i1000],color=!ct.red,psym=-1
oplot,ener,siga[j,*,i10000],color=!ct.navy,psym=-2
oplot,glow_ener,glow_siga[j,*,glow_i10000],color=!ct.red,psym=-1

print,siga[j,0,i20],glow_siga[j,0,glow_i20]
print,siga[j,0,i100],glow_siga[j,0,glow_i100]
print,siga[j,0,i1000],glow_siga[j,0,glow_i1000]
print,siga[j,0,i10000],glow_siga[j,0,glow_i10000]
print,max(abs((siga[j,*,i20]-glow_siga[j,*,glow_i20]))),!c
print,max(abs((siga[j,*,i100]-glow_siga[j,*,glow_i100]))),!c
print,max(abs((siga[j,*,i1000]-glow_siga[j,*,glow_i10000]))),!c
print,max(abs((siga[j,*,i1000]-glow_siga[j,*,glow_i10000]))),!c
for i=0,24 do print,i, ener[i],siga[j,i,i20],glow_siga[j,i,glow_i20],siga[j,i,i100],glow_siga[j,i,glow_i100]

window,11
plot,ener,sec[j,*,i100],/xlog,xr=[1e-1,1e4],/ylog,/nodata,yr=[1e-24,1e-14]
oplot,ener,sec[j,*,i100],color=!ct.navy
oplot,glow_ener,glow_sec[j,*,glow_i100],color=!ct.red
oplot,ener,sec[j,*,i1000],color=!ct.navy
oplot,glow_ener,glow_sec[j,*,glow_i1000],color=!ct.red
oplot,ener,sec[j,*,i10000],color=!ct.navy
oplot,glow_ener,glow_sec[j,*,glow_i10000],color=!ct.red

window,12
plot,ener,sigloss[j,*],/xlog,xr=[1e-1,1e5],/ylog,/nodata,yr=[1e-24,1e-14]
oplot,ener,sigloss[j,*],color=!ct.navy
oplot,glow_ener,glow_sigloss[j,*],color=!ct.red


;stop
return
end

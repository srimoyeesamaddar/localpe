
estudy=185.5

dum=min(abs(e190-estudy),i190)
dum=min(abs(e3600-estudy),i3600)

window,7
ct=3
cs=2

plot,e3600,sec3600[0,*,i3600],charsize=cs,charthic=ct,xthick=ct,ythick=ct,$
xtitle='Energy',ytitle='Cross Section (cm!u2!n)',back=!snoe.ct.fg,color=!snoe.ct.bg,$
/nodata,xr=[120,200],yr=[0,8e-17]


for j=2,2 do oplot,e3600,sec3600[j,*,i3600],thick=ct,color=!snoe.ct.navy,psym=10
for i=5,5 do for j=2,2 do oplot,e3600,sigix3600[i,j,*],color=!snoe.ct.navy,psym=10

for j=2,2 do  oplot,e190,sec190[j,*,i190],thick=ct,color=!snoe.ct.red,psym=10

xyouts,.8,.8,strtrim(long(estudy),2),/norm,charsize=cs,charthick=ct,color=!snoe.ct.navy

print,total(sec190[2,*,i190]*d190),total(sec3600[2,*,i3600]*d3600)
end
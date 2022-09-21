;cd,'/home/padma/togo/Rocket programs/'
;.run read_ephotn2.pro
;get_lun,unit
file='ephot_xn2.txt'
openr,lun,file,/get_lun

s='a string'
for i=1,4 do readf,lun,s


format='$(f8.2,2x,f8.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f9.5,1x,f9.5)'

maxn=3000.

n=0

;a=(b=(c=(d=(e=(f=(j1=(k1=(k2=(k3=maxn)))))))))
 wv1=(wv2=(xn2=(an2=(bn2=(cn2=(fn2=(dissn2=(tionn2=(tabsn2=dblarr(maxn))))))))))
while not eof(lun) do begin
 readf,lun,a,b,c,d,e,f,j1,k1,k2,k3,format=format
 wv1(n)=a
 wv2(n)=b
 xn2(n)=c
 an2(n)=d
 bn2(n)=e
 cn2(n)=f
 fn2(n)=j1
 dissn2(n)=k1
 tionn2(n)=k2
 tabsn2(n)=k3
 n=n+1
;print,n

endwhile



free_lun,lun
nn=n-1
wv1=wv1(0:nn)
wv2=wv2(0:nn)
xn2=xn2(0:nn)
an2=an2(0:nn)
bn2=bn2(0:nn)
cn2=cn2(0:nn)
fn2=fn2(0:nn)
dissn2=dissn2(0:nn)
tionn2=tionn2(0:nn)
tabsn2=tabsn2(0:nn)

aa=wv1
bb=wv2
probn2=dblarr(6,123)
probn2(0,*)=xn2
probn2(1,*)=an2
probn2(2,*)=bn2
probn2(3,*)=cn2
probn2(4,*)=fn2
probn2(5,*)=dissn2
sigin2=tionn2
sigan2=tabsn2

savefile='read_ephotn2'+'.sav'
save,file=savefile,aa,bb,probn2,sigin2,sigan2
end

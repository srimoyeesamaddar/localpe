 
;cd,'/home/padma/togo/Rocket programs/'
;.run read_ephotoO2.pro
;get_lun,unit
file='ephoto_xo2.dat'
openr,lun,file,/get_lun

s='a string'
for i=1,4 do readf,lun,s


format='$(f8.2,2x,f8.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f9.5,1x,f9.5)'

maxn=3000.

n=0

;a=(b=(c=(d=(e=(f=(j1=(k1=(k2=(k3=maxn)))))))))
 aa=(bb=(bro2=(st1o2=(st2o2=(st3o2=(st4o2=(st5o2=(tiono2=(tabso2=dblarr(maxn))))))))))

;probo2=dblarr(6,123
while not eof(lun) do begin
 readf,lun,a,b,c,d,e,f,j1,k1,k2,k3,format=format
 aa(n)=a
 bb(n)=b
 bro2(n)=c
 st1o2(n)=d
 st2o2(n)=e
 st3o2(n)=f
 st4o2(n)=j1
 st5o2(n)=k1
 tiono2(n)=k2
 tabso2(n)=k3
 n=n+1
;print,n

endwhile



free_lun,lun
nn=n-1
aa=aa(0:nn)
bb=bb(0:nn)
bro2=bro2(0:nn)
st1o2=st1o2(0:nn)
st2o2=st2o2(0:nn)
st3o2= st3o2(0:nn)
st4o2= st4o2(0:nn)
st5o2= st5o2(0:nn)
tiono2=tiono2(0:nn)
tabso2=tabso2(0:nn)

probo2=dblarr(6,123)
probo2(0,*)=bro2
probo2(1,*)=st1o2
probo2(2,*)=st2o2
probo2(3,*)=st3o2
probo2(4,*)=st4o2
probo2(5,*)=st5o2
sigio2=tiono2
sigao2=tabso2

savefile='read_ephotoO2'+'.sav'
save,file=savefile,aa,bb,probo2,sigio2,sigao2

end
 

;cd,'/home/padma/togo/Rocket programs/'
;.run read_ephotoO.pro
;get_lun,unit
file='ephoto_xo.dat'
openr,lun,file,/get_lun

s='a string'
for i=1,4 do readf,lun,s


format='$(f8.2,2x,f8.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f9.5,1x,f9.5)'

maxn=3000.

n=0

;a=(b=(c=(d=(e=(f=(j1=(k1=(k2=(k3=maxn)))))))))
 wv1o=(wv2o=(st1o=(st2o=(st3o=(st4o=(st5o=(st6o=(tiono=(tabso=dblarr(maxn))))))))))
while not eof(lun) do begin
 readf,lun,a,b,c,d,e,f,j1,k1,k2,k3,format=format
 wv1o(n)=a
 wv2o(n)=b
 st1o(n)=c
 st2o(n)=d
 st3o(n)=e
 st4o(n)=f
 st5o(n)=j1
 st6o(n)=k1
 tiono(n)=k2
 tabso(n)=k3
 n=n+1
;print,n

endwhile



free_lun,lun
nn=n-1
wv1o=wv1o(0:nn)
wv2o=wv2o(0:nn)
st1o=st1o(0:nn)
st2o=st2o(0:nn)
st3o= st3o(0:nn)
st4o= st4o(0:nn)
st5o= st5o(0:nn)
st6o=st6o(0:nn)
tiono=tiono(0:nn)
tabso=tabso(0:nn)

aa=wv1o
bb=wv2o
probo=dblarr(6,123)
probo(0,*)=st1o
probo(1,*)=st2o
probo(2,*)=st3o
probo(3,*)=st4o
probo(4,*)=st5o
probo(5,*)=st6o
sigio=tiono
sigao=tabso

savefile='read_ephotoO'+'.sav'
save,file=savefile,aa,bb,probo,sigio,sigao
end
 

;Main photoelectron calling function

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport  
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
    
;............Model constants....................
nbins=180
nmaj=3
;lmax=139
lmax=156
;lmax=123
;lmax=22 ;Solomon and Qian 2005
lon=0.0
;altitude bins
z0=40.
zz=findgen(161)*2.+z0  
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]


;..............Inputs..................


;idate=1999266
idate=20031104
lat=0.

;f107=137.8
;f107a=157.077

f107=70.
f107a=70.
ap=15
;binsize=178 

;.....................................................
floc='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\'
;ploc=
;fname='Solomon&Qian2005_allrun'
fname='mjuly'
;.....................................................

tic

localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi
 
endif

;DEFAULT solar flask mask inputs

maskk=0  ;deafult values ;make maskk=1 when using one flux bin at a time
mask=0

;Solar flux inputs

    localpe_ssflux,maskk=maskk,mask

    dum=min(abs(ener-20.),ia20)
    dum=min(abs(ener-100.),ia100)
    dum=min(abs(ener-1000.),ia1000)

    n_hours=14
    glow20=fltarr(n_hours)
    glow100=fltarr(n_hours)
    glow1000=fltarr(n_hours)
    ace20=fltarr(n_hours)
    ace100=fltarr(n_hours)
    ace1000=fltarr(n_hours)
    asza=fltarr(n_hours)

    ;for i=0,n_hours-1 do begin
    ;  for i=0, 3 do begin
    i=6.
    hour=float(i)+6.
    utsec=3600.*hour
    localpe_neutralatm
    asza[i]=sza
    localpe_photoionz
  

    localpe_approxeden,zz,f107,eden,etemp
    ;eden=eden/8.0
    localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo


    ;photoelectron flux transport calculation
    ace_etransport
    

;.....................................
;z1=110.
;z2=250.
;dum1l=min(abs(local_zz-z1),ilz1)
;print,dum1l
;dum1l=min(abs(local_zz-z2),ilz2)
;print,dum1l
;cs=2
;ct=3

;window,i
;plot,local_ener,local_pespec[ilz1,*],xr=[1.,1000.],/yl,yr=[1e2,1e10],/xl,psym=10,/nodata,$
;charsize=cs,charthick=ct,xthick=ct,ythick=ct,color=!snoe.ct.bg,back=!snoe.ct.fg
;
;oplot,local_ener,local_pespec[ilz1,*],psym=10,color=!ct.navy,thick=ct*2
;oplot,local_ener,local_pespec[ilz2,*],psym=10,color=!ct.navy,thick=ct
;oplot,local_ener,tflux[ilz1,*],psym=10,color=!ct.blue,thick=ct*2
;oplot,local_ener,tflux[ilz2,*],psym=10,color=!ct.blue,thick=ct


;.........................................................................................

toc

; Energy conservation tests for photoionisation and photoelectrons


;the height integrated photoionisation energy should add up to the solar flux, 
;atleast for the lower wavlength bins 

photoi_ht=dblarr(lmax)  ;height integrated photoionisation energy

for j=0, nmaj-1 do begin
  for l=0, lmax-1 do begin
       photoi_ht(l)=photoi_ht(l)+total(photoi_wv(*,j,l)*dz*1.e5) ;height integrated ionisation rate   
  endfor
endfor


;...........................................................................................

;Plot variables...............

colors=['red','dark green','dodger blue','yellow']
plt_sym=0
lstyl=[2,1,0]
pos1=[0.1,0.5,0.45,0.95]
pos2=[0.1,0.1,0.45,0.45]
pos3=[0.47,0.5,0.95,0.95]
pos4=[0.47,0.1,0.95,0.45]
thck=5.
leg_loc=[0.90, 0.90]

head_tilt='Solomon & Qian 2005 ' 
;head_tilt='X28 '
sav_loc="C:\Users\Srimoyee\Desktop\agu_plots\"


;Electron Impact Ionisation plot 

xr=[1.e-5,1.e5]
yr=[min(zz),max(zz)]
xt='Electron impact ionization (Pe)!C !C(cm!U-3!N s!U-1!N) '  ;!C Black- local, Red-transport
yt='Altitude (km)'
tilt=head_tilt+'Eiionz'
maj=['O','O2','N2']
leg_loc=[0.90, 0.85]

w1 = WINDOW(WINDOW_TITLE=tilt, $
             DIMENSIONS=[800,800])

;local
ei1=plot(eiionz_transp[0,*],zz,name=maj[0],$
  xtitle=xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  thick=thck,color=colors[0],symbol=plt_sym,/current)

ei2=plot(eiionz_transp[1,*],zz,name=maj[1],$
  thick=thck,color=colors[1],symbol=plt_sym,/overplot)

ei3=plot(eiionz_transp[2,*],zz,name=maj[2],$
  thick=thck,color=colors[2],symbol=plt_sym,/overplot)


;;transport
;eit1=plot(eiionz_transp[0,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;eit2=plot(eiionz_transp[1,*],zz,$
;  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)
;
;eit3=plot(eiionz_transp[2,*],zz,$
;  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg1= LEGEND(TARGET=[ei1,ei2,ei3], POSITION=leg_loc,/normal)
w1.Save, sav_loc+"EX28_90_eiionz.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;photoionisation


xr=[1.e-5,1.e4]
yr=[min(zz),max(zz)]
xt='Photoionisation (Pi)!C !C(cm!U-3!N s!U-1!N)'
yt='Altitude (km)'
tilt='SQ Photoi'
leg_loc=[0.92,0.85]

w2= WINDOW(WINDOW_TITLE=tilt, $
           DIMENSIONS=[800,800])


p1=plot(photoi[0,*],zz,name=maj[0],$
              xtitle=xt,ytitle=yt,$
              xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
              thick=thck,color=colors[0],symbol=plt_sym,/current)



p2=plot(photoi[1,*],zz,name=maj[1],$
       thick=thck,color=colors[1],symbol=plt_sym,/overplot)

p3=plot(photoi[2,*],zz,name=maj[2],$
        thick=thck,color=colors[2],symbol=plt_sym,/overplot)


leg2 = LEGEND(TARGET=[p1,p2,p3], POSITION=leg_loc,/normal)

w2.Save, sav_loc+"X28_90_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;.............................................................................
;Ratio of tflux(transport) and pespec(local) wrt altitude

dif=tflux/pespec

xr=[1.e-3,1.e0]
yr=[min(zz),max(zz)]
xt='trasport flux/local flux'
yt='Altitude'
tilt=head_tilt+'tflux/pespec vs alt'

w3= WINDOW(WINDOW_TITLE=tilt, $
            DIMENSIONS=[800,800])

ener_pos=[29,47,170]

rz1=plot(dif[*,ener_pos[0]],zz,name=strtrim(string(ener[ener_pos[0]]),1)+'eV',$
      xtitle=xt,ytitle=yt,$
      xstyle=1,ystyle=1,yrange=yr_alt,xrange=xr,xlog=1,$;
      color=colors[0],thick=thck,/current)

rz2=plot(dif[*,ener_pos[1]],zz,name=strtrim(string(ener[ener_pos[1]]),1)+'eV',$
      color=colors[1],thick=thck,/overplot)

rz3=plot(dif[*,ener_pos[2]],zz,name=strtrim(string(ener[ener_pos[2]]),1)+'eV',$
      color=colors[2],thick=thck,/overplot)


leg3= LEGEND(TARGET=[rz1,rz2,rz3], POSITION=leg_loc,/normal)

;Ratio of tflux(transport) and pespec(local) wrt energy



yr=[0.0,5.]
xr=[10.,max(ener)]
yt='trasport flux/local flux'
xt='energy (eV)'
tilt=head_tilt+'tflux/pespec vs energy'

alt_pos=[10,47,100]


w4= WINDOW(WINDOW_TITLE=tilt, $
            DIMENSIONS=[800,800])


re1=plot(ener,dif[alt_pos[0],*],name=strtrim(string(zz[alt_pos[0]]),1)+'km',$
          xtitle=xt,ytitle=yt,$
          xstyle=1,ystyle=1,xrange=xr,xlog=1,yrange=yr,$;
          color=colors[0],thick=thck,/current)


re2=plot(ener,dif[alt_pos[1],*],name=strtrim(string(zz[alt_pos[1]]),1)+'km',$
         color=colors[1],thick=thck,/overplot)

re3=plot(ener,dif[alt_pos[2],*],name=strtrim(string(zz[alt_pos[2]]),1)+'km',$
        color=colors[2],thick=thck,/overplot)

leg4 = LEGEND(TARGET=[re1,re2,re3], POSITION=leg_loc,/normal)



;...................................................................................
;Pe/Pi ratio

xr=[0.1,1000.]
yr=[min(zz),max(zz)]
xt='Pe/Pi  '
yt='Altitude (km)'
tilt=head_tilt+'Pe/Pi'


w5 = WINDOW(WINDOW_TITLE=tilt, $
              DIMENSIONS=[800,800])


pp1=plot(eiionz_transp[0,*]/photoi[0,*],zz,name=maj[0],$
    xtitle=xt,ytitle=yt,$
    xstyle=1,ystyle=1,yrange=yr,xlog=1,xrange=xr,$;
    thick=thck,color=colors[0],symbol=plt_sym,/current)



pp2=plot(eiionz_transp[1,*]/photoi[1,*],zz,name=maj[1],$
  thick=thck,color=colors[1],symbol=plt_sym,/overplot)

pp3=plot(eiionz_transp[2,*]/photoi[2,*],zz,name=maj[2],$
  thick=thck,color=colors[2],symbol=plt_sym,/overplot)


;;transport
;ppt1=plot(eiionz_transp[0,*]/photoi[0,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;ppt2=plot(eiionz_transp[1,*]/photoi[1,*],zz,$
;  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)
;
;ppt3=plot(eiionz_transp[2,*]/photoi[2,*],zz,$
;  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg5 = LEGEND(TARGET=[pp1,pp2,pp3], POSITION=[0.90,0.85],/normal)
w5.Save, sav_loc+"X28_90_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;.................................................................................
;Electron impact ionisation flux for one altitde, for both local and transport
xr=[0.1,1.e5]
xt='energy (eV) for !C'
yt='Electron impact ionisation flux'
tilt=head_tilt+'tflux & pespec'
leg_loc1=[0.1,0.9]
leg_loc2=[0.55,0.9]
leg_loc3=[0.95,0.9]


w6 = WINDOW(WINDOW_TITLE=tilt, $
  DIMENSIONS=[800,800])

eif1=plot(ener,tflux[alt_pos[0],*],name='transport',$
      xrange=xr,histogram=1,xlog=1,ylog=1,$
      xtitle=xt+strtrim(string(zz[alt_pos[0]]),1)+'km',ytitle=yt,/current,layout=[3,1,1])
 
eif2=plot(ener,pespec[alt_pos[0],*],name='local',$
         histogram=1,color=colors[0],/overplot)
         
leg6 = LEGEND(TARGET=[eif1,eif2], POSITION=leg_loc1,/normal)
        
         
eif3=plot(ener,tflux[alt_pos[1],*],name='transport',$
         xrange=xr,histogram=1,xlog=1,ylog=1,$
         xtitle=xt+strtrim(string(zz[alt_pos[1]]),1)+'km',ytitle=yt,/current,layout=[3,1,2])

eif4=plot(ener,pespec[alt_pos[1],*],name='local',$
         histogram=1,color=colors[0],/overplot)

leg7 = LEGEND(TARGET=[eif3,eif4], POSITION=leg_loc2,/normal)


eif5=plot(ener,tflux[alt_pos[2],*],name='transport',$
         xrange=xr,histogram=1,xlog=1,ylog=1,$
         xtitle=xt+strtrim(string(zz[alt_pos[2]]),1)+'km',ytitle=yt,/current,layout=[3,1,3])

eif6=plot(ener,pespec[alt_pos[2],*],name='local',$
         histogram=1,color=colors[0],/overplot)

leg8 = LEGEND(TARGET=[eif5,eif6], POSITION=leg_loc3,/normal)

         
;Electron impact ionsation and column density of N2

xr=[0.1,1.e5]
xt='Column Density of N2'
yt='Altitude (km)'
tilt=head_tilt+'col dens and eiionz @ 1e4 eV'
leg_loc=[0.9,0.45]

w7=Window(window_title=tilt,$
          dimension=[800,800])
z1=plot(zcol(2,*),zz,$
        xtitle=xt,ytitle=yt,$
        xlog=1,xstyle=1,ystyle=1,$
        /current,layout=[1,2,1])




ener_pos=(where(ener ge 1.e4 and ener le 2.e4))[2]

xt='Electron Impact ionisation of N2 @'+strtrim(string(ener[ener_pos]))+'eV'
ein21=plot(eiionz_wv(*,2,ener_pos),zz,name='Local',$
      xtitle=xt, ytitle=yt,$
      xlog=1,xstyle=1,ystyle=1,$
      /current,layout=[1,2,2])

ein22=plot(eiionz_transp_wv(*,2,ener_pos),zz,name='Transport',$
           color=colors[0],/overplot)


leg9 = LEGEND(TARGET=[ein21,ein22], POSITION=leg_loc,/normal)

;Mean free path calculation
mfp= dblarr(nmaj+1,jmax,nbins)
for i=0,jmax-1 do begin
  for j=0,nbins-1 do begin
   mfp[0,i,j]=1.e-5/total(zmaj[0,i]*(sigix[*,0,j]+sigex[*,0,j]))
   mfp[1,i,j]=1.e-5/total(zmaj[1,i]*(sigix[*,1,j]+sigex[*,1,j]))
   mfp[2,i,j]=1.e-5/total(zmaj[2,i]*(sigix[*,2,j]+sigex[*,2,j]))
  endfor
endfor
mfp[3,*,*]=total(mfp,1) ;total mean free path of all species

xr=[1.,1.e7]
yr=[min(zz),max(zz)]

xt='Mean free path (km)for '
yt='Altitude (km)'
tilt=head_tilt+'mean free path'
leg_loc1=[0.1,0.9]
leg_loc2=[0.55,0.9]
leg_loc3=[0.95,0.9]
ener_pos=[29,47,170]

w8=window(window_title=tilt,$
          dimension=[800,800])
          
m1=plot(mfp[0,*,ener_pos[0]],zz,name=maj[0],$
        xtitle=xt+strtrim(string(ener[ener_pos[0]]),1)+'eV',ytitle=yt,$
        yrange=yr,xlog=1,xstyle=1,ystyle=1,$
        color=colors[0],/current,layout=[3,1,1])
        
m2=plot(mfp[1,*,ener_pos[0]],zz,name=maj[1],$
        thick=thck,$
        color=colors[1],/overplot)
        
m3=plot(mfp[2,*,ener_pos[0]],zz,name=maj[2],$
        thick=thck,$
        color=colors[2],/overplot)
        
tm1=plot(mfp[3,*,ener_pos[0]],zz,name='total',$
        thick=thck,$
        color=colors[3],/overplot)

leg10=legend(target=[m1,m2,m3,tm1],position=leg_loc1,/normal)   
;w8.Save, "C:\Users\Srimoyee\Desktop\agu_plots\test.jpeg", BORDER=10, RESOLUTION=300, /TRANSPARENT

;xr=[1.e-3,1.e4]
m4=plot(mfp[0,*,ener_pos[1]],zz,name=maj[0],$
  xtitle=xt+strtrim(string(ener[ener_pos[1]]),1)+'eV',ytitle=yt,$
  xlog=1,xrange=xr,yrange=yr,xstyle=1,ystyle=1,$
  color=colors[0],/current,layout=[3,1,2])

m5=plot(mfp[1,*,ener_pos[1]],zz,name=maj[1],$
  thick=thck,$
  color=colors[1],/overplot)

m6=plot(mfp[2,*,ener_pos[1]],zz,name=maj[2],$
  thick=thck,$
  color=colors[2],/overplot)

tm2=plot(mfp[3,*,ener_pos[1]],zz,name='total',$
    thick=thck,$
    color=colors[3],/overplot)


leg11=legend(target=[m4,m5,m6,tm2],position=leg_loc2,/normal)

;xr=[1.e-3,1.e4]
m7=plot(mfp[0,*,ener_pos[2]],zz,name=maj[0],$
  xtitle=xt+strtrim(string(ener[ener_pos[2]]),1)+'eV',ytitle=yt,$
  xlog=1,xrange=xr,yrange=yr,xstyle=1,ystyle=1,$
  color=colors[0],/current,layout=[3,1,3])

m8=plot(mfp[1,*,ener_pos[2]],zz,name=maj[1],$
  thick=thck,$
  color=colors[1],/overplot)

m9=plot(mfp[2,*,ener_pos[2]],zz,name=maj[2],$
  thick=thck,$
  color=colors[2],/overplot)


tm3=plot(mfp[3,*,ener_pos[2]],zz,name='total',$
    thick=thck,$
    color=colors[3],/overplot)

leg12=legend(target=[m7,m8,m9,tm3],position=leg_loc3,/normal)     




w9=window(window_title='Solar Flux',$
  dimension=[800,800])

s=plot((wv1+wv2)/2., ssflux, $
  ylog=1,ystyle=1,xlog=1,xstyle=1,$
  xrange=[0.5,1.75e3],yrange=[1.e4,1.e12],$
  thick=thck,histogram=1,color='blue',/current,$ ;, position=[0.18,0.1,0.95,0.85]
  xtitle='Wavlength ('+'!3' + STRING(197B) + '!X'+')' ,$
  ytitle='Solar Flux (Photons cm!U-2!N s!U-1!N)')
  
  w9.Save, sav_loc+"X28_90_flux.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
  ;idl save file containing all the data for the plots
save, eiionz,eiionz_transp,eiionz_wv,eiionz_transp_wv,photoi,tflux,pespec,ener,mfp,zz,wv1,wv2,ssflux,zcol,$
  filename=floc+fname+'.sav'


;....................Plots...................................................................
;loc='C:\Users\Srimoyee\Desktop\nrl_files\plots2\'
;fname=['X28.ps','X28.pdf']
;gspath='C:\Program Files\gs\gs9.27'



;cgps_open, filename=loc+fname[0]

cgdisplay,1000,1000, /free
cgplot,(wv1+wv2)/2., ssflux, $
       /yl,/ys,/xl,/xs,xr=[0.5,1.75e3],yr=[1.e0,1.e12],$
       thick=5,psym=10,aspect=1.,color='red',charthick=2,$ ;, position=[0.18,0.1,0.95,0.85]
       xtitle='!C !C Wavlength ($\Angstrom$)',$
       ytitle='Solar Flux (Photons cm$\up-2$ s$\up-1$)'
 cgplot, (wv1+wv2)/2.,photoi_ht,color='dodger blue',/overplot       



cgdisplay,1000,1000,/free
cgplot, photoi(0,*), zz ,$ ;0 for O, 1 for O2, 2 for N2
       /xs,/ys,/xl,yr=[40.,360.],xr=[1.,1.e4],$;,xr=[1.e-8,1.e4],
       thick=5,charthick=2,aspect=1.,color='red',$
       xtitle='Photoionization rates(cm$\up-3$ s$\up-1$)',$
       ytitle='Altitude (km)'
cgplot, photoi(1,*), zz,/overplot,$
       color='charcoal',thick=5
cgplot,photoi(2,*), zz,/overplot,$
      color='dodger blue',thick=5
cglegend, titles=['O','O$\down2$','N$\down2$'],$
      psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
      alignment=1,/box,Location=[0.9, 0.9]



cgdisplay, 1000,1000, /free
cgplot,eiionz(0,*), zz ,$
      /xs,/ys,/xl,yr=[40.,250.],xr=[1.,1.e4],$ ;xr=[1.e-6,1.e4],
      color='red', thick=5.,linestyle=2,charthick=2, $ 
      xtitle='Phototelectron ionization rates (cm$\up-3$ s$\up-1$)',$
      ytitle='Altitude (km)',$
      title='Photoelectron rate with localpe'
cgplot,eiionz(1,*), zz,/overplot,$
       color='charcoal',thick=5.,linestyle=2
cgplot,eiionz(2,*), zz,/overplot,$
      color='dodger blue',thick=5.,linestyle=2
cglegend, titles=['O','O$\down2$','N$\down2$'],$
      psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
      alignment=1,/box,Location=[0.9, 0.9]


cgdisplay, 1000,1000, /free
    cgplot,eiionz_transp(0,*), zz ,$
      /xs,/ys,/xl,yr=[40.,250.],xr=[1.,1.e4],$ ;xr=[1.e-6,1.e4],
      color='red', thick=5.,linestyle=2,charthick=2, $
      xtitle='Phototelectron ionization rates (cm$\up-3$ s$\up-1$)',$
      ytitle='Altitude (km)',$
      title='Photoelectron rate with transport'
    cgplot,eiionz_transp(1,*), zz,/overplot,$
      color='charcoal',thick=5.,linestyle=2
    cgplot,eiionz_transp(2,*), zz,/overplot,$
      color='dodger blue',thick=5.,linestyle=2
    cglegend, titles=['O','O$\down2$','N$\down2$'],$
      psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
      alignment=1,/box,Location=[0.9, 0.9]


cgdisplay,1000,1000, /free
cgplot,eiionz(0,*)/ photoi(0,*), zz ,$
       /xs,/ys,/xl,yr=[50.,400.],xr=[0.01,1.e3],$
       color='red', charthick=2, thick=5.,$
       xtitle='Photoelectron/Photoionization',$
       ytitle='Altitude (km)',$
       title='Pe/Pi local calculation'
cgplot,eiionz(1,*)/ photoi(1,*), zz,/overplot,$
      color='charcoal',thick=5.
cgplot, eiionz(2,*)/photoi(2,*), zz,/overplot,$
      color='dodger blue',thick=5.
cglegend,titles=['O','O$\down2$','N$\down2$'],$
         psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
         alignment=1,/box,Location=[0.9, 0.9]
         
cgdisplay,1000,1000, /free
cgplot,eiionz_transp(0,*)/ photoi(0,*), zz ,$
       /xs,/ys,/xl,yr=[50.,400.],xr=[0.01,1.e3],$
       color='red', charthick=2, thick=5.,$
       xtitle='Photoelectron/Photoionization',$
       ytitle='Altitude (km)',$
       title='Pe/Pi transport calculation'
cgplot,eiionz_transp(1,*)/ photoi(1,*), zz,/overplot,$
       color='charcoal',thick=5.
cgplot, eiionz_transp(2,*)/photoi(2,*), zz,/overplot,$
        color='dodger blue',thick=5.
cglegend,titles=['O','O$\down2$','N$\down2$'],$
         psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
         alignment=1,/box,Location=[0.9, 0.9]


cgdisplay,1000,1000, /free
cgplot,eiionz_transp(2,*)/photoi(2,*), zz,$
      /xs,/ys,/xl,yr=[80.,200.],xr=[0.1,1.e3],$
     xtitle='N2 Photoelectron/Photoionization',ytitle='Altitude (km)',$
     color='charcoal',charthick=2, thick=3.
cgplot,[8.,6.,4.,1.8,1.,0.6,0.5,.45,.35,.22,.2],findgen(11)*10.+100.,/overplot,$
        color='charcoal',thick=3.,psym=-2,linestyle=1
cgplot, [11.8,10.6,8.,5.,3.,0.7,0.3,.18],[95.,100.,105.,110.,130.,140.,170.,200.],/overplot,$
       color='charcoal',thick=3.,psym=-4,linestyle=2
cglegend, titles=['VT','Buon92','TFR93'],$
        psyms=[-15,-2,-4], symcolors=['charcoal','charcoal','charcoal'],$
      alignment=1,/box,Location=[0.9, 0.9]

;cgps2pdf,loc+fname[0],$
;         loc+fname[1],$
;         gs_path=gspath
;cgps_close,/pdf


;.............Data files..................................


;openw,lun_O,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_O.dat',/get_lun
;printf, lun_O,'Altitude (km)','Photoionisation ','Photoelectron(local) ','Photoelectron(transport) ',$
;              'Pe/Pi(local) ','Pe/Pi(transport) '
;for i=0,n_elements(zz)-1 do $
;    printf,lun_O,zz(i),'  ',photoi(0,i),'  ',eiionz(0,i),'  ',eiionz_v2(0,i),'  ',$
;           eiionz(0,i)/ photoi(0,i),'  ',eiionz_v2(0,i)/ photoi(0,i)
;
;
;openw,lun_O2,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_O2.dat',/get_lun
;printf, lun_O2,'Altitude (km)','Photoionisation ','Photoelectron(local) ','Photoelectron(transport) ',$
;              'Pe/Pi(local) ','Pe/Pi(transport) '
;for i=0,n_elements(zz)-1 do $
;    printf,lun_O2,zz(i),'  ',photoi(1,i),'  ',eiionz(1,i),'  ',eiionz_v2(1,i),'  ',$
;           eiionz(1,i)/ photoi(1,i),'  ',eiionz_v2(1,i)/ photoi(1,i)
;
;openw,lun_N2,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_N2.dat',/get_lun
;printf, lun_N2,'Altitude (km)','Photoionisation ','Photoelectron(local) ','Photoelectron(transport) ',$
;              'Pe/Pi(local) ','Pe/Pi(transport) '
;for i=0,n_elements(zz)-1 do $
;    printf,lun_N2,zz(i),'  ',photoi(2,i),'  ',eiionz(2,i),'  ',eiionz_v2(2,i),'  ',$
;           eiionz(2,i)/ photoi(2,i),'  ',eiionz_v2(2,i)/ photoi(2,i)
;    
;    
;free_lun,lun_O,lun_O2,lun_N2
;strn='eiionz- photoelectron rate (local calculation);eiionz_v2- photoelectron rate(transport) ; indices 0-O, 1-O2, 2-N2'
;save,strn,ssflux,wv1,wv2,zz,photoi,eiionz,eiionz_v2,$
; filename='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_ionz.sav'


;ssflux2=ssflux
;
;save,ssflux2,wv1,wv2,filename='C:\Users\Srimoyee\Desktop\nov04.sav'
  
  
  
  
  
;  cgdisplay, 1000,1000, /free ,title='New Photoelectron rate with tansport '
;  cgplot,eiionz_v2(0,*), zz ,color='red', thick=5.,linestyle=2, $ ;0 for O, 1 for O2, 2 for N2,
;    /xs,/ys,yr=[40.,400.],/xl,xr=[1.,2.e4],charthick=2,$
;    xtitle='Phototelectron rates (photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)'
;  cgplot,  eiionz_v2(1,*), zz,color='charcoal',/overplot,thick=5.,linestyle=2
;  cgplot,  eiionz_v2(2,*), zz,color='dodger blue',/overplot,thick=5.,linestyle=2
;  cglegend, titles=['O','O$\down2$','N$\down2$'],psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;



;  cgdisplay, 1000,1000, /free, title='Pe/Pi comparisons'
;  cgplot,eiionz(0,*)/ photoi(0,*), zz ,color='red', thick=5., $ ;0 for O, 1 for O2, 2 for N2,
;    /xs,/ys,yr=[40.,150.],/xl,xr=[0.1,1.e3],charthick=2,$
;    xtitle='Pe/Pi !C local(bold) transport (dashed)',$
;     ytitle='Altitude (km)'
;  ;cgplot,  eiionz(0,*), zz,color='red',/overplot,thick=5.,linestyle=2
;  cgplot,  eiionz(1,*)/ photoi(1,*), zz,color='charcoal',/overplot,thick=5.
;  cgplot,  eiionz(2,*)/ photoi(2,*), zz,color='dodger blue',/overplot,thick=5.
;  
;
;
;  cgplot,eiionz_v2(0,*)/ photoi(0,*), zz ,color='red', thick=5.,linestyle=2,/overplot
;  cgplot,  eiionz_v2(1,*)/ photoi(1,*), zz,color='charcoal',/overplot,thick=5.,linestyle=2
;  cgplot,  eiionz_v2(2,*)/ photoi(2,*), zz,color='dodger blue',/overplot,thick=5.,linestyle=2
;  cglegend, titles=['O','O$\down2$','N$\down2$'],psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;  cglegend, titles=['O','O$\down2$','N$\down2$'],psyms=[15,15,15], symcolors=['red','charcoal','dodger blue'],$
;    alignment=1,/box,Location=[0.9, 0.9]
    
   
   
   
; tau plots
;plot variables 
thck=3


;#1

;i_init=140 ; initial wavelength index to be plotted
;i_fin=156 ; final wavelength index to be plotted
;col_ind=3
;tau_pltvec=[]; array for storing the plot objects

;max_tau=max(tau(*,i_init:i_fin-1)) ;1.e4;
;min_tau=min(tau(*,i_init:i_fin-1)) ;0.1;


;  tiltnam=strtrim(string(wv1(i_init)),1)+'-'+strtrim(string(wv2(i_fin-1)),1)
;  
;  tau_w = WINDOW(WINDOW_TITLE="Optical depth vs altitude "+tiltnam+" Angstrom", $
;                 DIMENSIONS=[800,800],BACKGROUND_COLOR='light gray')
;  nam=strtrim(string(wv1(i_init)),1)+'-'+strtrim(string(wv2(i_init)),1)
;
;  tau_plt=plot(tau(*,i_init),zz,name=nam,$
;        xtitle='Optical depth',ytitle='Altitude (km)',$
;        xr=[min_tau,max_tau],yr=[zz[0],zz[-1]],xlog=1,ylog=0,xstyle=1,ystyle=1,$  ;
;        color=!color.(col_ind),thick=thck,/current,BACKGROUND_COLOR='light gray')
;
;  tau_pltvec=[tau_pltvec,tau_plt]  
;  
;  for i=i_init+1, i_fin-1 do begin
;         
;         nam=strtrim(string(wv1(i)),1)+'-'+strtrim(string(wv2(i)),1)
;         col_ind+=3
;         tau_plt=plot(tau(*,i),zz,name=nam,$
;                       color=!color.(col_ind),thick=thck,/overplot)
;         tau_pltvec=[tau_pltvec,tau_plt]  
;
;  endfor
;     
;; Add the legend.
;leg1 = LEGEND(TARGET=tau_pltvec[0:9],POSITION=[0.55,0.90],/AUTO_TEXT_COLOR,transparency=100)
;leg2 = LEGEND(TARGET=tau_pltvec[10:-1],POSITION=[0.90,0.90],/AUTO_TEXT_COLOR,transparency=100)

end
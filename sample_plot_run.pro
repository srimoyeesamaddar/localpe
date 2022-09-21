

;plots of 156 runs of photoionisations and photoelectrons for 156 wavelength runs for X28 spectra
; and 22 runs for Solomon&Qian 2005

loc_sav='C:\Users\srimo\Desktop\nrl_files\sav_files\'
;fname='Solomon&Qian2005' ;change this each time
fname='Mjuly' 

filename0=loc_sav+fname+'.sav'
;file has:
; photoi_run,eiionz_run,eiionzt_run,pepi,pepi_t,zz,ssflux_run,wv1,wv2,tau
restore, filename0


;pepi[where(~finite(pepi),/null)]=0.
;pepi_t[where(~finite(pepi_t),/null)]=0.
;4 bands in Solomon&Qian 2005
;32-70 A
;70-155 A
;155- 224 A
;290-320 A

ii_band1=where(wv1 eq 32.)
ii_band2=where(wv1 eq 70.)
ii_band3=where(wv1 eq 155.)
ii_band4=where(wv1 eq 290.)


;ii_band1=where(wv1 eq 30.)  ;NRL
;ii_band2=where(wv1 eq 70.)
;ii_band3=where(wv1 eq 150.)
;ii_band4=where(wv1 eq 300.)
;
;

;Plot variables...............

colors=['red','black','blue']
pos1=[0.1,0.5,0.45,0.95]
pos2=[0.1,0.1,0.45,0.45]
pos3=[0.47,0.5,0.95,0.95]
pos4=[0.47,0.1,0.95,0.45]

xr=[0.1,100.]
yr_alt=[min(zz),max(zz)]
yr_tau1=[tau[where(zz eq yr_alt[0],/null),ii_band1], $
        tau[where(zz eq yr_alt[1],/null),ii_band1] ]
        
yr_tau2=[tau[where(zz eq yr_alt[0],/null),ii_band2], $
        tau[where(zz eq yr_alt[1],/null),ii_band2] ]
        
yr_tau3=[tau[where(zz eq yr_alt[0],/null),ii_band3], $
        tau[where(zz eq yr_alt[1],/null),ii_band3] ]

yr_tau4=[tau[where(zz eq yr_alt[0],/null),ii_band4], $
        tau[where(zz eq yr_alt[1],/null),ii_band4] ]
;yr_opt=[4.,0.]

wv1_x=strtrim(string(wv1[ii_band1]),1)+'-'+strtrim(string(wv2[ii_band1]),1)
wv2_x=strtrim(string(wv1[ii_band2]),1)+'-'+strtrim(string(wv2[ii_band2]),1)
wv3_x=strtrim(string(wv1[ii_band3]),1)+'-'+strtrim(string(wv2[ii_band3]),1)
wv4_x=strtrim(string(wv1[ii_band4]),1)+'-'+strtrim(string(wv2[ii_band4]),1)

xt='Pe/Pi '
yt_tau='Optical Depth'
yt_alt='Altitude'
maj=['O- ','O2- ','N2- ']
tilt='Solomon & Qian 2005'
thck=3
lstyl=[2,1,0]

leg_loc=[0.90, 0.90]



;Plots  

;# 1st band  32-70 A
w_1 = WINDOW(WINDOW_TITLE=tilt, $
             DIMENSIONS=[800,800])
             


p_band11=plot(pepi[0,*,ii_band1],zz,name=maj[0],$
          xtitle=xt+wv1_x,ytitle=yt_alt,$
          xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_alt,$
          linestyle=lstyl[0],thick=thck,/current,layout=[2,2,1])



p_band12=plot(pepi[1,*,ii_band1],zz,name=maj[1],$
        linestyle=lstyl[1],thick=thck,/overplot)
        
p_band13=plot(pepi[2,*,ii_band1],zz,name=maj[2],$
        linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band11=plot(pepi_t[0,*,ii_band1],zz,$
        linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band12=plot(pepi_t[1,*,ii_band1],zz,$
        linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band13=plot(pepi_t[2,*,ii_band1],zz,$
        linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg1 = LEGEND(TARGET=[p_band11,p_band12,p_band13], POSITION=[0.45,0.95],/normal)



;#2 70-155 A



p_band21=plot(pepi[0,*,ii_band2],zz,name=maj[0],$
  xtitle=xt+wv2_x,ytitle=yt_alt,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_alt,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,2])



p_band22=plot(pepi[1,*,ii_band2],zz,name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band23=plot(pepi[2,*,ii_band2],zz,name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band21=plot(pepi_t[0,*,ii_band2],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band22=plot(pepi_t[1,*,ii_band2],zz,$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band23=plot(pepi_t[2,*,ii_band2],zz,$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg2 = LEGEND(TARGET=[p_band21,p_band22,p_band23], POSITION=[0.95,0.95],/normal)


;#3 155-224 A



p_band31=plot(pepi[0,*,ii_band3],zz,name=maj[0],$
  xtitle=xt+wv3_x,ytitle=yt_alt,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_alt,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,3])



p_band32=plot(pepi[1,*,ii_band3],zz,name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band33=plot(pepi[2,*,ii_band3],zz,name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band31=plot(pepi_t[0,*,ii_band3],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band32=plot(pepi_t[1,*,ii_band3],zz,$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band33=plot(pepi_t[2,*,ii_band3],zz,$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg3 = LEGEND(TARGET=[p_band31,p_band32,p_band33], POSITION=[0.45,0.45],/normal)



;#2 290- 320 A



p_band41=plot(pepi[0,*,ii_band4],zz,name=maj[0],$
  xtitle=xt+wv4_x,ytitle=yt_alt,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_alt,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,4])



p_band42=plot(pepi[1,*,ii_band4],zz,name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band43=plot(pepi[2,*,ii_band4],zz,name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band41=plot(pepi_t[0,*,ii_band4],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band42=plot(pepi_t[1,*,ii_band4],zz,$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band43=plot(pepi_t[2,*,ii_band4],zz,$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg4 = LEGEND(TARGET=[p_band41,p_band42,p_band43], POSITION=[0.95,0.45],/normal)


;...................................................................................

;# 1st band  32-70 A
w_1 = WINDOW(WINDOW_TITLE=tilt, $
  DIMENSIONS=[800,800])



p_band11=plot(pepi[0,*,ii_band1],tau(*,ii_band1),name=maj[0],$
  xtitle=xt+wv1_x,ytitle=yt_tau,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_tau2,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,1])



p_band12=plot(pepi[1,*,ii_band1],tau(*,ii_band1),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band13=plot(pepi[2,*,ii_band1],tau(*,ii_band1),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band11=plot(pepi_t[0,*,ii_band1],tau(*,ii_band1),$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band12=plot(pepi_t[1,*,ii_band1],tau(*,ii_band1),$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band13=plot(pepi_t[2,*,ii_band1],tau(*,ii_band1),$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg1 = LEGEND(TARGET=[p_band11,p_band12,p_band13], POSITION=[0.45,0.95],/normal)



;#2 70-155 A



p_band21=plot(pepi[0,*,ii_band2],tau(*,ii_band2),name=maj[0],$
  xtitle=xt+wv2_x,ytitle=yt_tau,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_tau2,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,2])



p_band22=plot(pepi[1,*,ii_band2],tau(*,ii_band2),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band23=plot(pepi[2,*,ii_band2],tau(*,ii_band2),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band21=plot(pepi_t[0,*,ii_band2],tau(*,ii_band2),$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band22=plot(pepi_t[1,*,ii_band2],tau(*,ii_band2),$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band23=plot(pepi_t[2,*,ii_band2],tau(*,ii_band2),$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg2 = LEGEND(TARGET=[p_band21,p_band22,p_band23], POSITION=[0.95,0.95],/normal)


;#3 155-224 A



p_band31=plot(pepi[0,*,ii_band3],tau(*,ii_band3),name=maj[0],$
  xtitle=xt+wv3_x,ytitle=yt_tau,$
  xlog=1,xstyle=1,ystyle=1,xrange=xr,yrange=yr_tau3,$
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,3])



p_band32=plot(pepi[1,*,ii_band3],tau(*,ii_band3),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band33=plot(pepi[2,*,ii_band3],tau(*,ii_band3),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band31=plot(pepi_t[0,*,ii_band3],tau(*,ii_band3),$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band32=plot(pepi_t[1,*,ii_band3],tau(*,ii_band3),$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band33=plot(pepi_t[2,*,ii_band3],tau(*,ii_band3),$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg3 = LEGEND(TARGET=[p_band31,p_band32,p_band33], POSITION=[0.45,0.45],/normal)



;#2 290- 320 A



p_band41=plot(pepi[0,*,ii_band4],tau(*,ii_band4),name=maj[0],$
  xtitle=xt+wv4_x,ytitle=yt_tau,$
  xlog=1,xstyle=1,ystyle=1,xrange=[0.1,100.],yrange=[max(tau(*,ii_band4)),min(tau(*,ii_band4))],ylog=1,$;
  linestyle=lstyl[0],thick=thck,/current,layout=[2,2,4])



p_band42=plot(pepi[1,*,ii_band4],tau(*,ii_band4),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band43=plot(pepi[2,*,ii_band4],tau(*,ii_band4),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band41=plot(pepi_t[0,*,ii_band4],tau(*,ii_band4),$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band42=plot(pepi_t[1,*,ii_band4],tau(*,ii_band4),$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band43=plot(pepi_t[2,*,ii_band4],tau(*,ii_band4),$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg4 = LEGEND(TARGET=[p_band41,p_band42,p_band43], POSITION=[0.95,0.45],/normal)

;..............................................................................
;..............................................................................

;#2 290- 320 A

w_15 = WINDOW(WINDOW_TITLE=tilt, $
  DIMENSIONS=[800,800])


p_band51=plot(eiionz_run[0,*,ii_band4],zz(*,ii_band4),name=maj[0],$
  xtitle='Electron impact ionization'+wv4_x+'!C Black- local, Red-transport',ytitle=yt_alt,$
  xstyle=1,ystyle=1,yrange=[40.,360.],$;
  linestyle=lstyl[0],thick=thck,/current)



p_band52=plot(eiionz_run[1,*,ii_band4],zz(*,ii_band4),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band53=plot(eiionz_run[2,*,ii_band4],zz(*,ii_band4),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


;transport pepi
pt_band511=plot(eiionzt_run[0,*,ii_band4],zz(*,ii_band4),$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

pt_band522=plot(eiionzt_run[1,*,ii_band4],zz(*,ii_band4),$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

pt_band533=plot(eiionzt_run[2,*,ii_band4],zz(*,ii_band4),$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg4 = LEGEND(TARGET=[p_band51,p_band52,p_band53], POSITION=[0.95,0.45],/normal)







w_16 = WINDOW(WINDOW_TITLE=tilt, $
  DIMENSIONS=[800,800])


p_band61=plot(photoi_run[0,*,ii_band4],zz(*,ii_band4),name=maj[0],$
  xtitle='Photoionization'+wv4_x,ytitle=yt_alt,$
  xstyle=1,ystyle=1,yrange=[40.,360.],$;
  linestyle=lstyl[0],thick=thck,/current)



p_band62=plot(photoi_run[1,*,ii_band4],zz(*,ii_band4),name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

p_band63=plot(photoi_run[2,*,ii_band4],zz(*,ii_band4),name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


  leg4 = LEGEND(TARGET=[p_band61,p_band62,p_band63], POSITION=[0.95,0.45],/normal)


  w_17 = WINDOW(WINDOW_TITLE=tilt, $
    DIMENSIONS=[800,800])


  p_band71=plot(eiionzt_run[0,*,ii_band4],zz(*,ii_band4),name=maj[0],$
    xtitle='Transport eiionz'+wv4_x,ytitle=yt_alt,$
    xstyle=1,ystyle=1,yrange=[40.,360.],$;
    linestyle=lstyl[0],thick=thck,/current,color='red')



  p_band72=plot(eiionzt_run[1,*,ii_band4],zz(*,ii_band4),name=maj[1],$
    linestyle=lstyl[1],thick=thck,/overplot,color='red')

  p_band73=plot(eiionzt_run[2,*,ii_band4],zz(*,ii_band4),name=maj[2],$
    linestyle=lstyl[2],thick=thck,/overplot,color='red')


  leg7= LEGEND(TARGET=[p_band71,p_band72,p_band73], POSITION=[0.95,0.45],/normal)

close,/all
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------

;Comparisons of Solomon &Qian 2005 and NRL with all solar bins included 


fname='Solomon&Qian2005_allrun' ;change this each time
sav_loc="C:\Users\srimo\Desktop\agu_plots\"
filename1=loc_sav+fname+'.sav'
;file has:
;eiionz,eiionz_transp,eiionz_wv,eiionz_transp_wv,photoi,tflux,pespec,ener,mpf,zz
restore, filename1


eiionz1=eiionz
eiionz_transp1=eiionz_transp
eiionz_wv1=eiionz_wv
eiionz_transp_wv1=eiionz_transp_wv
photoi1=photoi
tflux1=tflux
pespec1=pespec
mfp1=mfp
wv11=wv1
wv21=wv2
zcol1=zcol

fname='X28_allrun'
filename2=loc_sav+fname+'.sav'
restore,filename2

eiionz2=eiionz
eiionz_transp2=eiionz_transp
eiionz_wv2=eiionz_wv
eiionz_transp_wv2=eiionz_transp_wv
photoi2=photoi
tflux2=tflux
pespec2=pespec
mfp2=mfp
wv12=wv1
wv22=wv2
zcol2=zcol

;Comparisons of 2 two different spectral binnings


;.................Plots.................................................


;Plot variables...............

colors=['red','dark green','dodger blue']
lstyl=[2,1,0]
pos1=[0.1,0.5,0.45,0.95]
pos2=[0.1,0.1,0.45,0.45]
pos3=[0.47,0.5,0.95,0.95]
pos4=[0.47,0.1,0.95,0.45]
thck=5
leg_loc=[0.90, 0.90]




;Electron Impact Ionisation plot


yr=[min(zz),max(zz)]
leg_loc=[0.4,0.9]
xt='Electron impact ionization !C Black- S&Q05, Red-NRL X28'
yt='Altitude (km)'
tilt='Local Electron impact ionization'
maj=['O','O2','N2','Total']

;local

w1 = WINDOW(WINDOW_TITLE=tilt, $
  DIMENSIONS=[800,800])


ei1=plot(eiionz1[0,*],zz,name=maj[0],$
  xtitle=xt,ytitle=yt,title=tilt,$
  xstyle=1,ystyle=1,yrange=yr,xlog=1,$;
  linestyle=lstyl[0],thick=thck,/current)

ei2=plot(eiionz1[1,*],zz,name=maj[1],$
  linestyle=lstyl[1],thick=thck,/overplot)

ei3=plot(eiionz1[2,*],zz,name=maj[2],$
  linestyle=lstyl[2],thick=thck,/overplot)


ei4=plot(eiionz2[0,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

ei5=plot(eiionz2[1,*],zz,$
  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)

ei6=plot(eiionz2[2,*],zz,$
  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)


leg1= LEGEND(TARGET=[ei1,ei2,ei3], POSITION=leg_loc,/normal)

; transport---------------------

yr=[min(zz),max(zz)]
xr=[1.e-5,1.e5]
xt='Electron impact ionization (Pe) !C!C (cm!U-3!N s!U-1!N)'
yt='Altitude'
maj=['O','O2','N2']
leg_loc=[0.92,0.85]
wtilt='Transport Electron impact ionization'

w21 = WINDOW(WINDOW_TITLE=wtilt, $
  DIMENSIONS=[800,800])


eit1=plot(eiionz_transp1[0,*],zz,name=maj[0],$
  xtitle='O '+xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  color=colors[0],thick=thck,/current)
eit2=plot(eiionz_transp2[0,*],zz,$
    linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

w21.Save, sav_loc+"comp_O_eiionz.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

w22 = WINDOW(WINDOW_TITLE=wtilt, $
    DIMENSIONS=[800,800])


eit3=plot(eiionz_transp1[1,*],zz,name=maj[1],$
    xtitle='O!L2!N '+xt,ytitle=yt,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    color=colors[1],thick=thck,/current)
eit4=plot(eiionz_transp2[1,*],zz,name=maj[1],$
     linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)

w22.Save, sav_loc+"comp_O2_eiionz.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

w23 = WINDOW(WINDOW_TITLE=wtilt, $
    DIMENSIONS=[800,800])
eit5=plot(eiionz_transp1[2,*],zz,name=maj[2],$
    xtitle='N!L2!N '+xt,ytitle=yt,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    color=colors[2],thick=thck,/current)
eit6=plot(eiionz_transp2[2,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
w23.Save, sav_loc+"comp_N2_eiionz.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;leg2= LEGEND(TARGET=[eit1,eit2,eit3], POSITION=leg_loc,/normal)
 
; ###############################################################################
 w211= WINDOW(WINDOW_TITLE=wtilt, $
  DIMENSIONS=[800,800])


eit1=plot(eiionz_transp1[0,*],zz,name=maj[0],$
  xtitle=xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  color=colors[0],thick=thck,/current)
eit2=plot(eiionz_transp2[0,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)



eit3=plot(eiionz_transp1[1,*],zz,name=maj[1],$
 color=colors[1],thick=thck,/overplot)
eit4=plot(eiionz_transp2[1,*],zz,$
  linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)


eit5=plot(eiionz_transp1[2,*],zz,name=maj[2],$
  color=colors[2],thick=thck,/overplot)
eit6=plot(eiionz_transp2[2,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
  
  leg211= LEGEND(TARGET=[eit1,eit3,eit5], POSITION=leg_loc,/normal)
w211.Save, sav_loc+"comp_eiionz.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT


 
 
; ################################################################################
 
;photoionisation----------------------------


xr=[1.e-5,1.e5]
yr=[min(zz),max(zz)]
xt='Photoionisation (Pi)(cm!U-3!N s!U-1!N) '
yt='Altitude (km)'
tilt= 'Photoi'



w31 = WINDOW(WINDOW_TITLE=wtilt, $
  DIMENSIONS=[800,800])


p1=plot(photoi1[0,*],zz,name=maj[0],$
  xtitle='O '+xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  color=colors[0],thick=thck,/current)
p2=plot(photoi2[0,*],zz,$
    linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

w31.Save, sav_loc+"comp_O_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

w32 = WINDOW(WINDOW_TITLE=wtilt, $
    DIMENSIONS=[800,800])


p3=plot(photoi1[1,*],zz,name=maj[1],$
    xtitle='O!L2!N '+xt,ytitle=yt,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    color=colors[1],thick=thck,/current)
p4=plot(photoi2[1,*],zz,name=maj[1],$
     linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)

w32.Save, sav_loc+"comp_O2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

w33 = WINDOW(WINDOW_TITLE=wtilt, $
    DIMENSIONS=[800,800])
p5=plot(photoi1[2,*],zz,name=maj[2],$
    xtitle='N!L2!N '+xt,ytitle=yt,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
    color=colors[2],thick=thck,/current)
p6=plot(photoi2[2,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
w33.Save, sav_loc+"comp_N2_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

;#######################################################################################

w311 = WINDOW(WINDOW_TITLE=wtilt, $
  DIMENSIONS=[800,800])


p1=plot(photoi1[0,*],zz,name=maj[0],$
  xtitle=xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
  color=colors[0],thick=thck,/current)
  
p2=plot(photoi2[0,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)

p3=plot(photoi1[1,*],zz,name=maj[1],$
   color=colors[1],thick=thck,/overplot)

p4=plot(photoi2[1,*],zz,$
  linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)


p5=plot(photoi1[2,*],zz,name=maj[2],$
    color=colors[2],thick=thck,/overplot)
    
p6=plot(photoi2[2,*],zz,$
  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
  
  l311=legend(target=[p1,p3,p5],position=[0.93,0.85],/normal)
w311.Save, sav_loc+"comp_photoi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT




restore,'C:\Users\srimo\Desktop\nrl_files\sav_files\minuslymanbeta.sav'
mphotoi=photoi

w411 = WINDOW(WINDOW_TITLE=wtilt, $
  DIMENSIONS=[800,800])


mp1=plot(photoi2[1,*],zz,name='With Lyman beta',$
  xtitle='O!L2!N'+xt,ytitle=yt,$
  xstyle=1,ystyle=1,yrange=[80,max(zz)],xrange=[0.1,1.e4],xlog=1,$;
  linestyle=2,color=colors[1],thick=thck,/current)

mp2=plot(mphotoi[1,*],zz,name='Without Lyman beta',$
  thick=thck,color=colors[1],/overplot)
  
  lmp=legend(target=[mp1,mp2],position=[0.9,0.85],/normal)
w411.Save, sav_loc+"mlbeta.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT





;;######################################################################################
;
;;.............................................................................
;;Ratio of tflux(transport) and pespec(local) wrt altitude
;
;dif1=tflux1/pespec1
;dif2=tflux2/pespec2
;
;xr=[1.e-3,5.e0]
;yr=[min(zz),max(zz)]
;xt='trasport flux/local flux !C black- S&Q05, red- NRL X28'
;yt='Altitude'
;tilt='tflux/pespec vs alt'
;
;w4= WINDOW(WINDOW_TITLE=tilt, $
;     DIMENSIONS=[800,800])
;
;ener_pos=[29,47,170]
;
;rz1=plot(dif1[*,ener_pos[0]],zz,name=strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;  xtitle=xt,ytitle=yt,$
;  xstyle=1,ystyle=1,yrange=yr_alt,xrange=xr,xlog=1,$;
;  linestyle=lstyl[0],thick=thck,/current)
;
;rz2=plot(dif1[*,ener_pos[1]],zz,name=strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;    linestyle=lstyl[1],thick=thck,/overplot)
;
;rz3=plot(dif1[*,ener_pos[2]],zz,name=strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;    linestyle=lstyl[2],thick=thck,/overplot)
;
;
;rz4=plot(dif2[*,ener_pos[0]],zz,name=strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;    linestyle=lstyl[1],color=colors[0],thick=thck,/overplot)
;
;
;rz5=plot(dif2[*,ener_pos[1]],zz,name=strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;    linestyle=lstyl[1],color=colors[0],thick=thck,/overplot)
;
;rz6=plot(dif2[*,ener_pos[2]],zz,name=strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;    linestyle=lstyl[2],color=colors[0],thick=thck,/overplot)
;
;leg4= LEGEND(TARGET=[rz1,rz2,rz3], POSITION=leg_loc,/normal)
;
;
;
;
;;Ratio of tflux(transport) and pespec(local) wrt energy
;
;yr=[0.0,5.]
;xr=[10.,max(ener)]
;leg_loc=[0.9,0.9]
;yt='trasport flux/local flux'
;xt='energy (eV)!C black-S&Q05 red- NRL X28'
;tilt='tflux/pespec vs energy'
;
;alt_pos=[10,47,100]
;
;
;w5= WINDOW(WINDOW_TITLE=tilt, $
;  DIMENSIONS=[800,800])
;
;
;re1=plot(ener,dif1[alt_pos[0],*],name=strtrim(string(zz[alt_pos[0]]),1)+'km',$
;  xtitle=xt,ytitle=yt,$
;  xstyle=1,ystyle=1,xrange=xr,xlog=1,yrange=yr,$;
;  linestyl=lstyl[0],thick=thck,/current)
;
;
;re2=plot(ener,dif1[alt_pos[1],*],name=strtrim(string(zz[alt_pos[1]]),1)+'km',$
;  linestyl=lstyl[1],thick=thck,/overplot)
;
;re3=plot(ener,dif1[alt_pos[2],*],name=strtrim(string(zz[alt_pos[2]]),1)+'km',$
;  linestyl=lstyl[2],thick=thck,/overplot)
;  
;  
;re4=plot(ener,dif2[alt_pos[0],*],name=strtrim(string(zz[alt_pos[0]]),1)+'km',$
;    linestyl=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;re5=plot(ener,dif2[alt_pos[1],*],name=strtrim(string(zz[alt_pos[1]]),1)+'km',$
;    linestyl=lstyl[1],thick=thck,color=colors[0],/overplot)
;
;re6=plot(ener,dif2[alt_pos[2],*],name=strtrim(string(zz[alt_pos[2]]),1)+'km',$
;    linestyl=lstyl[2],thick=thck,color=colors[0],/overplot)  
;
;leg5 = LEGEND(TARGET=[re1,re2,re3], POSITION=leg_loc,/normal)
;
;
;
;;...................................................................................
;;Pe/Pi ratio-----------------------------------------------------
;
;xr=[0.01,1000.]
;yr=[min(zz),150.]
;leg_loc=[0.9,0.85]
;xt='Pe/Pi '
;yt='Altitude (km)'
;tilt='Local Pe/Pi'
;
;
;;local
;w6 = WINDOW(WINDOW_TITLE=tilt, $
;  DIMENSIONS=[800,800])
;
;
;pp1=plot(eiionz1[0,*]/photoi1[0,*],zz,name=maj[0],$
;  xtitle=xt,ytitle=yt,title=tilt,$
;  xstyle=1,ystyle=1,yrange=yr,xlog=1,xrange=xr,$;
;  linestyle=lstyl[0],thick=thck,/current)
;
;
;
;pp2=plot(eiionz1[1,*]/photoi1[1,*],zz,name=maj[1],$
;  linestyle=lstyl[1],thick=thck,/overplot)
;
;pp3=plot(eiionz1[2,*]/photoi1[2,*],zz,name=maj[2],$
;  linestyle=lstyl[2],thick=thck,/overplot)
;
;
;
;pp4=plot(eiionz2[0,*]/photoi2[0,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;pp5=plot(eiionz2[1,*]/photoi2[1,*],zz,$
;  linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)
;
;pp6=plot(eiionz2[2,*]/photoi2[2,*],zz,$
;  linestyle=lstyl[2],thick=1.,color=colors[0],/overplot)
;
;
;leg6 = LEGEND(TARGET=[pp1,pp2,pp3], POSITION=leg_loc,/normal)
;
;
;
;;Transport
;tilt='Transport Pe/Pi'
;
;
;w71 = WINDOW(WINDOW_TITLE=wtilt, $
;  DIMENSIONS=[800,800])
;
;
;pp1=plot(eiionz_transp1[0,*]/photoi1[0,*],zz,name=maj[0],$
;  xtitle='O '+xt,ytitle=yt,$
;  xstyle=1,ystyle=1,yrange=[70.,150.],xrange=xr,xlog=1,$;
;  color=colors[0],thick=thck,/current)
;pp2=plot(eiionz_transp2[0,*]/photoi2[0,*],zz,$
;    linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;
; w71.Save, sav_loc+"comp_O_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
; 
;w72 = WINDOW(WINDOW_TITLE=wtilt, $
;    DIMENSIONS=[800,800])
;
;
;pp3=plot(eiionz_transp1[1,*]/photoi1[1,*],zz,name=maj[1],$
;    xtitle='O!L2!N '+xt,ytitle=yt,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    color=colors[1],thick=thck,/current)
;pp4=plot(eiionz_transp2[1,*]/photoi2[1,*],zz,name=maj[1],$
;     linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)
;
;w72.Save, sav_loc+"comp_O2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;w73 = WINDOW(WINDOW_TITLE=wtilt, $
;    DIMENSIONS=[800,800])
;pp5=plot(eiionz_transp1[2,*]/photoi1[2,*],zz,name=maj[2],$
;    xtitle='N!L2!N '+xt,ytitle=yt,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    color=colors[2],thick=thck,/current)
;pp6=plot(eiionz_transp2[2,*]/photoi2[2,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
;w73.Save, sav_loc+"comp_N2_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;;############################################################################################
;xr=[0.01,1000.]
;yr=[min(zz),max(zz)]
;leg_loc=[0.9,0.85]
;xt='Pe/Pi '
;yt='Altitude (km)'
;tilt='Pe/Pi'
;
;
;
;w711 = WINDOW(WINDOW_TITLE=tilt, $
;  DIMENSIONS=[800,800])
;
;
;pp1=plot(eiionz_transp1[0,*]/photoi1[0,*],zz,name=maj[0],$
;  xtitle=xt,ytitle=yt,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  color=colors[0],thick=thck,/current)
;
;pp2=plot(eiionz_transp2[0,*]/photoi2[0,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;
;pp3=plot(eiionz_transp1[1,*]/photoi1[1,*],zz,name=maj[1],$
;      color=colors[1],thick=thck,/overplot)
;
;pp4=plot(eiionz_transp2[1,*]/photoi2[1,*],zz,$
;  linestyle=lstyl[0],color=colors[1],thick=thck,/overplot)
;
;
;pp5=plot(eiionz_transp1[2,*]/photoi1[2,*],zz,name=maj[2],$
;  color=colors[2],thick=thck,/overplot)
;  
;pp6=plot(eiionz_transp2[2,*]/photoi2[2,*],zz,$
;  linestyle=lstyl[0],thick=thck,color=colors[2],/overplot)
;  
;  l711=legend(target=[pp1,pp3,pp5],position=[0.9,0.85],/normal)
;  
;w711.Save, sav_loc+"comp_pepi.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;
;
;
;
;
;
;
;
;
;
;
;
;;########################################################################################################
;;Electron impact ionisation flux for one altitude, for both local and transport
;xr=[0.1,1.e5]
;xt='energy (eV) for !C black- SQ05 , red- NRL X28'
;yt='Electron impact ionisation flux'
;tilt='tflux & pespec'
;leg_loc1=[0.1,0.9]
;leg_loc2=[0.55,0.9]
;leg_loc3=[0.95,0.9]
;
;
;w6 = WINDOW(WINDOW_TITLE=tilt, $
;  DIMENSIONS=[800,800])
;
;eif1=plot(ener,tflux1[alt_pos[0],*],name='transport',$
;  xrange=xr,histogram=1,xlog=1,ylog=1,$
;  xtitle=xt+strtrim(string(zz[alt_pos[0]]),1)+'km',ytitle=yt,$
;  linestyle=lstyl[0],/current,thick=thck,layout=[3,1,1])
;
;eif2=plot(ener,pespec1[alt_pos[0],*],name='local',$
;  histogram=1,linestyle=lstyl[1],thick=thck,/overplot)
;  
;eif3=plot(ener,tflux2[alt_pos[0],*],name='transport',$
;    histogram=1,linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;eif4=plot(ener,pespec2[alt_pos[0],*],name='local',$
;    histogram=1,linestyle=lstyl[1],thick=thck,color=colors[1],/overplot)
;
;leg6 = LEGEND(TARGET=[eif1,eif2], POSITION=leg_loc1,/normal)
;
;
;eif5=plot(ener,tflux1[alt_pos[1],*],name='transport',$
;  xrange=xr,histogram=1,xlog=1,ylog=1,$
;  xtitle=xt+strtrim(string(zz[alt_pos[1]]),1)+'km',ytitle=yt,$
;  linestyle=lstyl[0],/current,thick=thck,layout=[3,1,2])
;
;eif6=plot(ener,pespec1[alt_pos[1],*],name='local',$
;  histogram=1,linestyle=lstyl[1],thick=thck,/overplot)
;
;eif7=plot(ener,tflux2[alt_pos[1],*],name='transport',$
;    histogram=1,linestyle=lstyl[0],thick=thck,color=colors[0],/overplot)
;
;eif8=plot(ener,pespec2[alt_pos[1],*],name='local',$
;    histogram=1,linestyle=lstyl[1],thick=thck,color=colors[0],/overplot)
;
;leg7 = LEGEND(TARGET=[eif5,eif6], POSITION=leg_loc2,/normal)
;
;
;eif9=plot(ener,tflux1[alt_pos[2],*],name='transport',$
;  xrange=xr,histogram=1,xlog=1,ylog=1,$
;  xtitle=xt+strtrim(string(zz[alt_pos[2]]),1)+'km',ytitle=yt,$
;  linestyle=lstyl[0],thick=thck,/current,layout=[3,1,3])
;
;eif10=plot(ener,pespec1[alt_pos[2],*],name='local',$
;  histogram=1,linestyle=lstyl[0],thick=thck,/overplot)
;
;eif11=plot(ener,tflux2[alt_pos[2],*],name='transport',$
;    histogram=1,linestyle=lstyl[0],color=colors[0],thick=thck,/overplot)
;
;  eif12=plot(ener,pespec[alt_pos[2],*],name='local',$
;    histogram=1,linestyle=lstyl[1],color=colors[0],thick=thck,/overplot)
;
;leg8 = LEGEND(TARGET=[eif9,eif10], POSITION=leg_loc3,/normal)
;
;
;; Column density of N2
;
;xr=[0.1,1.e5]
;xt='Column Density of N2 !C black- SQ05 , red- NRL X28'
;yt='Altitude (km)'
;tilt='col dens' 
;leg_loc=[0.9,0.45]
;
;w7=Window(window_title=tilt,$
;         dimension=[800,800])
;         
;z1=plot(zcol1(2,*),zz,$
;   xtitle=xt,ytitle=yt,$
;   xlog=1,xstyle=1,ystyle=1,$
;   /current)
;
;z2=plot(zcol2(2,*),zz,$
;    color=colors[0],linestyle=lstyl[0],/overplot)
;    
;
;;Electron impact ionsation of N2
;tilt='eeionz of N2 @ 1e4 eV'
;leg_loc1=[0.4,0.9]
;leg_loc2=[0.4,0.45]
;ener_pos=(where(ener ge 1.e4 and ener le 2.e4))[2]
;xt='Electron Impact ionisation of N2 @'+strtrim(string(ener[ener_pos]))+'eV'
;
;
;w8=Window(window_title=tilt,$
;  dimension=[800,800])
;
;ein21=plot(eiionz_wv1(*,2,ener_pos),zz,name='SQ05 local',$
;           xtitle=xt, ytitle=yt,$
;          xlog=1,xstyle=1,ystyle=1,$
;         /current,layout=[1,2,1])
;
;ein22=plot(eiionz_wv2(*,2,ener_pos),zz,name='NRL local',$
;           color=colors[0],/overplot)
;
;
;leg9=legend(target=[ein21,ein22],position=leg_loc2,/normal)
;
;ein23=plot(eiionz_transp_wv1(*,2,ener_pos),zz,name='SQ05 transp',$
;            xtitle=xt, ytitle=yt,$
;            xlog=1,xstyle=1,ystyle=1,$
;            /current,layout=[1,2,2])
;
;
;ein24=plot(eiionz_transp_wv2(*,2,ener_pos),zz,name='NRL transp',$
;            color=colors[0],/overplot)
;
;
;leg10 = LEGEND(TARGET=[ein23,ein24], POSITION=leg_loc1,/normal)
;
;;Mean free path calculation
;
;xr=[1.,1.e7]
;yr=[min(zz),max(zz)]
;xt='Mean free path (km)'
;yt='Altitude (km)'
;tilt='Black- SQ05 red- NRL !C'
;wtilt='Mean free path'
;leg_loc1=[0.1,0.9]
;leg_loc2=[0.55,0.9]
;leg_loc3=[0.95,0.9]
;ener_pos=[29,47,170]
;
;w9=window(window_title=wtilt,$
;           dimension=[800,800])
;
;m1=plot(mfp1[0,*,ener_pos[0]],zz,$
;  xtitle=maj[0]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,1])
;
;m2=plot(mfp2[0,*,ener_pos[0]],zz,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;    
;m3=plot(mfp1[1,*,ener_pos[0]],zz,$
;    xtitle=maj[1]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;    yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;    /current,layout=[2,2,2])
;
;m4=plot(mfp2[1,*,ener_pos[0]],zz,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;m5=plot(mfp1[2,*,ener_pos[0]],zz,$
;    xtitle=maj[2]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;    yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;    /current,layout=[2,2,3])
;
;m6=plot(mfp2[2,*,ener_pos[0]],zz,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;    
;;m7=plot(mfp1[3,*,ener_pos[0]],zz,$
;;    xtitle=maj[3]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[0]]),1)+'eV',$
;;    yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;;    /current,layout=[2,2,4])
;;
;;m8=plot(mfp2[3,*,ener_pos[0]],zz,$
;;    thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;;    
;    
;;----------------------------------------------------------------------------------------
;w10=window(window_title=wtilt,$
;   dimension=[800,800])
;
;m9=plot(mfp1[0,*,ener_pos[1]],zz,$
;  xtitle=maj[0]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,1])
;
;m10=plot(mfp2[0,*,ener_pos[1]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;m11=plot(mfp1[1,*,ener_pos[1]],zz,$
;  xtitle=maj[1]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,2])
;
;m12=plot(mfp2[1,*,ener_pos[1]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;m13=plot(mfp1[2,*,ener_pos[1]],zz,$
;  xtitle=maj[2]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,3])
;
;m14=plot(mfp2[2,*,ener_pos[1]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;;m15=plot(mfp1[3,*,ener_pos[1]],zz,$
;;  xtitle=maj[3]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[1]]),1)+'eV',$
;;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;;  /current,layout=[2,2,4])
;;
;;m16=plot(mfp2[3,*,ener_pos[1]],zz,$
;;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
; 
;;--------------------------------------------------------------------------------------
;
;w11=window(window_title=wilt,$
;  dimension=[800,800])
;
;m17=plot(mfp1[0,*,ener_pos[2]],zz,$
;  xtitle=maj[0]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,1])
;
;m18=plot(mfp2[0,*,ener_pos[2]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;m19=plot(mfp1[1,*,ener_pos[2]],zz,$
;  xtitle=maj[1]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,2])
;
;m20=plot(mfp2[1,*,ener_pos[2]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;m21=plot(mfp1[2,*,ener_pos[2]],zz,$
;  xtitle=maj[2]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;  /current,layout=[2,2,3])
;
;m22=plot(mfp2[2,*,ener_pos[2]],zz,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;
;;m23=plot(mfp1[3,*,ener_pos[2]],zz,$
;;  xtitle=maj[3]+xt,ytitle=yt,title=tilt+strtrim(string(ener[ener_pos[2]]),1)+'eV',$
;;  yrange=yr,xlog=1,xstyle=1,ystyle=1,$
;;  /current,layout=[2,2,4])
;;
;;m24=plot(mfp2[3,*,ener_pos[2]],zz,$
;;  thick=thck,color=colors[0],linestyle=lstyl[0],/overplot)
;;
;
;
;    
;;    ----------------------------------------------------
;;-----------------------------------------------------------
;;DELETE THIS PART LATER
;;#2 NRl Spectrum
;
;filename='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\x28_custombins_myver.sav'
;restore, filename
;nrlwv1=wavelength_low
;nrlwv2=wavelength_high
;nrlssflux=spectrum[90,*] ;for x_28custombins_myver.sav file only
;
;
;; #3 Solomon and Qian 2005
;filename='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\Solomon_Qian2005.sav'
;;file has:
;;A_fac, fref, wavelength_low, wavelength_high
;restore, filename
;F107=70.
;F107a=70.
;P107 = (F107+F107A)/2.
;n_wvl=n_elements(a_fac)
;ssflux= fltarr(n_wvl)
;for L=0,n_wvl-1 do begin
;          SSFLUX(L) = fref(L) * (1. + A_fac(L)*(P107-80.))
;          IF (SSFLUX(L) LT 0.8*fref(L)) then SSFLUX(L) = 0.8*fref(L)
;endfor
;
;euvwv1=wavelength_low
;euvwv2=wavelength_high
;euvssflux=ssflux
;
;
;ww1=window(window_title='Solar Flux',$
;  dimension=[800,800])
;
;s1=plot(euvwv1, euvssflux, $
;  ylog=1,ystyle=1,xlog=1,xstyle=1,$
;  xrange=[0.5,1.75e3],yrange=[1.e1,1.e12],$
;  thick=thck,histogram=1,color='blue',/current,$ ;, position=[0.18,0.1,0.95,0.85]
;  xtitle='Wavlength ('+'!3' + STRING(197B) + '!X'+')' ,$
;  ytitle='Solar Flux (Photons cm!U-2!N s!U-1!N)',$
;  name='EUVAC reference flux')
;
;s2=plot(nrlwv1,nrlssflux,$
;         thick=thck,histogram=1,color='red',$
;         name='NRL flux',/overplot)
;
;ll1=legend(target=[s1,s2],position=[0.55,0.85],/normal)
;         
;ww1.Save, sav_loc+"comp_flux.jpeg", BORDER=20, RESOLUTION=600, /TRANSPARENT
;
;;lyman beta 
;save, photoi,eiionz,eiionz_transp,zz,ssflux,wv1,wv2,$
;  filename='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\'+fname+'.sav'
;  
;  
; lw= window(window_title='Lyman beta ionisation',dimension=[800,800])
;  
;  lb1=plot(photoi[1,*],zz ,$
;    ystyle=1,xlog=1,xstyle=1,$
;    xrange=[1.e-2,1.e4],yrange=[min(zz),max(zz)],$
;    thick=thck,color=colors[1],/current,$ ;, position=[0.18,0.1,0.95,0.85]
;    xtitle='Photoinisation of O!L2!N due to Lyman beta !C!C (cm!U-3!N s!U-1!N)' ,$
;    ytitle='Altitude (km)')
;
;  lw.Save, sav_loc+"lymanbeta_photO2.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;  
;  
;  close,/all
;  
;  restore,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_334.sav'
;  thck=5
; 
;  colors=['deep pink','magenta','medium purple','medium orchid','violet',$
;           'dark blue','indigo','purple','blue violet','dodger blue','royal blue']
;;  p334=intarr(11)
;  ww10=window(window_title="Spectra",dimension=[800,800])
;  
;  p334= plot(wv1,ssflux_fi(0,*),$
;               ylog=1,ystyle=1,xlog=1,xstyle=1,$
;               xrange=[0.4,1.75e3],yrange=[1.e1,1.e12],$
;               thick=thck,histogram=1,color=colors[0],/current,$ ;, position=[0.18,0.1,0.95,0.85]
;               xtitle='Wavlength ('+'!3' + STRING(197B) + '!X'+')' ,$
;               ytitle='Solar Flux (Photons cm!U-2!N s!U-1!N)',$
;               name='X28 spectrum 1')
;               
;  for i=1,10 do begin
;    
;    
;    p334= plot(wv1,ssflux_fi(i,*),name='X 28 spectrum'+ strtrim(string(i+1),1),$
;           thick=thck,histogram=1,color=colors[i],/overplot)
;    
;  endfor
;;  ll10=legend(target=p334,position=[0.45,0.85],/normal)
;  
;  ww10.save,sav_loc+"solarflux_334.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;  
;;  Opepi334=intarr(11)
;  ww11=window(window_title="O PEPI",dimension=[800,800])
;
;  Opepi334= plot(pepi(0,0,*),zz,$
;    ystyle=1,xlog=1,xstyle=1,$
;    xrange=[0.01,1.0e3],yrange=[min(zz),max(zz)],$
;    thick=thck,color=colors[0],/current,$ ;, position=[0.18,0.1,0.95,0.85]
;    xtitle='Pe/Pi for O' ,$
;    ytitle='Altitude (km)')
;
;  for i=1,10 do begin
;
;
;    Opepi334= plot(pepi(0,i,*),zz,$
;      thick=thck,color=colors[i],/overplot)
;
;  endfor
;;  ll11=legend(target=Opepi334,position=[0.45,0.85],/normal)
;  ww11.save,sav_loc+"Opepi_334.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;  
;;  O2pepi334=intarr(11)
;  ww12=window(window_title="O2 PEPI",dimension=[800,800])
;
;  O2pepi334= plot(pepi(1,0,*),zz,$
;    ystyle=1,xlog=1,xstyle=1,$
;    xrange=[0.01,1.0e3],yrange=[min(zz),max(zz)],$
;    thick=thck,color=colors[0],/current,$ ;, position=[0.18,0.1,0.95,0.85]
;    xtitle='Pe/Pi for O!L2!N' ,$
;    ytitle='Altitude (km)' )
;
;  for i=1,10 do begin
;
;
;    O2pepi334= plot(pepi(1,i,*),zz,$
;      thick=thck,color=colors[i],/overplot)
;
;  endfor
;;  ll12=legend(target=O2pepi334,position=[0.45,0.85],/normal)
;
;ww12.save,sav_loc+"O2pepi_334.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;
;;  N2pepi334=intarr(11)
;  ww13=window(window_title="N2 PEPI",dimension=[800,800])
;
;  N2pepi334= plot(pepi(2,0,*),zz,$
;    ystyle=1,xlog=1,xstyle=1,$
;    xrange=[0.01,1.0e3],yrange=[min(zz),max(zz)],$
;    thick=thck,color=colors[0],/current,$ ;, position=[0.18,0.1,0.95,0.85]
;    xtitle='Pe/Pi for N!L2!N' ,$
;    ytitle='Altitude (km)')
;  for i=1,10 do begin
;
;
;    N2pepi334= plot(pepi(2,i,*),zz,$
;      thick=thck,color=colors[i],/overplot)
;
;  endfor
;;  ll12=legend(target=N2pepi334,position=[0.45,0.85],/normal)
;
;
;ww13.save,sav_loc+"N2pepi_334.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT
;restore,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\lymanbeta_nrl_334.sav'
;
;wlb=window(window_title="O2 Lyman beta photoionisation",dimension=[800,800])
;
;lb= plot(photoi_lbeta(1,0,*),zz,$
;  ystyle=1,xlog=1,xstyle=1,$
;  xrange=[0.01,1.0e4],yrange=[min(zz),max(zz)],$
;  thick=thck,color=colors[0],linestyle=4,/current,$ ;, position=[0.18,0.1,0.95,0.85]
;  xtitle='Photoionisation of O!L2!N due to Lyman beta !C!C (cm!U-3!N s!-1!N)' ,$
;  ytitle='Altitude (km)')
;for i=1,10 do begin
;
;
;  lb= plot(photoi_lbeta(1,i,*),zz,$
;    thick=thck,color=colors[i],/overplot)
;
;endfor
;;  ll12=legend(target=N2pepi334,position=[0.45,0.85],/normal)


close,/all

end
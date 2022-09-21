;pro average_ener,N2_exctrate1,ener,avg_ener_3,avg_ener_2,avg_ener_1,N2_rate_1
;  nei=(size(N2_exctrate1))[1]
;  jmax=(size(N2_exctrate1))[2]
;  nbins=(size(N2_exctrate1))[3]
;  lmax=19 ;(size(N2_exctrate1))[4]
;
;
;
;  N2_tot_3 = fltarr(nei,jmax,lmax)
;  N2_tot_2 = fltarr(nei,jmax)
;  N2_tot_1 = fltarr(jmax)
;
;  N2_rate_3 = fltarr(nei,jmax,lmax)
;  N2_rate_2 = fltarr(nei,jmax)
;  N2_rate_1 = fltarr(jmax)
;
;  avg_ener_3 = fltarr(nei,jmax,lmax)
;  avg_ener_2 = fltarr(nei,jmax)
;  avg_ener_1 = fltarr(jmax)
;
;
;
;
;  for l=0,lmax-1 do begin  ;wavelength bins
;
;    for e=0,nbins-1 do begin
;      for j=0,jmax-1 do begin  ;altitude
;
;        for s=0,nei-1 do begin  ;states
;
;          N2_tot_3[s,j,l]=N2_tot_3[s,j,l]+N2_exctrate1[s,j,e,l]*ener[e]
;          N2_tot_2[s,j]=N2_tot_2[s,j]+N2_exctrate1[s,j,e,l]*ener[e]
;          N2_tot_1[j]=N2_tot_1[j]+N2_exctrate1[s,j,e,l]*ener[e]
;
;
;          N2_rate_3[s,j,l]=N2_rate_3[s,j,l]+N2_exctrate1[s,j,e,l]
;          N2_rate_2[s,j]=N2_rate_2[s,j]+N2_exctrate1[s,j,e,l]
;          N2_rate_1[j]=N2_rate_1[j]+N2_exctrate1[s,j,e,l]
;
;
;
;
;        endfor
;      endfor
;    endfor
;
;  endfor
;
;  avg_ener_3=N2_tot_3/N2_rate_3
;  avg_ener_2=N2_tot_2/N2_rate_2
;  avg_ener_1=N2_tot_1/N2_rate_1
;
;
;
;
;end
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;pro max_ener,N2_exctrate1,ener,maxener,maxener_states,sdener,sdener_states,N2_rate_1,N2_rate_2
;  nei=(size(N2_exctrate1))[1]
;  jmax=(size(N2_exctrate1))[2]
;  nbins=(size(N2_exctrate1))[3]
;  lmax=19 ;(size(N2_exctrate1))[4]
;
;
;
;  
;  N2_rate_1 = fltarr(nei,jmax,nbins)
;  N2_rate_2 = fltarr(jmax,nbins)
;  
;  N2_tot_1 = fltarr(nei,jmax,nbins)
;  N2_tot_2 = fltarr(jmax,nbins)
;
;
;  maxener = fltarr(jmax)
;  maxener_states = fltarr(nei,jmax)
;  
;  sdener = fltarr(jmax)
;  sdener_states = fltarr(nei,jmax)
;
;  
;;  N2_rate_1=reform(total(N2_exctrate1,4))
;;  N2_rate_2=reform(total(N2_rate_1,1))
;;  
;  for l=0,lmax-1 do begin  ;wavelength bins
;
;    for e=0,nbins-1 do begin
;      for j=0,jmax-1 do begin  ;altitude
;
;        for s=0,nei-1 do begin  ;states
;
;          
;
;           N2_rate_1[s,j,e]=N2_rate_1[s,j,e]+N2_exctrate1[s,j,e,l]
;           N2_rate_2[j,e]=N2_rate_2[j,e]+N2_exctrate1[s,j,e,l]
;          
;           N2_tot_1[s,j,e]=N2_tot_1[s,j,e]+N2_exctrate1[s,j,e,l]*ener[e]
;           N2_tot_2[j,e]=N2_tot_2[j,e]+N2_exctrate1[s,j,e,l]*ener[e]
;          
;           
;
;
;
;
;        endfor
;      endfor
;    endfor
;
;  endfor
;
;
; 
;    for j=0,jmax-1 do begin  ;altitude
;
;         mi1=max(reform(N2_rate_2[j,*]),ind1)
;         maxener[j]=ener[ind1]
;         
;         N2_tot=fltarr(nbins)
;         N2_tot[0]=N2_tot_2[j,0]
;         
;          for e=1,nbins-1 do N2_tot[e]=total(N2_tot_2[j,0:e],2)
;;          plot,ener,N2_tot,/xl
;;          stop
;          
;         jj=where(  ( (total(N2_tot_2[j,*])*0.95)-N2_tot )  ge 0,/null)
;;         print,ener[jj[-1]]
;;         stop
;         sdener[j]=ener[jj[-1]]
;  
;    endfor
;    
;    
;                    
;    for j=0,jmax-1 do begin  ;altitude
;
;      for s=0,nei-1 do begin  ;states
;        mi2=max(reform(N2_rate_1[s,j,*]),ind2)
;        maxener_states[s,j]=ener[ind2]
;
;        N2_tot=fltarr(nbins)
;        N2_tot[0]=N2_tot_1[s,j,0]
;        
;        for e=1,nbins-1 do N2_tot[e]=total(N2_tot_1[s,j,0:e],3)
;        
;        jj=where(  ( (total(N2_tot_1[s,j,*])*0.95)-N2_tot )  ge 0,/null)
;        
;        
;        sdener_states[s,j]=ener[jj[-1]]
;
;
;
;      endfor
;    endfor
;                
;
;
;
;end
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;; Testing the energy deposition  due to electron impact
;;NOTE:Plot code for sample_localpe_v6.pro
;;____________________________________________________________________________________________________________________________________________________________________
;; N2 Excitation States
;;                    predissociation
;;0  A 3Σu+   6.169
;;1  B 3Пg    7.353
;;2  W 3∆u    7.362
;;3  B’ 3Σu-  8.165
;;4  a’ 1Σu-  8.399
;;5  a 1Пg    8.549     12%
;;6  C 3Пu    11.032    50%
;;7  E 3Σg+   11.875
;;8  a’’1Σg+  12.255
;;9  b 1Пu    12.500    95%
;;10  c’4 1Σu+ 12.935   10%
;;11  b’ 1Σu+  12.854   84%
;;15  w 1∆u    8.890
;;16  c 1Пu    12.910   99%
;
;;States added to account for dissociation due to electron impact- Strickland 1996
;;12  15.8 eV peak 16.4 eV  100%
;;13  17.3 eV peak  17.4eV 100%
;;14  triplet manifold  11eV 100%
;
;;17  VUV       23.7 eV 100%
;;18  Ryd atoms 40 eV 100%
;;____________________________________________________________________________________________________________________________________________________________________
;
;
;sav_loc="/home/srimoyee/Desktop/nrl_files/sav_files/"
;
;fname=[sav_loc+'N2exctrate1OLD_Jun2021.sav',$              ;0 flare,old,transp ; file has: N2_exctrate1,tflux_all,zz,ener,del,wv1,wv2,sza,lat,lon,tau
;
;       sav_loc+'N2exctrate1TEST_Jun2021.sav',$             ;1 flare,new,transp
;       
;       sav_loc+'N2exctrate1localOLD_Jun2021.sav',$         ;2 flare,old,local
;       sav_loc+'N2exctrate1localTEST_Jun2021.sav',$        ;3 flare,new,local
;       
;       
;       sav_loc+'N2exctrate1OLD_Jun2021_QUIET.sav',$         ;4 quiet,old,transp
;       sav_loc+'N2exctrate1_Jun2021_QUIET.sav',$            ;5 quiet,new,transp
;       
;       sav_loc+'N2exctrate1localOLD_Jun2021_QUIET.sav',$    ;6 quiet,old,local
;       sav_loc+'N2exctrate1local_Jun2021_QUIET.sav',$       ;7 quiet,new,local
;       
;       sav_loc+'N2crosssecOLd_Jun2021.sav',$                ; 8 old cross-section, file has:sigix,sigex,ener
;       sav_loc+'N2cross_secJun2021.sav']                    ; 9 new cross-section
;
;plt_loc="/home/srimoyee/Desktop/nrl_files/newNRL_plots/"
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;restore,fname[8]
;
;sigex_OLD=sigex
;totdisscrss_OLD= sigex_OLD[4,2,*]+sigex_OLD[5,2,*]+sigex_OLD[6,2,*]
;
;restore,fname[9]
;
;totdisscrss_NEW= 0.12*sigex[5,2,*]+0.5*sigex[6,2,*]+0.95*sigex[9,2,*]+0.10*sigex[10,2,*]+0.84*sigex[11,2,*]+sigex[12,2,*]+sigex[13,2,*]+sigex[14,2,*]+$
;                0.99*sigex[16,2,*]+sigex[17,2,*]+sigex[18,2,*] 
;
;
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;;#1 Flare Old Transp
;
;restore,fname[0]
;
;nei=3
;jmax=(size(N2_exctrate1))[2]
;nbins=(size(N2_exctrate1))[3]
;lmax=19 
;
;N2_exctrate1_FOT=N2_exctrate1
;N2_dissrate=fltarr(nei,jmax,nbins,lmax)
;
;tflux_FOT=fltarr(jmax,nbins)
;tflux_FOT=total(reform(tflux_all[*,*,0:18]),3)
;
;
;N2_dissrate[0:2,*,*,*]=N2_exctrate1[4:6,*,*,0:18]
;
;average_ener,N2_dissrate,ener,avg_ener_3_FOT,avg_ener_2_FOT,avg_ener_1_FOT,N2dissrate_FOT
;
;max_ener,N2_dissrate,ener,maxener_FOT,maxener_states_FOT,sdener_FOT,sdener_states_FOT,N2_rate_1_FOT,N2_rate_2_FOT
;
;
;;#2 Flare New Transp
;
;
;restore,fname[1]
;
;nei=11
;jmax=(size(N2_exctrate1))[2]
;nbins=(size(N2_exctrate1))[3]
;lmax=19
;
;N2_exctrate1_FNT=N2_exctrate1
;
;N2_dissrate=fltarr(nei,jmax,nbins,lmax)
;
;N2_dissrate[0:1,*,*,*]=N2_exctrate1[5:6,*,*,0:18]
;N2_dissrate[2:7,*,*,*]=N2_exctrate1[9:14,*,*,0:18]
;N2_dissrate[8:10,*,*,*]=N2_exctrate1[16:18,*,*,0:18]
;
;tflux_FNT=fltarr(jmax,nbins)
;tflux_FNT=total(reform(tflux_all[*,*,0:18]),3)
;
;
;average_ener,N2_dissrate,ener,avg_ener_3_FNT,avg_ener_2_FNT,avg_ener_1_FNT,N2dissrate_FNT
;max_ener,N2_dissrate,ener,maxener_FNT,maxener_states_FNT,sdener_FNT,sdener_states_FNT,N2_rate_1_FNT,N2_rate_2_FNT
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;;#4 Quiet Old Transp
;
;restore,fname[4]
;
;nei=3
;jmax=(size(N2_exctrate1))[2]
;nbins=(size(N2_exctrate1))[3]
;lmax=19
;
;N2_exctrate1_QOT=N2_exctrate1
;N2_dissrate=fltarr(nei,jmax,nbins,lmax)
;
;tflux_QOT=fltarr(jmax,nbins)
;tflux_QOT=total(reform(tflux_all[*,*,0:18]),3)
;
;
;N2_dissrate[0:2,*,*,*]=N2_exctrate1[4:6,*,*,0:18]
;
;average_ener,N2_dissrate,ener,avg_ener_3_QOT,avg_ener_2_QOT,avg_ener_1_QOT,N2dissrate_QOT
;
;max_ener,N2_dissrate,ener,maxener_QOT,maxener_states_QOT,sdener_QOT,sdener_states_QOT,N2_rate_1_QOT,N2_rate_2_QOT
;
;
;;#5 Quiet New Transp
;
;
;restore,fname[5]
;
;nei=11
;jmax=(size(N2_exctrate1))[2]
;nbins=(size(N2_exctrate1))[3]
;lmax=19
;
;N2_exctrate1_QNT=N2_exctrate1
;
;N2_dissrate=fltarr(nei,jmax,nbins,lmax)
;
;N2_dissrate[0:1,*,*,*]=N2_exctrate1[5:6,*,*,0:18]
;N2_dissrate[2:7,*,*,*]=N2_exctrate1[9:14,*,*,0:18]
;N2_dissrate[8:10,*,*,*]=N2_exctrate1[16:18,*,*,0:18]
;
;tflux_QNT=fltarr(jmax,nbins)
;tflux_QNT=total(reform(tflux_all[*,*,0:18]),3)
;
;
;average_ener,N2_dissrate,ener,avg_ener_3_QNT,avg_ener_2_QNT,avg_ener_1_QNT,N2dissrate_QNT
;max_ener,N2_dissrate,ener,maxener_QNT,maxener_states_QNT,sdener_QNT,sdener_states_QNT,N2_rate_1_QNT,N2_rate_2_QNT
;
;;____________________________________________________________________________________________________________________________________________________________________
;
;;Universal Plot variables
;xg=0 ;gridstyles
;yg=0
;xtl=1.  ;ticklen
;ytl=1.
;xsg=1
;ysg=1
;xstl=1.
;ystl=1.
;fntsz=15
;
;colors=['red','dodger blue','spring green','purple','gray','deep pink','firebrick','medium blue','orange','fuchsia','dark khaki','dark slate gray','gold']
;lstyl=[0,2,3,4]
;postop=[0.1,0.5,0.45,0.95]
;posbot=[0.1,0.1,0.45,0.45]
;
;postleft=[0.1,0.1,0.5,0.90]
;postright=[0.55,0.1,0.95,0.90]
;
;thck=5.
;leg_loc=[0.47, 0.85]
;leg_loc_left=[0.37, 0.87]
;maj=['O ','O!L2!N ','N!L2!N ']
;
;
;
;leg_nam=['a!U1!N $\Pi$ !Lg!N', 'C!U3!N $\Pi$ !Lu!N','b!U1!N $\Pi$ !Lu!N',"c'!L4!N !U1!N $\Sigma$ !Lu!N !U+!N",$
;  "b'!U1!N $\Sigma$ !Lu!N !U+!N",'15.8 eV','17.3 eV','Triplet manifold','c!U1!N $\Pi$ !Lu!N','VUV',$
;  'Ryd atoms','Dissociative Ionisation','Dissociative excitation']
;
;leg_nam2=['Flare PE MODEL-with Glow Cross-section','Flare PE Model with New Cross-section',$
;           'Quiet PE MODEL-with Glow Cross-section','Quiet PE Model with New Cross-section']
;
;  ;____________________________________________________________________________________________________________________________________________________________________
;
;
;  ;N2 Electron Impact Dissociation percentage contribution plot
;
;
;
;  xt='Average Energy (eV)'
;  yt='Altitude (km)'
;
;  
;  ;xr=[]
;   yr=[min(zz),max(zz)]
;  titl='Average energy of dissociation of N!L2!N '
;
;
;  w1 = WINDOW(DIMENSIONS=[1000,800])
;
;  d1a=plot(avg_ener_1_FOT,zz,name=leg_nam2[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  d1b=plot(avg_ener_1_FNT,zz,name=leg_nam2[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;
;  leg1 = LEGEND(TARGET=[d1a,d1b], POSITION=leg_loc,/normal)
;
;  ;_____________________________________________________________________________________________________
;  
;  w2 = WINDOW(DIMENSIONS=[1000,800])
;
;  d1a=plot(avg_ener_2_FNT[0,*],zz,name=leg_nam[0],$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[0],linestyle=lstyl[0],/current)
;
;  d1b=plot(avg_ener_2_FNT[1,*],zz,name=leg_nam[1],$
;    thick=thck,color=colors[11],linestyle=lstyl[0],/overplot)
;
;  d1c=plot(avg_ener_2_FNT[2,*],zz,name=leg_nam[2],$
;    thick=thck,color=colors[2],linestyle=lstyl[0],/overplot)
;
;  d1d=plot(avg_ener_2_FNT[3,*],zz,name=leg_nam[3],$
;    thick=thck,color=colors[3],linestyle=lstyl[0],/overplot)
;
;  d1e=plot(avg_ener_2_FNT[4,*],zz,name=leg_nam[4],$
;    thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;  d1f=plot(avg_ener_2_FNT[5,*],zz,name=leg_nam[5],$
;    thick=thck,color=colors[5],linestyle=lstyl[0],/overplot)
;
;  d1g=plot(avg_ener_2_FNT[6,*],zz,name=leg_nam[6],$
;    thick=thck,color=colors[6],linestyle=lstyl[0],/overplot)
;
;  d1h=plot(avg_ener_2_FNT[7,*],zz,name=leg_nam[7],$
;    thick=thck,color=colors[7],linestyle=lstyl[0],/overplot)
;
;  d1i=plot(avg_ener_2_FNT[8,*],zz,name=leg_nam[8],$
;    thick=thck,color=colors[8],linestyle=lstyl[0],/overplot)
;
;  d1j=plot(avg_ener_2_FNT[9,*],zz,name=leg_nam[9],$
;    thick=thck,color=colors[9],linestyle=lstyl[0],/overplot)
;
;  d1k=plot(avg_ener_2_FNT[10,*],zz,name=leg_nam[10],$
;    thick=thck,color=colors[10],linestyle=lstyl[0],/overplot)
;
;  
;  leg1 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f,d1g,d1h,d1i,d1j,d1k], POSITION=leg_loc,/normal)
;O
;  ;_____________________________________________________________________________________________________
;
;  titl='Average energy of dissociation of N!L2!N (GLOW cross-section)'
;  w3 = WINDOW(DIMENSIONS=[1000,800])
;
;  d1a=plot(avg_ener_2_FOT[0,*],zz,name='!U1!N $\Pi$ !Lu!N',$
;    xtitle=xt,ytitle=yt,title=titl,$
;    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;    thick=thck,color=colors[2],linestyle=lstyl[0],/current)
;
;  d1b=plot(avg_ener_2_FOT[1,*],zz,name="b'",$
;    thick=thck,color=colors[4],linestyle=lstyl[0],/overplot)
;
;  d1c=plot(avg_ener_2_FOT[2,*],zz,name='Ryds',$
;    thick=thck,color=colors[10],linestyle=lstyl[0],/overplot)
;
;leg3 = LEGEND(TARGET=[d1a,d1b,d1c], POSITION=leg_loc,/normal)
;
;;_____________________________________________________________________________________________________
;
;
;
;
;xt='Energy (eV)'
;yt='Altitude (km)'
;
;
;xr=[0.1,1.e4]
;yr=[min(zz),max(zz)]
;titl='Energy of maximum dissociation of N!L2!N '
;
;
;w4 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(maxener_FOT,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(maxener_FNT,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(maxener_QOT,zz,name=leg_nam2[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(maxener_QNT,zz,name=leg_nam2[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;
;leg4 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;;_____________________________________________________________________________________________________
;
;xt='Energy (eV)'
;yt='Altitude (km)'
;
;;#1
;;xr=[]
;yr=[min(zz),max(zz)]
;titl='Energy of maximum dissociation of N!L2!N !C(solid-Quiet,dashed-flare) '
;
;
;w5 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(maxener_states_FNT[0,*],zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(maxener_states_FNT[1,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[11],linestyle=lstyl[1],/overplot)
;
;d1c=plot(maxener_states_FNT[2,*],zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],linestyle=lstyl[1],/overplot)
;
;d1d=plot(maxener_states_FNT[3,*],zz,name=leg_nam[3],$
;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;
;d1e=plot(maxener_states_FNT[4,*],zz,name=leg_nam[4],$
;  thick=thck,color=colors[4],linestyle=lstyl[1],/overplot)
;
;d1f=plot(maxener_states_FNT[5,*],zz,name=leg_nam[5],$
;  thick=thck,color=colors[5],linestyle=lstyl[1],/overplot)
;
;d1g=plot(maxener_states_FNT[6,*],zz,name=leg_nam[6],$
;  thick=thck,color=colors[6],linestyle=lstyl[1],/overplot)
;
;d1h=plot(maxener_states_FNT[7,*],zz,name=leg_nam[7],$
;  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)
;
;d1i=plot(maxener_states_FNT[8,*],zz,name=leg_nam[8],$
;  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;
;d1j=plot(maxener_states_FNT[9,*],zz,name=leg_nam[9],$
;  thick=thck,color=colors[9],linestyle=lstyl[1],/overplot)
;
;d1k=plot(maxener_states_FNT[10,*],zz,name=leg_nam[10],$
;  thick=thck,color=colors[10],linestyle=lstyl[1],/overplot)
;
;
;leg5 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f,d1g,d1h,d1i,d1j,d1k], POSITION=leg_loc,/normal)
;
;d1a=plot(maxener_states_QNT[0,*],zz,name=leg_nam[0],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1b=plot(maxener_states_QNT[1,*],zz,name=leg_nam[1],$
;  thick=2,color=colors[11],linestyle=lstyl[0],/overplot)
;
;d1c=plot(maxener_states_QNT[2,*],zz,name=leg_nam[2],$
;  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)
;
;d1d=plot(maxener_states_QNT[3,*],zz,name=leg_nam[3],$
;  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;
;d1e=plot(maxener_states_QNT[4,*],zz,name=leg_nam[4],$
;  thick=2,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d1f=plot(maxener_states_QNT[5,*],zz,name=leg_nam[5],$
;  thick=2,color=colors[5],linestyle=lstyl[0],/overplot)
;
;d1g=plot(maxener_states_QNT[6,*],zz,name=leg_nam[6],$
;  thick=2,color=colors[6],linestyle=lstyl[0],/overplot)
;
;d1h=plot(maxener_states_QNT[7,*],zz,name=leg_nam[7],$
;  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
;
;d1i=plot(maxener_states_QNT[8,*],zz,name=leg_nam[8],$
;  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
;d1j=plot(maxener_states_QNT[9,*],zz,name=leg_nam[9],$
;  thick=2,color=colors[9],linestyle=lstyl[0],/overplot)
;
;d1k=plot(maxener_states_QNT[10,*],zz,name=leg_nam[10],$
;  thick=2,color=colors[10],linestyle=lstyl[0],/overplot)
;
;
;;_____________________________________________________________________________________________________
;
;titl='Energy of maximum dissociation of N!L2!N (GLOW cross-section)!C (solid-Quiet,dashed-Flare)'
;w6 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(maxener_states_FOT[0,*],zz,name='!U1!N $\Pi$ !Lu!N',$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[2],linestyle=lstyl[1],font_size=fntsz,/current)
;totdisscrss_OLD
;d1b=plot(maxener_states_FOT[1,*],zz,name="b'",$
;  thick=thck,color=colors[4],linestyle=lstyl[1],/overplot)
;
;d1c=plot(maxener_states_FOT[2,*],zz,name='Ryds',$
;  thick=thck,color=colors[10],linestyle=lstyl[1],/overplot)
;
;leg6 = LEGEND(TARGET=[d1a,d1b,d1c], POSITION=leg_loc,/normal)
;
;d1a=plot(maxener_states_QOT[0,*],zz,name='!U1!N $\Pi$ !Lu!N',$
;  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)
;
;d1b=plot(maxener_states_QOT[1,*],zz,name="b'",$
;  thick=2,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d1c=plot(maxener_states_QOT[2,*],zz,name='Ryds',$
;  thick=2,color=colors[10],linestyle=lstyl[0],/overplot)
;
;
;;_____________________________________________________________________________________________________
;
;;p1 = surface( transpose(N2_rate_1_OLD[0,*,*]), ener,zz,title='OLD 0',xlog=1)
;;p2 = surface( transpose(N2_rate_1_OLD[1,*,*]), ener,zz,title='OLD 1',xlog=1)
;;p3 = surface( transpose(N2_rate_1_OLD[2,*,*]), ener,zz,title='OLD 2',xlog=1) 
;;
;;q1 = surface( transpose(N2_rate_1_NEW[0,*,*]), ener,zz,title='NEW 0',xlog=1)
;;q2 = surface( transpose(N2_rate_1_NEW[1,*,*]), ener,zz,title='NEW 1',xlog=1)
;;q3 = surface( transpose(N2_rate_1_NEW[2,*,*]), ener,zz,title='NEW 2',xlog=1)
;;q4 = surface( transpose(N2_rate_1_NEW[3,*,*]), ener,zz,title='NEW 3',xlog=1)
;;q5 = surface( transpose(N2_rate_1_NEW[4,*,*]), ener,zz,title='NEW 4',xlog=1)
;;q6 = surface( transpose(N2_rate_1_NEW[5,*,*]), ener,zz,title='NEW 5',xlog=1)
;;q7 = surface( transpose(N2_rate_1_NEW[6,*,*]), ener,zz,title='NEW 6',xlog=1)
;;q8 = surface( transpose(N2_rate_1_NEW[7,*,*]), ener,zz,title='NEW 7',xlog=1)
;;q9 = surface( transpose(N2_rate_1_NEW[8,*,*]), ener,zz,title='NEW 8',xlog=1)
;;q10 = surface( transpose(N2_rate_1_NEW[9,*,*]), ener,zz,title='NEW 9',xlog=1)
;;q11 = surface( transpose(N2_rate_1_NEW[10,*,*]), ener,zz,title='NEW 10',xlog=1)
;
;xt=' Energy (eV)'
;yt='Photoelectron Dissociation rate (cm!U-3!N s!U-1!N)'
;
;xr=[1.,1.e4]
;;yr=[0.,2.]
;titl='Dissociation rate of N!L2!N at altitiude '+string(zz[17])+'km'
;
;
;w10a = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(ener,N2_rate_2_FOT[17,*],name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;  
;d1b=plot(ener,N2_rate_2_FNT[17,*],name=leg_nam2[1],$
;      thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;    
;d1c=plot(ener,N2_rate_2_QOT[17,*],name=leg_nam2[2],$
;      thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,N2_rate_2_QNT[17,*],name=leg_nam2[3],$
;      thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;
;
;
;titl='Dissociation rate of N!L2!N at altitiude '+string(zz[60])+'km'
;yr=[0.,1200.]
;w10b = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(ener,N2_rate_2_FOT[60,*],name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;  
;d1b=plot(ener,N2_rate_2_FNT[60,*],name=leg_nam2[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(ener,N2_rate_2_QOT[60,*],name=leg_nam2[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,N2_rate_2_QNT[60,*],name=leg_nam2[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;
;
;
;titl='Dissociation rate of N!L2!N at altitiude '+string(zz[150])+'km'
;yr=[0.,100.]
;w10c = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(ener,N2_rate_2_FOT[150,*],name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=[],xrange=xr,xlog=1,ylog=1,$;
;;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;  
;d1b=plot(ener,N2_rate_2_FNT[150,*],name=leg_nam2[1],$
;    thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(ener,N2_rate_2_QOT[150,*],name=leg_nam2[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(ener,N2_rate_2_QNT[150,*],name=leg_nam2[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;
;;_____________________________________________________________________________________________________
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;xr=[10,1.e4]
;yr=[min(zz),max(zz)]
;titl='Energy of 95% dissociation of N!L2!N !C (solid-Quiet,dashed-Flare)'
;
;
;w7 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(sdener_states_FNT[0,*],zz,name=leg_nam[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(sdener_states_FNT[1,*],zz,name=leg_nam[1],$
;  thick=thck,color=colors[11],linestyle=lstyl[1],/overplot)
;
;d1c=plot(sdener_states_FNT[2,*],zz,name=leg_nam[2],$
;  thick=thck,color=colors[2],linestyle=lstyl[1],/overplot)
;
;d1d=plot(sdener_states_FNT[3,*],zz,name=leg_nam[3],$
;  thick=thck,color=colors[3],linestyle=lstyl[1],/overplot)
;
;d1e=plot(sdener_states_FNT[4,*],zz,name=leg_nam[4],$
;  thick=thck,color=colors[4],linestyle=lstyl[1],/overplot)
;
;d1f=plot(sdener_states_FNT[5,*],zz,name=leg_nam[5],$
;  thick=thck,color=colors[5],linestyle=lstyl[1],/overplot)
;
;d1g=plot(sdener_states_FNT[6,*],zz,name=leg_nam[6],$
;  thick=thck,color=colors[6],linestyle=lstyl[1],/overplot)
;
;d1h=plot(sdener_states_FNT[7,*],zz,name=leg_nam[7],$
;  thick=thck,color=colors[7],linestyle=lstyl[1],/overplot)
;
;d1i=plot(sdener_states_FNT[8,*],zz,name=leg_nam[8],$
;  thick=thck,color=colors[8],linestyle=lstyl[1],/overplot)
;
;d1j=plot(sdener_states_FNT[9,*],zz,name=leg_nam[9],$
;  thick=thck,color=colors[9],linestyle=lstyl[1],/overplot)
;
;d1k=plot(sdener_states_FNT[10,*],zz,name=leg_nam[10],$
;  thick=thck,color=colors[10],linestyle=lstyl[1],/overplot)
;
;
;leg7 = LEGEND(TARGET=[d1a,d1b,d1c,d1d,d1e,d1f,d1g,d1h,d1i,d1j,d1k], POSITION=leg_loc,/normal)
;
;d2a=plot(sdener_states_QNT[0,*],zz,name=leg_nam[0],$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;
;d2b=plot(sdener_states_QNT[1,*],zz,name=leg_nam[1],$
;  thick=2,color=colors[11],linestyle=lstyl[0],/overplot)
;
;d2c=plot(sdener_states_QNT[2,*],zz,name=leg_nam[2],$
;  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)
;
;d2d=plot(sdener_states_QNT[3,*],zz,name=leg_nam[3],$
;  thick=2,color=colors[3],linestyle=lstyl[0],/overplot)
;
;d2e=plot(sdener_states_QNT[4,*],zz,name=leg_nam[4],$
;  thick=2,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d2f=plot(sdener_states_QNT[5,*],zz,name=leg_nam[5],$
;  thick=2,color=colors[5],linestyle=lstyl[0],/overplot)
;
;d2g=plot(sdener_states_QNT[6,*],zz,name=leg_nam[6],$
;  thick=2,color=colors[6],linestyle=lstyl[0],/overplot)
;
;d2h=plot(sdener_states_QNT[7,*],zz,name=leg_nam[7],$
;  thick=2,color=colors[7],linestyle=lstyl[0],/overplot)
;
;d2i=plot(sdener_states_QNT[8,*],zz,name=leg_nam[8],$
;  thick=2,color=colors[8],linestyle=lstyl[0],/overplot)
;
;d2j=plot(sdener_states_QNT[9,*],zz,name=leg_nam[9],$
;  thick=2,color=colors[9],linestyle=lstyl[0],/overplot)
;
;d2k=plot(sdener_states_QNT[10,*],zz,name=leg_nam[10],$
;  thick=2,color=colors[10],linestyle=lstyl[0],/overplot)
;
;;_____________________________________________________________________________________________________
;
;titl='Energy of 95% dissociation of N!L2!N (GLOW cross-section)!C (solid-Quiet,dashed-Flare)'
;w8 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(sdener_states_FOT[0,*],zz,name='!U1!N $\Pi$ !Lu!N',$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[2],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(sdener_states_FOT[1,*],zz,name="b'",$
;  thick=thck,color=colors[4],linestyle=lstyl[1],/overplot)
;
;d1c=plot(sdener_states_FOT[2,*],zz,name='Ryds',$
;  thick=thck,color=colors[10],linestyle=lstyl[1],/overplot)
;
;leg8 = LEGEND(TARGET=[d1a,d1b,d1c], POSITION=leg_loc,/normal)
;
;d2a=plot(sdener_states_QOT[0,*],zz,name='!U1!N $\Pi$ !Lu!N',$
;  thick=2,color=colors[2],linestyle=lstyl[0],/overplot)
;
;
;d2b=plot(sdener_states_QOT[1,*],zz,name="b'",$
;  thick=2,color=colors[4],linestyle=lstyl[0],/overplot)
;
;d2c=plot(sdener_states_QOT[2,*],zz,name='Ryds',$
;  thick=2,color=colors[10],linestyle=lstyl[0],/overplot)
;
;
;;_____________________________________________________________________________________________________
;
;xt=' Energy (eV)'
;yt='Altitude (km)'
;
;
;;xr=[0.1,1.e4]
;yr=[min(zz),max(zz)]
;titl='Energy of 95% dissociation of N!L2!N '
;
;
;w4 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(sdener_FOT,zz,name=leg_nam2[0],$
;  xtitle=xt,ytitle=yt,title=titl,$
;  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
;  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(sdener_FNT,zz,name=leg_nam2[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;
;d1c=plot(sdener_QOT,zz,name=leg_nam2[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(sdener_QNT,zz,name=leg_nam2[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg4 = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;;_____________________________________________________________________________________________________
;
;
;xt='Energy (eV)'
;xr=[0.1,1.e4]
;leg_nam3=['Photoelectron Flux(GLOW,Flare)','Photoelectron Flux(NEW,Flare)','Photoelectron Flux(GLOW,Quiet)','Photoelectron Flux(NEW,Quiet)',$
;          'Cross-section(GLOW)','Cross-section(NEW)']
;
;
;titl='Dissociation of N!L2!N at altitiude '+string(zz[60])+'km'
;w11 = WINDOW(DIMENSIONS=[1000,800])
;
;
;plot_margin = [0.15, 0.25, 0.15, 0.15]
;
;p1 = plot(ener, reform(tflux_FOT[60,*])*del,name=leg_nam3[0], $
;          xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N )",title=titl,$
;          axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
;;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;          thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;          
;          
;p1b=plot(ener,reform(tflux_FNT[60,*])*del,name=leg_nam3[1],$
;            thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;p1c=plot(ener,reform(tflux_QOT[60,*])*del,name=leg_nam3[2],$
;            thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;p1d=plot(ener,reform(tflux_QNT[60,*])*del,name=leg_nam3[3],$
;            thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;p2 = plot(ener, totdisscrss_OLD,name=leg_nam3[4] ,$
;          axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;          thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,totdisscrss_NEW,name=leg_nam3[5],$
;              thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;            
;
;
;a2 = axis('y',target=p2, $
;          location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;          tickfont_size=fntsz,$
;          textpos=1, $                           ; text faces outward
;          tickdir=1, $                           ; ticks face inward
;          title='Cross-section (cm!U2!N)')
;
;l11a=legend(target=[p1,p1b,p1c])
;l11b=legend(target=[p1d,p2,p2b])
;
;
;
;titl='Dissociation of N!L2!N at altitiude '+string(zz[17])+'km'
;w12 = WINDOW(DIMENSIONS=[1000,800])
;
;
;
;
;p1 = plot(ener, reform(tflux_FOT[17,*])*del,name=leg_nam3[0], $
;  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N )",title=titl,$
;  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;p1b=plot(ener,reform(tflux_FNT[17,*])*del,name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;  
;p1c=plot(ener,reform(tflux_QOT[17,*])*del,name=leg_nam3[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;p1d=plot(ener,reform(tflux_QNT[17,*])*del,name=leg_nam3[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;p2 = plot(ener, totdisscrss_OLD,name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,totdisscrss_NEW,name=leg_nam3[5],$
;  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;a2 = axis('y',target=p2, $
;  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;  tickfont_size=fntsz,$
;  textpos=1, $                           ; text faces outward
;  tickdir=1, $                           ; ticks face inward
;  title='Cross-section (cm!U2!N)')
;
;l12a=legend(target=[p1,p1b,p1c])
;l12b=legend(target=[p1d,p2,p2b])
;
;
;
;
;titl='Dissociation of N!L2!N at altitiude '+string(zz[150])+'km'
;w13 = WINDOW(DIMENSIONS=[1000,800])
;
;
;
;
;p1 = plot(ener, reform(tflux_FOT[150,*])*del,name=leg_nam3[0], $
;  xtitle=xt,ytitle="Photoelectron Flux (cm!U-2!N s!U-1!N )",title=titl,$
;  axis_style=1,yrange=[],ylog=1,xrange=xr,xlog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;p1b=plot(ener,reform(tflux_FNT[150,*])*del,name=leg_nam3[1],$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;  
;p1c=plot(ener,reform(tflux_QOT[150,*])*del,name=leg_nam3[2],$
;    thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;p1d=plot(ener,reform(tflux_QNT[150,*])*del,name=leg_nam3[3],$
;    thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;p2 = plot(ener, totdisscrss_OLD,name=leg_nam3[4] ,$
;  axis_style=0,yrange=[],xrange=xr,xlog=1,ylog=1,margin=plot_margin,$
;  ;          xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[0],font_size=fntsz,/current)
;
;p2b=plot(ener,totdisscrss_NEW,name=leg_nam3[5],$
;  thick=thck,color=colors[1],linestyle=lstyl[0],/overplot)
;
;
;
;a2 = axis('y',target=p2, $
;  location=[max(p2.xrange),0,0], $   ; right axis, data coordinates
;  tickfont_size=fntsz,$
;  textpos=1, $                           ; text faces outward
;  tickdir=1, $                           ; ticks face inward
;  title='Cross-section (cm!U2!N)')
;
;l13a=legend(target=[p1,p1b,p1c])
;l13b=legend(target=[p1d,p2,p2b])
;
;;_____________________________________________________________________________________________________
;
;titl='Total Dissociation rate of N!L2!N 
;yt=' Altitude (km)'
;xt='Photoelectron Dissociation rate (cm!U-3!N s!U-1!N)'
;;yr=[0.,100.]
;w14 = WINDOW(DIMENSIONS=[1000,800])
;
;d1a=plot(N2dissrate_FOT,zz,name='GLOW,Flare',$
;  xtitle=xt,ytitle=yt,title=titl,position=postleft,$
;  xstyle=1,ystyle=1,yrange=[],xrange=[],xlog=1,ylog=0,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,/current)
;
;d1b=plot(N2dissrate_FNT,zz,name='New,Flare',$
;  thick=thck,color=colors[1],linestyle=lstyl[1],/overplot)
;
;d1c=plot(N2dissrate_QOT,zz,name='GLOW,Quiet',$
;  thick=2,color=colors[0],linestyle=lstyl[0],/overplot)
;
;d1d=plot(N2dissrate_QNT,zz,name='New,Flare',$
;  thick=2,color=colors[1],linestyle=lstyl[0],/overplot)
;
;leg10c = LEGEND(TARGET=[d1a,d1b,d1c,d1d], POSITION=leg_loc,/normal)
;
;titl='Ratio of change in  dissociation rate of N!L2!N
;xt='Ratio (GLOW/NEW)'
;e1a=plot(N2dissrate_FOT/N2dissrate_FNT,zz,name='Flare',$
;  xtitle=xt,ytitle=[],title=titl,position=postright,$
;  xstyle=1,ystyle=1,yrange=[],xrange=[1.,2.],xlog=0,ylog=0,$;
;    xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;  thick=thck,color=colors[-2],linestyle=lstyl[1],font_size=fntsz,/current)
;
;
;e1b=plot(N2dissrate_QOT/N2dissrate_QNT,zz,name='Quiet',$
;  thick=2,color=colors[-2],linestyle=lstyl[0],/overplot)
;
;
;leg10c = LEGEND(TARGET=[e1a,e1b], POSITION=leg_loc,/normal)
;
;
;end
;
;
;

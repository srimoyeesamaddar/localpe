; This version tests the SQ'05 +FE line provided in 
; M5_stanbads_myver.pro and X9_stanbands_myver.pro
;testing for a total of 4 time periods
;m5 time=0, time= 316. x9 time=0, time= 720

;____________________________________________________________________________________________________________________________________________________________________

;_______________________________________________Main photoelectron calling function__________________________________________________________________________________

;Electron impact ionisation cross-sections
@localpe_setup_exsect
@ace_etransport  
@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
    
;________________________________________________Model constants and Inputs_____________________________________________________________________________________________________

nbins=176
lmax=46 ;32 46 22
nmaj=3

lat=5.    ;x9 flare
lon=0.0

;
;lat=5.    ;x9 flare 60 degrees
;lon=60.0

;lat=75.    ;x9 flare  85 degrees
;lon=85.0



;
;lat=20.  ;m5 flare low sza
;lon=105.0

;lat=35.  ;m5 flare  60 degrees sza
;lon=40.0
;
;lat=20.  ;m5 flare   85 degrees sza
;lon=15.



;altitude bins
z0=40.
zz=findgen(161)*1.+z0  
;zz=findgen(17)*10.+z0  
;zz=findgen(400)*0.5+z0
zz=findgen(171)*5.+80.  
jmax=n_elements(zz)
dz=(zz(1:*)-zz(0:-2))
dz=[dz[0],dz]

;idate=2002089 ;SQ'05
;idate =2016205  ;m5 flare
idate=2017249   ;x9 flare
;idate=2017253 ;x8 flare

floc='/home/srimoyee/Desktop/nrl_files/sav_files/'

;________________________________________________F10.7 and ap values _____________________________________________________________________________________________________


;restore,floc+'f107datafile.sav'
;;this file was created from f107_values.pro from
;;data available from
;;https://lasp.colorado.edu/lisird/data/penticton_radio_flux/
;; file has: month, day, year, hour, minute, second,yyyyddd,f107d,f107_81
;f107_ii=where(yyyyddd eq idate,count)
;if count gt 0 then begin
;  f107= f107d[f107_ii[0]]
;  f107a=f107_81[f107_ii[0]]
;  print, f107,f107a
;endif else begin
;  f107=70.
;  f107a=70.
;  print,'F107 values not found for the given date, taking default value of (f107, f107-81day) ',f107, f107a
;endelse


;OR USE YOUR OWN VALUES
f107=85.
f107a=85.

ap=15
;_______________________________________________________________________________________________________________________________________________________________________

;output variables - for calculating parameters in each bin

;pepi=fltarr(nmaj,lmax)
;tau1h=fltarr(lmax)
;ii_check=intarr(12)


nstate=20 ;10
pstate=3;10 ;equal to nst

Sp_exctrate1_transp=fltarr(nstate,nmaj,jmax,nbins,lmax)
Sp_ionzrate1_transp=fltarr(nstate,nmaj,jmax,nbins,lmax)
Sp_exctrate1_local=fltarr(nstate,nmaj,jmax,nbins,lmax)
Sp_ionzrate1_local=fltarr(nstate,nmaj,jmax,nbins,lmax)
Sp_photoi=fltarr(nmaj,jmax,lmax)
Sp_photoki=fltarr(nmaj,pstate,jmax,lmax)
tflux_all=fltarr(jmax,nbins,lmax)
pespec_all=fltarr(jmax,nbins,lmax)
eiionz_wv=fltarr(nmaj,jmax,lmax)

min_tau_arr=fltarr(lmax)
zz_min=fltarr(lmax)
;
;............................................................................
;wv_lo= fltarr(lmax)
;wv_hi= fltarr(lmax)
;
;;O
;O_pepi_tot= fltarr(lmax)
;O_pepi_4S= fltarr(lmax)
;O_pepi_2D= fltarr(lmax)
;O_pepi_2P= fltarr(lmax)
;
;
;
;
;;O2
;O2_pepi_tot=  fltarr(lmax)
;O2_pepi_ion=  fltarr(lmax)
;O2_pepi_dion=  fltarr(lmax)
;O2_pepi_dissoc=  fltarr(lmax)
;
;
;;N2
;N2_pepi_tot=  fltarr(lmax)
;N2_pepi_ion=  fltarr(lmax)
;N2_pepi_dion=  fltarr(lmax)
;N2_pepi_dissoc=  fltarr(lmax)

;.................................................................................


;tic
;_______________________________________________________________________________________________________________________________________________________________________


localpe_setup_pxsect
if fix((total(size(first_exsect))) eq 0) then begin
  localpe_setup_exsect
  ace_setup_elasticx,ener,nbins,nmaj,sigs,pe,pi
 
endif

   
    dum=min(abs(ener-20.),ia20)
    dum=min(abs(ener-100.),ia100)
    dum=min(abs(ener-1000.),ia1000)
;____________________________________________________________If you need to loop over time in a day________________________________________________________________________________
;    n_hours=14
;    glow20=fltarr(n_hours)
;    glow100=fltarr(n_hours)
;    glow1000=fltarr(n_hours)
;    ace20=fltarr(n_hours)
;    ace100=fltarr(n_hours)
;    ace1000=fltarr(n_hours)
;    asza=fltarr(n_hours)




    ;for i=0,n_hours-1 do begin
;     i=6  ;i th hour
;     hour= 5  ; m5 flare UT time
     hour =12  ; x9 flare UT time
;    hour=float(i)+6.
    utsec=3600.*hour
    
              
    localpe_neutralatm
;    asza[i]=sza

    localpe_approxeden;,zz,f107,f107a,eden,etemp,idate,lat,lon,utsec
    
;______________________________________________________________________  Solar flux__________________________________________________________________________________________

;______________________________________________________________________________________________________________________________________________________________________-



;for hi=0, lmax-1 do begin  ;solar flux loop   lmax-1
     
      localpe_ssflux;,hi

      localpe_photoionz
  
      localpe,idate,lat,hour,local_zz,local_ener,local_pespec,local_del;,/demo

      ace_etransport
    
;     Printing out parameters for optical depth tau=1
;      ii=where(reform(tau(*,hi)) ge 1.,/null)
;      min_tau=min(tau(ii,hi),min_i)
;      
;      
;      min_dif=min(abs(reform(tau(*,hi))-1.),min_i)
;      print, ii[min_i]
;      ii_check[hi]=ii[min_i];

      ;print, min_i
      
;eiionzk_transp=fltarr(nei,nmaj,jmax)
     
;      print,'Optical Depth: ' ,min_tau
;      print,'Height of unit optical depth:',hi,'.',zz[ii[min_i]],'Wavelength',wv1[hi],'-',wv2[hi]
;     print,wv1[hi],'-',wv2[hi];,zz[ii[min_i]]
;     min_tau_arr[hi]=min_tau
;     zz_min[hi]=zz[ii[min_i]]
;   
;   min_tau_arr[hi]=tau[min_i,hi]
;zz_min[hi]=zz[min_i]
;;
;   
;   print,wv1[hi],wv2[hi],zz[ii[min_i]],eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]];,$   ;O
;           eiionzk_transp[0,0,ii[min_i]]/photoi[0,ii[min_i]],$
;           eiionzk_transp[1,0,ii[min_i]]/photoi[0,ii[min_i]],$
;           eiionzk_transp[2,0,ii[min_i]]/photoi[0,ii[min_i]]

   
   
   
   
;     print,wv1[hi],wv2[hi],zz[ii[min_i]],eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]],$   ;O2
;                eiionzk_transp[0,1,ii[min_i]]/photoi[1,ii[min_i]],$
;                eiionzk_transp[1,1,ii[min_i]]/photoi[1,ii[min_i]],$
;                total(exct_transp[1:5,1,ii[min_i]])/photoi[1,ii[min_i]],$    
;               ( total(exct_transp[1:5,1,ii[min_i]],1)+eiionzk_transp[1,1,ii[min_i]])/photoi[1,ii[min_i]]
     
;      print,wv1[hi],wv2[hi],zz[ii[min_i]],eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]],$    ;N2
;            eiionzk_transp[0,2,ii[min_i]]/photoi[2,ii[min_i]],$
;            eiionzk_transp[1,2,ii[min_i]]/photoi[2,ii[min_i]],$
;             
;           (exct_transp[5,2,ii[min_i]]+exct_transp[6,2,ii[min_i]]+exct_transp[9,2,ii[min_i]]+exct_transp[10,2,ii[min_i]]+$
;              exct_transp[11,2,ii[min_i]]+exct_transp[12,2,ii[min_i]]+exct_transp[13,2,ii[min_i]]+exct_transp[14,2,ii[min_i]]+$
;              exct_transp[16,2,ii[min_i]]+exct_transp[17,2,ii[min_i]]+exct_transp[18,2,ii[min_i]])/photoi[2,ii[min_i]],$
;            
;            (exct_transp[5,2,ii[min_i]]+exct_transp[6,2,ii[min_i]]+exct_transp[9,2,ii[min_i]]+exct_transp[10,2,ii[min_i]]+$
;              exct_transp[11,2,ii[min_i]]+exct_transp[12,2,ii[min_i]]+exct_transp[13,2,ii[min_i]]+exct_transp[14,2,ii[min_i]]+$
;              exct_transp[16,2,ii[min_i]]+exct_transp[17,2,ii[min_i]]+exct_transp[18,2,ii[min_i]]+ eiionzk_transp[1,2,ii[min_i]])/photoi[2,ii[min_i]]

    ;........................................................................................................................        
            
;            wv_lo[hi]= wv1[hi]
;            wv_hi[hi]= wv2[hi]
;            
;            ;O
;            O_pepi_tot[hi]= eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]   
;            O_pepi_4S[hi]= eiionzk_transp[0,0,ii[min_i]]/photoi[0,ii[min_i]]
;            O_pepi_2D[hi]= eiionzk_transp[1,0,ii[min_i]]/photoi[0,ii[min_i]]
;            O_pepi_2P[hi]= eiionzk_transp[2,0,ii[min_i]]/photoi[0,ii[min_i]]
;            
;            
;
;            
;            ;O2
;            O2_pepi_tot[hi]= eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]  
;            O2_pepi_ion[hi]= total(eiionzk_transp[0:3,1,ii[min_i]])/photoi[1,ii[min_i]]
;            O2_pepi_dion[hi]= total(eiionzk_transp[4:6,1,ii[min_i]])/photoi[1,ii[min_i]]
;            O2_pepi_dissoc[hi]=   total(exct_transp[3:5,1,ii[min_i]])/photoi[1,ii[min_i]]
;            
;
;            ;N2
;            N2_pepi_tot[hi]= eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]    
;            N2_pepi_ion[hi]= total(eiionzk_transp[0:2,2,ii[min_i]])/photoi[2,ii[min_i]]
;            N2_pepi_dion[hi]= total(eiionzk_transp[3:5,2,ii[min_i]])/photoi[2,ii[min_i]]
;            N2_pepi_dissoc[hi]= total(exct_transp[4:6,2,ii[min_i]])/photoi[2,ii[min_i]]
            
;            ............................................................................................................................

;            print,wv1[hi],wv2[hi], total(exct_transp[3:5,1,ii[min_i]])/photoi[1,ii[min_i]]   ;O2
;            print,wv1[hi],wv2[hi], total(exct_transp[4:6,2,ii[min_i]])/photoi[2,ii[min_i]]   ;N2
            
            
;      print,wv1[hi],'-',wv2[hi];,eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]],eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]],eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]
;      print,'O:',eiionz[0,ii[min_i]]/photoi[0,ii[min_i]],eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]
;      print,'O2:',eiionz[1,ii[min_i]]/photoi[1,ii[min_i]],eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]
;      print,'N2:',eiionz[2,ii[min_i]]/photoi[2,ii[min_i]],eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]
;      
;         
       
;      pepi[0,hi]=eiionz_transp[0,ii[min_i]]/photoi[0,ii[min_i]]
;      pepi[1,hi]=eiionz_transp[1,ii[min_i]]/photoi[1,ii[min_i]]
;      pepi[2,hi]=eiionz_transp[2,ii[min_i]]/photoi[2,ii[min_i]]
;      tau1h[hi]=zz[ii[min_i]]
      
         
;      min_difO=min(abs(reform(tau_wv(0,*,i))-1.),Oi)
;      min_difO2=min(abs(reform(tau_wv(1,*,i))-1.),O2i)
;      min_difN2=min(abs(reform(tau_wv(2,*,i))-1.),N2i)
;         
;      print,'O:',eiionz[0,Oi]/photoi[0,Oi]
;      print,'O2:',eiionz[1,O2i]/photoi[1,O2i]
;      print,'N2:',eiionz[2,N2i]/photoi[2,N2i]
;      print,""
;         print,hi
;         


;      N2_wv_prediss[*,hi] = exct_transp[5,2,*]+exct_transp[6,2,*]+exct_transp[9,2,*]+exct_transp[10,2,*]+$
;                            exct_transp[11,2,*]+exct_transp[12,2,*]+exct_transp[13,2,*]+exct_transp[14,2,*]+$
;                            exct_transp[16,2,*]+exct_transp[17,2,*]+exct_transp[18,2,*]
;            
     
 

;           N2_disionzrate = N2_disionzrate+reform( eiionzk_transp[3,2,*]+eiionzk_transp[4,2,*]+eiionzk_transp[5,2,*])
;            N2_wv_prediss[*,hi] =exct_transp[4,2,*]+exct_transp[5,2,*]+exct_transp[6,2,*]



;       Sp_exctrate1_transp[*,*,*,*,hi]=  exct_transp1
;       Sp_ionzrate1_transp[*,*,*,*,hi]=  eiionz_transp1
;       Sp_exctrate1_local[*,*,*,*,hi]=   exct_local1
;       Sp_ionzrate1_local[*,*,*,*,hi]=   eiionz_local1
;       Sp_photoi[*,*,hi]= photoi
;       Sp_photoki[*,*,*,hi]=photoki[*,*,*,hi]
;       tflux_all[*,*,hi]=tflux
;       pespec_all[*,*,hi]=pespec
;;      
   
  ;solar flux loop end
;    
;pepi[where(~finite(pepi),/null)]=0.

;  eiionz_wv[*,*,hi]=eiionz_transp

;toc
;____________________________________________________________________________________________________________________________________________________________________

;; Energy conservation tests for photoionisation and photoelectrons
;
;
;;the height integrated photoionisation energy should add up to the solar flux,
;;atleast for the lower wavlength bins
;
;photoi_ht=dblarr(lmax)  ;height integrated photoionisation energy
;
;;photoi_wv:Array[alt, species, wavelength]
;photoi_stot=total(photoi_wv,2)  ;total photoization in all species
;
;  for l=0, lmax-1 do begin
;    photoi_ht(l)=total(reform(photoi_stot(*,l)*dz*1.e5)) ;height integrated ionisation rate
;  endfor
;h=6.62607004e-34 
;c=3.e8
;E=h*c/(0.5*1.e-10*(wv1+wv2))   ;J cm-2 s-1
;photoi_ht=photoi_ht*E;*6.242e18
;ssflux_ht=ssflux*E
;
;print,total(photoi_ht),total(ssflux_ht)
leg_nam=['A (6.17 eV) ',    'BW (8.16 eV)', "B'(11.03 eV)",     'C (8.40 eV)',     "aa'w (12.85 eV)",   '1Pu (14.00 eV)',   "b'(13.75 eV)",       'Ryds (1.85 eV)',  'vib',' ']

colors=['red','dark green','dodger blue','purple',"dark goldenrod","midnight blue","orchid","dim gray"]
thck=3
;
;w1 = WINDOW(DIMENSIONS=get_screen_size())
;
;p1=plot(ener,sigex(0,2,*),name=leg_nam[0],$
;  xtitle='Energy (eV)',ytitle='Cross-section (cm!U-2!N)',title='N!L2!N',$
;  xstyle=1,ystyle=1,xrange=[0.,300.],xlog=0,$;yrange=yr,,
;  thick=thck,color=colors[0],/current)
;
;p2=plot(ener,sigex(1,2,*),name=leg_nam[1],$
;  thick=thck,color=colors[1],/overplot)
;
;p3=plot(ener,sigex(2,2,*),name=leg_nam[2],$
;  thick=thck,color=colors[2],/overplot)
;
;  p4=plot(ener,sigex(3,2,*),name=leg_nam[3],$
;    thick=thck,color=colors[3],/overplot)
;
;  p5=plot(ener,sigex(4,2,*)/1.e-16,name=leg_nam[4],$
;    thick=thck,color=colors[4],/overplot)
;
;  p6=plot(ener,sigex(5,2,*)/1.e-16,name=leg_nam[5],$
;    thick=thck,color=colors[5],/overplot)
;
;  p7=plot(ener,sigex(6,2,*)/1.e-16,name=leg_nam[6],$
;    thick=thck,color=colors[6],/overplot)
;
;  p8=plot(ener,sigex(7,2,*),name=leg_nam[7],$
;    thick=thck,color=colors[7],/overplot)
;
;  leg1 = LEGEND(TARGET=[p1,p2,p3,p4,p5,p6,p7,p8])
  
;  p5=plot(ener,(sigex(3,1,*)+sigex(4,1,*)+sigex(5,1,*))/1.e-18,name=leg_nam[4],$
;    xtitle='Energy (eV)',ytitle='Cross-section (10!U-18!N cm!U2!N)',title='O!L2!N',$
;    xstyle=1,ystyle=1,xrange=[0.,1000.],yrange=[0.,200.],xlog=0,$;
;    thick=thck,color=colors[3])

;  
;save,ener, sigex, filename=floc+'ele_exct_crss.sav'  
;file='/home/srimoyee/Desktop/nrl_files/sav_files/O_elecimpactionz_cross.csv'
;WRITE_CSV, file ,ener,  reform(sigix(0,0,*)),reform(sigix(1,0,*)),reform(sigix(2,0,*))


;;_____________________________________________________________________________________________________________________________
;
;  ;energy conservation test  for transport code
;  
;  eiionz_ht= fltarr(nbins)
;  eiionz_stot=total(eiionz_transp_wv,2)  ;species total   ; eiionz_wv:Array[alt, species, wavelength]
;  
;  for i=0, nbins-1 do $
;    eiionz_ht[i]=total(reform(eiionz_stot(*,i)*dz*1.e5)) ;height integrated electron ionisation rate 
;  
;  eiionz_ht=eiionz_ht;*1.60218e-19*ener  ;J cm-2 s-1
;;  
;;  eiionz_tot=fltarr(lmax)  ;putting the electron ionisation energy in solar flux bins
;;  prim_tot=fltarr(jmax,lmax)
;;  
;;  for l=0,lmax-1 do begin
;;    ii=where(((12397./ener ge wv1[l])and (12397./ener lt wv2[l])),cnt,/null)
;;    if cnt gt 0 then begin $
;;         eiionz_tot[l]=total(eiionz_ht(ii)) 
;;         prim_tot[*,l]=total(primary(*,ii))
;;    endif
;;  endfor
;  
;;  for l=0, lmax-1 do begin
;;    ii=where(((ener ge wv1[l])and (ener lt wv2[l])),cnt,/null)
;;    prim_tot[*,l]=total(primary(*,ii))
;;    
;;    
;;  endfor
;  tflux_ht=fltarr(nbins)
;  tflux_ht= total(tflux,1)
;;  for i=0, nbins-1 do $
;;      tflux_ht(i)=total(reform(tflux(*,i)*del(i)))
;;  
;
;pespec_ht=fltarr(nbins)
;pespec_ht= total(pespec,1)
;
;  prim_ht=fltarr(nbins)
;  for i=0, nbins-1 do $
;      prim_ht(i)=total(reform(primary(*,i)*dz*1.e5))
;
;;____________________________________________________________________________________________________________________________________________________________________
;  
;;Universal Plot variables
;;xg=0 ;gridstyles
;;yg=0
;;xtl=1.  ;ticklen
;;ytl=1.
;;xsg=-1
;;ysg=-1
;;xstl=0.5
;;ystl=0.5
;
;xg=0 ;gridstyles
;yg=0
;xtl=1.  ;ticklen
;ytl=1.
;xsg=1
;ysg=1
;xstl=1.
;ystl=1.
;
;colors=['red','dark green','dodger blue','purple',"dark goldenrod","midnight blue","orchid","dim gray"]
;lstyl=[0,2,3,4]
;fnt_sz=15
;posleft=[0.05,0.15,0.33,0.9]
;posmid=[0.39,0.15,0.68,0.9]
;posright=[0.73,0.15,0.98,0.9]
;thck=3.
;ln_thck=5
;sym_thck=2
;sym_incr=10
;leg_loc=[0.74, 0.90]
;
;leg_loc_left=[0.04, 0.9]
;leg_loc_right=[0.98, 0.9]
;leg_transp;maj=['O ','O!L2!N ','N!L2!N ']
;
;sym=["*","+","tu"]
;
;wv=[0.5,1.0,1.5,2.0,2.5,3.0,4,5,6,8,10,14,18,32,70,155,224,290,320]
;
; p=[ 51649.5 ,     830008.,  5.97678e+06,  7.51428e+06,  1.31009e+07,  5.28521e+07,  7.96660e+07,  1.34595e+08,  4.33168e+08,  6.12483e+08,  3.82449e+09,  6.04859e+09,$
;     3.23430e+09,  2.80280e+09,  5.55741e+09,  8.52751e+09,  1.53705e+10,  1.50474e+10,  1.34524e+10]
;     
;t=[4.24520e+14,  1.18620e+15,  3.10727e+15,  1.05363e+15,  1.89404e+14,  2.38248e+13,  7.94692e+11,  5.04372e+11,  4.53308e+11,  1.79552e+11,  6.67363e+11,  9.89528e+11,$
;   1.04037e+12,  3.54115e+11,  5.91034e+11,  7.19905e+11,  1.15663e+12,  1.03665e+12,  9.30721e+11]
;
;e=[1.16925e+07,  1.12056e+08,  5.65808e+08,  5.49042e+08,  7.80189e+08,  2.50292e+09,  2.90646e+09,  3.99769e+09,  1.02079e+10,  1.10950e+10,  5.25242e+10,  6.08774e+10,$
;    4.49472e+10,  1.92646e+10,  1.50645e+10,  8.81441e+09,  7.89190e+09,  4.66137e+09,  1.00993e+09]
;      
;    w5 = WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())
;
;    N2_pi1=plot(wv,e,name='Height Integrated Photoelectron flux',$
;      xtitle='Wavelength ($\AA$)',ytitle='photons cm!U-2!N s!U-1!N ',$
;      xstyle=1,ystyle=1,yrange=[1.e4,1.e16],xrange=[0.1,1.e3],xlog=1,ylog=1,histogram=1,$;
;      xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
;      thick=thck,color=colors[1],font_size=fnt_sz,/current)
;
;    N2_pi2=plot(wv,t,name='Tflux',histogram=1,$
;      thick=thck,color=colors[3],/overplot)
;    ;
;    N2_pi3=plot(wv,p,name='Primary Electron Flux',histogram=1,$
;      thick=thck,color=colors[4],/overplot)
;      
;      leg1 = LEGEND(TARGET=[N2_pi1,N2_pi2,N2_pi3], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)
;
;    cgdisplay,2000,2000, /free
;  
;    cgplot,ener,eiionz_ht,color='red',thick=5,/xs,/ys,/yl,/xl,yr=[1.e2,1.e10],psym=-46,$
;      ytitle='photons cm-2 s-1 ',xtitle='Electron Energy (eV)',xr=[0.01,1.e5];,output='/home/srimoyee/Desktop/pic_sept30_2020/pic_0.jpeg'
;    cgplot, ener,tflux_ht*del,color='dodger blue',/overplot,thick=5
;    cgplot, ener,prim_ht,color='green',/overplot,thick=5
;    cglegend, titles=['Height Integrated Photoelectron flux','Tflux','Primary'],$
;      psyms=[15], symcolors=['red','dodger blue','green'],$
;      alignment=1,/box,Location=[0.93, 0.9]
;
;;print,hi,wv1[hi],wv2[hi],total(prim_ht),total(tflux_ht*del),total(eiionz_ht)
;p[hi-3]=total(prim_ht)
;t[hi-3]=total(tflux_ht*del)
;e[hi-3]=total(eiionz_ht)
; endfor  ;COMMENT
; header= ["Energy",'ABW', "B'",  'C',  "aa'w",   '1Pu', "b'","Ryds"]
; file='/home/srimoyee/Desktop/nrl_files/sav_files/N2_elecimpact_cross3.csv'
; WRITE_CSV, file , ener, sigex(0,2,*),sigex(1,2,*),sigex(2,2,*),sigex(3,2,*),sigex(4,2,*),sigex(5,2,*),sigex(6,2,*),header=header
;
;
; header= ["Energy",'X', "A",  'B',  "D",   'C', "40 eV"]
; file='/home/srimoyee/Desktop/nrl_files/sav_files/N2_ionzimpact_cross.csv'
; WRITE_CSV, file , ener, sigix(0,2,*),sigix(1,2,*),sigix(2,2,*),sigix(3,2,*),sigix(4,2,*),sigix(5,2,*),header=header

;header= ["Energy","X2Pig",    "a4Piu", "A2Piu",    "b4Sigmag",  "B2Sigmag", "c4Sigmau",  "37eV"]
;file='/home/srimoyee/Desktop/nrl_files/sav_files/O2_ionzcross.csv'
;WRITE_CSV, file , ener, sigix(0,1,*),sigix(1,1,*),sigix(2,1,*), sigix(3,1,*),sigix(4,1,*),sigix(5,1,*),sigix(6,1,*),header=header

;
;print,p
;print,t
;print,e
;            
;            ;O
;            O_pepi_tot[where(~finite(O_pepi_tot),/null)]=0.   
;            O_pepi_4S[where(~finite(O_pepi_4S),/null)]=0.
;            O_pepi_2D[where(~finite(O_pepi_2D),/null)]=0.
;            O_pepi_2P[where(~finite(O_pepi_2P),/null)]=0.
;            
;            
;
;            
;            ;O2
;            O2_pepi_tot[where(~finite(O2_pepi_tot),/null)]= 0.  
;            O2_pepi_ion[where(~finite(O2_pepi_ion),/null)]= 0.
;            O2_pepi_dion[where(~finite(O2_pepi_dion),/null)]= 0.
;            O2_pepi_dissoc[where(~finite(O2_pepi_dissoc),/null)]= 0.
;            
;
;            ;N2
;            N2_pepi_tot[where(~finite(N2_pepi_tot),/null)]= 0.  
;            N2_pepi_ion[where(~finite(N2_pepi_ion),/null)]= 0.
;            N2_pepi_dion[where(~finite(N2_pepi_dion),/null)]= 0.
;            N2_pepi_dissoc[where(~finite(N2_pepi_dissoc),/null)]= 0.
;            
;____________________________________________________________________________________________________________________________________________________________________

; output- for running each solar flux bin at a time

;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_m5time0.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_m5time316.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_x9time0.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'NRL_pepitab_x9time720.sav'



;sza_deg=sza*!radeg
;save, wv_lo,wv_hi,sza_deg,O_pepi_tot,O_pepi_4S,O_pepi_2D,O_pepi_2P,$
;      O2_pepi_tot,O2_pepi_ion,O2_pepi_dion,O2_pepi_dissoc,$
;      N2_pepi_tot,N2_pepi_ion,N2_pepi_dion,N2_pepi_dissoc,$
;      filename=floc+'NRL_pepi_final.sav'




;save, photoi, eiionz_transp,pepi, zz, ssflux, wv1,wv2,zcol,eden,etemp, filename=floc+'SQ05_zcol_neut3.sav'
;save,  ssflux, wv1,wv2, filename=floc+'NEWFLAREflux.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_org.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_iri.sav'
;save,pepi,tau1h,wv1,wv2, filename=floc+'SQ05_pepitab_neut1.sav'
;save,pepi,tau1h,wv1,wv2,sza,filename=floc+'SQ05_sza6_4.sav'

;save,eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'SQpeakm5_wv_szalo_test.sav'
;save,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'SQM5flare_Jan2021.sav'

;save,eionz_wv_loc,eionz_wv_transp ,photi_wv ,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'peakm5_oct2020.sav'

;save,eionz_wv_transp ,photi_wv ,N2_wv_prediss,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'N2predissrate_May2021.sav'
;save,N2_exctrate,N2_disionzrate,zz,wv1,wv2,sza,lat,lon,tau,filename =floc+'N2exctrate_Jun2021.sav'

wv_lo=wv1
wv_hi=wv2
sza_deg=sza*!radeg
;save,Sp_exctrate1_transp,Sp_ionzrate1_transp,Sp_exctrate1_local,Sp_ionzrate1_local,Sp_photoi,Sp_photoki,$
;      tflux_all,pespec_all,zz,ener,del,wv1,wv2,sza,lat,lon,tau,$
;      zmaj,sigix,sigex,filename =floc+'SQX9_NC_17JUL2022.sav'
;;      
;save,sigix,sigex,ener,filename =floc+'N2crosssecOLd_Jun2021.sav'
;save,min_tau_arr,zz_min,wv1,wv2,filename =floc+'UnitOptdepth_Sept2022.sav'
save,zz,photoi,sza_deg,filename =floc+'Piratesza0_Sept2022.sav'
 
 
;output- for running the entire solar flux in one go
;
;pepi=fltarr(nmaj,lmax)
;pepi=eiionz_transp/photoi
;pepi[where(~finite(pepi),/null)]=0.

;save, photoi, photoi_wv,eiionz_transp, zz, ssflux, flux,wv1,wv2, zcol,filename=floc+'peakm5.sav'

;save, photoi, photoi_wv,eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'peakx9.sav'
;
;save,eionz_wv, zz, ssflux, flux,wv1,wv2,filename=floc+'peakm5_eiionz_v2.sav'
;save, eionz_wv, zz, ssflux, flux,wv1,wv2, filename=floc+'peakx9_eiionz_v2.sav'

;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'m5_time0.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'m5_time316.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'x9_time0.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv1,wv2, filename=floc+'x9_pt5-4.sav'

;save, photoi, eiionz_transp, zz, ssflux, wv_lo,wv_hi,sza_deg,filename=floc+'NRL_m5.sav'
;save, photoi, eiionz_transp, zz, ssflux, wv_lo,wv_hi,sza_deg, filename=floc+'NRL_x9.sav'

;save,  photoi_wv, tau,zz,lat,lon,sza_deg,wv1,wv2,filename=floc+'NRLX8flare_Aug2021_1.sav'
;save,  eiionz_wv,filename=floc+'NRLX8flare_Aug2021_2.sav'
;___________________________________________________________________________________________________
;z=zz*1.e5
;am=[16., 32., 28.] ; amu for O, O2 anf N2
;kb= 1.38064852e-16   ; Boltzmann's constant ergs per k
;
;re=6378.e5 ; radius of the Earth in cm
;g=980.  ;acceleration due to gravity in cm^2/sec
;gr=g*(re/(re+z))^2 ; account for small change in g with altitude
;
;amave=(.2*32 + .79*28 + .01*40.) ; average mass of atmosphere, used for lower altitudes
;gperamu=1.66054*1.e-24  ;coverting amu to grams
;
;
;;    Assumption: at low altitudes, every constituent follows
;;                a single scale height, above that, they follow
;;                their own, assume change occurs at 100km
;
;
;ilow=z le 100.*1e5
;ihi=1-ilow
;NEW
;H=fltarr(3,n_elements(zz))
;for i=0,2 do begin
;   amass=amave*ilow + am[i]*ihi     ;mass for every height
;
;   H[i,*]=(kb*zt/(amass*gperamu*gr))/1.e5 
;   
;
;endfor
;___________________________________________________________________________________________________

;Theoritical calculation of scale height H
 
;z_0=0.
;H=7.
;n0=1.90120e19
;n0_1=1.90120e19;2.7e19
;chi=0
;abs_n2=[2.48E-05,1.47E-04,4.57E-04,1.02E-03,1.90E-03,4.00E-03,8.55E-03,$
;  1.55E-02,3.12E-02,6.32E-02,0.14,0.3]*1e-18
;  
;;zmax=[43.,55.5,63.5,69.1,73.3,78.7,84.,88.1,93,98,103.5,109.]  ;km
;; 
;;;nz=n0*exp((zmax-z_0)/H)
;;  
;;  
;;n0_2=[1.63198e+19,3.01390e+19,5.72532e+19,1.13414e+20,1.85562e+20,2.44977e+20,6.30649e+20,1.83762e+21,1.57760e+21,$
;;    9.80887e+20,1.19087e+21,5.36025e+20]
;zmax_z1=fltarr(12)
;;zmax_z2=fltarr(12)
;
;
;;
;zmax_z=z_0+H*alog(H*1.e5*n0*abs_n2/cos(chi))
;;
;for i=0,11 do begin
;  zmax_z1[i]=z_0+H[2,ii_check[i]]*alog(H[2,ii_check[i]]*1.e5*n0_1*abs_n2[i])
;  print, H[2,ii_check[i]]
;;  zmax_z2[i]=z_0+H[2,ii_check[i]]*alog(H[2,ii_check[i]]*1.e5*n0_2[i]*abs_n2[i])
;
;endfor
;;zmax=[45.,58,65.5,71,75,79,83.5,87,90.5,94.5,99.,103.]
;;n_0=fltarr(12)
;;for i=0,11 do $
;;   n_0[i]=(1/(H[2,ii_check[i]]*1.e5*abs_n2[i]))*exp((zmax[i]-z_0)/H[2,ii_check[i]])
;
;
;


;h=[28.821991430024,$
;31.696343433164,$
;28.8580941239989,$
;34.1378393351801,$
;34.8455653159535,$
;35.5617817390518,$
;36.3108504755989,$
;36.8286440192398,$
;37.33480886947,$
;37.8383934021851,$
;38.2284248802575,$
;38.4911364573313,$
;1.72404580832459,$
;39.1205633172018,$
;37.1962908237929,$
;32.7503365012033,$
;43.1166835822067,$
;68.428042820609,$
;23.1923986506764]
;
;phisto = PLOT(wv1[0:18], h , /CURRENT,  $
;
;TITLE='Percentage change in Dissocitaion Pe/Pi', XTITLE='Wavelength (angstrom)', YTITLE='Percentage change', $
;   COLOR='red',xrange=[0.45,540.],thick=5,xlog=1,ylog=0,yr=[0.,80],histogram=1);,$
;  symbol='*', SYM_COLOR = "blue",SYM_FILLED = 1,SYM_FILL_COLOR = 0,sym_thick=5)
;;  
;  
;  ; Create the first BARPLOT
;  b1 = BARPLOT(wv1[0:18], h, INDEX=0, NBARS=2, $
;    FILL_COLOR='green', YRANGE=[0, 100], YMINOR=0, $
;    YTITLE='Percentage change ', XTITLE='Wavelength (angstrom)', $
;    TITLE = 'Percentage change in Dissocitaion Pe/Pi' )
;
;____________________________________________________________________________________________________

;____________________________________________________________________________________________________________________________________________________________________

end
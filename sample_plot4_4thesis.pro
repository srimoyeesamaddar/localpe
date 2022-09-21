;plotting photons cross-section and also making tables

const=1;1.e-16
;For Table
restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav'
;file has: wave_gcm1,wave_gcm2,meanpeakx9,peakm5

;New Files
probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states.sav';THIS ONE WAS USED 17JUL2022
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
restore, probstate_crss
restore, ionabs_crss

;print,'O'
;for i=0, n_elements(wave_gcm1)-1 do $
;  print,wave_gcm1[i]*0.1,wave_gcm2[i]*0.1,sigab_o[i],sigi_o[i],o_prob_state[1,i],o_prob_state[2,i],$
;        o_prob_state[3,i],o_prob_state[4,i],o_prob_state[5,i],o_prob_state[6,i]

;print,'O2'
;for i=0, n_elements(wave_gcm1)-1 do $
;     print,wave_gcm1[i]*0.1,wave_gcm2[i]*0.1,sigab_o2[i],sigi_o2[i],$
;          o2_prob_state[1,i],o2_prob_state[2,i],$
;        o2_prob_state[3,i],o2_prob_state[4,i],o2_prob_state[5,i],o2_prob_state[6,i],$
;        o2_prob_state[7,i],o2_prob_state[8,i],o2_prob_state[9,i]


print,'N2'
for i=0, n_elements(wave_gcm1)-1 do $
;     print,wave_gcm1[i]*0.1,wave_gcm2[i]*0.1,sigab_n2[i],sigi_n2[i];
      print,wave_gcm1[i]*0.1,wave_gcm2[i]*0.1,n2_prob_state[1,i],n2_prob_state[2,i],$
        n2_prob_state[3,i],n2_prob_state[4,i],n2_prob_state[5,i],n2_prob_state[6,i],$
        n2_prob_state[7,i],n2_prob_state[8,i],n2_prob_state[9,i],n2_prob_state[10,i]


;____________________________________________________________________________________________________________________________________________________________________

;Universal Plot variables
xg=0 ;gridstyles
yg=0
xtl=1.  ;ticklen
ytl=1.
xsg=1
ysg=1
xstl=1.
ystl=1.
fntsz=14
lgdsz=12
fntstyl=1 ;Bold 0-Normal
transp=100

colors=['red','dodger blue','spring green','fuchsia','yellow','deep pink','firebrick','medium blue','orange','purple','dark khaki','blue','gold']
lstyl=[2,0,3,4]
postop=[0.1,0.5,0.45,0.95]
posbot=[0.1,0.1,0.45,0.45]

postleft=[0.1,0.1,0.5,0.90]
postright=[0.55,0.1,0.95,0.90]

thck=5.
leg_loc=[0.47, 0.85]
leg_loc_left=[0.37, 0.87]
maj=['O ','O!L2!N ','N!L2!N ']

;____________________________________________________________________________________________________________________________________________________________________

;Base Cross-section file
restore,'/home/srimoyee/Desktop/nrl_files/sav_files/ion&abs_cross_files.sav'
;file has:
;wave_o,wave_o2,wave_n2,totionz_o,totionz_o2,totionz_n2,totabs_o,totabs_o2,totabs_n2

;Base States
restore,'/home/srimoyee/Desktop/nrl_files/sav_files/baseline_states.sav'
;File has:
;O_crss_states,O2_crss_states,N2_crss_states

;O: wavelength,4So,2Do,2Po,4P,2P,K

;O2:  wavelength,X,a+A,b,B,2pi+c,2sig,33 eV;2,4sig;K

;N2: wavelength,X,A,B,C, F, G+E,HP,H, N2++,K



oplus_states=['!U4!NS!UO!N','!U2!ND!UO!N','!U2!NP!UO!N','!U4!NP','!U2!NP','(1s)!U-1!N']
n2plus_states=['$X^2\Sigma_g^+$','$A^2\Pi_u$','$B^2\Sigma_u^+$','$C^2\Sigma_u^+$',$
'$F^2\Sigma_g^+$','$G^2\Sigma_g^+$+$E^2\Sigma_u^+$',"$H'^2\Sigma_g^+$",'H ','N!L2!N!U++!N','K Eject']

o2plus_states=['$X^2\Pi_g$','$a^4\Pi_u$+$A^2\Pi_u$','$b^4\Sigma_g^-$','$B^2\Sigma_g^-$',$
  '$^2\Pi_u$+$c^4\Sigma_u^-$','$^2\Sigma_u^+$',"33eV",'$^2,^4\Sigma_g^-$','K eject']
plotmargin=[0.2,0.2,0.2,0.2]
xt='Wavelength (nm)'
yt='Cross-section (cm!U2!N)';10!U-16!N 


xr=[1,100]
yr=[1.e-21,1e-16]
w7a = WINDOW(DIMENSIONS=[450,450])
  titl='(N!L2!N!U+!N)'+n2plus_states[0]+'Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[1,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7b = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,1e-16]
  titl='(N!L2!N!U+!N)'+n2plus_states[1]+'Ionization!C!CCross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[2,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7c = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,1e-17]
  titl='(N!L2!N!U+!N)'+n2plus_states[2]+'Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[3,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7d = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,2e-18]
  titl='(N!L2!N!U+!N)'+n2plus_states[3]+'Ionization!C!CCross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[4,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7e = WINDOW(DIMENSIONS=[450,450])
  titl='(N!L2!N!U+!N)'+n2plus_states[4]+'Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[5,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7f = WINDOW(DIMENSIONS=[450,450])
  titl='(N!L2!N!U+!N) '+n2plus_states[5]+'!C!C Ionization Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[6,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7g = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,3.e-18]
  titl='(N!L2!N!U+!N) '+n2plus_states[6]+'Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[7,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7h = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,1.e-18]
  titl='(N!L2!N!U+!N) '+n2plus_states[7]+'Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[8,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7i = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,1.e-18]
  titl='(N!L2!N!U+!N) '+n2plus_states[8]+' Ionization!C!C Cross-section '
  a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[9,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w7j = WINDOW(DIMENSIONS=[450,450])
yr=[0.,1.2]
xr=[1.,10]
const=1.e-18
yt='Cross-section (10!U-18!Ncm!U2!N)'
titl='(N!L2!N!U+!N) '+n2plus_states[9]+' Ionization!C!C Cross-section '
a1= plot(n2_crss_states[0,*]*0.1,n2_crss_states[10,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,ylog=0,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)


const=1.
xr=[1,200]
yr=[1.e-21,1.e-16]
yt='Cross-section (10!U-18!Ncm!U2!N)'

w8a = WINDOW(DIMENSIONS=[450,450])
  titl='(O!L2!N!U+!N) '+o2plus_states[0]+' Ionization!C!C Cross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[1,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8b = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,2.e-17]
titl='(O!L2!N!U+!N) '+o2plus_states[1]+'!C!CIonization Cross-section '
a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[2,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8c = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-21,2.e-17]
titl='(O!L2!N!U+!N) '+o2plus_states[2]+' Ionization!C!CCross-section '
a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[3,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8d = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,1.e-17]
titl='(O!L2!N!U+!N)'+o2plus_states[3]+' Ionization!C!CCross-section '
a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[4,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8e = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,1.e-17]
titl='(O!L2!N!U+!N) '+o2plus_states[4]+'!C!CIonization Cross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[5,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8f = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,2.e-18]
  titl='(O!L2!N!U+!N) '+o2plus_states[5]+'!C!CIonization Cross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[6,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8g = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,2.e-18]
  titl='(O!L2!N!U+!N) '+o2plus_states[6]+' Ionization!C!CCross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[7,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8h = WINDOW(DIMENSIONS=[450,450])
  titl='(O!L2!N!U+!N) '+o2plus_states[7]+' Ionization!C!CCross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[8,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w8a = WINDOW(DIMENSIONS=[450,450])
const=1.e-18
xr=[0.,10.]
yr=[0.,1.2]
  titl='(O!L2!N!U+!N) '+o2plus_states[8]+' Ionization!C!CCross-section '
  a1= plot(o2_crss_states[0,*]*0.1,o2_crss_states[9,*]/const,$
    xtitle=xt,ytitle=yt,title=titl,$
    xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,ylog=0,$;
    thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)


const=1.
xr=[1,100]
yr=[1.e-20,1.e-17]


w9a = WINDOW(DIMENSIONS=[450,450])
            titl='(O!U+!N) '+oplus_states[0]+' Ionization!C!CCross-section '
            a1= plot(o_crss_states[0,*]*0.1,o_crss_states[1,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w9b = WINDOW(DIMENSIONS=[450,450])
          titl='(O!U+!N) '+oplus_states[1]+' Ionization!C!CCross-section '
          a1= plot(o_crss_states[0,*]*0.1,o_crss_states[2,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w9c = WINDOW(DIMENSIONS=[450,450])
          titl='(O!U+!N) '+oplus_states[2]+' Ionization!C!CCross-section '
          a1= plot(o_crss_states[0,*]*0.1,o_crss_states[3,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w9d = WINDOW(DIMENSIONS=[450,450])
yr=[1.e-20,1.e-18]
titl='(O!U+!N) '+oplus_states[3]+' Ionization!C!CCross-section '
a1= plot(o_crss_states[0,*]*0.1,o_crss_states[4,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w9e = WINDOW(DIMENSIONS=[450,450])
          titl='(O!U+!N)'+oplus_states[4]+' Ionization!C!CCross-section '
          a1= plot(o_crss_states[0,*]*0.1,o_crss_states[5,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w9f = WINDOW(DIMENSIONS=[450,450])
const=1.e-18
xr=[0.,10]
yr=[0.,0.6]
titl='(O!U+!N) '+oplus_states[5]+' Ionization!C!CCross-section '
a1= plot(o_crss_states[0,*]*0.1,o_crss_states[6,*]/const,$
            xtitle=xt,ytitle=yt,title=titl,$
            xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,ylog=0,$;
            thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)






;____________________________________________________________________________________________________________________________________________________________________

const=1


yr=[1.e-22,1.e-15]
xr=[0.05,500]


w1 = WINDOW(DIMENSIONS=[450,450])
titl='Absorption Cross-section of N!L2!N '
a1= plot(wave_n2*0.1,totabs_N2/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
;         xtickvalues=[0.1,1,10,100,1000],xtickname=['0.1','1','10','100','1000'],$
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w2 = WINDOW(DIMENSIONS=[450,450])
       titl='Absorption Cross-section of O!L2!N '
       a1= plot(wave_o2*0.1,totabs_o2/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w3 = WINDOW(DIMENSIONS=[450,450])
       titl='Absorption Cross-section of O '
       a1= plot(wave_o*0.1,totabs_o/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w4 = WINDOW(DIMENSIONS=[450,450])
       titl='Ionization Cross-section of N!L2!N '
       a1= plot(wave_n2*0.1,totionz_N2/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w5 = WINDOW(DIMENSIONS=[450,450])
       titl='Ionization Cross-section of O!L2!N '
       a1= plot(wave_o2*0.1,totionz_o2/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)

w6 = WINDOW(DIMENSIONS=[450,450])
       titl='Ionization Cross-section of O '
       a1= plot(wave_o*0.1,totionz_o/const,$
         xtitle=xt,ytitle=yt,title=titl,$
         xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,ylog=1,$;
         thick=thck,color=colors[0],linestyle=lstyl[1],font_size=fntsz,font_style=fntstyl,margin=plotmargin,/current)



;____________________________________________________________________________________________________________________________________________________________________  
end

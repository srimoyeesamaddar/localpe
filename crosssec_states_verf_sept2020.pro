;Model Cross-section files

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_conway.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2

restore, probstate_crss
restore, ionabs_crss

NRL_wv=spec_wave


nmaj=3
lmax= n_elements(spec_wave)

O_prob_state_nrl=  O_prob_state
O2_prob_state_nrl= O2_prob_state
N2_prob_state_nrl= N2_prob_state

sigionx_nrl=fltarr(nmaj,lmax)
sigabs_nrl=fltarr(nmaj,lmax)

sigabs_nrl(0,*)=sigab_o
sigabs_nrl(1,*)=sigab_o2
sigabs_nrl(2,*)=sigab_n2

sigionx_nrl(0,*)=sigi_o
sigionx_nrl(1,*)=sigi_o2
sigionx_nrl(2,*)=sigi_n2



;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________
;NRL witih all states

probstate_crss= '/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_conway_allstates.sav'
restore, probstate_crss


O_prob_state_nrlall=  O_prob_state
O2_prob_state_nrlall= O2_prob_state
N2_prob_state_nrlall= N2_prob_state

; O data
;Column 1: wavelength
;Column 2:4So
;Column 3: 2Do
;Column 4: 2Po
;Column 5: 4P
;Column 6: 2P
;Column 7: K
;
;O2 data
;Column 1: wavelength
;Column 2:   X
;Column 3:   a+A
;Column 4:    b
;Column 5:  B
;Column 6:   2pi+c
;Column 7:  2sig
;Column 8:   33 eV
;Column 9: 2,4sig
;Column 10:   K


;N2 data
;Column 1: wavelength
;Column 2: X
;Column 3: A
;Column 4: B
;Column 5: C
;Column 6: F
;Column 7: G+E
;Column 8: HP
;Column 9: H
;Column 10: N2++
;Column 11: K






;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________


;  SQ cross-section files

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_SQ05.sav'
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec_SQ05.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2
;
restore, probstate_crss
restore, ionabs_crss

SQ_wv=spec_wave

lmax= n_elements(spec_wave)

O_prob_state_SQ=  O_prob_state
O2_prob_state_SQ= O2_prob_state
N2_prob_state_SQ= N2_prob_state


sigionx_SQ=fltarr(nmaj,lmax)
sigabs_SQ=fltarr(nmaj,lmax)

sigabs_SQ(0,*)=sigab_o
sigabs_SQ(1,*)=sigab_o2
sigabs_SQ(2,*)=sigab_n2

sigionx_SQ(0,*)=sigi_o
sigionx_SQ(1,*)=sigi_o2
sigionx_SQ(2,*)=sigi_n2





;__________________________________________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________________________________________


plt_loc="/home/srimoyee/Desktop/nrl_files/jul12_plots/x9_"

;Universal Plot variables
;xg=0 ;gridstyles
;yg=0
;xtl=1.  ;ticklen
;ytl=1.
;xsg=-1
;ysg=-1
;xstl=0.5
;ystl=0.5

xg=0 ;gridstyles
yg=0
xtl=1.  ;ticklen
ytl=1.
xsg=1
ysg=1
xstl=1.
ystl=1.

colors=['red','dark green','dodger blue','purple',"dark goldenrod","midnight blue","orchid","dim gray"]
lstyl=[0,2,3,4]
fnt_sz=15
posleft=[0.05,0.15,0.33,0.9]
posmid=[0.39,0.15,0.68,0.9]
posright=[0.73,0.15,0.98,0.9]
thck=3.
ln_thck=5
sym_thck=3
sym_incr=10
leg_loc=[0.25, 0.87]

leg_loc_left=[0.25, 0.87]
leg_loc_right=[0.98, 0.9]
leg_transp=50
maj=['O ','O!L2!N ','N!L2!N ']

sym=["*","+","tu"]

;___________________________________________________________________________________________________________________________________________________________________________________

;Photoionisation plot

O_nam=['NRL 4S','NRL 2D','NRL 2P','SQ 4S','SQ 2D','SQ 2P']
O2_nam=['NRL O!L2!N!U+!N','NRL D.I','SQ O!L2!N!U+!N','SQ D.I','NRL_allstate O!L2!N!U+!N','NRL_allstate D.I']
N2_nam=['NRL N!L2!N!U+!N','NRL D.I','SQ N!L2!N!U+!N','SQ D.I','NRL_allstate N!L2!N!U+!N','NRL_allstate D.I']


xt='Wavelength ($\AA$)'
yt='Cross-section (cm!U2!N)'
titl='Cross-sections'

; yr=[0.,1.1]
; xr=[0.1,32.]
; yr=[0.,2.e-19]





;O

xr=[0.1,20.]
w1= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O_crss1=plot(NRL_wv,O_prob_state_nrl[1,*],name=O_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;yrange=yr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O_crss2=plot(NRL_wv,O_prob_state_nrl[2,*],name=O_nam[1],$
  thick=thck,color=colors[1],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss3=plot(NRL_wv,O_prob_state_nrl[3,*],name=O_nam[2],$
  thick=thck,color=colors[2],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss4=plot(SQ_wv,O_prob_state_SQ[1,*],name=O_nam[3],$
  thick=thck,color=colors[0],symbol=sym[1],sym_thick=sym_thck,/overplot)

O_crss5=plot(SQ_wv,O_prob_state_SQ[2,*],name=O_nam[4],$
  thick=thck,color=colors[1],symbol=sym[1],sym_thick=sym_thck,/overplot)

O_crss6=plot(SQ_wv,O_prob_state_SQ[3,*],name=O_nam[5],$
  thick=thck,color=colors[2],symbol=sym[1],sym_thick=sym_thck,/overplot)





leg1 = LEGEND(TARGET=[O_crss1,O_crss2,O_crss3,O_crss4,O_crss5,O_crss6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)



;__________________________________________________________________________________________________________________________________________________________________________________

;O

xr=[0.1,20.]
w4= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O_crss1=plot(NRL_wv,O_prob_state_nrl[1,*],name=O_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;yrange=yr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O_crss2=plot(SQ_wv,O_prob_state_SQ[1,*],name=O_nam[3],$
  thick=thck,color=colors[0],symbol=sym[1],sym_thick=sym_thck,/overplot)

leg4 = LEGEND(TARGET=[O_crss1,O_crss2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)




w5= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O_crss1=plot(NRL_wv,O_prob_state_nrl[2,*],name=O_nam[1],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;yrange=yr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O_crss2=plot(SQ_wv,O_prob_state_SQ[2,*],name=O_nam[4],$
  thick=thck,color=colors[0],symbol=sym[1],sym_thick=sym_thck,/overplot)

leg5 = LEGEND(TARGET=[O_crss1,O_crss2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)



w6= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O_crss1=plot(NRL_wv,O_prob_state_nrl[3,*],name=O_nam[2],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;yrange=yr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O_crss2=plot(SQ_wv,O_prob_state_SQ[3,*],name=O_nam[5],$
  thick=thck,color=colors[0],symbol=sym[1],sym_thick=sym_thck,/overplot)

leg6 = LEGEND(TARGET=[O_crss1,O_crss2], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)

;__________________________________________________________________________________________________________________________________________________________________________________

;O  all states

O_nam=['NRL 4S!U0!N','NRL 2D!U0!N','NRL 2P!U0!N','NRL 4P','NRL 2P','NRL K']


xr=[0.1,20.]
w7= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O_crss1=plot(NRL_wv,O_prob_state_nrlall[1,*],name=O_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[0]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;yrange=yr,
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O_crss2=plot(NRL_wv,O_prob_state_nrlall[2,*],name=O_nam[1],$
  thick=thck,color=colors[1],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss3=plot(NRL_wv,O_prob_state_nrlall[3,*],name=O_nam[2],$
  thick=thck,color=colors[2],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss4=plot(NRL_wv,O_prob_state_nrlall[4,*],name=O_nam[3],$
  thick=thck,color=colors[4],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss5=plot(NRL_wv,O_prob_state_nrlall[5,*],name=O_nam[4],$
  thick=thck,color=colors[5],symbol=sym[2],sym_thick=sym_thck,/overplot)

O_crss6=plot(NRL_wv,O_prob_state_nrlall[6,*],name=O_nam[5],$
  thick=thck,color=colors[6],symbol=sym[2],sym_thick=sym_thck,/overplot)





leg7 = LEGEND(TARGET=[O_crss1,O_crss2,O_crss3,O_crss4,O_crss5,O_crss6], POSITION=leg_loc,font_size=fnt_sz,transparency=leg_transp,/normal)














;__________________________________________________________________________________________________________________________________________________________________________________



  yr=[-0.1,1.1]
;O2

w2= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

O2_crss1=plot(NRL_wv,O2_prob_state_NRL[1,*],name=O2_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[1]+titl,$
  xstyle=1,ystyle=1,xrange=xr,yrange=yr,xlog=0,ylog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

O2_crss2=plot(NRL_wv,O2_prob_state_NRL[2,*],name=O2_nam[1],$
  thick=thck,color=colors[1],symbol=sym[2],sym_thick=sym_thck,/overplot)

O2_crss3=plot(SQ_wv,O2_prob_state_SQ[1,*],name=O2_nam[2],$
  thick=thck,color=colors[2],symbol=sym[1],sym_thick=sym_thck,/overplot)

O2_crss4=plot(SQ_wv,O2_prob_state_SQ[2,*],name=O2_nam[3],$
  thick=thck,color=colors[3],symbol=sym[1],sym_thick=sym_thck,/overplot)

 O2_crss5=plot(nrL_wv,  total(O2_prob_state_NRLALL[1:3,*],1),name=O2_nam[4],$
   thick=thck,color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

 O2_crss6=plot(NRL_wv,total(O2_prob_state_NRLALL[4:-1,*],1),name=O2_nam[5],$
   thick=thck,color=colors[1],symbol=sym[0],sym_thick=sym_thck,/overplot)


leg2 = LEGEND(TARGET=[O2_crss1,O2_crss2,O2_crss3,O2_crss4,O2_crss5,O2_crss6], POSITION=leg_loc_left,font_size=fnt_sz,transparency=leg_transp,/normal)




;__________________________________________________________________________________________________________________________________________________________________________________

;xr=[0.1,20.]

;N2

w3= WINDOW(WINDOW_TITLE=titl,DIMENSIONS=get_screen_size())

N2_crss1=plot(NRL_wv,N2_prob_state_NRL[1,*],name=N2_nam[0],$
  xtitle=xt,ytitle=yt,title=maj[2]+titl,$
  xstyle=1,ystyle=1,xrange=xr,xlog=0,ylog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  thick=thck,color=colors[0],symbol=sym[2],sym_thick=sym_thck,font_size=fnt_sz,/current)

N2_crss2=plot(NRL_wv,N2_prob_state_NRL[2,*],name=N2_nam[1],$
  thick=thck,color=colors[1],symbol=sym[2],sym_thick=sym_thck,/overplot)

N2_crss3=plot(SQ_wv,N2_prob_state_SQ[1,*],name=N2_nam[2],$
  thick=thck,color=colors[0],symbol=sym[1],sym_thick=sym_thck,/overplot)

N2_crss4=plot(SQ_wv,N2_prob_state_SQ[2,*],name=N2_nam[3],$
  thick=thck,color=colors[1],symbol=sym[1],sym_thick=sym_thck,/overplot)

 N2_crss5=plot(nrl_wv,total(N2_prob_state_nrlall[1:3,*],1),name=N2_nam[4],$
   thick=thck,color=colors[0],symbol=sym[0],sym_thick=sym_thck,/overplot)

 N2_crss6=plot(nrl_wv,total(N2_prob_state_nrlall[4:-1,*],1),name=N2_nam[5],$
   thick=thck,color=colors[1],symbol=sym[0],sym_thick=sym_thck,/overplot)


leg3 = LEGEND(TARGET=[N2_crss1,N2_crss2,N2_crss3,N2_crss4,N2_crss5,N2_crss6], POSITION=leg_loc_left,font_size=fnt_sz,transparency=leg_transp,/normal)






end
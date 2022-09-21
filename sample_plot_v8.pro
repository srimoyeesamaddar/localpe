; Plotting contours of photoiionization for different solar zenith angle and altiudes
; code in sample_localpe_v8


floc = '/home/srimoyee/Desktop/nrl_files/sav_files/'

file_m5 = floc+"peak_m5_sza.sav"
file_SQm5 = floc+"SQ_peak_m5_sza.sav"

file_x9 =floc+"peak_x9_sza.sav"
file_SQx9 = floc +"SQ_peak_x9_sza.sav"


;file has :
;    photoi_wv_sza,photoi_sza,eiionz_sza,eiionz_transp_sza,sza_arr,lotlon_arr, zz, wv1,wv2, ssflux

restore, file_m5 
m5_photoi_wv = photoi_wv_sza
nrl_wv1= wv1
nrl_wv2= wv2
m5_sza= sza_arr

restore, file_SQm5 
SQm5_photoi_wv = photoi_wv_sza
SQ_wv1= wv1
SQ_wv2= wv2

restore, file_x9
x9_photoi_wv = photoi_wv_sza
x9_sza= sza_arr

restore, file_SQx9
SQx9_photoi_wv = photoi_wv_sza

close,/all
;__________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________

; m5 flare- 1st 12 bins of NRL and last 19 bins of SQ'05
tot_m5_photoi= total(m5_photoi_wv[*,*,0:11,*],3)+ total(SQm5_photoi_wv[*,*,3:21,*],3)

;x9 flare -  1st 12 bins of NRL and last 19 bins of SQ'05
tot_x9_photoi= total(x9_photoi_wv[*,*,0:11,*],3)+ total(SQx9_photoi_wv[*,*,3:21,*],3)




;__________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________

;Universal Plot variables
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
posleft=[0.05,0.25,0.45,0.9]
posright=[0.55,0.25,0.95,0.9]
postot= [0.15,0.25,0.9,0.9]

thck=5.
ln_thck=5
sym_thck=5
leg_loc=[0.92, 0.88]

leg_loc_left=[0.05, 0.9]
leg_loc_right=[0.9, 0.9]

maj=['O ','O!L2!N ','N!L2!N ']
leg_nam=['NRL (0.5-18 $\AA$)','SQ (0.5-18 $\AA$)']

  sym=["*","+","tu"]

ct=colortable(26)
poscb=[0.15,0.1,0.85,0.15]
poscbleft=[0.15,0.1,0.45,0.15]
poscbright=[0.55,0.1,0.85,0.15]
;__________________________________________________________________________________________________________________________


cb_titl='Photoionization (Pi)(cm!U-3!N s!U-1!N) '
yt='Altitude (km)'
xt= 'Solar Zenith Angle (degree)'
;xr=[min(m5_sza),max(m5_sza)]*!radeg
xr=[min([m5_sza,x9_sza])*!radeg,89.]
yr=[min(zz),210.]

; O photoiionisation


dr=[  min( [min( tot_m5_photoi[*,0,*]), min( tot_x9_photoi[*,0,*])] ) , $

      max( [max( tot_m5_photoi[*,0,*]), max( tot_x9_photoi[*,0,*])] ) ]
      
c_val= findgen(20)*150.+190.


w1 = WINDOW(WINDOW_TITLE=tilt,DIMENSIONS=get_screen_size())

titl='NRL M5 spectrum -'+ maj[0]+' Photoionisation'

cp11_o= contour(transpose(reform(tot_m5_photoi[*,0,*])),m5_sza*!radeg,zz,$
                rgb_table=ct,axis_style=2,$
                xtitle=xt,ytitle=yt,title=titl,$
                yrange=yr,xrange=xr,min_value=dr[0],max_value=dr[1],c_value=c_val,$  
                xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
                position=posleft,font_size=fnt_sz,c_thick=thck,/current)



titl='NRL X9 spectrum -'+ maj[0]+' Photoionisation'


cp21_o= contour(transpose(reform(tot_x9_photoi[*,0,*])),x9_sza*!radeg,zz,$
                rgb_table=ct,axis_style=2,$
                xtitle=xt,ytitle=yt,title=titl,$
                yrange=yr,xrange=xr,min_value=dr[0],max_value=dr[1],c_value=c_val,$
                xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
                position=posright,font_size=fnt_sz,c_thick=thck,/current)


cb_o= colorbar( rgb_table=ct, orientation = 0, position=poscb,font_size=fnt_sz,title=cb_titl,range=dr,tickvalues=c_val,/border)



;__________________________________________________________________________________________________________________________


; O2 photoiionisation


dr=[  min( [min( tot_m5_photoi[*,1,*]), min( tot_x9_photoi[*,1,*])] ) , $

  max( [max( tot_m5_photoi[*,1,*]), max( tot_x9_photoi[*,1,*])] ) ]

;c_val= findgen(28)*400.+500.
c_val= findgen(11)*1000.+400.

w2 = WINDOW(WINDOW_TITLE=tilt,DIMENSIONS=get_screen_size())

titl='NRL M5 spectrum -'+ maj[1]+' Photoionisation'

cp1_o2= contour(transpose(reform(tot_m5_photoi[*,1,*])),m5_sza*!radeg,zz,$
            rgb_table=ct,axis_style=2,$
            xtitle=xt,ytitle=yt,title=titl,$
            yrange=yr,xrange=xr,min_value=0.,max_value=dr[1],c_value=c_val,$
            xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
            position=posleft,font_size=fnt_sz,c_thick=thck,/current)



titl='NRL X9 spectrum -'+ maj[1]+' Photoionisation'

cp2_o2= contour(transpose(reform(tot_x9_photoi[*,1,*])),x9_sza*!radeg,zz,$
            rgb_table=ct,axis_style=2,$
            xtitle=xt,ytitle=yt,title=titl,$
            yrange=yr,xrange=xr,min_value=0.,max_value=dr[1],c_value=c_val,$;
            xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
            position=posright,font_size=fnt_sz,c_thick=thck,/current)

cb_o2= colorbar( rgb_table=ct, orientation = 0, position=poscb,font_size=fnt_sz,title=cb_titl,range=[0.,dr[1],tickvalues=c_val,/border)
;  
;;__________________________________________________________________________________________________________________________

; N2 photoiionisation

w3 = WINDOW(WINDOW_TITLE=tilt,DIMENSIONS=get_screen_size())
titl='NRL M5 spectrum -'+ maj[2]+' Photoionisation'

cp1_n2= contour(transpose(reform(tot_m5_photoi[*,2,*])),m5_sza*!radeg,zz,$
  rgb_table=ct,axis_style=2,$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  position=posleft,font_size=fnt_sz,c_thick=thck,/current)

cb1_n2= colorbar(target= cp1_n2, orientation = 0, position=poscbleft,font_size=fnt_sz,title=cb_titl,/border)

titl='NRL X9 spectrum -'+ maj[2]+' Photoionisation'

cp2_n2= contour(transpose(reform(tot_x9_photoi[*,2,*])),x9_sza*!radeg,zz,$
  rgb_table=ct,axis_style=2,$
  xtitle=xt,ytitle=yt,title=titl,$
  xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=0,$;
  xgridstyle=xg,ygridstyle=yg,xticklen=xtl,yticklen=ytl,xsubgridstyle=xsg,ysubgridstyle=ysg,xsubticklen=xstl,ysubticklen=ystl,$
  position=posright,font_size=fnt_sz,c_thick=thck,/current)

cb2_n2= colorbar(target= cp2_n2, orientation = 0, position=poscbright,font_size=fnt_sz,title=cb_titl,/border)













;__________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________


  end
pro sample_localpe_plots_srimoyee

;Comparing the photoionisation & photoelectrons for Fenn, Leiden & Leiden 1nm
 
file1='C:\Users\srimo\Desktop\nrl_files\sav_files\photoionz_fenncrss.sav'
;file has :
;photoi_fenncrss,eiionz_fenncrss,zz_fenncrss

file2='C:\Users\srimo\Desktop\nrl_files\sav_files\photoionz_leiden.sav'
;file has :
;photoi_leiden,eiionz_leiden,zz_leiden

file3='C:\Users\srimo\Desktop\nrl_files\sav_files\photoionz_leiden1nm.sav'
;file has :
;photoi_leiden1nm,eiionz_leiden1nm,zz_leiden1nm

restore, file1
restore, file2
restore, file3

;Plots for Photoionisation & photoelectrons for Fennely, Leiden and leiden 1nm cross-sections- Spectrum 2...........

;cgps_open, filename='C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_comp.ps'
;
;;#1
;cgdisplay,1000,1000, /free
;cgplot, photoi_fenncrss(0,*), zz_fenncrss ,color='red',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;  xtitle='O Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;  charthick=2, thick=3
;cgplot, photoi_leiden(0,*), zz_leiden,color='magenta',/overplot,thick=3,psym=-46
;cgplot, photoi_leiden1nm(0,*), zz_leiden1nm,color='orange',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['red','magenta','orange'],$
;  alignment=1,/box,Location=[0.9, 0.9]
; 
;
;;#2
;cgdisplay,1000,1000, /free
;cgplot, photoi_fenncrss(1,*), zz_fenncrss ,color='violet',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;    xtitle='O$\down2$ Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;    charthick=2, thick=3
;cgplot, photoi_leiden(1,*), zz_leiden,color='dodger blue',/overplot,thick=3
;cgplot, photoi_leiden1nm(1,*), zz_leiden1nm,color='purple',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,15,15], symcolors=['violet','dodger blue','purple'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;;#3
;cgdisplay,1000,1000, /free
;cgplot, photoi_fenncrss(2,*), zz_fenncrss ,color='sienna',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;    xtitle='N$\down2$ Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;    charthick=2, thick=3
;cgplot, photoi_leiden(2,*), zz_leiden,color='sandy brown',/overplot,thick=3, psym=-46
;cgplot, photoi_leiden1nm(2,*), zz_leiden1nm,color='coral',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['sienna','sandy brown','coral'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;
;
;;#4
;cgdisplay,1000,1000, /free
;cgplot, eiionz_fenncrss(0,*), zz_fenncrss ,color='red',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;    xtitle='O Photoelectron flux (photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;    charthick=2, thick=3
;cgplot, eiionz_leiden(0,*), zz_leiden,color='magenta',/overplot,thick=3,psym=-46
;cgplot, eiionz_leiden1nm(0,*), zz_leiden1nm,color='orange',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['red','magenta','orange'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;
;;#5
;cgdisplay,1000,1000, /free
;cgplot, eiionz_fenncrss(1,*), zz_fenncrss ,color='violet',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;    xtitle='O$\down2$ Photoelectron flux (photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;    charthick=2, thick=3
;cgplot, eiionz_leiden(1,*), zz_leiden,color='dodger blue',/overplot,thick=3, psym=-46
;cgplot, eiionz_leiden1nm(1,*), zz_leiden1nm,color='purple',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['violet','dodger blue','purple'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;;#6
;cgdisplay,1000,1000, /free
;cgplot, eiionz_fenncrss(2,*), zz_fenncrss ,color='sienna',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
;    xtitle='N$\down2$ Photoielectron flux (photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
;    charthick=2, thick=3
;cgplot, eiionz_leiden(2,*), zz_leiden,color='sandy brown',/overplot,thick=3,psy=-46
;cgplot, eiionz_leiden1nm(2,*), zz_leiden1nm,color='coral',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['sienna','sandy brown','coral'],$
;    alignment=1,/box,Location=[0.9, 0.9]
;
;
;cgps2pdf,'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_comp.ps',$
;    'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_comp.pdf', gs_path='C:\Program Files\gs\gs9.27'
;cgps_close,/pdf
;



;..........................................................................................
;Pe/Pi for Fennely, Leiden and Leiden 1nm cross sections -Spectrum2

;cgps_open, filename='C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_ratio.ps'

;;#7
;cgdisplay,1000,1000, /free
;cgplot,eiionz_fenncrss(0,*)/ photoi_fenncrss(0,*), zz_fenncrss ,color='red',/xs,/ys,/xl,yr=[40.,300.],xr=[0.1,1.e3],$
;  xtitle='O Photoelectron/Photoionization',ytitle='Altitude (km)',$
;  charthick=2, thick=3
;cgplot,eiionz_leiden(0,*)/ photoi_leiden(0,*), zz_leiden,color='magenta',/overplot,thick=3,psym=-46
;cgplot, eiionz_leiden1nm(0,*)/photoi_leiden1nm(0,*), zz_leiden1nm,color='orange',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['red','magenta','orange'],$
;  alignment=1,/box,Location=[0.9, 0.9]
;
;
;;#8
;cgdisplay,1000,1000, /free
;cgplot, eiionz_fenncrss(1,*)/ photoi_fenncrss(1,*), zz_fenncrss ,color='violet',/xs,/ys,/xl,yr=[40.,300.],xr=[0.1,1.e3],$
;  xtitle='O$\down2$ Photoelectron/Photoionization',ytitle='Altitude (km)',$
;  charthick=2, thick=3
;cgplot, eiionz_leiden(1,*)/ photoi_leiden(1,*),zz_leiden,color='dodger blue',/overplot,thick=3
;cgplot, eiionz_leiden1nm(1,*)/photoi_leiden1nm(1,*), zz_leiden1nm,color='purple',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,15,15], symcolors=['violet','dodger blue','purple'],$
;  alignment=1,/box,Location=[0.9, 0.9]
;
;;#9
;cgdisplay,1000,1000, /free
;cgplot, eiionz_fenncrss(2,*)/ photoi_fenncrss(2,*), zz_fenncrss ,color='sienna',/xs,/ys,/xl,yr=[40.,300.],xr=[0.1,1.e3],$
;  xtitle='N$\down2$ Photoelectron/Photoionization',ytitle='Altitude (km)',$
;  charthick=2, thick=3
;cgplot, eiionz_leiden(2,*)/ photoi_leiden(2,*), zz_leiden,color='sandy brown',/overplot,thick=3, psym=-46
;cgplot, eiionz_leiden1nm(2,*)/photoi_leiden1nm(2,*), zz_leiden1nm,color='coral',/overplot,thick=3
;cglegend, titles=['Fennelly','Leiden','Leiden 1nm'],psyms=[15,-46,15], symcolors=['sienna','sandy brown','coral'],$
;  alignment=1,/box,Location=[0.9, 0.9]
;

;cgps2pdf,'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_ratio.ps',$
;  'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photo_eiion_ratio.pdf', gs_path='C:\Program Files\gs\gs9.27'
;cgps_close,/pdf

;...............................................................................

; Comparing Photoionisation for two different spectra Spectrum 2 and Spectrum 5


file4='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\photoionz_fenncrss_sp5.sav'
;file has :
;photoi_fenncrss_sp5,eiionz_fenncrss_sp5,zz_fenncrss_sp5

file5='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\photoionz_leiden_sp5.sav'
;file has :
;photoi_leiden_sp5,eiionz_leiden_sp5,zz_leiden_sp5

file6='C:\Users\Srimoyee\Desktop\nrl_files\sav_files\photoionz_leiden1nm_sp5.sav'
;file has :
;photoi_leiden1nm_sp5,eiionz_leiden1nm_sp5,zz_leiden1nm_sp5

restore, file4
restore, file5
restore, file6

;#10 Fennely Cross-sections
;Photoionisation O

;cgps_open, filename='C:\Users\Srimoyee\Desktop\nrl_files\plots2\photoeiion_speccomp_fenn.ps'

cgdisplay,1000,1000, /free
cgplot, photoi_fenncrss(0,*), zz_fenncrss ,color='red',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
  xtitle='O Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
  charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
cgplot, photoi_fenncrss_sp5(0,*), zz_fenncrss_sp5,color='magenta',/overplot,thick=3,psym=-46
cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['red','magenta'],$
  alignment=1,/box,Location=[0.9, 0.9]

cgplot,photoi_fenncrss(0,*)/photoi_fenncrss_sp5(0,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$ ;
       xtitle='O Photoionization ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
       charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'
;O2
cgdisplay,1000,1000, /free
cgplot, photoi_fenncrss(1,*), zz_fenncrss ,color='violet',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
    xtitle='O$\down2$ Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
cgplot, photoi_fenncrss_sp5(1,*), zz_fenncrss_sp5,color='dodger blue',/overplot,thick=3,psym=-46
cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['violet','dodger blue'],$
    alignment=1,/box,Location=[0.9, 0.9]

cgplot,photoi_fenncrss(1,*)/photoi_fenncrss_sp5(1,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$
    xtitle='O2 Photoionization ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'

;N2
cgdisplay,1000,1000, /free
cgplot, photoi_fenncrss(2,*), zz_fenncrss ,color='sienna',/xs,/ys,yr=[40.,300.],$
    xtitle='N$\down2$ Photoionization(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
cgplot, photoi_fenncrss_sp5(2,*), zz_fenncrss_sp5,color='sandy brown',/overplot,thick=3,psym=-46
cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['sienna','sandy brown'],$
    alignment=1,/box,Location=[0.9, 0.9]

cgplot,photoi_fenncrss(2,*)/photoi_fenncrss_sp5(2,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$
    xtitle='N2 Photoionization ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'

;Photoelectron O
cgdisplay,1000,1000, /free
  cgplot, eiionz_fenncrss(0,*), zz_fenncrss ,color='red',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
    xtitle='O Photoelectron flux(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
  cgplot, eiionz_fenncrss_sp5(0,*), zz_fenncrss_sp5,color='magenta',/overplot,thick=3,psym=-46
  cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['red','magenta'],$
    alignment=1,/box,Location=[0.9, 0.9]
    
cgplot,eiionz_fenncrss(0,*)/eiionz_fenncrss_sp5(0,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$
    xtitle='O Photoelectron ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'

; O2
cgdisplay,1000,1000, /free
cgplot, eiionz_fenncrss(1,*), zz_fenncrss ,color='violet',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
    xtitle='O$\down2$ Photoelectron flux(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
cgplot, eiionz_fenncrss_sp5(1,*), zz_fenncrss_sp5,color='dodger blue',/overplot,thick=3,psym=-46
cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['violet','dodger blue'],$
    alignment=1,/box,Location=[0.9, 0.9]

cgplot,eiionz_fenncrss(1,*)/eiionz_fenncrss_sp5(1,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$
    xtitle='O2 Photoelectron ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'


;N2
cgdisplay,1000,1000, /free
cgplot, eiionz_fenncrss(2,*), zz_fenncrss ,color='sienna',/xs,/ys,yr=[40.,300.],/xl,xr=[1.,2.e4],$
    xtitle='N$\down2$ Photoelectron flux(photons cm$\up-3$ s$\up-1$)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.55,0.9,0.95]
cgplot, eiionz_fenncrss_sp5(2,*), zz_fenncrss_sp5,color='sandy brown',/overplot,thick=3,psym=-46
cglegend, titles=['Spectrum 2','Spectrum 5'],psyms=[15,-46], symcolors=['sienna','sandy brown'],$
    alignment=1,/box,Location=[0.9, 0.9]

cgplot,eiionz_fenncrss(1,*)/eiionz_fenncrss_sp5(1,*), zz_fenncrss,/xs,/ys,yr=[40.,300.],/xl,xr=[0.1,1.e4],$
    xtitle='N2 Photoelectron ratio (spec 2/spec 5)',ytitle='Altitude (km)',$
    charthick=2, thick=3,position=[0.18,0.1,0.9,0.43],/noerase,color='deep pink'

;cgps2pdf,'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photoeiion_speccomp_fenn.ps',$
;    'C:\Users\Srimoyee\Desktop\nrl_files\plots2\photoeiion_speccomp_fenn.pdf', gs_path='C:\Program Files\gs\gs9.27'
;cgps_close,/pdf

restore,'C:\Users\Srimoyee\Desktop\july23.sav'
restore,'C:\Users\Srimoyee\Desktop\nov04.sav'


cgps_open, filename='C:\Users\Srimoyee\Desktop\nrl_spec.ps'

cgdisplay,1000,1000, /free
cgplot,(wv11+wv22)/2., ssflux1, thick=5,xtitle='!C !C Wavlength ($\Angstrom$)',$
       ytitle='Solar Flux (Photons cm$\up-2$ s$\up-1 $ $\Angstrom$ $\up-1 $)', /yl,/ys,/xl,/xs,psym=10,$
       color='red',charthick=2, position=[0.18,0.1,0.95,0.85],xr=[1.,1000.],yr=[1.e6,1.e12]

cgplot,(wv1+wv2)/2., ssflux2,/overplot,thick=5,color='dodger blue',psym=10

cglegend, titles=['Spectrum1 (July 23, 2016)','Spectrum2 (Nov 4, 2003)'],$
          psyms=[15,15], symcolors=['red','dodger blue'],$
          alignment=1,/box,Location=[0.95, 0.85]


cgps2pdf,'C:\Users\Srimoyee\Desktop\nrl_spec.ps',$
        'C:\Users\Srimoyee\Desktop\nrl_spec.pdf', gs_path='C:\Program Files\gs\gs9.27'
cgps_close,/pdf






close, /all

end
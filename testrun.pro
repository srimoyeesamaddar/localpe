  
  
  restore,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\mjuly.sav'
  
  
  eiionz_transp1=eiionz_transp
  thck=5
 photoi1=photoi
  wv11=wv1
  wv21=wv2

ssflux1=ssflux
  
  restore,'C:\Users\Srimoyee\Desktop\nrl_files\sav_files\X28_allrun.sav'
  
  
  eiionz_transp2=eiionz_transp
  
  photoi2=photoi
  ssflux2=ssflux
  wv12=wv1
  wv22=wv2
  
  xr=[0.01,1000.]
  yr=[min(zz),max(zz)]
  leg_loc=[0.9,0.85]
  xt='Pe/Pi '
  yt='Altitude (km)'
  tilt='Local Pe/Pi'
  
  w1 = WINDOW(WINDOW_TITLE=tilt, $
      DIMENSIONS=[800,800])
    
    
    pp1=plot(eiionz_transp1[0,*]/photoi1[0,*],zz,name='X28 Nov 4, 2003',$
      xtitle='Pe/Pi for O',ytitle='Altitude (km)',$
      xstyle=1,ystyle=1,yrange=[70.,max(zz)],xrange=xr,xlog=1,$;
      color='red',thick=thck,/current)
    
    pp2=plot(eiionz_transp2[0,*]/photoi2[0,*],zz,name='M5 July 23, 2016',$
             thick=thck,color='blue',/overplot)
             
        sav_loc='C:\Users\Srimoyee\Desktop\agu_plots\'     
           l1=legend(target=[pp1,pp2],position=[0.9,0.85],/normal)

           w1.Save, sav_loc+"test1.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

           
           
           w2 = WINDOW(WINDOW_TITLE=tilt, $
             DIMENSIONS=[800,800])
    
    
    pp3=plot(eiionz_transp1[1,*]/photoi1[1,*],zz,name='X28 Nov 4, 2003',$
      xtitle='Pe/Pi for O!L2!N',ytitle='Altitude (km)',$
      xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
          color='red',thick=thck,/current)
    
    pp4=plot(eiionz_transp2[1,*]/photoi2[1,*],zz,name='M5 July 23, 2016',$
      color='blue',thick=thck,/overplot)
    
     l2=legend(target=[pp3,pp4],position=[0.9,0.85],/normal)

           w2.Save, sav_loc+"test2.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

    w3 = WINDOW(WINDOW_TITLE=tilt, $
      DIMENSIONS=[800,800])
    
    pp5=plot(eiionz_transp1[2,*]/photoi1[2,*],zz,name='X28 Nov 4, 2003',$
      xtitle='Pe/Pi for N!L2!N',ytitle='Altitude (km)',$
      xstyle=1,ystyle=1,yrange=yr,xrange=xr,xlog=1,$;
      color='red',thick=thck,/current)
    
    pp6=plot(eiionz_transp2[2,*]/photoi2[2,*],zz,name='M5 July 23, 2016',$
      thick=thck,color='blue',/overplot)
      
      l3=legend(target=[pp5,pp6],position=[0.9,0.85],/normal)

    w3.Save, sav_loc+"test3.jpeg", BORDER=10, RESOLUTION=600, /TRANSPARENT

    
     
    
 
    
    
    
    w4=window(window_title='Solar Flux',$
        dimension=[800,800])
      
      s1=plot(wv11, ssflux1, $
        ylog=1,ystyle=1,xlog=1,xstyle=1,$
        xrange=[0.5,1.75e3],yrange=[1.e1,1.e12],$
        thick=thck,histogram=1,color='red',/current,$ ;, position=[0.18,0.1,0.95,0.85]
        xtitle='Wavlength ('+'!3' + STRING(197B) + '!X'+')' ,$
        ytitle='Solar Flux (Photons cm!U-2!N s!U-1!N)',$
        name='X28 Nov 4, 2003')
      
      s2=plot(wv12,ssflux2,$
               thick=thck,histogram=1,color='blue',$
               name='M5 July 23, 2016',/overplot)
      
      ll1=legend(target=[s1,s2],position=[0.55,0.85],/normal)
      
      w4.Save, sav_loc+"test4.jpeg", BORDER=20, RESOLUTION=600, /TRANSPARENT
      
    
    end
    
  
  
 
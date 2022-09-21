;Using the Chapman function fit [Smith and Smith [JGR 77, 3592, 1972]

Function chapman ,chi, z, t, i

  ;     Universal Constants
 

  am=[16., 32., 28.] ; amu for O, O2 anf N2
  kb= 1.38064852e-16   ; Boltzmann's constant ergs per k

  re=6378.e5 ; radius of the Earth in cm
  g=980.  ;acceleration due to gravity in cm^2/sec
  gr=g*(re/(re+z))^2 ; account for small change in g with altitude

  amave=(.2*32 + .79*28 + .01*40.) ; average mass of atmosphere, used for lower altitudes
  gperamu=1.66054*1.e-24  ;coverting amu to grams


  ;    Assumption: at low altitudes, every constituent follows
  ;                a single scale height, above that, they follow
  ;                their own, assume change occurs at 100km


  ilow=z le 110.*1e5
;  ilow=z le 115.*1e5


  ihi=1-ilow
  amass=amave*ilow + am[i]*ihi     ;mass for every height

  H=kb*t/(amass*gperamu*gr) ; scale height  CHANGED NAME FROM hn TO H
  
;  H[ilow]=7.e5   ;GIVING CONSTANT SCALE HEIGHT TO ZZ<= 115KM
  Xp=(re+z)/H              ;CHANGED NAME FROM hg TO Xp
  y=sqrt(0.5*Xp*(cos(chi)^2))    ;CHANGED NAME FROM hf TO y NOTE: chi should be in radians


  ;    Constants for calculating error function as given in paper

  a=1.0606963
  b=0.55643831
  c=1.0619896
  d=1.7245609
  e=0.56498823
  f=0.06651874

  erfcy=fltarr(n_elements(y))

  ;    Calculating the fractional part of the error function

  for ii=0, n_elements(y)-1 do begin

    if ((y[ii] ge 0.) && (y[ii] lt 8.)) then $
      erfcy[ii] = (a+b*y[ii]) / (c+d*y[ii]+y[ii]^2) else $
      if ((y[ii] ge 8.) && (y[ii] le 100.)) then $
      erfcy[ii]=e/(f+y[ii]) $
    else begin
      print,'y out of range in error function calculation in ace_rcolum code'
      stop
    endelse

  endfor

  ;    Calculating the error function

  chap=sqrt(0.5*!pi*Xp)*erfcy

  return,chap
end
;____________________________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________________________

;Vertical Column Density

Pro vcd,zz_cm,zmaj,jmax,nmaj,zvcd    
  
  zvcd=fltarr(nmaj,jmax) ; array of vertical column densities

  ; the following integrates number density down in altitude
  ; and includes a factor to take into account that the atmosphere
  ; varies exponentially within any bin

  for i=0,nmaj-1 do begin
    zvcd[i,jmax-1] =   zmaj[i,jmax-1] $
      * (zz_cm[jmax-1]-zz_cm[jmax-2]) $
      / alog(zmaj[i,jmax-2]/zmaj[i,jmax-1])
    for j=jmax-2,0,-1 do begin
      ;                  if zmaj[i,j] ne 0. then begin  ;added this check: Srimoyee
      rat = zmaj[i,j+1] / zmaj[i,j]
      zvcd[i,j] =   zvcd[i,j+1] $
        + zmaj[i,j] * (zz_cm[j]-zz_cm[j+1]) / alog(rat) * (1.-rat)
      ;                  endif
    endfor
  endfor
  ; adding the column densities for lower altitudes to be constant when the densities at lower altitudes become zero
  for i=0, nmaj-1 do begin
    zvcd_ind=where(~finite(reform(zvcd[i,*])) , ct_i)

    if ct_i gt 0 then begin  ;infinte column density
      min_ind=min(where(finite(reform(zvcd[i,*])) ,/null),min_i)
      zvcd[i,zvcd_ind]=zvcd[i,min_ind]
    endif

  endfor

  return
end

;____________________________________________________________________________________________________________________________________________________________________________________
;____________________________________________________________________________________________________________________________________________________________________________________
 
            
; Stan Solomon, 1988, 1991
; Modified for IDL, 12/93 Scott Bailey
; Updated for !a environment and allowed for array processing, SMB 1/2013
;
; Calculates the column density ZCOL for each species ZMAJ above height
; ZZ at zenith angle CHI.  Uses a Chapman function fit [Smith and Smith,
; JGR 77, 3592, 1972].  If CHI is less than 90 degrees, column
; densities are calculated directly; if CHI is greater than 90 degrees
; the column density at grazing height for 90 degrees is calculated and
; doubled and the column density above ZZ(J) is subtracted; if CHI is
; such that the grazing height is less than the radius of the earth the
; column densities are set to 'infinity', i.e., 1.0E30.  Densities
; supplied in array ZMAJ are used in the calculation except where
; grazing height is below the lowest level specified, in this case
; values are interpolated logarithmically from a US Standard Atmosphere
; at sea level, the tropopause, the stratopause, and the mesopause.

      Pro ace_rcolum ,chi, z, zmaj, tn, zcol, zvcd, jmax, nmaj

      zz_cm=z*1e5 ; (calcultaions are done in cm.)
      zcol=fltarr(nmaj,jmax)
      zvcd=fltarr(nmaj,jmax)
      nm=3
      nu=4
      zcus=fltarr(nm,nu)
      zcg=fltarr(nmaj)
      re=6.37e8 ; cm
       ; us standard atmosphere
      zus=[0., 1.5e6, 5.e6, 9.e6]& tnus=[288., 217., 271., 187.]
      zcus=zcus+[   8.00e17, 4.54e24, 1.69e25, $
                    8.00e17, 5.46e23, 2.03e24, $
                    8.00e17, 3.63e21, 1.35e22, $
                    7.80e17, 8.48e18, 3.16e19]

;____________________________________________________________________________________________________________________________________

; start by getting vertical column densities
      vcd ,zz_cm, zmaj, jmax,nmaj,zvcd
      
;____________________________________________________________________________________________________________________________________

;       #1
   if chi le !pi/2.0 then begin

     for i=0,nmaj-1 do zcol[i,*]=zvcd[i,*]*chapman(chi,zz_cm,tn,i)

   endif else $


 if ((chi gt !pi/2.0) and (chi lt 2.)) then begin  ;#2

            for j=0,jmax-1 do begin  ;altitude loop
                grflag =0  ; this is a flag that tracks the grazing conditon; executes the Chapman function only when flag=1
                ghrg=(re+zz_cm(j))*sin(chi)
                ghz=ghrg-re

               if (ghz le 0.) then begin
                            zcol(*,j) = 1.0e30
;                            print,'enetered loop1'
               endif   else $  ;begin $


              if (ghz ge zz_cm(0)) then begin     ; grazing >zz[0] ;***
                     grflag =1
;                     print,'entered loop2'
                      for jjjg=0,j-2 do begin
                             if (zz_cm(jjjg) le ghz and zz_cm(jjjg+1) gt ghz) then jg=jjjg
                      endfor
                      tng = tn(jg)+(tn(jg+1)-tn(jg))*(ghz-zz_cm(jg))/(zz_cm(jg+1)-zz_cm(jg))
                      for i=0,nmaj-1 do $
                            zcg(i) = zvcd(i,jg) * (zvcd(i,jg+1) / zvcd(i,jg)) ^ $
                                    ((ghz-zz_cm(jg)) / (zz_cm(jg+1)-zz_cm(jg)))
              endif else begin ;grazing< zz[0]
                      grflag =1                
;                      print,'entered loop3'
                      for jjjg=0,2 do begin
                            if (zus(jjjg) lt ghz and zus(jjjg+1) gt ghz) then jg=jjjg
                      endfor
                      tng = tnus(jg) + (tnus(jg+1)-tnus(jg))*(ghz-zus(jg))/(zus(jg+1)-zus(jg))
                      for i=0,nmaj-1 do $
                            zcg(i) = zcus(i,jg) * (zcus(i,jg+1) / zcus(i,jg)) ^ $
                              ((ghz-zus(jg)) / (zus(jg+1)-zus(jg)))
         
              endelse ;***
       

;                          grflag=1  ;grazing >0 **** Not using the US standard atmoshere
;                           tng= INTERPOL(tn,zz, ghz)  ;taking the altitude ranges from zz[0] to zz[j+4] for original curve
;                           for i=0, nmaj-1 do   zcg[i]= INTERPOL(reform(alog(zvcd[i,*])),zz, ghz)
;                           
;       
;                      endelse ;****

       if (grflag eq 1) then begin
         for i=0,nmaj-1 do begin
             zcol(i,j) = 2. * zcg(i) * chapman(!pi/2.,ghz,tng,i) $
                         - zvcd(i,j) * chapman(chi,zz_cm(j),tn(j),i)
         endfor
      endif

     endfor ; altitude loop;#2
     
     
   endif else begin ; (chi ge 2.)  #3
     zcol=zcol+1.0e30

   endelse

         
;_____________________________________________________________________________________________________________________________________________
      
return
End

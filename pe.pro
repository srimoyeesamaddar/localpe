;cd,'/home/padma/togo/Rocket programs/'
;.run pe.pro

    ;C PROBO, PROBO2, PROBN2; branching ratio data arrays
    filename='read_ephotn2.sav'  ;from read_ephotn2.pro
    restore,file=filename ; gives aa,bb,probn2,sigin2,sigan2
    ;probn2(0,*)=x: branching ratios
    ;probn2(1,*)=an2
    ;probn2(2,*)=bn2
    ;probn2(3,*)=cn2
    ;probn2(4,*)=fn2
    ;probn2(5,*)=dissn2

    filename='read_ephotoO2.sav'  ;from read_ephotoO2.pro
    restore,file=filename  ; gives aa,bb,probo2,sigio2,sigao2

    ;probo2(0,*)=x: branching ratios
    ;probo2(1,*)=a+A
    ;probo2(2,*)=b
    ;probo2(3,*)=diss
    ;probo2(4:5,*)=all zeroes

    filename='read_ephotoO.sav'  ;from read_ephotn2.pro
    restore,file=filename ; gives aa,bb,probo,sigio,sigao

    ;probo(0,*)=branching ratios
    ;probo(1,*)=4s
    ;probo(2,*)=2Do
    ;probo(3,*)=2Po
    ;probo(4,*)=4Pe
    ;probo(5,*)=2Pe

    ;lmax=123
    ;JMAX=120


    ;glow.h
    ;C SIGABS  photoabsorption cross sections, O, O2, N2; cm2
    ;C SIGIONx  photoionization cross sections, O, O2, N2; cm2
    sigabs=dblarr(3,123) ;0 for O, 1 for O2, 2 for N2, 123 for wavelengths
    sigionx=dblarr(3,123) ;making this x so that it is not the same as the function sigion
    for l=0,122 do begin
    sigabs(0,*)=sigao(*)*1e-18
    sigabs(1,*)=sigao2(*)*1e-18
    sigabs(2,*)=sigan2(*)*1e-18
    sigionx(0,*)=sigio(*)*1e-18
    sigionx(1,*)=sigio2(*)*1e-18
    sigionx(2,*)=sigin2(*)*1e-18
    endfor
    
    nst=6                       
                                ; NST     number of states produced by photoionization/dissociation
                                ;TPOT    ionization potentials for each species, state; eV
    tpot=[[13.61, 16.93, 18.63, 28.50, 40.00,  0.00],$
     [12.07, 16.10, 18.20, 20.00,  0.00,  0.00],$
     [15.60, 16.70, 18.80, 30.00, 34.80, 25.00]]
  tpot=transpose(tpot)

                                ;C PROB    branching ratios for each state, species, and wavelength bin:
                                ;C         O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
                                ;C         O2+ states: X, a+A, b, dissoc.
                                ;C         N2+ states: X, A, B, C, F, dissoc.

    prob=dblarr(6,3,123) ;states species and energies
    l=0
    for l=0,122 do begin
     for k=0,5 do begin
      prob(k,0,l)=probo(k,l)
      prob(k,1,l)=probo2(k,l) 
      prob(k,2,l)=probn2(k,l)
     endfor
    endfor


 ; solar irradiance file

idate=2003351.
utsec=3600.*12.
sun=get_index()
dum=min(abs(sun.tyd-idate),index_idate)
if dum ne 0. then begin
   print,"You're screwed"
stop
endif
f107=sun.ten7(index_idate)
smf107=smooth(sun.ten7,81)
f107a=smf107(index_idate)

; solar irradiance file
;iscale = 1 ; EUVAC
;if keyword_set(see) or keyword_set(ut_see) then begin ; SEE data
;   iscale=2
;   if keyword_set(see) then see_euvac_glow,idate,f107,f107a
   ;if keyword_set(see_ut) then see_euvac_glow,ida
   ; if keyword_set(ut_see) then see_euvac_glow,ida
;    ;if keyword_set(ut_see) then see_euvac_glow,idate,f107,f107a,ut=ut
;   print,'We are now getting the solar spectrum...'

; solar irradiance file
;qte,f107,f107a,ut=ut
;endif

see_euvac_glow,idate,f107,f107a
restore,'exsect_cs.sav'

openr,lun,'ssflux_user.dat',/get_lun
s='a string'
readf,lun,s
lmax=123                        ; number of wavelength intervals for solar flux
wv1=fltarr(lmax)
wv2=fltarr(lmax)
ssflux=fltarr(lmax)
for i=0,lmax-1 do begin
   readf,lun,a,b,c
   wv1(i)=a
   wv2(i)=b
   ssflux(i)=c
endfor
free_lun,lun

;ssflux[33]=ssflux[33]*50.;;;;;;;;;;;;;;;;;;;;;;;;;fix this;;;;;;;;;;;;;;;;;;;;;
;idate=2003070.;changed apr4, 2011, jy;NOTE: must use date from 1999 if
;comparison with euvac spectrum, via plot_euvac_pespec.pro, is desired
;idate=1999340

;lat=30.;changed apr 4 2011, jy
lat=0.
lon=0.

dz=2.
zz=findgen(jmax)*dz + 100. ;altitude bin in km, jmax:number of altitude levels
z=zz*1e5                   ;altitude bin in cm
nrlmsis,  idate, utsec, lat, lon, zz, t_alt, t_exo, nd
zmaj=fltarr(3,jmax)
zmaj(0,*)=nd.o                  ;number density of o
zmaj(1,*)=nd.o2                 ;number density of o2
zmaj(2,*)=nd.n2                 ;number density of n2
zt=t_alt
save,file='pevar.sav',/var

END

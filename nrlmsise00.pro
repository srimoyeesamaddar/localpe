pro nrlmsise00, date, ut, glat, glon, ft7, ft7a, fap, alt, t_alt, t_exo, nd, md, SI=si_yesno, DRAG=drag_yesno

;Attempts to get Picone, Hedin, and Drob's, nrlmsise-00 into
;Linux idl.  Trying to replicate B. Knapp's similar feat for msis90.
;Does not allow for including the 7-array 3-day average aps, only
;daily.  Much code copied straight from Knapp's interface.  Thanks.

;Zachary Mullen
;Summer, 2004.  HAO, UCAR


display = 'more '
image = 'nrlmsise00.so'
ext = '_pass__'

;print usage stuff
  if n_params() eq 0  then begin
     print," "
     print," nrlmsise00.pro is an IDL interface to the Fortran code of"
     print," Al Hedin, Mike Picone, and Doug Drob's two-dimensional"
     print," thermospheric model , which provides estimates of temperature"
     print," and of number density of seven prevalant species (H, He, N,"
     print," O, N2, O2, and Ar).  It uses the 10.7 cm radio flux (for the"
     print,"  previous day, and a 81-day smooth centered on the requested"
     print," day) and the Ap geomagnetic indices as input."
     print," "
     print," Much of this program is copied and adapted from Barry Knapp's"
     print,"  similar interface for MSIS90."
     print," "
     print," Usage:"
     print," "
     print,"     nelmsise00, date, ut, glat, glon, ft7, ft7a, fap, alt, "
     print,"             t_alt, t_exo, nd, md, /SI, /DRAG"
     print," "
     return
 endif


;check inputs array sizes, etc.
ft7=FLOAT(ft7)
ft7a=FLOAT(ft7a)
fap=FLOAT(fap)
glat=FLOAT(glat)
glon=FLOAT(glon)

; Convert inputs to correct types and ranges (???)
  ldate = y2toy4( date )
  yyddd = long( ldate mod 100000L )
  fsec = float( ut )
  flat = float( glat )
  flon = ( (float( glon )+180.) mod 360. ) - 180.
  if flon lt -180. then $
     flon = flon + 360. $
  else if flon ge  180. then $
     flon = flon - 360.


szalt = size ( alt) 
nz = n_elements(szalt)
if szalt[0] eq 0 then falt = [ float( alt ) ] else falt = float( alt )

; Are we computing only temperatures, or also number densities?
  if n_params() lt 8 then req_code = 0L else req_code = 48L

;switch stuff
sw = fltarr(25)
for i=0, 24 do begin & sw(i)=1 & endfor
;print, sw
entry = 'tretrv'+ext
dummy = call_external( image, entry, sw )
;print, sw

for i=0, 24 do begin & sw(i)=1 & endfor
;reset switch stuff
entry = 'tselec'+ext
dummy = call_external( image, entry, sw) 
;print, sw

;local solar time
stl = ut/3600 + flon/15

;form of outputs?
  if szalt[0] eq 0 then nalt=1 else nalt=szalt[1]
  md = fltarr( nalt )
  t_alt = fltarr( nalt )
  t_exo = 0.
  nd_type = { _nd_, h:0., he:0., n:0., o:0., n2:0., o2:0., ar:0. }
  nd = replicate( nd_type, nalt )


;Prepare gtd7 variables
d = fltarr(9)
t = fltarr(2)


;call si or cgs?
if keyword_set( si_yesno) then si=-1L else si=0L
 entry = 'meters'+ext
 dummy=call_external( image, entry, si )


;Call gtd7
 entry= 'gtd7'+ext
 for j = 0,nalt-1 do begin
 dummy= call_external( image, entry, $
                       date, ut, falt(j), flat, flon, stl, ft7a, ft7, $
                       fap, req_code, d, t, /f_value)

;     dnew=reform(d)
     t_alt[j] = t[1]
     nd[j].h  = d[6]
     nd[j].he = d[0]
     nd[j].n  = d[7]
     nd[j].o  = d[1]
     nd[j].n2 = d[2]
     nd[j].o2 = d[3]
     nd[j].ar = d[4]
     md[j]    = d[5]
;     print, d
    
;call gtd7d, if desired, for inclusion of anomalous oxygen for drag
    if keyword_set(drag_yesno) then begin
    entry= 'gtd7d'+ext
    dummy=call_external( image, entry, $
                       date, ut, falt(j), flat, flon, stl, ft7a, ft7, $
                       fap, req_code, d, t)
    md[j]     = d[5]
 
    endif
    endfor
t_exo=t[0]
  if szalt[0] eq 0 then begin
     t_alt = t_alt[0]
     md = md[0]
  endif

;output tests
;plot, nd.h, alt, /XLOG, XRANGE=[1e5, 1e25], YRANGE=[0, 800]
;oplot, nd.he, alt
;oplot, nd.n, alt
;oplot, nd.o, alt
;oplot, nd.n2, alt
;oplot, nd.o2, alt
;oplot, nd.ar, alt
;plot, md, alt, /XLOG

  return

  end

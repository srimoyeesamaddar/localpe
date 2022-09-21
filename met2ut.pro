;***************************************************************************************
;+
; NAME: met2ut.pro
;
; LOCATION IN LIBRARY:
;
; TYPE:
;
; PURPOSE:
;   Converts MET values from missions including ARGOS, STS-39 (UVLIM)
;   mission into a data structure with numerous useful UT time information.
;
; Special Note:  UVLIM MET is generally stated in deciseconds!
;
; CATEGORY:
;
; CALLING SEQUENCE: met2ut,met,mission_name
;
; INPUT ARGUMENTS: Variable Name         Type                    Description
;
; OUTPUTS: Variable Name         Type                    Description
;   return       structure timestruct   stucture filled with universal time values
; INPUT KEYWORDS: these are included for backwards compatability
;   argos     argos P91 mission, this is the default
;   sts39     sts-39 UVLIM mission
;   STRV1D       CERTO-PLUS mission
;   f16          DMSP F16 mission (SSULI)
;   f17          DMSP F17 mission (SSULI)
;   f18          DMSP F17 mission (SSULI)
;       COSMIC   FORMOSAT3/COSMIC (TIP)
;       ANDERR   STS-116 ANDE Risk Reduction
;       ANDE     STS-127 ANDE
;       RAIDS    HTV-1 RAIDS to ISS
;
; OUTPUT KEYWORDS:
;
; COMMON BLOCKS:
;
; CALLED ROUTINES:
;   fjulian
;   jul2cal
;   mmddyyyyhhmm_d
;
; SIDE EFFECTS:
;   may define system variable !launchtimes by running launchtimesinit if not defined
;
; RESTRICTIONS: This IDL routine is to be distributed only by the Naval
;          Research Laboratory, Code 7607.
;
; NOTES:
;
; MODIFICATION HISTORY:
;     Version 0.9: Scott Alan Budzien, NRL Code 7607, 07/07/99  met2ut
;     Version 1.1: Clyde Fortna, SFA inc. @NRL Code 7607, 07/17/2002 split off launchtimes structure
;-
;***************************************************************************************
function met2ut,met,mission_name,argos=argos,sts39=sts39,STRV1D=strv1d,F16=f16,F17=f17,F18=f18,$
         cosmic=COSMIC,anderr=ANDERR,ande=ANDE,raids=RAIDS
;
dayname = ['Sun','Mon','Tue','Wed','Thu','Fri','Sat']
;
; Backward compatibility
IF (n_params() eq 1) THEN BEGIN ; set the mission string from parameters
    if      keyword_set(sts39) then mission_name='STS-39' $
    else if keyword_set(strv1d) then mission_name='STRV-1D' $
    else if keyword_set(f16) then mission_name='F16' $
    else if keyword_set(cosmic) then mission_name='COSMIC' $
    else if keyword_set(f17) then mission_name='F17' $
    else if keyword_set(f18) then mission_name='F18' $
    else if keyword_set(anderr) then mission_name='ANDERR' $
    else if keyword_set(ande) then mission_name='ANDE' $
    else if keyword_set(raids) then mission_name='RAIDS' $
    else mission_name='ARGOS' ; default
ENDIF
DEFSYSV, '!launchtimes', EXISTS = exists
IF NOT exists THEN launchtimesinit ; initialize system value if not already there
s=(where(!launchtimes.mission_name eq mission_name))(0)
if s lt 0 then begin
    print, 'mission_name must be one of:'
    print, !launchtimes.mission_name
    STOP
    return,-1
endif
    met1 = met/(!launchtimes(s).ticksperday/86400d0) ; value of MET in seconds
    launch_sec = !launchtimes(s).gpssecofweek
    launch_week = !launchtimes(s).gpsweek
    launch_jd = fjulian(!launchtimes(s).year,!launchtimes(s).month, $
                        !launchtimes(s).day + (launch_sec mod 86400)/86400d0 )
;
sec = (launch_sec + met1) mod 604800d0
week = fix((launch_sec + met1)/604800d0 + launch_week)
jd = launch_jd + double(met1)/86400d0
caldat = jul2cal(long(jd),jd-long(jd))
fd = caldat.dy-fix(caldat.dy)
doy = mmddyyyyhhmm_d(caldat.mo,fix(caldat.dy),caldat.yr,0,0,0)
;
out = {timestruct,year:0,month:0B,day:0B,hour:0B,min:0B,sec:0.,$
       secondofday:0d0,secondofweek:0d0,dayofweek:0B,dayname:'',$
       dayofyear:0d0,julian:0d0,met:0L}
if n_elements(met) gt 1 then out = replicate(out,n_elements(met))
out.year = caldat.yr
out.month = caldat.mo
out.day = caldat.dy
out.hour = fd*24.
out.min = fix(fd*1440.) mod 60
out.secondofday = fd*86400d0
out.sec = out.secondofday mod 60.
out.dayofyear = doy
out.secondofweek = (((launch_sec + met1(*)) mod 604800d0) + 604800d0) mod 604800d0 ;avoid negative values
out.dayofweek = byte(out.secondofweek/86400d0)+1
out.dayname = dayname(out.dayofweek-1)
out.julian = jd(*)
out.met = met(*)
return,out
end

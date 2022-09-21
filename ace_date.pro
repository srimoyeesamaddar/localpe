


function ace_date,indate,$
from_yyyyddd=from_yyyyddd,from_ymd=from_ymd,from_tyr=from_tyr,from_snoedom=from_snoedom,from_aimdom=from_aimdom,from_scottdom=from_scottdom,from_today=from_today,$
to_yyyyddd=to_yyyyddd,to_ymd=to_ymd,to_tyr=to_tyr,to_snoedom=to_snoedom,to_aimdom=to_aimdom,to_scottdom=to_scottdom

; Note ths is a beta version and while there are no known problems it has not been heavily vetted!!!!!!!!
; 
; SMB, January 1, 2019 finally getting around to writing this
; should convert dates from any given format to any desired format

; note, the first day of any year or month is 1, not zero - but this is not true for "day of mission" or dom which all start at zero
; this approach allows the dom variables to work in arrays

day0_snoe =[1998,3,11]
day0_aim  =[2007,4,25]
day0_scott=[1966,3,11]

; first put it in a julian date (variable=infuldate) and then convert to whatever is requested

;;;;;;;;;; from_yyyyddd
if keyword_set(from_yyyyddd) then begin

year= long(indate)/1000l
doy = long(indate) - long(year)*1000l

doy_firstofmonth=long( julday(findgen(12)+1,fltarr(12)+1,fltarr(12)+year,fltarr(12),fltarr(12),fltarr(12))-julday(1,1,year,0,0,0) ) + 1L

diff=doy-doy_firstofmonth
ind=where(diff ge 0,n_ind)
if n_ind eq 0 then stop
month=ind[n_ind-1] + 1
dayofmonth=diff[ind[n_ind-1]] + 1

injuldate = julday(month,dayofmonth,year,0.,0.,0.)

endif ; from_yyyyddd

;;;;;;;;;; from_ymd
if keyword_set(from_ymd) then begin

ymd=indate
injuldate = julday(ymd[1],ymd[2],ymd[0],0.,0.,0.)

endif ; from_ymd

;;;;;;;;;; from_tyr
if keyword_set(from_tyr) then begin

year= double(long(indate))
if (long(indate) mod 4) eq 0 then yearlen=366. else yearlen=365.
fracyear=(indate-year)
injuldate=julday(1,1,year,0.,0.,0.)+fracyear*yearlen

endif ; from_tyr

;;;;;;;;;; from_snoedom
if keyword_set(from_snoedom) then begin

findyearsmonths,indate,day0_snoe,jday
injuldate=jday

endif ; from_snoedom

;;;;;;;;;; from_aimdom
if keyword_set(from_aimdom) then begin

findyearsmonths,indate,day0_aim,jday
injuldate=jday

endif ; from_aimdom

;;;;;;;;;; from_scottdom
if keyword_set(from_scottdom) then begin

findyearsmonths,indate,day0_scott,jday
injuldate=jday

endif ; from_scottdom


;;;;;;;;; from_today
if keyword_set(from_today) then begin

injuldate=systime(/julian)

endif ; from_today

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Now set up the output - outdate


;;;;;;;;;;;to_yyyyddd
if keyword_set(to_yyyyddd) then begin

caldat,injuldate,month,dayofmonth,year
doy=long(injuldate - julday(1,1,year,0.,0.,0.)) + 1L
outdate=long(year)*1000L + doy

endif ; to_yyyyddd

;;;;;;;;;;;to_ymd
if keyword_set(to_ymd) then begin

caldat,injuldate,month,dayofmonth,year
ymd=[year,month,dayofmonth]
outdate=ymd

endif ; to_ymd

;;;;;;;;;;;to_tyr
if keyword_set(to_tyr) then begin

caldat,injuldate,month,dayofmonth,year
if (year mod 4) eq 0 then yearlen=366. else yearlen=365.
fracyear=(injuldate - julday(1,1,year,0.,0.,0.))/yearlen
outdate=double(year)+fracyear
endif ; to_tyr

;;;;;;;;;;;to_snoedom
if keyword_set(to_snoedom) then begin

ndays = injuldate - julday(day0_snoe[1],day0_snoe[2],day0_snoe[0],0.,0.,0.)
outdate=long(ndays)

endif ; to_snoedom

;;;;;;;;;;;to_aimdom
if keyword_set(to_aimdom) then begin

ndays = injuldate - julday(day0_aim[1],day0_aim[2],day0_aim[0],0.,0.,0.)
outdate=long(ndays)

endif ; to_aimdom

;;;;;;;;;;;to_scottdom
if keyword_set(to_scottdom) then begin

ndays = injuldate - julday(day0_scott[1],day0_scott[2],day0_scott[0],0.,0.,0.)
outdate=long(ndays)

endif ; to_scottdom


return,outdate
end


;;;;;;;;;;Required Procedures.....


pro findyearsmonths,ndays,ymd,jday

; find the julian day of the day that is ndays from ymd

firstdom=julday(ymd[1],ymd[2],ymd[0],0.,0.,0.)
; get approx years
n_years=fix(ndays)/365 + 2

;get year
firstdayyear=julday(fltarr(n_years)+1,fltarr(n_years)+1,findgen(n_years)+1+ymd[2],0.,0.,0.)
diff=ndays - (firstdayyear-firstdom)
ind=where(diff ge 0,n_ind)
year=ind[0]+1+ymd[2]

;get doy, month, day
if n_ind ge 1 then begin
  doy=diff[ind[0]]
  doy2md,year,doy,md
  month=md[0]
  dayofmonth=md[1]
endif else begin
  if long(year) mod 4 eq 0 then dt=double(ndays)/366. else dt=double(ndays)/365.
  doy=firstdom+dt
  doy2md,year,doy,md
  month=md[0]
  dayofmonth=md[1]
endelse

jday=julday(month,dayofmonth,year,0.,0.,0.)-1

return
end

pro doy2md,year,doy,md

doy_firstofmonth=long( julday(findgen(12)+1,fltarr(12)+1,fltarr(12)+year,fltarr(12),fltarr(12),fltarr(12))-julday(1,1,year,0,0,0) ) + 1L

diff=doy-doy_firstofmonth
ind=where(diff ge 0,n_ind)
if n_ind eq 0 then stop
month=ind[n_ind-1] + 1
dayofmonth=diff[ind[n_ind-1]] + 1

md=[month,dayofmonth]

return
end





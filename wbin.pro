pro wbin,ener,del,lo,hi,w


enerlo=ener-del/2.0
enerhi=ener+del/2.0


ilo=max(where(lo ge enerlo,n_ilo))>0
ihi=min(where(hi le enerhi,n_ihi))<(n_elements(ener)-1)
if ihi eq -1 then ihi=n_elements(ener)-1

fraclo=(enerhi(ilo)-lo)/del(ilo)
frachi=(hi-enerlo(ihi))/del(ihi)

w=ener-ener

if (lo ne hi) then begin
    w(ilo)=fraclo
    w(ihi)=frachi
    if ihi-ilo ge 2 then w(ilo+1:ihi-1)=1.0
    w=w<1.0

endif else $
  w(ihi)=1. ;else begin

;if (lo eq hi) then w(ihi)=1. else begin
;     print,'Bin width should be zero-Something is wrong!! Check your code!!'
;     stop  
;endelse

  


return
end



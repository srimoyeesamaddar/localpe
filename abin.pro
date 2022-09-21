pro abin,arr,val,loc,frac

; this program finds the first location in arr that is greater than val
; it returns that location and the weight for that bin relative to the next bin lower

test=arr gt val
dum=min(reverse(test),iloc)
loc=n_elements(test)-iloc ; first array element with value larger than val

ind=where(arr gt val)
testloc=ind[0]
;if testloc ne loc then print,'You are a bonehead!',loc,testloc
loc=testloc
if dum ne 1 then begin
frac=1.0-abs(arr(loc)-val)/(arr(loc)-arr(loc-1))
endif else begin
loc=0;-1
dm1=0.0
dp1=0.0
endelse

return
end

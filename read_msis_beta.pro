;pro read_msis_beta

msis_filename='C:\Users\Srimoyee\Desktop\nrl_files\msis1.97\msis_idl_interface.so'  ; sharable idl executable file


iyd=1999201
sec=0.
alt=[55.,60.,65.,70.,75.,80.,85.,90.,95.,100.]
glat=0.
glong=0.
stl=0
f107a=100
f107=100
ap=5
mass=1
d=fltarr(9)
t=fltarr(2)

;fortran subroutine
;gtd8(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

for i=0, n_elements(alt)-1 do begin
   spawn,[msis_filename],result,/noshell  ;msis_gtd8_interface is the executable in fortran

;   S = CALL_EXTERNAL(msis_filename, $
;       'msis_gtd8_interface', iyd,sec,alt,$
;        glat,glong,stl,f107a,f107,ap,mass,d,t) ;msis_gtd8_interface is the fortran wrapper
  
print,d(1)

  
endfor
  
  


end

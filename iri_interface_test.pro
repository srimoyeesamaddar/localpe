;Using the new IRI model to calculate the electron temperature and electron density values
;___________________________________________________________________________________________________________________________________________________
;Making some inputs for testing

z0=60.
zz=findgen(161)*2.+z0

alati=0.
along=0.

iyyyy=fix(2016)
mmdd=fix(0723)


hour=12.    ;Universa time in hours 0-24 hrs
utsec=3600.*hour  ;Universal time in seconds
;localhour=(((utsec/240. + lon)/15. + 24. ) mod 24.)
localhour=(utsec/3600. + along/15.)  mod 24.  ;My version NEEDS TO BE CHECKED!!
;dhour=hour+25.                             ;LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
dhour=localhour

f107a=70.
f107=70.
;___________________________________________________________________________________________________________________________________________________


;#1 make a file with the inputs to the IRI model

openw, lun_iriinp,'/home/srimoyee/Desktop/nrl_files/sav_files/iri_inputs.dat',/get_lun
printf,lun_iriinp,'alati', 'along', 'iyyyy','mmdd','iut','dhour',$
                    'f107a','f107','height'
printf,lun_iriinp,n_elements(zz)
printf,lun_iriinp,alati
printf,lun_iriinp,along
printf,lun_iriinp,iyyyy
printf,lun_iriinp,mmdd
printf,lun_iriinp,dhour
printf,lun_iriinp,f107a
printf,lun_iriinp,f107
printf,lun_iriinp,zz

close, lun_iriinp
;__________________________________________________________________________________________

;#2 make the iri executable file

spawn,$
  'f77 -o iri_idl_interface.exe iri_idl_interface.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for'
;__________________________________________________________________________________________

;#3 call the fortran wrapper using spawn command

iri_filename=['/home/srimoyee/Desktop/iri/iri_idl_interface.exe']  ; fortran executable file/wrapper. 
spawn,iri_filename,result,/noshell

;__________________________________________________________________________________________

;#4 open the output file generated by fortran IRI to get the varibles

iri_outdata=dblarr(2,n_elements(zz))  ;2 columns: eden and etemp

iri_outfile='/home/srimoyee/Desktop/nrl_files/sav_files/iri_outputs.dat'
openr, lun_iri, iri_outfile, /get_lun
skip_lun, lun_iri, 1, /lines ; 1 header line
readf, lun_iri , iri_outdata
free_lun, lun_iri

;__________________________________________________________________________________________

;#5 save the file data into the local variables and sav file

floc='/home/srimoyee/Desktop/nrl_files/sav_files/'
;floc='/home/srimoyee/Desktop/'

iri_eden=reform(iri_outdata(0,*))
iri_etemp=reform(iri_outdata(1,*))

save, zz,iri_eden,iri_etemp, filename=floc+'iri_edentemp_test.sav'

end
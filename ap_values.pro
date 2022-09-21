



;idate =2016205  ;m5 flare
idate=2017249   ;x9 flare


str_idate=strtrim(string(fix(idate/1000)),1)
ap_file= "/home/srimoyee/Desktop/ap_files/KP_AP/"+str_idate;+".txt"
openr, lun_file, ap_file, /get_lun

;ip_mnth=07 ;m5
;ip_day=23

ip_mnth=09 ;x9
ip_day=06


; hour= 5  ; m5 flare UT time
 hour =12  ; x9 flare UT time
ip_hr= hour
;_________________________________________________________________________________________________________________

ap_UT= fltarr(9)

if (ip_hr eq 24) then ip_hr=0

format= '$(I2,I2,I2,I4,I2,I2,I2,I2,I2,I2,I2,I2,I2,I3,I3,I3,I3,I3,I3,I3,I3,I3,I3)';,F3.1,I1,I3,F5.1, I1
;a1=fix(0)
;a2=fix(0)
;a3=fix(0)
;a4=fix(0)
;a5=fix(0)
;a6=fix(0)
;a7=fix(0)
;a8=fix(0)
;a9=fix(0)
;a10=fix(0)
;a11=fix(0)
;a12=fix(0)
;a13=fix(0)
;a14=fix(0)
;a15=fix(0)
;a16=fix(0)
;a17=fix(0)
;a18=fix(0)
;a19=fix(0)
;a20=fix(0)
;a21=fix(0)
;a22=fix(0)
;a23=0.0
;a24=fix(0)
;a25=fix(0)
;a26=0.0
;a27=fix(0)
while ~eof(lun_file) do begin
  
  readf, lun_file, format, $
         a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23;,a24,a25,a26,a27,a28
        
          if ((a2 eq ip_mnth) and (a3 eq ip_day))then begin
            ap_UT[0] = a15
            ap_UT[1] = a16
            ap_UT[2] = a17
            ap_UT[3] = a18
            ap_UT[4] = a19
            ap_UT[5] = a20
            ap_UT[6] = a21
            ap_UT[7] = a22
            ap_UT[8]=  a23
            break
            
          endif
        
endwhile

op_ap=ap_UT[fix(ip_hr/3)]

close,lun_file

;1- 2     I2    YEAR
;3- 4     I2    MONTH
;5- 6     I2    DAY
;
;7-10     I4    BARTELS SOLAR ROTATION NUMBER--a sequence of 27-day intervals
;               counted continuously from February 8, 1832.
;11-12     I2    NUMBER OF DAY within the Bartels 27-day cycle.
;
;13-14     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 0000 - 0300 UT.
;15-16     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 0300 - 0600 UT.
;17-18     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 0600 - 0900 UT.
;19-20     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 0900 - 1200 UT.
;21-22     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 1200 - 1500 UT.
;23-24     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 1500 - 1800 UT.
;25-26     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 1800 - 2100 UT.
;27-28     I2    Kp or PLANETARY 3-HOUR RANGE INDEX for 2100 - 2400 UT.
;29-31     I3    SUM of the eight Kp indices for the day expressed to the near-
;                est third of a unit.
;
;32-34     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 0000 - 0300 UT.
;35-37     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 0300 - 0600 UT.
;38-40     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 0600 - 0900 UT.
;41-43     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 0900 - 1200 UT.
;44-46     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 1200 - 1500 UT.
;47-49     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 1500 - 1800 UT.
;50-52     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 1800 - 2100 UT.
;53-55     I3    ap or PLANETARY EQUIVALENT AMPLITUDE for 2100 - 2400 UT.
;56-58     I3    Ap or PLANETARY EQUIVALENT DAILY AMPLITUDE--the arithmetic mean
;                of the day's eight ap values.
;
;59-61     F3.1  Cp or PLANETARY DAILY CHARACTER FIGURE--a qualitative estimate
;                of overall level of magnetic activity for the day determined
;                from the sum of the eight ap amplitudes.  Cp ranges, in steps
;                of one-tenth, from 0 (quiet) to 2.5 (highly disturbed).
;
;62-62     I1    C9--a conversion of the 0-to-2.5 range of the Cp index to one
;                digit between 0 and 9.
;
;63-65     I3    INTERNATIONAL SUNSPOT NUMBER.  Records contain the Zurich num-
;                ber through December 31, 1980, and the International Brus-
;                sels number thereafter.
;
;66-70     F5.1  OTTAWA 10.7-CM SOLAR RADIO FLUX ADJUSTED TO 1 AU--measured at
;                1700 UT daily and expressed in units of 10 to the -22 Watts/
;                meter sq/hertz.  Observations began on February 14, 1947.
;                From that date through December 31, 1973, the fluxes given
;                here don't reflect the revisions Ottawa made in 1966. NOTE:
;                If a solar radio burst is in progress during the observation
;                the pre-noon or afternoon value is used (as indicated by a
;                flux qualifier value of 1 in column 71.
;
;71-71     I1    FLUX QUALIFIER.  "0" indicates flux required no adjustment;
;                "1" indicates flux required adjustment for burst in progress
;                at time of measurement; "2" indicates a flux approximated by
;                either interpolation or extrapolation; and "3" indicates no
;                observation.





end
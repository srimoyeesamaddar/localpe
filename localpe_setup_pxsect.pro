pro localpe_setup_pxsect

@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
; see this file for definitions
 
  nst=10 ; NST     number of states produced by photoionization/dissociation, According to Conway 1988
;  nst=3
;___________________________________________________________________________________________________________________________________________________________________  
;  Input cross-section files
; 
; probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states_SQ05.sav'
;;  file has : O_prob_state,O2_prob_state,N2_prob_state
; ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec_SQ05.sav'
;;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2 

;____________________________________________________________________________________________

probstate_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_prob_states.sav';THIS ONE WAS USED 17JUL2022
;  file has : O_prob_state,O2_prob_state,N2_prob_state
ionabs_crss='/home/srimoyee/Desktop/nrl_files/sav_files/new_crosssec.sav'
;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2

;probstate_crss='/home/srimoyee/Desktop/new_prob_states.sav'
;;  file has : O_prob_state,O2_prob_state,N2_prob_state
;ionabs_crss='/home/srimoyee/Desktop/new_crosssec.sav'
;;  file has: spec_wave,sigi_o,sigi_o2,sigi_n2,sigab_o,sigab_o2,sigab_n2

;____________________________________________________________________________________________

restore, probstate_crss
restore, ionabs_crss
;o_prob_state[1,0:5]=0.39314517
;o_prob_state[2,0:5]=0.38104841 
;o_prob_state[3,0:5]=0.22580646


;o_prob_state[1,0:5]=0.0
;o_prob_state[2,0:5]=1.0 
;
;n2_prob_state[1,0:5]=0.039999999 
;n2_prob_state[2,0:5]=0.95999998
;
;
;
;sigi_o[0:5]= 2.2816000e-21
;sigab_o[0:5]= 2.3000001e-21 
;
;sigi_o2[0:5]= 4.4999999e-21 
;sigab_o2[0:5]= 4.4999999e-21
;
;
;sigi_n2[0:5]= 2.5000001e-21
;sigab_n2[0:5]= 2.5000001e-21 

lmax= n_elements(spec_wave)

;____________________________________________________________________________________
;cal=[2.48E-05,1.47E-04,4.57E-04,1.02E-03,1.90E-03,4.00E-03,8.55E-03,$
;  1.55E-02,3.12E-02,6.32E-02,0.14,0.3]*1e-18
;
;sigi_n2[0:11]=cal
;sigab_n2[0:11]=cal
;____________________________________________________________________________________

;___________________________________________________________________________________________________________________________________________________________________

;;     TPOT    ionization potentials for each species, state    ; eV 

;    tpot=[[13.61, 16.93, 18.63, 28.50, 40.00,  0.00],$    ; Dr Bailey's compilation 
;     [12.07, 16.10, 18.20, 20.00,  0.00,  0.00],$
;     [15.60, 16.70, 18.80, 30.00, 34.80, 25.00]]
;


      tpot=[[13.6, 16.9, 18.6, 28.50, 40.00, 531.70,0.0,0.0,0.0,0.0],$  ;adding extra potentials for N2 and O2 for extra states
        [12.1, 16.10, 18.20, 20.3,23.2,27.2, 33.0,39.8,531.70,0.0],$      ;Conway 1988
        [15.60,16.70, 18.80, 25.3,29.0,33.40,36.80,37.8, 43.6,400.]] ;Srimoyee was here !!


;  
;  tpot=[[13.61, 16.93, 18.63],$    ; For NRL with two states O2+/N2+ and Frag
;         [12.07, 18.20, 0.0],$
;         [15.60,25.3, 0.0]]
    
    
;
;tpot=[[13.61, 16.93, 18.63],$    ; For SQ'05 and Fe line
;       [12.07, 20.3, 0.0],$
;       [15.60,25.3, 0.0]]  

;  ___________________________________________________________________________________________________________________________



di=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0],$              ;dissociation energy arraged in tpot fashion for use in photoionisation code
  [0.0, 0.0, 0.0, 5.115,5.115,5.115,5.115,5.115,0.0,0.0],$      ;tpot states are from Conway 1988
  [0.0,0.0, 0.0, 9.759,9.759,9.759,9.759,9.759, 0.0,0.0]]       ;keeping the dissociation energy for Auger State to be zero. Will modify it during photoionisation calculation


  ;  ___________________________________________________________________________________________________________________________


  ;  flag states as Ionisation, D.I or Auger
  ;  1-ionisation, 2- D.I, 3- Auger
  
    flag_states=[[1, 1, 1, 1, 1, 3,0.0,0.0,0.0,0.0],$              
                 [1, 1, 1, 2,2,2,2,2,3,0.0],$      
                 [1,1, 1, 2,2,2,2,2, 1,3]]     

  ;  ___________________________________________________________________________________________________________________________


tpot=transpose(tpot)
di=transpose(di)
flag_states=transpose(flag_states)

;  ___________________________________________________________________________________________________________________________
;From Conway 1988

auger_energy=[531.70,531.70,400.]
auger_wvln=12397./auger_energy; O, O2 and N2


;___________________________________________________________________________________________________________________________________________________________________

;     ..............................................................   
;;    C PROBO, PROBO2, PROBN2; branching ratio data arrays............Dr. Bailey's file
;    filename='read_ephotn2.sav'  ;from read_ephotn2.pro
;    restore,file=filename ; gives aa,bb,probn2,sigin2,sigan2
;    filename='read_ephotoO2.sav'  ;from read_ephotoO2.pro
;    restore,file=filename  ; gives aa,bb,probo2,sigio2,sigao2
;    filename='read_ephotoO.sav'  ;from read_ephotn2.pro
;    restore,file=filename ; gives aa,bb,probo,sigio,sigao
;.....................................................................
;.....................................................................
    
; SIGABS  photoabsorption cross sections,0- O, 1-O2, 2-N2; cm2
; SIGIONx  photoionization cross sections,0- O,1- O2,2- N2; cm2

sigionx=fltarr(nmaj,lmax)
sigabs=fltarr(nmaj,lmax)
     
sigabs(0,*)=sigab_o
sigabs(1,*)=sigab_o2
sigabs(2,*)=sigab_n2

sigionx(0,*)=sigi_o
sigionx(1,*)=sigi_o2
sigionx(2,*)=sigi_n2

     

;PROB    branching ratios for each state, species, and wavelength bin:

prob=fltarr(nst,nmaj,lmax) ;states species and energies, 0-O, 1-O2, 2- N2


;0-O    
n_cols=(size(O_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
    for k=1,n_cols-1 do $
        prob(k-1,0,l)=O_prob_state(k,l)
    endfor

;1-O2    
n_cols=(size(O2_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
      for k=1,n_cols-1 do $
        prob(k-1,1,l)=O2_prob_state(k,l)
    endfor


;2-N2    
n_cols=(size(N2_prob_state))[1]  ;number of states

for l=0,lmax-1 do begin
      for k=1,n_cols-1 do $
        prob(k-1,2,l)=N2_prob_state(k,l)
    endfor

;prob[*,0,12]=prob[*,0,11]
;prob[*,1,12]=prob[*,1,11]
;prob[*,2,12]=prob[*,2,11]

first_pxsect=1

return
END





;pro localpe_setup_pxsect
;
;  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
;  ; see this file for definitions
;
;  ;  auger_energy=[500.,500.,360.]
;  auger_wvln=[23.,23.,30.]
;  auger_energy=12397./auger_wvln
;  nst=6  ; NST     number of states produced by photoionization/dissociation
;  ;     ..............................................................
;  ;;    C PROBO, PROBO2, PROBN2; branching ratio data arrays
;  ;    filename='read_ephotn2.sav'  ;from read_ephotn2.pro
;  ;    restore,file=filename ; gives aa,bb,probn2,sigin2,sigan2
;  ;    filename='read_ephotoO2.sav'  ;from read_ephotoO2.pro
;  ;    restore,file=filename  ; gives aa,bb,probo2,sigio2,sigao2
;  ;    filename='read_ephotoO.sav'  ;from read_ephotn2.pro
;  ;    restore,file=filename ; gives aa,bb,probo,sigio,sigao
;  ;;;    .............................................................
;  ;
;  ;    Srimoyee.......Replacing the above section with my read code
;  ;    Column_0: Low Wavelength Bin
;  ;    Column_1: High Wavelength Bin(A)]
;  ;    Column_2: X
;  ;    Column_3: A
;  ;    Column_4: B
;  ;    Column_5: C
;  ;    Column_6: F
;  ;    Column_7: Diss
;  ;    Column_8: TotIon
;  ;    Column_9: TotAbs
;
;  format="$((1x,f7.2,3x,f7.2,6(3x,f4.2),1x,f9.5,1x,f9.5))"
;  n_rows= lmax
;  n_columns=10
;
;  data_n2=  fltarr(n_columns,n_rows)
;  data_o2=fltarr(n_columns,n_rows)
;  data_o=fltarr(n_columns,n_rows)
;
;  file_n2='C:\Users\Srimoyee\Desktop\bailey_bins\newbins_n2.dat'
;  openr, lun_n2, file_n2, /get_lun
;  skip_lun, lun_n2, 4, /lines ;4 comment lines
;  readf, lun_n2 , data_n2,format=format
;  free_lun, lun_n2
;
;  file_o2='C:\Users\Srimoyee\Desktop\bailey_bins\newbins_o2.dat'
;  openr, lun_o2, file_o2, /get_lun
;  skip_lun, lun_o2, 4, /lines ;4 comment lines
;  readf, lun_o2 , data_o2,format=format
;  free_lun, lun_o2
;
;  file_o='C:\Users\Srimoyee\Desktop\bailey_bins\newbins_o.dat'
;  openr, lun_o, file_o, /get_lun
;  skip_lun, lun_o, 4, /lines ;4 comment lines
;  readf, lun_o , data_o,format=format
;  free_lun, lun_o
;
;  sigao=data_o[9,*]
;  sigao2=data_o2[9,*]
;  sigan2=data_n2[9,*]
;  sigio=data_o[8,*]
;  sigio2=data_o2[8,*]
;  sigin2=data_n2[8,*]
;
;  probn2=dblarr(nst,lmax)
;  probo2=dblarr(nst,lmax)
;  probo=dblarr(nst,lmax)
;
;  probn2[0,*]=data_n2[2,*]
;  probn2[1,*]=data_n2[3,*]
;  probn2[2,*]=data_n2[4,*]
;  probn2[3,*]=data_n2[5,*]
;  probn2[4,*]=data_n2[6,*]
;  probn2[5,*]=data_n2[7,*]
;
;  probo2[0,*]=data_o2[2,*]
;  probo2[1,*]=data_o2[3,*]
;  probo2[2,*]=data_o2[4,*]
;  probo2[3,*]=data_o2[5,*]
;  probo2[4,*]=data_o2[6,*]
;  probo2[5,*]=data_o2[7,*]
;
;  probo[0,*]=data_o[2,*]
;  probo[1,*]=data_o[3,*]
;  probo[2,*]=data_o[4,*]
;  probo[3,*]=data_o[5,*]
;  probo[4,*]=data_o[6,*]
;  probo[5,*]=data_o[7,*]
;
;  ;    ................................................................
;
;
;
;
;  ;aa is the lower wavelength bins
;  ;bb is the upper wavelength bins
;  wave_n2=(data_n2[0,*]+data_n2[1,*])/2.
;  wave_o2=(data_o2[0,*]+data_o2[1,*])/2.
;  wave_o=(data_o[0,*]+data_o[1,*])/2.
;
;  ;C SIGABS  photoabsorption cross sections, O, O2, N2; cm2
;  ;C SIGIONx  photoionization cross sections, O, O2, N2; cm2
;  sigabs=dblarr(nmaj,lmax) ;0 for O, 1 for O2, 2 for N2, 123 for wavelengths
;  sigionx=dblarr(nmaj,lmax) ;making this x so that it is not the same as the function sigion
;;      for l=0,lmax-1 do begin
;      sigabs(0,*)=sigao(*)*1e-18
;      sigabs(1,*)=sigao2(*)*1e-18
;      sigabs(2,*)=sigan2(*)*1e-18
;     sigionx(0,*)=sigio(*)*1e-18
;      sigionx(1,*)=sigio2(*)*1e-18
;      sigionx(2,*)=sigin2(*)*1e-18
;  ;    endfor
;
; 
;  ;TPOT    ionization potentials for each species, state; eV
;  tpot=[[13.61, 16.93, 18.63, 28.50, 40.00,  0.00],$
;    [12.07, 16.10, 18.20, 20.00,  0.00,  0.00],$
;    [15.60, 16.70, 18.80, 30.00, 34.80, 25.00]]
;  tpot=transpose(tpot)
;
;  ;C PROB    branching ratios for each state, species, and wavelength bin:
;  ;C         O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
;  ;C         O2+ states: X, a+A, b, dissoc.
;  ;C         N2+ states: X, A, B, C, F, dissoc.
;
;  prob=dblarr(nst,nmaj,lmax) ;states species and energies
;  l=0
;  for l=0,lmax-1 do begin
;    for k=0,nst-1 do begin
;      prob(k,0,l)=probo(k,l)
;      prob(k,1,l)=probo2(k,l)
;      prob(k,2,l)=probn2(k,l)
;    endfor
;  endfor
;
;  first_pxsect=1
;
;  return
;END
;


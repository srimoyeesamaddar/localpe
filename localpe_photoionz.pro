pro localpe_photoionz

  ; this procedure calculates the photoionization rates and fills the primary photoelectron
  ; spectrum array

  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
  ; see this file for definitions


  ;__________________________________________________________________________________________________________________________

  ; photoionization
  tau=fltarr(jmax,lmax)
  tau_wv=fltarr(nmaj,jmax,lmax)

  flux=fltarr(jmax,lmax)
  photoi=fltarr(nmaj,jmax)
  photoki=fltarr(nmaj,nst,jmax,lmax)

  photoi_wv=fltarr(nmaj,jmax,lmax)

  primary=dblarr(jmax,nbins) ; primary photoelectron spectrum
  ;__________________________________________________________________________________________________________________________


  ; start with photoionization and creating the primary electron spectrum at each altitude
  ;i-altitude, j-species, k-states, l- wavelength

  for i=jmax-1,0,-1 do begin  ;altitude


    tau(i,*)=zcol(0,i)*sigabs(0,*) + zcol(1,i)*sigabs(1,*) + zcol(2,i)*sigabs(2,*); optical depth
    flux(i,*)=ssflux*exp(-1.*tau(i,*))  ; solar flux at each altitude

    tau_wv(0,i,*)=zcol(0,i)*sigabs(0,*)
    tau_wv(1,i,*)=zcol(1,i)*sigabs(1,*)
    tau_wv(2,i,*)= zcol(2,i)*sigabs(2,*)


    for j=0,nmaj-1 do begin     ;species



      for k=0,nst-1 do begin    ;states


        for l=0,lmax-1 do begin  ;wavelength



          ephoton1=12397./wv1(l)
          ephoton2=12397./wv2(l)
          ;                 ephoton=12397./(0.5*(wv1(l)+wv2(l)))

          ep=tpot(j,k)
          if ((wv2[l] le auger_wvln[j])) then ep=auger_energy[j];+tpot(j,k)  ; auger ionization is when photon can break k-shell electron

          eelec1=ephoton1-ep
          eelec2=ephoton2-ep

          if (eelec1 gt 0.)  and (eelec2 gt 0.)then begin
            ;                    print, j,k,l,eelec1, eelec2

            kion=flux(i,l)*zmaj(j,i)*sigionx(j,l)*prob(k,j,l)      ; # ionizations of species j, alt i, state k, by wvln l
            photoki(j,k,i,l)= kion                         ; photoionization by species, state, altitude
            photoi(j,i)=photoi(j,i)+kion

            photoi_wv(j,i,l)=photoi_wv(j,i,l)+kion

            ;__________________________________________________________________________________________________________________________

            ;                       Primary electron energy calculation

            wbin,ener,del,eelec2,eelec1,w
            if max(w) eq 0. then w(0)=1.0

            ;                        print, ener[where(w ne 0)]
            ;                       Normal primary electron energy
            if (eelec1 eq eelec2) then $
              primary(i,0:nbins-1)=(primary(i,0:nbins-1)) + kion*w(0:nbins-1) $
            else  $
              primary(i,0:nbins-1)=(primary(i,0:nbins-1)) + kion*(w(0:nbins-1))*(del(0:nbins-1))/(eelec1-eelec2) ; fill the primary spectrum


            ;stop

                                 if (eelec1-eelec2) lt 0. then print,'Youve got a problem',i,j,k,l   ; check for any problems
            if kion lt 0. then print,'kion is negative',i,j,k,l
            if (min(w) lt 0.) then stop ;print,'your program sux',eelec1,eelec2

            ;                     Auger ionisation check

            if (wv2[l] le auger_wvln[j]) and (k eq 0)  then begin ; if auger ionization, a second electron is released

              ;              
;                                            dum=min(abs(ener-(auger_energy[j]-tpot(j,k))),ebin) ;Auger electron energy
              dum=min(abs(ener-(auger_energy[j])),ebin)
              ;                             print, 'Auger',ener[ebin]

              primary(i,ebin)=primary(i,ebin) + flux(i,l)*zmaj(j,i)*sigionx(j,l)
            endif
            ;
          endif ; eelec1 and eelec2 > 0.
          ;                   print,'________________________________________________________________________________________________________________-'

        endfor ; wavelength

      endfor ; state

      ;      photoi(j,i)=total(photoki(j,*,i)) ; photoionization by species and altitude

    endfor ; species


  endfor ; altitude

  ;__________________________________________________________________________________________________________________________

  ;primary=primary>1e-12
  if keyword_set(demo) then begin ; in demo mode, look only at 1 wavelength from solar flux, default bin 33, 300-310A
    primary=primary*0.
    primary[*,33]=100.
  endif





  return
end



;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;************************************************************************************************************************************************************
;**********************************************************TRYING MY NEW CODE HERE**************************************************************************
;

;
;pro localpe_photoionz
;
;  ; this procedure calculates the photoionization rates and fills the primary photoelectron
;  ; spectrum array
;
;  @ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
;  ; see this file for definitions
;
;
;  ;__________________________________________________________________________________________________________________________
;
;  ; photoionization
;  tau=fltarr(jmax,lmax)
;  tau_wv=fltarr(nmaj,jmax,lmax)
;
;  flux=fltarr(jmax,lmax)
;  photoi=fltarr(nmaj,jmax)
;  photoki=fltarr(nmaj,nst,jmax,lmax)
;
;  photoi_wv=fltarr(jmax,nmaj,lmax)
;
;  primary=dblarr(jmax,nbins) ; primary photoelectron spectrum
;  ;__________________________________________________________________________________________________________________________
;
;
;
;  ; start with photoionization and creating the primary electron spectrum at each altitude
;  ;i-altitude, j-species, k-states, l- wavelength
;
;  for i=jmax-1,0,-1 do begin  ;altitude
;
;
;    tau(i,*)=zcol(0,i)*sigabs(0,*) + zcol(1,i)*sigabs(1,*) + zcol(2,i)*sigabs(2,*); optical depth
;    flux(i,*)=ssflux*exp(-1.*tau(i,*))  ; solar flux at each altitude
;
;    tau_wv(0,i,*)=zcol(0,i)*sigabs(0,*)
;    tau_wv(1,i,*)=zcol(1,i)*sigabs(1,*)
;    tau_wv(2,i,*)= zcol(2,i)*sigabs(2,*)
;
;
;    for j=0,nmaj-1 do begin     ;species
;
;
;
;      for k=0,nst-1 do begin    ;states
;
;
;        for l=0,lmax-1 do begin  ;wavelength
;
;
;          ephoton1=12397./wv1(l)
;          ephoton2=12397./wv2(l)
;
;          ep=tpot(j,k)
;          eelec1=ephoton1-ep-di[j,k]  ;subtracts the dissociation energyy from the D.I if any dissociated states exist; Auger states dissociation treated separately
;          eelec2=ephoton2-ep-di[j,k]  ;Auger dissociation energy is kept zero here; added later for Auger electron
;
;          if (eelec1 gt 0.)  and (eelec2 gt 0.)then begin
;
;
;
;            kion=flux(i,l)*zmaj(j,i)*sigionx(j,l)*prob(k,j,l)      ; # ionizations of species j, alt i, state k, by wvln l
;
;            photoki(j,k,i,l)= kion                         ; photoionization by species, state, altitude
;            photoi(j,i)=photoi(j,i)+kion
;            photoi_wv(i,j,l)=photoi_wv(i,j,l)+kion
;
;;           __________________________________________________________________________________________________________________________________
;
;;           Primary electron energy calculation
;
;            wbin,ener,del,eelec2,eelec1,w
;            if max(w) eq 0. then w(0)=1.0
;
;            if (eelec1 eq eelec2) then $
;              primary(i,0:nbins-1)=(primary(i,0:nbins-1)) + kion*w(0:nbins-1) $
;            else  $
;              primary(i,0:nbins-1)=(primary(i,0:nbins-1)) + kion*(w(0:nbins-1))*(del(0:nbins-1))/(eelec1-eelec2) ; fill the primary spectrum
;
;
;            ;stop
;
;            ;                     if (eelec1-eelec2) lt 0. then print,'Youve got a problem',i,j,k,l   ; check for any problems
;            if kion lt 0. then print,'kion is negative',i,j,k,l
;            if (min(w) lt 0.) then stop ;print,'your program sux',eelec1,eelec2
;
;
;;__________________________________________________________________________________________________________________________________
;
;;                        Auger States
;
;            if (wv1[l] le auger_wvln[j]) and (flag_states[j,k] eq 3)  then begin ; if auger ionization, a second electron is released: flag_states = 3 is Auger state defined in pxsect
;
;
;              e2ener=  tpot(j,k)-(tpot(j,0)+di[j,0])
;              dum=min(abs(ener-e2ener), ebin)
;              primary(i,ebin)=primary(i,ebin) + flux(i,l)*zmaj(j,i)*sigionx(j,l)*prob(k,j,l)
;
;            endif
;
;
;
;
;
;          endif ; eelec1 and eelec2 > 0.
;
;
;
;        endfor ; wavelength
;
;      endfor ; state
;
;
;    endfor ; species
;
;
;  endfor ; altitude
;
;
;
;  END


;
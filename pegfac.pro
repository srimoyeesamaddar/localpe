;cd,'/mnt/snoesci/snoe/other/photoelectrons/Programs_labels/photoelectrons'
;.run pegfac.pro
;Note: Changed sigion to sigionx

restore,'pevar.sav'  ; from pe.pro


;JMAX    number of altitude levels
;LMAX    number of wavelength intervals for solar flux
;NMAJ    number of major species
;NEI     number of states produced by electron impact
;NBINS   number of energetic electron energy bins
;NST     number of states produced by photoionization/dissociation
;TAU     optical depth, dimensionless
;SIGABS  photoabsorption cross sections, O, O2, N2; cm2
;ZCOL    slant column density for species O, O2, N2, altitude; cm-2
;WV1   wavelength array, upper bound; Angstroms
;WV2   wavelength array, lower bound; Angstroms
;ZMAJ    density array for species O, O2, N2, altitude; cm-3
; PROB    branching ratios for each state, species, and wavelength bin:
;         O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
;         O2+ states: X, a+A, b, dissoc.
;         N2+ states: X, A, B, C, F, dissoc.
;TPOT    ionization potentials for each species, state; eV


; photoionization
tau=fltarr(jmax,lmax)
zcol=fltarr(nmaj,jmax)
flux=fltarr(jmax,lmax)
photoi=fltarr(nmaj,jmax)
photoki=fltarr(nmaj,nei,jmax)

primary=fltarr(jmax,nbins) ; primary photoelectron spectrum

for i=jmax-1,0,-1 do begin

   for j=0,2 do zcol(j,i)=total(zmaj(j,i:*))*dz*1e5
   tau(i,*)=zcol(0,i)*sigabs(0,*) + zcol(1,i)*sigabs(1,*) + zcol(2,i)*sigabs(2,*)
   flux(i,*)=ssflux*exp(-1.*tau(i,*))  ; solar flux at each altitude

   for j=0,nmaj-1 do begin
      for k=0,nst-1 do begin
         for l=0,lmax-1 do begin

            ephoton1=12400./wv1(l)
            ephoton2=12400./wv2(l)
            eelec1=ephoton1-tpot(j,k)
            eelec2=ephoton2-tpot(j,k)

            if (eelec1 gt 0.) and (eelec2 gt 0.) then begin

            kion=flux(i,l)*zmaj(j,i)*sigionx(j,l)*prob(k,j,l)      ; how many ionizations of species j, alt alt i, state k, by wvln l

           ; then print,'kion is negative',i,j,k,l
            if (min(w) lt 0.) then print,'your program sux',eelec1,eelec2
            photoki(j,k,i)=photoki(j,k,i)+kion

            wbin,ener,del,eelec2,eelec1,w
            if max(w) eq 0. then w(0)=1.0
            primary(i,*)=primary(i,*) + kion*w/(eelec1-eelec2)
            if (eelec1-eelec2) lt 0. then print,'Youve got a problem',i,j,k,l
            if kion lt 0. then print,'kion is negative',i,j,k,l
            if (min(w) lt 0.) then stop ;print,'your program sux',eelec1,eelec2

         endif

         endfor

      endfor
            
      photoi(j,i)=total(photoki(j,*,i))

   endfor
endfor


;;;;
pespec=fltarr(jmax,nbins)
;iz=30
for iz=0,jmax-1 do begin
foundfirst=0
   for iener=nbins-1,0,-1 do begin

      if not( foundfirst eq 0 and primary[iz,iener] eq 0.0) then begin

         if foundfirst eq 0 then begin
            foundfirst=1

; now look at highest energy production where production comes only from
; primary

            prod=primary[iz,iener]
            loss=0.
            
            for imaj=0,nmaj-1 do loss=loss+total(zmaj[imaj,iz]*siga[imaj,*,iener])

            if loss gt 0. then pespec(iz,iener) = prod / loss

         endif else begin

            prod=primary[iz,iener]
            
            for imaj=0,nmaj-1 do begin
               for iupper=nbins-1,iener,-1 do begin
                  prod=prod+zmaj[imaj,iz]*pespec[iz,iupper]*sec[imaj,iener,iupper]
               endfor
            endfor

            loss=0.
            for imaj=0,nmaj-1 do loss=loss+total(zmaj(imaj,iz)*siga(imaj,*,iener))

            if loss gt 0. then pespec(iz,iener) = prod / loss

         endelse
      endif

   endfor                       ; loop over energy
   endfor;loop over altitude

sign2a=sigex(0,2,*)+sigex(1,2,*)+sigex(2,2,*)
pn2a150=fltarr(1)
for i=0,nbins-1 do begin
pn2a150=pn2a150+sign2a(0,0,i)*pespec(25,i)*zmaj(2,25)*del(i)
endfor
print,pn2a150
end



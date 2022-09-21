pro ace_etransport

@ace_common_blocks.prg ; common blocks used for all atmospheric chemistry and energetics (ace) software
; see this file for definitions
        
common impitck, alpha, beta, gamma, psi, $
     delz, del2, dela, delp,$
     delm, dels, den, fac

toaflux=fltarr(nbins)
fac=0.0    
avemu=0.5
dip=23.34
sindip=sin(dip/!radeg)
rmusin=1./(sindip*avemu)
phiinf=toaflux/avemu

; setup arrays
upflux=fltarr(jmax,nbins)
downflux=fltarr(jmax,nbins)
prod=fltarr(jmax)
produp=fltarr(jmax,nbins)
proddown=fltarr(jmax,nbins)
l2thermals=fltarr(jmax,nbins)
ddcasc=fltarr(nbins)
ddsec=fltarr(nbins)
; Calcualte delta z's in cm
zzcm=zz*1e5
;_________________________________
;delz=fltarr(jmax)
;delz[1:*]=zzcm[1:*]-zzcm
;delz[0]=zzcm[1]-zz[0]
;________________________________
delz=(zz(1:*)-zz(0:-2))*1.e5
delz=[delz[0],delz]
;________________________________
del2=fltarr(jmax)
dela=fltarr(jmax)
delp=fltarr(jmax)
delm=fltarr(jmax)
dels=fltarr(jmax)

for i=0,jmax-2 do begin
        del2(i) = delz(i)+delz(i+1)
        dela(i) = del2(i)/2.
        delp(i) = dela(i)*delz(i+1)
        delm(i) = dela(i)*delz(i)
        dels(i) = delz(i)*delz(i+1)
 endfor
 
del2(jmax-1) = del2(jmax-2)
delp(jmax-1) = delp(jmax-2)
delm(jmax-1) = delm(jmax-2)
dels(jmax-1) = dels(jmax-2)

ienertest=112
; energy loop
firstbin=nbins-1
indp=where(primary[jmax-1,*] gt 0.,n_indp)
if n_indp gt 0 then firstbin=indp[n_elements(indp)-1]

for iener=firstbin,0,-1 do begin

    t1=fltarr(jmax)
    t2=fltarr(jmax)
    secprodup=fltarr(jmax,nbins)
    cascadeup=fltarr(jmax,nbins)
    secproddown=fltarr(jmax,nbins)
    cascadedown=fltarr(jmax,nbins)

    ; make alt arrays of primary production, loss to thermal electrons, and collision terms
    for iz=0,jmax-1 do begin

           prod(iz)=(primary(iz,iener) * rmusin / del[iener])>1e-30  ; primary production (photoionization) in per eV units
           
           ee=8.618e-5*etemp[iz]
           dag=ener[iener]-ener[(iener-1)>0]
           d1=(ener[iener]-ee)>0.
           d2=(ener[iener]-0.53*ee)
           if d2 gt 0. then begin
                l2thermals[iz,iener]=(  eden[iz] * 3.37e-12 * ((d1/d2)^2.36)   /  (ener[iener]^0.94)  /  ((eden[iz]^0.03)>1.0)  )>0. ; loss to thermals
           endif else l2thermals[iz,iener]=0.0
       
           if finite(l2thermals[iz,iener]) eq 0b then print,'thermal problem: ',eden[iz],d1,d2,ener[iener],ee,dag

           for inmaj=0,nmaj-1 do begin
                t1(iz)=t1(iz) + zmaj[inmaj,iz]* sigs[inmaj,iener] * pe[inmaj,iener]
                t2(iz)=t2(iz) + zmaj[inmaj,iz]*(sigs[inmaj,iener] * pe[inmaj,iener]  + sigloss[inmaj,iener])
           endfor 
       
    endfor 

    if iener gt 0 then dag=ener[iener]-ener[iener-1] else dag=del[0]      
    l2thermals[*,iener]=l2thermals[*,iener] * rmusin / dag
    t1=t1*rmusin
    t2=t2*rmusin + l2thermals[*,iener] 

    ; we now calculate the terms of the ODE that are solved; see SMB dissertation Chapter 3
    ; including Equations 3-1 to 3-11

    ; we'll need the following derivatives and arrays:
    dt1dz=deriv(zzcm,t1)
    dt2dz=deriv(zzcm,t2)
    dproddz=deriv(zzcm,prod>1e-30)
    dproddndz=deriv(zzcm,proddown[*,iener]>1e-30)
    
    ;if iener eq ienertest then stop

    phiup=fltarr(jmax)
    phidown=fltarr(jmax)
    alpha=fltarr(jmax)
    beta=fltarr(jmax)
    gamma=fltarr(jmax)
    proddown=proddown>1e-30
    den=fltarr(jmax)

    ; fill the arrays for each term of the ODE
    alpha[1:jmax-2]=(-1.0/t1[1:jmax-2])*dt1dz[1:jmax-2]
    beta[1:jmax-2]=(t1[1:jmax-2]^2 - t2[1:jmax-2]^2 + t2[1:jmax-2]/t1[1:jmax-2]*dt1dz[1:jmax-2]-dt2dz[1:jmax-2])
    signbeta=beta/abs(beta)


    beta=signbeta*(abs(beta)>1e-20)

    beta[where(~finite(beta),/null)]=0.

    beta[0]=0
    beta[jmax-1]=0.
    gamma[1:jmax-2]=prod[1:jmax-2]/2.0 *(-t1[1:jmax-2]-t2[1:jmax-2]-alpha[1:jmax-2] - dproddz[1:jmax-2]/prod[1:jmax-1]) + $
              proddown[1:jmax-2,iener] * (-1.*alpha[1:jmax-2]-t2[1:jmax-2] - dproddndz[1:jmax-2]/proddown[1:jmax-2,iener]) - $
              produp[1:jmax-2,iener]*t1[1:jmax-2]
    psi=fltarr(jmax)+1.0
  
;if iener eq ienertest then begin ; debugging help...
;print,'ener,sigs,pe,pi,sigloss'
;print,ener[iener],ener[iener],sigs[2,iener],pe[2,iener],pi[2,iener],total(sigloss[*,iener],1),format='$(2f10.2,4e12.2)'
;print,'ener,zz,eden,etemp,zmaj1.2.3'
;for ia=70,80 do print,ener[iener],zz[ia],eden[ia],etemp[ia],zt[ia],zmaj[0,ia],zmaj[1,ia],zmaj[2,ia],format='$(2f10.2,6e12.2)'
;print,'ener,zz,alpha,beta,gamma,t1,t2,thermals'
;for ia=80,90 do print,ener[iener],zz[ia],alpha[ia],beta[ia],gamma[ia],t1[ia],t2[ia],l2thermals[ia,iener],format='$(2f10.2,6e12.2)'
;for ia=70,80 do print,ener[iener],zz[ia],proddown[ia,iener],beta[ia],gamma[ia],t1[ia],t2[ia],l2thermals[ia,iener],format='$(2f10.2,6e12.2)'
;stop
;endif

;stop

      ; solve the equation for phidown for this energy
      phidown[1] = gamma[1]/beta[1]
      den[0]=phidown[1]
      fluxj=phiinf[iener]
      impitck,jmax,fluxj
      phidown=den

      ; now apply the boundary condition at the bottom wich is the phiup=phidown
      ; this works since the lower altitude should be low enough that any flux is essentially
      ; zero
      phiup[0]=phidown[0]

      ; now that we have phidown we integrate up in altitude to get phiup
      r1=(t1*phidown + (prod + 2.0*produp[*,iener]) / 2.0)  / t2
      taue=(t2*delz)<60.
      expt2=exp(-1.*taue)
      for iz=1,jmax-1 do phiup[iz]=r1[iz] + (phiup[iz-1]-r1[iz])*expt2[iz] 
 
 ;if iener eq ienertest then begin ; debugging help...
 ;if min(den) lt 0. then begin
 ;print,'got a neg value at ener=',ener[iener],iener
 ;print,'zz. den, phidn, phiup',ener[iener]
;for ia=80,90 do print,zz[ia],den[ia],phidown[ia],phiup[ia],r1(ia),taue(ia),expt2(ia),format='$(f10.2,6e12.2)'
;stop
;endif

      ; we now have phiup and phidown for this energy bin!!!
      upflux[*,iener]=phiup * avemu
      downflux[*,iener]=phidown * avemu

      ; now that we've solved for the current energy bin, we calculate the impact of the electrons in that bin in all
      ; lower bins: degradation of the primary (cascade) and production of secondaries

       dcasc=fltarr(jmax)
       dsec=fltarr(jmax)
       for iz=0,jmax-1 do begin
            for ilower=0,iener-1 do begin
                icasc=iener-1-ilower
                for imaj=0,nmaj-1 do begin
                      dummycascup =       zmaj[imaj,iz]*phidown[iz]*siga[imaj,ilower,iener]*pi[imaj,iener] + $
                                                      zmaj[imaj,iz]*phiup[iz]     *siga[imaj,ilower,iener]*(1.0-pi[imaj,iener])
                      dummycascdown = zmaj[imaj,iz]*phiup[iz]    *siga[imaj,ilower,iener]*pi[imaj,iener]  + $
                                                      zmaj[imaj,iz]*phidown[iz]*siga[imaj,ilower,iener]*(1.0-pi[imaj,iener])
                      dcasc[iz]=dcasc[iz]+dummycascdown
                      cascadeup[iz,icasc]     = cascadeup[iz,icasc]      + dummycascdown
                      cascadedown[iz,icasc]= cascadedown[iz,icasc] + dummycascup
                      produp[iz,icasc]      = produp[iz,icasc]      + dummycascup      * rmusin
                      proddown[iz,icasc] = proddown[iz,icasc] + dummycascdown * rmusin    
                     ;if dummycascdown lt 0. then stop
                      if icasc eq ienertest and iz eq 85 then ddcasc[iener]=ddcasc[iener]+dummycascdown * rmusin  
               
                endfor  ; imaj
            endfor    ;  ilower  
             
            for imaj=0,nmaj-1 do begin
               for ilower=iener-1,0,-1 do begin
                   dummysec = zmaj[imaj,iz]*(phiup[iz]+phidown[iz])*sec[imaj,ilower,iener]*0.5
                   secprodup[iz,ilower]=secprodup[iz,ilower]+dummysec
                   dsec[iz]=dsec[iz]+dummysec               
                   secproddown[iz,ilower]=secprodup[iz,ilower]
                   produp[iz,ilower]      = produp[iz,ilower]      + dummysec * rmusin
                   proddown[iz,ilower] = proddown[iz,ilower] + dummysec * rmusin     
                   if ilower eq ienertest and iz eq 85 then ddsec[iener]=ddsec[iener]+dummysec * rmusin  
            endfor
           endfor  
               
           if iener ne 0 then begin
                  kk=iener-1
                  produp[iz,kk]      = produp[iz,kk]      + l2thermals[iz,iener] * phiup[iz]       * del[iener]/del[kk]
                  proddown[iz,kk] = proddown[iz,kk] + l2thermals[iz,iener] * phidown[iz] * del[iener]/del[kk]
           endif
 
       endfor    ; iz
;if iener eq ienertest+9 then begin ; debugging help...
 ;print,'ener,thermals,sec,casc
;for ia=15,105 do print,zz[ia],ener[iener],l2thermals[ia,iener],dsec[ia],dcasc[ia];,format='$(2f10.2,5e12.2)'
;print,'ener,dsec,dcasc'
;for i=0,nbins-1 do print,i,ener[i],ddsec[i],ddcasc[i]
;endif

 
endfor ; iener





;Trying to add the eiionz here..............................
Rs=695510.
R=149.6e6
constant=1.;!pi;2.*!pi*((Rs/R)^2.)
tflux=(upflux+downflux)*2.*constant
tflux[where(~finite(tflux),/null)]=0. 

; calculate electron impact ionization rates
;eiionz,eiionzk
eiionzk_transp=fltarr(nei,nmaj,jmax)
eiionz_transp=fltarr(nmaj,jmax) ; electron impact ionization  of O, O2, and N2
eiionz_transp1=fltarr(nei,nmaj,jmax,nbins) 


exct_transp=fltarr(nei,nmaj,jmax) ;electron impact excitation
exct_transp1=fltarr(nei,nmaj,jmax,nbins) ;electron impact excitation


;Making a copy of sigex to test predissociation of N2**
sigex_copy=sigex
eii=[where(ener le 10,/null)]
sigex_copy(5,2,eii)=0.
sigex_copy(5,2,*)=0.12*sigex_copy(5,2,*)
sigex_copy(6,2,*)=0.50*sigex_copy(6,2,*)
sigex_copy(9,2,*)=0.95*sigex_copy(9,2,*)
sigex_copy(10,2,*)=0.10*sigex_copy(10,2,*)
sigex_copy(11,2,*)=0.84*sigex_copy(11,2,*)
sigex_copy(16,2,*)=0.99*sigex_copy(16,2,*)


;Zero out predissociation of O2 below dissociation limit- AA'c 5.583 eV and SR 7.083 eV**
eii=[where(ener le 5.583,/null)]
sigex_copy(1,1,eii)=0.

eii=[where(ener le 7.083,/null)]
sigex_copy(2,1,eii)=0.


;;Making a copy of sigex to test predissociation of N2**-GLOW cross-sections
;sigex_copy=sigex
;sigex_copy(2,2,*)=0.50*sigex_copy(2,2,*)
;sigex_copy(5,2,*)=0.84*sigex_copy(5,2,*)
;
;;Zero out predissociation cross-section of O2 below dissociation limit- AA'c 5.583 eV **- GLOW cross-sections
;
;eii=[where(ener le 5.583,/null)]
;sigex_copy(2,1,eii)=0.



for i=0,jmax-1 do begin  ;altitude
    
    for j=0,nmaj-1 do begin  ;species
    
    for l=0, nbins-1 do begin
     

    ionz=0.
    
    for k=0,nei-1 do begin  ;states

      ; need to loop over ion states / ionizatoin processes to get total ionization rate

      eiionzk_transp[k,j,i]= eiionzk_transp[k,j,i]+(tflux[i,l]*del[l]*sigix[k,j,l]*zmaj[j,i])

      eiionz_transp1[k,j,i,l]=  eiionz_transp1[k,j,i,l]+(tflux[i,l]*del[l]*sigix[k,j,l]*zmaj[j,i])
     
      eiionz_transp[j,i]=eiionz_transp[j,i]+ ( tflux[i,l]*del[l]*sigix[k,j,l]*zmaj[j,i]  )
      
;      ionz=ionz+eiionzk_transp[k,j,i]

;     need to loop over ion states / ionizatoin processes to get total ionization rate

      exct_transp[k,j,i]= exct_transp[k,j,i]+(tflux[i,l]*del[l]*sigex_copy[k,j,l]*zmaj[j,i])
;      exct_transp[k,j,i]= exct_transp[k,j,i]+(tflux[i,l]*del[l]*sigex[k,j,l]*zmaj[j,i])
      
      exct_transp1[k,j,i,l]=  exct_transp1[k,j,i,l]+(tflux[i,l]*del[l]*sigex_copy[k,j,l]*zmaj[j,i])
;      exct_transp1[k,j,i,l]=  exct_transp1[k,j,i,l]+(tflux[i,l]*del[l]*sigex[k,j,l]*zmaj[j,i])
      
      
    endfor  ; k loop

;    eiionz_transp[j,i]=ionz

  endfor ; l loop

endfor ; j loop


endfor  ; i loop


;.....................................................................................

return
end

; Procedure impitck  solves a parabolic differential equation by implicit
; crank-nicholson method
; this is taken directly from GLOW (etrans.f) and translated to IDL
; see GLOW documentation for further details

Pro impitck,jmax,fluxj

k=fltarr(jmax)
l=fltarr(jmax)
a=fltarr(jmax)
b=fltarr(jmax)
c=fltarr(jmax)
d=fltarr(jmax)
      
common impitck, alpha, beta, gamma, psi, $
     delz, del2, dela, delp,$
     delm, dels, den, fac

      i1 = jmax - 2

        a[0:i1] = psi[0:i1] / delp[0:i1] + alpha[0:i1] / del2[0:i1]
        b[0:i1] = -2. * psi[0:i1] / dels[0:i1] + beta[0:i1]
        c[0:i1] = psi[0:i1] / delm[0:i1] - alpha[0:i1] / del2[0:i1]
        d[0:i1] = gamma[0:i1]
   
      k(1) = (d(1) - c(1)*den(0)) / b(1)
      l(1) = a(1) / b(1)
      
for i = 2, i1 do begin
        dem = b(i) - c(i) * l(i-1)
        k(i) = (d(i) - c(i)*k(i-1)) / dem
        l(i) = a(i) / dem
endfor

      den(i1) = (k(i1) - l(i1)*fluxj) / (1. + l(i1)*fac)
      den(jmax-1) = den(i1)
         
for kk = 0, jmax-4 do begin
        jk = i1 - kk
        den(jk) = k(jk) - l(jk) * den(jk + 1)
endfor
   
      return
      end



!!!
!!! Authors:  Douglas P. Drob and John Emmert, NRL 7632
!!!
!!! =================================================
!!! MSIS temperature vertical profile functions
!!! =================================================
    
module tfn  

    use constants,only:  nl
    
    type tnparm
        sequence
        real(8)              :: cf(0:nl)  ! Spline coefficients
        
        real(8)              :: tzetaF    ! Tn at zetaF
        
        real(8)              :: tzetaA    ! Tn at zetaA (reference altitude for O1, H1)
        real(8)              :: dlntdzA   ! log-temperature gradient at zetaA (km^-1)
        
        real(8)              :: lndtotF   ! ln total number density at zetaF (m^-3)

        real(8)              :: tex
        real(8)              :: tgb0
        real(8)              :: tb0
        real(8)              :: sigma
        real(8)              :: sigmasq
        real(8)              :: b       ! b = 1-tb0/tex
        
        real(8)              :: beta(0:nl)  ! 1st integration coefficients on k=5 splines 
        real(8)              :: gamma(0:nl)  ! 2nd integration coefficients on k=6 splines 
        real(8)              :: cVs         ! 1st integration constant (spline portion)
        real(8)              :: cVb         ! 1st integration constant (Bates portion)
        real(8)              :: cWs         ! 2nd integration constant (spline portion)
        real(8)              :: cWb         ! 2nd integration constant (Bates portion)
        
        real(8)              :: VzetaF      ! 1st indefinite integral at zetaF
        real(8)              :: VzetaA      ! 1st indefinite integral at zetaA
        real(8)              :: WzetaA      ! 2nd indefinite integral at zetaA
        real(8)              :: Vzeta0      ! 1st indefinite integral at zeta=0 (needed for pressure calculation)
        
    end type tnparm

    contains
    
! ==========================================================
! Compute the vertical temperature and species independent
! profile parameters
! ==========================================================

    subroutine tfnparm(gf,tpro)

        use constants
        use msis,only: smod,sfluxmod,geomag,utdep,TN,PR

        implicit none
        real(8),intent(in)       :: gf(0:maxnbf-1)    
        type(tnparm),intent(out) :: tpro 
    
        integer(4)               :: ix
        real(8)                  :: bc(3)
        
        real(8), external        :: dilog
                        

        ! Unconstrained spline coefficients
        do ix = 0,itb0-1
            tpro%cf(ix) = dot_product(TN%beta(0:mbf,ix),gf(0:mbf))
            if (smod(ix)) then
                tpro%cf(ix) = tpro%cf(ix) + sfluxmod(ix,gf,TN,1d0/TN%beta(0,ix))    !sfluxmod adds F10.7 modulation of tides
            endif
        enddo

        ! Bates temperature profile parameters
        tpro%tex = dot_product(TN%beta(0:mbf,itex),gf(0:mbf))
        !if (smod(itex)) tpro%tex = tpro%tex + sfluxmod(itex,gf,TN,1d0/TN%beta(0,itex))         
        tpro%tex = tpro%tex + sfluxmod(itex,gf,TN,1d0/TN%beta(0,itex))         
        tpro%tex = tpro%tex + geomag(TN%beta(cmag:cmag+nmag-1,itex),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
        tpro%tex = tpro%tex + utdep(TN%beta(cut:cut+nut-1,itex),gf(cut:cut+8))
        
        tpro%tgb0 = dot_product(TN%beta(0:mbf,itgb0),gf(0:mbf))
        if (smod(itgb0)) tpro%tgb0 = tpro%tgb0 + sfluxmod(itgb0,gf,TN,1d0/TN%beta(0,itgb0))         
        tpro%tgb0 = tpro%tgb0 + geomag(TN%beta(cmag:cmag+nmag-1,itgb0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
        
        tpro%tb0 = dot_product(TN%beta(0:mbf,itb0),gf(0:mbf))
        if (smod(itb0)) tpro%tb0 = tpro%tb0 + sfluxmod(itb0,gf,TN,1d0/TN%beta(0,itb0))         
        tpro%tb0 = tpro%tb0 + geomag(TN%beta(cmag:cmag+nmag-1,itb0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))

        tpro%sigma = tpro%tgb0/(tpro%tex-tpro%tb0)

        ! Spline coefficients constrained so that profile is C2 continuous
        bc(1) = 1.0d0/tpro%tb0
        bc(2) = -tpro%tgb0/(tpro%tb0*tpro%tb0)
        bc(3) = -bc(2)*(tpro%sigma + 2.0d0*tpro%tgb0/tpro%tb0)    
        tpro%cf(itb0:itex) = matmul(bc, c2tn)

        ! Reference temperature at zetaF (70 km)
        tpro%tzetaF = 1.0d0 / dot_product(tpro%cf(izFx:izFx+2),S4zetaF)

        ! Reference temperature and gradient at zetaA (85 km)
        tpro%tzetaA = 1.0d0 / dot_product(tpro%cf(izAx:izAx+2),S4zetaA)
        tpro%dlntdzA = -dot_product(tpro%cf(izAx:izAx+2),wghtAxdz) * tpro%tzetaA
        
        ! Calculate spline coefficients for first and second integrals 
        tpro%beta(0) = tpro%cf(0)*wbeta(0)
        do ix = 1, nl
            tpro%beta(ix) = tpro%beta(ix-1) + tpro%cf(ix)*wbeta(ix)
        enddo
        tpro%gamma(0) = tpro%beta(0)*wgamma(0)
        do ix = 1, nl
            tpro%gamma(ix) = tpro%gamma(ix-1) + tpro%beta(ix)*wgamma(ix)
        enddo
        
        ! Integration terms and constants
        tpro%b = 1 - tpro%tb0 / tpro%tex
        tpro%sigmasq = tpro%sigma * tpro%sigma
        tpro%cVS = -dot_product(tpro%beta(itb0-1:itb0+2),S5zetaB)
        tpro%cWS = -dot_product(tpro%gamma(itb0-2:itb0+2),S6zetaB)
        tpro%cVB = -log(1-tpro%b) / (tpro%sigma * tpro%tex)
        tpro%cWB = -dilog(tpro%b) / (tpro%sigmasq * tpro%tex)
        tpro%VzetaF = dot_product(tpro%beta(izfx-1:izfx+2),S5zetaF) + tpro%cVS
        tpro%VzetaA = dot_product(tpro%beta(izax-1:izax+2),S5zetaA) + tpro%cVS
        tpro%WzetaA = dot_product(tpro%gamma(izax-2:izax+2),S6zetaA) + tpro%cVS*(zetaA-zetaB) + tpro%cWS
        tpro%Vzeta0 = dot_product(tpro%beta(0:2),S5zeta0) + tpro%cVS

        ! Compute total number density at zetaF
        tpro%lndtotF = lnP0 - fh*(tpro%VzetaF - tpro%Vzeta0) - log(boltz*tpro%TzetaF)

        return
    
    end subroutine tfnparm
    
    
!!! =================================================
!!! MSIS temperature function
!!! =================================================

    real(8) function tfnx(z,iz,wght,tpro)

      use constants

      implicit none
      real(8),intent(in)      :: z
      integer,intent(in)      :: iz
      real(8),intent(in)      :: wght(-3:0)
      type(tnparm),intent(in) :: tpro
      integer                 :: i, j
  
      if (z .lt. zetaB) then 
         ! pure and matched B-splines
         i = max(iz-3,0)
         if (iz .lt. 3) then
             j = -iz
         else
             j = -3
         endif
         tfnx = 1.0d0 / dot_product(tpro%cf(i:iz),wght(j:0))
      else
         ! Thermosphere (pure Bates profile region)
         tfnx = tpro%tex - (tpro%tex - tpro%tb0)*exp(-tpro%sigma * (z - zetaB))
      endif

      return

    end function tfnx
    
end module tfn
    

!!!
!!! Authors:  Douglas P. Drob and John Emmert, NRL 7632
!!!
!!! =================================================
!!! MSIS density and composition vertical profiles
!!! =================================================

module dfn

    use constants, only : nl, nsplO1, nsplNO
    type dnparm
        sequence

        real(8)              :: lnPhiF      ! (Except O, H) Natural log of mixing ratio at zetaF (70 km), before chemical and dynamical corrections are applied (ln m^-3) (global term only)
        real(8)              :: lndref      ! Natural log of number density at reference height

        real(8)              :: zetaM       ! "Turbopause Height": Height of midpoint of effective mass transition (km)
        real(8)              :: HML         ! Scale height of lower portion of effective mass profile (km)
        real(8)              :: HMU         ! Scale height of upper portion of effective mass profile (km)

        real(8)              :: C           ! Chapman term coefficient
        real(8)              :: zetaC       ! Chapman term reference height (km)
        real(8)              :: HC          ! Chapman term scale height (km)
        
        real(8)              :: R           ! Chemical/dynamical term coefficient
        real(8)              :: zetaR       ! Chemical/dynamical term reference height (km)
        real(8)              :: HR          ! Chemical/dynamical term scale height (km)
        
        real(8)              :: cf(0:nsplO1+1)! Merged spline coefficients (for chemistry-dominated region of O1, NO, and (eventually), H, N)

        real(8)              :: zref        ! Reference height for hydrostatic integral and ideal gas terms
        real(8)              :: Mi(0:4)     ! Effective mass at nodes of piecewise mass profile (derived from zetaM, HML, HMU)
        real(8)              :: zetaMi(0:4) ! Height of nodes of piecewise mass profile (derived from zetaM, HML, HMU)
        real(8)              :: aMi(0:4) = 0d0    ! Slopes of piecewise mass profile segments (derived from zetaM, HML, HMU)
        real(8)              :: WMi(0:4) = 0d0   ! 2nd indefinite integral of 1/T at mass profile nodes
        real(8)              :: XMi(0:4) = 0d0   ! Cumulative adjustment to M/T integral due to changing effective mass
        real(8)              :: Izref      ! Indefinite hydrostatic integral at reference height
        real(8)              :: Tref      ! Temperature at reference height (for ideal gas law term)
        real(8)              :: zmin      ! Minimum height of profile (missing values below)
        real(8)              :: zhyd      ! Hydrostatic terms needed above this height

        integer(8)           :: massid

    end type dnparm

    contains

!==========================================================
!  A subroutine to evaluate all composition parameters that
!  are independent of altitude, but not location, time, solar
!  flux, and/or mass number.  These only need to be updated
!  if gf (and implicitly tpro) have changed.
!==========================================================

    subroutine dfnparm(massid,gf,tpro,dpro)

        use constants
        use tfn, only : tnparm
        use msis,only: eta,sfluxmod,geomag,utdep,TN,PR,N2,O2,O1,HE,H1,AR,N1,OA,NO

        implicit none
        integer,intent(in)       :: massid
        real(8),intent(in)       :: gf(0:maxnbf-1)
        type(tnparm),intent(in)  :: tpro
        type(dnparm),intent(out) :: dpro
        
        integer     :: izf, i, i1, iz
        real(8)     :: Cterm
        real(8)     :: bc(2)
        real(8)     :: hbetaL,hbetaU
        real(8)     :: delM, delz
        real(8)     :: Wi   !2nd indefinite integral at a piecewise mass profile node
    	real(8)     :: Si(-5:0,2:6) !Array of b-spline values at a mass profile node
        real(8)     :: Mzref     !Effective mass at reference altitude

        real(8), external  :: dilog

    ! ==================================================================
    ! Compute the species dependent profile parameters
    ! ==================================================================

        dpro%massid = massid

        select case(massid)

        ! Molecular Nitrogen ----------------------
        case(2)
            ! Reference mixing ratio
            dpro%lnPhiF = log(volmix(massid))
            dpro%lndref = tpro%lndtotF + dpro%lnPhiF
            dpro%zref = zetaF
            dpro%zmin = -1d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(N2%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(N2%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(N2%beta(0:mbf,3),gf(0:mbf))
            ! Photochemical correction
            dpro%R     = dot_product(N2%beta(0:mbf,7),gf(0:mbf))
            dpro%zetaR = dot_product(N2%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(N2%beta(0:mbf,9),gf(0:mbf))

        ! Molecular Oxygen ------------------------
        case(3)
            ! Reference mixing ratio
            dpro%lnPhiF = log(volmix(massid))
            dpro%lndref = tpro%lndtotF + dpro%lnPhiF
            dpro%zref = zetaF
            dpro%zmin = -1d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(O2%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(O2%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(O2%beta(0:mbf,3),gf(0:mbf))
            ! Photochemical correction
            dpro%R     = dot_product(O2%beta(0:mbf,7),gf(0:mbf))
            dpro%R     = dpro%R + geomag(O2%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%zetaR = dot_product(O2%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(O2%beta(0:mbf,9),gf(0:mbf))

        ! Atomic Oxygen --------------------------
        case(4)
            ! Reference number density
            dpro%lnPhiF = 0d0
            dpro%lndref = dot_product(O1%beta(0:mbf,0),gf(0:mbf))
            dpro%zref = zetarefO1
            dpro%zmin = nodesO1(3)
            dpro%zhyd = zetarefO1
            ! Effective mass
            dpro%zetaM = dot_product(O1%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(O1%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(O1%beta(0:mbf,3),gf(0:mbf))
            ! Chapman correction
            dpro%C     = dot_product(O1%beta(0:mbf,4),gf(0:mbf))
            dpro%zetaC = dot_product(O1%beta(0:mbf,5),gf(0:mbf))
            dpro%HC    = dot_product(O1%beta(0:mbf,6),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(O1%beta(0:mbf,7),gf(0:mbf))
            dpro%R     = dpro%R + sfluxmod(7,gf,O1,0d0)         
            dpro%R     = dpro%R + geomag(O1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%R     = dpro%R + utdep(O1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
            dpro%zetaR = dot_product(O1%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(O1%beta(0:mbf,9),gf(0:mbf))
            ! Unconstrained splines
            do izf = 0,nsplO1-1
                dpro%cf(izf) = dot_product(O1%beta(0:mbf,izf+10),gf(0:mbf))
            enddo
            !Constrained splines calculated after case statement

        ! Helium ----------------------
        case(5)
            ! Reference mixing ratio
            dpro%lnPhiF = log(volmix(massid))
            dpro%lndref = tpro%lndtotF + dpro%lnPhiF
            dpro%zref = zetaF
            dpro%zmin = -1d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(HE%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(HE%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(HE%beta(0:mbf,3),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(HE%beta(0:mbf,7),gf(0:mbf))
            dpro%R     = dpro%R + sfluxmod(7,gf,HE,1d0)         
            dpro%R     = dpro%R + geomag(HE%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%R     = dpro%R + utdep(HE%beta(cut:cut+nut-1,7),gf(cut:cut+8))
            dpro%zetaR = dot_product(HE%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(HE%beta(0:mbf,9),gf(0:mbf))

        ! Atomic Hydrogen ----------------------
        case(6)
            ! Reference number density
            dpro%lnPhiF = 0d0
            dpro%lndref = dot_product(H1%beta(0:mbf,0),gf(0:mbf))
            dpro%zref = zetaA
            dpro%zmin = 75d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(H1%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(H1%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(H1%beta(0:mbf,3),gf(0:mbf))
            ! Chapman correction
            dpro%C     = dot_product(H1%beta(0:mbf,4),gf(0:mbf))
            dpro%zetaC = dot_product(H1%beta(0:mbf,5),gf(0:mbf))
            dpro%HC    = dot_product(H1%beta(0:mbf,6),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(H1%beta(0:mbf,7),gf(0:mbf))
            dpro%R     = dpro%R + sfluxmod(7,gf,H1,0d0)        
            dpro%R     = dpro%R + geomag(H1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%R     = dpro%R + utdep(H1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
            dpro%zetaR = dot_product(H1%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(H1%beta(0:mbf,9),gf(0:mbf))

        ! Argon ----------------------
        case(7)
            ! Reference mixing ratio
            dpro%lnPhiF = log(volmix(massid))
            dpro%lndref = tpro%lndtotF + dpro%lnPhiF
            dpro%zref = zetaF
            dpro%zmin = -1d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(AR%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(AR%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(AR%beta(0:mbf,3),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(AR%beta(0:mbf,7),gf(0:mbf))
            dpro%R     = dpro%R + geomag(AR%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%R     = dpro%R + utdep(AR%beta(cut:cut+nut-1,7),gf(cut:cut+8))
            dpro%zetaR = dot_product(AR%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(AR%beta(0:mbf,9),gf(0:mbf))

!        ! Atomic Nitrogen ----------------------
        case(8)
            ! Reference number density
            dpro%lnPhiF = 0d0
            dpro%lndref = dot_product(N1%beta(0:mbf,0),gf(0:mbf))
            dpro%lndref = dpro%lndref + sfluxmod(0,gf,N1,0d0)         
            dpro%lndref = dpro%lndref + geomag(N1%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%lndref = dpro%lndref + utdep(N1%beta(cut:cut+nut-1,0),gf(cut:cut+8))
            dpro%zref = zetaF
            dpro%zmin = 110d0
            dpro%zhyd = zetaF
            ! Effective mass
            dpro%zetaM = dot_product(N1%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(N1%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(N1%beta(0:mbf,3),gf(0:mbf))
            ! Chapman correction
            dpro%C     = dot_product(N1%beta(0:mbf,4),gf(0:mbf))
            dpro%zetaC = dot_product(N1%beta(0:mbf,5),gf(0:mbf))
            dpro%HC    = dot_product(N1%beta(0:mbf,6),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(N1%beta(0:mbf,7),gf(0:mbf))
            dpro%zetaR = dot_product(N1%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(N1%beta(0:mbf,9),gf(0:mbf))

!        ! Anomalous Oxygen ----------------------
        case(9)
            dpro%lndref = dot_product(OA%beta(0:mbf,0),gf(0:mbf))
            dpro%lndref = dpro%lndref + geomag(OA%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
            dpro%zref = zetarefOA
            dpro%zmin = 110d0
            dpro%zhyd = 0d0
            dpro%C     = OA%beta(0,4)
            dpro%zetaC = OA%beta(0,5)
            dpro%HC    = OA%beta(0,6)
            return              !No further parameters needed for legacy anomalous oxygen profile

!        ! Nitric Oxide ----------------------
        case(10)
            ! Reference number density
            dpro%lnPhiF = 0d0
            dpro%lndref = dot_product(NO%beta(0:mbf,0),gf(0:mbf))
            dpro%zref = zetarefNO
            dpro%zmin = nodesNO(3)
            dpro%zhyd = zetarefNO
            ! Effective mass
            dpro%zetaM = dot_product(NO%beta(0:mbf,1),gf(0:mbf))
            dpro%HML   = dot_product(NO%beta(0:mbf,2),gf(0:mbf))
            dpro%HMU   = dot_product(NO%beta(0:mbf,3),gf(0:mbf))
            ! Chapman correction
            dpro%C     = dot_product(NO%beta(0:mbf,4),gf(0:mbf))
            dpro%zetaC = dot_product(NO%beta(0:mbf,5),gf(0:mbf))
            dpro%HC    = dot_product(NO%beta(0:mbf,6),gf(0:mbf))
            ! Dynamical correction
            dpro%R     = dot_product(NO%beta(0:mbf,7),gf(0:mbf))
            dpro%zetaR = dot_product(NO%beta(0:mbf,8),gf(0:mbf))
            dpro%HR    = dot_product(NO%beta(0:mbf,9),gf(0:mbf))
            ! Unconstrained splines
            do izf = 0,nsplNO-1
                dpro%cf(izf) = dot_product(NO%beta(0:mbf,izf+10),gf(0:mbf))
            enddo
            !Constrained splines calculated after case statement

! Failsafe -----   ---------------------------
        case default
            stop 'Species not yet implemented'

        endselect
        
        !Compute piecewise mass profile values and integration terms
        dpro%zetaMi(0) = dpro%zetaM - 2d0*dpro%HML
        dpro%zetaMi(1) = dpro%zetaM - dpro%HML
        dpro%zetaMi(2) = dpro%zetaM
        dpro%zetaMi(3) = dpro%zetaM + dpro%HMU
        dpro%zetaMi(4) = dpro%zetaM + 2d0*dpro%HMU
        dpro%Mi(0) = mbar
        dpro%Mi(4) = amu(massid)
        dpro%Mi(2) = (dpro%Mi(0) + dpro%Mi(4)) / 2d0
        delM = tanh1 * (dpro%Mi(4) - dpro%Mi(0)) / 2d0
        dpro%Mi(1) = dpro%Mi(2) - delM
        dpro%Mi(3) = dpro%Mi(2) + delM
        do i = 0, 4
            i1 = i + 1
            if (i .lt. 4) dpro%aMi(i) = (dpro%Mi(i1) - dpro%Mi(i)) / (dpro%zetaMi(i1) - dpro%zetaMi(i))
            delz = dpro%zetaMi(i) - zetaB
            if (dpro%zetaMi(i) .lt. zetaB) then
                call bspline(dpro%zetaMi(i),nodes,nd+2,6,eta,Si,iz)
                dpro%WMi(i) = dot_product(tpro%gamma(iz-5:iz),Si(:,6)) + tpro%cVS*delz + tpro%cWS
            else
                dpro%WMi(i) = (0.5d0*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
                              + tpro%cVB*delz + tpro%cWB
            endif
        end do
        dpro%XMi(0) = -dpro%aMi(0) * dpro%WMi(0)
        do i = 1, 3
            dpro%XMi(i) = dpro%XMi(i-1) - dpro%WMi(i) * (dpro%aMi(i) - dpro%aMi(i-1))
        end do
        dpro%XMi(4) = dpro%XMi(3) + dpro%WMi(4) * dpro%aMi(3)

        !Calculate hydrostatic integral at reference height, and copy temperature
        if (dpro%zref .eq. zetaF) then
            Mzref = mbar
            dpro%Tref = tpro%TzetaF
            dpro%Izref = mbar * tpro%VzetaF
        else if (dpro%zref .eq. zetaB) then
            Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
            dpro%Tref = tpro%Tb0
            dpro%Izref = 0.0d0
            if ((zetaB .gt. dpro%zetaMi(0)) .and. (zetaB .lt. dpro%zetaMi(4))) then
                i = 0
                do i1 = 1, 3
                    if (zetaB .lt. dpro%zetaMi(i1)) then
                        exit
                    else
                        i = i1
                    endif
                enddo
                dpro%Izref = dpro%Izref -  dpro%XMi(i)
            else
                dpro%Izref = dpro%Izref - dpro%XMi(4)                
            endif
        else if (dpro%zref .eq. zetaA) then
            Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
            dpro%Tref = tpro%TzetaA
            dpro%Izref = Mzref * tpro%VzetaA
            if ((zetaA .gt. dpro%zetaMi(0)) .and. (zetaA .lt. dpro%zetaMi(4))) then
                i = 0
                do i1 = 1, 3
                    if (zetaA .lt. dpro%zetaMi(i1)) then
                        exit
                    else
                        i = i1
                    endif
                enddo
                dpro%Izref = dpro%Izref - (dpro%aMi(i)*tpro%WzetaA + dpro%XMi(i))
            else
                dpro%Izref = dpro%Izref - dpro%XMi(4)                
            endif
        else 
            stop 'Integrals at reference height not available'
        endif

        !C1 constraint for O1 at 85 km
        if (massid .eq. 4) then
            Cterm = dpro%C*exp(-(dpro%zref - dpro%zetaC)/dpro%HC)
            bc(1) = dpro%lndref - Cterm - dpro%cf(7)*c1o1adj(1)      !Reference density, Chapman term, and subtraction of last unconstrained spline(7)
            bc(2) = -Mzref*g0divR/tpro%tzetaA &   !Gradient of hydrostatic term
                    -tpro%dlntdzA &  !Gradient of ideal gas law term
                    +Cterm/dpro%HC & !Gradient of Chapman term
                    -dpro%cf(7)*c1o1adj(2)  !Subtraction of gradient of last unconstrained spline(7)
            ! Compute coefficients for constrained splines
            dpro%cf(8:9) = matmul(bc,c1o1)
        endif
        
        !C1 constraint for NO at 122.5 km
        if (massid .eq. 10) then
            Cterm = dpro%C*exp(-(dpro%zref - dpro%zetaC)/dpro%HC)
            bc(1) = dpro%lndref - Cterm - dpro%cf(7)*c1noadj(1)      !Reference density, Chapman term, and subtraction of last unconstrained spline(7)
            bc(2) = -Mzref*g0divR/tpro%tb0 &   !Gradient of hydrostatic term
                    -tpro%tgb0/tpro%tb0 &  !Gradient of ideal gas law term
                    +Cterm/dpro%HC & !Gradient of Chapman term
                    -dpro%cf(7)*c1noadj(2)  !Subtraction of gradient of last unconstrained spline(7)
            ! Compute coefficients for constrained splines
            dpro%cf(8:9) = matmul(bc,c1no)
        endif
        
    return

    end subroutine dfnparm

! =======================================================================
! Evaluate density from the profile parameters
! =======================================================================

    real(8) function dfnx(z,tnz,lndtotz,Vz,Wz,tpro,dpro)

        use constants
        use msis, only      : etaO1, etaNO
        use tfn,only        : tnparm

        implicit none
        real(8),intent(in)      :: z
        real(8),intent(in)      :: tnz,lndtotz   !Temperature, total number density at input z 
        real(8),intent(in)      :: Vz,Wz   !First and second indefinite integrals of 1/T at z
        type(tnparm),intent(in) :: tpro
        type(dnparm),intent(in) :: dpro

        integer(4)          :: i,i1,iz
        real(8)             :: Mz
        real(8)             :: Sz(-5:0,2:6)
        real(8)             :: Ihyd   !Hydrostatic definite integral


        ! Below minimum height of profile
        if (z .lt. dpro%zmin) then
            dfnx = 1.0d-32
            return
        endif
        
        ! Anomalous Oxygen (legacy MSISE-00 formulation)
        if (dpro%massid .eq. 9) then
            dfnx = dpro%lndref - (z - dpro%zref)/HOA - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC)
            dfnx = exp(dfnx)
            return               !No further calculation needed for anomalous oxygen
        endif
        
        ! Below height where hydrostatic terms are needed
        if (z .lt. dpro%zhyd) then
            select case(dpro%massid)
            case(2,3,5,7)               !For N2, O2, He, and Ar, apply mixing ratios and exit
                dfnx = exp(lndtotz + dpro%lnPhiF)
                return
            case(4)                     !For O, evaluate splines
                call bspline(z,nodesO1,ndO1,4,etaO1,Sz,iz)
                dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
                return
            case(10)                     !For NO, evaluate splines
                call bspline(z,nodesNO,ndNO,4,etaNO,Sz,iz)
                dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
                return
            endselect
        endif
        
        ! Calculate hydrostatic term and apply to reference density
        Mz = pwmp(z,dpro%zetaMi,dpro%Mi,dpro%aMi)
        Ihyd = Mz * Vz - dpro%Izref
        if ((z .gt. dpro%zetaMi(0)) .and. (z .lt. dpro%zetaMi(4))) then
            i = 0
            do i1 = 1, 3
                if (z .lt. dpro%zetaMi(i1)) then
                    exit
                else
                    i = i1
                endif
            enddo
            Ihyd = Ihyd - (dpro%aMi(i)*Wz + dpro%XMi(i))
        else if (z .ge. dpro%zetaMi(4)) then
            Ihyd = Ihyd - dpro%XMi(4)                
        endif
        dfnx = dpro%lndref - Ihyd * g0divR
        
        ! Apply Chapman and logistic corrections
        select case(dpro%massid)
        case(2,3,5,7)               !For N2, O2, He, and Ar, logistic correction only
            dfnx = dfnx + dpro%R*(1+tanh((z-dpro%zetaR)/(2d0*dpro%HR)))/2d0
        case(4,6,8,10)                   !For O, H, N, NO Chapman and logistic corrections
            dfnx = dfnx - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC) &
                        + dpro%R*(1+tanh((z-dpro%zetaR)/(2d0*dpro%HR)))/2d0
        endselect

        ! Apply ideal gas law                
        dfnx = exp(dfnx) * dpro%Tref/tnz

        return

        end function dfnx

        ! -----------------------------------------------------
        ! Piecewise mass profile interpolation
        ! -----------------------------------------------------

        real(8) function pwmp(z,zm,m,dmdz)
            real(8), intent(in)  :: z
            real(8), intent(in)  :: zm(0:4)
            real(8), intent(in)  :: m(0:4)
            real(8), intent(in)  :: dmdz(0:3)
            integer              :: irng  !Index of piecwise interval
            integer              :: inode

            ! Most probable case
            if (z .ge. zm(4)) then
                pwmp = m(4)
                return
            endif
            ! Second most probable case
            if (z .le. zm(0)) then
                pwmp = m(0)
                return
            endif
            ! None of the above
            do inode = 0,3
              if (z .lt. zm(inode+1)) then
                pwmp = m(inode) + dmdz(inode)*(z - zm(inode))
                return
              endif
            enddo
            ! If we are here this is a problem
            stop 'Error in pwmp'

        end function pwmp

end module dfn


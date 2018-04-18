MODULE MOXTHProperties
  !+-------------------------------------------------------------------+
  !|                                                                   |
  !|  MOXTHproperties.f90                                              |
  !|                                                                   |
  !|  Data from NUREG/CR-6150 (EGG-2720) Volume 1V, MATPRO Library of  |
  !|  Material Properties, INEL 1993                                   |
  !|                                                                   |
  !+-------------------------------------------------------------------+
  USE DefineKinds, ONLY: Ikind, Dkind
  IMPLICIT NONE
  !
  REAL (Dkind), PARAMETER :: K1_U=296.7, K2_U=2.43D-2, K3_U=8.745D7,   &
       &   Theta_U=535.285, ED_U=1.577D5, K1_P=3.95D-4, K2_P=3.95D-4,       &
       &   K3_P=3.860D7, Theta_P=571.000, ED_P=1.967D5, R=8.3143,           &
       &   kB=1.38D-23, L1_U=1.0D-5, L2_U=3.0D-3, L3_U=4.0D-2,              &
       &   Gd_U=6.9D-20,L1_P=9.0D-6, L2_P=2.7D-3, L3_P=7.0D-2,              &
       &   Gd_P=7.0D-20
  !
CONTAINS
  !
  ! define T/H properties of Mixed Uranium-Plutonium Oxide functions
  !
  REAL (Dkind) FUNCTION Tm_MOX(FCOMP,FBu)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, FBu
    !
    ! FCOMP = weight fraction of PuO2 in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! Tm_MOX = Temp. solid MOX begins to melt (T-liquidus) (Deg K)
    Tm_MOX = 3.11315D3 + FCOMP*(-3.21660D2 + FCOMP*(-1.448518D2))    &
         &         -3.2D-3*FBu
    !
    RETURN
  END FUNCTION Tm_MOX

  REAL (Dkind) FUNCTION Tf_MOX(FCOMP,FBu)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, FBu
    !
    ! FCOMP = weight fraction of PuO2 in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! Tf_MOX = Temp. liquid MOX begins to freeze (T-solidus) (Deg K)
    Tf_MOX = 3.11315D3 + FCOMP*(-5.41395D2 + FCOMP*(7.468390D1))     &
         &           -3.2D-3*FBu
    !
    RETURN
  END FUNCTION Tf_MOX

  REAL (Dkind) FUNCTION Cp_MOX(FCOMP,Y,FBu,X,T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, Y, FBu, T, T1, FCPUO2, FCPPUO2, X
    !
    ! FCOMP = weight fraction of PuO2 in MOX
    ! Y = oxygen to metal ratio in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! X is fraction melted
    !
    ! Temperature dependant function for MOX specific heat
    ! in units of J/(Kg DegK) for Valid Range: 30K < T < 4000 Deg K
    !
    IF(T .LT. Tm_MOX(FCOMP,FBU)) THEN
       T1 = T
    ELSE
       T1 = Tm_MOX(FCOMP,FBU)
    END IF
    !  values if less than melt temperature of solid MOX (Deg K)
    FCPUO2 = K1_U*(Theta_U**2)*EXP(Theta_U/T1)/                      &
         &       (T1*T1*(EXP(Theta_U/T1)-1.0D0)**2) + K2_U*T1 +           &
         &       ((Y*K3_U*ED_U)/(2.0D0*R*T1*T1))*EXP(-ED_U/(R*T1))
    FCPPUO2 = K1_P*(Theta_P**2)*EXP(Theta_P/T1)/                     &
         &       (T1*T1*(EXP(Theta_P/T1)-1.0D0)**2) + K2_P*T1 +           &
         &       ((Y*K3_P*ED_P)/(2.0D0*R*T1*T1))*EXP(-ED_P/(R*T1))
    IF(T .GE. Tf_MOX(FCOMP,FBu)) THEN
       ! greater than freeze temperature of liquid MOX (deg K)
       FCPUO2 = 503.0D0
       FCPPUO2 = 503.0D0
    END IF
    IF(T .LT. Tm_MOX(FCOMP,FBu) .OR. T .GE. Tf_MOX(FCOMP,FBu)) THEN
       Cp_MOX = FCPUO2*(1.0D0-FCOMP)+FCPPUO2*FCOMP
    ELSE
       ! between melt and freeze temperatures some melting
       Cp_MOX = (1.0D0 -X)*(FCPUO2*(1.0D0-FCOMP)+FCPPUO2*FCOMP)       &
            &        + X*503.0D0
    END IF
    !
    RETURN
  END FUNCTION Cp_MOX

  REAL (Dkind) FUNCTION H_MOX(FCOMP,Y,FBu,X,T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, Y, FBu, T, T1, HUO2, HPUO2, X, HMIX
    !
    ! FCOMP = weight fraction of PuO2 in MOX
    ! Y = oxygen to metal ratio in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! X is fraction melted
    !
    ! Temperature dependant function for MOX specific enthalpy
    ! in units of J/Kg for Valid Range: 30K < T < 4000 Deg K
    !
    IF(T .LT. Tm_MOX(FCOMP,FBU)) THEN
       T1 = T
    ELSE
       T1 = Tm_MOX(FCOMP,FBU)
    END IF
    HUO2 = K1_U*Theta_U/(EXP(Theta_U/T1)-1.0D0) + K2_U*T1*T1/2.0D0 + &
         &       (Y/2.0D0)*(K3_U*EXP(-ED_U/(R*T1)))
    HPUO2 = K1_P*Theta_P/(EXP(Theta_P/T1)-1.0D0)+ K2_P*T1*T1/2.0D0 + &
         &       (Y/2.0D0)*(K3_P*EXP(-ED_P/(R*T1)))
    HMIX = HUO2*(1.0D0-FCOMP)+HPUO2*FCOMP
    IF(T .LT. Tm_MOX(FCOMP,FBu)) THEN
       H_MOX = HMIX
    ELSE IF((T-Tm_MOX(FCOMP,FBu)) .LT. 1.0D-2) THEN
       H_MOX = HMIX + X*2.75D5
    ELSE IF(T .GT. Tm_MOX(FCOMP,FBu)) THEN
       H_MOX = HMIX +2.75D5+503.0D0*(T-Tm_MOX(FCOMP,FBu))
    END IF
    !
    RETURN
  END FUNCTION H_MOX

  REAL (Dkind) FUNCTION Epsilon_MOX(FCOMP,FBu,X,T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, FBu, X, T, Epsilon_UO2, Epsilon_PUO2, DelT
    !
    ! X = fraction of UO2/PUO2 melted
    ! FCOMP = weight fraction of PuO2 in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! T = Temperature (Deg K)
    ! Temperature dependant function for MOX Thermal
    ! Expansion strain in m/m  Valid Range: 300 < T < 3500K
    IF(T .LT. Tm_MOX(FCOMP,FBU)) THEN
       Epsilon_UO2 = L1_U*T - L2_U + L3_U*EXP(-Gd_U/(kB*T))
       Epsilon_PUO2 = L1_P*T - L2_P + L3_P*EXP(-Gd_P/(kB*T))
    ELSE IF((T-Tm_MOX(FCOMP,FBu)) .LT. 1.0D-2) THEN
       Epsilon_UO2 = L1_U*T - L2_U + L3_U*EXP(-Gd_U/(kB*T))           &
            &             + 4.3D-2*X
       Epsilon_PUO2 = L1_P*T - L2_P + L3_P*EXP(-Gd_P/(kB*T))          &
            &             + 4.3D-2*X
    ELSE IF(T .GT. Tm_MOX(FCOMP,FBu)) THEN
       DelT = Tm_MOX(FCOMP,FBU)-Tf_MOX(FCOMP,FBU)
       Epsilon_UO2 = L1_U*T - L2_U + L3_U*EXP(-Gd_U/(kB*T))           &
            + 4.3D-2 +3.6D-5*(T-(Tm_MOX(FCOMP,FBU)+DelT))
       Epsilon_PUO2 = L1_P*T - L2_P + L3_P*EXP(-Gd_P/(kB*T))          &
            + 4.3D-2 +3.6D-5*(T-(Tm_MOX(FCOMP,FBU)+DelT))
    END IF
    Epsilon_MOX = Epsilon_UO2*(1.0D0-FCOMP)+Epsilon_PUO2*FCOMP +     &
         &            X*4.3D-2
    !
    RETURN
  END FUNCTION  Epsilon_MOX

  REAL (Dkind) FUNCTION Rho_UO2(FCOMP,FBu,X,T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP,FBu, X, T
    !
    ! T = Temperature (Deg K)
    ! Epsilon_UO2 is thermal strain in UO2
    ! Temperature dependant function for UO2 theoretical density
    ! density units are MT/m^3 = g/cc =10^3 Kg/m^3
    Rho_UO2 = 10.980D+3*(1.0D0 - 3.0D0*Epsilon_MOX(FCOMP,FBu,X,T))
    !
    RETURN
  END FUNCTION  Rho_UO2

  REAL (Dkind) FUNCTION Cv_MOX(FCOMP,T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) FCOMP, T, FCVUO2, FCVPUO2
    !
    ! FCOMP = weight fraction of PuO2 in MOX
    !
    ! Temperature dependant function for MOX specific heat at
    ! constant volume along saturation linein units of J/(Kg DegK)
    ! for Valid Range: 30K < T < 4000 Deg K
    FCVUO2 = K1_U*(Theta_U**2)*EXP(Theta_U/T)/                       &
         &       (T*T*(EXP(Theta_U/T)-1.0D0)**2)
    FCVPUO2 = K1_P*(Theta_P**2)*EXP(Theta_P/T)/                      &
         &       (T*T*(EXP(Theta_P/T)-1.0D0)**2)
    !
    Cv_MOX = FCVUO2*(1.0D0-FCOMP)+FCVPUO2*FCOMP
    !
    RETURN
  END FUNCTION Cv_MOX

  REAL (Dkind) FUNCTION K_MOX(D, FCOMP,FBu, Y, X, T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) D, FCOMP, FBu, Y, X, T, AA, BB, TRM1, TRM2,        &
         &            TRM2a, TRM2b
    !
    ! D = fration of theoretical density
    ! FCOMP = weight fraction of PuO2 in MOX
    ! FBu = fuel burnup (MWd/tU)
    ! Y = oxygen to metal ratio in MOX
    ! X = fraction melted
    ! T = Temperature (Deg K)
    !
    ! Temperature dependant function for MOX Thermal
    ! Conductivity in W/(m-Deg.K). Valid Range:  T < 2300.K
    !
    IF (T .LT. 1364.0D0) THEN
       TRM1 = D/(1.0D0+(6.5D0-4.69D-3*T)*(1.0D0-D))
    ELSE IF (T .LT. 1834.0D0) THEN
       TRM1= 1.0D0+((D/(1.0D0+(6.5D0-4.69D-3*T)*(1.0D0-D)))-1.0D0)*   &
            &     (T-1834.0D0)/(1364.0D0-1834.0D0)
    ELSE
       TRM1 =1.0D0
    END IF
    AA = 0.399D0+12.6D0*(2.0D0-Y)
    BB = 6.867D-2*(1.0D0-0.6238D0*FCOMP)
    !

    IF (T .LT. 1800.0D0) THEN
       TRM2 = Cv_MOX(FCOMP,T)/                                        &
            &      ((AA+BB*T)*(1.0D0+3.0D0*Epsilon_MOX(FCOMP,FBu,X,T)))
    ELSE IF (T .LT. 2300.0D0) THEN
       TRM2a= Cv_MOX(FCOMP,1800.0D0)/((AA+BB*1800.0D0)*               &
            &      (1.0D0+3.0D0*Epsilon_MOX(FCOMP,FBu,X,1800.0D0)))
       TRM2b= Cv_MOX(FCOMP,2050.0D0)/((AA+BB*2050.0D0)*               &
            &      (1.0D0+3.0D0*Epsilon_MOX(FCOMP,FBu,X,2050.0D0)))
       TRM2= TRM2a+(TRM2b-TRM2a)*(T-1800.0D0)/(500.0D0)
    ELSE
       TRM2= Cv_MOX(FCOMP,2050.0D0)/((AA+BB*2050.0D0)*               &
            &      (1.0D0+3.0D0*Epsilon_MOX(FCOMP,FBu,X,2050.0D0)))
    END IF
    !
    K_MOX = TRM1*TRM2 + (5.2997D-3)*T*EXP(-13358.0D0/T)*             &
         &      (1.0D0 + 1.69D-1*((13358.0D0/T)+2.0D0)**2)
    !
    RETURN
  END FUNCTION K_MOX
  !
END MODULE MOXTHProperties

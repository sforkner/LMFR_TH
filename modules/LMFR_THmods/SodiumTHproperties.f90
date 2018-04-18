MODULE SodiumTHProperties
  ! +----------------------------------------------------------------+
  ! |                                                                |
  ! |  Module SodiumTHProperties.f90                                 |
  ! |                                                                |
  ! |  Data from:                                                    |
  ! |  ANL/RE-95/2 "Thermodynamic and Transport Properties of Sodium |
  ! |  Liquid and Vapor", J. K. Fink and L. Leibowitz, January, 1995.|
  ! |                                                                |
  ! |                                                                |
  ! +----------------------------------------------------------------+
  USE DefineKinds, ONLY: Ikind, Dkind
  IMPLICIT NONE
  !
  ! Sodium critical temperature (Deg K), critical Pressure (MPa),
  ! critical density (kg/m**3), Solid sodium melting temperature
  ! (Deg K) and liquid sodium boiling Point (Deg K) at 1 Atm pressure
  REAL (Dkind), PARAMETER :: Tc_Na=2503.7, Pc_Na=25.64,              &
       &                     Rhoc_Na=219.0, Tm_Na=371.0, Tb_Na=1154.7
  !
CONTAINS
  !
  ! define T/H properties of Sodium coolant functions
  !
  REAL (Dkind) FUNCTION Hfg_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, DHg
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium enthalpy
    ! of vaporization (kJ/kg)
    Hfg_Na = 393.37D0*(1.0D0 - T/Tc_Na) +                          &
         &        4398.6*(1.0D0 - T/Tc_Na)**0.29302D0
    !
    RETURN
  END FUNCTION Hfg_Na

  REAL (Dkind) FUNCTION Hf_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium enthalpy
    ! relative to zero for solid sodium at 298.15 K in units of
    ! kJ/Kg for Valid Range: 371K < T < 2503K
    IF (Tm_Na .LE. T .AND. T .LE. 2000.0D0) THEN
       Hf_Na = -365.77D0 + T*(1.6582D0 + T*(-4.2395D4 + T*(1.4847D7)))&
            &        + 2992.6D0/T
    ELSE IF (2000.0D0 .LT. T .AND. T .LE. Tc_Na) THEN
       Hf_Na = 2128.4D0+0.86496D0*T - Hfg_Na(T)/2.0D0
    ELSE
       STOP ' Sodium Temperature out valid of range********'
    END IF
    !
    RETURN
  END FUNCTION Hf_Na

  REAL (Dkind) FUNCTION Hg_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Sodium vapor enthalpy
    ! relative to zero for solid sodium at 298.15 K in units of
    ! kJ/Kg for Valid Range: 371K < T < 2503K
    IF (Tm_Na .LE. T .AND. T .LE. 2000.0D0) THEN
       Hg_Na = Hf_Na(T) + Hfg_Na(T)
    ELSE IF (2000.0D0 .LT. T .AND. T .LE. Tc_Na) THEN
       Hg_Na = 2128.4D0+0.86496D0*T + Hfg_Na(T)/2.0D0
    ELSE
       STOP ' Sodium Temperature out valid of range********'
    END IF
    !
    RETURN
  END FUNCTION Hg_Na

  REAL (Dkind) FUNCTION Rhol_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium density
    ! in units of Kg/m**3 valid range: 371K < T < 2503.7K
    Rhol_Na = Rhoc_Na + 275.32D0*(1.0D0 - T/Tc_Na)                   &
         &                 + 511.58D0*SQRT(1.0D0 - T/Tc_Na)
    !
    RETURN
  END FUNCTION Rhol_Na

  REAL (Dkind) FUNCTION GammaSigma_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium temperature
    ! derivative with pressure along saturation curve
    GammaSigma_Na = (1.263373D4/(T*T) -0.4672D0/T)                  &
         &               *EXP(11.9463D0 - 1.263373D4/T - 0.4672D0*LOG(T))
    !
    RETURN
  END FUNCTION GammaSigma_Na

  REAL (Dkind) FUNCTION Rhog_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Sodium vapor density
    ! in units of Kg/m**3 valid range: 371K < T < 2503.7K
    Rhog_Na = 1.0D0/(Hfg_Na(T)/(T*GammaSigma_Na(T))                 &
         &             +  1.0D0/Rhol_Na(T))
    !
    RETURN
  END FUNCTION Rhog_Na

  REAL (Dkind) FUNCTION DHDT_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium temperature
    ! derivative od enthalpy with temperaturealong saturation curve
    IF (T .GE. Tm_Na .AND. T .LE. 2000.D0) THEN
       DHDT_Na =  1.6582D0 + T*(-8.4790D-4 + T*(4.4541D-7))           &
            &         - 2.9926D3/(T*T)
    ELSE IF (T .GE. 2000.D0 .AND. T.LE. Tc_Na) THEN
       DHDT_Na = 0.86496D0-0.5D0*(-0.157115D0-0.514789D0/             &
            &         (1.0D0-3.99409D-4*T)**0.70698D0)
    ELSE
       STOP ' Sodium Temperature out valid of range********'
    END IF
    !
    RETURN
  END FUNCTION DHDT_Na

  REAL (Dkind) FUNCTION CSigma_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for derivative Sodium enthalpy
    ! along saturation curve
    CSigma_Na = DHDT_Na(T) - GammaSigma_Na(T)/Rhol_Na(T)
    !
    RETURN
  END FUNCTION CSigma_Na

  REAL (Dkind) FUNCTION AlphaSigma_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, Anumer, Adenom
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium coefficeint of
    ! Thermal Expansion along saturation curve
    Anumer = 0.109965D0 + 0.102165D0/(1.0D0-3.99409D-4*T)**0.5
    Adenom = 494.32D0+511.58D0*(1.0D0-3.99409D-4*T)**0.5 -  0.109965D0*T
    AlphaSigma_Na = Anumer/Adenom
    !
    RETURN
  END FUNCTION AlphaSigma_Na

  REAL (Dkind) FUNCTION BetaS_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, ThetaS_Na
    REAL (Dkind), PARAMETER :: BetaSM=1.717D-4
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Sodium adiabatic
    ! compressibility along saturation curve
    ThetaS_Na = (T-Tm_Na)/(Tc_Na-Tm_Na)
    BetaS_Na = BetaSM*(1.0D0+ThetaS_Na/3.2682D0)/(1.0D0-ThetaS_Na)
    !
    RETURN
  END FUNCTION BetaS_Na

  REAL (Dkind) FUNCTION BetaT_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, Anumer, Adenom
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Sodium isothermal
    ! compressibility along saturation curve
    Anumer = BetaS_Na(T)*CSigma_Na(T)+(AlphaSigma_Na(T)/Rhol_Na(T))* &
         &        (AlphaSigma_Na(T)+BetaS_Na(T)*GammaSigma_Na(T))
    Adenom = CSigma_Na(T)-(GammaSigma_Na(T)/Rhol_Na(T))* &
         &        (AlphaSigma_Na(T)+BetaS_Na(T)*GammaSigma_Na(T))
    BetaT_Na = Anumer/Adenom
    !
    RETURN
  END FUNCTION BetaT_Na

  REAL (Dkind) FUNCTION AlphaP_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, ThetaS_Na
    REAL (Dkind), PARAMETER :: BetaSM=1.717D-4
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for saturate Sodium liquid
    ! thermal expansion coefficient
    AlphaP_Na = AlphaSigma_Na(T)+ BetaT_Na(T)*GammaSigma_Na(T)
    !
    RETURN
  END FUNCTION AlphaP_Na

  REAL (Dkind) FUNCTION Cp_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Liquid Sodium specific heat
    ! capacity at constant pressure in units kJ/(Kg DegK)
    ! valid range:371K < T < 2500K
    Cp_Na = CSigma_Na(T)+(T*AlphaP_Na(T)*GammaSigma_Na(T))/Rhol_Na(T)
    !
    RETURN
  END FUNCTION Cp_Na

  REAL (Dkind) FUNCTION OldCp_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Liquid Sodium specific heat
    ! capacity at constant pressure in units kJ/(Kg DegK)
    ! valid range:371K < T < 2500K
    OldCp_Na = 1.6582D0 + T*(-8.4790D-4 + T*(4.4541D-7)) -          &
         &          2.9926D3/(T*T)
    !
    RETURN
  END FUNCTION OldCp_Na

  REAL (Dkind) FUNCTION Pvsat_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Sodium vapor pressure (MPa)
    ! in equilibrium with saturated liquid
    Pvsat_Na = EXP(11.9463D0 - 1.263373D4/T - 0.4672D0*LOG(T))
    !
    RETURN
  END FUNCTION Pvsat_Na

  REAL (Dkind) FUNCTION v_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium speed of
    ! sound (m/s)
    IF (T .GE. Tm_Na .AND. T .LE. 1773.0D0) THEN
       v_Na = 2660.7D0 + T*(-0.37667D0 + T*(-9.0356D-5))
    ELSE IF (T .GT. 1773.0D0 .AND. T .LE. Tc_Na) THEN
       v_Na = 1000.0D0/(Rhol_Na(T)*BetaS_Na(T))**0.5
    ELSE
       STOP ' Sodium Temperature out valid of range********'
    END IF
    !
    RETURN
  END FUNCTION v_Na

  REAL (Dkind) FUNCTION SurfTen_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg L)
    ! Temperature dependant function for Liquid Sodium surface tension
    ! in units mN/m
    ! valid range:371K < T < 1600K
    SurfTen_Na = 240.5D0*(1.0D0 - T/Tc_Na)**1.126
    !
    RETURN
  END FUNCTION SurfTen_Na

  REAL (Dkind) FUNCTION K_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium Thermal
    ! Conductivity in W/(m-Deg.K). Valid Range: 371K < T < 1500K
    K_Na = 124.67D0 + T*(-0.11381D0 + T*(5.5226D-5 + T*(-1.1842D-8)))
    !
    RETURN
  END FUNCTION K_Na

  REAL (Dkind) FUNCTION Vis_Na(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for liquid Sodium Dynamic
    ! Viscosity in (Pa s). Valid Range: 371K < T < 2500K
    Vis_Na = EXP(-6.4406D0 -0.3958D0*LOG(T) +556.835D0/T)
    !
    RETURN
  END FUNCTION Vis_Na

  REAL (Dkind) FUNCTION Tsat_Na()
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    !
    ! Liquid sodium boiling temperature (Deg K)
    Tsat_Na = Tb_Na
    !
    RETURN
  END FUNCTION Tsat_Na
  !
END MODULE SodiumTHProperties

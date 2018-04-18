MODULE SS304THProperties
  !+-------------------------------------------------------------------+
  !|                                                                   |
  !|  SS304THproperties.f90                                            |
  !|                                                                   |
  !|  Data from NUREG/CR-6150 (EGG-2720) Volume 1V, MATPRO Library of  |
  !|  Material Properties, INEL 1993                                   |
  !|                                                                   |
  !+-------------------------------------------------------------------+
  USE DefineKinds, ONLY: Ikind, Dkind
  IMPLICIT NONE
  !
CONTAINS
  !
  ! define T/H properties of Stainless Steel 304 functions
  !
  REAL (Dkind) FUNCTION Cp_SS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Stainless Steel 304 specific
    ! heat in units of J/(Kg DegK) for Valid Range: 300K < T < 1671K
    IF (T .LT. 1558.0D0) THEN
       Cp_SS = 326.0D0 +T*(2.98D-1 + T*(-9.56D-5 ))
    ELSE IF (T .GE. 1558.0D0) THEN
       Cp_SS = 558.228D0
    END IF
    !
    RETURN
  END FUNCTION Cp_SS

  REAL (Dkind) FUNCTION H_SS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Stainless Steel 304 enthalpy
    ! relative to zero for U at 300K in units of J/Kg for Valid
    ! Range: 371K < T < 1800K
    IF (T .LT. 1558.0D0) THEN
       H_SS = T*(326.0D0 +T*(1.49D-1 + T*(-3.187D-5)))
    ELSE IF (1558.0D0 .LT. T .AND. T .LT. 1671.0D0) THEN
       H_SS = -1.206610D+5 + T*(558.228D0)
    ELSE IF (1671.0D0 .LE. T .AND. T .LE. 1727.0D0) THEN
       H_SS = -8.4756610D+5 + T*(5558.228D0)
    ELSE
       H_SS =  1.593390D+5 + T*(558.228D0)
    END IF
    !
    RETURN
  END FUNCTION H_SS

  REAL (Dkind) FUNCTION K_SS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Stainless Steel 304 Thermal
    ! Conductivity in W/(m-Deg.K). Valid Range: 300 < T < 1800.0K
    IF (T .LT. 1671.0D0) THEN
       K_SS = 7.58D0 + T*(1.89D-2)
    ELSE IF (1671.0D0 .LE. T .AND. T .LT. 1727.0D0) THEN
       K_SS = 610.9393D0 + T*(-3.421767D-1)
    ELSE
       K_SS =  20.0D0
    END IF
    !
    RETURN
  END FUNCTION K_SS

  REAL (Dkind) FUNCTION Epsilon_SS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Stainless Steel 304 Thermal
    ! Expansion strain in m/m  Valid Range: 300 < T < 1800.0K
    IF (T .LT. 1671.0D0) THEN
       Epsilon_SS = T*(1.57D-5 + T*(1.69D-9))
    ELSE IF (1671.0D0 .LE. T .AND. T .LT. 1727.0D0) THEN
       Epsilon_SS = -2.9866340D-1 + T*(1.972573D-4)
    ELSE
       Epsilon_ss = 4.20D-2
    END IF
    !
    RETURN
  END FUNCTION  Epsilon_SS

  REAL (Dkind) FUNCTION Rho_SS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind), PARAMETER :: SS300=7.90D+3
    REAL (Dkind) T
    !
    ! T = Temperature (Deg K)
    ! Temperature dependant function for Stainless Steel Thermal
    ! Expansion strain in m/m  Valid Range: 300 < T < 1800.0K
    ! density is value of Stainless Steel density at 300K (i.e.
    ! 7.9MT/m**3)divided by the cubic power of one plus thermal
    ! strain at temperature T (strain is zero at 300K)
    ! density units are MT/m^3 = g/cc =10^3 Kg/m^3
    Rho_SS = SS300/(1.0D0 + Epsilon_SS(T))**3
    !
    RETURN
  END FUNCTION  Rho_SS

  REAL (Dkind) FUNCTION Tmelt_SS()
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    !
    ! melting point of solid Stainless Steel temperature (Deg K)
    Tmelt_SS = 1727.0D0
    !
    RETURN
  END FUNCTION Tmelt_SS

END MODULE SS304THProperties

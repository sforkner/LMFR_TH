MODULE THProperties
  ! +----------------------------------------------------------------+
  ! |                                                                |
  ! |  Module THProperties.f90                                       |
  ! |                                                                |
  ! +----------------------------------------------------------------+
  USE DefineKinds, ONLY: Ikind, Dkind
  USE SodiumTHProperties, ONLY: K_Na, VIS_Na, Rhol_Na, Cp_Na, Tsat_Na
  USE UmetalTHProperties, ONLY: K_U
  USE SS304THProperties, ONLY:  K_SS
  USE MOXTHProperties, ONLY: K_MOX
  IMPLICIT NONE
  !
CONTAINS
  !
  ! define T/H properties of coolant functions
  !
  REAL (Dkind) FUNCTION KFU(T)
    USE DefineKinds, ONLY: Ikind, Dkind, Lkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    LOGICAL (Lkind) UO2
    !
    !UO2 = .TRUE.
    UO2 = .FALSE.
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for fuel conductivity
    IF (UO2) THEN
       ! K_MOX is W/(m-Deg K) convert to BTU/(Hr ft-Deg F) by
       ! multiplying by 0.578176
       KFU = K_MOX(0.96D0, 0.0D0,0.0D0, 2.0D0,0.0D0, TK)*0.578176D0
    ELSE
       ! K_U is W/(m-Deg K) convert to BTU/(Hr ft-Deg F) by
       ! multiplying by 0.578176
       KFU = K_U(TK)*0.578176D0
    END IF
    !
    RETURN
  END FUNCTION KFU

  REAL (Dkind) FUNCTION KCL(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for clad conductivity
    ! K_SS is W/(m-Deg K) convert to BTU/(Hr ft-Deg F) by
    ! multiplying by 0.578176
    KCL = K_SS(TK)*0.578176D0
    !
    RETURN
  END FUNCTION KCL

  REAL (Dkind) FUNCTION KBC(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for coolant conductivity
    ! K_Na is W/(m-Deg K) convert to BTU/(Hr ft-Deg F) by
    ! multiplying by 0.578176
    KBC = K_Na(TK)*0.578176D0
    !
    RETURN
  END FUNCTION KBC

  REAL (Dkind) FUNCTION VIS(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for Coolant viscosity
    ! VIS_Na is (Pa s)convert to lb/(ft Hr) by multipling
    ! by 0.67192D0*3600.0D0
    VIS = VIS_Na(TK)*0.67192D0*3600.0D0
    !
    RETURN
  END FUNCTION VIS

  REAL (Dkind) FUNCTION RHO(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for Coolant density
    ! Rho_Na is Kg/m**3 convert to lb/ft**3 by multipling
    ! by 0.062428
    RHO = Rhol_Na(TK)*6.2428D-2
    !
    RETURN
  END FUNCTION RHO

  REAL (Dkind) FUNCTION CP(T)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) T, TK
    !
    ! T = Temperature (Deg F)
    ! TK = Temperature (Deg K)
    TK = (T-32.0D0)*(5.0D0/9.0D0)+273.15D0
    ! Temperature dependant function for Coolant specific heat
    ! Cp_Na is kJ/(Kg DegK) convert to BTU/(Lb DegF) by multiplying
    ! by 0.23884597
    CP = Cp_Na(TK)*0.23884597D0
    !
    RETURN
  END FUNCTION CP

  REAL (Dkind) FUNCTION TSAT(P)
    USE DefineKinds, ONLY: Ikind, Dkind
    IMPLICIT NONE
    REAL (Dkind) P
    !
    ! P = local pressure (psia)
    ! Tsat_Na returns temperature in Deg K convert to Deg F
    ! Coolant saturation temperature as function of pressure
    ! Tsat_Na returns temperature in Deg K convert to Deg F
    ! and is not a function of pressure
    TSAT = (Tsat_Na()-273.15D0)*(9.0D0/5.0D0)+32.0D0
    !
    RETURN
  END FUNCTION TSAT
  !
END MODULE THProperties

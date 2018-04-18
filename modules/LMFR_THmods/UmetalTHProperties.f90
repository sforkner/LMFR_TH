Module   UmetalTHProperties
  !+-------------------------------------------------------------------+
  !|                                                                   |
  !|  UmetalTHproperties.f90                                           |
  !|                                                                   |
  !|  Data from NUREG/CR-6150 (EGG-2720) Volume 1V, MATPRO Library of  |
  !|  Material Properties, INEL 1993                                   |
  !|                                                                   |
  !+-------------------------------------------------------------------+
  Use DefineKinds, Only: Ikind, Dkind
  Implicit None
  Real (Dkind), Parameter :: T_Uab=938.0, H_Uab=12500.0, T_Ubg=1049.0, &
  &  H_Ubg=20060.0, T_Ugl=1405.6, H_Ugl=82350.0, T_Um=1132.3
  !
  Contains
    !
    ! define T/H properties of Metalic Uranium functions
    !
    Real (Dkind) Function Epsilon_U(T)
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      Real (Dkind) T
      !
      ! T = Temperature (Deg K)
      ! Temperature dependant function for metalic Uranium Thermal
      ! Expansion strain in m/m  Valid Range: 300 < T < 1132.3K
      If (T .LT. 942.0D0) Then
        Epsilon_U = (-3.0033D-1 +T*(7.1847D-4 + T*(1.0498D-6)))/1.0D2
      Else If (942.0D0 .LE. T .AND. T .LT. 1045.0D0) Then
        Epsilon_U = (-2.8340D-1 + T*(1.9809D-3))/1.0D2
      Else If (1045.0D0 .LE. T .AND.T .LE. T_Um) Then
        Epsilon_U = (-2.7120D-1 + T*(2.2298D-3))/1.0D2
      Else
        STOP '*STOP* Metallic Uranium Temp out of range in Epsilon_U'
      End If
      Return
    End Function  Epsilon_U

    Real (Dkind) Function Rho_U(T)
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      Real (Dkind), Parameter :: U300=19.05D+3
      Real (Dkind) T
      !
      ! T = Temperature (Deg K)
      ! Temperature dependant function for metalic Uranium Thermal
      ! Expansion strain in m/m.  Density is value of metallic
      ! uranium density at 300K (i.e.19.05MT/m**3) corrected by the
      ! thermal strain Epsilon_U(T) at temperature T (strain is zero
      ! at 300K). Density units are MT/m^3 = g/cc =10^3 Kg/m^3
      !
      Rho_U = U300/(1.0D0 + Epsilon_U(T))**3
      !
      Return
    End Function  Rho_U

    Real (Dkind) Function Cp_U(T)
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      Real (Dkind) T
      !
      ! T = Temperature (Deg K)
      ! Temperature dependant function for metalic Uranium specific
      ! heat in units of J/(Kg DegK) for Valid Range: 371K < T < 1173K
      If (T .LT. T_Uab) Then
        Cp_U = 104.82D0 +T*(5.3686D-3 + T*(10.1823D-5))
      Else If (T_Uab .LE. T .AND. T .LT. T_Ubg) Then
        Cp_U = 176.41311D0
      Else If (T .GE. T_Ubg) Then
        Cp_U = 156.80756D0
      End If
      !
      Return
    End Function Cp_U

    Real (Dkind) Function H_U(T)
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      Real (Dkind) T
      !
      ! T = Temperature (Deg K)
      ! Temperature dependant function for metalic Uranium enthalpy
      ! relative to zero for U at 300K in units of J/Kg for Valid
      ! Range: 371K < T < 2000K
      If (T .LT. T_Uab) Then
        H_U = -3.255468D4 +T*(1.0466D2 + T*(2.685D-3 + T*(3.389D-5)))
      Else If (T_Uab .LE. T .AND. T .LT. T_Ubg) Then
        H_U = -5.1876776D4 + T*(1.7092D2)
      Else If (T_Ubg .LE. T .AND. T .LT. T_Ugl) Then
        H_U = -2.0567496D-4 + T*(1.602D2)
      Else
        H_U = 6.177850D5 + T*(1.602D2)
      End If
      !
      Return
    End Function H_U

    Real (Dkind) Function K_U(T)
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      Real (Dkind) T
      !
      ! T = Temperature (Deg K)
      ! Temperature dependant function for metalic Uranium Thermal
      ! Conductivity in W/(m-Deg.K). Valid Range:  T < T_Ugl
      !K_U = 20.457D0 + T*(1.2047D-2 + T*(-5.7368D-6))
      K_U = 20.457D0 + T*(1.2047D-2 + T*(5.7368D-6))
      !
      Return
    End Function K_U

    Real (Dkind) Function Tmelt_U()
      Use DefineKinds, Only: Ikind, Dkind
      Implicit None
      !
      ! melting point of metalic Uranium temperature (Deg K)
      Tmelt_U = T_Ugl
      !
      Return
    End Function Tmelt_U

End Module UmetalTHProperties

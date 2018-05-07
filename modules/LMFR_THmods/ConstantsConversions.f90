MODULE  ConstantsConversions
  !------------------------------------------------------------------
  ! Module ConstantsConversions Description:
  !   Define Globally used Constants (Parameters) and Conversions
  !  (Parameters and Functions).
  !
  ! Language: Fortran 90+
  !
  ! Author   : Samuel L. Forkner
  ! Company  : Signal Mountain Software
  ! Date     : 2018-04-18 15:58
  !------------------------------------------------------------------
  USE  DefineKinds, ONLY: Ikind, Dkind
  IMPLICIT NONE
  !
  ! PI (ratio of circle circumference to diameter)
  REAL (kind=Dkind), PARAMETER :: PI = 3.141592653589793238462_Dkind
  ! Degrees Fahrenheit at zero Degrees Celsius
  REAL (kind=Dkind), PARAMETER :: F_0C = 32.0_Dkind
  ! Ratio a degree Celsius (or Kelvin) to Degree Fahrenheit
  REAL (kind=Dkind), PARAMETER :: CperF = 0.555555555555556_Dkind
  ! Ratio a degree Fahrenheit to Degree Celsius (or Kelvin)
  REAL (kind=Dkind), PARAMETER :: FperC = 1.80000000000000_Dkind
  ! Degrees Kelvin at zero Degrees Celsius
  REAL (kind=Dkind), PARAMETER :: K_0C = 273.15_Dkind
  ! Density Conversion from Kg/m**3 to Lb/Ft**3
  REAL (kind=Dkind), PARAMETER :: DensityConvert = 6.2427960576145D-2
  ! Mass Conversion Kg to Lb
  REAL (kind=Dkind), PARAMETER :: MassConvert = 2.20462262184878_Dkind
  ! Specific Heat Conversion (kJ/(Kg deg K)) to (BTU/(lbm deg F))
  REAL (kind=Dkind), PARAMETER :: SpecificHeatC = 0.238845896627496_Dkind
  ! Thermal Conductivity Conversion from (w/(m deg K) to BTU/(Hr Ft deg F))
  REAL (kind=Dkind), PARAMETER :: TconductivityC = 0.577789316542998_Dkind
  ! Dynamic Viscosity conversion from Kg/(m-s) to lbm/(ft-Hr)
  REAL (kind=Dkind), PARAMETER :: DynViscosityC = 2419.08831050222_Dkind
  ! Pressure Atmospheres to psia
  REAL (kind=Dkind), PARAMETER :: PressureC = 14.6959487755134_Dkind
  ! Acceleration due to Gravity (g) ft/(sec**2)
  REAL (kind=Dkind), PARAMETER :: AccelGravity = 32.1740_Dkind
  ! Avagrado's Number (atoms or molecules)/(Mole)
  REAL (kind=Dkind), PARAMETER :: Avagrado = 6.02214D23
  !
CONTAINS

  REAL(Dkind) FUNCTION  F2Kconvert(temperature)
    !------------------------------------------------------------------
    ! Function F2Kconvert temperature
    !  Converts temperature From Degree Fahrenheit to Degree Kelvin
    !
    ! Language: Fortran 90+
    !
    ! Author   : Samuel L. Forkner
    ! Company  : Signal Mountain Software
    ! Date     : 2018/04/18 @ 17:40 PM
    !------------------------------------------------------------------
    IMPLICIT NONE


    !  Function arguments
    REAL(DKind), INTENT(IN) :: temperature
    !
    F2Kconvert =(temperature - F_0C) * CperF + K_0C

  END FUNCTION F2Kconvert

  REAL(Dkind) FUNCTION  K2Fconvert(temperature)
    !------------------------------------------------------------------
    ! Function K2Fconvert temperature
    !  Converts temperature From Degree Kelvin to Degree Fahrenheit
    !
    ! Language: Fortran 90+
    !
    ! Author   : Samuel L. Forkner
    ! Company  : Signal Mountain Software
    ! Date     : 2018/04/18 @ 17:40 PM
    !------------------------------------------------------------------
    IMPLICIT NONE


    !  Function arguments
    REAL(DKind), INTENT(IN) :: temperature
    !
    K2Fconvert =(temperature - K_0C) * FperC + F_0C

  END FUNCTION K2Fconvert

END MODULE ConstantsConversions

PROGRAM  LMFR_TH
  !------------------------------------------------------------------
  ! Program LMFR_TH Description:
  !   Program for Liquid Metal cooled Fast Reactor (LMFR)
  !   Thermal-Hydraulic Analysis. Current only single channel for
  !   PRISM reactor with U/PU metal fuel. Eventually will expand to
  !   other designs and full core analysis.
  !
  ! Methods & References:
  !  <methods & references>
  !
  ! Language: Fortran 90+
  !
  ! Author   : Samuel L. Forkner
  ! Company  : Signal Mountain Software
  ! Date     : 2018/04/16 @ 9:12 AM
  !------------------------------------------------------------------
  USE DefineKinds, ONLY: Ikind, Dkind
  USE CaseInfo, ONLY: InitializeFileInfo, ReadInput
  USE Calculations, ONLY: Calc1
  !
  IMPLICIT NONE
  !
  ! initialize the file information
  CALL InitializeFileInfo
  !
  ! Read Input Data
  !
  CALL ReadInput
  !
  ! do Calculations
  !
  CALL Calc1
  !
  STOP
END PROGRAM LMFR_TH

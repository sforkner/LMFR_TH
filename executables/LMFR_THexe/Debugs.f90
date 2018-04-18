MODULE  Debugs
  !------------------------------------------------------------------
  ! Module Debugs Description:
  !   Used for debug Edits
  !
  ! Language: Fortran 90+
  !
  ! Author   : Samuel L. Forkner
  ! Company  : Signal Mountain Software
  ! Date     : 2018/04/16 @ 20:36 PM
  !------------------------------------------------------------------
  USE DefineKinds, ONLY: Ikind, Dkind, Ckind
  USE CaseInfo, ONLY:  TIN
  USE SodiumTHProperties, ONLY: Cp_Na, K_Na
  USE SS304THProperties, ONLY:  Cp_SS, K_SS
  USE UmetalTHProperties, ONLY:  Cp_U, K_U
  USE MOXTHProperties, ONLY: K_MOX, Cp_MOX
  USE THProperties, ONLY: KFU, KCL, KBC, VIS, RHO, CP, TSAT
  !
  IMPLICIT NONE
  REAL (Dkind) TKIN
  !
CONTAINS
  !
  SUBROUTINE InletDebug
    WRITE(2,"(' RHO(TIN)=',G15.5)") RHO(TIN)
    WRITE(2,"(' VIS(TIN)=',G15.5)") VIS(TIN)
    WRITE(2,"(' CP(TIN)=',G15.5)") CP(TIN)
    WRITE(2,"(' KFU(TIN)=',G15.5)") KFU(TIN)
    WRITE(2,"(' KCL(TIN)=',G15.5)") KCL(TIN)
    WRITE(2,"(' KBC(TIN)=',G15.5)") KBC(TIN)
    TKIN = (TIN-32.0D0)*(5.0D0/9.0D0)+273.15D0
    WRITE(2,"(' TKIN (Deg K)=',G15.5)") TKIN
    WRITE(2,"(' Cp_Na(TKIN)=',G15.5)") Cp_Na(TKIN)
    WRITE(2,"(' Cp_SS(TKIN)=',G15.5)") Cp_SS(TKIN)
    WRITE(2,"(' Cp_U(TKIN)=',G15.5)") Cp_U(TKIN)
    WRITE(2,"(' Cp_MOX(TKIN)=',G15.5)")   &
         &        Cp_MOX(0.0D0,2.0D0,0.0D0,0.0D0,TKIN)
    WRITE(2,"(' K_Na(TIKN)=',G15.5)") K_Na(TKIN)
    WRITE(2,"(' K_U(TKIN)=',G15.5)") K_U(TKIN)
    WRITE(2,"(' K_MOX(TKIN)=',G15.5)")   &
         &       K_MOX(0.96D0, 0.0D0,0.0D0, 2.0D0, 0.0D0,TKIN)
    WRITE(2,"(' K_SS(TKIN)=',G15.5)") K_SS(TKIN)
    !
    RETURN
  END  SUBROUTINE InletDebug
  !
END MODULE Debugs

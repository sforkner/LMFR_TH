MODULE Calculations
  !+-------------------------------------------------------------------+
  !|                                                                   |
  !|  Module Calculations.f90                                          |
  !|                                                                   |
  !+-------------------------------------------------------------------+
  !
  USE DefineKinds, ONLY: Ikind, Dkind
  USE ConstantsConversions, ONLY: GSUBC, PI, SECPERHR, SQINPERSQFT
  !
  USE CaseInfo, ONLY:  DFU, DCL, S, HP, H, GIN, TIN, PIN, QTP0, HE,     &
       &                    noutput, DWW
  USE THProperties, ONLY: KFU, KCL, KBC, VIS, RHO, CP, TSAT
  USE Debugs, ONLY: InletDebug
  !
  IMPLICIT NONE
  !
  INTEGER (Ikind), PARAMETER :: INCREMENTS = 25
  REAL (Dkind) AFL, PW, DEQ, DR, SRAT, FRAT, Y, FWWM, MDOT, DELL,      &
       & UB, RE, FCIRC, F, DELPF, DELPE, P1, TB1, ZC, QTPAV, TB2, TBAV, P2, &
       & PAV, TSATAV, PR, PE, NU, PE0, HTC, K, TSAV, KCL1, TF, TCLAV, KCL2, &
       & KFU1, T0, TFAV, KFU2, T14, T12, T34, PEX, TEX, POVL, T0max,        &
       & TCLAVMax, TBmax, DTsubcmin, DELPH
  !
CONTAINS
  SUBROUTINE Calc1
    INTEGER (Ikind)  i, j
    !
    ! Geometry Calculations
    !
    ! wire wrap diameter (DWW) assumed to be gap between pins
    ! DWW = S - DCL
    ! flow area
    AFL = COS(PI/6.0)*S*S/2.0 - (DCL*DCL + DWW*DWW)*PI/8.0
    ! wetted perimeter
    PW = (DCL + DWW)*PI/2.0
    ! equivalent diameter
    DEQ = 4.0*AFL/PW
    ! diameter ratio wire wrap to equivalent channel
    DR = DWW/DEQ
    ! pitch to rod diameter ratio
    SRAT = S/DCL
    ! ratio of friction factor in smooth triangular coolant channel
    ! to circular tube friction factor
    FRAT = -3.4571D0 + SRAT*(8.0528D0 + SRAT*(-4.7047D0 +            &
         &                  SRAT*0.91621D0))
    ! ratio of friction factor in wire-wrapped coolant channel to
    ! that in smooth triangular channel
    Y = 0.18598D0 + DR*(-1.3640D0 + DR*2.8944D0)
    FWWM = 1.195D0 - Y*(HP - 0.75D0)/0.25D0
    ! mass flow rate
    MDOT = GIN*AFL
    !  axial calculation increment
    DELL = H/INCREMENTS
    !
    !***debug
    CALL InletDebug
    WRITE(2,*) '  '
    WRITE(2,"(' AFL=',G15.5)") AFL
    WRITE(2,"(' DEQ=',G15.5)") DEQ
    WRITE(2,"(' DR=',G15.5)") DR
    WRITE(2,"(' SRAT=',G15.5)") SRAT
    WRITE(2,"(' FRAT=',G15.5)") FRAT
    WRITE(2,"(' FWWM=',G15.5)") FWWM
    WRITE(2,"(' MDOT=',G15.5)") MDOT
    WRITE(2,*) '  '
    !***debug
    !
    ! initialize max and min vlues storages
    T0max = 0.0
    TCLAVMax = 0.0
    TBmax = 0.0
    DTsubcmin = 100000.0
    !
    ! Start calculations of pressure loss at core inlet
    ! Loop over all increments
    DO i=1,INCREMENTS
       !
       IF (i.EQ.1) THEN
          ! inlet loss coefficient = 7.5 a guess
          K = 7.5
          !
          ! calculate velocity, friction factor Reynols No, delta Ps
          CALL DeltaP(TIN,K,4.54D0)
          !
          !
          ! Calculations of conditions at beginning of first increment
          !
          ! Bottom of increment pressure
          P1 = PIN -(DELPF+DELPE+DELPH)/SQINPERSQFT
          ! Bottom of increment coolant temperature
          TB1 = TIN
          !
          ! mid 1st increment axial coordinate
          ZC =  -H/2.0 + DELL/2.0
          !
          WRITE(noutput,15)   UB, RE, FCIRC, (DELPF/SQINPERSQFT),      &
               &                   (DELPE/SQINPERSQFT),(DELPH/SQINPERSQFT)
15        FORMAT(/,'  At Inlet of The Core',/,                         &
               &        '  Channel Ave. Velocity (Ft/Sec):', T60,F8.2,/,    &
               &        '  Reynolds Number (Dimensionless):',T58,F10.1,/,   &
               &        '  Friction Factor (Dimensionless):',T60,F8.4,/,    &
               &        '  tube friction Pressure Loss (psi):',T60,F8.4,/,  &
               &        '  Inlet Local Pressure Loss (psi):',T60,F8.4,/,    &
               &        '  Inlet Elevation Pressure Loss (psi):',T60,F8.4,/,&
               &        '                                                 ')
       ELSE
          ! other local losses are zero
          K = 0.0
          !  increment inlet temp & pressure are outlet of previous
          TB1 = TB2
          P1 = P2
          !
          ! axial coordinate at center of increment
          ZC = ZC + DELL
       END IF
       !
       IF(i.EQ.  1.OR.i.EQ. 51.OR.i.EQ.101.OR.i.EQ.151.OR.i.EQ.201    &
            &          .OR.i.EQ.251.OR.i.EQ.301.OR.i.EQ.351)THEN
          !
          ! Print out Column Headings
          !
          WRITE(noutput,25)
25        FORMAT(///,' Inc.    Z-Loc     TFc      TFs    T-Clad',       &
               &          '  T-Bulk    P-Avg    T-Sat'    &
               &          '     HT-Coeff    Delta P   Delta P   Delta P',/   &
               &          '  No.       ft    Deg F    Deg F    Deg F    '    &
               &          'Deg F     psia    Deg F BTU/Hr-Ft**2'             &
               &          '      psi       psi      psi  ')
       END IF
       !
       ! Calculate heat input to increment
       QTPAV = QTP0*COS(PI*ZC/HE)
       !
       ! calcuate top of increment coolant temperature
       !
       ! make two passes on CP
       TBAV = TB1
       TB2 = TB1+(QTPAV*(PI/4.0)*DFU*DFU*DELL/(MDOT*CP(TBAV)))/2.0
       !
       ! update increment average coolant temperature
       TBAV = (TB1+TB2)/2.0
       !  2nd pass use new TBAV
       TB2 = TB1+(QTPAV*(PI/4.0)*DFU*DFU*DELL/(MDOT*CP(TBAV)))/2.0
       ! update increment average coolant temperature
       TBAV = (TB1+TB2)/2.0
       !
       ! calculate increment delta p based on TBAV, local loss coefficient=0.0
       CALL DeltaP(TBAV,K,1.0D0*DELL)
       ! Top of increment pressure
       P2 = P1 -(DELPF+DELPE+DELPH)/SQINPERSQFT
       ! average pressure in increment
       PAV = (P1+P2)/2.0
       ! average saturation temperature for increment
       TSATAV = TSAT(PAV*SQINPERSQFT)
       ! Prandtl number
       PR= CP(TBAV)*VIS(TBAV)/KBC(TBAV)
       ! Peclet number
       PE = RE*PR
       ! Nusselt Number
       NU = -88.0 + SRAT*(130.8 - 41.67*SRAT)
       PE0 =336.0*SRAT-110.0
       IF(PE.GT.PE0) THEN
          NU = NU+0.025*(PE-PE0)**0.8
       END IF
       ! consider using this W correlation for NU
       ! NU = 4.0+0.33*(SRAT**3.8)*(PE/100)**0.86 + 0.16*SRAT**5
       ! for 1.1 <=SRAT<=1.5 and 10<PE<5000
       !**********
       ! inteveral heat transfer coefficient
       HTC = NU*KBC(TBAV)/DEQ
       !
       ! calculate temperature at outside surface of clad
       TSAV = TBAV+QTPAV*DFU*DFU/(4.0*DCL*HTC)
       !
       ! check for possibility of surface boiling
       IF(TSAV.GE.TSATAV) THEN
          ! possible boiling at clad surface-violates code assumption
          WRITE(noutput,35) i, TSAV,TSATAV
35        FORMAT('0 ***************Surface Boiling*****************',/,&
               &      '  Increment=',I3,' Coolant Temp=',F9.2,              &
               &      '  Sat. Temp=',F9.2,/)
          STOP
       END IF
       !
       ! calculate temperature at outside surface of fuel pellet
       ! iterate until converged
       !
       KCL1 = KCL(TSAV)
       DO
          TF = TSAV+(QTPAV*DFU*DFU/(8.0*KCL1))*LOG(DCL/DFU)
          TCLAV = (TSAV+TF)/2.0
          KCL2 = KCL(TCLAV)
          IF(ABS((KCL1-KCL2)/KCL2).LT.0.0001) EXIT
          KCL1 = KCL2
       END DO
       !
       ! calculate temperature at center of fuel pellet
       ! iterate until converged
       !
       KFU1 = KFU(TF)
       DO
          T0 = TF+QTPAV*DFU*DFU/(16.0*KFU1)
          TFAV = (T0+TF)/2.0
          KFU2 = KFU(TFAV)
          IF(ABS((KFU1-KFU2)/KFU2).LT.0.0001) EXIT
          KFU1 = KFU2
       END DO
       !
       ! Calculate fuel temperature at 1/4, 1/2 and 3/4 radius
       !
       ! T14  = T0-QTPAV*DFU*DFU/(256.0*KFU1)
       T14  = T0-QTPAV*(DFU**2.2)/(256.0*KFU1)
       T12  = T0-QTPAV*DFU*DFU/(64.0*KFU1)
       T34  = T0-9.0*QTPAV*DFU*DFU/(256.0*KFU1)
       !
       ! print out results for this increment
       !
       WRITE(noutput,45) I, ZC, T0, TF, TSAV, TBAV,    &
            &    PAV, TSATAV, HTC, (DELPF/SQINPERSQFT),     &
            &   (DELPE/SQINPERSQFT), (DELPH/SQINPERSQFT)
45     FORMAT(1H ,I4,F9.4,6F9.1,F13.1,3F10.4)
       !
       ! store extreme values
       !
       IF(T0.GE.T0max) T0max = T0
       IF(TCLAV.GT.TCLAVMax) TCLAVMax = TCLAV
       IF(TB2.GT.TBmax) TBmax = TB2
       IF((TSATAV-TB2).LT.DTsubcmin) DTsubcmin = (TSATAV-TB2)
       !
       ! print exit conditions
       IF (i.EQ.INCREMENTS) THEN
          ! exit loss coefficient = 3.0 a guess
          K = 3.0
          !
          ! calculate velocity, friction factor Reynols No, delta Ps
          CALL DeltaP(TB2,K,6.75D0)
          !
          WRITE(noutput,55)   UB, RE, FCIRC, (DELPF/SQINPERSQFT),      &
               &                   (DELPE/SQINPERSQFT), (DELPH/SQINPERSQFT)
55        FORMAT('0 At Outlet of The Core',/,                          &
               &      '  Channel Ave. Velocity (Ft/Sec):', T60,F8.2,/,      &
               &      '  Reynolds Number (Dimensionless):',T58,F10.1,/,     &
               &      '  Friction Factor (Dimensionless):',T60,F8.4,/,      &
               &      '  tube friction Pressure Loss (psi):',T60,F8.4,/,    &
               &      '  Local Pressure Loss (psi):',T60,F8.4,/,            &
               &      '  Elevation Pressure Loss (psi):',T60,F8.4,/,        &
               &      '                                                 ')
          PEX = P2 - (DELPF+DELPE+DELPH)/SQINPERSQFT
          TEX = TB2
          WRITE(noutput,65) TEX, (TEX-TIN), PEX, (PIN-PEX), T0max,     &
               &                 TCLAVmax, TBmax, DTsubcmin
65        FORMAT(' Exit Temperature(Deg F):',T60,F8.2,/,               &
               &      ' Core Temperature Rise(Deg F):',T60,F8.2,/,          &
               &      ' Exit Pressure(psia):',T60,F8.2,/,                   &
               &      ' Core Pressure Drop (psi):',T60,F8.4,/,              &
               &      ' Max Centerline Fuel Temp.(Deg F)',T60,F8.2,/,       &
               &      ' Max Avg Clad Temperature(Deg F)',T60,F8.2,/,        &
               &      ' Max Bulk Coolant Temp (Deg F)',T60,F8.2,/,          &
               &      ' Min Coolant Subcooling (Deg F)',T60,F8.2,//)
          !
       END IF
       !
    END DO
    !
    ! calculate rod average KW/FT
    POVL = (QTP0*HE*DFU*DFU*SIN(PI*H/(2.0*HE)))/(2.0*H*3412.0)
    !
    WRITE(noutput,75) POVL
75  FORMAT(/,' Average Power Per Unit Length (KW/FT):',T58,F10.3,//)
    !
    RETURN
    !
  END SUBROUTINE Calc1
  SUBROUTINE DeltaP(T,KL,L)
    !+-------------------------------------------------------------------+
    !|                                                                   |
    !|  DeltaP.F90                                                       |
    !|                                                                   |
    !|  Computes dressure drop components for increment based on         |
    !|  specified local temperature (T) and loss coefficient KL.         |
    !|                                                                   |
    !|  T is temperature to use for properties                           |
    !|  KL is the local (geometry change) loss coeficient                |
    !|  L is the increment length for tube friction                      |
    !|                                                                   |
    !+-------------------------------------------------------------------+
    IMPLICIT NONE
    !
    REAL (Dkind) T, KL, L
    !
    ! average velocity in flow channel
    UB = (Gin/SECPERHR)/RHO(T)
    ! Reynolds Number
    RE = UB*DEQ*RHO(T)/(VIS(T)/SECPERHR)
    ! friction factor in circular tubes
    FCIRC = 0.184D0/RE**0.2D0
    ! friction factor
    F = FCIRC*FRAT*FWWM
    ! tube friction pressure drop
    DELPF = F*L*RHO(T)*UB*UB/(2.0*DEQ*GSUBC)
    !  local loss coefficient (geometry change) pressure drop
    DELPE = KL*RHO(T)*UB*UB/(2.0*GSUBC)
    ! elevation head loss (psft)
    DELPH = RHO(T)*L
    !
  END SUBROUTINE DeltaP

END MODULE Calculations

MODULE CaseInfo
  ! +----------------------------------------------------------------+
  ! |                                                                |
  ! |  Module CaseInfo.f90                                           |
  ! |                                                                |
  ! +----------------------------------------------------------------+
  USE DefineKinds, ONLY: FileInfo, Ikind, Lkind, Dkind, FileInfo
  IMPLICIT NONE
  !
  INTEGER (Ikind)  ninput, noutput
  REAL (Dkind) DFU, DCL, S, HP, H, DWW, GIN, TIN, PIN, QTP0, HE,       &
       &            DFUinch, DCLinch, Sinch, HPinch, Hinch, HEinch, DWWinch
  !
  INTEGER (Ikind), PRIVATE :: i, io
  CHARACTER (Len=132) card, title
  CHARACTER (Len=1), PARAMETER :: exclam = "!", quote = "'"
  CHARACTER (Len=1) :: Yes
  REAL (Dkind), PARAMETER :: PI = 3.141592654, INCHESPERFT = 12.0,     &
       &  GSUBC = 32.2, SQINPERSQFT =144.0, SECPERHR = 3600.0
  TYPE (FileInfo), DIMENSION(2) :: Files
  !
  !  DFU       = Fuel pellet outer diameter (feet)
  !  DFUinch   = Fuel pellet outer diameter (inches)
  !  DCL       = Fuel rod clad outer diameter (feet)
  !  DCLinch   = Fuel rod clad outer diameter (inches)
  !  S         = Fuel rod pitch (feet)
  !  Sinch     = Fuel rod pitch (inches)
  !  HP        = Wire wrap pitch (feet)
  !  HPinch    = Wire wrap pitch (inches)
  !  H         = Active core height (feet)
  !  Hinch     = Active core height (inches)
  !  DWW       = Wire Wrap outer diameter (feet)
  !  DWWinch   = Wire Wrap outer diameter (inches)
  !  GIN       = Inlet mass velocity Lb/(Hr-Ft**2)
  !  TIN       = Core inlet temperature (Deg F)
  !  PIN       = Core inlet pressure (psia)
  !  QTP0      = Volumetric thermal source strength BTU/(Hr-Ft**3)
  !  HE        = Extrapolated core height (feet)
  !  HEinch    = Extrapolated core height (inches)
  !
CONTAINS
  !
  ! initilize files info and read in names of input & Output files
  !
  SUBROUTINE InitializeFileInfo
    !
    Files%name = "                                                      &
         &                                                              "
    Files%unitno = (/(i,i=1,2)/)
    Files%opened = .FALSE.
    !
    DO
       WRITE(*,'(A29)',Advance='NO') '  Input name of Print file: '
       READ(*,*) files(2)%name
       WRITE(*,*) TRIM(files(2)%name)
       WRITE(*,*)
       OPEN (UNIT=files(2)%unitno,FILE=files(2)%name,STATUS="NEW",IOSTAT=io)
       IF(io /= 0) THEN
          WRITE(*,"(A44)",Advance="NO") " File exists! Do you want to reuse (Y/N)? "
          READ(*,'(A1)') Yes
          WRITE(*,*)
          IF (Yes == 'y' .OR. Yes == 'Y') THEN
             OPEN (UNIT=files(2)%unitno,FILE=files(2)%name,STATUS='REPLACE')
             EXIT
          ELSE
             CYCLE
          ENDIF
       ELSE
          EXIT
       ENDIF
    END DO
    files(2)%opened = .TRUE.
    noutput = files(2)%unitno
    !
    DO
       WRITE(*,'(A31)',Advance='NO') '  Input name of Input Deck: '
       READ(*,'(A)') files(1)%name
       WRITE(*,*) TRIM(files(1)%name)
       WRITE(*,*)
       OPEN (UNIT=files(1)%unitno,FILE=files(1)%name,STATUS='OLD',IOSTAT=io)
       IF(io /= 0) THEN
          WRITE(*,'(A35)') ' File Does not exist! Try again.'
          CYCLE
        ELSE
          EXIT
       ENDIF
    END DO
    files(1)%opened = .TRUE.
    ninput = files(1)%unitno
    !
  END SUBROUTINE InitializeFileInfo
  !
  ! Read input data and list on output
  !
  SUBROUTINE ReadInput
    !
    ! put code information in print file
    !
    WRITE(noutput,'(//,7X,A)') 'PRISM LMFR Thermal/Hydraulics Code v0.1'
    WRITE(noutput,'(/,A)') ' '
    !
    ! Read in Title for case
    !
    READ(ninput,'(A132)',END=999) title
    WRITE(noutput,'(/,A132,/)') title
    WRITE(noutput,*) '  Output to written to file: ', files(2)%name
    WRITE(noutput,*)  '  Input deck from file: ', files(1)%name
    !
    ! read input data
    !
    WRITE(noutput,'(//,3X,A50)') ' Begin Reading Case Input Card 1'
    card(1:1) = quote
    DO WHILE (card(1:1).EQ.quote.OR.card(1:1).EQ.exclam)
       READ(ninput,'(A132)',END=999) card
       WRITE(noutput,'(1H ,A132)') card
    ENDDO
    BACKSPACE ninput
    READ(ninput,*,END=999)  DFUinch, DCLinch, Sinch, HPinch, Hinch,  &
         &                       DWWinch
    WRITE(noutput,'(1H ,3X,A50)') ' End Reading Case Input Card 1'
    !
    WRITE(noutput,'(//,3X,A50)') ' Begin Reading Case Input Card 2'
    card(1:1) = quote
    DO WHILE (card(1:1).EQ.quote.OR.card(1:1).EQ.exclam)
       READ(ninput,'(A132)',END=999) card
       WRITE(noutput,'(1H ,A132)') card
    ENDDO
    BACKSPACE ninput
    READ(ninput,*,END=999)  GIN, TIN, PIN
    WRITE(noutput,'(1H ,3X,A50)') ' End Reading Case Input Card 2'
    !
    WRITE(noutput,'(//,3X,A50)') ' Begin Reading Case Input Card 3'
    card(1:1) = quote
    DO WHILE (card(1:1).EQ.quote.OR.card(1:1).EQ.exclam)
       READ(ninput,'(A132)',END=999) card
       WRITE(noutput,'(1H ,A132)') card
    ENDDO
    BACKSPACE ninput
    READ(ninput,*,END=999)  QTP0, HEinch
    WRITE(noutput,'(1H ,3X,A50)') ' End Reading Case Input Card 3'
    !
    ! all input read -- close input file unit
    CLOSE(UNIT=ninput)
    files(1)%opened = .FALSE.
    !
    ! Write input data to output file
    !
    WRITE(noutput,15) DFUinch, DCLinch, Sinch, HPinch, DWWinch,      &
         &                 Hinch, GIN, TIN, PIN, QTP0, HEinch
15  FORMAT(//,' Fuel pellet outer diameter (inches):',T60,F8.4,/,    &
         &         ' Fuel rod clad outer diameter (inches):',T60,F8.4,/,  &
         &         ' Fuel rod pitch (inches):',T60,F8.4,/,                &
         &         ' Wire wrap pitch (inches):',T60,F8.4,/,               &
         &         ' Wire wrap diameter (inches):',T60,F8.4,/,            &
         &         ' Active core height (inches):',T60,F8.4,/             &
         &         ' Inlet mass velocity Lb/(Hr-Ft**2):',T55,1PG13.4,/    &
         &         ' Core inlet temperature (Deg F):',T60,0PF8.2,/        &
         &         ' Core inlet pressure (psia):',T60,F8.2,/              &
         &         ' Volumetric thermal source strength BTU/(Hr-Ft**3):', &
         &           T55,1PG13.4,/                                        &
         &         ' Extrapolated core height (Inches):',T60,0PF8.2,/,    &
         &         '   End Of Input',//)
    !
    ! convert DFU, S, DCL, H, HE & HP form units of inches to feet
    !
    DFU = DFUinch / INCHESPERFT
    S   = Sinch / INCHESPERFT
    DCL = DCLinch / INCHESPERFT
    H   = Hinch / INCHESPERFT
    HE  = HEinch / INCHESPERFT
    HP  = HPinch / INCHESPERFT
    DWW = DWWinch / INCHESPERFT
    !
    RETURN
    !
999 WRITE(noutput,*) ' *==> End of file on Unit 1  <==*'
    STOP
    !
  END SUBROUTINE ReadInput
  !
END MODULE CaseInfo

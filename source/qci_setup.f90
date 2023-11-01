MODULE QCISETUP

   CONTAINS
      SUBROUTINE QCI_INIT(PARAMETERFILE, ALIGNT)
         USE QCIKEYS, ONLY: NATOMS, QCIAMBERT
         USE QCICONSTRAINTS, ONLY: CREATE_CONSTRAINTS
         USE QCIPERMDIST, ONLY: ALLOC_QCIPERM, INIT_PERMALLOW(NATOMS)
         USE MOD_FREEZE, ONLY: GET_FROZEN_ATOMS
         USE CHIRALITY, ONLY: FIND_CHIRAL_CENTRES
         USE REPULSION, ONLY: NREPCURR, ALLOC_REP_VARS
         IMPLICIT NONE
         CHARACTER(30), INTENT(IN) :: PARAMETERFILE
         LOGICAL, INTENT(IN) :: ALIGNT

         ! parse settings
         CALL PARSE_SETTINGS(PARAMETERFILE)
         ! starting permutational setup perm.allow
         CALL ALLOC_QCIPERM()
         CALL INIT_PERMALLOW(NATOMS)
         ! align endpoints
         IF (ALIGNT) CALL ALIGN_ENDPOINTS()
         ! get frozen atoms setup
         IF (QCIFREEZET) CALL GET_FROZEN_ATOMS()
         ! get constraints
         CALL CREATE_CONSTRAINTS()
         ! setting up repulsions
         ! we use NATOMS as initial size here for the number of repulsions
         NREPCURR = NATOMS
         CALL ALLOC_REP_VARS(NREPCURR)
         ! get chiral information for AMBER
         IF (QCIAMBERT) CALL FIND_CHIRAL_CENTRES()
      END SUBROUTINE QCI_INIT

      SUBROUTINE PARSE_SETTINGS(PARAMETERFILE)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         IMPLICIT NONE
         CHARACTER(30), INTENT(IN) :: PARAMETERFILE
         INTEGER :: PARAMUNIT, IOS
         LOGICAL :: YESNO, EOFT
         CHARACTER(25) :: ENTRY, VALUE 
            
         ! open parameter file
         INQUIRE(FILE=PARAMETERFILE, EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(*,*) " parse_settings> Cannot locate QCI settings file ", PARAMETERFILE
            STOP
         END IF
         PARAMUNIT = GETUNIT()
         OPEN(PARAMUNIT, FILE=PARAMETERFILE, STATUS="OLD")
            
         !loop over lines in file
         EOFT = .FALSE.
         DO WHILE (.NOT. EOFT)
            READ(PARAMUNIT,*,IOSTAT=IOS) ENTRY, VALUE
            IF (IOS.GT.0) THEN
               EOFT = .TRUE.     
            ELSE
               CALL SETKEYS(ENTRY, VALUE)
            ENDIF
         END DO
         CLOSE(PARAMUNIT)
      END SUBROUTINE PARSE_SETTINGS

      SUBROUTINE ALIGN_ENDPOINTS()
         USE QCIPERMDIST, ONLY: LOPERMDIST

         CALL LOPERMDIST(XFINAL,XSTART,NATOMS,.FALSE.,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,0,NMOVE,NEWPERM)
      END SUBROUTINE ALIGN_ENDPOINTS

      SUBROUTINE SETKEYS(ENTRY, VAL)
         USE QCIKEYS
         USE AMBER_CONSTRAINTS, ONLY: AMBERCONSTRFILE, TOPFILENAME
         IMPLICIT NONE
         CHARACTER(25), INTENT(IN) :: ENTRY, VAL

         IF (ENTRY.EQ."COMMENT") THEN
            RETURN
         ELSE IF (ENTRY.EQ."QCIDOBACK") THEN
            QCIDOBACK = .TRUE.           
         ELSE IF (ENTRY.EQ."QCIRESTART") THEN
            QCIRESTART = .TRUE.
         ELSE IF (ENTRY.EQ."QCIREADGUESS") THEN
            QCIREADGUESS = .TRUE.            
         ELSE IF (ENTRY.EQ."FREEZEFILE") THEN
            QCIFREEZET = .TRUE.            
            FREEZEFILE = VAL
         ELSE IF (ENTRY.EQ."NMINUNFROZEN") THEN
            READ(VAL, *) NMINUNFROZEN
         ELSE IF (ENTRY.EQ."FREEZEFILE") THEN
            FREEZEFILE = VAL
         ELSE IF (ENTRY.EQ."QCIFREEZE") THEN
            QCIFREEZET = .TRUE.
            READ(VAL, *) QCIFREEZETOL
         ELSE IF (ENTRY.EQ."QCICONSEP") THEN
            READ(VAL, *) QCICONSEP
         ELSE IF (ENTRY.EQ."QCICONSTRAINTTOL") THEN
            READ(VAL, *) QCICONSTRAINTTOL
         ELSE IF (ENTRY.EQ."QCICONCUT") THEN
            READ(VAL, *) QCICONCUT
         ELSE IF (ENTRY.EQ."GEOMFILE") THEN
            GEOMFILE = VAL
         ELSE IF (ENTRY.EQ."CONSTRFILE") THEN
            CONSTRFILE = VAL
         ELSE IF (ENTRY.EQ."QCIPERMCHECKINT") THEN
            READ(VAL, *) QCIPERMCHECKINT
         ELSE IF (ENTRY.EQ."QCIPERMCHECK") THEN
            IF (VAL.EQ."T") QCIPERMCHECK = .TRUE.
         ELSE IF (ENTRY.EQ."QCILPERMDIST") THEN
            IF (VAL.EQ."T") QCILPERMDIST = .TRUE.
         ELSE IF (ENTRY.EQ."QCIPERMT") THEN
            IF (VAL.EQ."T") QCIPERMT = .TRUE.
         ELSE IF (ENTRY.EQ."QCIPERMCUT") THEN
            READ(VAL, *) QCIPERMCUT
         ELSE IF (ENTRY.EQ."AMBERCONSTRFILE") THEN
            AMBERCONSTRFILE = VAL
         ELSE IF (ENTRY.EQ."TOPFILENAME") THEN
            TOPFILENAME = VAL
         ELSE IF (ENTRY.EQ."HIRECONSTRFILE") THEN
            HIRECONSTRFILE = VAL
         ELSE IF (ENTRY.EQ."HIRETOPFILE") THEN
            HIRETOPFILE = VAL
         ELSE IF (ENTRY.EQ."SBMCONTACTFILE") THEN
            SBMCONTACTFILE = VAL
         ELSE IF (ENTRY.EQ."QCIMODE") THEN
            IF (VAL.EQ."AMBER") THEN
               QCIAMBERT = .TRUE.
            ELSE IF (VAL.EQ."HIRE") THEN
               QCIHIRET = .TRUE.
            ELSE IF (VAL.EQ."SBM") THEN
               QCISBT = .TRUE.
            ELSE 
               WRITE(*,*) " setkeys> QCI mode ", VAL, " is not a valid function"
            END IF
         ELSE
            WRITE(*,*) " setkeys> Cannot find setting ", ENTRY, " - will skip this entry"
         END IF
      END SUBROUTINE SETKEYS

END MODULE QCISETUP
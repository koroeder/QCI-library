MODULE MOD_FREEZE

   CONTAINS

      SUBROUTINE GET_FROZEN_ATOMS()
         USE QCIKEYS, ONLY: NATOMS, QCIFROZEN, NQCIFROZEN, FREEZE, QCIFREEZETOL, NMINUNFROZEN
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE QCIPREC
         IMPLICIT NONE
         REAL(KIND=REAL64) :: DSORTED(NATOMS)
         INTEGER :: IDSORTED(NATOMS)
         REAL(KIND=REAL64) :: FREEZETOL2, DIST
         INTEGER :: J1, J2, J3
         
         CALL ALLOC_FREEZE(NATOMS)
         CALL READ_FROZEN_ATOMS()
         QCIFROZEN(1:NATOMS) = .FALSE.
         NQCIFROZEN = 0

         IDSORTED(1:NATOMS)=-1
         DSORTED(1:NATOMS)=1.0D100

         FREEZETOL2 = QCIFREEZETOL**2

         DO J1=1,NATOMS
            DIST = (XSTART(3*(J1-1)+1)-XFINAL(3*(J1-1)+1))**2 &
                 + (XSTART(3*(J1-1)+2)-XFINAL(3*(J1-1)+2))**2 &
                 + (XSTART(3*(J1-1)+3)-XFINAL(3*(J1-1)+3))**2
            IF ((DIST.LT.FREEZETOL2).OR.FREEZE(J1)) THEN
               NQCIFROZEN  = NQCIFROZEN + 1
               QCIFROZEN(J1) = .TRUE.
               ! make sure we keep all forced ones later
               IF (FREEZE(J1)) THEN
                  DIST = -DIST
               END IF
            END IF

            DO J2=1,J1
               IF (DIST.LT.DSORTED(J2)) THEN
                  DO J3=J1,J2+1,-1
                     DSORTED(J3) = DSORTED(J3-1)
                     IDSORTED(J3) = IDSORTED(J3-1)
                  END DO
                  DSORTED(J2) = DIST
                  IDSORTED(J2) = J1
                  EXIT
               END IF
            END DO
         END DO
         !check the minimum number of atoms is unfrozen 
         IF ((NATOMS-NQCIFROZEN.LT.NMINUNFROZEN)) THEN
            DO J1=NATOMS, NATOMS-NQCIFROZEN+1,-1
               QCIFROZEN(IDSORTED(J1)) = .FALSE.
            END DO
            NQCIFROZEN = MAX(0,NATOMS-NMINUNFROZEN)
         END IF
      END SUBROUTINE GET_FROZEN_ATOMS

      SUBROUTINE READ_FROZEN_ATOMS()
         USE QCIKEYS, ONLY: FREEZEFILE, FREEZE, NATOMS
         USE QCIFILEHANDLER, ONLY: FILE_LENGTH, GETUNIT
         IMPLICIT NONE
         LOGICAL :: YESNO, EOFT
         INTEGER :: FREEZEUNIT, IDX, IOS

         FREEZE(1:NATOMS) = .FALSE.
         INQUIRE(FILE=FREEZEFILE, EXIST=YESNO)
         IF (.NOT.YESNO) RETURN !if we dont have this file, no need to read anything
         FREEZEUNIT = GETUNIT()
         OPEN(FREEZEUNIT,FILE=FREEZEFILE,STATUS='OLD')
         EOFT = .FALSE.
         DO WHILE (.NOT.EOFT)
            READ(FREEZEUNIT,*,IOSTAT=IOS) IDX
            IF (IOS.GT.0) THEN
               EOFT = .TRUE.
               EXIT
            ELSE
               FREEZE(IDX) = .TRUE.
            END IF
         END DO
         CLOSE(FREEZEUNIT)
      END SUBROUTINE READ_FROZEN_ATOMS

      SUBROUTINE ADD_CONSTR_AND_REP_FROZEN_ATOMS(NBEST)
         USE QCIPREC, ONLY: REAL64
         USE QCIKEYS, ONLY: NATOMS, QCIREPCUT, QCIFROZEN, QCIINTREPMINSEP
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, ATOMACTIVE, NCONSTRAINTON
         USE REPULSION, ONLY: NREPCURR, NREPULSIVE, REPI, REPJ, REPCUT, &
                              DOUBLE_ALLOC_REP
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NBEST
         INTEGER :: I, J, NI, NJ, CURRCI, CURRCJ
         LOGICAL :: SKIPREP
         REAL(KIND = REAL64) :: DSTART, DFINAL, DMIN

         ! start with finding active constraints
         NI = CONI(NBEST)
         NJ = CONJ(NBEST)
         DO I=1,NCONSTRAINT
            ! we are only interested in active constraints
            IF (.NOT.CONACTIVE(I)) CYCLE
            IF (((CONI(I).EQ.NI).OR.(CONI(I).EQ.NJ)).AND.(ATOMACTIVE(CONJ(I))).OR. &
                ((CONJ(I).EQ.NI).OR.(CONJ(I).EQ.NJ)).AND.(ATOMACTIVE(CONI(I)))) THEN
               CONACTIVE(I) = .TRUE.
               NCONSTRAINTON = NCONSTRAINTON + 1
               ! WRITE(*,*) "add_c+r_frozen> Turning on constraint ", I, " for atoms ", &
               ! CONI(I), " and ", CONJ(I)
            END IF
         END DO

         ! now find potential repulsion terms
         ! we need to run two cycles - one for NI and one for NJ
         DO I=1,NATOMS
            ! we only want active atoms
            IF (.NOT.ATOMACTIVE(I)) CYCLE
            ! we have a minimum separation for atoms in sequence
            IF (ABS(I-NI).LE.QCIINTREPMINSEP) CYCLE
            ! ignore pairs of frozen atoms
            IF (QCIFROZEN(I).AND.QCIFROZEN(NI)) CYCLE
            ! make sure we are not adding a repulsion to an existing cosntraint
            SKIPREP = .FALSE.
            DO J=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J)) CYCLE
               CURRCI = CONI(J)
               CURRCJ = CONJ(J)
               IF (((CURRCI.EQ.I).AND.(CURRCJ.EQ.NI)).OR. & 
                   ((CURRCJ.EQ.I).AND.(CURRCI.EQ.NI))) THEN
                  SKIPREP = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT.SKIPREP) THEN
               CALL DISTANCE_TWOATOMS(NATOMS,XSTART,NI,I,DSTART)
               CALL DISTANCE_TWOATOMS(NATOMS,XFINAL,NI,I,DFINAL)
               DMIN = MIN(DSTART,DFINAL)
               DMIN = MIN(DMIN-1.0D-3,QCIREPCUT)
               NREPULSIVE = NREPULSIVE + 1
               IF (NREPULSIVE.GT.NREPCURR) CALL DOUBLE_ALLOC_REP()
               REPI(NREPULSIVE) = I
               REPJ(NREPULSIVE) = NI
               REPCUT(NREPULSIVE) = DMIN
            END IF
         END DO

         ! repeat for NJ
         DO I=1,NATOMS
            ! we only want active atoms
            IF (.NOT.ATOMACTIVE(I)) CYCLE
            ! we have a minimum separation for atoms in sequence
            IF (ABS(I-NJ).LE.QCIINTREPMINSEP) CYCLE
            ! ignore pairs of frozen atoms
            IF (QCIFROZEN(I).AND.QCIFROZEN(NJ)) CYCLE
            ! make sure we are not adding a repulsion to an existing cosntraint
            SKIPREP = .FALSE.
            DO J=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J)) CYCLE
               CURRCI = CONI(J)
               CURRCJ = CONJ(J)
               IF (((CURRCI.EQ.I).AND.(CURRCJ.EQ.NJ)).OR. & 
                     ((CURRCJ.EQ.I).AND.(CURRCI.EQ.NJ))) THEN
                  SKIPREP = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT.SKIPREP) THEN
               CALL DISTANCE_TWOATOMS(NATOMS,XSTART,NJ,I,DSTART)
               CALL DISTANCE_TWOATOMS(NATOMS,XFINAL,NJ,I,DSTART)
               DMIN = MIN(DSTART,DFINAL)
               DMIN = MIN(DMIN-1.0D-3,QCIREPCUT)
               NREPULSIVE = NREPULSIVE + 1
               IF (NREPULSIVE.GT.NREPCURR) CALL DOUBLE_ALLOC_REP()
               REPI(NREPULSIVE) = I
               REPJ(NREPULSIVE) = NJ
               REPCUT(NREPULSIVE) = DMIN
            END IF
         END DO
         
      END SUBROUTINE ADD_CONSTR_AND_REP_FROZEN_ATOMS

      SUBROUTINE ALLOC_FREEZE(NATOMS)
         USE QCIKEYS, ONLY: QCIFROZEN, FREEZE
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         CALL DEALLOC_FREEZE()
         ALLOCATE(QCIFROZEN(NATOMS))
         ALLOCATE(FREEZE(NATOMS))
      END SUBROUTINE ALLOC_FREEZE

      SUBROUTINE DEALLOC_FREEZE()
         USE QCIKEYS, ONLY: QCIFROZEN, FREEZE
         IF (ALLOCATED(QCIFROZEN)) DEALLOCATE(QCIFROZEN)
         IF (ALLOCATED(FREEZE)) DEALLOCATE(FREEZE)
      END SUBROUTINE DEALLOC_FREEZE


END MODULE MOD_FREEZE
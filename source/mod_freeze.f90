MODULE MOD_FREEZE

   CONTAINS

      SUBROUTINE GET_FROZEN_ATOMS()
         USE QCIKEYS, ONLY: NATOMS, QCIFROZEN, NCQIFROZEN, FREEZE, QCIFREEZETOL, XSTART, XFINAL, NMINUNFROZEN
         USE QCIPREC
         IMPLICIT NONE
         REAL(KIND=REAL64) :: DSORTED(NATOMS)
         INTEGER :: IDSORTED(NATOMS)
         REAL(KIND=REAL64) :: FREEZETOL2, DIST
         
         CALL ALLOC_FREEZE(NATOMS)
         CALL READ_FROZEN_ATOMS()
         QCIFROZEN(1:NATOMS) = .FALSE.
         NQCIFREEZE = 0

         IDSORTED(1:NATOMS)=-1
         DSORTED(1:NATOMS)=1.0D100

         FREEZETOL2 = QCIFREEZETOL**2

         DO J1=1,NATOMS
            DIST = (XSTART(3*(J1-1)+1)-XFINAL(3*(J1-1)+1))**2 &
                 + (XSTART(3*(J1-1)+2)-XFINAL(3*(J1-1)+2))**2 &
                 + (XSTART(3*(J1-1)+3)-XFINAL(3*(J1-1)+3))**2
            IF ((DIST.LT.FREEZETOL2).OR.FREEZE(J1)) THEN
               NCQIFROZEN  = NCQIFROZEN + 1
               QCIFROZEN(J1) = .TRUE.
               ! make sure we keep all forced ones later
               IF (FREEZE(J1)) THEN
                  DIST = -DIST
               END IF
            END IF

            DO J2=1,J1
               IF (DIST.LT.DSORTED(J2)) THEN
                  DO J3=J1,J2+1,-1
                     DMOVED(J3) = DMOVED(J3-1)
                     IDSORTED(J3) = IDSORTED(J3-1)
                  END DO
                  DSORTED(J2) = DIST
                  IDSORTED(J2) = J1
                  EXIT
               END IF
            END DO
         END DO
         !check the minimum number of atoms is unfrozen 
         IF ((NATOMS-NCQIFROZEN.LT.NMINUNFROZEN)) THEN
            DO J1=NATOMS, NATOMS-NQCIFROZEN+1,-1
               QCIFROZEN(IDSORTED(J1)) = .FALSE.
            END DO
            NQCIFROZEN = MAX(0,NATOMS-NMINUNFROZEN)
         END IF
      END SUBROUTINE GET_FROZEN_ATOMS

      SUBROUTINE READ_FROZEN_ATOMS()
         USE QCIKEYS, ONLY: FREEZEFILE, FREEZE
         USE QCIFILEHANDLER, ONLY: FILE_LENGTH, GETUNIT
         IMPLICIT NONE
         LOGICAL :: YESNO, EOFT
         INTEGER :: FREEZEUNIT, IDX, IOS

         FREEZE(1:NATOMS) = .FALSE.
         INQUIRE(FILE=FREEZEFILE, EXIST=YESNO)
         IF (.NOT.YESNO) RETURN !if we dont have this file, no need to read anything
         FREEZEUNIT = GETUNIT()
         OPEN(FREEZEUNIT,FILE=FREEZEFILE,STATUS='OLD'))
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
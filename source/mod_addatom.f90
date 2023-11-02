MODULE ADDATOM
   USE QCIPREC
   IMPLICIT NONE
   LOGICAL :: BBDONE = .FALSE.

   CONTAINS

      SUBROUTINE ADDNEXTATOM()
         USE QCIKEYS, ONLY: QCIDOBACK
         IMPLICIT NONE
         
         IF (QCIDOBACK.AND.(.NOT.BBDONE)) CALL CHECK_BBLIST()

      END SUBROUTINE ADDNEXTATOM

      SUBROUTINE CHECK_BBLIST()
         USE QCIKEYS, ONLY: NATOMS
         DO J1=1,NATOMS
            IF (AABACK(J1)) THEN
               IF (.NOT.ATOMACTIVE(J1)) THEN
                  ! if there are inactive BB atoms, return
                  RETURN
               END IF
            ENDIF
         ENDDO
         BBDONE = .TRUE.
      END SUBROUTINE CHECK_BBLIST


END MODULE ADDATOM
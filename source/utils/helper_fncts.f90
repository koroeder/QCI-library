MODULE HELPER_FNCTS
   USE QCIPREC
   IMPLICIT NONE

   CONTAINS
      SUBROUTINE DISTANCE_TWOATOMS(NATOMS, X, IDX1, IDX2, DIST)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)
         INTEGER, INTENT(IN) :: IDX1, IDX2
         REAL(KIND = REAL64), INTENT(OUT) :: DIST
         REAL(KIND = REAL64) :: XYZ1(3), XYZ2(3)
         INTEGER :: J
         DO J = 1,3
            XYZ1(J) = X(3*(IDX1-1)+J)
            XYZ2(J) = X(3*(IDX2-1)+J)
         END DO
         DIST = SQRT((XYZ1(1)-XYZ2(1))**2 + (XYZ1(2)-XYZ2(2))**2 + (XYZ1(3)-XYZ2(3))**2)
      END SUBROUTINE DISTANCE_TWOATOMS

      SUBROUTINE DISTANCE_ATOM_DIFF_IMAGES(NATOMS, X1, X2, IDX, DIST)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X1(3*NATOMS), X2(3*NATOMS)
         INTEGER, INTENT(IN) :: IDX
         REAL(KIND = REAL64), INTENT(OUT) :: DIST
         INTEGER :: I

         DIST = 0.0D0
         DO I = 1,3
            DIST = DIST + (X1(3*(IDX-1)+I) - X2(3*(IDX-1)+I))**2          
         END DO
         DIST = SQRT(DIST)
      END SUBROUTINE DISTANCE_ATOM_DIFF_IMAGES

      ! double size of REAL64 array preserving old data and initialising the new parts of the array
      SUBROUTINE DOUBLE_REAL64_ARRAY(NSIZE,ARRAY,INITVAL)
         IMPLICIT NONE
         INTEGER :: NSIZE
         REAL(KIND = REAL64), ALLOCATABLE :: ARRAY(:)
         REAL(KIND = REAL64) :: INITVAL
         REAL(KIND = REAL64) :: TEMP_DATA(NSIZE)

         ! if the array is not allocated, we just allocate it here
         IF (.NOT.ALLOCATED(ARRAY)) THEN
            ALLOCATE(ARRAY(2*NSIZE))
            ARRAY(1:2*NSIZE) = INITVAL
         ! otherwise first save the data, and then reallocate
         ELSE
            TEMP_DATA(1:NSIZE) = ARRAY(1:NSIZE)
            DEALLOCATE(ARRAY)
            ALLOCATE(ARRAY(2*NSIZE))
            ARRAY(1:NSIZE) = TEMP_DATA(1:NSIZE)
            ARRAY(NSIZE+1:2*NSIZE) = INITVAL
         END IF
      END SUBROUTINE DOUBLE_REAL64_ARRAY

      ! double size of INT array preserving old data and initialising the new parts of the array
      SUBROUTINE DOUBLE_INT_ARRAY(NSIZE,ARRAY,INITVAL)
         IMPLICIT NONE
         INTEGER :: NSIZE
         INTEGER, ALLOCATABLE :: ARRAY(:)
         INTEGER :: INITVAL
         INTEGER :: TEMP_DATA(NSIZE)

         ! if the array is not allocated, we just allocate it here
         IF (.NOT.ALLOCATED(ARRAY)) THEN
            ALLOCATE(ARRAY(2*NSIZE))
            ARRAY(1:2*NSIZE) = INITVAL
         ! otherwise first save the data, and then reallocate
         ELSE
            TEMP_DATA(1:NSIZE) = ARRAY(1:NSIZE)
            DEALLOCATE(ARRAY)
            ALLOCATE(ARRAY(2*NSIZE))
            ARRAY(1:NSIZE) = TEMP_DATA(1:NSIZE)
            ARRAY(NSIZE+1:2*NSIZE) = INITVAL
         END IF
      END SUBROUTINE DOUBLE_INT_ARRAY      

      SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
         CHARACTER(*), INTENT(IN) :: LINE
         INTEGER, INTENT(IN) :: NWORDS
         CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
         INTEGER:: J1,START_IND,END_IND,J2
         CHARACTER(25) :: WORD
         START_IND=0
         END_IND=0
         J1=1
         J2=0
         DO WHILE(J1.LE.LEN(LINE))
            IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
               START_IND=J1
            ENDIF
            IF (START_IND.GT.0) THEN
               IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
               IF (J1.EQ.LEN(LINE)) END_IND=J1
               IF (END_IND.GT.0) THEN
                  J2=J2+1
                  WORD=LINE(START_IND:END_IND)
                  WORDSOUT(J2)=TRIM(WORD)
                  START_IND=0
                  END_IND=0
               ENDIF
            ENDIF
            J1=J1+1
         ENDDO
      END SUBROUTINE

END MODULE HELPER_FNCTS
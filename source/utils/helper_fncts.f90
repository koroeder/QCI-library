MODULE HELPER_FNCTS
   USE QCIPREC
   IMPLICIT NONE

   CONTAINS

      !> Euclidean norm of a vector
      FUNCTION EUC_NORM(V)
         REAL(KIND = REAL64) :: EUC_NORM
         REAL(KIND = REAL64), INTENT(IN) :: V(3)
         
         EUC_NORM = DSQRT(DOTP(3,V,V))
      END FUNCTION EUC_NORM

      !> Function to calculate dot product for any dimension of vector
      FUNCTION DOTP(NSIZE,V1,V2)
         REAL(KIND = REAL64) :: DOTP     ! dot product
         INTEGER, INTENT(IN) :: NSIZE    ! input size for vectors
         REAL(KIND = REAL64), INTENT(IN) :: V1(NSIZE), V2(NSIZE) !the vectors to be used
         INTEGER :: I !iteration variable
         
         DOTP = 0.0D0
         !use element-wise addition by iterating from 1 to NSIZE
         DO I=1,NSIZE
            DOTP = DOTP + V1(I)*V2(I)
         END DO
      END FUNCTION DOTP

      !> subroutine to get norm and normed vector
      SUBROUTINE NORM_VEC(V,VN,NORM)
         REAL(KIND = REAL64), INTENT(IN) :: V(3)
         REAL(KIND = REAL64), INTENT(OUT) :: VN(3)
         REAL(KIND = REAL64), INTENT(OUT) :: NORM

         NORM = EUC_NORM(V)
         VN(1:3) = V(1:3)/NORM
      END SUBROUTINE NORM_VEC

      FUNCTION CROSS_PROD(V1,V2) RESULT(A)
         REAL(KIND = REAL64) :: V1(3), V2(3)
         REAL(KIND = REAL64) :: A(3)

         A(1) = V1(2) * V2(3) - V1(3) * V2(2)
         A(2) = V1(3) * V2(1) - V1(1) * V2(3)
         A(3) = V1(1) * V2(2) - V1(2) * V2(1)

      END FUNCTION CROSS_PROD

      !Gram-Schmidt projection for vector V1 using V2
      FUNCTION GS_PROJECTION(V1,V2) RESULT(U)
         REAL(KIND = REAL64) :: V1(3), V2(3), U(3)

         U(1:3) = DOTP(3,V1,V2)/DOTP(3,V2,V2)*V2(1:3)
      END FUNCTION GS_PROJECTION

      REAL(KIND = REAL64) FUNCTION ANGLE(COORDS) RESULT(ANG)
         REAL(KIND = REAL64) :: COORDS(9)
         REAL(KIND = REAL64) :: VEC1(3), VEC2(3)
         REAL(KIND = REAL64) :: DP, CPMAG

         VEC1(1:3) = COORDS(1:3) - COORDS(4:6)
         VEC2(1:3) = COORDS(7:9) - COORDS(4:6)

         DP = DOTP(3,VEC1,VEC2)

         CPMAG = EUC_NORM(CROSS_PROD(VEC1,VEC2))

         ANG = ATAN2(CPMAG,DP)
      END FUNCTION ANGLE

      REAL(KIND = REAL64) FUNCTION DIHEDRAL(COORDS) RESULT(DIH)
         REAL(KIND = REAL64) :: COORDS(12)
         REAL(KIND = REAL64) :: VECS(9)
         REAL(KIND = REAL64) :: B1xB2(3), B2xB3(3), X(3), B2NORM(3)
         INTEGER :: I

         DO I=1,3
            VECS(I)   = COORDS(I+3) - COORDS(I)
            VECS(I+3) = COORDS(I+6) - COORDS(I+3)
            VECS(I+6) = COORDS(I+9) - COORDS(I+6)
         END DO
         B1xB2 = CROSS_PROD(VECS(1:3), VECS(4:6))
         B2xB3 = CROSS_PROD(VECS(4:6), VECS(7:9))
         X = CROSS_PROD(B1xB2, B2xB3)
         B2NORM = VECS(4:6)/EUC_NORM(VECS(4:6))

         ! calculate angle according to formula (Blondel and Karplus, 1996):
         ! phi = atan2( ([b1 x b2] x [b2 x b3]) . (b2/|b2|), [b1 x b2] . [b2 x b3] )
         DIH = ATAN2(DOTP(3,X,B2NORM), DOTP(3,B1xB2,B2xB3))
      END FUNCTION DIHEDRAL
    
      SUBROUTINE DISTANCE_SIMPLE(X1, X2, DIST)
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: X1(3), X2(3)
         REAL(KIND = REAL64), INTENT(OUT) :: DIST
         INTEGER :: I

         DIST = 0.0D0
         DO I=1,3
            DIST = DIST + (X1(I) - X2(I))**2
         END DO
         DIST = SQRT(DIST)
      END SUBROUTINE DISTANCE_SIMPLE

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

         ! WRITE(*,*) " read_line> ", LINE
         WORDSOUT(1:NWORDS) = ""
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
      END SUBROUTINE READ_LINE

END MODULE HELPER_FNCTS
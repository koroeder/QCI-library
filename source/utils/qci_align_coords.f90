MODULE QCIMINDIST
   USE QCIPREC   
   IMPLICIT NONE

   CONTAINS
      SUBROUTINE ALIGNXBTOA(XA, XB, NSIZE)
         USE QCIKEYS, ONLY: NATOMS
         USE QCIPREC
         IMPLICIT NONE 
         REAL(KIND = REAL64), INTENT(IN) :: XA(3*NATOMS)
         REAL(KIND = REAL64), INTENT(INOUT) :: XB(3*NATOMS)
         INTEGER, INTENT(IN) :: NSIZE
         REAL(KIND = REAL64) :: RA(3*NSIZE), RB(3*NSIZE)
         REAL(KIND = REAL64) :: CXA(3), CXB(3), RMAT(3,3)
         REAL(KIND = REAL64) :: DIST, R0(3), R1
         INTEGER :: I, J, K

         RA(1:3*NSIZE) = XA(1:3*NSIZE)
         RB(1:3*NSIZE) = XB(1:3*NSIZE)
         !centre RA
         CALL FIND_ORIGIN(NSIZE,RA,CXA)
         CALL MOVE_COORDS(NSIZE,RA,CXA)
         !centre RB
         CALL FIND_ORIGIN(NSIZE,RB,CXB)
         CALL MOVE_COORDS(NSIZE,RB,CXB)

         !align coordinates
         CALL FIND_ALIGNMENT(NSIZE, RB, RA, DIST, RMAT)

         DO I=1,NSIZE
            DO J=1,3
               R0(J) = RB(3*(I-1)+1)
            END DO
            DO J=1,3
               R1=0.0D0
               DO K=1,3
                  R1 = R1 + RMAT(J,K)*R0(K)
               ENDDO
               RB(3*(I-1)+J) = R1
            ENDDO
         ENDDO     

         ! move coordinates back to original frame position
         CALL MOVE_COORDS(NSIZE,RB,-CXB)
         ! write coordinates back to XB to be returned
         XB(1:3*NSIZE) = RB(1:3*NSIZE)

      END SUBROUTINE ALIGNXBTOA

      SUBROUTINE FIND_ORIGIN(NATOMS,X,CX)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND = REAL64), INTENT(OUT) :: CX(3)
         INTEGER :: I, J

         CX(1:3) = 0.0D0

         DO I=1,NATOMS
            DO J = 1,3
               CX(J) = CX(J) + X(3*(I-1) + J)
            END DO
         END DO
         CX(1:3) = CX(1:3)/DBLE(NATOMS)
      END SUBROUTINE FIND_ORIGIN

      SUBROUTINE MOVE_COORDS(NATOMS,X,CX)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
         REAL(KIND = REAL64), INTENT(IN) :: CX(3)
         INTEGER :: I, J, IDX
         
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1)+J
               X(IDX) = X(IDX) - CX(J)
            END DO
         END DO
      END SUBROUTINE MOVE_COORDS 
      
      SUBROUTINE FIND_ALIGNMENT(NATOMS, X, REFX, DIST, RMAT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS) 
         REAL(KIND = REAL64), INTENT(IN) :: REFX(3*NATOMS) 
         REAL(KIND = REAL64), INTENT(OUT) :: DIST
         REAL(KIND = REAL64), INTENT(OUT) :: RMAT(3,3)   
         INTEGER, PARAMETER :: LWORK=12       
         REAL(KIND = REAL64) :: CX(3)
         REAL(KIND = REAL64) :: QMAT(4,4) , XM, YM, ZM, XP, YP, ZP, DIAG(4), TEMPA(LWORK), MINV
         REAL(KIND = REAL64) :: Q1, Q2, Q3, Q4
         INTEGER :: I, J, JMIN, INFO

         JMIN = -1

         ! Analystic method based on quaternions and general angle-axis
         ! See: Kearsley, Acta Cryst. A, 45, 208-210, 1989
         !      Griffiths, Niblett and Wales, JCTC, 13, 4914-1931, 2017
         QMAT(1:4,1:4)=0.0D0
         DO I=1,NATOMS
            XM=REFX(3*(I-1)+1)-X(3*(I-1)+1)
            YM=REFX(3*(I-1)+2)-X(3*(I-1)+2)
            ZM=REFX(3*(I-1)+3)-X(3*(I-1)+3)
            XP=REFX(3*(I-1)+1)+X(3*(I-1)+1)
            YP=REFX(3*(I-1)+2)+X(3*(I-1)+2)
            ZP=REFX(3*(I-1)+3)+X(3*(I-1)+3)
            QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
            QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
            QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
            QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
            QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
            QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
            QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
            QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
            QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
            QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
         ENDDO
         QMAT(2,1)=QMAT(1,2)
         QMAT(3,1)=QMAT(1,3)
         QMAT(3,2)=QMAT(2,3)
         QMAT(4,1)=QMAT(1,4)
         QMAT(4,2)=QMAT(2,4)
         QMAT(4,3)=QMAT(3,4)
         !Eigendecomposition fo the quarternion
         CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,LWORK,INFO)
         MINV=1.0D100
         DO J=1,4
            IF (DIAG(J).LT.MINV) THEN
               JMIN=J
               MINV=DIAG(J)
            ENDIF
         ENDDO
         IF (MINV.LT.0.0D0) THEN
            IF (ABS(MINV).LT.1.0D-6) THEN
               MINV=0.0D0
            ELSE
               WRITE(*,'(A,G20.10,A)') 'newmindist> WARNING MINV is ',MINV,' change to absolute value'
               MINV=-MINV
            ENDIF
         ENDIF
         ! This is the Euclidean distance!
         DIST=SQRT(MINV)
         IF (JMIN.EQ.-1) THEN
            WRITE(*,*) "QMAT: ", QMAT
            WRITE(*,*) "MINV: ", MINV, " DIST: ", DIST
            WRITE(*,*) "DIAG: ", DIAG
         END IF
         ! Get the rotational matrix
         Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
         RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
         RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
         RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
         RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
         RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
         RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
         RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
         RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
         RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2             
      END SUBROUTINE FIND_ALIGNMENT
END MODULE QCIMINDIST
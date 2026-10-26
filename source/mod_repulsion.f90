MODULE REPULSION
   USE QCIPREC
   USE QCIKEYS
   IMPLICIT NONE
   INTEGER :: NREPCURR = 0
   INTEGER, ALLOCATABLE :: REPI(:), REPJ(:)
   INTEGER, ALLOCATABLE :: NREPI(:), NREPJ(:) 
   REAL(KIND = REAL64), ALLOCATABLE :: REPCUT(:), NREPCUT(:)
   
   INTEGER :: NNREPULSIVE, NREPULSIVE
   CONTAINS
      SUBROUTINE ALLOC_REP_VARS(NSIZE)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NSIZE 

         CALL DEALLOC_REP_VARS()
         ALLOCATE(REPI(NSIZE), REPJ(NSIZE), NREPI(NSIZE), NREPJ(NSIZE), &
                  REPCUT(NSIZE), NREPCUT(NSIZE))
      END SUBROUTINE ALLOC_REP_VARS


      SUBROUTINE DEALLOC_REP_VARS()
         IF (ALLOCATED(REPI)) DEALLOCATE(REPI)
         IF (ALLOCATED(REPJ)) DEALLOCATE(REPJ)
         IF (ALLOCATED(NREPI)) DEALLOCATE(NREPI)
         IF (ALLOCATED(NREPJ)) DEALLOCATE(NREPJ)
         IF (ALLOCATED(REPCUT)) DEALLOCATE(REPCUT)
         IF (ALLOCATED(NREPCUT)) DEALLOCATE(NREPCUT)        
      END SUBROUTINE DEALLOC_REP_VARS


      SUBROUTINE DOUBLE_ALLOC_REP()
         USE HELPER_FNCTS, ONLY: DOUBLE_INT_ARRAY, DOUBLE_REAL64_ARRAY
         IMPLICIT NONE
         ! call routines to increase array sizer (they preserve data)
         ! input is the old array size (NREPCURR) and the initialisation value for the new parts of the array
         CALL DOUBLE_INT_ARRAY(NREPCURR, REPI, -1)
         CALL DOUBLE_INT_ARRAY(NREPCURR, REPJ, -1)
         CALL DOUBLE_INT_ARRAY(NREPCURR, NREPI, -1)
         CALL DOUBLE_INT_ARRAY(NREPCURR, NREPJ, -1)
         CALL DOUBLE_REAL64_ARRAY(NREPCURR, REPCUT, -1.0D0)
         CALL DOUBLE_REAL64_ARRAY(NREPCURR, NREPCUT, -1.0D0)   
         !finally increase NREPCURR
         NREPCURR = 2*NREPCURR
      END SUBROUTINE DOUBLE_ALLOC_REP

      !SUBROUTINE CHECKREP(XYZ, NNSTART, NSTART)
      !   USE QCIKEYS, ONLY: NATOMS, NIMAGES
      !   USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
      !   IMPLICIT NONE
      !   REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))
      !   INTEGER, INTENT(IN) :: NSTART, NNSTART
      !   INTEGER :: J, K, NI, NJ
      !   REAL(KIND = REAL64) :: COMPARE, LDIST

      !   NNREPULSIVE = NNSTART

         ! loop over part of the repulsion list
      !   DO J = NSTART, NREPULSIVE
      !      COMPARE = (CHECKREPCUTOFF*REPCUT(J))
      !      NI = REPI(J)
      !      NJ = REPJ(J)
      !      ! now check for the distance between atoms in all images
      !      DO K=1,NIMAGES+2
      !         CALL DISTANCE_TWOATOMS(NATOMS, XYZ((K-1)*3*NATOMS+1:3*K*NATOMS) , NI, NJ, LDIST)
      !         IF (LDIST.LT.COMPARE) THEN
      !            NNREPULSIVE = NNREPULSIVE + 1
      !            NREPI(NNREPULSIVE) = NI
      !            NREPJ(NNREPULSIVE) = NJ
      !            NREPCUT(NNREPULSIVE) = REPCUT(J)
      !            EXIT
      !         END IF
      !      END DO
      !   END DO
      !   WRITE(*,'(A,2I8)') ' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE
      !END SUBROUTINE CHECKREP

      SUBROUTINE CHECKREP(XYZ, NNSTART, NSTART)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES, QCIINTREPMINSEP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))
         INTEGER, INTENT(IN) :: NSTART, NNSTART
         INTEGER :: J, K, OI, OJ
         LOGICAL :: NEAR, CHECKSEG
         REAL(KIND = REAL64) :: COMPARE2, G1(3), G2(3), DG(3), GT(3)
         REAL(KIND = REAL64) :: DSQ2, DGSQ, T
 
         NNREPULSIVE = NNSTART
 
         DO J = NSTART, NREPULSIVE
            COMPARE2 = (CHECKREPCUTOFF*REPCUT(J))**2
            ! GET_REPULSION_E never checks internal minima for these pairs
            CHECKSEG = (ABS(REPI(J)-REPJ(J)).GE.QCIINTREPMINSEP)
            OI = 3*(REPI(J)-1)
            OJ = 3*(REPJ(J)-1)
 
!            ! image 1
            G2(1:3) = XYZ(OI+1:OI+3) - XYZ(OJ+1:OJ+3)
            DSQ2 = DOT_PRODUCT(G2,G2)
            NEAR = (DSQ2.LT.COMPARE2)
 
            ! images 2..NIMAGES+2, plus the segment ending at each
            K = 2
           DO WHILE ((.NOT.NEAR).AND.(K.LE.NIMAGES+2))
               G1(1:3) = G2(1:3)
               G2(1:3) = XYZ((K-1)*3*NATOMS+OI+1:(K-1)*3*NATOMS+OI+3) &
                       - XYZ((K-1)*3*NATOMS+OJ+1:(K-1)*3*NATOMS+OJ+3)
               DSQ2 = DOT_PRODUCT(G2,G2)
               IF (DSQ2.LT.COMPARE2) THEN
                  NEAR = .TRUE.
               ELSE IF (CHECKSEG) THEN
                  ! closest approach on G1 + t*(G2-G1), t in (0,1)
                  DG(1:3) = G2(1:3) - G1(1:3)
                  DGSQ = DOT_PRODUCT(DG,DG)
                  IF (DGSQ.GT.0.0D0) THEN
                     T = -DOT_PRODUCT(G1,DG)/DGSQ
                    IF ((T.GT.0.0D0).AND.(T.LT.1.0D0)) THEN
                        GT(1:3) = G1(1:3) + T*DG(1:3)
                        NEAR = (DOT_PRODUCT(GT,GT).LT.COMPARE2)
                    END IF
                  END IF
               END IF
               K = K + 1
            END DO
 
            IF (NEAR) THEN
               NNREPULSIVE = NNREPULSIVE + 1
               NREPI(NNREPULSIVE) = REPI(J)
               NREPJ(NNREPULSIVE) = REPJ(J)
               NREPCUT(NNREPULSIVE) = REPCUT(J)
            END IF
         END DO
         WRITE(*,'(A,2I8)') ' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE
      END SUBROUTINE CHECKREP

END MODULE REPULSION

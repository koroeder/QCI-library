MODULE SBM_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: SBM_NCONST = 0
   INTEGER, ALLOCATABLE :: SBM_CONI(:), SBM_CONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: SBM_CONDISTREF(:)
   REAL(KIND = REAL64), ALLOCATABLE :: SBM_CONCUT(:)
   CHARACTER(LEN=30) :: SBMCONTACTFILE = "contacts.sbm"
   CONTAINS
      ! for SBM go model
      SUBROUTINE SBMMODEL_QCI_CONSTRAINTS(NATOMS)
         USE QCI_KEYS, ONLY: XSTART, XFINAL
         USE QCIFILEHANDLER, ONLY: GETUNIT, FILE_LENGTH
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY ! counter for constraints
         REAL(KIND = REAL64) :: DF, DS ! distance in final and initial structure
         LOGICAL :: YESNO ! do we have a contact file
         INTEGER :: CONUNIT ! unit for opening file
         INTEGER :: NADDCONSTR ! number of additional constraints
         INTEGER :: IDX1, IDX2 !indices read in from file

         ! check for additional constraints in file
         INQUIRE(FILE=SBMCONTACTFILE, EXIST=YESNO)
         IF (YESNO) THEN
            NADDCONSTR = FILE_LENGTH(SBMCONTACTFILE)
         ELSE
            NADDCONSTR = 0
         END IF
         ! allocate arrays
         SBM_NCONST = 2*NATOMS - 3 + NADDCONSTR
         CALL ALLOC_SBM_CONST()

         NDUMMY = 0
         ! add bond constraints - we have NATOMS - 1 bonds, all between adjacent atom indices
         DO J1=1,NATOMS-1
            NDUMMY = NDUMMY + 1
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, J1, J1+1, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, J1, J1+1, DF)
            SBM_CONI(NDUMMY) = J1
            SBM_CONJ(NDUMMY) = J1 + 1
            SBM_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            SBM_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
         END DO

         ! add angle constraints - we have NATOMS - 2 angles, between i and i+2
         DO J1=1,NATOMS-2
            NDUMMY = NDUMMY + 1
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, J1, J1+2, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, J1, J1+2, DF)
            SBM_CONI(NDUMMY) = J1
            SBM_CONJ(NDUMMY) = J1 + 2
            SBM_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            SBM_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
         END DO

         ! now add contacts in file provided
         IF (YESNO) THEN
            CONUNIT = GETUNIT()
            OPEN(CONUNIT,FILE=SBMCONTACTFILE,STATUS='OLD')
            DO J1=1,NADDCONSTR
               READ(CONUNIT,*) IDX1, IDX2
               NDUMMY = NDUMMY + 1
               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF) 
               SBM_CONI(NDUMMY) = IDX1
               SBM_CONJ(NDUMMY) = IDX2
               SBM_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
               SBM_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0                             
            END DO
         END IF
         WRITE(*,*) " sbm_constraints> Identified ", SBM_NCONST, " constraints"
         WRITE(*,*) "                  Bonds: ", NATOMS-1, ", angles: ", NATOMS-2, ", additional constraints: ", NADDCONSTR
      END SUBROUTINE SBMMODEL_QCI_CONSTRAINTS

      SUBROUTINE ALLOC_SBM_CONST()
         CALL DEALLOC_SBM_CONST()
         ALLOCATE(SBM_CONI(SBM_NCONST))
         ALLOCATE(SBM_CONJ(SBM_NCONST))
         ALLOCATE(SBM_CONDISTREF(SBM_NCONST))
         ALLOCATE(SBM_CONCUT(SBM_NCONST))
      END SUBROUTINE ALLOC_SBM_CONST

      SUBROUTINE DEALLOC_SBM_CONST()
         IF (ALLOCATED(SBM_CONI)) DEALLOCATE(SBM_CONI)
         IF (ALLOCATED(SBM_CONJ)) DEALLOCATE(SBM_CONJ)
         IF (ALLOCATED(SBM_CONDISTREF)) DEALLOCATE(SBM_CONDISTREF)
         IF (ALLOCATED(SBM_CONCUT)) DEALLOCATE(SBM_CONCUT)
      END SUBROUTINE DEALLOC_SBM_CONST

      SUBROUTINE DISTANCE_TWOATOMS(NATOMS, X, IDX1, IDX2, DIST)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64) :: X(3*NATOMS)
         INTEGER, INTENT(IN) :: IDX1, IDX2
         REAL(KIND = REAL64) :: DIST
         REAL(KIND = REAL64) :: XYZ1(3), XYZ2(3)
         INTEGER :: J
         DO J = 1,3
            XYZ1(J) = X(3*(IDX1-1)+J)
            XYZ2(J) = X(3*(IDX2-1)+J)
         END DO
         DIST = SQRT((XYZ1(1)-XYZ2(1))**2 + (XYZ1(2)-XYZ2(2))**2 + (XYZ1(3)-XYZ2(3))**2)
      END SUBROUTINE DISTANCE_TWOATOMS

END MODULE SBM_CONSTRAINTS
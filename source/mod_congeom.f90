MODULE CONGEOM
   USE QCIPREC, ONLY: REAL64
   INTEGER :: NGEOMCONST
   INTEGER, ALLOCATABLE :: GEOMCONI(:), GEOMCONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: GEOMCONDISTREF(:), GEOMCONCUT(:)
   LOGICAL :: USEENDPOINTS

   INTEGER :: NADDCONSTR
   INTEGER, ALLOCATABLE :: FILE_CONI(:), FILE_CONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: FILE_CONDISTREF(:), FILE_CONCUT(:)

   CONTAINS

      ! allocate arrays
      SUBROUTINE ALLOC_CONGEOM()
         CALL DEALLOC_CONGEOM()
         ALLOCATE(GEOMCONI(NGEOMCONST))
         ALLOCATE(GEOMCONJ(NGEOMCONST))
         ALLOCATE(GEOMCONDISTREF(NGEOMCONST))
         ALLOCATE(GEOMCONCUT(NGEOMCONST))
      END SUBROUTINE ALLOC_CONGEOM

      ! deallocate arrays
      SUBROUTINE DEALLOC_CONGEOM()
         IF (ALLOCATED(GEOMCONI)) DEALLOCATE(GEOMCONI)
         IF (ALLOCATED(GEOMCONJ)) DEALLOCATE(GEOMCONJ)
         IF (ALLOCATED(GEOMCONDISTREF)) DEALLOCATE(GEOMCONDISTREF)
         IF (ALLOCATED(GEOMCONCUT)) DEALLOCATE(GEOMCONCUT)
      END SUBROUTINE DEALLOC_CONGEOM
      
      SUBROUTINE READ_GEOMS()
         USE QCIKEYS, ONLY: NATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONGEOM, CONGEOM, GEOMFILE
         USE QCIFILEHANDLER, ONLY: FILE_LENGTH, GETUNIT
         IMPLICIT NONE 
         LOGICAL :: YESNO
         LOGICAL :: NDUMMY
         INTEGER :: FLENGTH, GEOMUNIT
         INTEGER :: I, J

         USEENDPOINTS = .FALSE.
         INQUIRE(FILE=GEOMFILE, EXIST=YESNO)
         IF (YESNO) THEN
            FLENGTH = FILE_LENGTH(GEOMFILE)
            IF (FLENGTH.LT.2*NATOMS) THEN
               USEENDPOINTS = .TRUE.
               WRITE(*,*) " read_geoms> Not enough geometries found in file, derive constraints from endpoints"
            END IF 
         ELSE
            USEENDPOINTS = .TRUE.
            WRITE(*,*) " read_geoms> Geometry file ", ADJUSTL(TRIM(GEOMFILE)), " not found - derive constraints from endpoints"
            RETURN
         END IF 

         !open topology
         GEOMUNIT = GETUNIT()
         OPEN(GEOMUNIT,FILE=GEOMFILE,STATUS='OLD')
         IF (MOD(FLENGTH,NATOMS).EQ.0) THEN
            NCONGEOM = FLENGTH/NATOMS
            ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
            DO J=1,NCONGEOM
               READ(GEOMUNIT,*) (CONGEOM(J,I), I=1,3*NATOMS)
            END DO
         ELSE
            NCONGEOM = (FLENGTH-1)/NATOMS
            ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
            READ(GEOMUNIT,*)
            DO J=1,NCONGEOM
               READ(GEOMUNIT,*) (CONGEOM(J,I), I=1,3*NATOMS)
            END DO           
         END IF
         CLOSE(GEOMUNIT)
      END SUBROUTINE READ_GEOMS

      ! create constraints from endpoints
      SUBROUTINE CREATE_FROM_ENDPOINTS(MODIFYTOL)
         USE QCIKEYS, ONLY:  NATOMS
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE QCIPREC
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: QCICONSEP, QCICONSTRAINTTOL, NCONGEOM, QCICONCUT
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: MODIFYTOL
         INTEGER :: J1, J2, J3, MAXCOUNT
         REAL(KIND = REAL64) :: DS, DF
         LOGICAL :: ACCEPTCONST
         INTEGER :: NDUMMY
         INTEGER :: DCONI(2*QCICONSEP*NATOMS), DCONJ(2*QCICONSEP*NATOMS)
         REAL(KIND = REAL64) :: DCONDISTREF(2*QCICONSEP*NATOMS), DCONCUT(2*QCICONSEP*NATOMS)
         REAL(KIND = REAL64) :: LOCALTOL

         NDUMMY = 0

         LOCALTOL = MODIFYTOL * QCICONSTRAINTTOL
         DO J1=1,NATOMS-1
            MAXCOUNT = MIN(J1+QCICONSEP,NATOMS)
            DO J2=J1+1,MAXCOUNT
               ACCEPTCONST = .TRUE.
               DO J3=1,NCONGEOM
                  CALL DISTANCE_TWOATOMS(NATOMS, XSTART, J1, J2, DS)
                  CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, J1, J2, DF)
                  ! if the distance is too big in any image, do not use this as constraint
                  IF ((DS.GT.QCICONCUT).OR.(DF.GT.QCICONCUT)) THEN
                     ACCEPTCONST = .FALSE.
                     EXIT
                  END IF 
                  ! we do not want the variance between geometries to be too big
                  IF (ABS(DF-DS).GT.LOCALTOL) THEN
                     ACCEPTCONST = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (ACCEPTCONST) THEN
                  NDUMMY = NDUMMY + 1
                  DCONI(NDUMMY)=J2
                  DCONJ(NDUMMY)=J3
                  DCONDISTREF(NDUMMY)=(DF+DS)/2.0D0 
                  DCONCUT(NDUMMY)=ABS(DF-DS)/2.0D0
               END IF
            END DO
         END DO
         NGEOMCONST = NDUMMY
         CALL ALLOC_CONGEOM()
         GEOMCONI(1:NGEOMCONST) = DCONI(1:NGEOMCONST)
         GEOMCONJ(1:NGEOMCONST) = DCONJ(1:NGEOMCONST)
         GEOMCONDISTREF(1:NGEOMCONST) = DCONDISTREF(1:NGEOMCONST)
         GEOMCONCUT(1:NGEOMCONST) = DCONCUT(1:NGEOMCONST)           
      END SUBROUTINE CREATE_FROM_ENDPOINTS

      ! create constraints from input geometries
      SUBROUTINE CREATE_FROM_GEOMETRIES(MODIFYTOL)
         USE QCIKEYS, ONLY: NATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: QCICONSEP, QCICONSTRAINTTOL, NCONGEOM, CONGEOM, QCICONCUT
         USE QCIPREC
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: MODIFYTOL
         INTEGER :: J1, J2, J3, MAXCOUNT
         REAL(KIND = REAL64) :: DSMIN, DSMAX, DS, DSMEAN
         LOGICAL :: ACCEPTCONST
         INTEGER :: NDUMMY
         INTEGER :: DCONI(2*QCICONSEP*NATOMS), DCONJ(2*QCICONSEP*NATOMS)
         REAL(KIND = REAL64) :: X(3*NATOMS)
         REAL(KIND = REAL64) :: DCONDISTREF(2*QCICONSEP*NATOMS), DCONCUT(2*QCICONSEP*NATOMS)
         REAL(KIND = REAL64) :: LOCALTOL

         NDUMMY = 0

         LOCALTOL = MODIFYTOL * QCICONSTRAINTTOL
         DO J1=1,NATOMS-1
            MAXCOUNT = MIN(J1+QCICONSEP,NATOMS)
            DO J2=J1+1,MAXCOUNT
               DSMIN = 1.0D100
               DSMAX = -1.0D100
               DSMEAN = 0.0D0
               ACCEPTCONST = .TRUE.
               DO J3=1,NCONGEOM
                  X(1:3*NATOMS) = CONGEOM(J3,1:3*NATOMS)
                  CALL DISTANCE_TWOATOMS(NATOMS, X, J1, J2, DS)
                  ! if the distance is too big in any image, do not use this as constraint
                  IF (DS.GT.QCICONCUT) THEN
                     ACCEPTCONST = .FALSE.
                     EXIT
                  END IF 
                  IF (DS.GT.DSMAX) DSMAX = DS
                  IF (DS.LT.DSMIN) DSMIN = DS
                  ! we do not want the variance between geometries to be too big
                  IF (ABS(DSMAX-DSMIN).GT.LOCALTOL) THEN
                     ACCEPTCONST = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (ACCEPTCONST) THEN
                  NDUMMY = NDUMMY + 1
                  DCONI(NDUMMY)=J2
                  DCONJ(NDUMMY)=J3
                  DCONDISTREF(NDUMMY)=(DSMAX+DSMIN)/2.0D0 
                  DCONCUT(NDUMMY)=(DSMAX-DSMIN)/2.0D0
               END IF
            END DO
         END DO
         NGEOMCONST = NDUMMY
         CALL ALLOC_CONGEOM()
         GEOMCONI(1:NGEOMCONST) = DCONI(1:NGEOMCONST)
         GEOMCONJ(1:NGEOMCONST) = DCONJ(1:NGEOMCONST)
         GEOMCONDISTREF(1:NGEOMCONST) = DCONDISTREF(1:NGEOMCONST)
         GEOMCONCUT(1:NGEOMCONST) = DCONCUT(1:NGEOMCONST)         
      END SUBROUTINE CREATE_FROM_GEOMETRIES
    
      ! add constraints list
      SUBROUTINE ADD_CONSTRAINT_LIST()
         USE QCIKEYS, ONLY: NATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: CONSTRFILE
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE QCIFILEHANDLER, ONLY: FILE_LENGTH, GETUNIT
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE QCIPREC
         IMPLICIT NONE
         REAL(KIND = REAL64) :: DS, DF
         LOGICAL :: YESNO
         INTEGER :: CONUNIT, NDUMMY, J1, IDX1, IDX2

         INQUIRE(FILE=CONSTRFILE, EXIST=YESNO)
         IF (YESNO) THEN
            NADDCONSTR = FILE_LENGTH(CONSTRFILE)
         ELSE
            NADDCONSTR = 0
            RETURN
         END IF
         CALL ALLOC_ADDCONST()
         ! now add additional constraints from file provided
         IF (YESNO) THEN
            CONUNIT = GETUNIT()
            OPEN(CONUNIT,FILE=CONSTRFILE,STATUS='OLD')
            NDUMMY = 0
            DO J1=1,NADDCONSTR
               READ(CONUNIT,*) IDX1, IDX2
               NDUMMY = NDUMMY + 1
               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF) 
               IF (IDX1.LT.IDX2) THEN
                  FILE_CONI(NDUMMY) = IDX1
                  FILE_CONJ(NDUMMY) = IDX2
               ELSE
                  FILE_CONI(NDUMMY) = IDX2
                  FILE_CONJ(NDUMMY) = IDX1
               END IF
               FILE_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
               FILE_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0                             
            END DO
            CLOSE(CONUNIT)
         END IF
      END SUBROUTINE ADD_CONSTRAINT_LIST

      SUBROUTINE ALLOC_ADDCONST()
         CALL DEALLOC_ADDCONST()
         ALLOCATE(FILE_CONI(NADDCONSTR))
         ALLOCATE(FILE_CONJ(NADDCONSTR))
         ALLOCATE(FILE_CONDISTREF(NADDCONSTR))
         ALLOCATE(FILE_CONCUT(NADDCONSTR))
      END SUBROUTINE ALLOC_ADDCONST

      SUBROUTINE DEALLOC_ADDCONST()
         IF (ALLOCATED(FILE_CONI)) DEALLOCATE(FILE_CONI)
         IF (ALLOCATED(FILE_CONJ)) DEALLOCATE(FILE_CONJ)
         IF (ALLOCATED(FILE_CONDISTREF)) DEALLOCATE(FILE_CONDISTREF)
         IF (ALLOCATED(FILE_CONCUT)) DEALLOCATE(FILE_CONCUT)  
      END SUBROUTINE DEALLOC_ADDCONST
END MODULE CONGEOM
MODULE QCI_LINEAR
   USE QCIPREC
   IMPLICIT NONE
   ! number of atoms for QCIlinear
   INTEGER :: NQCILINEAR = 0
   ! cutoff for QCIlinear treatment
   REAL(KIND=REAL64) :: LINEARCUT = 0.05D0
   ! list of linear atoms
   INTEGER, ALLOCATABLE :: LINEARATOMS(:)
   ! file name linear atoms
   CHARACTER(25) :: LINEARFILE = "QCIlinear"

   CONTAINS

      SUBROUTINE GET_LINEAR_ATOMS()
         USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET
         USE QCIFILEHANDLER, ONLY: FILE_LENGTH, GETUNIT        
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_ATOM_DIFF_IMAGES
         USE AMBER_CONSTRAINTS, ONLY: AMBERBB => BACKBONE
         USE HIRE_CONSTRAINTS, ONLY: HIREBB => BACKBONE
         INTEGER :: NDUMMY, DUMMY
         INTEGER :: LINEART(NATOMS)
         REAL(KIND=REAL64) :: DIST
         INTEGER :: J1
         INTEGER :: LINUNIT
         LOGICAL :: YESNO

         LINEART(1:NATOMS) = 0
         INQUIRE(FILE=LINEARFILE, EXIST=YESNO)
         IF (YESNO) THEN
            WRITE(*,*) " get_linear_atoms> Reading in linear atoms from file"
            LINUNIT = GETUNIT()
            OPEN(LINUNIT,FILE=LINEARFILE,STATUS='OLD')
            READ(LINUNIT, '(I6)') NDUMMY
            DO J1=1,NDUMMY
               READ(LINUNIT, '(I6)') DUMMY
               LINEART(DUMMY) = 1
            END DO
            CLOSE(LINUNIT)
         END IF

         IF (LINEARBBT) THEN
            IF (QCIAMBERT.OR.QCIHIRET) THEN
               DO J1=1,NATOMS
                  IF (ISBBATOM(J1)) LINEART(J1) = 1
               END DO
            ELSE
               WRITE(*,*) "WARNING: Linear backbone interpolation set, but neither AMBER not HiRE are used"
            END IF
         END IF

         DO J1=1,NATOMS
            CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, XSTART, XFINAL, J1, DIST)
            IF (DIST.LT.LINEARCUT) THEN
               LINEART(J1) = 1
            END IF
         END DO

         NQCILINEAR = SUM(LINEART)
         CALL ALLOC_QCI_LINEAR()
         INLINLIST(1:NATOMS) = .FALSE.
         DUMMY=0
         DO J1=1,NATOMS
            IF (LINEART(J1).EQ.1) THEN
               DUMMY = DUMMY + 1
               LINEARATOMS(DUMMY) = J1
               INLINLIST(J1) = .TRUE.
            END IF
         END DO
         WRITE(*,*) " linear list: ", LINEARATOMS(1:DUMMY)
      END SUBROUTINE GET_LINEAR_ATOMS
   
      SUBROUTINE ALLOC_QCI_LINEAR()
         USE QCIKEYS, ONLY: NATOMS, INLINLIST
         CALL DEALLOC_QCI_LINEAR
         ALLOCATE(LINEARATOMS(NQCILINEAR))
         ALLOCATE(INLINLIST(NATOMS))
      END SUBROUTINE ALLOC_QCI_LINEAR

      SUBROUTINE DEALLOC_QCI_LINEAR()
         USE QCIKEYS, ONLY: INLINLIST
         IF (ALLOCATED(LINEARATOMS)) DEALLOCATE(LINEARATOMS)
         IF (ALLOCATED(INLINLIST)) DEALLOCATE(INLINLIST)
      END SUBROUTINE DEALLOC_QCI_LINEAR

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !IDEA FOR QUASI-RIGID BODY LINEAR LIST
      !DO NOT USE 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !SUBROUTINE DETECT_LINEAR()
      !   USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET       
      !   USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
      !   USE HELPER_FNCTS, ONLY: DISTANCE_ATOM_DIFF_IMAGES
      !   USE QCI_CONSTRAINT_KEYS
      !   USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL

      !   IMPLICIT NONE

       !  REAL(KIND=REAL64), PARAMETER :: TOLERANCE
       !  INTEGER :: J1 

       !  CHECK_RIGIDITY = .TRUE.
       !  DO J1 = 1, NCONSTRAINT
       !     A = CONI(J1)
       !     B = CONJ(J1)
       !  ! Check if either atom is in this group
       !  IF (.NOT. ANY(ATOM_LIST(1:NATOMS) == A) .AND. &
       !      .NOT. ANY(ATOM_LIST(1:NATOMS) == B)) THEN
       !  CYCLE
       !  END IF
       !  
       !  CHARACTER(LEN=4), INTENT(IN) :: ATNAME1, ATNAME2
       !  INTEGER, INTENT(IN) :: RESID
       !  INTEGER :: ID1, ID2 1
      !
      !   CALL GET_ATOMID(ATNAME1,RESID,ID1)
      !   CALL GET_ATOMID(ATNAME2,RESID,ID2)

         ! Bond length in start image
!         DX = START_XYZ(B,1) - START_XYZ(A,1)
!         DY = START_XYZ(B,2) - START_XYZ(A,2)
!         DZ = START_XYZ(B,3) - START_XYZ(A,3)
!         D_START = SQRT(DX*DX + DY*DY + DZ*DZ)
!    
!    ! Bond length in end image
!    DX = END_XYZ(B,1) - END_XYZ(A,1)
!    DY = END_XYZ(B,2) - END_XYZ(A,2)
!    DZ = END_XYZ(B,3) - END_XYZ(A,3)
!    D_END = SQRT(DX*DX + DY*DY + DZ*DZ)
!    
!    IF (ABS(D_END - D_START) > TOLERANCE) THEN
!      CHECK_RIGIDITY = .FALSE.
!      RETURN
!    END IF
!  END DO

         
     
     
 !     END SUBROUTINE


END MODULE QCI_LINEAR
MODULE HIRE_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: NBACKBONE
   INTEGER :: NRES
   INTEGER, ALLOCATABLE :: BACKBONE(:)
   INTEGER :: HIRE_NCONST = 0
   INTEGER, ALLOCATABLE :: HIRE_CONI(:), HIRE_CONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: HIRE_CONDISTREF(:)
   REAL(KIND = REAL64), ALLOCATABLE :: HIRE_CONCUT(:)
   CHARACTER(LEN=30) :: HIRECONSTRFILE = "constraintfile"
   CHARACTER(LEN=25) :: HIRETOPFILE = "parameters.top"
   INTEGER, ALLOCATABLE :: BONDS(:,:)
   INTEGER, ALLOCATABLE :: ANGLES(:,:) 
   INTEGER :: NBOND = 0
   INTEGER :: NANGLE = 0
   INTEGER, ALLOCATABLE :: RESFINAL(:), RESSTART(:)
   CHARACTER(LEN=4), ALLOCATABLE :: HIRE_NAMES(:), RESNAMES(:)
 
   CONTAINS
      ! from HiRE topology
      SUBROUTINE HIRE_QCI_CONSTRAINTS()
         USE QCI_KEYS, ONLY: XSTART, XFINAL
         USE QCIFILEHANDLER, ONLY: GETUNIT, FILE_LENGTH
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY ! counter for constraints
         REAL(KIND = REAL64) :: DF, DS ! distance in final and initial structure
         LOGICAL :: YESNO ! do we have a contact file
         INTEGER :: CONUNIT ! unit for opening file
         INTEGER :: NADDCONSTR ! number of additional constraints
 
         INTEGER :: IDX1, IDX2 !indices for atoms

         ! parse topology
         CALL READ_TOPOLOGY()
         ! check for additional constraints in file
         INQUIRE(FILE=HIRECONSTRFILE, EXIST=YESNO)
         IF (YESNO) THEN
            NADDCONSTR = FILE_LENGTH(HIRECONTACTFILE)
         ELSE
            NADDCONSTR = 0
         END IF

         ! allocate arrays
         HIRE_NCONST = NBOND + NANGLE + NADDCONSTR
         ! allocate the required arrays
         CALL ALLOC_HIRE_CONSTR()
         NDUMMY = 0
         ! add bond constraints
         DO J1=1,NBOND
            NDUMMY = NDUMMY + 1
            IDX1 = BONDS(J1,1)/3+1
            IDX2 = BONDS(J1,2)/3+1
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF)
            IF (IDX1.LT.IDX2) THEN
               HIRE_CONI(NDUMMY) = IDX1
               HIRE_CONJ(NDUMMY) = IDX2
            ELSE
               HIRE_CONI(NDUMMY) = IDX2
               HIRE_CONJ(NDUMMY) = IDX1
            END IF
            HIRE_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            HIRE_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
         END DO
         ! add angle constraints
         DO J1=1,NANGLE
            NDUMMY = NDUMMY + 1
            IDX1 = ANGLES(J1,1)/3+1
            IDX2 = ANGLES(J1,2)/3+1
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF)
            IF (IDX1.LT.IDX2) THEN
               HIRE_CONI(NDUMMY) = IDX1
               HIRE_CONJ(NDUMMY) = IDX2
            ELSE
               HIRE_CONI(NDUMMY) = IDX2
               HIRE_CONJ(NDUMMY) = IDX1
            END IF
            HIRE_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            HIRE_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
         END DO
         ! now add additional constraints from file provided
         IF (YESNO) THEN
            CONUNIT = GETUNIT()
            OPEN(CONUNIT,FILE=AMBERCONSTRFILE,STATUS='OLD')
            DO J1=1,NADDCONSTR
               READ(CONUNIT,*) IDX1, IDX2
               NDUMMY = NDUMMY + 1
               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF) 
               IF (IDX1.LT.IDX2) THEN
                  HIRE_CONI(NDUMMY) = IDX1
                  HIRE_CONJ(NDUMMY) = IDX2
               ELSE
                  HIRE_CONI(NDUMMY) = IDX2
                  HIRE_CONJ(NDUMMY) = IDX1
               END IF
               HIRE_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
               HIRE_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0                             
            END DO
            CLOSE(CONUNIT)
         END IF
         DEALLOCATE(BONDS,ANGLES)
         WRITE(*,*) " HiRE_constraints> Identified ", HIRE_NCONST, " constraints"
         WRITE(*,*) "                   Bonds: ", NBOND, ", angles: ", NANGLE, ", additional constraints: ", NADDCONSTR

      END SUBROUTINE HIRE_QCI_CONSTRAINTS

      SUBROUTINE GET_BACKBONE_HIRE(NATOMS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY
         INTEGER :: DUMMY_BB(NATOMS)
         CHARACTER(LEN=4) :: ATNAME
         INTEGER :: J1
         DO J1=1,NATOMS
            ATNAME = ADJUSTL(TRIM(HIRE_NAMES(J1)))
            IF ((ATNAME.EQ.'C').OR.(ATNAME.EQ.'O').OR.(ATNAME.EQ.'P').OR.(ATNAME.EQ.'R4')) THEN
               NDUMMY = NDUMMY + 1
               DUMMY_BB(NDUMMY) = J1
            END IF
         END DO

         NBACKBONE = NDUMMY
         ALLOCATE(BACKBONE(NBACKBONE))
         BACKBONE(1:NBACKBONE) = DUMMY_BB(1:NBACKBONE)
      END SUBROUTINE GET_BACKBONE_HIRE


      ! allocate HiRE constraints array
      SUBROUTINE ALLOC_HIRE_CONSTR()
         CALL DEALLOC_HIRE_CONSTR()
         ALLOCATE(HIRE_CONI(HIRE_NCONST))
         ALLOCATE(HIRE_CONJ(HIRE_NCONST))
         ALLOCATE(HIRE_CONDISTREF(HIRE_NCONST))
         ALLOCATE(HIRE_CONCUT(HIRE_NCONST))
      END SUBROUTINE ALLOC_HIRE_CONSTR

      ! dealloc HIRE constraints arrays
      SUBROUTINE DEALLOC_HIRE_CONSTR()
         IF (ALLOCATED(HIRE_CONI)) DEALLOCATE(HIRE_CONI)
         IF (ALLOCATED(HIRE_CONJ)) DEALLOCATE(HIRE_CONJ)
         IF (ALLOCATED(HIRE_CONDISTREF)) DEALLOCATE(HIRE_CONDISTREF)
         IF (ALLOCATED(HIRE_CONCUT)) DEALLOCATE(HIRE_CONCUT)
      END SUBROUTINE DEALLOC_HIRE_CONSTR

      SUBROUTINE READ_TOPOLOGY()
         USE QCIFILEHANDLER, ONLY: GETUNIT
         USE HELPER_FNCTS, ONLY: READ_LINE
         IMPLICIT NONE

         LOGICAL :: YESNO, ENDFILET
         INTEGER :: TOPUNIT
         INTEGER :: NPARTICLES, NTYPEP, NDIHS, NRES, NUMBND, NUMANG, NPTRA, &
                     NCHAINS, NOPT
         INTEGER :: IEND, J
         INTEGER, PARAMETER :: NWORDS = 15
         CHARACTER(LEN=30), DIMENSION(NWORDS) :: WORDSLINE
         CHARACTER(LEN=250) :: THISLINE
         CHARACTER(LEN=30) :: CATNAME
      
         INTEGER, ALLOCATABLE :: DUMMY(:,:), DUMMY1(:), DUMMY5(:), DUMMY6(:), DUMMY7(:)
         REAL(KIND = REAL64), ALLOCATABLE :: DUMMY1R(:), DUMMY2(:), DUMMY3(:), DUMMY4(:)

         !check topology exists
         INQUIRE(FILE=HIRETOPFILE, EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(*,*) " qci_hire_constr> Cannot locate file ", HIRETOPFILE
            STOP
         END IF
         !open topology
         TOPUNIT = GETUNIT()
         OPEN(TOPUNIT,FILE=HIRETOPFILE,STATUS='OLD')
      
         !First line is topology name - ignore it
         READ(TOPUNIT, *)
         !Now we iterate, the format is:
         ! SECTION CATNAME
         ! entries ...
         ENDFILET = .FALSE.
         DO WHILE (.NOT.(ENDFILET))
            READ(TOPUNIT, '(A)', IOSTAT=IEND) THISLINE
            IF (IEND.LT.0) THEN
               ENDFILET = .TRUE.
               CYCLE
            ENDIF
            CALL READ_LINE(THISLINE,NWORDS,WORDSLINE)
            CATNAME = WORDSLINE(2)
            SELECT CASE (CATNAME)
               CASE("DEFINITIONS")
                  READ(TOPUNIT,'(12I6)') NPARTICLES, NTYPEP, NBOND, NANGLE, &
                                          NDIHS, NRES, NUMBND, NUMANG, NPTRA, &
                                          NCHAINS 
                  NOPT = 3 * NPARTICLES
                  !allocation of arrays from bond, angles and dihedral modules
                  ALLOCATE(BONDS(NBOND,2))
                  ALLOCATE(ANGLES(NANGLE,2))
                  ALLOCATE(RESSTART(NRES),RESFINAL(NRES),RESNAMES(NRES))
                  ALLOCATE(HIRE_NAMES(NPARTICLES))
                  !dummies not needed
                  ALLOCATE(DUMMY(NCHAINS,2))
                  ALLOCATE(DUMMY1(NPARTICLES),DUMMY2(NUMBND),DUMMY3(NUMANG),DUMMY4(NPTRA))
                  ALLOCATE(DUMMY1R(NPARTICLES))
                  ALLOCATE(DUMMY5(NBONDS),DUMMY6(NANGLES),DUMMY7(NDIHS))
               CASE("PARTICLE_NAMES")
                  READ(TOPUNIT,'(20A4)') (HIRE_NAMES(J), J=1,NPARTICLES)
               CASE("RESIDUE_LABELS")
                  READ(TOPUNIT,'(20A4)') (RESNAMES(J), J=1,NRES)
               CASE("RESIDUE_POINTER")
                  READ(TOPUNIT,'(12I6)') (RESSTART(J),RESFINAL(J), J=1,NRES)
               CASE("CHAIN_POINTER")
                  READ(TOPUNIT,'(12I6)') &
                           (DUMMY(J,1),DUMMY(J,2), J=1,NCHAINS)
               CASE("PARTICLE_MASSES")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY1R(J), J=1,NPARTICLES)
               CASE("PARTICLE_TYPE")
                  READ(TOPUNIT,'(12I6)') (DUMMY1(J), J=1,NPARTICLES)
               CASE("CHARGES")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY1R(J), J=1,NPARTICLES)             
               CASE("BOND_FORCE_CONSTANT")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY2(J), J=1,NUMBND)
               CASE("BOND_EQUIL_VALUE") 
                  READ(TOPUNIT,'(5E16.8)') (DUMMY2(J), J=1,NUMBND)
               CASE("ANGLE_FORCE_CONSTANT")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY3(J), J=1,NUMANG)
               CASE("ANGLE_EQUIL_VALUE")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY3(J), J=1,NUMANG) 
               CASE("DIHEDRAL_FORCE_CONSTANT")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY4(J), J=1,NPTRA) 
               CASE("DIHEDRAL_PERIODICITY")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY4(J), J=1,NPTRA) 
               CASE("DIHEDRAL_PHASE")
                  READ(TOPUNIT,'(5E16.8)') (DUMMY4(J), J=1,NPTRA) 
               CASE("BONDS")
                  READ(TOPUNIT,'(12I6)') (BONDS(J,1), BONDS(J,2), DUMMY5(J), J=1,NBONDS)
               CASE("ANGLES")
                  READ(TOPUNIT,'(12I6)') &
                              (ANGLES(J,1), DUMMY6(J), ANGLES(J,2), DUMMY6(J), J=1,NANGLES)
               CASE("DIHEDRALS") 
                  READ(TOPUNIT,'(12I6)') &
                              (DUMMY7(J), DUMMY7(J), DUMMY7(J), DUMMY7(J), DUMMY7(J), J=1,NDIHS)
            END SELECT
         ENDDO
         CLOSE(TOPUNIT)
         DEALLOCATE(DUMMY,DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY1R,DUMMY5,DUMMY6,DUMMY7)
      END SUBROUTINE READ_TOPOLOGY

END MODULE HIRE_CONSTRAINTS
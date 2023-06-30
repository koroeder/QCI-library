MODULE AMBER_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: AMBER_NCONST = 0
   INTEGER, ALLOCATABLE :: AMBER_CONI(:), AMBER_CONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: AMBER_CONDISTREF(:)
   REAL(KIND = REAL64), ALLOCATABLE :: AMBER_CONCUT(:)
   CHARACTER(LEN=30) :: AMBERCONSTRFILE = "constraintfile"
   CHARACTER(LEN=25) :: TOPFILENAME = "coords.prmtop"
   INTEGER, ALLOCATABLE :: BONDS(:,:)
   INTEGER, ALLOCATABLE :: ANGLES(:,:)
   INTEGER :: NBOND = 0
   INTEGER :: NANGLE = 0

   CONTAINS
      ! from AMBER
      SUBROUTINE AMBER_QCI_CONSTRAINTS()
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
         CALL PARSE_TOPOLOGY()
         ! check for additional constraints in file
         INQUIRE(FILE=AMBERCONSTRFILE, EXIST=YESNO)
         IF (YESNO) THEN
            NADDCONSTR = FILE_LENGTH(AMBERCONTACTFILE)
         ELSE
            NADDCONSTR = 0
         END IF
         ! allocate arrays
         AMBER_NCONST = NBOND + NANGLE + NADDCONSTR
         ! allocate the required arrays
         CALL ALLOC_AMBER_CONSTR()
         NDUMMY = 0
         ! add bond constraints
         DO J1=1,NBOND
            NDUMMY = NDUMMY + 1
            IDX1 = BONDS(J1,1)
            IDX2 = BONDS(J1,2)
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF)
            IF (IDX1.LT.IDX2) THEN
               AMBER_CONI(NDUMMY) = IDX1
               AMBER_CONJ(NDUMMY) = IDX2
            ELSE
               AMBER_CONI(NDUMMY) = IDX2
               AMBER_CONJ(NDUMMY) = IDX1
            END IF
            AMBER_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            AMBER_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
         END DO
         ! add angle constraints
         DO J1=1,NBOND
            NDUMMY = NDUMMY + 1
            IDX1 = ANGLES(J1,1)
            IDX2 = ANGLES(J1,2)
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, IDX1, IDX2, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, IDX1, IDX2, DF)
            IF (IDX1.LT.IDX2) THEN
               AMBER_CONI(NDUMMY) = IDX1
               AMBER_CONJ(NDUMMY) = IDX2
            ELSE
               AMBER_CONI(NDUMMY) = IDX2
               AMBER_CONJ(NDUMMY) = IDX1
            END IF
            AMBER_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
            AMBER_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0
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
                  AMBER_CONI(NDUMMY) = IDX1
                  AMBER_CONJ(NDUMMY) = IDX2
               ELSE
                  AMBER_CONI(NDUMMY) = IDX2
                  AMBER_CONJ(NDUMMY) = IDX1
               END IF
               AMBER_CONDISTREF(NDUMMY) = (DF+DS)/2.0D0
               AMBER_CONCUT(NDUMMY) = ABS(DF-DS)/2.0D0                             
            END DO
            CLOSE(CONUNIT)
         END IF
         DEALLOCATE(BONDS,ANGLES)
         WRITE(*,*) " amber_constraints> Identified ", AMBER_NCONST, " constraints"
         WRITE(*,*) "                    Bonds: ", NBOND, ", angles: ", NANGLE, ", additional constraints: ", NADDCONSTR

      END SUBROUTINE AMBER_QCI_CONSTRAINTS

      ! allocate amber constraints array
      SUBROUTINE ALLOC_AMBER_CONSTR()
         CALL DEALLOC_AMBER_CONST()
         ALLOCATE(AMBER_CONI(AMBER_NCONST))
         ALLOCATE(AMBER_CONJ(AMBER_NCONST))
         ALLOCATE(AMBER_CONDISTREF(AMBER_NCONST))
         ALLOCATE(AMBER_CONCUT(AMBER_NCONST))
      END SUBROUTINE ALLOC_AMBER_CONSTR

      ! dealloc amber constraints arrays
      SUBROUTINE DEALLOC_AMBER_CONSTR()
         IF (ALLOCATED(AMBER_CONI)) DEALLOCATE(AMBER_CONI)
         IF (ALLOCATED(AMBER_CONJ)) DEALLOCATE(AMBER_CONJ)
         IF (ALLOCATED(AMBER_CONDISTREF)) DEALLOCATE(AMBER_CONDISTREF)
         IF (ALLOCATED(AMBER_CONCUT)) DEALLOCATE(AMBER_CONCUT)
      END SUBROUTINE DEALLOC_AMBER_CONSTR

      ! parse topology
      SUBROUTINE PARSE_TOPOLOGY()
         USE HELPER_FNCTS, ONLY: READ_LINE
         USE QCIFILEHANDLER, ONLY: GETUNIT, FILE_LENGTH
         IMPLICIT NONE
         INTEGER :: TOPUNIT !topology unit
         INTEGER :: TOPLENGTH !topology file length
         INTEGER :: LINECOUNTER !current line being read
         LOGICAL :: YESNO
         CHARACTER(100) :: ENTRY !string for line from topology
         INTEGER , PARAMETER :: NWORDS=20 !maximum number of entries per line
         CHARACTER(25) :: ENTRIES(NWORDS)='' !array of entries per line
         INTEGER :: NBONDH, NBONDA, NANGH, NANGA
         INTEGER :: NLINES
         INTEGER, ALLOCATABLE :: INDICES(:)
         INTEGER, PARAMETER :: NENT_NBONDH = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NBONDA = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NANGH = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NANGA = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NPERENT_BOND = 3 ! integers per bond entry
         INTEGER, PARAMETER :: NPERENT_ANG = 4 ! integers per angle entry 


         !check topology exists
         INQUIRE(FILE=TOPFILENAME, EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(*,*) " qci_amber_constr> Cannot locate file ", TOPFILENAME
            STOP
         END IF
         !open topology
         TOPUNIT = GETUNIT()
         TOPLENGTH = FILE_LENGTH(TOPFILENAME)
         LINECOUNTER = 0

         OPEN(TOPUNIT,FILE=TOPFILENAME,STATUS='OLD')
         !parse topology file
         DO
            !check whether we are at the end of the file
            LINECOUNTER = LINECOUNTER + 1
            IF (LINECOUNTER.GT.TOPLENGTH) EXIT
            ! read line and parse it
            READ(TOPUNIT,'(A)') ENTRY
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES) 
            IF (ENTRIES(2).EQ.'POINTERS') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               READ(TOPUNIT,'(A)') ENTRY
               LINECOUNTER = LINECOUNTER + 1
               CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
               READ(ENTRIES(3),'(I8)') NBONDH
               READ(ENTRIES(4),'(I8)') NBONDA
               NBOND = NBONDH + NBONDA
               WRITE(MYUNIT,'(A,I8)') 'readtopology> Number of bonds:',NBOND
               ALLOCATE(BONDS(NBOND,2))
               READ(ENTRIES(3),'(I8)') NANGH
               READ(ENTRIES(4),'(I8)') NANGA
               NANGLE = NANGH + NANGA
               WRITE(MYUNIT,'(A,I8)') 'readtopology> Number of angles:',NANGLE
               ALLOCATE(ANGLES(NANGLE,2))
            END IF
            
            IF (ENTRIES(2).EQ. 'BONDS_INC_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               CALL GET_NLINES(NBONDH, NPERENT_BOND, NENT_NBONDH, NLINES)
               ALLOCATE(INDICES(NENT_NBONDH*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,NENT_NBONDH
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*NENT_BONDH+J2) 
                  END DO
               END DO
               DO J1=1,NBONDH
                  IDX = (J1-1)*NPERENT_BOND
                  BONDS(IDX,1) = INDICES(IDX+1)/3+1
                  BONDS(IDX,2) = INDICES(IDX+2)/3+1
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'BONDS_WITHOUT_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               CALL GET_NLINES(NBONDA, NPERENT_BOND, NENT_NBONDA, NLINES)
               ALLOCATE(INDICES(NENT_NBONDA*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,NENT_NBONDA
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*NENT_BONDA+J2) 
                  END DO
               END DO
               DO J1=1,NBONDA
                  IDX = (J1-1)*NPERENT_BOND
                  BONDS(NBONDH+IDX,1) = INDICES(IDX+1)/3+1
                  BONDS(NBONDH+IDX,2) = INDICES(IDX+2)/3+1
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'ANGLES_INC_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               CALL GET_NLINES(NANGH, NPERENT_ANG, NENT_NANGH, NLINES)
               ALLOCATE(INDICES(NENT_NANGH*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,NENT_NANGH
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*NENT_ANGH+J2) 
                  END DO
               END DO
               DO J1=1,NANGH
                  IDX = (J1-1)*NPERENT_ANG
                  ANGLES(IDX,1) = INDICES(IDX+1)/3+1 ! atom i 
                  ANGLES(IDX,2) = INDICES(IDX+3)/3+1 ! atom k
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'ANGLES_WITHOUT_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               CALL GET_NLINES(NANGA, NPERENT_ANG, NENT_NANGA, NLINES)
               ALLOCATE(INDICES(NENT_NANGA*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,NENT_NANGA
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*NENT_ANGA+J2) 
                  END DO
               END DO
               DO J1=1,NANGA
                  IDX = (J1-1)*NPERENT_ANG
                  ANGLES(NANGH+IDX,1) = INDICES(IDX+1)/3+1 ! atom i 
                  ANGLES(NANGH+IDX,2) = INDICES(IDX+3)/3+1 ! atom k
               END DO
               DEALLOCATE(INDICES)
            END IF

         END DO
         CLOSE(TOPUNIT)

      END SUBROUTINE PARSE_TOPOLOGY

      SUBROUTINE GET_NLINES(NENTRIES, PERENTRY, ENTRIESPERLINE, NLINES)
         INTEGER, INTENT(IN) :: NENTRIES
         INTEGER, INTENT(IN) :: PERENTRY
         INTEGER, INTENT(IN) :: ENTRIESPERLINE
         INTEGER, INTENT(OUT) :: NLINES
         INTEGER :: MODVAL

         MODVAL = MOD(NENTRIES*PERENTRY,ENTRIESPERLINE)
         NLINES = (NENTRIES*PERENTRY-MODVAL)/ENTRIESPERLINE
         IF (MODVAL.NE.0) NLINES = NLINES + 1
      END SUBROUTINE GET_NLINES

END MODULE AMBER_CONSTRAINTS
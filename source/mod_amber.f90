MODULE AMBER_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: NBACKBONE
   INTEGER :: NRES
   INTEGER, ALLOCATABLE :: BACKBONE(:)
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
   INTEGER, ALLOCATABLE :: RESFINAL
   CHARACTER(LEN=4), ALLOCATABLE :: AMBER_NAMES(:), RESNAMES(:)

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
         INTEGER :: NBIOCONSTR ! biological constraints - cis-trans and planarity 
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

      SUBROUTINE GET_BACKBONE(NATOMS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY
         INTEGER :: DUMMY_BB(NATOMS)
         CHARACTER(LEN=4) :: ATNAME
         INTEGER :: J1

         NDUMMY = 0
         DO J1=1,NATOMS
            ATNAME = ADJUSTL(TRIM(AMBER_NAMES(J1)))
            IF ((ATNAME.EQ.'C').OR.(ATNAME.EQ.'O').OR.(ATNAME.EQ.'N').OR.(ATNAME.EQ.'H').OR.(ATNAME.EQ.'CA') &
                .OR.(ATNAME.EQ.'HA').OR.(ATNAME.EQ.'CB').OR.(ATNAME.EQ.'HA2').OR.(ATNAME.EQ.'HA3').OR.(ATNAME.EQ.'P') &
                .OR.(ATNAME.EQ.'OP1').OR.(ATNAME.EQ.'OP2').OR.(ATNAME.EQ."O5'").OR.(ATNAME.EQ."C5'").OR.(ATNAME.EQ."H5'") &
                .OR.(ATNAME.EQ."H5''").OR.(ATNAME.EQ."C4'").OR.(ATNAME.EQ."H4'").OR.(ATNAME.EQ."O4'").OR.(ATNAME.EQ."C3'") &
                .OR.(ATNAME.EQ."H3'").OR.(ATNAME.EQ."C2'").OR.(ATNAME.EQ."O3'").OR.(ATNAME.EQ."HO3'").OR.(ATNAME.EQ."HO5'")) THEN
               NDUMMY = NDUMMY + 1
               DUMMY_BB(NDUMMY) = J1
            END IF
         END DO

         NBACKBONE = NDUMMY
         ALLOCATE(BACKBONE(NBACKBONE))
         BACKBONE(1:NBACKBONE) = DUMMY_BB(1:NBACKBONE)
      END SUBROUTINE GET_BACKBONE




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
         INTEGER :: J1, J2
         INTEGER, ALLOCATABLE :: INDICES(:)
         INTEGER, PARAMETER :: NENT_NBONDH = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NBONDA = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NANGH = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NENT_NANGA = 10 ! entries per line in the topology
         INTEGER, PARAMETER :: NPERENT_BOND = 3 ! integers per bond entry
         INTEGER, PARAMETER :: NPERENT_ANG = 4 ! integers per angle entry 
         CHARACTER(4) :: NAMES_CURR(20)

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
               !ignore format identifier after flag
               READ(TOPUNIT,*)                             
               LINECOUNTER = LINECOUNTER + 1
               ! the first line contains the information we ar eintereste din 
               READ(TOPUNIT,'(A)') ENTRY
               LINECOUNTER = LINECOUNTER + 1
               CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
               READ(ENTRIES(3),'(I8)') NBONDH
               READ(ENTRIES(4),'(I8)') NBONDA
               NBOND = NBONDH + NBONDA
               WRITE(MYUNIT,'(A,I8)') 'readtopology> Number of bonds:',NBOND
               ALLOCATE(BONDS(NBOND,2))
               READ(ENTRIES(5),'(I8)') NANGH
               READ(ENTRIES(6),'(I8)') NANGA
               NANGLE = NANGH + NANGA
               WRITE(MYUNIT,'(A,I8)') 'readtopology> Number of angles:',NANGLE
               ALLOCATE(ANGLES(NANGLE,2))
               READ(TOPUNIT,'(A)') ENTRY
               LINECOUNTER = LINECOUNTER + 1
               CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
               READ(ENTRIES(2),'(I8)') NRES
               ALLOCATE(RESNAMES(NRES),RES_START(NRES),RES_END(NRES),RESTYPE(NRES))
               RESTYPE(1:NRES) = 0
               ALLOCATE(AMBER_NAMES(NATOMS))
            END IF
            
            IF (ENTRIES(2).EQ. 'ATOM_NAME') THEN
               !ignore format identifier after flag
               READ(TOPUNIT,*)                             
               LINECOUNTER = LINECOUNTER + 1
               NLINES = NATOMS/20
               IF (NATOMS.GT.NLINES*20) NLINES=NLINES+1
               NDUMMY = 1
               DO J1=1,NLINES
                  READ(TOPUNIT,'(20 A4)') NAMES_CURR
                  LINECOUNTER = LINECOUNTER + 1
                  DO J2=1,20
                     AMBER_NAMES(NDUMMY)= NAMES_CURR(J4)
                     NDUMMY=NDUMMY+1
                     IF(NDUMMY.GT.NATOMS) EXIT
                  END DO
               END DO
            END IF

            IF (ENTRIES(2).EQ. 'RESIDUE_LABEL') THEN
               !ignore format identifier after flag
               READ(TOPUNIT,*)                             
               LINECOUNTER = LINECOUNTER + 1
               NLINES = NRES/20
               IF (NRES.GT.NLINES*20) NLINES=NLINES+1
               NDUMMY=1
               DO J1=1,NLINES
                  READ(TOPUNIT,'(20 A4)') NAMES_CURR
                  LINECOUNTER = LINECOUNTER + 1
                  DO J2=1,20
                     RESNAMES(NDUMMY) = NAMES_CURR(J4)
                     NDUMMY = NDUMMY+1
                     IF(NDUMMY.GT.NRES) EXIT
                  END DO
               END DO
            END IF

            IF (ENTRIES(2).EQ. 'RESIDUE_POINTER') THEN
               !ignore format identifier after flag
               READ(TOPUNIT,*)                             
               LINECOUNTER = LINECOUNTER + 1
               NLINES=NRES/10   
               IF(NRES.GT.NLINES*10) NLINES=NLINES+1 
               NDUMMY=1
               DO J1=1,LINES !go through all lines  
                  READ(TOPUNIT,'(A)') ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES) 
                  J2=1
                  DO WHILE(J2.LE.10)
                     READ(ENTRIES(J2),'(I8)') INTDUM
                     RES_START(NDUMMY) = INTDUM
                     IF (NDUMMY.GT.1) RES_END(NDUMMY-1) = INTDUM - 1 
                     J2=J2+1
                     NDUMMY = NDUMMY+1               
                     IF(NDUMMY.GT.NRES) EXIT
                  ENDDO
                  RES_END(NDUMMY-1) = NATOMS
               ENDDO
            ENDIF

            IF (ENTRIES(2).EQ. 'BONDS_INC_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               NLINES=(3*NBONDH)/10
               IF (NLINES*10.LT.NBONDH*3) NLINES=NLINES+1
               ALLOCATE(INDICES(10*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,10
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*10+J2) 
                  END DO
               END DO
               DO J1=1,NBONDH
                  IDX = (J1-1)*3
                  BONDS(IDX,1) = INDICES(IDX+1)/3+1
                  BONDS(IDX,2) = INDICES(IDX+2)/3+1
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'BONDS_WITHOUT_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               NLINES=(3*NBONDA)/10
               IF (NLINES*10.LT.NBONDA*3) NLINES=NLINES+1
               ALLOCATE(INDICES(10*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,NENT_NBONDA
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*10+J2) 
                  END DO
               END DO
               DO J1=1,NBONDA
                  IDX = (J1-1)*3
                  BONDS(NBONDH+IDX,1) = INDICES(IDX+1)/3+1
                  BONDS(NBONDH+IDX,2) = INDICES(IDX+2)/3+1
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'ANGLES_INC_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               NLINES=(3*NANGH)/10
               IF (NLINES*10.LT.NANGH*3) NLINES=NLINES+1
               ALLOCATE(INDICES(10*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,10
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*10+J2) 
                  END DO
               END DO
               DO J1=1,NANGH
                  IDX = (J1-1)*3
                  ANGLES(IDX,1) = INDICES(IDX+1)/3+1 ! atom i 
                  ANGLES(IDX,2) = INDICES(IDX+3)/3+1 ! atom k
               END DO
               DEALLOCATE(INDICES)
            END IF

            IF (ENTRIES(2).EQ. 'ANGLES_WITHOUT_HYDROGEN') THEN
               READ(TOPUNIT,*)                             !ignore format identifier after flag
               LINECOUNTER = LINECOUNTER + 1
               NLINES=(3*NANGH)/10
               IF (NLINES*10.LT.NANGA*3) NLINES=NLINES+1
               ALLOCATE(INDICES(10*NLINES))
               DO J1=1,NLINES
                  READ(TOPUNIT,*) ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
                  DO J2=1,10
                     READ(ENTRIES(J2), '(I8)') INDICES((J1-1)*10+J2) 
                  END DO
               END DO
               DO J1=1,NANGA
                  IDX = (J1-1)*3
                  ANGLES(NANGH+IDX,1) = INDICES(IDX+1)/3+1 ! atom i 
                  ANGLES(NANGH+IDX,2) = INDICES(IDX+3)/3+1 ! atom k
               END DO
               DEALLOCATE(INDICES)
            END IF

         END DO
         CLOSE(TOPUNIT)
         DO J1=1,NRES
            CALL CHECK_RES(J1,AAT,DNAT,RNAT)
            IF (AAT) THEN
               RESTYPE(J1) = 1
            ELSE IF (DNAT) THEN
               RESTYPE(J1) = 2
            ELSE IF (RNAT) THEN
               RESTYPE(J1) = 3
            ENDIF
         ENDDO
      END SUBROUTINE PARSE_TOPOLOGY

        ! check if residue is an amino acid - brute force ...
      SUBROUTINE CHECK_RES(RESID,AAT,DNAT,RNAT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         LOGICAL, INTENT(OUT) :: AAT, DNAT, RNAT
         CHARACTER(4) :: DNAME
         LOGICAL :: TERTEST
         
         DNAME = ADJUSTL(TRIM(RESNAMES(RESID)))
         AAT = .FALSE.
         RNAT = .FALSE.
         DNAT = .FALSE.
         TERTEST = .FALSE.
60       CONTINUE
         IF ((DNAME.EQ."ALA").OR.(DNAME.EQ."ARG").OR.(DNAME.EQ."ASH").OR.  &
            (DNAME.EQ."ASN").OR.(DNAME.EQ."ASP").OR.(DNAME.EQ."CYM").OR.  &
            (DNAME.EQ."CYS").OR.(DNAME.EQ."CYX").OR.(DNAME.EQ."GLH").OR.  &
            (DNAME.EQ."GLN").OR.(DNAME.EQ."GLU").OR.(DNAME.EQ."GLY").OR.  &
            (DNAME.EQ."HID").OR.(DNAME.EQ."HIE").OR.(DNAME.EQ."HIP").OR.  &
            (DNAME.EQ."HYP").OR.(DNAME.EQ."ILE").OR.(DNAME.EQ."LEU").OR.  &
            (DNAME.EQ."LYN").OR.(DNAME.EQ."LYS").OR.(DNAME.EQ."MET").OR.  &
            (DNAME.EQ."PHE").OR.(DNAME.EQ."PRO").OR.(DNAME.EQ."SER").OR.  &
            (DNAME.EQ."THR").OR.(DNAME.EQ."TRP").OR.(DNAME.EQ."TYR").OR.  &
            (DNAME.EQ."VAL")) THEN
            AAT = .TRUE.
         ELSE IF ((DNAME.EQ."A").OR.(DNAME.EQ."A3").OR.(DNAME.EQ."A5").OR.(DNAME.EQ."AN").OR.  &
                  (DNAME.EQ."C").OR.(DNAME.EQ."C3").OR.(DNAME.EQ."C5").OR.(DNAME.EQ."CN").OR.  &
                  (DNAME.EQ."G").OR.(DNAME.EQ."G3").OR.(DNAME.EQ."G5").OR.(DNAME.EQ."GN").OR.  &
                  (DNAME.EQ."U").OR.(DNAME.EQ."U3").OR.(DNAME.EQ."U5").OR.(DNAME.EQ."UN")) THEN
            RNAT = .TRUE.
         ELSE IF ((DNAME.EQ."DA").OR.(DNAME.EQ."DA3").OR.(DNAME.EQ."DA5").OR.(DNAME.EQ."DAN").OR.  &
                  (DNAME.EQ."DC").OR.(DNAME.EQ."DC3").OR.(DNAME.EQ."DC5").OR.(DNAME.EQ."DCN").OR.  &
                  (DNAME.EQ."DG").OR.(DNAME.EQ."DG3").OR.(DNAME.EQ."DG5").OR.(DNAME.EQ."DGN").OR.  &
                  (DNAME.EQ."DT").OR.(DNAME.EQ."DT3").OR.(DNAME.EQ."DT5").OR.(DNAME.EQ."DTN")) THEN
            DNAT = .TRUE.
         ENDIF
         !make sure it is not a terminal residue (we technically should never encounter one)
         IF (((.NOT.AAT).AND.(.NOT.TERTEST)).AND.((DNAME(1:1).EQ."N").OR.(DNAME(1:1).EQ."C"))) THEN
            TERTEST=.TRUE.
            DNAME = DNAME(2:4)
            GOTO 60
         ENDIF
      END SUBROUTINE CHECK_RES

      !get atom id from name for given residue
      !returns 0 if atom does not exist in residue
      SUBROUTINE GET_ATOMID(ATNAME,RESID,ATOMID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         CHARACTER(*), INTENT(IN) :: ATNAME
         INTEGER, INTENT(OUT) :: ATOMID
         INTEGER :: FIRST, LAST, J1
         
         ATOMID = 0
         FIRST=RES_START(RESID)
         LAST=RES_END(RESID)
         DO J1=FIRST,LAST
            IF (ADJUSTL(TRIM(ATNAMES(J1))).EQ.ADJUSTL(TRIM(ATNAME))) THEN
               ATOMID = J1
               EXIT
            ENDIF
         ENDDO
         RETURN
      END SUBROUTINE GET_ATOMID

END MODULE AMBER_CONSTRAINTS
MODULE AMBER_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
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
   INTEGER, ALLOCATABLE :: BIOCONSTR(:,:)   
   INTEGER, ALLOCATABLE :: ELEMENT(:)
   INTEGER :: NBOND = 0
   INTEGER :: NANGLE = 0
   INTEGER :: NBIOCONSTR = 0 ! biological constraints - cis-trans and planarity
   INTEGER, ALLOCATABLE :: RES_START(:), RES_END(:)
   CHARACTER(LEN=4), ALLOCATABLE :: AMBER_NAMES(:), RESNAMES(:), RESTYPE(:)

   CONTAINS
      ! from AMBER
      SUBROUTINE AMBER_QCI_CONSTRAINTS(NATOMS)
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE QCIFILEHANDLER, ONLY: GETUNIT, FILE_LENGTH
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY ! counter for constraints
         REAL(KIND = REAL64) :: DF, DS ! distance in final and initial structure
         LOGICAL :: YESNO ! do we have a contact file
         INTEGER :: CONUNIT ! unit for opening file
         INTEGER :: NADDCONSTR ! number of additional constraints
         INTEGER :: J1
         INTEGER :: IDX1, IDX2 !indices for atoms

         ! parse topology
         CALL PARSE_TOPOLOGY()
         ! create atomstores list
         CALL CREATE_ATOMS2RES()

         ! check for additional constraints in file
         INQUIRE(FILE=AMBERCONSTRFILE, EXIST=YESNO)
         IF (YESNO) THEN
            NADDCONSTR = FILE_LENGTH(AMBERCONSTRFILE)
         ELSE
            NADDCONSTR = 0
         END IF
         CALL GET_BIOCONSTR()

         ! allocate arrays
         AMBER_NCONST = NBOND + NANGLE + NBIOCONSTR + NADDCONSTR
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
         DO J1=1,NANGLE
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
         ! add bioconstraints
         DO J1=1,NBIOCONSTR
            NDUMMY = NDUMMY + 1
            IDX1 = BIOCONSTR(J1,1)
            IDX2 = BIOCONSTR(J1,2)
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
         WRITE(*,*) " amber_constraints> Identified ", AMBER_NCONST, " constraints"
         WRITE(*,*) "                    Bonds: ", NBOND, ", angles: ", NANGLE, ", additional constraints: ", NADDCONSTR

      END SUBROUTINE AMBER_QCI_CONSTRAINTS

      SUBROUTINE CREATE_ATOMS2RES()
         USE QCIKEYS, ONLY: ATOMS2RES, NATOMS
         IMPLICIT NONE
         INTEGER :: I, START, END
         
         IF (.NOT.ALLOCATED(ATOMS2RES)) DEALLOCATE(ATOMS2RES)
         ALLOCATE(ATOMS2RES(NATOMS))

         DO I=1,NRES
            START = RES_START(I)
            END = RES_END(I)
            ATOMS2RES(START:END) = I
         END DO

      END SUBROUTINE CREATE_ATOMS2RES


      SUBROUTINE AMBER_QCI_DEALLOCATE()
         IF (ALLOCATED(BACKBONE)) DEALLOCATE(BACKBONE)
         IF (ALLOCATED(BIOCONSTR)) DEALLOCATE(BIOCONSTR)
         IF (ALLOCATED(BONDS)) DEALLOCATE(BONDS)
         IF (ALLOCATED(ANGLES)) DEALLOCATE(ANGLES)
         IF (ALLOCATED(RESNAMES)) DEALLOCATE(RESNAMES)
         IF (ALLOCATED(RES_START)) DEALLOCATE(RES_START)
         IF (ALLOCATED(RES_END)) DEALLOCATE(RES_END)
         IF (ALLOCATED(AMBER_NAMES)) DEALLOCATE(AMBER_NAMES)
         IF (ALLOCATED(ELEMENT)) DEALLOCATE(ELEMENT)
         IF (ALLOCATED(RESTYPE)) DEALLOCATE(RESTYPE)
      END SUBROUTINE AMBER_QCI_DEALLOCATE

      SUBROUTINE GET_BACKBONE(NATOMS)
         USE QCIKEYS, ONLY: NBACKBONE, ISBBATOM
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER :: NDUMMY
         INTEGER :: DUMMY_BB(NATOMS)
         CHARACTER(LEN=4) :: ATNAME
         INTEGER :: J1

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
         DO J1=1,NBACKBONE
            ISBBATOM(BACKBONE(J1)) = .TRUE.
         END DO
      END SUBROUTINE GET_BACKBONE

      SUBROUTINE GET_BIOCONSTR()
         IMPLICIT NONE
        
         INTEGER :: NDUMMY, J1, OPOS1, CPOS1, CPOS2, HPOS2, ATOMID
         LOGICAL :: AAT, RNAT, DNAT, CAPT, AALIST(NRES), ISTER(NRES), ISCAP(NRES)
         INTEGER, ALLOCATABLE :: DUMMYC(:,:)


         ! the current largest number of additional constraints per residue is 8 and we have two additional constraints for peptide bonds,
         ! and seven constraints for chiral centres and sugar and base in nucleic acids - the worst case seems like 16 constraints per residue
         ! we use 20 for this array size, and should have plenty of space
         ALLOCATE(BIOCONSTR(20*NRES,2))
         BIOCONSTR(:,:) = -1
         NDUMMY = 0
         AALIST(1:NRES) = .FALSE.
         ISTER(1:NRES) = .FALSE.
         ISCAP(1:NRES) = .FALSE.

         ! this is the first round to add constraints for planarity, the sugars in NAs etc.
         DO J1=1,NRES
            !find what kind of residue we have
            CALL CHECK_RES(J1, AAT, DNAT, RNAT, CAPT)
            IF (CAPT) THEN
               ISCAP(J1) = .TRUE.
               AALIST(J1) = .TRUE.
            ! for amino acids check whether this is a temrinal residue and then find specific constraints
            ELSE IF (AAT) THEN
               CALL GET_ATOMID("OXT",J1,ATOMID)
               IF (ATOMID.NE.0) THEN
                  ISTER(J1) = .TRUE.
               ELSE
                  CALL GET_ATOMID("H1",J1,ATOMID)
                  IF (ATOMID.NE.0) THEN
                     ISTER(J1) = .TRUE.
                  END IF
               END IF
               AALIST(J1) = .TRUE.
               CALL GET_AA_CONSTR(J1)
            ELSE IF (RNAT.OR.DNAT) THEN
               CALL GET_NA_CONSTR(J1)
            END IF
         END DO

         ! peptide bonds
         DO J1=1,NRES-1
            IF (AALIST(J1)) THEN
               ! chekc we are not looking at capping groups
               IF ((.NOT.ISCAP(J1)).AND.(.NOT.ISCAP(J1+1))) THEN
                  IF (.NOT.ISTER(J1)) THEN
                     CALL GET_ATOMID("O",J1,OPOS1)
                     CALL GET_ATOMID("CA",J1,CPOS1)
                     CALL GET_ATOMID("H",J1+1,HPOS2)
                     CALL GET_ATOMID("CA",J1+1,CPOS2)
                     ! O(i) - H(i+1)
                     IF ((OPOS1.GT.0).AND.(HPOS2.GT.0)) THEN
                        NBIOCONSTR = NBIOCONSTR + 1
                        BIOCONSTR(NBIOCONSTR,1) = OPOS1
                        BIOCONSTR(NBIOCONSTR,2) = HPOS2
                     END IF
                     ! CA(i) - CA(i+1)
                     IF ((CPOS1.GT.0).AND.(CPOS2.GT.0)) THEN
                        NBIOCONSTR = NBIOCONSTR + 1
                        BIOCONSTR(NBIOCONSTR,1) = CPOS1
                        BIOCONSTR(NBIOCONSTR,2) = CPOS2
                     END IF
                     ! CA(i) - H(i+1)
                     IF ((CPOS1.GT.0).AND.(HPOS2.GT.0)) THEN
                        NBIOCONSTR = NBIOCONSTR + 1
                        BIOCONSTR(NBIOCONSTR,1) = CPOS1
                        BIOCONSTR(NBIOCONSTR,2) = HPOS2
                     END IF
                  END IF
               END IF
            END IF
         END DO

         ALLOCATE(DUMMYC(NBIOCONSTR,2))
         DUMMYC(1:NBIOCONSTR,1:2) = BIOCONSTR(1:NBIOCONSTR,1:2)
         DEALLOCATE(BIOCONSTR)
         ALLOCATE(BIOCONSTR(NBIOCONSTR,2))
         BIOCONSTR(1:NBIOCONSTR,1:2) = DUMMYC(1:NBIOCONSTR,1:2)

      END SUBROUTINE GET_BIOCONSTR

      SUBROUTINE GET_AA_CONSTR(RESID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID

         CHARACTER(LEN=4) :: RES

         RES = RESNAMES(RESID)
         IF ((RES.EQ.'ARG').OR.(RES.EQ.'CARG').OR.(RES.EQ.'NARG')) THEN
            CALL ADD_CONSTRAINT("HH12 ","HH22",RESID)
            CALL ADD_CONSTRAINT("HH12 ","HH21",RESID)
            CALL ADD_CONSTRAINT("HH11 ","HH22",RESID)
            CALL ADD_CONSTRAINT("HH11 ","HH21",RESID)
            CALL ADD_CONSTRAINT("HH11 ","NE  ",RESID)
            CALL ADD_CONSTRAINT("HH12 ","NE  ",RESID)
            CALL ADD_CONSTRAINT("HH21 ","NE  ",RESID)
            CALL ADD_CONSTRAINT("HH22 ","NE  ",RESID)
         ELSE IF ((RES.EQ.'HID').OR.(RES.EQ.'HIE').OR.(RES.EQ.'HIS').OR.(RES.EQ.'HIP').OR. &
                  (RES.EQ.'CHID').OR.(RES.EQ.'CHIE').OR.(RES.EQ.'CHIS').OR.(RES.EQ.'CHIP').OR. &
                  (RES.EQ.'NHID').OR.(RES.EQ.'NHIE').OR.(RES.EQ.'NHIS').OR.(RES.EQ.'NHIP')) THEN
            CALL ADD_CONSTRAINT("HD1 ","NE2 ",RESID)
            CALL ADD_CONSTRAINT("HE1 ","CD2 ",RESID)
            CALL ADD_CONSTRAINT("HE2 ","CG  ",RESID)
            CALL ADD_CONSTRAINT("HD2 ","ND1 ",RESID)
         ELSE IF ((RES.EQ.'PHE').OR.(RES.EQ.'CPHE').OR.(RES.EQ.'NPHE')) THEN
            CALL ADD_CONSTRAINT("HD2 ","CE1 ",RESID)
            CALL ADD_CONSTRAINT("HE2 ","CD1 ",RESID)
            CALL ADD_CONSTRAINT("HZ  ","CG  ",RESID)
            CALL ADD_CONSTRAINT("HE1 ","CD2 ",RESID)
            CALL ADD_CONSTRAINT("HD1 ","CE2 ",RESID)
            CALL ADD_CONSTRAINT("CB  ","CZ  ",RESID)
         ELSE IF ((RES.EQ.'TRP').OR.(RES.EQ.'CTRP').OR.(RES.EQ.'NTRP')) THEN
            CALL ADD_CONSTRAINT("HD1 ","CD2 ",RESID)
            CALL ADD_CONSTRAINT("HE1 ","CG  ",RESID)
            CALL ADD_CONSTRAINT("HZ2 ","CE3 ",RESID)
            CALL ADD_CONSTRAINT("HH2 ","CD2 ",RESID)
            CALL ADD_CONSTRAINT("HZ3 ","CE2 ",RESID)
            CALL ADD_CONSTRAINT("HE3 ","CZ2 ",RESID)
         ELSE IF ((RES.EQ.'TYR').OR.(RES.EQ.'CTYR').OR.(RES.EQ.'NTYR')) THEN
            CALL ADD_CONSTRAINT("HD1 ","CE2 ",RESID)
            CALL ADD_CONSTRAINT("HE1 ","CD2 ",RESID)
            CALL ADD_CONSTRAINT("OH  ","CG  ",RESID)
            CALL ADD_CONSTRAINT("HE2 ","CD1 ",RESID)
            CALL ADD_CONSTRAINT("HD2 ","CE1 ",RESID)
            CALL ADD_CONSTRAINT("CB  ","CZ  ",RESID)
         END IF 
      END SUBROUTINE GET_AA_CONSTR

      SUBROUTINE GET_NA_CONSTR(RESID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID

         CHARACTER(LEN=4) :: NAME

         NAME = RESNAMES(RESID)
         IF ((NAME.EQ.'A').OR.(NAME.EQ.'A3').OR.(NAME.EQ.'A5').OR.  &
             (NAME.EQ.'AN').OR.(NAME.EQ.'DA').OR.(NAME.EQ.'DA3').OR.  &
             (NAME.EQ.'DA5').OR.(NAME.EQ.'DAN')) THEN
            
            CALL ADD_CONSTRAINT("N1  ","N9  ",RESID)        
            CALL ADD_CONSTRAINT("C2  ","N7  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","C6  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","N6  ",RESID) 
            CALL ADD_CONSTRAINT("C4  ","N1  ",RESID) 
            CALL ADD_CONSTRAINT("C5  ","H2  ",RESID) 
            CALL ADD_CONSTRAINT("C8  ","C4  ",RESID) 
            CALL ADD_CONSTRAINT("C8  ","C5  ",RESID) 

         ELSE IF ((NAME.EQ.'C').OR.(NAME.EQ.'C3').OR.(NAME.EQ.'C5').OR.  &
                  (NAME.EQ.'CN').OR.(NAME.EQ.'DC').OR.(NAME.EQ.'DC3').OR.  &
                  (NAME.EQ.'DC5').OR.(NAME.EQ.'DCN')) THEN

            CALL ADD_CONSTRAINT("N1  ","N4  ",RESID) 
            CALL ADD_CONSTRAINT("C2  ","H5  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","H6  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","C6  ",RESID) 
            CALL ADD_CONSTRAINT("C4  ","N1  ",RESID) 
            CALL ADD_CONSTRAINT("C5  ","O2  ",RESID) 
        
         ELSE IF ((NAME.EQ.'G').OR.(NAME.EQ.'G3').OR.(NAME.EQ.'G5').OR.  &
                  (NAME.EQ.'GN').OR.(NAME.EQ.'DG').OR.(NAME.EQ.'DG3').OR.  &
                  (NAME.EQ.'DG5').OR.(NAME.EQ.'DGN')) THEN
        
            CALL ADD_CONSTRAINT("N1  ","N9  ",RESID)        
            CALL ADD_CONSTRAINT("C2  ","N7  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","C6  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","O6  ",RESID) 
            CALL ADD_CONSTRAINT("C4  ","H1  ",RESID) 
            CALL ADD_CONSTRAINT("C5  ","N2  ",RESID) 
            CALL ADD_CONSTRAINT("C8  ","C4  ",RESID) 
            CALL ADD_CONSTRAINT("C8  ","C5  ",RESID) 

         ELSE IF ((NAME.EQ.'U').OR.(NAME.EQ.'U3').OR.(NAME.EQ.'U5').OR.  &
                  (NAME.EQ.'UN')) THEN

            CALL ADD_CONSTRAINT("N1  ","O4  ",RESID) 
            CALL ADD_CONSTRAINT("C2  ","H5  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","H6  ",RESID) 
            CALL ADD_CONSTRAINT("H3  ","C6  ",RESID) 
            CALL ADD_CONSTRAINT("C4  ","N1  ",RESID) 
            CALL ADD_CONSTRAINT("C5  ","O2  ",RESID)
                  
         ELSE IF ((NAME.EQ.'DT').OR.(NAME.EQ.'DT3').OR.  &
                  (NAME.EQ.'DT5').OR.(NAME.EQ.'DTN')) THEN

            CALL ADD_CONSTRAINT("N1  ","O4  ",RESID) 
            CALL ADD_CONSTRAINT("C2  ","C7  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","H6  ",RESID) 
            CALL ADD_CONSTRAINT("N3  ","C6  ",RESID) 
            CALL ADD_CONSTRAINT("C4  ","N1  ",RESID) 
            CALL ADD_CONSTRAINT("C5  ","O2  ",RESID) 

         ENDIF

         !constrain the sugar to keep its shape
         CALL ADD_CONSTRAINT("O4' ","C2' ",RESID) 
         CALL ADD_CONSTRAINT("C1' ","C3' ",RESID)         

         !constrain chiral centres
         CALL ADD_CONSTRAINT("O2' ","O3' ",RESID)
         CALL ADD_CONSTRAINT("H2''","O3' ",RESID)
         CALL ADD_CONSTRAINT("H3' ","C5' ",RESID)
         CALL ADD_CONSTRAINT("H3' ","H2' ",RESID)
         CALL ADD_CONSTRAINT("O3''","H4' ",RESID)

      END SUBROUTINE GET_NA_CONSTR

      SUBROUTINE ADD_CONSTRAINT(ATNAME1,ATNAME2,RESID)
         CHARACTER(LEN=4), INTENT(IN) :: ATNAME1, ATNAME2
         INTEGER, INTENT(IN) :: RESID
         INTEGER :: ID1, ID2 

         CALL GET_ATOMID(ATNAME1,RESID,ID1)
         CALL GET_ATOMID(ATNAME2,RESID,ID2)
         IF ((ID1.GT.0).AND.(ID2.GT.0)) THEN
            NBIOCONSTR = NBIOCONSTR + 1
            BIOCONSTR(NBIOCONSTR,1) = ID1
            BIOCONSTR(NBIOCONSTR,2) = ID2
         END IF
      END SUBROUTINE ADD_CONSTRAINT

      ! allocate amber constraints array
      SUBROUTINE ALLOC_AMBER_CONSTR()
         CALL DEALLOC_AMBER_CONSTR()
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
         USE QCIKEYS, ONLY: NATOMS
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
         INTEGER :: J1, J2, IDX, INTDUM, NDUMMY
         INTEGER, ALLOCATABLE :: INDICES(:)
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
               WRITE(*,'(A,I8)') ' readtopology> Number of bonds:',NBOND
               ALLOCATE(BONDS(NBOND,2))
               READ(ENTRIES(5),'(I8)') NANGH
               READ(ENTRIES(6),'(I8)') NANGA
               NANGLE = NANGH + NANGA
               WRITE(*,'(A,I8)') ' readtopology> Number of angles:',NANGLE
               ALLOCATE(ANGLES(NANGLE,2))
               READ(TOPUNIT,'(A)') ENTRY
               LINECOUNTER = LINECOUNTER + 1
               CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
               READ(ENTRIES(2),'(I8)') NRES
               ALLOCATE(RESNAMES(NRES),RES_START(NRES),RES_END(NRES))
               ALLOCATE(AMBER_NAMES(NATOMS),ELEMENT(NATOMS),RESTYPE(NRES))
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
                     AMBER_NAMES(NDUMMY)= NAMES_CURR(J2)
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
                     RESNAMES(NDUMMY) = NAMES_CURR(J2)
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
               DO J1=1,NLINES !go through all lines  
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

            IF (ENTRIES(2).EQ. 'ATOMIC_NUMBER') THEN  
               !ignore format identifier after flag
               READ(TOPUNIT,*)                             
               LINECOUNTER = LINECOUNTER + 1
               NLINES=NATOMS/10   
               IF(NATOMS.GT.NLINES*10) NLINES=NLINES+1 
               NDUMMY=1
               DO J1=1,NLINES !go through all lines  
                  READ(TOPUNIT,'(A)') ENTRY
                  LINECOUNTER = LINECOUNTER + 1
                  CALL READ_LINE(ENTRY,NWORDS,ENTRIES) 
                  J2=1
                  DO WHILE(J2.LE.10)
                     READ(ENTRIES(J2),'(I8)') INTDUM
                     ELEMENT(NDUMMY) = INTDUM
                     NDUMMY = NDUMMY+1               
                     IF(NDUMMY.GT.NATOMS) EXIT
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
                  DO J2=1,10
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
      END SUBROUTINE PARSE_TOPOLOGY

        ! check if residue is an amino acid - brute force ...
      SUBROUTINE CHECK_RES(RESID,AAT,DNAT,RNAT,CAPT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         LOGICAL, INTENT(OUT) :: AAT, DNAT, RNAT, CAPT
         CHARACTER(4) :: DNAME
         LOGICAL :: TERTEST
         
         DNAME = ADJUSTL(TRIM(RESNAMES(RESID)))
         AAT = .FALSE.
         RNAT = .FALSE.
         DNAT = .FALSE.
         CAPT = .FALSE.
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
            RESTYPE(RESID) = "AA"
         ELSE IF ((DNAME.EQ."A").OR.(DNAME.EQ."A3").OR.(DNAME.EQ."A5").OR.(DNAME.EQ."AN").OR.  &
                  (DNAME.EQ."C").OR.(DNAME.EQ."C3").OR.(DNAME.EQ."C5").OR.(DNAME.EQ."CN").OR.  &
                  (DNAME.EQ."G").OR.(DNAME.EQ."G3").OR.(DNAME.EQ."G5").OR.(DNAME.EQ."GN").OR.  &
                  (DNAME.EQ."U").OR.(DNAME.EQ."U3").OR.(DNAME.EQ."U5").OR.(DNAME.EQ."UN")) THEN
            RNAT = .TRUE.
            RESTYPE(RESID) = "RNA"           
         ELSE IF ((DNAME.EQ."DA").OR.(DNAME.EQ."DA3").OR.(DNAME.EQ."DA5").OR.(DNAME.EQ."DAN").OR.  &
                  (DNAME.EQ."DC").OR.(DNAME.EQ."DC3").OR.(DNAME.EQ."DC5").OR.(DNAME.EQ."DCN").OR.  &
                  (DNAME.EQ."DG").OR.(DNAME.EQ."DG3").OR.(DNAME.EQ."DG5").OR.(DNAME.EQ."DGN").OR.  &
                  (DNAME.EQ."DT").OR.(DNAME.EQ."DT3").OR.(DNAME.EQ."DT5").OR.(DNAME.EQ."DTN")) THEN
            DNAT = .TRUE.
            RESTYPE(RESID) = "DNA"            
         ELSE IF ((DNAME.EQ."ACE").OR.(DNAME.EQ."NHE").OR.(DNAME.EQ."NME")) THEN
            CAPT = .TRUE. 
         ENDIF
         !make sure it is not a terminal residue (we technically should never encounter one)
         IF (((.NOT.AAT).AND.(.NOT.CAPT).AND.(.NOT.TERTEST)).AND.((DNAME(1:1).EQ."N").OR.(DNAME(1:1).EQ."C"))) THEN
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
            IF (ADJUSTL(TRIM(AMBER_NAMES(J1))).EQ.ADJUSTL(TRIM(ATNAME))) THEN
               ATOMID = J1
               EXIT
            ENDIF
         ENDDO
         RETURN
      END SUBROUTINE GET_ATOMID

END MODULE AMBER_CONSTRAINTS
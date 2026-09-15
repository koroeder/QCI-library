!> Integration tests for MODULE CHIRALITY.
!!
!! Build (after applying the fixes):
!!   gfortran -cpp -DSWAP_HAS_SUCCESS -std=f2008 -fcheck=all -O0 \
!!            -o test_chirality \
!!            test_stubs.f90 comparelist.f90 chirality.f90 \
!!            test_fixtures.f90 test_chirality.f90
!!
!! Omit -DSWAP_HAS_SUCCESS to build against a SWITCH_SMALL_GROUPS that still
!! has the original three-argument signature; the geometry tests still run.
!!
!! Each test names the specific defect it covers, so a failure points straight
!! at the patch that regressed.
PROGRAM TEST_CHIRALITY_MAIN
   USE QCIPREC
   USE QCIKEYS,            ONLY: NATOMS, NIMAGES, DEBUG, ATOMS2RES
   USE AMBER_CONSTRAINTS,  ONLY: NBOND, BONDS, ELEMENT, AMBER_NAMES, RESTYPE
   USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
   USE CHIRALITY
   USE TEST_FIXTURES
   IMPLICIT NONE

   INTEGER :: NPASS, NFAIL

   NPASS = 0
   NFAIL = 0
   DEBUG = .FALSE.

   WRITE(*,'(A)') "=============================================="
   WRITE(*,'(A)') " CHIRALITY integration tests"
   WRITE(*,'(A)') "=============================================="

   CALL T_DETECT_ALA(NPASS,NFAIL)
   CALL T_DETECT_THR(NPASS,NFAIL)
   CALL T_DETECT_ILE(NPASS,NFAIL)
   CALL T_DETECT_SYNTH_NC3(NPASS,NFAIL)
   CALL T_DETECT_RNA(NPASS,NFAIL)
   CALL T_DETECT_DNA(NPASS,NFAIL)
   CALL T_RESTYPE_GUARD(NPASS,NFAIL)
   CALL T_GHOST_ATOMS(NPASS,NFAIL)
   CALL T_SWAP_GROUPS_ALA(NPASS,NFAIL)
   CALL T_SWAP_GROUPS_THR(NPASS,NFAIL)
   CALL T_PRIORITY_IS_PERMUTATION(NPASS,NFAIL)
   CALL T_ASSIGNMENT_SR(NPASS,NFAIL)
   CALL T_SWITCH_GEOMETRY(NPASS,NFAIL)
   CALL T_CHIRALITY_CHECK_BAND(NPASS,NFAIL)

   WRITE(*,'(A)') "----------------------------------------------"
   WRITE(*,'(A,I5,A,I5)') " passed: ", NPASS, "   failed: ", NFAIL
   WRITE(*,'(A)') "=============================================="
   IF (NFAIL.GT.0) THEN
      WRITE(*,'(A)') " RESULT: FAIL"
      STOP 1
   ELSE
      WRITE(*,'(A)') " RESULT: PASS"
   END IF

CONTAINS

   SUBROUTINE REPORT(NAME,OK,NPASS,NFAIL,DETAIL)
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      LOGICAL, INTENT(IN) :: OK
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: DETAIL
      IF (OK) THEN
         NPASS = NPASS + 1
         WRITE(*,'(A,A)') "  ok   : ", TRIM(NAME)
      ELSE
         NFAIL = NFAIL + 1
         WRITE(*,'(A,A)') "  FAIL : ", TRIM(NAME)
         IF (PRESENT(DETAIL)) WRITE(*,'(A,A)') "         ", TRIM(DETAIL)
      END IF
   END SUBROUTINE REPORT

   !> .TRUE. if atom NAME is one of the detected chiral centres.
   LOGICAL FUNCTION IS_CENTRE(NAME)
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      INTEGER :: J1, ID
      ID = ATID(NAME)
      IS_CENTRE = .FALSE.
      DO J1=1,NCHIRAL
         IF (CHIR_INFO(J1,1).EQ.ID) IS_CENTRE = .TRUE.
      END DO
   END FUNCTION IS_CENTRE

   INTEGER FUNCTION CENTRE_INDEX(NAME)
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      INTEGER :: J1, ID
      ID = ATID(NAME)
      CENTRE_INDEX = -1
      DO J1=1,NCHIRAL
         IF (CHIR_INFO(J1,1).EQ.ID) CENTRE_INDEX = J1
      END DO
   END FUNCTION CENTRE_INDEX

   SUBROUTINE RESET()
      CALL DEALLOC_CHIR_INTERNALS()
      NCHIRAL = -1
   END SUBROUTINE RESET

   !====================================================================
   ! detection
   !====================================================================

   SUBROUTINE T_DETECT_ALA(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- ALA: single chiral centre"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED()
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("ALA has exactly 1 chiral centre", NCHIRAL.EQ.1, NPASS,NFAIL,MSG)
      CALL REPORT("ALA centre is CA", IS_CENTRE("CA"), NPASS,NFAIL)
      ! CB carries three hydrogens -> DISCOUNT_H must reject it
      CALL REPORT("ALA methyl CB not a centre", .NOT.IS_CENTRE("CB"), NPASS,NFAIL)
   END SUBROUTINE T_DETECT_ALA

   SUBROUTINE T_DETECT_THR(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- THR: two chiral centres (EQ23 pair branch on CB)"
      CALL RESET(); CALL FIXTURE_THR_CAPPED()
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("THR has exactly 2 chiral centres", NCHIRAL.EQ.2, NPASS,NFAIL,MSG)
      CALL REPORT("THR CA is a centre", IS_CENTRE("CA"), NPASS,NFAIL)
      CALL REPORT("THR CB is a centre", IS_CENTRE("CB"), NPASS,NFAIL)
   END SUBROUTINE T_DETECT_THR

   !> ILE CB has substituent elements C,C,C,H -> TRIPLE with EQ12.AND.EQ23.
   !! Covers the first half of the FIND_CHIRAL_CENTRES TRIPLE branch.
   SUBROUTINE T_DETECT_ILE(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- ILE: TRIPLE branch, EQ12.AND.EQ23 (C,C,C,H)"
      CALL RESET(); CALL FIXTURE_ILE_CAPPED()
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("ILE has exactly 2 chiral centres", NCHIRAL.EQ.2, NPASS,NFAIL,MSG)
      CALL REPORT("ILE CA is a centre", IS_CENTRE("CA"), NPASS,NFAIL)
      CALL REPORT("ILE CB is a centre", IS_CENTRE("CB"), NPASS,NFAIL)
   END SUBROUTINE T_DETECT_ILE

   !> Synthetic N,C,C,C centre -> TRIPLE with EQ23.AND.EQ34.
   !! This is the branch that compared LIST1..LIST3 instead of LIST2..LIST4
   !! and never compared substituent 4 at all.  With the original code the
   !! three carbon groups are ranked without ever being compared against each
   !! other correctly, so this test is the direct regression guard.
   SUBROUTINE T_DETECT_SYNTH_NC3(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- synthetic amine: TRIPLE branch, EQ23.AND.EQ34 (N,C,C,C)"
      CALL RESET(); CALL FIXTURE_SYNTH_NC3()
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("synthetic CX detected as chiral", IS_CENTRE("CX"), NPASS,NFAIL,MSG)
      CALL REPORT("synthetic: no spurious extra centres", NCHIRAL.EQ.1, NPASS,NFAIL,MSG)
   END SUBROUTINE T_DETECT_SYNTH_NC3

   !> RNA: C1', C2', C3', C4' are all chiral.  C3' is only resolvable via the
   !! hardcoded AMBER_NAMES=="C3'" branch, because the 2'-OH makes the
   !! first-shell lists of C2' and C4' identical.
   SUBROUTINE T_DETECT_RNA(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- RNA adenosine: four ribose centres"
      CALL RESET(); CALL FIXTURE_RNA_ADENOSINE(.TRUE.,"RNA")
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("RNA has exactly 4 chiral centres", NCHIRAL.EQ.4, NPASS,NFAIL,MSG)
      CALL REPORT("RNA C1' is a centre", IS_CENTRE("C1'"), NPASS,NFAIL)
      CALL REPORT("RNA C2' is a centre", IS_CENTRE("C2'"), NPASS,NFAIL)
      CALL REPORT("RNA C3' is a centre", IS_CENTRE("C3'"), NPASS,NFAIL)
      CALL REPORT("RNA C4' is a centre", IS_CENTRE("C4'"), NPASS,NFAIL)
   END SUBROUTINE T_DETECT_RNA

   !> Deoxyribose: C2' carries two hydrogens and must be rejected, and C3' must
   !! resolve WITHOUT the special case.
   SUBROUTINE T_DETECT_DNA(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      CHARACTER(LEN=80) :: MSG
      WRITE(*,'(A)') " -- DNA deoxyadenosine: three ribose centres"
      CALL RESET(); CALL FIXTURE_RNA_ADENOSINE(.FALSE.,"DNA")
      CALL FIND_CHIRAL_CENTRES()
      WRITE(MSG,'(A,I3)') "NCHIRAL=",NCHIRAL
      CALL REPORT("DNA has exactly 3 chiral centres", NCHIRAL.EQ.3, NPASS,NFAIL,MSG)
      CALL REPORT("DNA C2' rejected (two hydrogens)", .NOT.IS_CENTRE("C2'"), NPASS,NFAIL)
      CALL REPORT("DNA C3' resolved without special case", IS_CENTRE("C3'"), NPASS,NFAIL)
   END SUBROUTINE T_DETECT_DNA

   !> The C3' escape hatch is gated on RESTYPE=="RNA".  Feed it a 2'-OH sugar
   !! labelled DNA: the branch must NOT fire, so C3' is dropped.  If this test
   !! starts passing C3' as chiral, the RESTYPE guard has been lost.
   SUBROUTINE T_RESTYPE_GUARD(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      WRITE(*,'(A)') " -- RESTYPE guard on the C3' special case"
      CALL RESET(); CALL FIXTURE_RNA_ADENOSINE(.TRUE.,"DNA")
      CALL FIND_CHIRAL_CENTRES()
      CALL REPORT("ribose+DNA label: C3' not resolved", .NOT.IS_CENTRE("C3'"), NPASS,NFAIL)
   END SUBROUTINE T_RESTYPE_GUARD

   !====================================================================
   ! ghost atoms
   !====================================================================

   !> In adenine, C5 is bonded to two other sp2 carbons (C4 and C6).  Without
   !! the EXIT after a C=C match in FIND_DOUBLE_BONDS it is assigned two double
   !! bonds and therefore two ghost atoms, corrupting its priority list.
   SUBROUTINE T_GHOST_ATOMS(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      INTEGER :: J1, J2, NG, NSP2C, WORST, WORSTID
      LOGICAL :: OKONE, OKSLOT
      CHARACTER(LEN=100) :: MSG

      WRITE(*,'(A)') " -- ghost atoms on the adenine ring"
      CALL RESET(); CALL FIXTURE_RNA_ADENOSINE(.TRUE.,"RNA")
      CALL FIND_CHIRAL_CENTRES()

      ! no atom may carry more than one ghost
      WORST = 0; WORSTID = -1
      DO J1=1,NATOMS
         NG = 0
         DO J2=1,6
            IF (BONDEDATS(J1,J2).EQ.0) NG = NG + 1
         END DO
         IF (NG.GT.WORST) THEN
            WORST = NG; WORSTID = J1
         END IF
      END DO
      WRITE(MSG,'(A,I3,A,A)') "max ghosts=",WORST," on atom ", &
                              TRIM(MERGE(AMBER_NAMES(MAX(WORSTID,1)),"none",WORSTID.GT.0))
      CALL REPORT("no atom carries two ghost atoms", WORST.LE.1, NPASS,NFAIL,MSG)

      ! every sp2 carbon must carry exactly one ghost
      OKONE = .TRUE.
      NSP2C = 0
      DO J1=1,NATOMS
         IF ((ELEMENT(J1).NE.6).OR.(NBONDED(J1).NE.3)) CYCLE
         NSP2C = NSP2C + 1
         NG = 0
         DO J2=1,6
            IF (BONDEDATS(J1,J2).EQ.0) NG = NG + 1
         END DO
         IF (NG.NE.1) THEN
            OKONE = .FALSE.
            WRITE(*,'(A,A,A,I2)') "         sp2 carbon ",TRIM(AMBER_NAMES(J1))," has ghosts: ",NG
         END IF
      END DO
      WRITE(MSG,'(A,I3)') "sp2 carbons examined=",NSP2C
      CALL REPORT("every sp2 carbon has exactly one ghost", OKONE, NPASS,NFAIL,MSG)
      CALL REPORT("adenine sp2 carbons found", NSP2C.GE.5, NPASS,NFAIL,MSG)

      ! ghosting must not overrun the six slots
      OKSLOT = .TRUE.
      DO J1=1,NATOMS
         NG = 0
         DO J2=1,6
            IF (BONDEDATS(J1,J2).NE.-1) NG = NG + 1
         END DO
         IF (NG.GT.6) OKSLOT = .FALSE.
      END DO
      CALL REPORT("bonded+ghost entries fit in six slots", OKSLOT, NPASS,NFAIL)
   END SUBROUTINE T_GHOST_ATOMS

   !====================================================================
   ! swap groups
   !====================================================================

   !> Checks GET_ATTACHED_SIZE / FIND_SWAPPABLE_GROUPS on ALA.  Guards three
   !! separate defects: the ghost id 0 being used as an array index, the
   !! chiral centre itself ending up inside its own swap group, and groups
   !! being accepted whose atoms were never collected.
   SUBROUTINE T_SWAP_GROUPS_ALA(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      INTEGER :: CI, J1, ID, NHA, NCB, NBAD
      LOGICAL :: HASCA, HASN, HASC, GOTHA, GOTCB
      CHARACTER(LEN=100) :: MSG

      WRITE(*,'(A)') " -- ALA swap groups"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED()
      CALL FIND_CHIRAL_CENTRES()
      CI = CENTRE_INDEX("CA")
      CALL REPORT("ALA CA found", CI.GT.0, NPASS,NFAIL)
      IF (CI.LE.0) RETURN

      CALL REPORT("ALA CA is swappable", POTENTIAL_SWAPT(CI), NPASS,NFAIL)

      ! the chiral centre must never appear in either group
      HASCA = .FALSE.; HASN = .FALSE.; HASC = .FALSE.
      NHA = 0; NCB = 0; NBAD = 0
      GOTHA = .FALSE.; GOTCB = .FALSE.
      DO J1=1,NSWAPCUT
         ID = SWAPGROUPS(CI,1,J1)
         IF (ID.EQ.ATID("CA")) HASCA = .TRUE.
         IF (ID.EQ.ATID("N"))  HASN  = .TRUE.
         IF (ID.EQ.ATID("C"))  HASC  = .TRUE.
         IF (ID.EQ.ATID("HA")) GOTHA = .TRUE.
         IF (ID.EQ.ATID("CB")) GOTCB = .TRUE.
         IF (ID.GT.0) NHA = NHA + 1
         IF (ID.EQ.0) NBAD = NBAD + 1
         ID = SWAPGROUPS(CI,2,J1)
         IF (ID.EQ.ATID("CA")) HASCA = .TRUE.
         IF (ID.EQ.ATID("N"))  HASN  = .TRUE.
         IF (ID.EQ.ATID("C"))  HASC  = .TRUE.
         IF (ID.EQ.ATID("HA")) GOTHA = .TRUE.
         IF (ID.EQ.ATID("CB")) GOTCB = .TRUE.
         IF (ID.GT.0) NCB = NCB + 1
         IF (ID.EQ.0) NBAD = NBAD + 1
      END DO

      WRITE(MSG,'(A,4I4,A,4I4)') "grp1=",SWAPGROUPS(CI,1,1:NSWAPCUT), &
                                 "  grp2=",SWAPGROUPS(CI,2,1:NSWAPCUT)
      CALL REPORT("chiral centre not inside its own swap group", .NOT.HASCA, NPASS,NFAIL,MSG)
      CALL REPORT("no ghost id 0 in swap groups", NBAD.EQ.0, NPASS,NFAIL,MSG)
      CALL REPORT("swap groups are HA and CB", GOTHA.AND.GOTCB, NPASS,NFAIL,MSG)
      CALL REPORT("backbone N not selected (group too large)", .NOT.HASN, NPASS,NFAIL,MSG)
      CALL REPORT("backbone C not selected (group too large)", .NOT.HASC, NPASS,NFAIL,MSG)
      CALL REPORT("group sizes are 1 and 4", &
                  ((NHA.EQ.1).AND.(NCB.EQ.4)).OR.((NHA.EQ.4).AND.(NCB.EQ.1)), &
                  NPASS,NFAIL,MSG)
   END SUBROUTINE T_SWAP_GROUPS_ALA

   !> THR CB: the hydroxyl group {OG1,HG1} must be collected complete.
   !! A group whose atoms are only partially collected would leave HG1 behind
   !! when the group is rotated, stretching the O-H bond to nonsense.
   SUBROUTINE T_SWAP_GROUPS_THR(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      INTEGER :: CI, J1, ID
      LOGICAL :: GOTOG, GOTHG, GOTHB, HASCB
      CHARACTER(LEN=100) :: MSG

      WRITE(*,'(A)') " -- THR CB swap groups (hydroxyl)"
      CALL RESET(); CALL FIXTURE_THR_CAPPED()
      CALL FIND_CHIRAL_CENTRES()
      CI = CENTRE_INDEX("CB")
      CALL REPORT("THR CB found", CI.GT.0, NPASS,NFAIL)
      IF (CI.LE.0) RETURN
      CALL REPORT("THR CB is swappable", POTENTIAL_SWAPT(CI), NPASS,NFAIL)

      GOTOG=.FALSE.; GOTHG=.FALSE.; GOTHB=.FALSE.; HASCB=.FALSE.
      DO J1=1,NSWAPCUT
         ID = SWAPGROUPS(CI,1,J1)
         IF (ID.EQ.ATID("OG1")) GOTOG=.TRUE.
         IF (ID.EQ.ATID("HG1")) GOTHG=.TRUE.
         IF (ID.EQ.ATID("HB"))  GOTHB=.TRUE.
         IF (ID.EQ.ATID("CB"))  HASCB=.TRUE.
         ID = SWAPGROUPS(CI,2,J1)
         IF (ID.EQ.ATID("OG1")) GOTOG=.TRUE.
         IF (ID.EQ.ATID("HG1")) GOTHG=.TRUE.
         IF (ID.EQ.ATID("HB"))  GOTHB=.TRUE.
         IF (ID.EQ.ATID("CB"))  HASCB=.TRUE.
      END DO
      WRITE(MSG,'(A,4I4,A,4I4)') "grp1=",SWAPGROUPS(CI,1,1:NSWAPCUT), &
                                 "  grp2=",SWAPGROUPS(CI,2,1:NSWAPCUT)
      CALL REPORT("smallest groups are HB and OG1", GOTHB.AND.GOTOG, NPASS,NFAIL,MSG)
      CALL REPORT("hydroxyl H collected with its O", GOTHG, NPASS,NFAIL,MSG)
      CALL REPORT("centre CB not inside its own group", .NOT.HASCB, NPASS,NFAIL,MSG)
   END SUBROUTINE T_SWAP_GROUPS_THR

   !====================================================================
   ! priority integrity
   !====================================================================

   !> CHIR_INFO(J,2:5) is written through FINALPRIO as an index.  If FINALPRIO
   !! is ever not a permutation of 1..4, one slot is written twice and another
   !! keeps whatever was in the uninitialised scratch array - which then gets
   !! used as an XYZ index.  Check every centre of every fixture.
   SUBROUTINE T_PRIORITY_IS_PERMUTATION(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      WRITE(*,'(A)') " -- CHIR_INFO integrity across all fixtures"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED();  CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("ALA",NPASS,NFAIL)
      CALL RESET(); CALL FIXTURE_THR_CAPPED();  CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("THR",NPASS,NFAIL)
      CALL RESET(); CALL FIXTURE_ILE_CAPPED();  CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("ILE",NPASS,NFAIL)
      CALL RESET(); CALL FIXTURE_SYNTH_NC3();   CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("NC3",NPASS,NFAIL)
      CALL RESET(); CALL FIXTURE_HIS_CAPPED();  CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("HIS",NPASS,NFAIL)
      CALL RESET(); CALL FIXTURE_RNA_ADENOSINE(.TRUE.,"RNA"); CALL FIND_CHIRAL_CENTRES()
      CALL CHECK_CHIR_INFO("RNA",NPASS,NFAIL)
   END SUBROUTINE T_PRIORITY_IS_PERMUTATION

   SUBROUTINE CHECK_CHIR_INFO(TAG,NPASS,NFAIL)
      CHARACTER(LEN=*), INTENT(IN) :: TAG
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      INTEGER :: J1, J2, J3, ID
      LOGICAL :: OKRANGE, OKUNIQ, OKBOND
      CHARACTER(LEN=100) :: MSG

      OKRANGE = .TRUE.; OKUNIQ = .TRUE.; OKBOND = .TRUE.
      DO J1=1,NCHIRAL
         DO J2=1,5
            ID = CHIR_INFO(J1,J2)
            IF ((ID.LT.1).OR.(ID.GT.NATOMS)) OKRANGE = .FALSE.
         END DO
         ! the four neighbours must be distinct
         DO J2=2,5
            DO J3=J2+1,5
               IF (CHIR_INFO(J1,J2).EQ.CHIR_INFO(J1,J3)) OKUNIQ = .FALSE.
            END DO
         END DO
         ! and each must actually be bonded to the centre
         DO J2=2,5
            IF (.NOT.ARE_BONDED(CHIR_INFO(J1,1),CHIR_INFO(J1,J2))) OKBOND = .FALSE.
         END DO
      END DO
      WRITE(MSG,'(A,A,A,I3)') "fixture ",TAG,"  NCHIRAL=",NCHIRAL
      CALL REPORT(TAG//": CHIR_INFO ids in range",   OKRANGE, NPASS,NFAIL,MSG)
      CALL REPORT(TAG//": four distinct neighbours", OKUNIQ,  NPASS,NFAIL,MSG)
      CALL REPORT(TAG//": neighbours bonded to centre", OKBOND, NPASS,NFAIL,MSG)
   END SUBROUTINE CHECK_CHIR_INFO

   LOGICAL FUNCTION ARE_BONDED(A,B)
      INTEGER, INTENT(IN) :: A, B
      INTEGER :: J1
      ARE_BONDED = .FALSE.
      DO J1=1,NBOND
         IF (((BONDS(J1,1).EQ.A).AND.(BONDS(J1,2).EQ.B)).OR. &
             ((BONDS(J1,1).EQ.B).AND.(BONDS(J1,2).EQ.A))) ARE_BONDED = .TRUE.
      END DO
   END FUNCTION ARE_BONDED

   !====================================================================
   ! geometry
   !====================================================================

   !> A mirrored copy must get the opposite S/R answer, and the answer must be
   !! stable under a rigid rotation of the whole molecule.
   SUBROUTINE T_ASSIGNMENT_SR(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      REAL(KIND=REAL64), ALLOCATABLE :: XYZ(:)
      REAL(KIND=REAL64) :: NB(12), CE(3)
      LOGICAL :: SR1, SR2
      INTEGER :: CI

      WRITE(*,'(A)') " -- ASSIGNMENT_SR"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED()
      NIMAGES = 0
      CALL FIND_CHIRAL_CENTRES()
      CI = CENTRE_INDEX("CA")
      ALLOCATE(XYZ(3*NATOMS*(NIMAGES+2)))
      XYZ = 0.0D0
      CALL BUILD_COORDS_ALA(XYZ,1)
      CALL MIRROR_IMAGE(XYZ,1,2)

      CALL GATHER(XYZ,1,CI,NB,CE)
      SR1 = ASSIGNMENT_SR(NB,CE)
      CALL GATHER(XYZ,2,CI,NB,CE)
      SR2 = ASSIGNMENT_SR(NB,CE)

      CALL REPORT("mirror image flips the S/R answer", SR1.NEQV.SR2, NPASS,NFAIL)
      DEALLOCATE(XYZ)
   END SUBROUTINE T_ASSIGNMENT_SR

   SUBROUTINE GATHER(XYZ,IMG,CI,NB,CE)
      REAL(KIND=REAL64), INTENT(IN) :: XYZ(:)
      INTEGER, INTENT(IN) :: IMG, CI
      REAL(KIND=REAL64), INTENT(OUT) :: NB(12), CE(3)
      INTEGER :: J1, ID, OFFS
      OFFS = 3*NATOMS*(IMG-1)
      ID = CHIR_INFO(CI,1)
      CE(1) = XYZ(OFFS+3*ID-2); CE(2) = XYZ(OFFS+3*ID-1); CE(3) = XYZ(OFFS+3*ID)
      DO J1=1,4
         ID = CHIR_INFO(CI,J1+1)
         NB(3*J1-2) = XYZ(OFFS+3*ID-2)
         NB(3*J1-1) = XYZ(OFFS+3*ID-1)
         NB(3*J1)   = XYZ(OFFS+3*ID)
      END DO
   END SUBROUTINE GATHER

   !> The core geometry test.  SWITCH_SMALL_GROUPS must be a rigid rotation of
   !! each group, so:
   !!   (a) every moved atom keeps its distance to the chiral centre
   !!   (b) every intra-group bond length is unchanged
   !!   (c) the handedness is actually inverted
   !! (a) and (b) fail if the rotation axis is not normalised or if the
   !! components of the rotated vector are computed from partially updated
   !! values.  (c) is what the routine is for.
   SUBROUTINE T_SWITCH_GEOMETRY(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      REAL(KIND=REAL64), ALLOCATABLE :: XYZ(:), XSAVE(:)
      REAL(KIND=REAL64) :: NB(12), CE(3), P(3), Q(3), C0(3)
      REAL(KIND=REAL64) :: D0, D1, WORSTR, WORSTB
      LOGICAL :: SRBEFORE, SRAFTER, SRREF
      INTEGER :: CI, J1, J2, ID, IDB, OFFS
      CHARACTER(LEN=100) :: MSG
      REAL(KIND=REAL64), PARAMETER :: TOL = 1.0D-9
#ifdef SWAP_HAS_SUCCESS
      LOGICAL :: SWAPOK
#endif

      WRITE(*,'(A)') " -- SWITCH_SMALL_GROUPS rigidity and inversion"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED()
      NIMAGES = 0
      CALL FIND_CHIRAL_CENTRES()
      CI = CENTRE_INDEX("CA")
      ALLOCATE(XYZ(3*NATOMS*(NIMAGES+2)), XSAVE(3*NATOMS*(NIMAGES+2)))
      XYZ = 0.0D0
      CALL BUILD_COORDS_ALA(XYZ,1)
      CALL MIRROR_IMAGE(XYZ,1,2)
      XSAVE = XYZ

      CALL GATHER(XYZ,1,CI,NB,CE)
      SRREF = ASSIGNMENT_SR(NB,CE)
      CALL GATHER(XYZ,2,CI,NB,CE)
      SRBEFORE = ASSIGNMENT_SR(NB,CE)
      CALL REPORT("setup: image 2 starts inverted", SRBEFORE.NEQV.SRREF, NPASS,NFAIL)

#ifdef SWAP_HAS_SUCCESS
      CALL SWITCH_SMALL_GROUPS(XYZ,CI,2,SWAPOK)
      CALL REPORT("switch reported success", SWAPOK, NPASS,NFAIL)
#else
      CALL SWITCH_SMALL_GROUPS(XYZ,CI,2)
#endif

      OFFS = 3*NATOMS
      CALL GETXYZ(XYZ,OFFS,CHIR_INFO(CI,1),C0)

      ! (a) distance to the centre preserved for every moved atom
      WORSTR = 0.0D0
      DO J2=1,2
         DO J1=1,NSWAPCUT
            ID = SWAPGROUPS(CI,J2,J1)
            IF (ID.LE.0) CYCLE
            CALL GETXYZ(XSAVE,OFFS,ID,P)
            CALL GETXYZ(XYZ,  OFFS,ID,Q)
            D0 = DIST(P,C0); D1 = DIST(Q,C0)
            WORSTR = MAX(WORSTR,ABS(D1-D0))
         END DO
      END DO
      WRITE(MSG,'(A,ES12.4)') "max |delta r| to centre = ",WORSTR
      CALL REPORT("rotation preserves distance to centre", WORSTR.LT.TOL, NPASS,NFAIL,MSG)

      ! (b) intra-group bond lengths preserved
      WORSTB = 0.0D0
      DO J2=1,2
         DO J1=1,NSWAPCUT
            ID = SWAPGROUPS(CI,J2,J1)
            IF (ID.LE.0) CYCLE
            DO IDB=1,NATOMS
               IF (.NOT.ARE_BONDED(ID,IDB)) CYCLE
               IF (.NOT.IN_GROUP(CI,J2,IDB)) CYCLE
               CALL GETXYZ(XSAVE,OFFS,ID,P);  CALL GETXYZ(XSAVE,OFFS,IDB,Q)
               D0 = DIST(P,Q)
               CALL GETXYZ(XYZ,OFFS,ID,P);    CALL GETXYZ(XYZ,OFFS,IDB,Q)
               D1 = DIST(P,Q)
               WORSTB = MAX(WORSTB,ABS(D1-D0))
            END DO
         END DO
      END DO
      WRITE(MSG,'(A,ES12.4)') "max |delta bond| = ",WORSTB
      CALL REPORT("rotation preserves intra-group bond lengths", WORSTB.LT.TOL, NPASS,NFAIL,MSG)

      ! (c) handedness restored
      CALL GATHER(XYZ,2,CI,NB,CE)
      SRAFTER = ASSIGNMENT_SR(NB,CE)
      CALL REPORT("switch inverts the centre", SRAFTER.NEQV.SRBEFORE, NPASS,NFAIL)
      CALL REPORT("switch restores reference handedness", SRAFTER.EQV.SRREF, NPASS,NFAIL)

      ! the chiral centre itself must not have moved
      CALL GETXYZ(XSAVE,OFFS,CHIR_INFO(CI,1),P)
      CALL REPORT("chiral centre itself did not move", DIST(P,C0).LT.TOL, NPASS,NFAIL)

      DEALLOCATE(XYZ,XSAVE)
   END SUBROUTINE T_SWITCH_GEOMETRY

   LOGICAL FUNCTION IN_GROUP(CI,G,ID)
      INTEGER, INTENT(IN) :: CI, G, ID
      INTEGER :: J1
      IN_GROUP = .FALSE.
      DO J1=1,NSWAPCUT
         IF (SWAPGROUPS(CI,G,J1).EQ.ID) IN_GROUP = .TRUE.
      END DO
   END FUNCTION IN_GROUP

   REAL(KIND=REAL64) FUNCTION DIST(A,B)
      REAL(KIND=REAL64), INTENT(IN) :: A(3), B(3)
      DIST = SQRT(SUM((A-B)**2))
   END FUNCTION DIST

   !> End-to-end: a two-image band whose second image is inverted must come
   !! out of CHIRALITY_CHECK with both images the same handedness.  This is
   !! the assertion the module never makes about itself - without it a repair
   !! that silently fails just prints the same warning on every image.
   SUBROUTINE T_CHIRALITY_CHECK_BAND(NPASS,NFAIL)
      INTEGER, INTENT(INOUT) :: NPASS, NFAIL
      REAL(KIND=REAL64), ALLOCATABLE :: XYZ(:)
      REAL(KIND=REAL64) :: NB(12), CE(3)
      LOGICAL :: SR1, SR2
      INTEGER :: CI

      WRITE(*,'(A)') " -- CHIRALITY_CHECK end to end"
      CALL RESET(); CALL FIXTURE_ALA_CAPPED()
      NIMAGES = 0
      CALL FIND_CHIRAL_CENTRES()
      CI = CENTRE_INDEX("CA")
      ALLOCATE(XYZ(3*NATOMS*(NIMAGES+2)))
      XYZ = 0.0D0
      CALL BUILD_COORDS_ALA(XYZ,1)
      CALL MIRROR_IMAGE(XYZ,1,2)
      ATOMACTIVE(1:NATOMS) = .TRUE.

      CALL CHIRALITY_CHECK(XYZ)

      CALL GATHER(XYZ,1,CI,NB,CE); SR1 = ASSIGNMENT_SR(NB,CE)
      CALL GATHER(XYZ,2,CI,NB,CE); SR2 = ASSIGNMENT_SR(NB,CE)
      CALL REPORT("band is consistent after the check", SR1.EQV.SR2, NPASS,NFAIL)
      DEALLOCATE(XYZ)
   END SUBROUTINE T_CHIRALITY_CHECK_BAND

END PROGRAM TEST_CHIRALITY_MAIN
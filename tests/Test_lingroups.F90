!====================================================================
! Unit tests for QCI_LINEAR's group-construction logic
! (REDUCE_MULTI_HINGE_GROUP / GROUP_LINEAR_ATOMS and helpers).
!====================================================================
MODULE TEST_ASSERT
   IMPLICIT NONE
   INTEGER :: NPASS = 0
   INTEGER :: NFAIL = 0
CONTAINS
   SUBROUTINE ASSERT_TRUE(COND, MSG)
      LOGICAL, INTENT(IN) :: COND
      CHARACTER(*), INTENT(IN) :: MSG
      IF (COND) THEN
         NPASS = NPASS + 1
         PRINT '(A,A)', "  [PASS] ", MSG
      ELSE
         NFAIL = NFAIL + 1
         PRINT '(A,A)', "  [FAIL] ", MSG
      END IF
   END SUBROUTINE ASSERT_TRUE

   SUBROUTINE ASSERT_EQ_INT(ACTUAL, EXPECTED, MSG)
      INTEGER, INTENT(IN) :: ACTUAL, EXPECTED
      CHARACTER(*), INTENT(IN) :: MSG
      CHARACTER(200) :: FULLMSG
      IF (ACTUAL == EXPECTED) THEN
         NPASS = NPASS + 1
         WRITE(FULLMSG,'(A,A,I0,A)') MSG, " (=", ACTUAL, ")"
         PRINT '(A,A)', "  [PASS] ", TRIM(FULLMSG)
      ELSE
         NFAIL = NFAIL + 1
         WRITE(FULLMSG,'(A,A,I0,A,I0,A)') MSG, " (expected ", EXPECTED, ", got ", ACTUAL, ")"
         PRINT '(A,A)', "  [FAIL] ", TRIM(FULLMSG)
      END IF
   END SUBROUTINE ASSERT_EQ_INT

   SUBROUTINE REPORT_SUMMARY()
      PRINT *, ""
      PRINT '(A,I0,A,I0,A)', "===== ", NPASS, " passed, ", NFAIL, " failed ====="
   END SUBROUTINE REPORT_SUMMARY
END MODULE TEST_ASSERT


PROGRAM TEST_LINEAR_GROUPS
   USE QCIKEYS, ONLY: NATOMS
   USE QCI_CONSTRAINT_KEYS, ONLY: BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
   USE TEST_GRAPH_BUILDER
   USE QCI_LINEAR_CORE
   USE TEST_ASSERT
   IMPLICIT NONE

   LOGICAL :: LINATOM(40)

   PRINT *, "=========================================="
   PRINT *, "TEST 1: chain already single-anchor"
   PRINT *, "=========================================="
   CALL TEST_1_SIMPLE_CHAIN()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 2: isolated fragment (0 anchors) dropped"
   PRINT *, "=========================================="
   CALL TEST_2_ISOLATED_FRAGMENT()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 3: valid single-anchor group below MIN_SIZE dropped"
   PRINT *, "=========================================="
   CALL TEST_3_TOO_SMALL()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 4: chain anchored at both ends -> unsalvageable, drop all"
   PRINT *, "=========================================="
   CALL TEST_4_CHAIN_BOTH_ENDS_ANCHORED()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 5: hub with two clean branches -> splits into two groups"
   PRINT *, "=========================================="
   CALL TEST_5_FORCED_SPLIT()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 6: ring, anchors on opposite atoms -> unsalvageable, drop all"
   PRINT *, "=========================================="
   CALL TEST_6_RING_TWO_ANCHORS()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 7: two independent DFS groups in one run"
   PRINT *, "        (regression test: GROUPID reset bug)"
   PRINT *, "=========================================="
   CALL TEST_7_MULTIPLE_DFS_GROUPS()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 8: dropped atoms leave no stale group id"
   PRINT *, "        (regression test: ATOM2LINGROUP alloc bug)"
   PRINT *, "=========================================="
   CALL TEST_8_ATOM2LINGROUP_SENTINEL()

   PRINT *, ""
   PRINT *, "=========================================="
   PRINT *, "TEST 9: larger irregular topology, invariants only"
   PRINT *, "=========================================="
   CALL TEST_9_LARGER_STRESS()

   CALL REPORT_SUMMARY()
   IF (NFAIL > 0) STOP 1

CONTAINS

   !> Recomputes anchors for LINEAR_GROUPS(GROUP_ID,:) against the full
   !! bond graph and asserts exactly one -- the core single-hinge
   !! contract every committed group must satisfy, independent of any
   !! particular growth/tie-break path that produced it.
   SUBROUTINE ASSERT_SINGLE_HINGE(GROUP_ID, MSG)
      INTEGER, INTENT(IN) :: GROUP_ID
      CHARACTER(*), INTENT(IN) :: MSG
      LOGICAL :: IN_SET(NATOMS)
      INTEGER :: ANCHORS(NATOMS), N_ANCHORS

      CALL RMH_MEMBERSHIP(LINEAR_GROUPS(GROUP_ID,1:NINGROUP(GROUP_ID)), NINGROUP(GROUP_ID), IN_SET)
      CALL RMH_COMPUTE_ANCHORS(LINEAR_GROUPS(GROUP_ID,1:NINGROUP(GROUP_ID)), NINGROUP(GROUP_ID), &
                                IN_SET, ANCHORS, N_ANCHORS)
      CALL ASSERT_EQ_INT(N_ANCHORS, 1, MSG)
   END SUBROUTINE ASSERT_SINGLE_HINGE


   SUBROUTINE TEST_1_SIMPLE_CHAIN()
      INTEGER :: CHAIN(7)
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(8, 4)
      CHAIN = [(J1, J1=1,7)]
      CALL ADD_CHAIN(CHAIN, 7)
      CALL ADD_BOND(1, 8)       ! atom 1 -- external anchor 8

      LINATOM(1:8) = .FALSE.
      LINATOM(1:7) = .TRUE.     ! atom 8 is the external anchor, not a linear atom

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:8))

      CALL ASSERT_EQ_INT(NLINGROUPS, 1, "exactly one group formed")
      IF (NLINGROUPS >= 1) THEN
         CALL ASSERT_EQ_INT(NINGROUP(1), 7, "group keeps all 7 chain atoms")
         CALL ASSERT_SINGLE_HINGE(1, "group has exactly one external anchor")
      END IF
      CALL ASSERT_TRUE(ALL(LINATOM(1:7)), "all chain atoms remain marked linear")
   END SUBROUTINE TEST_1_SIMPLE_CHAIN


   SUBROUTINE TEST_2_ISOLATED_FRAGMENT()
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(6, 4)
      DO J1 = 1, 5
         CALL ADD_BOND(J1, J1+1)
      END DO
      CALL ADD_BOND(6, 1)       ! close the ring -- no external bonds at all

      LINATOM(1:6) = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:6))

      CALL ASSERT_EQ_INT(NLINGROUPS, 0, "no groups formed (nothing to hinge on)")
      CALL ASSERT_TRUE(.NOT. ANY(LINATOM(1:6)), "all 6 atoms dropped from the linear set")
   END SUBROUTINE TEST_2_ISOLATED_FRAGMENT


   SUBROUTINE TEST_3_TOO_SMALL()
      INTEGER :: CHAIN(3)
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(4, 4)
      CHAIN = [(J1, J1=1,3)]
      CALL ADD_CHAIN(CHAIN, 3)
      CALL ADD_BOND(3, 4)       ! atom 3 -- external anchor 4

      LINATOM(1:4) = .FALSE.
      LINATOM(1:3) = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:4))

      CALL ASSERT_EQ_INT(NLINGROUPS, 0, "structurally-valid 3-atom group still dropped (< MIN_LINEAR_GROUP_SIZE)")
      CALL ASSERT_TRUE(.NOT. ANY(LINATOM(1:3)), "all 3 atoms dropped")
   END SUBROUTINE TEST_3_TOO_SMALL


   !> A straight chain with a REAL external anchor at BOTH true ends is
   !! a rigid body pinned at two independent points -- it cannot rotate
   !! as a single-hinge unit no matter where you cut it. Every possible
   !! non-empty subset of this chain touches exactly 2 distinct atoms
   !! outside itself (the anchor on whichever end it reaches, plus the
   !! chain atom just beyond wherever it stops on the other side), so
   !! the correct behaviour is to drop the whole pool, not to keep a
   !! "dominant" N-1 piece that only LOOKS single-hinge because the
   !! excluded neighbour was never re-checked against the full graph.
   SUBROUTINE TEST_4_CHAIN_BOTH_ENDS_ANCHORED()
      INTEGER :: CHAIN(10)
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(12, 4)
      CHAIN = [(J1, J1=1,10)]
      CALL ADD_CHAIN(CHAIN, 10)
      CALL ADD_BOND(1, 11)      ! anchor A at one end
      CALL ADD_BOND(10, 12)     ! anchor B at the other end

      LINATOM(1:12) = .FALSE.
      LINATOM(1:10) = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:12))

      CALL ASSERT_EQ_INT(NLINGROUPS, 0, "no valid single-hinge subset exists; whole pool dropped")
      CALL ASSERT_TRUE(.NOT. ANY(LINATOM(1:10)), "all 10 chain atoms dropped")
   END SUBROUTINE TEST_4_CHAIN_BOTH_ENDS_ANCHORED


   !> Hub H bonded to two REAL external anchors of its own (EA, EB),
   !! plus two "clean" branches (P and Q, 5 atoms each) that carry no
   !! anchor of their own. H itself can never be part of any group (it
   !! touches 2 real anchors directly), but each clean branch CAN be
   !! cleanly cut off at its single bond to H, giving each branch its
   !! own valid single-hinge group hinged on H. This is a genuine
   !! one-pool-splits-into-two-groups case, and also demonstrates a
   !! hinge that is an ordinary pool atom (H) rather than a real
   !! external anchor -- H itself ends up dropped once both branches
   !! are peeled off.
   SUBROUTINE TEST_5_FORCED_SPLIT()
      INTEGER :: PCHAIN(5), QCHAIN(5)
      INTEGER :: J1
      INTEGER :: NDROPPED, NCOMMITTED

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(13, 5)
      ! pool atoms: P1-5 = 1-5, H = 6, Q1-5 = 7-11 ; anchors: EA=12, EB=13 (both on H)
      PCHAIN = [(J1, J1=1,5)]
      CALL ADD_CHAIN(PCHAIN, 5)
      CALL ADD_BOND(6, 1)       ! H -- P1
      QCHAIN = [(J1, J1=7,11)]
      CALL ADD_CHAIN(QCHAIN, 5)
      CALL ADD_BOND(6, 7)       ! H -- Q1
      CALL ADD_BOND(6, 12)      ! EA on H
      CALL ADD_BOND(6, 13)      ! EB on H

      LINATOM(1:13) = .FALSE.
      LINATOM(1:11) = .TRUE.    ! atoms 1-11 are the pool; 12-13 are external anchors

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:13))

      CALL ASSERT_EQ_INT(NLINGROUPS, 2, "pool resolves into exactly two groups")

      IF (NLINGROUPS == 2) THEN
         CALL ASSERT_EQ_INT(NINGROUP(1), 5, "first group has 5 atoms")
         CALL ASSERT_EQ_INT(NINGROUP(2), 5, "second group has 5 atoms")
         CALL ASSERT_SINGLE_HINGE(1, "first group has exactly one external anchor")
         CALL ASSERT_SINGLE_HINGE(2, "second group has exactly one external anchor")

         NCOMMITTED = COUNT(ATOM2LINGROUP(1:11) > 0)
         NDROPPED   = COUNT((.NOT. LINATOM(1:11)) .AND. ATOM2LINGROUP(1:11) == -1)
         CALL ASSERT_EQ_INT(NCOMMITTED, 10, "10 of the 11 pool atoms ended up committed (both branches)")
         CALL ASSERT_EQ_INT(NDROPPED, 1, "hub H itself was dropped (touches 2 real anchors)")
         CALL ASSERT_TRUE(ATOM2LINGROUP(6) == -1, "H specifically is the dropped atom")
      END IF
   END SUBROUTINE TEST_5_FORCED_SPLIT


   !> 8-membered ring with REAL anchors on two opposite atoms. Just
   !! like the two-ended chain, a ring pinned at 2 independent points
   !! is over-constrained: any arc you keep still touches both the
   !! anchor at whichever end it reaches AND the excluded ring atom at
   !! the other end (going around either direction doesn't help, since
   !! a ring offers no THIRD path that avoids both anchors). Correct
   !! behaviour is to drop the whole ring.
   SUBROUTINE TEST_6_RING_TWO_ANCHORS()
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(10, 4)
      DO J1 = 1, 7
         CALL ADD_BOND(J1, J1+1)
      END DO
      CALL ADD_BOND(8, 1)       ! close the 8-ring
      CALL ADD_BOND(1, 9)       ! EA on atom 1
      CALL ADD_BOND(5, 10)      ! EB on atom 5 (opposite side of the ring)

      LINATOM(1:10) = .FALSE.
      LINATOM(1:8)  = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:10))

      CALL ASSERT_EQ_INT(NLINGROUPS, 0, "no valid single-hinge subset exists; whole ring dropped")
      CALL ASSERT_TRUE(.NOT. ANY(LINATOM(1:8)), "all 8 ring atoms dropped")
   END SUBROUTINE TEST_6_RING_TWO_ANCHORS


   !> Two fully disconnected 6-atom chains, each already single-anchor.
   !! GROUP_LINEAR_ATOMS must find 2 separate DFS groups and commit
   !! BOTH correctly. Before the fix, GROUPID was left holding the DFS
   !! group *count* (2) when REDUCE_MULTI_HINGE_GROUP's commit loop
   !! started, so real commits landed at slots 3 and 4 while slots 1
   !! and 2 stayed as empty phantom groups (NINGROUP==0) -- this test
   !! fails immediately on that bug via the NINGROUP(1) check.
   SUBROUTINE TEST_7_MULTIPLE_DFS_GROUPS()
      INTEGER :: CHAIN1(6), CHAIN2(6)
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(14, 4)
      CHAIN1 = [(J1, J1=1,6)]
      CALL ADD_CHAIN(CHAIN1, 6)
      CALL ADD_BOND(6, 7)        ! anchor for fragment 1

      CHAIN2 = [(J1, J1=8,13)]
      CALL ADD_CHAIN(CHAIN2, 6)
      CALL ADD_BOND(13, 14)      ! anchor for fragment 2

      LINATOM(1:14) = .FALSE.
      LINATOM(1:6)  = .TRUE.
      LINATOM(8:13) = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:14))

      CALL ASSERT_EQ_INT(NLINGROUPS, 2, "two independent groups formed")
      IF (NLINGROUPS == 2) THEN
         CALL ASSERT_EQ_INT(NINGROUP(1), 6, "slot 1 actually holds fragment 1 (not a phantom empty group)")
         CALL ASSERT_EQ_INT(NINGROUP(2), 6, "slot 2 actually holds fragment 2 (not a phantom empty group)")
         CALL ASSERT_EQ_INT(NQCILINEAR, 12, "total linear atom count across both groups is 12")
      END IF
   END SUBROUTINE TEST_7_MULTIPLE_DFS_GROUPS


   !> Mixed run: one fragment big enough to commit, one fragment too
   !! small to survive. Confirms ATOM2LINGROUP (allocated lazily on
   !! first commit) correctly distinguishes committed atoms from
   !! dropped ones, and that indexing it never crashes even though it
   !! is allocated *after* the dropped fragment is already processed.
   SUBROUTINE TEST_8_ATOM2LINGROUP_SENTINEL()
      INTEGER :: SMALLCHAIN(3), BIGCHAIN(6)
      INTEGER :: J1

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(11, 4)
      ! fragment A: too small, processed FIRST (lower atom indices),
      ! so no group has been committed yet when it's handled --
      ! exercises RMH_ENSURE_GROUP_CAPACITY/ATOM2LINGROUP allocation
      ! not yet having happened at all.
      SMALLCHAIN = [(J1, J1=1,3)]
      CALL ADD_CHAIN(SMALLCHAIN, 3)
      CALL ADD_BOND(3, 4)

      ! fragment B: big enough, processed second
      BIGCHAIN = [(J1, J1=5,10)]
      CALL ADD_CHAIN(BIGCHAIN, 6)
      CALL ADD_BOND(10, 11)

      LINATOM(1:11) = .FALSE.
      LINATOM(1:3)  = .TRUE.
      LINATOM(5:10) = .TRUE.

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:11))

      CALL ASSERT_EQ_INT(NLINGROUPS, 1, "only fragment B survives")
      CALL ASSERT_TRUE(ALLOCATED(ATOM2LINGROUP), "ATOM2LINGROUP got allocated once a group was committed")
      IF (ALLOCATED(ATOM2LINGROUP)) THEN
         CALL ASSERT_TRUE(ALL(ATOM2LINGROUP(1:3) == -1), "dropped fragment A atoms are unassigned (-1)")
         CALL ASSERT_TRUE(ALL(ATOM2LINGROUP(5:10) > 0), "surviving fragment B atoms are all assigned a group")
      END IF
   END SUBROUTINE TEST_8_ATOM2LINGROUP_SENTINEL

   !> A 15-atom spine anchored at both ends, with a clean 6-atom branch
   !! (no anchor) hanging off the middle, a 6-atom branch with its own
   !! anchor hanging off elsewhere, and a small unanchored ring fused
   !! onto the spine for extra structural variety. Deliberately NOT
   !! hand-solved -- too easy to make the same kind of manual error
   !! that caused Tests 4-6 to need correcting in the first place.
   !! Instead this only checks the invariants that MUST hold for any
   !! correct output, atom-by-atom, against the full bond graph.
   SUBROUTINE TEST_9_LARGER_STRESS()
      INTEGER :: SPINE(15), BRANCHA(6), BRANCHB(6), RING(3)
      INTEGER :: J1, GID
      INTEGER :: NCOMMITTED, NDROPPED, NUNACCOUNTED
      LOGICAL :: NO_ATOM_LEFT_DANGLING
      LOGICAL :: ALL_GROUPS_BIG_ENOUGH
      LOGICAL :: ALL_GROUPS_SINGLE_HINGE

      CALL DEALLOC_LINEAR_GROUP()
      CALL INIT_GRAPH(33, 5)

      SPINE = [(J1, J1=1,15)]
      CALL ADD_CHAIN(SPINE, 15)

      BRANCHA = [(J1, J1=16,21)]
      CALL ADD_BOND(4, 16)
      CALL ADD_CHAIN(BRANCHA, 6)
      CALL ADD_BOND(21, 33)      ! EA on branch-A tip

      BRANCHB = [(J1, J1=22,27)]
      CALL ADD_BOND(10, 22)
      CALL ADD_CHAIN(BRANCHB, 6) ! branch B tip (27) has no anchor -- clean

      RING = [28, 29, 30]
      CALL ADD_BOND(7, 28)
      CALL ADD_CHAIN(RING, 3)
      CALL ADD_BOND(30, 7)       ! closes a small ring fused onto the spine

      CALL ADD_BOND(1, 31)       ! E_START on spine atom 1
      CALL ADD_BOND(15, 32)      ! E_END on spine atom 15

      LINATOM(1:33)  = .FALSE.
      LINATOM(1:30)  = .TRUE.    ! atoms 1-30 are the pool; 31-33 are external anchors

      CALL GROUP_LINEAR_ATOMS(LINATOM(1:33))

      PRINT '(A,I0,A)', "  (", NLINGROUPS, " group(s) formed -- not hand-predicted, checking invariants)"

      ALL_GROUPS_BIG_ENOUGH  = .TRUE.
      ALL_GROUPS_SINGLE_HINGE = .TRUE.
      DO GID = 1, NLINGROUPS
         IF (NINGROUP(GID) < MIN_LINEAR_GROUP_SIZE) ALL_GROUPS_BIG_ENOUGH = .FALSE.
      END DO
      CALL ASSERT_TRUE(ALL_GROUPS_BIG_ENOUGH, "every committed group meets MIN_LINEAR_GROUP_SIZE")

      DO GID = 1, NLINGROUPS
         CALL ASSERT_SINGLE_HINGE(GID, "committed group is genuinely single-hinge (full-graph recheck)")
      END DO

      IF (ALLOCATED(ATOM2LINGROUP)) THEN
         NCOMMITTED = COUNT(ATOM2LINGROUP(1:30) > 0)
      ELSE
         NCOMMITTED = 0
      END IF
      NDROPPED = COUNT(.NOT. LINATOM(1:30))
      CALL ASSERT_EQ_INT(NCOMMITTED + NDROPPED, 30, "every pool atom is either committed or dropped, none missing")

      NO_ATOM_LEFT_DANGLING = .TRUE.
      DO J1 = 1, 30
         IF (LINATOM(J1)) THEN
            IF (.NOT. ALLOCATED(ATOM2LINGROUP)) THEN
               NO_ATOM_LEFT_DANGLING = .FALSE.
            ELSE IF (ATOM2LINGROUP(J1) <= 0) THEN
               NO_ATOM_LEFT_DANGLING = .FALSE.
            END IF
         END IF
      END DO
      CALL ASSERT_TRUE(NO_ATOM_LEFT_DANGLING, "no atom is still marked linear without a real group assignment")

      NUNACCOUNTED = COUNT(NINGROUP(1:MAX(NLINGROUPS,1)) < 0)   ! sanity: NINGROUP never negative
      CALL ASSERT_TRUE(NUNACCOUNTED == 0, "group sizes are never negative/uninitialised")
   END SUBROUTINE TEST_9_LARGER_STRESS

END PROGRAM TEST_LINEAR_GROUPS
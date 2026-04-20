!! @TODO CONTINUE here. 
!! External anchors does not give correct results 

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
   ! list of linear groups
   INTEGER, ALLOCATABLE :: LINEAR_GROUPS(:,:)
   INTEGER, ALLOCATABLE :: ATOM2LINGROUP(:)
   INTEGER, ALLOCATABLE :: NINGROUP(:)          ! Number of atoms in each group
   INTEGER :: NLINGROUPS = 0
   CONTAINS

      SUBROUTINE GET_LINEAR_ATOMS()
         USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET, QCILINEART
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

      SUBROUTINE ALLOC_LINEAR_GROUP(NGROUPS, MAXGROUPSIZE)
         USE QCIKEYS, ONLY: NATOMS
         INTEGER, INTENT(IN) :: NGROUPS, MAXGROUPSIZE
         ALLOCATE(LINEAR_GROUPS(NGROUPS,MAXGROUPSIZE))
         ALLOCATE(ATOM2LINGROUP(NATOMS))
         ALLOCATE(NINGROUP(NGROUPS))
      END SUBROUTINE ALLOC_LINEAR_GROUP

      SUBROUTINE DEALLOC_LINEAR_GROUP()
         USE QCI_CONSTRAINT_KEYS, ONLY: NBONDS, BOND_LIST, BONDS_PER_ATOM_LIST ,N_BONDS_PER_ATOM
         IF (ALLOCATED(LINEAR_GROUPS)) DEALLOCATE(LINEAR_GROUPS)
         IF (ALLOCATED(ATOM2LINGROUP)) DEALLOCATE(ATOM2LINGROUP)
         IF (ALLOCATED(NINGROUP)) DEALLOCATE(NINGROUP) 
         IF (ALLOCATED(BOND_LIST)) DEALLOCATE(BOND_LIST)
         IF (ALLOCATED(BONDS_PER_ATOM_LIST)) DEALLOCATE(BONDS_PER_ATOM_LIST)
         IF (ALLOCATED(N_BONDS_PER_ATOM)) DEALLOCATE(N_BONDS_PER_ATOM)
      END SUBROUTINE DEALLOC_LINEAR_GROUP
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !quasi-rigid body linear groups
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE GET_LINEAR_GROUPS()   
         USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET       
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS, ANGLE, DIHEDRAL, DIHEDRAL_DIFF
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL, NCONPERATOM, &
                                       CONLIST, BOND_LIST, NBONDS, MAX_BONDS_PER_ATOM, &
                                        BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
         USE QCICONSTRAINTS, ONLY: GET_NBONDS_PER_ATOM
         USE AMBER_CONSTRAINTS, ONLY: INGROUP, GROUPLOOKUP, NPLACINGGROUPS, SIZEPLACINGGROUPS
         IMPLICIT NONE

         REAL(KIND=REAL64), PARAMETER :: TOLERANCE = 0.1D0
         REAL(KIND=REAL64) :: DS, DF
         LOGICAL :: LINATOM(NATOMS) 
         INTEGER :: NLINHERE
         INTEGER :: A, B, C, D, J1, J2, J3, J4, K

         INTEGER :: CURRENTGROUP(NATOMS)      ! Group ID for each linear atom, -1 if not assigned
         INTEGER :: GROUPS(NATOMS, NATOMS)    ! List of atoms in each group
         INTEGER :: GROUPID                   ! Total number of groups
         INTEGER :: VISITED(NATOMS)           ! For DFS traversal
         INTEGER :: STACK(NATOMS)             ! Stack for DFS
         INTEGER :: STACKPTR
         INTEGER :: ATOM, NEIGHBOR
         INTEGER :: NINGROUP_TEMP(3*NATOMS)
         INTEGER :: MAXGROUPSIZE
         INTEGER :: NGROUPS                     !number of linear groups
                
        
         ! ========== NEW: Attachment point tracking ==========
         INTEGER :: EXTERNAL_ANCHORS(NATOMS)  ! distinct external atoms bonded to this group 
         INTEGER :: N_EXTERNAL_ANCHORS        ! number of distinct external anchor atoms
         INTEGER :: HINGE_ATOM                ! unique external anchor (if exists)

  
         LOGICAL :: DEBUG_MODE = .TRUE.       ! Set to .TRUE. for verbose output
         
         INTEGER :: EXTERNAL_PER_ATOM(NATOMS)
         INTEGER :: N_BACKBONE_ATOMS
         INTEGER :: BACKBONE_ATOM
         LOGICAL :: IS_VALID_GROUP
         INTEGER :: TEMP_EXTERNAL_COUNT
         INTEGER :: NEIGHBOR_ATOM
      
         ! ========== Compare Amber groups with Linear groups ==========
         INTEGER :: AMBER_ATOM, AMBER_GROUP
         INTEGER :: NATOMS_THIS_GROUP
         INTEGER :: ATOMS_IN_BOTH(NATOMS)
         INTEGER :: ATOMS_ONLY_IN_AMBER(NATOMS)
         INTEGER :: N_IN_BOTH, N_ONLY_AMBER

         

        
         LINATOM(:) = .FALSE.
         

         !Get bonds arrays we need - cannot use contrainsts here
         CALL GET_NBONDS_PER_ATOM()
         

         !Part 1: Find linear atoms based on bonds, angles, and dihedrals constraints
         CALL FIND_LINEAR_ATOMS(LINATOM)  

     
         !Part 2: Group LINEAR atoms by bond connectivity 
         
         CURRENTGROUP(:)   = -1
         GROUPS(:,:)       = -1
         NINGROUP_TEMP(:)  = 0
         NGROUPS           = 0

         VISITED(:) = 0

         DO J1 = 1, NATOMS

            ! Start DFS only from unvisited LINEAR atoms
            IF (.NOT. LINATOM(J1)) CYCLE
            IF (VISITED(J1) == 1) CYCLE

            ! New group
            NGROUPS = NGROUPS + 1
            GROUPID = NGROUPS

            STACKPTR   = 1
            STACK(1)   = J1
            VISITED(J1)= 1

            ! -------- DFS --------
            DO WHILE (STACKPTR > 0)

               ATOM = STACK(STACKPTR)
               STACKPTR = STACKPTR - 1

               ! Add atom to current group
               NINGROUP_TEMP(GROUPID) = NINGROUP_TEMP(GROUPID) + 1
               GROUPS(GROUPID, NINGROUP_TEMP(GROUPID)) = ATOM

               ! Traverse bonded neighbors
               DO J2 = 1, N_BONDS_PER_ATOM(ATOM)

                  NEIGHBOR = BONDS_PER_ATOM_LIST(ATOM, J2)

                  ! Skip non-linear atoms
                  IF (.NOT. LINATOM(NEIGHBOR)) CYCLE

                  ! Skip already visited atoms
                  IF (VISITED(NEIGHBOR) == 1) CYCLE

                  ! Visit neighbor
                  VISITED(NEIGHBOR) = 1
                  STACKPTR = STACKPTR + 1

                  IF (STACKPTR > NATOMS) THEN
                     PRINT *, "ERROR: DFS stack overflow"
                     STOP
                  END IF

                  STACK(STACKPTR) = NEIGHBOR

               END DO
            END DO
            ! -------- end DFS --------

            ! Assign group ID *after* full traversal
            DO J2 = 1, NINGROUP_TEMP(GROUPID)
               ATOM = GROUPS(GROUPID, J2)
               CURRENTGROUP(ATOM) = GROUPID
            END DO

         END DO

      ! ========== PART 3: Validate attachment points and filter groups ==========
      
      ! First pass: identify single-atom groups and unmark them
      DO J1 = 1, NGROUPS
         IF (NINGROUP_TEMP(J1).LE.1) THEN
            DO J2 = 1, NINGROUP_TEMP(J1)
               ATOM = GROUPS(J1, J2)
               LINATOM(ATOM) = .FALSE.
               !IF (DEBUG_MODE) PRINT *, "Unmarking atom ", ATOM, " (single-atom group)"
            END DO
         END IF
      END DO

      ! Allocate output arrays
      MAXGROUPSIZE = MAXVAL(NINGROUP_TEMP)
      
      !Should this be done later?
      CALL ALLOC_LINEAR_GROUP(NGROUPS, MAXGROUPSIZE)

      ! Initialize output arrays
      ATOM2LINGROUP(:) = -1
      LINEAR_GROUPS(:,:) = -1
      NINGROUP(:) = 0
     
      NQCILINEAR = 0
      NLINGROUPS = 0

  
      ! ========== PART 3A: Identify distinct external anchor atoms ==========
      ! ========== PART 3B: Validate group topology =========================

      GROUPID = 0

      DO J1 = 1, NGROUPS

         ! Skip trivial groups - shouldn't happen.
         IF (NINGROUP_TEMP(J1) <= 1) CYCLE

         ! Reset external anchor tracking
         EXTERNAL_ANCHORS(:) = -1
         N_EXTERNAL_ANCHORS  = 0
         HINGE_ATOM          = -1
         IS_VALID_GROUP      = .FALSE.

         ! Loop over atoms in this group
         DO J2 = 1, NINGROUP_TEMP(J1)
            ATOM = GROUPS(J1, J2)

            ! Loop over bonded neighbors
            DO J3 = 1, N_BONDS_PER_ATOM(ATOM)

               NEIGHBOR_ATOM = BONDS_PER_ATOM_LIST(ATOM, J3)

               ! Validate neighbor index
               IF ( (NEIGHBOR_ATOM.LT.1).OR. (NEIGHBOR_ATOM.GT.NATOMS) ) THEN
                  PRINT *, "ERROR: Invalid neighbor index ", NEIGHBOR_ATOM
                  CYCLE
               END IF

               ! Bond crosses the group boundary
               ! If neighbor is not in the same groupp, it's an external anchor
               IF (CURRENTGROUP(NEIGHBOR_ATOM).NE.CURRENTGROUP(ATOM)) THEN

                  ! Record only DISTINCT external atoms
                  IF (.NOT. IS_IN_LIST(NEIGHBOR_ATOM, EXTERNAL_ANCHORS, N_EXTERNAL_ANCHORS)) THEN

                     N_EXTERNAL_ANCHORS = N_EXTERNAL_ANCHORS + 1
                     EXTERNAL_ANCHORS(N_EXTERNAL_ANCHORS) = NEIGHBOR_ATOM

                     ! First external anchor becomes the hinge
                     IF (N_EXTERNAL_ANCHORS == 1) THEN
                        HINGE_ATOM = NEIGHBOR_ATOM   
                     END IF
                  END IF

               END IF

            END DO   ! J3
         END DO      ! J2

         ! ========== PART 3B: Validate group ==========
         ! A valid group interacts with the rest of the molecule
         ! through exactly ONE distinct external atom (hinge)
    

         IF (N_EXTERNAL_ANCHORS == 1) THEN
            IS_VALID_GROUP = .TRUE.
            IF (DEBUG_MODE) THEN
               PRINT *, "Accepting group ", J1
               PRINT *, "  Number of distinct external anchors: ", &
                        N_EXTERNAL_ANCHORS
               PRINT *, "  External anchor (hinge): ", HINGE_ATOM
               PRINT *, "  Group atoms: ", &
                        (GROUPS(J1,K), K=1,NINGROUP_TEMP(J1))
            END IF
         ELSE
            IF(DEBUG_MODE) THEN
               PRINT *, "Group ", J1, " has ", N_EXTERNAL_ANCHORS, " external anchors. Attempting to reduce group..."
               PRINT *, "  External anchors: ", &
                        (EXTERNAL_ANCHORS(K), K=1,N_EXTERNAL_ANCHORS)
               PRINT *, "  Group atoms before reduction: ", &
                        (GROUPS(J1,K), K=1,NINGROUP_TEMP(J1))
            END IF
            
            
            !Try to create a valid group by reducing the curent group
            !CALL REDUCE_MULTI_HINGE_GROUP(LINATOM, GROUPS(J1,:), NINGROUP_TEMP(J1), IS_VALID_GROUP, CURRENTGROUP)        

         END IF

         
         ! ========== PART 3C: Copy valid groups to output arrays ==========
 

         IF (IS_VALID_GROUP) THEN
            GROUPID = GROUPID + 1
            NINGROUP(GROUPID) = NINGROUP_TEMP(J1)
           
            !Copy atoms from temporary storage to final arrays
            DO J2 = 1, NINGROUP_TEMP(J1)
               ATOM = GROUPS(J1, J2)
               LINEAR_GROUPS(GROUPID, J2) = ATOM
               ATOM2LINGROUP(ATOM) = GROUPID
               NQCILINEAR = NQCILINEAR + 1
            END DO
            
            !Print group information
            PRINT *, "Valid Linear Group - mew group id: ", GROUPID, ": ", NINGROUP(GROUPID), " atoms"
            PRINT *, "  Hinge atom (attachment point): ", HINGE_ATOM
            PRINT *, "  Group atoms: ", (LINEAR_GROUPS(GROUPID, J2), J2=1, NINGROUP(GROUPID))
            
         ELSE
            ! Unmark atoms in invalid groups
            DO J2 = 1, NINGROUP_TEMP(J1)
               ATOM = GROUPS(J1, J2)
               LINATOM(ATOM) = .FALSE.
            END DO
         END IF
      END DO

      ! Store final count of valid groups
      NLINGROUPS = GROUPID

      PRINT *, "========================================================"
      ! ========== Validation ==========
      PRINT *, ""
      PRINT *, "========== LINEAR GROUP DETECTION SUMMARY =========="
      PRINT *, "Total linear groups found: ", NLINGROUPS
      PRINT *, "Total linear atoms: ", NQCILINEAR
      PRINT *, "Linear atoms remaining: ", COUNT(LINATOM(:))

      IF (COUNT(LINATOM(:)) /= NQCILINEAR) THEN
         PRINT *, "WARNING: Mismatch between LINATOM count and NQCILINEAR"
         PRINT *, "LINATOM count: ", COUNT(LINATOM(:))
         PRINT *, "NQCILINEAR: ", NQCILINEAR
      END IF
      
      PRINT *, "====================================================="
   END SUBROUTINE GET_LINEAR_GROUPS

   
   SUBROUTINE GET_LIN_ROT_TRANSLATION(GROUP, GROUPSIZE, CXS, CXF, QSTART, QFINAL)
      USE QCIKEYS, ONLY: NATOMS
      USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
      USE QCIMINDIST, ONLY: FIND_ORIGIN, MOVE_COORDS, FIND_ALIGNMENT
      USE QUATERNIONS

      INTEGER, INTENT(IN) :: GROUP !< group id
      INTEGER, INTENT(IN) :: GROUPSIZE 
      REAL(KIND=REAL64), INTENT(OUT) :: CXS(3), CXF(3) !< centre of coordinates 
      REAL(KIND=REAL64), INTENT(OUT) :: QSTART(4), QFINAL(4) !< unit quaternions for SLERP
      REAL(KIND=REAL64) :: RS(3*GROUPSIZE), RF(3*GROUPSIZE) !group coordinates
      REAL(KIND=REAL64), DIMENSION(3,3) :: RMAT !rotational matrix 
      REAL(KIND=REAL64) :: DIST
      INTEGER :: J1, ATOMID

      DO J1=1, GROUPSIZE
         ATOMID = LINEAR_GROUPS(GROUP,J1)
         RS(3*(J1-1)+1:3*(J1-1)+3) = XSTART(3*(ATOMID-1)+1:3*(ATOMID-1)+3)
         RF(3*(J1-1)+1:3*(J1-1)+3) = XFINAL(3*(ATOMID-1)+1:3*(ATOMID-1)+3)
      END DO
      
      !centre RA
      CALL FIND_ORIGIN(GROUPSIZE,RS,CXS)
      CALL MOVE_COORDS(GROUPSIZE,RS,CXS)
      !centre RB around origin
      CALL FIND_ORIGIN(GROUPSIZE,RF,CXF)
      CALL MOVE_COORDS(GROUPSIZE,RF,CXF)

      !WRITE(*,*) "get_lin_rot_translation> GROUPSIZE ", GROUPSIZE, "RS ", RS, " RF ", RF
      !align coordinates
      CALL FIND_ALIGNMENT(GROUPSIZE, RS, RF, DIST, RMAT)

      QSTART = [1.0D0, 0.0D0, 0.0D0, 0.0D0]

      CALL MATRIX_TO_QUATERNION(RMAT,QFINAL)
      
   END SUBROUTINE GET_LIN_ROT_TRANSLATION

   SUBROUTINE FIND_LINEAR_ATOMS(LINATOM)
       USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET       
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS, ANGLE, DIHEDRAL, DIHEDRAL_DIFF
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL, NCONPERATOM, &
                                       CONLIST, BOND_LIST, NBONDS, MAX_BONDS_PER_ATOM, &
                                        BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
         USE QCICONSTRAINTS, ONLY: GET_NBONDS_PER_ATOM
         IMPLICIT NONE

         LOGICAL, INTENT(OUT) :: LINATOM(NATOMS) 

         REAL(KIND=REAL64), PARAMETER :: TOLERANCE = 0.1D0
         REAL(KIND=REAL64) :: DS, DF
         
         INTEGER :: NLINHERE
         INTEGER :: A, B, C, D, J1, J2, J3, J4, K
    
          ! Grouping variables
         INTEGER :: CURRENTGROUP(NATOMS)      ! Group ID for each linear atom
         INTEGER :: GROUPS(NATOMS, NATOMS)    ! List of atoms in each group
         INTEGER :: GROUPID
         INTEGER :: VISITED(NATOMS)           ! For DFS traversal
         INTEGER :: STACK(NATOMS)             ! Stack for DFS
         INTEGER :: STACKPTR
         INTEGER :: ATOM, NEIGHBOR
         INTEGER :: NINGROUP_TEMP(3*NATOMS)
         INTEGER :: MAXGROUPSIZE 
         INTEGER :: NGROUPS   !number of linear groups
         
         LOGICAL :: IS_RIGID_BODY(NATOMS)
         REAL(KIND=REAL64) :: ANGLE_COORDS(9)
         REAL(KIND=REAL64) :: ANGLE_DEVIATION, DIHEDRAL_DEVIATION
         REAL(KIND=REAL64) :: ANGLE_START, ANGLE_FINAL
         REAL(KIND=REAL64) :: ANGLE_TOLERANCE = 0.08726646259971647 !5.0D0 degrees  !0.1745329251994329 !
         REAL(KIND=REAL64) :: DIHEDRAL_TOLERANCE = 0.17453292519943295  !10.0D0  ! degrees
         REAL(KIND=REAL64) :: DIH_COORDS(12)
         REAL(KIND=REAL64) :: DIH_START, DIH_FINAL
         
         LINATOM(:) = .FALSE.
         IS_RIGID_BODY(:) = .TRUE.

         !Get bonds arrays we need
         !CALL GET_NBONDS_PER_ATOM()
        
          ! ========== PART 1: Detect linear atoms ==========
         DO J1 = 1, NATOMS
            NLINHERE = 0
         
            DO J2 = 1, N_BONDS_PER_ATOM(J1) !NCONPERATOM(J1)
               A = J1
               B = BONDS_PER_ATOM_LIST(A,J2) !CONLIST(A, J2)

               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, A, B, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, A, B, DF)

               IF(DABS(DS - DF) < TOLERANCE) NLINHERE = NLINHERE + 1
            END DO

            ! Mark atom as linear if all constraints have unchanged distances
            !IF(NLINHERE == NCONPERATOM(J1)) LINATOM(J1) = .TRUE.
            IF(NLINHERE == N_BONDS_PER_ATOM(J1)) LINATOM(J1) = .TRUE.
            
          END DO
    
          
         ! Check if atoms in a group maintain angles
         DO J1 = 1, NATOMS
            IF (.NOT. LINATOM(J1)) CYCLE
            
            ! For each atom, check angles to neighbors
            DO J2 = 1, N_BONDS_PER_ATOM(J1) !NCONPERATOM(J1)
               DO J3 = J2+1, N_BONDS_PER_ATOM(J1) !NCONPERATOM(J1)
                     !A = CONLIST(J1, J2)
                     A = BONDS_PER_ATOM_LIST(J1, J2)
                     B = J1
                     C = BONDS_PER_ATOM_LIST(J1, J3)
                     !C = CONLIST(J1, J3)
                     
                     ANGLE_COORDS(1:3) = XSTART(3*(A-1)+1:3*(A-1)+3)
                     ANGLE_COORDS(4:6) = XSTART(3*(B-1)+1:3*(B-1)+3)
                     ANGLE_COORDS(7:9) = XSTART(3*(C-1)+1:3*(C-1)+3)

                     ANGLE_START = ANGLE(ANGLE_COORDS)
                     
                     ANGLE_COORDS(1:3) = XFINAL(3*(A-1)+1:3*(A-1)+3)
                     ANGLE_COORDS(4:6) = XFINAL(3*(B-1)+1:3*(B-1)+3)
                     ANGLE_COORDS(7:9) = XFINAL(3*(C-1)+1:3*(C-1)+3)


                     ANGLE_FINAL = ANGLE(ANGLE_COORDS)
                     
                     ANGLE_DEVIATION = DABS(ANGLE_START - ANGLE_FINAL)
                     
                     IF (ANGLE_DEVIATION > ANGLE_TOLERANCE) THEN
                        IS_RIGID_BODY(J1) = .FALSE.
                        EXIT
                     END IF
               END DO
               IF (.NOT. IS_RIGID_BODY(J1)) EXIT
            END DO
         END DO

         !dihedrals check
         DO J1 = 1, NATOMS
            IF (.NOT. LINATOM(J1)) CYCLE

            ! B = central atom
            B = J1

            DO J2 = 1, N_BONDS_PER_ATOM(B)  !NCONPERATOM(B)
               A = BONDS_PER_ATOM_LIST(B,J2) !CONLIST(B, J2)
               IF (.NOT. LINATOM(A)) CYCLE

               DO J3 = J2+1, N_BONDS_PER_ATOM(B)  !NCONPERATOM(B)
                  C = BONDS_PER_ATOM_LIST(B,J3) !CONLIST(B, J3)
                  IF (.NOT. LINATOM(C)) CYCLE

                  ! Look for D bonded to C, excluding B
                  DO J4 = 1,  N_BONDS_PER_ATOM(C) !NCONPERATOM(C)
                     D = BONDS_PER_ATOM_LIST(C, J4)  !CONLIST(C, J4)
                     IF (D == B) CYCLE
                     IF (.NOT. LINATOM(D)) CYCLE

                     ! --- START geometry ---
                     DIH_COORDS(1:3)   = XSTART(3*(A-1)+1:3*(A-1)+3)
                     DIH_COORDS(4:6)   = XSTART(3*(B-1)+1:3*(B-1)+3)
                     DIH_COORDS(7:9)   = XSTART(3*(C-1)+1:3*(C-1)+3)
                     DIH_COORDS(10:12) = XSTART(3*(D-1)+1:3*(D-1)+3)

                     DIH_START = DIHEDRAL(DIH_COORDS)

                     ! --- FINAL geometry ---
                     DIH_COORDS(1:3)   = XFINAL(3*(A-1)+1:3*(A-1)+3)
                     DIH_COORDS(4:6)   = XFINAL(3*(B-1)+1:3*(B-1)+3)
                     DIH_COORDS(7:9)   = XFINAL(3*(C-1)+1:3*(C-1)+3)
                     DIH_COORDS(10:12) = XFINAL(3*(D-1)+1:3*(D-1)+3)

                     DIH_FINAL = DIHEDRAL(DIH_COORDS)

                     DIHEDRAL_DEVIATION = DIHEDRAL_DIFF(DIH_START,DIH_FINAL)
                     

                     IF (DIHEDRAL_DEVIATION.GT.DIHEDRAL_TOLERANCE) THEN
                        IS_RIGID_BODY(B) = .FALSE.
                        EXIT
                     END IF
                  END DO
                  IF (.NOT. IS_RIGID_BODY(B)) EXIT
               END DO
               IF (.NOT. IS_RIGID_BODY(B)) EXIT
            END DO
         END DO


         ! Update LINATOM to only include truly rigid atoms
         LINATOM(:) = LINATOM(:) .AND. IS_RIGID_BODY(:)

   END SUBROUTINE FIND_LINEAR_ATOMS
   
   
   SUBROUTINE COMPARE_AMBER_LINEAR_GROUPS()
      USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET       
            USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
            USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS, ANGLE, DIHEDRAL, DIHEDRAL_DIFF
            USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL, NCONPERATOM, &
                                          CONLIST, BOND_LIST, NBONDS, MAX_BONDS_PER_ATOM, &
                                          BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
            USE QCICONSTRAINTS, ONLY: GET_NBONDS_PER_ATOM
            USE AMBER_CONSTRAINTS, ONLY: INGROUP, GROUPLOOKUP, NPLACINGGROUPS, SIZEPLACINGGROUPS
      IMPLICIT NONE

      INTEGER :: AMBER_ATOM, AMBER_GROUP
      INTEGER :: NATOMS_THIS_GROUP
      INTEGER :: ATOMS_IN_BOTH(NATOMS)
      INTEGER :: ATOMS_ONLY_IN_AMBER(NATOMS)
      INTEGER :: N_IN_BOTH, N_ONLY_AMBER
      INTEGER :: J1, K

      PRINT *, ""
      PRINT *, "========== AMBER vs LINEAR GROUP COMPARISON =========="

      ! Handle case with no linear groups
      IF (NLINGROUPS == 0) THEN
         PRINT *, "No linear groups detected."
         DO AMBER_GROUP = 1, NPLACINGGROUPS
               NATOMS_THIS_GROUP = SIZEPLACINGGROUPS(AMBER_GROUP)
               PRINT *, ""
               PRINT *, "Amber Group ", AMBER_GROUP, " (", NATOMS_THIS_GROUP, " atoms)"
               PRINT *, "  Atoms ALSO in linear groups: NONE (no linear groups exist)"
               PRINT *, "  Atoms ONLY in Amber (not in linear groups): ALL"
               PRINT *, "  -> NO OVERLAP: No linear groups to compare"
         END DO
         PRINT *, "========================================================"
         RETURN
      END IF

      ! Handle case with no Amber groups
      IF (NPLACINGGROUPS == 0) THEN
         PRINT *, "No Amber placing groups defined."
         PRINT *, "========================================================"
         RETURN
      END IF

      ! Loop over each Amber placing group
      DO AMBER_GROUP = 1, NPLACINGGROUPS
         
         NATOMS_THIS_GROUP = SIZEPLACINGGROUPS(AMBER_GROUP)
         N_IN_BOTH = 0
         N_ONLY_AMBER = 0
         ATOMS_IN_BOTH(:) = 0
         ATOMS_ONLY_IN_AMBER(:) = 0
         
         ! Scan all atoms to find those belonging to this Amber group
         DO J1 = 1, NATOMS
               IF (GROUPLOOKUP(J1) == AMBER_GROUP) THEN
                  AMBER_ATOM = J1
                  
                  ! Check if this atom is in a linear group
                  IF (ATOM2LINGROUP(AMBER_ATOM) > 0) THEN
                     N_IN_BOTH = N_IN_BOTH + 1
                     ATOMS_IN_BOTH(N_IN_BOTH) = AMBER_ATOM
                  ELSE
                     N_ONLY_AMBER = N_ONLY_AMBER + 1
                     ATOMS_ONLY_IN_AMBER(N_ONLY_AMBER) = AMBER_ATOM
                  END IF
               END IF
         END DO
         
         ! Print comparison for this Amber group
         PRINT *, ""
         PRINT *, "Amber Group ", AMBER_GROUP, " (", NATOMS_THIS_GROUP, " atoms)"
         
         IF (N_IN_BOTH > 0) THEN
               PRINT *, "  Atoms ALSO in linear groups (", N_IN_BOTH, "): ", &
                        (ATOMS_IN_BOTH(K), K=1,N_IN_BOTH)
         ELSE
               PRINT *, "  Atoms ALSO in linear groups: NONE"
         END IF
         
         IF (N_ONLY_AMBER > 0) THEN
               PRINT *, "  Atoms ONLY in Amber (not in linear groups) (", N_ONLY_AMBER, "): ", &
                        (ATOMS_ONLY_IN_AMBER(K), K=1,N_ONLY_AMBER)
         ELSE
               PRINT *, "  Atoms ONLY in Amber: NONE"
         END IF
         
         ! Summary
         IF (N_ONLY_AMBER > 0 .AND. N_IN_BOTH > 0) THEN
               PRINT *, "  -> MIXED: Some atoms in linear groups, some not"
         ELSE IF (N_ONLY_AMBER == 0 .AND. N_IN_BOTH == NATOMS_THIS_GROUP) THEN
               PRINT *, "  -> FULL OVERLAP: All Amber atoms in linear groups"
         ELSE IF (N_ONLY_AMBER == NATOMS_THIS_GROUP) THEN
               PRINT *, "  -> NO OVERLAP: No Amber atoms in linear groups"
         END IF

      END DO

      PRINT *, "========================================================"

   END SUBROUTINE COMPARE_AMBER_LINEAR_GROUPS


   !> Allocate and save bonds
   SUBROUTINE SAVE_BONDS(BONDS_IN, NBOND_IN)
      
      USE QCI_CONSTRAINT_KEYS, ONLY: NBONDS, BOND_LIST
      IMPLICIT NONE

      ! === Inputs ===
      INTEGER, INTENT(IN) :: NBOND_IN
      INTEGER, INTENT(IN) :: BONDS_IN(NBOND_IN,2)

      ! === Local ===
      INTEGER :: I

      ! === Sanity checks ===
      IF (NBOND_IN <= 0) THEN
         PRINT *, "SAVE_BONDS: nbonds <= 0, nothing to store"
         RETURN
      END IF

      ! === Deallocate old storage if present ===
      IF (ALLOCATED(BOND_LIST)) THEN
         DEALLOCATE(BOND_LIST)
      END IF

      ! === Allocate new storage ===
      ALLOCATE(BOND_LIST(2, NBOND_IN))
      NBONDS = NBOND_IN

      DO I=1, NBONDS
         BOND_LIST( 1, I) = BONDS_IN(I,1)
         BOND_LIST( 2, I) = BONDS_IN(I,2)
      END DO

   END SUBROUTINE SAVE_BONDS
  

   !> Currently deallocation done in DEALLOC_LINEAR_GROUP, so this doesn't need to be called
   SUBROUTINE DEALLOC_BONDS()
      USE QCI_CONSTRAINT_KEYS, ONLY: NBONDS, BOND_LIST
      IF(ALLOCATED(BOND_LIST)) DEALLOCATE(BOND_LIST)
   END SUBROUTINE DEALLOC_BONDS


   SUBROUTINE REDUCE_MULTI_HINGE_GROUP(LINEARATOMS, GROUP_ATOMS, N_GROUP_ATOMS, SUCCESS, CURRENT_GROUP,DEBUG_MODE)

      USE QCIKEYS, ONLY: NATOMS 
      USE QCI_CONSTRAINT_KEYS, ONLY: BOND_LIST, NBONDS, MAX_BONDS_PER_ATOM, &
                                       BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
      USE QCICONSTRAINTS, ONLY: GET_NBONDS_PER_ATOM

      IMPLICIT NONE

      ! Input/Output
      LOGICAL, INTENT(INOUT) :: LINEARATOMS(:)
      INTEGER, INTENT(INOUT) :: GROUP_ATOMS(:)
      INTEGER, INTENT(INOUT) :: N_GROUP_ATOMS
      INTEGER, INTENT(INOUT) :: CURRENT_GROUP(:)
      LOGICAL, INTENT(OUT)   :: SUCCESS
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_MODE


      ! Local variables
      INTEGER :: WORKING_ATOMS(NATOMS)
      INTEGER :: N_WORKING
      INTEGER :: EXTERNAL_ANCHORS(NATOMS)
      INTEGER :: N_EXTERNAL_ANCHORS
      INTEGER :: INTERNAL_DEGREE(NATOMS)
      INTEGER :: ATOM, NEIGHBOR, J1, J2, J3
      INTEGER :: ATOM_TO_REMOVE
      INTEGER :: MIN_INTERNAL
      INTEGER :: MIN_SIZE
      LOGICAL :: IS_VALID
      LOGICAL :: LOCAL_DEBUG

      ! Parameters
      INTEGER, PARAMETER :: MIN_GROUP_SIZE = 3  ! Don't reduce below this

      ! Check debug mode
      LOCAL_DEBUG = .TRUE.
      IF (PRESENT(DEBUG_MODE)) LOCAL_DEBUG = DEBUG_MODE

      ! Initialize
      SUCCESS = .FALSE.
      
      ! Copy initial group to working array
      WORKING_ATOMS(1:N_GROUP_ATOMS) = GROUP_ATOMS(1:N_GROUP_ATOMS)
      N_WORKING = N_GROUP_ATOMS

      IF (LOCAL_DEBUG) THEN
         PRINT *, ""
         PRINT *, "=== Starting reduction for group with ", N_WORKING, " atoms ==="
         PRINT *, "Initial atoms: ", (WORKING_ATOMS(J1), J1=1,N_WORKING)
      END IF

      ! ==========================================================================
      ! Main reduction loop
      ! ==========================================================================
      REDUCTION_LOOP: DO WHILE (N_WORKING.GT.MIN_GROUP_SIZE)

         ! Step 1: Find external anchors for current working group
         CALL FIND_EXTERNAL_ANCHORS(WORKING_ATOMS, N_WORKING, &
                                    EXTERNAL_ANCHORS, N_EXTERNAL_ANCHORS, CURRENT_GROUP)

         IF (LOCAL_DEBUG) THEN
            PRINT *, ""
            PRINT *, "Current working group size: ", N_WORKING
            PRINT *, "Number of external anchors: ", N_EXTERNAL_ANCHORS
            PRINT *, "External anchors: ", (EXTERNAL_ANCHORS(J1), J1=1,N_EXTERNAL_ANCHORS)
         END IF

         ! Step 2: Check if we have valid single-hinge group
         IF (N_EXTERNAL_ANCHORS.EQ.1) THEN
            SUCCESS = .TRUE.
            
            N_GROUP_ATOMS = N_WORKING
            GROUP_ATOMS(1:N_WORKING) = WORKING_ATOMS(1:N_WORKING)
            
            IF (LOCAL_DEBUG) THEN
               PRINT *, ""
               PRINT *, "*** SUCCESS: Reduced to single-hinge group ***"
               PRINT *, "Final reduced atoms: ", (GROUP_ATOMS(J1), J1=1,N_GROUP_ATOMS)
               PRINT *, "Hinge atom: ", EXTERNAL_ANCHORS(1)
            END IF
            
            RETURN
         END IF

         ! Step 3: Check if we've hit minimum size
         IF (N_WORKING.LE.MIN_GROUP_SIZE) THEN
            IF (LOCAL_DEBUG) THEN
               PRINT *, ""
               PRINT *, "*** FAILED: Reached minimum group size ", MIN_GROUP_SIZE, " ***"
            END IF
            EXIT REDUCTION_LOOP
         END IF

         ! Step 4: Calculate internal degree for each atom in working group
         INTERNAL_DEGREE(:) = 0
         
         DO J1 = 1, N_WORKING
            ATOM = WORKING_ATOMS(J1)
            
            ! Count connections to other atoms in the working group
            DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
               NEIGHBOR = BONDS_PER_ATOM_LIST(ATOM, J2)
               
               ! Check if neighbor is also in the working group
               IF (IS_IN_LIST(NEIGHBOR, WORKING_ATOMS, N_WORKING)) THEN
                  INTERNAL_DEGREE(J1) = INTERNAL_DEGREE(J1) + 1
               END IF
            END DO
         END DO

         IF (LOCAL_DEBUG) THEN
            PRINT *, "Internal degrees:"
            DO J1 = 1, N_WORKING
               ATOM = WORKING_ATOMS(J1)
               PRINT *, "  Atom ", ATOM, ": ", INTERNAL_DEGREE(J1), " internal connections"
            END DO
         END IF

         ! Step 5: Find atom with minimum internal degree to remove
         MIN_INTERNAL = 100
         ATOM_TO_REMOVE = -1
         

         DO J1 = 1, N_WORKING
            
            !If we find an atom with zero internal connections, remove it immediately (can't disconnect the group)
            IF (INTERNAL_DEGREE(J1).EQ.0) THEN
               ATOM_TO_REMOVE = J1
               MIN_INTERNAL = INTERNAL_DEGREE(J1)
               EXIT
            END IF
            !If the atom has minimum number of internal connections, it is candidate for removal. 
            ATOM = WORKING_ATOMS(J1)
            
            DO J2 = 1, N_EXTERNAL_ANCHORS
            ! If atom is connected to external anchor, it is candidate for removal
               
               IF(IS_IN_LIST(EXTERNAL_ANCHORS(J2), BONDS_PER_ATOM_LIST(ATOM,:), N_BONDS_PER_ATOM(ATOM))) THEN
                  WRITE(*,*) "reduce_multi_hinge_group> Atom ", WORKING_ATOMS(J1), " is connected to external anchor ", EXTERNAL_ANCHORS(J2), &
                                 " with internal degree ", INTERNAL_DEGREE(J1)
                  WRITE(*,*) "J1: ", ATOM, "BONDS_PER_ATOM_LIST(J1): ", BONDS_PER_ATOM_LIST(ATOM,:)
                  IF (INTERNAL_DEGREE(J1).LT.MIN_INTERNAL) THEN
                     ATOM_TO_REMOVE = J1
                     MIN_INTERNAL = INTERNAL_DEGREE(J1) !INTERNAL_DEGREE(ATOM)
                  END iF
               END IF
            END DO
           
         END DO

         IF (LOCAL_DEBUG) THEN
            PRINT *, "Removing atom ",WORKING_ATOMS(ATOM_TO_REMOVE), &
                     " with ", MIN_INTERNAL, " internal connections"
         END IF

         LINEARATOMS(WORKING_ATOMS(ATOM_TO_REMOVE)) = .FALSE.
         CURRENT_GROUP(WORKING_ATOMS(ATOM_TO_REMOVE)) = -1
         ! Step 6: Remove atom from working group (shift remaining atoms)
         DO J1 = ATOM_TO_REMOVE, N_WORKING - 1
            WORKING_ATOMS(J1) = WORKING_ATOMS(J1 + 1)
         END DO
         N_WORKING = N_WORKING - 1
         

      END DO REDUCTION_LOOP

      ! ==========================================================================
      ! Reduction failed - return what we have and mark as invalid)
      ! ==========================================================================
      N_GROUP_ATOMS = N_WORKING
      GROUP_ATOMS(1:N_WORKING) = WORKING_ATOMS(1:N_WORKING)
      SUCCESS = .FALSE.

      IF (LOCAL_DEBUG) THEN
         PRINT *, ""
         PRINT *, "=== Reduction complete - FAILED ==="
         PRINT *, "Final atoms: ", (GROUP_ATOMS(J1), J1=1,N_GROUP_ATOMS)
      END IF

   END SUBROUTINE REDUCE_MULTI_HINGE_GROUP


   SUBROUTINE FIND_EXTERNAL_ANCHORS(GROUP_ATOMS, N_GROUP_ATOMS, EXTERNAL_ANCHORS, N_EXTERNAL_ANCHORS, CURRENT_GROUP)
      
      USE QCI_CONSTRAINT_KEYS, ONLY: BOND_LIST, NBONDS, MAX_BONDS_PER_ATOM, &
                                          BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM

      IMPLICIT NONE
      
      INTEGER, INTENT(IN)  :: GROUP_ATOMS(:)
      INTEGER, INTENT(IN)  :: N_GROUP_ATOMS
      INTEGER, INTENT(OUT) :: EXTERNAL_ANCHORS(:)
      INTEGER, INTENT(OUT) :: N_EXTERNAL_ANCHORS
      INTEGER, INTENT(IN)  :: CURRENT_GROUP(:)
      INTEGER :: ATOM, NEIGHBOR, J1, J2

      EXTERNAL_ANCHORS(:) = -1
      N_EXTERNAL_ANCHORS  = 0

      DO J1 = 1, N_GROUP_ATOMS
         ATOM = GROUP_ATOMS(J1)

         DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
            NEIGHBOR = BONDS_PER_ATOM_LIST(ATOM, J2)

            ! Check if neighbor belongs to a different group
            IF (CURRENT_GROUP(NEIGHBOR).NE.CURRENT_GROUP(ATOM)) THEN
               ! This is an external anchor candidate
               IF (.NOT. IS_IN_LIST(NEIGHBOR, EXTERNAL_ANCHORS, N_EXTERNAL_ANCHORS)) THEN
                  N_EXTERNAL_ANCHORS = N_EXTERNAL_ANCHORS + 1
                  EXTERNAL_ANCHORS(N_EXTERNAL_ANCHORS) = NEIGHBOR
               END IF
            END IF
         END DO !J2
      END DO !J1

   END SUBROUTINE FIND_EXTERNAL_ANCHORS

   !> helper function to check if a value (integer) is in a list of length N
   LOGICAL FUNCTION IS_IN_LIST(VALUE, LIST, N)
      INTEGER, INTENT(IN) :: VALUE
      INTEGER, INTENT(IN) :: LIST(:)
      INTEGER, INTENT(IN) :: N
      INTEGER :: I

      IS_IN_LIST = .FALSE.
      DO I = 1, N
         IF (LIST(I).EQ.VALUE) THEN
            IS_IN_LIST = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION IS_IN_LIST

END MODULE QCI_LINEAR
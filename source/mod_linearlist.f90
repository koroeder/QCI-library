!! @TODO CONTINUE here. 
!! External anchors does not give correct results 

MODULE QCI_LINEAR
   USE QCIPREC
   IMPLICIT NONE

   !------------------old variables-----------------------------------!
   INTEGER :: NQCILINEAR = 0 !< number of atoms for QCIlinear 
   INTEGER, ALLOCATABLE :: LINEARATOMS(:) !< list of linear atoms
   CHARACTER(25) :: LINEARFILE = "QCIlinear" !< file name linear atoms
    REAL(KIND=REAL64) :: LINEARCUT = 0.05D0 !< user defined value, cutoff for linear atoms 
   !------------------------------------------------------------------!


   INTEGER, ALLOCATABLE :: LINEAR_GROUPS(:,:)  !< list of linear groups
   INTEGER, ALLOCATABLE :: ATOM2LINGROUP(:)
   INTEGER, ALLOCATABLE :: NINGROUP(:)         !< Number of atoms in each group
   INTEGER :: NLINGROUPS = 0                   !< Number of linear groups
   
  
   REAL(KIND=REAL64), PARAMETER :: TOLERANCE = 0.1D0               !< bond length change tol between start and finish 
   REAL(KIND=REAL64), PARAMETER :: ANGLE_TOLERANCE = 0.08726646259971647     !< 5.0D0 degrees  
   REAL(KIND=REAL64), PARAMETER :: DIHEDRAL_TOLERANCE = 0.17453292519943295  !< 10.0D0 degrees
   INTEGER, PARAMETER :: MIN_LINEAR_GROUP_SIZE = 5
   
   CONTAINS

      !> Get linear atoms from a file. 
      !! TO BE REMOVED
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

         !DO J1=1,NATOMS
         !   CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, XSTART, XFINAL, J1, DIST)
         !   IF (DIST.LT.LINEARCUT) THEN
         !      LINEART(J1) = 1
         !   END IF
         !END DO

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
      
      !> Old function. For removal 
      SUBROUTINE ALLOC_QCI_LINEAR()
         USE QCIKEYS, ONLY: NATOMS, INLINLIST
         CALL DEALLOC_QCI_LINEAR
         ALLOCATE(LINEARATOMS(NQCILINEAR))
         ALLOCATE(INLINLIST(NATOMS))
      END SUBROUTINE ALLOC_QCI_LINEAR

      !> Old function. For removal 
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

         REAL(KIND=REAL64), PARAMETER :: TOLERANCE = 0.01D0
         REAL(KIND=REAL64) :: DS, DF
         LOGICAL :: LINATOM(NATOMS) 
         INTEGER :: NLINHERE
         INTEGER :: A, B, C, D, J1, J2, J3, J4, K

         INTEGER :: CURRENTGROUP(NATOMS)      !< Group ID for each linear atom, -1 if not assigned
         INTEGER :: GROUPS(NATOMS, NATOMS)    !< List of atoms in each group
         INTEGER :: GROUPID                   !< Total number of groups
         INTEGER :: VISITED(NATOMS)           !< For DFS traversal
         INTEGER :: STACK(NATOMS)             !< Stack for DFS
         INTEGER :: STACKPTR
         INTEGER :: ATOM, NEIGHBOR
         INTEGER :: NINGROUP_TEMP(3*NATOMS)
         INTEGER :: MAXGROUPSIZE
         INTEGER :: NGROUPS                     !< number of linear groups
                
        
         INTEGER :: EXTERNAL_ANCHORS(NATOMS)  !< distinct external atoms bonded to this group 
         INTEGER :: N_EXTERNAL_ANCHORS        !< number of distinct external anchor atoms
         INTEGER :: HINGE_ATOM                !< unique external anchor (if exists)

  
         LOGICAL :: DEBUG_MODE = .FALSE.       !< Set to .TRUE. for verbose output
         
         INTEGER :: EXTERNAL_PER_ATOM(NATOMS)
         INTEGER :: N_BACKBONE_ATOMS
         INTEGER :: BACKBONE_ATOM
         LOGICAL :: IS_VALID_GROUP
         INTEGER :: TEMP_EXTERNAL_COUNT
         INTEGER :: NEIGHBOR_ATOM
         INTEGER :: N_NEW_GROUPS
      
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

         !Part 3 - Validate groups

      GROUPID = 0
      DO J1 = 1, NGROUPS
         IF (NINGROUP_TEMP(J1) <= 1) CYCLE
            CALL REDUCE_MULTI_HINGE_GROUP(LINATOM, GROUPS(J1,1:NINGROUP_TEMP(J1)), NINGROUP_TEMP(J1), &
                                          GROUPID, N_NEW_GROUPS, DEBUG_MODE)
            END DO
      NLINGROUPS = GROUPID
         NQCILINEAR = SUM(NINGROUP(1:NLINGROUPS))
      

   END SUBROUTINE GET_LINEAR_GROUPS

   
   SUBROUTINE GET_LIN_ROT_TRANSLATION(GROUP, GROUPSIZE, CXS, CXF, QSTART, QFINAL)
      USE QCIKEYS, ONLY: NATOMS
      USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
      USE align, ONLY: move_to_origin
      use lpermdist, only: get_rot_matrix
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
      call move_to_origin(GROUPSIZE,RS,CXS)
      !centre RB around origin
      CALL move_to_origin(GROUPSIZE,RF,CXF)


      !WRITE(*,*) "get_lin_rot_translation> GROUPSIZE ", GROUPSIZE, "RS ", RS, " RF ", RF
      !align coordinates
      call get_rot_matrix(GROUPSIZE, RF, RS, RMAT, DIST)

      QSTART = [1.0D0, 0.0D0, 0.0D0, 0.0D0]

      CALL MATRIX_TO_QUATERNION(RMAT,QFINAL)
      
   END SUBROUTINE GET_LIN_ROT_TRANSLATION

   !> Linear atoms here are defined as atoms whose bonds, angles and dihedrals change less than 
   !! the tolarence between start & finish.  
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

         
         REAL(KIND=REAL64) :: DS, DF
         
         INTEGER :: NLINHERE
         INTEGER :: A, B, C, D, J1, J2, J3, J4, K
    
          ! Grouping variables
         INTEGER :: CURRENTGROUP(NATOMS)      ! Group ID for each linear atom
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

         REAL(KIND=REAL64) :: DIH_COORDS(12)
         REAL(KIND=REAL64) :: DIH_START, DIH_FINAL
         
         LINATOM(:) = .FALSE.
         IS_RIGID_BODY(:) = .TRUE.
      
          ! ========== PART 1: Detect linear atoms ==========
         DO J1 = 1, NATOMS
            NLINHERE = 0
         
            DO J2 = 1, N_BONDS_PER_ATOM(J1) !NCONPERATOM(J1)
               A = J1
               B = BONDS_PER_ATOM_LIST(A,J2) !CONLIST(A, J2)

               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, A, B, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, A, B, DF)

               IF(DABS(DS - DF).LT.TOLERANCE) NLINHERE = NLINHERE + 1
            END DO

            ! Mark atom as linear if all constraints have unchanged distances
            IF(NLINHERE == N_BONDS_PER_ATOM(J1)) LINATOM(J1) = .TRUE.
            
          END DO
    
          
         ! Check if atoms in a group maintain angles
         DO J1 = 1, NATOMS
            IF (.NOT. LINATOM(J1)) CYCLE
            
            ! For each atom, check angles to neighbors
            DO J2 = 1, N_BONDS_PER_ATOM(J1) 
               DO J3 = J2+1, N_BONDS_PER_ATOM(J1) 
                    
                     A = BONDS_PER_ATOM_LIST(J1, J2)
                     B = J1
                     C = BONDS_PER_ATOM_LIST(J1, J3)
                                          
                     ANGLE_COORDS(1:3) = XSTART(3*(A-1)+1:3*(A-1)+3)
                     ANGLE_COORDS(4:6) = XSTART(3*(B-1)+1:3*(B-1)+3)
                     ANGLE_COORDS(7:9) = XSTART(3*(C-1)+1:3*(C-1)+3)

                     ANGLE_START = ANGLE(ANGLE_COORDS)
                     
                     ANGLE_COORDS(1:3) = XFINAL(3*(A-1)+1:3*(A-1)+3)
                     ANGLE_COORDS(4:6) = XFINAL(3*(B-1)+1:3*(B-1)+3)
                     ANGLE_COORDS(7:9) = XFINAL(3*(C-1)+1:3*(C-1)+3)


                     ANGLE_FINAL = ANGLE(ANGLE_COORDS)
                     
                     ANGLE_DEVIATION = DABS(ANGLE_START - ANGLE_FINAL)
                     
                     IF (ANGLE_DEVIATION.GT.ANGLE_TOLERANCE) THEN
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

            DO J2 = 1, N_BONDS_PER_ATOM(B)  
               A = BONDS_PER_ATOM_LIST(B,J2) 
               IF (.NOT. LINATOM(A)) CYCLE

               DO J3 = J2+1, N_BONDS_PER_ATOM(B)  
                  C = BONDS_PER_ATOM_LIST(B,J3) 
                  IF (.NOT. LINATOM(C)) CYCLE

                  ! Look for D bonded to C, excluding B
                  DO J4 = 1,  N_BONDS_PER_ATOM(C) 
                     D = BONDS_PER_ATOM_LIST(C, J4) 
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
      IF (NLINGROUPS.EQ.0) THEN
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
      IF (NPLACINGGROUPS.EQ.0) THEN
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
               IF (GROUPLOOKUP(J1).EQ.AMBER_GROUP) THEN
                  AMBER_ATOM = J1
                  ! Check if this atom is in a linear group
                  IF (ATOM2LINGROUP(AMBER_ATOM).GT.0) THEN
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
         
         IF (N_IN_BOTH.GT.0) THEN
               PRINT *, "  Atoms ALSO in linear groups (", N_IN_BOTH, "): ", &
                        (ATOMS_IN_BOTH(K), K=1,N_IN_BOTH)
         ELSE
               PRINT *, "  Atoms ALSO in linear groups: NONE"
         END IF
         
         IF (N_ONLY_AMBER.GT.0) THEN
               PRINT *, "  Atoms ONLY in Amber (not in linear groups) (", N_ONLY_AMBER, "): ", &
                        (ATOMS_ONLY_IN_AMBER(K), K=1,N_ONLY_AMBER)
         ELSE
               PRINT *, "  Atoms ONLY in Amber: NONE"
         END IF
         
         ! Summary
         IF ( (N_ONLY_AMBER.GT.0) .AND. (N_IN_BOTH.GT.0) ) THEN
               PRINT *, "  -> MIXED: Some atoms in linear groups, some not"
         ELSE IF ( (N_ONLY_AMBER.EQ.0) .AND. (N_IN_BOTH.EQ.NATOMS_THIS_GROUP)) THEN
               PRINT *, "  -> FULL OVERLAP: All Amber atoms in linear groups"
         ELSE IF (N_ONLY_AMBER == NATOMS_THIS_GROUP) THEN
               PRINT *, "  -> NO OVERLAP: No Amber atoms in linear groups"
         END IF

      END DO

      PRINT *, "========================================================"

   END SUBROUTINE COMPARE_AMBER_LINEAR_GROUPS



  

   !> Currently deallocation done in DEALLOC_LINEAR_GROUP, so this doesn't need to be called
   SUBROUTINE DEALLOC_BONDS()
      USE QCI_CONSTRAINT_KEYS, ONLY: NBONDS, BOND_LIST
      IF(ALLOCATED(BOND_LIST)) DEALLOCATE(BOND_LIST)
   END SUBROUTINE DEALLOC_BONDS


  

   

   !====================================================================
   ! Resolves a connected group of linear atoms that touches more than
   ! one external anchor into one or more valid sub-groups, by
   ! repeatedly CONSTRUCTING the largest possible single-anchor group
   ! from the current pool (growing outward from each candidate anchor
   ! and keeping the best), committing it, and recursing on whatever
   ! atoms are left. Atoms that end up in no group of at least
   ! MIN_LINEAR_GROUP_SIZE atoms are dropped (LINATOM_MASK cleared).
   !====================================================================
   SUBROUTINE REDUCE_MULTI_HINGE_GROUP(LINATOM_MASK, GROUP_ATOMS_IN, N_GROUP_ATOMS_IN, &
                                        GROUPID, N_NEW_GROUPS, DEBUG_MODE)
      USE QCIKEYS, ONLY: NATOMS 
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: LINATOM_MASK(:)
      INTEGER, INTENT(IN)    :: GROUP_ATOMS_IN(:)
      INTEGER, INTENT(IN)    :: N_GROUP_ATOMS_IN
      INTEGER, INTENT(INOUT) :: GROUPID
      INTEGER, INTENT(OUT)   :: N_NEW_GROUPS
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_MODE

      ! explicit work-list of pending pools (stack); worst case one
      ! pool per input atom, e.g. a star graph shattering into singles
      INTEGER :: POOL_STACK(N_GROUP_ATOMS_IN, N_GROUP_ATOMS_IN)
      INTEGER :: POOL_STACK_SIZE(N_GROUP_ATOMS_IN)
      INTEGER :: NSTACK

      INTEGER :: CUR_ATOMS(N_GROUP_ATOMS_IN), CUR_SIZE
      LOGICAL :: IN_SET(NATOMS)
      INTEGER :: ANCHORS(N_GROUP_ATOMS_IN), N_ANCHORS
      INTEGER :: TRY_ATOMS(N_GROUP_ATOMS_IN), TRY_SIZE
      INTEGER :: BEST_ATOMS(N_GROUP_ATOMS_IN), BEST_SIZE
      INTEGER :: REMAINDER(N_GROUP_ATOMS_IN), REM_SIZE
      INTEGER :: COMP_ATOMS(N_GROUP_ATOMS_IN, N_GROUP_ATOMS_IN)
      INTEGER :: COMP_SIZES(N_GROUP_ATOMS_IN)
      INTEGER :: NCOMP
      LOGICAL :: SUCCESS, LOCAL_DEBUG
      INTEGER :: J1

      LOCAL_DEBUG = .FALSE.
      IF (PRESENT(DEBUG_MODE)) LOCAL_DEBUG = DEBUG_MODE

      N_NEW_GROUPS = 0
      NSTACK = 1
      POOL_STACK_SIZE(1) = N_GROUP_ATOMS_IN
      POOL_STACK(1, 1:N_GROUP_ATOMS_IN) = GROUP_ATOMS_IN(1:N_GROUP_ATOMS_IN)

      IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> starting with ", N_GROUP_ATOMS_IN, " atoms"

      DO WHILE (NSTACK > 0)

         CUR_SIZE = POOL_STACK_SIZE(NSTACK)
         CUR_ATOMS(1:CUR_SIZE) = POOL_STACK(NSTACK, 1:CUR_SIZE)
         NSTACK = NSTACK - 1

         IF (CUR_SIZE < MIN_LINEAR_GROUP_SIZE) THEN
            IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> pool of ", CUR_SIZE, " too small, dropping"
            CALL RMH_DROP_ATOMS(LINATOM_MASK, CUR_ATOMS, CUR_SIZE)
            CYCLE
      END IF

         CALL RMH_MEMBERSHIP(CUR_ATOMS, CUR_SIZE, IN_SET)
         CALL RMH_COMPUTE_ANCHORS(CUR_ATOMS, CUR_SIZE, IN_SET, ANCHORS, N_ANCHORS)

         IF (N_ANCHORS == 1) THEN
            CALL RMH_COMMIT_GROUP(CUR_ATOMS, CUR_SIZE, LINATOM_MASK, GROUPID)
            N_NEW_GROUPS = N_NEW_GROUPS + 1
            IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> committed group of size ", &
                                        CUR_SIZE, ", hinge atom ", ANCHORS(1)
            CYCLE
         END IF

         IF (N_ANCHORS == 0) THEN
            IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> dropping isolated fragment of size ", CUR_SIZE
            CALL RMH_DROP_ATOMS(LINATOM_MASK, CUR_ATOMS, CUR_SIZE)
            CYCLE
         END IF

         ! N_ANCHORS >= 2: construct the largest single-anchor group by
         ! growing from each candidate anchor and keeping the best
         BEST_SIZE = 0
         DO J1 = 1, N_ANCHORS
            CALL RMH_GROW_FROM_ANCHOR(CUR_ATOMS, CUR_SIZE, IN_SET, ANCHORS(J1), TRY_ATOMS, TRY_SIZE)
            IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> anchor ", ANCHORS(J1), &
                                        " grows to ", TRY_SIZE, " atoms"
            IF (TRY_SIZE > BEST_SIZE) THEN
               BEST_SIZE = TRY_SIZE
               BEST_ATOMS(1:TRY_SIZE) = TRY_ATOMS(1:TRY_SIZE)
            END IF
         END DO

         IF (BEST_SIZE >= MIN_LINEAR_GROUP_SIZE) THEN

            CALL RMH_COMMIT_GROUP(BEST_ATOMS, BEST_SIZE, LINATOM_MASK, GROUPID)
            N_NEW_GROUPS = N_NEW_GROUPS + 1
            IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> constructed group of size ", BEST_SIZE

            ! what's left may fall apart into several disconnected
            ! pieces once BEST_ATOMS is removed -- resolve each
            ! independently
            CALL RMH_REMOVE_ATOMS(CUR_ATOMS, CUR_SIZE, BEST_ATOMS, BEST_SIZE, REMAINDER, REM_SIZE)

            IF (REM_SIZE > 0) THEN
               CALL RMH_CONNECTED_COMPONENTS(REMAINDER, REM_SIZE, COMP_ATOMS, COMP_SIZES, NCOMP)
               DO J1 = 1, NCOMP
                  NSTACK = NSTACK + 1
                  POOL_STACK_SIZE(NSTACK) = COMP_SIZES(J1)
                  POOL_STACK(NSTACK, 1:COMP_SIZES(J1)) = COMP_ATOMS(J1, 1:COMP_SIZES(J1))
               END DO
         END IF

         ELSE
            ! construction couldn't reach the minimum from any single
            ! anchor -- last-resort atom-by-atom trim before giving up
            CALL RMH_GREEDY_TRIM(CUR_ATOMS, CUR_SIZE, MIN_LINEAR_GROUP_SIZE, SUCCESS, LOCAL_DEBUG)
            IF (SUCCESS) THEN
               CALL RMH_COMMIT_GROUP(CUR_ATOMS, CUR_SIZE, LINATOM_MASK, GROUPID)
               N_NEW_GROUPS = N_NEW_GROUPS + 1
               IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> fallback trim salvaged group of size ", CUR_SIZE
            ELSE
               IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> pool unsalvageable, dropping ", CUR_SIZE, " atoms"
               CALL RMH_DROP_ATOMS(LINATOM_MASK, CUR_ATOMS, CUR_SIZE)
            END IF
         END IF

      END DO

      IF (LOCAL_DEBUG) PRINT *, "reduce_multi_hinge_group> input of ", N_GROUP_ATOMS_IN, &
                                  " atoms resolved into ", N_NEW_GROUPS, " valid sub-group(s)"

   END SUBROUTINE REDUCE_MULTI_HINGE_GROUP


   !--------------------------------------------------------------------
   ! Private helpers (prefixed RMH_)
   !--------------------------------------------------------------------

   !> IN_SET(atom) = .TRUE. iff atom appears in ATOMS(1:NAT)
   SUBROUTINE RMH_MEMBERSHIP(ATOMS, NAT, IN_SET)
      USE QCIKEYS, ONLY: NATOMS
      INTEGER, INTENT(IN)  :: ATOMS(:), NAT
      LOGICAL, INTENT(OUT) :: IN_SET(NATOMS)
      INTEGER :: J1
      IN_SET(:) = .FALSE.
      DO J1 = 1, NAT
         IN_SET(ATOMS(J1)) = .TRUE.
      END DO
   END SUBROUTINE RMH_MEMBERSHIP

   !> Distinct atoms bonded into ATOMS(1:NAT) from outside the set
   !! (checked against the full bond graph).
   SUBROUTINE RMH_COMPUTE_ANCHORS(ATOMS, NAT, IN_SET, ANCHORS, N_ANCHORS)
      USE QCI_CONSTRAINT_KEYS, ONLY: BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
      INTEGER, INTENT(IN)  :: ATOMS(:), NAT
      LOGICAL, INTENT(IN)  :: IN_SET(:)
      INTEGER, INTENT(OUT) :: ANCHORS(:)
      INTEGER, INTENT(OUT) :: N_ANCHORS
      INTEGER :: J1, J2, ATOM, NB

      N_ANCHORS = 0
      DO J1 = 1, NAT
         ATOM = ATOMS(J1)
         DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
            NB = BONDS_PER_ATOM_LIST(ATOM, J2)
            IF (IN_SET(NB)) CYCLE
            IF (.NOT. IS_IN_LIST(NB, ANCHORS, N_ANCHORS)) THEN
               N_ANCHORS = N_ANCHORS + 1
               ANCHORS(N_ANCHORS) = NB
               END IF
            END DO
         END DO
   END SUBROUTINE RMH_COMPUTE_ANCHORS

   !> Grows the maximal connected subset of POOL_ATOMS reachable from
   !! atoms bonded to TARGET_ANCHOR, refusing to cross any pool atom
   !! that itself has an external bond to a DIFFERENT anchor (such an
   !! atom would give the resulting group a second anchor, so it and
   !! everything only reachable through it are excluded).
   SUBROUTINE RMH_GROW_FROM_ANCHOR(POOL_ATOMS, POOL_SIZE, IN_SET, TARGET_ANCHOR, OUT_ATOMS, OUT_SIZE)
      USE QCIKEYS, ONLY: NATOMS
      USE QCI_CONSTRAINT_KEYS, ONLY: BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
      INTEGER, INTENT(IN)  :: POOL_ATOMS(:), POOL_SIZE
      LOGICAL, INTENT(IN)  :: IN_SET(:)
      INTEGER, INTENT(IN)  :: TARGET_ANCHOR
      INTEGER, INTENT(OUT) :: OUT_ATOMS(:)
      INTEGER, INTENT(OUT) :: OUT_SIZE

      LOGICAL :: BLOCKED(POOL_SIZE), VISITED(POOL_SIZE)
      INTEGER :: ATOM_TO_LOCAL(NATOMS)
      INTEGER :: QUEUE(POOL_SIZE), QHEAD, QTAIL, CUR, CUR_ATOM
      INTEGER :: J1, J2, NB, NB_LOCAL, ATOM

      ATOM_TO_LOCAL(:) = 0
      DO J1 = 1, POOL_SIZE
         ATOM_TO_LOCAL(POOL_ATOMS(J1)) = J1
      END DO

      ! an atom is blocked if it has an external bond to any anchor
      ! OTHER than the one we're growing towards
      BLOCKED(:) = .FALSE.
      DO J1 = 1, POOL_SIZE
         ATOM = POOL_ATOMS(J1)
         DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
            NB = BONDS_PER_ATOM_LIST(ATOM, J2)
            IF (IN_SET(NB)) CYCLE
            IF (NB /= TARGET_ANCHOR) THEN
               BLOCKED(J1) = .TRUE.
               EXIT
            END IF
         END DO
      END DO

      VISITED(:) = .FALSE.
      QHEAD = 1; QTAIL = 0
      OUT_SIZE = 0

      ! seed: unblocked pool atoms directly bonded to TARGET_ANCHOR
      DO J1 = 1, POOL_SIZE
         IF (BLOCKED(J1)) CYCLE
         ATOM = POOL_ATOMS(J1)
         DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
            IF (BONDS_PER_ATOM_LIST(ATOM,J2) == TARGET_ANCHOR) THEN
               VISITED(J1) = .TRUE.
               QTAIL = QTAIL + 1
               QUEUE(QTAIL) = J1
               EXIT
               END IF
            END DO
      END DO

      DO WHILE (QHEAD <= QTAIL)
         CUR = QUEUE(QHEAD); QHEAD = QHEAD + 1
         CUR_ATOM = POOL_ATOMS(CUR)
         OUT_SIZE = OUT_SIZE + 1
         OUT_ATOMS(OUT_SIZE) = CUR_ATOM

         DO J2 = 1, N_BONDS_PER_ATOM(CUR_ATOM)
            NB = BONDS_PER_ATOM_LIST(CUR_ATOM, J2)
            NB_LOCAL = ATOM_TO_LOCAL(NB)
            IF (NB_LOCAL == 0) CYCLE       ! not in pool
            IF (BLOCKED(NB_LOCAL)) CYCLE
            IF (VISITED(NB_LOCAL)) CYCLE
            VISITED(NB_LOCAL) = .TRUE.
            QTAIL = QTAIL + 1
            QUEUE(QTAIL) = NB_LOCAL
         END DO
         END DO

   END SUBROUTINE RMH_GROW_FROM_ANCHOR

   !> OUT_ATOMS = ATOMS(1:NAT) with every atom in TO_REMOVE(1:N_REMOVE) removed.
   SUBROUTINE RMH_REMOVE_ATOMS(ATOMS, NAT, TO_REMOVE, N_REMOVE, OUT_ATOMS, N_OUT)
      USE QCIKEYS, ONLY: NATOMS
      INTEGER, INTENT(IN)  :: ATOMS(:), NAT
      INTEGER, INTENT(IN)  :: TO_REMOVE(:), N_REMOVE
      INTEGER, INTENT(OUT) :: OUT_ATOMS(:)
      INTEGER, INTENT(OUT) :: N_OUT
      LOGICAL :: REMOVE_MASK(NATOMS)
      INTEGER :: J1

      CALL RMH_MEMBERSHIP(TO_REMOVE, N_REMOVE, REMOVE_MASK)
      N_OUT = 0
      DO J1 = 1, NAT
         IF (.NOT. REMOVE_MASK(ATOMS(J1))) THEN
            N_OUT = N_OUT + 1
            OUT_ATOMS(N_OUT) = ATOMS(J1)
         END IF
      END DO
   END SUBROUTINE RMH_REMOVE_ATOMS

   !> Splits ATOMS(1:NAT) into its connected components using only
   !! internal bonds (both endpoints within ATOMS). COMP_ATOMS(k,:) /
   !! COMP_SIZES(k) hold component k.
   SUBROUTINE RMH_CONNECTED_COMPONENTS(ATOMS, NAT, COMP_ATOMS, COMP_SIZES, NCOMP)
      USE QCIKEYS, ONLY: NATOMS
      USE QCI_CONSTRAINT_KEYS, ONLY: BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
      INTEGER, INTENT(IN)  :: ATOMS(:), NAT
      INTEGER, INTENT(OUT) :: COMP_ATOMS(:,:)
      INTEGER, INTENT(OUT) :: COMP_SIZES(:)
      INTEGER, INTENT(OUT) :: NCOMP

      LOGICAL :: IN_SET(NATOMS), VISITED(NAT)
      INTEGER :: ATOM_TO_LOCAL(NATOMS)
      INTEGER :: QUEUE(NAT), QHEAD, QTAIL, CUR, NB, NB_LOCAL
      INTEGER :: J1, J2

      CALL RMH_MEMBERSHIP(ATOMS, NAT, IN_SET)
      ATOM_TO_LOCAL(:) = 0
      DO J1 = 1, NAT
         ATOM_TO_LOCAL(ATOMS(J1)) = J1
         END DO

      VISITED(:) = .FALSE.
      NCOMP = 0

      DO J1 = 1, NAT
         IF (VISITED(J1)) CYCLE
         NCOMP = NCOMP + 1
         COMP_SIZES(NCOMP) = 0

         VISITED(J1) = .TRUE.
         QHEAD = 1; QTAIL = 1
         QUEUE(1) = J1

         DO WHILE (QHEAD <= QTAIL)
            CUR = QUEUE(QHEAD); QHEAD = QHEAD + 1
            COMP_SIZES(NCOMP) = COMP_SIZES(NCOMP) + 1
            COMP_ATOMS(NCOMP, COMP_SIZES(NCOMP)) = ATOMS(CUR)

            DO J2 = 1, N_BONDS_PER_ATOM(ATOMS(CUR))
               NB = BONDS_PER_ATOM_LIST(ATOMS(CUR), J2)
               IF (.NOT. IN_SET(NB)) CYCLE
               NB_LOCAL = ATOM_TO_LOCAL(NB)
               IF (NB_LOCAL == 0) CYCLE
               IF (VISITED(NB_LOCAL)) CYCLE
               VISITED(NB_LOCAL) = .TRUE.
               QTAIL = QTAIL + 1
               QUEUE(QTAIL) = NB_LOCAL
            END DO
         END DO
      END DO
   END SUBROUTINE RMH_CONNECTED_COMPONENTS

   !> Rare-case safety net for pools where no single-anchor group of
   !! sufficient size could be grown directly (e.g. every contact point
   !! with every anchor is itself a multi-anchor junction atom): trims
   !! one atom at a time -- preferring the lowest-internal-degree atom
   !! that touches an anchor -- re-splitting to the largest remaining
   !! component after each removal, until a single anchor remains or
   !! the pool drops below MIN_SIZE.
   SUBROUTINE RMH_GREEDY_TRIM(ATOMS, NAT, MIN_SIZE, SUCCESS, DEBUG_MODE)
      USE QCIKEYS, ONLY: NATOMS
      USE QCI_CONSTRAINT_KEYS, ONLY: BONDS_PER_ATOM_LIST, N_BONDS_PER_ATOM
      INTEGER, INTENT(INOUT) :: ATOMS(:)
      INTEGER, INTENT(INOUT) :: NAT
      INTEGER, INTENT(IN)    :: MIN_SIZE
      LOGICAL, INTENT(OUT)   :: SUCCESS
      LOGICAL, INTENT(IN)    :: DEBUG_MODE

      LOGICAL :: IN_SET(NATOMS), TOUCHES_ANCHOR
      INTEGER :: ANCHORS(NAT), N_ANCHORS
      INTEGER :: INTERNAL_DEGREE(NAT)
      INTEGER :: COMP_ATOMS(NAT, NAT), COMP_SIZES(NAT), NCOMP, BEST_COMP
      INTEGER :: ATOM, NEIGHBOR, J1, J2, J3
      INTEGER :: REMOVE_IDX, BEST_DEGREE

      SUCCESS = .FALSE.

      DO WHILE (NAT >= MIN_SIZE)

         CALL RMH_MEMBERSHIP(ATOMS, NAT, IN_SET)
         CALL RMH_COMPUTE_ANCHORS(ATOMS, NAT, IN_SET, ANCHORS, N_ANCHORS)

         IF (N_ANCHORS == 1) THEN
            SUCCESS = .TRUE.
            RETURN
         END IF
         IF (N_ANCHORS == 0 .OR. NAT == MIN_SIZE) THEN
            SUCCESS = .FALSE.
            RETURN
      END IF

         DO J1 = 1, NAT
            ATOM = ATOMS(J1)
            INTERNAL_DEGREE(J1) = 0
            DO J2 = 1, N_BONDS_PER_ATOM(ATOM)
               IF (IN_SET(BONDS_PER_ATOM_LIST(ATOM,J2))) INTERNAL_DEGREE(J1) = INTERNAL_DEGREE(J1) + 1
            END DO
         END DO

         REMOVE_IDX = -1
         BEST_DEGREE = HUGE(1)
         DO J1 = 1, NAT
            ATOM = ATOMS(J1)
            TOUCHES_ANCHOR = .FALSE.
            DO J3 = 1, N_BONDS_PER_ATOM(ATOM)
               NEIGHBOR = BONDS_PER_ATOM_LIST(ATOM, J3)
               IF (.NOT. IN_SET(NEIGHBOR)) THEN
                  TOUCHES_ANCHOR = .TRUE.
                  EXIT
               END IF
            END DO
            IF (TOUCHES_ANCHOR .AND. INTERNAL_DEGREE(J1) < BEST_DEGREE) THEN
               REMOVE_IDX = J1
               BEST_DEGREE = INTERNAL_DEGREE(J1)
            END IF
         END DO

         IF (REMOVE_IDX == -1) THEN
            SUCCESS = .FALSE.
            RETURN
         END IF

         IF (DEBUG_MODE) PRINT *, "rmh_greedy_trim> removing atom ", ATOMS(REMOVE_IDX)

         DO J1 = REMOVE_IDX, NAT - 1
            ATOMS(J1) = ATOMS(J1 + 1)
         END DO
         NAT = NAT - 1

         CALL RMH_CONNECTED_COMPONENTS(ATOMS, NAT, COMP_ATOMS, COMP_SIZES, NCOMP)
         BEST_COMP = 1
         DO J1 = 2, NCOMP
            IF (COMP_SIZES(J1) > COMP_SIZES(BEST_COMP)) BEST_COMP = J1
         END DO
         NAT = COMP_SIZES(BEST_COMP)
         ATOMS(1:NAT) = COMP_ATOMS(BEST_COMP, 1:NAT)

      END DO

      SUCCESS = .FALSE.

   END SUBROUTINE RMH_GREEDY_TRIM

   !> Commits ATOMS(1:NAT) as a finished group, growing LINEAR_GROUPS /
   !! NINGROUP on demand so a single input group splitting into several
   !! output groups is never limited by the caller's original sizing.
   SUBROUTINE RMH_COMMIT_GROUP(ATOMS, NAT, LINATOM_MASK, GROUPID)
      INTEGER, INTENT(IN)    :: ATOMS(:), NAT
      LOGICAL, INTENT(INOUT) :: LINATOM_MASK(:)
      INTEGER, INTENT(INOUT) :: GROUPID
      INTEGER :: J1

      GROUPID = GROUPID + 1
      CALL RMH_ENSURE_GROUP_CAPACITY(GROUPID, NAT)

      NINGROUP(GROUPID) = NAT
      DO J1 = 1, NAT
         LINEAR_GROUPS(GROUPID, J1) = ATOMS(J1)
         ATOM2LINGROUP(ATOMS(J1)) = GROUPID
         LINATOM_MASK(ATOMS(J1)) = .TRUE.
      END DO
   END SUBROUTINE RMH_COMMIT_GROUP

   SUBROUTINE RMH_DROP_ATOMS(LINATOM_MASK, ATOMS, NAT)
      LOGICAL, INTENT(INOUT) :: LINATOM_MASK(:)
      INTEGER, INTENT(IN)    :: ATOMS(:), NAT
      INTEGER :: J1
      DO J1 = 1, NAT
         LINATOM_MASK(ATOMS(J1)) = .FALSE.
      END DO
   END SUBROUTINE RMH_DROP_ATOMS

   !> Grows module-level LINEAR_GROUPS / NINGROUP (doubling) if GROUPID
   !! or GROUPSIZE_NEEDED exceeds the current allocation.
   SUBROUTINE RMH_ENSURE_GROUP_CAPACITY(GROUPID, GROUPSIZE_NEEDED)
      USE QCIKEYS, ONLY:NATOMS
      INTEGER, INTENT(IN) :: GROUPID, GROUPSIZE_NEEDED
      INTEGER, ALLOCATABLE :: TMP_GROUPS(:,:)
      INTEGER, ALLOCATABLE :: TMP_NINGROUP(:)
      INTEGER :: OLD_NGROUPS, OLD_MAXSIZE, NEW_NGROUPS, NEW_MAXSIZE


      IF (.NOT. ALLOCATED(ATOM2LINGROUP)) THEN
         ALLOCATE(ATOM2LINGROUP(NATOMS))
         ATOM2LINGROUP(:) = -1
               END IF

      OLD_NGROUPS = 0; OLD_MAXSIZE = 0
      IF (ALLOCATED(LINEAR_GROUPS)) THEN
         OLD_NGROUPS = SIZE(LINEAR_GROUPS, 1)
         OLD_MAXSIZE = SIZE(LINEAR_GROUPS, 2)
            END IF

      IF (GROUPID <= OLD_NGROUPS .AND. GROUPSIZE_NEEDED <= OLD_MAXSIZE) RETURN

      NEW_NGROUPS = MAX(GROUPID, OLD_NGROUPS * 2, 1)
      NEW_MAXSIZE = MAX(GROUPSIZE_NEEDED, OLD_MAXSIZE)

      ALLOCATE(TMP_GROUPS(NEW_NGROUPS, NEW_MAXSIZE)); TMP_GROUPS(:,:) = -1
      ALLOCATE(TMP_NINGROUP(NEW_NGROUPS));            TMP_NINGROUP(:) = 0

      IF (OLD_NGROUPS > 0) THEN
         TMP_GROUPS(1:OLD_NGROUPS, 1:OLD_MAXSIZE) = LINEAR_GROUPS(1:OLD_NGROUPS, 1:OLD_MAXSIZE)
         TMP_NINGROUP(1:OLD_NGROUPS) = NINGROUP(1:OLD_NGROUPS)
         DEALLOCATE(LINEAR_GROUPS)
         DEALLOCATE(NINGROUP)
      END IF

      CALL MOVE_ALLOC(TMP_GROUPS, LINEAR_GROUPS)
      CALL MOVE_ALLOC(TMP_NINGROUP, NINGROUP)
   END SUBROUTINE RMH_ENSURE_GROUP_CAPACITY

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
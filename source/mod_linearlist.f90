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
         IF (ALLOCATED(LINEAR_GROUPS)) DEALLOCATE(LINEAR_GROUPS)
         IF (ALLOCATED(ATOM2LINGROUP)) DEALLOCATE(ATOM2LINGROUP)
         IF (ALLOCATED(NINGROUP)) DEALLOCATE(NINGROUP)
      END SUBROUTINE DEALLOC_LINEAR_GROUP
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !IDEA FOR QUASI-RIGID BODY LINEAR LIST
      !DO NOT USE 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE DETECT_LINEAR()   
         USE QCIKEYS, ONLY: NATOMS, INLINLIST, LINEARBBT, ISBBATOM, QCIAMBERT, QCIHIRET       
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL, NCONPERATOM, CONLIST

         IMPLICIT NONE

         REAL(KIND=REAL64), PARAMETER :: TOLERANCE = 0.5D-1
         REAL(KIND=REAL64) :: DS, DF
         LOGICAL :: LINATOM(NATOMS) 
         INTEGER :: NLINHERE
         INTEGER :: A, B, J1, J2
    
          ! Grouping variables
         INTEGER :: CURRENTGROUP(NATOMS)      ! Group ID for each linear atom
         INTEGER :: GROUPS(NATOMS, NATOMS)    ! List of atoms in each group
                            ! Total number of groups
         INTEGER :: GROUPID
         INTEGER :: VISITED(NATOMS)           ! For DFS traversal
         INTEGER :: STACK(NATOMS)             ! Stack for DFS
         INTEGER :: STACKPTR
         INTEGER :: ATOM, NEIGHBOR
         INTEGER :: NINGROUP_TEMP(3*NATOMS)
         INTEGER :: MAXGROUPSIZE
                    !number of linear groups
         INTEGER :: NGROUPS
         LINATOM(:) = .FALSE.
    
         ! ========== PART 1: Detect linear atoms ==========
         DO J1 = 1, NATOMS
            NLINHERE = 0
        
            DO J2 = 1, NCONPERATOM(J1)
               A = J1
               B = CONLIST(A, J2)

               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, A, B, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, A, B, DF)

               IF(DABS(DS - DF) < TOLERANCE) NLINHERE = NLINHERE + 1
            END DO

            ! Mark atom as linear if all constraints have unchanged distances
            IF(NLINHERE == NCONPERATOM(J1)) LINATOM(J1) = .TRUE.
         END DO
    
         ! ========== PART 2: Group LINEAR atoms by constraint connectivity ==========
         CURRENTGROUP(:) = -1
         NINGROUP_TEMP(:) = 0
         GROUPS(:, :) = -1
         NGROUPS = 0
    
         ! Iterate through all atoms
         DO J1 = 1, NATOMS
            ! Skip if atom is not linear or already assigned to a group
            IF(.NOT. LINATOM(J1) .OR. CURRENTGROUP(J1) /= -1) CYCLE
                  
            ! Start a new group
            NGROUPS = NGROUPS + 1
            GROUPID = NGROUPS
        
            ! Use depth-first search (DFS) to find all connected LINEAR atoms
            VISITED(:) = 0
            STACKPTR = 1
            STACK(1) = J1
            VISITED(J1) = 1
        
            ! Process stack until empty
            DO WHILE(STACKPTR > 0)
                ATOM = STACK(STACKPTR)
               STACKPTR = STACKPTR - 1
               
               ! Assign atom to current group
               CURRENTGROUP(ATOM) = GROUPID
               NINGROUP_TEMP(GROUPID) = NINGROUP_TEMP(GROUPID) + 1
               
               ! Check bounds to avoid array overflow
               IF(NINGROUP_TEMP(GROUPID) > NATOMS) THEN
                  PRINT *, "ERROR: Group size exceeds NATOMS"
                  STOP
               END IF
               
               GROUPS(GROUPID, NINGROUP_TEMP(GROUPID)) = ATOM
               
               ! Find all neighbors (atoms connected via constraints)
               DO J2 = 1, NCONPERATOM(ATOM)
                  NEIGHBOR = CONLIST(ATOM, J2)
                  
                  ! Validate neighbor index
                  IF(NEIGHBOR < 1 .OR. NEIGHBOR > NATOMS) THEN
                     PRINT *, "ERROR: Invalid neighbor index ", NEIGHBOR
                     CYCLE
                  END IF
                  
                  ! Only process if neighbor is LINEAR
                  IF(.NOT. LINATOM(NEIGHBOR)) CYCLE
                  
                  ! If neighbor already in a group, merge groups
                  IF(CURRENTGROUP(NEIGHBOR) /= -1) THEN
                     IF(CURRENTGROUP(NEIGHBOR) /= GROUPID) THEN
                           CALL MERGE_GROUPS(CURRENTGROUP, GROUPS, NINGROUP, &
                                            GROUPID, CURRENTGROUP(NEIGHBOR))
                           !CONTINUE
                     END IF
                  ELSE IF(VISITED(NEIGHBOR) == 0) THEN
                     ! Neighbor not visited, add to stack
                     VISITED(NEIGHBOR) = 1
                     STACKPTR = STACKPTR + 1
                     
                     ! Check stack bounds
                     IF(STACKPTR > NATOMS) THEN
                           PRINT *, "ERROR: Stack overflow"
                           STOP
                     END IF
                     
                     STACK(STACKPTR) = NEIGHBOR
                  END IF
               END DO
            END DO
         END DO
    
         ! ========== Optional: Print results ==========
          DO GROUPID = 1, NGROUPS
              IF(NINGROUP_TEMP(GROUPID) > 0) THEN
                  PRINT *, "Linear Group ", GROUPID, ": ", NINGROUP_TEMP(GROUPID), " atoms"
                  DO J1 = 1, NINGROUP_TEMP(GROUPID)
                      ATOM = GROUPS(GROUPID, J1)
                      PRINT *, "  Atom ", ATOM
                  END DO
              END IF
          END DO

          MAXGROUPSIZE = MAXVAL(NINGROUP_TEMP)
          CALL ALLOC_LINEAR_GROUP(NGROUPS,MAXGROUPSIZE)
          ATOM2LINGROUP(:) = -1
          LINEAR_GROUPS(:,:) = -1
          NINGROUP(:) = 0
          GROUPID = 0
          DO J1 = 1, NGROUPS
            NINGROUP(J1) = NINGROUP_TEMP(J1)
            DO J2=1, NINGROUP_TEMP(J1)
                  LINEAR_GROUPS(J1,J2) = GROUPS(J1, J2)
                  ATOM2LINGROUP(GROUPS(J1,J2)) = J1 
            END DO
         END DO

   END SUBROUTINE DETECT_LINEAR

   SUBROUTINE MERGE_GROUPS(CURRENTGROUP, GROUPS, NINGROUP, GROUP1, GROUP2)
           
      USE QCIKEYS, ONLY: NATOMS
     
      INTEGER, INTENT(INOUT) :: CURRENTGROUP(:), GROUPS(:,:), NINGROUP(:)
      INTEGER, INTENT(IN) :: GROUP1, GROUP2
      INTEGER :: K, ATOM
        
      ! Move all atoms from GROUP2 to GROUP1
      DO K = 1, NINGROUP(GROUP2)
         ATOM = GROUPS(GROUP2, K)
         CURRENTGROUP(ATOM) = GROUP1
         NINGROUP(GROUP1) = NINGROUP(GROUP1) + 1
         GROUPS(GROUP1, NINGROUP(GROUP1)) = ATOM
      END DO
        
      ! Clear GROUP2
      NINGROUP(GROUP2) = 0
      GROUPS(GROUP2, :) = -1
   END SUBROUTINE MERGE_GROUPS

   SUBROUTINE GET_LIN_ROT_TRANSLATION(GROUP, GROUPSIZE, CXS, CXF, QSTART, QFINAL)
      USE QCIKEYS, ONLY: NATOMS
      USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
      USE QCIMINDIST, ONLY: FIND_ORIGIN, MOVE_COORDS, FIND_ALIGNMENT
      USE QUATERNIONS

      INTEGER, INTENT(IN) :: GROUP, GROUPSIZE !< group id
      REAL(KIND=REAL64), INTENT(OUT) :: CXS(3), CXF(3) !< centre of coordinates 
      REAL(KIND=REAL64), INTENT(OUT) :: QSTART(4), QFINAL(4) !< unit quaternions for SLERP
      REAL(KIND=REAL64) :: RS(3*GROUPSIZE), RF(3*GROUPSIZE) !group coordinates
      REAL(KIND=REAL64), DIMENSION(3,3) :: RMAT !rotational matrix 
      REAL(KIND=REAL64) :: DIST
      INTEGER :: J1, ATOMID

      DO J1=1, GROUPSIZE
         ATOMID = LINEAR_GROUPS(GROUP,J1)
         RS(3*(J1-1)+1:3*(J1-1)+2) = XSTART(3*(ATOMID-1)+1:3*(ATOMID-1)+2)
         RF(3*(J1-1)+1:3*(J1-1)+2) = XFINAL(3*(ATOMID-1)+1:3*(ATOMID-1)+2)
      END DO
      
      !centre RA
      CALL FIND_ORIGIN(GROUPSIZE,RS,CXS)
      CALL MOVE_COORDS(GROUPSIZE,RS,CXS)
      !centre RB around origin
      CALL FIND_ORIGIN(GROUPSIZE,RF,CXF)
      CALL MOVE_COORDS(GROUPSIZE,RF,CXF)

      !align coordinates
      CALL FIND_ALIGNMENT(GROUPSIZE, RS, RF, DIST, RMAT)

      QSTART = [1.0D0, 0.0D0, 0.0D0, 0.0D0]

      CALL MATRIX_TO_QUATERNION(RMAT,QFINAL)
      
   END SUBROUTINE GET_LIN_ROT_TRANSLATION

END MODULE QCI_LINEAR
MODULE QCICONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: NCONSTRAINT = 0
   INTEGER, ALLOCATABLE :: CONI(:), CONJ(:), CONION(:), CONJON(:)
   REAL(KIND = REAL64), ALLOCATABLE :: CONDISTREF(:)
   REAL(KIND = REAL64), ALLOCATABLE :: CONCUT(:)
   INTEGER, ALLOCATABLE :: NCONPERATOM(:), CONLIST(:,:)
   INTEGER :: MAXCONSTRAINTS

   CONTAINS

      SUBROUTINE SET_CONCUTABS()
         USE QCIKEYS, ONLY: CONCUTABS, CONCUTABS2, CONCUTABSSAVE, CONCUTABS2SAVE
         IMPLICIT NONE
         CONCUTABSSAVE = CONCUTABS
         CONCUTABS2SAVE = CONCUTABS2
      END SUBROUTINE SET_CONCUTABS

      ! create constraints and check them
      SUBROUTINE CREATE_CONSTRAINTS()
         USE QCIKEYS, ONLY: NATOMS, QCIAMBERT, QCIHIRET, QCISBMT
         USE SBM_CONSTRAINTS, ONLY: SBMMODEL_QCI_CONSTRAINTS, DEALLOC_SBM_CONST, SBM_NCONST, &
                                    SBM_CONI, SBM_CONJ, SBM_CONDISTREF, SBM_CONCUT
         USE AMBER_CONSTRAINTS, ONLY: AMBER_QCI_CONSTRAINTS, DEALLOC_AMBER_CONSTR, GET_BACKBONE, AMBER_NCONST, &
                                      AMBER_CONI, AMBER_CONJ, AMBER_CONDISTREF, AMBER_CONCUT
         USE HIRE_CONSTRAINTS, ONLY: HIRE_QCI_CONSTRAINTS, DEALLOC_HIRE_CONSTR, GET_BACKBONE_HIRE, HIRE_NCONST, &
                                     HIRE_CONI, HIRE_CONJ, HIRE_CONDISTREF, HIRE_CONCUT
         IMPLICIT NONE
         LOGICAL :: PERCT

         !call the routines that get us the constraints
         IF (QCIAMBERT) THEN
            ! get constraints
            CALL AMBER_QCI_CONSTRAINTS(NATOMS)
            ! get backbone information
            IF (QCIBBT) THEN
               CALL GET_BACKBONE(NATOMS)
            END IF
            ! set the global number of constraints
            NCONSTRAINT = AMBER_NCONST
            ! allocate the arrays and copy the data from AMBER module
            CALL ALLOC_CONSTR()
            CONI(1:NCONSTRAINT) = AMBER_CONI(1:AMBER_NCONST)
            CONJ(1:NCONSTRAINT) = AMBER_CONJ(1:AMBER_NCONST)
            CONDISTREF(1:NCONSTRAINT) = AMBER_CONDISTREF(1:AMBER_NCONST)
            CONCUT(1:NCONSTRAINT) = AMBER_CONCUT(1:AMBER_NCONST)
            ! deallocate information in AMBER module about constraints
            CALL DEALLOC_AMBER_CONSTR()
         ELSE IF (QCIHIRET) THEN
            CALL HIRE_QCI_CONSTRAINTS(NATOMS)
            IF (QCIBBT) THEN
               CALL GET_BACKBONE_HIRE(NATOMS)
            END IF
            NCONSTRAINT = HIRE_NCONST
            CALL ALLOC_CONSTR()
            CONI(1:NCONSTRAINT) = HIRE_CONI(1:HIRE_NCONST)
            CONJ(1:NCONSTRAINT) = HIRE_CONJ(1:HIRE_NCONST)
            CONDISTREF(1:NCONSTRAINT) = HIRE_CONDISTREF(1:HIRE_NCONST)
            CONCUT(1:NCONSTRAINT) = HIRE_CONCUT(1:HIRE_NCONST)
            CALL DEALLOC_HIRE_CONSTR()
         ELSE IF (QCISBMT) THEN
            CALL SBMMODEL_QCI_CONSTRAINTS(NATOMS)
            NCONSTRAINT = SBM_NCONST
            CALL ALLOC_CONSTR()
            CONI(1:NCONSTRAINT) = SBM_CONI(1:SBM_NCONST)
            CONJ(1:NCONSTRAINT) = SBM_CONJ(1:SBM_NCONST)
            CONDISTREF(1:NCONSTRAINT) = SBM_CONDISTREF(1:SBM_NCONST)
            CONCUT(1:NCONSTRAINT) = SBM_CONCUT(1:SBM_NCONST)
            CALL DEALLOC_HIRE_CONSTR()
         ELSE
            ! geometries and endpoints and potential additional constraints
            CALL GET_GEOMCONSTRAINTS()
         END IF
        
         ! check percolation 
         CALL CHECK_DUPLICATES()
         CALL CHECK_PERCOLATION(NATOMS, PERCT)
         IF (.NOT.PERCT) THEN
            WRITE(*,*) " create_con> Constraints are not including all atoms - STOP"
            STOP
         END IF
         !finally, get list of constraints and maximum number of constraints per atom
         CALL GET_NCON_PERATOM()
      END SUBROUTINE CREATE_CONSTRAINTS

      SUBROUTINE GET_GEOMCONSTRAINTS()
         USE CONGEOM
         IMPLICIT NONE
         LOGICAL :: PERCT
         INTEGER :: COUNTER
         REAL(KIND = REAL64), PARAMETER :: TOLMOD = 1.1
         INTEGER, PARAMETER :: MAXCYCLES = 5

         CALL READ_GEOMS()
         PERCT = .FALSE.
         COUNTER = 0
         READADDLIST = .FALSE.
         DO WHILE (.NOT.PERCT)
            COUNTER = COUNTER + 1
            IF (COUNTER.GT.MAXCYCLES) THEN
               WRITE(*,*) " get_geomconstraints> Canot not find a percolating network - STOP"
               STOP
            END IF
            IF (USEENDPOINTS) THEN
               CALL CREATE_FROM_ENDPOINTS(TOLMOD**(COUNTER-1))
            ELSE
               CALL CREATE_FROM_GEOMETRIES(TOLMOD**(COUNTER-1))
            END IF
            IF (.NOT.READADDLIST) THEN
               CALL ADD_CONSTRAINT_LIST()
               READADDLIST = .TRUE.
            END IF
            NCONSTRAINT = NGEOMCONST + NADDCONSTR
            CALL ALLOC_CONSTR()
            CONI(1:NGEOMCONST) = GEOMCONI(1:NGEOMCONST)
            CONJ(1:NGEOMCONST) = GEOMCONJ(1:NGEOMCONST)
            CONCUT(1:NGEOMCONST) = GEOMCONCUT(1:NGEOMCONST)
            CONDISTREF(1:NGEOMCONST) = GEOMCONDISTREF(1:NGEOMCONST)
            IF (NADDCONSTR.GT.0) THEN
               CONI((NGEOMCONST+1):(NGEOMCONST+NADDCONSTR)) = FILE_CONI(1:NADDCONSTR)
               CONJ((NGEOMCONST+1):(NGEOMCONST+NADDCONSTR)) = FILE_CONJ(1:NADDCONSTR)
               CONDISTREF((NGEOMCONST+1):(NGEOMCONST+NADDCONSTR)) = FILE_CONDISTREF(1:NADDCONSTR)
               CONCUT((NGEOMCONST+1):(NGEOMCONST+NADDCONSTR)) = FILE_CONCUT(1:NADDCONSTR)
            END IF
            CALL CHECK_PERCOLATION(NATOMS, PERCT)
         END DO
         CALL DEALLOC_CONGEOM()
         CALL DEALLOC_ADDCONST()
      END SUBROUTINE GET_GEOMCONSTRAINTS


      SUBROUTINE ALLOC_CONSTR()
         USE QCI_KEYS, ONLY:NATOMS
         CALL DEALLOC_CONSTR()
         ALLOCATE(CONI(NCONSTRAINT))
         ALLOCATE(CONJ(NCONSTRAINT))
         ALLOCATE(CONDISTREF(NCONSTRAINT))
         ALLOCATE(CONCUT(NCONSTRAINT))
         ALLOCATE(NCONPERATOM(NATOMS))        
      END SUBROUTINE ALLOC_CONSTR

      SUBROUTINE DEALLOC_CONSTR()
         IF (ALLOCATED(CONI)) DEALLOCATE(CONI)
         IF (ALLOCATED(CONJ)) DEALLOCATE(CONJ)
         IF (ALLOCATED(CONDISTREF)) DEALLOCATE(CONDISTREF)
         IF (ALLOCATED(CONCUT)) DEALLOCATE(CONCUT)
         IF (ALLOCATED(NCONPERATOM)) DEALLOCATE(NCONPERATOM)
         IF (ALLOCATED(CONLIST)) DEALLOCATE(CONLIST)
      END SUBROUTINE DEALLOC_CONSTR    

      SUBROUTINE CHECK_DUPLICATES()
         IMPLICIT NONE
         INTEGER :: NEWNCONST
         INTEGER :: NDUPL
         INTEGER :: NEWCONI(NCONSTRAINT)
         INTEGER :: NEWCONJ(NCONSTRAINT)
         REAL(KIND = REAL64) :: NEWCONDISTREF(NCONSTRAINT)
         REAL(KIND = REAL64) :: NEWCONCUT(NCONSTRAINT)                
         REAL(KIND = REAL64), PARAMETER :: TOL=1.0D-5

         NDUPL = 0
         NEWNCONST = 0

         DO J1=1,NCONSTRAINT
            I = CONI(J1)
            J = CONJ(J1)
            DUPLICATE = .FALSE.
            DO J2=1,NEWNCONST
               IF (((NEWCONI(J2).EQ.I).AND.(NEWCONJ(J2).EQ.J)).OR. &
                   ((NEWCONI(J2).EQ.J).AND.(NEWCONJ(J2).EQ.I))) THEN
                  DUPLICATE = .TRUE.
                  NDUPL = NDUPL + 1
                  IF (((NEWCONDISTREF(J2)-CONDISTREF(J1)).GT.TOL).OR.
                      ((NEWCONUT(J2)-CONCUT(J1)).GT.TOL)) THEN
                     WRITE(*,*) " WARNING> Duplicate has different constraints:"
                     WRITE(*,*) " CONDISTREF: ", NEWCONDISTREF(J2), CONDISTREF(J1)
                     WRITE(*,*) " CONCUT: ", NEWCONCUT(J2), CONCUT(J1)
                  END IF
                  EXIT
               END IF
            END DO
            IF (.NOT.DUPLICATE) THEN
               ! sanity check so we have valid constraints
               IF (I.NE.J) THEN
                  NEWNCONST = NEWNCONST + 1
                  ! sort i and j
                  IF (I.LT.J) THEN
                     NEWCONI(NEWNCONST) = I
                     NEWCONJ(NEWNCONST) = J
                  ELSE IF (J.LT.I) THEN
                     NEWCONI(NEWNCONST) = J
                     NEWCONJ(NEWNCONST) = I
                  END IF                  
                  NEWCONDISTREF(NEWNCONST) = CONDISTREF(J1)
                  NEWCONCUT(NEWNCONST) = CONCUT(J1)
               END IF
            END IF
         END DO

         IF (NDUPL.GT.0) THEN
            DEALLOCATE(CONI,CONJ,CONDISTREF,CONCUT)
            NCONSTRAINT = NEWNCONST
            ALLOCATE(CONI(NCONSTRAINT))
            ALLOCATE(CONJ(NCONSTRAINT))            
            ALLOCATE(CONDISTREF(NCONSTRAINT))
            ALLOCATE(CONCUT(NCONSTRAINT))
            CONI(1:NCONSTRAINT) = NEWCONI(1:NEWNCONST)
            CONJ(1:NCONSTRAINT) = NEWCONJ(1:NEWNCONST)
            CONDISTREF(1:NCONSTRAINT) = NEWCONDISTREF(1:NEWNCONST)
            CONCUT(1:NCONSTRAINT) = NEWCONCUT(1:NEWNCONST)
         END IF

      END SUBROUTINE CHECK_DUPLICATES

      SUBROUTINE CHECK_PERCOLATION(NATOMS, PERCT)
         USE QCIKEYS, ONLY: NQCIFROZEN, QCIFROZEN
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         LOGICAL, INTENT(OUT) :: PERCT
         INTEGER :: CURRENTGROUP(NATOMS) !current group for each atom
         INTEGER :: NINGROUP(NATOMS) !number of atoms in group
         INTEGER :: GROUPS(NATOMS,NATOMS) !list of atoms currently in each group
         INTEGER :: J1, IDAT, ID1, ID2

         CURRENTGROUP(1:NATOMS) = -1
         NINGROUP(1:NATOMS) = 0
         GROUPS(1:NATOMS,1:NATOMS) = -1

         IF (NQCIFROZEN.GT.0) THEN
            DO J1=1,NATOMS
               IF (QCIFROZEN(J1)) THEN
                  CURRENTGROUP(J1) = 1
                  NINGROUP(1) = NINGROUP(1) + 1
                  GROUPS(1,NINGROUP(1)) = J1
               END IF
            END DO
         END IF

         DO J1=1,NCONSTRAINT
            ID1 = CONI(J1)
            ID2 = CONJ(J1)
            GROUP1 = CURRENTGROUP(ID1)
            GROUP2 = CURRENTGROUP(ID2)
            ! if neither of them is assigned yet, make a new group
            IF ((GROUP1.EQ.-1).AND.(GROUP2.EQ.-1)) THEN
               DO J2=1,NATOMS
                  IF (NINGROUP(J2).EQ.0) THEN
                     GROUPS(J2,1) = ID1
                     GROUPS(J2,2) = ID2  
                     CURRENTGROUP(ID1) = J2
                     CURRENTGROUP(ID2) = J2
                     NINGROUP(J2) = 2
                     EXIT
                  ELSE
                     CYCLE
                  END IF                
               END DO
            ! if they are in the same group
            ELSE IF (GROUP1.EQ.GROUP2) THEN
               CYCLE
            ! if one of them is in a grop, simply add the other one to that group
            ELSE IF ((GROUP1.NE.-1).AND.(GROUP2.EQ.-1)) THEN
               NINGROUP(GROUP1) = NINGROUP(GROUP1) + 1
               CURRENTGROUP(ID2) = GROUP1
               GROUPS(GROUP1,NINGROUP(GROUP1)) = ID2
            ELSE IF ((GROUP1.EQ.-1).AND.(GROUP2.NE.-1)) THEN
               NINGROUP(GROUP2) = NINGROUP(GROUP2) + 1
               CURRENTGROUP(ID1) = GROUP2
               GROUPS(GROUP2,NINGROUP(GROUP2)) = ID1
            ! the atoms are in different groups - merge the groups
            ELSE            
               ! we merge from higher id to lower id
               IF (GROUP1.LT.GROUP2) THEN
                  DO J2=1,NINGROUP(GROUP2)
                     NINGROUP(GROUP1) = NINGROUP(GROUP1) + 1
                     CURRENTGROUP(GROUPS(GROUP2,J2)) = GROUP1
                     GROUPS(GROUP1,NINGROUP(GROUP1)) = GROUPS(GROUP2,J2)
                     GROUPS(GROUP2,J2) = -1
                  END DO
                  NINGROUP(GROUP2) = 0
               ELSE IF (GROUP1.GT.GROUP2) THEN
                  DO J2=1,NINGROUP(GROUP1)
                     NINGROUP(GROUP2) = NINGROUP(GROUP2) + 1
                     CURRENTGROUP(GROUPS(GROUP1,J2)) = GROUP2
                     GROUPS(GROUP2,NINGROUP(GROUP2)) = GROUPS(GROUP1,J2)
                     GROUPS(GROUP1,J2) = -1
                  END DO
                  NINGROUP(GROUP1) = 0
               END IF
            END IF
         END DO

         IF (NINGROUP(1).NE.NATOMS) THEN
            WRITE(*,*) " check_perc> constraints do not form a fully percolating network"
            PERCT = .FALSE.
         ELSE
            PERCT = .TRUE.
         END IF
      END SUBROUTINE CHECK_PERCOLATION

      SUBROUTINE GET_NCON_PERATOM()
         USE QCI_KLEYS, ONLY: NATOMS
         IMPLICIT NONE
         INTEGER :: J1, JMAX
         NCONPERATOM(1:NATOMS) = 0
         JMAX = 0
         DO J1=1,NCONSTRAINT
            NCONPERATOM(CONI(J1))=NCONPERATOM(CONI(J1))+1
            NCONPERATOM(CONJ(J1))=NCONPERATOM(CONJ(J1))+1
         ENDDO
         MAXCONSTRAINTS=-1
         DO J1=1,NATOMS
            IF (NCONATOM(J1).GT.MAXCONSTRAINTS) THEN
               MAXCONSTRAINTS=NCONPERATOM(J1)
               JMAX=J1
            ENDIF
         ENDDO
         WRITE(*,'(A,I6,A,I6)') ' getnconperatom> maximum constraints ',MAXCONSTRAINTS,' for atom ',JMAX

         ALLOCATE(CONLIST(NATOMS,MAXCONSTRAINTS))
         CONLIST(1:NATOMS,1:MAXCONSTRAINTS)=0
         NCONPERATOM(1:NATOMS)=0
         DO J1=1,NCONSTRAINT
            NCONPERATOM(CONI(J1))=NCONPERATOM(CONI(J1))+1
            NCONPERATOM(CONJ(J1))=NCONPERATOM(CONJ(J1))+1
            CONLIST(CONI(J1),NCONPERATOM(CONI(J1)))=CONJ(J1)
            CONLIST(CONJ(J1),NCONPERATOM(CONJ(J1)))=CONI(J1)
         ENDDO

      END SUBROUTINE GET_NCON_PERATOM

END MODULE QCICONSTRAINTS
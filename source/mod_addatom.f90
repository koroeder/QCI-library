MODULE ADDATOM
   USE QCIPREC
   IMPLICIT NONE
   LOGICAL :: BBDONE = .FALSE.
   INTEGER, ALLOCATABLE :: NCONTOACTIVE(:)
   !limits on finding atoms for local axis set
   REAL(KIND=REAL64) :: QCIDISTCUT = 10.0D0
   INTEGER :: QCIATOMSEP = 15

   CONTAINS

      SUBROUTINE ALLOC_ADDATOM()
         USE QCIKEYS, ONLY: NATOMS
         CALL DEALLOC_ADDATOM()
         ALLOCATE(NCONTOACTIVE(NATOMS))
      END SUBROUTINE ALLOC_ADDATOM

      SUBROUTINE DEALLOC_ADDATOM()
         IF (ALLOCATED(NCONTOACTIVE)) DEALLOCATE(NCONTOACTIVE)
      END SUBROUTINE DEALLOC_ADDATOM

      SUBROUTINE ADDATOM()
         USE QCIKEYS, ONLY: QCIDOBACK, NATOMS, DEBUG, QCILINEART
         USE REPULSION, ONLY: NNREPULSIVE, NREPULSIVE
         USE QCI_LINEAR, ONLY: NQCILINEAR, INLINLIST
         USE CONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE
         IMPLICIT NONE
         INTEGER :: NTOADD, NADDED                          !number to be added and number already added
         INTEGER :: NNREPSAVE, NREPSAVE                     !variables for saving repulsion list
         LOGICAL :: MORETOADD                               !logical switch to stay in main loop
         REAL(KIND = REAL64) :: INVDTOACTIVE(1:NACTIVE)     !inverse distances of constraints to active atoms
         INTEGER :: NMAXCON                                 !current largest number of constraints to active set from any inactive atom
         LOGICAL :: CHOSENACID                              !switch whether we currently adding a residue in its entirety
         INTEGER :: ACID                                    !id of residue currently added in full
         INTEGER :: NEXTATOM, NCONTOACT                     !atom to be added and the number of constraints it has to the active set
         REAL(KIND=REAL64) :: SHORTESTCON                   !shortest distance of constraint to active set for NEXTATOM
         REAL(KIND=REAL64) :: BESTDIST(NATOMS)              !list of sorted average distance to newatom
         INTEGER :: BESTIDX, NDISTNEWATOM                   !associated ids and total number found
         REAL(KIND=REAL64) :: BESTCONDIST(NATOMS)           !list of sorted average constraint distances to newatom
         INTEGER :: BESTCONIDX, NCONNEWATOM                 !associated ids and total number found
      
         ! setup book keeping
         NTOADD = 1
         IF (QCILINEART) NTOADD = NQCILINEAR - 2
         NADDED = 0
         CHOSENACID = .FALSE.
         ! check whether we have more backbone atoms to add
         IF (QCIDOBACK.AND.(.NOT.BBDONE)) CALL CHECK_BBLIST()

         !Save current repulsion to speed up checks later
         NNREPSAVE=NNREPULSIVE
         NREPSAVE=NREPULSIVE

         !Set variable for tracking whether we completed adding atoms
         MORETOADD = .TRUE.
         DO WHILE (MORETOADD)
            !get list of atoms by number of constraints to current active atoms
            CALL CREATE_NCONTOACTIVE_LIST(INVDTOACTIVE,NMAXCON)
            !find the next atom to be added
            !the priorities are: 1. QCILINEAR, 2. CHOOSEACID, 3. QCIDOBACK, and then the highest number of constraints to the active set
            CALL FIND_NEXT_ATOM(CHOSENACID,ACID,NEXTATOM,NCONTOACT,SHORTESTCON)
            !set chosenacid if needed
            IF ((QCIADDACIDT.AND.(.NOT.QCIDOBACK)).OR.QCIDOBACKALL) THEN
               IF (.NOT.CHOSENACID) THEN
                  ACID = ATOMS2RES(NEWATOM)
                  CHOOSEACID = .TRUE.
               END IF
            END IF
            
            WRITE(*,*) " addatom> Adding atom ", NEWATOM, ", which has ", NCONTOACT, " constraints to active set out of maximum ", NMAXCON
            WRITE(*,*) "          Shortest distance constraint to active set is: ", SHORTESTCON

            ! The interpolation for the new atom relies on a local axis system formed by three atoms.
            ! We look for a sorted list, according to how well the end point distace is preserved.
            ! We sort by shortest average distance to avoid distant atoms having accidentally well preserved distance.
            CALL GET_ATOMS_BY_DISTANCE(NEWATOM,NDISTNEWATOM,BESTDIST,BESTIDX)
            ! Now update the constraints including the new atom 
            CALL UPDATE_CONSTRAINTS(NEWATOM,NCONNEWATOM,BESTCONDIST,BESTCONIDX)
            !update the repulsions - TOO continue here (line 64 in sub_doatom and line 2845)
            CALL UPDATE_REPULSIONS(NEWATOM)

            !activate new atom
            ATOMACTIVE(NEWATOM)=.TRUE.
            NACTIVE=NACTIVE+1 
            TURNONORDER(NACTIVE)=NEWATOM
            !check consistency
            CALL CHECK_NACTIVE()


         END DO


      END SUBROUTINE ADDATOM

      SUBROUTINE CHECK_NACTIVE()
         USE INTERPOLATION_KEYS, ONLY: NACTIVE, ATOMACTIVE
         USE QCIKEYS, ONLY: NATOMS
         USE MOD_TERMINATE, ONLY: INT_ERR_TERMINATE
         IMPLICIT NONE
         INTEGER :: NDUMMY, J1

         NDUMMY = 0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         END DO
         IF (NDUMMY.NE.NACTIVE) THEN
            WRITE(*,*) " check_active> Inconsistent number of active atoms"
            WRITE(*,*) "               ", NDUMMY, " atoms active, should be: ", NACTIVE
            CALL INT_ERR_TERMINATE()
         END IF
      END SUBROUTINE CHECK_NACTIVE

      SUBROUTINE UPDATE_REPULSIONS(NEWATOM)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, QCIREPCUT
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE  
         USE REPULSION
         USE CONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ

         IMPLICIT NONE  
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER :: J1
         LOGICAL :: ISCONSTRAINED  
         REAL(KIND = REAL64) :: DF, DS , DMIN      
         
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE !ignore atoms that are not active
            IF (ABS(J1-NEWATOM).LE.SEPREPULSION) CYCLE !no repulsions if atoms are close in sequence
            ISCONSTRAINED = .FALSE.
            DO J2=1,NCONSTRAINT
               IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.NEWATOM)).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.NEWATOM))) THEN
                  ISCONSTRAINED = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT.ISCONSTRAINED) THEN
               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, NEWATOM, J1, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, NEWATOM, J1, DF)
               DMIN=MIN(DS,DF)
               DMIN=MIN(DMIN,QCIREPCUT)
               NREPULSIVE=NREPULSIVE+1
               IF (NREPULSIVE.GT.NREPMAX) CALL DOUBLE_ALLOC_REP()
               REPI(NREPULSIVE)=J1
               REPJ(NREPULSIVE)=NEWATOM
               REPCUT(NREPULSIVE)=DMIN
            END IF
         END DO

      END SUBROUTINE UPDATE_REPULSIONS


      SUBROUTINE UPDATE_CONSTRAINTS(NEWATOM,NCONNEWATOM,BESTCONDIST,BESTCONIDX)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, MAXCONUSE
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, ATOMACTIVE
         USE CONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREF
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         REAL(KIND=REAL64), INTENT(OUT) :: BESTCONDIST(NATOMS)  !sorted list of average distances
         INTEGER, INTENT(OUT)  :: BESTCONIDX(NATOMS)            !associated list of atom ids
         INTEGER, INTENT(OUT) :: NCONNEWATOM                    !number of atoms found that are within the cut offs and active     
         INTEGER :: J1, ATOM1, ATOM2

         BESTCONDIST(1:NATOMS) = 1.0D100
         BESTCONIDX(1:NATOMS) = -1
         NCONNEWATOM = 0

         DO J1=1,NCONSTRAINT
            IF (CONACTIVE(J1)) CYCLE
            ATOM1 = CONI(J1)
            ATOM2 = CONJ(J1)
            !if either atom1 is the new atom and atom2 is active or vice versa, this is a new possible constraint
            IF (((ATOM1.EQ.NEWATOM).AND.(ATOMACTIVE(ATOM2))).OR.((ATOM2.EQ.NEWATOM).AND.(ATOMACTIVE(ATOM1)))) THEN
               NCONNEWATOM = NCONNEWATOM + 1
               DO J2=1,NCONNEWATOM
                  IF (CONDISTREF(J1).LT.BESTCONDIST(J2)) THEN
                     DO J3=NCONNEWATOM,J2+1,-1
                        BESTCONDIST(J3)=BESTCONDIST(J3-1)
                        BESTCONIDX(J3)=BESTCONIDX(J3-1)
                     END DO
                     BESTCONDIST(J2) = CONDISTREF(J1)
                     IF (ATOM1.EQ.NEWATOM) BESTCONIDX(J2) = ATOM2
                     IF (ATOM2.EQ.NEWATOM) BESTCONIDX(J2) = ATOM1
                     EXIT
                  END IF
               END DO
            END IF              
         END DO
         
         IF (DEBUG) THEN
            WRITE(*,*) " get_constraints_by_dist> Constraints including new atom by average reference distance:"
            WRITE(*,'(10I6)') BESTDIST(1:MIN(10,NCONNEWATOM))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10G12.4)') BESTIDX(1:MIN(10,NCONNEWATOM))
         END IF

         ! Turning on constraints identified
         DO J1=1,MIN(MAXCONUSE,NCONNEWATOM)
            DO J2=1,NCONSTRAINT
               ATOM1 = CONI(J2)
               ATOM2 = CONJ(J2)
               IF (((ATOM1.EQ.NEWATOM).AND.(ATOM2.EQ.BESTCONIDX(J1))).OR.((ATOM2.EQ.NEWATOM).AND.(ATOM1.EQ.BESTCONIDX(J1)))) THEN
                  CONACTIVE(J2) = .TRUE.
                  IF (DEBUG) THEN
                     WRITE(*,'(A,I6,A,2I6)') ' addatom> Turning on constraint ',J2,' for atoms ', ATOM1, ATOM2
                  END IF
               END IF
            END DO
         END DO



      END SUBROUTINE UPDATE_CONSTRAINTS

      SUBROUTINE GET_ATOMS_BY_DISTANCE(NEWATOM,NDISTNEWATOM,BESTDIST,BESTIDX)
         USE QCIKEYS, ONLY: NATOMS, DEBUG
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         REAL(KIND=REAL64), INTENT(OUT) :: BESTDIST(NATOMS)  !sorted list of average distances
         INTEGER, INTENT(OUT)  :: BESTIDX(NATOMS)            !associated list of atom ids
         INTEGER, INTENT(OUT) :: NDISTNEWATOM                !number of atoms found that are within the cut offs and active
         REAL(KIND = REAL64) : : DF, DS, AVD
         INTEGER :: J1, J2, STARTIDX, ENDIDX

         BESTDIST(1:NATOMS) = 1.0D100
         BESTIDX(1:NATOMS) = -1
         NDISTNEWATOM = 0
         
         !we limit the search in sequence of atoms, so we ignore distant atoms
         STARTIDX = MAX(1,NEWATOM-QCIATOMSEP)
         ENDIDX = MIN(NATOMS,NEWATOM+QCIATOMSEP)

         DO J1=STARTIDX,ENDIDX
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            CALL DISTANCE_TWOATOMS(NATOMS, XSTART, NEWATOM, J1, DS)
            CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, NEWATOM, J1, DF)
            !ignore atoms if they are separated by a large distance
            IF ((DS.GT.QCIDISTCUT).OR.(DF.GT.QCIDISTCUT)) CYCLE
            AVD=0.5*(DS+DF)
            NDISTNEWATOM = NDISTNEWATOM + 1
            ! go through sorting loop to mke sure we maintain an ordered list
            DO J2=1,NDISTNEWATOM
               IF (AVD.LT.BESTDIST(J2)) THEN
                  DO J3=NDISTNEWATOM,J2+1,-1
                     BESTDIST(J3)=BESTDIST(J3-1)
                     BESTIDX(J3)=BESTIDX(J3-1)
                  END DO
                  BESTDIST(J2) = AVD
                  BESTIDX(J2) = J1
                  EXIT
               END IF
            END DO
         END DO

         IF (DEBUG) THEN
            WRITE(*,*) " get_atoms_by_dist> Closest atoms to new atom by average endpoint distance:"
            WRITE(*,'(10I6)') BESTDIST(1:MIN(10,NDISTNEWATOM))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10G12.4)') BESTIDX(1:MIN(10,NDISTNEWATOM))
         END IF

      END SUBROUTINE GET_ATOMS_BY_DISTANCE

      SUBROUTINE FIND_NEXT_ATOM(CHOSENACID,ACID,NEXTATOM,NCONTOACT,SHORTESTCON)
         USE QCIKEYS, ONLY: QCILINEAR, INLINLIST, ATOMS2RES, QCIDOBACK, ISBBATOM
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE, CONACTIVE
         USE CONSTRAINTS, ONLY: NCONSTRAINT, CONDISTREF
         USE MOD_TERMINATE, ONLY: INT_ERR_TERMINATE
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ACID
         LOGICAL, INTENT(IN) :: CHOSENACID
         INTEGER, INTENT(OUT) :: NEXTATOM
         INTEGER, INTENT(OUT) :: NCONTOACT
         REAL(KIND=REAL64), INTENT(OUT)  :: SHORTESTCON
         INTEGER :: J1, IDXACTIVE, IDXINACTIVE

         !resetting all dummy variuables
         NEWATOM = 0
         NCONTOACT = 0
         SHORTESTCON = 1.0D100

         DO J1=1,NCONSTRAINT
            !ignore active constraints
            IF (CONACTIVE(J1)) CYCLE
            !set the active and inactive atoms in the constraint and cycle if neither or both (shouldn't happen) are active
            IF (ATOMACTIVE(CONI(J1)).AND.(.NOT.ATOMACTIVE(CONJ(J1)))) THEN
               IDXACTIVE = CONI(J1)
               IDXINACTIVE = CONJ(J1)
            ELSE IF (ATOMACTIVE(CONJ(J1)).AND.(.NOT.ATOMACTIVE(CONI(J1)))) THEN
               IDXACTIVE = CONJ(J1)
               IDXINACTIVE = CONI(J1) 
            ELSE IF (ATOMACTIVE(CONJ(J1)).AND.ATOMACTIVE(CONI(J1))) THEN
               WRITE(*,*) " find_next_atom> WARNING: atoms ", CONI(J1), " and ", CONJ(J1), " are active, but the constraint", J1," between them is not!"
               CYCLE
            ELSE 
               CYCLE
            END IF 
            ! Now test for various options to get the next atom
            ! 1. Are we using QCIlinear, and is the inactive atom in the list?
            IF (QCILINEAR.AND.(.NOT.INLINLIST(IDXINACTIVE))) CYCLE
            ! 2. Is CHOSENACID set, and is the inactive atom in the residue to be added?
            IF (CHOSENACID.AND.(.NOT.(ATOMS2RES(IDXINACTIVE).EQ.ACID))) CYCLE
            ! 3. Are we adding backbone atoms, and is the inactive atom a backbone atom?
            IF (QCIDOBACK.AND.(.NOT.BBDONE).AND.(.NOT.ISBBATOM(IDXINACTIVE))) CYCLE

            IF (NCONTOACTIVE(IDXINACTIVE).GE.NCONTOACT) THEN
               IF (CONDISTREF(J1).LT.SHORTESTCON) THEN
                  SHORTESTCON = CONDISTREF(J1)
                  NEWATOM = IDXINACTIVE
                  NCONTOACT = NCONTOACTIVE(IDXINACTIVE)
               END IF
            END IF
         END DO

         !sanity check that we have found an atom to be added - if not we terminate
         IF (NEWATOM*NCONTOACT.EQ.0) THEN
            WRITE(*,*) " find_next_atom> Error - new active atom not set, NEWATOM: ", NEWATOM, " NCONTOACT: ", NCONTOACT
            CALL INT_ERR_TEMRINATE()
         END IF


      END SUBROUTINE FIND_NEXT_ATOM


      SUBROUTINE CREATE_NCONTOACTIVE_LIST(INVDTOACTIVE,NBEST)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE, CONACTIVE
         USE QCICONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREF
         INTEGER, INTENT(OUT) :: NBEST ! maximum number of constraintson inactive atom
         REAL(KIND = REAL64), INTENT(OUT) :: INVDTOACTIVE(1:NACTIVE)  !inverse distance
         INTEGER :: I, ATOM1, ATOM2
         REAL(KIND = REAL64) :: INVDIST 

         NBEST = 0
         NCONTOACTIVE(1:NATOMS) = 0
         INVDTOACTIVE(1:NATOMS) = 0.0D0
         
         DO I=1,NCONSTRAINT
            IF (CONACTIVE(I)) CYCLE
            ATOM1 = CONI(I)
            ATOM2 = CONJ(I)
            INVDIST = 1.0D0/CONDISTREF(I)
            !check if the first atom is active and the second atom is inactive
            IF (ATOMACTIVE(ATOM1).AND.(.NOT.ATOMACTIVE(ATOM2))) THEN
               NCONTOACTIVE(ATOM2) = NCONTOACTIVE(ATOM2) + 1
               IF (INVDIST.GT.INVDTOACTIVE(ATOM2)) INVDTOACTIVE(ATOM2)=INVDIST
             !check if the second atom is active and the first atom is inactive  
            ELSE IF (ATOMACTIVE(ATOM2).AND.(.NOT.ATOMACTIVE(ATOM1))) THEN
               NCONTOACTIVE(ATOM1) = NCONTOACTIVE(ATOM1) + 1
               IF (INVDIST.GT.INVDTOACTIVE(ATOM1)) INVDTOACTIVE(ATOM1)=INVDIST               
            END IF
            !update current best
            IF (NCONTOACTIVE(ATOM1).GT.BEST) NBEST = NCONTOACTIVE(ATOM1)
            IF (NCONTOACTIVE(ATOM2).GT.BEST) NBEST = NCONTOACTIVE(ATOM2)
         END DO
      END SUBROUTINE CREATE_NCONTOACTIVE_LIST

      SUBROUTINE CHECK_BBLIST()
         USE QCIKEYS, ONLY: NATOMS
         DO J1=1,NATOMS
            IF (AABACK(J1)) THEN
               IF (.NOT.ATOMACTIVE(J1)) THEN
                  ! if there are inactive BB atoms, return
                  RETURN
               END IF
            ENDIF
         ENDDO
         BBDONE = .TRUE.
         IF (DEBUG) WRITE(*,*) " check_bblist> All backbone atoms are active"
      END SUBROUTINE CHECK_BBLIST


END MODULE ADDATOM
MODULE ADDINGATOM
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

      !debugging function
      SUBROUTINE CHECK_DIFF_FINAL()
         USE MOD_INTCOORDS, ONLY: XYZ, XFINAL
         USE QCIKEYS, ONLY: NATOMS, NIMAGES

         REAL(KIND = REAL64) :: DIFF(3*NATOMS)

         DIFF(1:3*NATOMS) = XYZ((3*NATOMS)*(NIMAGES+1)+1:(3*NATOMS)*(NIMAGES+2)) - XFINAL(1:3*NATOMS) 

         WRITE(*,*) " COMPARING FINISH COORDS: ", SUM(DIFF)
      END SUBROUTINE CHECK_DIFF_FINAL

      SUBROUTINE ADDATOM()
         USE MOD_INTCOORDS, ONLY: XYZ, EEE, GGG, RMS
         USE QCIKEYS, ONLY: QCIDOBACK, QCIADDACIDT, NATOMS, NIMAGES, DEBUG, QCILINEART, INLINLIST, &
                            QCITRILATERATION, QCIDOBACKALL, ATOMS2RES, ISBBATOM, CHECKCHIRAL, QCIUSEGROUPS
         USE REPULSION, ONLY: NNREPULSIVE, NREPULSIVE, CHECKREP
         USE CONSTR_E_GRAD, ONLY: CONGRAD
         USE QCI_LINEAR, ONLY: NQCILINEAR
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, NACTIVE, TURNONORDER, ATOMACTIVE
         USE MOD_INTCOORDS, ONLY: WRITE_ACTIVE_BAND
         USE CHIRALITY, ONLY: GET_ACTIVE_CHIRAL_CENTRES, NCHIRAL, CHIRALITY_CHECK, CHECK_SINGLE_CHIRAL_CENTRE
         USE AMBER_CONSTRAINTS, ONLY: GROUPLOOKUP, CURRENT_GROUP, CURRENTLY_ADDING_GROUP, INGROUP, PLACINGGROUPS, SIZEPLACINGGROUPS
         IMPLICIT NONE
         INTEGER :: NTOADD, NADDED                          !number to be added and number already added
         INTEGER :: NNREPSAVE, NREPSAVE                     !variables for saving repulsion list
         LOGICAL :: MORETOADD                               !logical switch to stay in main loop
         REAL(KIND = REAL64) :: INVDTOACTIVE(1:NATOMS)      !inverse distances of constraints to active atoms
         INTEGER :: NMAXCON                                 !current largest number of constraints to active set from any inactive atom
         LOGICAL :: CHOSENACID                              !switch whether we currently adding a residue in its entirety
         INTEGER :: ACID                                    !id of residue currently added in full
         INTEGER :: NEXTATOM, NCONTOACT                     !atom to be added and the number of constraints it has to the active set
         REAL(KIND=REAL64) :: SHORTESTCON                   !shortest distance of constraint to active set for NEXTATOM
         REAL(KIND=REAL64) :: LOCALDIST(4)
         INTEGER :: NLOCAL, LOCALIDX(4)
         LOGICAL :: ADDEDTHISCYCLE                          !have we added an atom this cycle?
         REAL(KIND = REAL64), PARAMETER :: FRAC = 1.0D0     !fraction for linear interpolation
         INTEGER :: THISIMAGE, NEWATOMOFFSET, CONONEOFFSET, ENDPOINT !offsets used to simplify linear interpolation
         REAL(KIND = REAL64) :: STARTWEIGHT, ENDWEIGHT      !weights of endpoints in linear interpolation
         REAL(KIND = REAL64) :: ELIN, ECON, EDIST, ETOTAL   !linear, constraint and distance measure-based energies
         REAL(KIND = REAL64) :: XLIN(NIMAGES,3), XCON(NIMAGES,3), XDIST(NIMAGES,3) !saved coordinates for interpolations
         INTEGER :: J1, IDX1, IDX2, IDX3, IDX4
         CHARACTER(LEN=6) :: ATOMSTRING
         INTEGER :: INIT_NCHIRACTIVE, FINAL_NCHIRACTIVE, CHIRALCENTRE
         LOGICAL :: INIT_ACTIVE_CHIR_CENTRES(NCHIRAL), FINAL_ACTIVE_CHIR_CENTRES(NCHIRAL)
         LOGICAL :: ALLADDED
      
         ! setup book keeping
         NTOADD = 1
         IF (QCILINEART) NTOADD = MIN(NQCILINEAR - 2,NQCILINEAR-NACTIVE)
         NADDED = 0
         CHOSENACID = .FALSE.
         ! check whether we have more backbone atoms to add
         IF (QCIDOBACK.AND.(.NOT.BBDONE)) CALL CHECK_BBLIST()

         !Save current repulsion to speed up checks later
         NNREPSAVE=NNREPULSIVE
         NREPSAVE=NREPULSIVE

         ! call check on number of chiral centres
         IF (CHECKCHIRAL) CALL GET_ACTIVE_CHIRAL_CENTRES(INIT_NCHIRACTIVE,INIT_ACTIVE_CHIR_CENTRES)

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
                  ACID = ATOMS2RES(NEXTATOM)
                  CHOSENACID = .TRUE.
               END IF
            END IF
            
            WRITE(*,'(A,I6,A,I4,A,I4)') "  addatom> Adding atom ", NEXTATOM, ", which has ", NCONTOACT, " constraints to active set out of maximum ", NMAXCON
            WRITE(*,'(A,F8.4)') "           Shortest distance constraint to active set is: ", SHORTESTCON

            ! We update the cosntraints and get a list of constraint and closest atoms to construct local 
            ! if we have four atoms, we use all four to build the local axis system, if we have three, we use three, otherwise we go linear
            CALL GET_ATOMS_FOR_LOCAL_AXIS(NEXTATOM,NLOCAL,LOCALIDX,LOCALDIST)
            !update the repulsions
            CALL UPDATE_REPULSIONS(NEXTATOM)

            !activate new atom
            ATOMACTIVE(NEXTATOM)=.TRUE.
            NACTIVE=NACTIVE+1 
            TURNONORDER(NACTIVE)=NEXTATOM
            !check consistency
            CALL CHECK_NACTIVE()

            !check whether we have a new group or completed the old one
            IF (QCIUSEGROUPS) THEN
               !if we have an active group, check whether we added all atoms
               IF (CURRENTLY_ADDING_GROUP) THEN
                  ALLADDED = .TRUE.
                  DO J1=1,SIZEPLACINGGROUPS(CURRENT_GROUP)
                     IF (.NOT.ATOMACTIVE(PLACINGGROUPS(CURRENT_GROUP,J1))) THEN
                        ALLADDED = .FALSE.
                        EXIT
                     END IF
                  END DO
                  !if so, we are done with the group and turn the group addition option opff for now
                  IF (ALLADDED) THEN
                     CURRENT_GROUP = 0
                     CURRENTLY_ADDING_GROUP = .FALSE.
                  END IF
               ! if we are not actively adding a group, we have to check whether we started a new group
               ELSE
                  IF (INGROUP(NEXTATOM)) THEN
                     CURRENTLY_ADDING_GROUP = .TRUE.
                     CURRENT_GROUP = GROUPLOOKUP(NEXTATOM)
                  END IF
               END IF
            END IF
            

            !now we need to actually add the atom
            ADDEDTHISCYCLE = .FALSE.
            ELIN = 1.0D100
            ECON = 1.0D100

            !if we have three or more constraints, we use them and construct a local axis system
            IF (NLOCAL.GE.3) THEN
               ADDEDTHISCYCLE = .TRUE.
               !use vectors within the local coordinate frame
               !CALL PLACE_ATOM(NEXTATOM,NLOCAL,LOCALIDX)
               !build internal coordinates for interpolation
               CALL PLACE_INTERNALS(NEXTATOM,NLOCAL,LOCALIDX)
               IF (QCITRILATERATION) THEN
                  CALL TRILATERATE_ATOMS(NEXTATOM,LOCALIDX,LOCALDIST)
               END IF
               ! before we continue check repulsion neighbour list
               CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
               ! call congrad routine
               CALL CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)
               ECON = ETOTAL
               DO J1=1,NIMAGES
                  NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                  XCON(J1,1:3) = XYZ((NEWATOMOFFSET+1):(NEWATOMOFFSET+3)) 
               END DO
      
            END IF
            !if we haven't suceeded in adding the atom or it is QCIlinear, go for a linear interpolation for new atom based on the tightest constraint
            IF ((.NOT.ADDEDTHISCYCLE).OR.QCILINEART) THEN
               DO J1=1,NIMAGES
                  ! the images are technically running from 2 to NIMAGES+1, but the offset is (J1-1)*3*NATOMs, 
                  ! so we can just use a shifted range from 1 to NIMAGES
                  THISIMAGE = J1*3*NATOMS
                  ENDPOINT = 3*NATOMS*(NIMAGES+1)
                  NEWATOMOFFSET = 3*(NEXTATOM-1)
                  !there is always at least one constraint!
                  CONONEOFFSET = 3*(LOCALIDX(1)-1)
                  STARTWEIGHT = FRAC*(NIMAGES+1-J1)/(NIMAGES+1) 
                  ENDWEIGHT = 1.0D0*J1/(NIMAGES+1) !need to keep this 1.0D0* in here to convert the type correctly
                  !X coordinate
                  XYZ(THISIMAGE+NEWATOMOFFSET+1) = XYZ(THISIMAGE+CONONEOFFSET+1) + STARTWEIGHT*(XYZ(NEWATOMOFFSET+1)-XYZ(CONONEOFFSET+1)) + &
                                                   ENDWEIGHT*(XYZ(ENDPOINT+NEWATOMOFFSET+1)-XYZ(ENDPOINT+CONONEOFFSET+1))
                  !Y coordinate
                  XYZ(THISIMAGE+NEWATOMOFFSET+2) = XYZ(THISIMAGE+CONONEOFFSET+2) + STARTWEIGHT*(XYZ(NEWATOMOFFSET+2)-XYZ(CONONEOFFSET+2)) + &
                                                   ENDWEIGHT*(XYZ(ENDPOINT+NEWATOMOFFSET+2)-XYZ(ENDPOINT+CONONEOFFSET+2))
                  !Z coordinate
                  XYZ(THISIMAGE+NEWATOMOFFSET+3) = XYZ(THISIMAGE+CONONEOFFSET+3) + STARTWEIGHT*(XYZ(NEWATOMOFFSET+3)-XYZ(CONONEOFFSET+3)) + &
                                                   ENDWEIGHT*(XYZ(ENDPOINT+NEWATOMOFFSET+3)-XYZ(ENDPOINT+CONONEOFFSET+3))    
               END DO                                                                                

               CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
               ! call congrad routine
               CALL CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)
               ELIN = ETOTAL
               DO J1=1,NIMAGES
                  NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                  XLIN(J1,1:3) = XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) 
               END DO
            END IF   

            !select which interpolation is used based on energy if multiple are used (we shouldn't have this case except for QCIlinear)
            WRITE(*,*) "ECON, ELIST: ", ECON, ELIN
            IF (QCILINEART) THEN
               WRITE(*,*) " addatom> Using linear interpolation for new atom ", NEXTATOM, " from linear list"
            ELSE 
               IF ((ELIN.LT.ECON)) THEN
                  WRITE(*,*) " addatom> Using linear interpolation for new atom ", NEXTATOM
               ELSE
                  WRITE(*,*) " addatom> Using interpolation from preserved distances for new atom ", NEXTATOM
                  ETOTAL = ECON
                  DO J1=1,NIMAGES
                     NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                     XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) = XCON(J1,1:3)
                  END DO
               END IF
            END IF

            NADDED = NADDED + 1
            IF (NADDED.LT.NTOADD) CYCLE !if we already know we need to add more, skip ahead and add the next atom

            !we should not have more atoms to add unless we hit certain criteria
            MORETOADD = .FALSE.
            ! if we add one residue at a time, check whether we added all of them
            IF (QCIADDACIDT.AND.(.NOT.QCIDOBACK)) THEN
               IF (NLOCAL.LT.3) THEN
                  WRITE(*,*) " addatom> Added atom constraints ", NLOCAL, " is less than 3 - not attempting to add another atom."
               ELSE
                  ! iterate over all atoms and check whether they are in the residue and not active
                  DO J1=1,NATOMS
                     IF ((ATOMS2RES(J1).EQ.ACID).AND.(.NOT.(ATOMACTIVE(J1)))) THEN
                        !we found another atom to add, so we set moretoadd to true and leave this loop
                        MORETOADD = .TRUE.
                        EXIT
                     END IF
                  END DO
               END IF
            END IF

            IF (QCIDOBACKALL.AND.ISBBATOM(NEXTATOM)) THEN
               DO J1=1,NATOMS
                  IF ((ATOMS2RES(J1).EQ.ACID).AND.(ISBBATOM(J1)).AND.(.NOT.(ATOMACTIVE(J1)))) THEN
                     MORETOADD = .TRUE.
                     EXIT
                  END IF
               END DO           
            END IF
            !this is the end of the add atom loop - the loop will continue if MORETOADD is set to TRUE, otherwise we leave the loop
            WRITE(ATOMSTRING,'(I6)') NEXTATOM
            CALL WRITE_ACTIVE_BAND("int.active.addedatom_"//TRIM(ADJUSTL(ATOMSTRING))//".xyz")
         END DO

         !check number of chiral centres after addition
         IF (CHECKCHIRAL.AND..NOT.CURRENTLY_ADDING_GROUP) THEN
            CALL GET_ACTIVE_CHIRAL_CENTRES(FINAL_NCHIRACTIVE,FINAL_ACTIVE_CHIR_CENTRES)
         
            IF (INIT_NCHIRACTIVE.NE.FINAL_NCHIRACTIVE) THEN
               IF (INIT_NCHIRACTIVE+1.EQ.FINAL_NCHIRACTIVE) THEN
                  DO J1=1,NCHIRAL
                     IF (.NOT.INIT_ACTIVE_CHIR_CENTRES(J1).AND.FINAL_ACTIVE_CHIR_CENTRES(J1)) THEN
                        CHIRALCENTRE=J1
                     END IF
                  END DO
                  WRITE(*,*) " addatom> New chiral centre has been activated - checking this centre"
                  CALL CHECK_SINGLE_CHIRAL_CENTRE(CHIRALCENTRE,XYZ)
               ELSE
                  WRITE(*,*) " addatom> Multiple chiral centres activated - checking entire band"
                  CALL CHIRALITY_CHECK(XYZ)
               END IF
            END IF
         END IF

         CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
         ! call congrad routine
         CALL CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)
         !we are done with QCIlinear, so set it to false
         QCILINEART = .FALSE.
         WRITE(*,*) "FINISHED LINEAR INTERPOLATION"
         CALL WRITE_ACTIVE_BAND("after_linear.xyz")


      END SUBROUTINE ADDATOM

      SUBROUTINE GET_ATOMS_FOR_LOCAL_AXIS(NEWATOM,NLOCAL,LOCALIDX,LOCALDIST)
         USE QCIKEYS, ONLY: NATOMS
         INTEGER , INTENT(IN) :: NEWATOM
         REAL(KIND=REAL64), INTENT(OUT) :: LOCALDIST(4)
         INTEGER , INTENT(OUT) :: NLOCAL, LOCALIDX(4)
         REAL(KIND=REAL64) :: BESTDIST(NATOMS)              !list of sorted average distance to newatom
         INTEGER :: BESTIDX(NATOMS), NDISTNEWATOM           !associated ids and total number found
         REAL(KIND=REAL64) :: BESTCONDIST(NATOMS)           !list of sorted average constraint distances to newatom
         INTEGER :: BESTCONIDX(NATOMS), NCONNEWATOM         !associated ids and total number found
         REAL(KIND=REAL64) :: SECCONDIST(NATOMS)            !list of sorted average constraint distances to newatom
         INTEGER :: SECCONIDX(NATOMS), NSECCONSTR           !associated ids and total number found

         INTEGER :: NCONST, NDISTS, NCOUNT, NSEC
         INTEGER :: IDX1, IDX2, IDX3, IDX4, J

         ! The interpolation for the new atom relies on a local axis system formed by three atoms.
         ! We look for a sorted list, according to how well the end point distace is preserved.
         ! We sort by shortest average distance to avoid distant atoms having accidentally well preserved distance.
         CALL GET_ATOMS_BY_DISTANCE(NEWATOM,NDISTNEWATOM,BESTDIST,BESTIDX)
         ! Now update the constraints including the new atom, and get list of constraints ordered by distance
         CALL UPDATE_CONSTRAINTS(NEWATOM,NCONNEWATOM,BESTCONDIST,BESTCONIDX)

         NCONST = MIN(4,NCONNEWATOM)
         NDISTS = MIN(4,NDISTNEWATOM)
         NCOUNT = 0
         LOCALDIST(1:4) = 0.0D0
         LOCALIDX(1:4) = -1
         IF (NCONST.NE.0) THEN
            DO J = 1,NCONST
               LOCALDIST(J) = BESTCONDIST(J)
               LOCALIDX(J) = BESTCONIDX(J)
               NCOUNT = NCOUNT + 1
            END DO
         END IF
         !if we don't have enough constraints - get atoms constraint to the constraints sorted by distance
         IF (NCOUNT.LT.4) THEN
            CALL GET_SECONDARY_CONSTRAINTS(NEWATOM,NCONNEWATOM,BESTCONIDX,NSECCONSTR,SECCONDIST,SECCONIDX)
            IF (NSECCONSTR.NE.0) THEN
               NSEC=MIN(4,NSECCONSTR)
               IF (NCOUNT.EQ.0) THEN
                  DO J = 1,NSEC
                     LOCALDIST(J) = SECCONDIST(J)
                     LOCALIDX(J) = SECCONIDX(J)
                     NCOUNT = NCOUNT + 1
                     IF (NCOUNT.EQ.4) EXIT    
                  END DO
               ELSE IF (NCOUNT.EQ.1) THEN
                  IDX1 = LOCALIDX(1)
                  DO J = 1,NSEC
                     IDX2 = SECCONIDX(J)
                     IF (IDX1.NE.IDX2) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = SECCONDIST(J)
                        LOCALIDX(NCOUNT) = SECCONIDX(J)  
                     END IF   
                     IF (NCOUNT.EQ.4) EXIT                  
                  END DO
               ELSE IF (NCOUNT.EQ.2) THEN
                  IDX1 = LOCALIDX(1)
                  IDX2 = LOCALIDX(2)
                  DO J = 1,NSEC
                     IDX3 = SECCONIDX(J)
                     IF ((IDX1.NE.IDX3).AND.(IDX2.NE.IDX3)) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = SECCONDIST(J)
                        LOCALIDX(NCOUNT) = SECCONIDX(J)  
                     END IF 
                     IF (NCOUNT.EQ.4) EXIT                
                  END DO
               ELSE IF (NCOUNT.EQ.3) THEN
                  IDX1 = LOCALIDX(1)
                  IDX2 = LOCALIDX(2)
                  IDX3 = LOCALIDX(3)
                  DO J = 1,NSEC
                     IDX4 = SECCONIDX(J)
                     IF ((IDX1.NE.IDX4).AND.(IDX2.NE.IDX4).AND.(IDX3.NE.IDX4)) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = SECCONDIST(J)
                        LOCALIDX(NCOUNT) = SECCONIDX(J)  
                     END IF    
                     IF (NCOUNT.EQ.4) EXIT                 
                  END DO
               END IF
            END IF
         END IF

         IF (NLOCAL.LT.4) THEN
            IF (NDISTS.NE.0) THEN
               IF (NCOUNT.EQ.0) THEN
                  DO J = 1,NDISTS
                     LOCALDIST(J) = BESTDIST(J)
                     LOCALIDX(J) = BESTIDX(J)
                     NCOUNT = NCOUNT + 1
                     IF (NCOUNT.EQ.4) EXIT    
                  END DO
               ELSE IF (NCOUNT.EQ.1) THEN
                  IDX1 = LOCALIDX(1)
                  DO J = 1,NDISTS
                     IDX2 = BESTIDX(J)
                     IF (IDX1.NE.IDX2) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = BESTDIST(J)
                        LOCALIDX(NCOUNT) = BESTIDX(J)  
                     END IF   
                     IF (NCOUNT.EQ.4) EXIT                  
                  END DO
               ELSE IF (NCOUNT.EQ.2) THEN
                  IDX1 = LOCALIDX(1)
                  IDX2 = LOCALIDX(2)
                  DO J = 1,NDISTS
                     IDX3 = BESTIDX(J)
                     IF ((IDX1.NE.IDX3).AND.(IDX2.NE.IDX3)) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = BESTDIST(J)
                        LOCALIDX(NCOUNT) = BESTIDX(J)  
                     END IF 
                     IF (NCOUNT.EQ.4) EXIT                
                  END DO
               ELSE IF (NCOUNT.EQ.3) THEN
                  IDX1 = LOCALIDX(1)
                  IDX2 = LOCALIDX(2)
                  IDX3 = LOCALIDX(3)
                  DO J = 1,NDISTS
                     IDX4 = BESTIDX(J)
                     IF ((IDX1.NE.IDX4).AND.(IDX2.NE.IDX4).AND.(IDX3.NE.IDX4)) THEN
                        NCOUNT = NCOUNT + 1
                        LOCALDIST(NCOUNT) = BESTDIST(J)
                        LOCALIDX(NCOUNT) = BESTIDX(J)  
                     END IF    
                     IF (NCOUNT.EQ.4) EXIT                 
                  END DO
               END IF
            END IF
         END IF
         NLOCAL = NCOUNT
      END SUBROUTINE GET_ATOMS_FOR_LOCAL_AXIS


      SUBROUTINE TRILATERATE_ATOMS(NEWATOM,CONIDXLIST,CONDISTLIST)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ         
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER, INTENT(IN) :: CONIDXLIST(4)        
         REAL(KIND=REAL64), INTENT(IN) :: CONDISTLIST(4)
         REAL(KIND=REAL64) :: P1(3), P2(3), P3(3), R1, R2, R3
         REAL(KIND=REAL64) :: SOL1(3), SOL2(3), PREV(3), D1SQ, D2SQ
         INTEGER :: IMAGEOFFSET
         INTEGER :: N1, N2, N3, IDX1, IDX2, IDX3, J1
         LOGICAL :: FTEST

         WRITE(*,*) " trilaterate_atom> Adding atom via trilateration"

         !set initial guess to the best three constrained atoms
         N1=1; N2=2; N3=3
         IDX1 = CONIDXLIST(N1); IDX2 = CONIDXLIST(N2); IDX3 = CONIDXLIST(N3)

         DO J1=2,NIMAGES+1
            IMAGEOFFSET = (J1-1)*3*NATOMS
            P1(1:3)=XYZ((IMAGEOFFSET+3*(IDX1-1)+1):(IMAGEOFFSET+3*(IDX1-1)+3))
            P2(1:3)=XYZ((IMAGEOFFSET+3*(IDX2-1)+1):(IMAGEOFFSET+3*(IDX2-1)+3))
            P3(1:3)=XYZ((IMAGEOFFSET+3*(IDX3-1)+1):(IMAGEOFFSET+3*(IDX3-1)+3))           
            R1=CONDISTLIST(N1)
            R2=CONDISTLIST(N2)
            R3=CONDISTLIST(N3) 
            CALL TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
            IF (.NOT.FTEST) THEN
               PREV(1:3) = XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3))
               D1SQ = (SOL1(1)-PREV(1))**2 + (SOL1(2)-PREV(2))**2 + (SOL1(3)-PREV(3))**2
               D2SQ = (SOL2(1)-PREV(1))**2 + (SOL2(2)-PREV(2))**2 + (SOL2(3)-PREV(3))**2
               IF (D1SQ.LT.D2SQ) THEN
                  XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3)) = SOL1(1:3)
               ELSE
                  XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3)) = SOL2(1:3)
               END IF
            END IF
         END DO
      END SUBROUTINE TRILATERATE_ATOMS

      !Intersection of three spheres with centres P1...P3 with radii R1...R3
      SUBROUTINE TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, NORM_VEC
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: P1(3), P2(3), P3(3), R1, R2, R3
         REAL(KIND=REAL64), INTENT(OUT) :: SOL1(3), SOL2(3)
         LOGICAL, INTENT(OUT) :: FTEST
         REAL(KIND=REAL64) :: V1(3), V2(3), V3(3), EX(3), EY(3), EZ(3)
         REAL(KIND=REAL64) :: NORM, NORM2, DOTX2, DOTY2, X, Y, Z, TEMP

         FTEST=.FALSE.
         !get vectors between centres
         V1(1:3) = P2(1:3) - P1(1:3)  
         V2(1:3) = P3(1:3) - P1(1:3)  
         !get normalised version of V1 
         CALL NORM_VEC(V1,EX,NORM)     

         DOTX2 = EX(1)*V2(1) + EX(2)*V2(2) + EX(3)*V2(3) 
         V3(1:3) = V2(1:3) - DOTX2*V1(1:3)
         CALL NORM_VEC(V3,EY,NORM2)
         EZ = CROSS_PROD(EX,EY)

         DOTY2 = EY(1)*V2(1) + EY(2)*V2(2) + EY(3)*V2(3) 

         X=(R1*R1 - R2*R2 + NORM*NORM) / (2.0D0*NORM)
         Y=(R1*R1 - R3*R3 -2.0D0*DOTX2*X + DOTX2*DOTX2 + DOTY2*DOTY2) / (2.0D0*DOTY2)
         TEMP = R1*R1 - X*X - Y*Y
   
         
         IF (TEMP.LT.0.0D0) THEN
            FTEST=.TRUE.
         ELSE
            FTEST=.FALSE.
            Z=SQRT(TEMP)
            SOL1(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) + Z*EZ(1:3)
            SOL2(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) - Z*EZ(1:3)
         END IF
      END SUBROUTINE TRILATERATION


      SUBROUTINE PLACE_INTERNALS(NEWATOM,NLOCAL,CONIDXLIST)
         USE QCIFILEHANDLER, ONLY: FILE_OPEN
         USE HELPER_FNCTS, ONLY: DOTP, EUC_NORM, NORM_VEC, DIHEDRAL, ANGLE, DISTANCE_SIMPLE, CROSS_PROD
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER, INTENT(IN) :: NLOCAL
         INTEGER, INTENT(IN) :: CONIDXLIST(4)
         
         REAL(KIND=REAL64) :: D1, D2, ANG1, ANG2, DIH1, DIH2
         REAL(KIND=REAL64) :: COORDS(12)
         REAL(KIND=REAL64) :: ALPHA
         REAL(KIND=REAL64) :: DNEW, ANGNEW, DIHNEW
         REAL(KIND=REAL64) :: C(3), V1(3), V2(3), V3(3), POS(3)
         REAL(KIND=REAL64) :: SINTH, COSTH, SINPHI, COSPHI

         REAL(KIND=REAL64) :: BC(3), BA(3), A(3), B(3), D(3), N(3), M(3), NORM

         INTEGER :: IDX1, IDX2, IDX3, IMAGEOFFSET
         INTEGER :: J1, I

         !!DEBUG - check the order of things for the initial and the final calculation - are we applying things to the correct atoms?
         IDX1 = CONIDXLIST(1); IDX2 = CONIDXLIST(2); IDX3 = CONIDXLIST(3)
         WRITE(*,'(A,I8)') " place_internal> New atom: ", NEWATOM

         ! start coordinates
         COORDS(1:3) = XYZ((3*(IDX3-1)+1):((3*(IDX3-1)+3)))
         COORDS(4:6) = XYZ((3*(IDX2-1)+1):((3*(IDX2-1)+3)))
         COORDS(7:9) = XYZ((3*(IDX1-1)+1):((3*(IDX1-1)+3)))
         COORDS(10:12) = XYZ((3*(NEWATOM-1)+1):((3*(NEWATOM-1)+3)))

         ! get distance, angle and dihedral
         CALL DISTANCE_SIMPLE(COORDS(7:9), COORDS(10:12), D1)
         ANG1 = ANGLE(COORDS(4:12))
         DIH1 = DIHEDRAL(COORDS)

         !finish coordinates
         IMAGEOFFSET = (3*NATOMS)*(NIMAGES+1)
         COORDS(1:3) = XYZ(IMAGEOFFSET+(3*(IDX3-1)+1):IMAGEOFFSET+((3*(IDX3-1)+3)))
         COORDS(4:6) = XYZ(IMAGEOFFSET+(3*(IDX2-1)+1):IMAGEOFFSET+((3*(IDX2-1)+3)))
         COORDS(7:9) = XYZ(IMAGEOFFSET+(3*(IDX1-1)+1):IMAGEOFFSET+((3*(IDX1-1)+3)))
         COORDS(10:12) = XYZ(IMAGEOFFSET+(3*(NEWATOM-1)+1):IMAGEOFFSET+((3*(NEWATOM-1)+3)))

         ! get distance, angle and dihedral
         CALL DISTANCE_SIMPLE(COORDS(7:9), COORDS(10:12), D2)
         ANG2 = ANGLE(COORDS(4:12))
         DIH2 = DIHEDRAL(COORDS)

         WRITE(*,*) "<><>Start: ", D1, ANG1, DIH1
         WRITE(*,*) "<><>Final: ", D2, ANG2, DIH2

         ! iterate over all images and place the new atom
         DO J1=2,NIMAGES+1
            IMAGEOFFSET = (J1-1)*3*NATOMS
            ! get linear interpolation in internal coordinates
            ALPHA = 1.0D0/J1
            DNEW = (1-ALPHA)*D1 + ALPHA*D2
            ANGNEW = (1-ALPHA)*ANG1 + ALPHA*ANG2
            DIHNEW = (1-ALPHA)*DIH1 + ALPHA*DIH2
            
            !!translate back into Cartesians
            ! orthogonal basis 
            A(1:3) = XYZ((IMAGEOFFSET+3*(IDX3-1)+1):(IMAGEOFFSET+3*(IDX3-1)+3))
            B(1:3) = XYZ((IMAGEOFFSET+3*(IDX2-1)+1):(IMAGEOFFSET+3*(IDX2-1)+3))
            C(1:3) = XYZ((IMAGEOFFSET+3*(IDX1-1)+1):(IMAGEOFFSET+3*(IDX1-1)+3))
            BC(1:3) = C(1:3) - B(1:3)
            BA(1:3) = A(1:3) - B(1:3)

            CALL NORM_VEC(BC,V1,NORM)
            CALL NORM_VEC(CROSS_PROD(BA,BC),V2,NORM)
            CALL NORM_VEC(CROSS_PROD(BC,V2),V3,NORM)

            COSTH = COS(ANGNEW)
            SINTH = SIN(ANGNEW)

            COSPHI = COS(DIHNEW)
            SINPHI = SIN(DIHNEW)

            DO I=1,3
               POS(I) = C(I) + DNEW*(-V1(I)*COSTH + V2(I)*SINTH*COSPHI + V3(I)*SINTH*SINPHI)
            END DO
            XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3)) =  POS(1:3)
         END DO
      END SUBROUTINE PLACE_INTERNALS


      SUBROUTINE PLACE_ATOM(NEWATOM,NLOCAL,CONIDXLIST)
         USE QCIFILEHANDLER, ONLY: FILE_OPEN
         USE HELPER_FNCTS, ONLY: DOTP, EUC_NORM, NORM_VEC
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER, INTENT(IN) :: NLOCAL
         INTEGER, INTENT(IN) :: CONIDXLIST(4)

         REAL(KIND=REAL64) :: VEC1(3), VEC2(3), B1(3), B2(3), B3(3), POS1(3), NEW(3), NEW2(3)
         REAL(KIND=REAL64) :: C(3), D(3), ROT(3), ANGLE, SCALED(3)
         REAL(KIND=REAL64) :: C1, C2, C3, D1, D2, D3, E1, E2, E3, SCALE, NORM
         INTEGER :: IMAGEOFFSET
         INTEGER :: IDX1, IDX2, IDX3, IDX4
         INTEGER :: J1

         IDX1 = CONIDXLIST(1); IDX2 = CONIDXLIST(2); IDX3 = CONIDXLIST(3); IDX4 = CONIDXLIST(4)
         
         WRITE(*,'(A,I8)') " place_atom> New atom: ", NEWATOM

         !Setting up local axis system for the start image
         IF (NLOCAL.EQ.3) THEN
            CALL GET_LOCAL_AXIS(IDX1,IDX2,IDX3,1,B1,B2,B3)
            WRITE(*,'(A)') " place_atom> Using three atoms for placement"
            IF (DEBUG) THEN
               WRITE(*,*) " place_atom> Use the following three active atoms to place new atom"
               WRITE(*,*) "             New atom: ", NEWATOM, "Closest active atoms: ", CONIDXLIST(1:3)
            END IF
         ELSE
            CALL GET_LOCAL_AXIS2(IDX1,IDX2,IDX3,IDX4,1,B1,B2,B3)
            WRITE(*,'(A)') " place_atom> Using four atoms for placement"
            IF (DEBUG) THEN
               WRITE(*,*) " place_atom> Use the following four active atoms to place new atom"
               WRITE(*,*) "             New atom: ", NEWATOM, "Closest active atoms: ", CONIDXLIST(1:4)
            END IF
         END IF
         !VEC1 is pointing from IDX1 to NEWATOM
         VEC1(1:3) = XYZ((3*(NEWATOM-1)+1):((3*(NEWATOM-1)+3))) - XYZ((3*(IDX1-1)+1):((3*(IDX1-1)+3)))
         !get relative coordinates of NEWATOM in B1,B2,B3 axis system
         C1 = DOTP(3,VEC1,B1)
         C2 = DOTP(3,VEC1,B2)
         C3 = DOTP(3,VEC1,B3)

         !alternative scheme using rotations
         !C = C1*B1(1:3) + C2*B2(1:3) + C3*B3(1:3)
         !CALL GET_ROT_AXIS(B1,B2,B3,ROT)

         !Setting up local axis system for the final image
         IF (NLOCAL.EQ.3) THEN
            CALL GET_LOCAL_AXIS(IDX1,IDX2,IDX3,NIMAGES+2,B1,B2,B3)
         ELSE
            CALL GET_LOCAL_AXIS2(IDX1,IDX2,IDX3,IDX4,NIMAGES+2,B1,B2,B3)
         END IF
         !VEC2 is pointing from IDX1 to NEWATOM in the last image
         IMAGEOFFSET = (3*NATOMS)*(NIMAGES+1)
         VEC2(1:3) = XYZ(IMAGEOFFSET+(3*(NEWATOM-1)+1):IMAGEOFFSET+((3*(NEWATOM-1)+3))) - &
                     XYZ(IMAGEOFFSET+(3*(IDX1-1)+1):(IMAGEOFFSET+(3*(IDX1-1)+3)))  
         !get relative coordinates of NEWATOM in B1,B2,B3 axis system
         D1 = DOTP(3,VEC2,B1)
         D2 = DOTP(3,VEC2,B2)
         D3 = DOTP(3,VEC2,B3)

         !D = D1*B1(1:3) + D2*B2(1:3) + D3*B3(1:3)
         !CALL GET_ROT_ANGLE(C,D,ROT,ANGLE)
         !WRITE(*,'(A,3F12.7,A,F12.7)') " Rotational axis: ", ROT, " angle from C to D: ", ANGLE

         SCALE = 0.5*(EUC_NORM(VEC1) + EUC_NORM(VEC2))
         !iterate over images and place NEWATOM
         DO J1=2,NIMAGES+1
            !get B1,B2,B3 for current image 
            IF (NLOCAL.EQ.3) THEN
               CALL GET_LOCAL_AXIS(IDX1,IDX2,IDX3,J1,B1,B2,B3)
            ELSE
               CALL GET_LOCAL_AXIS2(IDX1,IDX2,IDX3,IDX4,J1,B1,B2,B3)
            END IF         
            IMAGEOFFSET = (J1-1)*3*NATOMS
            !position of reference atom 1
            POS1(1:3) = XYZ((IMAGEOFFSET+3*(IDX1-1)+1):(IMAGEOFFSET+3*(IDX1-1)+3))
            !place new atom using fractional coordinates relative to first and last image depoending on image number
            E1 = (1.0D0*(NIMAGES+2-J1)/(NIMAGES+1))*C1 + (1.0D0*(J1-1)/(NIMAGES+1))*D1
            E2 = (1.0D0*(NIMAGES+2-J1)/(NIMAGES+1))*C2 + (1.0D0*(J1-1)/(NIMAGES+1))*D2
            E3 = (1.0D0*(NIMAGES+2-J1)/(NIMAGES+1))*C3 + (1.0D0*(J1-1)/(NIMAGES+1))*D3

            !alternative procedure that doesn't work as well
            !C = C1*B1(1:3) + C2*B2(1:3) + C3*B3(1:3)
            !D = D1*B1(1:3) + D2*B2(1:3) + D3*B3(1:3)
            !CALL GET_ROT_AXIS(B1,B2,B3,ROT)
            !CALL APPLY_ROT(C,ROT,ANGLE*(1.0D0*J1/(NIMAGES+2)),NEW2)
            !CALL NORM_VEC(NEW2,NEW,NORM)
            
            !use mixed coefficients and scale the distance properly (avoids artificial shortening)
            CALL NORM_VEC(E1*B1(1:3) + E2*B2(1:3) + E3*B3(1:3),NEW,NORM)
            NEW(1:3) = SCALE*NEW(1:3)
            XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3)) =  POS1(1:3) + NEW(1:3)
            
         END DO
      END SUBROUTINE PLACE_ATOM

      SUBROUTINE GET_LOCAL_AXIS(IDX1,IDX2,IDX3,IMAGE,B1,B2,B3)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         USE HELPER_FNCTS, ONLY: NORM_VEC, CROSS_PROD
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IDX1, IDX2, IDX3
         INTEGER, INTENT(IN) :: IMAGE
         REAL(KIND=REAL64), INTENT(OUT) :: B1(3), B2(3), B3(3)
         REAL(KIND=REAL64) :: VEC1(3), VEC2(3), VEC3(3), NORM, DOT12
         INTEGER :: IMAGEOFFSET

         IMAGEOFFSET = (IMAGE-1)*3*NATOMS
         ! VEC1 is pointing from IDX1 to IDX2
         VEC1(1:3) = XYZ((IMAGEOFFSET+3*(IDX2-1)+1):((IMAGEOFFSET+3*(IDX2-1)+3))) - XYZ((IMAGEOFFSET+3*(IDX1-1)+1):((IMAGEOFFSET+3*(IDX1-1)+3))) 
         !VEC2 is pointing from IDX1 to IDX3
         VEC2(1:3) = XYZ((IMAGEOFFSET+3*(IDX3-1)+1):((IMAGEOFFSET+3*(IDX3-1)+3))) - XYZ((IMAGEOFFSET+3*(IDX1-1)+1):((IMAGEOFFSET+3*(IDX1-1)+3))) 

         !B1 (first base vector) is the normed VEC1
         CALL NORM_VEC(VEC1,B1,NORM)
         !to get the second base vector (B2) we use the orthogonal component of VEC2 to VEC1
         DOT12 = VEC1(1)*VEC2(1) + VEC1(2)*VEC2(2) + VEC1(3)*VEC2(3)
         VEC2(1:3) = VEC2(1:3) - DOT12*B1(1:3)
         CALL NORM_VEC(VEC2,B2,NORM)
         !The final base vector is the cross product of B1 and B2
         B3 = CROSS_PROD(B1,B2)
      END SUBROUTINE GET_LOCAL_AXIS

      SUBROUTINE GET_LOCAL_AXIS2(IDX1,IDX2,IDX3,IDX4,IMAGE,B1,B2,B3)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         USE HELPER_FNCTS, ONLY: NORM_VEC, CROSS_PROD, GS_PROJECTION, DOTP
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IDX1, IDX2, IDX3, IDX4
         INTEGER, INTENT(IN) :: IMAGE
         REAL(KIND=REAL64), INTENT(OUT) :: B1(3), B2(3), B3(3)
         REAL(KIND=REAL64) :: VEC1(3), VEC2(3), VEC3(3), NORM, DOT12, NORM1, NORM2, CP1(3), CP2(3), VOL
         REAL(KIND=REAL64), PARAMETER :: VOLCUT=0.1D0
         INTEGER :: IMAGEOFFSET

         IMAGEOFFSET = (IMAGE-1)*3*NATOMS
         ! VEC1 is pointing from IDX1 to IDX2
         VEC1(1:3) = XYZ((IMAGEOFFSET+3*(IDX2-1)+1):((IMAGEOFFSET+3*(IDX2-1)+3))) - XYZ((IMAGEOFFSET+3*(IDX1-1)+1):((IMAGEOFFSET+3*(IDX1-1)+3))) 
         !VEC2 is pointing from IDX1 to IDX3
         VEC2(1:3) = XYZ((IMAGEOFFSET+3*(IDX3-1)+1):((IMAGEOFFSET+3*(IDX3-1)+3))) - XYZ((IMAGEOFFSET+3*(IDX1-1)+1):((IMAGEOFFSET+3*(IDX1-1)+3))) 
         !VEC3 is pointing from IDX1 to IDX4
         VEC3(1:3) = XYZ((IMAGEOFFSET+3*(IDX4-1)+1):((IMAGEOFFSET+3*(IDX4-1)+3))) - XYZ((IMAGEOFFSET+3*(IDX1-1)+1):((IMAGEOFFSET+3*(IDX1-1)+3))) 
         !debugging output
         CALL NORM_VEC(VEC1,B1,NORM)
         CALL NORM_VEC(VEC2,B2,NORM)
         CALL NORM_VEC(VEC3,B3,NORM)
         CALL NORM_VEC(CROSS_PROD(B1,B2),CP1,NORM1)
         CALL NORM_VEC(CROSS_PROD(B1,B3),CP2,NORM2)
         VOL = DOTP(3,B1,CROSS_PROD(B2,B3))
         
         !B1 (first base vector) is the normed VEC1
         CALL NORM_VEC(VEC1,B1,NORM)
         !to get the second base vector (B2) we use the orthogonal component of VEC2 to VEC1
         CALL NORM_VEC(VEC2-GS_PROJECTION(VEC2,VEC1),B2,NORM)
         !The final base vector is the orthogonal component of VEC3 to VEC1 and VEC2
         IF (VOL.GT.VOLCUT) THEN
            !WRITE(*,*) "local_axis> Local volume indicates this is not a planer set of axis, using Gram-Schmidt projection"
            CALL NORM_VEC(VEC3-GS_PROJECTION(VEC3,VEC1)-GS_PROJECTION(VEC3,VEC2),B3,NORM)
         ELSE
            !WRITE(*,*) "local_axis> Local volume indicates this is a planer set of axis, using normal vector"
            CALL NORM_VEC(CROSS_PROD(B1,B2),B3,NORM)
         END IF
      END SUBROUTINE GET_LOCAL_AXIS2     

      SUBROUTINE GET_ROT_ANGLE(C,D,ROT,ANGLE)
         USE HELPER_FNCTS, ONLY: GS_PROJECTION, NORM_VEC, DOTP, CROSS_PROD, EUC_NORM
         REAL(KIND = REAL64), INTENT(IN) :: C(3), D(3), ROT(3)
         REAL(KIND = REAL64), INTENT(OUT) :: ANGLE
         REAL(KIND = REAL64) :: CP(3), DP(3)
         REAL(KIND = REAL64) :: NORM, COSTH, SINTH

         CALL NORM_VEC(C-GS_PROJECTION(C,ROT), CP, NORM)
         CALL NORM_VEC(D-GS_PROJECTION(D,ROT), DP, NORM)
         
         COSTH = DOTP(3,CP,DP)
         SINTH = EUC_NORM(CROSS_PROD(CP,DP))

         ANGLE = ATAN2(SINTH,COSTH)
      END SUBROUTINE GET_ROT_ANGLE

      SUBROUTINE APPLY_ROT(C,ROT,ANGLE,NEW)
         REAL(KIND = REAL64), INTENT(IN) :: C(3), ROT(3), ANGLE
         REAL(KIND = REAL64), INTENT(OUT) :: NEW(3)
         REAL(KIND = REAL64) :: CX, CY, CZ, RX, RY, RZ
         REAL(KIND = REAL64) :: COSTH, SINTH, DP

         CX = C(1); CY = C(2); CZ = C(3)
         RX = ROT(1); RY = ROT(2); RZ = ROT(3)
         COSTH = DCOS(ANGLE)
         SINTH = DSIN(ANGLE)

         DP = CX*RX + CY*RY + CZ*RZ

         NEW(1) = RX*DP*(1.0D0-COSTH) + CX*COSTH + (-RZ*CY + RY*CZ)*SINTH
         NEW(2) = RY*DP*(1.0D0-COSTH) + CY*COSTH + ( RZ*CX - RX*CZ)*SINTH
         NEW(3) = RZ*DP*(1.0D0-COSTH) + CZ*COSTH + (-RY*CX + RX*CY)*SINTH
      END SUBROUTINE APPLY_ROT
      
      SUBROUTINE SCALE_NEW_POSITION(NEW,C,D,SCALED)
         USE HELPER_FNCTS, ONLY: EUC_NORM, NORM_VEC
         REAL(KIND = REAL64), INTENT(IN) :: C(3), D(3), NEW(3)
         REAL(KIND = REAL64), INTENT(OUT) :: SCALED(3)
         REAL(KIND = REAL64) :: NORM, NORM1, NORM2

         CALL NORM_VEC(NEW,SCALED,NORM)

         NORM1 = EUC_NORM(C)
         NORM2 = EUC_NORM(D)

         SCALED = SCALED * 0.5 * (NORM1+NORM2)
      END SUBROUTINE SCALE_NEW_POSITION

      SUBROUTINE GET_ROT_AXIS(B1,B2,B3,ROT)
         USE HELPER_FNCTS, ONLY: NORM_VEC
         REAL(KIND = REAL64), INTENT(IN) :: B1(3), B2(3), B3(3)
         REAL(KIND = REAL64), INTENT(OUT) :: ROT(3)
         REAL(KIND = REAL64) :: NORM

         CALL NORM_VEC(B1+B2+B3,ROT,NORM)
      END SUBROUTINE GET_ROT_AXIS

      SUBROUTINE CHECK_NACTIVE()
         USE INTERPOLATION_KEYS, ONLY: NACTIVE, ATOMACTIVE
         USE QCIKEYS, ONLY: NATOMS
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
         USE QCIKEYS, ONLY: NATOMS, DEBUG, QCIREPCUT, QCIINTREPMINSEP
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE  
         USE REPULSION
         USe MOD_INTCOORDS, ONLY: XSTART, XFINAL
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         IMPLICIT NONE  
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER :: J1, J2
         LOGICAL :: ISCONSTRAINED  
         REAL(KIND = REAL64) :: DF, DS , DMIN      
         
         WRITE(*,*) " Update_repulsion> Before update: ", NREPULSIVE

         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE !ignore atoms that are not active
            IF (ABS(J1-NEWATOM).LE.QCIINTREPMINSEP) CYCLE !no repulsions if atoms are close in sequence
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
               IF (NREPULSIVE.GT.NREPCURR) CALL DOUBLE_ALLOC_REP()
               REPI(NREPULSIVE)=J1
               REPJ(NREPULSIVE)=NEWATOM
               REPCUT(NREPULSIVE)=DMIN
               !WRITE(*,*) " new repulsion: ",NREPULSIVE, " repi: ", REPI(NREPULSIVE), " repj: ", REPJ(NREPULSIVE), " repcut: ", REPCUT(NREPULSIVE)
            END IF
         END DO
         WRITE(*,*) " Update_repulsion> After update: ", NREPULSIVE

      END SUBROUTINE UPDATE_REPULSIONS


      SUBROUTINE UPDATE_CONSTRAINTS(NEWATOM,NCONNEWATOM,BESTCONDIST,BESTCONIDX)
         USE QCIKEYS, ONLY: NATOMS, DEBUG
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, ATOMACTIVE, NCONSTRAINTON
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREF, MAXCONUSE
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         REAL(KIND=REAL64), INTENT(OUT) :: BESTCONDIST(NATOMS)  !sorted list of average distances
         INTEGER, INTENT(OUT)  :: BESTCONIDX(NATOMS)            !associated list of atom ids
         INTEGER, INTENT(OUT) :: NCONNEWATOM                    !number of atoms found that are within the cut offs and active     
         INTEGER :: J1, J2, J3, ATOM1, ATOM2

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
               WRITE(*,*) " New constraint: ", ATOM1, ATOM2, "new atom: ", NEWATOM, " distance: ", CONDISTREF(J1)
            END IF              
         END DO

         
         IF (DEBUG) THEN
            WRITE(*,*) " get_constraints_by_dist> Constraints including new atom by average reference distance:"
            WRITE(*,'(10G12.4)') BESTCONDIST(1:MIN(10,NCONNEWATOM))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10I6)') BESTCONIDX(1:MIN(10,NCONNEWATOM))
         END IF

         ! Turning on constraints identified
         DO J1=1,MIN(MAXCONUSE,NCONNEWATOM)
            DO J2=1,NCONSTRAINT
               ATOM1 = CONI(J2)
               ATOM2 = CONJ(J2)
               IF (((ATOM1.EQ.NEWATOM).AND.(ATOM2.EQ.BESTCONIDX(J1))).OR.((ATOM2.EQ.NEWATOM).AND.(ATOM1.EQ.BESTCONIDX(J1)))) THEN
                  CONACTIVE(J2) = .TRUE.
                  NCONSTRAINTON = NCONSTRAINTON + 1
                  IF (DEBUG) THEN
                     WRITE(*,'(A,I6,A,2I6)') ' addatom> Turning on constraint ',J2,' for atoms ', ATOM1, ATOM2
                  END IF
               END IF
            END DO
         END DO
      END SUBROUTINE UPDATE_CONSTRAINTS

      SUBROUTINE GET_SECONDARY_CONSTRAINTS(NEWATOM,NCONSTR,CONSTIDX,NSECCONSTR,SECCONDIST,SECCONIDX)
         USE QCIKEYS, ONLY: NATOMS, DEBUG
         USE AMBER_CONSTRAINTS, ONLY: NBOND, BONDS
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, ATOMACTIVE
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ
         USE HELPER_FNCTS, ONLY: DISTANCE_TWOATOMS
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM                         !new atom id
         INTEGER, INTENT(IN) :: CONSTIDX(NATOMS)                !atoms ids for primary constraints
         INTEGER, INTENT(IN) :: NCONSTR                         !number of primary constraints        
         REAL(KIND=REAL64), INTENT(OUT) :: SECCONDIST(NATOMS)   !sorted list of average distances for secondary constraints
         INTEGER, INTENT(OUT) :: SECCONIDX(NATOMS)              !associated list of atom ids
         INTEGER, INTENT(OUT) :: NSECCONSTR                     !number of atoms found for secondary constraints
         INTEGER :: PRIMARYIDX, IDX1, IDX2
         REAL(KIND = REAL64) :: DF, DS, AVD, DEV
         LOGICAL :: SECCONSTATOMS(NATOMS)
         INTEGER :: J1, J2, J3
         INTEGER :: NDUMMY
         LOGICAL :: USE_BOND_DIST
         INTEGER :: BOND_DIST(NATOMS)
         LOGICAL :: BONDED_ACTIVE(NATOMS)
         INTEGER :: NBACTIVE

         !TODO: bonded bits can be removed again ...

         IF (NBOND.GT.0) THEN
            !set every atom to be not connected
            USE_BOND_DIST = .TRUE.
            BONDED_ACTIVE(1:NATOMS) = .FALSE.
            BOND_DIST(1:NATOMS) = -1
            !activate new atom
            BONDED_ACTIVE(NEWATOM) = .TRUE.
            BOND_DIST(NEWATOM) = 0
            NBACTIVE = 1
            DO WHILE (NBACTIVE.LT.NATOMS)
               DO J1=1,NBOND
                  IDX1 = BONDS(J1,1)
                  IDX2 = BONDS(J1,2)
                  IF (BONDED_ACTIVE(IDX1)) THEN
                     IF (BONDED_ACTIVE(IDX2)) THEN
                        IF (BOND_DIST(IDX1).GT.BOND_DIST(IDX2)+1) THEN
                           BOND_DIST(IDX1) = BOND_DIST(IDX2)+1
                        ELSE IF (BOND_DIST(IDX2).GT.BOND_DIST(IDX1)+1) THEN
                           BOND_DIST(IDX2) = BOND_DIST(IDX1)+1 
                        END IF 
                     ELSE
                        BOND_DIST(IDX2) = BOND_DIST(IDX1)+1 
                        BONDED_ACTIVE(IDX2) = .TRUE.
                        NBACTIVE = NBACTIVE + 1
                     END IF
                  ELSE IF (BONDED_ACTIVE(IDX2)) THEN
                     BOND_DIST(IDX1) = BOND_DIST(IDX2)+1 
                     BONDED_ACTIVE(IDX1) = .TRUE.                     
                     NBACTIVE = NBACTIVE + 1
                  END IF
               END DO
            END DO
         ELSE
            USE_BOND_DIST = .FALSE.
         END IF
         SECCONSTATOMS(1:NATOMS) = .FALSE.
         DO J1=1,NCONSTR
            !get id of primary atom
            PRIMARYIDX = CONSTIDX(J1)
            DO J2=1,NCONSTRAINT
               !we only consider active constraints
               IF (.NOT.CONACTIVE(J2)) CYCLE
               IDX1 = CONI(J2)
               IDX2 = CONJ(J2)
               !ignore constraints that don't contain the primary constraint atom
               IF ((IDX1.NE.PRIMARYIDX).AND.(IDX2.NE.PRIMARYIDX)) CYCLE
               !ignore the constraint that contains the primary and new atom
               IF (((IDX1.EQ.PRIMARYIDX).AND.(IDX2.EQ.NEWATOM)).OR.((IDX2.EQ.PRIMARYIDX).AND.(IDX1.EQ.NEWATOM))) CYCLE
               ! if we get to this point, we have an active constraint with the primary atom that is not the constraint to NEWATOM
               IF (IDX1.EQ.PRIMARYIDX) THEN
                  SECCONSTATOMS(IDX2) = .TRUE.
               ELSE
                  SECCONSTATOMS(IDX1) = .TRUE.
               END IF
            END DO
         END DO

         NSECCONSTR = 0
         SECCONDIST(1:NATOMS) = 1.0D100
         SECCONIDX(1:NATOMS) = -1

         DO J1=1,NATOMS
            !ignore atoms not selected in the previous loop
            IF (.NOT.SECCONSTATOMS(J1)) CYCLE
            !add a sanity check for active atoms
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            NSECCONSTR = NSECCONSTR + 1
            IF (USE_BOND_DIST) THEN
               AVD = 1.0D0*BOND_DIST(J1)
            ELSE
               CALL DISTANCE_TWOATOMS(NATOMS, XSTART, NEWATOM, J1, DS)
               CALL DISTANCE_TWOATOMS(NATOMS, XFINAL, NEWATOM, J1, DF)
               AVD=0.5*(DS+DF)
               DEV=0.5*ABS(DS-DF)
            END IF
            ! go through sorting loop to mke sure we maintain an ordered list
            DO J2=1,NSECCONSTR
               IF (AVD.LT.SECCONDIST(J2)) THEN
                  DO J3=NSECCONSTR,J2+1,-1
                     SECCONDIST(J3)=SECCONDIST(J3-1)
                     SECCONIDX(J3)=SECCONIDX(J3-1)
                  END DO
                  SECCONDIST(J2) = AVD
                  SECCONIDX(J2) = J1
                  EXIT
               END IF
            END DO
         END DO

         IF (DEBUG) THEN
            WRITE(*,*) " get_secondary_constraint> Closest secondary constraint atoms to new atom by average endpoint distance:"
            WRITE(*,'(10G12.4)') SECCONDIST(1:MIN(10,NSECCONSTR))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10I6)') SECCONIDX(1:MIN(10,NSECCONSTR))          
         END IF
      END SUBROUTINE GET_SECONDARY_CONSTRAINTS

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
         REAL(KIND = REAL64) :: DF, DS, AVD, DEV
         INTEGER :: J1, J2, J3, STARTIDX, ENDIDX
         REAL(KIND = REAL64) :: DISTDEV(NATOMS)

         BESTDIST(1:NATOMS) = 1.0D100
         DISTDEV(1:NATOMS) = -1.0D0
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
            DEV=0.5*ABS(DS-DF)
            NDISTNEWATOM = NDISTNEWATOM + 1
            ! go through sorting loop to mke sure we maintain an ordered list
            DO J2=1,NDISTNEWATOM
               IF (AVD.LT.BESTDIST(J2)) THEN
                  DO J3=NDISTNEWATOM,J2+1,-1
                     BESTDIST(J3)=BESTDIST(J3-1)
                     DISTDEV(J3)=DISTDEV(J3-1)
                     BESTIDX(J3)=BESTIDX(J3-1)
                  END DO
                  BESTDIST(J2) = AVD
                  DISTDEV(J2) = DEV
                  BESTIDX(J2) = J1
                  EXIT
               END IF
            END DO
         END DO

         IF (DEBUG) THEN
            WRITE(*,*) " get_atoms_by_dist> Closest atoms to new atom by average endpoint distance:"
            WRITE(*,'(10G12.4)') BESTDIST(1:MIN(10,NDISTNEWATOM))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10I6)') BESTIDX(1:MIN(10,NDISTNEWATOM))
            WRITE(*,*) "                    Distance deviations:"
            WRITE(*,'(10G12.4)') DISTDEV(1:MIN(10,NDISTNEWATOM))            
         END IF
      END SUBROUTINE GET_ATOMS_BY_DISTANCE

      SUBROUTINE FIND_NEXT_ATOM(CHOSENACID,ACID,NEWATOM,NCONTOACT,SHORTESTCON)
         USE QCIKEYS, ONLY: QCILINEART, INLINLIST, ATOMS2RES, QCIDOBACK, ISBBATOM, QCIUSEGROUPS
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE, CONACTIVE, NACTIVE
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONDISTREF, CONI, CONJ
         USE AMBER_CONSTRAINTS, ONLY: CURRENTLY_ADDING_GROUP, CHECK_IN_GROUP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ACID
         LOGICAL, INTENT(IN) :: CHOSENACID
         INTEGER, INTENT(OUT) :: NEWATOM
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
            IF (QCILINEART.AND.(.NOT.INLINLIST(IDXINACTIVE))) CYCLE
            ! 2. Is CHOSENACID set, and is the inactive atom in the residue to be added?
            IF (CHOSENACID.AND.(.NOT.(ATOMS2RES(IDXINACTIVE).EQ.ACID))) CYCLE
            ! 3. Are we adding backbone atoms, and is the inactive atom a backbone atom?
            IF (QCIDOBACK.AND.(.NOT.BBDONE).AND.(.NOT.ISBBATOM(IDXINACTIVE))) CYCLE
            ! 4. Is the atom in the current group?
            IF (QCIUSEGROUPS.AND.CURRENTLY_ADDING_GROUP) THEN
               IF (.NOT.CHECK_IN_GROUP(IDXINACTIVE)) CYCLE
            END IF

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
            WRITE(*,*) NCONTOACTIVE
            WRITE(*,*) "Nactive: ", NACTIVE
            WRITE(*,*) "QCI linear: ", QCILINEART
            WRITE(*,*) NEWATOM, NCONTOACT
            WRITE(*,*) " find_next_atom> Error - new active atom not set, NEWATOM: ", NEWATOM, " NCONTOACT: ", NCONTOACT
            CALL INT_ERR_TERMINATE()
         END IF
      END SUBROUTINE FIND_NEXT_ATOM


      SUBROUTINE CREATE_NCONTOACTIVE_LIST(INVDTOACTIVE,NBEST)
         USE QCIKEYS, ONLY: NATOMS
         USE INTERPOLATION_KEYS, ONLY: NACTIVE, ATOMACTIVE, CONACTIVE
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREF
         INTEGER, INTENT(OUT) :: NBEST ! maximum number of constraintson inactive atom
         REAL(KIND = REAL64), INTENT(OUT) :: INVDTOACTIVE(1:NATOMS)  !inverse distance
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
            IF (NCONTOACTIVE(ATOM1).GT.NBEST) NBEST = NCONTOACTIVE(ATOM1)
            IF (NCONTOACTIVE(ATOM2).GT.NBEST) NBEST = NCONTOACTIVE(ATOM2)
         END DO
         
      END SUBROUTINE CREATE_NCONTOACTIVE_LIST

      SUBROUTINE CHECK_BBLIST()
         USE QCIKEYS, ONLY: NATOMS, DEBUG, ISBBATOM
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         IMPLICIT NONE
         INTEGER :: J1
         DO J1=1,NATOMS
            IF (ISBBATOM(J1)) THEN
               IF (.NOT.ATOMACTIVE(J1)) THEN
                  ! if there are inactive BB atoms, return
                  RETURN
               END IF
            ENDIF
         ENDDO
         BBDONE = .TRUE.
         IF (DEBUG) WRITE(*,*) " check_bblist> All backbone atoms are active"
      END SUBROUTINE CHECK_BBLIST


END MODULE ADDINGATOM
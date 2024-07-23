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

      SUBROUTINE ADDATOM()
         USE MOD_INTCOORDS, ONLY: XYZ, EEE, GGG, RMS
         USE QCIKEYS, ONLY: QCIDOBACK, QCIADDACIDT, NATOMS, NIMAGES, DEBUG, QCILINEART, INLINLIST, &
                            QCITRILATERATION, QCIDOBACKALL, ATOMS2RES, ISBBATOM, CHECKCONINT
         USE REPULSION, ONLY: NNREPULSIVE, NREPULSIVE, CHECKREP
         USE CONSTR_E_GRAD, ONLY: CONGRAD1, CONGRAD2
         USE QCI_LINEAR, ONLY: NQCILINEAR
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, NACTIVE, TURNONORDER, ATOMACTIVE
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
         REAL(KIND=REAL64) :: BESTDIST(NATOMS)              !list of sorted average distance to newatom
         INTEGER :: BESTIDX(NATOMS), NDISTNEWATOM           !associated ids and total number found
         REAL(KIND=REAL64) :: BESTCONDIST(NATOMS)           !list of sorted average constraint distances to newatom
         INTEGER :: BESTCONIDX(NATOMS), NCONNEWATOM         !associated ids and total number found
         LOGICAL :: ADDEDTHISCYCLE                          !have we added an atom this cycle?
         REAL(KIND = REAL64), PARAMETER :: FRAC = 1.0D0     !fraction for linear interpolation
         INTEGER :: THISIMAGE, NEWATOMOFFSET, CONONEOFFSET, ENDPOINT !offsets used to simplify linear interpolation
         REAL(KIND = REAL64) :: STARTWEIGHT, ENDWEIGHT      !weights of endpoints in linear interpolation
         REAL(KIND = REAL64) :: ELIN, ECON, EDIST, ETOTAL   !linear, constraint and distance measure-based energies
         REAL(KIND = REAL64) :: XLIN(NIMAGES,3), XCON(NIMAGES,3), XDIST(NIMAGES,3) !saved coordinates for interpolations
         INTEGER :: J1
      
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
                  ACID = ATOMS2RES(NEXTATOM)
                  CHOSENACID = .TRUE.
               END IF
            END IF
            
            WRITE(*,'(A,I6,A,I4,A,I4)') "  addatom> Adding atom ", NEXTATOM, ", which has ", NCONTOACT, " constraints to active set out of maximum ", NMAXCON
            WRITE(*,'(A,F8.4)') "           Shortest distance constraint to active set is: ", SHORTESTCON

            ! The interpolation for the new atom relies on a local axis system formed by three atoms.
            ! We look for a sorted list, according to how well the end point distace is preserved.
            ! We sort by shortest average distance to avoid distant atoms having accidentally well preserved distance.
            CALL GET_ATOMS_BY_DISTANCE(NEXTATOM,NDISTNEWATOM,BESTDIST,BESTIDX)
            ! Now update the constraints including the new atom, and get list of constraints ordered by distance
            CALL UPDATE_CONSTRAINTS(NEXTATOM,NCONNEWATOM,BESTCONDIST,BESTCONIDX)
            !update the repulsions
            CALL UPDATE_REPULSIONS(NEXTATOM)

            !activate new atom
            ATOMACTIVE(NEXTATOM)=.TRUE.
            NACTIVE=NACTIVE+1 
            TURNONORDER(NACTIVE)=NEXTATOM
            !check consistency
            CALL CHECK_NACTIVE()

            !now we need to actually add the atom
            ADDEDTHISCYCLE = .FALSE.
            EDIST = 1.0D100
            ELIN = 1.0D100
            ECON = 1.0D100
            !if we have three or more constraints, we use them and construct a local axis system
            IF (NCONNEWATOM.GE.3) THEN
               ADDEDTHISCYCLE = .TRUE.
               CALL PLACE_ATOM(NEXTATOM,BESTCONIDX)
               IF (QCITRILATERATION) THEN
                  CALL TRILATERATE_ATOMS(NEXTATOM,BESTCONIDX,BESTCONDIST)
               END IF
               ! before we continue check repulsion neighbour list
               CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
               ! call congrad routine
               IF (CHECKCONINT) THEN
                  CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
               ELSE
                  CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
               END IF
               ECON = ETOTAL
               DO J1=1,NIMAGES
                  NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                  XCON(J1,1:3) = XYZ((NEWATOMOFFSET+1):(NEWATOMOFFSET+3)) 
               END DO
            END IF
            !if we don't have enough constraints, use the closest active atoms instead to construct the axis system
            IF ((.NOT.ADDEDTHISCYCLE).AND.(NDISTNEWATOM.GE.3)) THEN
               ADDEDTHISCYCLE = .TRUE.
               CALL PLACE_ATOM(NEXTATOM,BESTIDX)
               IF (QCITRILATERATION) THEN
                  CALL TRILATERATE_ATOMS(NEXTATOM,BESTIDX,BESTDIST)
               END IF
               ! before we continue check repulsion neighbour list
               !TODO: check checkrep here and in the original version match up
               CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
               ! call congrad routine
               IF (CHECKCONINT) THEN
                  CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
               ELSE
                  CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
               END IF
               EDIST = ETOTAL
               DO J1=1,NIMAGES
                  NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                  XDIST(J1,1:3) = XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) 
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
                  CONONEOFFSET = 3*(BESTCONIDX(1)-1)
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
               IF (CHECKCONINT) THEN
                  CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
               ELSE
                  CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
               END IF
               ELIN = ETOTAL
               DO J1=1,NIMAGES
                  NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                  XLIN(J1,1:3) = XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) 
               END DO
            END IF
            
            !select which interpolation is used based on energy if multiple are used (we shouldn't have this case except for QCIlinear)
            IF (QCILINEART) THEN
               WRITE(*,*) " addatom> Using linear interpolation for new atom ", NEXTATOM, " from linear list"
            ELSE 
               IF ((ELIN.LT.ECON).AND.(ELIN.LT.EDIST)) THEN
                  WRITE(*,*) " addatom> Using linear interpolation for new atom ", NEXTATOM
               ELSE IF ((ECON.LT.ELIN).AND.(ECON.LT.EDIST)) THEN
                  WRITE(*,*) " addatom> Using interpolation from preserved distances for new atom ", NEXTATOM
                  ETOTAL = ECON
                  DO J1=1,NIMAGES
                     NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                     XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) = XCON(J1,1:3)
                  END DO
               ELSE IF ((EDIST.LT.ELIN).AND.(EDIST.LT.ECON)) THEN
                  WRITE(*,*) " addatom> Using interpolation from closest atoms for new atom ", NEXTATOM
                  ETOTAL = EDIST
                  DO J1=1,NIMAGES
                     NEWATOMOFFSET = J1*3*NATOMS + 3*(NEXTATOM-1)
                     XYZ(NEWATOMOFFSET+1:NEWATOMOFFSET+3) = XDIST(J1,1:3)
                  END DO
               END IF
            END IF

            NADDED = NADDED + 1
            IF (NADDED.LT.NTOADD) CYCLE !if we already know we need to add more, skip ahead and add the next atom

            !we should not have more atoms to add unless we hit certain criteria
            MORETOADD = .FALSE.
            ! if we add one residue at a time, check whether we added all of them
            IF (QCIADDACIDT.AND.(.NOT.QCIDOBACK)) THEN
               IF (NCONNEWATOM.LT.3) THEN
                  WRITE(*,*) " addatom> Addad atom constraints ", NCONNEWATOM, " is less than 3 - not attempting to add another atom."
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
         END DO

         CALL CHECKREP(XYZ,NNREPSAVE,NREPSAVE+1)
         ! call congrad routine
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
         ELSE
            CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
         END IF
         !we are done with QCIlinear, so set it to false
         QCILINEART = .FALSE.
      END SUBROUTINE ADDATOM

      SUBROUTINE TRILATERATE_ATOMS(NEWATOM,CONIDXLIST,CONDISTLIST)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ         
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER, INTENT(IN) :: CONIDXLIST(NATOMS)        
         REAL(KIND=REAL64), INTENT(IN) :: CONDISTLIST(NATOMS)
         REAL(KIND=REAL64) :: P1(3), P2(3), P3(3), R1, R2, R3
         REAL(KIND=REAL64) :: SOL1(3), SOL2(3), PREV(3), D1SQ, D2SQ
         INTEGER :: IMAGEOFFSET
         INTEGER :: N1, N2, N3, IDX1, IDX2, IDX3, J1
         LOGICAL :: FTEST

         !set initial guess to the best three constrained atoms
         N1=1; N2=2; N3=3
         IDX1 = CONIDXLIST(N1); IDX2 = CONIDXLIST(N2); IDX3 = CONIDXLIST(N3)

         DO J1=2,NIMAGES+1
            IMAGEOFFSET = (J1-1)*3*NIMAGES
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

      SUBROUTINE PLACE_ATOM(NEWATOM,CONIDXLIST)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NEWATOM
         INTEGER, INTENT(IN) :: CONIDXLIST(NATOMS)

         REAL(KIND=REAL64) :: VEC(3), B1(3), B2(3), B3(3), POS1(3)
         REAL(KIND=REAL64) :: C1, C2, C3
         INTEGER :: IMAGEOFFSET
         INTEGER :: N1, N2, N3, IDX1, IDX2, IDX3
         INTEGER :: J1

         !set initial guess to the best three constrained atoms
         N1=1; N2=2; N3=3
         IDX1 = CONIDXLIST(N1); IDX2 = CONIDXLIST(N2); IDX3 = CONIDXLIST(N3)
         IF (DEBUG) THEN
            WRITE(*,*) " place_atom> Use the closest three constrained active atoms as initial guess"
            WRITE(*,*) "             New atom: ", NEWATOM, "Closest active atoms: ", CONIDXLIST(N1:N3)
         END IF

         !Setting up local axis system
         !VEC is pointing from IDX1 to NEWATOM
         VEC(1:3) = XYZ((3*(NEWATOM-1)+1):((3*(NEWATOM-1)+3))) - XYZ((3*(IDX1-1)+1):((3*(IDX1-1)+3)))
         CALL GET_LOCAL_AXIS(IDX1,IDX2,IDX3,1,B1,B2,B3)
         !get relative coordinates of NEWATOM in B1,B2,B3 axis system
         C1 = DOT_PRODUCT(VEC,B1)
         C2 = DOT_PRODUCT(VEC,B2)
         C3 = DOT_PRODUCT(VEC,B3)

         !iterate over images and place NEWATOM
         DO J1=2,NIMAGES+1
            !get B1,B2,B3 for current image
            CALL GET_LOCAL_AXIS(IDX1,IDX2,IDX3,J1,B1,B2,B3)
            IMAGEOFFSET = (J1-1)*3*NATOMS
            !position of reference atom 1
            POS1(1:3) = XYZ((IMAGEOFFSET+3*(IDX1-1)+1):(IMAGEOFFSET+3*(IDX1-1)+3))
            !place new atom usig fractional coordinates from start image
            XYZ((IMAGEOFFSET+3*(NEWATOM-1)+1):(IMAGEOFFSET+3*(NEWATOM-1)+3)) =  POS1(1:3) + C1*B1(1:3) + C2*B2(1:3) + C3*B3(1:3) 
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
         !to get the second base vector (B2) we use the orthogonal component of VEC2 to B1
         DOT12 = B1(1)*VEC2(1) + B1(2)*VEC2(2) + B1(3)*VEC2(3)
         VEC2(1:3) = VEC2(1:3) - DOT12*VEC1(1:3)
         CALL NORM_VEC(VEC2,B2,NORM)
         !The final base vector is the cross product of B1 and B2
         B3 = CROSS_PROD(B1,B2)
      END SUBROUTINE GET_LOCAL_AXIS

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
            END IF
         END DO

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
         REAL(KIND = REAL64) :: DF, DS, AVD
         INTEGER :: J1, J2, J3, STARTIDX, ENDIDX

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
            WRITE(*,'(10G12.4)') BESTDIST(1:1:MIN(10,NDISTNEWATOM))
            WRITE(*,*) "                    Average distances:"
            WRITE(*,'(10I6)') BESTIDX(1:1:MIN(10,NDISTNEWATOM))
         END IF
      END SUBROUTINE GET_ATOMS_BY_DISTANCE

      SUBROUTINE FIND_NEXT_ATOM(CHOSENACID,ACID,NEWATOM,NCONTOACT,SHORTESTCON)
         USE QCIKEYS, ONLY: QCILINEART, INLINLIST, ATOMS2RES, QCIDOBACK, ISBBATOM
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE, CONACTIVE, NACTIVE
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONDISTREF, CONI, CONJ
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
            WRITE(*,*) "Nactive; ", NACTIVE
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
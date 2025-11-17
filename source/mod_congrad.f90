! Module containing the congrad routines
!
! version 1 (congrad1) - tests for internal minimum in repulsions only
! version 2 (congrad2) - tests for internal minimum in repulsion and constraints
MODULE CONSTR_E_GRAD
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: MAXCONIMAGE, MAXREPIMAGE, MAXSPRIMAGE      ! image with highest constraint / repulsion
   INTEGER :: MAXCONSTR, MAXREP                          ! index for worst constraint / repulsion
   REAL(KIND=REAL64) :: CONVERGECONTEST, CONVERGEREPTEST ! energy of that term
   REAL(KIND=REAL64) :: FCONTEST, FREPTEST
   REAL(KIND=REAL64) :: EMAXSPR
   REAL(KIND=REAL64) :: FCONMAX, FREPMAX                 ! maximum gradient
   INTEGER :: CALLN = 0
   CONTAINS

      SUBROUTINE CONGRAD(ETOTAL, XYZ, GGG, EEE, RMS)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, CHECKCONINT
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2)             ! energy for each image
         REAL(KIND = REAL64), INTENT(OUT) :: ETOTAL                    ! overall energy
         REAL(KIND = REAL64), INTENT(OUT) :: RMS                       ! total force

         ! call correct congrad routine
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
         ELSE
            CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
         END IF
      END SUBROUTINE CONGRAD


      SUBROUTINE CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, QCICONSTRREP, KINT, QCIFREEZET, QCIFROZEN, INTCONSTRAINTDEL, &
                            USEDIHEDRALCONST
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2)             ! energy for each image
         REAL(KIND = REAL64), INTENT(OUT) :: ETOTAL                    ! overall energy
         REAL(KIND = REAL64), INTENT(OUT) :: RMS                       ! total force
        
         REAL(KIND = REAL64) :: ECON, EREP, ESPR, EDIH      ! QUERY: should these be globals in keys?
         REAL(KIND = REAL64) :: EEEC(NIMAGES+2), GGGC(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEER(NIMAGES+2), GGGR(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEES(NIMAGES+2), GGGS(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEED(NIMAGES+2), GGGD(3*NATOMS*(NIMAGES+2))
         INTEGER :: J1, J2

         CALLN = CALLN + 1

         ! initiate some variables
         EEE(1:NIMAGES+2)=0.0D0; EEEC(1:NIMAGES+2)=0.0D0; EEER(1:NIMAGES+2)=0.0D0; EEES(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0; GGGC(1:(3*NATOMS)*(NIMAGES+2))=0.0D0; GGGR(1:(3*NATOMS)*(NIMAGES+2))=0.0D0; GGGS(1:(3*NATOMS)*(NIMAGES+2))=0.0D0
         ECON = 0.0D0; EREP = 0.0D0; ESPR = 0.0D0


         ! QUERY: what is INTCONSTRAINTDEL? seems like a scaling for the potential
         IF (.NOT.(INTCONSTRAINTDEL.EQ.0.0D0)) THEN
            CALL GET_CONSTRAINT_E_NOINTERNAL(XYZ,GGGC,EEEC,ECON)
         END IF
         IF (.NOT.(QCICONSTRREP.EQ.0.0D0)) THEN
            CALL GET_REPULSION_E(XYZ,GGGR,EEER,EREP)
         END IF
         IF (.NOT.(KINT.EQ.0.0D0)) THEN
            CALL GET_SPRING_E(XYZ, GGGS, EEES, ESPR)
         END IF
         IF (USEDIHEDRALCONST) THEN
            CALL GET_DIH_CON_E(XYZ,GGGD,EEED,EDIH)
         END IF

         ! add all contributions
         EEE = EEEC + EEER + EEES + EEED
         GGG = GGGC + GGGR + GGGS + GGGD
         ! freeze atoms that should be frozen
         IF (QCIFREEZET) THEN
            DO J1=2,NIMAGES+1
               DO J2=1,NATOMS
                  IF (QCIFROZEN(J2)) THEN
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
                  END IF
               END DO
            END DO
         END IF

         ! Set gradients to zero for start and finish images.
         GGG(1:(3*NATOMS))=0.0D0
         GGG((NIMAGES+1)*(3*NATOMS)+1:(NIMAGES+2)*(3*NATOMS))=0.0D0

         ! get the RMS force
         RMS=0.0D0
         DO J1=2,NIMAGES+1
            DO J2=1,3*NATOMS
               RMS = RMS + GGG((3*NATOMS)*(J1-1)+J2)**2
            END DO
         END DO
         RMS = SQRT(RMS/(3*NATOMS*NIMAGES))
         ETOTAL = SUM(EEE(2:NIMAGES+1))
         WRITE(*,*) " congrad> E total: ", ETOTAL, "RMS: ", RMS, " E rep: ", SUM(EEER), " E constr: ", SUM(EEEC)
         WRITE(*,*) "                                              E spring: ", SUM(EEES), " E dih: ", SUM(EEED)
         WRITE(*,*) " congrad> FCONTEST: ", FCONTEST, " FREPTEST: ", FREPTEST
      END SUBROUTINE CONGRAD1


      SUBROUTINE CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, KINT, QCIFREEZET, QCIFROZEN, QCICONSTRREP, INTCONSTRAINTDEL, &
                            USEDIHEDRALCONST
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2)             ! energy for each image
         REAL(KIND = REAL64), INTENT(OUT) :: ETOTAL                    ! overall energy
         REAL(KIND = REAL64), INTENT(OUT) :: RMS                       ! total force
        
         REAL(KIND = REAL64) :: ECON, EREP, ESPR, EDIH      ! QUERY: should these be globals in keys?
         REAL(KIND = REAL64) :: EEEC(NIMAGES+2), GGGC(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEER(NIMAGES+2), GGGR(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEES(NIMAGES+2), GGGS(3*NATOMS*(NIMAGES+2))
         REAL(KIND = REAL64) :: EEED(NIMAGES+2), GGGD(3*NATOMS*(NIMAGES+2))
         INTEGER :: J1, J2


         ! initiate some variables
         EEE(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0
         ECON = 0.0D0; EREP = 0.0D0; ESPR = 0.0D0

         ! QUERY: what is INTCONSTRAINTDEL? seems like a scaling for the potential
         IF (.NOT.(INTCONSTRAINTDEL.EQ.0.0D0)) THEN
            !use the constraint energy including the internal extrema
            CALL GET_CONSTRAINT_E(XYZ,GGGC,EEEC,ECON)
         END IF
         IF (.NOT.(QCICONSTRREP.EQ.0.0D0)) THEN
            CALL GET_REPULSION_E2(XYZ,GGGR,EEER,EREP)
         END IF
         IF (.NOT.(KINT.EQ.0.0D0)) THEN
            CALL GET_SPRING_E(XYZ, GGGS, EEES, ESPR)
         END IF
         IF (USEDIHEDRALCONST) THEN
            CALL GET_DIH_CON_E(XYZ,GGGD,EEED,EDIH)
         END IF

         ! add all contributions
         EEE = EEEC + EEER + EEES + EEED
         GGG = GGGC + GGGR + GGGS + GGGD

         ! freeze atoms that should be frozen
         IF (QCIFREEZET) THEN
            DO J1=2,NIMAGES+1
               DO J2=1,NATOMS
                  IF (QCIFROZEN(J2)) THEN
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
                     GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
                  END IF
               END DO
            END DO
         END IF

         ! Set gradients to zero for start and finish images.
         GGG(1:(3*NATOMS))=0.0D0
         GGG((NIMAGES+1)*(3*NATOMS)+1:(NIMAGES+2)*(3*NATOMS))=0.0D0

         ! get the RMS force
         RMS=0.0D0
         DO J1=2,NIMAGES+1
            DO J2=1,3*NATOMS
               RMS = RMS + GGG((3*NATOMS)*(J1-1)+J2)**2
            END DO
         END DO
         RMS = SQRT(RMS/(3*NATOMS*NIMAGES))
         ETOTAL = SUM(EEE(2:NIMAGES+1))
      END SUBROUTINE CONGRAD2    

      SUBROUTINE GET_CONSTRAINT_E_NOINTERNAL(XYZ,GGG,EEE,ECON)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), ECON       ! energy for constraints
         INTEGER :: J1, J2, I
         INTEGER :: NI1, NJ1                         ! indices for atoms A and B in XYZ and GGG
         REAL(KIND=REAL64) :: CCLOCAL, CCLOCAL2      ! local concut
         REAL(KIND=REAL64) :: XA(3), XB(3)           ! coordinates for atoms A and B 
         REAL(KIND=REAL64) :: DIST, DUMMY, DUMMY2    ! distance and related measure
         REAL(KIND=REAL64) :: G2(3), GRADAB(3)
         REAL(KIND=REAL64) :: EMAX, FMIN, FMAX
         INTEGER :: IMAX, JMAX

         EMAX = -(HUGE(1.0D0))
         FMAX = -(HUGE(1.0D0))
         FMIN = HUGE(1.0D0)
         IMAX = -1
         JMAX = -1
         EEE(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0
         ECON = 0.0D0
         DO J2=1,NCONSTRAINT
            ! only active constraints contribute
            IF (.NOT.CONACTIVE(J2)) CYCLE
            ! get constraint cut off for this contraint
            CALL GET_CCLOCAL(J2,CCLOCAL)
            !!!!!Debugging WRITE(*,*) " Constraint: ", J2, " reference: ", CONDISTREFLOCAL(J2), " cut: ", CCLOCAL
            ! go through all images
            DO J1=2,NIMAGES+1
               NI1=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
               NJ1=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
               ! get coordinates for atoms involved in constraint
               DO I=1,3
                  XA(I) = XYZ(NI1+I)
                  XB(I) = XYZ(NJ1+I)
               END DO
               ! get distance and various related measures
               CALL DISTANCE_SIMPLE(XA, XB, DIST)
               DUMMY = DIST - CONDISTREFLOCAL(J2)
               ! now check whether we are beyond the constraint cutoff
               IF (DUMMY.GT.CCLOCAL) THEN
                  DUMMY2 = DUMMY**2
                  CCLOCAL2 = CCLOCAL**2  
                  ! calculate gradient and energy
                  G2 = (XA-XB)/DIST
                  GRADAB(1:3) = 2*(DUMMY2-CCLOCAL2)*DUMMY*G2(1:3)/(CCLOCAL2**2)
                  DUMMY = (DUMMY2-CCLOCAL2)**2/(2.0D0*CCLOCAL2**2)
                  EEE(J1) = EEE(J1) + DUMMY
                  ECON = ECON + DUMMY
                  !!!!!!debugging WRITE(*,*) " For image ", J1, " dist: ", DIST, " dummy: ", DIST - CONDISTREFLOCAL(J2), " energy: ", DUMMY 
                  ! save largest contribution to ECON
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX = J1
                     JMAX = J2
                     EMAX = DUMMY
                  END IF
                  GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+GRADAB(1:3)
                  GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-GRADAB(1:3)
                  DUMMY2=MINVAL(GRADAB)
                  IF (DUMMY2.LT.FMIN) FMIN=DUMMY2
                  DUMMY2=MAXVAL(GRADAB)
                  IF (DUMMY2.GT.FMAX) FMAX=DUMMY2
               END IF
            END DO
         END DO
         !QUERY: FIs FMAX definitely for the same JMAX and IMAX as EMAX?
         IF (-FMIN.GT.FMAX) FMAX=-FMIN
         FCONTEST=FMAX
         CONVERGECONTEST=EMAX
         MAXCONIMAGE = JMAX
         MAXCONSTR = IMAX

         IF (JMAX.GT.0) THEN
            WRITE(*,*) ' congrad> Highest constraint for image ',IMAX, ', con ',JMAX, ', atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX
         ENDIF
      END SUBROUTINE GET_CONSTRAINT_E_NOINTERNAL

      SUBROUTINE GET_DIH_CON_E(XYZ,GGG,EEE,EDIH)
         USE DIHEDRAL_CONSTRAINTS, ONLY: DIHEDRAL, DIHEDRALS, NDIH, REFDIH, ALLDIHACTIVE, DIHACTIVE, &
                                         CHECK_DIH_ACTIVE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), EDIH       ! energy for constraints  

         INTEGER :: A, B, C, D, N
         REAL(KIND = REAL64) :: XA(3), XB(3), XC(3), XD(4) !coordinates of atoms in dihedral
         REAL(KIND = REAL64) :: FA(3), FB(3), FC(3), FD(3) !returned gradient for individual atoms
         REAL(KIND = REAL64) :: PHIREF, THISE

         EDIH = 0.0D0
         GGG = 0.0D0
         EEE = 0.0D0

         ! if not all dihedral constraints are activated, update the list
         IF (.NOT.ALLDIHACTIVE) CALL CHECK_DIH_ACTIVE()

         DO J=1,NDIH
            !cycle if the dihedral constraint is inactive
            IF (.NOT.DIHACTIVE(J)) CYCLE
            !look up atoms in dihedral
            A = DIHEDRALS(J,1)
            B = DIHEDRALS(J,2)
            C = DIHEDRALS(J,3)
            D = DIHEDRALS(J,4)
            PHIREF = REF(DIH(J))
            DO I=2,NIMAGES+1
               !reference for image we are in (The x ccoord of the first atom of the current image is N+1)
               N = 3*NATOMS*(I-1)
               !extract relevant coordinates
               XA(1:3) = XYZ(N+3*A-2:N+3*A)
               XB(1:3) = XYZ(N+3*B-2:N+3*B)
               XC(1:3) = XYZ(N+3*C-2:N+3*C)
               XD(1:3) = XYZ(N+3*D-2:N+3*D)
               !call routine to compute dihedral and get gradient
               CALL DIHEDRAL(XA, XB, XC, XD, PHIREF, THISE, FA, FB, FC, FD)
               !add results to appropriate variables
               EDIH = EDIH + THISE
               EEE(I) = EEE(I) + THISE
               GGG(N+3*A-2:N+3*A) = GGG(N+3*A-2:N+3*A) - FA(1:3)
               GGG(N+3*B-2:N+3*B) = GGG(N+3*B-2:N+3*B) - FB(1:3)
               GGG(N+3*C-2:N+3*C) = GGG(N+3*C-2:N+3*C) - FC(1:3)
               GGG(N+3*D-2:N+3*D) = GGG(N+3*D-2:N+3*D) - FD(1:3)
            END DO
         END DO
      END SUBROUTINE GET_DIH_CON_E

      SUBROUTINE GET_CONSTRAINT_E(XYZ,GGG,EEE,ECON)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, INTMINFAC, CHECKCONINT, INTCONSTRAINTDEL, &
                            CONACTINACT, USECONACTINACT
         USE QCI_CONSTRAINT_KEYS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREFLOCAL
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE, ATOMACTIVE
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE, DOTP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), ECON       ! energy for constraints         
         INTEGER :: J1, J2       
         INTEGER :: NI1, NI2, NJ1, NJ2                                 ! indices for atom I and J in images 1 and 2
         REAL(KIND = REAL64) :: G1(3), G2(3), DSQ1, DSQ2, DINTMIN, DP_G12 ! dummy variables
         REAL(KIND = REAL64) :: CONSTGRAD(3), DUMMY
         REAL(KIND = REAL64) :: CCLOCAL
         LOGICAL :: NOINT                                            ! do we hae an internal minimum
         REAL(KIND = REAL64) :: DINT, DSQI, G1INT(3), G2INT(3)       ! information for internal minimum contribution
         REAL(KIND = REAL64) , PARAMETER :: DINTTEST = 1.0D-50       ! cutoff for DINTMIN to test for internal min
         REAL(KIND = REAL64) :: D1, D2, D12                         ! distances in image 1 and 2         
         REAL(KIND=REAL64) :: EMAX, FMIN, FMAX
         REAL(KIND=REAL64) :: LOCALCONFACTOR
         INTEGER :: IMAX, JMAX

         EEE(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0
         ECON = 0.0D0
         EMAX = -(HUGE(1.0D0))
         FMAX = -(HUGE(1.0D0))
         FMIN = HUGE(1.0D0)
         IMAX = -1
         JMAX = -1

         !  Constraint energy and forces.
         !
         ! For J1 we consider the line segment between image J1-1 and J1.
         ! There are NIMAGES+1 line segments in total, with an energy contribution
         ! and corresponding gradient terms for each. 
         ! A and B refer to atoms, 1 and 2 to images J1-1 and J1 corresponding to J1-2 and J1-1 below.

         DO J2=1,NCONSTRAINT
            IF (.NOT.USECONACTINACT) THEN
               ! only active constraints contribute
               IF (.NOT.CONACTIVE(J2)) CYCLE
            END IF

            IF (ATOMACTIVE(CONI(J2)).AND.ATOMACTIVE(CONJ(J2))) THEN
               !both atoms active:
               LOCALCONFACTOR = 1.0D0
            ELSE IF ((ATOMACTIVE(CONI(J2)).AND.(.NOT.ATOMACTIVE(CONJ(J2)))).OR.((.NOT.ATOMACTIVE(CONI(J2))).AND.ATOMACTIVE(CONJ(J2)))) THEN
               ! one atom active
               LOCALCONFACTOR = CONACTINACT
            ELSE
               CYCLE
            END IF
            ! get constraint cut off for this contraint
            CALL GET_CCLOCAL(J2,CCLOCAL)            
            ! go through all images 
            ! QUERY: here we differ from the no internal minimum routine - why are we going up to NIMAGES+2
            DO J1=2,NIMAGES+2
               NI1=(3*NATOMS)*(J1-2)+3*(CONI(J2)-1)
               NI2=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
               NJ1=(3*NATOMS)*(J1-2)+3*(CONJ(J2)-1)
               NJ2=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
         
               G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3) !vector from j to i in image 1
               G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3) !vector from j to i in image 2

               ! squared distance between atoms in image 1 (theta = pi/2)
               DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
               ! squared distance between atoms in image 2 (theta = 0)
               DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
               DP_G12 = DOTP(3,G1,G2)
               DINTMIN = DSQ1+DSQ2-2.0D0*DP_G12

               ! Convert derivatives of distance^2 to derivative of distance.
               ! We have cancelled a factor of two above and below
               D1 = SQRT(DSQ1); D2 = SQRT(DSQ2)
               G1(1:3) = G1(1:3)/D1; G2(1:3) = G2(1:3)/D2

               IF ((.NOT.CHECKCONINT).OR.(DINTMIN.LT.DINTTEST)) THEN
                  NOINT = .TRUE.
               ELSE
                  CALL INTMIN_CONSTRAINT(G1,G2,DSQ1,DSQ2,DP_G12,DINTMIN,NOINT,DSQI,DINT,G1INT,G2INT)
               END IF
               !TODO: CHECK FORMUATION - ARE THESE THE SAME?!?!?
               ! Need to include both D2 and D1 contributions if they are both outside tolerance.
               ! Otherwise we get discontinuities if they are very close and swap over.

               ! terms for image J1 - non-zero derivatives only for J1. D2 is the distance for image J1.
               ! these are the constraint energies as in get_constraint_e_nointernal
               DUMMY = D2-CONDISTREFLOCAL(J2)
               IF ((DUMMY.GT.CCLOCAL).AND.(J1.LT.NIMAGES+2)) THEN  
                  !CONSTGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
                  !DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
                  CONSTGRAD(1:3)=LOCALCONFACTOR*2*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
                  DUMMY=LOCALCONFACTOR*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)                  
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX=J1
                     JMAX=J2
                     EMAX=DUMMY
                  ENDIF
                  EEE(J1)=EEE(J1)+DUMMY
                  ECON=ECON+DUMMY
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+CONSTGRAD(1:3)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-CONSTGRAD(1:3)
               END IF
               ! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
               ! terms for image J1-1 - non-zero derivatives only for J1-1. D1 is the distance for image J1-1.
               ! here we have the internal extremum contribution
               !QUERY: Does this need to be absolute??? Or should it be the relative value? Not sure it makes sense to have this as the absolute value?
               IF (CHECKCONINT.AND.(.NOT.NOINT).AND.(ABS(DINT-CONDISTREFLOCAL(J2)).GT.CCLOCAL)) THEN
                  DUMMY=DINT-CONDISTREFLOCAL(J2)  
                  !CONSTGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1INT(1:3)
                  CONSTGRAD(1:3)=2*INTMINFAC*LOCALCONFACTOR*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1INT(1:3)
                  GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+CONSTGRAD(1:3)
                  GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-CONSTGRAD(1:3)
                  !CONSTGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2INT(1:3)
                  !DUMMY=INTMINFAC*INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
                  CONSTGRAD(1:3)=2*INTMINFAC*LOCALCONFACTOR*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2INT(1:3)
                  DUMMY=INTMINFAC*LOCALCONFACTOR*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
                  ECON=ECON+DUMMY
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX=J1
                     JMAX=J2
                     EMAX=DUMMY
                  ENDIF
                  IF (J1.EQ.2) THEN
                     EEE(J1)=EEE(J1)+DUMMY
                  ELSE IF (J1.LT.NIMAGES+2) THEN
                     EEE(J1)=EEE(J1)+DUMMY/2.0D0
                     EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
                  ELSE IF (J1.EQ.NIMAGES+2) THEN
                     EEE(J1-1)=EEE(J1-1)+DUMMY
                  ENDIF
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+CONSTGRAD(1:3)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-CONSTGRAD(1:3)
               ENDIF
            END DO
         END DO
         CONVERGECONTEST=EMAX/INTCONSTRAINTDEL
         IF (-FMIN.GT.FMAX) FMAX=-FMIN
         FCONTEST=FMAX
         MAXCONIMAGE = JMAX
         MAXCONSTR = IMAX

         IF (JMAX.GT.0) THEN
            WRITE(*,*) ' congrad> Highest constraint for image ',IMAX, ', con ',JMAX, ', atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX
         ENDIF
      END SUBROUTINE GET_CONSTRAINT_E

      SUBROUTINE GET_REPULSION_E(XYZ,GGG,EEE,EREP)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, INTMINFAC, QCIINTREPMINSEP
         USE REPULSION, ONLY: NNREPULSIVE, NREPI, NREPJ, NREPCUT
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE, DOTP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), EREP       ! energy for repulsions
         
         INTEGER :: J1, J2, OFFSET1, OFFSET2
         REAL(KIND = REAL64) :: X1(3*NATOMS), X2(3*NATOMS)           ! coordinates for image 1 and 2
         REAL(KIND = REAL64) :: GLOCAL1(3*NATOMS), GLOCAL2(3*NATOMS) ! gradients for each image
         REAL(KIND = REAL64) :: EREP1, EREP2                         ! repulsion energy contributions
         INTEGER :: NI, NJ                                           ! indices for atom I and J
         REAL(KIND = REAL64) :: G1(3), G2(3), DSQ1, DSQ2             ! dummy variables
         REAL(KIND = REAL64) :: DCUT, DINTMIN, DP_G12                ! values used to test for internal minimum
         REAL(KIND = REAL64) , PARAMETER :: DINTTEST = 1.0D-50       ! cutoff for DINTMIN to test for internal min
         REAL(KIND = REAL64) :: D1, D2, D12                         ! distances in image 1 and 2
         REAL(KIND = REAL64) :: RPLOCAL, RPLOCAL2, RPLOCALINV 
         REAL(KIND = REAL64) :: REPGRAD(3)
         LOGICAL :: NOINT                                            ! do we hae an internal minimum
         REAL(KIND = REAL64) :: DINT, DSQI, G1INT(3), G2INT(3)       ! information for internal minimum contribution
         REAL(KIND=REAL64) :: EMAX, FMIN, FMAX
         REAL(KIND=REAL64) :: DUMMY, DUMMY2
         INTEGER :: IMAX, JMAX

         EMAX = -(HUGE(1.0D0))
         FMAX = -(HUGE(1.0D0))
         FMIN = HUGE(1.0D0)
         EEE(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0         
         EREP = 0.0D0
         IMAX = -1
         JMAX = -1

         DO J1=2,NIMAGES+1
            ! get coordinates for images
            OFFSET2=(3*NATOMS)*(J1-1)
            OFFSET1=(3*NATOMS)*(J1-2)
            X1(1:3*NATOMS)=XYZ(OFFSET1+1:OFFSET1+3*NATOMS)
            X2(1:3*NATOMS)=XYZ(OFFSET2+1:OFFSET2+3*NATOMS) 
            ! initialise some local vairables
            GLOCAL1(1:3*NATOMS) = 0.0D0; GLOCAL2(1:3*NATOMS) = 0.0D0
            EREP1 = 0.0D0; EREP2 = 0.0D0
            ! iterate over all repulsions
            DO J2 = 1,NNREPULSIVE
               NI=3*(NREPI(J2)-1)
               NJ=3*(NREPJ(J2)-1)             
               G1(1:3)=X1(NI+1:NI+3)-X1(NJ+1:NJ+3) !vector from j to i in image 1
               G2(1:3)=X2(NI+1:NI+3)-X2(NJ+1:NJ+3) !vector from j to i in image 2
               ! squared distance between atoms in image 1 (theta = pi/2)
               DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
               ! squared distance between atoms in image 2 (theta = 0)
               DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
               DCUT = NREPCUT(J2)**2
               ! don't look for an internal minimum if both repulsions outside cutoff
               IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE
               ! don't check for internal minimum in distance - atoms too close for chain crossing.
               IF (ABS(NREPI(J2)-NREPJ(J2)).LT.QCIINTREPMINSEP) THEN
                  DINTMIN = 0.0D0
               ELSE 
                  DP_G12 = DOTP(3,G1,G2)
                  DINTMIN = DSQ1+DSQ2-2.0D0*DP_G12
               END IF

               ! Convert derivatives of distance^2 to derivative of distance.
               ! We have cancelled a factor of two above and below
               D1 = SQRT(DSQ1); D2 = SQRT(DSQ2)
               G1(1:3) = G1(1:3)/D1; G2(1:3) = G2(1:3)/D2

               ! if the denominator in the d^2 is approx zero, we do not need to check for an internal minimum
               IF (DINTMIN.LT.DINTTEST) THEN
                  NOINT = .TRUE.
               ELSE
                  CALL INTMIN_REPULSION(G1,G2,DSQ1,DSQ2,DP_G12,DINTMIN,NOINT,DSQI,DINT,G1INT,G2INT)
               END IF

               RPLOCAL = NREPCUT(J2)
               RPLOCAL2 = RPLOCAL**2
               RPLOCALINV = 1.0D0/RPLOCAL

               IF (D2.LT.RPLOCAL) THEN
                  DUMMY = RPLOCAL2/DSQ2 + 2.0D0*D2*RPLOCALINV - 3.0D0
                  !WRITE(*,*) " EREP1: ", DUMMY, DSQ2, D2, RPLOCAL2, RPLOCALINV
                  EREP1 = EREP1 + DUMMY
                  EREP = EREP + DUMMY
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX = J1
                     JMAX = J2
                     EMAX = DUMMY
                  END IF
                  DUMMY=-2.0D0*(RPLOCAL2/(D2*DSQ2)-RPLOCALINV)
                  REPGRAD(1:3) = DUMMY*G2(1:3)
                  GLOCAL2(NI+1:NI+3)=GLOCAL2(NI+1:NI+3)+REPGRAD(1:3)
                  GLOCAL2(NJ+1:NJ+3)=GLOCAL2(NJ+1:NJ+3)-REPGRAD(1:3)
               END IF  
               ! For internal minima we are counting edges. 
               ! Edge J1 is between images J1-1 and J1, starting from J1=2.
               ! Energy contributions are shared evenly, except for
               ! edge 1, which was assigned to image 2, and edge NIMAGES+1, which
               ! was assigned to image NIMAGES+1. 
               DUMMY = 0.0D0
               IF ((.NOT.NOINT).AND.(DINT.LT.RPLOCAL).AND.(J1.NE.2)) THEN
                  D12 = DSQI !from call to find internal minimum
                  DUMMY=INTMINFAC*(RPLOCAL2/D12+2.0D0*DINT*RPLOCALINV-3.0D0)
                  !WRITE(*,*) " EREP2: ", DUMMY, DSQI, DINT, RPLOCAL2, RPLOCALINV
                  EREP2=EREP2+DUMMY
                  EREP=EREP+DUMMY
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX=J1
                     JMAX=J2
                     EMAX=DUMMY
                  ENDIF
                  DUMMY=-2.0D0*(RPLOCAL2/(DINT*D12)-RPLOCALINV)
                  ! Gradient contributions for image J1-1
                  REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
                  GLOCAL1(NI+1:NI+3)=GLOCAL1(NI+1:NI+3)+REPGRAD(1:3)
                  GLOCAL1(NJ+1:NJ+3)=GLOCAL1(NJ+1:NJ+3)-REPGRAD(1:3)
                  DUMMY2=MINVAL(REPGRAD)
                  ! Gradient contributions for image J1
                  REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
                  GLOCAL2(NI+1:NI+3)=GLOCAL2(NI+1:NI+3)+REPGRAD(1:3)
                  GLOCAL2(NJ+1:NJ+3)=GLOCAL2(NJ+1:NJ+3)-REPGRAD(1:3)
               END IF
               !WRITE(*,*) " EREP1: ", EREP1, " EREP2: ", EREP2
            END DO
            GGG(OFFSET1+1:OFFSET1+3*NATOMS)=GGG(OFFSET1+1:OFFSET1+3*NATOMS)+GLOCAL1(1:3*NATOMS)
            GGG(OFFSET2+1:OFFSET2+3*NATOMS)=GGG(OFFSET2+1:OFFSET2+3*NATOMS)+GLOCAL2(1:3*NATOMS)
            EEE(J1)=EEE(J1)+EREP1
            IF (J1.EQ.2) THEN
               EEE(J1)=EEE(J1)+EREP2
            ELSE
               EEE(J1)=EEE(J1)+EREP2/2.0D0
               EEE(J1-1)=EEE(J1-1)+EREP2/2.0D0
            ENDIF
         END DO
         FMIN=MINVAL(GGG(3*NATOMS+1:3*NATOMS*(NIMAGES+1)))
         FMAX=MAXVAL(GGG(3*NATOMS+1:3*NATOMS*(NIMAGES+1)))
         IF (-FMIN.GT.FMAX) FMAX=-FMIN
         FREPTEST=FMAX
         CONVERGEREPTEST=EMAX
         MAXCONIMAGE = JMAX
         MAXCONSTR = IMAX
      END SUBROUTINE GET_REPULSION_E

      ! should be the same as GET_REPULSION_E, but the iteration inverts the order - outer loops is repulsions
      SUBROUTINE GET_REPULSION_E2(XYZ,GGG,EEE,EREP)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, INTMINFAC, QCICONSTRREP
         USE REPULSION, ONLY: NNREPULSIVE, NREPI, NREPJ, NREPCUT
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE, DOTP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), EREP       ! energy for repulsions
         
         INTEGER :: J1, J2, OFFSET1, OFFSET2
         REAL(KIND = REAL64) :: X1(3*NATOMS), X2(3*NATOMS)           ! coordinates for image 1 and 2
         REAL(KIND = REAL64) :: GLOCAL1(3*NATOMS), GLOCAL2(3*NATOMS) ! gradients for each image
         REAL(KIND = REAL64) :: EREP1, EREP2                         ! repulsion energy contributions
         INTEGER :: NI1, NI2, NJ1, NJ2                               ! indices for atom I and J in images 1 and 2
         REAL(KIND = REAL64) :: G1(3), G2(3), DSQ1, DSQ2             ! dummy variables
         REAL(KIND = REAL64) :: DCUT, DINTMIN, DP_G12                ! values used to test for internal minimum
         REAL(KIND = REAL64) , PARAMETER :: DINTTEST = 1.0D-50       ! cutoff for DINTMIN to test for internal min
         REAL(KIND = REAL64) :: D1, D2, D12                         ! distances in image 1 and 2
         REAL(KIND = REAL64) :: RPLOCAL, INTCONST, INTCONSTINV
         REAL(KIND = REAL64) :: REPGRAD(3)
         LOGICAL :: NOINT                                            ! do we hae an internal minimum
         REAL(KIND = REAL64) :: DINT, DSQI, G1INT(3), G2INT(3)       ! information for internal minimum contribution
         REAL(KIND=REAL64) :: EMAX, FMIN, FMAX
         REAL(KIND=REAL64) :: DUMMY
         INTEGER :: IMAX, JMAX

         WRITE(*,*) " E repulsion 2"
         EMAX = -(HUGE(1.0D0))
         FMAX = -(HUGE(1.0D0))
         FMIN = HUGE(1.0D0)
         EEE(1:NIMAGES+2)=0.0D0
         GGG(1:(3*NATOMS)*(NIMAGES+2))=0.0D0         
         EREP = 0.0D0
         IMAX = -1
         JMAX = -1

         DO J2=1,NNREPULSIVE
            RPLOCAL = NREPCUT(J2)
            INTCONST = RPLOCAL**3
            INTCONSTINV = 1.0D0/INTCONST

            DO J1=2,NIMAGES+2
               NI1=(3*NATOMS)*(J1-2)+3*(NREPI(J2)-1)
               NI2=(3*NATOMS)*(J1-1)+3*(NREPI(J2)-1)
               NJ1=(3*NATOMS)*(J1-2)+3*(NREPJ(J2)-1)
               NJ2=(3*NATOMS)*(J1-1)+3*(NREPJ(J2)-1)
         
               G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3) !vector from j to i in image 1
               G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3) !vector from j to i in image 2

               ! squared distance between atoms in image 1 (theta = pi/2)
               DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
               ! squared distance between atoms in image 2 (theta = 0)
               DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
               DCUT=NREPCUT(J2)**2
               ! don't look for an internal minimum if both repulsions outside cutoff
               IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE 
               !QUERY: in the other repulsion routine we use an additional cutoff with QCIINTREPMINSEP - why not here?
               DP_G12 = DOTP(3,G1,G2)
               DINTMIN = DSQ1+DSQ2-2.0D0*DP_G12

               ! Convert derivatives of distance^2 to derivative of distance.
               ! We have cancelled a factor of two above and below
               D1 = SQRT(DSQ1); D2 = SQRT(DSQ2)
               G1(1:3) = G1(1:3)/D1; G2(1:3) = G2(1:3)/D2

               ! if the denominator in the d^2 is approx zero, we do not need to check for an internal minimum
               IF (DINTMIN.LT.DINTTEST) THEN
                  NOINT = .TRUE.
               ELSE
                  CALL INTMIN_REPULSION(G1,G2,DSQ1,DSQ2,DP_G12,DINTMIN,NOINT,DSQI,DINT,G1INT,G2INT)
               END IF
               ! terms for image J1 - non-zero derivatives only for J1
               IF ((D2.LT.RPLOCAL).AND.(J1.LT.NIMAGES+2)) THEN
                  DUMMY=QCICONSTRREP*(1.0D0/DSQ2+(2.0D0*D2-3.0D0*RPLOCAL)*INTCONSTINV)
                  EEE(J1)=EEE(J1)+DUMMY
                  EREP=EREP+DUMMY
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX=J1
                     JMAX=J2
                     EMAX=DUMMY
                  ENDIF
                  DUMMY=-2.0D0*QCICONSTRREP*(1.0D0/(D2*DSQ2)-INTCONSTINV)
                  REPGRAD(1:3)=DUMMY*G2(1:3)
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
               END IF
               DUMMY=0.0D0
               IF ((.NOT.NOINT).AND.(DINT.LT.RPLOCAL)) THEN
                  DUMMY=INTMINFAC*QCICONSTRREP*(1.0D0/DSQI+(2.0D0*DINT-3.0D0*RPLOCAL)*INTCONSTINV)
                  EREP=EREP+DUMMY
                  IF (DUMMY.GT.EMAX) THEN
                     IMAX=J1
                     JMAX=J2
                     EMAX=DUMMY
                  ENDIF
                  IF (J1.EQ.2) THEN
                     EEE(J1)=EEE(J1)+DUMMY         
                  ELSE IF (J1.LT.NIMAGES+2) THEN
                     EEE(J1)=EEE(J1)+DUMMY/2.0D0
                     EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
                  ELSE IF (J1.EQ.NIMAGES+2) THEN
                     EEE(J1-1)=EEE(J1-1)+DUMMY
                  ENDIF
                  DUMMY=-2.0D0*QCICONSTRREP*(1.0D0/(DINT*DSQI)-INTCONSTINV)
                  REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
                  GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
                  GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
                  REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
               ENDIF
            END DO
         END DO
         FMIN=MINVAL(GGG(3*NATOMS+1:3*NATOMS*(NIMAGES+1)))
         FMAX=MAXVAL(GGG(3*NATOMS+1:3*NATOMS*(NIMAGES+1)))
         IF (-FMIN.GT.FMAX) FMAX=-FMIN
         FREPTEST=FMAX
         CONVERGEREPTEST=EMAX/QCICONSTRREP
         MAXCONIMAGE = JMAX
         MAXCONSTR = IMAX
      END SUBROUTINE GET_REPULSION_E2

      SUBROUTINE GET_SPRING_E(XYZ, GGG, EEE, ESPR)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS, KINT, KINTSCALED, QCIADJUSTKT, QCISPRINGACTIVET
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGES+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGES+2), ESPR      ! energy for repulsions
         
         INTEGER :: J1, J2, NI1, NI2
         REAL(KIND = REAL64) :: DPLUS, DUMMY, EMAX, SPGRAD(3)
         REAL(KIND = REAL64) :: DVEC(NIMAGES+1)
         INTEGER :: IMAX

         ESPR = 0.0D0
         EMAX = -(HUGE(1.0D0))
         IMAX = -1

         DO J1=1,NIMAGES+1
            NI1 = (3*NATOMS)*(J1-1)
            NI2 = (3*NATOMS)*J1
            DPLUS = 0.0D0
            ! get energy for springs
            DO J2=1,NATOMS
               IF ((.NOT.QCISPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
                  DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
        &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
        &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
               ENDIF               
            END DO
            DVEC(J1) = SQRT(DPLUS)
            DUMMY = KINT*0.5D0*DPLUS/KINTSCALED
            IF (DUMMY.GT.EMAX) THEN
               IMAX=J1
               EMAX=DUMMY
            ENDIF
            ESPR = ESPR + DUMMY
            ! get gradient
            DUMMY=KINT/KINTSCALED
            DO J2=1,NATOMS
               SPGRAD = 0.0D0
               IF ((.NOT.QCISPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
                  SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
                  GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
                  GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
               ENDIF
            ENDDO
         END DO
         MAXSPRIMAGE = IMAX
         EMAXSPR = EMAX
         IF (QCIADJUSTKT) CALL GET_AV_DEV(DVEC)
      END SUBROUTINE GET_SPRING_E

      SUBROUTINE GET_AV_DEV(DVEC)
         USE QCIKEYS, ONLY: NIMAGES, QCIAVDEV
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: DVEC(1:NIMAGES+1)
         REAL(KIND = REAL64) :: SEPARATION, DEVIATION(1:NIMAGES+1)

         SEPARATION = SUM(DVEC(1:NIMAGES+1))
         DEVIATION(1:NIMAGES+1)=ABS(100*((NIMAGES+1)*DVEC(1:NIMAGES+1)/SEPARATION-1.0D0))
         QCIAVDEV=SUM(DEVIATION)/(NIMAGES+1)
      END SUBROUTINE GET_AV_DEV

      SUBROUTINE GET_CCLOCAL(CONID,CCLOCAL)
         USE QCI_CONSTRAINT_KEYS, ONLY: CONCUTLOCAL, CONCUTABST, CONCUTABS, &
                                        CONCUTFRACT, CONCUTFRAC, CONDISTREFLOCAL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CONID
         REAL(KIND=REAL64), INTENT(OUT) :: CCLOCAL
         !QUERY: this seems odd - why are they applied consecutively?
         CCLOCAL=CONCUTLOCAL(CONID)
         IF (CONCUTABST) THEN
            CCLOCAL=CCLOCAL+CONCUTABS
         ELSE IF (CONCUTFRACT) THEN
            CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(CONID)
         END IF
      END SUBROUTINE GET_CCLOCAL

      SUBROUTINE INTMIN_REPULSION(G1,G2,DSQ1,DSQ2,DP_G12,DINTMIN,NOINT,DSQI,DINT,G1INT,G2INT)
         USE QCIKEYS, ONLY: QCIREPCUT
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: G1(3), G2(3)
         REAL(KIND = REAL64), INTENT(IN) :: DSQ1, DSQ2
         REAL(KIND = REAL64), INTENT(IN) :: DP_G12, DINTMIN
         LOGICAL, INTENT(OUT) :: NOINT
         REAL(KIND = REAL64), INTENT(OUT) :: DSQI, DINT
         REAL(KIND = REAL64), INTENT(OUT) :: G1INT(3), G2INT(3)
         REAL(KIND = REAL64) :: DUMMY, DUMMY2, DP_G12_SQ
         REAL(KIND = REAL64), PARAMETER :: LARGEDIST = 1.0D10

         !initliase all output variables
         NOINT = .TRUE.
         DSQI = LARGEDIST; DINT = LARGEDIST
         G1INT(1:3)=0.0D0; G2INT(1:3)=0.0D0

         ! are we having an internal minimum?
         DUMMY = (DSQ1-DP_G12)/DINTMIN
         IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) THEN
            NOINT=.FALSE.
            DP_G12_SQ = DP_G12**2
            DUMMY2 = DP_G12_SQ - DSQ1*DSQ2
            DSQI = MAX(DUMMY2/DINTMIN,0.0D0)
            DINT = SQRT(DSQI)
            IF (DINT.LE.0.0D0) THEN
               NOINT = .TRUE.
            ELSE IF (DINT.LE.QCIREPCUT) THEN
               DUMMY = DINT*DINTMIN**2
               ! to convert derivatives of distance^2 to derivative of distance.
               G1INT(1:3)= (DUMMY2*(G1(1:3) - G2(1:3)) + DINTMIN*(G1(1:3)*DSQ2 -G2(1:3)*DP_G12))/DUMMY
               G2INT(1:3)= (DUMMY2*(G2(1:3) - G1(1:3)) + DINTMIN*(G2(1:3)*DSQ1 -G1(1:3)*DP_G12))/DUMMY
            END IF              
         END IF
      END SUBROUTINE INTMIN_REPULSION

      ! QUERY: these functions are identical apart from use of the ondition for repulsions - is that condition required?
      SUBROUTINE INTMIN_CONSTRAINT(G1,G2,DSQ1,DSQ2,DP_G12,DINTMIN,NOINT,DSQI,DINT,G1INT,G2INT)
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: G1(3), G2(3)
         REAL(KIND = REAL64), INTENT(IN) :: DSQ1, DSQ2
         REAL(KIND = REAL64), INTENT(IN) :: DP_G12, DINTMIN
         LOGICAL, INTENT(OUT) :: NOINT
         REAL(KIND = REAL64), INTENT(OUT) :: DSQI, DINT
         REAL(KIND = REAL64), INTENT(OUT) :: G1INT(3), G2INT(3)
         REAL(KIND = REAL64) :: DUMMY, DUMMY2, DP_G12_SQ
         REAL(KIND = REAL64), PARAMETER :: LARGEDIST = 1.0D10

         !initliase all output variables
         NOINT = .TRUE.
         DSQI = LARGEDIST; DINT = LARGEDIST
         G1INT(1:3)=0.0D0; G2INT(1:3)=0.0D0

         ! are we having an internal minimum?
         DUMMY = (DSQ1-DP_G12)/DINTMIN
         IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) THEN
            NOINT=.FALSE.
            DP_G12_SQ = DP_G12**2
            DUMMY2 = DP_G12_SQ - DSQ1*DSQ2
            DSQI = MAX(DUMMY2/DINTMIN,0.0D0)
            DINT = SQRT(DSQI)
            IF (DINT.LE.0.0D0) THEN
               NOINT = .TRUE.
            ELSE
               DUMMY = DINT*DINTMIN**2
               ! to convert derivatives of distance^2 to derivative of distance.
               G1INT(1:3)= (DUMMY2*(G1(1:3) - G2(1:3)) + DINTMIN*(G1(1:3)*DSQ2 -G2(1:3)*DP_G12))/DUMMY
               G2INT(1:3)= (DUMMY2*(G2(1:3) - G1(1:3)) + DINTMIN*(G2(1:3)*DSQ1 -G1(1:3)*DP_G12))/DUMMY
            END IF              
         END IF
      END SUBROUTINE INTMIN_CONSTRAINT

END MODULE CONSTR_E_GRAD
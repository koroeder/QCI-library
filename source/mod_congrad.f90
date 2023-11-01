! Module containing the congrad routines
!
! version 1 (congrad1) - tests for internal minimum in repulsions only
! version 2 (congrad2) - tests for internal minimum in repulsion and constraints
MODULE CONSTR_E_GRAD
   USE QCIPREC
   IMPLICIT NONE
   INTEGER :: MAXCONIMAGE, MAXREPIMAGE                   ! image with highest constraint / repulsion
   INTEGER :: MAXCONSTR, MAXREP                          ! index for worst constraint / repulsion
   REAL(KIND=REAL64) :: CONVERGECONTEST, CONVERGEREPTEST ! energy of that term
   REAL(KIND=REAL64) :: FCONMAX, FREPMAX                 ! maximum gradient
   CONTAINS

      SUBROUTINE CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
         USE QCIKEYS, ONLY: NIMAGE, NATOMS
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGE+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGE+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGE+2)             ! energy for each image
         REAL(KIND = REAL64), INTENT(OUT) :: ETOTAL                    ! overall energy
         REAL(KIND = REAL64), INTENT(OUT) :: RMS                       ! total force
        
         REAL(KIND = REAL64) :: ECON, EREP      ! QUERY: should these be globals in keys?
         REAL(KIND = REAL64) :: EEEC(NIMAGE+2), GGGC(3*NATOMS*(NIMAGE+2))
         REAL(KIND = REAL64) :: EEER(NIMAGE+2), GGGR(3*NATOMS*(NIMAGE+2))

         ! initiate some variables
         EEE(1:INTIMAGE+2)=0.0D0
         GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
         ECON = 0.0D0; EREP = 0.0D0

         ! QUERY: what is INTCONSTRAINTDEL? seems like a scaling for the potential
         IF (.NOT.(INTCONSTRAINTDEL.EQ.0.0D0)) THEN
            CALL GET_CONSTRAINT_E_NOINTERNAL(XYZ,GGGC,EEEC,ECON)

         END IF

      END SUBROUTINE CONGRAD1


      SUBROUTINE CONGRAD2()

      END SUBROUTINE CONGRAD2    

      SUBROUTINE GET_CONSTRAINT_E_NOINTERNAL(XYZ,GGG,EEE,ECON)
         USE QCIKEYS, ONLY: NIMAGE, NATOMS
         USE QCICONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ
         USE INTERPOLATION_KEYS, ONLY: CONACTIVE
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGE+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGE+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGE+2), ECON       ! energy for constraints
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
         JMAX = -1
         EEE(1:INTIMAGE+2)=0.0D0
         GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
         ECON = 0.0D0

         DO J2=1,NCONSTRAINT
            ! only active constraints contribute
            IF (.NOT.CONACTIVE(J2)) CYCLE
            ! get constraint cut off for this contraint
            CALL GET_CCLOCAL(CCLOCAL)
            ! go through all images
            DO J1=2,INTIMAGE+1
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
               DUMMY2 = DUMMY**2
               CCLOCAL2 = CCLOCAL**2
               ! now check whether we are beyond the constraint cutoff
               IF (DUMMY2.GT.CCLOCAL2) THEN
                  ! calculate gradient and energy
                  G2 = (XA-XB)/D2
                  GRADAB(1:3) = 2*(DUMMY2-CCLOCAL2)*DUMMY*G2(1:3)/(CCLOCAL2**2)
                  DUMMY = (DUMMY2-CCLOCAL2)**2/(2.0D0*CCLOCAL2**2)
                  EEE(J1) = EEE(J1) + DUMMY
                  ECON = ECON + DUMMY
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

         ! IF (JMAX.GT.0) THEN
            ! WRITE(*,*) ' congrad> Highest constraint for image ',IMAX, ', con ',JMAX, ', atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX
         ! ENDIF
      END SUBROUTINE GET_CONSTRAINT_E_NOINTERNAL

      SUBROUTINE GET_REPULSION_E(XYZ,GGG,EEE,EREP)
         USE QCIKEYS, ONLY: NIMAGE, NATOMS
         USE REPULSION, ONLY: NNREPULSIVE, NREPI, NREPJ, NREPCUT
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGE+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGE+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGE+2), EREP       ! energy for repulsions
         
         INTEGER :: J1, J2, OFFSET1, OFFSET2
         REAL(KIND = REAL64) :: X1(3*NATOMS), X2(3*NATOMS)           ! coordinates for image 1 and 2
         REAL(KIND = REAL64) :: GLOCAL1(3*NATOMS), GLOCAL2(3*NATOMS) ! gradients for each image
         REAL(KIND = REAL64) :: EREP1, EREP2                         ! repulsion energy contributions
         INTEGER :: NI1, NI2, NJ1, NJ2                               ! indices for atom I and J in images 1 and 2
         REAL(KIND = REAL64) :: G1(3), G2(3), DSQ1, DSQ2             ! dummy variables
         REAL(KIND = REAL64) :: DCUT, DINTMIN, DP_G12                ! values used to test for internal minimum
         REAL(KIND = REAL64) , PARAMETER :: DINTTEST = 1.0D-50       ! cutoff for DINTMIN to test for internal min
         REAL(KIND = REAL64) :: D1, D2                               ! distances in image 1 and 2
         LOGICAL :: NOINT                                            ! no internal minimum

         DO J1=2,NIMAGE+1
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
               NI1=3*(NREPI(J2)-1)
               NI2=3*(NREPI(J2)-1)
               NJ1=3*(NREPJ(J2)-1)
               NJ2=3*(NREPJ(J2)-1)               
               G1(1:3)=X1(NI1+1:NI1+3)-X1(NJ1+1:NJ1+3) !vector from j to i in image 1
               G2(1:3)=X2(NI2+1:NI2+3)-X2(NJ2+1:NJ2+3) !vector from j to i in image 2
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
                  DP_G12 = DOT_PRODUCT(G1,G2)
                  DINTMIN = DSQ1+DSQ2-2.0D0*DP_G12
               END IF

               ! if the denominator in the d^2 is approx zero, we do not need to check for an internal minimum
               IF (DINTMIN.LT.DINTTEST) THEN
                  NOINT = .TRUE.
                  D1 = SQRT(DSQ1); D2 = SQRT(DSQ1)
                  G1(1:3) = G1(1:3)/D1; G2(1:3) = G2(1:3)/D2
               ELSE
                  !TODO: write MINMAX routine
                  CALL MINMAXD2R()
               END IF
               ! TODO: continue line 242
            END DO 
         END DO


      END SUBROUTINE GET_REPULSION_E


      SUBROUTINE GET_CCLOCAL(CCLOCAL)
         USE QCICONSTRAINTS, ONLY: CONCUTLOCAL, CONCUTABST, CONCUTABS, &
                                   CONCUTFRACT, CONCUTFRAC, CONDISTREFLOCAL
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(OUT) :: CCLOCAL
         !QUERY: this seems odd - why are they applied consecutively?
         CCLOCAL=CONCUTLOCAL(J2)
         IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
         IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
      END SUBROUTINE GET_CCLOCAL
END MODULE CONSTR_E_GRAD
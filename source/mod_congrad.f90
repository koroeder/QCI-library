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
         USE REPULSION, ONLY: 
         USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGE+2))   ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: GGG(3*NATOMS*(NIMAGE+2))  ! gradient for each atom in each image
         REAL(KIND = REAL64), INTENT(OUT) :: EEE(NIMAGE+2), EREP       ! energy for repulsions
         
         INTEGER :: J1, OFFSET1, OFFSET2
         REAL(KIND=REAL64) :: X1(3*NATOMS), X2(3*NATOMS)           ! coordinates for image 1 and 2


         DO J1=2,NIMAGE+1
            ! get coordinates for images
            OFFSET2=(3*NATOMS)*(J1-1)
            OFFSET1=(3*NATOMS)*(J1-2)
            X1(1:3*NATOMS)=XYZ(OFFSET1+1:OFFSET1+3*NATOMS)
            X2(1:3*NATOMS)=XYZ(OFFSET2+1:OFFSET2+3*NATOMS) 
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
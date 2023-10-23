!
! This version of congrad tests for an internal minimum in the
! constraint distances as well as the repulsions.
!
SUBROUTINE CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,FREEZENODEST, INTIMAGE, ATOMACTIVE, KINT, IMSEPMAX, &
  &            INTFREEZET, INTFROZEN, REPI, REPJ, CONCUT, CONCUTLOCAL, &
  &            CONCUTABS, CONCUTABST, CONCUTFRAC, CONCUTFRACT, INTMINFAC, INTSPRINGACTIVET, CHECKCONINT, JMAXCON, &
  &            NCONOFF, CONOFFTRIED, INTCONMAX, KINTENDS, ECON, EREP, ESPRING, &
  &  CONVERGECONTEST, CONVERGEREPTEST, KINTSCALED
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NCONINT(INTIMAGE+2),NREPINT(INTIMAGE+2),JMAX,IMAX,JMAXNOFF,IMAXNOFF
DOUBLE PRECISION :: ETOTAL, RMS, EMAX, EMAXNOFF
INTEGER JJMAX(INTIMAGE+2)
DOUBLE PRECISION  EEMAX(INTIMAGE+2)
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),CONEINT(INTIMAGE+2),RMSIMAGE(INTIMAGE+2)
LOGICAL NOINT, LPRINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2)
LOGICAL IMGFREEZE(INTIMAGE), PRINTE
DOUBLE PRECISION DPLUS, SPGRAD(3), CCLOCAL, DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq
DOUBLE PRECISION CONDMAX, CONREFMAX, CONCUTMAX

PRINTE=.FALSE.
111 CONTINUE

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
NCONINT(1:INTIMAGE+2)=0
NREPINT(1:INTIMAGE+2)=0
REPEINT(1:INTIMAGE+2)=0.0D0
CONEINT(1:INTIMAGE+2)=0.0D0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
LPRINT=.TRUE.
LPRINT=.FALSE.
!
!  Constraint energy and forces.
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 1 and 2 to images J1-1 and J1 corresponding to J1-2 and J1-1 below.
!
! IMGFREEZE(1:INTIMAGE) refers to the images excluding end points!
!
EMAX=-1.0D200
EMAXNOFF=-1.0D200
EEMAX(1:INTIMAGE+2)=-1.0D200
JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
JMAXNOFF=-1
IMAXNOFF=-1
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
   CCLOCAL=CONCUTLOCAL(J2)
   IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
   IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
!  IF (INTFROZEN(CONI(J2)).AND.INTFROZEN(CONJ(J2))) THEN
!     WRITE(*, '(A,I6,A,2I6)') ' congrad2> ERROR *** constraint ',J2,' between frozen atoms ',CONI(J2),CONJ(J2)
!     STOP
!  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
   DO J1=2,INTIMAGE+2
      IF (FREEZENODEST) THEN ! IMGFREEZE is not allocated otherwise!
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=2 and IMGFREEZE(1)=T cycle'
               CYCLE
            ENDIF
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=INTIMAGE+2 and IMGFREEZE(INTIMAGE)=T cycle'
               CYCLE
            ENDIF
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) THEN
!              IF (J2.EQ.1) PRINT '(A,I6,A)','J1=',J1,' IMGFREEZE(J1-2)=T and IMGFREEZE(J1-1)=T cycle'
               CYCLE
            ENDIF
         ENDIF
      ENDIF
      NI1=(3*NATOMS)*(J1-2)+3*(CONI(J2)-1)
      NI2=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(CONJ(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)

      G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3)
      G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3)
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
!
! Is there an internal extremum?
!
      r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b

      IF ((.NOT.CHECKCONINT).OR.(r1apr2bmr2amr1bsq.LT.1.0D-50)) THEN
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2(D2,D1,DINT,DSQ2,DSQ1,G1,G2,G1INT,G2INT,NOINT,DEBUG,r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!
! Need to include both D2 and D1 contributions if they are both outside tolerance.
! Otherwise we get discontinuities if they are very close and swap over.
!
!     CONCUT=CONCUTFRAC*CONDISTREF(J2)
!
! terms for image J1 - non-zero derivatives only for J1. D2 is the distance for image J1.
!
!     IF (LPRINT) WRITE(*, '(A,I6,5G15.5)') &
! &       'J1,D2,D1,DINT,MIN diff,CONCUT=',J1,D2,D1,DINT,ABS(D2-CONDISTREFLOCAL(J2)),CCLOCAL
      IF ((ABS(D2-CONDISTREFLOCAL(J2)).GT.CCLOCAL).AND.(J1.LT.INTIMAGE+2)) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=D2
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
         IF (DUMMY.GT.EMAXNOFF) THEN
            IF (.NOT.CONOFFTRIED(J2)) THEN
               IMAXNOFF=J1
               JMAXNOFF=J2
               EMAXNOFF=DUMMY
            ENDIF
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         EEE(J1)=EEE(J1)+DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         ECON=ECON      +DUMMY
         IF (LPRINT) WRITE(*, '(A,4I6,G15.5)') 'min J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
  &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
! !
! ! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
! !
! ! terms for image J1-1 - non-zero derivatives only for J1-1. D1 is the distance for image J1-1.
! !
! !     IF (LPRINT) WRITE(*, '(A,I6,5G15.5)') &
! ! &       'J1,D2,D1,DINT,MAX diff,CCLOCAL=',J1,D2,D1,DINT,ABS(D1-CONDISTREFLOCAL(J2)),CCLOCAL
!       IF ((ABS(D1-CONDISTREFLOCAL(J2)).GT.CCLOCAL).AND.(J1.GT.2)) THEN  
!          DUMMY=D1-CONDISTREFLOCAL(J2)  
!          REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1(1:3)
!          DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
! !        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
! !           WRITE(*, '(A,2I6,2L5,G20.10)') 'A CONI,CONJ,INTFROZEN(CONI),INTFROZEN(CONJ),DUMMY=', &
! ! &                                       CONI(J2),CONJ(J2),INTFROZEN(CONI(J2)),INTFROZEN(CONJ(J2)),DUMMY
! !        ENDIF
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             JMAX=J2
!             EMAX=DUMMY
!          ENDIF
!          IF (DUMMY.GT.EMAXNOFF) THEN
!             IF (.NOT.CONOFFTRIED(J2)) THEN
!                IMAXNOFF=J1
!                JMAXNOFF=J2
!                EMAXNOFF=DUMMY
!             ENDIF
!          ENDIF
!          IF (DUMMY.GT.EEMAX(J1-1)) THEN
!             JJMAX(J1-1)=J2
!             EEMAX(J1-1)=DUMMY
!          ENDIF
!          EEE(J1-1)=EEE(J1-1)+DUMMY
!          CONE(J1-1)=CONE(J1-1)+DUMMY
!          ECON=ECON      +DUMMY
!          IF (LPRINT) WRITE(*, '(A,4I6,G15.5)') 'max J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
!   &         SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
!          GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!          GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
!       ENDIF
!
      IF (CHECKCONINT.AND.(.NOT.NOINT).AND.(ABS(DINT-CONDISTREFLOCAL(J2)).GT.CCLOCAL)) THEN
         DUMMY=DINT-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1INT(1:3)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2INT(1:3)
         DUMMY=INTMINFAC*INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
!        WRITE(*,'(A,3I7,9F13.5)') 'J1,NI1,NJ1,INTMINFAC,INTCONSTRAINTDEL,DUMMY,GGG=',J1,NI1,NJ1,INTMINFAC,INTCONSTRAINTDEL, &
! &            DUMMY, GGG(NI1+1:NI1+3),GGG(NJ1+1:NJ1+3)   
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'CONI,CONJ,INTFROZEN(CONI),INTFROZEN(CONJ),DUMMY=', &
! &                                       CONI(J2),CONJ(J2),INTFROZEN(CONI(J2)),INTFROZEN(CONJ(J2)),DUMMY
!        ENDIF
         ECON=ECON+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=DINT
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
         IF (DUMMY.GT.EMAXNOFF) THEN
            IF (.NOT.CONOFFTRIED(J2)) THEN
               IMAXNOFF=J1
               JMAXNOFF=J2
               EMAXNOFF=DUMMY
            ENDIF
         ENDIF
         IF (DUMMY.GT.EEMAX(J1-1)) THEN
            JJMAX(J1-1)=J2
            EEMAX(J1-1)=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            CONEINT(J1)=CONEINT(J1)+DUMMY
            NCONINT(J1)=NCONINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            CONEINT(J1)=CONEINT(J1)+DUMMY/2.0D0
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY/2.0D0
            NCONINT(J1)=NCONINT(J1)+1
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ENDIF
!        WRITE(*, '(A,4I6,G15.5)') 'in2 J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
!        WRITE(*, '(A,3I7,9G13.5)') 'J1,NI2,NJ2,INTMINFAC,INTCONSTRAINTDEL,DUMMY,GGG=',J1,NI2,NJ2,INTMINFAC,INTCONSTRAINTDEL, & 
! &            DUMMY, GGG(NI2+1:NI2+3),GGG(NJ2+1:NJ2+3)   
!        WRITE(*,'(A,2G20.10)') 'in intmin block EEE(J1),EEE(J1-1)=',EEE(J1),EEE(J1-1)
      ENDIF
   ENDDO
ENDDO
IF (JMAX.GT.0) THEN
   IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,2I6,A,G20.10,A,3G20.10)') ' congrad> Highest constraint contribution for any image in image ',IMAX, &
 & ' constraint ',JMAX, &
 &                              ' atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX,' d,ref,cutoff=',CONDMAX,CONREFMAX,CONCUTMAX

ENDIF
CONVERGECONTEST=EMAX/INTCONSTRAINTDEL
! IF (JMAXNOFF.GT.0) THEN
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad2> Highest constraint contribution never turned off for any image in image ',IMAXNOFF, &
!  & ' constraint ',JMAXNOFF, &
!  &                              ' atoms ',CONI(JMAX),CONJ(JMAX),' off=',NCONOFF
! ELSEIF (JMAX.GT.0) THEN
!    JMAXNOFF=JMAX
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad2> *** Using highest constraint contribution for any image, setting NCONOFF to 0'
!    CONOFFTRIED(1:INTCONMAX)=.FALSE.
!    NCONOFF=0
! ENDIF
JMAXCON=JMAXNOFF

! INTCONST=INTCONSTRAINREPCUT**13

EMAX=-1.0D200
EEMAX(1:INTIMAGE+2)=-1.0D200
JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
   DO J1=2,INTIMAGE+2
      IF (FREEZENODEST) THEN
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
!     IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
!        WRITE(*, '(A,I6,A,2I6)') ' congrad2> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
!        STOP
!     ENDIF
      NI1=(3*NATOMS)*(J1-2)+3*(NREPI(J2)-1)
      NI2=(3*NATOMS)*(J1-1)+3*(NREPI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(NREPJ(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(NREPJ(J2)-1)

      G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3) 
      G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3) 
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      DCUT=NREPCUT(J2)**2
      IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE ! don't look for an internal minimum if both repulsions outside cutoff
      r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
!
! Is there an internal extremum?
!
      r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b

      IF (r1apr2bmr2amr1bsq.LT.1.0D-50) THEN
!        D1=1.0D100; D2=1.0D100; 
         NOINT=.TRUE.  
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!??????????????????????????????????????????????????????????????????????????????
!
!  WHY ARE WE DOING both D2 and D1? Isn't this double counting? See CONGRAD routine above
!
!
! Skip image INTIMAGE+2 - no non-zero gradients on other images and no energy contributions.
!
      IF ((D2.LT.NREPCUT(J2)).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
!        D12=DSQ2**6
         D12=DSQ2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'R1 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
!        ENDIF
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G2(1:3)
!        WRITE(*, '(A,4I6,G15.5)') 'min J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
      DUMMY=0.0D0
!
! SURELY THIS BLOCK JUST DUPLICATES THE TERMS FOR THE D2 BLOCK?
!
! !
! ! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
! !
! !     IF ((D1.LT.INTCONSTRAINREPCUT).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
!       IF ((D1.LT.NREPCUT(J2)).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
! !        D12=DSQ1**6
!          D12=DSQ1
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*NREPCUT(J2))/INTCONST)
!          DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D1-3.0D0*NREPCUT(J2))/INTCONST)
!          EEE(J1-1)=EEE(J1-1)+DUMMY
! !        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
! !           WRITE(*, '(A,2I6,2L5,G20.10)') 'R2 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! ! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
! !        ENDIF
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             JMAX=J2
!             EMAX=DUMMY
!          ENDIF
!          IF (DUMMY.GT.EEMAX(J1-1)) THEN
!             JJMAX(J1-1)=J2
!             EEMAX(J1-1)=DUMMY
!          ENDIF
!          REPE(J1-1)=REPE(J1-1)+DUMMY
!          EREP=EREP+DUMMY
! !        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
!          DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
!          REPGRAD(1:3)=DUMMY*G1(1:3)
! !        WRITE(*, '(A,4I6,G15.5)') 'max J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! ! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
!          GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!          GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
!       ENDIF
      DUMMY=0.0D0
!     IF ((.NOT.NOINT).AND.(DINT.LT.INTCONSTRAINREPCUT)) THEN
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
         EREP=EREP+DUMMY
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'R3 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
!        ENDIF
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1-1)) THEN
            JJMAX(J1-1)=J2
            EEMAX(J1-1)=DUMMY
         ENDIF
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            REPEINT(J1)=REPEINT(J1)+DUMMY
            NREPINT(J1)=NREPINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            REPEINT(J1)=REPEINT(J1)+DUMMY/2.0D0
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY/2.0D0
            NREPINT(J1)=NREPINT(J1)+1
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ENDIF
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO
IF (JMAX.GT.0) THEN
   WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad2> Highest repulsive  contribution for any image in image ',IMAX, &
 &  ' pair index ', &
 &                                JMAX,' atoms ',NREPI(JMAX),NREPJ(JMAX)
ENDIF
CONVERGEREPTEST=EMAX/INTCONSTRAINTREP
!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
EMAX=0.0D0
IMAX=0
IF (KINT.NE.0.0D0) THEN
   DO J1=1,INTIMAGE+1 ! sum over edges from J1 to J1+1
      NI1=(3*NATOMS)*(J1-1)
      NI2=(3*NATOMS)*J1
!
!  Edge between J1 and J1+1
!
      DPLUS=0.0D0
      DO J2=1,NATOMS
         IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
            DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
         ENDIF
!        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
      ENDDO
      DPLUS=SQRT(DPLUS)
!     IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2/KINTSCALED
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            EMAX=DUMMY
         ENDIF
!        DUMMY=KINT*0.5D0*DPLUS**2
         ESPRING=ESPRING+DUMMY
         IF (DUMMY.GT.EEMAX(J1+1)) THEN
            EEMAX(J1+1)=DUMMY
         ENDIF

!        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
!        DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DUMMY=KINT/KINTSCALED
         DO J2=1,NATOMS
            IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
               SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
               GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
               GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
            ENDIF
         ENDDO
!     ENDIF
   ENDDO
ENDIF
! IF (KINTENDS.GT.0.0D0) THEN
! !
! ! Extra terms for the two fixed end points.
! !
!    DO J1=2,INTIMAGE+1 ! sum over images
!       NI1=(3*NATOMS)*(J1-1)
!       NI2=(3*NATOMS)*(INTIMAGE+1)
! !
! !  Spring between 1 and J1
! !
!       DPLUS=0.0D0
!       DO J2=1,NATOMS
!          IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!             DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(3*(J2-1)+1))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(3*(J2-1)+2))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(3*(J2-1)+3))**2
!          ENDIF
! !        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
!       ENDDO
!       DPLUS=SQRT(DPLUS)
! !     IF (DPLUS.GT.IMSEPMAX) THEN
! !        DUMMY=KINTENDS*0.5D0*(DPLUS-IMSEPMAX)**2
!          DUMMY=KINTENDS*0.5D0*DPLUS**2
!          DUMMY=DUMMY*(1.0D0/(J1-1))**2 ! (INTIMAGE/(J1-1))**2
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             EMAX=DUMMY
!          ENDIF
! !        DUMMY=KINTENDS*0.5D0*DPLUS**2
!          ESPRING=ESPRING+DUMMY
!          IF (DUMMY.GT.EEMAX(J1)) THEN
!             EEMAX(J1)=DUMMY
!          ENDIF
! 
! !        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
! !        DUMMY=KINTENDS*(DPLUS-IMSEPMAX)/DPLUS
!          DUMMY=KINTENDS
!          DO J2=1,NATOMS
!             IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!                SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(3*(J2-1)+1:3*(J2-1)+3))
!                GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
!             ENDIF
!          ENDDO
! !     ENDIF
! !
! !  Spring between INTIMAGE+2 and J1
! !
!       DPLUS=0.0D0
!       DO J2=1,NATOMS
!          IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!             DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
!          ENDIF
! !        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
!       ENDDO
!       DPLUS=SQRT(DPLUS)
! !     IF (DPLUS.GT.IMSEPMAX) THEN
! !        DUMMY=KINTENDS*0.5D0*(DPLUS-IMSEPMAX)**2
!          DUMMY=KINTENDS*0.5D0*DPLUS**2
!          DUMMY=DUMMY*(INTIMAGE/(INTIMAGE+2-J1))**2
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             EMAX=DUMMY
!          ENDIF
! !        DUMMY=KINTENDS*0.5D0*DPLUS**2
!          ESPRING=ESPRING+DUMMY
!          IF (DUMMY.GT.EEMAX(J1)) THEN
!             EEMAX(J1)=DUMMY
!          ENDIF
! 
! !        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
! !        DUMMY=KINTENDS*(DPLUS-IMSEPMAX)/DPLUS
!          DUMMY=KINTENDS
!          DO J2=1,NATOMS
!             IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!                SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
!                GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
!             ENDIF
!          ENDDO
! !     ENDIF
!    ENDDO
! 
! ENDIF
WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad2> Highest spring  contribution for any image in image ',IMAX
IF (DEBUG) WRITE(*, '(A,3G20.10)') 'congrad2> ECON,EREP,ESPRING=',ECON,EREP,ESPRING
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients on locally frozen atoms to zero.
!
IF (INTFREEZET) THEN
   DO J1=2,INTIMAGE+1
      DO J2=1,NATOMS
         IF (INTFROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
IF (INTIMAGE.GT.0) THEN
   GGG(1:(3*NATOMS))=0.0D0
   GGG((INTIMAGE+1)*(3*NATOMS)+1:(INTIMAGE+2)*(3*NATOMS))=0.0D0
ENDIF
RMS=0.0D0
RMSIMAGE(1:INTIMAGE+2)=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,(3*NATOMS)
      RMSIMAGE(J1)=RMSIMAGE(J1)+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
   RMS=RMS+RMSIMAGE(J1)
   IF (LPRINT) WRITE(*, '(A,I6,2G20.10)') ' congrad2> J1,EEE,RMSIMAGE=',J1,EEE(J1),RMSIMAGE(J1)
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!

ETOTAL=0.0D0
MAXINT=-1.0D100
MININT=1.0D100
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
!  IF (DEBUG) PRINT '(A,I6,A,4G15.5)',' congrad2> con/rep and con/rep int image ', &
! &      J1,' ',CONE(J1),REPE(J1),CONEINT(J1),REPEINT(J1)
   IF (CONEINT(J1)+REPEINT(J1).LT.MININT) THEN
      MININT=CONEINT(J1)+REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (CONEINT(J1)+REPEINT(J1).GT.MAXINT) THEN
      MAXINT=CONEINT(J1)+REPEINT(J1)
      NMAXINT=J1
   ENDIF
ENDDO
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> largest  internal energy=',MAXINT,' for image ',NMAXINT
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> smallest internal energy=',MININT,' for image ',NMININT
IF (INTIMAGE.EQ.0) ETOTAL=EEE(1)+EEE(2)

! IF ((RMS.LT.1.0D-50).AND.(ETOTAL.GT.0.1D0).AND.(INTIMAGE.GT.0.0D0)) THEN
!    WRITE(*, '(A,2G20.10)') 'ETOTAL,RMS=',ETOTAL,RMS
!    WRITE(*, '(A,G20.10)') 'ECON=',ECON
!    WRITE(*, '(A,G20.10)') 'EREP=',EREP
!    IF (PRINTE) STOP
!    PRINTE=.TRUE.
!    GOTO 111
! ENDIF

END SUBROUTINE CONGRAD2
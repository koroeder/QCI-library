SUBROUTINE CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,INTIMAGE, KINT, IMSEPMAX, ATOMACTIVE, JMAXCON, &
  &            INTFREEZET, INTFROZEN, CONCUTLOCAL, CONCUT, CONCUTABST, CONCUTABS, CONCUTFRACT, CONCUTFRAC, &
  &  FREEZENODEST, INTSPRINGACTIVET, INTMINFAC, NCONOFF, CONOFFLIST, CONOFFTRIED, INTCONMAX, ECON, EREP, ESPRING, &
  &  CONVERGECONTEST, CONVERGEREPTEST, FCONTEST, FREPTEST, QCIINTREPMINSEP, QCIAVDEV, KINTSCALED
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
USE PORFUNCS
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,ISTAT,MYUNIT,JMAX,IMAX,JMAXNOFF,IMAXNOFF,NMAXINT,NMININT
DOUBLE PRECISION :: ETOTAL, RMS, EMAX, EMAXNOFF, FMAX, FMIN, SEPARATION
INTEGER OFFSET1, OFFSET2
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), D12, DSQ2, DSQ1, DSQI
LOGICAL NOINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2), CCLOCAL
DOUBLE PRECISION GGGR((3*NATOMS)*(INTIMAGE+2)), EREP1,EREP2
LOGICAL IMGFREEZE(INTIMAGE)
DOUBLE PRECISION DPLUS, SPGRAD(3), DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq,CUTMAX,DISTMAX
DOUBLE PRECISION CONDMAX, CONREFMAX, CONCUTMAX, DUMMY2, CCLOCAL2, DVEC(INTIMAGE+1), DEVIATION(INTIMAGE+1)
DOUBLE PRECISION GLOCAL1(3*NATOMS), GLOCAL2(3*NATOMS), XYZ1(3*NATOMS), XYZ2(3*NATOMS)

EEE(1:INTIMAGE+2)=0.0D0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
NMAXINT=0; NMININT=0 ! not actually used
MYUNIT=6

EMAX=-1.0D200
FMAX=-1.0D200
FMIN=1.0D200
! EEMAX(1:INTIMAGE+2)=-1.0D200
! JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
EMAXNOFF=-1.0D200
JMAXNOFF=-1
IMAXNOFF=-1
IF (INTCONSTRAINTDEL.EQ.0.0D0) GOTO 531
!
!  Constraint energy and forces.
!
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
      CCLOCAL=CONCUTLOCAL(J2)
      IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
      IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 2 refers to image J1.
!
   DO J1=2,INTIMAGE+1
!  DO J1=1,INTIMAGE+2  ! checking for zero!
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
      NI1=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      DUMMY=D2-CONDISTREFLOCAL(J2)  
      DUMMY2=DUMMY**2
      CCLOCAL2=CCLOCAL**2
!
! Reduced form for penalty function and gradient - multiply through by 1/(CCLOCAL**2*INTCONSTRAINTDEL)
! Should save CCLOCAL**2
! Could save DUMMY**2-CCLOCAL**2 outside IF block and base the test on difference of squares
!
      IF (DUMMY2.GT.CCLOCAL2) THEN 
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
!
!        REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
!        DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
!
         REPGRAD(1:3)=2*(DUMMY2-CCLOCAL2)*DUMMY*G2(1:3)/CCLOCAL2**2
         DUMMY=(DUMMY2-CCLOCAL2)**2/(2.0D0*CCLOCAL2**2)
         EEE(J1)=EEE(J1)  +DUMMY
         ECON=ECON        +DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=D2
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
!        IF (DUMMY.GT.EMAXNOFF) THEN
!           IF (.NOT.CONOFFTRIED(J2)) THEN
!              IMAXNOFF=J1
!              JMAXNOFF=J2
!              EMAXNOFF=DUMMY
!           ENDIF
!        ENDIF
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         DUMMY2=MINVAL(REPGRAD)
         IF (DUMMY2.LT.FMIN) FMIN=DUMMY2
         DUMMY2=MAXVAL(REPGRAD)
         IF (DUMMY2.GT.FMAX) FMAX=DUMMY2
      ENDIF
   ENDDO
ENDDO
IF (-FMIN.GT.FMAX) FMAX=-FMIN
FCONTEST=FMAX
IF (JMAX.GT.0) THEN
   IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,2I6,A,G15.5,A,3G15.5,A,G15.5)') ' congrad> Highest constraint for image ',IMAX, &
 & ' con ',JMAX, ' atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX,' d,ref,cutoff=',CONDMAX,CONREFMAX,CONCUTMAX, &
 & ' max grad=',FMAX

ENDIF
CONVERGECONTEST=EMAX

! IF (JMAXNOFF.GT.0) THEN
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,G20.10,A,I8)') ' congrad> Highest constraint contribution never turned off for any image in image ',IMAXNOFF, &
!  & ' constraint ',JMAXNOFF, &
!  &                              ' atoms ',CONI(JMAXNOFF),CONJ(JMAXNOFF),' value=',EMAXNOFF,' off=',NCONOFF
! ELSEIF (JMAX.GT.0) THEN
!    JMAXNOFF=JMAX
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad> *** Using highest constraint contribution for any image, setting NCONOFF to 0'
!    CONOFFTRIED(1:INTCONMAX)=.FALSE.
!    NCONOFF=0
! ENDIF
! JMAXCON=JMAXNOFF
531 CONTINUE

! GGG(1:(3*NATOMS))=0.0D0                            ! can delete when loop range above changes
! GGG((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=0.0D0 ! can delete when loop range above changes

! INTCONST=INTCONSTRAINREPCUT**13

EMAX=-1.0D200
FMAX=-1.0D200
FMIN=1.0D200
JMAX=-1
IMAX=-1

GGGR(1:3*NATOMS*(INTIMAGE+2))=0.0D0
IF (INTCONSTRAINTREP.EQ.0.0D0) GOTO 654
DO J1=2,INTIMAGE+1
! DO J1=2,INTIMAGE+2 ! we don't do anything for INTIMAGE+2 any more
!  DO J1=1,INTIMAGE+2 ! can change when zero energies are confirmed for end images
   IF (FREEZENODEST) THEN
      IF (J1.EQ.2) THEN
         IF (IMGFREEZE(1)) CYCLE
!     ELSE IF (J1.EQ.INTIMAGE+2) THEN
!        IF (IMGFREEZE(INTIMAGE)) CYCLE
      ELSE
         IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
      ENDIF
   ENDIF
   OFFSET2=(3*NATOMS)*(J1-1)
   OFFSET1=(3*NATOMS)*(J1-2)
   XYZ1(1:3*NATOMS)=XYZ(OFFSET1+1:OFFSET1+3*NATOMS)
   XYZ2(1:3*NATOMS)=XYZ(OFFSET2+1:OFFSET2+3*NATOMS)
   GLOCAL1(1:3*NATOMS)=0.0D0
   GLOCAL2(1:3*NATOMS)=0.0D0
   EREP1=0.0D0
   EREP2=0.0D0
   DO J2=1,NNREPULSIVE
! !  INTCONST=NREPCUT(J2)**13
!    INTCONST=NREPCUT(J2)**3
!      IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
!!        WRITE(*, '(A,I6,A,2I6)') ' congrad> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
!         STOP
!      ENDIF

!     NI1=OFFSET1+3*(NREPI(J2)-1)
!     NI2=OFFSET2+3*(NREPI(J2)-1)
!     NJ1=OFFSET1+3*(NREPJ(J2)-1)
!     NJ2=OFFSET2+3*(NREPJ(J2)-1)

      NI1=3*(NREPI(J2)-1)
      NI2=3*(NREPI(J2)-1)
      NJ1=3*(NREPJ(J2)-1)
      NJ2=3*(NREPJ(J2)-1)

!     G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3)
!     G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3)

      G1(1:3)=XYZ1(NI1+1:NI1+3)-XYZ1(NJ1+1:NJ1+3)
      G2(1:3)=XYZ2(NI2+1:NI2+3)-XYZ2(NJ2+1:NJ2+3)
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
      IF (ABS(NREPI(J2)-NREPJ(J2)).LT.QCIINTREPMINSEP) THEN ! don't check for internal minimum in distance - atoms too close for chain crossing.
         r1apr2bmr2amr1bsq=0.0D0
      ELSE
!        WRITE(*,'(A,I6,A,I6,A,2G20.10,A,G20.10)') 'distances in ',J1-1,' and ',J1,' are ',SQRT(DSQ1),SQRT(DSQ2),' cutoff=',NREPCUT(J2)
!        WRITE(*,'(A,I6,A,2I6,A,6G15.7)') 'image ',J1-1,' atoms ',NREPI(J2),NREPJ(J2),' coords ',XYZ(NI1+1:NI1+3),XYZ(NJ1+1:NJ1+3)
!        WRITE(*,'(A,I6,A,2I6,A,6G15.7)') 'image ',J1  ,' atoms ',NREPI(J2),NREPJ(J2),' coords ',XYZ(NI2+1:NI2+3),XYZ(NJ2+1:NJ2+3)
         r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
         r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b
      ENDIF
!
! Is the denominator of the d^2 internal minimum term close to zero?
!
      IF (r1apr2bmr2amr1bsq.LT.1.0D-50) THEN
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!
! Skip image INTIMAGE+2 - no non-zero gradients on other images and no energy contributions.
!
! Multiply energy and gradient by NREPCUT**2/INTCONSTRAINTREP to put values on a common scale.
!
!     IF ((D2.LT.NREPCUT(J2)).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
      IF (D2.LT.NREPCUT(J2)) THEN ! terms for image J1 - non-zero derivatives only for J1
!        D12=DSQ2**6
         D12=DSQ2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)*NREPCUT(J2)**2
         DUMMY=NREPCUT(J2)**2/D12+2.0D0*D2/NREPCUT(J2)-3.0D0  
!        EEE(J1)=EEE(J1)+DUMMY
         EREP1=EREP1+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            CUTMAX=NREPCUT(J2)
            DISTMAX=D2
            EMAX=DUMMY
         ENDIF
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*(1.0D0/(D2*D12)-1.0D0/INTCONST)*NREPCUT(J2)**2

         DUMMY=-2.0D0*(NREPCUT(J2)**2/(D2*D12)-1.0D0/NREPCUT(J2))
         REPGRAD(1:3)=DUMMY*G2(1:3)
!        GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
!        GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
         GLOCAL2(NI2+1:NI2+3)=GLOCAL2(NI2+1:NI2+3)+REPGRAD(1:3)
         GLOCAL2(NJ2+1:NJ2+3)=GLOCAL2(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!
! For internal minima we are counting edges. 
! Edge J1 is between images J1-1 and J1, starting from J1=2.
! Energy contributions are shared evenly, except for
! edge 1, which was assigned to image 2, and edge INTIMAGE+1, which
! was assigned to image INTIMAGE+1. Gradients are set to zero for the end images.
! 20/11/17 - changed to turn off internal minimum checks involving the fixed endpoints. DJW
!
      DUMMY=0.0D0
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2)).AND.(J1.NE.2)) THEN
!        D12=DSQI**6
         D12=DSQI
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTMINFAC*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)*NREPCUT(J2)**2
         DUMMY=INTMINFAC*(NREPCUT(J2)**2/D12+2.0D0*DINT/NREPCUT(J2)-3.0D0)
!        IF (J1.EQ.2) THEN
!           EEE(J1)=EEE(J1)+DUMMY
!        ELSE ! IF (J1.LT.INTIMAGE+2) THEN
!           EEE(J1)=EEE(J1)+DUMMY/2.0D0
!           EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
!        ENDIF
         EREP2=EREP2+DUMMY
         EREP=EREP+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            DISTMAX=DINT
            CUTMAX=NREPCUT(J2)
         ENDIF
! !        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*(1.0D0/(DINT*D12)-1.0D0/INTCONST)*NREPCUT(J2)**2
         DUMMY=-2.0D0*(NREPCUT(J2)**2/(DINT*D12)-1.0D0/NREPCUT(J2))
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1-1
!

!        GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!        GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         GLOCAL1(NI1+1:NI1+3)=GLOCAL1(NI1+1:NI1+3)+REPGRAD(1:3)
         GLOCAL1(NJ1+1:NJ1+3)=GLOCAL1(NJ1+1:NJ1+3)-REPGRAD(1:3)
         DUMMY2=MINVAL(REPGRAD)
!
! Gradient contributions for image J1
!
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
!        GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
!        GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
         GLOCAL2(NI2+1:NI2+3)=GLOCAL2(NI2+1:NI2+3)+REPGRAD(1:3)
         GLOCAL2(NJ2+1:NJ2+3)=GLOCAL2(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
   GGGR(OFFSET1+1:OFFSET1+3*NATOMS)=GGGR(OFFSET1+1:OFFSET1+3*NATOMS)+GLOCAL1(1:3*NATOMS)
   GGGR(OFFSET2+1:OFFSET2+3*NATOMS)=GGGR(OFFSET2+1:OFFSET2+3*NATOMS)+GLOCAL2(1:3*NATOMS)
   EEE(J1)=EEE(J1)+EREP1
   IF (J1.EQ.2) THEN
      EEE(J1)=EEE(J1)+EREP2
   ELSE ! IF (J1.LT.INTIMAGE+2) THEN
      EEE(J1)=EEE(J1)+EREP2/2.0D0
      EEE(J1-1)=EEE(J1-1)+EREP2/2.0D0
   ENDIF
ENDDO
FMIN=MINVAL(GGGR(3*NATOMS+1:3*NATOMS*(INTIMAGE+1)))
FMAX=MAXVAL(GGGR(3*NATOMS+1:3*NATOMS*(INTIMAGE+1)))
IF (-FMIN.GT.FMAX) FMAX=-FMIN
FREPTEST=FMAX

GGG(1:3*NATOMS*(INTIMAGE+2))=GGG(1:3*NATOMS*(INTIMAGE+2))+GGGR(1:3*NATOMS*(INTIMAGE+2))

IF (JMAX.GT.0) THEN
   IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,2I6,A,G15.5,A,2G15.5,A,G15.5)') ' congrad> Highest repulsion  for image ',IMAX, &  
 &  ' ind ',JMAX,' atoms ',NREPI(JMAX),NREPJ(JMAX),' value=',EMAX,' d,cutoff=',DISTMAX,CUTMAX, &
 &  ' max grad=',FMAX
ENDIF
CONVERGEREPTEST=EMAX
654 CONTINUE

!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
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
      ENDDO
      DPLUS=SQRT(DPLUS)
      DVEC(J1)=DPLUS
!     IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2/KINTSCALED
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            EMAX=DUMMY
         ENDIF
         ESPRING=ESPRING+DUMMY
!        IF (J1.EQ.1) THEN
!           EEE(2)=EEE(2)+DUMMY
!        ELSEIF (J1.EQ.INTIMAGE+1) THEN
!           EEE(INTIMAGE+1)=EEE(INTIMAGE+1)+DUMMY
!        ELSE
!           EEE(J1)=EEE(J1)+DUMMY/2.0D0
!           EEE(J1+1)=EEE(J1+1)+DUMMY/2.0D0
!        ENDIF
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
SEPARATION=SUM(DVEC(1:INTIMAGE+1))
DEVIATION(1:INTIMAGE+1)=ABS(100*((INTIMAGE+1)*DVEC(1:INTIMAGE+1)/SEPARATION-1.0D0))
QCIAVDEV=SUM(DEVIATION)/(INTIMAGE+1)
IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad> Highest spring  contribution for any image in image ',IMAX
IF (DEBUG) THEN
!  WRITE(*, '(A,3G20.10)') ' congrad> ECON,EREP,ESPRING=',ECON,EREP,ESPRING
!  WRITE(*,'(A)') '   edge         gap                deviation      '
!  WRITE(*,'(I6,3X,G20.10,G20.10)') (J1,DVEC(J1),DEVIATION(J1),J1=1,INTIMAGE+1)
   WRITE(*, '(A,2G20.10)') ' congrad> mean gap and mean deviation=',SEPARATION/(INTIMAGE+1),QCIAVDEV
ENDIF
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
GGG(1:(3*NATOMS))=0.0D0
GGG((INTIMAGE+1)*(3*NATOMS)+1:(INTIMAGE+2)*(3*NATOMS))=0.0D0
RMS=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,(3*NATOMS)
      RMS=RMS+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!
ETOTAL=0.0D0
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
ENDDO

END SUBROUTINE CONGRAD











!
! Call this version for additional repulsive terms between constraints.
!
SUBROUTINE CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,INTIMAGE, KINT, IMSEPMAX, ATOMACTIVE, QCINOREPINT, &
  &            INTFREEZET, INTFROZEN, CONCUTLOCAL, CONCUT, CONCUTABST, CONCUTABS, CONCUTFRACT, CONCUTFRAC, &
  &  INTMINFAC, FREEZENODEST, INTSPRINGACTIVET, QCIADDREP, QCIADDREPCUT, QCIADDREPEPS, QCIBONDS, KINTSCALED
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
USE PORFUNCS
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NREPINT(INTIMAGE+2),ISTAT,J3,J4,J5,NINTMIN,NINTMIN2,MYUNIT
DOUBLE PRECISION :: ECON, EREP, ETOTAL, RMS
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),RMSIM(INTIMAGE+2)
LOGICAL NOINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2), CCLOCAL
LOGICAL IMGFREEZE(INTIMAGE)
DOUBLE PRECISION DPLUS, ESPRING, SPGRAD(3), X1, X2, Y1, Y2, Z1, Z2, DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
REPEINT(1:INTIMAGE+2)=0.0D0
NREPINT(1:INTIMAGE+2)=0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
MYUNIT=6

IF (QCIADDREP.LT.1) THEN
   WRITE(*,'(A,I6)') 'congrad3> ERROR congrad3 called with no QCIADDREP=',QCIADDREP
   STOP
ENDIF
!
!  Constraint energy and forces.
!
OPEN(UNIT=852,FILE='test.xyz',STATUS='UNKNOWN')
INTCONST=QCIADDREPCUT**3
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
   CCLOCAL=CONCUTLOCAL(J2)
   IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
   IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
   DO J1=2,INTIMAGE+1
      NI1=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF (ABS(D2-CONDISTREFLOCAL(J2)).GT.CCLOCAL) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
         EEE(J1)=EEE(J1)  +DUMMY
         ECON=ECON        +DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
      ENDIF
!     WRITE(MYUNIT,'(A,2I6,5G20.10)') 'J1,J2,D2,CONDISTREFLOCAL,CCLOCAL,EEE,CONE=',J1,J2,D2,CONDISTREFLOCAL(J2),CCLOCAL,EEE(J1),CONE(J1)
      IF (J2.GT.QCIBONDS) CYCLE
      DO J3=J2+1,QCIBONDS
         IF (.NOT.CONACTIVE(J3)) CYCLE
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         NI2=(3*NATOMS)*(J1-1)+3*(CONI(J3)-1)
         NJ2=(3*NATOMS)*(J1-1)+3*(CONJ(J3)-1)
         DO J4=1,QCIADDREP
            X1=(J4*XYZ(NI1+1)+(QCIADDREP+1-J4)*XYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*XYZ(NI1+2)+(QCIADDREP+1-J4)*XYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*XYZ(NI1+3)+(QCIADDREP+1-J4)*XYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*XYZ(NI2+1)+(QCIADDREP+1-J5)*XYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*XYZ(NI2+2)+(QCIADDREP+1-J5)*XYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*XYZ(NI2+3)+(QCIADDREP+1-J5)*XYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.QCIADDREPCUT) THEN 
                  D12=D2**2
                  DUMMY=QCIADDREPEPS*(1.0D0/D12+(2.0D0*D2-3.0D0*QCIADDREPCUT)/INTCONST)
                  EEE(J1)=EEE(J1)+DUMMY
                  REPE(J1)=REPE(J1)+DUMMY
                  EREP=EREP+DUMMY
!                 WRITE(*,'(A,4I6,A,2I6,A,2G20.10)') 'congrad3> atoms: ',CONI(J2),CONJ(J2),CONI(J3),CONJ(J3), &
! &                     ' sites ',J4,J5,' dist,erep=',D2,DUMMY   
                  DUMMY=-2.0D0*QCIADDREPEPS*(1.0D0/(D2*D12)-1.0D0/INTCONST)
                  G2(1)=(X1-X2)/D2
                  G2(2)=(Y1-Y2)/D2
                  G2(3)=(Z1-Z2)/D2
                  REPGRAD(1:3)=DUMMY*G2(1:3)
                  GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+J4*REPGRAD(1:3)/(QCIADDREP+1.0D0) ! forces on the four atoms involved in the two constraints
                  GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)+(QCIADDREP+1-J4)*REPGRAD(1:3)/(QCIADDREP+1.0D0)
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)-J5*REPGRAD(1:3)/(QCIADDREP+1.0D0)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-(QCIADDREP+1-J5)*REPGRAD(1:3)/(QCIADDREP+1.0D0)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
CLOSE(852)

GGG(1:(3*NATOMS))=0.0D0                                      ! can delete when loop range above changes
GGG((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=0.0D0 ! can delete when loop range above changes

! INTCONST=INTCONSTRAINREPCUT**13

DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
   DO J1=2,INTIMAGE+2
!  DO J1=1,INTIMAGE+2 ! can change when zero energies are confirmed for end images
      IF (FREEZENODEST) THEN
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
      IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
         WRITE(*, '(A,I6,A,2I6)') ' congrad3> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
         STOP
      ENDIF
!     WRITE(*,'(A,2I8,6G20.10)') 'congrad3> B J1,J2,GGG(1:6)=',J1,J2,GGG(1:6)
      NI2=(3*NATOMS)*(J1-1)+3*(NREPI(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(NREPJ(J2)-1)
      R2AX=XYZ(NI2+1); R2AY=XYZ(NI2+2); R2AZ=XYZ(NI2+3)
      R2BX=XYZ(NJ2+1); R2BY=XYZ(NJ2+2); R2BZ=XYZ(NJ2+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF (D2.LT.NREPCUT(J2)) THEN ! term for image J1
!        D12=D2**12
         D12=D2**2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=DUMMY*G2(1:3)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!     WRITE(MYUNIT,'(A,2I6,4G20.10)') 'J1,J2,D2,NREPCUT,EEE,REPE=',J1,J2,D2,NREPCUT(J2),EEE(J1),REPE(J1)
!
! For internal minima we are counting edges. 
! Edge J1 is between images J1-1 and J1, starting from J1=2.
! Energy contributions are shared evenly, except for
! edge 1, which is assigned to image 2, and edge INTIMAGE+1, which
! is assigned to image INTIMAGE+1. Gradients are set to zero for
! the end images.
!
      IF (J1.EQ.1) CYCLE
      IF (QCINOREPINT) CYCLE
      NI1=(3*NATOMS)*(J1-2)+3*(NREPI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(NREPJ(J2)-1)

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

!     IF (r2ax**2+r2ay**2+r2az**2+r2bx**2+r2by**2+r2bz**2-2*(r2ax*r2bx+r2ay*r2by+r2az*r2bz).EQ.0.0D0) THEN
      IF (r1apr2bmr2amr1bsq.LT.1.0D-10) THEN
!        D1=1.0D100; D2=1.0D100;
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq) 
         IF (.NOT.NOINT) THEN
!           WRITE(*,'(A,I6,A,I6,A,2I6,A,2G20.10)') 'congrad3> internal minimum images ',J1-1,' and ',J1,' atoms: ',NREPI(J2),NREPJ(J2), &
! &                        ' distance,cutoff=',DINT,NREPCUT(J2)
            NINTMIN=NINTMIN+1
         ENDIF
      ENDIF
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
         NINTMIN2=NINTMIN2+1
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
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
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1-1
!
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1
!
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO
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
      IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2/KINTSCALED
         ESPRING=ESPRING+DUMMY
!        DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DUMMY=KINT/KINTSCALED
         DO J2=1,NATOMS
            IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
               SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
               GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
               GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
            ENDIF
         ENDDO
      ENDIF
   ENDDO
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
   RMSIM(J1)=0.0D0
   DO J2=1,(3*NATOMS)
      RMS=RMS+GGG((3*NATOMS)*(J1-1)+J2)**2
      RMSIM(J1)=RMSIM(J1)+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
   RMSIM(J1)=SQRT(RMSIM(J1)/(3*NATOMS))
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
!  WRITE(*, '(A,I6,A,3G20.10)') ' congrad3> con/rep/RMS image ',J1,' ',CONE(J1),REPE(J1),RMSIM(J1)
   IF (REPEINT(J1).LT.MININT) THEN
      MININT=REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (REPE(J1).GT.MAXINT) THEN
      MAXINT=REPE(J1)
      NMAXINT=J1
   ENDIF
ENDDO

END SUBROUTINE CONGRAD3
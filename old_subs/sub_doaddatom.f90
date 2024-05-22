!SUBROUTINE DOADDATOM(NCONSTRAINT,NTRIES,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE,AABACK,BACKDONE)
SUBROUTINE DOADDATOM(NCONSTRAINT,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE,AABACK,BACKDONE,INLIST)
USE KEY, ONLY : CONACTIVE, CONI, CONJ, ATOMACTIVE, CONDISTREF, REPI, REPJ, REPCUT, INTREPSEP,  &
  &             INTCONSTRAINREPCUT, NREPULSIVE, NREPMAX, MAXCONUSE, CHECKCONINT, &
  &             FREEZENODEST, NNREPULSIVE, INTFROZEN, ATOMSTORES, QCIADDACIDT, &
  &             NREPULSIVEFIX, REPIFIX, REPJFIX, REPCUTFIX, NREPI, NREPJ, NREPCUT, MAXNACTIVE, &
  &             NCONSTRAINTFIX, CONIFIX, CONJFIX, INTCONCUT, INTCONSEP, QCIRADSHIFTT, QCIRADSHIFT, &
  &             CONOFFTRIED, QCIADDREP, DOBACK, DOBACKALL, QCITRILAT, QCILINEARLIST, QCILINEARN, QCILINEARATOMS
USE COMMONS, ONLY: NATOMS, DEBUG
IMPLICIT NONE
INTEGER INTIMAGE
INTEGER NBEST, NBEST2, NCONTOACTIVE(NATOMS),  NCONSTRAINT, J2, NEWATOM,  CONLIST(NATOMS), N1, N2, N3, &
  &     NTOADD, NADDED, NMININT, NMAXINT, TURNONORDER(NATOMS), NDUMMY, J1, J3, NITERDONE, NCONFORNEWATOM, NACTIVE!, NTRIES(NATOMS)
DOUBLE PRECISION DUMMY, DUMMY2, CONDIST(NATOMS), DMIN!, RANDOM, DPRAND
INTEGER ACID!, NDFORNEWATOM, BESTPRESERVEDN(NATOMS)
DOUBLE PRECISION BESTCLOSESTD(NATOMS), INVDTOACTIVE(NATOMS)!, BESTPRESERVEDD(NATOMS)
LOGICAL IMGFREEZE(INTIMAGE), CHOSENACID, AABACK(NATOMS), BACKDONE, FTEST, IDONE, INLIST(NATOMS)
DOUBLE PRECISION P1(3), P2(3), P3(3), SOL1(3), SOL2(3), R1, R2, R3
DOUBLE PRECISION C1, C2, C3, VEC1(3), VEC2(3), VEC3(3), ESAVED, ESAVEC, ESAVE0
DOUBLE PRECISION C4, C5, C6, VEC4(3), VEC5(3), VEC6(3), F1, F2
INTEGER NCFORNEWATOM, BESTCLOSESTN(NATOMS), NNREPSAVE, NREPSAVE
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), XSAVED(3,INTIMAGE+2), XSAVEC(3,INTIMAGE+2), XSAVE0(3,INTIMAGE+2),FRAC,&
  &              RMS,EEE(INTIMAGE+2),GGG((3*NATOMS)*(INTIMAGE+2)),ETOTAL,DS,DF,D1SQ,D2SQ!, DNORM, RAN1
INTEGER  REPITEMP(NREPMAX), REPJTEMP(NREPMAX)
DOUBLE PRECISION REPCUTTEMP(NREPMAX)


NTOADD=1
IF (QCILINEARLIST) NTOADD=QCILINEARN-2 ! we should already have two done
NADDED=0
CHOSENACID=.FALSE.

IF (DOBACK.AND.(.NOT.BACKDONE)) THEN
   DO J1=1,NATOMS
      IF (AABACK(J1)) THEN
         IF (.NOT.ATOMACTIVE(J1)) GOTO 763
      ENDIF
   ENDDO
   IF (DEBUG) WRITE(*,'(A,I6,A)') ' intlbfgs> All backbone atoms are active'
   BACKDONE=.TRUE.
ENDIF

763   CONTINUE

!
! ATOMSTORES(J1) is the residue for atom J1
!
! AMBER12_GET_RESDATA needs a data type in amber12 interface and the number of residues
! defined in amber12interface.f90
!

!
! Save current number of repulsions and number that are active to speed up the
! calls to CHECKREP
!
NNREPSAVE=NNREPULSIVE
NREPSAVE=NREPULSIVE
542   CONTINUE
!     DUMMY=1.0D100

!
! Initial guess for new active atom position. This is crucial for success in INTCONSTRAINT schemes!
!
      ESAVED=1.0D100
      ESAVE0=1.0D100
      ESAVEC=1.0D100
      FTEST=.TRUE.
      IDONE=.FALSE.
      IF (NCONFORNEWATOM.GE.3) THEN
         IDONE=.TRUE.
!
! Move the new atom consistently in the local environment of its three nearest actively constrained atoms.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
!        N1=NCONFORNEWATOM-2; N2=NCONFORNEWATOM-1; N3=NCONFORNEWATOM
         N1=1; N2=2; N3=3
         IF (DEBUG) WRITE(*,'(A,3I6)') ' intlbfgs> initial guess from closest three constrained active atoms, ',CONLIST(N1:N3)
!        IF (DEBUG) WRITE(*,'(A,3I6)') ' intlbfgs> initial guess from furthest three constrained active atoms, ',CONLIST(N1:N3)
         VEC1(1:3)=XYZ(3*(CONLIST(N2)-1)+1:3*(CONLIST(N2)-1)+3)-XYZ(3*(CONLIST(N1)-1)+1:3*(CONLIST(N1)-1)+3)
         DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
         IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
         VEC2(1:3)=XYZ(3*(CONLIST(N3)-1)+1:3*(CONLIST(N3)-1)+3)-XYZ(3*(CONLIST(N1)-1)+1:3*(CONLIST(N1)-1)+3)
         DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
         VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
         DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
         IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
         VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
         VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
         VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
         C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC1(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC1(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC1(3)
         C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC2(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC2(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC2(3)
         C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC3(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC3(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC3(3)

         DO J1=2,INTIMAGE+1
            VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+3) &
  &                  -XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+3) &
  &                  -XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)+C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)

!
! Alternative analytical solution from intersection of three spheres by trilateration
!
            IF (QCITRILAT) THEN
               P1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
               P2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+3)
               P3(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+3)
               R1=CONDIST(N1)
               R2=CONDIST(N2)
               R3=CONDIST(N3)
               CALL TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
               IF (FTEST) THEN
!                 WRITE(*,'(A,I8)') ' intlbfgs> WARNING *** no trilateration solution for image ',J1
               ELSE

!
! Try minmum distance from previous solution
!
                  D1SQ=(SOL1(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                   +(SOL1(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                   +(SOL1(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
                  D2SQ=(SOL2(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                   +(SOL2(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                   +(SOL2(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
!                 WRITE(*,'(A,2F20.10)') 'D1SQ,D2SQ=',D1SQ,D2SQ
                  IF (D1SQ.LT.D2SQ) THEN
                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL1(1:3)
                  ELSE
                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL2(1:3)
                  ENDIF
               ENDIF
            ENDIF

         ENDDO
         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
         IF (QCIADDREP.GT.0) THEN
            CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSEIF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         ESAVE0=ETOTAL
         DO J1=2,INTIMAGE+1
            XSAVE0(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
         ENDDO
      ENDIF

      IF ((.NOT.IDONE).AND.(NCFORNEWATOM.GE.3)) THEN
         IDONE=.TRUE.
!
! Choose three atoms from the BESTCLOSEST list at random with bias towards the
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
         DUMMY=0.0D0

         N1=1; N2=2; N3=3
         IF (DEBUG) WRITE(*,'(A,3I6,A)') ' intlbfgs> choosing positions ',N1,N2,N3,' in closest list'
!
! Define axes for starting point.
!
         VEC1(1:3)=XYZ(3*(BESTCLOSESTN(N2)-1)+1:3*(BESTCLOSESTN(N2)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
         DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
         IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
         VEC2(1:3)=XYZ(3*(BESTCLOSESTN(N3)-1)+1:3*(BESTCLOSESTN(N3)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
         DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
         VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
         DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
         IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
         VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
         VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
         VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
         C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC1(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC1(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC1(3)
         C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC2(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC2(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC2(3)
         C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC3(1)+ &
  &         (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC3(2)+ &
  &         (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC3(3)
!
! Define axes for end point.
!
         VEC4(1:3)=XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N2)-1)+1:3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N2)-1)+3) &
  &               -XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+1:3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+3)
         DUMMY=SQRT(VEC4(1)**2+VEC4(2)**2+VEC4(3)**2)
         IF (DUMMY.NE.0.0D0) VEC4(1:3)=VEC4(1:3)/DUMMY
         VEC5(1:3)=XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N3)-1)+1:3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N3)-1)+3) &
  &               -XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+1:3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+3)
         DUMMY=VEC4(1)*VEC5(1)+VEC4(2)*VEC5(2)+VEC4(3)*VEC5(3)
         VEC5(1:3)=VEC5(1:3)-DUMMY*VEC4(1:3)
         DUMMY=SQRT(VEC5(1)**2+VEC5(2)**2+VEC5(3)**2)
         IF (DUMMY.NE.0.0D0) VEC5(1:3)=VEC5(1:3)/DUMMY
         VEC6(1)= VEC4(2)*VEC5(3)-VEC4(3)*VEC5(2)
         VEC6(2)=-VEC4(1)*VEC5(3)+VEC4(3)*VEC5(1)
         VEC6(3)= VEC4(1)*VEC5(2)-VEC4(2)*VEC5(1)
         C4=(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+1))*VEC4(1)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+2))*VEC4(2)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+3))*VEC4(3)
         C5=(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+1))*VEC5(1)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+2))*VEC5(2)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+3))*VEC5(3)
         C6=(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+1))*VEC6(1)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+2))*VEC6(2)+ &
  &         (XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(BESTCLOSESTN(N1)-1)+3))*VEC6(3)
!        WRITE(*,'(A,6G20.10)') 'C1,C4,C2,C5,C3,C6=',C1,C4,C2,C5,C3,C6

         DO J1=2,INTIMAGE+1
!
! Analogous axes in image J1.
!
            VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3) &
  &                  -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3) &
  &                  -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            F1=(INTIMAGE+2.0D0-J1)/(INTIMAGE+1.0D0) ! fractions to use from the two end points.
            F2=(J1-1.0D0)/(INTIMAGE+1.0D0)

            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &         XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)+ &
  &                (C1*F1+C4*F2)*VEC1(1:3)+(C2*F1+C5*F2)*VEC2(1:3)+(C3*F1+C6*F2)*VEC3(1:3)

!
! Alternative analytical solution from intersection of three spheres by trilateration
!
            IF (QCITRILAT) THEN
               P1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               P2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3)
               P3(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3)
               R1=BESTCLOSESTD(N1)
               R2=BESTCLOSESTD(N2)
               R3=BESTCLOSESTD(N3)
               CALL TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
               IF (FTEST) THEN
!              WRITE(*,'(A,I8)') ' intlbfgs> WARNING *** no trilateration solution for image ',J1
               ELSE

!
! Try minmum distance from previous solution
!
                  D1SQ=(SOL1(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                   +(SOL1(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                   +(SOL1(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
                  D2SQ=(SOL2(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                   +(SOL2(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                   +(SOL2(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
!              WRITE(*,'(A,2F20.10)') 'D1SQ,D2SQ=',D1SQ,D2SQ
                  IF (D1SQ.LT.D2SQ) THEN
                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL1(1:3)
                  ELSE
                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL2(1:3)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
         IF (QCIADDREP.GT.0) THEN
            CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSEIF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         ESAVEC=ETOTAL
         DO J1=2,INTIMAGE+1
            XSAVEC(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
         ENDDO
      ENDIF
!
! Standard linear interpolation, with constraint distance scaled by FRAC.
! Works for FRAC as small as 0.1 with repulsion turned off.
! We use an appropriately weighted displacement from atom CONLIST(1) using the displacements
! in the two end points.
!
      ETOTAL=1.0D100
      IF ((.NOT.IDONE).OR.QCILINEARLIST) THEN
         FRAC=1.0D0
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1)  &
 &         +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))/(INTIMAGE+1) &
 &+(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+1))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+2)  &
 &         +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))/(INTIMAGE+1) &
 &+(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+2))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)  &
 &         +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))/(INTIMAGE+1) &
 &+(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+3))/(INTIMAGE+1)
         ENDDO
         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
         IF (QCIADDREP.GT.0) THEN
            CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSEIF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ENDIF

      IF (DEBUG) WRITE(*,'(A,4G15.5)') ' intlbfgs> energies for constrained, preserved, closest, and linear schemes=', &
  &           ESAVE0,ESAVED,ESAVEC,ETOTAL
 
      IF (QCILINEARLIST) THEN
         IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> using linear interpolation for atom in linearlist'
      ELSEIF ((ETOTAL.LT.ESAVEC).AND.(ETOTAL.LT.ESAVED).AND.(ETOTAL.LT.ESAVE0)) THEN
         IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from linear interpolation'
      ELSE IF ((ESAVEC.LT.ESAVED).AND.(ESAVEC.LT.ESAVE0)) THEN
         IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using closest atoms'
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVEC(1:3,J1)
         ENDDO
         ETOTAL=ESAVEC
      ELSE IF (ESAVED.LT.ESAVE0) THEN
         IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using preserved distances'
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVED(1:3,J1)
         ENDDO
         ETOTAL=ESAVED
      ELSE 
         IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> interpolation using closest constraints'
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVE0(1:3,J1)
         ENDDO
         ETOTAL=ESAVE0
      ENDIF

      NADDED=NADDED+1
      IF (NADDED.LT.NTOADD) GOTO 542
!
! Check whether we've added all atoms in the amino acid corresponding to the new atom. If not, go back to the top
! and choose the next candidate.
!
      IF (QCIADDACIDT.AND.(.NOT.DOBACK)) THEN
         IF (NCONFORNEWATOM.LT.3) THEN
            WRITE(*,'(A,I6,A)') 'doaddatom> Added atom constraints ',NCONFORNEWATOM,' is < 3; not attempting to add another atom now'
         ELSE
            DO J1=1,NATOMS
               IF ((ATOMSTORES(J1).EQ.ACID).AND.(.NOT.(ATOMACTIVE(J1)))) GOTO 542
            ENDDO
            WRITE(*,'(A,I6,A)') 'doaddatom> All atoms of residue ',ACID,' are active'
         ENDIF
      ENDIF
!
! Check whether we've added all atoms in the backbone of the amino acid corresponding to the new atom. If not, go back to the top
! and choose the next candidate.
!
      IF (DOBACKALL.AND.(AABACK(NEWATOM))) THEN
         DO J1=1,NATOMS
            IF ((ATOMSTORES(J1).EQ.ACID).AND.(.NOT.(ATOMACTIVE(J1))).AND.AABACK(J1)) GOTO 542
         ENDDO
         WRITE(*,'(A,I6,A)') 'doaddatom> All backbone atoms of residue ',ACID,' are active'
      ENDIF

      IF (QCIRADSHIFTT) THEN
         WRITE(*,'(A,F15.5)') ' intlbfgs> Applying radial shift for unconstrained atoms of ',QCIRADSHIFT
         WRITE(*,'(20I6)') CONLIST(1:NCONFORNEWATOM)
         DO J1=2,INTIMAGE+1
            scaleloop: DO J2=1,NATOMS
               IF (.NOT.ATOMACTIVE(J2)) CYCLE scaleloop
               IF (J2.EQ.NEWATOM) CYCLE scaleloop
               DO J3=1,NCONFORNEWATOM
                  IF (CONLIST(J3).EQ.J2) CYCLE scaleloop
               ENDDO
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)- &
   &                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)*QCIRADSHIFT/DUMMY
               XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)= &
   &           XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)+VEC1(1:3)

            ENDDO scaleloop
         ENDDO
      ENDIF
!
! Turn frozen images off for new added atom.
!
!     IF (DEBUG) WRITE(*,'(A)') ' intlbfgs> turning off frozen images'
!     IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
      CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
!
! need a new gradient since the active atom has changed !
!
      IF (QCIADDREP.GT.0) THEN
         CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSEIF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      QCILINEARLIST=.FALSE.
      IF (ALLOCATED(QCILINEARATOMS)) THEN
         DEALLOCATE(QCILINEARATOMS)
      ENDIF
   RETURN

END SUBROUTINE DOADDATOM
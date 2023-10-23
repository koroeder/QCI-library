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
      NBEST=0
      NCONTOACTIVE(1:NATOMS)=0
      INVDTOACTIVE(1:NATOMS)=0.0D0
      DO J2=1,NCONSTRAINT
         IF (CONACTIVE(J2)) CYCLE     ! count new, inactive constraints
         IF (CONOFFTRIED(J2)) CYCLE   ! if we've tried turning it off, it must actually be active. Don't try again.
         IF (ATOMACTIVE(CONI(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONJ(J2))) THEN
               NCONTOACTIVE(CONJ(J2))=NCONTOACTIVE(CONJ(J2))+1
               IF (1.0D0/CONDISTREF(J2).GT.INVDTOACTIVE(CONJ(J2))) INVDTOACTIVE(CONJ(J2))=1.0D0/CONDISTREF(J2)
!              INVDTOACTIVE(CONJ(J2))=INVDTOACTIVE(CONJ(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (ATOMACTIVE(CONJ(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONI(J2))) THEN
               NCONTOACTIVE(CONI(J2))=NCONTOACTIVE(CONI(J2))+1
!              INVDTOACTIVE(CONI(J2))=INVDTOACTIVE(CONI(J2))+1.0D0/CONDISTREF(J2)
               IF (1.0D0/CONDISTREF(J2).GT.INVDTOACTIVE(CONI(J2))) INVDTOACTIVE(CONI(J2))=1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (NCONTOACTIVE(CONI(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONI(J2))
         ENDIF
         IF (NCONTOACTIVE(CONJ(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONJ(J2))
         ENDIF
      ENDDO

!
!  Choose NEWATOM deterministically. Take the inactive atom with the shortest constrained distance
!  subject to maximum constraints.
!  If QCIADDACIDT try adding all atoms in a give residue until we hit < 3 constraints. 
!  However, we must allow such atoms if there is nothing else to try.
!  In this case CHOSENACID will be false the first time through, and we should not set it true.
!  If CHOSENACID is true, but < 3 constraints, add that atom but don't try another.
!  Should enable us to add one atom at a time later in the run where necessary.
!
!  If QCILINEARLIST is true we are going to add all the linear interpolation atoms first.
!  This option overrides CHOSENACID.
!  QCILINEARLIST is set to false once we have added all the atoms in INLIST
!
      DUMMY2=1.0D100
      NBEST2=0
     
      DO J1=1,NCONSTRAINT
         IF (CONACTIVE(J1)) CYCLE

         IF (ATOMACTIVE(CONJ(J1))) THEN
            IF (CHOSENACID.AND.(.NOT.(ATOMSTORES(CONI(J1)).EQ.ACID)).AND.(.NOT.QCILINEARLIST)) THEN
            ELSE
               IF (.NOT.ATOMACTIVE(CONI(J1))) THEN
                  IF (DOBACK.AND.(.NOT.BACKDONE).AND.(.NOT.AABACK(CONI(J1))).AND.(.NOT.QCILINEARLIST)) THEN
                  ELSEIF (QCILINEARLIST.AND.(.NOT.INLIST(CONI(J1)))) THEN
                  ELSE
                     IF (NCONTOACTIVE(CONI(J1)).GE.NBEST2) THEN
                        IF (CONDISTREF(J1).LT.DUMMY2) THEN
                           DUMMY2=CONDISTREF(J1)
                           NEWATOM=CONI(J1)
                           NBEST2=NCONTOACTIVE(CONI(J1))
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ATOMACTIVE(CONI(J1))) THEN
            IF (CHOSENACID.AND.(.NOT.(ATOMSTORES(CONJ(J1)).EQ.ACID)).AND.(.NOT.QCILINEARLIST)) THEN
            ELSE
               IF (.NOT.ATOMACTIVE(CONJ(J1))) THEN
                  IF (DOBACK.AND.(.NOT.BACKDONE).AND.(.NOT.AABACK(CONJ(J1))).AND.(.NOT.QCILINEARLIST)) THEN
                  ELSEIF (QCILINEARLIST.AND.(.NOT.INLIST(CONJ(J1)))) THEN
                  ELSE
                     IF (NCONTOACTIVE(CONJ(J1)).GE.NBEST2) THEN
                        IF (CONDISTREF(J1).LT.DUMMY2) THEN
                           DUMMY2=CONDISTREF(J1)
                           NEWATOM=CONJ(J1)
                           NBEST2=NCONTOACTIVE(CONJ(J1))
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF (NEWATOM*NBEST2.EQ.0) THEN ! sanity check
         WRITE(*,'(A,2I6)') ' intlbfgs> ERROR *** new active atom not set NEWATOM,NBEST2=', NEWATOM,NBEST2
         STOP
      ENDIF
      IF (DEBUG) WRITE(*,'(3(A,I8),I8,A,F15.5)') ' intlbfgs> Choosing new active atom ',NEWATOM,' new constraints=', &
  &            NCONTOACTIVE(NEWATOM),' maximum constraints available and possible=',NBEST,NBEST2,' shortest constraint=',DUMMY2

      IF (QCIADDACIDT.AND.(.NOT.CHOSENACID).AND.(.NOT.DOBACK)) THEN
         ACID=ATOMSTORES(NEWATOM)
         CHOSENACID=.TRUE.
      ENDIF
      IF ((.NOT.CHOSENACID).AND.DOBACKALL) THEN
         ACID=ATOMSTORES(NEWATOM)
         CHOSENACID=.TRUE.
      ENDIF
          
!
!  We need a sorted list of up to 3 active atoms, sorted according to how well the
!  end point distance is preserved, even if they don't satisfy the constraint 
!  condition. We want three atoms to use for a local axis system in the interpolation.
!
!  Try sorting on the shortest average distances in the endpoint structures instead, to avoid
!  problems with distant atoms acidentally having a well-preserved distance.
!

      IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.

      NCFORNEWATOM=0
      BESTCLOSESTD(1:NATOMS)=1.0D100
      DO J1=1,NATOMS
         IF (ABS(J1-NEWATOM).GT.INTCONSEP) CYCLE
         IF (.NOT.ATOMACTIVE(J1)) CYCLE
         DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &             +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &             +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
         DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &             +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &             +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
         IF (DS.GT.INTCONCUT) CYCLE
         IF (DF.GT.INTCONCUT) CYCLE
         DUMMY=(DS+DF)/2.0D0
         NCFORNEWATOM=NCFORNEWATOM+1
         DO J2=1,NCFORNEWATOM
            IF (DUMMY.LT.BESTCLOSESTD(J2)) THEN
               DO J3=NCFORNEWATOM,J2+1,-1
                  BESTCLOSESTD(J3)=BESTCLOSESTD(J3-1)
                  BESTCLOSESTN(J3)=BESTCLOSESTN(J3-1)
               ENDDO
               BESTCLOSESTD(J2)=DUMMY
               BESTCLOSESTN(J2)=J1
               GOTO 659
            ENDIF
         ENDDO
659      CONTINUE
      ENDDO
      IF (DEBUG) THEN
         WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' closest average distances in endpoints:'
         WRITE(*,'(20I6)') BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
         WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> sorted average distances:'
         WRITE(*,'(10G12.4)') BESTCLOSESTD(1:MIN(10,NCFORNEWATOM))
      ENDIF
!
!  Maintain a sorted list of active atoms that are constrained to the new atom, sorted
!  according to their distance.
!
      NCONFORNEWATOM=0
      CONDIST(1:NATOMS)=1.0D100
      IF (DEBUG) WRITE(*,'(3(A,I6))') ' intlbfgs> New active atom is number ',NEWATOM,' total=',NACTIVE+1, &
 &                     ' steps=',NITERDONE
      DO J1=1,NCONSTRAINT
         IF (CONACTIVE(J1)) CYCLE
         IF ((CONI(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONI(J1)))) THEN  
              NCONFORNEWATOM=NCONFORNEWATOM+1
!             CONACTIVE(J1)=.TRUE.
!             NITSTART(J1)=NITERDONE
!             NCONSTRAINTON=NCONSTRAINTON+1
! !
! ! The ...ON variables are not actually used in congrad.f90.
! !
!             CONDISTREFLOCALON(NCONSTRAINTON)=CONDISTREFLOCAL(J1)
!             CONDISTREFON(NCONSTRAINTON)=CONDISTREF(J1)
!             CONION(NCONSTRAINTON)=CONI(J1)
!             CONJON(NCONSTRAINTON)=CONJ(J1)
! 
!             IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
            IF (NCONFORNEWATOM.EQ.1) THEN
               CONDIST(1)=CONDISTREF(J1)
               IF (CONI(J1).EQ.NEWATOM) CONLIST(1)=CONJ(J1)
               IF (CONJ(J1).EQ.NEWATOM) CONLIST(1)=CONI(J1)
            ENDIF
            DO J2=1,NCONFORNEWATOM-1
               IF (CONDISTREF(J1).LT.CONDIST(J2)) THEN
!                 WRITE(*,'(A,I6,G12.4,I6,G12.4)') 'J1,CONDISTREF < J2,CONDIST: ',J1,CONDISTREF(J1),J2,CONDIST(J2)
                  DO J3=NCONFORNEWATOM,J2+1,-1
!                    WRITE(*,'(A,I6,A,I6,A,G12.4)') ' moving dist and list from ',J3-1,' to ',J3,' CONDIST=',CONDIST(J3-1)
                     CONDIST(J3)=CONDIST(J3-1)
                     CONLIST(J3)=CONLIST(J3-1)
                  ENDDO
                  CONDIST(J2)=CONDISTREF(J1)
!                 WRITE(*,'(A,I6,A,G12.4)') ' setting condist element ',J2,' to ',CONDISTREF(J1)
                  IF (CONI(J1).EQ.NEWATOM) CONLIST(J2)=CONJ(J1)
                  IF (CONJ(J1).EQ.NEWATOM) CONLIST(J2)=CONI(J1)
!                 WRITE(*,'(A,I6,A,G12.4)') ' setting conlist element ',J2,' to ',CONLIST(J2)
                  GOTO 654
               ENDIF
            ENDDO 
            CONDIST(NCONFORNEWATOM)=CONDISTREF(J1)
!           WRITE(*,'(A,I6,A,G12.4)') ' setting condist element ',NCONFORNEWATOM,' to ',CONDISTREF(J1)
            IF (CONI(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONJ(J1)
            IF (CONJ(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONI(J1)
!           WRITE(*,'(A,I6,A,G12.4)') ' setting conlist element ',NCONFORNEWATOM,' to ',CONLIST(NCONFORNEWATOM)
654       CONTINUE
         ENDIF
      ENDDO 
      IF (DEBUG) THEN
         WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' is constrained to ',NCONFORNEWATOM, &
  &                                    ' other active atoms:'
         WRITE(*,'(20I6)') CONLIST(1:NCONFORNEWATOM)
         WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> sorted distances:'
         WRITE(*,'(10G12.4)') CONDIST(1:NCONFORNEWATOM)
      ENDIF
      DO J1=1,MIN(MAXCONUSE,NCONFORNEWATOM)
         DO J2=1,NCONSTRAINT
            IF ((CONI(J2).EQ.NEWATOM).AND.(CONJ(J2).EQ.CONLIST(J1))) THEN
                  CONACTIVE(J2)=.TRUE.
                  IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
            ELSE IF ((CONJ(J2).EQ.NEWATOM).AND.(CONI(J2).EQ.CONLIST(J1))) THEN
                  CONACTIVE(J2)=.TRUE.
                  IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
            ENDIF
         ENDDO
      ENDDO

      DO J1=1,NATOMS
         IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
         IF (ABS(J1-NEWATOM).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
         DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
!           IF (.NOT.CONACTIVE(J2)) CYCLE ! repulsions for inactive constraints 
            IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.NEWATOM)).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.NEWATOM))) GOTO 543
         ENDDO
         DMIN=1.0D100
         DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
            DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &                 (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &                 (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
            IF (DF.LT.DMIN) DMIN=DF
         ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
         DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
         NREPULSIVE=NREPULSIVE+1
         IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
         REPI(NREPULSIVE)=J1
         REPJ(NREPULSIVE)=NEWATOM
         REPCUT(NREPULSIVE)=DMIN
!        IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',NEWATOM,' with atom ',J1, &
! &                                                ' cutoff=',DMIN
543      CONTINUE
      ENDDO
      ATOMACTIVE(NEWATOM)=.TRUE.
      NACTIVE=NACTIVE+1
      IF (MAXNACTIVE.EQ.0) MAXNACTIVE=NATOMS
!
! Freeze atoms that became active more than NACTIVE-MAXNACTIVE events ago.
! For example, with MAXNACTIVE=5 and 40 active atoms, we would freeze those 
! turned on first, second, up to the 35th in the TURNONORDER list.
!
      IF (NACTIVE.GT.MAXNACTIVE) THEN
         NDUMMY=TURNONORDER(NACTIVE-MAXNACTIVE)
         IF (INTFROZEN(NDUMMY)) THEN
            IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' doaddatom> Not turning off frozen active atom ',NDUMMY,' already frozen'
         ELSE
            IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' doaddatom> Freezing active atom ',NDUMMY
            INTFROZEN(NDUMMY)=.TRUE.
!
! Turn off constraints and repulsions between frozen atoms to save time and remove fixed components of
! energy and gradient.
!
            DO J2=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J2)) CYCLE
               IF (INTFROZEN(CONI(J2)).AND.INTFROZEN(CONJ(J2))) THEN
                  CONACTIVE(J2)=.FALSE.
                  WRITE(*,'(A,I6,A,2I6)') 'doaddatom> turning off constraint ',J2,' between atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO

            J2=0
!
! We can't use the FIX variables because they are no longer being set if we employ the 
! builtin AMBER topology.
!

            DO J1=1,NREPULSIVE
               IF (INTFROZEN(REPI(J1)).AND.INTFROZEN(REPJ(J1))) THEN
                  WRITE(*,'(A,I6,A,2I6)') 'doaddatom> turning off repulsion ',J1,' between atoms ',REPI(J1),' and ',REPJ(J1) 
                  CYCLE
               ENDIF
               IF (ATOMACTIVE(REPI(J1)).AND.ATOMACTIVE(REPJ(J1))) THEN
                  DO J3=1,NCONSTRAINT
!                    IF (.NOT.CONACTIVE(J3)) CYCLE ! no repulsions for any constraints
                     IF ((CONI(J3).EQ.REPI(J1)).AND.(CONJ(J3).EQ.REPJ(J1))) GOTO 962 ! no repulsion
                     IF ((CONI(J3).EQ.REPJ(J1)).AND.(CONJ(J3).EQ.REPI(J1))) GOTO 962 ! no repulsion
                  ENDDO
                  J2=J2+1
                  REPITEMP(J2)=REPI(J1)
                  REPJTEMP(J2)=REPJ(J1)
                  REPCUTTEMP(J2)=REPCUT(J1)
962               CONTINUE
               ENDIF
            ENDDO
            NREPULSIVE=J2
            WRITE(*,'(A,I6,A)') ' doaddatom> After allowing for frozen atoms there are ',NREPULSIVE,' possible repulsions'
            REPI(1:NREPULSIVE)=REPITEMP(1:NREPULSIVE)
            REPJ(1:NREPULSIVE)=REPJTEMP(1:NREPULSIVE)
            NREPI(1:NREPULSIVE)=REPI(1:NREPULSIVE)
            NREPJ(1:NREPULSIVE)=REPJ(1:NREPULSIVE)
            NNREPULSIVE=NREPULSIVE
            REPCUT(1:NREPULSIVE)=REPCUTTEMP(1:NREPULSIVE)
            NREPCUT(1:NREPULSIVE)=REPCUT(1:NREPULSIVE)
         ENDIF
      ENDIF

      NDUMMY=0
      DO J1=1,NATOMS
         IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
      ENDDO
      IF (NDUMMY.NE.NACTIVE) THEN
         WRITE(*,'(A,I6)') ' doaddatom> ERROR *** inconsistency in number of active atoms. ',NDUMMY,' should be ',NACTIVE
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) WRITE(*,'(A,I6)') ' active atom ',J1
         ENDDO
         STOP
      ENDIF

      TURNONORDER(NACTIVE)=NEWATOM
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
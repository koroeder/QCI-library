SUBROUTINE CHECKPERC(LXYZ,LINTCONSTRAINTTOL,NQCIFREEZE,NCPFIT)
USE KEY, ONLY : ATOMACTIVE, NCONSTRAINT, INTFROZEN, CONI, CONJ, CONDISTREF, INTCONMAX, INTCONSTRAINTTOL, &
  &             INTCONSEP, NCONGEOM, CONGEOM, CONIFIX, CONJFIX, CONDISTREFFIX, INTCONCUT, &
  &             NCONSTRAINTFIX, BULKT, TWOD, RIGIDBODY, CONDATT, CONCUT, CONCUTFIX, &
  &             BONDS, QCIAMBERT,QCISBT, QCIADDREP, QCIADDREPCUT, QCIBONDS, QCISECOND, &
  &             QCIHIRET, NBONDSOPEP, BONDSOPEP, HIRET
USE COMMONS, ONLY: NATOMS, DEBUG, PARAM1, PARAM2, PARAM3
USE HIRE_OPEP_INTERFACE_MOD, ONLY: OPEP_BOND_INFO
USE HIRE_INTERFACE, ONLY: GET_NBONDS, HIRE_BOND_INFO
IMPLICIT NONE
INTEGER NDIST1(NATOMS), NCYCLE, DMIN1, DMAX1, NUNCON1, J1, J2, J3, NQCIFREEZE, J4, NCPFIT, LUNIT, GETUNIT
INTEGER NI1, NJ1, NI2, NJ2, J5, ATOM1, ATOM2
DOUBLE PRECISION LINTCONSTRAINTTOL, MAXCONDIST, MINCONDIST, DS, DF, LXYZ((3*NATOMS)*2)
DOUBLE PRECISION DSMIN, DSMAX, DSMEAN, D, DIST2, RMAT(3,3), X1, Y1, Z1, X2, Y2, Z2, DMIN, D2
LOGICAL CHANGED, LDEBUG, CONFILET
LOGICAL :: CALLED=.FALSE.
SAVE CALLED
!for QCIAMBER
INTEGER NBOND, NDUMMY,NANGLE
LOGICAL CONTACTFILE
!QCISB

LINTCONSTRAINTTOL=INTCONSTRAINTTOL
IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
!
! Fixed constraints based on congeom file entries
! Just need to adjust the list based on any frozen atoms. We
! want to exclude any constraints between two frozen atoms 
! from the list, because subsequent code depends on this.
!

IF (NCONGEOM.GE.2) THEN
   IF (CALLED.OR.CONDATT) THEN
      J2=0
      DO J1=1,NCONSTRAINTFIX
!
! If called with two minima check that CONCUTFIX is large enough to
! accommodate the separation of the two atoms in both minima.
!
         IF (NCPFIT.EQ.2) THEN
            DF=MAX(ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ(3*(CONIFIX(J1)-1)+1)-LXYZ(3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+2)-LXYZ(3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+3)-LXYZ(3*(CONJFIX(J1)-1)+3))**2)),&
                   ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+1)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+2)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+3)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+3))**2)))
            IF (DF.GT.CONCUTFIX(J1)) THEN
               IF (ABS(DF-CONCUTFIX(J1)).GT.1.0D-6) &
  &                WRITE(*,'(A,2I5,3(A,G15.5))') ' checkperc> Increasing con cutoff atoms ', &
  &                CONIFIX(J1),CONJFIX(J1),' from ',CONCUTFIX(J1),' to ',DF,' ref=',CONDISTREFFIX(J1)
               CONCUTFIX(J1)=DF
            ENDIF
         ENDIF
         IF (INTFROZEN(CONIFIX(J1)).AND.INTFROZEN(CONJFIX(J1))) CYCLE
         J2=J2+1
         CONI(J2)=CONIFIX(J1)
         CONJ(J2)=CONJFIX(J1)
         CONDISTREF(J2)=CONDISTREFFIX(J1)
         CONCUT(J2)=CONCUTFIX(J1)
      ENDDO
      NCONSTRAINT=J2
      WRITE(*,'(A,I6,A)') ' checkperc> After allowing for frozen atoms there are ',NCONSTRAINT,' constraints'
      RETURN 
   ELSE
!
! Put reference minima in optimal permutational alignment with reference minimum one.
!
      DO J2=2,NCONGEOM
         LDEBUG=.FALSE.
         CALL MINPERMDIST(CONGEOM(1,1:3*NATOMS),CONGEOM(J2,1:3*NATOMS),NATOMS,LDEBUG, &
  &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
      ENDDO
   ENDIF
   ALLOCATE(CONIFIX(INTCONMAX),CONJFIX(INTCONMAX),CONCUTFIX(INTCONMAX),CONDISTREFFIX(INTCONMAX))
ENDIF

INQUIRE(FILE='constraintfile',EXIST=CONFILET)

51   NCONSTRAINT=0 
MAXCONDIST=-1.0D0
MINCONDIST=1.0D100

IF ((QCISBT).AND.(QCIAMBERT)) THEN
WRITE(*,'(A,2I6,A,I6)')  'intlbfgs> Both QCISB and QCIAMB set to true!'
STOP
ENDIF

!
! This block should be a subroutine...
!
IF (QCISBT) THEN 
  IF (DEBUG) WRITE(*,'(A,2I6,A,I6)')  'intlbfgs> QCISB set to true.'
 !sn513: Adding constraints for SBM Go models.From QCIAMBERT.
  NBOND=NATOMS-1
  NANGLE=NATOMS-2
  IF (ALLOCATED(BONDS)) DEALLOCATE(BONDS)
  ALLOCATE(BONDS(NBOND,2))
  WRITE(*,*) 'intlbfgs> Number of bonds (QCISB):',NBOND
  WRITE(*,*) 'intlbfgs> Number of angles(QCISB):',NANGLE
! Pairs of bonds
  NDUMMY=0
  DO J1=1,NBOND
   BONDS(J1,1)=J1
   BONDS(J1,2)=J1+1
   NDUMMY=NDUMMY+1
  ENDDO
 !Loop through bonds.
 DO J2=1,NBOND                !loop through all bonds and add them to constraint list
      IF (INTFROZEN(BONDS(J2,1)).AND.INTFROZEN(BONDS(J2,2))) CYCLE ! no constraints between intfrozen atoms
      NCONSTRAINT=NCONSTRAINT+1
      IF (DEBUG) WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> BONDS: Adding constraint for atoms ',BONDS(J2,1),BONDS(J2,2), &
  &                     '  total=',NCONSTRAINT
      DS=SQRT((LXYZ(3*(BONDS(J2,1)-1)+1)-LXYZ(3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+2)-LXYZ(3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+3)-LXYZ(3*(BONDS(J2,2)-1)+3))**2)
      DF=SQRT((LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+1)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+2)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+3)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+3))**2)
      IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
      CONI(NCONSTRAINT)=MIN(BONDS(J2,1),BONDS(J2,2))
      CONJ(NCONSTRAINT)=MAX(BONDS(J2,1),BONDS(J2,2))
      CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
      CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
      IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
      IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
     IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs>constrain distance for ',CONI(NCONSTRAINT), &
 &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
 &            ' # bond constraints=',NCONSTRAINT,'QCISB'
   ENDDO
   QCIBONDS=NCONSTRAINT


  write(*,*)'NCONSTRAINT=',NCONSTRAINT,'NCONMAX=',INTCONMAX
  IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') 'intlbfgs> END OF QCISB, BONDS'
  ! Add constraints for second-nearest neighbours - should correspond to bond angles
!
   DO J2=1,NBOND
      inloop1: DO J3=J2+1,NBOND
        IF (BONDS(J2,1).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,1).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,1)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,1)
        ELSE
           CYCLE inloop1
        ENDIF
        WRITE(*,*) ATOM1,ATOM2
        IF (INTFROZEN(ATOM1).AND.INTFROZEN(ATOM2)) CYCLE ! no constraints between intfrozen atom
        NCONSTRAINT=NCONSTRAINT+1
       IF (DEBUG) WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> ANGLES:Adding constraint for second neighbours ',ATOM1,ATOM2, &
  &                     '  total=',NCONSTRAINT
         DS=SQRT((LXYZ(3*(ATOM1-1)+1)-LXYZ(3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*(ATOM1-1)+2)-LXYZ(3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*(ATOM1-1)+3)-LXYZ(3*(ATOM2-1)+3))**2)
         DF=SQRT((LXYZ(3*NATOMS+3*(ATOM1-1)+1)-LXYZ(3*NATOMS+3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+2)-LXYZ(3*NATOMS+3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+3)-LXYZ(3*NATOMS+3*(ATOM2-1)+3))**2)
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=MIN(ATOM1,ATOM2)
         CONJ(NCONSTRAINT)=MAX(ATOM1,ATOM2)
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
!        WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,2I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
!        &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
!        &            ' # second neighbour constraints, total=',QCISECOND,NCONSTRAINT
!
        NDUMMY=NCONSTRAINT
      ENDDO inloop1
   ENDDO
   QCISECOND=NCONSTRAINT-QCIBONDS
   WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> First and second neighbour constraints: ',QCIBONDS,QCISECOND,' total: ',NCONSTRAINT
   IF (DEBUG) write(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') 'intlbfgs> END OF QCISB, ANGLES'
!STOP

!Read contacts to be constrained from constraintfile
!   contactfile='contacts.sbm'
! Note: Contacts are by definition not identical to bond or angle constraints passed on earlier. 
   INQUIRE(FILE='contacts.sbm',EXIST=CONTACTFILE)
   write(*,*) 'Reading constraints from pairs in',CONTACTFILE
   IF (CONTACTFILE) THEN
    IF(DEBUG) WRITE(*,*)'intlbfgs> FOUND contacts.sbm file. Adding contacts to constraintlist.'
    OPEN(LUNIT,FILE='contacts.sbm',STATUS='OLD')
    DO 
       READ(LUNIT,*,END=801) J2,J3
       write(*,*) J2,J3
       IF (J3-J2.GT.INTCONSEP) CYCLE
       IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE !no constraints between frozen atoms
       NCONSTRAINT=NCONSTRAINT+1
       DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2)
       DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2)
       IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
       CONI(NCONSTRAINT)=J2
       CONJ(NCONSTRAINT)=J3
       CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
       CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
       IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
       IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
       WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
  &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &            ' # bond constraints=',NCONSTRAINT
    ENDDO
      801 continue
    WRITE(*,'(A,I6,2(A,G15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT, &
  &                               ' shortest=',MINCONDIST,' longest=',MAXCONDIST
    CLOSE(LUNIT)
   ENDIF
ENDIF

IF (QCIAMBERT.OR.QCIHIRET) THEN
   IF (QCIAMBERT.AND.QCIHIRET) THEN
      WRITE(*,*) " intlbfgs> both QCIAMBER and QCIHIRE set"
      STOP
   ENDIF
   IF (QCIAMBERT) THEN        
      CALL TOPOLOGY_READER(NBOND)
   ELSE
      IF (HIRET)  THEN
         CALL GET_NBONDS(NBOND)
         IF (ALLOCATED(BONDS)) DEALLOCATE(BONDS)
         ALLOCATE(BONDS(NBOND,2))
         CALL HIRE_BOND_INFO(NBOND,BONDS)
      ELSE  
         CALL OPEP_BOND_INFO(NATOMS)
         NBOND = NBONDSOPEP
         IF (ALLOCATED(BONDS)) DEALLOCATE(BONDS)
         ALLOCATE(BONDS(NBOND,2))
         BONDS(1:NBOND,1:2) = BONDSOPEP(1:NBOND,1:2)
      ENDIF
   ENDIF 
!
!  kr366> assume we use two endpoints and topology for amber constraints
!  get number of bonds and bonds from topology
!  loop through all bonds and add them to constraint list
!
   DO J2=1,NBOND                !loop through all bonds and add them to constraint list
      IF (INTFROZEN(BONDS(J2,1)).AND.INTFROZEN(BONDS(J2,2))) CYCLE ! no constraints between intfrozen atoms
      NCONSTRAINT=NCONSTRAINT+1
      IF (DEBUG) WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> Adding constraint for atoms ',BONDS(J2,1),BONDS(J2,2), &
  &                     '  total=',NCONSTRAINT
      DS=SQRT((LXYZ(3*(BONDS(J2,1)-1)+1)-LXYZ(3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+2)-LXYZ(3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+3)-LXYZ(3*(BONDS(J2,2)-1)+3))**2)
      DF=SQRT((LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+1)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+2)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+3)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+3))**2)
      IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
      CONI(NCONSTRAINT)=MIN(BONDS(J2,1),BONDS(J2,2))
      CONJ(NCONSTRAINT)=MAX(BONDS(J2,1),BONDS(J2,2))
      CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
      CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
      IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
      IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
     IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
  &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &            ' # bond constraints=',NCONSTRAINT

   ENDDO
   QCIBONDS=NCONSTRAINT
!
! Add constraints for second-nearest neighbours - should correspond to bond angles
!
   DO J2=1,NBOND
      inloop: DO J3=J2+1,NBOND
        IF (BONDS(J2,1).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,1).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,1)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,1)
        ELSE
           CYCLE inloop
        ENDIF
        IF (INTFROZEN(ATOM1).AND.INTFROZEN(ATOM2)) CYCLE ! no constraints between intfrozen atoms
        NCONSTRAINT=NCONSTRAINT+1
        IF (DEBUG) WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> Adding constraint for second neighbours ',ATOM1,ATOM2, &
  &             '  total=',NCONSTRAINT
         DS=SQRT((LXYZ(3*(ATOM1-1)+1)-LXYZ(3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*(ATOM1-1)+2)-LXYZ(3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*(ATOM1-1)+3)-LXYZ(3*(ATOM2-1)+3))**2)
         DF=SQRT((LXYZ(3*NATOMS+3*(ATOM1-1)+1)-LXYZ(3*NATOMS+3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+2)-LXYZ(3*NATOMS+3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+3)-LXYZ(3*NATOMS+3*(ATOM2-1)+3))**2)
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=MIN(ATOM1,ATOM2)
         CONJ(NCONSTRAINT)=MAX(ATOM1,ATOM2)
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
        NDUMMY=NCONSTRAINT

      ENDDO inloop
   ENDDO
   QCISECOND=NCONSTRAINT-QCIBONDS
   WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> First and second neighbour constraints: ',QCIBONDS,QCISECOND,' total: ',NCONSTRAINT
   IF (CONFILET) THEN
      LUNIT=GETUNIT()
      OPEN(LUNIT,FILE='constraintfile',STATUS='OLD')
!
!  Additional amber constraints, e.g. cis/trans
!
      DO
         READ(LUNIT,*,END=534)  J2, J3
!
! Forbid constraints corresponding to atoms distant in sequence. Set INTCONSEP to number of sites to
! turn this off
!
         IF (J3-J2.GT.INTCONSEP) CYCLE
         DO J4=1,NCONSTRAINT  ! check for duplicates - would mess up interpolation based on closest constrained atoms
            IF ((J2.EQ.CONI(J4)).AND.(J3.EQ.CONJ(J4))) GOTO 679
            IF ((J2.EQ.CONJ(J4)).AND.(J3.EQ.CONI(J4))) GOTO 679
         ENDDO
         IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
         NCONSTRAINT=NCONSTRAINT+1
         DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2)
         DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2)
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> extra constraint distance for ',CONI(NCONSTRAINT), &
  &                     CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &                  ' # constraints=',NCONSTRAINT
679      CONTINUE
      ENDDO
534   CONTINUE
      CLOSE(LUNIT)
      IF (NCONSTRAINT-NDUMMY.GT.0) WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Extra distance constraints: ',NCONSTRAINT-NDUMMY
       WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT,' shortest=',MINCONDIST,&
  &       ' longest=',MAXCONDIST
       CLOSE(LUNIT)
      ENDIF
ELSE IF (CONFILET) THEN 
    LUNIT=GETUNIT()
    OPEN(LUNIT,FILE='constraintfile',STATUS='OLD')
!
!  Add constraint for this distance to the list.
!
    DO 
       READ(LUNIT,*,END=531)  J2, J3
!
! Forbid constraints corresponding to atoms distant in sequence. Set INTCONSEP to number of sites to 
! turn this off
!
       IF (J3-J2.GT.INTCONSEP) CYCLE 
       DO J4=1,NCONSTRAINT  ! check for duplicates - would mess up interpolation based on closest constrained atoms
          IF ((J2.EQ.CONI(J4)).AND.(J3.EQ.CONJ(J4))) GOTO 678
          IF ((J2.EQ.CONJ(J4)).AND.(J3.EQ.CONI(J4))) GOTO 678
       ENDDO
       IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
       NCONSTRAINT=NCONSTRAINT+1
!      WRITE(*,'(A,2I6,A,I6)') ' intlbfgs> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
       DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &           +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &           +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2) 
       DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &           +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &           +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2) 
       IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
       CONI(NCONSTRAINT)=J2
       CONJ(NCONSTRAINT)=J3
       CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
       CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
       IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
       IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
       WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
  &                 CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &                ' # constraints=',NCONSTRAINT
678    CONTINUE
    ENDDO
531 CONTINUE
    WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT, &
  &                               ' shortest=',MINCONDIST,' longest=',MAXCONDIST
    CLOSE(LUNIT)

ELSE IF (NCONGEOM.LT.2) THEN 
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS

         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
         DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2) 
         IF (DS.GT.INTCONCUT) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
         DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2) 
         IF (DF.GT.INTCONCUT) CYCLE ! don't allow constraints if either endpoint separation is too large DJW

         IF (ABS(DS-DF).LT.LINTCONSTRAINTTOL) THEN
!
!  Add constraint for this distance to the list.
!
            NCONSTRAINT=NCONSTRAINT+1
            IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
            CONI(NCONSTRAINT)=J2
            CONJ(NCONSTRAINT)=J3
            CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
            CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
            IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
            IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)

         ENDIF
      ENDDO
   ENDDO
   IF (DEBUG) WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT, &
  &                                     ' shortest=',MINCONDIST,' longest=',MAXCONDIST
ELSE
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS
         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         DSMIN=1.0D100
         DSMAX=-1.0D100
         DSMEAN=0.0D0
         DO J4=1,NCONGEOM
            DS=SQRT((CONGEOM(J4,3*(J2-1)+1)-CONGEOM(J4,3*(J3-1)+1))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+2)-CONGEOM(J4,3*(J3-1)+2))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+3)-CONGEOM(J4,3*(J3-1)+3))**2) 
            IF (DS.GT.DSMAX) DSMAX=DS
            IF (DS.LT.DSMIN) DSMIN=DS
            IF ((J4.GT.1).AND.(ABS(DSMIN-DSMAX).GT.LINTCONSTRAINTTOL)) GOTO 753 ! unconstrained
            IF (DS.GT.INTCONCUT) GOTO 753 ! don't allow constraints if any image separation is too large DJW
            DSMEAN=DSMEAN+DS
         ENDDO
!
!  Add constraint for this distance to the list if we make it to here.
!
         NCONSTRAINT=NCONSTRAINT+1
         WRITE(*,'(A,2I6,A,I6)') 'checkperc> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DSMAX+DSMIN)/2.0D0 
         CONCUT(NCONSTRAINT)=(DSMAX-DSMIN)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) WRITE(*,'(A,2I5,A,2F10.4,A,F12.4,A,I8)') &
  &                       ' checkperc> constrain atoms ',CONI(NCONSTRAINT), &
  &                       CONJ(NCONSTRAINT),' max, min ',DSMAX,DSMIN, &
  &                       ' cutoff=',CONCUT(NCONSTRAINT),' constraints=',NCONSTRAINT
753      CONTINUE
      ENDDO
   ENDDO
   CONIFIX(1:NCONSTRAINT)=CONI(1:NCONSTRAINT)
   CONJFIX(1:NCONSTRAINT)=CONJ(1:NCONSTRAINT)
   CONDISTREFFIX(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
   CONCUTFIX(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)
   NCONSTRAINTFIX=NCONSTRAINT
ENDIF

IF (QCIADDREP.GT.0) THEN
   DMIN=1.0D100
   DO J2=1,QCIBONDS
!
! end point 1
!
      NI1=3*(CONI(J2)-1)
      NJ1=3*(CONJ(J2)-1)
      DO J3=J2+1,QCIBONDS
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
!
! end point 1
!
         NI2=3*(CONI(J3)-1)
         NJ2=3*(CONJ(J3)-1)
         DO J4=1,QCIADDREP
            X1=(J4*LXYZ(NI1+1)+(QCIADDREP+1-J4)*LXYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*LXYZ(NI1+2)+(QCIADDREP+1-J4)*LXYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*LXYZ(NI1+3)+(QCIADDREP+1-J4)*LXYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*LXYZ(NI2+1)+(QCIADDREP+1-J5)*LXYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*LXYZ(NI2+2)+(QCIADDREP+1-J5)*LXYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*LXYZ(NI2+3)+(QCIADDREP+1-J5)*LXYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.DMIN) DMIN=D2
!              WRITE(*,'(A,2I6,A,4I6,A,2I6,A,F20.10)') ' intlbfgs> start constraints ',J2,J3,' atoms ', &
! &                                CONI(J2),CONJ(J2),CONI(J3),CONJ(J3),' J4,J5 ',J4,J5,' distance=',D2
           ENDDO
         ENDDO
      ENDDO
!
! end point 2
!
      NI1=3*(CONI(J2)-1)+3*NATOMS
      NJ1=3*(CONJ(J2)-1)+3*NATOMS
      DO J3=J2+1,QCIBONDS
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
!
! end point 2
!
     NI2=3*(CONI(J3)-1)
         NI2=3*(CONI(J3)-1)+3*NATOMS
         NJ2=3*(CONJ(J3)-1)+3*NATOMS
         DO J4=1,QCIADDREP
            X1=(J4*LXYZ(NI1+1)+(QCIADDREP+1-J4)*LXYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*LXYZ(NI1+2)+(QCIADDREP+1-J4)*LXYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*LXYZ(NI1+3)+(QCIADDREP+1-J4)*LXYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*LXYZ(NI2+1)+(QCIADDREP+1-J5)*LXYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*LXYZ(NI2+2)+(QCIADDREP+1-J5)*LXYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*LXYZ(NI2+3)+(QCIADDREP+1-J5)*LXYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.DMIN) DMIN=D2
!              WRITE(*,'(A,2I6,A,4I6,A,2I6,A,F20.10)') ' intlbfgs> finish constraints ',J2,J3,' atoms ', &
! &                                CONI(J2),CONJ(J2),CONI(J3),CONJ(J3),' J4,J5 ',J4,J5,' distance=',D2
           ENDDO
         ENDDO
      ENDDO
   ENDDO
   WRITE(*,'(A,F20.10,A,F20.10)') ' intlbfgs> minimum decoration distance=',DMIN,' compared with cutoff ',QCIADDREPCUT
   QCIADDREPCUT=MIN(DMIN-1.0D-3,QCIADDREPCUT)
   WRITE(*,'(A,F20.10)') ' intlbfgs> cutoff after setup is ',QCIADDREPCUT
ENDIF
!
! Check that we have a percolating constraint network. If not, increase the tolerance and try again!
! Calculate minimum number of steps of each atom from number 1 or any frozen atom.
!
NDIST1(1:NATOMS)=1000000
IF (NQCIFREEZE.EQ.0) THEN
   NDIST1(1)=0
ELSE
   DO J1=1,NATOMS
      IF (INTFROZEN(J1)) NDIST1(J1)=0
   ENDDO
ENDIF
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN1=100000
DMAX1=0
NUNCON1=0
DO J1=1,NATOMS
   IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
   DO J2=1,NCONSTRAINT
      IF (CONI(J2).EQ.J1) THEN
         IF (NDIST1(CONJ(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONJ(J2))+1
         ENDIF
      ELSE IF (CONJ(J2).EQ.J1) THEN
         IF (NDIST1(CONI(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONI(J2))+1
         ENDIF
      ENDIF
   ENDDO
   IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
   IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
   IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
ENDDO
IF (CHANGED) GOTO 5
  IF (DEBUG) WRITE(*,'(3(A,I8))') ' checkperc> steps to atom 1 converged in ',NCYCLE-1, &
    &               ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
IF (NUNCON1.GT.0) THEN
   LINTCONSTRAINTTOL=LINTCONSTRAINTTOL*1.1D0
   IF (DEBUG) WRITE(*,'(A,F15.5)') ' checkperc> increasing the local constraint tolerance parameter to ',LINTCONSTRAINTTOL
   IF (LINTCONSTRAINTTOL.GT.100.0D0) THEN
      WRITE(*,'(A,G20.10)') 'checkperc> likely ERROR *** LINTCONSTRAINTTOL=',LINTCONSTRAINTTOL
      STOP
   ENDIF
   GOTO 51
ENDIF
! IF (DEBUG) WRITE(*,'(A,F15.5)') ' checkperc> Final constraint tolerance parameter ',LINTCONSTRAINTTOL

! WRITE(*,'(A,I6,3(A,F15.5))') ' checkperc> Total distance constraints=',NCONSTRAINT, &
!   &                    ' shortest=',MINCONDIST,' longest=',MAXCONDIST,' tolerance=',LINTCONSTRAINTTOL

CALLED=.TRUE.

END SUBROUTINE CHECKPERC
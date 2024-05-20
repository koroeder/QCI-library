MODULE ADDATOM
   USE QCIPREC
   IMPLICIT NONE
   LOGICAL :: BBDONE = .FALSE.

   CONTAINS

      SUBROUTINE ADDATOM()
         USE QCIKEYS, ONLY: QCIDOBACK, NATOMS, DEBUG
         USE REPULSION, ONLY: NNREPULSIVE, NREPULSIVE
         IMPLICIT NONE
         INTEGER :: NTOADD, NADDED
         INTEGER :: NNREPSAVE, NREPSAVE
         LOGICAL :: MORETOADD
         INTEGER :: NCONTOACTIVE(1:NATOMS)
         REAL(KIND = REAL64) :: INVDTOACTIVE(1:NACTIVE)

         ! setup book keeping
         NTOADD = 1
         NADDED = 0
         ! check whether we have more backbone atoms to add
         IF (QCIDOBACK.AND.(.NOT.BBDONE)) CALL CHECK_BBLIST()

         !Save current repulsion to speed up checks later
         NNREPSAVE=NNREPULSIVE
         NREPSAVE=NREPULSIVE

         !Set variable for tracking whether we completed adding atoms
         MORETOADD = .TRUE.
         DO WHILE (MORETOADD)
            !get list of atoms by number of constraints to current active atoms
            CALL CREATE_NCONTOACTIVE_LIST(NCONTOACTIVE,INVDTOACTIVE)

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

         END DO !line 702


      END SUBROUTINE ADDATOM

      SUBROUTINE CREATE_NCONTOACTIVE_LIST(NCONTOACTIVE,INVDTOACTIVE,NBEST)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE, CONACTIVE
         USE QCICONSTRAINTS, ONLY: NCONSTRAINT, CONI, CONJ, CONDISTREF
         INTEGER, INTENT(OUT) :: NCONTOACTIVE(1:NATOMS) !list of constraints
         INTEGER, INTENT(OUT) :: NBEST ! maximum number of constraintson inactive atom
         REAL(KIND = REAL64), INTENT(OUT) :: INVDTOACTIVE(1:NACTIVE)  !inverse distance
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
            IF (NCONTOACTIVE(ATOM1).GT.BEST) NBEST = NCONTOACTIVE(ATOM1)
            IF (NCONTOACTIVE(ATOM2).GT.BEST) NBEST = NCONTOACTIVE(ATOM2)
         END DO
      END SUBROUTINE CREATE_NCONTOACTIVE_LIST

      SUBROUTINE CHECK_BBLIST()
         USE QCIKEYS, ONLY: NATOMS
         DO J1=1,NATOMS
            IF (AABACK(J1)) THEN
               IF (.NOT.ATOMACTIVE(J1)) THEN
                  ! if there are inactive BB atoms, return
                  RETURN
               END IF
            ENDIF
         ENDDO
         BBDONE = .TRUE.
         IF (DEBUG) WRITE(*,*) " check_bblist> All backbone atoms are active"
      END SUBROUTINE CHECK_BBLIST


END MODULE ADDATOM
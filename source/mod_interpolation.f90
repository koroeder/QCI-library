MODULE QCIINTERPOLATION
   USE INTERPOLATION_KEYS
   IMPLICIT NONE

   CONTAINS
      SUBROUTINE RUN_QCI_INTERPOLATION()
         USE MOD_INTCOORDS, ONLY: INITIATE_INTERPOLATION_BAND, GET_DISTANCES_CONSTRAINTS
         USE QCI_KEYS, ONLY: QCIREADGUESS
         USE CONSTRAINTS
         IMPLICIT NONE
         INTEGER :: NBEST


         !initiate vairables and band
         CALL ALLOC_INTERPOLATION_VARS()
         CALL INITIALISE_INTERPOLATION_VARS()
         CALL INITIATE_INTERPOLATION_BAND()

         ! are we restarting a simulation?
         IF (QCIRESTART) THEN
            !TODO: add option to restart interpolation
            CONTINUE
         ! alternatively, are we reading in a guess?
         ELSE IF (QCIREADGUESS) THEN
            !TODO: add option to restart interpolation
            CONTINUE
         END IF

         ! get constraint with smallest distance between endpoints (respecting QCIDOBACK)
         CALL GET_DISTANCES_CONSTRAINTS(NBEST)

         ! Turning first constraint on
         CONACTIVE(NBEST)=.TRUE.
         ATOMACTIVE(CONI(NBEST))=.TRUE.
         ATOMACTIVE(CONJ(NBEST))=.TRUE.
         IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',NBEST,' for atoms ',CONI(NBEST),CONJ(NBEST)
         IF (.NOT.QCIFROZEN(CONI(NBEST))) THEN
            TURNONORDER(NACTIVE+1)=CONI(NBEST)
            NACTIVE=NACTIVE+1
         ENDIF
         IF (.NOT.QCIFROZEN(CONJ(NBEST))) THEN
            TURNONORDER(NACTIVE+1)=CONJ(NBEST)
            NACTIVE=NACTIVE+1
         ENDIF
         NTRIES(CONI(NBEST))=1
         NTRIES(CONJ(NBEST))=1
         NREPULSIVE=0
         NCONSTRAINTON=1
         CONION(1)=CONI(NBEST)
         CONJON(1)=CONJ(NBEST)

         ! add constraints and repulsions for all frozen atoms

         ! call congrad routine

         ! now enter main loop and add atom by atom going through congrad routines as we go along

      END SUBROUTINE RUN_QCI_INTERPOLATION



      SUBROUTINE INITIALISE_INTERPOLATION_VARS()
         CONACTIVE(:) = .FALSE.
         ATOMACTIVE(1:NATOMS) = .FALSE.
         NTRIES(1:NATOMS) = 0
         NTURNONRODER(1:NATOMS) = 0
      END SUBROUTINE INITIALISE_INTERPOLATION_VARS

END MODULE QCIINTERPOLATION
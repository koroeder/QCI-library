!Module that controls the interpolation coordinates and related information such as energy and gradient
!TODO add PREV here
!TODO CONTINUE HERE:
!move constraint variables into constraint_keys module to avoid circular dependencies
MODULE MOD_INTCOORDS
   USE QCIPREC
   USE QCIKEYS, ONLY: NATOMS, NIMAGES, DEBUG
   USE INTERPOLATION_KEYS
   IMPLICIT NONE 


   INTEGER :: NCONLARGEST, NCONSMALLEST
   REAL(KIND = REAL64) :: DCONLARGEST, DCONSMALLEST
   CONTAINS

      SUBROUTINE SETUP_ENDPOINTS(NATS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATS
         NATOMS = NATS
         ALLOCATE(XSTART(3*NATOMS),XFINAL(3*NATOMS))
      END SUBROUTINE SETUP_ENDPOINTS

      SUBROUTINE INITIATE_INTERPOLATION_BAND()
         USE QCIFILEHANDLER, ONLY: GETUNIT
         USE QCIKEYS, ONLY: QCIFROZEN, NQCIFROZEN
         IMPLICIT NONE
         INTEGER :: J1
         CHARACTER(LEN=25) :: LINFILE = "int_linear.xyz"

         CALL ALLOC_INTCOORDS()
         CALL ALLOC_PREVCOORDS()
         XYZ(1:(3*NATOMS))=XSTART(1:3*NATOMS)
         XYZ((3*NATOMS)*(NIMAGES+1)+1:(3*NATOMS)*(NIMAGES+2))=XFINAL(1:3*NATOMS)
         DO J1=1,NIMAGES+2
            XYZ((J1-1)*(3*NATOMS)+1:J1*(3*NATOMS))=((NIMAGES+2-J1)*XSTART(1:3*NATOMS)+(J1-1)*XFINAL(1:3*NATOMS))/(NIMAGES+1)
         ENDDO
         XPREV(:)  = XYZ(:)

         !Deal with any frozen atoms
         IF (NQCIFROZEN.GT.0) THEN
            DO J1=1,NATOMS
               IF (QCIFROZEN(J1)) THEN
                  ATOMACTIVE(J1)=.TRUE.
                  NACTIVE=NACTIVE+1
                  TURNONORDER(NACTIVE)=J1
                  NTRIES(J1)=1
               END IF
            END DO
         END IF

         IF (DEBUG) CALL WRITE_BAND(LINFILE)

      END SUBROUTINE INITIATE_INTERPOLATION_BAND

      SUBROUTINE GET_DISTANCES_CONSTRAINTS(NBEST)
         USE QCI_CONSTRAINT_KEYS
         USE QCIKEYS, ONLY: QCIDOBACK, ISBBATOM, QCIFROZEN, QCILINEART, INLINLIST
         USE HELPER_FNCTS, ONLY: DISTANCE_ATOM_DIFF_IMAGES
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: NBEST
         INTEGER :: J1, BACKBEST
         REAL(KIND = REAL64) :: DF, D1, D2, CURRLARGEST, CURRSMALLEST, CURRBBDIST

         BACKBEST = -1
         CURRLARGEST = -1.0D100
         CURRSMALLEST = 1.0D100
         CURRBBDIST = 1.0D100

         !adjust inlinlist to account for frozen atoms, these are active by default, so we don't need them in the linear interpolation
         DO J1=1,NATOMS
            IF (QCIFROZEN(J1)) INLINLIST(J1) = .FALSE.
         END DO

         DO J1=1,NCONSTRAINT
            IF (QCILINEART.AND.(.NOT.INLINLIST(CONI(J1)).OR.(.NOT.INLINLIST(CONJ(J1))))) CYCLE            
            ! we want to collect the change in atom positions between the endpoints
            CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, XSTART, XFINAL, CONI(J1), D1)
            CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, XSTART, XFINAL, CONJ(J1), D2)
            DF = D1 + D2
            IF (DF.LT.CURRSMALLEST) THEN
               CURRSMALLEST = DF
               NCONSMALLEST = J1
            END IF
            IF (DF.LT.CURRBBDIST) THEN
               IF (QCIDOBACK.AND.(ISBBATOM(CONI(J1))).AND.(ISBBATOM(CONJ(J1)))) THEN
                  CURRBBDIST = DF
                  BACKBEST = J1
               END IF
            END IF
            IF (DF.GT.CURRLARGEST) THEN
               CURRLARGEST = DF
               NCONLARGEST = J1
            END IF
         END DO
         IF (QCIDOBACK.AND.(BACKBEST.GT.0)) THEN
            NCONSMALLEST=BACKBEST  ! ensures NBEST is set if there are frozen atoms and DOBACK is set 
         END IF
         NBEST = NCONSMALLEST
         IF (DEBUG) THEN
            WRITE(*,*) ' get_dists_constr> Smallest overall motion for constraint ',NBEST, ' atoms ', &
                       CONI(NBEST),CONJ(NBEST),' distance=', CURRSMALLEST
            WRITE(*,*) ' get_dists_constr> Largest overall motion for constraint  ',NCONLARGEST,' atoms ', &
                       CONI(NCONLARGEST),CONJ(NCONLARGEST),' distance=',CURRLARGEST
         END IF

      END SUBROUTINE GET_DISTANCES_CONSTRAINTS

      SUBROUTINE READGUESS()
         USE QCIKEYS, ONLY: GUESSFILE
         USE QCIFILEHANDLER, ONLY: GETUNIT
         IMPLICIT NONE
         INTEGER :: XUNIT
         LOGICAL :: YESNO
         CHARACTER(LEN=2) :: SDUMMY
         INTEGER :: J1, J2, NDUMMY
         CHARACTER(LEN=25) :: XYZFILE = "int_after_readguess.xyz"

         INQUIRE(FILE=GUESSFILE, EXIST=YESNO)
         IF (YESNO) THEN
            XUNIT = GETUNIT()
            OPEN(XUNIT, FILE=GUESSFILE, STATUS='OLD')
         ELSE
            WRITE(*,*) " File ", GUESSFILE , "does not exist - cannot read in guess for interpolation band - STOP"
            STOP
         END IF
         DO J2=1,NIMAGES+2
            READ(XUNIT,*) NDUMMY
            READ(XUNIT,*) 
            DO J1=1,NATOMS
               READ(XUNIT,*) SDUMMY,XYZ(3*NATOMS*(J2-1)+3*(J1-1)+1),XYZ(3*NATOMS*(J2-1)+3*(J1-1)+2),XYZ(3*NATOMS*(J2-1)+3*(J1-1)+3)
            END DO
         END DO
         IF (DEBUG) CALL WRITE_BAND(XYZFILE)
      END SUBROUTINE READGUESS

      SUBROUTINE WRITE_BAND(FNAME)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         IMPLICIT NONE
         CHARACTER(25), INTENT(IN) :: FNAME
         INTEGER :: LUNIT
         INTEGER :: J2, J3

         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE=FNAME,STATUS='replace')
         DO J2=1,NIMAGES+2
            WRITE(LUNIT,'(I6)') NATOMS
            DO J3=1,NATOMS
                  WRITE(LUNIT,'(A5,1X,3F20.10)') 'LA   ', XYZ((J2-1)*3*NATOMS+3*(J3-1)+1),  &
                               XYZ((J2-1)*3*NATOMS+3*(J3-1)+2), XYZ((J2-1)*3*NATOMS+3*(J3-1)+3)
            ENDDO
         ENDDO
         CLOSE(LUNIT)
      END SUBROUTINE WRITE_BAND

      SUBROUTINE WRITE_PROFILE(FNAME, E)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         USE QCIKEYS, ONLY: NIMAGES, DEBUG
         IMPLICIT NONE
         CHARACTER(25), INTENT(IN) :: FNAME
         REAL(KIND = REAL64), INTENT(IN) :: E(NIMAGES+2)
         INTEGER :: LUNIT
         INTEGER :: I

         LUNIT=GETUNIT()

         OPEN(UNIT=LUNIT,FILE=FNAME,STATUS='replace')
         DO I=1,NIMAGES+2
            WRITE(LUNIT,'(G24.13)') E(I)
         END DO

         CLOSE(LUNIT)

         IF (DEBUG) WRITE(*,*) " write_profile> Interpolation energy profile written to ", FNAME
      END SUBROUTINE WRITE_PROFILE
END MODULE MOD_INTCOORDS
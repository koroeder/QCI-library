MODULE QCIINTERPOLATION
   USE INTERPOLATION_KEYS
   USE MOD_INTCOORDS
   USE PREC
   IMPLICIT NONE

   CONTAINS
      SUBROUTINE RUN_QCI_INTERPOLATION()
         USE QCIKEYS, ONLY: QCIREADGUESS, QCIRESTART, QCIFREEZET, MAXITER, MAXGRADCOMP, MAXCONE, &
                            CONCUTABS, CONCUTABSINC, CHECKCHIRAL
         USE MOD_FREEZE, ONLY: ADD_CONSTR_AND_REP_FROZEN_ATOMS
         USE CONSTR_E_GRAD, ONLY: CONGRAD1, CONGRAD2, CONVERGECONTEST, CONVERGEREPTEST, &
                                  FCONMAX, FREPMAX
         USE QCIPERMDIST, ONLY: CHECK_COMMON_CONSTR
         USE QCICONSTRAINTS
         IMPLICIT NONE
         INTEGER :: NBEST, NITERDONE
         LOGICAL :: QCICONVT 
         REAL(KIND = REAL64) :: GLAST(NDIMS), XLAST(NDIMS), ELAST
         CHARACTER(25) :: XYZFILE = "int.xyz"
         CHARACTER(25) :: EEEFILE = "int.EofS"
         CHARACTER(25) :: RESETXYZFILE = "QCIreset.int.xyz"
         CHARACTER(25) :: RESETEEEFILE = "QCIreset.int.EofS"     
         REAL(KIND=REAL64), PARAMETER :: INCREASETOL = 1.1D0   


         !initiate variables for the interpolation, including image density set nimage
         CALL ALLOC_INTERPOLATION_VARS()
         CALL INITIALISE_INTERPOLATION_VARS()
         ! allocate the coordinate, energy and gradient variables for the band and
         ! initiate the interpolation band
         CALL INITIATE_INTERPOLATION_BAND()

         ! are we reading in a guess?
         IF (QCIREADGUESS) THEN
            CALL READGUESS()
         END IF

         ! get constraint with smallest distance between endpoints (respecting QCIDOBACK)
         CALL GET_DISTANCES_CONSTRAINTS(NBEST)

         ! get common constraints for atoms in permutational groups
         IF (QCIPERMT) CALL CHECK_COMMON_CONSTR()

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
         IF (QCIFREEZET) THEN
            CALL ADD_CONSTR_AND_REP_FROZEN_ATOMS(NBEST)
         END IF

         ! before we continue check repulsion neighbour list -> line 983 in old routine
         CALL CHECKREP(XYZ,0,1)

         ! call congrad routine
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
         ELSE
            CALL CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
         END IF

         !scale gradient if necessary
         IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
         !save gradient and coordinates (we use the pointer to the images here)
         GLAST(1:DIMS) = G(1:DIMS)
         XLAST(1:DIMS) = X(1:DIMS)
         ELAST = ETOTAL 

         !call check for cold fusion
         CALL CHECK_FOR_COLDFUSION(ETOTAL)

         !save concut settings
         CALL SET_CONCUTABS()

         NITERDONE = 0
         NLASTGOODE = 0
         QCICONVT = .FALSE.
         ! -> starts line 1070 in old routine
         ! now enter main loop and add atom by atom going through congrad routines as we go along
         DO WHILE (NITERDONE.LT.MAXITER)
            NITERDONE = NITERDONE + 1

            IF (QCIRESET) THEN
               !check when we had the last good energy
               IF ((NITERDONE-NLASTGOODE.GT.QCIRESETINT1)) THEN 
                  ! save interpolation if we debug
                  IF (DEBUG) THEN
                     CALL WRITE_BAND(RESETXYZFILE)
                     CALL WRITE_PROFILE(RESETEEEFILE, EEE)
                  END IF

                  IF (MAX(CONVERGECONTEST,CONVERGEREPTEST).GT.MAXCONE) MAXCONE=MAXCONE*INCREASETOL
                  IF (MAX(FCONMAX,FREPMAX).GT.QCIRMSTOL) QCIRMSTOL=QCIRMSTOL*INCREASETOL
                  CONCUTABS=CONCUTABS+0.1D0
                  WRITE(*,*) " QCI-interp> Interpolation seems to be stuck. Increasing convergence thresholds."
                  WRITE(*,*) " MAXECON = ",MAXCONE, ", QCIRMSTOL = ", QCIRMSTOL, ", CONCUTABS = ", CONCUTABS
                  CONCUTABSINC=.TRUE.
                  NCONCUTABSINC=NITERDONE
                  NLASTGOODE=NITERDONE
               ELSEIF (CONCUTABSINC) THEN
                  IF (NITERDONE-NCONCUTABSINC.GT.QCIRESETINT1) THEN ! reset CONCUTABS
                     CONCUTABS=CONCUTABSSAVE
                     !TODO: what is this line below doing?
                     !->?? IF (CCABSPHASE2) CONCUTABS=CONCUTABSSAVE2
                     CONCUTABSINC=.FALSE.
                     WRITE(*,*) " QCI-interp> Interpolation seems to be stuck. Resetting concutabs"
                     WRITE(*,*) " MAXECON = ",MAXCONE, ", QCIRMSTOL = ", QCIRMSTOL, ", CONCUTABS = ", CONCUTABS
                  ENDIF
               ENDIF
            ENDIF

            ! Checking the permutational alignment. Maintain a list of the permutable groups where all
            ! members are active. See if we have any new complete groups. MUST update NDUMMY
            ! counter to step through permutable atom list.
            IF (QCIPERMT.AND.(MOD(NITERDONE-1,QCIPERMCHECKINT).EQ.0)) THEN
               IF (CHECKCHIRAL) THEN
                  CALL CHIRALITY_CHECK()
               END IF



            END IF
            !end of permutational checks of band


            ! spring constant dynamic adjustment
            IF (QCIADJUSTKT.AND.MOD(NITERDONE,QCIADJUSTKFRQ).EQ.0) THEN
               IF (QCIAVDEV.GT.QCIADJUSTKTOL) THEN
                  KINT=MIN(KINT*QCIKADJUSTFRAC,QCIKINTMAX)
               ELSE IF (QCIAVDEV.LT.QCIADJUSTKTOL) THEN
                  KINT=MAX(KINT/QCIKADJUSTFRAC,QCIKINTMIN)
               END IF
            END IF

            ! TODO: CONOFFLIST needs to be needed

            !need to add atom now
            CALL DOATOM()
            ! HERE TO CONTINUE:
            ! need to implement the doatom routine next and then continue with the main loop

         END DO
      END SUBROUTINE RUN_QCI_INTERPOLATION

   


      SUBROUTINE INITIALISE_INTERPOLATION_VARS()
         USE QCIKEYS, ONLY: USEIMAGEDENSITY, NIMAGES, E2E_DIST, IMAGEDENSITY, MAXINTIMAGE
         IMPLICIT NONE

         CONACTIVE(:) = .FALSE.
         ATOMACTIVE(1:NATOMS) = .FALSE.
         NTRIES(1:NATOMS) = 0
         TURNONRODER(1:NATOMS) = 0

         IF (USEIMAGEDENSITY) THEN
            NIMAGES=MIN(IMAGEDENSITY*E2E_DIST,1.0D0*MAXINTIMAGE)
         END IF
      END SUBROUTINE INITIALISE_INTERPOLATION_VARS

      ! Subroutine to sclae excessive gradient components
      SUBROUTINE SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: DIMS
         REAL(KIND=REAL64), INTENT(IN) :: MAXGRADCOMP
         REAL(KIND=REAL64), INTENT(INOUT) :: G(DIMS), RMS
         REAL(KIND=REAL64):: RMSSAVE
         INTEGER :: J1

         RMSSAVE=RMS
         RMS=0.0D0
         DO J1=1,DIMS
            IF (ABS(G(J1)).GT.MAXGRADCOMP) G(J1)=SIGN(MAXGRADCOMP,G(J1)) ! should be MAXGRADCOMP with the sign of G(J1)
            RMS=RMS+G(J1)**2
         ENDDO
         RMS=SQRT(RMS/REAL(DIMS))
         WRITE(*,'(A,3G20.10)') 'scalegrad> RMS, RMSSAVE, MAXGRADCOMP=',RMS,RMSSAVE,MAXGRADCOMP
      END SUBROUTINE SCALEGRAD

      SUBROUTINE CHECK_FOR_COLDFUSION(ECURRENT)
         USE MOD_TERMINATE, ONLY: INT_ERR_TERMINATE
         USE QCIKEYS, ONLY: NIMAGES, DEBUG, COLDFUSIONLIMIT
         USE MOD_INTCOORDS, ONLY: WRITE_BAND
         IMPLICIT NONE
         CHARACTER(25) :: FNAME="intcoords.aborted.xyz"
         REAL(KIND=REAL64) :: EPERIMAGE

         EPERIMAGE  = ECURRENT/REAL(NIMAGES)
         IF (EPERIMAGE.LT.COLDFUSIONLIMIT) THEN
            WRITE(*,*) " check_int> Cold fusion diagnosed - aborting interpolation"
            WRITE(*,*) "            E per image: ", EPERIMAGE, " Limit: ", COLDFUSIONLIMIT
            WRITE(*,*) "            Current interpolation band written to ", FNAME
            CALL WRITE_BAND(FNAME)
            CALL INT_ERR_TERMINATE()
         ELSE
            IF (DEBUG) THEN
               WRITE(*,*) " check_int> No cold fusion diagnosed, per image energy: ", EPERIMAGE
            END IF
         END IF

      END SUBROUTINE CHECK_FOR_COLDFUSION
END MODULE QCIINTERPOLATION
MODULE QCIINTERPOLATION
   USE INTERPOLATION_KEYS
   USE MOD_INTCOORDS
   USE PREC
   IMPLICIT NONE
   REAL(KIND = REAL64), ALLOCATABLE :: GTMP(:), DIAG(:), STP(:), SEARCHSTEP(:,:), GDIF(:,:), STEPIMAGE(:)
   
   CONTAINS
      SUBROUTINE RUN_QCI_INTERPOLATION()
         USE QCIKEYS, ONLY: QCIREADGUESS, QCIRESTART, QCIFREEZET, MAXITER, MAXGRADCOMP, MAXCONE, &
                            CONCUTABS, CONCUTABSINC, CHECKCHIRAL
         USE MOD_FREEZE, ONLY: ADD_CONSTR_AND_REP_FROZEN_ATOMS
         USE CONSTR_E_GRAD, ONLY: CONGRAD1, CONGRAD2, CONVERGECONTEST, CONVERGEREPTEST, &
                                  FCONMAX, FREPMAX
         USE QCIPERMDIST, ONLY: CHECK_COMMON_CONSTR, UPDATE_ACTIVE_PERMGROUPS
         USE QCICONSTRAINTS
         USE CHIRALITY, ONLY: ASSIGNMENT_SR
         IMPLICIT NONE
         INTEGER :: NBEST, NITERDONE, FIRSTATOM
         INTEGER :: J1
         LOGICAL :: QCICONVT 
         REAL(KIND = REAL64) :: GLAST(NDIMS), XLAST(NDIMS), ELAST
         CHARACTER(25) :: XYZFILE = "int.xyz"
         CHARACTER(25) :: EEEFILE = "int.EofS"
         CHARACTER(25) :: RESETXYZFILE = "QCIreset.int.xyz"
         CHARACTER(25) :: RESETEEEFILE = "QCIreset.int.EofS"     
         REAL(KIND = REAL64), PARAMETER :: INCREASETOL = 1.1D0   
         INTEGER :: NITERUSE, POINT, NPT ! variables needed for lbfgs steps
         REAL(KIND = REAL64) :: STPMIN ! minimum step-size
         REAL(KIND = REAL64) :: CURRMAXSEP, CURRMINSEP ! current minimum and maximum image separation
         INTEGER :: IDXMIN, IDXMAX ! id for minimum and maximum distance images

         !initiate variables for the interpolation, including image density set nimage
         CALL ALLOC_INTERPOLATION_VARS()
         CALL INITIALISE_INTERPOLATION_VARS()
         ! allocate the coordinate, energy and gradient variables for the band and
         ! initiate the interpolation band
         CALL INITIATE_INTERPOLATION_BAND()

         !initialise step taking variables

         ! are we reading in a guess?
         IF (QCIREADGUESS) THEN
            CALL READGUESS()
         END IF

         ! get constraint with smallest distance between endpoints (respecting QCILINEAR and QCIDOBACK)
         CALL GET_DISTANCES_CONSTRAINTS(NBEST)

         ! get common constraints for atoms in permutational groups
         IF (QCIPERMT) CALL CHECK_COMMON_CONSTR()

         ! Turning first constraint on
         CONACTIVE(NBEST)=.TRUE.
         ATOMACTIVE(CONI(NBEST))=.TRUE.
         ATOMACTIVE(CONJ(NBEST))=.TRUE.
         IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' QCIinterp> Turning on constraint ',NBEST,' for atoms ',CONI(NBEST),CONJ(NBEST)
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

         ! before we continue check repulsion neighbour list
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
         NITERUSE = 1
         NLASTGOODE = 0
         QCICONVT = .FALSE.

         ! now enter main loop and add atom by atom going through congrad routines as we go along
         DO WHILE (NITERDONE.LT.MAXITER)
            NITERDONE = NITERDONE + 1
            WRITE(*,*)
            WRITE(*,*) ">>> QCI interpolation cycle ", NITERDONE
            WRITE(*,*) "    Number of active atoms: ", NACTIVE, "   Number of active constraints: ", NCONSTRAINTON
            WRITE(*,*)
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
                  WRITE(*,*) " QCIinterp> Interpolation seems to be stuck. Increasing convergence thresholds."
                  WRITE(*,*) "            MAXECON = ",MAXCONE, ", QCIRMSTOL = ", QCIRMSTOL, ", CONCUTABS = ", CONCUTABS
                  CONCUTABSINC=.TRUE.
                  NCONCUTABSINC=NITERDONE
                  NLASTGOODE=NITERDONE
               ELSEIF (CONCUTABSINC) THEN
                  IF (NITERDONE-NCONCUTABSINC.GT.QCIRESETINT1) THEN ! reset CONCUTABS
                     CONCUTABS=CONCUTABSSAVE
                     !TODO: what is this line below doing?
                     !->?? IF (CCABSPHASE2) CONCUTABS=CONCUTABSSAVE2
                     CONCUTABSINC=.FALSE.
                     WRITE(*,*) " QCIinterp> Interpolation seems to be stuck. Resetting concutabs"
                     WRITE(*,*) "            MAXECON = ",MAXCONE, ", QCIRMSTOL = ", QCIRMSTOL, ", CONCUTABS = ", CONCUTABS
                  ENDIF
               ENDIF
            ENDIF

            ! Checking the permutational alignment. Maintain a list of the permutable groups where all
            ! members are active. See if we have any new complete groups. MUST update NDUMMY
            ! counter to step through permutable atom list.
            IF (QCIPERMT.AND.(MOD(NITERDONE-1,QCIPERMCHECKINT).EQ.0)) THEN
               IF (CHECKCHIRAL) THEN
                  WRITE(*,*) " QCIinterp> Checking chirality across band"
                  CALL CHIRALITY_CHECK() !line 1152
                  !TODO: need to complete this subroutine - discussion how we replace/fix issues with chirality
               END IF
               !update active permutational groups
               WRITE(*,*) " QCIinterp> Updating active permutational groups"
               CALL UPDATE_ACTIVE_PERMGROUPS()

               ! Checking all active groups across the band - do we have the best alignement?
               FIRSTATOM = 1
               WRITE(*,*) " QCIinterp> Checking permutational alignment across the band"
               DO J1=1,NPERMGROUP
                  IF (GROUPACTIVE(J1)) THEN
                     !check permutational consistency forward
                     CALL CHECK_PERM_BAND(J1, FIRSTATOM, .FALSE.)
                     !check permutational consistency in reverse
                     CALL CHECK_PERM_BAND(J1, FIRSTATOM, .TRUE.)
                  END IF
                  FIRSTATOM = FIRSTATOM + NPERMSIZE(J1)
               END DO
            END IF
            !end of permutational checks of band

            ! spring constant dynamic adjustment
            IF (QCIADJUSTKT.AND.MOD(NITERDONE,QCIADJUSTKFRQ).EQ.0) THEN
               IF (QCIAVDEV.GT.QCIADJUSTKTOL) THEN
                  WRITE(*,*) " QCIinterp> Lowering spring constant from ", KINT, " to ", MIN(KINT*QCIKADJUSTFRAC,QCIKINTMAX)
                  KINT=MIN(KINT*QCIKADJUSTFRAC,QCIKINTMAX)
               ELSE IF (QCIAVDEV.LT.QCIADJUSTKTOL) THEN
                  KINT=MAX(KINT/QCIKADJUSTFRAC,QCIKINTMIN)
                  WRITE(*,*) " QCIinterp> Increasing spring constant from ", KINT, " to ", MAX(KINT*QCIKADJUSTFRAC,QCIKINTMAX)
               END IF
            END IF

            !if not all atoms are active, we add an atom now to the active set and hence the images
            IF (ADDATOM.AND.(NACTIVE.LT.NATOMS)) THEN
               WRITE(*,*) " QCIinterp> Adding the next atom to active set" 
               CALL ADDATOM()
               !scale gradient if necessary
               IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
               NLASTGOODE=NITERDONE
            END IF

            GTMP(1:D)=0.0D0
            ! the variables needed for step taking are either module variables in this module or saved in mod_intcoords
            CALL MAKESTEP(NITERUSE,NPT,POINT)

            IF ((DOT_PRODUCT(G,GTMP)/MAX(1.0-100,SQRT(DOT_PRODUCT(G,G))*SQRT(DOT_PRODUCT(GTMP,GTMP)))).GT.0.0D0) THEN
               IF (DEBUG) WRITE(*,*) 'Search direction has positive projection onto gradient - reversing step'
               GTMP(1:DIMS)=-GTMP(1:DIMS)
               SEARCHSTEP(POINT,1:DIMS)=GTMP(1:DIMS)
            END IF
            GTMP(1:D)=G(1:D)
            ! Take the minimum scale factor for all images for LBFGS step to avoid discontinuities
            STPMIN = 1.0D0
            DO J1=1,NIMAGES
               STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,(3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1), &
                                                SEARCHSTEP(POINT,(3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1)))
               IF (STEPIMAGE(J1).GT.MAXQCIBFGS) THEN
                  STP((3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1) = MAXQCIBFGS/STEPIMAGE(J1)
                  STPMIN=MIN(STPMIN,STP((3*NATOMS)*(J1-1)+1))
               END IF
            END DO
            STP(1:D) = STPMIN

            !continue at line 1575 with the removing and adding image blocks
            !check for image addition or removal
            IF (REMOVEIMAGE.OR.(MOD(NITERDONE,QCIIMAGECHECK).EQ.0)) THEN
               MOREIMAGES=.TRUE.
               DO WHILE(MOREIMAGES)
                  CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
                  IF ((.NOT.REMOVEIMAGE).AND.((CURRMAXSEP.GT.IMSEPMAX).AND.(NIMAGES.LT.MAXNIMAGES))) THEN
                     !TODO CONTINUE LINE 1632
                     CALL ADD_IMAGE()
                  END IF
                  IF (REMOVEIMAGE.OR.((CURRMINSEP.LT.IMSEPMIN).AND.(INTIMAGE.GT.1))) THEN
                     CALL REMOVE_IMAGE()
                  END IF
               END DO
            END IF

         END DO
      END SUBROUTINE RUN_QCI_INTERPOLATION

      SUBROUTINE GET_IMAGE_SEPARATION(DMIN,DMAX,JMIN,JMAX)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         USE HELPER_FNCTS, ONLY: DISTANCE_ATOM_DIFF_IMAGES
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(OUT) :: DMIN, DMAX
         INTEGER, INTENT(OUT) :: JMIN, JMAX
         REAL(KIND = REAL64) :: DISTATOM, DISTTOTAL, X1(3*NATOMS), X2(3*NATOMS)
         INTEGER :: J1, J2

         DMIN = HUGE(1.0D0)
         DMAX = -1.0D0

         DO J1=1,NIMAGES+1
            DISTTOTAL = 0.0D0
            X1(1:3*NATOMS) = XYZ((J1-1)*3*NATOMS+1:J1*3*NATOMS)
            X2(1:3*NATOMS) = XYZ(J1*3*NATOMS+1:(J1+1)*3*NATOMS)
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J2)) THEN
                  CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, X1, X2, IDX, DISTATOM)
                  DISTTOTAL = DISTTOTAL + DISTATOM
               END IF
            END DO
            IF (DISTTOTAL.GT.DMAX) THEN
               DMAX = DISTTOTAL
               JMAX = J1
            END IF
            IF (DISTTOTAL.LT.DMIN) THEN
               DMIN = DISTTOTAL
               JMIN = J1
            END IF
         END DO
      END SUBROUTINE GET_IMAGE_SEPARATION


      SUBROUTINE ALLOC_STEPTAKING()
         USE QCIKEYS, ONLY: NIMAGES, NATOMS,MUPDATE
         IMPLICIT NONE
         CALL DEALLOC_STEPTAKING()
         ALLOCATE(DIAG(3*NATOMS*NIMAGES))
         ALLOCATE(GTMP(3*NATOMS*NIMAGES))
         ALLOCATE(GDIF(0:MUPDATE,(3*NATOMS)*NIMAGES))
         ALLOCATE(STP(3*NATOMS*NIMAGES))
         ALLOCATE(SEARCHSTEP(0:MUPDATE,(3*NATOMS)*NIMAGES))
         ALLCOATE(STEPIMAGE(NIMAGES))
      END SUBROUTINE ALLOC_STEPTAKING

      SUBROUTINE DEALLOC_STEPTAKING()
         IF (ALLOCATED(DIAG)) DEALLOCATE(DIAG)
         IF (ALLOCATED(GTMP)) DEALLOCATE(GTMP)
         IF (ALLOCATED(GDIF)) DEALLOCATE(GDIF)
         IF (ALLOCATED(STP)) DEALLOCATE(STP)
         IF (ALLOCATED(SEARCHSTEP)) DEALLOCATE(SEARCHSTEP)
         IF (ALLOCATED(STEPIMAGE)) DEALLOCATE(STEPIMAGE)
      END SUBROUTINE DEALLOC_STEPTAKING


      SUBROUTINE INITIALISE_INTERPOLATION_VARS()
         USE QCIKEYS, ONLY: USEIMAGEDENSITY, NIMAGES, E2E_DIST, IMAGEDENSITY, MAXINTIMAGE
         USE ADDATOM, ONLY: ALLOC_ADDATOM
         IMPLICIT NONE

         CONACTIVE(:) = .FALSE.
         ATOMACTIVE(1:NATOMS) = .FALSE.
         NTRIES(1:NATOMS) = 0
         TURNONRODER(1:NATOMS) = 0

         IF (USEIMAGEDENSITY) THEN
            NIMAGES=MIN(IMAGEDENSITY*E2E_DIST,1.0D0*MAXINTIMAGE)
         END IF

         CALL ALLOC_ADDATOM()
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
         REAL(KIND=REAL64), INTENT(IN) :: ECURRENT
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

      SUBROUTINE MAKESTEP(NITERDONE,NPT,POINT)
         USE QCIKEYS, ONLY: MUPDATE, DGUESS, NATOMS, NIMAGES
         USE MOD_INTCOORDS, ONLY: G, DIMS
         IMPLICIT NONE         
         INTEGER, INTENT(IN) :: NITERDONE
         INTEGER, INTENT(IN) :: NPT
         INTEGER, INTENT(OUT) :: POINT
         REAL(KIND = REAL64) :: GNORM
         REAL(KIND = REAL64) :: YS, YY, YR, SQ, BETA
         REAL(KIND = REAL64) :: RHO1(MUPDATE), ALPHA(MUPDATE)
         INTEGER :: BOUND, CP, I

         !if it is the first step, we use a cautious guess
         IF (NITERDONE.EQ.1) THEN
            POINT = 0
            DIAG(1:NDIMS) = DGUESS
            SEARCHSTEP(0,1:NDIMS) = -DGUESS*G(1:NDIMS)
            GTMP(1:DIMS) = SEARCHSTEP(0,1:NDIMS)
            GNORM =  MAX(SQRT(DOT_PRODUCT(G(1:DIMS),G(1:DIMS))),1.0D-100)
            STP(1:DIMS) = MIN(1.0D0/GNORM, GNORM)
            RETURN
         END IF

         IF (NITERDONE.GT.MUPDATE) THEN
            BOUND = MUPDATE
         ELSE
            BOUND = NITERDONE - 1
         END IF
         YS=DOT_PRODUCT(GDIF(NPT/DIMS,:),SEARCHSTEP(NPT/DIMS,:))
         IF (YS==0.0D0) YS=1.0D0

         ! Update estimate of diagonal inverse Hessian elements.
         YY=DOT_PRODUCT(GDIF(NPT/DIMS,:),GDIF(NPT/DIMS,:))
         IF (YY==0.0D0) YY=1.0D0
         DIAG(1) = YS/YY

         ! Compute -H*G using the formula given in: Nocedal, J. 1980, Mathematics of Computation, Vol.35, No.151, pp. 773-782
         IF (POINT.EQ.0) THEN
            CP = MUPDATE
         ELSE
            CP = POINT
         END IF

         RHO1(CP) = 1.0D0/YS
         GTMP(1:DIMS) = -G(1:DIMS)
         CP = POINT

         DO I=1,BOUND
            CP = CP - 1
            IF (CP.EQ.-1) CP = MUPDATE-1
            SQ = DOT_PRODUCT(SEARCHSTEP(CP,1:DIMS),GTMP(1:DIMS))
            ALPHA(CP+1) = RHO1(CP+1) * SQ
            GTMP(1:DIMS) = -ALPHA(CP+1)*GDIF(CP,1:DIMS) + GTMP(1:DIMS)
         END DO

         GTMP(1:DIMS)=DIAG(1)*GTMP(1:DIMS)

         DO I=1,BOUND
            YR = DOT_PRODUCT(GDIF(CP,1:DIMS),GTMP)
            BETA= RHO1(CP+1)*YR
            BETA= ALPHA(CP+1)-BETA
            GTMP(1:DIMS) = BETA*SEARCHSTEP(CP,1:DIMS) + GTMP(1:DIMS)
            CP=CP+1
            IF (CP.EQ.MUPDATE) CP=0
         END DO

         STP(1:DIMS) = 1.0D0
         SEARCHSTEP(POINT,1:DIMS)=GTMP(1:DIMS)
      END SUBROUTINE MAKESTEP

END MODULE QCIINTERPOLATION
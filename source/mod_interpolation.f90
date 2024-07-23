MODULE QCIINTERPOLATION
   USE INTERPOLATION_KEYS
   USE MOD_INTCOORDS
   USE QCIPREC
   IMPLICIT NONE
   LOGICAL :: QCICOMPLETE

   CONTAINS
      SUBROUTINE RUN_QCI_INTERPOLATION()
         USE QCIKEYS, ONLY: QCIREADGUESS, QCIFREEZET, MAXITER, MAXGRADCOMP, MAXCONE, &
                            CHECKCHIRAL, CHECKCONINT, CHECKREPINTERVAL, IMSEPMIN, IMSEPMAX, &
                            KINT, MAXERISE, MAXINTIMAGE, MAXQCIBFGS, MUPDATE, DGUESS, &
                            DUMPQCIXYZFRQS, DUMPQCIXYZ, QCIADJUSTKFRQ, QCIADJUSTKT, QCIAVDEV, &
                            QCIKINTMIN, QCIKINTMAX, QCIADJUSTKFRAC, QCIADJUSTKTOL, QCIIMAGECHECK, &
                            QCIPERMCHECKINT, QCIPERMT, QCIRMSTOL, QCIFROZEN, QCIRESET, QCIRESETINT1
         USE MOD_FREEZE, ONLY: ADD_CONSTR_AND_REP_FROZEN_ATOMS
         USE CONSTR_E_GRAD, ONLY: CONGRAD1, CONGRAD2, CONVERGECONTEST, CONVERGEREPTEST, &
                                  FCONMAX, FREPMAX
         USE QCIPERMDIST, ONLY: NPERMGROUP, CHECK_COMMON_CONSTR, UPDATE_ACTIVE_PERMGROUPS, GROUPACTIVE, &
                                NPERMSIZE, CHECK_PERM_BAND
         USE QCI_CONSTRAINT_KEYS
         USE CHIRALITY, ONLY: ASSIGNMENT_SR, CHIRALITY_CHECK
         USE ADDINGATOM, ONLY: ADDATOM
         USE REPULSION, ONLY: NREPULSIVE, CHECKREP
         USE ADDREMOVE_IMAGES, ONLY: ADD_IMAGE, REMOVE_IMAGE
         IMPLICIT NONE
         INTEGER :: NBEST, NITERDONE, FIRSTATOM, NCONCUTABSINC, NDECREASE, NFAIL, NLASTGOODE
         INTEGER :: J1
         LOGICAL :: QCICONVT 
         LOGICAL :: ADDATOMT
         LOGICAL :: ACCEPTEDSTEP
         LOGICAL :: MOREIMAGES
         CHARACTER(25) :: XYZFILE = "int.xyz"
         CHARACTER(25) :: EEEFILE = "int.EofS"
         CHARACTER(25) :: RESETXYZFILE = "QCIreset.int.xyz"
         CHARACTER(25) :: RESETEEEFILE = "QCIreset.int.EofS"     
         CHARACTER(6) :: ITERSTRING
         REAL(KIND = REAL64) :: ETOTAL, EPREV
         REAL(KIND = REAL64) :: CONCUTABSSAVE
         INTEGER :: EXITSTATUS
         REAL(KIND = REAL64), PARAMETER :: INCREASETOL = 1.1D0   
         REAL(KIND = REAL64), PARAMETER :: STPREDUCTION = 10.0D0
         INTEGER :: NITERUSE, POINT, NPT ! variables needed for lbfgs steps
         REAL(KIND = REAL64) :: RHO1(MUPDATE), ALPHA(MUPDATE)
         REAL(KIND = REAL64) :: STPMIN ! minimum step-size
         REAL(KIND = REAL64) :: CURRMAXSEP, CURRMINSEP ! current minimum and maximum image separation
         INTEGER :: IDXMIN, IDXMAX ! id for minimum and maximum distance images

         LOGICAL :: CHIRCHECK !debugging variable


         CHIRCHECK = .FALSE.

         !initiate variables for the interpolation, including image density set nimage
         CALL ALLOC_INTERPOLATION_VARS()
         CALL ALLOC_STEPTAKING()
         !RHO1 = 1.0D0
         !ALPHA = 1.0D0
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
         GPREV(:) = GGG(:)
         XPREV(:) = XYZ(:)
         EPREV = ETOTAL 

         !call check for cold fusion
         CALL CHECK_FOR_COLDFUSION(ETOTAL)

         !save concut settings
         CONCUTABSSAVE = CONCUTABS

         NITERDONE = 0
         NITERUSE = 1
         NPT = 0
         NLASTGOODE = 0
         QCICONVT = .FALSE.
         ADDATOMT = .TRUE.
         NCONCUTABSINC = 0
         CONCUTABSINC=.FALSE.

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
                     CONCUTABSINC=.FALSE.
                     WRITE(*,*) " QCIinterp> Interpolation seems to be stuck. Resetting concutabs"
                     WRITE(*,*) "            MAXECON = ",MAXCONE, ", QCIRMSTOL = ", QCIRMSTOL, ", CONCUTABS = ", CONCUTABS
                  ENDIF
               ENDIF
            ENDIF

            ! Checking the permutational alignment. Maintain a list of the permutable groups where all
            ! members are active. See if we have any new complete groups. MUST update NDUMMY
            ! counter to step through permutable atom list.
            IF (QCIPERMT.AND.(MOD(NITERDONE,QCIPERMCHECKINT).EQ.0)) THEN
               IF (CHECKCHIRAL) THEN
                  WRITE(*,*) " QCIinterp> Checking chirality across band"
                  CALL CHIRALITY_CHECK(XYZ) 
               END IF
               WRITE(*,*) "check images after chirality check"
               CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)

               !update active permutational groups
               WRITE(*,*) " QCIinterp> Updating active permutational groups"
               CALL UPDATE_ACTIVE_PERMGROUPS()
               WRITE(*,*) "check images after updating permutational groups"
               CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)

               ! Checking all active groups across the band - do we have the best alignement?
               FIRSTATOM = 1
               WRITE(*,*) " QCIinterp> Checking permutational alignment across the band, NPERMGROUPS: ", NPERMGROUP
               DO J1=1,NPERMGROUP
                  IF (GROUPACTIVE(J1)) THEN
                     WRITE(*,*) "check for permutational group: ", J1
                     !check permutational consistency forward
                     CALL CHECK_PERM_BAND(J1, FIRSTATOM, .FALSE.)
                     WRITE(*,*) "completed forwards pass"
                     !check permutational consistency in reverse
                     CALL CHECK_PERM_BAND(J1, FIRSTATOM, .TRUE.)
                     WRITE(*,*) "completed backwards pass"
                  END IF
                  FIRSTATOM = FIRSTATOM + NPERMSIZE(J1)
               END DO
               WRITE(*,*) "check images after updating permutational check"
               CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
               CHIRCHECK=.TRUE.
            END IF
            !end of permutational checks of band

            ! spring constant dynamic adjustment
            IF (QCIADJUSTKT.AND.MOD(NITERDONE,QCIADJUSTKFRQ).EQ.0) THEN
               IF (QCIAVDEV.GT.QCIADJUSTKTOL) THEN
                  WRITE(*,*) " QCIinterp> Lowering spring constant from ", KINT, " to ", MAX(KINT/QCIADJUSTKFRAC,QCIKINTMIN)
                  KINT=MAX(KINT/QCIADJUSTKFRAC,QCIKINTMIN)
               ELSE IF (QCIAVDEV.LT.QCIADJUSTKTOL) THEN
                  KINT=MIN(KINT*QCIADJUSTKFRAC,QCIKINTMAX)
                  WRITE(*,*) " QCIinterp> Increasing spring constant from ", KINT, " to ", MIN(KINT*QCIADJUSTKFRAC,QCIKINTMAX)
               END IF
            END IF

            !if not all atoms are active, we add an atom now to the active set and hence the images
            IF (ADDATOMT.AND.(NACTIVE.LT.NATOMS)) THEN
               WRITE(*,*) " QCIinterp> Adding the next atom to active set" 
               CALL ADDATOM()
               !scale gradient if necessary
               IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
               NLASTGOODE=NITERDONE
               IF (CHIRCHECK) THEN
                  WRITE(*,*) "check images after adding atom"
                  CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
               END IF
            END IF

            GTMP(1:DIMS)=0.0D0
            ! the variables needed for step taking are either module variables in this module or saved in mod_intcoords
            CALL MAKESTEP(NITERUSE,NPT,POINT,RHO1,ALPHA)

            IF (CHIRCHECK) THEN
               WRITE(*,*) "check images after make step"
               CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
            END IF

            IF ((DOT_PRODUCT(G,GTMP)/MAX(1.0-100,SQRT(DOT_PRODUCT(G,G))*SQRT(DOT_PRODUCT(GTMP,GTMP)))).GT.0.0D0) THEN
               IF (DEBUG) WRITE(*,*) ' QCIinterp - Search direction has positive projection onto gradient - reversing step'
               GTMP(1:DIMS)=-GTMP(1:DIMS)
               SEARCHSTEP(POINT,1:DIMS)=GTMP(1:DIMS)
            END IF

            GTMP(1:DIMS)=G(1:DIMS)
            ! Take the minimum scale factor for all images for LBFGS step to avoid discontinuities
            STPMIN = 1.0D0
            DO J1=1,NIMAGES
               STEPIMAGE(J1) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,(3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1), &
                                                SEARCHSTEP(POINT,(3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1)))
               IF (STEPIMAGE(J1).GT.MAXQCIBFGS) THEN
                  STP((3*NATOMS)*(J1-1)+1:(3*NATOMS)*J1) = MAXQCIBFGS/STEPIMAGE(J1)
                  STPMIN=MIN(STPMIN,STP((3*NATOMS)*(J1-1)+1))
               END IF
            END DO
            STP(1:DIMS) = STPMIN

            ACCEPTEDSTEP = .FALSE.
            NDECREASE = 0

            IF (CHIRCHECK) THEN
               WRITE(*,*) "check images after applying search step"
               CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
            END IF

            DO WHILE(.NOT.ACCEPTEDSTEP)


               !apply our step to our coordinates
               X(1:DIMS) = X(1:DIMS) + STP(1:DIMS)*SEARCHSTEP(POINT,1:DIMS)

               IF (CHIRCHECK) THEN
                  WRITE(*,*) "STP: ", STP
                  WRITE(*,*) "SEARCHSTEP: ", SEARCHSTEP
                  WRITE(*,*) "X: ", X
                  WRITE(*,*) "check images after updating coordinates"
                  CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
                  STOP
               END IF

               !adding or removing images if required
               IF (MOD(NITERDONE,QCIIMAGECHECK).EQ.0) THEN
                  MOREIMAGES=.TRUE.
                  DO WHILE(MOREIMAGES)
                     MOREIMAGES = .FALSE.
                     CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
                     IF ((CURRMAXSEP.GT.IMSEPMAX).AND.(NIMAGES.LT.MAXINTIMAGE)) THEN
                        WRITE(*,*) " QCIinterp> Adding image between images ", IDXMAX, " and ", IDXMAX+1
                        CALL ADD_IMAGE(IDXMAX,ETOTAL,RMS)
                        NITERUSE = 0
                        !scale gradient if necessary
                        IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
                        NLASTGOODE=NITERDONE
                        MOREIMAGES = .TRUE. !if we add an atom we will stay in the loop and check whether we should add more images
                     END IF
                     IF ((CURRMINSEP.LT.IMSEPMIN).AND.(NIMAGES.GT.1)) THEN
                        IF (IDXMIN.EQ.1) IDXMIN = 2
                        WRITE(*,*) " QCIinterp> Removing image ", IDXMIN
                        CALL REMOVE_IMAGE(IDXMAX,ETOTAL,RMS)
                        NITERUSE = 0
                        !scale gradient if necessary
                        IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
                        NLASTGOODE=NITERDONE
                        MOREIMAGES = .TRUE. !if we add an atom we will stay in the loop and check whether we should add more images
                     END IF
                     !line 1798 gets us to here
                     IF (.NOT.MOREIMAGES) WRITE(*,*) " QCIinterp> Not adding or removing further images"
                  END DO
               ELSE
                  CALL GET_IMAGE_SEPARATION(CURRMINSEP,CURRMAXSEP,IDXMIN,IDXMAX)
               END IF

               !get energies and gradient and check whether we are making progress
               
               !at some intervals check repulsion neighbour list
               IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREP(XYZ,0,1)

               ! call congrad routine
               IF (CHECKCONINT) THEN
                  CALL CONGRAD2(ETOTAL, XYZ, GGG, EEE, RMS)
               ELSE
                  CALL  CONGRAD1(ETOTAL, XYZ, GGG, EEE, RMS)
               END IF
               !TODO: check with old routine: There is a weirder energy call here that I am not sure I udnerstand:
               !CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)

               IF ((ETOTAL-EPREV.LT.MAXERISE).OR.ADDATOMT) THEN
                  EPREV = ETOTAL
                  GPREV(:) = GGG(:)
                  XPREV(:) = XYZ(:)
                  ACCEPTEDSTEP = .TRUE.
               ELSE
                  NDECREASE = NDECREASE + 1
                  !TODO: add NDECREASE and a parameter variable for its limit 
                  IF (NDECREASE.GT.5) THEN
                     NFAIL = NFAIL + 1
                     XYZ(:) = XPREV(:)
                     GGG(:) = GPREV(:)
                     WRITE(*,*) " QCIinterp> WARNING - LBFGS cannot find a lower energy, NFAIL=",NFAIL
                     ACCEPTEDSTEP = .TRUE. !we failed to many times, so for now we accept failure and leave the loop
                  ELSE
                     XYZ(:) = XPREV(:)
                     GGG(:) = GPREV(:)
                     !TODO: add stepreduction as parameter variable at 10.0D0
                     STP(1:DIMS) = STP(1:DIMS)/STPREDUCTION  
                     WRITE(*,*) " QCIinterp> Energy increased from ", EPREV, "to", ETOTAL, "; decreasing step size"                 
                  END IF
               END IF
            END DO ! end of loop for accepting the step size
            !scale gradient if necessary
            IF (MAXGRADCOMP.GT.0.0D0) CALL SCALEGRAD(DIMS,G,RMS,MAXGRADCOMP)
            CALL CHECK_FOR_COLDFUSION(ETOTAL)
            CALL GET_STATISTIC_INTERP()

            ! set DGUESS to a reasonable guess
            DGUESS=DIAG(1)
            !get exit status
            CALL SET_EXIT_STATUS(NITERDONE,EXITSTATUS)
            !first check if we should add an atom
            ADDATOMT = .FALSE.
            IF ((EXITSTATUS.GT.0).AND.(NACTIVE.LT.NATOMS)) THEN
               ADDATOMT = .TRUE. 
            END IF
            !if we don't add an atom, check if we converged
            IF (.NOT.(ADDATOMT).AND.(EXITSTATUS.EQ.1)) THEN
               !we have converged and there is no more atom to add - we can leave the main loop
               EXIT
            END IF

            ! Compute the new step and gradient change
            NPT=POINT*DIMS
            SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
            GDIF(POINT,:)=G-GTMP

            POINT=POINT+1
            IF (POINT.EQ.MUPDATE) POINT=0
            NITERUSE=NITERUSE+1

            IF (DUMPQCIXYZ.AND.(MOD(NITERDONE,DUMPQCIXYZFRQS).EQ.0)) THEN
               WRITE(ITERSTRING,'(I6)') NITERDONE 
               CALL WRITE_BAND(XYZFILE//"."//ADJUSTL(TRIM(ITERSTRING)))
               CALL WRITE_PROFILE(EEEFILE//"."//ADJUSTL(TRIM(ITERSTRING)),EEE)
            END IF

         END DO
         ! if the exit status is not 1, we left without convergence 
         IF (NACTIVE.LT.NATOMS) THEN
            WRITE(*,*) " QCIinterp> QCI did not activate all atoms during interpolation. Final number of active atoms: ", NACTIVE
            CALL INT_ERR_TERMINATE()
         END IF
         IF (EXITSTATUS.EQ.1) THEN
            WRITE(*,*) " QCIinterp> Converged after ", NITERDONE," steps, energy/image=",ETOTAL/NIMAGES, &
                                 ' RMS=',RMS,' images=',NIMAGES
            QCICOMPLETE = .TRUE.
         ELSE
            WRITE(*,*) " QCIinterp> Not converged after ", NITERDONE," steps, energy/image=",ETOTAL/NIMAGES, &
                                 ' RMS=',RMS,' images=',NIMAGES, ", but all atoms were activated."
            QCICOMPLETE = .FALSE.                     
         END IF
         CALL WRITE_BAND(XYZFILE)
         CALL WRITE_PROFILE(EEEFILE,EEE)
         WRITE(*,*) " QCIinterp> Leaving interpolation"
      END SUBROUTINE RUN_QCI_INTERPOLATION


      SUBROUTINE SET_EXIT_STATUS(NITERDONE,EXITSTATUS)
         USE CONSTR_E_GRAD, ONLY: CONVERGECONTEST, CONVERGEREPTEST, FCONTEST, FREPTEST
         USE QCIKEYS, ONLY: MAXCONE, QCIRMSTOL, INTADDATOM
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NITERDONE
         INTEGER, INTENT(OUT) :: EXITSTATUS

         EXITSTATUS = 0
         !is the simulation converged?
         IF ((FCONTEST.LT.QCIRMSTOL).AND.(FREPTEST.LT.QCIRMSTOL).AND.(CONVERGECONTEST.LT.MAXCONE).AND.(CONVERGEREPTEST.LT.MAXCONE).AND.(NITERDONE.GT.1)) THEN
            EXITSTATUS = 1
         END IF

         !should we add more atoms?
         IF (MOD(NITERDONE,INTADDATOM).EQ.0) EXITSTATUS = 2
      END SUBROUTINE SET_EXIT_STATUS

      SUBROUTINE GET_STATISTIC_INTERP()
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         IMPLICIT NONE
         REAL(KIND = REAL64) :: MINE, MAXE, SUME, SUME2, MAXRMS, THISE, THISRMS, SIGMAE
         INTEGER :: I, J, JMAX, JMIN

         MINE = 1.0D100
         MAXE = -1.0D100
         MAXRMS = -1.0D0
         SUME = 0.0D0
         SUME2 = 0.0D0

         DO I=2,NIMAGES+1
            THISE = EEE(I)
            SUME = SUME + THISE
            SUME2 = SUME2 + THISE**2
            !WRITE(*,*) "Image: ", I, " energy: ", THISE
            IF (THISE.GT.MAXE) THEN
               MAXE = THISE
               JMAX = I
            END IF
            IF (THISE.LT.MINE) THEN
               MINE = THISE
               JMIN = I
            END IF
            THISRMS = 0.0D0
            DO J=1,3*NATOMS
               THISRMS = THISRMS + GGG(3*NATOMS*(I-1)+J)**2
            END DO
            IF (THISRMS.GT.MAXRMS) THEN
               MAXRMS = THISRMS
            END IF
            !WRITE(*,*) "Image: ", I, " rms: ", THISRMS
         END DO

         MAXRMS = SQRT(MAXRMS/(3*NACTIVE))
         SUME = SUME/NIMAGES
         SUME2 = SUME2/(NIMAGES-1)
         SIGMAE = SQRT(MAX(SUME2-SUME**2,1.0D-100))
         WRITE(*,'(A,I6,A,F20.5,A,F20.5,A)') " get_interp_stat> The highest image ", JMAX, " with energy ", MAXE, " is ", ABS(MAXE-SUME)/SIGMAE, " sigma from the mean"
         WRITE(*,'(A,F20.5,A,F20.5)') "                  The average energy per image is ", SUME, " with variance of ", SUME2
         WRITE(*,*)

      END SUBROUTINE GET_STATISTIC_INTERP

      SUBROUTINE GET_IMAGE_SEPARATION(DMIN,DMAX,JMIN,JMAX)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         USE HELPER_FNCTS, ONLY: DISTANCE_ATOM_DIFF_IMAGES
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(OUT) :: DMIN, DMAX
         INTEGER, INTENT(OUT) :: JMIN, JMAX
         REAL(KIND = REAL64) :: ADMAX
         REAL(KIND = REAL64) :: DISTATOM, DISTTOTAL, X1(3*NATOMS), X2(3*NATOMS)
         INTEGER :: J1, J2, JAMAX_IMG, JAMAX_ATOM

         DMIN = HUGE(1.0D0)
         DMAX = -1.0D0
         ADMAX = -1.0D0

         DO J1=1,NIMAGES+1
            DISTTOTAL = 0.0D0
            X1(1:3*NATOMS) = XYZ((J1-1)*3*NATOMS+1:J1*3*NATOMS)
            X2(1:3*NATOMS) = XYZ(J1*3*NATOMS+1:(J1+1)*3*NATOMS)
            DO J2=1,NATOMS
               IF (ATOMACTIVE(J2)) THEN
                  CALL DISTANCE_ATOM_DIFF_IMAGES(NATOMS, X1, X2, J2, DISTATOM)
                  DISTTOTAL = DISTTOTAL + DISTATOM
                  ! WRITE(*,*) "images", J1, J1+1, "atom ", J2, " distance: ", DISTATOM 
                  IF (DISTATOM.GT.ADMAX) THEN
                     ADMAX = DISTATOM
                     JAMAX_ATOM = J2
                     JAMAX_IMG = J1
                  END IF
               END IF
            END DO 
            WRITE(*,*) "Images ", J1, " and ", J1+1, " - separation: ", DISTTOTAL          
            IF (DISTTOTAL.GT.DMAX) THEN
               DMAX = DISTTOTAL
               JMAX = J1
            END IF
            IF (DISTTOTAL.LT.DMIN) THEN
               DMIN = DISTTOTAL
               JMIN = J1
            END IF
         END DO
         WRITE(*,'(A,F15.5)') " get_image_separation> The largest distance between images is ", DMAX
         WRITE(*,'(A,F15.5)') "                       The smallest distance between images is ", DMIN
         WRITE(*,'(A,I6,A,I4,A,I4)') "                       The largest distance by atom is for atom ",JAMAX_ATOM," between images", JAMAX_IMG," and ", JAMAX_IMG+1
      END SUBROUTINE GET_IMAGE_SEPARATION


      SUBROUTINE INITIALISE_INTERPOLATION_VARS()
         USE QCIKEYS, ONLY: USEIMAGEDENSITY, NIMAGES, E2E_DIST, IMAGEDENSITY, MAXINTIMAGE
         USE ADDINGATOM, ONLY: ALLOC_ADDATOM
         IMPLICIT NONE

         CONACTIVE(:) = .FALSE.
         ATOMACTIVE(1:NATOMS) = .FALSE.
         NTRIES(1:NATOMS) = 0
         TURNONORDER(1:NATOMS) = 0

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

      SUBROUTINE MAKESTEP(NITERDONE,NPT,POINT, RHO1, ALPHA)
         USE QCIKEYS, ONLY: MUPDATE, DGUESS, NATOMS, NIMAGES
         USE MOD_INTCOORDS, ONLY: G, DIMS
         IMPLICIT NONE         
         INTEGER, INTENT(IN) :: NITERDONE
         INTEGER, INTENT(IN) :: NPT
         INTEGER, INTENT(OUT) :: POINT
         REAL(KIND = REAL64), INTENT(INOUT) :: RHO1(MUPDATE), ALPHA(MUPDATE)
         REAL(KIND = REAL64) :: GNORM
         REAL(KIND = REAL64) :: YS, YY, YR, SQ, BETA
         INTEGER :: BOUND, CP, I

         !if it is the first step, we use a cautious guess
         IF (NITERDONE.EQ.1) THEN
            POINT = 0
            DIAG(1:DIMS) = DGUESS
            SEARCHSTEP(0,1:DIMS) = -DGUESS*G(1:DIMS)
            GTMP(1:DIMS) = SEARCHSTEP(0,1:DIMS)
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
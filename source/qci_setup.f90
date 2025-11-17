MODULE QCISETUP
   USE QCIPREC
   CONTAINS
      SUBROUTINE QCI_INIT(PARAMETERFILE, ALIGNT)
         USE QCIKEYS, ONLY: NATOMS, QCIAMBERT, QCIFREEZET, USEIMAGEDENSITY, QCIPERMT
         USE QCICONSTRAINTS, ONLY: CREATE_CONSTRAINTS
         USE QCIPERMDIST, ONLY: INIT_PERMALLOW
         USE MOD_FREEZE, ONLY: GET_FROZEN_ATOMS
         USE CHIRALITY, ONLY: FIND_CHIRAL_CENTRES
         USE REPULSION, ONLY: NREPCURR, ALLOC_REP_VARS
         USE QCI_LINEAR, ONLY: GET_LINEAR_ATOMS
         IMPLICIT NONE
         CHARACTER(30), INTENT(IN) :: PARAMETERFILE
         LOGICAL, INTENT(IN) :: ALIGNT

         ! parse settings
         CALL PARSE_SETTINGS(PARAMETERFILE)
         ! starting permutational setup perm.allow
         IF (QCIPERMT) CALL INIT_PERMALLOW(NATOMS)
         ! align endpoints
         IF (ALIGNT.OR.USEIMAGEDENSITY) CALL ALIGN_ENDPOINTS()
         ! get frozen atoms setup
         IF (QCIFREEZET) CALL GET_FROZEN_ATOMS()
         ! get constraints
         CALL CREATE_CONSTRAINTS()
         ! get atoms for linear interpolation
         CALL GET_LINEAR_ATOMS()
         ! setting up repulsions
         ! we use NATOMS as initial size here for the number of repulsions
         NREPCURR = NATOMS
         CALL ALLOC_REP_VARS(NREPCURR)
         ! get chiral information for AMBER
         IF (QCIAMBERT) CALL FIND_CHIRAL_CENTRES()
         ! set up atom names for output
         CALL SET_ELEMENTS()
      END SUBROUTINE QCI_INIT

      SUBROUTINE SET_ELEMENTS()
         USE QCIKEYS, ONLY: NAMES, NATOMS, ALLOC_QCIKEYS
         USE AMBER_CONSTRAINTS, ONLY: ELEMENT
         IMPLICIT NONE
         
         INTEGER :: LUNIT
         INTEGER :: J1, J2, J3

         CALL ALLOC_QCIKEYS()

         NAMES(1:NATOMS) = 'LA   '

         IF (ALLOCATED(ELEMENT)) THEN
            DO J1=1,NATOMS
               IF (ELEMENT(J1).EQ.1) THEN
                  NAMES(J1) = 'H    '
               ELSE IF (ELEMENT(J1).EQ.6) THEN
                  NAMES(J1) = 'C    '
               ELSE IF (ELEMENT(J1).EQ.7) THEN
                  NAMES(J1) = 'N    '
               ELSE IF (ELEMENT(J1).EQ.8) THEN
                  NAMES(J1) = 'O    '
               ELSE IF (ELEMENT(J1).EQ.15) THEN
                  NAMES(J1) = 'P    '
               ELSE IF (ELEMENT(J1).EQ.16) THEN
                  NAMES(J1) = 'S    '                 
               END IF
            END DO
         END IF
      END SUBROUTINE SET_ELEMENTS

      SUBROUTINE PARSE_SETTINGS(PARAMETERFILE)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         USE HELPER_FNCTS, ONLY: READ_LINE
         IMPLICIT NONE
         CHARACTER(30), INTENT(IN) :: PARAMETERFILE
         INTEGER :: PARAMUNIT, IOS
         LOGICAL :: YESNO, EOFT
         INTEGER, PARAMETER :: NWORDS = 20
         CHARACTER(25) :: ENTRIES(NWORDS)=''
         CHARACTER(200) :: LINE
            
         ! open parameter file
         INQUIRE(FILE=PARAMETERFILE, EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(*,*) " parse_settings> Cannot locate QCI settings file ", PARAMETERFILE
            STOP
         END IF
         PARAMUNIT = GETUNIT()
         OPEN(PARAMUNIT, FILE=PARAMETERFILE, STATUS="OLD")
            
         !loop over lines in file
         EOFT = .FALSE.
         DO WHILE (.NOT. EOFT)
            READ(PARAMUNIT,'(A)',IOSTAT=IOS) LINE
            IF (IOS.NE.0) THEN
               EOFT = .TRUE.     
            ELSE
               CALL READ_LINE(LINE,NWORDS,ENTRIES)
               IF (ENTRIES(1).NE."") CALL SETKEYS(ENTRIES(1), ENTRIES(2))
            ENDIF
         END DO
         CLOSE(PARAMUNIT)
      END SUBROUTINE PARSE_SETTINGS

      SUBROUTINE ALIGN_ENDPOINTS()
         USE QCIPERMDIST, ONLY: LOPERMDIST
         USE QCIMINDIST, ONLY: ALIGNXBTOA
         USE QCIKEYS, ONLY : NATOMS, E2E_DIST, QCIPERMT, DEBUG
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         IMPLICIT NONE
         REAL(KIND=REAL64) :: DIST2, RMATBEST(3,3)
         INTEGER :: NMOVE, NEWPERM(NATOMS)

         !QUERY: will lopermdist actually change coordinates or do we need a wrapper to do so?
         IF (QCIPERMT) THEN
            CALL LOPERMDIST(XFINAL,XSTART,E2E_DIST,DIST2,RMATBEST,0,NMOVE,NEWPERM)
            WRITE(*,*) " align_endpoints> Distance between endpoints is ", E2E_DIST
         ELSE
            CALL ALIGNXBTOA(XSTART, XFINAL, NATOMS)
         END IF
         IF (DEBUG) THEN
            OPEN(UNIT=55,FILE="start.aligned")
            WRITE(55,'(3F20.7)') XSTART
            CLOSE(55)
            OPEN(UNIT=55,FILE="finish.aligned")
            WRITE(55,'(3F20.7)') XFINAL
            CLOSE(55)         
         END IF
      END SUBROUTINE ALIGN_ENDPOINTS

      SUBROUTINE SETKEYS(ENTRY, VAL)
         USE QCIKEYS
         USE QCI_CONSTRAINT_KEYS, ONLY: MAXCONUSE, QCICONSEP, QCICONSTRAINTTOL, QCICONCUT, GEOMFILE, CONSTRFILE, &
                                        CONCUTABST, CONCUTFRACT, CONCUTABSINC, CONCUTABS, CONCUTFRAC 
         USE AMBER_CONSTRAINTS, ONLY: AMBERCONSTRFILE, TOPFILENAME
         USE HIRE_CONSTRAINTS, ONLY: HIRECONSTRFILE, HIRETOPFILE 
         USE SBM_CONSTRAINTS, ONLY: SBMCONTACTFILE
         USE QCI_LINEAR, ONLY: LINEARCUT, LINEARFILE
         USE REPULSION, ONLY: CHECKREPCUTOFF
         USE ADDINGATOM, ONLY: QCIDISTCUT, QCIATOMSEP
         USE DIHEDRAL_CONSTRAINTS, ONLY: KDIH
         IMPLICIT NONE
         CHARACTER(25), INTENT(IN) :: ENTRY, VAL

         IF (ENTRY.EQ."COMMENT") THEN
            RETURN
         ! basic set up
         ! enable debug
         ELSE IF (ENTRY.EQ."DEBUG") THEN
            DEBUG = .TRUE.
         ! number of images to start with
         ELSE IF (ENTRY.EQ."NIMAGES") THEN
            READ(VAL, *) NIMAGES
         ! maximum number of images
         ELSE IF (ENTRY.EQ."MAXINTIMAGES") THEN
            READ(VAL, *) MAXINTIMAGE
         ! maximum number of iterations for band optimisation
         ELSE IF (ENTRY.EQ."MAXITERATIONS") THEN
            READ(VAL, *) MAXITER
         ! which potential is used?
         ELSE IF (ENTRY.EQ."QCIMODE") THEN
            IF (VAL.EQ."AMBER") THEN
               QCIAMBERT = .TRUE.
               CHECKCHIRAL = .TRUE.
            ELSE IF (VAL.EQ."HIRE") THEN
               QCIHIRET = .TRUE.
            ELSE IF (VAL.EQ."SBM") THEN
               QCISBMT = .TRUE.
            ELSE IF (VAL.EQ."GEOMETRY") THEN
               QCIGEOMT = .TRUE.
            ELSE 
               WRITE(*,*) " setkeys> QCI mode ", VAL, " is not a valid function"
            END IF
         ! options for these potentials
         ELSE IF (ENTRY.EQ."AMBERCONSTRFILE") THEN
            AMBERCONSTRFILE = VAL
         ELSE IF (ENTRY.EQ."DETECTBASEPAIRS") THEN
            BASEPAIRDETECTION = .TRUE.
         ELSE IF (ENTRY.EQ."TOPFILENAME") THEN
            TOPFILENAME = VAL
         ELSE IF (ENTRY.EQ."HIRECONSTRFILE") THEN
            HIRECONSTRFILE = VAL
         ELSE IF (ENTRY.EQ."HIRETOPFILE") THEN
            HIRETOPFILE = VAL
         ELSE IF (ENTRY.EQ."SBMCONTACTFILE") THEN
            SBMCONTACTFILE = VAL
         ! Strategies of how atoms are activated
         ELSE IF (ENTRY.EQ."QCIDOBACK") THEN
            QCIDOBACK = .TRUE. 
         ELSE IF (ENTRY.EQ."QCIDOBACKALL") THEN
            QCIDOBACKALL = .TRUE.             
         ELSE IF (ENTRY.EQ."QCIADDACID") THEN
            QCIADDACIDT = .TRUE.      
         ELSE IF (ENTRY.EQ."GEOMFILE") THEN
            GEOMFILE = VAL
         ELSE IF (ENTRY.EQ."ADDATOMINTERVAL") THEN
            READ(VAL,*) INTADDATOM
         !limits on finding atoms for local axis set
         ELSE IF (ENTRY.EQ."DISTCUTADDATOM") THEN
            READ(VAL,*) QCIDISTCUT
         ELSE IF (ENTRY.EQ."MAXSEPADDATOM") THEN
            READ(VAL,*) QCIATOMSEP
         !use minimisation after adding atom
         ELSE IF (ENTRY.EQ."MINIMISEAFTERADD") THEN
            OPTIMISEAFTERADDITION = .TRUE.
            READ(VAL, *) NMINAFTERADD
         ! use trilateration
         ELSE IF (ENTRY.EQ."TRILATERATE") THEN
            QCITRILATERATION = .TRUE.
            USEINTERNALST = .FALSE.
         ! use internal corodinates for inteprolation
         ELSE IF (ENTRY.EQ."USEINTERNALS") THEN
            USEINTERNALST = .TRUE.
         ! use four atom basis for interpolation
         ELSE IF (ENTRY.EQ."USEFOURATOMS") THEN
            USEFOURATOMST = .TRUE.
         ! use linear atoms
         ELSE IF (ENTRY.EQ."QCILINEAR") THEN
            QCILINEART = .TRUE.
            LINEARFILE = VAL
         ELSE IF (ENTRY.EQ."LINEARCUTOFF") THEN
            READ(VAL, *) LINEARCUT
         ELSE IF (ENTRY.EQ."LINEARBB") THEN
            LINEARBBT = .TRUE.           
         !reading a guess for the band
         ELSE IF (ENTRY.EQ."QCIREADGUESS") THEN
            QCIREADGUESS = .TRUE.   
            GUESSFILE = VAL
         ! permutational variables
         ELSE IF (ENTRY.EQ."QCIPERMCHECK") THEN
            QCIPERMT = .TRUE.
            READ(VAL, *) QCIPERMCHECKINT
         ELSE IF (ENTRY.EQ."QCIPERMCUT") THEN
            READ(VAL, *) QCIPERMCUT
         ! use of frozen atoms
         ELSE IF (ENTRY.EQ."FREEZEFILE") THEN
            QCIFREEZET = .TRUE.            
            FREEZEFILE = VAL
         ELSE IF (ENTRY.EQ."NMINUNFROZEN") THEN
            READ(VAL, *) NMINUNFROZEN
         ELSE IF (ENTRY.EQ."FREEZEFILE") THEN
            FREEZEFILE = VAL
         ELSE IF (ENTRY.EQ."QCIFREEZE") THEN
            QCIFREEZET = .TRUE.
            READ(VAL, *) QCIFREEZETOL
         ! check for internal minima in constraints?
         ELSE IF (ENTRY.EQ."CHECKINTMINCONSTR") THEN
            CHECKCONINT = .TRUE.
         ! scaling factor for internal minima
         ELSE IF (ENTRY.EQ."INTMINFACTOR") THEN 
            READ(VAL,*) INTMINFAC 
         ! minimum distacne in sequence for internal minima
         ELSE IF (ENTRY.EQ."REPINTMINSEP") THEN 
            READ(VAL,*) QCIINTREPMINSEP      
         ! ????
         ELSE IF (ENTRY.EQ."INTCONSTRAINTDEL") THEN 
            READ(VAL,*) INTCONSTRAINTDEL   
         ! image control
         ELSE IF (ENTRY.EQ."USEIMAGEDENSITY") THEN  
            USEIMAGEDENSITY = .TRUE.
            READ(VAL,*) IMAGEDENSITY
         ELSE IF (ENTRY.EQ."MAXIMAGESEPARATION") THEN
            READ(VAL,*) IMSEPMAX
         ELSE IF (ENTRY.EQ."MINIMAGESEPARATION") THEN
            READ(VAL,*) IMSEPMIN
         ! spring constant
         ELSE IF (ENTRY.EQ."KSPRING") THEN 
            READ(VAL,*) KINT
         ELSE IF (ENTRY.EQ."KSPRINGSCALING") THEN
            READ(VAL,*) KINTSCALED
         ELSE IF (ENTRY.EQ."ADJUSTSPRING") THEN  
            QCIADJUSTKT = .TRUE.
            READ(VAL,*) QCIADJUSTKFRQ
         ELSE IF (ENTRY.EQ."KSPACINGDEV") THEN
            READ(VAL,*) QCIADJUSTKTOL
         ELSE IF (ENTRY.EQ."KMIN") THEN
            READ(VAL,*) QCIKINTMIN
         ELSE IF (ENTRY.EQ."KMAX") THEN
            READ(VAL,*) QCIKINTMAX
         ELSE IF (ENTRY.EQ."KADJUSTFRAC") THEN
            READ(VAL,*) QCIADJUSTKFRAC            
         ! constraints
         ELSE IF (ENTRY.EQ."MAXCONUSE") THEN
            READ(VAL, *) MAXCONUSE        
         ELSE IF (ENTRY.EQ."QCICONSEP") THEN
            READ(VAL, *) QCICONSEP
         ELSE IF (ENTRY.EQ."QCICONSTRAINTTOL") THEN
            READ(VAL, *) QCICONSTRAINTTOL
         ELSE IF (ENTRY.EQ."QCICONCUT") THEN
            READ(VAL, *) QCICONCUT
         ELSE IF (ENTRY.EQ."CONCUTABS") THEN
            CONCUTABST = .TRUE.
            READ(VAL, *) CONCUTABS
         ELSE IF (ENTRY.EQ."CONCUTABSINC") THEN
            CONCUTABSINC = .TRUE.
         ELSE IF (ENTRY.EQ."CONCUTFRAC") THEN
            CONCUTFRACT = .TRUE.
            READ(VAL, *) CONCUTFRAC
         ELSE IF ( ENTRY.EQ."CONACTINACT") THEN
            USECONACTINACT = .TRUE.
            READ(VAL,*) CONACTINACT 
         ELSE IF (ENTRY.EQ."DIHEDRALCONSTR") THEN
            USEDIHEDRALCONST = .TRUE.
            READ(VAL,*) KDIH
         !repulsions
         ELSE IF (ENTRY.EQ."REPULSIONCUTOFF") THEN
            READ(VAL, *) QCIREPCUT
         ELSE IF (ENTRY.EQ."QCICONSTRREP") THEN
            READ(VAL, *) QCICONSTRREP
         !
         ELSE IF (ENTRY.EQ."CHECKREPCUTOFF") THEN
            READ(VAL,*) CHECKREPCUTOFF

         !maximum gradient component
         ELSE IF (ENTRY.EQ."MAXGRADCOMP") THEN
            READ(VAL, *) MAXGRADCOMP

         !settings for step taking
         ELSE IF (ENTRY.EQ."DIAGGUESS") THEN
            READ(VAL, *) DGUESS
         ELSE IF (ENTRY.EQ."UPDATES") THEN
            READ(VAL, *) MUPDATE
         ELSE IF (ENTRY.EQ."MAXQCIBFGS") THEN
            READ(VAL, *) MAXQCIBFGS
         ELSE IF (ENTRY.EQ."COLDFUSIONLIMIT") THEN
            READ(VAL, *) COLDFUSIONLIMIT

         !maximum constraint E - used for convergence
         ELSE IF (ENTRY.EQ."MAXCONSTRAINTE") THEN
            READ(VAL, *) MAXCONE
         !convergence for RMS
         ELSE IF (ENTRY.EQ."RMSTOLERANCE") THEN
            READ(VAL, *) QCIRMSTOL
         ELSE IF (ENTRY.EQ."MAXERISE") THEN
            READ(VAL, *) MAXERISE
         
         !QCI resetting
         ELSE IF (ENTRY.EQ."QCIRESET") THEN
            QCIRESET = .TRUE.
            READ(VAL, *) QCIRESETINT1
         !checking repulsions
         ELSE IF (ENTRY.EQ."CHECKREPINTERVAL") THEN
            READ(VAL, *) CHECKREPINTERVAL
         ! options for saving band
         ELSE IF (ENTRY.EQ."DUMPXYZ") THEN
            DUMPQCIXYZ = .TRUE.
            READ(VAL, *) DUMPQCIXYZFRQS
 
         ELSE
            WRITE(*,*) " setkeys> Cannot find setting ", ENTRY, " - will skip this entry"
         END IF
      END SUBROUTINE SETKEYS

END MODULE QCISETUP

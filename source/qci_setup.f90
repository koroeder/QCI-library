MODULE QCISETUP
   USE QCIPREC
   CONTAINS
      SUBROUTINE QCI_INIT(PARAMETERFILE, ALIGNT)
         USE QCIKEYS, ONLY: NATOMS, QCIAMBERT, QCIFREEZET, USEIMAGEDENSITY
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
         CALL INIT_PERMALLOW(NATOMS)
         ! align endpoints
         IF (ALIGNT.OR.USEIMAGEDENSITY) CALL ALIGN_ENDPOINTS()
         ! get atoms for linear interpolation
         CALL GET_LINEAR_ATOMS()
         ! get frozen atoms setup
         IF (QCIFREEZET) CALL GET_FROZEN_ATOMS()
         ! get constraints
         CALL CREATE_CONSTRAINTS()
         ! setting up repulsions
         ! we use NATOMS as initial size here for the number of repulsions
         NREPCURR = NATOMS
         CALL ALLOC_REP_VARS(NREPCURR)
         ! get chiral information for AMBER
         IF (QCIAMBERT) CALL FIND_CHIRAL_CENTRES()
      END SUBROUTINE QCI_INIT

      SUBROUTINE PARSE_SETTINGS(PARAMETERFILE)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         IMPLICIT NONE
         CHARACTER(30), INTENT(IN) :: PARAMETERFILE
         INTEGER :: PARAMUNIT, IOS
         LOGICAL :: YESNO, EOFT
         CHARACTER(25) :: ENTRY, VALUE 
            
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
            READ(PARAMUNIT,*,IOSTAT=IOS) ENTRY, VALUE
            IF (IOS.GT.0) THEN
               EOFT = .TRUE.     
            ELSE
               CALL SETKEYS(ENTRY, VALUE)
            ENDIF
         END DO
         CLOSE(PARAMUNIT)
      END SUBROUTINE PARSE_SETTINGS

      SUBROUTINE ALIGN_ENDPOINTS()
         USE QCIPERMDIST, ONLY: LOPERMDIST
         USE QCIKEYS, ONLY : NATOMS, E2E_DIST
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL
         IMPLICIT NONE
         REAL(KIND=REAL64) :: DIST2, RMATBEST(3,3)
         INTEGER :: NMOVE, NEWPERM(NATOMS)

         !QUERY: will lopermdist actually change coordinates or do we need a wrapper to do so?
         CALL LOPERMDIST(XFINAL,XSTART,E2E_DIST,DIST2,RMATBEST,0,NMOVE,NEWPERM)


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
         ! use trilateration
         ELSE IF (ENTRY.EQ."TRILATERATE") THEN
            QCITRILATERATION = .TRUE.
         ! use linear atoms
         ELSE IF (ENTRY.EQ."QCILINEAR") THEN
            QCILINEART = .TRUE.
            LINEARFILE = VAL
         ELSE IF (ENTRY.EQ."LINEARCUTOFF") THEN
            READ(VAL, *) LINEARCUT
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
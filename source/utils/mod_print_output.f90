!! This modules saves energies calculated by congrad for wiritng output

MODULE OUT_PRINT

    USE QCIPREC
    IMPLICIT NONE

    REAL(KIND = REAL64) :: SAVEECON, SAVEEREP, SAVEESPR, SAVEEDIH  
    REAL(KIND = REAL64) :: SAVEETOTAL                    ! overall energy
    REAL(KIND = REAL64) :: SAVERMS  
    REAL(KIND = REAL64) :: SAVEFCONTEST, SAVEFREPTEST, SAVEFDIHTEST, SAVEFSPRINGMAX
    REAL(KIND = REAL64) :: SAVECONVERGECONTEST, SAVECONVERGEREPTEST, SAVEEDIHTEST

    CONTAINS
    
    SUBROUTINE SAVE_OUT(ETOTAL, RMS, EREP, ECON, ESPR, EDIH, FCONTEST, FREPTEST, FDIHTEST, CONVERGECONTEST, CONVERGEREPTEST, CONVERGENCEDIHTEST,FSPRINGTEST)
        IMPLICIT NONE
        
        REAL(KIND = REAL64), INTENT(IN) :: ETOTAL  
        REAL(KIND = REAL64), INTENT(IN) :: EREP, ECON, ESPR, EDIH
        REAL(KIND = REAL64), INTENT(IN) :: FCONTEST, FREPTEST, FDIHTEST, FSPRINGTEST
        REAL(KIND = REAL64), INTENT(IN) :: CONVERGECONTEST, CONVERGEREPTEST, CONVERGENCEDIHTEST
        REAL(KIND = REAL64), INTENT(IN) :: RMS                       ! total force
        
        SAVEETOTAL = ETOTAL
        SAVEECON = ECON
        SAVEEREP = EREP
        SAVEEDIH = EDIH
        SAVEESPR  = ESPR

        SAVERMS = RMS

        SAVEFCONTEST = FCONTEST
        SAVEFREPTEST = FREPTEST
        SAVEFDIHTEST = FDIHTEST
        SAVECONVERGECONTEST = CONVERGECONTEST
        SAVECONVERGEREPTEST = CONVERGEREPTEST
        SAVEEDIHTEST = CONVERGENCEDIHTEST
        SAVEFSPRINGMAX = FSPRINGTEST

    END SUBROUTINE SAVE_OUT

    SUBROUTINE WRITE_CONGRADOUT()
        IMPLICIT NONE
         WRITE(*,*) " congrad> E total: ", SAVEETOTAL, "RMS: ", SAVERMS, " E rep: ", SAVEEREP, " E constr: ", SAVEECON
         WRITE(*,*) "                                              E spring: ", SAVEESPR, " E dih: ", SAVEEDIH
         WRITE(*,*) " congrad> FCONMAX: ", SAVEFCONTEST, " FREPMAX: ", SAVEFREPTEST, " FDIHMAX: ", SAVEFDIHTEST, "F_SPRING_MAX", SAVEFSPRINGMAX
         WRITE(*,*) " congrad> CONVERGECONTEST: ", SAVECONVERGECONTEST, " CONVERGEREPTEST: ", SAVECONVERGEREPTEST, " CONVERGEDIHTEST: ", SAVEEDIHTEST
    END SUBROUTINE WRITE_CONGRADOUT

    SUBROUTINE WRITE_QCI_KEYS()
        USE QCIKEYS
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
        
        INTEGER :: IU, I, N, NTRUE
        INTEGER, PARAMETER :: MAX_PRINT_ARRAY = 20

        ! Format strings for consistency
        CHARACTER(LEN=*), PARAMETER :: FMT_SEC  = '(1X,A)'
        CHARACTER(LEN=*), PARAMETER :: FMT_BOOL = '(3X,A, T50, "=", L1)'
        CHARACTER(LEN=*), PARAMETER :: FMT_INT  = '(3X,A, T50, "=", I0)'
        CHARACTER(LEN=*), PARAMETER :: FMT_REAL = '(3X,A, T50, "=", G0.6)'
        CHARACTER(LEN=*), PARAMETER :: FMT_CHAR = '(3X,A, T50, "=", A)'
   
        IU = OUTPUT_UNIT
        
        WRITE(IU,FMT_SEC) '============================================================'
        WRITE(IU,FMT_SEC) '           EFFECTIVE CONFIGURATION (QCI)'
        WRITE(IU,FMT_SEC) '============================================================'

        !============================================================================
        ! General control parameters
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- General control ---'
        WRITE(IU,FMT_BOOL) 'DEBUG', DEBUG
        WRITE(IU,FMT_INT)  'NATOMS', NATOMS
        WRITE(IU,FMT_INT)  'NIMAGES', NIMAGES
        WRITE(IU,FMT_INT)  'MAXITER', MAXITER
        WRITE(IU,FMT_INT)  'MAXINTIMAGE', MAXINTIMAGE

        !============================================================================
        ! QCI type flags
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- QCI type flags ---'
        WRITE(IU,FMT_BOOL) 'QCIAMBERT', QCIAMBERT
        WRITE(IU,FMT_BOOL) 'QCIHIRET', QCIHIRET
        WRITE(IU,FMT_BOOL) 'QCISBMT', QCISBMT
        WRITE(IU,FMT_BOOL) 'QCIGEOMT', QCIGEOMT

        !============================================================================
        ! Backbone settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Backbone settings ---'
        WRITE(IU,FMT_BOOL) 'QCIDOBACK', QCIDOBACK
        IF (ALLOCATED(ISBBATOM)) THEN
            N = SIZE(ISBBATOM)
            NTRUE = COUNT(ISBBATOM)
            WRITE(IU,FMT_INT) 'ISBBATOM (size)', N
            WRITE(IU,FMT_INT) 'ISBBATOM (.TRUE. count)', NTRUE
            IF (N <= MAX_PRINT_ARRAY) THEN
                WRITE(IU,'(3X,A45,1X,"=",1X,*(L1,1X))') 'ISBBATOM', (ISBBATOM(I), I=1,N)
            ELSE
                WRITE(IU,'(3X,A)') 'ISBBATOM: array too large to print in full.'
            END IF
        ELSE
            WRITE(IU,FMT_CHAR) 'ISBBATOM', 'not allocated'
        END IF
        WRITE(IU,FMT_INT)  'NBACKBONE', NBACKBONE
        WRITE(IU,FMT_BOOL) 'DETECTBBCROSSING', DETECTBBCROSSING
        WRITE(IU,FMT_INT)  'CHECKCROSSFREQ', CHECKCROSSFREQ

        !============================================================================
        ! Atom names
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Atom names ---'
        IF (ALLOCATED(NAMES)) THEN
            N = SIZE(NAMES)
            WRITE(IU,FMT_INT) 'NAMES (size)', N
            IF (N <= MAX_PRINT_ARRAY) THEN
                WRITE(IU,'(3X,A45,1X,"=",1X,*(A5,1X))') 'NAMES', (NAMES(I), I=1,N)
            ELSE
                WRITE(IU,'(3X,A45,1X,"=",1X,5(A5,1X),"...")') 'NAMES(1:5)', (NAMES(I), I=1,5)
            END IF
        ELSE
            WRITE(IU,FMT_CHAR) 'NAMES', 'not allocated'
        END IF

        !============================================================================
        ! Addition / optimisation logic
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Addition / optimisation logic ---'
        WRITE(IU,FMT_BOOL) 'QCIADDACIDT', QCIADDACIDT
        WRITE(IU,FMT_BOOL) 'QCIDOBACKALL', QCIDOBACKALL
        WRITE(IU,FMT_BOOL) 'OPTIMISEAFTERADDITION', OPTIMISEAFTERADDITION
        WRITE(IU,FMT_INT)  'NMINAFTERADD', NMINAFTERADD

        !============================================================================
        ! Interpolation methods
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Interpolation methods ---'
        WRITE(IU,FMT_BOOL) 'USEINTERNALST', USEINTERNALST
        WRITE(IU,FMT_BOOL) 'USEFOURATOMST', USEFOURATOMST
        WRITE(IU,FMT_BOOL) 'QCITRILATERATION', QCITRILATERATION
        WRITE(IU,FMT_BOOL) 'QCITRILATERATION2', QCITRILATERATION2

        !============================================================================
        ! Constraint / gradient settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Constraint / gradient settings ---'
        WRITE(IU,FMT_BOOL) 'CHECKCONINT', CHECKCONINT
        WRITE(IU,FMT_REAL) 'INTCONSTRAINTDEL', INTCONSTRAINTDEL
        WRITE(IU,FMT_REAL) 'MAXGRADCOMP', MAXGRADCOMP

        !============================================================================
        ! Guess file
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Guess file ---'
        WRITE(IU,FMT_BOOL) 'QCIREADGUESS', QCIREADGUESS
        WRITE(IU,FMT_CHAR) 'GUESSFILE', TRIM(GUESSFILE)

        !============================================================================
        ! Image density / separation
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Image density / separation ---'
        WRITE(IU,FMT_BOOL) 'USEIMAGEDENSITY', USEIMAGEDENSITY
        WRITE(IU,FMT_REAL) 'E2E_DIST', E2E_DIST
        WRITE(IU,FMT_REAL) 'IMAGEDENSITY', IMAGEDENSITY
        WRITE(IU,FMT_REAL) 'IMSEPMAX', IMSEPMAX
        WRITE(IU,FMT_REAL) 'IMSEPMIN', IMSEPMIN

        !============================================================================
        ! Linear atom / group settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Linear atom / group settings ---'
        WRITE(IU,FMT_BOOL) 'QCILINEART', QCILINEART
        IF (ALLOCATED(INLINLIST)) THEN
            N = SIZE(INLINLIST)
            NTRUE = COUNT(INLINLIST)
            WRITE(IU,FMT_INT) 'INLINLIST (size)', N
            WRITE(IU,FMT_INT) 'INLINLIST (.TRUE. count)', NTRUE
            IF (N <= MAX_PRINT_ARRAY) THEN
                WRITE(IU,'(3X,A45,1X,"=",1X,*(L1,1X))') 'INLINLIST', (INLINLIST(I), I=1,N)
            END IF
        ELSE
            WRITE(IU,FMT_CHAR) 'INLINLIST', 'not allocated'
        END IF
        WRITE(IU,FMT_BOOL) 'LINEARBBT', LINEARBBT
        WRITE(IU,FMT_BOOL) 'USELINGROUPS', USELINGROUPS

        !============================================================================
        ! Residue mapping
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Residue mapping ---'
        IF (ALLOCATED(ATOMS2RES)) THEN
            N = SIZE(ATOMS2RES)
            WRITE(IU,FMT_INT) 'ATOMS2RES (size)', N
            IF (N <= MAX_PRINT_ARRAY) THEN
                WRITE(IU,'(3X,A45,1X,"=",1X,*(I0,1X))') 'ATOMS2RES', (ATOMS2RES(I), I=1,N)
            END IF
        ELSE
            WRITE(IU,FMT_CHAR) 'ATOMS2RES', 'not allocated'
        END IF

        !============================================================================
        ! Frozen atoms
        !============================================================================
        !WRITE(IU,FMT_SEC) ''
        !WRITE(IU,FMT_SEC) '--- Frozen atoms ---'
        !WRITE(IU,FMT_BOOL) 'QCIFREEZET', QCIFREEZET
        !WRITE(IU,FMT_INT)  'NQCIFROZEN', NQCIFROZEN
        !WRITE(IU,FMT_INT)  'NMINUNFROZEN', NMINUNFROZEN
        !WRITE(IU,FMT_CHAR) 'FREEZEFILE', TRIM(FREEZEFILE)
        !WRITE(IU,FMT_REAL) 'QCIFREEZETOL', QCIFREEZETOL
        !IF (ALLOCATED(QCIFROZEN)) THEN
        !    N = SIZE(QCIFROZEN)
        !    NTRUE = COUNT(QCIFROZEN)
        !    WRITE(IU,FMT_INT) 'QCIFROZEN (size)', N
        !    WRITE(IU,FMT_INT) 'QCIFROZEN (.TRUE. count)', NTRUE
        !    IF (N <= MAX_PRINT_ARRAY) THEN
        !        WRITE(IU,'(3X,A45,1X,"=",1X,*(L1,1X))') 'QCIFROZEN', (QCIFROZEN(I), I=1,N)
        !    END IF
        !ELSE
        !    WRITE(IU,FMT_CHAR) 'QCIFROZEN', 'not allocated'
        !END IF
        !IF (ALLOCATED(FREEZE)) THEN
        !    N = SIZE(FREEZE)
        !    NTRUE = COUNT(FREEZE)
        !    WRITE(IU,FMT_INT) 'FREEZE (size)', N
        !    WRITE(IU,FMT_INT) 'FREEZE (.TRUE. count)', NTRUE
        !    IF (N <= MAX_PRINT_ARRAY) THEN
        !        WRITE(IU,'(3X,A45,1X,"=",1X,*(L1,1X))') 'FREEZE', (FREEZE(I), I=1,N)
        !    END IF
        !ELSE
        !    WRITE(IU,FMT_CHAR) 'FREEZE', 'not allocated'
        !END IF

        !============================================================================
        ! Image checking / internal minima
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Image checking / internal minima ---'
        WRITE(IU,FMT_INT)  'QCIIMAGECHECK', QCIIMAGECHECK
        WRITE(IU,FMT_REAL) 'INTMINFAC', INTMINFAC

        !============================================================================
        ! Repulsion settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Repulsion settings ---'
        WRITE(IU,FMT_REAL) 'QCIREPCUT', QCIREPCUT
        WRITE(IU,FMT_REAL) 'QCICONSTRREP', QCICONSTRREP
        WRITE(IU,FMT_INT)  'QCIINTREPMINSEP', QCIINTREPMINSEP
        WRITE(IU,FMT_REAL) 'K_REP', K_REP
        WRITE(IU,FMT_REAL) 'K_CONST', K_CONST

        !============================================================================
        ! Conactinact
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Conactinact ---'
        WRITE(IU,FMT_BOOL) 'USECONACTINACT', USECONACTINACT
        WRITE(IU,FMT_REAL) 'CONACTINACT', CONACTINACT

        !============================================================================
        ! Dihedral constraints
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Dihedral constraints ---'
        WRITE(IU,FMT_BOOL) 'USEDIHEDRALCONST', USEDIHEDRALCONST
        WRITE(IU,FMT_REAL) 'DIHDIFTOL', DIHDIFTOL

        !============================================================================
        ! Spring constants
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Spring constants ---'
        WRITE(IU,FMT_BOOL) 'QCISPRINGACTIVET', QCISPRINGACTIVET
        WRITE(IU,FMT_BOOL) 'QCIADJUSTKT', QCIADJUSTKT
        WRITE(IU,FMT_INT)  'QCIADJUSTKFRQ', QCIADJUSTKFRQ
        WRITE(IU,FMT_REAL) 'QCIADJUSTKTOL', QCIADJUSTKTOL
        WRITE(IU,FMT_REAL) 'QCIAVDEV', QCIAVDEV
        WRITE(IU,FMT_REAL) 'KINT', KINT
        WRITE(IU,FMT_REAL) 'KINTSCALED', KINTSCALED
        WRITE(IU,FMT_REAL) 'QCIKINTMIN', QCIKINTMIN
        WRITE(IU,FMT_REAL) 'QCIKINTMAX', QCIKINTMAX
        WRITE(IU,FMT_REAL) 'QCIADJUSTKFRAC', QCIADJUSTKFRAC

        !============================================================================
        ! BFGS / step settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- BFGS / step settings ---'
        WRITE(IU,FMT_REAL) 'DGUESS', DGUESS
        WRITE(IU,FMT_INT)  'MUPDATE', MUPDATE
        WRITE(IU,FMT_REAL) 'MAXQCIBFGS', MAXQCIBFGS

        !============================================================================
        ! Convergence / energy limits
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Convergence / energy limits ---'
        WRITE(IU,FMT_REAL) 'COLDFUSIONLIMIT', COLDFUSIONLIMIT
        WRITE(IU,FMT_REAL) 'MAXCONE', MAXCONE
        WRITE(IU,FMT_REAL) 'QCIRMSTOL', QCIRMSTOL
        WRITE(IU,FMT_REAL) 'MAXERISE', MAXERISE

        !============================================================================
        ! QCI resetting
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- QCI resetting ---'
        WRITE(IU,FMT_BOOL) 'QCIRESET', QCIRESET
        WRITE(IU,FMT_INT)  'QCIRESETINT1', QCIRESETINT1

        !============================================================================
        ! Misc flags / intervals
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Misc flags / intervals ---'
        WRITE(IU,FMT_INT)  'CHECKREPINTERVAL', CHECKREPINTERVAL
        WRITE(IU,FMT_BOOL) 'CHECKCHIRAL', CHECKCHIRAL
        WRITE(IU,FMT_BOOL) 'QCIUSEGROUPS', QCIUSEGROUPS
        WRITE(IU,FMT_BOOL) 'BASEPAIRDETECTION', BASEPAIRDETECTION
        WRITE(IU,FMT_INT)  'DUMPQCIXYZFRQS', DUMPQCIXYZFRQS
        WRITE(IU,FMT_BOOL) 'DUMPQCIXYZ', DUMPQCIXYZ

        !============================================================================
        ! Permutational settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Permutational settings ---'
        WRITE(IU,FMT_INT)  'QCIPERMCHECKINT', QCIPERMCHECKINT
        WRITE(IU,FMT_BOOL) 'QCIPERMT', QCIPERMT
        WRITE(IU,FMT_REAL) 'QCIPERMCUT', QCIPERMCUT
        WRITE(IU,FMT_REAL) 'ORBITTOL', ORBITTOL

        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '============================================================'
        WRITE(IU,FMT_SEC) '        END OF EFFECTIVE CONFIGURATION'
        WRITE(IU,FMT_SEC) '============================================================'



    END SUBROUTINE WRITE_QCI_KEYS

END MODULE OUT_PRINT

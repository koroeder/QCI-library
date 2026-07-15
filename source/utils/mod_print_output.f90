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
        USE INTERPOLATION_KEYS
        USE QCI_CONSTRAINT_KEYS
        USE AMBER_CONSTRAINTS, ONLY: TOPFILENAME, AMBERCONSTRFILE
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
        WRITE(IU,FMT_SEC) '--- General parameters ---'
        WRITE(IU,FMT_BOOL) 'DEBUG', DEBUG
        WRITE(IU,FMT_INT)  'NATOMS', NATOMS
        WRITE(IU,FMT_INT)  'NIMAGES', NIMAGES
        WRITE(IU,FMT_INT)  'MAXINTIMAGE', MAXINTIMAGE
        WRITE(IU,FMT_INT)  'MAXITERATIONS', MAXITER

        !============================================================================
        ! Guess file
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Guess file ---'
        WRITE(IU,FMT_BOOL) 'QCIREADGUESS', QCIREADGUESS
        WRITE(IU,FMT_CHAR) 'GUESSFILE', TRIM(GUESSFILE)
        
        !============================================================================
        ! Topology
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- QCI topology flags ---'
        WRITE(IU,FMT_BOOL) 'QCIAMBERT', QCIAMBERT
        WRITE(IU,FMT_BOOL) 'QCIHIRET', QCIHIRET
        WRITE(IU,FMT_BOOL) 'QCISBMT', QCISBMT
        WRITE(IU,FMT_BOOL) 'QCIGEOMT', QCIGEOMT

        IF (QCIAMBERT) THEN
            WRITE(IU,FMT_SEC) ''
            WRITE(IU,FMT_CHAR) 'TOPOLOGY', TRIM(TOPFILENAME)
            WRITE(IU,FMT_CHAR) 'CONSTRAINTS', TRIM(AMBERCONSTRFILE)
            WRITE(IU,FMT_BOOL) 'QCIUSEGROUPS', QCIUSEGROUPS
            WRITE(IU,FMT_BOOL) 'BASEPAIRDETECTION', BASEPAIRDETECTION
        END IF

        !Congeom only
        IF (QCIGEOMT) THEN
            WRITE(IU,FMT_SEC) ''
            WRITE(IU,FMT_REAL) 'QCICONSTRAINTTOL', QCICONSTRAINTTOL
        END IF

        !============================================================================
        ! Atom adding options
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Atom adding options ---'

        WRITE(IU,FMT_BOOL) 'USEINTERNALST', USEINTERNALST
        WRITE(IU,FMT_BOOL) 'USEFOURATOMST', USEFOURATOMST
        WRITE(IU,FMT_BOOL) 'QCITRILATERATION', QCITRILATERATION
        WRITE(IU,FMT_BOOL) 'QCITRILATERATION2', QCITRILATERATION2

        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_BOOL) 'QCIADDACIDT', QCIADDACIDT
        WRITE(IU,FMT_BOOL) 'QCIDOBACKALL', QCIDOBACKALL
        WRITE(IU,FMT_BOOL) 'OPTIMISEAFTERADDITION', OPTIMISEAFTERADDITION
        WRITE(IU,FMT_INT)  'NMINAFTERADD', NMINAFTERADD
        
        
        !============================================================================
        ! Backbone settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Backbone settings ---'
        WRITE(IU,FMT_BOOL) 'QCIDOBACK', QCIDOBACK
        IF (ALLOCATED(ISBBATOM).AND.DEBUG) THEN
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
        
        !============================================================================
        ! Penalty functions
        !============================================================================

        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Repulsion settings ---'
        WRITE(IU,FMT_REAL) 'INTMINFAC', INTMINFAC
        
        !============================================================================
        ! Repulsion settings
        !============================================================================
        
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Repulsion settings ---'
        WRITE(IU,FMT_REAL) 'K_REP', K_REP
        WRITE(IU,FMT_REAL) 'QCIREPCUT', QCIREPCUT
        WRITE(IU,FMT_REAL) 'CHECKREPCUTOFF', CHECKREPCUTOFF
        WRITE(IU,FMT_INT)  'CHECKREPINTERVAL', CHECKREPINTERVAL
        WRITE(IU,FMT_INT)  'QCIINTREPMINSEP', QCIINTREPMINSEP
        WRITE(IU,FMT_INT)  'CHECKREPINTERVAL', CHECKREPINTERVAL
        
        !============================================================================
        ! Constraint settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Constraint settings ---'
        
        WRITE(IU,FMT_REAL) 'K_CONT', K_CONST
        WRITE(IU,FMT_BOOL) 'CHECKINTMINCONSTR', CHECKCONINT
        WRITE(IU,FMT_REAL) 'CONCUTABS', CONCUTABS

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
        ! Linear atom / group settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Linear groups ---'
        WRITE(IU,FMT_BOOL) 'USELINGROUPS', USELINGROUPS
        
        !============================================================================
        ! Image control
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Image control ---'
        WRITE(IU,FMT_BOOL) 'USEIMAGEDENSITY', USEIMAGEDENSITY
        WRITE(IU,FMT_REAL) 'E2E_DIST', E2E_DIST
        WRITE(IU,FMT_REAL) 'IMAGEDENSITY', IMAGEDENSITY
        WRITE(IU,FMT_REAL) 'IMSEPMAX', IMSEPMAX
        WRITE(IU,FMT_REAL) 'IMSEPMIN', IMSEPMIN
        WRITE(IU,FMT_INT)  'QCIIMAGECHECK', QCIIMAGECHECK

        !============================================================================
        ! Permutational settings
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Permutations & chirality ---'
        WRITE(IU,FMT_INT)  'QCIPERMCHECKINT', QCIPERMCHECKINT
        WRITE(IU,FMT_BOOL) 'QCIPERMT', QCIPERMT
        WRITE(IU,FMT_REAL) 'QCIPERMCUT', QCIPERMCUT
        WRITE(IU,FMT_REAL) 'ORBITTOL', ORBITTOL

        WRITE(IU,FMT_BOOL) 'CHECKCHIRAL', CHECKCHIRAL
      
        !============================================================================
        ! Energy minimisation options
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Energy minimisation options ---'
        WRITE(IU,FMT_REAL) 'DGUESS', DGUESS
        WRITE(IU,FMT_INT)  'MUPDATE', MUPDATE
        WRITE(IU,FMT_REAL) 'MAXQCIBFGS', MAXQCIBFGS
        WRITE(IU,FMT_REAL) 'MAXGRADCOMP', MAXGRADCOMP

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
        ! Output control
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Output control ---' 

        WRITE(IU,FMT_INT)  'DUMPQCIXYZFRQS', DUMPQCIXYZFRQS
        WRITE(IU,FMT_BOOL) 'DUMPQCIXYZ', DUMPQCIXYZ

        
        !============================================================================
        ! Experimental / In development
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_BOOL) 'DETECTBBCROSSING', DETECTBBCROSSING
        WRITE(IU,FMT_INT)  'CHECKCROSSFREQ', CHECKCROSSFREQ

        
        !============================================================================
        ! Options which will be potentially removed in future 
        !============================================================================
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '--- Options which may be removed in future ---'
        WRITE(IU,FMT_BOOL) 'QCILINEART', QCILINEART
        IF (ALLOCATED(INLINLIST)) THEN
            N = SIZE(INLINLIST)
            NTRUE = COUNT(INLINLIST)
            WRITE(IU,FMT_INT) 'INLINLIST (size)', N
            WRITE(IU,FMT_INT) 'INLINLIST (.TRUE. count)', NTRUE
            !IF (N <= MAX_PRINT_ARRAY) THEN
            !    WRITE(IU,'(3X,A45,1X,"=",1X,*(L1,1X))') 'INLINLIST', (INLINLIST(I), I=1,N)
            !END IF
        ELSE
            WRITE(IU,FMT_CHAR) 'INLINLIST', 'not allocated'
        END IF
        WRITE(IU,FMT_BOOL) 'LINEARBBT', LINEARBBT
       
              
        WRITE(IU,FMT_SEC) ''
        WRITE(IU,FMT_SEC) '============================================================'
        WRITE(IU,FMT_SEC) '        END OF EFFECTIVE CONFIGURATION'
        WRITE(IU,FMT_SEC) '============================================================'


    END SUBROUTINE WRITE_QCI_KEYS

    SUBROUTINE WRITE_IMAGE_DIST(FILENAME)
             
        USE QCIKEYS, ONLY: NIMAGES
        USE INTERPOLATION_KEYS, ONLY: IMAGE_DIST

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        
        INTEGER :: I, UNIT_OUT

        ! OPEN FILE FOR WRITING
        OPEN(NEWUNIT=UNIT_OUT, FILE=FILENAME, STATUS='REPLACE', ACTION='WRITE')

        ! WRITE GNUPLOT SCRIPT HEADER
        WRITE(UNIT_OUT, '(A)') '#!/usr/bin/gnuplot'
        WRITE(UNIT_OUT, '(A)') 'set encoding iso_8859_1'
        WRITE(UNIT_OUT, '(A)') 'set xlabel "Image" font "Halvetica ,12" '
        WRITE(UNIT_OUT, '(A)') 'set ylabel "Distance-per-atom ({/Helverica=12 \305} )" font "Halvetica ,12"'
        WRITE(UNIT_OUT, '(A)') 'set grid'
        WRITE(UNIT_OUT, '(A)') 'set nokey'
        WRITE(UNIT_OUT, '(A)') ''
        WRITE(UNIT_OUT, '(A)') 'plot "-"  u 1:2 w p pt 6 '
        WRITE(UNIT_OUT, '(A)') '$DATA'

        ! WRITE DATA SECTION
        DO I = 1, NIMAGES+1
            WRITE(UNIT_OUT, '(I0, A, F10.4)') I, ' ', IMAGE_DIST(I)
        END DO

        ! WRITE GNUPLOT END MARKER
        WRITE(UNIT_OUT, '(A)') '$DATA'
        WRITE(UNIT_OUT, '(A)') ''

        ! CLOSE FILE
        CLOSE(UNIT_OUT)

        PRINT *, 'GNUPLOT SCRIPT GENERATED: ', TRIM(FILENAME)
  
    END SUBROUTINE WRITE_IMAGE_DIST

    SUBROUTINE WRITE_IMAGE_E(FILENAME, EEE)
        
        USE QCIKEYS, ONLY: NIMAGES
        
        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        REAL(KIND = REAL64) :: EEE(NIMAGES+2)

        INTEGER :: I, UNIT_OUT

        ! OPEN FILE FOR WRITING
        OPEN(NEWUNIT=UNIT_OUT, FILE=FILENAME, STATUS='REPLACE', ACTION='WRITE')

        ! WRITE GNUPLOT SCRIPT HEADER
        WRITE(UNIT_OUT, '(A)') '#!/usr/bin/gnuplot'
        WRITE(UNIT_OUT, '(A)') 'set encoding iso_8859_1'
        WRITE(UNIT_OUT, '(A)') 'set xlabel "Image" font "Halvetica ,12" '
        WRITE(UNIT_OUT, '(A)') 'set ylabel "Energy"  font "Halvetica ,12"'
        WRITE(UNIT_OUT, '(A)') 'set grid'
        WRITE(UNIT_OUT, '(A)') 'set nokey'
        WRITE(UNIT_OUT, '(A)') ''
        WRITE(UNIT_OUT, '(A)') 'plot "-"  u 1:2 w p pt 6 '
        WRITE(UNIT_OUT, '(A)') '$DATA'

        ! WRITE DATA SECTION
        DO I = 1, NIMAGES+2
            WRITE(UNIT_OUT, '(I0, A, F10.4)') I, ' ', EEE(I)
        END DO

        ! WRITE GNUPLOT END MARKER
        WRITE(UNIT_OUT, '(A)') '$DATA'
        WRITE(UNIT_OUT, '(A)') ''

        ! CLOSE FILE
        CLOSE(UNIT_OUT)

        PRINT *, 'GNUPLOT SCRIPT GENERATED: ', TRIM(FILENAME)
  
    END SUBROUTINE WRITE_IMAGE_E


END MODULE OUT_PRINT

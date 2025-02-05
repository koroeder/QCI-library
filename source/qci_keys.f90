MODULE QCIKEYS
   USE QCIPREC
   IMPLICIT NONE

   LOGICAL :: DEBUG = .FALSE.

   INTEGER :: NATOMS = 0

   INTEGER :: NIMAGES = 0

   INTEGER :: MAXITER = 1500

   INTEGER :: MAXINTIMAGE = 500

   LOGICAL :: QCIAMBERT = .FALSE.
   LOGICAL :: QCIHIRET = .FALSE.
   LOGICAL :: QCISBMT = .FALSE.
   LOGICAL :: QCIGEOMT = .FALSE.

   LOGICAL :: QCIDOBACK = .FALSE.    ! do backbone first
   LOGICAL, ALLOCATABLE :: ISBBATOM(:)
   INTEGER :: NBACKBONE = 0

   CHARACTER(5), ALLOCATABLE :: NAMES(:) 

   LOGICAL :: QCIADDACIDT = .FALSE.
   LOGICAL :: QCIDOBACKALL = .FALSE.

   INTEGER :: INTADDATOM = 1

   LOGICAL :: QCITRILATERATION = .FALSE.

   LOGICAL :: CHECKCONINT = .FALSE. !which congrad should be used?

   LOGICAL :: QCIREADGUESS = .FALSE. ! read guess for band
   CHARACTER(20) :: GUESSFILE = "int_guess.xyz"
   
   LOGICAL :: USEIMAGEDENSITY = .FALSE. ! base number of images on interpolation density
   REAL(KIND=REAL64) :: E2E_DIST = 0.D0 ! endpoint to endpoint distance after alignment
   REAL(KIND=REAL64) :: IMAGEDENSITY = 0.5 ! Image density per unit distance
   REAL(KIND = REAL64) :: IMSEPMAX=HUGE(1.0D0)
   REAL(KIND = REAL64) :: IMSEPMIN=-1.0D0

   ! using linearlist
   LOGICAL :: QCILINEART = .FALSE.
   ! logical list if atom is linear
   LOGICAL, ALLOCATABLE :: INLINLIST(:)
   !linear nterpolation for backbone
   LOGICAL :: LINEARBBT = .FALSE.

   ! Atom to residue mapping
   INTEGER, ALLOCATABLE :: ATOMS2RES(:)

   !frozen atoms 
   LOGICAL :: QCIFREEZET = .FALSE. ! Shall some atoms be frozen?
   INTEGER :: NQCIFROZEN = 0 ! total number of frozen atoms
   INTEGER :: NMINUNFROZEN = 0 ! minimum number unfrozen 
   LOGICAL, ALLOCATABLE :: QCIFROZEN(:)  ! frozen atoms in interpolation
   LOGICAL, ALLOCATABLE :: FREEZE(:) ! input of atoms to be frozen
   CHARACTER(LEN=30) :: FREEZEFILE = 'qci_frozen.dat' ! input file with atoms to be frozen -> used to populate FREEZE
   REAL(KIND=REAL64) :: QCIFREEZETOL = 1.0D-3 ! distance tolerance for atoms to be frozen

   INTEGER :: QCIIMAGECHECK = 10

   REAL(KIND=REAL64) :: INTMINFAC=1.0D0 !Scaling factor for internal minima

   !repulsions
   REAL(KIND=REAL64) :: QCIREPCUT = 1.0D-3
   REAL(KIND=REAL64) :: QCICONSTRREP=100.0D0
   INTEGER ::  QCIINTREPMINSEP=20 !Minimum separation in atom index for internal minimum check in repulsion

   !QUERY: what exactly is this
   REAL(KIND=REAL64) :: INTCONSTRAINTDEL=10.0D0


   !spring constants and adjustment
   !TODO: add initialisation and setting in qci setup
   LOGICAL :: QCISPRINGACTIVET = .TRUE. !QUERY: what is this?
   LOGICAL :: QCIADJUSTKT = .FALSE. ! adjust spring constant
   INTEGER :: QCIADJUSTKFRQ = 0
   REAL(KIND=REAL64) :: QCIADJUSTKTOL = 10.0D0
   REAL(KIND=REAL64) :: QCIAVDEV = 0.0D0
   REAL(KIND=REAL64) :: KINT = 1.0D0
   REAL(KIND=REAL64) :: KINTSCALED = 1.0D0
   REAL(KIND=REAL64) :: QCIKINTMIN = 1.0D-2
   REAL(KIND=REAL64) :: QCIKINTMAX = 1.0D2
   REAL(KIND=REAL64) :: QCIADJUSTKFRAC = 1.05D0

   !maximum gradient component
   REAL(KIND=REAL64) :: MAXGRADCOMP = -1.0

   !settings for step taking
   REAL(KIND = REAL64) :: DGUESS = 1.0D-3
   INTEGER :: MUPDATE = 4
   REAL(KIND = REAL64) :: MAXQCIBFGS = 0.2D0

   REAL(KIND=REAL64) :: COLDFUSIONLIMIT = -1.0D15

   !maximum constraint E - used for convergence
   REAL(KIND=REAL64) :: MAXCONE = 0.01D0
   !convergence for RMS
   REAL(KIND=REAL64) :: QCIRMSTOL=0.01D0

   REAL(KIND=REAL64) :: MAXERISE = 1.0D100

   !QCI resetting
   LOGICAL :: QCIRESET = .TRUE.
   INTEGER :: QCIRESETINT1 = 10

   INTEGER :: CHECKREPINTERVAL = 1

   !checking chirality
   LOGICAL :: CHECKCHIRAL = .FALSE.

   !using groups of atoms in AMBER
   LOGICAL :: QCIUSEGROUPS = .FALSE.

   INTEGER :: DUMPQCIXYZFRQS = 100
   LOGICAL :: DUMPQCIXYZ = .FALSE.


   !permutational stuff
   INTEGER :: QCIPERMCHECKINT = 100
   LOGICAL :: QCIPERMT = .FALSE.
   REAL(KIND=REAL64) :: QCIPERMCUT = 0.8D0

   ! for myorient
   REAL(KIND=REAL64) :: ORBITTOL = 0.3D0

   ! for future extensions if needed
   REAL(KIND=REAL64) :: BOXLX = 0.0D0, BOXLY = 0.0D0, BOXLZ = 0.0D0
   LOGICAL :: BULKT = .FALSE.
   LOGICAL :: TWOD = .FALSE.
   LOGICAL :: RIGID = .FALSE.
   LOGICAL :: STOCKT = .FALSE.

   CONTAINS
      SUBROUTINE ALLOC_QCIKEYS()
         IMPLICIT NONE
         CALL DEALLOC_QCIKEYS()
         ALLOCATE(NAMES(NATOMS))
      END SUBROUTINE ALLOC_QCIKEYS

      SUBROUTINE DEALLOC_QCIKEYS()
         IMPLICIT NONE
         IF (ALLOCATED(NAMES)) DEALLOCATE(NAMES)
      END SUBROUTINE DEALLOC_QCIKEYS

END MODULE QCIKEYS
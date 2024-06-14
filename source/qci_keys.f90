MODULE QCIKEYS
   USE QCIPREC
   IMPLICIT NONE

   LOGICAL :: DEBUG = .FALSE.

   LOGICAL :: STOCKT

   INTEGER :: NATOMS = 0

   INTEGER :: NIMAGES

   LOGICAL :: QCIAMBERT = .FALSE.
   LOGICAL :: QCIHIRET = .FALSE.
   LOGICAL :: QCISBMT = .FALSE.

   LOGICAL :: QCIDOBACK = .FALSE.    ! do backbone first
   LOGICAL, ALLOCATABLE :: ISBBATOM(:)
   INTEGER :: NBACKBONE = 0

   LOGICAL :: QCIADDACIDT = .FALSE.
   LOGICAL :: QCIDOBACKALL = .FALSE.

   LOGICAL :: QCITRILATERATION = .FALSE.

   LOGICAL :: CHECKCONINT = .FALSE. !which congrad should be used?

   LOGICAL :: QCIREADGUESS = .FALSE. ! read guess for band
   CHARACTER(20) :: GUESSFILE = "int_guess.xyz"
   
   LOGICAL :: USEIMAGEDENSITY = .FALSE. ! base number of images on interpolation density
   REAL(KIND=REAL64) :: E2E_DIST = 0.D0 ! endpoint to endpoint distance after alignment
   REAL(KIND=REAL64) :: IMAGEDENSITY = 0.5 ! Image density per unit distance
   INTEGER :: MAXNIMAGES = 250 ! maximum number of images 
   REAL(KIND = REAL64) :: IMSEPMAX=HUGE(1.0D0)
   REAL(KIND = REAL64) :: IMSEPMIN=-1.0D0

   ! using linearlist
   LOGICAL :: QCILINEART = .FALSE.
   ! logical list if atom is linear
   LOGICAL, ALLOCATABLE :: INLINLIST(:)

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
   REAL(KIND=REAL64) :: QCIADJUSTKTOL = 0.0D0
   REAL(KIND=REAL64) :: QCIAVDEV = 0.0D0
   REAL(KIND=REAL64) :: KINT
   REAL(KIND=REAL64) :: KINTSCALED
   REAL(KIND=REAL64) :: QCIKINTMIN, QCIKINTMAX
   REAL(KIND=REAL64) :: QCIADJUSTKFRAC

   !maximum gradient component
   REAL(KIND=REAL64) :: MAXGRADCOMP = -1.0

   !settings for step taking
   REAL(KIND = REAL64) :: DGUESS = 1.0D-3
   INTEGER :: MUPDATE = 4
   REAL(KIND = REAL64) :: MAXQCIBFGS = 0.2D0


   !maximum constraint E - used for convergence
   REAL(KIND=REAL64) :: MAXCONE = 0.01D0
   !convergence for RMS
   REAL(KIND=REAL64) :: QCIRMSTOL=0.01D0

   !QCI resetting
   LOGICAL :: QCIRESET = .TRUE.
   INTEGER :: QCIRESETINT1 = 10

   !checking chirality
   LOGICAL :: CHECKCHIRAL = .FALSE.

   !permutational stuff
   INTEGER :: QCIPERMCHECKINT = 100
   LOGICAL :: QCIPERMCHECK = .FALSE.
   LOGICAL :: QCILPERMDIST = .FALSE.
   LOGICAL :: QCIPERMT = .FALSE.
   REAL(KIND=REAL64) :: QCIPERMCUT = 0.8D0

   ! for myorient
   REAL(KIND=REAL64) :: ORBITTOL = 0.3D0
   LOGICAL :: PULLT = .FALSE.
END MODULE QCIKEYS
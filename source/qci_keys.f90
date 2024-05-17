MODULE QCIKEYS
   USE QCIPREC
   IMPLICIT NONE

   LOGICAL :: DEBUG = .FALSE.

   LOGICAL :: STOCKT

   INTEGER :: NATOMS = 0

   LOGICAL :: QCIAMBERT = .FALSE.
   LOGICAL :: QCIHIRET = .FALSE.
   LOGICAL :: QCISBT = .FALSE.

   LOGICAL :: QCIDOBACK = .FALSE.    ! do backbone first
   LOGICAL, ALLOCATABLE :: ISBBATOM(:)
   INTEGER :: NBACKBONE = 0

   LOGICAL :: QCIREADGUESS = .FALSE. ! read guess for band
   CHARACTER(20) :: GUESSFILE = "int_guess.xyz"
   
   LOGICAL :: USEIMAGEDENSITY = .FALSE. ! base number of images on interpolation density
   REAL(KIND=REAL64) :: E2E_DIST = 0.D0 ! endpoint to endpoint distance after alignment
   REAL(KIND=REAL64) :: IMAGEDENSITY = 0.5 ! Image density per unit distance
   INTEGER :: MAXINTIMAGE = 250 ! maximum number of images 

   !frozen atoms 
   LOGICAL :: QCIFREEZET = .FALSE. ! Shall some atoms be frozen?
   INTEGER :: NCQIFROZEN = 0 ! total number of frozen atoms
   INTEGER :: NMINUNFROZEN = 0 ! minimum number unfrozen 
   LOGICAL, ALLOCATABLE :: QCIFROZEN(:)  ! frozen atoms in interpolation
   LOGICAL, ALLOCATABLE :: FREEZE(:) ! input of atoms to be frozen
   CHARACTER(LEN=30) :: FREEZEFILE = 'qci_frozen.dat' ! input file with atoms to be frozen -> used to populate FREEZE
   REAL(KIND=REAL64) :: QCIFREEZETOL = 1.0D-3 ! distance tolerance for atoms to be frozen

   !constraint
   REAL(KIND=REAL64) :: QCICONSEP = 15
   REAL(KIND=REAL64) :: QCICONSTRAINTTOL = 0.1D0
   REAL(KIND=REAL64) :: QCICONCUT = 0.1D0   

   CHARACTER(LEN=30) :: GEOMFILE = "congeom.dat"
   INTEGER :: NCONGEOM = 0
   REAL(KIND=REAL64), ALLOCATABLE :: CONGEOM(:,:) 
   CHARACTER(LEN=30) :: CONSTRFILE = "constraintfile"

   !repulsions
   REAL(KIND=REAL64) :: QCIREPCUT = 1.0D-3

   !spring constants and adjustment
   !TODO: add initialisation and setting in qci setup
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

   !maximum constraint E - used for convergence
   REAL(KIND=REAL64) :: MAXCONE = 0.01D0
   !convergence for RMS
   REAL(KIND=REAL64) :: QCIRMSTOL=0.01D0

   !cut off for constraints
   LOGICAL :: CONCUTABSINC = .FALSE.
   REAL(KIND=REAL64) :: CONCUTABS = 0.7D0
   REAL(KIND=REAL64) :: CONCUTABS2 = 0.05D0

   !QCI resetting
   LOGICAL :: QCIRESET = .TRUE.
   INTEGER :: QCIRESETINT1 = 10


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
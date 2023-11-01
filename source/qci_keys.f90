MODULE QCIKEYS
   USE QCIPREC
   IMPLICIT NONE
   LOGICAL :: STOCKT

   INTEGER :: NATOMS = 0

   LOGICAL :: QCIAMBERT = .FALSE.
   LOGICAL :: QCIHIRET = .FALSE.
   LOGICAL :: QCISBT = .FALSE.

   LOGICAL :: QCIDOBACK = .FALSE.    ! do backbone first
   LOGICAL, ALLOCATABLE :: ISBBATOM(:)
   INTEGER :: NBACKBONE = 0

   LOGICAL :: QCIREADGUESS = .FALSE. ! read guess for band
   LOGICAL :: QCIRESTART = .FALSE.   ! restart QCI interpolation

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

   !repulsions
   REAL(KIND=REAL64) :: QCIREPCUT = 1.0D-3

   CHARACTER(LEN=30) :: GEOMFILE = "congeom.dat"
   INTEGER :: NCONGEOM = 0
   REAL(KIND=REAL64), ALLOCATABLE :: CONGEOM(:,:) 
   CHARACTER(LEN=30) :: CONSTRFILE = "constraintfile"

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
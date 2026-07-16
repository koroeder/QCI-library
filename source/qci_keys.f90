MODULE QCIKEYS
   USE QCIPREC
   IMPLICIT NONE

   LOGICAL :: DEBUG = .FALSE. !>Get extended output

   !---------------------General interpolation keys-------------------!

   INTEGER :: NATOMS = 0
   CHARACTER(5), ALLOCATABLE :: NAMES(:) 
   INTEGER, ALLOCATABLE :: ATOMS2RES(:) !< Atom to residue mapping
   INTEGER :: MAXITER = 1500

   LOGICAL :: QCIREADGUESS = .FALSE. !< read guess for band
   CHARACTER(20) :: GUESSFILE = "int_guess.xyz"
   
   !----------------------Topology keys -------------------------------!  

   LOGICAL :: QCIAMBERT = .FALSE.
   LOGICAL :: QCIHIRET = .FALSE.
   LOGICAL :: QCISBMT = .FALSE.
   LOGICAL :: QCIGEOMT = .FALSE.

   !------------------------- Amber -----------------------------------!
   
   !> using basepair detection
   LOGICAL :: BASEPAIRDETECTION = .FALSE.
   !> using groups of atoms in AMBER
   LOGICAL :: QCIUSEGROUPS = .FALSE.
   
   !----------------------Atom adding keys ----------------------------!
   
   !> use internal coordinates for local interpolation
   LOGICAL :: USEINTERNALST = .FALSE.
   !> use four atom basis for local interpolation
   LOGICAL :: USEFOURATOMST = .FALSE.
   !> use trilateration for interpolation
   LOGICAL :: QCITRILATERATION = .FALSE.
   !> use trilateration2 for interpolation
   LOGICAL :: QCITRILATERATION2 = .FALSE.

   LOGICAL :: QCIDOBACK = .FALSE.    !< do backbone first
   LOGICAL, ALLOCATABLE :: ISBBATOM(:)
   INTEGER :: NBACKBONE = 0

   !> Add entire residue at the time
   LOGICAL :: QCIADDACIDT = .FALSE.
   LOGICAL :: QCIDOBACKALL = .FALSE.

   !> use minimisations directly after adding atoms
   LOGICAL :: OPTIMISEAFTERADDITION = .FALSE.
   INTEGER :: NMINAFTERADD = 5

   !---------------------Penalty functions keys ----------------------!

   LOGICAL :: CHECKCONINT = .FALSE. !<which congrad should be used? TRUE=CONGRAD2, FALSE=CONGRAD1
   REAL(KIND=REAL64) :: INTMINFAC=1.0D0  !< Scaling factor for internal minima

   !--------------------- Repulsions ---------------------------------!
   
   REAL(KIND=REAL64) :: K_REP=0.2D0   !<scaling constant for repuslion
   REAL(KIND=REAL64) :: QCIREPCUT = 1.0D-3   !< Minimum repulsion cutoff
   INTEGER ::  QCIINTREPMINSEP=4  !< Minimum separation in atom index for internal minimum check in repulsion. Changed default value from 20 to 4 M.H
   REAL(KIND = REAL64) :: CHECKREPCUTOFF=1.25D0   !< factor used for checking repulsion neighbourhood 
   INTEGER :: CHECKREPINTERVAL = 1
   
   !--------------------- Constraints --------------------------------!
   
   !>scaling constant for constraint
   REAL(KIND=REAL64) :: K_CONST=1.0D0 
   !> Should we calculate active-inactive constraints
   LOGICAL :: USECONACTINACT = .FALSE. 
   !> conactinact settings - scale active-inactive atom interaction
   REAL(KIND = REAL64) :: CONACTINACT = 0.2D0
   
   !--------------------- Dihedrals ----------------------------------!

   !> use dihedral cosntraints for chiral atoms and planarity
   LOGICAL :: USEDIHEDRALCONST = .FALSE.
   !> tolerance for dihedral angle difference between start and final images. Preset 10 deg
   REAL(KIND=REAL64) :: DIHDIFTOL = 0.17453292519943295

   !-------------- Spring constants and adjustment -------------------!

   REAL(KIND=REAL64) :: SPRING_GRAD_CONV = 1.5D0
   LOGICAL :: QCISPRINGACTIVET = .TRUE. !< Should we apply spring only to the active atoms? (not in set-up)
   LOGICAL :: QCIADJUSTKT = .FALSE.     !< adjust spring constant
   INTEGER :: QCIADJUSTKFRQ = 0
   
   REAL(KIND=REAL64) :: QCIADJUSTKTOL = 10.0D0 !< tolerance for image spacing deviation (in percent)
   REAL(KIND=REAL64) :: QCIAVDEV = 0.0D0    !< Calculated deviation
   REAL(KIND=REAL64) :: KINT = 1.0D0        !< spring constant
   REAL(KIND=REAL64) :: KINTSCALED = 1.0D0  !< QUESTION Scaling for when we adjust spring constant?
   REAL(KIND=REAL64) :: QCIKINTMIN = 1.0D-2
   REAL(KIND=REAL64) :: QCIKINTMAX = 1.0D2
   REAL(KIND=REAL64) :: QCIADJUSTKFRAC = 1.05D0
   
   LOGICAL :: USEIMAGEDENSITY = .FALSE. !< base number of images on interpolation density
   REAL(KIND=REAL64) :: E2E_DIST = 0.D0 !< endpoint to endpoint distance after alignment
   REAL(KIND=REAL64) :: IMAGEDENSITY = 0.5D0 !< Image density per unit distance (unused atm)
   REAL(KIND = REAL64) :: IMSEPMAX=HUGE(1.0D0)
   REAL(KIND = REAL64) :: IMSEPMIN=-1.0D0

   !-----------------------Linear groups -----------------------------!
  
   LOGICAL :: USELINGROUPS = .FALSE.  !< use linear groups
  
   !----------------------- Image control ----------------------------!

   INTEGER :: NIMAGES = 0
   INTEGER :: MAXINTIMAGE = 500
   INTEGER :: QCIIMAGECHECK = 10
   
   !--------------------- Permutations and chirality -----------------!

   INTEGER :: QCIPERMCHECKINT = 100
   LOGICAL :: QCIPERMT = .FALSE.

   !> checking chirality (TRUE for Amber)
   LOGICAL :: CHECKCHIRAL = .FALSE.
  
   REAL(KIND=REAL64) :: QCIPERMCUT = 0.8D0  !< used in LOPERMDIST
   REAL(KIND=REAL64) :: ORBITTOL = 0.3D0  !< for myorient

   !----------------------- L-BFGS -----------------------------------!

   !> maximum gradient component
   REAL(KIND=REAL64) :: MAXGRADCOMP = -1.0D0

   !> settings for step taking
   REAL(KIND = REAL64) :: DGUESS = 1.0D-3
   INTEGER :: MUPDATE = 4
   REAL(KIND = REAL64) :: MAXQCIBFGS = 0.2D0

   REAL(KIND=REAL64) :: COLDFUSIONLIMIT = -1.0D15

   !> Convergence test for energy
   REAL(KIND=REAL64) :: MAXCONE = 0.01D0
   !> convergence test for RMS
   REAL(KIND=REAL64) :: QCIRMSTOL=0.01D0
   REAL(KIND=REAL64) :: MAXERISE = 1.0D100

   !> QCI resetting
   LOGICAL :: QCIRESET = .TRUE.
   INTEGER :: QCIRESETINT1 = 10

   !-------------------- Output control ------------------------------!

   INTEGER :: DUMPQCIXYZFRQS = 100
   LOGICAL :: DUMPQCIXYZ = .FALSE.
  
   !---------------------------In development ------------------------!
   
   LOGICAL :: DETECTBBCROSSING = .FALSE.
   INTEGER :: CHECKCROSSFREQ = 0
   
   !> for future extensions if needed
   REAL(KIND=REAL64) :: BOXLX = 0.0D0, BOXLY = 0.0D0, BOXLZ = 0.0D0
   LOGICAL :: BULKT = .FALSE.
   LOGICAL :: TWOD = .FALSE.
   LOGICAL :: RIGID = .FALSE.
   LOGICAL :: STOCKT = .FALSE.

   !----------------------- Candidates for removal--------------------!
  
   !> using linearlist
   LOGICAL :: QCILINEART = .FALSE.
   !> logical list if atom is linear
   LOGICAL, ALLOCATABLE :: INLINLIST(:)
   !> linear nterpolation for backbone
   LOGICAL :: LINEARBBT = .FALSE.
  
   !> frozen atoms 
   LOGICAL :: QCIFREEZET = .FALSE. !< Shall some atoms be frozen?
   INTEGER :: NQCIFROZEN = 0 !< total number of frozen atoms
   INTEGER :: NMINUNFROZEN = 0 !< minimum number unfrozen 
   LOGICAL, ALLOCATABLE :: QCIFROZEN(:)  !< frozen atoms in interpolation
   LOGICAL, ALLOCATABLE :: FREEZE(:) !< input of atoms to be frozen
   CHARACTER(LEN=30) :: FREEZEFILE = 'qci_frozen.dat' !< input file with atoms to be frozen -> used to populate FREEZE
   REAL(KIND=REAL64) :: QCIFREEZETOL = 1.0D-3 !< distance tolerance for atoms to be frozen

   !REAL(KIND=REAL64) :: INTCONSTRAINTDEL=10.0D0 
   
   !------------------------------------------------------------------!
   !------------------------ END of keys------------------------------!
   !------------------------------------------------------------------!

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
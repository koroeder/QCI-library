MODULE QCI_CONSTRAINT_KEYS
   USE QCIPREC
   IMPLICIT NONE
   !constraint
   INTEGER :: MAXCONUSE = 10 !< Max number of contraints per atom 
     

   !> only used for congeom
   REAL(KIND=REAL64) :: QCICONSTRAINTTOL = 0.1D0
   INTEGER :: QCICONSEP = 15 
   REAL(KIND=REAL64) :: QCICONCUT = 6.0D0 !< maximum distance between two atoms to create a constraint  
   CHARACTER(LEN=30) :: GEOMFILE = "congeom.dat"
   INTEGER :: NCONGEOM = 0
   REAL(KIND=REAL64), ALLOCATABLE :: CONGEOM(:,:) 
   
   CHARACTER(LEN=30) :: CONSTRFILE = "constraintfile"

   INTEGER :: NCONSTRAINT = 0
   INTEGER, ALLOCATABLE :: CONI(:), CONJ(:)
   REAL(KIND = REAL64), ALLOCATABLE :: CONDISTREF(:)
   REAL(KIND = REAL64), ALLOCATABLE :: CONDISTREFLOCAL(:) !<  mean(d_{AB})=(d^1_{AB} + d^M_{AB})/2 
   REAL(KIND = REAL64), ALLOCATABLE :: CONCUT(:), CONCUTLOCAL(:)
   INTEGER, ALLOCATABLE :: NCONPERATOM(:), CONLIST(:,:)
   INTEGER :: MAXCONSTRAINTS

   ! bonds array for linear groups and bond potential
   INTEGER, ALLOCATABLE :: BOND_LIST(:,:) !< BOND_LIST(atom_in_bond_id, bond_id)
   INTEGER, ALLOCATABLE :: N_BONDS_PER_ATOM(:)
   INTEGER, ALLOCATABLE :: BONDS_PER_ATOM_LIST(:,:)
   REAL(KIND=REAL64), ALLOCATABLE :: AVE_BL_LENGTH(:)
   INTEGER :: MAX_BONDS_PER_ATOM = 0
   INTEGER :: NBONDS = 0
  
   ! angles arrays for angle potential 
   INTEGER, ALLOCATABLE :: ANGLE_LIST(:,:) !< ANGLE_LIST(atom_id, angle_id)
   REAL(KIND=REAL64), ALLOCATABLE :: AVE_ANGLE(:) 
   INTEGER :: NANGLES = 0

   !> Should we use absolute value addition to adjust the cutoffs ?
   LOGICAL :: CONCUTABST = .FALSE.
   !> Adjustment cut off for constraints
   REAL(KIND=REAL64) :: CONCUTABS = 0.7D0

   !> Use fraction adjustment for constraints cutoffs? 
   LOGICAL :: CONCUTFRACT = .FALSE.
   !> Fraction for fractional adjustemnt of contraints cutoffs. 
   REAL(KIND=REAL64) :: CONCUTFRAC=0.1D0

   !> Should CONCUTABS be reset? Dynamically adjusted during the interpolation. 
   LOGICAL :: CONCUTABSINC = .FALSE. 
   
   CONTAINS

   SUBROUTINE DEALLOC_BONDS()
    IF (ALLOCATED(BOND_LIST)) DEALLOCATE(BOND_LIST)
    IF (ALLOCATED(N_BONDS_PER_ATOM)) DEALLOCATE(N_BONDS_PER_ATOM)
    IF (ALLOCATED(BONDS_PER_ATOM_LIST)) DEALLOCATE(BONDS_PER_ATOM_LIST)
    IF (ALLOCATED(ANGLE_LIST)) DEALLOCATE(ANGLE_LIST)
    IF (ALLOCATED(AVE_ANGLE)) DEALLOCATE(AVE_ANGLE)
    IF (ALLOCATED(AVE_BL_LENGTH)) DEALLOCATE(AVE_BL_LENGTH)
  END SUBROUTINE DEALLOC_BONDS

 END MODULE QCI_CONSTRAINT_KEYS
 
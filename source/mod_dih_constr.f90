MODULE DIHEDRAL_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   ! set value for pi
   REAL(KIND = REAL64), PARAMETER :: PI = 3.141592653589793D0
   ! number of dihedrals to be constraint
   INTEGER :: NDIH
   ! atoms in dihedrals (by id)
   INTEGER, ALLOCATABLE :: DIHEDRALS(:,:)
   ! reference dihedrals
   REAL(KIND = REAL64), ALLOCATABLE :: REFDIH(:), REFDIHS(:), REFDIHF(:)
   ! dihedral multiplicity
   INTEGER, ALLOCATABLE :: DIHMUL(:)
   ! regularisation of phiref
   REAL(KIND = REAL64), ALLOCATABLE :: S0(:), C0(:)
   ! activated all dihedrals?
   LOGICAL :: ALLDIHACTIVE = .FALSE.
   ! activated dihedrals?
   LOGICAL, ALLOCATABLE :: DIHACTIVE(:)
   !> spring restraint constant
   REAL(KIND = REAL64) :: KDIH !QUESTION should this be in qci keys?

   CONTAINS
      SUBROUTINE ALLOC_DIHVARS()
         CALL DEALLOC_DIHVARS
         ALLOCATE(DIHEDRALS(NDIH,4))
         ALLOCATE(REFDIH(NDIH), REFDIHS(NDIH), REFDIHF(NDIH))
         ALLOCATE(DIHACTIVE(NDIH))
         ALLOCATE(S0(NDIH),C0(NDIH))
         ALLOCATE(DIHMUL(NDIH))
      END SUBROUTINE ALLOC_DIHVARS

      SUBROUTINE DEALLOC_DIHVARS()
         IF (ALLOCATED(DIHEDRALS)) DEALLOCATE(DIHEDRALS)
         IF (ALLOCATED(REFDIH)) DEALLOCATE(REFDIH)
         IF (ALLOCATED(REFDIHS)) DEALLOCATE(REFDIHS)
         IF (ALLOCATED(REFDIHF)) DEALLOCATE(REFDIHF)
         IF (ALLOCATED(DIHACTIVE)) DEALLOCATE(DIHACTIVE)
         IF (ALLOCATED(S0)) DEALLOCATE(S0)
         IF (ALLOCATED(C0)) DEALLOCATE(C0)
         IF (ALLOCATED(DIHMUL)) DEALLOCATE(DIHMUL)
      END SUBROUTINE DEALLOC_DIHVARS

      SUBROUTINE SETUP_DIH_CONSTR()
         USE QCIKEYS, ONLY: NATOMS, DIHDIFTOL
         USE INTERPOLATION_KEYS, ONLY: XSTART, XFINAL
         USE AMBER_CONSTRAINTS, ONLY: RESNAMES, NRES, GET_ATOMID
         USE CHIRALITY, ONLY: NCHIRAL, CHIR_INFO

         REAL(KIND = REAL64), PARAMETER :: EPS3 = 1.0D-3
         REAL(KIND = REAL64), PARAMETER :: EPS6 = 1.0D-6
         INTEGER :: REFATOMS(NATOMS,4)
         INTEGER :: NCONS, AT1, AT2, AT3, AT4, I, J
         REAL(KIND=REAL64) :: THISDIHS, THISDIHF !phi_0 for start and final images
         REAL(KIND=REAL64) :: THISDIH, THISS0, THISC0

         NCONS = 0
         DO I=1,NRES
            IF ((RESNAMES(I).EQ."PHE").OR.(RESNAMES(I).EQ."CPHE").OR.(RESNAMES(I).EQ."NPHE")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CG",I,AT1)
               CALL GET_ATOMID("CD1",I,AT2)
               CALL GET_ATOMID("CE1",I,AT3)
               CALL GET_ATOMID("CZ",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE1",I,AT1)
               CALL GET_ATOMID("CZ",I,AT2)
               CALL GET_ATOMID("CE2",I,AT3)
               CALL GET_ATOMID("CD2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4     
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE2",I,AT1)
               CALL GET_ATOMID("CD2",I,AT2)
               CALL GET_ATOMID("CG",I,AT3)
               CALL GET_ATOMID("CD1",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4              
            ELSE IF ((RESNAMES(I).EQ."TYR").OR.(RESNAMES(I).EQ."CTYR").OR.(RESNAMES(I).EQ."NTYR")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CG",I,AT1)
               CALL GET_ATOMID("CD1",I,AT2)
               CALL GET_ATOMID("CE1",I,AT3)
               CALL GET_ATOMID("CZ",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE1",I,AT1)
               CALL GET_ATOMID("CZ",I,AT2)
               CALL GET_ATOMID("CE2",I,AT3)
               CALL GET_ATOMID("CD2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4     
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE2",I,AT1)
               CALL GET_ATOMID("CD2",I,AT2)
               CALL GET_ATOMID("CG",I,AT3)
               CALL GET_ATOMID("CD1",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
            ELSE IF ((RESNAMES(I).EQ."HIS").OR.(RESNAMES(I).EQ."CHIS").OR.(RESNAMES(I).EQ."NHIS")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CG",I,AT1)
               CALL GET_ATOMID("ND1",I,AT2)
               CALL GET_ATOMID("CE1",I,AT3)
               CALL GET_ATOMID("NE2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE1",I,AT1)
               CALL GET_ATOMID("NE2",I,AT2)
               CALL GET_ATOMID("CD2",I,AT3)
               CALL GET_ATOMID("CG",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4                
            ELSE IF ((RESNAMES(I).EQ."TRP").OR.(RESNAMES(I).EQ."CTRP").OR.(RESNAMES(I).EQ."NTRP")) THEN  
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CG",I,AT1)
               CALL GET_ATOMID("CD2",I,AT2)
               CALL GET_ATOMID("CE2",I,AT3)
               CALL GET_ATOMID("CZ2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("NE1",I,AT1)
               CALL GET_ATOMID("CE2",I,AT2)
               CALL GET_ATOMID("CD2",I,AT3)
               CALL GET_ATOMID("CE3",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4     
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE2",I,AT1)
               CALL GET_ATOMID("CZ2",I,AT2)
               CALL GET_ATOMID("CH2",I,AT3)
               CALL GET_ATOMID("CZ3",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               !fourth dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("CE2",I,AT1)
               CALL GET_ATOMID("NE1",I,AT2)
               CALL GET_ATOMID("CD1",I,AT3)
               CALL GET_ATOMID("CG",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
            ELSE IF ((RESNAMES(I).EQ."A").OR.(RESNAMES(I).EQ."A3").OR.(RESNAMES(I).EQ."A5").OR. &
               (RESNAMES(I).EQ."DA").OR.(RESNAMES(I).EQ."DA3").OR.(RESNAMES(I).EQ."DA5").OR. &
               (RESNAMES(I).EQ."FA").OR.(RESNAMES(I).EQ."FA3").OR.(RESNAMES(I).EQ."FA5")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N9",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !WRITE(*,*) "DIHEDRAL check: N9-C4-C5-C6: ", REFATOMS(NCONS,1:4)
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N7",I,AT1)
               CALL GET_ATOMID("C5",I,AT2)
               CALL GET_ATOMID("C4",I,AT3)
               CALL GET_ATOMID("N3",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4    
               !WRITE(*,*) "DIHEDRAL check: N7-C5-C4-N3: ", REFATOMS(NCONS,1:4)
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               !WRITE(*,*) "DIHEDRAL check: C5-C6-N1-C2: ", REFATOMS(NCONS,1:4)
               !fourth dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("N7",I,AT2)
               CALL GET_ATOMID("C8",I,AT3)
               CALL GET_ATOMID("N9",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !WRITE(*,*) "DIHEDRAL check: C5-N7-C8-N9: ", REFATOMS(NCONS,1:4)
            ELSE IF ((RESNAMES(I).EQ."C").OR.(RESNAMES(I).EQ."C3").OR.(RESNAMES(I).EQ."C5").OR. &
                  (RESNAMES(I).EQ."DC").OR.(RESNAMES(I).EQ."DC3").OR.(RESNAMES(I).EQ."DC5").OR. &
                  (RESNAMES(I).EQ."FC").OR.(RESNAMES(I).EQ."FC3").OR.(RESNAMES(I).EQ."FC5")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N1",I,AT1)
               CALL GET_ATOMID("C2",I,AT2)
               CALL GET_ATOMID("N3",I,AT3)
               CALL GET_ATOMID("C4",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N3",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4    
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               !might need fourth dihedral
            ELSE IF ((RESNAMES(I).EQ."G").OR.(RESNAMES(I).EQ."G3").OR.(RESNAMES(I).EQ."G5").OR. &
                     (RESNAMES(I).EQ."DG").OR.(RESNAMES(I).EQ."DG3").OR.(RESNAMES(I).EQ."DG5").OR. &
                     (RESNAMES(I).EQ."FG").OR.(RESNAMES(I).EQ."FG3").OR.(RESNAMES(I).EQ."FG5")) THEN
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N9",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N7",I,AT1)
               CALL GET_ATOMID("C5",I,AT2)
               CALL GET_ATOMID("C4",I,AT3)
               CALL GET_ATOMID("N3",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4    
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               !fourth dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("N7",I,AT2)
               CALL GET_ATOMID("C8",I,AT3)
               CALL GET_ATOMID("N9",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
            ELSE IF ((RESNAMES(I).EQ."U").OR.(RESNAMES(I).EQ."U3").OR.(RESNAMES(I).EQ."U5").OR. &
                     (RESNAMES(I).EQ."DT").OR.(RESNAMES(I).EQ."DT3").OR.(RESNAMES(I).EQ."DT5").OR. & 
                     (RESNAMES(I).EQ."FU").OR.(RESNAMES(I).EQ."FU3").OR.(RESNAMES(I).EQ."FU5").OR. &
                     (RESNAMES(I).EQ."FT").OR.(RESNAMES(I).EQ."FT3").OR.(RESNAMES(I).EQ."FT5")) THEN           
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N1",I,AT1)
               CALL GET_ATOMID("C2",I,AT2)
               CALL GET_ATOMID("N3",I,AT3)
               CALL GET_ATOMID("C4",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N3",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4    
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               !experimental add fourth dihedral
               NCONS = NCONS + 1
               CALL GET_ATOMID("C4",I,AT1)
               CALL GET_ATOMID("N3",I,AT2)
               CALL GET_ATOMID("C2",I,AT3)
               CALL GET_ATOMID("N1",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               
            END IF
         END DO
         
         !assign arrays for dihedrals
         NDIH = NCHIRAL + NCONS
         CALL ALLOC_DIHVARS()
         DIHACTIVE(1:NDIH) = .FALSE.
         REFDIH(1:NDIH) = 0.0D0
         REFDIHS(1:NDIH) = 0.0D0
         REFDIHF(1:NDIH) = 0.0D0
         
         DIHMUL(1:NDIH) = 0 
         
         DO J=1,NCHIRAL
            DIHEDRALS(J,1) = CHIR_INFO(J,1)
            DIHEDRALS(J,2) = CHIR_INFO(J,2)
            DIHEDRALS(J,3) = CHIR_INFO(J,3)
            DIHEDRALS(J,4) = CHIR_INFO(J,4)
            DIHMUL(J) = 1
         WRITE(*,*) "DIHEDRALS-CHIR_INFO:  J=", J, "atoms: ", DIHEDRALS(J,1:4)
         END DO
         DO J=1,NCONS
            DIHEDRALS(NCHIRAL+J,1) = REFATOMS(J,1)
            DIHEDRALS(NCHIRAL+J,2) = REFATOMS(J,2)
            DIHEDRALS(NCHIRAL+J,3) = REFATOMS(J,3)
            DIHEDRALS(NCHIRAL+J,4) = REFATOMS(J,4)            
            DIHMUL(NCHIRAL+J) = 2  
         WRITE(*,*) "DIHEDRALS-:  J=", J, "atoms: ",  DIHEDRALS(NCHIRAL+J,1:4)
         END DO
          
         !DO J=1,NCONS
         !   DIHEDRALS(J,1) = REFATOMS(J,1)
         !   DIHEDRALS(J,2) = REFATOMS(J,2)
         !   DIHEDRALS(J,3) = REFATOMS(J,3)
         !   DIHEDRALS(J,4) = REFATOMS(J,4) 
         !   DIHMUL(J) = 2           
         !WRITE(*,*) "DIHEDRALS-:  J=", J, "atoms: ",  DIHEDRALS(J,1:4)
         !END DO

         ! get reference angles
         DO J=1,NDIH
            AT1 = DIHEDRALS(J,1); AT2 = DIHEDRALS(J,2); AT3 = DIHEDRALS(J,3); AT4 = DIHEDRALS(J,4)
            !reference dihedrals in start and finish images 
            CALL COMPUTE_DIH(NATOMS, XSTART, AT1, AT2, AT3, AT4, THISDIHS)
            CALL COMPUTE_DIH(NATOMS, XFINAL, AT1, AT2, AT3, AT4, THISDIHF)

            REFDIHS(J) = THISDIHS 
            REFDIHF(J) = THISDIHF 
            
            !This should give a bisection of the angle 
            THISDIH = ATAN2(DSIN(THISDIHS)+DSIN(THISDIHF),DCOS(THISDIHS)+DCOS(THISDIHF))
            REFDIH(J) = THISDIH 
            
            ! now compute and store the regularised versions
            IF (DABS(THISDIH-PI).LE.EPS3) THISDIH = SIGN(PI,THISDIH)
            
            !IF (DIHMUL(J).EQ.2) THISDIH = THISDIH + PI
            

            THISC0 = DCOS(THISDIH)
            THISS0 = DSIN(THISDIH)
            IF (DABS(THISC0).LE.EPS6) THISC0 = 0.0D0
            IF (DABS(THISS0).LE.EPS6) THISS0 = 0.0D0

            C0(J) = THISC0
            S0(J) = THISS0
         END DO
      END SUBROUTINE SETUP_DIH_CONSTR

      SUBROUTINE COMPUTE_DIH(NATOMS, X, A, B, C, D, PHI)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, EUC_NORM, DOTP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS) !< all atoms in an image
         INTEGER, INTENT(IN) :: A, B, C, D
         REAL(KIND = REAL64), INTENT(OUT) :: PHI

         REAL(KIND = REAL64) :: RAB(3), RBC(3), RCD(3)
         REAL(KIND = REAL64) :: N1(3), N2(3), NORM1, NORM2, NORMBC, NPROD
         REAL(KIND = REAL64) :: COSPHI, SINPHI

         ! compute the vectors between atoms in order
         RAB(1:3) = X(3*B-2:3*B) - X(3*A-2:3*A)
         RBC(1:3) = X(3*C-2:3*C) - X(3*B-2:3*B)
         RCD(1:3) = X(3*D-2:3*D) - X(3*C-2:3*C)

         ! get normal vectors
         N1(1:3) = CROSS_PROD(RAB,RBC)
         N2(1:3) = CROSS_PROD(RBC,RCD)

         ! get norms
         NORM1 = EUC_NORM(N1)
         NORM2 = EUC_NORM(N2)
         NORMBC = EUC_NORM(RBC)
        
         N1(1:3) = N1(1:3)/NORM1
         N2(1:3) = N2(1:3)/NORM2
         
         NPROD = 1.0D0/(EUC_NORM(CROSS_PROD(N1,N2))) 
         
         !compute cos and sin components and get dihedral
         COSPHI = DOTP(3,N1,N2)     
         SINPHI = DOTP(3,CROSS_PROD(N1,N2),RBC)/NORMBC
        
         !ATAN2(y,x)
         PHI = ATAN2(SINPHI,COSPHI)
         
         !PHI = DSIGN(DACOS(COSPHI), DOTP(3,RBC,CROSS_PROD(N1,N2)))
      END SUBROUTINE COMPUTE_DIH

      SUBROUTINE CHECK_DIH_ACTIVE()
         USE INTERPOLATION_KEYS, ONLY: NACTIVE, ATOMACTIVE
         IMPLICIT NONE
         INTEGER :: I


         ALLDIHACTIVE = .TRUE.
         DO I=1,NDIH
            IF (.NOT.DIHACTIVE(I)) THEN
               IF (ATOMACTIVE(DIHEDRALS(I,1)).AND.ATOMACTIVE(DIHEDRALS(I,2)).AND. &
                   ATOMACTIVE(DIHEDRALS(I,3)).AND.ATOMACTIVE(DIHEDRALS(I,4))) THEN
                  DIHACTIVE(I) = .TRUE.
               ELSE
                  ALLDIHACTIVE = .FALSE.
               END IF
            END IF
         END DO
      END SUBROUTINE CHECK_DIH_ACTIVE

      !> Dihedral function based on cosine potential and dV/d(cos(phi)) derivation  
      SUBROUTINE DIHEDRAL(A, B, C, D, DIHREF, E, FA, FB, FC, FD)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, EUC_NORM, DOTP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: A(3), B(3), C(3), D(3)
         INTEGER, INTENT(IN) :: DIHREF
         REAL(KIND = REAL64), INTENT(OUT) :: E
         REAL(KIND = REAL64), INTENT(OUT) :: FA(3), FB(3), FC(3), FD(3)   

         REAL(KIND = REAL64) :: RAB(3), RBC(3), RCD(3)
         REAL(KIND = REAL64) :: N1(3), N2(3), NORM1, NORM2
         REAL(KIND = REAL64) :: PHI_REG, SINPHI, COSPHI, DF
         REAL(KIND = REAL64) :: DOT12, L1, L2, Z1, Z2, Z12, REGTERM1, CT0, CT1, SREG
         REAL(KIND = REAL64) :: FAB(3), FBC(3), FCD(3)
         REAL(KIND = REAL64) :: FREG1(3), FREG2(3)

         !cutoffs used for regularisation
         REAL(KIND = REAL64), PARAMETER :: EPS9 = 1.0D-9
         REAL(KIND = REAL64), PARAMETER :: EPS6 = 1.0D-6
         REAL(KIND = REAL64), PARAMETER :: EPS3 = 1.0D-3
         REAL(KIND = REAL64), PARAMETER :: P9999 = 0.9999

         ! we use the potential for each constraint of E=KDIH*[1+cos(PHI-PHIREF)]
         ! with KDIH as the force constant, PHI as the current dihedral and PHIREf the reference angle
         ! The normal factor of M in the cosine we drop as M=1 (we only want one minimum in the dihedral potential)
         ! This can be rewritten to
         ! E = Deps*KDIH*[1 + cos(PHI_REG)*C0 + sin(PHI_REG)*S0]
         ! The variables are as follows:
         ! PHI_REG = PI - SIGN(dotp(rBC,(N1xN2))) * arccos(DOTP(N1,N2))
         ! with rBC the vector between the central two atoms and N1 and N2 the two normals at B and C, respectively
         ! C0 and S0 are the cosine and sine of PHIREF_REG, i.e. the regularised version of the reference
         ! These are stored in arrays after the setup.
         ! Deps is a regularisation for small values (REGTERM1)

         ! compute the vectors between atoms in order
         RAB(1:3) = B(1:3) - A(1:3)
         RBC(1:3) = C(1:3) - B(1:3)
         RCD(1:3) = D(1:3) - C(1:3)

         ! get normal vectors
         N1(1:3) = CROSS_PROD(RAB,RBC)
         N2(1:3) = CROSS_PROD(RBC,RCD)  
         ! get norms
         NORM1 = EUC_NORM(N1)
         NORM2 = EUC_NORM(N2)
         DOT12 = DOTP(3,N1,N2)

         ! introduce regularisation step 1: step function to cutoff small values
         L1 = NORM1 + EPS9
         L2 = NORM2 + EPS9
         Z1 = 1.0D0/L1
         Z2 = 1.0D0/L2
         IF (EPS3.GT.L1) Z1 = 0.0D0
         IF (EPS3.GT.L2) Z2 = 0.0D0
         Z12 = Z1*Z2
         REGTERM1 = 0.0D0
         IF (Z12.NE.0.0D0) REGTERM1 = 1.0D0
                  
         ! second step of regularisation: restrict values to correct range (-1 to 1)
         CT0 = MIN(P9999,DOT12*Z12)
         CT1 = MAX(-P9999,CT0)

         ! compute regularised dihedral and the sine and cosine terms
         PHI_REG = PI - DSIGN(DACOS(CT1), DOTP(3,RBC,CROSS_PROD(N1,N2)))
         
         IF (DABS(PHI_REG-PI).LE.EPS3) PHI_REG = SIGN(PI,PHI_REG)
                 
         SINPHI = DSIN(PHI_REG)
         COSPHI = DCOS(PHI_REG)

         IF (DABS(SINPHI).LE.EPS6) SINPHI = 0.0D0
         IF (DABS(COSPHI).LE.EPS6) COSPHI = 0.0D0
     
         E = KDIH*(1+COSPHI*C0(DIHREF)+SINPHI*S0(DIHREF))*REGTERM1

         !use a regularised version of sine
         SREG = SINPHI + SIGN(1.0d-18,SINPHI)
         
         ! First part of gradient calculation
         ! dE/d(cos(phi)) 
         IF (ABS(SREG).LT.EPS6) THEN
         !  ! for small sine values we take the limit of the exact form
            DF = C0(DIHREF)*REGTERM1
         ELSE
            DF = (C0(DIHREF)*SINPHI - S0(DIHREF)*COSPHI)/SREG * REGTERM1
         END IF          
         
         ! Individual first derviatives with respect to the cartesian coords
         ! d(cos(phi))/d(R_xx)
         FREG1(1:3) = N2*Z12 - CT1*N1*Z1**2
         FREG2(1:3) = N1*Z12 - CT1*N2*Z2**2
         
         ! The minus gradient is accounted for in this step!
         ! Now, we have force terms!
         ! F_xx = - dE/dcos * dcos/d(R_xx)
         FAB(1:3) =  DF* CROSS_PROD(FREG1,RBC)
         FBC(1:3) = -DF*(CROSS_PROD(FREG1,RAB) + CROSS_PROD(FREG2, RCD))
         FCD(1:3) =  DF* CROSS_PROD(FREG2,RBC)

         !Note: sanity check - sum of forces on all atoms must be 0
         FA(1:3) =  FAB(1:3) 
         FB(1:3) = -FAB(1:3) + FBC(1:3) 
         FC(1:3) = -FCD(1:3) - FBC(1:3) 
         FD(1:3) =  FCD(1:3)   

      END SUBROUTINE DIHEDRAL

      !use of a harmonic constraint, which is not periodic -> not ideal here
      SUBROUTINE OLD_DIHEDRAL(A, B, C, D, PHIREF, E, FA, FB, FC, FD)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, EUC_NORM, DOTP
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: A(3), B(3), C(3), D(3)
         REAL(KIND = REAL64), INTENT(IN) :: PHIREF
         REAL(KIND = REAL64), INTENT(OUT) :: E
         REAL(KIND = REAL64), INTENT(OUT) :: FA(3), FB(3), FC(3), FD(3)

         REAL(KIND = REAL64) :: RAB(3), RBC(3), RCD(3)
         REAL(KIND = REAL64) :: N1(3), N2(3), NORM1, NORM2, NORMBC, NPROD, FE
         REAL(KIND = REAL64) :: COSPHI, SINPHI, PHI

         REAL(KIND = REAL64) :: DSINDN1(3), DSINDN2(3), DCOSDN1(3), DCOSDN2(3)
         REAL(KIND = REAL64) :: DPHIDN1(3), DPHIDN2(3)

         REAL(KIND = REAL64) :: DN1DA(3,3), DN1DB(3,3), DN1DC(3,3)
         REAL(KIND = REAL64) :: DN2DB(3,3), DN2DC(3,3), DN2DD(3,3)
         
         REAL(KIND = REAL64) :: X0(3), Y0(3), Z0(3)

         ! unit vectors references
         X0 = (/1.0, 0.0, 0.0/)
         Y0 = (/0.0, 1.0, 0.0/)
         Z0 = (/0.0, 0.0, 1.0/)

         ! compute the vectors between atoms in order
         RAB(1:3) = B(1:3) - A(1:3)
         RBC(1:3) = C(1:3) - B(1:3)
         RCD(1:3) = D(1:3) - C(1:3)

         ! get normal vectors
         N1(1:3) = CROSS_PROD(RAB,RBC)
         N2(1:3) = CROSS_PROD(RBC,RCD)

         ! get norms
         NORM1 = EUC_NORM(N1)
         NORM2 = EUC_NORM(N2)
         NORMBC = EUC_NORM(RBC)
         NPROD = 1.0D0/(NORM1*NORM2*NORMBC)
         N1(1:3) = N1(1:3)/NORM1
         N2(1:3) = N2(1:3)/NORM2

         !comput cos and sin components and get dihedral
         COSPHI = DOTP(3,N1,N2)
         SINPHI = DOTP(3,CROSS_PROD(N1,N2),RBC)*NPROD
         PHI = ATAN2(SINPHI,COSPHI)

         E = KDIH*(PHI-PHIREF)**2

         FE = 2.0D0*KDIH*(PHI-PHIREF)
         ! start the derivatives here
         !derivative of cos(phi) and sin(phi) w.r.t. n1 and n2
         DCOSDN1 = (N2 - COSPHI*N1)/NORM1
         DCOSDN2 = (N1 - COSPHI*N2)/NORM2
         DSINDN1 = CROSS_PROD(N2,RBC)*NPROD - (SINPHI/NORM1**2)*N1
         DSINDN2 = CROSS_PROD(RBC,N1)*NPROD - (SINPHI/NORM2**2)*N2

         !derivative of phi w.r.t. n1 and n2
         DPHIDN1 = (DSINDN1*COSPHI - DCOSDN1*SINPHI)/(NORM1*SINPHI)
         DPHIDN2 = (DSINDN2*COSPHI - DCOSDN2*SINPHI)/(NORM2*SINPHI)

         DN1DA(1,1:3) = -CROSS_PROD(X0,RBC)
         DN1DA(2,1:3) = -CROSS_PROD(Y0,RBC)
         DN1DA(3,1:3) = -CROSS_PROD(Z0,RBC)

         DN1DB(1,1:3) = CROSS_PROD(X0,RAB) - CROSS_PROD(X0,RBC)
         DN1DB(2,1:3) = CROSS_PROD(Y0,RAB) - CROSS_PROD(Y0,RBC)
         DN1DB(3,1:3) = CROSS_PROD(Z0,RAB) - CROSS_PROD(Z0,RBC)

         DN1DC(1,1:3) = CROSS_PROD(X0,RAB)
         DN1DC(2,1:3) = CROSS_PROD(Y0,RAB)
         DN1DC(3,1:3) = CROSS_PROD(Z0,RAB)

         DN2DB(1,1:3) = -CROSS_PROD(X0,RCD)
         DN2DB(2,1:3) = -CROSS_PROD(Y0,RCD)
         DN2DB(3,1:3) = -CROSS_PROD(Z0,RCD)

         DN2DC(1,1:3) = CROSS_PROD(X0,RBC) - CROSS_PROD(X0,RCD)
         DN2DC(2,1:3) = CROSS_PROD(Y0,RBC) - CROSS_PROD(Y0,RCD)
         DN2DC(3,1:3) = CROSS_PROD(Z0,RBC) - CROSS_PROD(Z0,RCD)

         DN2DD(1,1:3) = CROSS_PROD(X0,RBC)
         DN2DD(2,1:3) = CROSS_PROD(Y0,RBC)
         DN2DD(3,1:3) = CROSS_PROD(Z0,RBC)      

         FA(1) = DOTP(3,DPHIDN1, DN1DA(1,1:3))
         FA(2) = DOTP(3,DPHIDN1, DN1DA(2,1:3))
         FA(3) = DOTP(3,DPHIDN1, DN1DA(3,1:3))

         FB(1) = DOTP(3,DPHIDN1, DN1DB(1,1:3)) + DOTP(3,DPHIDN2,DN2DB(1,1:3))
         FB(2) = DOTP(3,DPHIDN1, DN1DB(2,1:3)) + DOTP(3,DPHIDN2,DN2DB(2,1:3))
         FB(3) = DOTP(3,DPHIDN1, DN1DB(3,1:3)) + DOTP(3,DPHIDN2,DN2DB(3,1:3))

         FC(1) = DOTP(3,DPHIDN1, DN1DC(1,1:3)) + DOTP(3,DPHIDN2,DN2DC(3,1:3))
         FC(2) = DOTP(3,DPHIDN1, DN1DC(2,1:3)) + DOTP(3,DPHIDN2,DN2DC(3,1:3))
         FC(3) = DOTP(3,DPHIDN1, DN1DC(3,1:3)) + DOTP(3,DPHIDN2,DN2DC(3,1:3))

         FD(1) = DOTP(3,DPHIDN2, DN2DD(1,1:3))
         FD(2) = DOTP(3,DPHIDN2, DN2DD(2,1:3))
         FD(3) = DOTP(3,DPHIDN2, DN2DD(3,1:3))

         FA(1:3) = FE*FA(1:3)
         FB(1:3) = FE*FB(1:3)
         FC(1:3) = FE*FC(1:3)
         FD(1:3) = FE*FD(1:3)
      END SUBROUTINE OLD_DIHEDRAL

      !> Following Blondel&Karplus (Journal of Computational Chemistry, Vol. 17, No. 9, 1132-1 141 (1996)
      ! for the calculation of dihedral foces. Chiral centres 
      SUBROUTINE DIHEDRAL2(I, A, B, C, D, DIHREF, E, FA, FB, FC, FD)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, EUC_NORM, DOTP
         USE QCIKEYS, ONLY: NIMAGES
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(IN) :: A(3), B(3), C(3), D(3)
         INTEGER, INTENT(IN) :: DIHREF
         INTEGER, INTENT(IN) :: I  !image n 
         REAL(KIND = REAL64), INTENT(OUT) :: E
         REAL(KIND = REAL64), INTENT(OUT) :: FA(3), FB(3), FC(3), FD(3)   

         REAL(KIND = REAL64) :: RAB(3), RBC(3), RCD(3)
         REAL(KIND = REAL64) :: NAB, NBC, NCD
         REAL(KIND = REAL64) :: N1(3), N2(3), NM1, NM2
         REAL(KIND = REAL64) :: PHI_REG, SINPHI, COSPHI, DF
         REAL(KIND = REAL64) :: DOT12, L1, L2, Z1, Z2, Z12, REGTERM1, CT0, CT1, SREG
         REAL(KIND = REAL64) :: FAB(3), FBC(3), FCD(3)
         REAL(KIND = REAL64) :: DIFF, THISREF, SINSUM, COSSUM
         INTEGER :: M
         !cutoffs used for regularisation
         REAL(KIND = REAL64), PARAMETER :: EPS9 = 1.0D-9
         REAL(KIND = REAL64), PARAMETER :: EPS6 = 1.0D-6
         REAL(KIND = REAL64), PARAMETER :: EPS3 = 1.0D-3
         REAL(KIND = REAL64), PARAMETER :: P9999 = 0.9999

         ! compute the vectors between atoms in order
         RAB(1:3) = B(1:3) - A(1:3)
         RBC(1:3) = C(1:3) - B(1:3)
         RCD(1:3) = D(1:3) - C(1:3)

         !Calculate their norms
         NAB = NORM2(RAB)
         NBC = NORM2(RBC)
         NCD = NORM2(RCD)

         ! get normal vectors
         N1(1:3) = CROSS_PROD(RAB,RBC)
         N2(1:3) = CROSS_PROD(RBC,RCD)  
         ! get norms
         NM1 = NORM2(N1)
         NM2 = NORM2(N2)
         DOT12 = DOT_PRODUCT(N1,N2) !DOTP(3,N1,N2)

         ! introduce regularisation step 1: step function to cutoff small values
         L1 = NM1 + EPS9
         L2 = NM2 + EPS9
         Z1 = 1.0D0/L1
         Z2 = 1.0D0/L2
         IF (EPS3.GT.L1) Z1 = 0.0D0
         IF (EPS3.GT.L2) Z2 = 0.0D0
         Z12 = Z1*Z2
         REGTERM1 = 0.0D0
         IF (Z12.NE.0.0D0) REGTERM1 = 1.0D0
         
         IF ((DOT12*Z12).GT.1.0D0)  WRITE(*,*) "WARNING DOTP > 1 something is very wrong with dihedrals"
         IF ((DOT12*Z12).LT.-1.0D0) WRITE(*,*) "WARNING DOTP < 1 something is very wrong with dihedrals"          
         
         ! second step of regularisation: restrict values to correct range (-1 to 1)
         CT0 = MIN(P9999,DOT12*Z12)
         CT1 = MAX(-P9999,CT0)
         CT1 = DOT12*Z12

         !This gives us correct dihedral angle, which corresponds to one given by VMD  
         PHI_REG = DSIGN(DACOS(CT1), DOTP(3,RBC,CROSS_PROD(N1,N2)))

         IF (DABS(PHI_REG-PI).LE.EPS3) PHI_REG = SIGN(PI,PHI_REG)
         IF (DABS(PHI_REG).LE.EPS3) PHI_REG = 0.0D0
         
         SINPHI = DSIN(PHI_REG)
         COSPHI = DCOS(PHI_REG)
         IF (DABS(SINPHI).LE.EPS6) SINPHI = 0.0D0
         IF (DABS(COSPHI).LE.EPS6) COSPHI = 0.0D0

         M = DIHMUL(DIHREF)
         
         SINSUM = ((NIMAGES+2-I)*DSIN(REFDIHS(DIHREF)) + (I-1)*DSIN(REFDIHF(DIHREF))) / (NIMAGES+1)
         COSSUM = ((NIMAGES+2-I)*DCOS(REFDIHS(DIHREF)) + (I-1)*COS(REFDIHF(DIHREF))) / (NIMAGES+1)
         THISREF = ATAN2(SINSUM, COSSUM)
         !DIFF = PHI_REG - REFDIH(DIHREF)
         !DIFF = ATAN2(DSIN(PHI_REG - THISREF), DCOS(PHI_REG - THISREF))
         
         !define M to choose which form of potential to use
         !M=1
         IF (M.EQ.1) THEN
            !We have chiral center and use quadratic potential
                    
            DIFF = ATAN2(DSIN(PHI_REG - THISREF), DCOS(PHI_REG - THISREF))
                      
            E = KDIH*DIFF**2  
            
            !First part of gradient calculation
            ! dE/d(cos(phi))          
            IF(DABS(DIFF).LT.EPS9) THEN
               DF = 0.0D0
            ELSE               
               DF = 2.0D0*KDIH*DIFF
            END IF
            
                
         ELSE IF (M.EQ.2) THEN
            !We have planar dihedral, use cosine with multiplicity 1   
            
            THISREF = THISREF + PI
            DIFF = ATAN2(DSIN(PHI_REG - THISREF), DCOS(PHI_REG - THISREF))
                        
            !Normal cosine dihedral      
            E = KDIH*(1.0D0+COSPHI*DCOS(THISREF)+SINPHI*DSIN(THISREF))*REGTERM1
           
            !First part of gradient calculation
            ! dE/d(phi)
            DF = -KDIH*(DCOS(THISREF)*SINPHI - DSIN(THISREF)*COSPHI) * REGTERM1
                       
         ELSE
            ! We should never be here!
            WRITE(*,*) "WARNING! We don't know what is this dihedral number: ", DIHREF
            DF = 0.0D0
         END IF
      
         
         ! Individual first derviatives with respect to the cartesian coords
         ! dphi/dr

         FAB = -NBC*Z1**2 * N1
         FBC = (DOT_PRODUCT(RAB,RBC)*Z1**2 / NBC)*N1 - (DOT_PRODUCT(RCD,RBC)*Z2**2 / NBC) * N2
         FCD = NBC*Z2**2 * N2
         
         FA = FAB * DF
         FB = (-FAB + FBC) * DF
         FC = (-FBC - FCD) * DF
         FD = FCD * DF 
                 
         !WRITE(*,*) "dihedral> DIHREF", DIHREF, " S0", S0(DIHREF), " C0 ", C0(DIHREF), " SINHPHI ", SINPHI, " COSPHI ", COSPHI
         !WRITE(*,*) "dihedral> E ", E, "PHI_0 ", THISREF , "PHI_REG ", PHI_REG, "DIFF= ", DIFF, "DF = ", DF
         !WRITE(*,*) "M", M
         !WRITE(*,*) "F_A ", FA 
         !WRITE(*,*) "F_B ", FB 
         !WRITE(*,*) "F_C ", FC 
         !WRITE(*,*) "F_D ", FD 

      END SUBROUTINE DIHEDRAL2    

END MODULE DIHEDRAL_CONSTRAINTS
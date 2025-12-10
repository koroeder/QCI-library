MODULE DIHEDRAL_CONSTRAINTS
   USE QCIPREC
   IMPLICIT NONE
   ! number of dihedrals to be constraint
   INTEGER :: NDIH
   ! atoms in dihedrals (by id)
   INTEGER, ALLOCATABLE :: DIHEDRALS(:,:)
   ! reference dihedrals
   REAL(KIND = REAL64), ALLOCATABLE :: REFDIH(:)
   ! activated all dihedrals?
   LOGICAL :: ALLDIHACTIVE = .FALSE.
   ! activated dihedrals?
   LOGICAL, ALLOCATABLE :: DIHACTIVE(:)
   ! spring restraint constant
   REAL(KIND = REAL64) :: KDIH

   CONTAINS
      SUBROUTINE ALLOC_DIHVARS()
         CALL DEALLOC_DIHVARS
         ALLOCATE(DIHEDRALS(NDIH,4))
         ALLOCATE(REFDIH(NDIH))
         ALLOCATE(DIHACTIVE(NDIH))
      END SUBROUTINE ALLOC_DIHVARS

      SUBROUTINE DEALLOC_DIHVARS()
         IF (ALLOCATED(DIHEDRALS)) DEALLOCATE(DIHEDRALS)
         IF (ALLOCATED(REFDIH)) DEALLOCATE(REFDIH)
         IF (ALLOCATED(DIHACTIVE)) DEALLOCATE(DIHACTIVE)
      END SUBROUTINE DEALLOC_DIHVARS

      SUBROUTINE SETUP_DIH_CONSTR()
         USE QCIKEYS, ONLY: NATOMS
         USE INTERPOLATION_KEYS, ONLY: XSTART
         USE AMBER_CONSTRAINTS, ONLY: RESNAMES, NRES, GET_ATOMID
         USE CHIRALITY, ONLY: NCHIRAL, CHIR_INFO

         INTEGER :: REFATOMS(NATOMS,4)
         INTEGER :: NCONS, AT1, AT2, AT3, AT4, I, J

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
               (RESNAMES(I).EQ."DA").OR.(RESNAMES(I).EQ."DA3").OR.(RESNAMES(I).EQ."DA5")) THEN
               WRITE(*,*) "Dihedrals for A"
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
               WRITE(*,*) REFATOMS(NCONS,1:4)  
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)
               !fourth dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("N7",I,AT2)
               CALL GET_ATOMID("C8",I,AT3)
               CALL GET_ATOMID("N9",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               WRITE(*,*) REFATOMS(NCONS,1:4)
            ELSE IF ((RESNAMES(I).EQ."C").OR.(RESNAMES(I).EQ."C3").OR.(RESNAMES(I).EQ."C5").OR. &
                  (RESNAMES(I).EQ."DC").OR.(RESNAMES(I).EQ."DC3").OR.(RESNAMES(I).EQ."DC5")) THEN
               WRITE(*,*) "Dihedrals for C"
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N1",I,AT1)
               CALL GET_ATOMID("C2",I,AT2)
               CALL GET_ATOMID("N3",I,AT3)
               CALL GET_ATOMID("C4",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               WRITE(*,*) REFATOMS(NCONS,1:4)
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N3",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4   
               WRITE(*,*) REFATOMS(NCONS,1:4)  
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)
            ELSE IF ((RESNAMES(I).EQ."G").OR.(RESNAMES(I).EQ."G3").OR.(RESNAMES(I).EQ."G5").OR. &
                     (RESNAMES(I).EQ."DG").OR.(RESNAMES(I).EQ."DG3").OR.(RESNAMES(I).EQ."DG5")) THEN
               WRITE(*,*) "Dihedrals for G"
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N9",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               WRITE(*,*) REFATOMS(NCONS,1:4)
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N7",I,AT1)
               CALL GET_ATOMID("C5",I,AT2)
               CALL GET_ATOMID("C4",I,AT3)
               CALL GET_ATOMID("N3",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               WRITE(*,*) REFATOMS(NCONS,1:4)     
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)
               !fourth dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("N7",I,AT2)
               CALL GET_ATOMID("C8",I,AT3)
               CALL GET_ATOMID("N9",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)
            ELSE IF ((RESNAMES(I).EQ."U").OR.(RESNAMES(I).EQ."U3").OR.(RESNAMES(I).EQ."U5").OR. &
                     (RESNAMES(I).EQ."DT").OR.(RESNAMES(I).EQ."DT3").OR.(RESNAMES(I).EQ."DT5")) THEN           
               WRITE(*,*) "Dihedrals for U/T"
               !first dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N1",I,AT1)
               CALL GET_ATOMID("C2",I,AT2)
               CALL GET_ATOMID("N3",I,AT3)
               CALL GET_ATOMID("C4",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4
               WRITE(*,*) REFATOMS(NCONS,1:4)
               !second dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("N3",I,AT1)
               CALL GET_ATOMID("C4",I,AT2)
               CALL GET_ATOMID("C5",I,AT3)
               CALL GET_ATOMID("C6",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)    
               !third dihedral 
               NCONS = NCONS + 1
               CALL GET_ATOMID("C5",I,AT1)
               CALL GET_ATOMID("C6",I,AT2)
               CALL GET_ATOMID("N1",I,AT3)
               CALL GET_ATOMID("C2",I,AT4)
               REFATOMS(NCONS,1) = AT1; REFATOMS(NCONS,2) = AT2; REFATOMS(NCONS,3) = AT3; REFATOMS(NCONS,4) = AT4 
               WRITE(*,*) REFATOMS(NCONS,1:4)
            END IF
         END DO
         WRITE(*,*) REFATOMS
         !assign arrays for dihedrals
         NDIH = NCHIRAL + NCONS
         WRITE(*,*) "NDIH,NCHIRAL,NCONS: ", NDIH, NCHIRAL, NCONS
         CALL ALLOC_DIHVARS()
         DIHACTIVE(1:NDIH) = .FALSE.
         DO J=1,NCHIRAL
            WRITE(*,*) "J: ", J
            DIHEDRALS(J,1) = CHIR_INFO(J,1)
            DIHEDRALS(J,2) = CHIR_INFO(J,2)
            DIHEDRALS(J,3) = CHIR_INFO(J,3)
            DIHEDRALS(J,4) = CHIR_INFO(J,4)
         END DO
         WRITE(*,*) "NCONS loop"
         DO J=1,NCONS
            WRITE(*,*) "J: ", NCHIRAL+J
            DIHEDRALS(NCHIRAL+J,1) = REFATOMS(J,1)
            DIHEDRALS(NCHIRAL+J,2) = REFATOMS(J,2)
            DIHEDRALS(NCHIRAL+J,3) = REFATOMS(J,3)
            DIHEDRALS(NCHIRAL+J,4) = REFATOMS(J,4)            
         END DO
         ! get reference angles
         DO J=1,NDIH
            AT1 = DIHEDRALS(J,1); AT2 = DIHEDRALS(J,2); AT3 = DIHEDRALS(J,3); AT4 = DIHEDRALS(J,4)
            CALL COMPUTE_DIH(NATOMS, XSTART, AT1, AT2, AT3, AT4, REFDIH(I))
         END DO
      END SUBROUTINE SETUP_DIH_CONSTR

      SUBROUTINE COMPUTE_DIH(NATOMS, X, A, B, C, D, PHI)
         USE HELPER_FNCTS, ONLY: CROSS_PROD, EUC_NORM, DOTP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)
         INTEGER, INTENT(IN) :: A, B, C, D
         REAL(KIND = REAL64), INTENT(OUT) :: PHI

         REAL(KIND = REAL64) :: RAB(3), RBC(3), RCD(3)
         REAL(KIND = REAL64) :: N1(3), N2(3), NORM1, NORM2, NORMBC, NPROD
         REAL(KIND = REAL64) :: COSPHI, SINPHI

         WRITE(*,*) "ABCD ", A, B, C, D
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
         NPROD = 1.0D0/(NORM1*NORM2*NORMBC)
         N1(1:3) = N1(1:3)/NORM1
         N2(1:3) = N2(1:3)/NORM2

         !comput cos and sin components and get dihedral
         COSPHI = DOTP(3,N1,N2)
         SINPHI = DOTP(3,CROSS_PROD(N1,N2),RBC)*NPROD
         PHI = ATAN2(SINPHI,COSPHI)
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

      SUBROUTINE DIHEDRAL(A, B, C, D, PHIREF, E, FA, FB, FC, FD)
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

         WRITE(*,*) ">>> dih_values phi - ", PHI, " ,ref: ", PHIREF

         E = KDIH*(PHI-PHIREF)**2

         FE = 2*KDIH*(PHI-PHIREF)
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
      END SUBROUTINE DIHEDRAL
END MODULE DIHEDRAL_CONSTRAINTS
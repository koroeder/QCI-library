!> Helper module to detect bond crossings
MODULE BOND_CROSSING_DETECTION
       
   USE QCIPREC
   IMPLICIT NONE
   REAL(KIND=REAL64), PARAMETER :: EPS = 1.0D-6
   
   CONTAINS

   !>  ------------------
   !> CHECK IF A POINT LIES INSIDE A TRIANGLE (3D).
   !> ASSUMES P IS COPLANAR WITH THE TRIANGLE.
   !> USES CROSS-PRODUCT SIGN TEST RELATIVE TO NORMAL N.
   !>------------------------------------------------------------------------------
   PURE FUNCTION POINT_IN_TRIANGLE(P, V0, V1, V2, N) RESULT(INSIDE)
    
      USE HELPER_FNCTS, ONLY: DOTP, CROSS_PROD
      
      REAL(KIND = REAL64), INTENT(IN):: P(3), V0(3), V1(3), V2(3), N(3)
      LOGICAL :: INSIDE
      REAL(KIND = REAL64) :: C(3), E(3), D
      LOGICAL :: POS, NEG
     

      POS = .FALSE.
      NEG = .FALSE.

      E = V1 - V0
      C = CROSS_PROD(E, P - V0)
      D = DOTP(3, C, N)
      IF (D.GT.EPS)  POS = .TRUE.
      IF (D.LT.-EPS) NEG = .TRUE.
      IF (POS.AND.NEG) THEN
         INSIDE = .FALSE.
         RETURN
      END IF

      E = V2 - V1
      C = CROSS_PROD(E, P - V1)
      D = DOTP(3,C, N)
      IF (D.GT.EPS)  POS = .TRUE.
      IF (D.LT.-EPS) NEG = .TRUE.
      IF (POS.AND.NEG) THEN
         INSIDE = .FALSE.
         RETURN
      END IF

      E = V0 - V2
      C = CROSS_PROD(E, P - V2)
      D = DOTP(3,C, N)
      IF (D.GT.EPS)  POS = .TRUE.
      IF (D.LT.-EPS) NEG = .TRUE.
      IF (POS.AND.NEG) THEN
         INSIDE = .FALSE.
         RETURN
      END IF

      INSIDE = .TRUE.

   END FUNCTION POINT_IN_TRIANGLE

   !>------------------------------------------------------------------------------
   !> SEGMENT-TRIANGLE INTERSECTION IN 3D.
   !> RETURNS .TRUE. IF THE OPEN SEGMENT (P1,P2) PIERCES THE TRIANGLE INTERIOR.
   !> COPLANAR SEGMENTS ARE IGNORED HERE; HANDLED BY THE CALLER.
   !>------------------------------------------------------------------------------
   PURE FUNCTION SEGMENT_TRIANGLE_INTERSECT(P1, P2, V0, V1, V2) RESULT(INTERSECT)
      USE HELPER_FNCTS, ONLY: DOTP, CROSS_PROD, EUC_NORM
      
      REAL(KIND = REAL64), INTENT(IN) :: P1(3), P2(3), V0(3), V1(3), V2(3)
      LOGICAL :: INTERSECT
      REAL(KIND = REAL64) :: N(3), E1(3), E2(3), DIR(3), W0(3)
      REAL(KIND = REAL64) :: NV, T, NORM_N, NORM_DIR
      REAL(KIND = REAL64) :: ISEC(3)

      INTERSECT = .FALSE.

      E1 = V1 - V0
      E2 = V2 - V0
      N  = CROSS_PROD(E1, E2)
      NORM_N = EUC_NORM(N)
      IF (NORM_N < EPS) RETURN          ! Degenerate triangle

      DIR = P2 - P1
      NORM_DIR = EUC_NORM(DIR)
      IF (NORM_DIR < EPS) RETURN        ! Zero-length segment
      
      NV = DOTP(3, N, DIR)
      
      ! Relative tolerance: reject only if genuinely parallel
      IF (ABS(NV) < EPS * NORM_N * NORM_DIR) RETURN

      W0 = P1 - V0
      T  = -DOTP(3, N, W0) / NV

      ! Closed segment [0,1] (allow tiny epsilon for round-off)
      IF (T < -EPS .OR. T > 1.0D0 + EPS) RETURN

      ISEC = P1 + T * DIR

      ! Normalise N for the point-in-triangle test
      N = N / NORM_N

      INTERSECT = POINT_IN_TRIANGLE(ISEC, V0, V1, V2, N)
   END FUNCTION SEGMENT_TRIANGLE_INTERSECT

  !>------------------------------------------------------------------------------
  !> PROJECT A TRIANGLE ONTO THE AXIS-ALIGNED PLANE WHERE THE NORMAL
  !> HAS THE SMALLEST COMPONENT.  RETURNS 2D COORDINATES.
  !>------------------------------------------------------------------------------
  PURE SUBROUTINE PROJECT_TRIANGLE(V0, V1, V2, N, X0, X1, X2)
    REAL(KIND = REAL64), INTENT(IN)  :: V0(3), V1(3), V2(3), N(3)
    REAL(KIND = REAL64), INTENT(OUT) :: X0(2), X1(2), X2(2)
    INTEGER :: I1, I2

    IF (ABS(N(1)).GE.ABS(N(2)) .AND. ABS(N(1)).GE.ABS(N(3))) THEN
       I1 = 2;  I2 = 3   ! DROP X
    ELSE IF (ABS(N(2)) .GE. ABS(N(3))) THEN
       I1 = 1;  I2 = 3   ! DROP Y
    ELSE
       I1 = 1;  I2 = 2   ! DROP Z
    END IF

    X0 = [V0(I1), V0(I2)]
    X1 = [V1(I1), V1(I2)]
    X2 = [V2(I1), V2(I2)]
  END SUBROUTINE PROJECT_TRIANGLE

  !>------------------------------------------------------------------------------
  !> 2D TRIANGLE-TRIANGLE INTERSECTION USING SEPARATING-AXIS THEOREM.
  !>------------------------------------------------------------------------------
   PURE FUNCTION TRIANGLE_INTERSECT_2D(V0, V1, V2, U0, U1, U2) RESULT(INTERSECT)
      
      REAL(KIND = REAL64), INTENT(IN) :: V0(2), V1(2), V2(2), U0(2), U1(2), U2(2)
      LOGICAL :: INTERSECT
      REAL(KIND = REAL64) :: AXES(2, 6), E(2), VAL
      REAL(KIND = REAL64) :: MINV, MAXV, MINU, MAXU
      INTEGER  :: I, K
      REAL(KIND = REAL64) :: PTSV(2,3), PTSU(2,3)

      PTSV(:,1) = V0;  PTSV(:,2) = V1;  PTSV(:,3) = V2
      PTSU(:,1) = U0;  PTSU(:,2) = U1;  PTSU(:,3) = U2

      ! 6 SEPARATING AXES: NORMAL TO EACH OF THE 3 EDGES OF EACH TRIANGLE
      AXES(:,1) = [  V1(2)-V0(2), -(V1(1)-V0(1)) ]   ! EDGE V0->V1
      AXES(:,2) = [  V2(2)-V1(2), -(V2(1)-V1(1)) ]   ! EDGE V1->V2
      AXES(:,3) = [  V0(2)-V2(2), -(V0(1)-V2(1)) ]   ! EDGE V2->V0
      AXES(:,4) = [  U1(2)-U0(2), -(U1(1)-U0(1)) ]   ! EDGE U0->U1
      AXES(:,5) = [  U2(2)-U1(2), -(U2(1)-U1(1)) ]   ! EDGE U1->U2
      AXES(:,6) = [  U0(2)-U2(2), -(U0(1)-U2(1)) ]   ! EDGE U2->U0

      DO I = 1, 6
         E = AXES(:, I)
         IF (ABS(E(1)) .LT. EPS .AND. ABS(E(2)) .LT. EPS) CYCLE

         MINV =  HUGE(1.0D0);  MAXV = -HUGE(1.0D0)
         MINU =  HUGE(1.0D0);  MAXU = -HUGE(1.0D0)

         DO K = 1, 3
            VAL = PTSV(1,K)*E(1) + PTSV(2,K)*E(2)
            MINV = MIN(MINV, VAL);  MAXV = MAX(MAXV, VAL)
            VAL = PTSU(1,K)*E(1) + PTSU(2,K)*E(2)
            MINU = MIN(MINU, VAL);  MAXU = MAX(MAXU, VAL)
         END DO

         IF (MAXV .LT. MINU - EPS .OR. MAXU .LT. MINV - EPS) THEN
            INTERSECT = .FALSE.
            RETURN
         END IF
      END DO

      INTERSECT = .TRUE.
   END FUNCTION TRIANGLE_INTERSECT_2D

   !>------------------------------------------------------------------------------
   !> 3D TRIANGLE-TRIANGLE INTERSECTION.
   !> FIRST TESTS ALL EDGE-TRIANGLE PAIRS, THEN FALLS BACK TO A 2D SAT TEST
   !> IF THE TRIANGLES ARE COPLANAR.
   !>------------------------------------------------------------------------------
   PURE FUNCTION TRIANGLE_INTERSECT_3D(V0, V1, V2, U0, U1, U2) RESULT(INTERSECT)
   USE HELPER_FNCTS, ONLY: DOTP, CROSS_PROD, EUC_NORM
      
   REAL(KIND = REAL64), INTENT(IN) :: V0(3), V1(3), V2(3), U0(3), U1(3), U2(3)
   LOGICAL :: INTERSECT
   REAL(KIND = REAL64) :: N(3), E1(3), E2(3)
   REAL(KIND = REAL64) :: D0, D1, D2, SCALE
   REAL(KIND = REAL64) :: V0_2D(2), V1_2D(2), V2_2D(2), U0_2D(2), U1_2D(2), U2_2D(2)

   INTERSECT = .FALSE.

   ! 1. EDGE-TRIANGLE TESTS
   IF (SEGMENT_TRIANGLE_INTERSECT(V0, V1, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
   IF (SEGMENT_TRIANGLE_INTERSECT(V1, V2, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
   IF (SEGMENT_TRIANGLE_INTERSECT(V2, V0, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF

   IF (SEGMENT_TRIANGLE_INTERSECT(U0, U1, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
   IF (SEGMENT_TRIANGLE_INTERSECT(U1, U2, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
   IF (SEGMENT_TRIANGLE_INTERSECT(U2, U0, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF

   ! 2. NO EDGE PIERCES.  CHECK IF COPLANAR (RELATIVE TOLERANCE).
   E1 = V1 - V0
   E2 = V2 - V0
   N = CROSS_PROD(E1, E2)
   N = N / EUC_NORM(N)

   D0 = ABS(DOTP(3, U0 - V0, N))
   D1 = ABS(DOTP(3, U1 - V0, N))
   D2 = ABS(DOTP(3, U2 - V0, N))

   ! Scale relative to the size of the triangles
   SCALE = MAX(EUC_NORM(V1-V0), EUC_NORM(V2-V0), EUC_NORM(V2-V1), &
               EUC_NORM(U1-U0), EUC_NORM(U2-U0), EUC_NORM(U2-U1))
   IF (SCALE < EPS) SCALE = 1.0D0

   IF (D0 > EPS*SCALE .OR. D1 > EPS*SCALE .OR. D2 > EPS*SCALE) THEN
      INTERSECT = .FALSE.
      RETURN
   END IF

   ! 3. COPLANAR: PROJECT TO 2D AND TEST OVERLAP
   CALL PROJECT_TRIANGLE(V0, V1, V2, N, V0_2D, V1_2D, V2_2D)
   CALL PROJECT_TRIANGLE(U0, U1, U2, N, U0_2D, U1_2D, U2_2D)
   INTERSECT = TRIANGLE_INTERSECT_2D(V0_2D, V1_2D, V2_2D, U0_2D, U1_2D, U2_2D)
END FUNCTION TRIANGLE_INTERSECT_3D
  !>------------------------------------------------------------------------------
  !> CHECK IF TWO QUADRILATERALS (EACH SPLIT INTO TWO TRIANGLES) INTERSECT.
  !> QUAD1: [P1, P2, P3, P4]  -> TRIANGLES (P1,P2,P3) AND (P1,P3,P4)
  !> QUAD2: [Q1, Q2, Q3, Q4]  -> TRIANGLES (Q1,Q2,Q3) AND (Q1,Q3,Q4)
  !>------------------------------------------------------------------------------
  PURE FUNCTION QUADRILATERAL_INTERSECT(P1, P2, P3, P4, Q1, Q2, Q3, Q4) RESULT(INTERSECT)
    REAL(KIND = REAL64), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3)
    REAL(KIND = REAL64), INTENT(IN) :: Q1(3), Q2(3), Q3(3), Q4(3)
    LOGICAL :: INTERSECT

    INTERSECT = .FALSE.

    ! 4 TRIANGLE-TRIANGLE TESTS
    IF (TRIANGLE_INTERSECT_3D(P1, P2, P3, Q1, Q2, Q3)) THEN; INTERSECT = .TRUE.; RETURN; END IF
    IF (TRIANGLE_INTERSECT_3D(P1, P2, P3, Q1, Q3, Q4)) THEN; INTERSECT = .TRUE.; RETURN; END IF
    IF (TRIANGLE_INTERSECT_3D(P1, P3, P4, Q1, Q2, Q3)) THEN; INTERSECT = .TRUE.; RETURN; END IF
    IF (TRIANGLE_INTERSECT_3D(P1, P3, P4, Q1, Q3, Q4)) THEN; INTERSECT = .TRUE.; RETURN; END IF
  END FUNCTION QUADRILATERAL_INTERSECT

   !>------------------------------------------------------------------------------
   !> DETECT BOND CROSSINGS BETWEEN CONSECUTIVE IMAGES.
   !>
   !> ONLY CHECKS BACKBONE BONDS (BBOND_LIST) AGAINST EACH OTHER.
   !> SPATIAL CULLING SKIPS PAIRS WHOSE CENTRES ARE > DIST_CUTOFF APART.
   !>
   !> INPUTS:
   !>   COORDS(3,NATOMS,NIMAGES) : COORDINATES FOR ALL IMAGES
   !>   NATOMS                   : NUMBER OF ATOMS
   !>   NIMAGES                  : NUMBER OF IMAGES
   !>   NBBONDS                  : NUMBER OF BACKBONE BONDS
   !>   BBOND_LIST(2,NBBONDS)    : BACKBONE BOND ATOM INDICES (1-BASED)
   !>   DIST_CUTOFF              : OPTIONAL SPATIAL CULLING RADIUS (ANGSTROM)
   !>
   !> OUTPUTS:
   !>   IS_CROSSED(NIMAGES-1)    : .TRUE. IF ANY BACKBONE PAIR CROSSES IMG->IMG+1
   !>   N_CROSSINGS              : TOTAL NUMBER OF CROSSING PAIRS (OPTIONAL)
   !>   CROSSING_PAIRS(NIMAGES-1): NUMBER OF CROSSING PAIRS PER TRANSITION (OPT)
   !>------------------------------------------------------------------------------
   SUBROUTINE DETECT_BOND_CROSSINGS(XYZ) 
 
      USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
      USE QCIKEYS, ONLY: NIMAGES, NATOMS, ISBBATOM
      USE QCI_CONSTRAINT_KEYS, ONLY: BOND_LIST, NBONDS
      
      REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2)) 
  
      REAL(KIND = REAL64) :: DIST_CUTOFF
     
      INTEGER :: N_CROSSINGS
      INTEGER :: CROSSING_PAIRS(NIMAGES-1)  
      
      LOGICAL :: IS_CROSSED(NIMAGES-1)
      

      INTEGER :: IMG, B1, B2, I1, J1, I2, J2, PAIR_COUNT, TOTAL_EVENTS
      INTEGER :: IMOFFSET, IMOFFSET1
      REAL(KIND = REAL64) :: P1(3), P2(3), P3(3), P4(3) !< points for intersect check
      REAL(KIND = REAL64) :: Q1(3), Q2(3), Q3(3), Q4(3)
      
      LOGICAL :: SHARE_ATOM, INTERSECT
      REAL(KIND = REAL64) :: RCUT, D
      REAL(KIND = REAL64) :: C1(3), C2(3)
      REAL(KIND = REAL64) :: DEFAULT_DIST_CUTOFF

      DEFAULT_DIST_CUTOFF = 100.0D0
      RCUT = DEFAULT_DIST_CUTOFF
     

      IS_CROSSED = .FALSE.
      TOTAL_EVENTS = 0
      CROSSING_PAIRS = 0

      DO IMG = 1, NIMAGES - 1
         PAIR_COUNT = 0

         DO B1 = 1, NBONDS-1
            
            
            
            I1 = BOND_LIST(1, B1)
            J1 = BOND_LIST(2, B1)

            !Loop over backbone bonds only
            IF (.NOT. (ISBBATOM(I1).AND.ISBBATOM(J1))) CYCLE

            IF (I1.LT.1 .OR. I1 .GT. NATOMS .OR. J1 .LT. 1 .OR. J1 .GT. NATOMS) CYCLE

            ! PRECOMPUTE CENTROID OF BOND B1 FOR THIS IMAGE PAIR

            IMOFFSET = 3*NATOMS*(IMG-1) 
            IMOFFSET1 = 3*NATOMS*IMG 

            C1(1) = XYZ(IMOFFSET+3*(I1-1)+1) + XYZ(IMOFFSET+3*(J1-1)+1) + XYZ(IMOFFSET1+3*(I1-1)+1) + XYZ(IMOFFSET1+3*(J1-1)+1)
            C1(2) = XYZ(IMOFFSET+3*(I1-1)+2) + XYZ(IMOFFSET+3*(J1-1)+2) + XYZ(IMOFFSET1+3*(I1-1)+2) + XYZ(IMOFFSET1+3*(J1-1)+2)
            C1(3) = XYZ(IMOFFSET+3*(I1-1)+3) + XYZ(IMOFFSET+3*(J1-1)+3) + XYZ(IMOFFSET1+3*(I1-1)+3) + XYZ(IMOFFSET1+3*(J1-1)+3)
            
            C1 = C1 * 0.250D0

            DO B2 = B1 + 1, NBONDS
                              
               I2 = BOND_LIST(1, B2)
               J2 = BOND_LIST(2, B2)

               !Loop over backbone atoms only
               IF (.NOT. (ISBBATOM(I2).AND.ISBBATOM(J2))) CYCLE

               IF (I2.LT.1 .OR. I2 .GT. NATOMS .OR. J2 .LT. 1 .OR. J2 .GT. NATOMS) CYCLE

               ! SKIP ADJACENT BACKBONE BONDS (SHARE AN ATOM)
               SHARE_ATOM = (I1.EQ.I2) .OR. (I1.EQ.J2) .OR. (J1.EQ.I2) .OR. (J1.EQ.J2)
               IF (SHARE_ATOM) CYCLE

               ! --- SPATIAL CULLING: SKIP IF CENTRES ARE TOO FAR APART ---
                              
               C2(1) = XYZ(IMOFFSET+3*(I2-1)+1) + XYZ(IMOFFSET+3*(J2-1)+1) + XYZ(IMOFFSET1+3*(I2-1)+1) + XYZ(IMOFFSET1+3*(J2-1)+1)
               C2(2) = XYZ(IMOFFSET+3*(I2-1)+2) + XYZ(IMOFFSET+3*(J2-1)+2) + XYZ(IMOFFSET1+3*(I2-1)+2) + XYZ(IMOFFSET1+3*(J2-1)+2)
               C2(3) = XYZ(IMOFFSET+3*(I2-1)+3) + XYZ(IMOFFSET+3*(J2-1)+3) + XYZ(IMOFFSET1+3*(I2-1)+3) + XYZ(IMOFFSET1+3*(J2-1)+3)
            
               C2 = C2 * 0.25D0
               
               CALL DISTANCE_SIMPLE(C1, C2, D)
               
               IF (D.GT.RCUT) CYCLE

               P1 = XYZ(IMOFFSET+3*(I1-1)+1:IMOFFSET+3*(I1-1)+3)
               P2 = XYZ(IMOFFSET+3*(J1-1)+1:IMOFFSET+3*(J1-1)+3)
               P4 = XYZ(IMOFFSET1+3*(I1-1)+1:IMOFFSET1+3*(I1-1)+3)
               P3 = XYZ(IMOFFSET1+3*(J1-1)+1:IMOFFSET1+3*(J1-1)+3)

               Q1 = XYZ(IMOFFSET+3*(I2-1)+1:IMOFFSET+3*(I2-1)+3)
               Q2 = XYZ(IMOFFSET+3*(J2-1)+1:IMOFFSET+3*(J2-1)+3)
               Q4 = XYZ(IMOFFSET1+3*(I2-1)+1:IMOFFSET1+3*(I2-1)+3)
               Q3 = XYZ(IMOFFSET1+3*(J2-1)+1:IMOFFSET1+3*(J2-1)+3)

               ! --- EXPENSIVE QUADRILATERAL-QUADRILATERAL INTERSECTION TEST ---
               INTERSECT = QUADRILATERAL_INTERSECT(P1, P2, P3, P4, Q1, Q2, Q3, Q4)
               !INTERSECT = PLANE_INTERSECT(P1, P2, P3, P4, Q1, Q2, Q3, Q4)

               IF (INTERSECT) THEN
                  IS_CROSSED(IMG) = .TRUE.
                  PAIR_COUNT = PAIR_COUNT + 1
               END IF

               

            END DO
         END DO

         !IF (PRESENT(CROSSING_PAIRS)) CROSSING_PAIRS(IMG) = PAIR_COUNT
         TOTAL_EVENTS = TOTAL_EVENTS + PAIR_COUNT  

      END DO

      !IF (PRESENT(N_CROSSINGS)) N_CROSSINGS = TOTAL_EVENTS

      !Write output
      WRITE(*,*) "Detect_bond_crossing> Check for backbone bond crossings "
      IF (TOTAL_EVENTS.EQ.0) THEN
         WRITE(*,*) "Detect_bond_crossing> Bond crossing not detected. "
      ELSEIF (TOTAL_EVENTS.GT.0) THEN
         WRITE(*,*) " WARNING! Detected ", TOTAL_EVENTS, " bond crossings. "
      END IF 

   END SUBROUTINE DETECT_BOND_CROSSINGS


  PURE FUNCTION PLANE_INTERSECT(P1, P2, P3, P4, Q1, Q2, Q3, Q4) RESULT(INTERSECT)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL64
    IMPLICIT NONE
    REAL(KIND = REAL64), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3)
    REAL(KIND = REAL64), INTENT(IN) :: Q1(3), Q2(3), Q3(3), Q4(3)
    LOGICAL :: INTERSECT

    REAL(KIND = REAL64), PARAMETER :: EPS = 1.0D-10
    REAL(KIND = REAL64) :: NP(3), NQ(3), D(3), A(3)
    REAL(KIND = REAL64) :: NPNORM, NQNORM, DNORM, DET
    REAL(KIND = REAL64) :: RHS1, RHS2, TMIN_P, TMAX_P, TMIN_Q, TMAX_Q
    INTEGER :: K, I, J
    LOGICAL :: VALID

    INTERSECT = .FALSE.

    !--------------------------------------------------------------------------
    ! 1. Plane normals from the first three vertices of each quad
    !--------------------------------------------------------------------------
    ! Newell's method: robust normal from all 4 vertices
      NP(1) = (P1(2)-P2(2))*(P1(3)+P2(3)) &
            + (P2(2)-P3(2))*(P2(3)+P3(3)) &
            + (P3(2)-P4(2))*(P3(3)+P4(3)) &
            + (P4(2)-P1(2))*(P4(3)+P1(3))

      NP(2) = (P1(3)-P2(3))*(P1(1)+P2(1)) &
            + (P2(3)-P3(3))*(P2(1)+P3(1)) &
            + (P3(3)-P4(3))*(P3(1)+P4(1)) &
            + (P4(3)-P1(3))*(P4(1)+P1(1))

      NP(3) = (P1(1)-P2(1))*(P1(2)+P2(2)) &
            + (P2(1)-P3(1))*(P2(2)+P3(2)) &
            + (P3(1)-P4(1))*(P3(2)+P4(2)) &
            + (P4(1)-P1(1))*(P4(2)+P1(2))

    NPNORM = NORM(NP)
    IF (NPNORM < EPS) RETURN                 ! Degenerate plane P

    ! Newell's method for quad Q (uses all four vertices Q1, Q2, Q3, Q4)
      NQ(1) = (Q1(2)-Q2(2))*(Q1(3)+Q2(3)) &
            + (Q2(2)-Q3(2))*(Q2(3)+Q3(3)) &
            + (Q3(2)-Q4(2))*(Q3(3)+Q4(3)) &
            + (Q4(2)-Q1(2))*(Q4(3)+Q1(3))

      NQ(2) = (Q1(3)-Q2(3))*(Q1(1)+Q2(1)) &
            + (Q2(3)-Q3(3))*(Q2(1)+Q3(1)) &
            + (Q3(3)-Q4(3))*(Q3(1)+Q4(1)) &
            + (Q4(3)-Q1(3))*(Q4(1)+Q1(1))

      NQ(3) = (Q1(1)-Q2(1))*(Q1(2)+Q2(2)) &
            + (Q2(1)-Q3(1))*(Q2(2)+Q3(2)) &
            + (Q3(1)-Q4(1))*(Q3(2)+Q4(2)) &
            + (Q4(1)-Q1(1))*(Q4(2)+Q1(2))
    NQNORM = NORM(NQ)
    IF (NQNORM < EPS) RETURN                 ! Degenerate plane Q

    NP = NP / NPNORM
    NQ = NQ / NQNORM

    !--------------------------------------------------------------------------
    ! 2. Direction of the planes' intersection line
    !--------------------------------------------------------------------------
    D = CROSS(NP, NQ)
    DNORM = NORM(D)

    IF (DNORM < EPS) THEN
        ! Planes are parallel.  Distinct or coplanar?
        IF (ABS(DOT(Q1 - P1, NP)) > EPS) RETURN   ! Distinct parallel → no hit

        ! Coplanar: project to 2D and run SAT
        K = MAXLOC([ABS(NP(1)), ABS(NP(2)), ABS(NP(3))], DIM=1)
        IF (K == 1) THEN; I = 2; J = 3
        ELSE IF (K == 2) THEN; I = 1; J = 3
        ELSE; I = 1; J = 2
        END IF

        INTERSECT = QUADS_OVERLAP_2D(P1,P2,P3,P4,Q1,Q2,Q3,Q4,I,J)
        RETURN
    END IF

    D = D / DNORM

    !--------------------------------------------------------------------------
    ! 3. One point A on the intersection line.
    !    Set the coordinate where |D| is largest to zero and solve the 2×2
    !    system given by the two plane equations.
    !--------------------------------------------------------------------------
    K = 1
    IF (ABS(D(2)) > ABS(D(1))) K = 2
    IF (ABS(D(3)) > ABS(D(K))) K = 3

    IF (K == 1) THEN; I = 2; J = 3
    ELSE IF (K == 2) THEN; I = 1; J = 3
    ELSE; I = 1; J = 2
    END IF

    RHS1 = DOT(NP, P1)
    RHS2 = DOT(NQ, Q1)
    DET = NP(I)*NQ(J) - NP(J)*NQ(I)          ! = D(K), guaranteed non-zero

    A = 0.0D0
    A(I) = (RHS1 * NQ(J) - NP(J) * RHS2) / DET
    A(J) = (NP(I) * RHS2 - RHS1 * NQ(I)) / DET

    !--------------------------------------------------------------------------
    ! 4. Clip the infinite line A+t*D against each convex quad.
    !    This gives the 1-D interval (in t) where the line is inside the quad.
    !--------------------------------------------------------------------------
    CALL LINE_CLIP_QUAD(A, D, P1, P2, P3, P4, NP, TMIN_P, TMAX_P, VALID)
    IF (.NOT. VALID) RETURN

    CALL LINE_CLIP_QUAD(A, D, Q1, Q2, Q3, Q4, NQ, TMIN_Q, TMAX_Q, VALID)
    IF (.NOT. VALID) RETURN

    !--------------------------------------------------------------------------
    ! 5. The quads intersect iff the two intervals on the line overlap.
    !--------------------------------------------------------------------------
    INTERSECT = (TMAX_P >= TMIN_Q - EPS) .AND. (TMAX_Q >= TMIN_P - EPS)

    !==========================================================================
    ! Internal helpers
    !==========================================================================
CONTAINS

    PURE FUNCTION CROSS(U, V) RESULT(W)
        REAL(KIND = REAL64), INTENT(IN) :: U(3), V(3)
        REAL(KIND = REAL64) :: W(3)
        W(1) = U(2)*V(3) - U(3)*V(2)
        W(2) = U(3)*V(1) - U(1)*V(3)
        W(3) = U(1)*V(2) - U(2)*V(1)
    END FUNCTION CROSS

    PURE FUNCTION DOT(U, V) RESULT(S)
        REAL(KIND = REAL64), INTENT(IN) :: U(3), V(3)
        REAL(KIND = REAL64) :: S
        S = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)
    END FUNCTION DOT

    PURE FUNCTION NORM(U) RESULT(L)
        REAL(KIND = REAL64), INTENT(IN) :: U(3)
        REAL(KIND = REAL64) :: L
        L = SQRT(DOT(U, U))
    END FUNCTION NORM

    PURE FUNCTION TRIPLE(U, V, W) RESULT(T)
        REAL(KIND = REAL64), INTENT(IN) :: U(3), V(3), W(3)
        REAL(KIND = REAL64) :: T
        T = DOT(CROSS(U, V), W)
    END FUNCTION TRIPLE

    !>--------------------------------------------------------------------------
    !> Clip an infinite line against a convex quadrilateral.
    !> The line is parameterised as X(t) = A + t*D and is assumed to lie in the
    !> plane of the quad (D . N == 0).  Returns the interval [TMIN,TMAX] for
    !> which X(t) is inside the quad.  VALID is .FALSE. if the interval is empty.
    !>--------------------------------------------------------------------------
    PURE SUBROUTINE LINE_CLIP_QUAD(A, D, V1, V2, V3, V4, N, TMIN, TMAX, VALID)
        REAL(KIND = REAL64), INTENT(IN) :: A(3), D(3), V1(3), V2(3), V3(3), V4(3), N(3)
        REAL(KIND = REAL64), INTENT(OUT) :: TMIN, TMAX
        LOGICAL, INTENT(OUT) :: VALID

        REAL(KIND = REAL64) :: C(3), VI(3), VJ(3), E(3)
        REAL(KIND = REAL64) :: AT, BT, S
        REAL(KIND = REAL64) :: TLOW, THIGH
        INTEGER :: EDGE

        C = (V1 + V2 + V3 + V4) * 0.25D0
        TLOW = -HUGE(1.0D0)
        THIGH = HUGE(1.0D0)
        VALID = .TRUE.

        DO EDGE = 1, 4
            SELECT CASE(EDGE)
            CASE(1); VI = V1; VJ = V2
            CASE(2); VI = V2; VJ = V3
            CASE(3); VI = V3; VJ = V4
            CASE(4); VI = V4; VJ = V1
            END SELECT

            E = VJ - VI
            S = TRIPLE(E, C - VI, N)          ! sign of interior side

            AT = TRIPLE(E, D, N)              ! coefficient of t
            BT = TRIPLE(E, A - VI, N)         ! constant term

            IF (ABS(AT) < EPS) THEN
                ! Line parallel to this edge
                IF (S * BT < -EPS) THEN
                    VALID = .FALSE.
                    RETURN
                END IF
            ELSE
                ! Constraint: S*(AT*t + BT) >= 0
                IF (S * AT > 0.0D0) THEN
                    TLOW = MAX(TLOW, -BT / AT)
                ELSE
                    THIGH = MIN(THIGH, -BT / AT)
                END IF
            END IF
        END DO

        IF (TLOW > THIGH + EPS) THEN
            VALID = .FALSE.
        ELSE
            TMIN = TLOW
            TMAX = THIGH
        END IF
    END SUBROUTINE LINE_CLIP_QUAD

    !>--------------------------------------------------------------------------
    !> 2-D separating-axis test for two convex quads.
    !> The coordinates I and J are the axes kept after projecting out the
    !> dominant normal component.
    !>--------------------------------------------------------------------------
    PURE FUNCTION QUADS_OVERLAP_2D(A1,A2,A3,A4,B1,B2,B3,B4,IX,IY) RESULT(OVERLAP)
        REAL(KIND = REAL64), INTENT(IN) :: A1(3),A2(3),A3(3),A4(3)
        REAL(KIND = REAL64), INTENT(IN) :: B1(3),B2(3),B3(3),B4(3)
        INTEGER, INTENT(IN) :: IX, IY
        LOGICAL :: OVERLAP

        REAL(KIND = REAL64) :: VA(2,4), VB(2,4), AX(2)
        REAL(KIND = REAL64) :: MINA, MAXA, MINB, MAXB, VAL, ALEN
        INTEGER :: I, J

        VA(:,1) = [A1(IX), A1(IY)]; VA(:,2) = [A2(IX), A2(IY)]
        VA(:,3) = [A3(IX), A3(IY)]; VA(:,4) = [A4(IX), A4(IY)]
        VB(:,1) = [B1(IX), B1(IY)]; VB(:,2) = [B2(IX), B2(IY)]
        VB(:,3) = [B3(IX), B3(IY)]; VB(:,4) = [B4(IX), B4(IY)]

        DO I = 1, 8
            IF (I <= 4) THEN
                J = I
                IF (J == 4) THEN
                    AX = [VA(2,1)-VA(2,4), -(VA(1,1)-VA(1,4))]
                ELSE
                    AX = [VA(2,J+1)-VA(2,J), -(VA(1,J+1)-VA(1,J))]
                END IF
            ELSE
                J = I - 4
                IF (J == 4) THEN
                    AX = [VB(2,1)-VB(2,4), -(VB(1,1)-VB(1,4))]
                ELSE
                    AX = [VB(2,J+1)-VB(2,J), -(VB(1,J+1)-VB(1,J))]
                END IF
            END IF

            ALEN = SQRT(AX(1)**2 + AX(2)**2)
            IF (ALEN < EPS) CYCLE
            AX = AX / ALEN

            MINA = HUGE(1.0D0); MAXA = -HUGE(1.0D0)
            MINB = HUGE(1.0D0); MAXB = -HUGE(1.0D0)
            DO J = 1, 4
                VAL = VA(1,J)*AX(1) + VA(2,J)*AX(2)
                MINA = MIN(MINA, VAL); MAXA = MAX(MAXA, VAL)
                VAL = VB(1,J)*AX(1) + VB(2,J)*AX(2)
                MINB = MIN(MINB, VAL); MAXB = MAX(MAXB, VAL)
            END DO

            IF (MAXA < MINB - EPS .OR. MAXB < MINA - EPS) THEN
                OVERLAP = .FALSE.
                RETURN
            END IF
        END DO

        OVERLAP = .TRUE.
    END FUNCTION QUADS_OVERLAP_2D

END FUNCTION PLANE_INTERSECT
END MODULE BOND_CROSSING_DETECTION
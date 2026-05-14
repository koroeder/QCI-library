!> Helper module to detect bond crossings
MODULE BOND_CROSSING_DETECTION
       
   USE QCIPREC
   IMPLICIT NONE
   REAL(KIND=REAL64), PARAMETER :: EPS = 1.0E-12_REAL64
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
      REAL(KIND=REAL64), PARAMETER :: EPS = 1.0E-12_REAL64

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
      REAL(KIND = REAL64) :: NV, T
      REAL(KIND = REAL64) :: ISEC(3)

      INTERSECT = .FALSE.

      E1 = V1 - V0
      E2 = V2 - V0
      N = CROSS_PROD(E1, E2)

      DIR = P2 - P1
      NV  = DOTP(3, N, DIR)
      
      ! SEGMENT PARALLEL TO TRIANGLE PLANE
      IF (ABS(NV).LT.EPS) RETURN

      W0 = P1 - V0
      T  = -DOTP(3, N, W0) / NV

      ! INTERSECTION OUTSIDE SEGMENT BOUNDS
      IF ( (T.LE.-EPS) .OR. (T.GE.1.0D0-EPS)) RETURN

      ISEC = P1 + T * DIR

      !Normalise N
      N = N / EUC_NORM(N)

      ! CHECK IF INTERSECTION POINT LIES INSIDE TRIANGLE
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
      REAL(KIND = REAL64) :: D0, D1, D2
      REAL(KIND = REAL64) :: V0_2D(2), V1_2D(2), V2_2D(2), U0_2D(2), U1_2D(2), U2_2D(2)

      INTERSECT = .FALSE.

      ! 1. TEST THE 3 EDGES OF TRIANGLE V AGAINST TRIANGLE U
      IF (SEGMENT_TRIANGLE_INTERSECT(V0, V1, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
      IF (SEGMENT_TRIANGLE_INTERSECT(V1, V2, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
      IF (SEGMENT_TRIANGLE_INTERSECT(V2, V0, U0, U1, U2)) THEN; INTERSECT = .TRUE.; RETURN; END IF

      ! 2. TEST THE 3 EDGES OF TRIANGLE U AGAINST TRIANGLE V
      IF (SEGMENT_TRIANGLE_INTERSECT(U0, U1, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
      IF (SEGMENT_TRIANGLE_INTERSECT(U1, U2, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF
      IF (SEGMENT_TRIANGLE_INTERSECT(U2, U0, V0, V1, V2)) THEN; INTERSECT = .TRUE.; RETURN; END IF

      ! 3. NO EDGE PIERCES THE OTHER TRIANGLE.  FOR NON-COPLANAR TRIANGLES THIS
      !    MEANS THEY ARE DISJOINT.  IF COPLANAR, ONE MAY LIE INSIDE THE OTHER.
      E1 = V1 - V0
      E2 = V2 - V0
      N = CROSS_PROD(E1, E2)

      D0 = ABS(DOTP(3, U0 - V0, N))
      D1 = ABS(DOTP(3, U1 - V0, N))
      D2 = ABS(DOTP(3, U2 - V0, N))

      IF (D0 .GT. EPS .OR. D1 .GT. EPS .OR. D2 .GT. EPS) THEN
         ! NOT COPLANAR -> NO INTERSECTION
         INTERSECT = .FALSE.
         RETURN
      END IF

      ! 4. COPLANAR: PROJECT TO 2D AND TEST OVERLAP
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
   !SUBROUTINE DETECT_BOND_CROSSINGS(XYZ, DIST_CUTOFF, &
   !                                   IS_CROSSED, N_CROSSINGS, CROSSING_PAIRS)
      USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE
      USE QCIKEYS, ONLY: NIMAGES, NATOMS, ISBBATOM
      USE QCI_CONSTRAINT_KEYS, ONLY: BOND_LIST, NBONDS
      
      REAL(KIND = REAL64), INTENT(IN) :: XYZ(3*NATOMS*(NIMAGES+2)) 
      !REAL(KIND = REAL64), INTENT(IN), OPTIONAL :: DIST_CUTOFF
      !LOGICAL, INTENT(OUT) :: IS_CROSSED(NIMAGES-1)
      !INTEGER, INTENT(OUT), OPTIONAL :: N_CROSSINGS
      !INTEGER, INTENT(OUT), OPTIONAL :: CROSSING_PAIRS(NIMAGES-1)
      
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

      DEFAULT_DIST_CUTOFF = 10.0D0

      RCUT = DEFAULT_DIST_CUTOFF
      !IF (PRESENT(DIST_CUTOFF)) RCUT = DIST_CUTOFF


      IS_CROSSED = .FALSE.
      TOTAL_EVENTS = 0
      !IF (PRESENT(CROSSING_PAIRS)) CROSSING_PAIRS = 0
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

            IMOFFSET = 3*NATOMS*(IMG-1) + 1
            IMOFFSET1 = 3*NATOMS*IMG + 1

            C1(1) = XYZ(IMOFFSET+3*(I1-1)+1) + XYZ(IMOFFSET+3*(J1-1)+1) + XYZ(IMOFFSET1+3*(I1-1)+1) + XYZ(IMOFFSET1+3*(J1-1)+1)
            C1(2) = XYZ(IMOFFSET+3*(I1-1)+2) + XYZ(IMOFFSET+3*(J1-1)+2) + XYZ(IMOFFSET1+3*(I1-1)+2) + XYZ(IMOFFSET1+3*(J1-1)+2)
            C1(3) = XYZ(IMOFFSET+3*(I1-1)+3) + XYZ(IMOFFSET+3*(J1-1)+2) + XYZ(IMOFFSET1+3*(I1-1)+3) + XYZ(IMOFFSET1+3*(J1-1)+3)
            
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
               P3 = XYZ(IMOFFSET1+3*(I1-1)+1:IMOFFSET1+3*(I1-1)+3)
               P4 = XYZ(IMOFFSET1+3*(J1-1)+1:IMOFFSET1+3*(J1-1)+3)

               Q1 = XYZ(IMOFFSET1+3*(I2-1)+1:IMOFFSET1+3*(I2-1)+3)
               Q2 = XYZ(IMOFFSET+3*(J2-1)+1:IMOFFSET+3*(J2-1)+3)
               Q3 = XYZ(IMOFFSET1+3*(I2-1)+1:IMOFFSET1+3*(I2-1)+3)
               Q4 = XYZ(IMOFFSET1+3*(J2-1)+1:IMOFFSET1+3*(J2-1)+3)

               ! --- EXPENSIVE QUADRILATERAL-QUADRILATERAL INTERSECTION TEST ---
               INTERSECT = QUADRILATERAL_INTERSECT(P1, P2, P3, P4, Q1, Q2, Q3, Q4)

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

  
END MODULE BOND_CROSSING_DETECTION
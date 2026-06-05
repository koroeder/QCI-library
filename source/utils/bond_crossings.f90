!> Helper module to detect bond crossings
MODULE BOND_CROSSING_DETECTION
       
   USE QCIPREC
   IMPLICIT NONE
   REAL(KIND=REAL64), PARAMETER :: EPS = 1.0D-6
   
   CONTAINS

  
   SUBROUTINE DETECT_BOND_CROSSINGS(XYZ) 
 
      USE HELPER_FNCTS, ONLY: DISTANCE_SIMPLE, EUC_NORM, DOTP, CROSS_PROD
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
      REAL(KIND = REAL64) :: Q1(3), Q2(3), R1(3), R2(3), T1(3), T2(3)
      REAL(KIND = REAL64) :: PQ1(3), PQ2(3), RT1(3), RT2(3)
      REAL(KIND = REAL64) :: D1(3), D2(3), S(3)
      REAL(KIND = REAL64) :: COS_ANGLE
      
      LOGICAL :: SHARE_ATOM, INTERSECT
      REAL(KIND = REAL64) :: RCUT, D
      REAL(KIND = REAL64) :: C1(3), C2(3)
      REAL(KIND = REAL64) :: DEFAULT_DIST_CUTOFF

      DEFAULT_DIST_CUTOFF = 100.0D0
      RCUT = DEFAULT_DIST_CUTOFF
     

      IS_CROSSED = .FALSE.
      TOTAL_EVENTS = 0
      CROSSING_PAIRS = 0

      WRITE(*,*) "hello "
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


               !We are looking at the bond between atoms P and Q at image 1 and 2 (P1-Q1, P2-Q2)
               !similar for atoms R and T 

               P1 = XYZ(IMOFFSET+3*(I1-1)+1:IMOFFSET+3*(I1-1)+3)  
               P2 = XYZ(IMOFFSET1+3*(I1-1)+1:IMOFFSET1+3*(I1-1)+3)

               Q1 = XYZ(IMOFFSET+3*(J1-1)+1:IMOFFSET+3*(J1-1)+3)
               Q2 = XYZ(IMOFFSET1+3*(J1-1)+1:IMOFFSET1+3*(J1-1)+3)

               R1 = XYZ(IMOFFSET+3*(I2-1)+1:IMOFFSET+3*(I2-1)+3)
               R2 = XYZ(IMOFFSET1+3*(I2-1)+1:IMOFFSET1+3*(I2-1)+3)

               T1 = XYZ(IMOFFSET+3*(J2-1)+1:IMOFFSET+3*(J2-1)+3)
               T2 = XYZ(IMOFFSET1+3*(J2-1)+1:IMOFFSET1+3*(J2-1)+3)

              
               !bond 1 vector
               PQ1 = Q1 - P1
               PQ2 = Q2 - P2 

               !bond 2 vector
               RT1 = T1 - R1
               RT2 = T2 - R2


               !distance vector between bond1 and bond 2 
               !CALL shortest_displacement(P1, Q1, R1, T1, D1)
               !CALL shortest_displacement(P2, Q2, R2, T2, D2)

                ! Vector between bond midpoints in image IMG
               D1 = 0.5D0*(R1 + T1) - 0.5D0*(P1 + Q1)
               
               ! Vector between bond midpoints in image IMG+1
               D2 = 0.5D0*(R2 + T2) - 0.5D0*(P2 + Q2)

               !normalise these vectors
               D1 = D1 / EUC_NORM(D1)
               D2 = D2 / EUC_NORM(D2)

               !calculate cos between D1 and D2 

               COS_ANGLE = DOTP(3, D1, D2)

               INTERSECT = .FALSE.
               IF (COS_ANGLE.LE.0D0) THEN
                    WRITE(*,*) "POTENTIAL BOND CROSSING!!!"
                    INTERSECT = vec_proj_in_segment(D1, P1, Q1)
               END IF

                    

               IF (INTERSECT) THEN
                  IS_CROSSED(IMG) = .TRUE.
                  PAIR_COUNT = PAIR_COUNT + 1
               END IF

               

            END DO
         END DO

         !IF (PRESENT(CROSSING_PAIRS)) CROSSING_PAIRS(IMG) = PAIR_COUNT
         TOTAL_EVENTS = TOTAL_EVENTS + PAIR_COUNT  

      END DO

      
      !Write output
      WRITE(*,*) "Detect_bond_crossing> Check for backbone bond crossings "
      IF (TOTAL_EVENTS.EQ.0) THEN
         WRITE(*,*) "Detect_bond_crossing> Bond crossing not detected. "
      ELSEIF (TOTAL_EVENTS.GT.0) THEN
         WRITE(*,*) " WARNING! Detected ", TOTAL_EVENTS, " bond crossings. "
      END IF 

   END SUBROUTINE DETECT_BOND_CROSSINGS


   pure logical function vec_proj_in_segment(D1, P1, Q1) result(in_seg)
        implicit none
        real(8), intent(in) :: D1(3), P1(3), Q1(3)
        real(8) :: eps
        real(8) :: PQ1(3), sq_len, t, epsilon

        ! Default tolerance for double-precision noise
        epsilon = 1.0d-12
    

        ! Bond vector PQ1 = Q1 - P1
        PQ1 = Q1 - P1
        sq_len = dot_product(PQ1, PQ1)

        ! Degenerate bond (zero length)
        if (sq_len == 0.0d0) then
            in_seg = .false.
            return
        end if

        ! Scalar projection factor: proj_vec = t * PQ1
        t = dot_product(D1, PQ1) / sq_len

        ! Check bounds with tolerance: [-eps, 1+eps]
        in_seg = (t >= -epsilon .and. t <= 1.0d0 + epsilon)

    end function vec_proj_in_segment

    pure subroutine shortest_displacement(p1, q1, r1, t1, disp, s, t, dist_sq)
        
        real(kind=REAL64), intent(in)  :: p1(3), q1(3), r1(3), t1(3)
        real(kind=REAL64), intent(out) :: disp(3)
        real(kind=REAL64), intent(out), optional :: s, t, dist_sq

        real(kind=REAL64) :: u(3), v(3), w0(3)
        real(kind=REAL64) :: a, b, c, d, e, DET   ! <-- DET instead of D
        real(kind=REAL64) :: sN, sD, tN, tD, sc, tc
        real(kind=REAL64) :: c1(3), c2(3)
        real(kind=REAL64), parameter :: SMALL_NUM = 1.0D-12

        ! Segment directions
        u = q1 - p1
        v = t1 - r1
        w0 = p1 - r1

        a = dot_product(u, u)
        b = dot_product(u, v)
        c = dot_product(v, v)
        d = dot_product(u, w0)
        e = dot_product(v, w0)

        DET = a*c - b*b   ! <-- renamed here too
        sD = DET
        tD = DET

        ! --- Find the closest points on the two infinite lines ---
        if (DET < SMALL_NUM) then
            ! Parallel (or nearly so): pick s = 0, find closest t
            sN = 0.0D0
            sD = 1.0D0
            tN = e
            tD = c
        else
            sN = b*e - c*d
            tN = a*e - b*d
            ! Clamp s to [0,1] and recompute t accordingly
            if (sN < 0.0D0) then
                sN = 0.0D0
                tN = e
                tD = c
            else if (sN > sD) then
                sN = sD
                tN = e + b
                tD = c
            end if
        end if

        ! --- Clamp t to [0,1] and recompute s if necessary ---
        if (tN < 0.0D0) then
            tN = 0.0D0
            if (-d < 0.0D0) then
                sN = 0.0D0
            else if (-d > a) then
                sN = sD
            else
                sN = -d
                sD = a
            end if
        else if (tN > tD) then
            tN = tD
            if ((-d + b) < 0.0D0) then
                sN = 0.0D0
            else if ((-d + b) > a) then
                sN = sD
            else
                sN = (-d + b)
                sD = a
            end if
        end if

        ! --- Final parameters ---
        if (abs(sD) > SMALL_NUM) then
            sc = sN / sD
        else
            sc = 0.0D0
        end if

        if (abs(tD) > SMALL_NUM) then
            tc = tN / tD
        else
            tc = 0.0D0
        end if

        ! --- Closest points and displacement vector ---
        c1 = p1 + sc * u
        c2 = r1 + tc * v
        disp = c2 - c1   ! vector from bond1 to bond2

        if (present(dist_sq)) dist_sq = dot_product(disp, disp)
        if (present(s)) s = sc
        if (present(t)) t = tc

    end subroutine shortest_displacement


END MODULE BOND_CROSSING_DETECTION
!> Helper module to detect bond crossings
MODULE BOND_CROSSING_DETECTION
       
   USE QCIPREC
   IMPLICIT NONE
   REAL(KIND=REAL64), PARAMETER :: EPS = 1.0D-6
   
   CONTAINS

  LOGICAL FUNCTION BONDS_CROSS(A0,B0,C0,D0,A1,B1,C1,D1)
         ! does bond AB pass through bond CD as atoms move linearly from image 0 to image 1?
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: A0(3),B0(3),C0(3),D0(3),A1(3),B1(3),C1(3),D1(3)
         REAL(KIND=REAL64) :: U0(3),U1(3),V0(3),V1(3),W0(3),W1(3),CF(0:3)
         REAL(KIND=REAL64) :: TB(4),QA,QB,QC,DISC,Q,TLO,THI,TM,FLO,FM,TMP
         INTEGER :: NB,J,IT

         U0 = B0-A0; U1 = (B1-A1)-U0
         V0 = C0-A0; V1 = (C1-A1)-V0
         W0 = D0-A0; W1 = (D1-A1)-W0
         ! coplanarity f(t) = (B-A).((C-A)x(D-A)) is cubic in t
         CF(0) = TRIPLE(U0,V0,W0)
         CF(1) = TRIPLE(U1,V0,W0)+TRIPLE(U0,V1,W0)+TRIPLE(U0,V0,W1)
         CF(2) = TRIPLE(U1,V1,W0)+TRIPLE(U1,V0,W1)+TRIPLE(U0,V1,W1)
         CF(3) = TRIPLE(U1,V1,W1)

         ! split [0,1] at stationary points of f so each piece has at most one root
         NB = 1; TB(1) = 0.0D0
         QA = 3.0D0*CF(3); QB = 2.0D0*CF(2); QC = CF(1)
         IF (QA.NE.0.0D0) THEN
            DISC = QB*QB-4.0D0*QA*QC
            IF (DISC.GE.0.0D0) THEN
               Q = -0.5D0*(QB+SIGN(SQRT(DISC),QB))
               IF (Q.NE.0.0D0) THEN
                  CALL ADDT(Q/QA); CALL ADDT(QC/Q)
               END IF
            END IF
         ELSE IF (QB.NE.0.0D0) THEN
            CALL ADDT(-QC/QB)
         END IF
         IF ((NB.EQ.3).AND.(TB(3).LT.TB(2))) THEN
            TMP = TB(2); TB(2) = TB(3); TB(3) = TMP
         END IF
         NB = NB+1; TB(NB) = 1.0D0

         ! bisect each sign change, then test whether the coplanar segments intersect
         BONDS_CROSS = .FALSE.
         DO J = 1, NB-1
            TLO = TB(J); THI = TB(J+1)
            FLO = FCUB(TLO)
            IF (FLO*FCUB(THI).GT.0.0D0) CYCLE
            DO IT = 1, 60
               TM = 0.5D0*(TLO+THI); FM = FCUB(TM)
               IF (FLO*FM.LE.0.0D0) THEN
                  THI = TM
               ELSE
                  TLO = TM; FLO = FM
               END IF
            END DO
            IF (STRADDLE(0.5D0*(TLO+THI))) THEN
               BONDS_CROSS = .TRUE.
               RETURN
            END IF
         END DO

      CONTAINS
         SUBROUTINE ADDT(T)
            REAL(KIND=REAL64), INTENT(IN) :: T
            IF ((T.GT.0.0D0).AND.(T.LT.1.0D0)) THEN
               NB = NB+1; TB(NB) = T
            END IF
         END SUBROUTINE ADDT
         REAL(KIND=REAL64) FUNCTION FCUB(T)
            REAL(KIND=REAL64), INTENT(IN) :: T
            FCUB = CF(0)+T*(CF(1)+T*(CF(2)+T*CF(3)))
         END FUNCTION FCUB
         REAL(KIND=REAL64) FUNCTION TRIPLE(X,Y,Z)
            REAL(KIND=REAL64), INTENT(IN) :: X(3),Y(3),Z(3)
            TRIPLE = DOT_PRODUCT(X,CROSS3(Y,Z))
         END FUNCTION TRIPLE
         FUNCTION CROSS3(X,Y)
            REAL(KIND=REAL64), INTENT(IN) :: X(3),Y(3)
            REAL(KIND=REAL64) :: CROSS3(3)
            CROSS3 = (/X(2)*Y(3)-X(3)*Y(2), X(3)*Y(1)-X(1)*Y(3), X(1)*Y(2)-X(2)*Y(1)/)
         END FUNCTION CROSS3
         ! coplanar segments intersect iff each straddles the other's line
         LOGICAL FUNCTION STRADDLE(T)
            REAL(KIND=REAL64), INTENT(IN) :: T
            REAL(KIND=REAL64) :: A(3),B(3),C(3),D(3)
            A = A0+T*(A1-A0); B = B0+T*(B1-B0)
            C = C0+T*(C1-C0); D = D0+T*(D1-D0)
            STRADDLE = (DOT_PRODUCT(CROSS3(B-A,C-A),CROSS3(B-A,D-A)).LE.0.0D0).AND. &
                       (DOT_PRODUCT(CROSS3(D-C,A-C),CROSS3(D-C,B-C)).LE.0.0D0)
         END FUNCTION STRADDLE
      END FUNCTION BONDS_CROSS
   


END MODULE BOND_CROSSING_DETECTION
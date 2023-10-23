SUBROUTINE TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
IMPLICIT NONE
DOUBLE PRECISION P1(3), P2(3), P3(3), R1, R2, R3, SOL1(3), SOL2(3), I
DOUBLE PRECISION  TEMP1(3), EX(3), EY(3), EZ(3), DUMMY, TEMP2(3), TEMP3(3), D, J, X, Y, Z, TEMP4
LOGICAL FTEST

! # Find the intersection of three spheres                 
! # P1,P2,P3 are the centers, r1,r2,r3 are the radii       
! # Implementaton based on Wikipedia Trilateration article.                              

FTEST=.FALSE.
TEMP1(1:3)=P2(1:3)-P1(1:3)
D=SQRT( TEMP1(1)**2+TEMP1(2)**2+TEMP1(3)**2 )
EX(1:3)=TEMP1(1:3)/D
TEMP2(1:3)=P3(1:3)-P1(1:3)
I=EX(1)*TEMP2(1)+EX(2)*TEMP2(2)+EX(3)*TEMP2(3)
TEMP3(1:3)=TEMP2(1:3)-I*EX(1:3)
DUMMY=SQRT( TEMP3(1)**2+TEMP3(2)**2+TEMP3(3)**2 )
EY(1:3)=TEMP3(1:3)/DUMMY
EZ(1)= EX(2)*EY(3)-EX(3)*EY(2)
EZ(2)=-EX(1)*EY(3)+EX(3)*EY(1)
EZ(3)= EX(1)*EY(2)-EX(2)*EY(1)
J=EY(1)*TEMP2(1)+EY(2)*TEMP2(2)+EY(3)*TEMP2(3)
X=(R1*R1 - R2*R2 + D*D) / (2.0D0*D)
Y=(R1*R1 - R3*R3 -2.0D0*I*X + I*I + J*J) / (2.0D0*J)
TEMP4=R1*R1 - X*X - Y*Y

! WRITE (*,'(A,9G15.5)') 'trilateration> p1, p2, p3: ',P1(1:3),P2(1:3),P3(1:3)
! WRITE (*,'(A,9G15.5)') 'trilateration> ex, ey, ez: ',EX(1:3),EY(1:3),EZ(1:3)
! WRITE (*,'(A,9G15.5)') 'trilateration> norms:      ',&
! & EX(1)**2+EX(2)**2+EX(3)**2,EY(1)**2+EY(2)**2+EY(3)**2,EZ(1)**2+EZ(2)**2+EZ(3)**2
! WRITE (*,'(A,9G15.5)') 'trilateration> r1, r2, r3: ',R1,R2,R3
! WRITE (*,'(A,9G15.5)') 'trilateration> X, Y, TEMP4:    ',X,Y,TEMP4

! PRINT *,'TEMP4=',TEMP4
!PRINT *,'TEMP4.LT.0.0D0=',TEMP4.LT.0.0D0

IF (TEMP4.LT.0.0D0) THEN
   FTEST=.TRUE.
   RETURN
ELSE
   FTEST=.FALSE.
   Z=SQRT(TEMP4)
   SOL1(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) + Z*EZ(1:3)
   SOL2(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) - Z*EZ(1:3)
!  PRINT *,'Z=',Z
ENDIF

END SUBROUTINE TRILATERATION 
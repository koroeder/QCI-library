SUBROUTINE MINMAXD2(D2,D1,DINT,DSQ2,DSQ1,G1,G2,G1INT,G2INT,NOINT,DEBUG,r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
USE KEY, ONLY : CHECKCONINT
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DINT
DOUBLE PRECISION G1(3),G2(3),G1INT(3),G2INT(3)
DOUBLE PRECISION DSQ2, DSQ1, DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq
DOUBLE PRECISION r1amr1bdr2amr2b, r1amr1bdr2amr2bsq, DUMMY, DUMMY2
LOGICAL NOINT, DEBUG

!
! Is there an internal extremum?
! THIS TEST IS NOT NEEDED IF CHECKCONINT IS FALSE. SHOULD SKIP.
!
! IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN ! now done in calling routine
!
DUMMY=(DSQ1-r1amr1bdr2amr2b)/r1apr2bmr2amr1bsq
NOINT=.TRUE.
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
G1INT(1:3)=0.0D0
G2INT(1:3)=0.0D0
D2=SQRT(DSQ2)
D1=SQRT(DSQ1)
DSQI=1.0D10
DINT=1.0D10
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
   DSQI=MAX(-DUMMY2/r1apr2bmr2amr1bsq,0.0D0)
   DUMMY=r1apr2bmr2amr1bsq**2
   DINT=SQRT(DSQI)
   IF (DINT.LE.0.0D0) THEN
      NOINT=.TRUE.
   ELSE
     DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
     DUMMY=DUMMY*DINT ! Convert derivatives of distance^2 to derivative of distance.
     G1INT(1:3)= (DUMMY2*(G1(1:3) - G2(1:3)) + r1apr2bmr2amr1bsq*(G1(1:3)*DSQ2 -G2(1:3)*r1amr1bdr2amr2b))/DUMMY
     G2INT(1:3)= (DUMMY2*(G2(1:3) - G1(1:3)) + r1apr2bmr2amr1bsq*(G2(1:3)*DSQ1 -G1(1:3)*r1amr1bdr2amr2b))/DUMMY
   ENDIF
ENDIF
!
! Convert derivatives of distance^2 to derivative of distance.
! We have cancelled a factor of two above and below!
!
G2(1:3)=G2(1:3)/D2
G1(1:3)=G1(1:3)/D1
! IF (.NOT.NOINT) THEN
!    G1INT(1:3)=G1INT(1:3)/DINT
!    G2INT(1:3)=G2INT(1:3)/DINT
! ENDIF

END SUBROUTINE MINMAXD2
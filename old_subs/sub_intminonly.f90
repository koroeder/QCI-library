SUBROUTINE INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,NOINT)
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,DUMMY
DOUBLE PRECISION DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq, r1amr1bdr2amr2b, r1amr1bdr2amr2bsq
LOGICAL NOINT
!
! Is there an internal extremum?
!
! PRINT '(A,4G20.10)','r1ax,r1bx,r2ax,r2bx=',r1ax,r1bx,r2ax,r2bx
! PRINT '(A,G20.10)','(r1ax-r1bx-r2ax+r2bx)**2=',(r1ax-r1bx-r2ax+r2bx)**2
! PRINT '(A,4G20.10)','r1ay,r1by,r2ay,r2by=',r1ay,r1by,r2ay,r2by
! PRINT '(A,G20.10)','(r1ay-r1by-r2ay+r2by)**2=',(r1ay-r1by-r2ay+r2by)**2
! PRINT '(A,4G20.10)','r1az,r1bz,r2az,r2bz=',r1az,r1bz,r2az,r2bz
! PRINT '(A,G20.10)','(r1az-r1bz-r2az+r2bz)**2=',(r1az-r1bz-r2az+r2bz)**2
r1apr2bmr2amr1bsq=(r1ax-r1bx-r2ax+r2bx)**2+(r1ay-r1by-r2ay+r2by)**2+(r1az-r1bz-r2az+r2bz)**2
NOINT=.TRUE.
DINT=1.0D100
IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN
   RETURN ! just to skip the internal solution
ELSE
   DUMMY=((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+ &
 &      (r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
ENDIF
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2b=(r1ax-r1bx)*(r2ax-r2bx)+(r1ay-r1by)*(r2ay-r2by)+(r1az-r1bz)*(r2az-r2bz)
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   r1amr1bsq=(r1ax - r1bx)**2 + (r1ay - r1by)**2 + (r1az - r1bz)**2
   r2amr2bsq=(r2ax - r2bx)**2 + (r2ay - r2by)**2 + (r2az - r2bz)**2
   DSQI=MAX((-r1amr1bdr2amr2bsq + r1amr1bsq*r2amr2bsq)/r1apr2bmr2amr1bsq,0.0D0)
   DINT=SQRT(DSQI)
ENDIF

RETURN

END SUBROUTINE INTMINONLY
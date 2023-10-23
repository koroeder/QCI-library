!
! Neighbour list for repulsions to reduce cost of constraint potential.
!
SUBROUTINE CHECKREP(INTIMAGE,XYZ,NOPT,NNSTART,NSTART)
USE KEY,ONLY : NREPI, NREPJ, NREPCUT, NNREPULSIVE, NREPULSIVE, REPI, REPJ, REPCUT, CHECKREPCUTOFF, &
  &            NNREPULSIVE, intconstraintrep
USE COMMONS, ONLY : DEBUG
USE PORFUNCS
IMPLICIT NONE
INTEGER JJ, KK, INTIMAGE, NOPT, NI, NJ, NNSTART, NSTART!, NI1, NJ1, NI2, NJ2
DOUBLE PRECISION LDIST, XYZ(NOPT*(INTIMAGE+2)),COMPARE
!DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN
!LOGICAL NOINT

IF (INTCONSTRAINTREP.EQ.0) THEN
   NNREPULSIVE=0
   RETURN
ENDIF

NNREPULSIVE=NNSTART
DO JJ=NSTART,NREPULSIVE
   COMPARE=(CHECKREPCUTOFF*REPCUT(JJ))**2
   NI=REPI(JJ)
   NJ=REPJ(JJ)
   DO KK=1,INTIMAGE+2 ! first check for standard distances within threshold
      LDIST=(XYZ((KK-1)*NOPT+3*(NI-1)+1)-XYZ((KK-1)*NOPT+3*(NJ-1)+1))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+2)-XYZ((KK-1)*NOPT+3*(NJ-1)+2))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+3)-XYZ((KK-1)*NOPT+3*(NJ-1)+3))**2
      IF (LDIST.LT.COMPARE) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=NI
         NREPJ(NNREPULSIVE)=NJ
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
!        IF ((REPI(JJ).EQ.2024).OR.(REPJ(JJ).EQ.2024)) THEN
!           WRITE(*,'(A)') 'checkrep> KK,JJ,ldist < compare : active'
!        ENDIF
         GOTO 246
      ENDIF
   ENDDO 
!
! We don't check for internal minima in repulsions in congrad now unless both distances are
! within threshold.
!

246 CONTINUE
ENDDO
IF (DEBUG) WRITE(*,'(A,2I8)') ' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE

END SUBROUTINE CHECKREP
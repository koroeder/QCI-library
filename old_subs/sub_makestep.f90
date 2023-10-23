SUBROUTINE MAKESTEP(NITERDONE,POINT,DIAG,INTIMAGE,SEARCHSTEP,G,GTMP,STP,GDIF,NPT,D,RHO1,ALPHA)
USE KEY, ONLY : INTMUPDATE, INTDGUESS
USE COMMONS, ONLY: NATOMS
IMPLICIT NONE
INTEGER NITERDONE, POINT, BOUND, NPT, D, CP, INTIMAGE, I
DOUBLE PRECISION DIAG(3*NATOMS*INTIMAGE),SEARCHSTEP(0:INTMUPDATE,(3*NATOMS)*INTIMAGE),G((3*NATOMS)*INTIMAGE), &
  &  GTMP(3*NATOMS*INTIMAGE), GNORM, STP(3*NATOMS*INTIMAGE), YS, GDIF(0:INTMUPDATE,(3*NATOMS)*INTIMAGE), YY, &
  &  SQ, YR, BETA
DOUBLE PRECISION, DIMENSION(INTMUPDATE)     :: RHO1,ALPHA
SAVE

MAIN: IF (NITERDONE==1) THEN
     POINT = 0
     DIAG(1:D)=INTDGUESS
     SEARCHSTEP(0,1:D)= -G(1:D)*INTDGUESS            ! NR STEP FOR DIAGONAL INVERSE HESSIAN
     GTMP(1:D)        = SEARCHSTEP(0,1:D)
     GNORM            = MAX(SQRT(DOT_PRODUCT(G(1:D),G(1:D))),1.0D-100)
     STP(1:D)         = MIN(1.0D0/GNORM, GNORM) ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
ELSE MAIN
     BOUND=NITERDONE-1
     IF (NITERDONE.GT.INTMUPDATE) BOUND=INTMUPDATE
     YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
     IF (YS==0.0D0) YS=1.0D0
    
! Update estimate of diagonal inverse Hessian elements.
! We divide by both YS and YY at different points, so they had better not be zero!

     YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
     IF (YY==0.0D0) YY=1.0D0
!    DIAG = ABS(YS/YY)
     DIAG(1) = YS/YY
      
! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, 
! "Updating quasi-Newton matrices with limited storage",
! Mathematics of Computation, Vol.35, No.151, pp. 773-782

     CP= POINT; IF (POINT==0) CP = INTMUPDATE
     RHO1(CP)=1.0D0/YS
     GTMP(1:D) = -G(1:D)
     CP= POINT 
                   
     DO I= 1,BOUND 
          CP = CP - 1; IF (CP == -1) CP = INTMUPDATE - 1
          SQ= DOT_PRODUCT( SEARCHSTEP(CP,1:D),GTMP(1:D) )
          ALPHA(CP+1) = RHO1(CP+1) * SQ
          GTMP(1:D)        = -ALPHA(CP+1)*GDIF(CP,1:D) + GTMP(1:D)
     ENDDO
              
     GTMP(1:D)=DIAG(1)*GTMP(1:D)

     DO I=1,BOUND
          YR= DOT_PRODUCT( GDIF(CP,1:D) , GTMP )
          BETA= RHO1(CP+1)*YR
          BETA= ALPHA(CP+1)-BETA
!         WRITE(*,'(A,I8,4G20.10)') 'makestep> I,YR,BETA,RHO1,ALPHA=',I,YR,BETA,RHO1(CP+1),ALPHA(CP+1)
          GTMP(1:D) = BETA*SEARCHSTEP(CP,1:D) + GTMP(1:D)
          CP=CP+1
!         IF (CP==M) CP=0
          IF (CP==INTMUPDATE) CP=0
     ENDDO
              
     STP(1:D) = 1.0D0
ENDIF MAIN

!  Store the new search direction
IF (NITERDONE.GT.1) SEARCHSTEP(POINT,1:D)=GTMP(1:D)

END SUBROUTINE MAKESTEP
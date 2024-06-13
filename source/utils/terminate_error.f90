! Subroutine that is called when an error has occurred
! The error and debugging information should be provided before this routine is called.
! It will not give any specific information on the error, but clean up and stop the programm
! Note: We do not call this routine, if we completed the interpolation as it uses a STOP statement.
SUBROUTINE INT_ERR_TERMINATE()
   USE MOD_TERMINATE, ONLY: FINISH_QCI
   IMPLICIT NONE

   WRITE(*,*) " terminateQCI> An error was encountered - QCI is terminated"
   WRITE(*,*) " terminateQCI> Deallocating variables"
   CALL FINISH_QCI()
   STOP
END SUBROUTINE INT_ERR_TEMRINATE
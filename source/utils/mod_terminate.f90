! Module dealing with the temrination of QCI
MODULE MOD_TERMINATE

   CONTAINS
      ! Clean up routine that deallocates everything that has been allocated at some point
      ! All DEALLOC routines use protected IF (ALLOCATED(...)) DEALLOCATE(...), so no issue should arise here
      ! New modules/additional variables should be tidied in the same way
      SUBROUTINE FINISH_QCI()
         USE INTERPOLATION_KEYS, ONLY: DEALLOC_INTERPOLATIONS_VARS, DEALLOC_STEPTAKING
         USE CHIRALITY, ONLY: DEALLOC_CHIR_INTERNALS
         USE CONGEOM, ONLY: DEALLOC_CONGEOM, DEALLOC_ADDCONST
         USE QCICONSTRAINTS, ONLY: DEALLOC_CONSTR
         USE MOD_FREEZE, ONLY: DEALLOC_FREEZE
         USE MOD_INTCOORDS, ONLY: DEALLOC_INTCOORDS
         USE QCI_LINEAR, ONLY: DEALLOC_QCI_LINEAR
         USE REPULSION, ONLY: DEALLOC_REP_VARS
         USE ADDATOM, ONLY: DEALLOC_ADDATOM
         USE MOD_INTCOORDS, ONLY: XSTART, XFINAL

         IMPLICIT NONE
         CALL DEALLOC_INTERPOLATIONS_VARS()
         CALL DEALLOC_STEPTAKING()
         CALL DEALLOC_CHIR_INTERNALS()
         CALL DEALLOC_CONGEOM()
         CALL DEALLOC_ADDCONST()
         CALL DEALLOC_CONSTR()
         CALL DEALLOC_FREEZE()
         CALL DEALLOC_INTCOORDS()
         CALL DEALLOC_QCI_LINEAR()
         CALL DEALLOC_REP_VARS()
         CALL DEALLOC_ADDATOM()
         ! XSTART and XFINAL are special cases, and are allocated outside of a specific ALLOC routine
         ! Hence they need to be deallocated here
         IF (ALLOCATED(XSTART)) DEALLOCATE(XSTART)
         IF (ALLOCATED(XFINAL)) DEALLOCATE(XFINAL)
      END SUBROUTINE FINISH_QCI

      ! Subroutine that is called when an error has occurred
      ! The error and debugging information should be provided before this routine is called.
      ! It will not give any specific information on the error, but clean up and stop the programm
      ! Note: We do not call this routine, if we completed the interpolation as it uses a STOP statement.
      SUBROUTINE INT_ERR_TERMINATE()
         IMPLICIT NONE

         WRITE(*,*) " terminateQCI> An error was encountered - QCI is terminated"
         WRITE(*,*) " terminateQCI> Deallocating variables"
         CALL FINISH_QCI()
         STOP

      END SUBROUTINE INT_ERR_TEMRINATE

END MODULE MOD_TERMINATE

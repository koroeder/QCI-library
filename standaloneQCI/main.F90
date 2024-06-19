PROGRAM QCI_STANDALONE
   USE QCI_INTERFACE
   IMPLICIT NONE
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)
   
   INTEGER :: NARGS, NATOMS, NIMAGES
   INTEGER, PARAMETER :: XUNIT = 11 
   INTEGER, PARAMETER :: STDOUT = 6 
   CHARACTER(LEN=10) ::NATSDUMMY
   CHARACTER(LEN=50) :: PARAMFILE !Name of parameter file
   REAL(KIND = REAL64), ALLOCATABLE :: XYZ(:), XS(:), XF(:)
   LOGICAL :: COMPLETED
   INTEGER :: I

   ! check number of arguments
   NARGS = COMMAND_ARGUMENT_COUNT()
   ! We expect two arguments, the number of atoms and the parameter file name
   IF (NARGS.EQ.2) THEN
      CALL GET_COMMAND_ARGUMENT(1, NATSDUMMY)
      READ(NATSDUMMY,*) NATOMS
      CALL GET_COMMAND_ARGUMENT(2, PARAMFILE) 
   ELSE
      WRITE(STDOUT,'(A,I4)') "Expecting two arguments, but got ", NARGS
      STOP
   END IF
   ! Read in coordinates
   ALLOCATE(XS(3*NATOMS),XF(3*NATOMS))
   OPEN(XUNIT, FILE="start")
   READ(XUNIT, *) (XS(I), I=1,3*NATOMS)
   CLOSE(XUNIT)
   OPEN(XUNIT, FILE="finish")
   READ(XUNIT, *) (XF(I), I=1,3*NATOMS)
   CLOSE(XUNIT)   
   WRITE(STDOUT,'(A)') "Read coordinates for endpoints"
   CALL PASS_DATA(NATOMS, XS, XF)
   WRITE(STDOUT,'(A)') "Initialising QCI ..."
   CALL QCI_INITIALISE(PARAMFILE, .TRUE.)
   WRITE(STDOUT,'(A)') " "
   WRITE(STDOUT,'(A)') "Running interpolation ..."
   FLUSH(STDOUT) 
   CALL QCI_INTERPOLATION()
   WRITE(STDOUT,'(A)') " "
   WRITE(STDOUT,'(A)') "Getting results ..."
   CALL GET_QCI_INFO(NIMAGES,COMPLETED)
   IF (COMPLETED) THEN
       WRITE(STDOUT,'(A)') "Quasi-continuous interpolation band converged"
   ELSE
       WRITE(STDOUT,'(A)') "Quasi-continuous interpolation band did not converge"   
   END IF   
   WRITE(STDOUT,'(A,I6,A)') "QCI band has ",NIMAGES, " images"
   ALLOCATE(XYZ((3*NATOMS)*(NIMAGES+2)))
   CALL GET_QCI_INTERPOLATION(NATOMS,NIMAGES,XYZ)
   !potentially do something here
   !...
   CALL QCI_TERMINATE()
END PROGRAM QCI_STANDALONE

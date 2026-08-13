PROGRAM QCI_STANDALONE
   USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: COMPILER_VERSION
   USE QCI_INTERFACE
   USE TIMER_MODULE
   IMPLICIT NONE
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)
   
   INTEGER :: NARGS, NATOMS, NIMAGES
   INTEGER :: NLINES
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
   IF (NARGS.EQ.1) THEN
      CALL GET_COMMAND_ARGUMENT(1, PARAMFILE) 
   ELSE IF (NARGS.EQ.2) THEN
      WRITE(*,'(A)')   " WARNING: Usage has changed, please change input script... "
      WRITE(*,'(A)')   " Usage: ./QCI <params_file> > output"
      CALL GET_COMMAND_ARGUMENT(1, NATSDUMMY)
      READ(NATSDUMMY,*) NATOMS
      CALL GET_COMMAND_ARGUMENT(2, PARAMFILE) 
   ELSE IF (NARGS.EQ.0) THEN
      WRITE(*,'(A)')   " QCI - Standalone Quasi-Continuous Interpolation"
      WRITE(*,'(A)')   " Usage: ./QCI <params_file> > output"
      WRITE(*,'(A,A,A,A)') ' Built on:    ', __DATE__, ' ', __TIME__
      WRITE(*,'(A,A)')   ' Compiler:    ', TRIM(COMPILER_VERSION())
      STOP
   ELSE
      WRITE(STDOUT,'(A,I4)') "Expecting two arguments, but got ", NARGS
      STOP
   END IF

   !Start the clock
   CALL TIMER_START()

   ! Read in coordinates
   
   !start
   OPEN(XUNIT, FILE="start")
   CALL COUNT_LINES(XUNIT, "start", NLINES)
   NATOMS = NLINES
   ALLOCATE(XS(3*NATOMS),XF(3*NATOMS))
   READ(XUNIT, *) (XS(I), I=1,3*NATOMS)
   CLOSE(XUNIT)
   
   !finish
   OPEN(XUNIT, FILE="finish")
   CALL COUNT_LINES(XUNIT, "finish", NLINES)
   IF (NATOMS.NE.NLINES) THEN
      WRITE(STDOUT, '(A)') "Number of atoms in start does not match number of atoms in finish."
      WRITE(STDOUT, '(A)') "Terminating ..."
      
      DEALLOCATE(XS,XF)
      STOP
   END IF

   READ(XUNIT, *) (XF(I), I=1,3*NATOMS)
   CLOSE(XUNIT)   

   WRITE(STDOUT,'(A)') "QCI - Standalone Quasi-Continious Interpolation"
   
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
   DEALLOCATE(XS,XF,XYZ)
   
   !Stop the timer and output time elapsed
   CALL TIMER_STOP()


   CONTAINS

   SUBROUTINE COUNT_LINES(unit, fname, n_lines)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: unit
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(OUT) :: n_lines
      INTEGER :: ios, ios2
      CHARACTER(LEN=512) :: line
      REAL(KIND=REAL64) :: x, y, z, extra

      n_lines = 0
      DO
         READ(unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0) EXIT
         n_lines = n_lines + 1

         ! Must have at least 3 numbers
         READ(line, *, IOSTAT=ios) x, y, z
         IF (ios /= 0) THEN
            WRITE(*,'(A,A,A,I0,A)') "Error in ", TRIM(fname), &
                 ": line ", n_lines, " has fewer than 3 numbers"
            STOP
         END IF

         ! Must not have more than 3
         READ(line, *, IOSTAT=ios2) x, y, z, extra
         IF (ios2 == 0) THEN
            WRITE(*,'(A,A,A,I0,A)') "Error in ", TRIM(fname), &
                 ": line ", n_lines, " has more than 3 numbers"
            STOP
         END IF
      END DO

      REWIND(unit)
   END SUBROUTINE COUNT_LINES

  

END PROGRAM QCI_STANDALONE

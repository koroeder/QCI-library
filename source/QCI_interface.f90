! This is the interface to the QCI library.
! External software should only interact with this interface and not with any other code.
! The interface subroutines should be used to pass data back and forth and to call routines.
! Three module variables control the flow.
! The correct order of use is 
! 1. a call to PASS_DATA
! 2. a call to QCI_INITIALISE
! 3. a call to QCI_INTERPOLATION
! 4. calls to GET_QCI_INFO and GET_QCI_INTERPOLATION
! The GET_QCI_INFO subroutine should be used to get the number of QCI images in the final interpolation
! This can then be used to allocate variables before the band coordinates are obtained through GET_QCI_INTERPOLATION
! For double precision here, the type R64 should be used.
! This type matches the type used throughout the library.
MODULE QCI_INTERFACE
   INTEGER, PARAMETER  :: R64 = SELECTED_REAL_KIND(15, 307) 
   LOGICAL :: SETUP_COMPLETE = .FALSE.
   LOGICAL :: INITIALISED = .FALSE.
   LOGICAL :: INTERPOLATION_COMPLETE = .FALSE.
   CONTAINS

      ! Subroutine to pass the number of atoms and end point coordinates to the QCI library
      ! The call to setup endpoints only allocates the starting point coordinate arrays and sets NATOMS
      SUBROUTINE PASS_DATA(NATS, XS, XF)
         USE MOD_INTCOORDS, ONLY: SETUP_ENDPOINTS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATS
         REAL(KIND = R64), INTENT(IN) :: XS(3*NATS), XF(3*NATS)

         CALL SETUP_ENDPOINTS(NATS)
         XSTART(1:3*NATS) = XS(1:3*NATS)
         XFINAL(1:3*NATS) = XF(1:3*NATS)
         SETUP_COMPLETE = .TRUE.
      END SUBROUTINE PASS_DATA

      ! Subroutine to set up the interpolation
      ! The two input variables are the name of the parameter file that is used to set up the interpoaltion
      ! and a logical switch whether the end points should be permutationally aligned.
      ! Likely, this option should always be used.
      ! The subroutine calls QCI_INIT, which handles the parsing of information and the correct setup of the interpolation.
      SUBROUTINE QCI_INITIALISE(PARAMFILE, ALIGNT)
         USE QCISETUP, ONLY: QCI_INIT
         IMPLICIT NONE
         CHARACTER(*), INTENT(IN) :: PARAMFILE
         LOGICAL, INTENT(IN) :: ALIGNT

         IF (.NOT.SETUP_COMPLETE) THEN
            WRITE(*,*) " qci_initialise> PASS_DATA has not been called - setup not complete, cannot initialise"
            RETURN
         END IF
         CALL QCI_INIT(PARAMFILE, ALIGNT)
         INITIALISED = .TRUE.
      END SUBROUTINE QCI_INITIALISE

      SUBROUTINE QCI_INTERPOLATION()
         USE QCIINTERPOLATION, ONLY: RUN_QCI_INTERPOLATION
         IMPLICIT NONE

         IF (.NOT.INITIALISED) THEN
            WRITE(*,*) " qci_interpoaltion> QCI has not been initialised correctly, cannot run interpolation"
            RETURN
         END IF

         CALL RUN_QCI_INTERPOLATION()
         INTERPOLATION_COMPLETE = .TRUE.
      END SUBROUTINE QCI_INTERPOLATION

      SUBROUTINE GET_QCI_INFO(NQCIIMAGES,COMPLETED)
         USE QCIKEYS, ONLY: NIMAGES
         USE QCIINTERPOLATION, ONLY: QCICOMPLETE
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: NQCIIMAGES
         LOGICAL, INTENT(OUT) :: COMPLETED

         IF (.NOT.INTERPOLATION_COMPLETE) THEN
            WRITE(*,*) " get_qci_info> QCI interpolation did not finish, cannot fetch data"
            RETURN
         END IF

         NQCIIMAGES = NIMAGES
         COMPLETED = QCICOMPLETE
      END SUBROUTINE GET_QCI_INFO

      SUBROUTINE GET_QCI_INTERPOLATION(XBAND)
         USE QCIKEYS, ONLY: NIMAGES, NATOMS        
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         REAL(KIND = R64), INTENT(OUT) :: XBAND((3*NATOMS)*(NIMAGES+2))

         IF (.NOT.INTERPOLATION_COMPLETE) THEN
            WRITE(*,*) " get_qci_interpolation> QCI interpolation did not finish, cannot fetch data"
            RETURN
         END IF

         XBAND = XYZ
      END SUBROUTINE GET_QCI_INTERPOLATION

      SUBROUTINE QCI_TERMINATE()
         USE MOD_TERMINATE, ONLY: FINISH_QCI

         CALL FINISH_QCI()

      END SUBROUTINE QCI_TERMINATE

END MODULE QCI_INTERFACE
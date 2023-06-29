MODULE QCIFILEHANDLER

   CONTAINS
      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)     
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
       
         FILE_UNIT = GETUNIT()
         IF (APPEND) THEN
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN', POSITION='APPEND')
         ELSE
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN')
         END IF   
      END SUBROUTINE FILE_OPEN

      INTEGER FUNCTION FILE_LENGTH(FILE_NAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         INTEGER                       :: FILE_UNIT
         INTEGER                       :: IO_STATUS

         CALL FILE_OPEN(FILE_NAME, FILE_UNIT, .FALSE.)
         FILE_LENGTH = 0
         DO
            READ(FILE_UNIT, *, IOSTAT=IO_STATUS)
            IF (IO_STATUS /= 0) EXIT
            FILE_LENGTH = FILE_LENGTH + 1
         ENDDO
         CLOSE(FILE_UNIT)
      END FUNCTION FILE_LENGTH 

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
         INTEGER :: UNITNUM
         ! start checking for available units > 103, to avoid system default units
         ! 100, 101 and 102 are stdin, stdout and stderr respectively.
         INUSE=.TRUE.
         UNITNUM=103

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT

END MODULE QCIFILEHANDLER
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

      function file_exists(filename) result(res)
         implicit none
         character(len=*),intent(in) :: filename
         logical                     :: res

         ! Check if the file exists
         inquire( file=trim(filename), exist=res )
      end function

      subroutine count_constraints(filename, nconstraints)
         use QCIKEYS, only: NATOMS
         implicit none
         character(len=*), intent(in)  :: filename
         integer,          intent(out) :: nconstraints

         integer, parameter :: line_len = 256
         character(len=line_len) :: line
         integer :: unit_no, ios, n, m, extra, line_num
         logical :: extra_present

         nconstraints = 0
         line_num     = 0

         open(newunit=unit_no, file=trim(filename), status='old', &
            action='read', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'ERROR> could not open file ', trim(filename)
            CALL INT_ERR_TERMINATE()
         end if

         do
            read(unit_no, '(A)', iostat=ios) line
            if (ios /= 0) exit          ! end of file (or read error)

            line_num = line_num + 1

            ! Guard against empty (or whitespace-only) lines
            if (len_trim(line) == 0) cycle

            ! Must contain exactly two integers "N M"
            read(line, *, iostat=ios) n, m
            if (ios /= 0) then
               write(*,*) 'ERROR> file' , trim(filename), ' line ', line_num, &
                           ' is not of the form "N M": ', trim(line)
               CALL INT_ERR_TERMINATE()
            end if

            ! Reject a line with a third number lurking on it
            extra_present = .false.
            read(line, *, iostat=ios) n, m, extra
            if (ios == 0) extra_present = .true.

            if (extra_present) then
               write(*,*) 'ERROR> file ', trim(filename), ' line ', line_num, &
                           ' has more than two numbers: ', trim(line)
               CALL INT_ERR_TERMINATE()
            end if

            ! Range check against NATOMS
            if (n < 1 .or. n > natoms .or. m < 1 .or. m > natoms) then
               write(*,*) 'ERROR> file ', trim(filename), ', line ', line_num, &
                           ' has atom index out of range (1..', natoms, &
                           '): ', trim(line)
               CALL INT_ERR_TERMINATE()
            end if

            nconstraints = nconstraints + 1
         end do

         close(unit_no)

      end subroutine count_constraints

END MODULE QCIFILEHANDLER
SUBROUTINE DOCHICHECK(NEIGHBOUR_COORDS,CENTRE_COORDS,XYZ,NUM_CHIRAL_CENTRES)
USE CHIRALITY
USE KEY, ONLY : ATOMACTIVE, INTIMAGE, ATOMSTORES, TWOD, RIGIDBODY, BULKT
USE COMMONS, ONLY: NATOMS, DEBUG
IMPLICIT NONE                
INTEGER J5,J3,J4,NCHANGE,J2, CHANGEAT(NATOMS), ACID
DOUBLE PRECISION NEIGHBOUR_COORDS(12), CENTRE_COORDS(3),XYZ((3*NATOMS)*(INTIMAGE+2)),COORDSA(3*NATOMS), COORDSB(3*NATOMS)
DOUBLE PRECISION DIST, DWORST, RMAT(3,3)
LOGICAL CHIRALSR, CHIRALSRP
INTEGER NUM_CHIRAL_CENTRES,ATOM_NUMBER
CHARACTER(LEN=5) ZUSE

IF (DEBUG) PRINT *,'DOING CHIRALCHECK NOW'
!       IF (DEBUG) WRITE(*,'(A)') ' intlbfgs> dump state before CHIRALCHECK index -4'
!       IF (DEBUG) CALL INTRWG2(NACTIVE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)

   IF (NUM_CHIRAL_CENTRES.EQ.-1) THEN
      PRINT '(A)', 'intlbfgs> ERROR *** NUM_CHIRAL_CENTRES is -1 *** forgotten QCIAMBER keyword?'
      STOP
   ENDIF

   chicheck: DO J5=1, NUM_CHIRAL_CENTRES
      ATOM_NUMBER=SR_ATOMS(J5, 1) 
!     WRITE(*,'(A,I6,A,I6,A,I6)') 'chiral centre ',J5,' is atom ',atom_number
      IF (.NOT.ATOMACTIVE(ATOM_NUMBER)) CYCLE chicheck
      DO J2=1,4
         IF (.NOT.ATOMACTIVE(SR_ATOMS(J5,J2))) CYCLE chicheck
      ENDDO

      DO J3=1,INTIMAGE+2

         CENTRE_COORDS(1)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+1)
         CENTRE_COORDS(2)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+2)
         CENTRE_COORDS(3)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+3)

         DO J4=1,4
            J2=sr_atoms(J5, J4 + 1)
            NEIGHBOUR_COORDS(3*(J4-1)+1)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+1)
            NEIGHBOUR_COORDS(3*(J4-1)+2)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+2)
            NEIGHBOUR_COORDS(3*(J4-1)+3)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+3)
         ENDDO

         CHIRALSR=CHIRALITY_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
!        WRITE(*,'(A,I6,I6,2L5)') 'image, atom, chirality, initial=',J3,atom_number,CHIRALSR,sr_states_initial(J5)
         IF (J3.EQ.1) CHIRALSRP=sr_states_initial(J5)
         IF (CHIRALSR.NEQV.CHIRALSRP) THEN
            WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> Atom ',atom_number,' image ',J3,&
     &          ' chirality CHANGED; use previous image coordinates'  
!            NLASTCHANGE=NITERDONE
!
! The idea is to change all the coordinates for active atoms in the same amino acid to the values
! in the frame before. We should align the replacements on the (sub)set they are replacing.
!
            ACID=ATOMSTORES(atom_number)

            NCHANGE=0
            CHANGEAT(1:NATOMS)=-1
            DO J4=1,NATOMS
               IF (ATOMSTORES(J4).NE.ACID) CYCLE
               IF (.NOT.ATOMACTIVE(J4)) CYCLE
               WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> Changing active atom ',J4,' image ',J3
               NCHANGE=NCHANGE+1
               CHANGEAT(NCHANGE)=J4
               COORDSA(3*(NCHANGE-1)+1:3*(NCHANGE-1)+3)=XYZ(3*NATOMS*(J3-1)+3*(J4-1)+1:3*NATOMS*(J3-1)+3*(J4-1)+3)
               COORDSB(3*(NCHANGE-1)+1:3*(NCHANGE-1)+3)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+1:3*NATOMS*(J3-2)+3*(J4-1)+3)

               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+1)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+1)
               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+2)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+2)
               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+3)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+3)

            ENDDO
!
! Align the replacement atoms with the original atoms as far as possible
!
            ZUSE='LA   '
            CALL NEWMINDIST(COORDSA,COORDSB,NCHANGE,DIST,BULKT,TWOD,ZUSE,.FALSE.,RIGIDBODY,DEBUG,RMAT,DWORST)

            DO J4=1,NCHANGE
               XYZ(3*NATOMS*(J3-1)+3*(CHANGEAT(J4)-1)+1:3*NATOMS*(J3-1)+3*(CHANGEAT(J4)-1)+3)=COORDSB(3*(J4-1)+1:3*(J4-1)+3)  
            ENDDO

         ENDIF
!        IF (J3.EQ.1) CHIRALSRP=CHIRALSR  ! just use result for fixed end point image 1
      ENDDO
   ENDDO chicheck

END SUBROUTINE DOCHICHECK
MODULE CHIRALITY
   USE QCIPREC
   USE QCIKEYS, ONLY: NATOMS
   USE AMBER_CONSTRAINTS, ONLY: NRES, NBOND, BONDS, RESNAMES, RES_START, RES_END, AMBER_NAMES, ELEMENT
   IMPLICIT NONE
   INTEGER, SAVE :: NCHIRAL = -1
   INTEGER :: NCSP2                                          !number of sp2 hybridised carbons
   INTEGER :: NDOUBLE                                        !number of double bonds
   INTEGER, ALLOCATABLE :: CSP2(:)                           !list of sp2 hybridised carbon
   INTEGER, ALLOCATABLE :: DOUBLEB(:,:)                                !array to store double bonds
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: CHIR_INFO   !array for chiral centres
   INTEGER, ALLOCATABLE, SAVE :: NBONDED(:)                  !number of bonds for each atom 
   INTEGER, ALLOCATABLE, SAVE :: NNEIGHBOURS(:)              !number of bonded atoms 
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: BONDEDATS   !atoms bonded by atomid (up to 6)
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: BONDEDELS   !atoms bonded by element (up to 6) 
   LOGICAL, ALLOCATABLE, SAVE :: POSSIBLE_CC(:)              !potential chiral centres
   INTEGER, ALLOCATABLE, SAVE :: POSSIBLE_IDS(:,:)           !ids of bonded atoms      

   CONTAINS
      SUBROUTINE DEALLOC_CHIR_INTERNALS()
         IF (ALLOCATED(DOUBLEB)) DEALLOCATE(DOUBLEB)
         IF (ALLOCATED(NBONDED)) DEALLOCATE(NBONDED)
         IF (ALLOCATED(NNEIGHBOURS)) DEALLOCATE(NNEIGHBOURS)
         IF (ALLOCATED(CHIR_INFO)) DEALLOCATE(CHIR_INFO) 
         IF (ALLOCATED(BONDEDATS)) DEALLOCATE(BONDEDATS)  
         IF (ALLOCATED(BONDEDELS)) DEALLOCATE(BONDEDELS)
         IF (ALLOCATED(POSSIBLE_CC)) DEALLOCATE(POSSIBLE_CC) 
         IF (ALLOCATED(POSSIBLE_IDS)) DEALLOCATE(POSSIBLE_IDS)  
         IF (ALLOCATED(CSP2)) DEALLOCATE(CSP2)
      END SUBROUTINE DEALLOC_CHIR_INTERNALS

      SUBROUTINE CHIRALITY_CHECK(XYZ)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES, ATOMS2RES
         USE AMBER_CONSTRAINTS, ONLY: RES_START, RES_END
         USE QCIMINDIST, ONLY: ALIGNXBTOA
         IMPLICIT NONE
         REAL(KIND = REAL64), INTENT(INOUT) :: XYZ((3*NATOMS)*(NIMAGES+2))
         REAL(KIND = REAL64) :: NEIGHBOUR_COORDS(12), CENTRE_COORDS(3) !intent needed?
         REAL(KIND = REAL64) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
         INTEGER :: CHIRALCENTRE
         INTEGER :: J1, J2, J3, J4, ATOMID, ACID
         INTEGER :: NCHANGE, CHANGEAT(NATOMS), THISIMAGE, PREVIMAGE, OFFSET, OFFSET2
         LOGICAL :: CENTREACTIVE
         LOGICAL :: THIS_SR, PREV_SR

         IF (DEBUG) WRITE(*,*) " chirality-check> Running check for chirality conservation across all images"

         IF (NCHIRAL.EQ.-1) THEN
            WRITE(*,*) " chirality-check> The number of chiral centres is -1. It looks like FIND_CHIRAL was not run. Is QCIAMBER set?"
            CALL INT_ERR_TERMINATE()
         ELSE IF (NCHIRAL.EQ.0) THEN
            WRITE(*,*) " chirality-check>  There are no chiral centres - check this is correct as CHECKCHIRAL is set to true."
            RETURN
         END IF

         PREV_SR = .FALSE.

         DO J1=1,NCHIRAL
            CENTREACTIVE = .TRUE.
            CHIRALCENTRE = CHIR_INFO(J1,1)
            !check whether the chiral centre and connected atoms are active
            IF (.NOT.ATOMACTIVE(CHIRALCENTRE)) CENTREACTIVE = .FALSE.
            DO J2=2,5
               IF (.NOT.ATOMACTIVE(CHIR_INFO(J1,J2))) CENTREACTIVE = .FALSE.
            END DO
            IF (.NOT.CENTREACTIVE) CYCLE
            !now apply the actual check            
            DO J3=1,NIMAGES+2
               CENTRE_COORDS(1) = XYZ(3*NATOMS*(J3-1)+3*CHIRALCENTRE-2)
               CENTRE_COORDS(2) = XYZ(3*NATOMS*(J3-1)+3*CHIRALCENTRE-1)
               CENTRE_COORDS(3) = XYZ(3*NATOMS*(J3-1)+3*CHIRALCENTRE)

               DO J4=1,4
                  ATOMID = CHIR_INFO(J1,J4+1)
                  NEIGHBOUR_COORDS(3*J4-2) = XYZ(3*NATOMS*(J3-1)+3*ATOMID-2)
                  NEIGHBOUR_COORDS(3*J4-1) = XYZ(3*NATOMS*(J3-1)+3*ATOMID-1)
                  NEIGHBOUR_COORDS(3*J4)   = XYZ(3*NATOMS*(J3-1)+3*ATOMID)
               END DO
               ! if this is the first image (i.e. the starting point), set SR for the current group
               IF (J3.EQ.1) THEN
                  PREV_SR = ASSIGNMENT_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
                  CYCLE
               END IF
               ! S/R assignment for current image
               THIS_SR = ASSIGNMENT_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
               IF (THIS_SR.NEQV.PREV_SR) THEN
                  WRITE(*,*) " chirality_check> Atom ", CHIRALCENTRE, " image ", J3, " chirality changed"
                  WRITE(*,*) "                  Using previous image coordinates."
                  ACID = ATOMS2RES(CHIRALCENTRE)
                  NCHANGE = 0
                  CHANGEAT(1:NATOMS) = -1
                  THISIMAGE = 3*NATOMS*(J3-1)
                  PREVIMAGE = 3*NATOMS*(J3-2)
                  DO J4=RES_START(ACID),RES_END(ACID)
                     IF (.NOT.ATOMACTIVE(J4)) CYCLE
                     ! WRITE(*,*) " chirality_check> Changing active atom ", J4, " in image ", J3
                     NCHANGE=NCHANGE+1
                     CHANGEAT(NCHANGE)=J4
                     OFFSET = 3*(NCHANGE-1)
                     OFFSET2 = 3*(J4-1)
                     COORDSA(OFFSET+1:OFFSET+3)=XYZ(THISIMAGE+OFFSET2+1:THISIMAGE+OFFSET2+3)
                     COORDSB(OFFSET+1:OFFSET+3)=XYZ(PREVIMAGE+OFFSET2+1:PREVIMAGE+OFFSET2+3)
                  END DO

                  ! Align the replacement atoms with the original atoms as far as possible
                  CALL ALIGNXBTOA(COORDSA, COORDSB, NCHANGE)
                  DO J4=1,NCHANGE
                     OFFSET = 3*(CHANGEAT(J4)-1)
                     OFFSET2 = 3*(J4-1)
                     XYZ(THISIMAGE+OFFSET+1:THISIMAGE+OFFSET+3)=COORDSB(OFFSET2+1:OFFSET2+3)  
                  ENDDO

               END IF
            END DO
         END DO
      END SUBROUTINE CHIRALITY_CHECK


      ! Works out whether a molecule is left- or right-handed (i.e. S or R)
      ! This is calculated by applying the CIP rules to the system and working out the
      ! dihedral angle between the atoms formed by:
      ! 1 - centre - 4 - 2
      LOGICAL FUNCTION ASSIGNMENT_SR(COORDS,CENTRE) RESULT(RIGHT_HANDED)
         USE HELPER_FNCTS, ONLY: DIHEDRAL
         IMPLICIT NONE
         REAL(KIND = REAL64) :: COORDS(12)
         REAL(KIND = REAL64) :: CENTRE(3)
         REAL(KIND = REAL64) :: DIHX(12)
         REAL(KIND = REAL64) :: ANGLE

         DIHX(1:3) = COORDS(1:3)
         DIHX(4:6) = CENTRE(1:3)
         DIHX(7:9) = COORDS(10:12)
         DIHX(10:12) = COORDS(4:6)

         ANGLE = DIHEDRAL(DIHX)

         RIGHT_HANDED = (ANGLE.GT.0.0D0)
      END FUNCTION ASSIGNMENT_SR


      SUBROUTINE FIND_CHIRAL_CENTRES()
         USE QCIKEYS, ONLY: ATOMS2RES
         USE COMPARELIST
         USE AMBER_CONSTRAINTS, ONLY: AMBER_NAMES, RESTYPE
         IMPLICIT NONE
         INTEGER :: J1, J2, ATOMID, J3, ATOMID2, VALENCY, GHOSTEL, J4
         INTEGER :: DUMMY_ELS(4), DUMMY_ELS_2(4), DUMMY_ATS(4), DUMMY_ATS_2(4)
         INTEGER :: PRIO4(4), PRIO3(3), FINALPRIO(4)
         INTEGER :: LIST1(4), LIST2(4), LIST3(4), LIST4(4)
         LOGICAL :: CHIRALT, EQ12, EQ23, EQ34, GHOSTS(4), TRIPLE, QUADRUPLE
         LOGICAL :: GT12T,ID12T,GT23T,ID23T,GT34T,ID34T,CHECKUNIQUE(4),IDFOUND
         INTEGER :: DUMMY_CHIR_INFO(NATOMS,5)     !array for chiral centres
         INTEGER :: C3AT1
         !get the number of atoms bonded to each atom 
         !anything with four bonds is considered a possible chiral centre
         CALL NBONDED_ATOMS()
         !need to add ghost atoms here, not later (too complicate otherwise)   
         CALL CREATE_GHOST_ATOMS()
         !sort list of bonded atoms by element (highest priority first)
         CALL SORT_BONDEDATS()
         !when we compare paths then we need a termination check to find the end of paths as well
         !remove any centre from the list with more than one Hydrogen attached
         CALL DISCOUNT_H()
         !now check all possible centres
         NCHIRAL=0      
         DO J1=1,NATOMS
            IF (POSSIBLE_CC(J1)) THEN
               CHIRALT = .FALSE.   
               !first copy the data for the directly attached atoms
               DUMMY_ELS(1:4) = BONDEDELS(J1,1:4)
               DUMMY_ATS(1:4) = BONDEDATS(J1,1:4)
               !first test: are all neighbours different elements? 
               CALL TEST_NEIGHBOURS(J1,CHIRALT)
               IF (CHIRALT) THEN
                  FINALPRIO(1:4) = (/1,2,3,4/)
               !now we need to compare paths, this means we need to compare the neighbouring lists for identical atoms  
               ELSE
                  !as the data is sorted, we have four cases:
                  !a) a pair of identical atoms        
                  !b) a triplet of identical atoms     
                  !c) a quadruplet of identical atoms   
                  !d) two pairs of identical atoms     
                  EQ12 = .FALSE.
                  EQ23 = .FALSE.
                  EQ34 = .FALSE.
                  IF (DUMMY_ELS(1).EQ.DUMMY_ELS(2)) EQ12=.TRUE.
                  IF (DUMMY_ELS(2).EQ.DUMMY_ELS(3)) EQ23=.TRUE.
                  IF (DUMMY_ELS(3).EQ.DUMMY_ELS(4)) EQ34=.TRUE.
                  ! before we go for comparisons, let's figure out how many comparisons we need
                  TRIPLE=.FALSE.
                  QUADRUPLE=.FALSE.
                  LIST1(1:4) = BONDEDELS(DUMMY_ATS(1),1:4)
                  LIST2(1:4) = BONDEDELS(DUMMY_ATS(2),1:4)
                  LIST3(1:4) = BONDEDELS(DUMMY_ATS(3),1:4)
                  LIST4(1:4) = BONDEDELS(DUMMY_ATS(4),1:4)
                  IF (EQ12.AND.EQ23.AND.EQ34) THEN
                     QUADRUPLE=.TRUE.
                  ELSE IF ((EQ12.AND.EQ23).OR.(EQ23.AND.EQ34)) THEN
                     TRIPLE=.TRUE.
                  ENDIF
                  IF ((.NOT.TRIPLE).AND.(.NOT.QUADRUPLE)) THEN
                     IF (EQ12) CALL COMPARE_2LISTS(LIST1,LIST2,4,ID12T,GT12T)                   
                     IF (EQ23) CALL COMPARE_2LISTS(LIST2,LIST3,4,ID23T,GT23T)                 
                     IF (EQ34) CALL COMPARE_2LISTS(LIST3,LIST4,4,ID34T,GT34T)                 
                  ELSE IF (TRIPLE) THEN
                     CALL COMPARE_3LISTS(LIST1,LIST2,LIST3,4,PRIO3)               
                  ELSE IF (QUADRUPLE) THEN
                     CALL COMPARE_4LISTS(LIST1,LIST2,LIST3,LIST4,4,PRIO4)
                  ENDIF
                  ! we now have the priority information based on the bound atoms 
                  ! of the directly attached atoms
                  ! there is one case in nucleic acids where this isn't sufficient
                  ! we deal with this later
                  ! now create a priority list for all four
                  ! for four identical atoms, just copy the result
                  IF (QUADRUPLE) THEN
                     FINALPRIO(1:4) = PRIO4(1:4)
                  ! for a triple assign priorities based on 
                  ELSE IF (TRIPLE) THEN
                     FINALPRIO(1:3) = PRIO3(1:3)             
                     IF (EQ12.AND.EQ23) THEN
                        FINALPRIO(4) = 4
                     ELSE
                        FINALPRIO(4) = 0
                        DO J2=1,4
                          FINALPRIO(J2) = FINALPRIO(J2) + 1
                        ENDDO                  
                     ENDIF
                  ELSE IF (EQ12.AND.EQ34) THEN
                     IF (ID12T) THEN
                        FINALPRIO(1:2) = (/1,1/)
                     ELSE
                        IF (GT12T) THEN
                           FINALPRIO(1:2) = (/1,2/)
                        ELSE
                           FINALPRIO(1:2) = (/2,1/)
                        ENDIF
                     ENDIF
                     IF (ID34T) THEN
                        FINALPRIO(3:4) = (/3,3/)
                     ELSE
                        IF (GT34T) THEN 
                           FINALPRIO(3:4) = (/3,4/)
                        ELSE
                           FINALPRIO(3:4) = (/4,3/) 
                        ENDIF                                         
                     ENDIF
                  ELSE
                     IF (EQ12) THEN
                        IF (ID12T) THEN
                           FINALPRIO(1:4) = (/1,1,2,3/)
                        ELSE
                           IF (GT12T) THEN
                              FINALPRIO(1:4) = (/1,2,3,4/)
                           ELSE
                              FINALPRIO(1:4) = (/2,1,3,4/)
                           ENDIF
                        ENDIF
                     ELSE IF (EQ23) THEN
                        IF (ID23T) THEN
                           FINALPRIO(1:4) = (/1,2,2,3/)
                        ELSE
                           IF (GT23T) THEN
                              FINALPRIO(1:4) = (/1,2,3,4/)
                           ELSE
                              FINALPRIO(1:4) = (/1,3,2,4/)
                           ENDIF
                        ENDIF                
                     ELSE IF (EQ34) THEN
                        IF (ID34T) THEN
                           FINALPRIO(1:4) = (/1,2,3,3/)
                        ELSE
                           IF (GT34T) THEN
                              FINALPRIO(1:4) = (/1,2,3,4/)
                           ELSE
                              FINALPRIO(1:4) = (/1,2,4,3/)
                           ENDIF
                        ENDIF     
                     ENDIF                              
                  ENDIF !end creating priority list
                  ! now we need to check if the list contains four different numbers.
                  ! if it does we have a chiral centre
                  CHECKUNIQUE(:) = .FALSE.
                  IDFOUND = .FALSE.
                  DO J3=1,4
                     !duplicate entry
                     IF (CHECKUNIQUE(FINALPRIO(J3))) THEN
                        IDFOUND=.TRUE.
                        EXIT
                     ELSE
                        CHECKUNIQUE(FINALPRIO(J3)) = .TRUE.
                     ENDIF
                  ENDDO
                  IF (.NOT.IDFOUND) THEN
                     CHIRALT=.TRUE.
                  ! in RNA we have one issue - for the C3' C2' and C4' 
                  ! can only be resolved with an additional list
                  ! in other nucleic acids we can, so if the atom name is C3', we need one more step
                  ELSE IF ((AMBER_NAMES(J1).EQ."C3'").AND.(RESTYPE(ATOMS2RES(J1)).EQ."RNA")) THEN
                     C3AT1 = DUMMY_ATS(2)
                     IF (AMBER_NAMES(C3AT1).EQ."C2'") THEN
                        FINALPRIO(3) = FINALPRIO(3) + 1
                        FINALPRIO(4) = FINALPRIO(4) + 1
                        CHIRALT=.TRUE.
                     ELSE IF (AMBER_NAMES(C3AT1).EQ."C4'") THEN 
                        FINALPRIO(2) = FINALPRIO(2) + 1
                        FINALPRIO(4) = FINALPRIO(4) + 1  
                        CHIRALT=.TRUE. 
                     ENDIF             
                  ENDIF
               ENDIF !end chiralt
               !now we should know if the centre is chiral
               IF (.NOT.CHIRALT) THEN
                  POSSIBLE_CC(J1) = .FALSE.
               ELSE
                  NCHIRAL = NCHIRAL + 1
                  DUMMY_CHIR_INFO(NCHIRAL,1) = J1
                  DO J3=1,4
                     DUMMY_CHIR_INFO(NCHIRAL,FINALPRIO(J3)+1)=DUMMY_ATS(J3)
                  ENDDO
               ENDIF
            ENDIF !possible chiral centre
         ENDDO
         ALLOCATE(CHIR_INFO(NCHIRAL,5))
         DO J1=1,NCHIRAL
            CHIR_INFO(J1,1) = DUMMY_CHIR_INFO(J1,1)
            CHIR_INFO(J1,2) = DUMMY_CHIR_INFO(J1,2)
            CHIR_INFO(J1,3) = DUMMY_CHIR_INFO(J1,3)
            CHIR_INFO(J1,4) = DUMMY_CHIR_INFO(J1,4)
            CHIR_INFO(J1,5) = DUMMY_CHIR_INFO(J1,5)
         ENDDO
         WRITE(*,*) " find_chiral_centres> Completed chiral centre detection, number of chiral centres: ", NCHIRAL
         RETURN    
       END SUBROUTINE FIND_CHIRAL_CENTRES

      SUBROUTINE SORT_BONDEDATS()
         IMPLICIT NONE
         INTEGER :: J1, J2, ENDL, MAXPOS, NDUMMY
         INTEGER :: DUMMY_ATS(6), DUMMY_ELS(6)
         INTEGER :: NEW_ATS(6), NEW_ELS(6)  
         
         ALLOCATE(NNEIGHBOURS(NATOMS))
         NNEIGHBOURS(1:NATOMS)=0
         DO J1=1,NATOMS
            NEW_ATS(1:6) = -1
            NEW_ELS(1:6) = -1
            DUMMY_ATS(1:6) = BONDEDATS(J1,1:6)
            DUMMY_ELS(1:6) = BONDEDELS(J1,1:6)
            J2=6
            !find number of entries that are defined (all other entries are -1)
            DO WHILE(J2.GT.0)
               IF (DUMMY_ATS(J2).EQ.-1) THEN
                  J2=J2-1
               ELSE
                  NNEIGHBOURS(J1) = J2
                  EXIT
               ENDIF
            ENDDO
            ENDL = J2
            NDUMMY=1
            !as we have a very short list to sort, we just do it with a basic sorting
            DO WHILE(ENDL.GT.0)
               MAXPOS = MAXLOC(DUMMY_ELS,DIM=1)
               NEW_ATS(NDUMMY) = DUMMY_ATS(MAXPOS)
               NEW_ELS(NDUMMY) = DUMMY_ELS(MAXPOS)
               DUMMY_ELS(MAXPOS) = -1
               ENDL=ENDL-1
               NDUMMY=NDUMMY+1          
            ENDDO
            BONDEDATS(J1,1:6) = NEW_ATS(1:6)
            BONDEDELS(J1,1:6) = NEW_ELS(1:6)
         ENDDO
         RETURN
       END SUBROUTINE SORT_BONDEDATS

      ! get sp2 carbons
      SUBROUTINE GET_SP2_CARBON()
         IMPLICIT NONE
         INTEGER :: LIST_DUMMY(NATOMS)
         INTEGER :: J1
         NCSP2 = 0
         DO J1=1,NATOMS
            IF ((ELEMENT(J1).EQ.6).AND.(NBONDED(J1).EQ.3)) THEN
               NCSP2=NCSP2+1
               LIST_DUMMY(NCSP2) = J1
            ENDIF
         ENDDO
         ALLOCATE(CSP2(NCSP2))
         CSP2(1:NCSP2)=LIST_DUMMY(1:NCSP2)
      END SUBROUTINE GET_SP2_CARBON

      ! find all double bonds (we find a set of double bonds that enables us to add ghost atoms for the priority calculation)
      SUBROUTINE FIND_DOUBLE_BONDS()
         IMPLICIT NONE
         LOGICAL :: PROCESSED(NCSP2)
         INTEGER :: I, J, BONDED2AT(3)
         INTEGER :: ATOM1, ATOM2

         PROCESSED(1:NCSP2) = .FALSE.
         NDOUBLE = 0
         DOUBLEB(1:NATOMS,1:2) = -1
         ! we iterate over the list of C sp2 atoms and find their binding partner
         DO I=1,NCSP2
            ! for C=C bonds we have two atoms in the list and have to avoid adding the bond twice
            IF (.NOT.PROCESSED(I)) THEN
               ! get a list of atoms bonded to this atom (we know there are three atoms)
               ATOM1 = CSP2(I)
               BONDED2AT(1:3) = BONDEDATS(Atom1,1:3)
               ! first we need to make sure we account for all C=C bonds (again making sure we do not double count!)
               DO J=I,NCSP2
                  IF (.NOT.PROCESSED(J)) THEN
                     ATOM2=CSP2(J)
                     IF ((BONDED2AT(1).EQ.ATOM2).OR.(BONDED2AT(2).EQ.ATOM2).OR.(BONDED2AT(3).EQ.ATOM2)) THEN
                         NDOUBLE=NDOUBLE+1
                         DOUBLEB(NDOUBLE,1) = ATOM1
                         DOUBLEB(NDOUBLE,2) = ATOM2
                         PROCESSED(I) = .TRUE.            
                         PROCESSED(J) = .TRUE.
                     ENDIF
                  END IF
               END DO
               ! at this stage we have either accounted for the dopuble bond or it is not a C=C bond,
               ! so we check if it is a C=O bond
               IF (.NOT.PROCESSED(I)) THEN
                  DO J=1,3
                     ATOM2 = BONDED2AT(J)
                     IF ((ELEMENT(ATOM2).EQ.8).AND.(NBONDED(ATOM2).EQ.1)) THEN
                        NDOUBLE=NDOUBLE+1
                        DOUBLEB(NDOUBLE,1) = ATOM1
                        DOUBLEB(NDOUBLE,2) = ATOM2
                        PROCESSED(I) = .TRUE.            
                        EXIT   !if we find a O sp2 - add the bond and leave the loop
                     ENDIF             
                  ENDDO
               END IF
               ! if is no C=C or C=O bond it must be C=N
               IF (.NOT.PROCESSED(I)) THEN 
                  DO J=1,3
                     ATOM2 = BONDED2AT(J)
                     IF (ELEMENT(ATOM2).EQ.7) THEN
                        NDOUBLE=NDOUBLE+1
                        DOUBLEB(NDOUBLE,1) = ATOM1
                        DOUBLEB(NDOUBLE,2) = ATOM2
                        PROCESSED(I) = .TRUE.            
                        EXIT   
                     ENDIF             
                  ENDDO
               ENDIF
            END IF
         END DO

      END SUBROUTINE FIND_DOUBLE_BONDS

      !add ghost atoms
      SUBROUTINE CREATE_GHOST_ATOMS()
         IMPLICIT NONE
         INTEGER :: J1, J2, ATOM1, ATOM2, ELEMENT1, ELEMENT2

         !first we find all C sp2 centres (all double bonds in amino acids and nucleotides contain C)
         CALL GET_SP2_CARBON()
         !we now find the double bonds such that we can add ghost atoms
         CALL FIND_DOUBLE_BONDS()
         !we now have a  list of double bonds
         !let's add the ghost atoms accordingly (element like bonded atom, atom id set to 0)
         DO J1=1,NDOUBLE
            ATOM1 = DOUBLEB(NDOUBLE,1)
            ATOM2 = DOUBLEB(NDOUBLE,2)
            ELEMENT1 = ELEMENT(ATOM1)
            ELEMENT2 = ELEMENT(ATOM2)
            DO J2=1,6
              IF (BONDEDATS(ATOM1,J2).EQ.-1) THEN
                 BONDEDATS(ATOM1,J2) = 0
                 BONDEDELS(ATOM1,J2) = ELEMENT2
                 EXIT
              ENDIF
            ENDDO
            DO J2=1,6
              IF (BONDEDATS(ATOM2,J2).EQ.-1) THEN
                 BONDEDATS(ATOM2,J2) = 0
                 BONDEDELS(ATOM2,J2) = ELEMENT1
                 EXIT
              ENDIF
            ENDDO   
         ENDDO
         DEALLOCATE(CSP2)
         RETURN
      END SUBROUTINE CREATE_GHOST_ATOMS

      !get the number of bonds for each atom and initialise the arrays for bonded atoms and their elements
      SUBROUTINE NBONDED_ATOMS()
         IMPLICIT NONE
         INTEGER :: J1, ATOMID1, ATOMID2
         
         ALLOCATE(DOUBLEB(NBOND,2))
         ALLOCATE(NBONDED(NATOMS))
         ALLOCATE(BONDEDATS(NATOMS,6))
         ALLOCATE(BONDEDELS(NATOMS,6))   
         ALLOCATE(POSSIBLE_CC(NATOMS))
         POSSIBLE_CC(1:NATOMS) = .FALSE.
         NBONDED(1:NATOMS) = 0
         BONDEDATS(1:NATOMS,1:6) = -1
         BONDEDELS(1:NATOMS,1:6) = -1    
         DO J1=1,NBOND
            ATOMID1 = BONDS(J1,1)
            ATOMID2 = BONDS(J1,2)
            NBONDED(ATOMID1) = NBONDED(ATOMID1) + 1 
            NBONDED(ATOMID2) = NBONDED(ATOMID2) + 1        
            BONDEDATS(ATOMID1,NBONDED(ATOMID1)) = ATOMID2
            BONDEDELS(ATOMID1,NBONDED(ATOMID1)) = ELEMENT(ATOMID2)        
            BONDEDATS(ATOMID2,NBONDED(ATOMID2)) = ATOMID1  
            BONDEDELS(ATOMID2,NBONDED(ATOMID2)) = ELEMENT(ATOMID1)     
         ENDDO   
         DO J1=1,NATOMS
            IF (NBONDED(J1).EQ.4) THEN
               POSSIBLE_CC(J1) = .TRUE.
            END IF
         ENDDO
         RETURN
      END SUBROUTINE NBONDED_ATOMS 

      ! check whether we have more than one hydrogens on a carbon - if so it can't be chiral
      SUBROUTINE DISCOUNT_H()
         IMPLICIT NONE
         INTEGER :: J1, J2, HCOUNT
         
         DO J1=1,NATOMS    
            IF (POSSIBLE_CC(J1)) THEN
               HCOUNT=0
               DO J2=1,4
                  IF (BONDEDELS(J1,J2).EQ.1) HCOUNT=HCOUNT+1
               ENDDO
               IF (HCOUNT.GT.1) POSSIBLE_CC(J1)=.FALSE.
            ENDIF
         ENDDO
         RETURN
       END SUBROUTINE DISCOUNT_H       

      !test whether the four neighbours are different elements
      SUBROUTINE TEST_NEIGHBOURS(ATOMID,CHIRALT) 
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ATOMID
         LOGICAL, INTENT(OUT) :: CHIRALT
         INTEGER :: N1, N2, N3, N4
     
         CHIRALT = .FALSE.    
         N1 = BONDEDELS(ATOMID,1) 
         N2 = BONDEDELS(ATOMID,2)    
         N3 = BONDEDELS(ATOMID,3)    
         N4 = BONDEDELS(ATOMID,4)
         !the lists were sorted before, so we can use a simple comparison
         IF ((N1.NE.N2).AND.(N2.NE.N3).AND.(N3.NE.N4)) CHIRALT = .TRUE.
      END SUBROUTINE TEST_NEIGHBOURS

      !check whether atoms are bonded
      SUBROUTINE CHECK_BOND(ATOM1,ATOM2,BONDEDT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ATOM1, ATOM2
         LOGICAL, INTENT(OUT) :: BONDEDT
         INTEGER :: J1
         
         BONDEDT = .FALSE.
         DO J1=1,NBOND
            IF (((BONDS(J1,1).EQ.ATOM1).AND.(BONDS(J1,2).EQ.ATOM2)).OR.  &
               ((BONDS(J1,1).EQ.ATOM2).AND.(BONDS(J1,2).EQ.ATOM1))) THEN
               BONDEDT = .TRUE.
               EXIT
            ENDIF
         ENDDO
      END SUBROUTINE CHECK_BOND 

END MODULE CHIRALITY
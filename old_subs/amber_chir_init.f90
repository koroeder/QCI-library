MODULE AMBER_CHIR_INIT
      
  IMPLICIT NONE
  INTEGER, SAVE :: NBOND                                    !number of bonds
  INTEGER, SAVE :: NATS                                     !number of atoms
  INTEGER, SAVE :: NCISTRANS                                !number of peptide bonds
  INTEGER, SAVE :: NCHIRAL                                  !number of chiral centres
  INTEGER, SAVE :: NRES                                     !number of residues
  INTEGER, SAVE :: NPRO                                     !number of PRO/HYP
  INTEGER, ALLOCATABLE, SAVE :: NBONDED(:)                  !number of bonds for each atom 
  INTEGER, ALLOCATABLE, SAVE :: NNEIGHBOURS(:)              !number of bonded atoms 
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: BONDS       !array for bonds
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: CT_INFO     !array for cis-trans dih
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: CHIR_INFO   !arreay for chiral centres
  CHARACTER(4), DIMENSION(:), ALLOCATABLE, SAVE :: ATNAMES  !atom names
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ELEMENT       !atomic number
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: BONDEDATS   !atoms bonded by atomid (up to 6)
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: BONDEDELS   !atoms bonded by element (up to 6)  
  INTEGER , ALLOCATABLE, SAVE ::  RES_START(:), RES_END(:)  !start/end indices
  CHARACTER(4), DIMENSION(:), ALLOCATABLE, SAVE :: RESNAMES !residue names
  CHARACTER(3), DIMENSION(:), ALLOCATABLE, SAVE :: TYPERES  !RNA,DNA,PRO(tein) or OTH(er) for each atom
  LOGICAL, ALLOCATABLE, SAVE :: PRO_RES(:)                  !list of prolines
  LOGICAL, ALLOCATABLE, SAVE :: IS_AAT(:)                   !are residues amino acids?
  LOGICAL, ALLOCATABLE, SAVE :: POSSIBLE_CC(:)              !potential chiral centres
  INTEGER, ALLOCATABLE, SAVE :: POSSIBLE_IDS(:,:)           !ids of bonded atoms
  CONTAINS
!add global deallocate
  SUBROUTINE SETUP_CT_CHIRAL()
     CALL CLEAR_VAR()           !make sure everything is deallocated
     CALL TOP_READER()          !get information from topology
     CALL PARSE_PROLINE()       !create some data for prolines
     CALL FIND_CISTRANS()       !get the peptide bonds
     CALL FIND_CHIRAL_CENTRES() !get chiral centres
     RETURN
  END SUBROUTINE SETUP_CT_CHIRAL

  SUBROUTINE CLEAR_VAR()
  
    IF (ALLOCATED(NBONDED)) DEALLOCATE(NBONDED)
    IF (ALLOCATED(NNEIGHBOURS)) DEALLOCATE(NNEIGHBOURS)
    IF (ALLOCATED(BONDS)) DEALLOCATE(BONDS)
    IF (ALLOCATED(CT_INFO)) DEALLOCATE(CT_INFO)
    IF (ALLOCATED(CHIR_INFO)) DEALLOCATE(CHIR_INFO) 
    IF (ALLOCATED(ATNAMES)) DEALLOCATE(ATNAMES)  
    IF (ALLOCATED(ELEMENT)) DEALLOCATE(ELEMENT) 
    IF (ALLOCATED(BONDEDATS)) DEALLOCATE(BONDEDATS)  
    IF (ALLOCATED(BONDEDELS)) DEALLOCATE(BONDEDELS)
    IF (ALLOCATED(RES_START)) DEALLOCATE(RES_START)
    IF (ALLOCATED(RES_END)) DEALLOCATE(RES_END)  
    IF (ALLOCATED(RESNAMES)) DEALLOCATE(RESNAMES) 
    IF (ALLOCATED(PRO_RES)) DEALLOCATE(PRO_RES)
    IF (ALLOCATED(IS_AAT)) DEALLOCATE(IS_AAT)
    IF (ALLOCATED(POSSIBLE_CC)) DEALLOCATE(POSSIBLE_CC) 
    IF (ALLOCATED(POSSIBLE_IDS)) DEALLOCATE(POSSIBLE_IDS)     
    RETURN
  END SUBROUTINE CLEAR_VAR

  SUBROUTINE FIND_CHIRAL_CENTRES()
    USE COMMONS, ONLY: DEBUG
    USE KEY, ONLY: MYUNIT
    IMPLICIT NONE
    INTEGER :: J1, J2, ATOMID, J3, ATOMID2, VALENCY, GHOSTEL, J4
    INTEGER :: DUMMY_ELS(4), DUMMY_ELS_2(4), DUMMY_ATS(4), DUMMY_ATS_2(4)
    INTEGER :: PRIO4(4), PRIO3(3), FINALPRIO(4)
    INTEGER :: LIST1(4), LIST2(4), LIST3(4), LIST4(4)
    LOGICAL :: CHIRALT, EQ12, EQ23, EQ34, GHOSTS(4), TRIPLE, QUADRUPLE
    LOGICAL :: GT12T,ID12T,GT23T,ID23T,GT34T,ID34T,CHECKUNIQUE(4),IDFOUND
    INTEGER :: DUMMY_CHIR_INFO(NATS,5)     !array for chiral centres
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
    DO J1=1,NATS
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
             ELSE IF ((ATNAMES(J1).EQ."C3'").AND.(TYPERES(J1).EQ."RNA")) THEN
                C3AT1 = DUMMY_ATS(2)
                IF (ATNAMES(C3AT1).EQ."C2'") THEN
                   FINALPRIO(3) = FINALPRIO(3) + 1
                   FINALPRIO(4) = FINALPRIO(4) + 1
                   CHIRALT=.TRUE.
                ELSE IF (ATNAMES(C3AT1).EQ."C4'") THEN 
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
    WRITE(*,'(A,I6,A)') " find_chiral> ", NCHIRAL, " chiral atoms found"
    RETURN    
  END SUBROUTINE FIND_CHIRAL_CENTRES

  !comparing four sorted lists of atoms
  SUBROUTINE COMPARE_4LISTS(LIST1,LIST2,LIST3,LIST4,LENGTH,PRIORITY)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: LENGTH                           !length of lists
     INTEGER, INTENT(IN) :: LIST1(LENGTH), LIST2(LENGTH)     !list of elements
     INTEGER, INTENT(IN) :: LIST3(LENGTH), LIST4(LENGTH)     !list of elements     
     INTEGER, INTENT(OUT) :: PRIORITY(4)                     !priorities of list (1 to 3)
     INTEGER :: P123(3), COMP4(3), J1, J2
     LOGICAL :: ID14T, GT14T, ID24T, GT24T, ID34T, GT34T

     CALL COMPARE_3LISTS(LIST1,LIST2,LIST3,LENGTH,P123)
     CALL COMPARE_2LISTS(LIST1,LIST4,LENGTH,ID14T,GT14T)
     CALL COMPARE_2LISTS(LIST2,LIST4,LENGTH,ID24T,GT24T)     
     CALL COMPARE_2LISTS(LIST3,LIST4,LENGTH,ID34T,GT34T)
     IF (ID14T) THEN
        COMP4(1) = 0
     ELSE IF (GT14T) THEN
        COMP4(1) = -1 !lower priority than one
     ELSE
        COMP4(1) = 1
     ENDIF
     IF (ID24T) THEN
        COMP4(2) = 0
     ELSE IF (GT24T) THEN
        COMP4(2) = -1 !lower priority than two
     ELSE
        COMP4(2) = 1
     ENDIF 
      IF (ID34T) THEN
        COMP4(3) = 0
     ELSE IF (GT34T) THEN
        COMP4(3) = -1 !lower priority than three
     ELSE
        COMP4(3) = 1
     ENDIF     
     !we know have a list of priorities for 1 to 3 and the comparison for 4
     !if comp4 has any zeros, we simply assign it the same priority
     DO J1=1,3
        IF (COMP4(J1).EQ.0) THEN
           PRIORITY(1:3) = P123(1:3)        
           PRIORITY(4) = P123(J1)
           RETURN
        ENDIF
     ENDDO
     ! the priorities need to be adjusted now to fit in the 4th list
     DO J2=1,3
        IF (COMP4(J2).EQ.1) THEN
           PRIORITY(J2) = P123(J2) + 1
        ELSE
           PRIORITY(J2) = P123(J2)
        ENDIF
     ENDDO
     !the possible values of the sum are -3,-1,1,3, substracting that from 5 and 
     !dividing by 2 gives the correct priority
     PRIORITY(4) = (5-SUM(COMP4(1:3)))/2
     RETURN
  END SUBROUTINE COMPARE_4LISTS

  !comparing three sorted lists of atoms
  SUBROUTINE COMPARE_3LISTS(LIST1,LIST2,LIST3,LENGTH,PRIORITY)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: LENGTH                           !length of lists
     INTEGER, INTENT(IN) :: LIST1(LENGTH), LIST2(LENGTH), LIST3(LENGTH)     !list of elements
     INTEGER, INTENT(OUT) :: PRIORITY(3)                     !priorities of list (1 to 3)
     LOGICAL :: ID12T, ID13T, ID23T, GT12T, GT13T, GT23T
     
     !return for priorities. 0 means can't assign (i.e. identical lists),
     !otherwise we assign 1 to 3 for highest to lowest priority
     PRIORITY(1:3) = 0
     !we simply use the comparison between two lists for all three combinations
     CALL COMPARE_2LISTS(LIST1,LIST2,LENGTH,ID12T,GT12T)
     CALL COMPARE_2LISTS(LIST1,LIST3,LENGTH,ID13T,GT13T)
     CALL COMPARE_2LISTS(LIST2,LIST3,LENGTH,ID23T,GT23T)
     !if list 1 and 2 are identical
     IF (ID12T) THEN
        ! if list 1 and 3 are also identical, we can't assign priorities here
        IF (ID13T) THEN
           PRIORITY(1:3) = (/1,1,1/)
           RETURN
        ELSE
           ! if 1 has a higher priority we know 3 has the lowest priority
           IF (GT13T) THEN
              PRIORITY(1:3) = (/1,1,2/)
              RETURN
           ! otherwise it has the highest
           ELSE
              PRIORITY(1:3) = (/2,2,1/) 
              RETURN
           ENDIF
        ENDIF
     !if list 1 and 3 are identical, we don't need to check for all three identical
     !this would have been caught before
     ELSE IF (ID13T) THEN 
        ! if 1 has a higher priority we know 2 has the lowest priority           
        IF (GT12T) THEN
           PRIORITY(1:3) = (/1,2,1/)
           RETURN
        ELSE
           PRIORITY(1:3) = (/2,1,2/)
           RETURN
        ENDIF
     !if list 1 and 3 are identical, we don't need to check for all three identical
     !this would have been caught before        
     ELSE IF (ID23T) THEN
        ! if 1 has a higher priority we know 1 has the highest priority           
        IF (GT12T) THEN
           PRIORITY(1:3) = (/1,2,2/)
           RETURN
        ELSE
           PRIORITY(1:3) = (/2,1,1/)
           RETURN
        ENDIF
     !none are identical
     ELSE
        !there are six possible orderings
        !1>2>3
        IF (GT12T.AND.GT23T) THEN
           PRIORITY(:) = (/1,2,3/)
        ELSE IF (GT12T.AND.(.NOT.GT23T)) THEN
           !1>3>2
           IF (GT13T) THEN
              PRIORITY(:) = (/1,3,2/)
           !3>1>2
           ELSE
              PRIORITY(:) = (/3,1,2/)
           ENDIF
        !2>1>3
        ELSE IF (.NOT.GT12T.AND.GT13T) THEN
           PRIORITY(:) = (/2,1,3/)
        ELSE IF (.NOT.GT12T.AND.(.NOT.GT13T)) THEN
           !2>3>1
           IF (GT23T) THEN
              PRIORITY(:) = (/2,3,1/)
           !3>2>1
           ELSE
              PRIORITY(:) = (/3,2,1/)
           ENDIF
        ENDIF             
     ENDIF
     RETURN
  END SUBROUTINE COMPARE_3LISTS

  !comparing two sorted lists of atoms
  SUBROUTINE COMPARE_2LISTS(LIST1,LIST2,LENGTH,IDENTICALT,HIGH1T)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: LENGTH                           !length of lists
     INTEGER, INTENT(IN) :: LIST1(LENGTH), LIST2(LENGTH)     !list of priorities
     LOGICAL, INTENT(OUT) :: IDENTICALT                      !are the lists identical
     LOGICAL, INTENT(OUT) :: HIGH1T                          !higher priority for list 1
     INTEGER :: J1
     
     IDENTICALT = .FALSE.
     DO J1=1,LENGTH
        IF (LIST1(J1).GT.LIST2(J1)) THEN
           HIGH1T = .TRUE.
           RETURN
        ELSE IF (LIST1(J1).LT.LIST2(J1)) THEN
           HIGH1T = .FALSE.
           RETURN
        ENDIF
     ENDDO
     IDENTICALT = .TRUE.
     HIGH1T = .FALSE. !set this to stop compiler complaining
     RETURN
  END SUBROUTINE COMPARE_2LISTS

  SUBROUTINE SORT_BONDEDATS()
    IMPLICIT NONE
    INTEGER :: J1, J2, ENDL, MAXPOS, NDUMMY
    INTEGER :: DUMMY_ATS(6), DUMMY_ELS(6)
    INTEGER :: NEW_ATS(6), NEW_ELS(6)  
    
    ALLOCATE(NNEIGHBOURS(NATS))
    NNEIGHBOURS(1:NATS)=0
    DO J1=1,NATS
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

  !add ghost atoms
  SUBROUTINE CREATE_GHOST_ATOMS()
    IMPLICIT NONE
    INTEGER :: J1, J2, J3, ATOM1, ATOM2, J4, ELEMENT1, ELEMENT2, J5
    INTEGER :: DOUBLEB(NBOND,2), NDOUBLE, NCSP2, BONDED2AT(3)
    INTEGER, ALLOCATABLE :: LIST_DUMMY(:), CSP2(:)
    LOGICAL, ALLOCATABLE :: GHOST_ADDED(:)
    !first we need to check the valency vs the number of bonds
    !in amino acids and nucleic acids all double bonds involve at least one C
    !so we first search for those
    NCSP2 = 0
    ALLOCATE(LIST_DUMMY(NATS))
    DO J1=1,NATS
       IF ((ELEMENT(J1).EQ.6).AND.(NBONDED(J1).EQ.3)) THEN
          NCSP2=NCSP2+1
          LIST_DUMMY(NCSP2) = J1
       ENDIF
    ENDDO
    ALLOCATE(CSP2(NCSP2),GHOST_ADDED(NCSP2))
    CSP2(1:NCSP2)=LIST_DUMMY(1:NCSP2)
    GHOST_ADDED(1:NCSP2)=.FALSE.
    DEALLOCATE(LIST_DUMMY)
    !we now go through the ids for every atom in the list and see if another sp2 matches
    !if it does we add them to the list for ghost atoms to be added
    NDOUBLE=0
    DO J2=1,NCSP2
       IF (.NOT.GHOST_ADDED(J2)) THEN
          ATOM1=CSP2(J2)
          BONDED2AT(1:3)=BONDEDATS(ATOM1,1:3)
          DO J3=J2,NCSP2
             IF (.NOT.GHOST_ADDED(J3)) THEN
                ATOM2=CSP2(J3)
                IF ((BONDED2AT(1).EQ.ATOM2).OR.(BONDED2AT(2).EQ.ATOM2).OR. &
                    (BONDED2AT(3).EQ.ATOM2)) THEN
                    NDOUBLE=NDOUBLE+1
                    DOUBLEB(NDOUBLE,1) = ATOM1
                    DOUBLEB(NDOUBLE,2) = ATOM2
                    GHOST_ADDED(J2) = .TRUE.            
                    GHOST_ADDED(J3) = .TRUE.
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !we now have added all double bonds between two carbons
    !next we add C=O and C=N
    !as we accoutned for C=C, we will check for C=O and otherwise add C=N
    DO J2=1,NCSP2
       IF (.NOT.GHOST_ADDED(J2)) THEN  
          ATOM1=CSP2(J2)
          BONDED2AT(1:3)=BONDEDATS(ATOM1,1:3)
          DO J3=1,3
             ATOM2 = BONDED2AT(J3)
             IF ((ELEMENT(ATOM2).EQ.8).AND.(NBONDED(ATOM2).EQ.1)) THEN
                NDOUBLE=NDOUBLE+1
                DOUBLEB(NDOUBLE,1) = ATOM1
                DOUBLEB(NDOUBLE,2) = ATOM2
                GHOST_ADDED(J2) = .TRUE.            
                EXIT   !if we find a O sp2 - add the bond and leave the loop
             ENDIF             
          ENDDO
          !only option left now is to find a C=N bond, valency is an issue here due to charges
          IF (.NOT.GHOST_ADDED(J2)) THEN 
             DO J3=1,3
                ATOM2 = BONDED2AT(J3)
                IF (ELEMENT(ATOM2).EQ.7) THEN
                   NDOUBLE=NDOUBLE+1
                   DOUBLEB(NDOUBLE,1) = ATOM1
                   DOUBLEB(NDOUBLE,2) = ATOM2
                   GHOST_ADDED(J2) = .TRUE.            
                   EXIT   
                ENDIF             
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    !we now have a  list of double bonds
    !let's add the ghost atoms accordingly (element like bonded atom, atom id set to 0)
    DO J4=1,NDOUBLE
       ATOM1 = DOUBLEB(NDOUBLE,1)
       ATOM2 = DOUBLEB(NDOUBLE,2)
       ELEMENT1 = ELEMENT(ATOM1)
       ELEMENT2 = ELEMENT(ATOM2)
       DO J5=1,6
         IF (BONDEDATS(ATOM1,J5).EQ.-1) THEN
            BONDEDATS(ATOM1,J5) = 0
            BONDEDELS(ATOM1,J5) = ELEMENT2
            EXIT
         ENDIF
       ENDDO
       DO J5=1,6
         IF (BONDEDATS(ATOM2,J5).EQ.-1) THEN
            BONDEDATS(ATOM2,J5) = 0
            BONDEDELS(ATOM2,J5) = ELEMENT1
            EXIT
         ENDIF
       ENDDO   
    ENDDO
    DEALLOCATE(GHOST_ADDED,CSP2)
    RETURN
  END SUBROUTINE CREATE_GHOST_ATOMS

  !sum up the number of bonds for each atom
  SUBROUTINE NBONDED_ATOMS()
    IMPLICIT NONE
    INTEGER :: J1, ATOMID1, ATOMID2
      
    ALLOCATE(NBONDED(NATS))
    ALLOCATE(BONDEDATS(NATS,6))
    ALLOCATE(BONDEDELS(NATS,6))   
    ALLOCATE(POSSIBLE_CC(NATS))
    POSSIBLE_CC(1:NATS) = .FALSE.
    NBONDED(1:NATS) = 0
    BONDEDATS(1:NATS,1:6) = -1
    BONDEDELS(1:NATS,1:6) = -1    
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
    DO J1=1,NATS
       IF (NBONDED(J1).EQ.4) POSSIBLE_CC(J1) = .TRUE.
    ENDDO
    RETURN
  END SUBROUTINE NBONDED_ATOMS 


  SUBROUTINE DISCOUNT_H()
    IMPLICIT NONE
    INTEGER :: J1, J2, HCOUNT
    
    DO J1=1,NATS    
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
    RETURN
  END SUBROUTINE TEST_NEIGHBOURS

  SUBROUTINE FIND_CISTRANS()
    USE KEY, ONLY: MYUNIT
    IMPLICIT NONE
    INTEGER J1, NPOS, OPOS, HPOS, CPOS
    LOGICAL :: CNT, OCT, NHT
    INTEGER :: DUMMY_CT_INFO(NRES,4)     !array for cis-trans dihedrals
    
    NCISTRANS=0
    !iterate over all residues
    DO J1=1,NRES-1
       !get the atom positions
       CALL GET_ATOMID("O",J1,OPOS)      
       CALL GET_ATOMID("C",J1,CPOS)   
       CALL GET_ATOMID("N",J1+1,NPOS)
       IF (PRO_RES(J1+1)) THEN !don't add a constraint for proline
          CYCLE
!          CALL GET_ATOMID("CD",J1+1,HPOS) !if we have proline use CD instead of H
       ELSE
          CALL GET_ATOMID("H",J1+1,HPOS)             
       ENDIF
       !if all of them exists, check bonding    
       IF ((OPOS.GT.0).AND.(NPOS.GT.0).AND.(HPOS.GT.0).AND.(CPOS.GT.0)) THEN
          CALL CHECK_BOND(CPOS,NPOS,CNT)
          CALL CHECK_BOND(NPOS,HPOS,NHT)
          CALL CHECK_BOND(OPOS,CPOS,OCT)
          IF (CNT.AND.NHT.AND.OCT) THEN
             IF (.NOT.(IS_AAT(J1).AND.IS_AAT(J1+1))) THEN
                WRITE(*,'(A,I6,A,I6)') " find_cistrans> WARNING! peptide bond found for non-standard residue(s): ", J1, "and", J1+1     
             ENDIF
             NCISTRANS = NCISTRANS+1
             DUMMY_CT_INFO(NCISTRANS,1) = OPOS
             DUMMY_CT_INFO(NCISTRANS,2) = CPOS   
             DUMMY_CT_INFO(NCISTRANS,3) = NPOS
             DUMMY_CT_INFO(NCISTRANS,4) = HPOS
          ENDIF
       ENDIF                                   
    ENDDO
    ALLOCATE(CT_INFO(NCISTRANS,4))
    DO J1=1,NCISTRANS
       CT_INFO(J1,1) = DUMMY_CT_INFO(J1,1)
       CT_INFO(J1,2) = DUMMY_CT_INFO(J1,2)       
       CT_INFO(J1,3) = DUMMY_CT_INFO(J1,3)       
       CT_INFO(J1,4) = DUMMY_CT_INFO(J1,4)       
    ENDDO
    WRITE(*,'(A,I6,A)') " find_cistrans> ", NCISTRANS, " peptide bonds detected"
    RETURN
  END SUBROUTINE

  !get atom id from name for given residue
  !returns 0 if atom does not exist in residue
  SUBROUTINE GET_ATOMID(ATNAME,RESID,ATOMID)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RESID
    CHARACTER(*), INTENT(IN) :: ATNAME
    INTEGER, INTENT(OUT) :: ATOMID
    INTEGER :: FIRST, LAST, J1
    
    ATOMID = 0
    FIRST=RES_START(RESID)
    LAST=RES_END(RESID)
    DO J1=FIRST,LAST
       IF (ADJUSTL(TRIM(ATNAMES(J1))).EQ.ADJUSTL(TRIM(ATNAME))) THEN
          ATOMID = J1
          EXIT
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE GET_ATOMID
  
  !check whether atoms are bonded
  SUBROUTINE CHECK_BOND(ATOM1,ATOM2,BONDEDT)
    USE COMMONS, ONLY: DEBUG
    USE KEY, ONLY: MYUNIT
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
    IF (DEBUG) THEN
       IF (BONDEDT) THEN
          WRITE(*,'(A,I6,A,I6,A)') " check_bond> Atom ", ATOM1 , &
             " and ", ATOM2 , " are bonded"
       ELSE
          WRITE(*,'(A,I6,A,I6,A)') " check_bond> Atom ", ATOM1 , &
             " and ", ATOM2 , " are not bonded"       
       ENDIF
    ENDIF 
    RETURN
  END SUBROUTINE CHECK_BOND  
   
  ! check if residue is terminal
  SUBROUTINE CHECK_TER(RESID,NTERT,CTERT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RESID
    LOGICAL, INTENT(OUT) :: NTERT,CTERT
    INTEGER :: FIRST, LAST, J1
    
    IF ((RESID.LT.1).OR.(RESID.GT.NRES)) THEN
       STOP (' check_termini> residue id out of range')
    ENDIF
    NTERT = .FALSE.
    CTERT = .FALSE.
    FIRST = RES_START(RESID)
    LAST = RES_END(RESID)
    DO J1=FIRST,LAST
       IF (ATNAMES(J1).EQ."OXT") THEN
          CTERT=.TRUE.
       ELSE IF (ATNAMES(J1).EQ."H1") THEN
          NTERT=.TRUE.
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE CHECK_TER
  
  ! check if residue is an amino acid - brute force ...
  SUBROUTINE CHECK_AA(RESID,AAT)
    USE COMMONS, ONLY: DEBUG
    USE KEY, ONLY: MYUNIT
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RESID
    LOGICAL, INTENT(OUT) :: AAT
    CHARACTER(4) :: DNAME
    LOGICAL :: TERTEST
    
    DNAME = ADJUSTL(TRIM(RESNAMES(RESID)))
    AAT = .FALSE.
    TERTEST = .FALSE.
60  CONTINUE
    IF ((DNAME.EQ."ALA").OR.(DNAME.EQ."ARG").OR.(DNAME.EQ."ASH").OR.  &
        (DNAME.EQ."ASN").OR.(DNAME.EQ."ASP").OR.(DNAME.EQ."CYM").OR.  &
        (DNAME.EQ."CYS").OR.(DNAME.EQ."CYX").OR.(DNAME.EQ."GLH").OR.  &
        (DNAME.EQ."GLN").OR.(DNAME.EQ."GLU").OR.(DNAME.EQ."GLY").OR.  &
        (DNAME.EQ."HID").OR.(DNAME.EQ."HIE").OR.(DNAME.EQ."HIP").OR.  &
        (DNAME.EQ."HYP").OR.(DNAME.EQ."ILE").OR.(DNAME.EQ."LEU").OR.  &
        (DNAME.EQ."LYN").OR.(DNAME.EQ."LYS").OR.(DNAME.EQ."MET").OR.  &
        (DNAME.EQ."PHE").OR.(DNAME.EQ."PRO").OR.(DNAME.EQ."SER").OR.  &
        (DNAME.EQ."THR").OR.(DNAME.EQ."TRP").OR.(DNAME.EQ."TYR").OR.  &
        (DNAME.EQ."VAL")) THEN
       AAT = .TRUE.
       IF(DEBUG) WRITE(*,'(A,I6,A)') " check_aa> Residue ",RESID," is an amino acid"
    ENDIF
    !make sure it is not a terminal residue (we technically should never encounter one)
    IF (((.NOT.AAT).AND.(.NOT.TERTEST)).AND.((DNAME(1:1).EQ."N").OR.(DNAME(1:1).EQ."C"))) THEN
       TERTEST=.TRUE.
       DNAME = DNAME(2:4)
       GOTO 60
    ENDIF
    RETURN
  END SUBROUTINE CHECK_AA
  
  !get residue type for every atom
  SUBROUTINE FILL_RESTYPE()
    INTEGER :: STARTAT,ENDAT
    INTEGER :: J1, J2
    CHARACTER(4) :: DNAME     
    
    TYPERES(:) = "OTH"
    DO J1=1,NRES
       STARTAT = RES_START(J1)
       ENDAT = RES_END(J1)
       DNAME = ADJUSTL(TRIM(RESNAMES(J1))) 
       IF ((DNAME.EQ."A").OR.(DNAME.EQ."A3").OR.(DNAME.EQ."A5").OR.(DNAME.EQ."AN").OR. &
           (DNAME.EQ."C").OR.(DNAME.EQ."C3").OR.(DNAME.EQ."C5").OR.(DNAME.EQ."CN").OR. &   
           (DNAME.EQ."G").OR.(DNAME.EQ."G3").OR.(DNAME.EQ."G5").OR.(DNAME.EQ."GN").OR. &   
           (DNAME.EQ."U").OR.(DNAME.EQ."U3").OR.(DNAME.EQ."U5").OR.(DNAME.EQ."UN")) THEN
          DO J2=STARTAT,ENDAT
             TYPERES(J2) = "RNA"
          ENDDO
       ELSE IF ((DNAME.EQ."DA").OR.(DNAME.EQ."DA3").OR.(DNAME.EQ."DA5").OR.(DNAME.EQ."DAN").OR. &
                (DNAME.EQ."DC").OR.(DNAME.EQ."DC3").OR.(DNAME.EQ."DC5").OR.(DNAME.EQ."DCN").OR. &   
                (DNAME.EQ."DG").OR.(DNAME.EQ."DG3").OR.(DNAME.EQ."DG5").OR.(DNAME.EQ."DGN").OR. &   
                (DNAME.EQ."DT").OR.(DNAME.EQ."DT3").OR.(DNAME.EQ."DU5").OR.(DNAME.EQ."DUN")) THEN
          DO J2=STARTAT,ENDAT
             TYPERES(J2) = "DNA"
          ENDDO
       ELSE                           
          !as this is called after is_aat is assigned, we simply can use this assignment for amino acids
          IF (IS_AAT(J1)) THEN
             DO J2=STARTAT,ENDAT
                TYPERES(J2) = "PRO"
             ENDDO         
          ENDIF
       ENDIF       
    ENDDO
    RETURN
  END SUBROUTINE FILL_RESTYPE
  
  
  
  !find prolines, as we need to use a different definition for the omega angle
  SUBROUTINE PARSE_PROLINE()
    USE COMMONS, ONLY: DEBUG
    USE KEY, ONLY: MYUNIT
    IMPLICIT NONE
    INTEGER :: J1
    CHARACTER(4) :: STRDUMMY
    
    NPRO=0
    !first iteration to find number of prolines
    DO J1=1,NRES
        STRDUMMY = ADJUSTL(TRIM(RESNAMES(J1)))
        IF ((STRDUMMY.EQ.'PRO').OR.(STRDUMMY.EQ.'CPRO').OR.(STRDUMMY.EQ.'NPRO') &
             .OR.(STRDUMMY.EQ.'HYP').OR.(STRDUMMY.EQ.'CHYP')                    &
             .OR.(STRDUMMY.EQ.'NHYP')) NPRO=NPRO+1
    ENDDO
    ALLOCATE(PRO_RES(NRES))
    PRO_RES(:) = .FALSE.
    !if there are prolines iterate over the listof residue names again to get indices
    IF (NPRO.GT.0) THEN
       DO J1=1,NRES
          STRDUMMY = ADJUSTL(TRIM(RESNAMES(J1)))
          IF ((STRDUMMY.EQ.'PRO').OR.(STRDUMMY.EQ.'CPRO').OR.(STRDUMMY.EQ.'NPRO') &
             .OR.(STRDUMMY.EQ.'HYP').OR.(STRDUMMY.EQ.'CHYP')                      &
             .OR.(STRDUMMY.EQ.'NHYP')) THEN
             PRO_RES(J1)=.TRUE.
             IF (DEBUG) THEN
                WRITE(*,'(A,I8)') ' cistrans_init> PRO/HYP detected, residue: ', J1
             ENDIF
          ENDIF
       ENDDO   
    ENDIF
    RETURN
  END SUBROUTINE PARSE_PROLINE

  SUBROUTINE TOP_READER()
    USE KEY,ONLY : MYUNIT
    USE PORFUNCS 
    IMPLICIT NONE
    CHARACTER(100) ENTRY
    INTEGER :: MYUNIT2, GETUNIT
    INTEGER :: J1, START_IND, END_IND, NBONDH, NBONDA, NENTRIES
    INTEGER :: HENTRIES, J3, J4, J5, NDUMMY, INTDUM, J6
    INTEGER , PARAMETER :: NWORDS=20
    CHARACTER(25) :: ENTRIES(NWORDS)=''
    CHARACTER(4) :: NAMES_CURR(20)
    LOGICAL :: AAT

    MYUNIT2=GETUNIT()
    OPEN(MYUNIT2,FILE='coords.prmtop',STATUS='OLD')
reading:DO
98    READ(MYUNIT2,'(A)',END=99) ENTRY
      CALL READ_LINE(ENTRY,NWORDS,ENTRIES)      !get all words in line
      IF (ENTRIES(2).EQ.'POINTERS') THEN        !get number of bonds
         READ(MYUNIT2,*)                        !ignore format identifier after flag
         READ(MYUNIT2,'(A)',END=99) ENTRY
         CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
         READ(ENTRIES(1),'(I8)') NATS
         READ(ENTRIES(3),'(I8)') NBONDH
         READ(ENTRIES(4),'(I8)') NBONDA
         READ(MYUNIT2,'(A)',END=99) ENTRY
         CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
         READ(ENTRIES(2),'(I8)') NRES       
         NBOND = NBONDH + NBONDA
         WRITE(*,'(A,I8)') ' chirality init> Number of bonds:',NBOND
         ALLOCATE(BONDS(NBOND,2))
         ALLOCATE(ATNAMES(NATS),ELEMENT(NATS),TYPERES(NATS))
         ALLOCATE(RESNAMES(NRES),RES_START(NRES),RES_END(NRES),IS_AAT(NRES))
         IS_AAT(1:NRES) = .FALSE.
      ENDIF
      IF (ENTRIES(2).EQ. 'ATOM_NAME') THEN      
         READ(MYUNIT2,*)                       !ignore format identifier after flag   
         NENTRIES=NATS/20
         IF(NATS.GT.NENTRIES*20) NENTRIES=NENTRIES+1
         NDUMMY=1
         DO J3=1,NENTRIES !go through all lines
            READ(MYUNIT2,'(20 A4)',END=99) NAMES_CURR
            DO J4=1,20
               ATNAMES(NDUMMY) = NAMES_CURR(J4)
               NDUMMY = NDUMMY+1
               IF(NDUMMY.GT.NATS) EXIT
            ENDDO
         ENDDO 
      ENDIF 
      IF (ENTRIES(2).EQ. 'RESIDUE_LABEL') THEN
         READ(MYUNIT2,*)                       !ignore format identifier after flag   
         NENTRIES=NRES/20
         IF(NRES.GT.NENTRIES*20) NENTRIES=NENTRIES+1
         NDUMMY=1
         DO J3=1,NENTRIES !go through all lines
            READ(MYUNIT2,'(20 A4)',END=99) NAMES_CURR
            DO J4=1,20
               RESNAMES(NDUMMY) = NAMES_CURR(J4)
               NDUMMY = NDUMMY+1
               IF(NDUMMY.GT.NRES) EXIT
            ENDDO
         ENDDO 
      ENDIF       
      IF (ENTRIES(2).EQ. 'RESIDUE_POINTER') THEN
         READ(MYUNIT2,*)                       !ignore format identifier after flag   
         NENTRIES=NRES/10   
         IF(NRES.GT.NENTRIES*10) NENTRIES=NENTRIES+1 
         NDUMMY=1
         DO J3=1,NENTRIES !go through all lines  
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES) 
            J4=1
            DO WHILE(J4.LE.10)
               READ(ENTRIES(J4),'(I8)') INTDUM
               RES_START(NDUMMY) = INTDUM
               IF (NDUMMY.GT.1) RES_END(NDUMMY-1) = INTDUM - 1 
               J4=J4+1
               NDUMMY = NDUMMY+1               
               IF(NDUMMY.GT.NRES) EXIT
            ENDDO
            RES_END(NDUMMY-1) = NATS
         ENDDO
      ENDIF
      IF (ENTRIES(2).EQ. 'ATOMIC_NUMBER') THEN  
         READ(MYUNIT2,*)                       !ignore format identifier after flag   
         NENTRIES=NATS/10   
         IF(NATS.GT.NENTRIES*10) NENTRIES=NENTRIES+1 
         NDUMMY=1
         DO J3=1,NENTRIES !go through all lines  
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES) 
            J4=1
            DO WHILE(J4.LE.10)
               READ(ENTRIES(J4),'(I8)') INTDUM
               ELEMENT(NDUMMY) = INTDUM   
               J4=J4+1                                                   
               NDUMMY = NDUMMY+1               
               IF(NDUMMY.GT.NATS) EXIT
            ENDDO
         ENDDO         
      ENDIF                            
      IF (ENTRIES(2).EQ. 'BONDS_INC_HYDROGEN') THEN
         READ(MYUNIT2,*)                             !ignore format identifier after flag
         HENTRIES=(NBONDH*3)/10
         HENTRIES=HENTRIES+((NBONDH*3)-(HENTRIES*10)) !number of lines of entries
         NDUMMY=1
         J5=1
         DO J3=1,HENTRIES                             !go through all lines
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
            J4=1
            DO WHILE(J4.LE.10)
               IF (NDUMMY.LE.NBONDH) THEN
                  IF (J5.EQ.1) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,1) = INTDUM/3+1        !atom1
                     J5=2
                  ELSE IF (J5.EQ.2) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,2) = INTDUM/3+1        !atom2
                     J5=3
                  ELSE
                     J5=1
                     NDUMMY=NDUMMY+1
                  ENDIF
               ELSE
                  GOTO 98
               ENDIF
               J4=J4+1
            ENDDO
         ENDDO
      ENDIF
      IF (ENTRIES(2).EQ. 'BONDS_WITHOUT_HYDROGEN') THEN
         READ(MYUNIT2,*)                             !ignore format identifier after flag
         HENTRIES=(NBONDA*3)/10
         HENTRIES=HENTRIES+((NBONDA*3)-(HENTRIES*10)) !number of lines of entries
         NDUMMY=NBONDH+1
         J5=1
         DO J3=1,HENTRIES                             !go through all lines
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
            J4=1
            DO WHILE(J4.LE.10)
               IF (NDUMMY.LE.(NBONDH+NBONDA)) THEN
                  IF (J5.EQ.1) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,1) = INTDUM/3+1
                     J5=2
                  ELSE IF (J5.EQ.2) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,2) = INTDUM/3+1
                     J5=3
                  ELSE
                     J5=1
                     NDUMMY=NDUMMY+1
                  ENDIF
               ELSE
                  GOTO 98
               ENDIF
               J4=J4+1
            ENDDO
         ENDDO
      ENDIF
    ENDDO reading
99  CLOSE(MYUNIT2)
    DO J6=1,NRES
       CALL CHECK_AA(J6,AAT)
       IS_AAT(J6) = AAT
    ENDDO
    CALL FILL_RESTYPE()
!  DO J6=1,NBOND
!     WRITE(MYUNIT,'(A,I8,A,I8)') 'readtopology> Bond between',BONDS(J6,1),' and',BONDS(J6,2)
!  ENDDO
  END SUBROUTINE TOP_READER


! Routine to read line and split by spaces
! For some parts of the topology using the format flags is better
  SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
      CHARACTER(*), INTENT(IN) :: LINE
      INTEGER, INTENT(IN) :: NWORDS
      CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
      INTEGER:: J1,START_IND,END_IND,J2
      CHARACTER(25) :: WORD
      START_IND=0
      END_IND=0
      J1=1
      J2=0
      DO WHILE(J1.LE.LEN(LINE))
          IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
             START_IND=J1
          ENDIF
          IF (START_IND.GT.0) THEN
             IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
             IF (J1.EQ.LEN(LINE)) END_IND=J1
             IF (END_IND.GT.0) THEN
                J2=J2+1
                WORD=LINE(START_IND:END_IND)
                WORDSOUT(J2)=TRIM(WORD)
                START_IND=0
                END_IND=0
             ENDIF
          ENDIF
          J1=J1+1
      ENDDO
  END SUBROUTINE READ_LINE

END MODULE AMBER_CHIR_INIT

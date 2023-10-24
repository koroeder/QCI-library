MODULE COMPARELIST
   IMPLICIT NONE
   CONTAINS
      !comparing two sorted lists of atoms
      SUBROUTINE COMPARE_2LISTS(LIST1,LIST2,LENGTH,IDENTICALT,HIGH1T)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: LENGTH                           !length of lists
         INTEGER, INTENT(IN) :: LIST1(LENGTH), LIST2(LENGTH)     !list of priorities
         LOGICAL, INTENT(OUT) :: IDENTICALT                      !are the lists identical
         LOGICAL, INTENT(OUT) :: HIGH1T                          !higher priority for list 1
         INTEGER :: J1
         
         IDENTICALT = .FALSE.
         HIGH1T = .FALSE.
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
         RETURN
      END SUBROUTINE COMPARE_2LISTS

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

END MODULE COMPARELIST
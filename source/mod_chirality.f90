MODULE CHIRALITY
   USE QCIPREC
   USE QCIKEYS, ONLY: NATOMS
   USE AMBER_CONSTRAINTS, ONLY: NRES, NBOND, BONDS, RESNAMES, RES_START, RES_END, AMBER_NAMES, ELEMENT
   IMPLICIT NONE
   INTEGER, SAVE :: NCHIRAL = -1
   !!DEBUG
   INTEGER, SAVE :: NCHIRALBANDS = 0
   CHARACTER(LEN=4) :: NBAND
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
   LOGICAL, ALLOCATABLE :: POTENTIAL_SWAPT(:)                ! Can we use a swap of atoms to change the chirality later?
   INTEGER, PARAMETER :: NSWAPCUT = 4                        ! Cut off for largest swap size - currently this is capped at methyl group size.  
   INTEGER, ALLOCATABLE :: SWAPGROUPS(:,:,:)                   ! Atom ids to be swapped (dim1: chiral centre, dim2: groups, dim3: atoms to be swapped)

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
         IF (ALLOCATED(POTENTIAL_SWAPT)) DEALLOCATE(POTENTIAL_SWAPT)
         IF (ALLOCATED(SWAPGROUPS)) DEALLOCATE(SWAPGROUPS)
      END SUBROUTINE DEALLOC_CHIR_INTERNALS

      SUBROUTINE GET_ACTIVE_CHIRAL_CENTRES(NACTCHIRAL,ACTIVE_CENTRES)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE QCIKEYS, ONLY: NATOMS
         INTEGER, INTENT(OUT) :: NACTCHIRAL
         LOGICAL, INTENT(OUT) :: ACTIVE_CENTRES(NCHIRAL)
         LOGICAL :: ACTIVE
         INTEGER :: J1, J2

         ACTIVE_CENTRES(1:NCHIRAL) = .FALSE.
         NACTCHIRAL = 0

         DO J1=1,NCHIRAL
            ACTIVE = .TRUE.
            DO J2=1,5
               IF (.NOT.ATOMACTIVE(CHIR_INFO(J1,J2))) ACTIVE = .FALSE.
            END DO
            IF (ACTIVE) THEN
               ACTIVE_CENTRES(J1) = .TRUE.
               NACTCHIRAL = NACTCHIRAL + 1
            END IF
         END DO
      END SUBROUTINE GET_ACTIVE_CHIRAL_CENTRES

      SUBROUTINE CHECK_SINGLE_CHIRAL_CENTRE(CHIRALGROUP,XYZ)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE QCIKEYS, ONLY: NATOMS, DEBUG, NIMAGES, ATOMS2RES
         USE AMBER_CONSTRAINTS, ONLY: NBOND, BONDS, ELEMENT
         INTEGER, INTENT(IN) :: CHIRALGROUP
         REAL(KIND = REAL64), INTENT(INOUT) :: XYZ((3*NATOMS)*(NIMAGES+2))
         INTEGER :: CHIRALCENTRE
         INTEGER :: J1, J2, J3, IDX1, IDX2
         REAL(KIND = REAL64) :: NEIGHBOUR_COORDS(12), CENTRE_COORDS(3)
         REAL(KIND = REAL64) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
         INTEGER :: ATOMID
         INTEGER :: NCHANGE, CHANGEAT(NATOMS), THISIMAGE, PREVIMAGE, OFFSET, OFFSET2
         LOGICAL :: CENTREACTIVE
         LOGICAL :: THIS_SR, PREV_SR
         INTEGER :: NCONSTPERATOM(4), SORTED_IDX(4), SORTED_NB(4), SORTED_ELS(4), NB, EL
         INTEGER :: SWITCH1, SWITCH2

         !we set the number initially to -1 to discount the bond to the chiral centre
         NCONSTPERATOM(1:4) = -1
         !get the number of bonds to active atoms for the four connected atoms
         DO J1=1,NBOND
            IDX1 = BONDS(J1,1)
            IDX2 = BONDS(J1,2)
            IF (ATOMACTIVE(IDX1).AND.ATOMACTIVE(IDX2)) THEN
               DO J2=1,4
                  ATOMID = CHIR_INFO(CHIRALGROUP,J2+1)
                  IF ((IDX1.EQ.ATOMID).OR.(IDX2.EQ.ATOMID)) THEN
                     NCONSTPERATOM(J2) = NCONSTPERATOM(J2) + 1
                  END IF
               END DO
            END IF
         END DO

         !get a priority list based on (1) number of bonds (the lower the better) and (2) the element
         SORTED_IDX(1:4) = -1
         SORTED_NB(1:4) = 100
         SORTED_ELS(1:4) = 100

         DO J1=1,4
            IDX1 = CHIR_INFO(CHIRALGROUP,J1+1)
            NB = NCONSTPERATOM(J1)
            EL = ELEMENT(IDX1)
            DO J2=1,4
               IF ((NB.LT.SORTED_NB(J2)).OR.((NB.EQ.SORTED_NB(J2)).AND.(EL.LT.SORTED_ELS(J2)))) THEN
                  DO J3=4,J2+1,-1
                     SORTED_NB(J3)=SORTED_NB(J3-1)
                     SORTED_ELS(J3)=SORTED_ELS(J3-1)
                     SORTED_IDX(J3)=SORTED_IDX(J3-1)
                  END DO
                  SORTED_NB(J2) = NB
                  SORTED_ELS(J2) = EL
                  SORTED_IDX(J2) = IDX1
                  EXIT
               END IF
            END DO
         END DO
         SWITCH1 = SORTED_IDX(1)
         SWITCH2 = SORTED_IDX(2)

         CHIRALCENTRE = CHIR_INFO(CHIRALGROUP,1)
         WRITE(*,*) " check_single_centre> Checking chirality for atom ", CHIRALCENTRE
         WRITE(*,*) " check_single_centre> Atoms attached: ", SORTED_IDX
         WRITE(*,*) " check_single_centre> Bonds to active atoms: ", SORTED_NB
         WRITE(*,*) " check_single_centre> Elements: ", SORTED_ELS
         PREV_SR = .FALSE.
         !apply the actual check            
         DO J1=1,NIMAGES+2
            CENTRE_COORDS(1) = XYZ(3*NATOMS*(J1-1)+3*CHIRALCENTRE-2)
            CENTRE_COORDS(2) = XYZ(3*NATOMS*(J1-1)+3*CHIRALCENTRE-1)
            CENTRE_COORDS(3) = XYZ(3*NATOMS*(J1-1)+3*CHIRALCENTRE)

            DO J2=1,4
               ATOMID = CHIR_INFO(CHIRALGROUP,J2+1)
               NEIGHBOUR_COORDS(3*J2-2) = XYZ(3*NATOMS*(J1-1)+3*ATOMID-2)
               NEIGHBOUR_COORDS(3*J2-1) = XYZ(3*NATOMS*(J1-1)+3*ATOMID-1)
               NEIGHBOUR_COORDS(3*J2)   = XYZ(3*NATOMS*(J1-1)+3*ATOMID)
            END DO

            ! if this is the first image (i.e. the starting point), set SR for the current group
            IF (J1.EQ.1) THEN
               PREV_SR = ASSIGNMENT_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
               CYCLE
            END IF
            ! S/R assignment for current image
            THIS_SR = ASSIGNMENT_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
            IF (THIS_SR.NEQV.PREV_SR) THEN
               WRITE(*,*) " check_single_centre> Atom ", CHIRALCENTRE, " image ", J1, " chirality changed"
               WRITE(*,*) " check_single_centre> Attempting to switch atom pair: ", SWITCH1, SWITCH2
               CALL SWITCH_PAIR(CHIRALCENTRE,SWITCH1,SWITCH2,XYZ(3*NATOMS*(J1-1)+1:3*NATOMS*J1))                  
            END IF
         END DO
      END SUBROUTINE CHECK_SINGLE_CHIRAL_CENTRE

      SUBROUTINE SWITCH_PAIR(CENTRE,IDX1,IDX2,COORDS)
         USE QCIKEYS, ONLY: NATOMS
         USE HELPER_FNCTS, ONLY: NORM_VEC
         INTEGER, INTENT(IN) :: CENTRE
         INTEGER, INTENT(IN) :: IDX1, IDX2
         REAL(KIND = REAL64), INTENT(INOUT) :: COORDS(3*NATOMS)
         REAL(KIND = REAL64) :: C(3), AT1(3), AT2(3)
         REAL(KIND = REAL64) :: LEN1, LEN2, VEC1(3), VEC2(3), N1(3), N2(3)

         !atom positions
         C(1:3) = COORDS(3*(CENTRE-1)+1:3*(CENTRE-1)+3)
         AT1(1:3) = COORDS(3*(IDX1-1)+1:3*(IDX1-1)+3)
         AT2(1:3) = COORDS(3*(IDX2-1)+1:3*(IDX2-1)+3)
         !position vectors
         VEC1(1:3) = AT1(1:3) - C(1:3)
         VEC2(1:3) = AT2(1:3) - C(1:3)
         !get normalised vectors
         CALL NORM_VEC(VEC1,N1,LEN1)
         CALL NORM_VEC(VEC2,N2,LEN2)
         !now switch their positions
         !for atom1 we place it at the distance of C to atom1 along n2, and vice versa for atom2
         COORDS(3*(IDX1-1)+1:3*(IDX1-1)+3) = C(1:3) + LEN1*N2(1:3)
         COORDS(3*(IDX2-1)+1:3*(IDX2-1)+3) = C(1:3) + LEN2*N1(1:3)
      END SUBROUTINE SWITCH_PAIR

      SUBROUTINE SWITCH_SMALL_GROUPS(XYZ,CENTRE,IMAGE)
         USE QCIKEYS, ONLY: NATOMS, NIMAGES
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE HELPER_FNCTS, ONLY: NORM_VEC, CROSS_PROD, DOTP, EUC_NORM
         REAL(KIND = REAL64), INTENT(INOUT) :: XYZ((3*NATOMS)*(NIMAGES+2))
         INTEGER, INTENT(IN) :: CENTRE, IMAGE
         INTEGER :: CHIRALCENTRE
         INTEGER :: ATOMS1(NSWAPCUT), ATOMS2(NSWAPCUT)
         REAL(KIND = REAL64) :: C(3), AT1(3), AT2(3), VEC1(3), VEC2(3)
         REAL(KIND = REAL64) :: N1(3), N2(3), LEN1, LEN2, ROTAX(3)
         REAL(KIND = REAL64) :: COSTH, SINTH, ANGLE, RX, RY, RZ, DP
         REAL(KIND = REAL64) :: ATOMX(3)
         INTEGER :: ATID, J1

         CHIRALCENTRE = CHIR_INFO(CENTRE,1)
         ATOMS1(1:NSWAPCUT) = SWAPGROUPS(CENTRE,1,1:NSWAPCUT)
         ATOMS2(1:NSWAPCUT) = SWAPGROUPS(CENTRE,2,1:NSWAPCUT)
         
         !get centre coordinates
         C(1) = XYZ(3*NATOMS*(IMAGE-1)+3*CHIRALCENTRE-2)
         C(2) = XYZ(3*NATOMS*(IMAGE-1)+3*CHIRALCENTRE-1)
         C(3) = XYZ(3*NATOMS*(IMAGE-1)+3*CHIRALCENTRE)

         !get bonded atoms
         AT1(1) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS1(1)-2)
         AT1(2) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS1(1)-1)
         AT1(3) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS1(1))
         AT2(1) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS2(1)-2)
         AT2(2) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS2(1)-1)
         AT2(3) = XYZ(3*NATOMS*(IMAGE-1)+3*ATOMS2(1))
         
         !position vectors
         VEC1(1:3) = AT1(1:3) - C(1:3)
         VEC2(1:3) = AT2(1:3) - C(1:3)  

         !get normalised vectors
         CALL NORM_VEC(VEC1,N1,LEN1)
         CALL NORM_VEC(VEC2,N2,LEN2)

         !get perpendicular axis for rotation
         ROTAX = CROSS_PROD(N1,N2)
         RX = ROTAX(1); RY = ROTAX(2); RZ = ROTAX(3)
         !get the rotation angle N1->N2
         COSTH = DOTP(3,N1,N2)
         SINTH = EUC_NORM(ROTAX)
         ANGLE = ATAN2(SINTH,COSTH)

         DO J1=1,NSWAPCUT
            ATID = ATOMS1(J1)
            IF (ATID.NE.-1) THEN
               !coordinates for atoms
               ATOMX(1) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID-2)
               ATOMX(2) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID-1)
               ATOMX(3) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID)
               ! recentering using chiral centre as reference
               ATOMX(1:3) = ATOMX(1:3) - C(1:3)
               ! apply rotation
               DP = DOTP(3,ATOMX,ROTAX)
               ATOMX(1) = RX*DP*(1.0D0-COSTH) + ATOMX(1)*COSTH + (-RZ*ATOMX(2) + RY*ATOMX(3))*SINTH
               ATOMX(2) = RX*DP*(1.0D0-COSTH) + ATOMX(2)*COSTH + ( RZ*ATOMX(1) - RX*ATOMX(3))*SINTH
               ATOMX(3) = RX*DP*(1.0D0-COSTH) + ATOMX(3)*COSTH + (-RY*ATOMX(1) + RX*ATOMX(2))*SINTH
               !write new coordinates back to XYZ array
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID-2) =  ATOMX(1) + C(1)
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID-1) =  ATOMX(2) + C(2)
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID) =  ATOMX(3) + C(3)
            END IF
         END DO

         !get the rotation angle N2->N1
         COSTH = DCOS(-ANGLE)
         SINTH = DSIN(-ANGLE)
         DO J1=1,NSWAPCUT
            ATID = ATOMS2(J1)
            IF (ATID.NE.-1) THEN
               !coordinates for atoms
               ATOMX(1) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID-2)
               ATOMX(2) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID-1)
               ATOMX(3) = XYZ(3*NATOMS*(IMAGE-1)+3*ATID)
               ! recentering using chiral centre as reference
               ATOMX(1:3) = ATOMX(1:3) - C(1:3)
               ! apply rotation
               DP = DOTP(3,ATOMX,ROTAX)
               ATOMX(1) = RX*DP*(1.0D0-COSTH) + ATOMX(1)*COSTH + (-RZ*ATOMX(2) + RY*ATOMX(3))*SINTH
               ATOMX(2) = RX*DP*(1.0D0-COSTH) + ATOMX(2)*COSTH + ( RZ*ATOMX(1) - RX*ATOMX(3))*SINTH
               ATOMX(3) = RX*DP*(1.0D0-COSTH) + ATOMX(3)*COSTH + (-RY*ATOMX(1) + RX*ATOMX(2))*SINTH
               !write new coordinates back to XYZ array
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID-2) =  ATOMX(1) + C(1)
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID-1) =  ATOMX(2) + C(2)
               XYZ(3*NATOMS*(IMAGE-1)+3*ATID) =  ATOMX(3) + C(3)
            END IF
         END DO
         WRITE(*,*) "Switchable group start"
      END SUBROUTINE SWITCH_SMALL_GROUPS


      SUBROUTINE CHIRALITY_CHECK(XYZ)
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE MOD_INTCOORDS, ONLY: WRITE_BAND
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
                  !!DEBUG
                  NCHIRALBANDS = NCHIRALBANDS + 1
                  WRITE(NBAND,'(I4)') NCHIRALBANDS
                  CALL WRITE_BAND("chirality_issue."//TRIM(ADJUSTL(NBAND))//"a.xyz")
                  
                  !!END DEBUG

                  IF (POTENTIAL_SWAPT(J1)) THEN
                     WRITE(*,*) "                  Located a swappable pair - fixing chirality by switching small groups."
                     CALL SWITCH_SMALL_GROUPS(XYZ,J1,J3)
                  ELSE
                     WRITE(*,*) "                  Not a swappable pair - Using previous image coordinates."
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
                  CALL WRITE_BAND("chirality_issue."//TRIM(ADJUSTL(NBAND))//"b.xyz")

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
            WRITE(*,*) "Centre: ", J1, " atom id: ", CHIR_INFO(J1,1), " atoms attached: ", CHIR_INFO(J1,2:5)
         ENDDO
         WRITE(*,*) " find_chiral_centres> Completed chiral centre detection, number of chiral centres: ", NCHIRAL
         
         CALL FIND_SWAPPABLE_GROUPS()
         RETURN    
      END SUBROUTINE FIND_CHIRAL_CENTRES

      SUBROUTINE FIND_SWAPPABLE_GROUPS()
         INTEGER :: J1, J2
         INTEGER :: GROUPSIZE(4)
         INTEGER :: GROUPIDS(4,NSWAPCUT)
         INTEGER :: NSMALL, SMALLEST1, SMALLEST2, S1ID, S2ID, NSWAPPABLE

         ! allocate the potential swap variables
         ALLOCATE(POTENTIAL_SWAPT(NCHIRAL),SWAPGROUPS(NCHIRAL,2,NSWAPCUT))

         NSWAPPABLE = 0
         POTENTIAL_SWAPT(1:NCHIRAL) = .FALSE.
         DO J1=1,NCHIRAL
            CALL GET_ATTACHED_SIZE(J1,GROUPSIZE,GROUPIDS)
            NSMALL = 0
            DO J2=1,4
               IF (GROUPSIZE(J2).LE.NSWAPCUT) NSMALL=NSMALL+1
            END DO
            IF (NSMALL.GE.2) THEN
               POTENTIAL_SWAPT(J1) = .TRUE.
               NSWAPPABLE = NSWAPPABLE + 1
               SMALLEST1 = NSWAPCUT + 1
               SMALLEST2= NSWAPCUT + 1
               S1ID = 0
               S2ID = 0
               DO J2=1,4
                  IF (GROUPSIZE(J2).LT.SMALLEST1) THEN
                     SMALLEST2 = SMALLEST1
                     S2ID = S1ID
                     SMALLEST1 = GROUPSIZE(J2)
                     S1ID = J2 
                  ELSE IF (GROUPSIZE(J2).LT.SMALLEST2) THEN
                     SMALLEST2 = GROUPSIZE(J2)
                     S2ID = J2 
                  END IF
               END DO
               SWAPGROUPS(J1,1,1:NSWAPCUT) = GROUPIDS(S1ID,1:NSWAPCUT)
               SWAPGROUPS(J1,2,1:NSWAPCUT) = GROUPIDS(S2ID,1:NSWAPCUT)
            END IF
         END DO
         WRITE(*,*) " find_swappable_groups> Found ", NSWAPPABLE, " chiral centres with at least two bonded groups smaller than NSWAPCUT=",NSWAPCUT
      END SUBROUTINE FIND_SWAPPABLE_GROUPS

      SUBROUTINE GET_ATTACHED_SIZE(CENTRE,GROUPSIZE,GROUPIDS)
         INTEGER, INTENT(IN) :: CENTRE
         INTEGER, INTENT(OUT) :: GROUPSIZE(4)
         INTEGER, INTENT(OUT) :: GROUPIDS(4,NSWAPCUT)
         INTEGER :: ATID, ATID2, J1, J2, J3
         INTEGER :: TOTAL_BONDED

         GROUPSIZE(1:4) = 0
         GROUPIDS(1:4,1:NSWAPCUT) = -1
         DO J1=1,4
            ATID = CHIR_INFO(CENTRE,J1+1)
            ! check if it is a hydrogen (There is no other atom here, as we already preselected chiral centres)
            IF (NBONDED(ATID).EQ.1) THEN
               GROUPSIZE(J1) = 1
               GROUPIDS(J1,1) = ATID
            ! otherwise, we need to actually find the group size
            ELSE
               TOTAL_BONDED = 1
               ! iterate over all the atoms bonded to the current atom under consideration
               ! we add the number of atoms bonded
               ! for H, we get +1, i.e. accounting for the atom
               ! for anything else, we get 2+, accounting for the atom+other bonded atoms
               DO J2=1,NBONDED(ATID)
                  ATID2 = BONDEDATS(ATID,J2)
                  TOTAL_BONDED = TOTAL_BONDED + NBONDED(ATID2)
               END DO
               !we are counting the chiral centre here as well, so we subtract the number of bonds from the chiral centre to correct for this 
               TOTAL_BONDED = TOTAL_BONDED - 4
               !note this test only works for small NSWAPCUT sizes!!!
               GROUPSIZE(J1) = TOTAL_BONDED
               IF (TOTAL_BONDED.LE.NSWAPCUT) THEN
                  GROUPIDS(J1,1) = ATID
                  DO J3=2,TOTAL_BONDED
                     GROUPIDS(J1,J3) = BONDEDATS(ATID,J3)  
                  END DO
               END IF
            END IF
         END DO
      END SUBROUTINE GET_ATTACHED_SIZE

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
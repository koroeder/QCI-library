MODULE QCIPERMDIST
   USE QCIPREC
   IMPLICIT NONE
   ! number of permutational groups
   INTEGER :: NPERMGROUP = 0
   ! back up for NPERMGROUP
   INTEGER :: NPERMGROUPBACK = 0
   ! maximum number of subsets (in OPTIM this is set to five as default but we actually fix it later to the actual size)
   INTEGER :: MAXNSETS = 5
   ! arrays to store permutational information
   INTEGER, ALLOCATABLE :: NPERMSIZE(:)
   INTEGER, ALLOCATABLE :: PERMGROUP(:)
   INTEGER, ALLOCATABLE :: NSETS(:)
   INTEGER, ALLOCATABLE :: SETS(:,:,:)
   ! backup arrays for the groups 
   INTEGER, ALLOCATABLE :: NPERMSIZEBACK(:)
   INTEGER, ALLOCATABLE :: PERMGROUPBACK(:)
   INTEGER, ALLOCATABLE :: NSETSBACK(:)
   INTEGER, ALLOCATABLE :: SETSBACK(:,:,:)
   ! same arrays for QCI - QUERY: are these really needed?
   INTEGER, ALLOCATABLE :: NPERMSIZEQCI(:)
   INTEGER, ALLOCATABLE :: PERMGROUPQCI(:)
   INTEGER, ALLOCATABLE :: NSETSQCI(:)
   INTEGER, ALLOCATABLE :: SETSQCI(:,:,:)   
   ! best permutation
   INTEGER, ALLOCATABLE :: BESTPERM(:) 

   INTEGER, ALLOCATABLE :: ORDERNUM(:)

   ! common constraints
   INTEGER :: NCOMMONMAX = -1
   INTEGER, ALLOCATABLE :: NCONCOMMON(:)
   INTEGER, ALLOCATABLE :: CONCOMMON(:,:)

   !variables in keys
   INTEGER :: LOCALPERMNEIGH = 11
   INTEGER :: LOCALPERMMAXSEP = 3

   LOGICAL :: PERMDISTINIT = .FALSE.
   LOGICAL :: PERMGUESS = .FALSE.
   LOGICAL :: LPERMOFF=.FALSE.

   LOGICAL :: PBETTER = .FALSE.
    
   REAL(KIND=REAL64) :: LOCALPERMCUT = 0.5D0
   REAL(KIND=REAL64) :: LOCALPERMCUT2 = 5.0D0
   ! is this ever used or initialised?
   REAL(KIND=REAL64) :: LOCALPERMCUTINC

   REAL(KIND=REAL64) :: PDISTANCE
   REAL(KIND=REAL64) :: NOPDISTANCE
    
   ! minperm variables
   ! Save the largest arrays between iterations to reduce allocations.
   ! cc, kk: Sparse matrix of distances
   INTEGER, ALLOCATABLE :: CC(:), KK(:)

   !start and end indices and active groups
   LOGICAL, ALLOCATABLE :: GROUPACTIVE(:)
   INTEGER, ALLOCATABLE :: STARTGROUP(:), ENDGROUP(:)


   CONTAINS

      SUBROUTINE ALLOC_QCIPERM(NATOMS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS

         CALL DEALLOC_QCIPERM()
         ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,NPERMGROUP,MAXNSETS))
         ! Backup arrays in case QCIPERM is used.
         ALLOCATE(NPERMSIZEBACK(3*NATOMS),PERMGROUPBACK(3*NATOMS),NSETSBACK(3*NATOMS),SETSBACK(NATOMS,NPERMGROUP,MAXNSETS))
         ! The above dimensions were fixed at NATOMS because:
         ! (a) Atoms were not allowed to appear in more than one group.
         ! (b) The maximum number of pair exchanges associated with a group is three.
         !
         ! However, for flexible water models we need to exchange all waters,
         ! and we can exchange H's attached to the same O. The dimension required
         ! becomes 3*NATOMS
         ALLOCATE(BESTPERM(NATOMS))
         !common constraints
         ALLOCATE(NCONCOMMON(NPERMGROUP))
         !indices
         ALLOCATE(STARTGROUP(NPERMGROUP),ENDGROUP(NPERMGROUP),GROUPACTIVE(NPERMGROUP))
      END SUBROUTINE ALLOC_QCIPERM

      SUBROUTINE DEALLOC_QCIPERM()
         IMPLICIT NONE
         IF (ALLOCATED(NPERMSIZE)) DEALLOCATE(NPERMSIZE)
         IF (ALLOCATED(PERMGROUP)) DEALLOCATE(PERMGROUP)
         IF (ALLOCATED(NSETS)) DEALLOCATE(NSETS)
         IF (ALLOCATED(SETS)) DEALLOCATE(SETS)
         IF (ALLOCATED(NPERMSIZEBACK)) DEALLOCATE(NPERMSIZEBACK)
         IF (ALLOCATED(PERMGROUPBACK)) DEALLOCATE(PERMGROUPBACK)
         IF (ALLOCATED(NSETSBACK)) DEALLOCATE(NSETSBACK)
         IF (ALLOCATED(SETSBACK)) DEALLOCATE(SETSBACK)
         IF (ALLOCATED(BESTPERM)) DEALLOCATE(BESTPERM)
         IF (ALLOCATED(NCONCOMMON)) DEALLOCATE(NCONCOMMON)
         IF (ALLOCATED(STARTGROUP)) DEALLOCATE(STARTGROUP)
         IF (ALLOCATED(ENDGROUP)) DEALLOCATE(ENDGROUP)
         IF (ALLOCATED(GROUPACTIVE)) DEALLOCATE(GROUPACTIVE)
      END SUBROUTINE DEALLOC_QCIPERM

      SUBROUTINE INIT_PERMALLOW(NATOMS)
         USE QCIFILEHANDLER, ONLY: GETUNIT
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         LOGICAL :: PERMFILE
         INTEGER :: NDUMMY, J1, J2, J3
         INTEGER :: PERMUNIT

         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         IF (PERMFILE) THEN
            PERMUNIT=GETUNIT()
            OPEN(UNIT=PERMUNIT,FILE='perm.allow',STATUS='OLD')
            READ(PERMUNIT,*) NPERMGROUP
            NPERMGROUPBACK=NPERMGROUP
            CALL ALLOC_QCIPERM(NATOMS)
            NPERMSIZE(:) = 0
            PERMGROUP(:) = 0
            NSETS(:) = 0
            SETS(:,:,:) = 0
            

            NDUMMY = 1
            DO J1=1,NPERMGROUP
               READ(PERMUNIT,*) NPERMSIZE(J1),NSETS(J1)
               ! Sanity checks!
               IF (NSETS(J1).GT.MAXNSETS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1), ' is > ', MAXNSETS
                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
               READ(PERMUNIT,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J1,J2), &
                                J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1), J2=1,NSETS(J1))
               STARTGROUP(J1)=NDUMMY
               GROUPACTIVE(J1)=.FALSE.              
               NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDGROUP(J1)=NDUMMY-1
            ENDDO
            CLOSE(PERMUNIT)
            MAXNSETS=SIZE(SETS,2)
         ELSE 
            WRITE(*,*) " init_permallow> Cannot find perm.allow file"
            STOP
         END IF
      END SUBROUTINE INIT_PERMALLOW

      SUBROUTINE UPDATE_ACTIVE_PERMGROUPS()
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         IMPLICIT NONE
         INTEGER :: NDUMMY
         INTEGER :: J1, J2
         LOGICAL :: CURRENTGROUP

         NDUMMY = 1
         DO J1=1,NPERMGROUP
            IF (.NOT.GROUPACTIVE(J1)) THEN
               CURRENTGROUP = .TRUE.
               DO J2=1,NPERMSIZE(J1)
                  IF (.NOT.ATOMACTIVE(PERMGROUP(NDUMMY+J2-1))) THEN
                     CURRENTGROUP = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (CURRENTGROUP) GROUPACTIVE(J1) = .TRUE.
            END IF
            NDUMMY = NDUMMY + NPERMSIZE(J1)
         END DO
      END SUBROUTINE UPDATE_ACTIVE_PERMGROUPS

      SUBROUTINE CHECK_PERM_BAND(PERMGROUPIDX, FIRSTATOM, REVERSET)
         USE QCIKEYS, ONLY: DEBUG, NATOMS, NIMAGES
         USE MOD_INTCOORDS, ONLY: XYZ
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: PERMGROUPIDX !number of permutational group
         INTEGER, INTENT(IN) :: FIRSTATOM !first atom id in permutational group
         LOGICAL, INTENT(IN) :: REVERSET  !stepping direction through images

         INTEGER :: J1, J2, IDX, FIRSTIMAGE, SECONDIMAGE
         INTEGER :: STARTIDX, ENDIDX, STEP
         REAL(KIND = REAL64) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
         INTEGER :: PERMP(NATOMS), NMOVEP
         REAL(KIND = REAL64) :: BOXLX,BOXLY,BOXLZ,RMATBEST(3,3)
         LOGICAL :: BULKT = .FALSE., TWOD = .FALSE.
         REAL(KIND=REAL64) :: DISTANCE, DIST2

         ! going forward through the images
         IF (.NOT.REVERSET) THEN
            STARTIDX = 1
            ENDIDX = NIMAGES
            STEP = 1
         ! going backward through the images
         ELSE
            STARTIDX = NIMAGES
            ENDIDX  = 1
            STEP = -1
         END IF

         DO J1=STARTIDX,ENDIDX,STEP 
            ! Test alignment with neighbouring image
            ! Including endpoints (J): (start) 1 - 2 - 3 - ... -  J-1 -   J  - J+1 - J+2  - ... - NIMAGES+1 - NIMAGES (finish)
            ! Excluding endpoints (J1):            1 - 2 - ... - J1-2 - J1-1 -  J1 - J1+1 - ... - NIMAGES
            ! We want to go forward and start from the starting image up to the last interpolation image,
            ! and reverse from the the finish coordinates up to the first interpolation image
            ! So for the forward direction we want to start with comparing J=1 to J=2 (start to the first interpolation image)
            ! up to comparing J=NIMAGES to J=NIMAGES+1 (second-to-last to last interpoaltion image).
            ! For the reverse direction we want to start with comparing J=NIMAGES+2 to J=NIMASGES+1 (finish to last interpolation image)
            ! down to comparing J=3 to J=2 (second to first interpolation image).
            ! We are running J1 from 1 to NIMAGES/NIMAGES to 1 and need to get the right coordinates for the above images.
            ! For the forward direction J=J1 and J=J1+1 for the two images, for the reverse it is J=J1+2 and J=J1+1.
            IF (STEP.EQ.1) THEN
               FIRSTIMAGE = J1
            ELSE
               FIRSTIMAGE = J1 + 2
            END IF 
            SECONDIMAGE = J1+1
            !coordinates for image 1
            COORDSB(1:3*NATOMS) = XYZ((3*NATOMS*(FIRSTIMAGE-1)+1):3*NATOMS*FIRSTIMAGE)
            !coordinates for image 2
            COORDSA(1:3*NATOMS) = XYZ((3*NATOMS*(SECONDIMAGE-1)+1):3*NATOMS*SECONDIMAGE)            

            CALL LOPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,RMATBEST,PERMGROUPIDX,NMOVEP,PERMP)
   
            IF (NMOVEP.GT.0)  THEN
               WRITE(*,*) ' check_perm_band> group ',PERMGROUPIDX,' alignment of images ',SECONDIMAGE,FIRSTIMAGE,' moves=',&
                          NMOVEP, ' permutations, distance =',DISTANCE    
            END IF

            !swap atoms if we found a better permutational alignment
            IF ((NMOVEP.GT.0).AND.PBETTER) THEN 
               COORDSA(1:3*NATOMS) = XYZ((3*NATOMS*(SECONDIMAGE-1)+1):3*NATOMS*SECONDIMAGE) 
               DO J2=1,NPERMSIZE(PERMGROUPIDX)
                  IDX = PERMGROUP(FIRSTATOM+J2-1)
                  IF (PERMP(IDX).NE.IDX) THEN
                     WRITE(*,*) ' check_perm_band> image ',SECONDIMAGE,' move atom ',PERMP(IDX),' to position ',IDX
                     COORDSA(3*(IDX-1)+1)=XYZ(3*NATOMS*(SECONDIMAGE-1)+3*(PERMP(IDX)-1)+1)
                     COORDSA(3*(IDX-1)+2)=XYZ(3*NATOMS*(SECONDIMAGE-1)+3*(PERMP(IDX)-1)+2)
                     COORDSA(3*(IDX-1)+3)=XYZ(3*NATOMS*(SECONDIMAGE-1)+3*(PERMP(IDX)-1)+3)
                  ENDIF
               ENDDO
               XYZ((3*NATOMS*(SECONDIMAGE-1)+1):3*NATOMS*SECONDIMAGE) = COORDSA(1:3*NATOMS)
            ENDIF
         END DO
      END SUBROUTINE CHECK_PERM_BAND

      SUBROUTINE CHECK_COMMON_CONSTR()
         USE QCICONSTRAINTS, ONLY: NCONPERATOM, CONLIST

         INTEGER :: J1, J2, J3, J4, PATOM1, PATOM2, PATOMTEST
         INTEGER :: NENTRY
         LOGICAL :: COMMONCONST, THISATOMFOUND

         NCONCOMMON(1:NPERMGROUP) = 0
         NENTRY = 1 ! counter for position of entries in permutational groups
         DO J1=1,NPERMGROUP
            PATOM1 = PERMGROUP(NENTRY)
            ! For each entry in constraint list of first permutable atom, 
            ! check if it exists for the second, 
            ! if so, check the third, etc.
            DO J2=1,NCONPERATOM(PATOM1)
               ! cycle over all atoms in the constraint list that are constraint with PATOM1
               PATOMTEST = CONLIST(PATOM1,J2)
               COMMONCONST = .TRUE.
               DO J3=2,NPERMSIZE(J1)
                  ! cycle over all atoms in permutational group
                  PATOM2=PERMGROUP(NENTRY+J3-1)
                  THISATOMFOUND = .FALSE.
                  DO J4=1,NCONPERATOM(PATOM2)
                     IF (CONLIST(PATOM2,J4).EQ.PATOMTEST) THEN
                        !we found a match to PATOM2
                        THISATOMFOUND = .TRUE.
                        EXIT
                     END IF
                  END DO
                  !upate commonconst based on whether we found patom2
                  COMMONCONST = COMMONCONST.AND.THISATOMFOUND
                  ! if the latest atom was not a common one, we can leave the loop
                  IF (.NOT.COMMONCONST) THEN
                     EXIT
                  END IF
               END DO
               ! if we found a common constraint, add it to our list
               IF (COMMONCONST) THEN
                  NCONCOMMON(J1) = NCONCOMMON(J1) + 1
               END IF
            END DO
            NENTRY=NENTRY+NPERMSIZE(J1)
            IF (NCONCOMMON(J1).GT.NCOMMONMAX) NCOMMONMAX=NCONCOMMON(J1)
         END DO

         ! allocate common constraint array
         IF (ALLOCATED(CONCOMMON)) DEALLOCATE(CONCOMMON)
         ALLOCATE(CONCOMMON(NPERMGROUP,NCOMMONMAX))

         ! reiterating over all groups, this time to save common constraints
         NENTRY = 1 ! counter for position of entries in permutational groups
         DO J1=1,NPERMGROUP 
            NCONCOMMON(J1)=0
            PATOM1 = PERMGROUP(NENTRY)
            ! For each entry in constraint list of first permutable atom, 
            ! check if it exists for the second, 
            ! if so, check the third, etc.
            DO J2=1,NCONPERATOM(PATOM1)
               ! cycle over all atoms in the constraint list that are constraint with PATOM1
               PATOMTEST = CONLIST(PATOM1,J2)
               COMMONCONST = .TRUE.
               DO J3=2,NPERMSIZE(J1)
                  ! cycle over all atoms in permutational group
                  PATOM2=PERMGROUP(NENTRY+J3-1)
                  THISATOMFOUND = .FALSE.
                  DO J4=1,NCONPERATOM(PATOM2)
                     IF (CONLIST(PATOM2,J4).EQ.PATOMTEST) THEN
                        !we found a match to PATOM2
                        THISATOMFOUND = .TRUE.
                        EXIT
                     END IF
                  END DO
                  !upate commonconst based on whether we found patom2
                  COMMONCONST = COMMONCONST.AND.THISATOMFOUND
                  ! if the latest atom was not a common one, we can leave the loop
                  IF (.NOT.COMMONCONST) THEN
                     EXIT
                  END IF
               END DO
               ! if we found a common constraint, add it to our list
               IF (COMMONCONST) THEN
                  NCONCOMMON(J1)=NCONCOMMON(J1)+1
                  CONCOMMON(J1,NCONCOMMON(J1))=PATOMTEST
               END IF
            END DO
            NENTRY=NENTRY+NPERMSIZE(J1)
         END DO

      END SUBROUTINE CHECK_COMMON_CONSTR

      !
      SUBROUTINE LOPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,RMATBEST,DOGROUP,NMOVE,NEWPERM)
         USE QCIKEYS, ONLY: NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, BULKT, TWOD, RIGID, STOCKT
         USE QCIPREC
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         USE QCIMINDIST, ONLY: FIND_ALIGNMENT
         IMPLICIT NONE
         ! input and output variables
         REAL(KIND = REAL64), INTENT(INOUT) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS) ! coordinates for structure A and B
         REAL(KIND = REAL64), INTENT(OUT) :: DISTANCE ! distance between A and B
         REAL(KIND = REAL64), INTENT(OUT) :: DIST2 ! distance squared between A and B
         REAL(KIND = REAL64), INTENT(OUT) :: RMATBEST(3,3) ! best rotational matrix
         INTEGER, INTENT(IN) :: DOGROUP ! do group number for QCI
         INTEGER, INTENT(OUT) :: NMOVE ! number of permutational moves
         INTEGER, INTENT(OUT) :: NEWPERM(NATOMS) ! new permutation of atoms

         ! local variables
         INTEGER :: MAXREGION
         REAL(KIND = REAL64) :: XA, YA, ZA, XB, YB, ZB
         REAL(KIND = REAL64) :: CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
         INTEGER :: NPERM
         LOGICAL :: PERMUTABLE(NATOMS), PERMUTABLE2(NATOMS) ! lists of permutable atoms
         INTEGER :: NDUMMY, NDUMMY2, J1, J2, J3, J4 !counters etc
         REAL(KIND = REAL64) :: DUMMY(3*NATOMS), DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), TEMPB(3*NATOMS) ! coordinates copies 
         REAL(KIND = REAL64) :: PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), SPDUMMYA(3*NATOMS), SPDUMMYB(3*NATOMS) ! coordinates copies used in the bipartite matching
         REAL(KIND = REAL64) :: XBEST(3*NATOMS)
         REAL(KIND = REAL64) :: DSUM, LDISTANCE, DWORST ! Sum of local distances, local distance for this group and distance from newmindist
         INTEGER :: PATOMS ! size of group currently under consideration
         INTEGER :: TRIED(NATOMS) ! atoms tried in current permutation set
         INTEGER :: LPERMBEST(NATOMS), LPERM(NATOMS), LPERMBESTATOM(NATOMS) ! tracking of local permutations
         INTEGER :: INGROUP(NATOMS) ! tracking what is in the group
         REAL(KIND = REAL64) :: DMEAN(NATOMS) !mean distance trackign used in bipartite matching
         INTEGER :: NDMEAN ! counter for dmean
         REAL(KIND = REAL64) :: XDUMMY, DA, DB ! average for DA, DB, which are the distances after permutation
         REAL(KIND = REAL64) :: ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3)
         LOGICAL :: DONE !are we finished in the matching routine?
         LOGICAL :: USEATOM
         INTEGER :: SORTLIST(NATOMS) ! tracking in bipartite routine
         INTEGER :: DLIST(NATOMS), NOTHER
         REAL(KIND = REAL64) :: LDBEST(NPERMGROUP), LDBESTATOM
         INTEGER :: NADDED
         INTEGER :: NCHOOSE1, NCHOOSE2, NCHOOSEB1, NCHOOSEB2
         INTEGER :: NORBIT1, NORBIT2, NORBITB1, NORBITB2
         REAL(KIND = REAL64) :: RMAT(3,3)
         INTEGER :: ALLPERM(NATOMS), SAVEPERM(NATOMS)
         REAL(KIND = REAL64) :: WORSTRAD

         !initialise various variables
         MAXREGION = 0
         CMXA=0.0D0;CMYA=0.0D0;CMZA=0.0D0;CMXB=0.0D0;CMYB=0.0D0;CMZB=0.0D0
         PERMUTABLE(1:NATOMS)=.FALSE.
         PERMUTABLE2(1:NATOMS)=.FALSE.
         NDUMMY=1

         ! Loops to assign permutable atoms for book keeping
         ! This block is needed to allow associated permutational atoms to ve recognised as permutable later.
         ! This is not necessary of the associated permutable atons appear in perm.allow in their own list,
         ! as for methy and amine hydrogens, but it is necessary for associated H atoms in phenylalanine,
         ! as these cannot permute without the associated C atoms, etc. The choice of primary and associated
         ! sets for phenylalanine is arbitrary.
         DO J1=1,NPERMGROUP
            DO J2=1,NPERMSIZE(J1)
               PERMUTABLE(PERMGROUP(NDUMMY+J2-1))=.TRUE.
               INGROUP(PERMGROUP(NDUMMY+J2-1))=J1
            ENDDO
            IF (NSETS(J1).GT.0) THEN
               DO J2=1,NPERMSIZE(J1)
                  DO J3=1,NSETS(J1)
                     PERMUTABLE2(SETS(PERMGROUP(NDUMMY+J2-1),J1,J3))=.TRUE.
                    ! PRINT '(A,I8,A)','lopermdist> Made atom ',SETS(PERMGROUP(NDUMMY+J2-1),J3),' permutable'       
                  ENDDO
               ENDDO
            ENDIF
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO
         ! Copy coordinates into dummy arrays
         DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
         DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)
         ! Initiate permutational tracking and other variables required for the bipartite matching routine below
         DO J1=1,NATOMS
            NEWPERM(J1)=J1
         ENDDO
         DSUM=0.0D0
         LOCALPERMNEIGH=MIN(LOCALPERMNEIGH,NATOMS)

         !  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
         !  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case,
         !  and is evaluated as a sum of local distances squared for permutable groups.
         !  We return to label 10 after every round of permutational/orientational alignment
         !  unless we have converged to the identity permutation.
         !
         !  The maximum number of pair exchanges associated with a group is two.
         NDUMMY = 1 ! Reset dummy counter

         DO J1=1,NPERMGROUP  ! Loop for bipartite matching
            ! are we doing this group?
            IF (DOGROUP.GT.0) THEN
               IF (J1.GT.DOGROUP) EXIT
               IF (J1.LT.DOGROUP) GOTO 864 ! for QCI test single groups - MUST increment NDUMMY on line 864!!
            ENDIF
            !get information for this group to be used
            PATOMS=NPERMSIZE(J1)
            LDBEST(J1)=1.0D100
            TRIED(1:NATOMS)=0
            DO J2=1,PATOMS
               LPERMBEST(J2)=J2
            ENDDO
            XA=0.0D0; YA=0.0D0; ZA=0.0D0
            XB=0.0D0; YB=0.0D0; ZB=0.0D0
            DO J2=1,PATOMS
              ! TRIED(NEWPERM(PERMGROUP(NDUMMY+J2-1)))=-1
               TRIED(PERMGROUP(NDUMMY+J2-1))=-1
               PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
               PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
               PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
              ! PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
              ! PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
              ! PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
               PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
               PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
               PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
               XA=XA+PDUMMYA(3*(J2-1)+1)
               YA=YA+PDUMMYA(3*(J2-1)+2)
               ZA=ZA+PDUMMYA(3*(J2-1)+3)
               XB=XB+PDUMMYB(3*(J2-1)+1)
               YB=YB+PDUMMYB(3*(J2-1)+2)
               ZB=ZB+PDUMMYB(3*(J2-1)+3)
            ENDDO
            XA=XA/PATOMS; YA=YA/PATOMS; ZA=ZA/PATOMS
            XB=XB/PATOMS; YB=YB/PATOMS; ZB=ZB/PATOMS
            SPDUMMYA(1:3*PATOMS)=PDUMMYA(1:3*PATOMS)
            SPDUMMYB(1:3*PATOMS)=PDUMMYB(1:3*PATOMS)  

            ! TRIED(J2) is 0 if atom J2 is eligible to be a neighbour, but has not
            ! yet been tried. It is -1 if it is ineligible, or has been tried and
            ! broke the alignment. It is +1 if it has been tried and did not break
            ! the alignment. It is -1 for atoms already in the set of permutable
            ! atoms in question. We add neighbours one at a time in order of 
            ! increasing distance from primary permutable set
            ! and test whether they break the alignment.
            
            ! initiate DMEAN and NDMEAN counter
            DMEAN(1:NATOMS)=10.0*LOCALPERMCUT2
            NDMEAN=0

            ! Make a sorted list of distance from the permuting atoms.
            ! DMEAN, SORTLIST, TRIED, PERMUTABLE, PERMUTABLE2 and DLIST entries refer to original
            ! atom labels. Use NEWPERM to find where they are in coordinate lists. 

            DO J2=1,NATOMS  ! outer1 loop
               USEATOM=.TRUE.
               IF (DOGROUP.GT.0) THEN
                  IF (.NOT.ATOMACTIVE(J2)) USEATOM=.FALSE. ! must not attempt to access atomactive unless DOGROUP.GT.0 - not allocated!
               ENDIF
               ! Don't allow members of the same permutational group
               ! to appear as reference neighbours.
               IF (TRIED(J2).EQ.-1) THEN
                  XDUMMY=1.0D9
                  CYCLE !go to the next iteration cycle (outer1)
               ELSE
                  IF (.NOT.USEATOM) THEN ! only use active atoms for QCI single group
                     XDUMMY=1.0D9
                  ELSE
                     DA = (XA-DUMMYA(3*(NEWPERM(J2)-1)+1))**2 + &
                          (YA-DUMMYA(3*(NEWPERM(J2)-1)+2))**2 + &
                          (ZA-DUMMYA(3*(NEWPERM(J2)-1)+3))**2
                     DB = (XB-DUMMYB(3*(J2-1)+1))**2 + &
                          (YB-DUMMYB(3*(J2-1)+2))**2 + &
                          (ZB-DUMMYB(3*(J2-1)+3))**2
                     XDUMMY=(SQRT(DA)+SQRT(DB))/2.0D0
                     IF (XDUMMY.GT.LOCALPERMCUT2) CYCLE !go to the next iteration cycle (outer1)
                  ENDIF
               ENDIF
               NDMEAN = NDMEAN + 1
               DO J3=1,NDMEAN ! NDMEAN is up to J2 loop loop1
                  DONE=.FALSE.
                  IF (XDUMMY.LT.DMEAN(J3)) THEN
                     ! Move the rest down.
                     DO J4=NDMEAN,J3+1,-1   !    J2,J3+1,-1
                        DMEAN(J4)=DMEAN(J4-1)
                        SORTLIST(J4)=SORTLIST(J4-1)
                     ENDDO
                     DMEAN(J3)=XDUMMY
                     SORTLIST(J3)=J2
                     DONE=.TRUE.
                  ENDIF
                  IF (DONE) EXIT
               ENDDO ! end of loop1
            ENDDO ! end of outer1 loop

            ! If the atoms in the permuting group are associated with obligatory permutations 
            ! then add all the members in these SETS. NDUMMY is already set for this group J1.
            ! This setup means that such atoms can appear in the list of candidates sorted by distance above.
            ! They will not be used with a tried value of 1.
            ! We still include everything in the cutoff, but it would be more efficient not to include these
            ! atoms in the sorted candidayes list

            ! J1 is the permutable group
            ! J2 is counting the sets of extra swaps for group J1
            ! J3 runs over the permutable atoms in group J1 - each associated set of permuting atoms has to be the same size
            NOTHER=0
            IF (NSETS(J1).GT.0) THEN
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1     
                     NOTHER=NOTHER+1
                     DLIST(NOTHER)=SETS(PERMGROUP(J3),J1,J2) 
                     TRIED(DLIST(NOTHER))=1
                  ENDDO
               ENDDO
            ENDIF 

71          CONTINUE
            PDUMMYA(1:3*PATOMS)=SPDUMMYA(1:3*PATOMS)
            PDUMMYB(1:3*PATOMS)=SPDUMMYB(1:3*PATOMS)
         
            LDBESTATOM=1.0D100
            NOTHER=0
            DO J2=1,NATOMS
               IF (TRIED(J2).EQ.1) THEN
                  NOTHER=NOTHER+1
                  DLIST(NOTHER)=J2
               ENDIF
            ENDDO

            DO J2=1,NATOMS ! loop outer2
               IF (DMEAN(J2).GT.LOCALPERMCUT2) GOTO 91
               IF (TRIED(SORTLIST(J2)).EQ.0) THEN
                  NOTHER=NOTHER+1
                  IF (NOTHER+PATOMS.GT.NATOMS) THEN
                     PRINT '(A,A,I6)', ' lopermdist> ERROR *** number of neighbours plus', &
                                       ' number of permutable atoms exceeds total for group ',J1
                     STOP
                  ENDIF
                  DLIST(NOTHER)=SORTLIST(J2)
                  EXIT !leave loop outer2
               ENDIF
            ENDDO ! end of outer2 loop         

            NADDED=1
            IF (PERMUTABLE(DLIST(NOTHER))) THEN
               NDUMMY2=1
               DO J2=1,INGROUP(DLIST(NOTHER))-1
                  NDUMMY2=NDUMMY2+NPERMSIZE(J2)
               ENDDO
               DO J2=1,NPERMSIZE(INGROUP(DLIST(NOTHER)))
                  IF (PERMGROUP(NDUMMY2+J2-1).EQ.DLIST(NOTHER-NADDED+1)) CYCLE
                  IF (TRIED(PERMGROUP(NDUMMY2+J2-1)).EQ.0) THEN
                     NOTHER=NOTHER+1
                     NADDED=NADDED+1
                     IF (NOTHER+PATOMS.GT.NATOMS) THEN
                        PRINT '(A,A,I6)', ' lopermdist> ERROR *** number of neighbours plus', & 
                              ' number of permutable atoms exceeds total for group ',J1
                        STOP
                     ENDIF
                     DLIST(NOTHER)=PERMGROUP(NDUMMY2+J2-1)
                  ELSE
                     PRINT '(A,I6,A)',' lopermdist> ERROR *** Partner atom ',DLIST(NOTHER),' has already been tried'
                     STOP
                  ENDIF
               ENDDO
            ENDIF

            DO J2=1,NOTHER
               PDUMMYA(3*(PATOMS+J2-1)+1)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+1)
               PDUMMYA(3*(PATOMS+J2-1)+2)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+2)
               PDUMMYA(3*(PATOMS+J2-1)+3)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+3)
              ! PDUMMYB(3*(PATOMS+J2-1)+1)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+1)
              ! PDUMMYB(3*(PATOMS+J2-1)+2)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+2)
              ! PDUMMYB(3*(PATOMS+J2-1)+3)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+3)
               PDUMMYB(3*(PATOMS+J2-1)+1)=DUMMYB(3*(DLIST(J2)-1)+1)
               PDUMMYB(3*(PATOMS+J2-1)+2)=DUMMYB(3*(DLIST(J2)-1)+2)
               PDUMMYB(3*(PATOMS+J2-1)+3)=DUMMYB(3*(DLIST(J2)-1)+3)
            ENDDO
         
            ! Save PDUMMYA and PDUMMYB for cycling over possible orbits in MYORIENT alignment.
            SPDUMMYA(3*PATOMS+1:3*(PATOMS+NOTHER))=PDUMMYA(3*PATOMS+1:3*(PATOMS+NOTHER))
            SPDUMMYB(3*PATOMS+1:3*(PATOMS+NOTHER))=PDUMMYB(3*PATOMS+1:3*(PATOMS+NOTHER))
            NCHOOSEB1=0
66          NCHOOSEB1=NCHOOSEB1+1
            NCHOOSEB2=0
31          NCHOOSEB2=NCHOOSEB2+1
            NCHOOSE1=0
65          NCHOOSE1=NCHOOSE1+1
            NCHOOSE2=0
30          NCHOOSE2=NCHOOSE2+1    

            ! Reset the coordinates of the PATOMS+NOTHER atoms in PDUMMYA and PDUMMYB
            ! to the subset of atoms from COORDSA and COORDSB.

            PDUMMYA(1:3*(PATOMS+NOTHER))=SPDUMMYA(1:3*(PATOMS+NOTHER))
            PDUMMYB(1:3*(PATOMS+NOTHER))=SPDUMMYB(1:3*(PATOMS+NOTHER))

            IF (.NOT.LPERMOFF) THEN ! don't reorient at all in this case, so that newmindist2 rotation refers to the original coordinates and we can rotate all of them
               CALL MYORIENT(PDUMMYA,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,PATOMS+NOTHER,DEBUG,ROTA,ROTINVA,STOCKT)
               PDUMMYA(1:3*(PATOMS+NOTHER))=DUMMY(1:3*(PATOMS+NOTHER))
               CALL MYORIENT(PDUMMYB,DUMMY,NORBITB1,NCHOOSEB1,NORBITB2,NCHOOSEB2,PATOMS+NOTHER,DEBUG,ROTB,ROTINVB,STOCKT)
               PDUMMYB(1:3*(PATOMS+NOTHER))=DUMMY(1:3*(PATOMS+NOTHER))
            ELSE
               NORBIT1=1
               NORBIT2=1
               NORBITB1=1
               NORBITB2=1
            ENDIF

            ! Optimimise permutational isomer for the standard orientation for the
            ! current choice of atoms from the possible orbits.
            !
            ! MINPERM does not change PDUMMYB and PDUMMYA.
            !
            ! Note that LDISTANCE is actually the distance squared. LDBEST also has dimensions of
            ! length squared.
            LDISTANCE=0.0D0
            IF (LPERMOFF) THEN
               CALL NEWMINDIST2(PDUMMYA,PDUMMYB,PATOMS+NOTHER,LDISTANCE,DEBUG,RMAT,CMXA,CMYA,CMZA,CMXB,CMYB,CMZB,DWORST)
               ! LDISTANCE=LDISTANCE**2 ! assumed to be squared below
               LDISTANCE=DWORST**2
               DO J2=1,PATOMS+NOTHER
                  LPERM(J2)=J2
               ENDDO
            ELSE
               CALL MINPERM(PATOMS+NOTHER, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD) 
            ENDIF

            DO J2=1,PATOMS
               IF (LPERM(J2).GT.PATOMS) THEN
                  LDISTANCE=1.0D300
                  EXIT
               ENDIF
            ENDDO

            DO J2=1,NOTHER
               IF (LPERM(PATOMS+J2).NE.PATOMS+J2) THEN
                  IF (.NOT.(PERMUTABLE(DLIST(J2)).OR.PERMUTABLE2(DLIST(J2)))) THEN
                     LDISTANCE=1.0D300
                  ENDIF
               ENDIF
            ENDDO

            ! Save the best permutation and local distance for this subset of atoms.
            ! NEWPERM and coordinates are only reset after all the cycles over orbits and NEWMINDIST.
            ! Hence we need to track a cumulative permutation and save the best current values.  
            IF (LDISTANCE.LT.LDBESTATOM) THEN
               LDBESTATOM=LDISTANCE
               LPERMBESTATOM(1:PATOMS)=LPERM(1:PATOMS)
            ENDIF

            IF (NCHOOSE2.LT.NORBIT2) GOTO 30
            IF (NCHOOSE1.LT.NORBIT1) GOTO 65
            IF (NCHOOSEB2.LT.NORBITB2) GOTO 31
            IF (NCHOOSEB1.LT.NORBITB1) GOTO 66           

            IF (SQRT(LDBESTATOM).GT.LOCALPERMCUT) THEN
               TRIED(DLIST(NOTHER))=-1
               IF (NADDED.GT.1) THEN
                  TRIED(DLIST(NOTHER-NADDED+1:NOTHER-1))=-1
               ENDIF
               GOTO 71
            ELSE
               TRIED(DLIST(NOTHER))=1
               IF (NADDED.GT.1) THEN
                  TRIED(DLIST(NOTHER-NADDED+1:NOTHER-1))=1
               ENDIF
               LDBEST(J1)=LDBESTATOM
               LPERMBEST(1:PATOMS)=LPERMBESTATOM(1:PATOMS)
            ENDIF

            ! Add the next eligible atom and try alignment again.
            ! Stop if we already have LOCALPERMNEIGH neighbours.
            IF (NOTHER.LT.LOCALPERMNEIGH) GOTO 71
    
91          CONTINUE ! jump here when there are no atoms left to try.   

            IF (LPERMOFF.AND.(PATOMS+NOTHER.GT.MAXREGION)) THEN
               MAXREGION=PATOMS+NOTHER
               PRINT '(A,I8,A)','lopermdist> with permutations off, setting size of maximum region to ',MAXREGION,' and saving orientation'  
               !
               ! align the whole system 
               ! RMAT refers to rotation of the second set of coordinates argument to newmindist2, which is PDUMMYB
               ! we have turned off myorient for LPERMOFF
               !
               ! Translate the B coordinates subset to the centre of coordinates of A.
               ! In newrotgeom cmxa, cmya, cmza are subtracted, the rotation is done, then the centre of coordinates is made coincident with A again
               !
               DO J2=1,NATOMS
                  TEMPB(3*(J2-1)+1)=COORDSB(3*(J2-1)+1)-CMXB+CMXA
                  TEMPB(3*(J2-1)+2)=COORDSB(3*(J2-1)+2)-CMYB+CMYA
                  TEMPB(3*(J2-1)+3)=COORDSB(3*(J2-1)+3)-CMZB+CMZA
               ENDDO
               CALL NEWROTGEOM(NATOMS,TEMPB,RMAT,CMXA,CMYA,CMZA)
            ENDIF

            ! We now have the best permutation for group J1 and standard orientations
            ! based upon all atoms belonging to the two possible orbits that appear
            ! for the standard alignment.
            LPERM(1:PATOMS)=LPERMBEST(1:PATOMS)

            ! Fill SAVEPERM with NEWPERM, which contains the current best permutation
            ! after the previous pass through J1
            SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)

            ! Update best permutation for atoms in subset J1, specified by PERMGROUP
            ! with offset NDUMMY (updated below after each pass through J1)
            DO J2=1,PATOMS
               SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERMBEST(J2)-1))
            ENDDO

            ! Update permutation of associated atoms, if any.
            ! We must do this as we go along, because these atoms could move in more than
            ! one permutational group now.
            IF (NSETS(J1).GT.0) THEN
               DO J2=1,PATOMS
                  DO J3=1,NSETS(J1)
                     SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J1,J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J1,J3)
                  ENDDO
               ENDDO
            ENDIF

            ! Save current optimal permutation in NEWPERM
            NEWPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
            DSUM=DSUM+SQRT(LDBEST(J1))

            ! Update NDUMMY, the cumulative offset for PERMGROUP
864         NDUMMY=NDUMMY+NPERMSIZE(J1)
         END DO ! End of loop for bipartite matching

         PBETTER=.FALSE.
         IF (DOGROUP.GT.0) THEN
            DISTANCE=SQRT(LDBEST(DOGROUP))
            NMOVE=0
            DO J2=1,NATOMS
               IF (NEWPERM(J2).NE.J2) THEN
                  NMOVE=NMOVE+1
               ENDIF
            ENDDO

            ! Test distance for COORDSA without permutation applied in NEWPERM
            DO J1=1,NATOMS
               DUMMYA(3*(J1-1)+1)=COORDSA(3*(J1-1)+1)
               DUMMYA(3*(J1-1)+2)=COORDSA(3*(J1-1)+2)
               DUMMYA(3*(J1-1)+3)=COORDSA(3*(J1-1)+3)
            ENDDO
            CALL NEWMINDIST2(COORDSB,DUMMYA,NATOMS,NOPDISTANCE,DEBUG,RMAT,CMXA,CMYA,CMZA,CMXB,CMYB,CMZB,DWORST)
            ! Test distance for COORDSA with permutation applied in NEWPERM
            DO J1=1,NATOMS
               DUMMYA(3*(J1-1)+1)=COORDSA(3*(NEWPERM(J1)-1)+1)
               DUMMYA(3*(J1-1)+2)=COORDSA(3*(NEWPERM(J1)-1)+2)
               DUMMYA(3*(J1-1)+3)=COORDSA(3*(NEWPERM(J1)-1)+3)
            ENDDO
            CALL NEWMINDIST2(COORDSB,DUMMYA,NATOMS,PDISTANCE,DEBUG,RMAT,CMXA,CMYA,CMZA,CMXB,CMYB,CMZB,DWORST)           
            IF (PDISTANCE.LT.NOPDISTANCE) PBETTER=.TRUE.
            RETURN
         ENDIF
         NMOVE=0
         DO J2=1,NATOMS
            IF (NEWPERM(J2).NE.J2) THEN
               NMOVE=NMOVE+1
            ENDIF
         ENDDO
         IF (DEBUG) PRINT '(A,I6)',' lopermdist> Total permutations for optimal alignment (will be applied)=',NMOVE

         ! NEWPERM(J1) is the atom that moves to position J1 to map COORDSA
         ! to the current best alignment. 
         ! This loop just appears to set SAVEPERM and ALLPERM equal to the current
         ! NEWPERM.

         ! Putting the ALLPERM(J1)=J1 into the second loop causes pgf90 to miscompile!!
         DO J1=1,NATOMS
            ALLPERM(J1)=J1
         ENDDO
         DO J1=1,NATOMS
            SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
         ENDDO
         ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
         ! At this point DUMMYA should not have changed from COORDSA, so we are
         ! putting COORDSA in DUMMY
         DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
         NPERM=0

         ! Update coordinates in DUMMYA to current best overall permutation using NEWPERM.
         ! We are doing this to operate with NEWPERMDIST in the next block.
         DO J3=1,NATOMS
            DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
            DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
            DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)
         
            IF (J3.NE.NEWPERM(J3)) THEN
               NPERM=NPERM+1
            ENDIF
         ENDDO

         DISTANCE=DSUM
         DIST2=DISTANCE**2
         ! Save current best overall distance, permuted version of COORDSA, and permutation.
         XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
         BESTPERM(1:NATOMS)=ALLPERM(1:NATOMS)
         ! At this point NEWPERM, ALLPERM, SAVEPERM, BESTPERM
         ! are all the same!
         CALL FIND_ALIGNMENT(NATOMS, DUMMYB, XBEST, DISTANCE, RMAT)
         IF (DEBUG) PRINT '(A,G20.10)',' lopermdist> after overall alignment distance=',DISTANCE
         RMATBEST(1:3,1:3)=RMAT(1:3,1:3)
         
         IF (LPERMOFF) THEN ! align COORDSB instead
            COORDSB(1:3*NATOMS)=TEMPB(1:3*NATOMS)
         ELSE
            COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!
         ENDIF
      END SUBROUTINE LOPERMDIST

      ! TODO add routines below
      ! Put permutational isomers into a standard orientation
      SUBROUTINE MYORIENT(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROT,ROTINV,STOCKT)
         USE QCIKEYS, ONLY: ORBITTOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND=REAL64), INTENT(INOUT) :: Q1(3*NATOMS)
         REAL(KIND=REAL64), INTENT(OUT) ::  T1(3*NATOMS)
         INTEGER, INTENT(OUT) :: NORBIT1, NORBIT2
         INTEGER, INTENT(IN) :: NCHOOSE1, NCHOOSE2
         LOGICAL, INTENT(IN) :: DEBUG, STOCKT
         REAL(KIND=REAL64), INTENT(OUT) :: ROT(3,3), ROTINV(3,3)

         REAL(KIND=REAL64) :: CUTOFF1
         REAL(KIND=REAL64) :: XVEC(3), YVEC(3), ZVEC(3) ! unit vectors for tracking
         REAL(KIND=REAL64) :: CMX, CMY, CMZ ! centre of mass
         INTEGER :: I, J1, JMAX1, J2, JMAX2 ! counters
         REAL(KIND=REAL64) :: DIST(NATOMS), DMAX, DMAX2
         REAL(KIND=REAL64) :: COST, SINT
         REAL(KIND=REAL64) :: RVEC(3), RDOTN, DUM1, PROJ
         REAL(KIND=REAL64) :: TXVEC(3), TYVEC(3), TZVEC(3)
         REAL(KIND=REAL64) :: T2(3*NATOMS)

         CUTOFF1=ORBITTOL

         ! Use unit vectors along the x, y, z axes to keep track of overall transformation matrix.
         XVEC(1)=1.0D0; XVEC(2)=0.0D0; XVEC(3)=0.0D0
         YVEC(1)=0.0D0; YVEC(2)=1.0D0; YVEC(3)=0.0D0
         ZVEC(1)=0.0D0; ZVEC(2)=0.0D0; ZVEC(3)=1.0D0

         CMX=0.0D0; CMY=0.0D0; CMZ=0.0D0

         ! Move centre of mass to the origin.
         DO I=1,NATOMS
            CMX=CMX+Q1(3*(I-1)+1)
            CMY=CMY+Q1(3*(I-1)+2)
            CMZ=CMZ+Q1(3*(I-1)+3)
         ENDDO
         CMX=CMX/NATOMS
         CMY=CMY/NATOMS
         CMZ=CMZ/NATOMS
         DO I=1,NATOMS
            Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
            Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
            Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
         ENDDO

         DMAX=-1.0D0
         NORBIT1=1

         JMAX1 = -1

         DO J1=1,NATOMS
            DIST(J1)=SQRT(Q1(3*(J1-1)+1)**2+Q1(3*(J1-1)+2)**2+Q1(3*(J1-1)+3)**2)
            IF (ABS(DIST(J1)-DMAX).LT.CUTOFF1) THEN
               NORBIT1=NORBIT1+1
               IF (NORBIT1.EQ.NCHOOSE1) THEN
                  JMAX1=J1
               ENDIF
            ELSE IF (DIST(J1).GT.DMAX) THEN
               DMAX=DIST(J1)
               NORBIT1=1
               JMAX1=J1
            ENDIF
         ENDDO 
         ! the following debug statement overcomes a runtime bug of uninitialised values - not sure why - likely a code optimisation somewhere
         IF (DEBUG) WRITE(*,*) " myorient> Atom chosen is ", J1, " with distance ", DMAX

         ! For tagged atoms the choice of the first atom matters if it belongs to an orbit of size > 1.
         IF ((ABS(Q1(3*(JMAX1-1)+1)).LT.1.0D-8).AND.(ABS(Q1(3*(JMAX1-1)+2)).LT.1.0D-8)) THEN
            IF (Q1(3*(JMAX1-1)+3).GT.0.0D0) THEN
               T1(1:3*NATOMS)=Q1(1:3*NATOMS)
            ELSE  ! rotate about the x axis DO NOT INVERT!!
               DO J1=1,NATOMS
                  T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)
                  T1(3*(J1-1)+2)=-Q1(3*(J1-1)+2)
                  T1(3*(J1-1)+3)=-Q1(3*(J1-1)+3)
               ENDDO
               YVEC(2)=-1.0D0
               ZVEC(3)=-1.0D0
            ENDIF
            GOTO 10
         ENDIF  

         ! For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!!
         COST=Q1(3*(JMAX1-1)+3)/DIST(JMAX1)
         SINT=SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)/DIST(JMAX1)

         ! Rotate atom JMAX1 onto the z axis.
         ! Rotate all the atoms through ANGLE about RVEC. Use rotation formula from Goldstein p. 165.
         RVEC(1)= Q1(3*(JMAX1-1)+2)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
         RVEC(2)=-Q1(3*(JMAX1-1)+1)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
         RVEC(3)=0.0D0

         !  ifort/64/2015/3/187 miscompiles the loop below and gives garbage without a print statement!!
         !  Now split into three loops to make ifort work.  
         DO J1=1,NATOMS
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            DUM1=RDOTN*(1.0D0-COST)
            T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)*COST + RVEC(1)*DUM1-(Q1(3*(J1-1)+2)*RVEC(3)-Q1(3*(J1-1)+3)*RVEC(2))*SINT
         ENDDO
         DO J1=1,NATOMS
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            DUM1=RDOTN*(1.0D0-COST)
            T1(3*(J1-1)+2)=Q1(3*(J1-1)+2)*COST + RVEC(2)*DUM1-(Q1(3*(J1-1)+3)*RVEC(1)-Q1(3*(J1-1)+1)*RVEC(3))*SINT
         ENDDO
         DO J1=1,NATOMS
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            DUM1=RDOTN*(1.0D0-COST)
            T1(3*(J1-1)+3)=Q1(3*(J1-1)+3)*COST + RVEC(3)*DUM1-(Q1(3*(J1-1)+1)*RVEC(2)-Q1(3*(J1-1)+2)*RVEC(1))*SINT
         ENDDO

         ! Track transformation of unit vectors.
         XVEC(1)=COST + RVEC(1)*RVEC(1)*(1-COST)
         XVEC(2)=       RVEC(2)*RVEC(1)*(1-COST)+RVEC(3)*SINT
         XVEC(3)=       RVEC(3)*RVEC(1)*(1-COST)-RVEC(2)*SINT
   
         YVEC(1)=       RVEC(1)*RVEC(2)*(1-COST)-RVEC(3)*SINT
         YVEC(2)=COST + RVEC(2)*RVEC(2)*(1-COST)
         YVEC(3)=       RVEC(3)*RVEC(2)*(1-COST)+RVEC(1)*SINT
   
         ZVEC(1)=       RVEC(1)*RVEC(3)*(1-COST)+RVEC(2)*SINT
         ZVEC(2)=       RVEC(2)*RVEC(3)*(1-COST)-RVEC(1)*SINT
         ZVEC(3)=COST + RVEC(3)*RVEC(3)*(1-COST) 

10       CONTINUE    
         !  Find the atom with the largest distance from the z axis.
         DMAX=-1.0D0
         DO J1=1,NATOMS
            DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
            IF (DIST(J1).GT.DMAX) DMAX=DIST(J1) 
         ENDDO
         DMAX2=-1.0D100     

         !  PROJ is the sum of the x components. TVEC arguments are dummys here.
         !  Use T2 as a dummy in order not to change T1 until we have decided which 
         !  atom to put in the xz plane.  
         TXVEC(1:3)=0.0D0; TYVEC(1:3)=0.0D0;TZVEC(1:3)=0.0D0
         DO J1=1,NATOMS
            IF (ABS(DIST(J1)-DMAX).LT.CUTOFF1) THEN
               T2(1:3*NATOMS)=T1(1:3*NATOMS)
               CALL ROTXZ(NATOMS,J1,T2,PROJ,DIST,TXVEC,TYVEC,TZVEC,CUTOFF1)
               IF (ABS(PROJ-DMAX2).LT.CUTOFF1) THEN
                  NORBIT2=NORBIT2+1
                  IF (NORBIT2.EQ.NCHOOSE2) THEN
                     JMAX2=J1
                  ENDIF
               ELSE IF (PROJ.GT.DMAX2) THEN
                  NORBIT2=1
                  DMAX2=PROJ
                  JMAX2=J1
               ENDIF
            ENDIF
         ENDDO
         ! Rotate it into the xz plane
         ! The correct distance from the z axis is used for each atom - they are saved in the vector DIST
         ! DMAX2 is a dummy here.
         CALL ROTXZ(NATOMS,JMAX2,T1,DMAX2,DIST,XVEC,YVEC,ZVEC,CUTOFF1)
         
         ROT(1:3,1)=XVEC(1:3); ROT(1:3,2)=YVEC(1:3); ROT(1:3,3)=ZVEC(1:3)
         ROTINV(1,1:3)=XVEC(1:3); ROTINV(2,1:3)=YVEC(1:3); ROTINV(3,1:3)=ZVEC(1:3) ! transpose
         
      END SUBROUTINE MYORIENT

      SUBROUTINE ROTXZ(NATOMS,JDO,T1,PROJ,DIST,XVEC,YVEC,ZVEC,CUTOFF1)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         INTEGER, INTENT(IN) :: JDO
         REAL(KIND=REAL64), INTENT(INOUT) :: T1(3*NATOMS)
         REAL(KIND=REAL64), INTENT(INOUT) :: PROJ, DIST(NATOMS)
         REAL(KIND=REAL64), INTENT(INOUT) :: XVEC(3), YVEC(3), ZVEC(3)
         REAL(KIND=REAL64), INTENT(IN) :: CUTOFF1

         INTEGER :: J1, J2, J3
         REAL(KIND=REAL64) ::  COST, SINT
         REAL(KIND=REAL64) :: TX, tY, TZ, RDOTN, RVEC(3)

         IF (ABS(T1(3*(JDO-1)+2)).LT.1.0D-8) THEN
            IF (T1(3*(JDO-1)+1).LT.0.0D0) THEN ! rotate about the z axis DO NOT INVERT!!
               DO J1=1,NATOMS
                  T1(3*(J1-1)+1)=-T1(3*(J1-1)+1)
                  T1(3*(J1-1)+2)=-T1(3*(J1-1)+2)
               ENDDO
               XVEC(1)=-XVEC(1); XVEC(2)=-XVEC(2)
               YVEC(1)=-YVEC(1); YVEC(2)=-YVEC(2)
               ZVEC(1)=-ZVEC(1); ZVEC(2)=-ZVEC(2)
            ENDIF
         ELSE   
            COST=T1(3*(JDO-1)+1)/DIST(JDO)
            SINT=T1(3*(JDO-1)+2)/DIST(JDO)         
            IF (ABS(COST**2+SINT**2-1.0D0).GT.2.0D-3) THEN
               WRITE(*, '(A,G20.10)') 'ERROR - in ROTXZ cos**2+sin**2=',COST**2+SINT**2
               STOP
            ENDIF

            DO J1=1,NATOMS
               J2=3*(J1-1)
               TX=T1(J2+1)*COST + T1(J2+2)*SINT
               TY=T1(J2+2)*COST - T1(J2+1)*SINT
               T1(J2+1)=TX
               T1(J2+2)=TY
            ENDDO

            RDOTN=XVEC(3)
            TX=XVEC(1)*COST+XVEC(2)*SINT
            TY=XVEC(2)*COST-XVEC(1)*SINT
            XVEC(1)=TX; XVEC(2)=TY
            RDOTN=YVEC(3)
            TX=YVEC(1)*COST + YVEC(2)*SINT
            TY=YVEC(2)*COST - YVEC(1)*SINT
            YVEC(1)=TX; YVEC(2)=TY
            RDOTN=ZVEC(3)
            TX=ZVEC(1)*COST + ZVEC(2)*SINT
            TY=ZVEC(2)*COST - ZVEC(1)*SINT
            ZVEC(1)=TX; ZVEC(2)=TY
            ! For possibly sloppy alignments with LOCALPERMDIST it may be important to
            ! increase the cutoff and use CUTOFF1 ! DJW 27/8/09
         ENDIF 

         PROJ=0.0D0 
         DO J1=1,NATOMS
            IF (T1(3*(J1-1)+3).GT.CUTOFF1) PROJ=PROJ+T1(3*(J1-1)+1)
         ENDDO       
      END SUBROUTINE ROTXZ


      SUBROUTINE NEWMINDIST2(RA,RB,NATOMS,DIST,DEBUG,RMAT,CMXA,CMYA,CMZA,CMXB,CMYB,CMZB,DWORST)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: RA(3*NATOMS), RB(3*NATOMS)
         REAL(KIND = REAL64), INTENT(OUT) :: DIST, DWORST
         REAL(KIND = REAL64), INTENT(OUT) :: RMAT(3,3)
         REAL(KIND = REAL64), INTENT(OUT) :: CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
         LOGICAL, INTENT(IN) :: DEBUG

         REAL(KIND = REAL64) :: XA(3*NATOMS), XB(3*NATOMS) ! copy of coordinates to be manipulated
         INTEGER :: J1, JMIN, J2 ! counters
         REAL(KIND = REAL64) :: QMAT(4,4), Q1, Q2, Q3, Q4, XM, YM, ZM, XP, YP, ZP ! quaternion
         REAL(KIND = REAL64) :: DIAG(4), TEMPA(9*NATOMS) ! arrays needed for diagonalisation
         INTEGER :: INFO ! exit status for dysev
         REAL(KIND = REAL64) :: MINV
         REAL(KIND = REAL64) :: DISTANCE, SCDUM

         !copy coordinates so we can manipulate them
         XA(1:3*NATOMS)=RA(1:3*NATOMS)
         XB(1:3*NATOMS)=RB(1:3*NATOMS)

         ! Move centre of coordinates of XA and XB to the origin
         CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
         DO J1=1,NATOMS
            CMXA=CMXA+XA(3*(J1-1)+1)
            CMYA=CMYA+XA(3*(J1-1)+2)
            CMZA=CMZA+XA(3*(J1-1)+3)
         ENDDO
         CMXA=CMXA/NATOMS; CMYA=CMYA/NATOMS; CMZA=CMZA/NATOMS
         DO J1=1,NATOMS
            XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-CMXA
            XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-CMYA
            XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-CMZA
         ENDDO
         CMXB=0.0D0; CMYB=0.0D0; CMZB=0.0D0
         DO J1=1,NATOMS
            CMXB=CMXB+XB(3*(J1-1)+1)
            CMYB=CMYB+XB(3*(J1-1)+2)
            CMZB=CMZB+XB(3*(J1-1)+3)
         ENDDO
         CMXB=CMXB/NATOMS; CMYB=CMYB/NATOMS; CMZB=CMZB/NATOMS
         DO J1=1,NATOMS
            XB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB
            XB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB
            XB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB
         ENDDO

         !  The formula below is not invariant to overall translation because XP, YP, ZP
         !  involve a sum of coordinates! We need to have XA and XB coordinate centres both 
         !  at the origin!!
         QMAT(1:4,1:4)=0.0D0
         DO J1=1,NATOMS
            XM=XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
            YM=XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
            ZM=XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
            XP=XA(3*(J1-1)+1)+XB(3*(J1-1)+1)
            YP=XA(3*(J1-1)+2)+XB(3*(J1-1)+2)
            ZP=XA(3*(J1-1)+3)+XB(3*(J1-1)+3)
         !  PRINT '(A,I8,6G18.8)','J1,XM,YM,ZM,XP,YP,ZP=',J1,XM,YM,ZM,XP,YP,ZP
            QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
            QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
            QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
            QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
            QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
            QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
            QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
            QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
            QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
            QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
         ENDDO
         QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3); QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)
         CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)
         IF (INFO.NE.0) PRINT '(A,I6,A)',' newmindist2> WARNING - INFO=',INFO,' in DSYEV'

         MINV=1.0D200
         JMIN=1
         DO J1=1,4
            IF (DIAG(J1).LT.MINV) THEN
            JMIN=J1
            MINV=DIAG(J1)
            ENDIF
         ENDDO
         IF (MINV.LT.0.0D0) THEN
            IF (ABS(MINV).LT.1.0D-6) THEN
            MINV=0.0D0
            ELSE
            PRINT '(A,G20.10,A)',' newmindist2> WARNING MINV is ',MINV,' change to absolute value'
            MINV=-MINV
            ENDIF
         ENDIF
         DIST=SQRT(MINV) 

         !  IF (DEBUG) PRINT '(A,G20.10,A,I6)',' newmindist2> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
         Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
         ! RMAT will contain the matrix that maps XB onto the best correspondence with XA
         RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
         RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
         RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
         RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
         RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
         RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
         RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
         RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
         RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2

         ! Test the rotation, check the distance, and find the biggest atomic displacement
         CALL NEWROTGEOM(NATOMS,XB,RMAT,0.0D0,0.0D0,0.0D0)
         DWORST=0.0D0
         DISTANCE=0.0D0
         DO J2=1,NATOMS
            SCDUM=(XA(3*(J2-1)+1)-XB(3*(J2-1)+1))**2+(XA(3*(J2-1)+2)-XB(3*(J2-1)+2))**2+(XA(3*(J2-1)+3)-XB(3*(J2-1)+3))**2
            DISTANCE=DISTANCE+SCDUM
            IF (SCDUM.GT.DWORST) DWORST=SCDUM
         ENDDO
         DISTANCE=SQRT(DISTANCE)
         DWORST=SQRT(DWORST)
      END SUBROUTINE NEWMINDIST2

      SUBROUTINE MINPERM(N, P, Q, SX, SY, SZ, PBC, PERM, DIST, WORSTDIST, WORSTRADIUS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N ! system size
         REAL(KIND = REAL64), INTENT(IN) :: P(3*N), Q(3*N) ! coordinates
         REAL(KIND = REAL64), INTENT(IN) :: SX, SY, SZ ! box length
         REAL(KIND = REAL64), INTENT(OUT) :: WORSTDIST, WORSTRADIUS
         LOGICAL, INTENT(IN) :: PBC ! periodic boundary conditions?
         INTEGER, INTENT(OUT) :: PERM(N) ! Permutation so that p(i) <--> q(perm(i))
         REAL(KIND = REAL64), INTENT(OUT) :: DIST ! Minimum attainable distance

         REAL(KIND = REAL64) :: S(3), DUMMY
         DOUBLE PRECISION, PARAMETER :: SCALE = 1.0D6 ! Precision
         INTEGER, PARAMETER :: MAXNEI = 60 ! Maximum number of closest neighbours
         !     Internal variables
         !     first:
         !       Sparse matrix of distances
         !     first(i):
         !       Beginning of row i in data,index vectors
         !     kk(first(i)..first(i+1)-1):
         !       Column indexes of existing elements in row i
         !     cc(first(i)..first(i+1)-1):
         !       Matrix elements of row i
         ! These integers might need to be be INT64 - but we try INT32 for now
         INTEGER :: FIRST(N+1), X(N), Y(N)
         INTEGER :: U(N), V(N), H
         INTEGER :: M, I, J, K, L, L2, T, A, I3, J3
         INTEGER :: N8, SZ8, D
         INTEGER :: NDONE, J1, J2

         ! allocate KK and CC
         IF (.NOT. ALLOCATED(KK)) ALLOCATE(KK(N*MAXNEI))
         IF (SIZE(KK) .NE. N*MAXNEI) THEN
             DEALLOCATE(KK)
             ALLOCATE(KK(N*MAXNEI))
         END IF
         IF (.NOT. ALLOCATED(CC)) ALLOCATE(CC(N*MAXNEI))
         IF (SIZE(CC) .NE. N*MAXNEI) THEN
             DEALLOCATE(CC)
             ALLOCATE(CC(N*MAXNEI))
         END IF

         S(1)=SX; S(2)=SY; S(3)=SZ
         M = MIN(MAXNEI,N)
         SZ8 = M*N
         N8 = N

         DO I=0,N
            FIRST(I+1) = I*M + 1
         ENDDO
   
         IF (M.EQ.N) THEN
            ! compute the full matrix
            DO I=1,N
               K = FIRST(I) - 1
               DO J=1,N
                  CC(K+J) = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC) * SCALE
                  KK(K+J) = J
               END DO
            END DO
         ELSE
            ! We need to store the distances of the maxnei closeest neighbors
            ! of each particle. The following builds a heap to keep track of
            ! the maxnei closest neighbours seen so far. It might be more
            ! efficient to use quick-select instead... (This is definately
            ! true in the limit of infinite systems.)
            DO I=1,N
               K=FIRST(I)-1
               DO J=1,M
                  CC(K+J) = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC) * SCALE
                  KK(K+J) = J
                  L = J
10                IF (L.LE.1) GOTO 11
                  L2 = L/2
                  IF (CC(K+L2).LT.CC(K+L)) THEN
                     H = CC(K+L2)
                     CC(K+L2) = CC(K+L)
                     CC(K+L) = H
                     T = KK(K+L2)
                     KK(K+L2) = KK(K+L)
                     KK(K+L) = T
                     L = L2
                     GOTO 10
                  ENDIF
11             ENDDO

               DO J=M+1,N
                  D = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC) * SCALE
                  IF (D.LT.CC(K+1)) THEN
                     CC(K+1) = D
                     KK(K+1) = J
                     L = 1
20                   L2 = 2*L
                     IF (L2+1.GT.M) GOTO 21
                     IF (CC(K+L2+1).GT.CC(K+L2)) THEN
                        A = K+L2+1
                     ELSE
                        A = K+L2
                     END IF
                     IF (CC(A).GT.CC(K+L)) THEN
                        H = CC(A)
                        CC(A) = CC(K+L)
                        CC(K+L) = H
                        T = KK(A)
                        KK(A) = KK(K+L)
                        KK(K+L) = T
                        L = A - K
                        GOTO 20
                     ENDIF
21                   IF (L2.LE.M) THEN ! split IF statements to avoid a segmentation fault
                        IF (CC(K+L2).GT.CC(K+L)) THEN
                           H = CC(K+L2)
                           CC(K+L2) = CC(K+L)
                           CC(K+L) = H
                           T = KK(K+L2)
                           KK(K+L2) = KK(K+L)
                           KK(K+L) = T
                        END IF
                     END IF
                  END IF
               END DO
            END DO
         ENDIF

         ! Call bipartite matching routine
         CALL JOVOSAP(N8, SZ8, CC, KK, FIRST, X, Y, U, V, H)
         ! If initial guess correct, deduce solution distance, which is not done in jovosap
         IF (H.LT.0) THEN
            H = 0
            DO I=1,n
               J = FIRST(I)
30             IF (J.GT.N*MAXNEI) THEN
                  DO J1=1,N
                     PERM(J1) = J1
                  END DO
                  RETURN
               END IF
               IF (KK(J).NE.X(I)) THEN
                  J = J + 1
                  GOTO 30
               END IF
               H = H + CC(J)
            END DO
         END IF

         DO I=1,N
            PERM(I) = X(I)
            IF (PERM(I).GT.N) PERM(I) = N
            IF (PERM(I).LT.1) PERM(I) = 1
         ENDDO

         DIST = DBLE(H)/SCALE
   
         WORSTDIST=-1.0D0
         DO I=1,N
           DUMMY=(p(3*(i-1)+1)-q(3*(perm(i)-1)+1))**2+(p(3*(i-1)+2)-q(3*(perm(i)-1)+2))**2+(p(3*(i-1)+3)-q(3*(perm(i)-1)+3))**2
            IF (DUMMY.GT.WORSTDIST) THEN
               WORSTDIST=DUMMY 
               WORSTRADIUS=p(3*(i-1)+1)**2+p(3*(i-1)+2)**2+p(3*(i-1)+3)**2
            ENDIF
         ENDDO
         WORSTDIST=SQRT(WORSTDIST)
         WORSTRADIUS=MAX(SQRT(WORSTRADIUS),1.0D0)
      END SUBROUTINE MINPERM

      SUBROUTINE NEWROTGEOM(NATOMS,COORDS,MYROTMAT,CX,CY,CZ)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND=REAL64), INTENT(INOUT) :: COORDS(3*NATOMS)
         REAL(KIND=REAL64), INTENT(IN) :: MYROTMAT(3,3), CX, CY, CZ

         INTEGER :: I, J, K
         REAL(KIND=REAL64) :: R1, R0(3)
            
         DO I=1,NATOMS
            R0(1)=COORDS(3*(I-1)+1)-CX
            R0(2)=COORDS(3*(I-1)+2)-CY
            R0(3)=COORDS(3*(I-1)+3)-CZ
            DO J=1,3
               R1=0.0D0
               DO K=1,3
                  R1=R1+MYROTMAT(J,K)*R0(K)
               ENDDO
               IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
               IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
               IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
            ENDDO
         ENDDO
      END SUBROUTINE NEWROTGEOM

      !     permdist is the distance or weight function. It is coded
!     separately for clarity. Just hope that the compiler
!     knows how to to do proper inlining!
!     Input
!       p,q: Coordinates
!       s  : Boxlengths (or dummy if open B.C.)
!       pbc: Periodic boundary conditions?

      PURE REAL(KIND=REAL64) FUNCTION PERMDIST(P, Q, S, PBC)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: P(3), Q(3) ! coordinates
         REAL(KIND=REAL64), INTENT(IN) :: S(3) ! box coordinates
         LOGICAL, INTENT(IN) :: PBC

         REAL(KIND=REAL64) :: T, D
         INTEGER :: I

         D = 0.0d0
         IF (PBC) THEN
            DO I = 1, 3
               IF (S(I).NE.0.0D0) THEN
                  T = Q(I) - P(I)
                  T = T - S(I)*ANINT(T/S(I))
                  D = D + T*T
               ENDIF
            END DO
         ELSE
            D = (Q(1) - P(1))**2+(Q(2) - P(2))**2+(Q(3) - P(3))**2
         END IF
         PERMDIST = D
      END FUNCTION PERMDIST

!     The following routine performs weighted bipartite matching for a sparse non-negative integer weight matrix.
!     "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems," Computing 38, 325-340, 1987
!     by R. Jonker and A. Volgenant, University of Amsterdam.

      SUBROUTINE JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N ! number of rows and columns
         INTEGER, INTENT(IN) :: SZ
         INTEGER, INTENT(INOUT) :: CC(SZ), KK(SZ), FIRST(N+1)
         INTEGER, INTENT(OUT) :: X(N)   ! X = COL ASSIGNED TO ROW
         INTEGER, INTENT(OUT) :: Y(N)   ! Y = ROW ASSIGNED TO COL
         INTEGER, INTENT(OUT) :: U(N)   ! U = DUAL ROW VARIABLE
         INTEGER, INTENT(OUT) :: V(N)   ! V = DUAL COLUMN VARIABLE
         INTEGER, INTENT(OUT) :: H   ! H = VALUE OF OPTIMAL SOLUTION

         INTEGER, PARAMETER :: BIGINT = HUGE(1)
         INTEGER :: J, J0, J1, I, I0, K, L, L0, T, T0, CNT
         INTEGER :: TODO(N), FREE(N), MIN, D(N)
         INTEGER :: V0, VJ
         LOGICAL :: OK(N)


         INTEGER*8 TD,DJ
         INTEGER*8 LAB(N)

         ! INITIALIZATION
         Y(1:N) = 0
         X(1:N) = 0
         TODO(1:N)=0
         H = -1
         J1 = 0
         V(1:N) = BIGINT
         DO I=1,N
            X(I) = 0
            DO T=FIRST(I),FIRST(I+1)-1
               J = KK(T)
               IF (CC(T).LT.V(J)) THEN
                  V(J) = CC(T)
                  Y(J) = I
               END IF
            ENDDO
         ENDDO

         DO J=1,N
            J0 = N-J+1
            I = Y(J0)
            IF (I.EQ.0) RETURN
            IF (X(I).NE.0) THEN
               X(I) = -ABS(X(I))
               Y(J0) = 0
            ELSE
               X(I) = J0
            END IF
         ENDDO

         L=0
         DO I=1,N
            IF (X(I).EQ.0) THEN
               L=L+1
               FREE(L) = I
            ELSE IF (X(I).LT.0) THEN
               X(I) = -X(I)
            ELSE
               J1 = X(I)
               MIN = BIGINT
               DO T=FIRST(I),FIRST(I+1)-1
                  J = KK(T)
                  IF (J.EQ.J1) CYCLE
                  IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
               ENDDO
               V(J1)=V(J1)-MIN
            END IF
         END DO

         ! Improve initial solution
         CNT = 0
         IF (L.EQ.0) RETURN
         
         DO WHILE ((L.GT.0).AND.(CNT.LT.2)) 
            L0 = L
            K = 1
            L = 0
            DO WHILE (K.LE.L0) 
               I = FREE(K)
               K = K+1
               V0 = BIGINT
               VJ = BIGINT
               DO T=FIRST(I),FIRST(I+1)-1
                  J=KK(T)
                  H=CC(T)-V(J)
                  IF (H.LT.VJ) THEN
                     IF (H.GE.V0) THEN
                        VJ=H
                        J1=J
                     ELSE
                        VJ=V0
                        V0=H
                        J1=J0
                        J0=J
                     END IF
                  END IF
               END DO

               I0=Y(J0)
               IF (V0.LT.VJ) THEN
                  V(J0) = V(J0) - VJ+V0               
               ELSE IF (I0.NE.0) THEN
                  J0=J1
                  I0=Y(J1)
               END IF
               IF (I0.NE.0) THEN
                  IF (V0.LT.VJ) THEN
                     K=K-1
                     FREE(K)=I0
                  ELSE
                     L=L+1
                     FREE(L)=I0
                  END IF
               END IF                 
               X(I)=J0
               Y(J0)=I
            END DO
            CNT=CNT+1
         END DO 

         ! AUGMENTATION PART
         L0=L
         DO L=1,L0
            OK(1:N) = .FALSE.
            D(1:N) = BIGINT
            MIN=BIGINT
            I0=FREE(L)
            TD=N
            DO T=FIRST(I0),FIRST(I0+1)-1
               J=KK(T)
               DJ=CC(T)-V(J)
               D(J)=DJ
               LAB(J)=I0
               IF (DJ.LE.MIN) THEN
                 IF (DJ.LT.MIN) THEN
                   MIN=DJ
                   K=1
                   TODO(1)=J
                 ELSE
                   K=K+1
                   TODO(K)=J
                 END IF
               END IF
            END DO
            DO H=1,K
               J=TODO(H)
               IF (J.EQ.0) RETURN
               IF (Y(J).EQ.0) GOTO 80
               OK(J)=.TRUE.
            END DO
            ! REPEAT UNTIL A FREE ROW HAS BEEN FOUND
60          IF (K.EQ.0) RETURN
            J0=TODO(K)
            K=K-1
            I=Y(J0)
            TODO(TD)=J0
            TD=TD-1
            T0=FIRST(I)
            T=T0-1
            T=T+1  !61
            DO WHILE (KK(T).NE.J0) 
               T=T+1
            END DO
            H=CC(T)-V(J0)-MIN
            DO T=T0,FIRST(I+1)-1
               J=KK(T)
               IF (.NOT. OK(J)) THEN
                  VJ=CC(T)-H-V(J)
                  IF (VJ.LT.D(J)) THEN
                     D(J)=VJ
                     LAB(J)=I
                     IF (VJ.EQ.MIN) THEN
                        IF (Y(J).EQ.0) GOTO 70
                        K=K+1
                        TODO(K)=J
                        OK(J)=.TRUE.
                     END IF
                  END IF
               END IF
            END DO
            IF (K.NE.0) GOTO 60
            MIN=BIGINT-1
            DO J=1,N
               IF (D(J).LE.MIN) THEN
                  IF (.NOT. OK(J)) THEN
                     IF (D(J).LT.MIN) THEN
                        MIN=D(J)
                        K=1
                        TODO(1)=J
                     ELSE
                        K=K+1
                        TODO(K)=J
                     END IF
                  END IF
               END IF
            END DO
            DO J0=1,K
               J=TODO(J0)
               IF (Y(J).EQ.0) GOTO 70
               OK(J)=.TRUE.
            END DO
            GOTO 60
70          IF (MIN.EQ.0) GOTO 80
            DO K=TD+1,N
               J0=TODO(K)
               V(J0)=V(J0)+D(J0)-MIN
            END DO
80          I=LAB(J)
            Y(J)=I
            K=J
            J=X(I)
            X(I)=K
            IF (I0.NE.I) GOTO 80
         END DO

         H=0
         DO I=1,N
            J=X(I)
            T=FIRST(I)
            IF (T.GT.SZ) THEN
               PRINT '(A,I6,A)','minperm> WARNING D - atom ',I,' not matched - maximum number of neighbours too small?'
               RETURN
            ENDIF
            DO WHILE (KK(T).NE.J)
               T=T+1
               IF (T.GT.SZ) THEN
                  PRINT '(A,I6,A)','minperm> WARNING D - atom ',I,' not matched - maximum number of neighbours too small?'
                  RETURN
               ENDIF 
            END DO 
            DJ=CC(T)
            U(I)=DJ-V(J)
            H=H+DJ
         END DO
   
      END SUBROUTINE JOVOSAP

END MODULE QCIPERMDIST 
MODULE QCIPERMDIST
   USE QCIPREC
   IMPLICIT NONE
   ! number of permutational groups
   INTEGER :: NPERMGROUP = 0
   ! back up for NPERMGROUP
   INTEGER :: NPERMGROUPBACK = 0
   ! maximum number of subsets (in OPTIM this is set to five as default but can be changed by keyword input)
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

   !variables in keys
   INTEGER :: LOCALPERMNEIGH = 11
   INTEGER :: LOCALPERMMAXSEP = 3

   LOGICAL :: PERMDIST = .FALSE.
   LOGICAL :: PERMDISTINIT = .FALSE.
   LOGICAL :: LOCALPERMDIST = .FALSE.
   LOGICAL :: LPERMDIST = .FALSE.
   LOGICAL :: PERMGUESS = .FALSE.
   LOGICAL :: LPERMOFF=.FALSE.

   LOGICAL :: PBETTER = .FALSE.
    
   REAL(KIND=REAL64) :: LOCALPERMCUT = 0.5D0
   REAL(KIND=REAL64) :: LOCALPERMCUT2 = 5.0D0
   ! is this ever used or initialised?
   REAL(KIND=REAL64) :: LOCALPERMCUTINC

   REAL(KIND=REAL64) :: PDISTANCE
   REAL(KIND=REAL64) :: NOPDISTANCE
    
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
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
            CLOSE(PERMUNIT)
         END IF
      END SUBROUTINE INIT_PERMALLOW

      !
      SUBROUTINE LOPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,DOGROUP,NMOVE,NEWPERM)
         USE QCIKEYS
         USE QCIPREC
         IMPLICIT NONE
         ! input and output variables
         INTEGER, INTENT(IN) :: NATOMS ! number of atoms
         REAL(KIND = REAL64), INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS) ! coordinates for structure A and B
         LOGICAL, INTENT(IN) :: DEBUG ! debugging switch
         REAL(KIND = REAL64), INTENT(IN) :: BOXLX, BOXLY, BOXLZ ! box coordinates
         LOGICAL, INTENT(IN) :: BULKT ! switch for bulk systems
         LOGICAL, INTENT(IN) :: TWOD ! switch for two dimenaional systems
         REAL(KIND = REAL64), INTENT(OUT) :: DISTANCE ! distance between A and B
         REAL(KIND = REAL64), INTENT(OUT) :: DIST2 ! distance squared between A and B
         LOGICAL, INTENT(IN) :: RIGID ! rigid body switch
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
            CALL NEWMINDIST(COORDSB,DUMMYA,NATOMS,NOPDISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT,DWORST)
            ! Test distance for COORDSA with permutation applied in NEWPERM
            DO J1=1,NATOMS
               DUMMYA(3*(J1-1)+1)=COORDSA(3*(NEWPERM(J1)-1)+1)
               DUMMYA(3*(J1-1)+2)=COORDSA(3*(NEWPERM(J1)-1)+2)
               DUMMYA(3*(J1-1)+3)=COORDSA(3*(NEWPERM(J1)-1)+3)
            ENDDO
            CALL NEWMINDIST(COORDSB,DUMMYA,NATOMS,PDISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT,DWORST)
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
         CALL NEWMINDIST(DUMMYB,XBEST,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT,DWORST)
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
         USE QCIKEYS
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
         INTEGER :: I, J1, J1MAX, J2, J2MAX ! counters
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

         ! We can only rotate around the z axis with a pulling potential!
         IF (PULLT) THEN
            T1(1:3*NATOMS)=Q1(1:3*NATOMS)
            GOTO 10
         ENDIF

         DMAX=-1.0D0
         NORBIT1=1

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


      SUBROUTINE NEWMINDIST2()

      END SUBROUTINE NEWMINDIST2

      SUBROUTINE MINPERM()

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


END MODULE QCIPERMDIST 
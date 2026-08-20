module lpermdist
   use qciprec
   implicit none
   contains
      !>  This routine uses local optimal alignment for each group of permutable atoms.
      !!  It is intended for use with CHARMM and AMBER.
      !!  Overall alignment is based on the transformation for the best preserved local group.
      !!
      !!  COORDSA becomes the optimal alignment of the optimal permutation
      !!  isomer with the coordinate permutations applied. DISTANCE is the residual square distance
      !!  for the best alignment with respect to permutation as well as
      !!  orientation and centre of mass.
      !!
      !!  The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
      !!  centre of coordinates of COORDSA will be the same as for COORDSB.
      !!
      subroutine lopermdist(coordsb,coordsa,distance,dist2,rmatbest,nmove,newperm, active_perm_only)
         use qcikeys, only: natoms !, debug
         use minperm_mod
         use perm_defs
         
         implicit none
         
         real(kind=real64), intent(in) :: coordsb(3*natoms) !< coordinates of one point preserved
         real(kind=real64), intent(inout) :: coordsa(3*natoms) !< coordinates of second point that will change
         real(kind=real64), intent(out) :: distance, dist2
         real(kind=real64), intent(out) :: rmatbest(3,3)
         integer, intent(out) :: nmove
         integer, intent(out) :: newperm(natoms)
         logical, intent(in) :: active_perm_only

         !arrays to track permutations
         integer :: lperm(natoms), lpermbest(natoms), lpermbestatom(natoms), updateperm(natoms)
         real(kind=real64) :: pdummya(3*natoms), pdummyb(3*natoms), spdummya(3*natoms), spdummyb(3*natoms), xbest(3*natoms)

         integer :: j1, j2, j3, j4, ndummy, ndummy2, idx1, idx2, nats, jmax1a, jmax2a, jmax1b, jmax2b
         !> local copies of the coordinates
         real(kind=real64) :: dummya(3*natoms), dummyb(3*natoms)
         !> sum of distances
         real(kind=real64) :: dsum
         !> number of atoms in current group
         integer :: patoms
         !> best distances (squared) for all groups and local one for finding it
         real(kind=real64) :: ldbest(npermgroup), ldbestatom, ldistance, worstrad, dworst
         !> array to track atoms we have tried so far
         integer :: tried(natoms)
         !> number of distances in list
         integer :: ndmean
         !> list of distances
         real(kind=real64) ::dmean(natoms)
         !> reference of atoms in distance list
         integer :: sortlist(natoms)
         !> centre for each group in a and b
         real(kind=real64) :: xa, ya, za, xb, yb, zb
         !> distances between atoms
         real(kind=real64) :: da, db, thisd
         !> logicals used to control list addition etc
         logical :: done, more2addt
         !> information about atoms in obligatory groups
         integer :: dlist(natoms), nother, nadded
         !> rotation matrices for standard alignments
         real(kind=real64) :: rota(3,3), rotinva(3,3), rotb(3,3), rotinvb(3,3)

         !  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
         !  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case,
         !  and is evaluated as a sum of local distances squared for permutable groups.

         !  We iterate over rounds of permutational/orientational alignment until we have converged to the identity permutation.
         !  The maximum number of pair exchanges associated with a group is two.

         ! Most arrays refer to the original atom order. For A, which is changed, we need to use NEWPERM to access current atom ids.

         ! copy coordinate data to local arrays
         dummya(1:3*natoms) = coordsa(1:3*natoms)
         dummyb(1:3*natoms) = coordsb(1:3*natoms)  

         !initialise newperm and set localpermneigh
         do j1=1,natoms
            newperm(j1)=j1
            lperm(j1)=-1
         enddo
         dsum=0.0d0
         localpermneigh=min(localpermneigh,natoms)

         ldbest(1:npermgroup) = 1.0d100

         !start iteration over all groups - we use dummy counter to access arrays with varying group sizes
         ndummy = 1
         do j1=1,npermgroup

            !Apply permutations to active groups only
            if (active_perm_only) then
               if (.not.groupactive(j1)) cycle
            end if
            

            !copy group size to local variable - we need it a lot
            patoms = npermsize(j1)
            !allocate some local variables, including lperm and lpermbest
            !initialise/reset variables for this group
            ! TRIED(j) is 0 if atom j is eligible to be a neighbour, but has not yet been tried.
            ! It is -1 if it is ineligible, or has been tried and broke the alignment. 
            ! It is +1 if it has been tried and did not break the alignment. 
            ! It is also -1 for atoms already in the set of permutable atoms in question. 
            lpermbest(1:natoms) = -1 !is this the correct initialisation?
            tried(1:natoms) = 0
            do j2=1,patoms
               lpermbest(j2)=j2
            end do
            !get the coordinates for atoms in the groups and the geometric centre
            xa = 0.0d0; ya=0.0d0; za = 0.0d0
            xb = 0.0d0; yb=0.0d0; zb = 0.0d0
            do j2=1,patoms
               tried(permgroup(ndummy+j2-1))=-1
               ! for coords from A, we need to consider potential applied permutations!
               idx1 = 3*(newperm(permgroup(ndummy+j2-1))-1)
               idx2 = 3*(j2-1)
               pdummya(idx2+1) = dummya(idx1+1)
               pdummya(idx2+2) = dummya(idx1+2)
               pdummya(idx2+3) = dummya(idx1+3)
               xa = xa + pdummya(idx2+1)
               ya = ya + pdummya(idx2+2)
               za = za + pdummya(idx2+3)
               ! for coords from B, we can just use the allocation we know from the perm groups
               idx1 = 3*(permgroup(ndummy+j2-1)-1)
               idx2 = 3*(j2-1)
               pdummyb(idx2+1) = dummyb(idx1+1)
               pdummyb(idx2+2) = dummyb(idx1+2)
               pdummyb(idx2+3) = dummyb(idx1+3)   
               xb = xb + pdummyb(idx2+1)
               yb = yb + pdummyb(idx2+2)
               zb = zb + pdummyb(idx2+3)            
            end do
            xa = xa/patoms; ya = ya/patoms; za = za/patoms
            xb = xb/patoms; yb = yb/patoms; zb = zb/patoms
            spdummya(1:3*patoms) = pdummya(1:3*patoms)
            spdummyb(1:3*patoms) = pdummyb(1:3*patoms)

            !initialise distance list 
            dmean(1:natoms) = 10.0*localpermcut2
            ndmean = 0
            sortlist(1:natoms) = -1
            !create a list of sorted distances to group
            do j2=1,natoms
               !compute distances
               thisd = 1.0d9 ! default have a large distance
               ! do not use members of the group as reference neighbours
               if (tried(j2).eq.-1) then
                  cycle
               else
                  !compute the actual distance
                  idx1 = 3*(newperm(j2)-1)
                  da = (xa - dummya(idx1+1))**2 + (ya - dummya(idx1+2))**2 + (za - dummya(idx1+3))**2
                  idx1 = 3*(j2-1)
                  db = (xb - dummyb(idx1+1))**2 + (yb - dummyb(idx1+2))**2 + (zb - dummyb(idx1+3))**2
                  thisd = 0.5d0*(sqrt(da)+sqrt(db))
                  if (thisd.gt.localpermcut2) cycle
               end if
               ndmean = ndmean + 1
               !now add it to the sorted list
               done = .false.
               do j3=1,ndmean
                  if (thisd.lt.dmean(j3)) then
                     !position found for the new atom
                     do j4=ndmean,j3+1,-1
                        dmean(j4) = dmean(j4-1)
                        sortlist(j4) = sortlist(j4-1)
                     end do
                     !add this distance
                     dmean(j3) = thisd
                     sortlist(j3) = j2
                     done = .true.
                  end if
                  if (done) exit
               end do
            end do ! end of creation of distance list

            ! If the atoms in the permuting group are associated with obligatory permutations 
            ! then add all the members in these SETS. NDUMMY is already set for this group J1.
            ! This setup means that such atoms can appear in the list of candidates sorted by distance above.
            ! They will not be used with a tried value of 1.
            ! We still include everything in the cutoff, but it would be more efficient not to include these
            ! atoms in the sorted candidates list.

            nother = 0
            if (nsets(j1).gt.0) then
               !there are multiple sets in group j1
               do j2 = 1,nsets(j1)
                  do j3=ndummy,ndummy+npermsize(j1)-1
                     !here we use the fact that exchangable groups must have the same size
                     nother = nother + 1
                     dlist(nother) = sets(permgroup(j3),j1,j2)
                     tried(dlist(nother)) = +1
                  end do
               end do
            end if

            !We are now at the end of the rather lengthy setup for the alignment
            ! In the following, we add neighbours one at a time in order of 
            ! increasing distance from primary permutable set and test whether they break the alignment.
            more2addt = .true.
            !this replaces the original goto statements
            do while (more2addt)
               !restore the coordinates from saved coords

               pdummya(1:3*patoms)=spdummya(1:3*patoms)
               pdummyb(1:3*patoms)=spdummyb(1:3*patoms)   
               !intialise used variables
               ldbestatom = 1.0d100
               nother = 0 
               !reset/update nother and dlist
               do j2 = 1,natoms
                  if (tried(j2).eq.1) then
                     nother = nother + 1
                     dlist(nother) = j2
                  end if
               end do
               !iterate over all atoms until we either find a new atom or run out of atoms
               do j2=1,natoms
                  !we do not have any atom left below the threshold
                  if (dmean(j2).gt.localpermcut2) then
                     more2addt = .false.
                     exit
                  end if
                  !atom has not been tried and is not in the original group
                  if (tried(sortlist(j2)).eq.0) then
                     nother = nother + 1
                     !sanity check
                     if (nother+patoms.gt.natoms) then
                        !call log_message(3, "Number of neighbours plus permutable atoms exceeds natoms")
                        write(*,*) "Number of neighbours plus permutable atoms exceeds natoms"
                        stop
                     end if
                     dlist(nother) = sortlist(j2)
                     exit
                  end if
               end do
               if (.not.more2addt) exit 

               ! we have found an atom to add, but this could also be permutable
               ! next test whether we can identify any additional permutable atoms 
               nadded = 1
               idx1 = dlist(nother)
               if (permutable(idx1)) then
                  !set up variable ndummy2 as reference to start of group for atom
                  ndummy2 = 1
                  !use iteration to increment ndummy2 by the sizes of all groups up to the relevant one
                  do j2=1,ingroup(idx1)-1
                     ndummy2 = ndummy2 + npermsize(j2)
                  end do
                  !now check the other atoms in the permutable group
                  do j2=1,npermsize(ingroup(idx1))
                     idx2 = ndummy2+j2-1
                     !skip idx1
                     if (permgroup(idx2).eq.idx1) cycle
                     !otherwise add the atom
                     if (tried(permgroup(idx2)).eq.0) then
                        nother = nother + 1
                        nadded = nadded + 1
                        !sanity check
                        if (nother+patoms.gt.natoms) then
                           !call log_message(3, " Number of neighbours plus permutable atoms exceeds natoms")
                           write(*,*) " Number of neighbours plus permutable atoms exceeds natoms"
                           stop
                        end if
                        dlist(nother) = permgroup(idx2)
                     else
                        !call log_message(3, " Partner atom has alreasy been tried")
                        write(*,*) " Partner atom has alreasy been tried"
                        stop                  
                     end if
                  end do
               end if !end test whether new atom is part of permutable group

               !we know have new atom(s) added and need to add them to pdummya and pdummyb
               do j2=1,nother
                  idx1 = 3*(patoms+j2-1)
                  idx2 = 3*(newperm(dlist(j2))-1)
                  pdummya(idx1+1) = dummya(idx2+1)
                  pdummya(idx1+2) = dummya(idx2+2)                 
                  pdummya(idx1+3) = dummya(idx2+3)
                  idx2 = 3*(dlist(j2)-1)    
                  pdummyb(idx1+1) = dummyb(idx2+1)
                  pdummyb(idx1+2) = dummyb(idx2+2)                 
                  pdummyb(idx1+3) = dummyb(idx2+3)                                
               end do
               !save coordinates as we will change pdummya and pdummyb
               spdummya(3*patoms+1:3*(patoms+nother))=pdummya(3*patoms+1:3*(patoms+nother))
               spdummyb(3*patoms+1:3*(patoms+nother))=pdummyb(3*patoms+1:3*(patoms+nother))
            
               ! we simplify the code here as we only work with biomolecules
               ! we want a standard orientation for the groups in A and B
               ! this orientation is given by:
               ! 1. move the CoX of the atoms in the group + neighbourhood to the origin
               ! 2. identify the most distant atom and put it on the z axis
               ! 3. identify the next most distant atom and put it in the xy plane
               ! the below routine is a simplified version of myorient in OPTIM
               nats = nother+patoms
               call standard_orient(nats,pdummya(1:3*nats),jmax1a,jmax2a,rota,rotinva)
               call standard_orient(nats,pdummyb(1:3*nats),jmax1b,jmax2b,rotb,rotinvb)
               write(*,*) " lpermdist> group ", j1, " for A: atom with z alignment is ", jmax1a, " and for xy alignment is ", jmax2a
               write(*,*) " lpermdist> group ", j1, " for B: atom with z alignment is ", jmax1b, " and for xy alignment is ", jmax2b 

               ! Optimimise permutational isomer for the standard orientation for the
               ! current choice of atoms from the possible orbits.
               ! MINPERM does not change PDUMMYB and PDUMMYA.
               call minperm(nats,pdummyb,pdummya, lperm, ldistance, dist2, worstrad)
               !check we are not mixing permutational group with outside atoms
               do j2=1,patoms
                  if (lperm(j2).gt.patoms) then
                     ldistance = 1.0d300
                     exit
                  end if
               end do

               !check we are not changing permutations of non-permutable atoms
               do j2=1,nother
                  if (lperm(patoms+j2).ne.patoms+j2) then
                     if ((.not.(permutable(dlist(j2)))).and.(.not.(permutable2(dlist(j2))))) then 
                        ldistance = 1.0d300
                        exit
                     end if
                  end if
               end do

               ! save the best permutation and local distances for this set of atoms
               if (ldistance.lt.ldbestatom) then
                  ldbestatom = ldistance
                  lpermbestatom(1:patoms) = lperm(1:patoms)
               end if

               ! test whether alignment is broken, if it is update tried and continue looping
               if (sqrt(ldbestatom).gt.localpermcut) then
                  tried(dlist(nother)) = -1
                  if (nadded.gt.1) then
                     tried(dlist(nother-nadded+1:nother-1))=-1
                  end if
               ! if the alignment is not broken update tried and only continue looping if we haven't found enough neighbours
               else
                  tried(dlist(nother)) = 1
                  if (nadded.gt.1) then
                     tried(dlist(nother-nadded+1:nother-1))=1
                  end if
                  ldbest(j1) = ldbestatom
                  lpermbest(1:patoms) = lpermbestatom(1:patoms)
                  if (nother.ge.localpermneigh) then
                     more2addt = .false.
                  end if    
               end if
            end do ! end of inner loop with more2addt

            ! at this point we should have the best permutation for group j1
            lperm(1:patoms) = lpermbest(1:patoms)

            !reset the updated array for perms
            updateperm(1:natoms) = newperm(1:natoms)

            !apply permutations to group
            do j2=1,patoms
               updateperm(permgroup(ndummy+j2-1)) = newperm(permgroup(ndummy+lperm(j2)-1))
            end do
            
            !update permutations of any associated atoms
            if (nsets(j1).gt.0) then
               do j2=1,patoms
                  do j3=1,nsets(j1)
                     updateperm(sets(permgroup(ndummy+j2-1),j1,j3)) = sets( newperm(permgroup(ndummy+lperm(j2)-1)) , j1, j3)
                  end do
               end do
            end if
            !update newperm
            newperm(1:natoms) = updateperm(1:natoms)
            ! update distance
            dsum = dsum + sqrt(ldbest(j1))
            !update ndummy to account for the size of the group
            ndummy = ndummy + npermsize(j1)
         end do ! end of loop to iterate over all perm groups

         !update the coordinates for dummya with the permutation
         !dummya is unchanged from coordsa, so we can use it as reference
         nmove = 0
         do j1=1,natoms
            if (newperm(j1).ne.j1) then
               idx1 = 3*(j1-1)
               idx2 = 3*(newperm(j1)-1)
               dummya(idx1+1) = coordsa(idx2+1)
               dummya(idx1+2) = coordsa(idx2+2)
               dummya(idx1+3) = coordsa(idx2+3)
               nmove = nmove + 1
               
            end if
         end do

         call mindist(dummyb,dummya,rmatbest,distance,dworst,xbest)
         dist2 = distance**2
         coordsa(1:3*natoms) = xbest(1:3*natoms)
         write(*,*) " lopermdist> Distance after alignment is ", distance
      end subroutine lopermdist
   
      subroutine standard_orient(nats,coords,jmax1,jmax2,rot,rotinv)
         
         implicit none
         
         integer, intent(in) :: nats !patoms + nother in lpermdist
         real(kind=real64), intent(inout) :: coords(3*nats) !coordinates
         integer, intent(out) :: jmax1, jmax2 
         real(kind=real64), intent(out) :: rot(3,3), rotinv(3,3)
         real(kind=real64) :: xvec(3), yvec(3), zvec(3), cox(3)
         integer :: j1, j2, idx
         real(kind=real64) :: dmax, dist(nats)
         real(kind=real64) :: newx(3*nats), finalx(3*nats)
         
         ! setup unit vectors for x, y and z
         xvec = (/ 1.0,0.0,0.0 /)
         yvec = (/ 0.0,1.0,0.0 /)
         zvec = (/ 0.0,0.0,1.0 /)

         ! determine centre of coordinates
         cox(1:3) = 0.0d0
         do j1=1,nats
            idx = 3*(j1-1)
            do j2=1,3
               cox(j2) = cox(j2) + coords(idx+j2)
            end do
         end do
         cox(1:3) = cox(1:3)/nats

         ! move centre of coordinates to origin (centre of mass with unit masses)
         do j1=1,nats
            idx = 3*(j1-1)
            do j2=1,3
               coords(idx+j2) = coords(idx+j2) - cox(j2)
            end do
         end do

         dmax = -1.0d0
         !orbital analysis based on mean recipricol pair distances
         ! we are making a key assumption here:
         ! in myorient, we test whether multiple atoms are at the same distance within a certian cutoff tolerance
         ! if there is symmetry (think LJ38 or similar clusters) this leads to multiple alignements to be tested
         ! for biomolecules, we assume that symmetry is broken, and hence we only need a single pass to find the best alignment
         ! i.e. norbit1 = 1 in all cases
         do j1 = 1,nats
            idx = 3*(j1-1)
            dist(j1) = sqrt(coords(idx+1)**2 + coords(idx+2)**2 + coords(idx+3)**2)
            if (dist(j1).gt.dmax) then
               dmax = dist(j1)
               jmax1 = j1
               !norbit1 = 1
            end if
         end do

         !rotate the atoms to the align jamx1 with the z axis
         call rotate_to_z_axis(nats,coords,newx,jmax1,dmax,xvec,yvec,zvec)

         ! now find the next atom to align to the xy plane -> again we simplifiy compared to myorient assuming there is a unique furthest atom
         ! note that the distances here are the xy distances only
         dmax = -1.0d0
         do j1=1, nats !natoms
            idx = 3*(j1-1)
            dist(j1) = sqrt(newx(idx+1)**2 + newx(idx+2)**2)
            if (dist(j1).gt.dmax) then
               dmax = dist(j1)
               jmax2 = j1
            end if
         end do

         !apply the actual rotation
         call rotxz(nats,newx,finalx,jmax2,dmax,xvec,yvec,zvec)

         coords(1:3*nats) = finalx(1:3*nats)

         !compute the actual rotation matrices
         rot(1:3,1) = xvec(1:3); rot(1:3,2) = yvec(1:3); rot(1:3,3) = zvec(1:3)
         rotinv(1,1:3) = xvec(1:3); rotinv(2,1:3) = yvec(1:3); rotinv(3,1:3) = zvec(1:3)
      end subroutine standard_orient

       SUBROUTINE CHECK_PERM_BAND( REVERSET)
         USE QCIKEYS, ONLY: DEBUG, NATOMS, NIMAGES, QCIPERMCUT
         USE MOD_INTCOORDS, ONLY: XYZ
         
         IMPLICIT NONE

         LOGICAL, INTENT(IN) :: REVERSET  !<stepping direction through images

         INTEGER :: J1, J2, IDX, FIRSTIMAGE, SECONDIMAGE
         INTEGER :: STARTIDX, ENDIDX, STEP
         REAL(KIND = REAL64) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
         INTEGER :: PERMP(NATOMS), NMOVEP
         REAL(KIND = REAL64) :: RMATBEST(3,3)
         REAL(KIND=REAL64) :: DISTANCE, DIST2

         !REAL(KIND = REAL64) :: SAVECOORDSA(3*NATOMS)
        
         REAL(KIND=REAL64) :: SAVELOCALPERMCUT
         !SAVELOCALPERMCUT=LOCALPERMCUT
         !LOCALPERMCUT=QCIPERMCUT !This was coppied from OPTIM 


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

         DO J1=STARTIDX,ENDIDX,STEP !cycling through images
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
               FIRSTIMAGE = J1     !1
            ELSE
               FIRSTIMAGE = J1 + 2 !NIMAGES+2
            END IF 
            SECONDIMAGE = J1+1     !2 or NIMAGES+1
            
            !coordinates for image 1
            COORDSB(1:3*NATOMS) = XYZ((3*NATOMS*(FIRSTIMAGE-1)+1):3*NATOMS*FIRSTIMAGE)
            !coordinates for image 2
            COORDSA(1:3*NATOMS) = XYZ((3*NATOMS*(SECONDIMAGE-1)+1):3*NATOMS*SECONDIMAGE)            
            
            CALL LOPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,RMATBEST,NMOVEP,PERMP, .TRUE.)      
            

            !IF (DEBUG.AND.NMOVEP.GT.0)  THEN
               WRITE(*,*) ' check_perm_band> alignment of images ',SECONDIMAGE,FIRSTIMAGE,' moves=',&
                          NMOVEP, ' permutations, distance =',DISTANCE 
            !END IF
               
            XYZ((3*NATOMS*(SECONDIMAGE-1)+1):3*NATOMS*SECONDIMAGE) = COORDSA(1:3*NATOMS)
            END DO
         !LOCALPERMCUT=SAVELOCALPERMCUT
         !WRITE(*,*) "Completed perm band check"
      END SUBROUTINE CHECK_PERM_BAND

      SUBROUTINE UPDATE_ACTIVE_PERMGROUPS()
         USE INTERPOLATION_KEYS, ONLY: ATOMACTIVE
         use perm_defs, only: groupactive, npermsize, npermgroup, permgroup
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

      subroutine rotate_to_z_axis(nats,x,newx,jmax,dmax,xvec,yvec,zvec)
         implicit none
         integer, intent(in) :: nats
         real(kind=real64), intent(in) :: x(3*nats)
         real(kind=real64), intent(out) :: newx(3*nats)
         integer, intent(in) :: jmax
         real(kind=real64), intent(in) :: dmax
         real(kind=real64), intent(inout) :: xvec(3), yvec(3), zvec(3)

         real(kind=real64) :: rvec(3), cost, sint, proj, rdotn
         integer :: idx, idx2, j1
         real(kind=real64), parameter :: thresh = 1.0d-8

         newx(1:3*nats) = 0.0d0
         idx = 3*(jmax-1)

         !check we don't have an unusual set of coordinates that are basically aligned
         if ((abs(x(idx+1)).lt.thresh).and.(abs(x(idx+2)).lt.thresh)) then
            !if z component is greater than zero we are already in the right orientation
            if (x(idx+3).gt.0.0d0) then
               newx(1:3*nats) = x(1:3*nats)
            ! otherwise rotate around the x axis (make sure not to invert by accident)
            else
               do j1=1,nats
                  idx2 = 3*(j1-1)
                  newx(idx2+1) = x(idx2+1)
                  newx(idx2+2) =-x(idx2+2)
                  newx(idx2+3) =-x(idx2+3)
               end do
               yvec(2) = -1.0d0
               zvec(3) = -1.0d0
            end if
            return !leave early as we got alignment
         end if

         !cos and sine based on x,y and z coordinates for jmax
         cost = x(idx+3)/dmax
         sint = sqrt(x(idx+1)**2 + x(idx+2)**2)/dmax

         !rotate jmax onto the z axis using rotation formula from Goldstein p. 165
         rvec(1) = x(idx+2)/sqrt(x(idx+1)**2 + x(idx+2)**2)
         rvec(2) =-x(idx+1)/sqrt(x(idx+1)**2 + x(idx+2)**2)
         rvec(3) = 0.0d0 

         !apply rotation
         do j1=1,nats
            idx = 3*(j1-1)
            rdotn = x(idx+1)*rvec(1) + x(idx+2)*rvec(2) + x(idx+3)*rvec(3)
            proj = rdotn * (1.0d0 - cost)
            newx(idx+1) = x(idx+1)*cost + rvec(1)*proj - (x(idx+2)*rvec(3) - x(idx+3)*rvec(2))*sint
            newx(idx+2) = x(idx+2)*cost + rvec(2)*proj - (x(idx+3)*rvec(1) - x(idx+1)*rvec(3))*sint
            newx(idx+3) = x(idx+3)*cost + rvec(3)*proj - (x(idx+1)*rvec(2) - x(idx+2)*rvec(1))*sint
         end do

         ! track unit vector transformation
         xvec(1) = cost + rvec(1)*rvec(1)*(1 - cost)
         xvec(2) =        rvec(2)*rvec(1)*(1 - cost) + rvec(3)*sint
         xvec(3) =        rvec(3)*rvec(1)*(1 - cost) - rvec(2)*sint

         yvec(1) =        rvec(1)*rvec(2)*(1 - cost) - rvec(3)*sint
         yvec(2) = cost + rvec(2)*rvec(2)*(1 - cost)
         yvec(3) =        rvec(3)*rvec(2)*(1 - cost) + rvec(1)*sint

         zvec(1) =        rvec(1)*rvec(3)*(1 - cost) + rvec(2)*sint
         zvec(2) =        rvec(2)*rvec(3)*(1 - cost) - rvec(1)*sint
         zvec(3) = cost + rvec(3)*rvec(3)*(1 - cost)
         
      end subroutine rotate_to_z_axis

      subroutine rotxz(nats,x,newx,jmax,dmax,xvec,yvec,zvec)
         implicit none
         integer, intent(in) :: nats
         real(kind=real64), intent(in) :: x(3*nats)
         real(kind=real64), intent(out) :: newx(3*nats)
         integer, intent(in) :: jmax
         real(kind=real64), intent(in) :: dmax
         real(kind=real64), intent(inout) :: xvec(3), yvec(3), zvec(3)

         real(kind=real64) :: cost, sint, tx, ty, rdotn
         integer :: idx, idx2, j1
         real(kind=real64), parameter :: thresh = 1.0d-8

         newx(1:3*nats) = 0.0d0
         idx = 3*(jmax-1)

         !check we are not already aligned
         if (abs(x(idx+2)).lt.thresh) then
            !already got the right alignment
            if (x(idx+1).gt.0) then
               newx(1:3*nats) = x(1:3*nats)
            !otherwise rotate about the z axis (avoid inversion!)
            else
               do j1=1,nats
                  idx2 = 3*(j1-1)
                  newx(idx2+1) =-x(idx2+1)
                  newx(idx2+2) =-x(idx2+2)
                  newx(idx2+3) = x(idx2+3)
               end do
               xvec(1) = -xvec(1); xvec(2) = -xvec(2)
               yvec(1) = -yvec(1); yvec(2) = -yvec(2)
               zvec(1) = -zvec(1); zvec(2) = -zvec(2)
            end if
            return
         end if

         ! otherwise apply the correct rotation in 2D
         cost = x(idx+1)/dmax
         sint = x(idx+2)/dmax

         do j1=1,nats
            idx2 = 3*(j1-1)
            newx(idx2+1) = x(idx2+1)*cost + x(idx2+2)*sint
            newx(idx2+2) = x(idx2+2)*cost - x(idx2+1)*sint
            newx(idx2+3) = x(idx2+3)
         end do

         !update unit vectors
         tx = xvec(1)*cost + xvec(2)*sint
         ty = xvec(2)*cost - xvec(1)*sint
         xvec(1) = tx
         xvec(2) = ty
         tx = yvec(1)*cost + yvec(2)*sint
         ty = yvec(2)*cost - yvec(1)*sint
         yvec(1) = tx
         yvec(2) = ty
         tx = zvec(1)*cost + zvec(2)*sint
         ty = zvec(2)*cost - zvec(1)*sint
         zvec(1) = tx
         zvec(2) = ty
      end subroutine rotxz

      subroutine mindist(xa,xb,rmat,distance,dworst,newxb)
         use qcikeys, only: natoms, debug
         implicit none
         real(kind=real64), intent(in) :: xa(3*natoms)
         real(kind=real64), intent(in) :: xb(3*natoms)
         real(kind=real64), intent(out) :: rmat(3,3)
         real(kind=real64), intent(out) :: distance
         real(kind=real64), intent(out) :: dworst
         real(kind=real64), intent(out) :: newxb(3*natoms)
         real(kind=real64) :: distbefore, thisdist
         integer :: idx, j
         real(kind=real64) :: cmxa(3), cmxb(3), cxa(3*natoms), cxb(3*natoms), rotxb(3*natoms)

         ! we need centre coordinates
         cmxa(1:3) = 0.0d0
         cmxb(1:3) = 0.0d0

         do idx=1,natoms
            j=3*(idx-1)
            cmxa(1:3) = cmxa(1:3) + xa(j+1:j+3)
            cmxb(1:3) = cmxb(1:3) + xb(j+1:j+3)
         end do

         cmxa(1:3) = cmxa(1:3)/natoms
         cmxb(1:3) = cmxb(1:3)/natoms

         do idx=1,natoms
            j=3*(idx-1)
            cxa(j+1:j+3) = xa(j+1:J+3) - cmxa(1:3)
            cxb(j+1:j+3) = xb(J+1:j+3) - cmxb(1:3)
         end do

         call get_rot_matrix(cxa,cxb,rmat,distbefore)
         call apply_rot(cxb,rmat,rotxb) 

         do idx=1,natoms
            j=3*(idx-1)
            newxb(j+1:j+3) = rotxb(j+1:j+3) + cmxa(1:3)
         end do

         dworst=0.0d0
         distance=0.0d0
         do idx=1,natoms
            j = 3*(idx-1)
            thisdist = (xa(j+1)-newxb(j+1))**2 + (xa(j+2)-newxb(j+2))**2 + (xa(j+3)-newxb(j+3))**2
            distance = distance + thisdist
            if (thisdist.gt.dworst) dworst=thisdist
         end do
         distance = sqrt(distance)
         dworst = sqrt(dworst)
         if (debug) then
            write(*,*) " mindist> Distance before rotation: ", distbefore, " | Distance after rotation: ", distance
            write(*,*) "          Largest distance for single atom: ", dworst
         end if
      end subroutine

      !assumes centred coordinates
      subroutine apply_rot(x,rmat,newx)
         use qcikeys, only: natoms
         implicit none
         real(kind=real64), intent(in) :: x(3*natoms)
         real(kind=real64), intent(out) :: newx(3*natoms)        
         real(kind=real64), intent(in) :: rmat(3,3)
         real(kind=real64) :: r0(3), r1
         integer :: i, j, k

         newx(1:3*natoms) = 0.0d0

         do i=1,natoms
            r0(1:3) = x(3*i-2:3*i)
            do j=1,3
               r1 = 0.0d0
               do k=1,3
                  r1 = r1 + rmat(j,k)*r0(k)
               end do
               newx(3*(i-1)+j) = r1
            end do
         end do
      end subroutine apply_rot 

      !> Find rotational matrix between xa and xb using quaternions
      !> assumes xa and xb are both centred at the origin
      !> New analytic method based on quaterions from Kearsley, Acta Cryst. A, 45, 208-210, 1989.
      subroutine get_rot_matrix(xa,xb,rmat,distance)
         use qcikeys, only: natoms
         !use logger, only: log_message
         implicit none
         real(kind=real64), intent(in) :: xa(3*natoms), xb(3*natoms)
         real(kind=real64), intent(out) :: rmat(3,3)
         real(kind=real64), intent(out) :: distance
         real(kind=real64) :: qmat(4,4)
         real(kind=real64) :: xm, ym, zm, xp, yp, zp
         integer :: idx, j, info
         real(kind=real64) :: tempa(9*natoms), diag(4)
         integer :: jmin
         real(kind=real64) :: minv, q1, q2, q3, q4
         real(kind=real64), parameter :: zerothresh = 1.0d-6

         qmat(1:4,1:4) = 0.0d0
         do idx = 1,natoms
            j=3*(idx-1)
            xm = xa(j+1) - xb(j+1)
            ym = xa(j+2) - xb(j+2)
            zm = xa(j+3) - xb(j+3)
            xp = xa(j+1) + xb(j+1)
            yp = xa(j+2) + xb(j+2)
            zp = xa(j+3) + xb(j+3)
            qmat(1,1) = qmat(1,1) + xm**2 + ym**2 + zm**2
            qmat(1,2) = qmat(1,2) + yp*zm - ym*zp
            qmat(1,3) = qmat(1,3) + xm*zp - xp*zm
            qmat(1,4) = qmat(1,4) + xp*ym - xm*yp
            qmat(2,2) = qmat(2,2) + yp**2 + zp**2 + xm**2
            qmat(2,3) = qmat(2,3) + xm*ym - xp*yp
            qmat(2,4) = qmat(2,4) + xm*zm - xp*zp
            qmat(3,3) = qmat(3,3) + xp**2 + zp**2 + ym**2
            qmat(3,4) = qmat(3,4) + ym*zm - yp*zp
            qmat(4,4) = qmat(4,4) + xp**2 + yp**2 + zm**2
         end do
         qmat(2,1)=qmat(1,2)
         qmat(3,1)=qmat(1,3)
         qmat(3,2)=qmat(2,3)
         qmat(4,1)=qmat(1,4)
         qmat(4,2)=qmat(2,4)
         qmat(4,3)=qmat(3,4)

         call dsyev('V','U',4,qmat,4,diag,tempa,9*natoms,info)
         if (info.ne.0) then
            !call log_message(2, "get_rot_matrix> non-zero exit status in DSYEV")
            write(*,*) "get_rot_matrix> non-zero exit status in DSYEV"
         end if

         ! find current distance and min value
         minv=1.0d200
         jmin=1
         do idx=1,4
            if (diag(idx).lt.minv) then
               jmin = idx
               minv = diag(idx)
            end if
         end do
         if (minv.lt.0.0d0) then
            if (abs(minv).lt.zerothresh) then
               minv = 0.0d0
            else
               minv = -minv
               !call log_message(2, " get_rot_matrix> minv is negative, change to absolute value")
               write(*,*) " get_rot_matrix> minv is negative, change to absolute value"
            end if
         end if

         distance = sqrt(minv)

         q1=qmat(1,jmin)
         q2=qmat(2,jmin)
         q3=qmat(3,jmin)
         q4=qmat(4,jmin)

         rmat(1,1) = q1**2 + q2**2 - q3**2 - q4**2
         rmat(1,2) = 2*(q2*q3 + q1*q4)
         rmat(1,3) = 2*(q2*q4 - q1*q3)
         rmat(2,1) = 2*(q2*q3 - q1*q4)
         rmat(2,2) = q1**2 + q3**2 - q2**2 - q4**2
         rmat(2,3) = 2*(q3*q4 + q1*q2)
         rmat(3,1) = 2*(q2*q4 + q1*q3)
         rmat(3,2) = 2*(q3*q4 - q1*q2)
         rmat(3,3) = q1**2 + q4**2 - q2**2 - q3**2
      end subroutine get_rot_matrix

end module lpermdist
